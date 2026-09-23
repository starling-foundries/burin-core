//! The polygon rule under a profile (nucleus rule).
//!
//! The cells of a polygon at resolution R are exactly the depth-R cells whose ellipsoidal
//! **nucleus** (the cell's reference point: the inverse projection of its planar centre) the
//! polygon contains, boundary excluded. The coverer descends from the B bases; prunes a node
//! whose conservative lon/lat bounding box misses the polygon; emits every descendant of a node
//! whose bounding box the polygon contains; otherwise recurses; at depth R applies the nucleus
//! test. Bounding boxes are conservative and can only fail to prune, so predicate ties at the
//! bbox stage cannot change the cell set; only the nucleus test decides.

use crate::error::{invalid, Result};
use crate::geo::Grid;
use crate::hash::{Ctx, Digest};
use crate::hierarchy::{Cid, Hierarchy};
use crate::profile::Profile;
use crate::tree::Tree;
use geo::{coord, Intersects, Line, MultiPolygon, Point, Polygon, Rect};
use rstar::{RTree, RTreeObject, AABB};
use serde_json::Value;
use std::sync::Arc;

pub const MAX_RESOLUTION: u32 = 15;

struct GridInfo {
    grid: Grid,
    polar_cap_lat: f64,
    /// Nucleus longitude of the equatorial bases O, P, Q, R.
    face_centre_lon: [f64; 4],
    /// Planar centre of the north and south polar squares (the poles), metres.
    pole: [(f64, f64); 2],
}

impl GridInfo {
    fn new(profile: &Profile) -> Result<GridInfo> {
        let grid = profile.grid();
        let cap = grid.boundary3(&[0])?.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
        let mut centres = [0.0; 4];
        for (i, c) in centres.iter_mut().enumerate() {
            *c = grid.nucleus(&[i as u32 + 1])?.0;
        }
        let w0 = grid.cell_width(0);
        let pole = [(grid.ul0[0].0 + w0 / 2.0, grid.ul0[0].1 - w0 / 2.0), (grid.ul0[5].0 + w0 / 2.0, grid.ul0[5].1 - w0 / 2.0)];
        Ok(GridInfo { grid, polar_cap_lat: cap, face_centre_lon: centres, pole })
    }

    /// A conservative `(lon_min, lat_min, lon_max, lat_max)` containing the whole cell.
    ///
    /// Equatorial cells: longitude depends only on x and latitude only on y, so the 8 boundary
    /// samples (corners included) bound the cell exactly. Polar cells: latitude is a function of
    /// the L∞ distance from the polar square's centre, so the extreme latitudes are at the
    /// farthest corner and at the nearest point of the planar rectangle; longitude is monotone
    /// along the perimeter as seen from the centre, so its extremes are at corners; the cap cell
    /// (centre inside) spans every longitude. Base cells reach the pole.
    fn cell_bbox(&self, suid: &[u32]) -> Result<(f64, f64, f64, f64)> {
        let cap = self.polar_cap_lat;
        let base = suid[0];
        if suid.len() == 1 {
            return Ok(match base {
                0 => (-180.0, cap, 180.0, 90.0),
                5 => (-180.0, -90.0, 180.0, -cap),
                b => {
                    let c = self.face_centre_lon[b as usize - 1];
                    let (mut lo, mut hi) = (c - 45.0, c + 45.0);
                    if lo < -180.0 {
                        lo += 360.0;
                        hi += 360.0;
                    }
                    (lo, -cap, hi, cap)
                }
            });
        }
        if (1..=4).contains(&base) {
            let pts = self.grid.boundary3(suid)?;
            let (lon_min, lon_max) = lon_range(pts.iter().map(|p| p.0));
            let lat_min = pts.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
            let lat_max = pts.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);
            return Ok((lon_min, lat_min, lon_max, lat_max));
        }
        let north = base == 0;
        let (x0, y0) = self.grid.ul_vertex(suid);
        let w = self.grid.cell_width((suid.len() - 1) as u32);
        let (x1, y1) = (x0 + w, y0 - w);
        let (cx, cy) = self.pole[if north { 0 } else { 1 }];
        if x0 < cx && cx < x1 && y1 < cy && cy < y0 {
            return Ok(if north { (-180.0, cap, 180.0, 90.0) } else { (-180.0, -90.0, 180.0, -cap) });
        }
        let corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)].iter().map(|&(x, y)| self.grid.inverse(x, y, None)).collect::<Result<Vec<_>>>()?;
        let (lon_min, lon_max) = lon_range(corners.iter().map(|p| p.0));
        let lat_far = if north { corners.iter().map(|p| p.1).fold(f64::INFINITY, f64::min) } else { corners.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max) };
        let (nx, ny) = (cx.max(x0).min(x1), cy.max(y1).min(y0));
        let lat_near = self.grid.inverse(nx, ny, None)?.1;
        Ok(if north {
            (lon_min, lat_far, lon_max, if lat_near > 60.0 { 90.0 } else { lat_near })
        } else {
            (lon_min, if lat_near < -60.0 { -90.0 } else { lat_near }, lon_max, lat_far)
        })
    }
}

/// A longitude interval containing every value; `hi > 180` signals an antimeridian wrap.
fn lon_range(lons: impl Iterator<Item = f64> + Clone) -> (f64, f64) {
    let lo = lons.clone().fold(f64::INFINITY, f64::min);
    let hi = lons.clone().fold(f64::NEG_INFINITY, f64::max);
    if hi - lo <= 180.0 {
        return (lo, hi);
    }
    let shifted = lons.map(|x| if x < 0.0 { x + 360.0 } else { x });
    let slo = shifted.clone().fold(f64::INFINITY, f64::min);
    let shi = shifted.fold(f64::NEG_INFINITY, f64::max);
    if shi - slo <= 180.0 {
        return (slo, shi);
    }
    (-180.0, 180.0)
}

fn rect(lon_lo: f64, lat_lo: f64, lon_hi: f64, lat_hi: f64) -> Rect {
    Rect::new(coord! { x: lon_lo, y: lat_lo }, coord! { x: lon_hi, y: lat_hi })
}

/// The bbox as one rectangle, or two when it wraps the antimeridian.
fn bbox_parts((lon_lo, lat_lo, lon_hi, lat_hi): (f64, f64, f64, f64)) -> Vec<Rect> {
    if lon_hi <= 180.0 {
        vec![rect(lon_lo, lat_lo, lon_hi, lat_hi)]
    } else {
        vec![rect(lon_lo, lat_lo, 180.0, lat_hi), rect(-180.0, lat_lo, lon_hi - 360.0, lat_hi)]
    }
}

/// One polygon edge in the R-tree, tagged with the part and ring it belongs to.
struct Edge {
    line: Line,
    part: u32,
    ring: u32, // 0 = exterior, k = k-th hole
}

impl RTreeObject for Edge {
    type Envelope = AABB<[f64; 2]>;
    fn envelope(&self) -> AABB<[f64; 2]> {
        let (a, b) = (self.line.start, self.line.end);
        AABB::from_corners([a.x.min(b.x), a.y.min(b.y)], [a.x.max(b.x), a.y.max(b.y)])
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Pos {
    Inside,
    Boundary,
    Outside,
}

/// The polygon with an R-tree over every edge. Point containment is an exact crossing test
/// (robust orientation, ring orientation irrelevant, holes honoured, boundary reported);
/// the bounding-box predicates the descent needs are built from corner containment and edge
/// queries, so every test is O(log n + k) rather than linear in the polygon's vertices.
struct Indexed {
    tree: RTree<Edge>,
    holes_per_part: Vec<u32>,
}

fn orient(a: geo::Coord, b: geo::Coord, p: geo::Coord) -> f64 {
    robust::orient2d(robust::Coord { x: a.x, y: a.y }, robust::Coord { x: b.x, y: b.y }, robust::Coord { x: p.x, y: p.y })
}

impl Indexed {
    fn new(polygon: &MultiPolygon) -> Indexed {
        let mut edges = Vec::new();
        let mut holes_per_part = Vec::with_capacity(polygon.0.len());
        for (pi, poly) in polygon.0.iter().enumerate() {
            holes_per_part.push(poly.interiors().len() as u32);
            for (ri, ring) in std::iter::once(poly.exterior()).chain(poly.interiors().iter()).enumerate() {
                for line in ring.lines() {
                    if line.start != line.end {
                        edges.push(Edge { line, part: pi as u32, ring: ri as u32 });
                    }
                }
            }
        }
        Indexed { tree: RTree::bulk_load(edges), holes_per_part }
    }

    /// Where `p` lies: a horizontal ray to +x, crossings counted per ring with exact orientation.
    fn position(&self, p: geo::Coord) -> Pos {
        let mut parity: std::collections::HashMap<(u32, u32), bool> = std::collections::HashMap::new();
        let ray = AABB::from_corners([p.x, p.y], [f64::INFINITY, p.y]);
        for e in self.tree.locate_in_envelope_intersecting(&ray) {
            let (a, b) = (e.line.start, e.line.end);
            let o = orient(a, b, p);
            if o == 0.0 && p.x >= a.x.min(b.x) && p.x <= a.x.max(b.x) && p.y >= a.y.min(b.y) && p.y <= a.y.max(b.y) {
                return Pos::Boundary;
            }
            // half-open rule on y so a vertex on the ray is counted exactly once
            let crosses = if a.y <= p.y {
                b.y > p.y && o > 0.0 // upward edge, p strictly left of it
            } else {
                b.y <= p.y && o < 0.0 // downward edge, p strictly right of it
            };
            if crosses {
                let v = parity.entry((e.part, e.ring)).or_insert(false);
                *v = !*v;
            }
        }
        for (part, &n_holes) in self.holes_per_part.iter().enumerate() {
            let part = part as u32;
            if !parity.get(&(part, 0)).copied().unwrap_or(false) {
                continue; // outside this part's exterior
            }
            let in_hole = (1..=n_holes).any(|h| parity.get(&(part, h)).copied().unwrap_or(false));
            if !in_hole {
                return Pos::Inside;
            }
        }
        Pos::Outside
    }

    fn any_edge_touches(&self, rect: &Rect) -> bool {
        let env = AABB::from_corners([rect.min().x, rect.min().y], [rect.max().x, rect.max().y]);
        self.tree.locate_in_envelope_intersecting(&env).any(|e| e.line.intersects(rect))
    }

    fn corners(rect: &Rect) -> [geo::Coord; 4] {
        let (lo, hi) = (rect.min(), rect.max());
        [coord! { x: lo.x, y: lo.y }, coord! { x: hi.x, y: lo.y }, coord! { x: hi.x, y: hi.y }, coord! { x: lo.x, y: hi.y }]
    }

    /// The rectangle meets the polygon (boundary contact counts).
    fn intersects(&self, rect: &Rect) -> bool {
        Indexed::corners(rect).iter().any(|c| self.position(*c) != Pos::Outside) || self.any_edge_touches(rect)
    }

    /// The rectangle lies in the polygon's interior (DE-9IM `T**FF*FF*`): every corner strictly
    /// inside and no polygon edge touching the closed rectangle.
    fn contains_properly(&self, rect: &Rect) -> bool {
        Indexed::corners(rect).iter().all(|c| self.position(*c) == Pos::Inside) && !self.any_edge_touches(rect)
    }

    /// Interior only: a point on the boundary is not contained.
    fn contains_point(&self, p: &Point) -> bool {
        self.position(p.0) == Pos::Inside
    }
}

/// Mixed-level cids: contained subtrees as one node, boundary leaves by the nucleus rule.
pub fn cover_nodes(polygon: &MultiPolygon, resolution: u32, profile: &Profile) -> Result<Vec<Cid>> {
    if resolution > MAX_RESOLUTION {
        return invalid(format!("resolution {resolution} out of range [0, {MAX_RESOLUTION}]"));
    }
    profile.validate()?;
    let h = profile.hierarchy();
    let info = GridInfo::new(profile)?;
    let poly = Indexed::new(polygon);
    let mut out = Vec::new();

    #[allow(clippy::too_many_arguments)]
    fn visit(suid: &mut Vec<u32>, cid: Cid, level: u32, res: u32, h: &Hierarchy, info: &GridInfo, poly: &Indexed, out: &mut Vec<Cid>) -> Result<()> {
        let parts = bbox_parts(info.cell_bbox(suid)?);
        if !parts.iter().any(|p| poly.intersects(p)) {
            return Ok(());
        }
        if level == res {
            let (lon, lat) = info.grid.nucleus(suid)?;
            if poly.contains_point(&Point::new(lon, lat)) {
                out.push(cid);
            }
            return Ok(());
        }
        // interior-only containment: a bbox touching the polygon boundary descends, so a nucleus
        // exactly on the boundary is judged by the leaf rule, never by the shortcut
        if parts.iter().all(|p| poly.contains_properly(p)) {
            out.push(cid);
            return Ok(());
        }
        for k in 0..h.a {
            suid.push(k);
            visit(suid, h.child(cid, k), level + 1, res, h, info, poly, out)?;
            suid.pop();
        }
        Ok(())
    }

    for b in 0..h.b {
        let mut suid = vec![b];
        visit(&mut suid, (h.a + b) as u64, 0, resolution, &h, &info, &poly, &mut out)?;
    }
    Ok(out)
}

/// Sorted depth-`resolution` cids whose nucleus `polygon` contains.
pub fn polyfill(polygon: &MultiPolygon, resolution: u32, profile: &Profile) -> Result<Vec<Cid>> {
    let h = profile.hierarchy();
    let mut out = Vec::new();
    for cid in cover_nodes(polygon, resolution, profile)? {
        out.extend(h.descendants(cid, resolution - h.level(cid)));
    }
    out.sort_unstable();
    out.dedup();
    Ok(out)
}

/// The canonical coverage tree of a polygon at `resolution` (depth = resolution). An empty
/// coverage (a polygon smaller than a cell) is a valid, empty tree.
pub fn cover(polygon: &MultiPolygon, resolution: u32, profile: &Profile, ctx: Arc<Ctx>) -> Result<Tree> {
    let nodes = cover_nodes(polygon, resolution, profile)?;
    Tree::from_cells(profile.hierarchy(), resolution, nodes, ctx)
}

/// 32-byte root of `cover(polygon, resolution, profile)`.
pub fn fingerprint_polygon(polygon: &MultiPolygon, resolution: u32, profile: &Profile) -> Result<Digest> {
    let ctx = Arc::new(profile.ctx(resolution)?);
    Ok(cover(polygon, resolution, profile, ctx)?.root())
}

fn ring_from_json(v: &Value) -> Result<geo::LineString> {
    let pts = v.as_array().ok_or_else(|| crate::Error::Invalid("ring must be an array".into()))?;
    let mut coords = Vec::with_capacity(pts.len());
    for p in pts {
        let xy = p.as_array().filter(|a| a.len() >= 2).ok_or_else(|| crate::Error::Invalid("position must be [lon, lat]".into()))?;
        let (x, y) = (xy[0].as_f64(), xy[1].as_f64());
        match (x, y) {
            (Some(x), Some(y)) => coords.push(coord! { x: x, y: y }),
            _ => return invalid("position must be numeric"),
        }
    }
    Ok(geo::LineString::new(coords))
}

fn polygon_from_json(rings: &Value) -> Result<Polygon> {
    let rs = rings.as_array().ok_or_else(|| crate::Error::Invalid("polygon must be an array of rings".into()))?;
    if rs.is_empty() {
        return invalid("polygon has no rings");
    }
    let exterior = ring_from_json(&rs[0])?;
    let holes = rs[1..].iter().map(ring_from_json).collect::<Result<Vec<_>>>()?;
    Ok(Polygon::new(exterior, holes))
}

/// A GeoJSON `Polygon`, `MultiPolygon`, or a `Feature` carrying one, as a `MultiPolygon`.
pub fn geometry_from_geojson(v: &Value) -> Result<MultiPolygon> {
    let ty = v.get("type").and_then(Value::as_str).unwrap_or("");
    match ty {
        "Feature" => geometry_from_geojson(v.get("geometry").ok_or_else(|| crate::Error::Invalid("feature without geometry".into()))?),
        "Polygon" => Ok(MultiPolygon::new(vec![polygon_from_json(v.get("coordinates").ok_or_else(|| crate::Error::Invalid("missing coordinates".into()))?)?])),
        "MultiPolygon" => {
            let polys = v.get("coordinates").and_then(Value::as_array).ok_or_else(|| crate::Error::Invalid("missing coordinates".into()))?;
            Ok(MultiPolygon::new(polys.iter().map(polygon_from_json).collect::<Result<Vec<_>>>()?))
        }
        other => invalid(format!("unsupported GeoJSON type {other:?}; need Polygon or MultiPolygon")),
    }
}

/// Convenience: fingerprint a GeoJSON geometry.
pub fn fingerprint_geojson(v: &Value, resolution: u32, profile: &Profile) -> Result<Digest> {
    fingerprint_polygon(&geometry_from_geojson(v)?, resolution, profile)
}
