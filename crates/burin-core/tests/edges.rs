//! Edge cases where two code paths or a numeric limit could quietly disagree: point lookup
//! across levels, at seams, corners and poles; lookup against the drawn polygon; the deepest
//! levels; the 64-bit limit; and the tree built from adversarial cell lists.

use burin_core::geo::{Grid, Shape};
use burin_core::hash::Ctx;
use burin_core::hierarchy::{Cid, Hierarchy, SPACE};
use burin_core::profile::Profile;
use burin_core::tree::Tree;
use burin_core::zone::cell_at;
use geo::{Contains, LineString, Point, Polygon};
use std::sync::Arc;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    /// An area-uniform point on the sphere, degrees.
    fn point(&mut self) -> (f64, f64) {
        (self.unit() * 360.0 - 180.0, libm::asin(self.unit() * 2.0 - 1.0).to_degrees())
    }
}

fn profiles() -> [Profile; 3] {
    [Profile::ogc(), Profile::burin_1(), Profile { lon_0_udeg: -123_456_789, north_square: 3, south_square: 1, ..Profile::default() }]
}

fn lookup(grid: &Grid, lon: f64, lat: f64, level: u32) -> Option<Cid> {
    grid.cell_from_lonlat(lon, lat, level).map(|(b, r, c)| cell_at(&SPACE, b, level, r, c).unwrap())
}

/// Points that sit on the hard places: every base-cell corner and edge midpoint and the corners of
/// cells at a few levels (projected back to lon/lat), the seams between base cells, the band edge.
fn awkward_points(grid: &Grid, profile: &Profile) -> Vec<(f64, f64)> {
    let mut pts = Vec::new();
    for level in 1..=4u32 {
        for c in (0..6u64).flat_map(|b| SPACE.descendants(9 + b, level)).step_by(7) {
            let path = SPACE.path(c).unwrap();
            let (x, y) = grid.ul_vertex(&path);
            let w = grid.cell_width(level);
            for (dx, dy) in [(0.0, 0.0), (w / 2.0, 0.0), (0.0, w / 2.0)] {
                if let Ok(p) = grid.inverse(x + dx, y - dy, None) {
                    pts.push(p);
                }
            }
        }
    }
    let lon_0 = profile.lon_0();
    for k in -4..=4 {
        for lat in [0.0, 30.0, -30.0, 41.9375, -41.9375, 41.94, -41.94, 89.9, -89.9] {
            let lon = ((lon_0 + 45.0 * k as f64 + 180.0).rem_euclid(360.0)) - 180.0;
            pts.push((lon, lat));
            pts.push((lon.next_up(), lat));
            pts.push((lon.next_down(), lat));
        }
    }
    pts
}

#[test]
fn a_point_keeps_its_ancestors_at_every_level() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    for profile in profiles() {
        let grid = profile.grid();
        let mut pts: Vec<(f64, f64)> = (0..20_000).map(|_| rng.point()).collect();
        pts.extend(awkward_points(&grid, &profile));
        for (lon, lat) in pts {
            let cells: Vec<Cid> = (0..=15).map(|l| lookup(&grid, lon, lat, l).unwrap()).collect();
            for l in 1..cells.len() {
                assert_eq!(cells[l] / 9, cells[l - 1], "({lon:?}, {lat:?}): level {l} is not inside level {} ({})", l - 1, profile.describe());
            }
        }
    }
}

#[test]
fn longitudes_wrap_poles_are_one_cell_and_bad_latitudes_are_refused() {
    let mut rng = Rng(7);
    for profile in profiles() {
        let grid = profile.grid();
        for _ in 0..2_000 {
            // lon a multiple of 1/64 degree, so lon ± 360 is exact
            let lon = (rng.next() % (360 * 64)) as f64 / 64.0 - 180.0;
            let lat = rng.point().1;
            let here = lookup(&grid, lon, lat, 15);
            assert_eq!(lookup(&grid, lon + 360.0, lat, 15), here, "({lon}, {lat}) + 360");
            assert_eq!(lookup(&grid, lon - 360.0, lat, 15), here, "({lon}, {lat}) - 360");
        }
        for level in [0u32, 1, 7, 15] {
            for (pole, base) in [(90.0, 0u32), (-90.0, 5u32)] {
                let cells: std::collections::BTreeSet<Cid> = (0..72).map(|k| lookup(&grid, k as f64 * 5.0 - 180.0, pole, level).unwrap()).collect();
                assert_eq!(cells.len(), 1, "pole {pole} at level {level} is split ({})", profile.describe());
                let path = SPACE.path(*cells.iter().next().unwrap()).unwrap();
                assert_eq!(path[0], base);
                assert_eq!(grid.shape(&path), Shape::Cap, "the pole lies in the cap cell");
            }
        }
        for (lon, lat) in [(0.0, 90.000_001), (0.0, -90.5), (10.0, 180.0), (f64::NAN, 0.0), (0.0, f64::NAN), (f64::INFINITY, 0.0), (0.0, f64::NEG_INFINITY)] {
            assert_eq!(grid.cell_from_lonlat(lon, lat, 5), None, "({lon}, {lat}) is not a point of the grid");
        }
    }
}

fn signed_area(ring: &[(f64, f64)]) -> f64 {
    ring.windows(2).map(|w| w[0].0 * w[1].1 - w[1].0 * w[0].1).sum::<f64>() / 2.0
}

#[test]
fn the_cell_a_point_is_assigned_contains_it_when_drawn() {
    let mut rng = Rng(0xdead_beef);
    for profile in profiles() {
        let grid = profile.grid();
        for _ in 0..3_000 {
            let (lon, lat) = rng.point();
            let level = 2 + (rng.next() % 9) as u32;
            let cid = lookup(&grid, lon, lat, level).unwrap();
            let path = SPACE.path(cid).unwrap();
            // how far inside its cell the point is, as a fraction of the cell width
            let (x, y) = grid.forward(lon, lat, None);
            let (ulx, uly) = grid.ul_vertex(&path);
            let w = grid.cell_width(level);
            let clearance = ((x - ulx).min(ulx + w - x).min(uly - y).min(y - (uly - w))) / w;
            let equatorial = (1..=4).contains(&path[0]);
            // equatorial edges are exact; polar edges are densified curves, so stay off them
            if clearance < if equatorial { 1e-9 } else { 0.05 } {
                continue;
            }
            let ring = grid.cell_polygon(&path, 9).unwrap();
            assert!(signed_area(&ring) > 0.0);
            let shifted = if ring.iter().any(|p| p.0 > 180.0) && lon < 0.0 { lon + 360.0 } else { lon };
            let poly = Polygon::new(LineString::from(ring), vec![]);
            assert!(poly.contains(&Point::new(shifted, lat)), "({lon}, {lat}) is assigned {path:?} but lies outside its polygon ({})", profile.describe());
        }
    }
}

#[test]
fn the_deepest_levels_round_trip_near_corners_and_caps() {
    for profile in profiles() {
        let grid = profile.grid();
        for level in 9..=15u32 {
            let mut paths = Vec::new();
            for b in 0..6u32 {
                for digit in [0u32, 2, 4, 6, 8] {
                    paths.push(std::iter::once(b).chain(std::iter::repeat_n(digit, level as usize)).collect::<Vec<u32>>());
                }
            }
            for path in paths {
                let cid = SPACE.cid(&path).unwrap();
                let (lon, lat) = grid.nucleus(&path).unwrap();
                assert_eq!(lookup(&grid, lon, lat, level), Some(cid), "{path:?} at level {level} ({})", profile.describe());
                let ring = grid.cell_polygon(&path, 5).unwrap();
                assert!(signed_area(&ring) > 0.0, "{path:?}: not counterclockwise");
                if grid.shape(&path) != Shape::Cap {
                    let lon = if ring.iter().any(|p| p.0 > 180.0) && lon < 0.0 { lon + 360.0 } else { lon };
                    assert!(Polygon::new(LineString::from(ring), vec![]).contains(&Point::new(lon, lat)), "{path:?} at level {level}");
                }
            }
        }
    }
}

#[test]
fn the_sixty_four_bit_limit_is_enforced_exactly() {
    assert_eq!(SPACE.max_level(), 18);
    assert_eq!(Hierarchy::new(2, 1).unwrap().max_level(), 62);
    let deepest: Vec<u32> = std::iter::once(5).chain(std::iter::repeat_n(8, 18)).collect();
    let last = SPACE.cid(&deepest).unwrap();
    assert_eq!(last, 15 * 9u64.pow(18) - 1, "the last cell at the deepest level");
    assert_eq!(SPACE.path(last).unwrap(), deepest);
    assert!(SPACE.cid(&[deepest.clone(), vec![0]].concat()).is_err(), "one level deeper overflows");
    assert!(cell_at(&SPACE, 0, 19, 0, 0).is_err());
    assert_eq!(cell_at(&SPACE, 5, 18, 3u64.pow(18) - 1, 3u64.pow(18) - 1).unwrap(), last);
    let ctx = |d| Arc::new(Ctx::sha256(9, 6, d));
    assert!(Tree::empty(SPACE, 19, ctx(19)).is_err());
    let t = Tree::from_cells(SPACE, 18, [last, SPACE.cid(&[0]).unwrap()], ctx(18)).unwrap();
    assert_eq!(t.cells(), vec![SPACE.cid(&[0]).unwrap(), last]);
    assert_eq!(t.leaf_count(), 9u64.pow(18) + 1);
    let one_by_one = Tree::empty(SPACE, 18, ctx(18)).unwrap().set_full(last).unwrap().set_full(SPACE.cid(&[0]).unwrap()).unwrap();
    assert_eq!(t.root(), one_by_one.root());
}

#[test]
fn the_bulk_builder_is_canonical_on_adversarial_lists() {
    let ctx = |d| Arc::new(Ctx::sha256(9, 6, d));
    let cid = |p: &[u32]| SPACE.cid(p).unwrap();
    let oracle = |d: u32, cells: &[Cid]| cells.iter().fold(Tree::empty(SPACE, d, ctx(d)).unwrap(), |t, &c| t.set_full(c).unwrap());
    let q4_children: Vec<Cid> = (0..9).map(|k| cid(&[3, 4, k])).collect();
    let cases: Vec<(&str, u32, Vec<Cid>, Vec<Cid>)> = vec![
        ("nine children are their parent", 3, q4_children.clone(), vec![cid(&[3, 4])]),
        ("eight children stay eight", 3, q4_children[..8].to_vec(), q4_children[..8].to_vec()),
        ("81 grandchildren are the base cell", 2, (0..81).map(|k| cid(&[3, k / 9, k % 9])).collect(), vec![cid(&[3])]),
        ("a parent swallows its children", 3, [q4_children.clone(), vec![cid(&[3, 4])], q4_children.clone()].concat(), vec![cid(&[3, 4])]),
        ("touching ranges across base cells stay apart", 2, vec![cid(&[3, 8, 8]), cid(&[4, 0, 0])], vec![cid(&[3, 8, 8]), cid(&[4, 0, 0])]),
        ("the planet from its 54 level-1 cells", 1, (0..54).map(|k| cid(&[k / 9, k % 9])).collect(), (0..6).map(|b| cid(&[b])).collect()),
        ("nothing", 4, vec![], vec![]),
    ];
    for (name, d, input, want) in cases {
        let bulk = Tree::from_cells(SPACE, d, input.clone(), ctx(d)).unwrap();
        assert_eq!(bulk.cells(), want, "{name}");
        assert_eq!(bulk.root(), oracle(d, &input).root(), "{name}: root differs from cell-by-cell insertion");
    }
    assert!(Tree::from_cells(SPACE, 1, [cid(&[3, 4, 5])], ctx(1)).is_err(), "a cell below the tree's depth");
    assert!(Tree::from_cells(SPACE, 3, [5u64], ctx(3)).is_err(), "not a cid");
}
