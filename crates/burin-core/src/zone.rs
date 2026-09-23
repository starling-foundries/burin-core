//! Zone topology: where a cell sits inside its base cell, the scanline order of its sub-zones,
//! its four edge neighbours, and whole-grid rasters with halos.
//!
//! ```text
//! position   (base, level, row, col)     row 0 at the top of the base cell, n = N^level
//! digit      d = N·row_i + col_i          one base-N digit of row and of col per level
//! scanline   index = row·N^k + col        a sub-zone at depth k inside its ancestor
//! raster     index = base·n² + row·n + col
//! ```
//! A cid's digits interleave its row and column, so every query here is integer arithmetic on
//! the cid, and none of it changes a root. The six base cells are laid out as rHEALPix lays them:
//! O, P, Q, R in one horizontal band that wraps, N above the band's column `north_square` and S
//! below column `south_square`. Crossing onto N or S away from its own column rotates the polar
//! square by the column offset in quarter turns. At each corner of a base cell only three base
//! cells meet, so the diagonal beyond a corner has no cell.

use crate::error::{invalid, Result};
use crate::hash::Ctx;
use crate::hierarchy::{Cid, Hierarchy};
use crate::profile::Profile;
use crate::tree::{make_branch, Node, Tree};
use std::sync::Arc;

/// The largest raster or halo index this module builds, in cells.
pub const MAX_RASTER_CELLS: u64 = 1 << 30;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Position {
    pub base: u32,
    pub level: u32,
    pub row: u64,
    pub col: u64,
}

/// The edge directions of the planar layout, in the order `neighbours` returns them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    Up,
    Right,
    Down,
    Left,
}

pub const DIRECTIONS: [Direction; 4] = [Direction::Up, Direction::Right, Direction::Down, Direction::Left];

impl Direction {
    pub fn name(self) -> &'static str {
        match self {
            Direction::Up => "up",
            Direction::Right => "right",
            Direction::Down => "down",
            Direction::Left => "left",
        }
    }

    fn step(self) -> (i64, i64) {
        match self {
            Direction::Up => (-1, 0),
            Direction::Right => (0, 1),
            Direction::Down => (1, 0),
            Direction::Left => (0, -1),
        }
    }
}

/// `N` for a square aperture `A = N²`.
pub fn side(h: &Hierarchy) -> Result<u64> {
    let n = (h.a as f64).sqrt().round() as u64;
    if n * n != h.a as u64 {
        return invalid(format!("aperture {} is not a square; cells have no row and column", h.a));
    }
    Ok(n)
}

fn span(h: &Hierarchy, level: u32) -> Result<u64> {
    side(h)?.checked_pow(level).ok_or_else(|| crate::Error::Invalid(format!("level {level} overflows")))
}

/// Where `cid` sits in its base cell.
pub fn position(h: &Hierarchy, cid: Cid) -> Result<Position> {
    let n = side(h)?;
    let path = h.path(cid)?;
    let (mut row, mut col) = (0u64, 0u64);
    for &d in &path[1..] {
        row = row * n + d as u64 / n;
        col = col * n + d as u64 % n;
    }
    Ok(Position { base: path[0], level: (path.len() - 1) as u32, row, col })
}

/// The cid of the cell at `level` containing `(lon, lat)` degrees, by the reference rule
/// (`Grid::cell_from_planar`). Non-finite coordinates, or a point outside the image, are an error.
pub fn cell_from_point(profile: &Profile, lon: f64, lat: f64, level: u32) -> Result<Cid> {
    profile.validate()?;
    match profile.grid().cell_from_lonlat(lon, lat, level) {
        Some((b, r, c)) => cell_at(&profile.hierarchy(), b, level, r, c),
        None => invalid(format!("({lon}, {lat}) is not a point of the grid")),
    }
}

/// The cid at `(row, col)` of base cell `base` at `level`.
pub fn cell_at(h: &Hierarchy, base: u32, level: u32, row: u64, col: u64) -> Result<Cid> {
    let n = side(h)?;
    let width = span(h, level)?;
    if base >= h.b {
        return invalid(format!("base {base} out of range [0, {})", h.b));
    }
    if row >= width || col >= width {
        return invalid(format!("({row}, {col}) outside a {width}x{width} base cell at level {level}"));
    }
    let mut c = (h.a + base) as u64;
    for j in (0..level).rev() {
        let p = n.pow(j);
        c = c * h.a as u64 + (row / p % n) * n + col / p % n;
    }
    Ok(c)
}

/// The scanline index of `cid` among the sub-zones of its ancestor `ancestor`.
pub fn subzone_index(h: &Hierarchy, ancestor: Cid, cid: Cid) -> Result<u64> {
    h.check(ancestor, None)?;
    h.check(cid, None)?;
    if !h.is_ancestor(ancestor, cid) {
        return invalid(format!("{ancestor} is not an ancestor of {cid}"));
    }
    let (pa, pc) = (position(h, ancestor)?, position(h, cid)?);
    let width = span(h, pc.level - pa.level)?;
    Ok((pc.row - pa.row * width) * width + (pc.col - pa.col * width))
}

/// The sub-zone of `ancestor` at relative `depth` with scanline index `index`.
pub fn subzone_at(h: &Hierarchy, ancestor: Cid, depth: u32, index: u64) -> Result<Cid> {
    let pa = position(h, ancestor)?;
    let width = span(h, depth)?;
    if index >= width.saturating_mul(width) {
        return invalid(format!("index {index} outside the {} sub-zones at depth {depth}", width.saturating_mul(width)));
    }
    let (row, col) = (pa.row.checked_mul(width), pa.col.checked_mul(width));
    match (row, col) {
        (Some(row), Some(col)) => cell_at(h, pa.base, pa.level + depth, row + index / width, col + index % width),
        _ => invalid(format!("depth {depth} below level {} overflows", pa.level)),
    }
}

/// Every sub-zone of `ancestor` at relative `depth`, in scanline order.
pub fn subzones(h: &Hierarchy, ancestor: Cid, depth: u32) -> Result<Vec<Cid>> {
    let width = span(h, depth)?;
    let count = width.checked_mul(width).filter(|&c| c <= MAX_RASTER_CELLS);
    let count = count.ok_or_else(|| crate::Error::Invalid(format!("depth {depth} has too many sub-zones")))?;
    (0..count).map(|i| subzone_at(h, ancestor, depth, i)).collect()
}

/// The rHEALPix arrangement of the six base cells under a profile.
#[derive(Debug, Clone, Copy)]
struct Layout {
    north: i64,
    south: i64,
}

impl Layout {
    fn of(profile: &Profile) -> Result<Layout> {
        profile.validate()?;
        if profile.n_base != 6 {
            return invalid(format!("the rHEALPix layout has 6 base cells, the profile has {}", profile.n_base));
        }
        Ok(Layout { north: profile.north_square as i64, south: profile.south_square as i64 })
    }

    /// The cell at `(r, c)` relative to base cell `base` of side `n`, where the point may lie
    /// outside the base cell by up to `n` in one coordinate; `None` beyond a corner.
    fn locate(&self, base: u32, n: i64, r: i64, c: i64) -> Option<(u32, i64, i64)> {
        let inside = |x: i64| (0..n).contains(&x);
        let cw = |(r, c): (i64, i64)| (c, n - 1 - r);
        let ccw = |(r, c): (i64, i64)| (n - 1 - c, r);
        match (inside(r), inside(c)) {
            (true, true) => return Some((base, r, c)),
            (false, false) => return None,
            _ => {}
        }
        let equatorial = |column: i64| 1 + column.rem_euclid(4) as u32;
        match base {
            1..=4 => {
                let column = base as i64 - 1;
                if !inside(c) {
                    return Some((equatorial(column + c.signum()), r, c.rem_euclid(n)));
                }
                if r < 0 {
                    let mut p = (n + r, c);
                    for _ in 0..(column - self.north).rem_euclid(4) {
                        p = ccw(p);
                    }
                    return Some((0, p.0, p.1));
                }
                let mut p = (r - n, c);
                for _ in 0..(column - self.south).rem_euclid(4) {
                    p = cw(p);
                }
                Some((5, p.0, p.1))
            }
            0 => {
                // edges counted from the one facing the band: bottom, right, top, left
                let k = if r >= n { 0 } else if c >= n { 1 } else if r < 0 { 2 } else { 3 };
                let mut p = (r, c);
                for _ in 0..k {
                    p = cw(p);
                }
                Some((equatorial(self.north + k), p.0 - n, p.1))
            }
            _ => {
                // edges counted from the one facing the band: top, right, bottom, left
                let k = if r < 0 { 0 } else if c >= n { 1 } else if r >= n { 2 } else { 3 };
                let mut p = (r, c);
                for _ in 0..k {
                    p = ccw(p);
                }
                Some((equatorial(self.south + k), p.0 + n, p.1))
            }
        }
    }
}

/// The four edge neighbours of `cid` in `DIRECTIONS` order: up, right, down, left in the
/// cell's own planar frame. Every cell has all four, across base-cell edges included.
pub fn neighbours(profile: &Profile, cid: Cid) -> Result<[Cid; 4]> {
    let layout = Layout::of(profile)?;
    let h = profile.hierarchy();
    let p = position(&h, cid)?;
    let n = span(&h, p.level)? as i64;
    let mut out = [0; 4];
    for (slot, dir) in out.iter_mut().zip(DIRECTIONS) {
        let (dr, dc) = dir.step();
        let (b, r, c) = layout.locate(p.base, n, p.row as i64 + dr, p.col as i64 + dc).expect("an edge step never leaves by a corner");
        *slot = cell_at(&h, b, p.level, r as u64, c as u64)?;
    }
    Ok(out)
}

/// A gather index that pads each base-cell raster at `level` with `width` halo cells on every
/// side. Entry `(base, i, j)` of the padded `(6, n + 2w, n + 2w)` array, flattened, is the raster
/// index of the cell drawn there, or -1 in the `w × w` blocks beyond a base cell's corners.
pub fn halo_index(profile: &Profile, level: u32, width: u32) -> Result<Vec<i64>> {
    let layout = Layout::of(profile)?;
    let h = profile.hierarchy();
    let n = span(&h, level)?;
    let w = width as u64;
    if w > n {
        return invalid(format!("a halo of {w} is wider than a base cell of {n} at level {level}"));
    }
    let m = n + 2 * w;
    let total = m.checked_mul(m).and_then(|x| x.checked_mul(6)).filter(|&t| t <= MAX_RASTER_CELLS);
    let total = total.ok_or_else(|| crate::Error::Invalid(format!("a halo index at level {level} is too large")))?;
    let (n, w, m) = (n as i64, w as i64, m as i64);
    let mut out = Vec::with_capacity(total as usize);
    for b in 0..6u32 {
        for i in 0..m {
            for j in 0..m {
                out.push(match layout.locate(b, n, i - w, j - w) {
                    Some((b2, r, c)) => b2 as i64 * n * n + r * n + c,
                    None => -1,
                });
            }
        }
    }
    Ok(out)
}

fn raster_len(h: &Hierarchy, level: u32) -> Result<(u64, usize)> {
    let n = span(h, level)?;
    let total = n.checked_mul(n).and_then(|x| x.checked_mul(h.b as u64)).filter(|&t| t <= MAX_RASTER_CELLS);
    let total = total.ok_or_else(|| crate::Error::Invalid(format!("a raster at level {level} is too large")))?;
    Ok((n, total as usize))
}

/// The raster index of each cid in the raster at its level; all cids must share one level.
pub fn raster_index(h: &Hierarchy, cids: &[Cid]) -> Result<Vec<u64>> {
    let Some(&first) = cids.first() else { return Ok(Vec::new()) };
    let level = h.check(first, None)?;
    let (n, _) = raster_len(h, level)?;
    cids.iter()
        .map(|&c| {
            let p = position(h, c)?;
            if p.level != level {
                return invalid(format!("cid {c} is at level {}, the first at {level}", p.level));
            }
            Ok(p.base as u64 * n * n + p.row * n + p.col)
        })
        .collect()
}

/// The cid at every position of the raster at `level`, in raster order.
pub fn raster_cells(h: &Hierarchy, level: u32) -> Result<Vec<Cid>> {
    let (n, total) = raster_len(h, level)?;
    let mut out = Vec::with_capacity(total);
    for b in 0..h.b {
        for row in 0..n {
            for col in 0..n {
                out.push(cell_at(h, b, level, row, col)?);
            }
        }
    }
    Ok(out)
}

/// The covered leaves of `tree` as one raster per base cell at the tree's depth, flattened
/// base-major then scanline: 1 where covered, 0 elsewhere.
pub fn raster(tree: &Tree) -> Result<Vec<u8>> {
    let h = tree.h;
    let (n, total) = raster_len(&h, tree.d)?;
    let mut out = vec![0u8; total];
    for c in tree.cells() {
        let p = position(&h, c)?;
        let block = span(&h, tree.d - p.level)?;
        let (r0, c0) = (p.row * block, p.col * block);
        for r in r0..r0 + block {
            let start = (p.base as u64 * n * n + r * n + c0) as usize;
            out[start..start + block as usize].fill(1);
        }
    }
    Ok(out)
}

/// The canonical tree of depth `d` whose covered leaves are the nonzero entries of a raster laid
/// out as `raster` returns it.
pub fn from_raster(h: Hierarchy, d: u32, cells: &[u8], ctx: Arc<Ctx>) -> Result<Tree> {
    let (n, total) = raster_len(&h, d)?;
    if cells.len() != total {
        return invalid(format!("raster has {} cells, depth {d} needs {total}", cells.len()));
    }
    let a = side(&h)?;
    let at = |base: u64, row: u64, col: u64| cells[(base * n * n + row * n + col) as usize] != 0;
    let bases = (0..h.b as u64).map(|b| block_node(&at, a, b, 0, 0, 0, d, &ctx)).collect();
    let empty = Tree::empty(h, d, ctx)?;
    Ok(Tree { bases, ..empty })
}

/// The canonical node of the block at `(level, row, col)` of base cell `base`, from a leaf test.
#[allow(clippy::too_many_arguments)]
pub(crate) fn block_node(at: &dyn Fn(u64, u64, u64) -> bool, a: u64, base: u64, level: u32, row: u64, col: u64, d: u32, ctx: &Ctx) -> Node {
    if level == d {
        return if at(base, row, col) { Node::Full } else { Node::Empty };
    }
    let children = (0..a * a).map(|k| block_node(at, a, base, level + 1, row * a + k / a, col * a + k % a, d, ctx)).collect();
    make_branch(children, d - level, ctx)
}
