//! The time line (SPEC §9–§12): ticks, the time hierarchy, the interval rule and Allen's relations.
//!
//! ```text
//! tick(t)   = floor((t − epoch_us) / tick_us) + 2^60       t in POSIX µs, 0 ≤ tick < 2^61
//! cell(r,k) = ticks [k·2^(61−r), (k+1)·2^(61−r)),   cid 2^(r+1) + k,   TIME = H(2, 1), depth 61
//! ```
//! Everything here is integer arithmetic, so it is exact on every platform.

use crate::error::{invalid, Result};
use crate::hash::{hasher_for, Ctx};
use crate::hierarchy::{Cid, Hierarchy};
use crate::profile::Profile;
use crate::tree::Tree;
use std::sync::Arc;

/// The time hierarchy: aperture 2 over one base cell, the whole line.
pub const TIME: Hierarchy = Hierarchy { a: 2, b: 1 };
/// The depth of every time tree; its leaves are ticks.
pub const TIME_DEPTH: u32 = 61;
/// The number of ticks on the line.
pub const TICKS: u64 = 1 << TIME_DEPTH;
/// The tick of the profile's epoch: the middle of the line.
pub const EPOCH_TICK: u64 = 1 << (TIME_DEPTH - 1);
/// The OGC-registered temporal CRS the line counts in.
pub const UNIX_TIME_CRS: &str = "https://www.opengis.net/def/crs/OGC/0/UnixTime";

fn tick_us(profile: &Profile) -> Result<i128> {
    if profile.tick_us == 0 {
        return invalid("tick_us must be >= 1");
    }
    Ok(profile.tick_us as i128)
}

/// The tick containing the POSIX instant `posix_us`, under the profile's tick and epoch.
pub fn tick(profile: &Profile, posix_us: i64) -> Result<u64> {
    let t = (posix_us as i128 - profile.epoch_us as i128).div_euclid(tick_us(profile)?) + EPOCH_TICK as i128;
    if !(0..TICKS as i128).contains(&t) {
        return invalid(format!("instant {posix_us} µs is off the time line of this profile"));
    }
    Ok(t as u64)
}

/// The POSIX instant at which `tick` starts; `tick_start(TICKS)` is where the line ends.
pub fn tick_start(profile: &Profile, tick: u64) -> Result<i64> {
    if tick > TICKS {
        return invalid(format!("tick {tick} is off the time line (0..2^61)"));
    }
    let t = (tick as i128 - EPOCH_TICK as i128) * tick_us(profile)? + profile.epoch_us as i128;
    i64::try_from(t).map_err(|_| crate::Error::Invalid(format!("tick {tick} starts outside the 64-bit POSIX range")))
}

/// The cid of cell `k` at `level`.
pub fn cell(level: u32, k: u64) -> Result<Cid> {
    if level > TIME_DEPTH || k >> level != 0 {
        return invalid(format!("no time cell {k} at level {level}"));
    }
    Ok((2u64 << level) + k)
}

/// The cid of the level-61 cell that is `tick`.
pub fn tick_cell(tick: u64) -> Result<Cid> {
    cell(TIME_DEPTH, tick)
}

/// `(level, k)` of a time cid.
pub fn level_index(cid: Cid) -> Result<(u32, u64)> {
    let level = TIME.check(cid, Some(TIME_DEPTH))?;
    Ok((level, cid - (2u64 << level)))
}

/// The ticks `[lo, hi)` a time cell spans.
pub fn cell_ticks(cid: Cid) -> Result<(u64, u64)> {
    let (level, k) = level_index(cid)?;
    let s = TIME_DEPTH - level;
    Ok((k << s, (k + 1) << s))
}

fn check_interval(lo: u64, hi: u64) -> Result<()> {
    if !(lo < hi && hi <= TICKS) {
        return invalid(format!("[{lo}, {hi}) is not a non-empty interval of the time line"));
    }
    Ok(())
}

/// The interval rule (§11): the compact cells of `[lo, hi)`, found by the HINT boundary walk,
/// in order of their first tick. At most 2·61 cells.
pub fn interval_cells(lo: u64, hi: u64) -> Result<Vec<Cid>> {
    check_interval(lo, hi)?;
    let (mut lo, mut hi, mut level) = (lo, hi, TIME_DEPTH);
    let (mut left, mut right) = (Vec::new(), Vec::new());
    while lo < hi {
        if lo & 1 == 1 {
            left.push((2u64 << level) + lo);
            lo += 1;
        }
        if hi & 1 == 1 {
            hi -= 1;
            right.push((2u64 << level) + hi);
        }
        (lo, hi) = (lo >> 1, hi >> 1);
        level = level.wrapping_sub(1); // only after the loop can it pass level 0
    }
    left.extend(right.into_iter().rev());
    Ok(left)
}

/// The most cells `coarse_cells` returns in one call.
pub const MAX_COARSE_CELLS: u64 = 1 << 28;

/// The level-`level` cells whose midpoint tick lies in `[lo, hi)` (§11, coarse cells).
pub fn coarse_cells(lo: u64, hi: u64, level: u32) -> Result<Vec<Cid>> {
    check_interval(lo, hi)?;
    if level > TIME_DEPTH {
        return invalid(format!("time level {level} is deeper than {TIME_DEPTH}"));
    }
    let s = TIME_DEPTH - level;
    let half = if s == 0 { 0 } else { 1u64 << (s - 1) };
    // cell k's midpoint is k·2^s + half; the cells wanted run from the first midpoint ≥ lo to the first ≥ hi
    let first = |t: u64| if t <= half { 0 } else { (t - half).div_ceil(1u64 << s) };
    let (k0, k1) = (first(lo), first(hi).min(1u64 << level));
    if k1.saturating_sub(k0) > MAX_COARSE_CELLS {
        return invalid(format!("{} cells at level {level}, more than MAX_COARSE_CELLS", k1 - k0));
    }
    Ok((k0..k1).map(|k| (2u64 << level) + k).collect())
}

/// A hashing context for time trees under the profile's hash.
pub fn ctx(profile: &Profile) -> Result<Ctx> {
    let h = hasher_for(&profile.hash).ok_or_else(|| crate::Error::Invalid(format!("unknown hash id {:?}", profile.hash)))?;
    Ok(Ctx::new(h, TIME.a, TIME.b, TIME_DEPTH))
}

/// The canonical time tree of a union of intervals.
pub fn intervals_tree(intervals: &[(u64, u64)], ctx: Arc<Ctx>) -> Result<Tree> {
    let mut cells = Vec::new();
    for &(lo, hi) in intervals {
        cells.extend(interval_cells(lo, hi)?);
    }
    Tree::from_cells(TIME, TIME_DEPTH, cells, ctx)
}

/// Allen's thirteen relations between two non-empty intervals (§12).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Allen {
    Before,
    Meets,
    Overlaps,
    Starts,
    During,
    Finishes,
    Equals,
    FinishedBy,
    Contains,
    StartedBy,
    OverlappedBy,
    MetBy,
    After,
}

impl Allen {
    /// Every relation, listed so that the inverse of `ALL[i]` is `ALL[12 − i]`.
    pub const ALL: [Allen; 13] = [
        Allen::Before,
        Allen::Meets,
        Allen::Overlaps,
        Allen::Starts,
        Allen::During,
        Allen::Finishes,
        Allen::Equals,
        Allen::FinishedBy,
        Allen::Contains,
        Allen::StartedBy,
        Allen::OverlappedBy,
        Allen::MetBy,
        Allen::After,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Allen::Before => "before",
            Allen::Meets => "meets",
            Allen::Overlaps => "overlaps",
            Allen::Starts => "starts",
            Allen::During => "during",
            Allen::Finishes => "finishes",
            Allen::Equals => "equals",
            Allen::FinishedBy => "finished-by",
            Allen::Contains => "contains",
            Allen::StartedBy => "started-by",
            Allen::OverlappedBy => "overlapped-by",
            Allen::MetBy => "met-by",
            Allen::After => "after",
        }
    }

    /// The relation of `b` to `a` when this is the relation of `a` to `b`.
    pub fn inverse(self) -> Allen {
        Allen::ALL[12 - Allen::ALL.iter().position(|&r| r == self).expect("listed")]
    }
}

/// How `[a.0, a.1)` stands to `[b.0, b.1)`. Both must be non-empty.
pub fn allen(a: (u64, u64), b: (u64, u64)) -> Result<Allen> {
    use std::cmp::Ordering::*;
    let ((a1, a2), (b1, b2)) = (a, b);
    if a1 >= a2 || b1 >= b2 {
        return invalid("Allen relations are between non-empty intervals");
    }
    Ok(match (a1.cmp(&b1), a2.cmp(&b2)) {
        (Equal, Equal) => Allen::Equals,
        (Equal, Less) => Allen::Starts,
        (Equal, Greater) => Allen::StartedBy,
        (Greater, Equal) => Allen::Finishes,
        (Less, Equal) => Allen::FinishedBy,
        (Greater, Less) => Allen::During,
        (Less, Greater) => Allen::Contains,
        (Less, Less) => match a2.cmp(&b1) {
            Less => Allen::Before,
            Equal => Allen::Meets,
            Greater => Allen::Overlaps,
        },
        (Greater, Greater) => match a1.cmp(&b2) {
            Greater => Allen::After,
            Equal => Allen::MetBy,
            Less => Allen::OverlappedBy,
        },
    })
}
