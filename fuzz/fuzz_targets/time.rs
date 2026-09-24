//! The time line from arbitrary bits: ticks, the interval walk, coarse cells and Allen's relations
//! refuse or answer, never panic, and what they answer is consistent.
#![no_main]
use burin_core::profile::Profile;
use burin_core::time::{self, allen, TICKS};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let word = |i: usize| data.get(8 * i..8 * i + 8).map(|b| u64::from_le_bytes(b.try_into().unwrap()));
    let (Some(a), Some(b)) = (word(0), word(1)) else { return };
    let p = Profile { tick_us: word(2).unwrap_or(1), epoch_us: word(3).unwrap_or(0) as i64, ..Profile::ogc() };
    if let Ok(k) = time::tick(&p, a as i64) {
        let start = time::tick_start(&p, k);
        if let Ok(s) = start {
            assert!(s <= a as i64, "a tick starts at or before its instants");
            assert_eq!(time::tick(&p, s).unwrap(), k);
        }
    }
    let (lo, hi) = (a.min(b) % (TICKS + 1), a.max(b) % (TICKS + 1));
    if let Ok(cells) = time::interval_cells(lo, hi) {
        let mut at = lo;
        for c in &cells {
            let (x, y) = time::cell_ticks(*c).unwrap();
            assert_eq!(x, at);
            at = y;
        }
        assert_eq!(at, hi);
        assert!(cells.len() <= 122);
    }
    if let Some(&level) = data.get(32) {
        let level = u32::from(level) % 62;
        if (hi.saturating_sub(lo)) >> (61 - level) < 100_000 {
            let _ = time::coarse_cells(lo, hi, level);
        }
    }
    if let (Some(c), Some(d)) = (word(2), word(3)) {
        if let (Ok(r), Ok(s)) = (allen((lo, hi), (c.min(d), c.max(d))), allen((c.min(d), c.max(d)), (lo, hi))) {
            assert_eq!(r.inverse(), s);
        }
    }
});
