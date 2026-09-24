//! Time (SPEC §9–§12) against `time.json`, which tools/gen_time_fixtures.py computes from the
//! specification's text alone, and the properties the rules promise.

use burin_core::hierarchy::{suid_to_cid, SPACE};
use burin_core::opening::{open_path, OpeningRecord};
use burin_core::profile::Profile;
use burin_core::setops::{holds, relation_transcript, verify_relation, SetRelation};
use burin_core::time::{self, allen, Allen, TICKS, TIME, TIME_DEPTH};
use burin_core::tree::Tree;
use serde_json::Value;
use std::path::PathBuf;
use std::sync::Arc;

fn fixture() -> Value {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/time.json");
    serde_json::from_str(&std::fs::read_to_string(p).expect("time.json")).unwrap()
}

fn ctx() -> Arc<burin_core::Ctx> {
    Arc::new(time::ctx(&Profile::ogc()).unwrap())
}

fn tree(intervals: &[(u64, u64)]) -> Tree {
    time::intervals_tree(intervals, ctx()).unwrap()
}

struct Rng(u64);
impl Rng {
    fn below(&mut self, n: u64) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 % n
    }
    /// A non-empty interval, often near the ends of the line or of a power of two.
    fn interval(&mut self) -> (u64, u64) {
        let mut point = || match self.below(4) {
            0 => self.below(64),
            1 => TICKS - self.below(64),
            2 => (1u64 << self.below(62)).saturating_sub(self.below(3)),
            _ => self.below(TICKS + 1),
        };
        let (a, b) = (point().min(TICKS), point().min(TICKS));
        if a == b {
            (a.min(TICKS - 1), a.min(TICKS - 1) + 1)
        } else {
            (a.min(b), a.max(b))
        }
    }
}

#[test]
fn ticks_and_their_starts_match_the_reference() {
    let fx = fixture();
    for (name, line) in fx["lines"].as_object().unwrap() {
        let p = Profile { tick_us: line["tick_us"].as_u64().unwrap(), epoch_us: line["epoch_us"].as_i64().unwrap(), ..Profile::ogc() };
        for row in line["ticks"].as_array().unwrap() {
            let t = row[0].as_i64().unwrap();
            assert_eq!(time::tick(&p, t).ok(), row[1].as_u64(), "{name}: tick({t})");
        }
        for row in line["starts"].as_array().unwrap() {
            let k = row[0].as_u64().unwrap();
            assert_eq!(time::tick_start(&p, k).ok(), row[1].as_i64(), "{name}: start of tick {k}");
        }
    }
}

#[test]
fn intervals_have_the_reference_cells_and_roots() {
    let fx = fixture();
    let empty = Tree::empty(TIME, TIME_DEPTH, ctx()).unwrap();
    assert_eq!(empty.root_hex(), fx["empty_root"].as_str().unwrap());
    assert_eq!(tree(&[(0, TICKS)]).root_hex(), fx["full_root"].as_str().unwrap());
    for case in fx["intervals"].as_array().unwrap() {
        let ivs: Vec<(u64, u64)> = case["intervals"].as_array().unwrap().iter().map(|i| (i[0].as_u64().unwrap(), i[1].as_u64().unwrap())).collect();
        let t = tree(&ivs);
        assert_eq!(t.root_hex(), case["root"].as_str().unwrap(), "{ivs:?}");
        let want: Vec<u64> = case["cells"].as_array().unwrap().iter().map(|c| c.as_u64().unwrap()).collect();
        let mut sorted = want.clone();
        sorted.sort_unstable(); // Tree::cells lists by id; the reference and the walk list by time
        assert_eq!(t.cells(), sorted, "{ivs:?}");
        if let [(lo, hi)] = ivs[..] {
            assert_eq!(time::interval_cells(lo, hi).unwrap(), want, "the walk of [{lo}, {hi})");
        }
    }
}

#[test]
fn the_walk_is_an_exact_minimal_cover() {
    let mut rng = Rng(0x7157);
    for _ in 0..20_000 {
        let (lo, hi) = rng.interval();
        let cells = time::interval_cells(lo, hi).unwrap();
        assert!(cells.len() <= 2 * TIME_DEPTH as usize);
        let mut at = lo;
        for &c in &cells {
            let (a, b) = time::cell_ticks(c).unwrap();
            assert_eq!(a, at, "[{lo}, {hi}): cells must tile in order");
            at = b;
        }
        assert_eq!(at, hi);
        for w in cells.windows(2) {
            let ((l0, k0), (l1, k1)) = (time::level_index(w[0]).unwrap(), time::level_index(w[1]).unwrap());
            assert!(!(l0 == l1 && k0 / 2 == k1 / 2), "[{lo}, {hi}): siblings {} and {} should be their parent", w[0], w[1]);
        }
    }
}

#[test]
fn a_root_does_not_depend_on_how_the_set_is_split() {
    let mut rng = Rng(0xa11e);
    for _ in 0..500 {
        let (lo, hi) = rng.interval();
        let whole = tree(&[(lo, hi)]);
        let mid = lo + rng.below(hi - lo);
        if mid > lo {
            assert_eq!(tree(&[(mid, hi), (lo, mid)]).root(), whole.root(), "[{lo}, {hi}) split at {mid}");
        }
        let (a, b) = (lo + rng.below(hi - lo), lo + rng.below(hi - lo));
        let overlapping = [(lo, a.max(b) + 1), (a.min(b), hi), (lo, hi)];
        assert_eq!(tree(&overlapping).root(), whole.root(), "[{lo}, {hi}) with overlapping parts");
    }
}

#[test]
fn coarse_cells_are_those_whose_midpoint_is_inside() {
    let mut rng = Rng(0xc0a5);
    for _ in 0..5_000 {
        let (lo, hi) = rng.interval();
        let level = rng.below(TIME_DEPTH as u64 + 1) as u32;
        let s = TIME_DEPTH - level;
        let mid = |k: u64| (k << s) + if s == 0 { 0 } else { 1 << (s - 1) };
        if (hi - lo) >> s > 100_000 {
            continue;
        }
        let cells = time::coarse_cells(lo, hi, level).unwrap();
        for &c in &cells {
            let (l, k) = time::level_index(c).unwrap();
            assert!(l == level && (lo..hi).contains(&mid(k)));
        }
        let (Some(&f), Some(&l)) = (cells.first(), cells.last()) else { continue };
        let (first, last) = (time::level_index(f).unwrap().1, time::level_index(l).unwrap().1);
        assert!(first == 0 || mid(first - 1) < lo, "the cell before is outside");
        assert!(last + 1 == 1 << level || mid(last + 1) >= hi, "the cell after is outside");
    }
    assert_eq!(time::coarse_cells(5, 9, TIME_DEPTH).unwrap(), (5..9).map(|t| time::tick_cell(t).unwrap()).collect::<Vec<_>>());
}

#[test]
fn allen_matches_the_reference_and_inverts() {
    let fx = fixture();
    let pair = |v: &Value| (v[0].as_u64().unwrap(), v[1].as_u64().unwrap());
    let mut seen = std::collections::BTreeSet::new();
    for row in fx["allen"].as_array().unwrap() {
        let (a, b) = (pair(&row[0]), pair(&row[1]));
        let r = allen(a, b).unwrap();
        assert_eq!(r.name(), row[2].as_str().unwrap(), "{a:?} vs {b:?}");
        assert_eq!(allen(b, a).unwrap(), r.inverse());
        seen.insert(r.name());
    }
    assert_eq!(seen.len(), 13, "the small cases reach every relation");
    assert!(Allen::ALL.iter().all(|r| r.inverse().inverse() == *r));
    assert!(allen((3, 3), (1, 2)).is_err() && allen((1, 2), (5, 4)).is_err());
}

#[test]
fn off_the_line_is_refused() {
    let p = Profile::ogc();
    assert!(time::interval_cells(4, 4).is_err());
    assert!(time::interval_cells(0, TICKS + 1).is_err());
    assert!(time::cell(3, 8).is_err() && time::cell(62, 0).is_err());
    assert!(time::tick_cell(TICKS).is_err());
    assert!(time::level_index(1 << 63).is_err(), "a level-62 id is below the time depth");
    assert!(time::tick(&Profile { tick_us: 0, ..p.clone() }, 0).is_err());
    assert!(time::coarse_cells(0, TICKS, TIME_DEPTH).is_err(), "2^61 cells are refused, not allocated");
}

#[test]
fn set_relations_have_evidence_that_verifies_exactly_when_they_hold() {
    let mut rng = Rng(0x5e7);
    let space = |cells: &[&str]| {
        let c = Arc::new(Profile::ogc().ctx(3).unwrap());
        Tree::from_cells(SPACE, 3, cells.iter().map(|s| suid_to_cid(s).unwrap()), c).unwrap()
    };
    let mut pairs = vec![
        (space(&["Q4"]), space(&["Q41", "Q418"])),
        (space(&["Q41"]), space(&["Q4", "N"])),
        (space(&["Q4", "N0"]), space(&["S"])),
        (space(&["P01"]), space(&["P01"])),
    ];
    for _ in 0..40 {
        let (a, b, c) = (rng.interval(), rng.interval(), rng.interval());
        pairs.push((tree(&[a]), tree(&[b, c])));
    }
    let (a, b) = (rng.interval(), rng.interval());
    pairs.push((tree(&[a]), tree(&[a, b])));
    let stranger = [7u8; 32];
    for (x, y) in &pairs {
        for rel in SetRelation::ALL {
            let proof = relation_transcript(x, y, rel).unwrap();
            let sent = proof.as_ref().map(|p| burin_core::SetOpProof::from_json(&p.to_json()).unwrap());
            assert_eq!(verify_relation(rel, &x.root(), &y.root(), sent.as_ref()), holds(x, y, rel).unwrap(), "{}", rel.name());
            if let Some(p) = &proof {
                // evidence about other sets, or offered for another relation, proves nothing
                assert!(!verify_relation(rel, &stranger, &y.root(), Some(p)));
                assert!(!verify_relation(rel, &x.root(), &stranger, Some(p)));
                assert!(!verify_relation(SetRelation::Equals, &x.root(), &y.root(), Some(p)));
            }
        }
    }
    assert!(!verify_relation(SetRelation::Disjoint, &pairs[0].0.root(), &pairs[0].1.root(), None));
}

#[test]
fn time_trees_open_like_spatial_ones_and_do_not_mix_with_them() {
    let t = tree(&[(1_000, 5_000), (7_000, 7_001)]);
    for (tick, inside) in [(1_000, true), (4_999, true), (5_000, false), (7_000, true), (0, false)] {
        let op = open_path(&t, time::tick_cell(tick).unwrap()).unwrap().unwrap();
        let rec = OpeningRecord::from_json(&OpeningRecord::new(&t, op).to_json()).unwrap();
        assert!(rec.verify());
        assert_eq!(rec.opening.terminal() == Some(burin_core::Claim::Full), inside, "tick {tick}");
    }
    let s = Tree::empty(SPACE, 3, Arc::new(Profile::ogc().ctx(3).unwrap())).unwrap();
    assert!(burin_core::union(&t, &s).is_err());
}
