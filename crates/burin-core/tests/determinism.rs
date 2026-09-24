//! burin against its own frozen answers, bit for bit: the cell containing each adversarial point
//! at every level, the forward projection, and the nuclei (`points_burin.json`, written by
//! `examples/freeze_points.rs`). No tolerance anywhere; CI runs this on every release platform.

mod support;

use burin_core::hierarchy::SPACE;
use serde_json::Value;
use std::path::PathBuf;
use support::goldens::{bits, lookup, profiles, self_goldens, unbits, DEEPEST};

fn committed() -> Value {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/points_burin.json");
    serde_json::from_str(&std::fs::read_to_string(p).expect("missing points_burin.json")).unwrap()
}

#[test]
fn the_goldens_are_what_this_build_computes() {
    let (now, then) = (self_goldens(), committed());
    for (name, _) in profiles() {
        for key in ["points", "forward", "nuclei"] {
            let (a, b) = (now["profiles"][name][key].as_array().unwrap(), then["profiles"][name][key].as_array().unwrap());
            assert_eq!(a.len(), b.len(), "{name}/{key}: a different number of rows");
            if let Some(i) = (0..a.len()).find(|&i| a[i] != b[i]) {
                panic!("{name}/{key}[{i}]: computed {} but frozen {} (rerun freeze_points only for an intended change)", a[i], b[i]);
            }
        }
    }
    assert_eq!(now, then);
}

#[test]
fn every_level_of_every_frozen_point_is_the_exact_ancestor() {
    let then = committed();
    for (name, profile) in profiles() {
        let grid = profile.grid();
        for row in then["profiles"][name]["points"].as_array().unwrap() {
            let (lon, lat, deepest) = (unbits(&row[0]), unbits(&row[1]), row[2].as_u64().unwrap());
            for level in 0..=DEEPEST {
                let want = (deepest != 0).then(|| deepest / 9u64.pow(DEEPEST - level));
                assert_eq!(lookup(&grid, lon, lat, level), want, "{name}: ({}, {}) at level {level}", bits(lon), bits(lat));
            }
        }
        for row in then["profiles"][name]["nuclei"].as_array().unwrap() {
            let cid = row[0].as_u64().unwrap();
            let (lon, lat) = grid.nucleus(&SPACE.path(cid).unwrap()).unwrap();
            assert_eq!((bits(lon), bits(lat)), (row[1].as_str().unwrap().to_string(), row[2].as_str().unwrap().to_string()), "{name}: nucleus of {cid}");
        }
    }
}
