//! Point → cell against rhealpixdggs-py (`points_*.json`), exactly: planar points on every kind of
//! edge and vertex; ellipsoidal points and the forward projection, where each difference from the
//! reference is listed in `reference_disagreements.json` to the ulp; and every nucleus lies in its cell.

use burin_core::hierarchy::{cid_to_suid, SPACE};
use burin_core::profile::Profile;
use burin_core::zone::{cell_at, cell_from_point};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

mod support;
use support::goldens::ulps;

fn fixture(name: &str) -> Value {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures").join(name);
    serde_json::from_str(&std::fs::read_to_string(&p).unwrap_or_else(|_| panic!("missing fixture {}", p.display()))).unwrap()
}

fn bits(v: &Value) -> f64 {
    f64::from_bits(u64::from_str_radix(v.as_str().unwrap(), 16).unwrap())
}

fn profile_of(fx: &Value) -> Profile {
    let p = &fx["profile"];
    Profile {
        lon_0_udeg: (p["lon_0"].as_f64().unwrap() * 1e6).round() as i64,
        north_square: p["ns"].as_u64().unwrap() as u8,
        south_square: p["ss"].as_u64().unwrap() as u8,
        ..Profile::default()
    }
}

const PROFILES: [&str; 3] = ["ogc", "burin1", "ns1ss2"];

#[test]
fn planar_points_match_the_reference_to_the_bit() {
    for name in PROFILES {
        let fx = fixture(&format!("points_{name}.json"));
        let (profile, grid) = (profile_of(&fx), profile_of(&fx).grid());
        for row in fx["planar"].as_array().unwrap() {
            let (x, y, level) = (bits(&row[0]), bits(&row[1]), row[2].as_u64().unwrap() as u32);
            let got = match grid.cell_from_planar(x, y, level) {
                Some((b, r, c)) => cid_to_suid(cell_at(&profile.hierarchy(), b, level, r, c).unwrap()).unwrap(),
                None => String::new(),
            };
            assert_eq!(got, row[3].as_str().unwrap(), "{name}: ({x:e}, {y:e}) at level {level}");
        }
    }
}

#[test]
fn ellipsoidal_points_match_the_reference_exactly_or_as_listed() {
    let listed_all = fixture("reference_disagreements.json");
    for name in PROFILES {
        let fx = fixture(&format!("points_{name}.json"));
        let profile = profile_of(&fx);
        let grid = profile.grid();
        let listed = &listed_all[name];
        let key = |r: &Value, n: usize| r.as_array().unwrap()[..n].iter().map(|v| v.to_string()).collect::<Vec<_>>().join("|");
        let lookups: BTreeMap<String, String> = listed["lookups"].as_array().unwrap().iter().map(|r| (key(r, 3), r[4].as_str().unwrap().to_string())).collect();
        let forward: BTreeMap<String, (u64, u64)> = listed["forward"].as_array().unwrap().iter().map(|r| (key(r, 2), (r[2].as_u64().unwrap(), r[3].as_u64().unwrap()))).collect();
        let (mut seen_lookups, mut seen_forward) = (BTreeSet::new(), BTreeSet::new());
        for row in fx["lonlat"].as_array().unwrap() {
            let (lon, lat, level) = (bits(&row[0]), bits(&row[1]), row[2].as_u64().unwrap() as u32);
            let want = row[3].as_str().unwrap();
            let got = cid_to_suid(cell_from_point(&profile, lon, lat, level).unwrap()).unwrap();
            match lookups.get(&key(row, 3)) {
                Some(listed_got) => {
                    assert_eq!(&got, listed_got, "{name}: ({lon}, {lat}) at level {level} is listed as {listed_got}");
                    seen_lookups.insert(key(row, 3));
                }
                None => assert_eq!(got, want, "{name}: ({lon}, {lat}) at level {level} differs from the reference and is not listed"),
            }
            if level == 0 {
                let (x, y) = grid.forward(lon, lat, None);
                let d = (ulps(x, bits(&row[4])), ulps(y, bits(&row[5])));
                match forward.get(&key(row, 2)) {
                    Some(&expected) => {
                        assert_eq!(d, expected, "{name}: forward({lon}, {lat}) differs by {d:?} ulps, listed as {expected:?}");
                        seen_forward.insert(key(row, 2));
                    }
                    None => assert_eq!(d, (0, 0), "{name}: forward({lon}, {lat}) differs from the reference and is not listed"),
                }
            }
        }
        assert_eq!(seen_lookups.len(), lookups.len(), "{name}: a listed lookup is no longer in the fixture or no longer differs");
        assert_eq!(seen_forward.len(), forward.len(), "{name}: a listed projection is no longer in the fixture or no longer differs");
    }
}

#[test]
fn every_nucleus_lies_in_its_cell() {
    for name in PROFILES {
        let fx = fixture(&format!("points_{name}.json"));
        let profile = profile_of(&fx);
        let grid = profile.grid();
        for level in 0..=8u32 {
            let cells: Vec<u64> = (0..6u64).flat_map(|b| SPACE.descendants(9 + b, level)).collect();
            for &c in cells.iter().step_by((cells.len() / 400).max(1)) {
                let (lon, lat) = grid.nucleus(&SPACE.path(c).unwrap()).unwrap();
                assert_eq!(cell_from_point(&profile, lon, lat, level).unwrap(), c, "{name}: nucleus of {}", cid_to_suid(c).unwrap());
            }
        }
    }
}

#[test]
fn points_off_the_grid_are_errors() {
    let p = Profile::ogc();
    for (lon, lat) in [(f64::NAN, 0.0), (0.0, f64::NAN), (f64::INFINITY, 10.0)] {
        assert!(cell_from_point(&p, lon, lat, 3).is_err());
    }
}
