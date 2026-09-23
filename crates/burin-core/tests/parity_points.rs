//! Point → cell against rhealpixdggs-py (`points_*.json`): exact on planar edge and vertex
//! points, and on ellipsoidal points away from cell edges; the forward projection to a stated
//! tolerance; and every cell's nucleus lies in that cell.

use burin_core::hierarchy::{cid_to_suid, suid_to_cid, SPACE};
use burin_core::profile::Profile;
use burin_core::zone::{cell_at, cell_from_point, neighbours};
use serde_json::Value;
use std::collections::BTreeSet;
use std::path::PathBuf;

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
/// Random points per level in the fixture (`gen_points` in tools/gen_fixtures.py).
const RANDOM: usize = 160;

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
fn ellipsoidal_points_match_the_reference() {
    for name in PROFILES {
        let fx = fixture(&format!("points_{name}.json"));
        let profile = profile_of(&fx);
        let grid = profile.grid();
        let (mut exact_xy, mut worst, mut near_edge) = (0usize, 0f64, 0usize);
        let rows = fx["lonlat"].as_array().unwrap();
        // rows run level by level over the same points; the first RANDOM of each block are random,
        // the rest are placed on edges, poles and seams on purpose
        let levels: BTreeSet<u64> = rows.iter().map(|r| r[2].as_u64().unwrap()).collect();
        let per_level = rows.len() / levels.len();
        for (i, row) in rows.iter().enumerate() {
            let random = i % per_level < RANDOM;
            let (lon, lat, level) = (bits(&row[0]), bits(&row[1]), row[2].as_u64().unwrap() as u32);
            let (rx, ry) = (bits(&row[4]), bits(&row[5]));
            let (x, y) = grid.forward(lon, lat, None);
            exact_xy += usize::from(x == rx && y == ry);
            worst = worst.max((x - rx).abs().max((y - ry).abs()) / rx.abs().max(ry.abs()).max(1.0));
            let want = suid_to_cid(row[3].as_str().unwrap()).unwrap();
            let got = cell_from_point(&profile, lon, lat, level).unwrap();
            let path = SPACE.path(want).unwrap();
            let (ulx, uly) = grid.ul_vertex(&path);
            let w = grid.cell_width(level);
            let clearance = (rx - ulx).min(ulx + w - rx).min(uly - ry).min(ry - (uly - w));
            if clearance > 1e-6 {
                assert_eq!(got, want, "{name}: ({lon}, {lat}) at level {level}, {clearance:.3e} m inside");
            } else {
                near_edge += usize::from(random);
                assert!(got == want || neighbours(&profile, want).unwrap().contains(&got), "{name}: ({lon}, {lat}) at level {level}");
            }
        }
        assert!(worst <= 1e-14, "{name}: forward projection differs by {worst:e} relative");
        assert!(exact_xy * 10 >= rows.len() * 9, "{name}: only {exact_xy}/{} forward projections are bit-identical", rows.len());
        eprintln!("{name}: forward bit-identical {exact_xy}/{}, worst relative {worst:e}; random points within 1 µm of an edge: {near_edge}", rows.len());
        assert!(near_edge * 100 <= rows.len(), "{name}: {near_edge} random points too near an edge to compare exactly");
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
