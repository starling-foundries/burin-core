//! Zone topology against both references: rhealpixdggs-py (`topology.json`, every polar
//! placement) and DGGAL, the OGC API - DGGS reference library (`dggal_rhealpix.json`).

use burin_core::hash::Ctx;
use burin_core::hierarchy::{cid_to_suid, suid_to_cid, Cid, Hierarchy, SPACE};
use burin_core::profile::Profile;
use burin_core::tree::Tree;
use burin_core::zone::{self, cell_at, halo_index, neighbours, position, subzone_at, subzone_index, subzones};
use serde_json::Value;
use std::collections::BTreeSet;
use std::path::PathBuf;
use std::sync::Arc;

fn fixture(name: &str) -> Value {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures").join(name);
    serde_json::from_str(&std::fs::read_to_string(&p).unwrap_or_else(|_| panic!("missing fixture {}", p.display()))).unwrap()
}

fn placement(ns: u8, ss: u8) -> Profile {
    Profile { north_square: ns, south_square: ss, ..Profile::default() }
}

fn cells_at(level: u32) -> Vec<Cid> {
    (0..6u64).flat_map(|b| SPACE.descendants(9 + b, level)).collect()
}

fn suid(c: Cid) -> String {
    cid_to_suid(c).unwrap()
}

fn strings(v: &Value) -> Vec<String> {
    v.as_array().unwrap().iter().map(|s| s.as_str().unwrap().to_string()).collect()
}

#[test]
fn position_and_cell_at_are_inverse() {
    for c in cells_at(3) {
        for cid in [c, SPACE.descendants(c, 4).start + 1234] {
            let p = position(&SPACE, cid).unwrap();
            assert_eq!(cell_at(&SPACE, p.base, p.level, p.row, p.col).unwrap(), cid);
        }
    }
    assert!(cell_at(&SPACE, 6, 1, 0, 0).is_err());
    assert!(cell_at(&SPACE, 0, 1, 3, 0).is_err());
}

#[test]
fn scanline_order_matches_the_planar_layout_under_every_profile() {
    let fx = fixture("topology.json");
    for (key, want) in fx["scanline"].as_object().unwrap() {
        let parts: Vec<&str> = key.split('/').collect();
        let got: Vec<String> = subzones(&SPACE, suid_to_cid(parts[1]).unwrap(), parts[2].parse().unwrap()).unwrap().into_iter().map(suid).collect();
        assert_eq!(got, strings(want), "{key}");
    }
}

#[test]
fn scanline_order_matches_dggal() {
    let fx = fixture("dggal_rhealpix.json");
    for (key, want) in fx["subzones"].as_object().unwrap() {
        let (parent, depth) = key.split_once('/').unwrap();
        let got: Vec<String> = subzones(&SPACE, suid_to_cid(parent).unwrap(), depth.parse().unwrap()).unwrap().into_iter().map(suid).collect();
        assert_eq!(got, strings(want), "{key}");
    }
}

#[test]
fn subzone_index_inverts_subzone_at() {
    for parent in ["N", "Q4", "S08", "R261"] {
        let p = suid_to_cid(parent).unwrap();
        for depth in 0..=3 {
            let subs = subzones(&SPACE, p, depth).unwrap();
            assert_eq!(subs.len() as u64, 9u64.pow(depth));
            assert_eq!(subs.iter().collect::<BTreeSet<_>>().len(), subs.len(), "{parent}/{depth} repeats a sub-zone");
            for (i, &s) in subs.iter().enumerate() {
                assert!(SPACE.is_ancestor(p, s));
                assert_eq!(subzone_index(&SPACE, p, s).unwrap(), i as u64);
                assert_eq!(subzone_at(&SPACE, p, depth, i as u64).unwrap(), s);
            }
            assert!(subzone_at(&SPACE, p, depth, 9u64.pow(depth)).is_err());
        }
    }
    let (q4, r0) = (suid_to_cid("Q4").unwrap(), suid_to_cid("R0").unwrap());
    assert!(subzone_index(&SPACE, q4, r0).is_err(), "a cell outside the ancestor has no index");
    assert!(subzone_index(&SPACE, suid_to_cid("Q45").unwrap(), q4).is_err(), "a descendant is not an ancestor");
}

#[test]
fn neighbours_match_the_reference_under_every_polar_placement() {
    let fx = fixture("topology.json");
    for ns in 0..4u8 {
        for ss in 0..4u8 {
            let p = placement(ns, ss);
            for row in fx["neighbours"][format!("{ns}{ss}")].as_array().unwrap() {
                let row = strings(row);
                let got: Vec<String> = neighbours(&p, suid_to_cid(&row[0]).unwrap()).unwrap().into_iter().map(suid).collect();
                assert_eq!(got, row[1..], "{} under ns={ns} ss={ss} (up, right, down, left)", row[0]);
            }
        }
    }
}

#[test]
fn neighbours_match_dggal() {
    let fx = fixture("dggal_rhealpix.json");
    let ogc = Profile::ogc();
    for (z, want) in fx["neighbours"].as_object().unwrap() {
        let got: BTreeSet<String> = neighbours(&ogc, suid_to_cid(z).unwrap()).unwrap().into_iter().map(suid).collect();
        assert_eq!(got, strings(want).into_iter().collect(), "{z}");
    }
}

#[test]
fn neighbours_are_symmetric_distinct_and_level_preserving() {
    for (ns, ss) in [(0, 0), (1, 2), (3, 3)] {
        let p = placement(ns, ss);
        for c in cells_at(3) {
            let nb = neighbours(&p, c).unwrap();
            assert_eq!(nb.iter().collect::<BTreeSet<_>>().len(), 4, "{} has a repeated neighbour", suid(c));
            for &n in &nb {
                assert_ne!(n, c);
                assert_eq!(SPACE.level(n), SPACE.level(c));
                assert!(neighbours(&p, n).unwrap().contains(&c), "{} -> {} is one-way (ns={ns} ss={ss})", suid(c), suid(n));
            }
        }
    }
}

fn raster_cid(level: u32, n: i64, idx: i64) -> Cid {
    let (b, rest) = (idx / (n * n), idx % (n * n));
    cell_at(&SPACE, b as u32, level, (rest / n) as u64, (rest % n) as u64).unwrap()
}

#[test]
fn a_width_one_halo_is_the_edge_neighbours() {
    for ns in 0..4u8 {
        for ss in 0..4u8 {
            let p = placement(ns, ss);
            for level in [0, 1, 2] {
                let idx = halo_index(&p, level, 1).unwrap();
                let n = 3i64.pow(level);
                let m = n + 2;
                assert_eq!(idx.len() as i64, 6 * m * m);
                assert_eq!(idx.iter().filter(|&&x| x < 0).count(), 24, "one missing diagonal per base-cell corner");
                for b in 0..6i64 {
                    let at = |i: i64, j: i64| idx[(b * m * m + i * m + j) as usize];
                    for r in 0..n {
                        for c in 0..n {
                            let own = raster_cid(level, n, at(r + 1, c + 1));
                            assert_eq!(own, cell_at(&SPACE, b as u32, level, r as u64, c as u64).unwrap());
                            let nb = neighbours(&p, own).unwrap();
                            let around = [at(r, c + 1), at(r + 1, c + 2), at(r + 2, c + 1), at(r + 1, c)];
                            for (k, &g) in around.iter().enumerate() {
                                assert_eq!(raster_cid(level, n, g), nb[k], "halo of {} (ns={ns} ss={ss})", suid(own));
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn a_wide_halo_continues_across_the_edge() {
    let p = placement(1, 2);
    let (level, width) = (2u32, 4u32);
    let idx = halo_index(&p, level, width).unwrap();
    let (n, w) = (3i64.pow(level), width as i64);
    let m = n + 2 * w;
    assert_eq!(idx.iter().filter(|&&x| x < 0).count() as i64, 24 * w * w);
    for b in 0..6i64 {
        let at = |i: i64, j: i64| idx[(b * m * m + i * m + j) as usize];
        // each strip from an edge cell outward: consecutive cells are edge neighbours, all distinct
        for t in 0..n {
            let strips: [Vec<i64>; 4] = [
                (0..=w).rev().map(|s| at(s, w + t)).collect(),
                (0..=w).map(|s| at(w + t, n + w - 1 + s)).collect(),
                (0..=w).map(|s| at(n + w - 1 + s, w + t)).collect(),
                (0..=w).rev().map(|s| at(w + t, s)).collect(),
            ];
            for strip in strips {
                let cells: Vec<Cid> = strip.iter().map(|&g| raster_cid(level, n, g)).collect();
                assert_eq!(cells.iter().collect::<BTreeSet<_>>().len(), cells.len());
                for pair in cells.windows(2) {
                    assert!(neighbours(&p, pair[0]).unwrap().contains(&pair[1]), "{} then {}", suid(pair[0]), suid(pair[1]));
                }
            }
        }
    }
    assert!(halo_index(&p, 1, 4).is_err(), "a halo wider than the base cell");
}

#[test]
fn raster_round_trips_and_counts_leaves() {
    let h = SPACE;
    let d = 3;
    let ctx = Arc::new(Ctx::sha256(9, 6, d));
    let mut x = 0x9e3779b97f4a7c15u64;
    let mut cells = Vec::new();
    for _ in 0..60 {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        let level = (x % 4) as u32;
        let mut path = vec![(x >> 8) as u32 % 6];
        for k in 0..level {
            path.push((x >> (12 + 4 * k)) as u32 % 9);
        }
        cells.push(h.cid(&path).unwrap());
    }
    let tree = Tree::from_cells(h, d, cells, ctx.clone()).unwrap();
    let r = zone::raster(&tree).unwrap();
    assert_eq!(r.len(), 6 * 27 * 27);
    assert_eq!(r.iter().filter(|&&v| v == 1).count() as u64, tree.leaf_count());
    for leaf in tree.leaves() {
        let p = position(&h, leaf).unwrap();
        assert_eq!(r[(p.base as u64 * 729 + p.row * 27 + p.col) as usize], 1);
    }
    let back = zone::from_raster(h, d, &r, ctx.clone()).unwrap();
    assert_eq!(back.root(), tree.root());
    assert_eq!(back.cells(), tree.cells(), "rebuilt tree is canonical");
    assert!(zone::from_raster(h, d, &r[1..], ctx).is_err());
}

#[test]
fn raster_index_and_raster_cells_are_inverse() {
    for level in 0..=3 {
        let cells = zone::raster_cells(&SPACE, level).unwrap();
        assert_eq!(cells.len() as u64, 6 * 9u64.pow(level));
        assert_eq!(cells.iter().collect::<BTreeSet<_>>().len(), cells.len());
        let idx = zone::raster_index(&SPACE, &cells).unwrap();
        assert!(idx.iter().enumerate().all(|(i, &x)| x == i as u64));
    }
    let mixed = [suid_to_cid("Q4").unwrap(), suid_to_cid("Q45").unwrap()];
    assert!(zone::raster_index(&SPACE, &mixed).is_err(), "cids of two levels have no common raster");
    assert!(zone::raster_index(&SPACE, &[]).unwrap().is_empty());
}

#[test]
fn a_non_square_aperture_has_no_rows() {
    let h = Hierarchy::new(8, 6).unwrap();
    assert!(position(&h, h.cid(&[0, 1]).unwrap()).is_err());
}
