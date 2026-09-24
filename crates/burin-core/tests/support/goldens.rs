//! The point and cell sets burin's own answers are frozen for (`points_burin.json`), and the
//! comparison with the reference fixtures that lists every disagreement exactly
//! (`reference_disagreements.json`). Shared by `examples/freeze_points.rs` and the tests.

use burin_core::geo::Grid;
use burin_core::hierarchy::{cid_to_suid, suid_to_cid, Cid, SPACE};
use burin_core::profile::Profile;
use burin_core::zone::cell_at;
use serde_json::{json, Value};
use std::path::Path;

/// The deepest level the point goldens record; every coarser level is its ancestor.
pub const DEEPEST: u32 = 15;

pub fn profiles() -> [(&'static str, Profile); 3] {
    [
        ("ogc", Profile::ogc()),
        ("burin1", Profile::burin_1()),
        ("ns1ss2", Profile { lon_0_udeg: 0, north_square: 1, south_square: 2, ..Profile::default() }),
    ]
}

pub fn bits(x: f64) -> String {
    format!("{:016x}", x.to_bits())
}

pub fn unbits(v: &Value) -> f64 {
    f64::from_bits(u64::from_str_radix(v.as_str().expect("hex bits"), 16).expect("hex bits"))
}

/// Distance in units in the last place, on the ordered line of finite doubles.
pub fn ulps(a: f64, b: f64) -> u64 {
    let key = |x: f64| {
        let i = x.to_bits() as i64;
        if i < 0 { i64::MIN - i } else { i }
    };
    key(a).abs_diff(key(b))
}

pub fn lookup(grid: &Grid, lon: f64, lat: f64, level: u32) -> Option<Cid> {
    grid.cell_from_lonlat(lon, lat, level).map(|(b, r, c)| cell_at(&SPACE, b, level, r, c).unwrap())
}

fn xorshift(state: &mut u64) -> f64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    (*state >> 11) as f64 / (1u64 << 53) as f64
}

fn random_points(n: usize, seed: u64) -> Vec<(f64, f64)> {
    let mut s = seed;
    (0..n).map(|_| (xorshift(&mut s) * 360.0 - 180.0, libm::asin(xorshift(&mut s) * 2.0 - 1.0).to_degrees())).collect()
}

fn seams(profile: &Profile) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for k in -4..=4 {
        let lon = ((profile.lon_0() + 45.0 * k as f64 + 180.0).rem_euclid(360.0)) - 180.0;
        for lat in [0.0, 30.0, -30.0, 41.9375, -41.9375, 41.94, -41.94, 89.9, -89.9, 90.0, -90.0] {
            out.push((lon, lat));
        }
    }
    out
}

/// Points on the hard places, each with its one-ulp neighbours: the corners and edge midpoints of
/// every cell at levels 1 and 2 (projected back to lon/lat), the seams between base cells at the
/// band edge and the poles, plus random points.
pub fn point_set(grid: &Grid, profile: &Profile) -> Vec<(f64, f64)> {
    let mut base = Vec::new();
    for level in 1..=2u32 {
        for c in (0..6u64).flat_map(|b| SPACE.descendants(9 + b, level)) {
            let path = SPACE.path(c).unwrap();
            let (x, y) = grid.ul_vertex(&path);
            let w = grid.cell_width(level);
            for (dx, dy) in [(0.0, 0.0), (w / 2.0, 0.0), (0.0, w / 2.0)] {
                if let Ok(p) = grid.inverse(x + dx, y - dy, None) {
                    base.push(p);
                }
            }
        }
    }
    base.extend(seams(profile));
    let mut out = Vec::with_capacity(base.len() * 5 + 400);
    for (lon, lat) in base {
        out.extend([(lon, lat), (lon.next_up(), lat), (lon.next_down(), lat), (lon, lat.next_up()), (lon, lat.next_down())]);
    }
    out.extend(random_points(400, 0x5eed_0001));
    out
}

/// Cells whose nucleus is frozen: every cell at levels 0-2, every ninth at level 3, and the cells
/// at levels 9-15 in each base cell's corners, centre and edges' midpoints.
pub fn nucleus_cells() -> Vec<Cid> {
    let mut out: Vec<Cid> = (0..=2).flat_map(|l| (0..6u64).flat_map(move |b| SPACE.descendants(9 + b, l))).collect();
    out.extend((0..6u64).flat_map(|b| SPACE.descendants(9 + b, 3)).step_by(9));
    for level in 9..=15usize {
        for b in 0..6u32 {
            for digit in [0u32, 2, 4, 6, 8] {
                out.push(SPACE.cid(&std::iter::once(b).chain(std::iter::repeat_n(digit, level)).collect::<Vec<_>>()).unwrap());
            }
        }
    }
    out
}

/// burin's own answers, bit for bit.
pub fn self_goldens() -> Value {
    let mut per_profile = serde_json::Map::new();
    for (name, profile) in profiles() {
        let grid = profile.grid();
        let points: Vec<Value> = point_set(&grid, &profile)
            .into_iter()
            .map(|(lon, lat)| json!([bits(lon), bits(lat), lookup(&grid, lon, lat, DEEPEST).unwrap_or(0)]))
            .collect();
        let forward: Vec<Value> = random_points(400, 0x5eed_0001)
            .into_iter()
            .chain(seams(&profile))
            .map(|(lon, lat)| {
                let (x, y) = grid.forward(lon, lat, None);
                json!([bits(lon), bits(lat), bits(x), bits(y)])
            })
            .collect();
        let nuclei: Vec<Value> = nucleus_cells()
            .into_iter()
            .map(|c| {
                let (lon, lat) = grid.nucleus(&SPACE.path(c).unwrap()).unwrap();
                json!([c, bits(lon), bits(lat)])
            })
            .collect();
        per_profile.insert(
            name.to_string(),
            json!({"lon_0_udeg": profile.lon_0_udeg, "north_square": profile.north_square, "south_square": profile.south_square,
                   "points": points, "forward": forward, "nuclei": nuclei}),
        );
    }
    json!({"deepest": DEEPEST, "profiles": per_profile})
}

fn read(dir: &Path, name: &str) -> Value {
    serde_json::from_str(&std::fs::read_to_string(dir.join(name)).unwrap_or_else(|_| panic!("missing fixture {name}"))).unwrap()
}

/// Every place burin's answer differs from the reference fixtures (rhealpixdggs-py 0.8.6): point
/// lookups, forward projections, and nuclei, with the distance in ulps.
pub fn reference_disagreements(dir: &Path) -> Value {
    let mut out = serde_json::Map::new();
    for (name, profile) in profiles() {
        let grid = profile.grid();
        let points_fx = read(dir, &format!("points_{name}.json"));
        let (mut lookups, mut forward) = (Vec::new(), Vec::new());
        for row in points_fx["lonlat"].as_array().unwrap() {
            let (lon, lat, level) = (unbits(&row[0]), unbits(&row[1]), row[2].as_u64().unwrap() as u32);
            let want = row[3].as_str().unwrap();
            let got = lookup(&grid, lon, lat, level).map(|c| cid_to_suid(c).unwrap()).unwrap_or_default();
            if got != want {
                lookups.push(json!([row[0], row[1], level, want, got]));
            }
            let (x, y) = grid.forward(lon, lat, None);
            let (rx, ry) = (unbits(&row[4]), unbits(&row[5]));
            if level == 0 && (x.to_bits() != rx.to_bits() || y.to_bits() != ry.to_bits()) {
                forward.push(json!([row[0], row[1], ulps(x, rx), ulps(y, ry)]));
            }
        }
        let cells_fx = read(dir, &format!("cells_{name}.json"));
        let mut nuclei = Vec::new();
        for c in cells_fx["cells"].as_array().unwrap().iter().chain(cells_fx["seam_nuclei"].as_array().unwrap()) {
            let suid = c["suid"].as_str().unwrap();
            let (lon, lat) = grid.nucleus(&SPACE.path(suid_to_cid(suid).unwrap()).unwrap()).unwrap();
            let (rlon, rlat) = (c["nucleus"][0].as_f64().unwrap(), c["nucleus"][1].as_f64().unwrap());
            if lon.to_bits() != rlon.to_bits() || lat.to_bits() != rlat.to_bits() {
                nuclei.push(json!([suid, ulps(lon, rlon), ulps(lat, rlat)]));
            }
        }
        nuclei.sort_by(|a, b| a[0].as_str().cmp(&b[0].as_str()));
        nuclei.dedup();
        out.insert(name.to_string(), json!({"lookups": lookups, "forward": forward, "nuclei": nuclei}));
    }
    Value::Object(out)
}
