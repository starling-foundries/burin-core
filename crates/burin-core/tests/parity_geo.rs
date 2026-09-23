//! Parity of the projection and cell layout against rhealpixdggs-py fixtures (set 1).

use burin_core::geo::ellipsoid::{auth_lat, auth_rad};
use burin_core::geo::healpix::*;
use burin_core::geo::rhealpix::*;
use burin_core::geo::{Ellipsoid, Grid, Region, Shape};
use burin_core::hierarchy::suid_to_path;
use serde_json::Value;
use std::path::PathBuf;

fn fixture(name: &str) -> Value {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures").join(name);
    let s = std::fs::read_to_string(&p).unwrap_or_else(|_| panic!("missing fixture {}", p.display()));
    serde_json::from_str(&s).unwrap()
}

fn g(v: &Value, k: &str) -> f64 {
    v[k].as_f64().unwrap_or_else(|| panic!("field {k} in {v}"))
}

/// Relative closeness with an absolute floor of `rel` (angles are O(1), so a cancellation to
/// 1e-16 is a last-ulp effect of an O(1) operand, not a relative error of 100%).
fn close(a: f64, b: f64, rel: f64, what: &str) {
    if a == b || (a.is_nan() && b.is_nan()) {
        return;
    }
    let scale = a.abs().max(b.abs()).max(1.0);
    let err = (a - b).abs() / scale;
    assert!(err <= rel, "{what}: {a:?} vs {b:?} (rel {err:.3e} > {rel:.1e})");
}

/// Strict relative closeness for planar metres.
fn close_rel(a: f64, b: f64, rel: f64, what: &str) {
    if a == b {
        return;
    }
    let err = (a - b).abs() / a.abs().max(b.abs());
    assert!(err <= rel, "{what}: {a:?} vs {b:?} (rel {err:.3e} > {rel:.1e})");
}

/// Longitudes compare on the circle: ±180 and the pole's −180 convention alias.
fn lon_close(a: f64, b: f64, rel: f64, what: &str) {
    let d = ((a - b).abs() + 180.0).rem_euclid(360.0) - 180.0;
    assert!(d.abs() <= rel * 360.0, "{what}: lon {a:?} vs {b:?}");
}

fn region_of(s: &str) -> Region {
    match s {
        "equatorial" => Region::Equatorial,
        "north_polar" => Region::NorthPolar,
        "south_polar" => Region::SouthPolar,
        _ => panic!("region {s}"),
    }
}

#[test]
fn auth_lat_parity() {
    let fx = fixture("auth_lat.json");
    let e = g(&fx, "e_wgs84");
    close(Ellipsoid::wgs84().e, e, 1e-15, "e");
    close(auth_rad(6378137.0, e), g(&fx, "auth_rad_wgs84"), 1e-15, "R_A");
    for r in fx["rows"].as_array().unwrap() {
        let got = auth_lat(g(r, "phi"), g(r, "e"), r["inverse"].as_bool().unwrap());
        close(got, g(r, "value"), 1e-14, &format!("auth_lat {r}"));
    }
}

#[test]
fn healpix_parity() {
    let fx = fixture("healpix.json");
    for r in fx["forward"].as_array().unwrap() {
        let e = g(r, "e");
        let (x, y) = if e == 0.0 { healpix_sphere(g(r, "lam"), g(r, "phi")) } else { healpix_ellipsoid(g(r, "lam"), g(r, "phi"), e) };
        close(x, g(r, "x"), 1e-12, &format!("fwd x {r}"));
        close(y, g(r, "y"), 1e-12, &format!("fwd y {r}"));
    }
    for r in fx["inverse"].as_array().unwrap() {
        let e = g(r, "e");
        let (lam, phi) = if e == 0.0 { healpix_sphere_inverse(g(r, "x"), g(r, "y")) } else { healpix_ellipsoid_inverse(g(r, "x"), g(r, "y"), e) }.unwrap();
        close(lam, g(r, "lam"), 1e-12, &format!("inv lam {r}"));
        close(phi, g(r, "phi"), 1e-12, &format!("inv phi {r}"));
    }
    for r in fx["poles"].as_array().unwrap() {
        let (lam, phi) = healpix_sphere_inverse(g(r, "x"), g(r, "y")).unwrap();
        if g(r, "y").abs() >= core::f64::consts::FRAC_PI_2 {
            // at the poles the convention (lon = -π) is exact, not approximate
            assert_eq!(lam, g(r, "lam"), "pole lam {r}");
            assert_eq!(phi, g(r, "phi"), "pole phi {r}");
        } else {
            close(lam, g(r, "lam"), 1e-12, &format!("near-pole lam {r}"));
            close(phi, g(r, "phi"), 1e-12, &format!("near-pole phi {r}"));
        }
    }
    for r in fx["image"].as_array().unwrap() {
        assert_eq!(in_healpix_image(g(r, "x"), g(r, "y")), r["inside"].as_bool().unwrap(), "image {r}");
    }
}

#[test]
fn rhealpix_parity() {
    let fx = fixture("rhealpix.json");
    for r in fx["triangle"].as_array().unwrap() {
        let (t, region) = triangle(g(r, "x"), g(r, "y"), r["ns"].as_i64().unwrap() as i32, r["ss"].as_i64().unwrap() as i32, r["inverse"].as_bool().unwrap());
        let want = r["t"].as_i64().map(|v| v as i32);
        assert_eq!(t, want, "triangle {r}");
        assert_eq!(region, region_of(r["region"].as_str().unwrap()), "region {r}");
    }
    for r in fx["combine"].as_array().unwrap() {
        let (x, y) = combine_triangles(g(r, "x"), g(r, "y"), r["ns"].as_i64().unwrap() as i32, r["ss"].as_i64().unwrap() as i32, r["inverse"].as_bool().unwrap());
        close(x, g(r, "rx"), 1e-13, &format!("combine x {r}"));
        close(y, g(r, "ry"), 1e-13, &format!("combine y {r}"));
    }
    for r in fx["image"].as_array().unwrap() {
        assert_eq!(in_rhealpix_image(g(r, "x"), g(r, "y"), r["ns"].as_i64().unwrap() as i32, r["ss"].as_i64().unwrap() as i32), r["inside"].as_bool().unwrap(), "image {r}");
    }
    let e = g(&fx, "e");
    for r in fx["ellipsoid"].as_array().unwrap() {
        let (ns, ss) = (r["ns"].as_i64().unwrap() as i32, r["ss"].as_i64().unwrap() as i32);
        let (x, y) = rhealpix_ellipsoid(g(r, "lam"), g(r, "phi"), e, ns, ss, None);
        close(x, g(r, "x"), 1e-12, &format!("ell x {r}"));
        close(y, g(r, "y"), 1e-12, &format!("ell y {r}"));
        let (lam, phi) = rhealpix_ellipsoid_inverse(g(r, "x"), g(r, "y"), e, ns, ss, None).unwrap();
        close(lam, g(r, "lam_back"), 1e-12, &format!("ell lam {r}"));
        close(phi, g(r, "phi_back"), 1e-12, &format!("ell phi {r}"));
    }
}

fn grid_for(p: &Value) -> Grid {
    Grid::new(Ellipsoid::wgs84(), g(p, "lon_0"), p["ns"].as_i64().unwrap() as i32, p["ss"].as_i64().unwrap() as i32, 3)
}

fn cells_parity(name: &str) {
    let fx = fixture(&format!("cells_{name}.json"));
    let grid = grid_for(&fx["profile"]);
    for (r, w) in fx["cell_width"].as_array().unwrap().iter().enumerate() {
        close_rel(grid.cell_width(r as u32), w.as_f64().unwrap(), 1e-15, &format!("width {r}"));
    }
    for (r, a) in fx["cell_area"].as_array().unwrap().iter().enumerate() {
        close_rel(grid.cell_area(r as u32), a.as_f64().unwrap(), 1e-15, &format!("area {r}"));
    }
    let mut n = 0;
    for c in fx["cells"].as_array().unwrap() {
        let suid = c["suid"].as_str().unwrap();
        let path = suid_to_path(suid).unwrap();
        let ul = grid.ul_vertex(&path);
        close_rel(ul.0, c["ul"][0].as_f64().unwrap(), 1e-15, &format!("{suid} ul.x"));
        close_rel(ul.1, c["ul"][1].as_f64().unwrap(), 1e-15, &format!("{suid} ul.y"));
        let np = grid.nucleus_planar(&path);
        close_rel(np.0, c["nucleus_planar"][0].as_f64().unwrap(), 1e-15, &format!("{suid} nucleus_planar.x"));
        close_rel(np.1, c["nucleus_planar"][1].as_f64().unwrap(), 1e-15, &format!("{suid} nucleus_planar.y"));
        let nu = grid.nucleus(&path).unwrap();
        lon_close(nu.0, c["nucleus"][0].as_f64().unwrap(), 1e-12, &format!("{suid} nucleus.lon"));
        close(nu.1, c["nucleus"][1].as_f64().unwrap(), 1e-12, &format!("{suid} nucleus.lat"));
        let mut b = grid.boundary3(&path).unwrap();
        b.sort_by(|p, q| p.partial_cmp(q).unwrap());
        let want: Vec<(f64, f64)> = c["boundary3"].as_array().unwrap().iter().map(|p| (p[0].as_f64().unwrap(), p[1].as_f64().unwrap())).collect();
        assert_eq!(b.len(), want.len(), "{suid} boundary count");
        for (i, (p, q)) in b.iter().zip(&want).enumerate() {
            lon_close(p.0, q.0, 1e-12, &format!("{suid} boundary[{i}].lon"));
            close(p.1, q.1, 1e-12, &format!("{suid} boundary[{i}].lat"));
        }
        assert_eq!(Grid::region(&path), region_of(c["region"].as_str().unwrap()), "{suid} region");
        let shape = match c["shape"].as_str().unwrap() {
            "quad" => Shape::Quad,
            "cap" => Shape::Cap,
            "dart" => Shape::Dart,
            "skew_quad" => Shape::SkewQuad,
            s => panic!("shape {s}"),
        };
        assert_eq!(grid.shape(&path), shape, "{suid} shape");
        n += 1;
    }
    assert!(n > 700, "too few cells in fixture {name}: {n}");
    // nuclei that fall exactly on the antimeridian seam (lon == -180.0 to the bit) are a tie
    // surface for split polygons; the two implementations must agree on them exactly
    let seam = fx["seam_nuclei"].as_array().unwrap();
    for c in seam {
        let path = suid_to_path(c["suid"].as_str().unwrap()).unwrap();
        let nu = grid.nucleus(&path).unwrap();
        assert_eq!(nu.0, -180.0, "{} seam nucleus lon", c["suid"]);
        close(nu.1, c["nucleus"][1].as_f64().unwrap(), 1e-12, &format!("{} seam nucleus lat", c["suid"]));
    }
    assert!(!seam.is_empty() || name == "ns1ss2", "expected some seam nuclei in {name}");
}

#[test]
fn cells_parity_ogc() {
    cells_parity("ogc");
}

#[test]
fn cells_parity_burin1() {
    cells_parity("burin1");
}

#[test]
fn cells_parity_ns1ss2() {
    cells_parity("ns1ss2");
}
