//! Polygon → canonical cell set parity against the nucleus-rule fixtures generated from
//! rhealpixdggs-py (set 1). Exact set equality.

use burin_core::hierarchy::cid_to_suid;
use burin_core::polyfill::{cover, geometry_from_geojson, polyfill};
use burin_core::profile::Profile;
use serde_json::Value;
use std::path::PathBuf;
use std::sync::Arc;

fn fixture(name: &str) -> Value {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures").join(name);
    let s = std::fs::read_to_string(&p).unwrap_or_else(|_| panic!("missing fixture {}", p.display()));
    serde_json::from_str(&s).unwrap()
}

fn run(name: &str, profile: &Profile) {
    let fx = fixture(&format!("polyfill_{name}.json"));
    assert_eq!(fx["profile"]["lon_0"].as_f64().unwrap(), profile.lon_0());
    let mut n = 0;
    for case in fx["cases"].as_array().unwrap() {
        let cname = case["name"].as_str().unwrap();
        let res = case["resolution"].as_u64().unwrap() as u32;
        let geom = geometry_from_geojson(&case["geometry"]).unwrap();
        let ctx = Arc::new(profile.ctx(res).unwrap());
        let tree = cover(&geom, res, profile, ctx).unwrap();
        let mut got: Vec<String> = tree.cells().iter().map(|&c| cid_to_suid(c).unwrap()).collect();
        got.sort_by(|a, b| (a.len(), a).cmp(&(b.len(), b)));
        let want: Vec<String> = case["canonical"].as_array().unwrap().iter().map(|v| v.as_str().unwrap().to_string()).collect();
        if got != want {
            let extra: Vec<_> = got.iter().filter(|s| !want.contains(s)).take(10).collect();
            let missing: Vec<_> = want.iter().filter(|s| !got.contains(s)).take(10).collect();
            panic!("{name}/{cname} @{res}: {} vs {} canonical cells; extra {extra:?} missing {missing:?}", got.len(), want.len());
        }
        assert_eq!(tree.leaf_count(), case["leaves"].as_u64().unwrap(), "{name}/{cname} leaf count");
        if res <= 8 {
            assert_eq!(polyfill(&geom, res, profile).unwrap().len() as u64, tree.leaf_count(), "{name}/{cname} polyfill");
        }
        n += 1;
    }
    assert!(n >= 20, "{name}: only {n} cases");
}

#[test]
fn polyfill_parity_ogc() {
    run("ogc", &Profile::ogc());
}

#[test]
fn polyfill_parity_burin1() {
    run("burin1", &Profile::burin_1());
}

#[test]
fn readme_style_fingerprint_is_stable_across_encodings() {
    // the same square as a Polygon and as a one-part MultiPolygon with a redundant vertex
    let p = Profile::ogc();
    let a = geometry_from_geojson(&serde_json::json!({"type":"Polygon","coordinates":[[[-0.13,51.5],[-0.10,51.5],[-0.10,51.52],[-0.13,51.52],[-0.13,51.5]]]})).unwrap();
    let b = geometry_from_geojson(&serde_json::json!({"type":"MultiPolygon","coordinates":[[[[-0.13,51.5],[-0.115,51.5],[-0.10,51.5],[-0.10,51.52],[-0.13,51.52],[-0.13,51.5]]]]})).unwrap();
    let fa = burin_core::polyfill::fingerprint_polygon(&a, 10, &p).unwrap();
    let fb = burin_core::polyfill::fingerprint_polygon(&b, 10, &p).unwrap();
    assert_eq!(fa, fb);
    assert_ne!(fa, burin_core::polyfill::fingerprint_polygon(&a, 9, &p).unwrap(), "resolution is part of the identity");
    assert_ne!(fa, burin_core::polyfill::fingerprint_polygon(&a, 10, &Profile::burin_1()).unwrap(), "profile is part of the identity");
}
