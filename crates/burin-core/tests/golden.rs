//! Fixture set 2: the frozen SHA-256 golden vectors (roots.json) must reproduce exactly.

use burin_core::hash::hex;
use burin_core::hierarchy::{suid_to_cid, SPACE};
use burin_core::opening::OpeningRecord;
use burin_core::polyfill::fingerprint_geojson;
use burin_core::profile::Profile;
use burin_core::setops::SetOpProof;
use burin_core::tree::Tree;
use serde_json::Value;
use std::path::PathBuf;
use std::sync::Arc;

fn fixture() -> Value {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/roots.json");
    serde_json::from_str(&std::fs::read_to_string(p).expect("run `cargo run --example freeze` first")).unwrap()
}

#[test]
fn ladders_and_commit_example() {
    let fx = fixture();
    let ctx = Profile::ogc().ctx(15).unwrap();
    for l in fx["ladders"].as_array().unwrap() {
        let d = l["d"].as_u64().unwrap() as u32;
        assert_eq!(hex(&ctx.ladders.empty(d)), l["empty"].as_str().unwrap());
        assert_eq!(hex(&ctx.ladders.full(d)), l["full"].as_str().unwrap());
    }
    let ce = &fx["commit_example"];
    let cells = ce["cells"].as_array().unwrap().iter().map(|s| suid_to_cid(s.as_str().unwrap()).unwrap());
    let t = Tree::from_cells(SPACE, ce["depth"].as_u64().unwrap() as u32, cells, Arc::new(Profile::ogc().ctx(4).unwrap())).unwrap();
    assert_eq!(t.root_hex(), ce["root"].as_str().unwrap());
    assert_eq!(Profile::ogc().id_hex(), fx["profile_ids"]["ogc"].as_str().unwrap());
    assert_eq!(Profile::burin_1().id_hex(), fx["profile_ids"]["burin1"].as_str().unwrap());
}

#[test]
fn polygon_roots() {
    let fx = fixture();
    for p in fx["polygons"].as_array().unwrap() {
        let profile = if p["profile"] == "ogc" { Profile::ogc() } else { Profile::burin_1() };
        let root = fingerprint_geojson(&p["geometry"], p["resolution"].as_u64().unwrap() as u32, &profile).unwrap();
        assert_eq!(hex(&root), p["root"].as_str().unwrap(), "{}", p["name"]);
    }
}

#[test]
fn openings_and_transcripts_verify_and_mutations_fail() {
    let fx = fixture();
    for o in fx["openings"].as_array().unwrap() {
        let rec = OpeningRecord::from_json(&o["record"]).unwrap();
        assert!(rec.verify(), "{}", o["case"]);
        assert_eq!(rec.cid(), Some(o["cid"].as_u64().unwrap()));
        let mut bad = o["record"].clone();
        bad["D"] = Value::from(rec.d + 1);
        assert!(!OpeningRecord::from_json(&bad).unwrap().verify());
        bad = o["record"].clone();
        bad["hash"] = Value::from("sha3");
        assert!(!OpeningRecord::from_json(&bad).unwrap().verify());
    }
    for t in fx["transcripts"].as_array().unwrap() {
        let p = SetOpProof::from_json(&t["proof"]).unwrap();
        assert!(p.verify(), "{} {} {}", t["a"], t["b"], t["op"]);
        let mut bad = t["proof"].clone();
        bad["op"] = Value::from(if t["op"] == "union" { "intersect" } else { "union" });
        assert!(!SetOpProof::from_json(&bad).unwrap().verify(), "the op is bound");
    }
}
