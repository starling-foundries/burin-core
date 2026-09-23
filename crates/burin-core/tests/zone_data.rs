//! DGGS-JSON zone data: round trips, sub-zone placement, and every refusal.

use burin_core::hierarchy::{suid_to_cid, Cid, SPACE};
use burin_core::profile::Profile;
use burin_core::setops::intersect;
use burin_core::tree::Tree;
use burin_core::zone::subzone_index;
use burin_core::zone_data::{dggrs_uri, from_dggs_json, to_dggs_json, Presence, RHEALPIX_DGGRS};
use serde_json::{json, Value};
use std::sync::Arc;

fn tree(d: u32, suids: &[&str]) -> Tree {
    let p = Profile::ogc();
    let cells: Vec<Cid> = suids.iter().map(|s| suid_to_cid(s).unwrap()).collect();
    Tree::from_cells(SPACE, d, cells, Arc::new(p.ctx(d).unwrap())).unwrap()
}

fn read(doc: &Value) -> burin_core::Result<Tree> {
    from_dggs_json(doc, &Profile::ogc(), None, None, Presence::NonNull)
}

#[test]
fn a_zone_round_trips_as_its_share_of_the_coverage() {
    let t = tree(4, &["Q4", "Q51", "Q520", "R0", "N8", "Q3888"]);
    let p = Profile::ogc();
    for zone in ["Q", "Q5", "Q52", "N", "R0", "S"] {
        let z = suid_to_cid(zone).unwrap();
        let doc = to_dggs_json(&t, &p, z, "coverage").unwrap();
        assert_eq!(doc["dggrs"], RHEALPIX_DGGRS);
        assert_eq!(doc["zoneId"], zone);
        let back = read(&doc).unwrap();
        let expected = intersect(&t, &tree(4, &[zone])).unwrap();
        assert_eq!(back.root(), expected.root(), "{zone}");
        assert_eq!(back.cells(), expected.cells(), "{zone}: rebuilt tree is canonical");
    }
}

#[test]
fn values_sit_at_their_scanline_index() {
    let t = tree(3, &["Q452"]);
    let zone = suid_to_cid("Q").unwrap();
    let doc = to_dggs_json(&t, &Profile::ogc(), zone, "coverage").unwrap();
    let data = doc["values"]["coverage"][0]["data"].as_array().unwrap();
    assert_eq!(data.len(), 729);
    let i = subzone_index(&SPACE, zone, suid_to_cid("Q452").unwrap()).unwrap() as usize;
    assert_eq!(data[i], json!(1));
    assert_eq!(data.iter().filter(|v| !v.is_null()).count(), 1);
    let full = to_dggs_json(&tree(2, &["Q"]), &Profile::ogc(), suid_to_cid("Q4").unwrap(), "c").unwrap();
    assert!(full["values"]["c"][0]["data"].as_array().unwrap().iter().all(|v| *v == json!(1)), "an ancestor covers every sub-zone");
}

fn doc(data: Vec<Value>) -> Value {
    let n = data.len();
    json!({"dggrs": RHEALPIX_DGGRS, "zoneId": "Q4", "depths": [1],
           "values": {"v": [{"depth": 1, "shape": {"count": n, "subZones": n}, "data": data}]}})
}

#[test]
fn presence_is_non_null_unless_asked_for_non_zero() {
    let d = doc(vec![json!(0), json!(2.5), Value::Null, json!(-1), json!(0.0), Value::Null, json!(1), json!(0), json!(7)]);
    let any = read(&d).unwrap();
    let nonzero = from_dggs_json(&d, &Profile::ogc(), Some("v"), Some(1), Presence::NonZero).unwrap();
    assert_eq!(any.leaf_count(), 7);
    assert_eq!(nonzero.leaf_count(), 4);
    assert!(any.covers(suid_to_cid("Q40").unwrap()).unwrap() && !nonzero.covers(suid_to_cid("Q40").unwrap()).unwrap());
}

#[test]
fn a_document_for_another_grid_is_refused() {
    let d = doc(vec![json!(1); 9]);
    let mut other = d.clone();
    other["dggrs"] = json!("https://www.opengis.net/def/dggrs/OGC/1.0/ISEA3H");
    assert!(read(&other).is_err());
    let mut missing = d.clone();
    missing.as_object_mut().unwrap().remove("dggrs");
    assert!(read(&missing).is_err());
    assert!(from_dggs_json(&d, &Profile::burin_1(), None, None, Presence::NonNull).is_err(), "lon_0 = 0 is not the registered grid");
    assert!(to_dggs_json(&tree(1, &["Q4"]), &Profile::burin_1(), suid_to_cid("Q").unwrap(), "c").is_err());
    assert_eq!(dggrs_uri(&Profile::ogc().with_hash("sha256")), Some(RHEALPIX_DGGRS));
    assert_eq!(dggrs_uri(&Profile { north_square: 1, ..Profile::ogc() }), None);
}

#[test]
fn malformed_zone_data_is_refused() {
    assert!(read(&doc(vec![json!(1); 8])).is_err(), "too few values");
    assert!(read(&doc(vec![json!(1); 10])).is_err(), "too many values");
    assert!(read(&doc(vec![json!("1"); 9])).is_err(), "strings are not zone data values");
    let mut shape = doc(vec![json!(1); 9]);
    shape["values"]["v"][0]["shape"]["subZones"] = json!(8);
    assert!(read(&shape).is_err());
    let mut dims = doc(vec![json!(1); 9]);
    dims["dimensions"] = json!([{"name": "time", "interval": ["2020-01-01", "2020-12-31"], "grid": {"cellsCount": 2}}]);
    assert!(read(&dims).is_err());
    let mut zone = doc(vec![json!(1); 9]);
    zone["zoneId"] = json!("X4");
    assert!(read(&zone).is_err());
}

#[test]
fn several_fields_or_depths_must_be_named() {
    let mut d = doc(vec![json!(1); 9]);
    d["values"]["w"] = d["values"]["v"].clone();
    assert!(read(&d).is_err());
    assert_eq!(from_dggs_json(&d, &Profile::ogc(), Some("w"), None, Presence::NonNull).unwrap().leaf_count(), 9);
    assert!(from_dggs_json(&d, &Profile::ogc(), Some("x"), None, Presence::NonNull).is_err());
    let mut deep = doc(vec![json!(1); 9]);
    let two: Vec<Value> = (0..81).map(|i| if i % 2 == 0 { json!(1) } else { Value::Null }).collect();
    deep["depths"] = json!([1, 2]);
    deep["values"]["v"].as_array_mut().unwrap().push(json!({"depth": 2, "shape": {"count": 81, "subZones": 81}, "data": two}));
    assert!(read(&deep).is_err());
    let t2 = from_dggs_json(&deep, &Profile::ogc(), None, Some(2), Presence::NonNull).unwrap();
    assert_eq!((t2.d, t2.leaf_count()), (3, 41));
}
