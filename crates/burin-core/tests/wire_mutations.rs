//! Wire records under every single edit of their JSON: each value replaced, removed, re-typed,
//! pushed out of range or re-spelled. A verifier must never panic, hang or allocate what a record
//! claims, and a mutant that still verifies must state something true: an opening of a cell is
//! the tree's own opening of that cell (swapping equal siblings relabels a proof to another cell
//! that really is covered), and a set-operation proof states the original operation and roots.

use burin_core::hierarchy::{suid_to_cid, Hierarchy, SPACE};
use burin_core::opening::{open_path, OpeningRecord};
use burin_core::profile::Profile;
use burin_core::setops::{prove, Op, SetOpProof};
use burin_core::tree::Tree;
use burin_core::zone_data::{from_dggs_json, to_dggs_json, Presence};
use serde_json::{json, Value};
use std::sync::Arc;

fn tree(d: u32, cells: &[&str]) -> Tree {
    let ctx = Arc::new(Profile::ogc().ctx(d).unwrap());
    Tree::from_cells(SPACE, d, cells.iter().map(|s| suid_to_cid(s).unwrap()), ctx).unwrap()
}

#[derive(Clone, Debug)]
enum Seg {
    Key(String),
    At(usize),
}

fn paths(v: &Value, here: &mut Vec<Seg>, out: &mut Vec<Vec<Seg>>) {
    out.push(here.clone());
    match v {
        Value::Object(m) => m.iter().for_each(|(k, x)| {
            here.push(Seg::Key(k.clone()));
            paths(x, here, out);
            here.pop();
        }),
        Value::Array(a) => a.iter().enumerate().for_each(|(i, x)| {
            here.push(Seg::At(i));
            paths(x, here, out);
            here.pop();
        }),
        _ => {}
    }
}

fn at<'a>(v: &'a mut Value, path: &[Seg]) -> &'a mut Value {
    path.iter().fold(v, |v, s| match s {
        Seg::Key(k) => &mut v[k.as_str()],
        Seg::At(i) => &mut v[*i],
    })
}

/// Every value the node at a path can be replaced by, plus `None` for removing it.
fn replacements(v: &Value) -> Vec<Option<Value>> {
    let mut out: Vec<Option<Value>> = vec![None, Some(Value::Null), Some(json!(true)), Some(json!(-1)), Some(json!(1.5)), Some(json!("x"))];
    out.extend([0u64, 1, 2, 4096, 4097, u32::MAX as u64, 1 << 32, u64::MAX].map(|n| Some(json!(n))));
    match v {
        Value::Number(n) => out.extend(n.as_u64().map(|n| Some(json!(n + (1 << 32))))),
        Value::String(s) if s.len() == 64 => {
            let mut b = s.clone().into_bytes();
            b[63] = if b[63] == b'0' { b'1' } else { b'0' };
            out.push(Some(Value::String(String::from_utf8(b).unwrap())));
            out.push(Some(Value::String(s.to_uppercase())));
            out.push(Some(Value::String(s[..62].to_string())));
        }
        Value::Array(a) if !a.is_empty() => {
            out.push(Some(Value::Array(vec![])));
            out.push(Some(Value::Array(a[1..].to_vec())));
            out.push(Some(Value::Array([a.clone(), vec![a[a.len() - 1].clone()]].concat())));
            out.push(Some(Value::Array(a.iter().rev().cloned().collect())));
        }
        _ => {}
    }
    out
}

fn mutants(original: &Value) -> Vec<Value> {
    let mut ps = Vec::new();
    paths(original, &mut Vec::new(), &mut ps);
    let mut out = Vec::new();
    for p in ps.iter().filter(|p| !p.is_empty()) {
        let target = at(&mut original.clone(), p).clone();
        for r in replacements(&target) {
            let mut m = original.clone();
            let (parent, last) = (at(&mut m, &p[..p.len() - 1]), &p[p.len() - 1]);
            match (r, last) {
                (Some(x), Seg::Key(k)) => parent[k.as_str()] = x,
                (Some(x), Seg::At(i)) => parent[*i] = x,
                (None, Seg::Key(k)) => {
                    parent.as_object_mut().unwrap().remove(k);
                }
                (None, Seg::At(i)) => {
                    parent.as_array_mut().unwrap().remove(*i);
                }
            }
            if m != *original {
                out.push(m);
            }
        }
    }
    out
}

fn check_all<R>(original: Value, parse: impl Fn(&Value) -> burin_core::Result<R>, verify: impl Fn(&R) -> bool, true_of: impl Fn(&R) -> bool) -> usize {
    assert!(verify(&parse(&original).unwrap()), "the unmutated record verifies");
    let ms = mutants(&original);
    for m in &ms {
        if let Ok(r) = parse(m) {
            assert!(!verify(&r) || true_of(&r), "a mutant verifies a false statement:\n{m}");
        }
    }
    ms.len()
}

#[test]
fn no_edit_of_an_opening_record_verifies() {
    let t = tree(3, &["Q4", "Q51", "N000", "S8"]);
    let mut n = 0;
    for cell in ["Q412", "Q518", "N000", "S820", "P000"] {
        let op = open_path(&t, suid_to_cid(cell).unwrap()).unwrap().unwrap();
        let own = |r: &OpeningRecord| r.cid().and_then(|c| open_path(&t, c).unwrap()).is_some_and(|o| OpeningRecord::new(&t, o) == *r);
        n += check_all(OpeningRecord::new(&t, op).to_json(), OpeningRecord::from_json, OpeningRecord::verify, own);
    }
    assert!(n > 2000, "only {n} mutants");
}

#[test]
fn no_edit_of_a_set_operation_proof_verifies() {
    let (a, b) = (tree(3, &["Q4", "Q51", "N0"]), tree(3, &["Q45", "Q5", "S"]));
    let mut n = 0;
    for op in [Op::Union, Op::Intersect, Op::Difference] {
        let p = prove(&a, &b, op).unwrap();
        let same = |r: &SetOpProof| (r.op, r.root_a, r.root_b, r.root_c) == (p.op, p.root_a, p.root_b, p.root_c);
        n += check_all(p.to_json(), SetOpProof::from_json, SetOpProof::verify, same);
    }
    assert!(n > 2000, "only {n} mutants");
}

#[test]
fn record_sizes_are_read_exactly_or_refused() {
    let t = tree(2, &["Q4"]);
    let rec = OpeningRecord::new(&t, open_path(&t, suid_to_cid("Q41").unwrap()).unwrap().unwrap()).to_json();
    let with = |k: &str, v: Value| {
        let mut r = rec.clone();
        r[k] = v;
        r
    };
    // a number past u32 is not read as its low bits: each field has one spelling
    for k in ["A", "B", "D"] {
        let n = rec[k].as_u64().unwrap();
        assert!(OpeningRecord::from_json(&with(k, json!(n + (1 << 32)))).is_err(), "{k} + 2^32");
    }
    // a hierarchy the verifier could not build a ladder for is refused before any hashing
    for a in [0u64, 1, 4097, u32::MAX as u64] {
        assert!(OpeningRecord::from_json(&with("A", json!(a))).is_err(), "A = {a}");
    }
    assert!(OpeningRecord::from_json(&with("D", json!(SPACE.max_level() + 1))).is_err());
    assert!(Hierarchy::new(4096, 4096 * 4095).is_ok_and(|h| h.max_level() >= 1));
    // a record built by hand with a degenerate hierarchy answers rather than looping
    let mut r = OpeningRecord::from_json(&rec).unwrap();
    for a in [0, 1] {
        r.a = a;
        assert!(!r.verify());
        let _ = r.cid();
    }
}

#[test]
fn zone_data_depths_are_read_exactly_or_refused() {
    let p = Profile::ogc();
    let t = tree(3, &["Q41", "Q418"]);
    let doc = to_dggs_json(&t, &p, suid_to_cid("Q4").unwrap(), "f").unwrap();
    assert_eq!(from_dggs_json(&doc, &p, None, None, Presence::NonNull).unwrap().root(), t.root());
    let mut wide = doc.clone();
    wide["depths"] = json!([2 + (1u64 << 32)]);
    wide["values"]["f"][0]["depth"] = json!(2 + (1u64 << 32));
    assert!(from_dggs_json(&wide, &p, None, None, Presence::NonNull).is_err());
    // a depth that fits the sub-zone limit but not below this zone
    let deep_zone = format!("Q{}", "4".repeat(SPACE.max_level() as usize));
    let mut deep = json!({"dggrs": doc["dggrs"], "zoneId": deep_zone, "depths": [1],
                          "values": {"f": [{"depth": 1, "data": vec![json!(1); 9]}]}});
    assert!(from_dggs_json(&deep, &p, None, None, Presence::NonNull).is_err());
    deep["zoneId"] = json!(format!("Q{}", "4".repeat(SPACE.max_level() as usize - 1)));
    assert!(from_dggs_json(&deep, &p, None, None, Presence::NonNull).is_ok());
}

#[test]
fn no_edit_of_zone_data_panics() {
    let p = Profile::ogc();
    let t = tree(2, &["Q41", "Q48"]);
    let doc = to_dggs_json(&t, &p, suid_to_cid("Q4").unwrap(), "f").unwrap();
    for m in mutants(&doc) {
        for presence in [Presence::NonNull, Presence::NonZero] {
            let _ = from_dggs_json(&m, &p, None, None, presence);
        }
    }
}
