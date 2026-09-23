//! Freeze fixture set 2: the golden SHA-256 roots, openings and transcripts this crate computes
//! for the set-1 polygon fixtures. Run `cargo run -p burin-core --example freeze`; the result is
//! committed and checked by `tests/golden.rs`. A change here is a wire-format bump.

use burin_core::hash::hex;
use burin_core::hierarchy::{suid_to_cid, SPACE};
use burin_core::opening::{open_path, OpeningRecord};
use burin_core::polyfill::{cover, geometry_from_geojson};
use burin_core::profile::Profile;
use burin_core::setops::{prove, Op};
use burin_core::tree::Tree;
use serde_json::{json, Value};
use std::path::PathBuf;
use std::sync::Arc;

fn main() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let ctx15 = Profile::ogc().ctx(15).unwrap();
    let ladders: Vec<Value> = (0..=15)
        .map(|d| json!({"d": d, "empty": hex(&ctx15.ladders.empty(d)), "full": hex(&ctx15.ladders.full(d))}))
        .collect();

    let commit_cells = ["Q453", "Q454", "R0"];
    let commit = Tree::from_cells(SPACE, 4, commit_cells.iter().map(|s| suid_to_cid(s).unwrap()), Arc::new(Profile::ogc().ctx(4).unwrap())).unwrap();

    let mut polygons = Vec::new();
    let mut openings = Vec::new();
    let mut transcripts = Vec::new();
    let mut trees: Vec<(String, Tree)> = Vec::new();
    for (name, profile) in [("ogc", Profile::ogc()), ("burin1", Profile::burin_1())] {
        let fx: Value = serde_json::from_str(&std::fs::read_to_string(dir.join(format!("polyfill_{name}.json"))).unwrap()).unwrap();
        for case in fx["cases"].as_array().unwrap() {
            let cname = case["name"].as_str().unwrap();
            let res = case["resolution"].as_u64().unwrap() as u32;
            let geom = geometry_from_geojson(&case["geometry"]).unwrap();
            let tree = cover(&geom, res, &profile, Arc::new(profile.ctx(res).unwrap())).unwrap();
            polygons.push(json!({"name": cname, "profile": name, "resolution": res, "geometry": case["geometry"], "root": tree.root_hex(), "canonical_count": tree.cells().len()}));
            // one membership and one non-membership opening per non-empty case
            if let Some(&first) = tree.cells().first() {
                let leaf = SPACE.descendants(first, res - SPACE.level(first)).next().unwrap();
                let op = open_path(&tree, leaf).unwrap().unwrap();
                openings.push(json!({"case": format!("{name}/{cname}"), "cid": leaf, "claim": "full", "record": OpeningRecord::new(&tree, op).to_json()}));
                // an absent leaf: walk the bases for an empty one, else a sibling
                let absent = (0..6u32).map(|b| SPACE.cid(&[b]).unwrap()).find(|&c| !tree.covers(c).unwrap() && !tree.node_at(c).unwrap().0.is_branch());
                if let Some(a) = absent {
                    let op = open_path(&tree, a).unwrap().unwrap();
                    openings.push(json!({"case": format!("{name}/{cname}"), "cid": a, "claim": "empty", "record": OpeningRecord::new(&tree, op).to_json()}));
                }
            }
            if res == 6 {
                trees.push((format!("{name}/{cname}"), tree));
            }
        }
    }
    // transcripts between pairs of res-6 fixtures under the same profile
    for i in 0..trees.len() {
        for j in (i + 1)..trees.len() {
            let (na, a) = &trees[i];
            let (nb, b) = &trees[j];
            if na.split('/').next() != nb.split('/').next() || a.d != b.d {
                continue;
            }
            for op in [Op::Union, Op::Intersect, Op::Difference] {
                let p = prove(a, b, op).unwrap();
                assert!(p.verify());
                transcripts.push(json!({"a": na, "b": nb, "op": op.name(), "root_c": hex(&p.root_c), "size": p.size(), "proof": p.to_json()}));
            }
            if transcripts.len() >= 12 {
                break;
            }
        }
        if transcripts.len() >= 12 {
            break;
        }
    }

    let out = json!({
        "hash": "sha256",
        "ladders": ladders,
        "commit_example": {"cells": commit_cells, "depth": 4, "root": commit.root_hex()},
        "profile_ids": {"ogc": Profile::ogc().id_hex(), "burin1": Profile::burin_1().id_hex()},
        "polygons": polygons,
        "openings": openings,
        "transcripts": transcripts,
    });
    std::fs::write(dir.join("roots.json"), serde_json::to_string_pretty(&out).unwrap() + "\n").unwrap();
    println!("wrote roots.json: {} polygons, {} openings, {} transcripts", polygons.len(), openings.len(), transcripts.len());
}
