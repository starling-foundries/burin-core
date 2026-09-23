//! OGC API - DGGS zone data (DGGS-JSON) for coverages: one field at one relative depth over the
//! sub-zones of one zone, listed in the DGGRS's sub-zone order (scanline for rHEALPix, SPEC §8).
//!
//! ```text
//! {"dggrs": uri, "zoneId": "Q4", "depths": [k],
//!  "values": {field: [{"depth": k, "shape": {"count": 9^k, "subZones": 9^k}, "data": [...]}]}}
//! ```
//! Writing puts 1 at each covered sub-zone and null elsewhere. Reading takes one field and one
//! depth and covers the sub-zones whose value is present (non-null, or non-zero if asked). A
//! document read against a profile whose grid is not the DGGRS it names is refused.

use crate::error::{invalid, Result};
use crate::hierarchy::{cid_to_suid, suid_to_cid, Cid, SPACE};
use crate::profile::{Profile, OGC_RHEALPIX};
use crate::tree::{make_branch, Node, Tree};
use crate::zone::{block_node, position, side, MAX_RASTER_CELLS};
use serde_json::{json, Value};
use std::sync::Arc;

/// The OGC-registered rHEALPix DGGRS.
pub const RHEALPIX_DGGRS: &str = "https://www.opengis.net/def/dggrs/OGC/1.0/rHEALPix";

/// Which values mark a sub-zone as covered when reading zone data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Presence {
    /// Any value that is not null.
    NonNull,
    /// Any value that is neither null nor zero.
    NonZero,
}

/// The registered DGGRS whose zones this profile's cells are, if there is one. The hash and the
/// time convention do not change what a zone is, so they are not compared.
pub fn dggrs_uri(profile: &Profile) -> Option<&'static str> {
    let g = |p: &Profile| (p.aperture, p.n_base, p.lon_0_udeg, p.north_square, p.south_square, p.a_um, p.inv_f_nano);
    (g(profile) == g(&OGC_RHEALPIX)).then_some(RHEALPIX_DGGRS)
}

fn registered(profile: &Profile) -> Result<&'static str> {
    profile.validate()?;
    match dggrs_uri(profile) {
        Some(uri) => Ok(uri),
        None => invalid(format!("profile {} is not a registered DGGRS; its zone ids would be read against another grid", profile.describe())),
    }
}

fn count_at(depth: u32) -> Result<u64> {
    let width = side(&SPACE)?.checked_pow(depth);
    let count = width.and_then(|w| w.checked_mul(w)).filter(|&c| c <= MAX_RASTER_CELLS);
    count.ok_or_else(|| crate::Error::Invalid(format!("depth {depth} has too many sub-zones")))
}

/// The coverage of `tree` over the sub-zones of `zone` at the tree's depth, as DGGS-JSON.
pub fn to_dggs_json(tree: &Tree, profile: &Profile, zone: Cid, field: &str) -> Result<Value> {
    let uri = registered(profile)?;
    if tree.h != SPACE || profile.hierarchy() != SPACE {
        return invalid("zone data needs the rHEALPix hierarchy");
    }
    let level = SPACE.check(zone, Some(tree.d))?;
    let depth = tree.d - level;
    let count = count_at(depth)?;
    let n = side(&SPACE)?;
    let width = n.pow(depth);
    let pz = position(&SPACE, zone)?;
    let mut data = vec![Value::Null; count as usize];
    for c in tree.cells() {
        if SPACE.is_ancestor(c, zone) {
            data.iter_mut().for_each(|v| *v = json!(1));
            break;
        }
        if !SPACE.is_ancestor(zone, c) {
            continue;
        }
        let pc = position(&SPACE, c)?;
        let block = n.pow(tree.d - pc.level);
        let (r0, c0) = (pc.row * block - pz.row * width, pc.col * block - pz.col * width);
        for r in r0..r0 + block {
            for v in &mut data[(r * width + c0) as usize..(r * width + c0 + block) as usize] {
                *v = json!(1);
            }
        }
    }
    Ok(json!({
        "dggrs": uri,
        "zoneId": cid_to_suid(zone)?,
        "depths": [depth],
        "values": {field: [{"depth": depth, "shape": {"count": count, "subZones": count}, "data": data}]},
    }))
}

fn pick<'a>(what: &str, items: Vec<(String, &'a Value)>, want: Option<String>) -> Result<&'a Value> {
    let names: Vec<String> = items.iter().map(|(k, _)| k.clone()).collect();
    match want {
        Some(w) => match items.into_iter().find(|(k, _)| *k == w) {
            Some((_, v)) => Ok(v),
            None => invalid(format!("no {what} {w:?}; the document has {names:?}")),
        },
        None if items.len() == 1 => Ok(items[0].1),
        None => invalid(format!("the document has several {what}s {names:?}; name one")),
    }
}

/// The coverage a DGGS-JSON document describes: the sub-zones of its zone, at one depth, whose
/// value in `field` is present. `field` and `depth` may be omitted when the document has one.
pub fn from_dggs_json(doc: &Value, profile: &Profile, field: Option<&str>, depth: Option<u32>, presence: Presence) -> Result<Tree> {
    let uri = registered(profile)?;
    if profile.hierarchy() != SPACE {
        return invalid("zone data needs the rHEALPix hierarchy");
    }
    match doc.get("dggrs").and_then(Value::as_str) {
        Some(d) if d == uri => {}
        other => return invalid(format!("document DGGRS {other:?} is not the profile's {uri:?}")),
    }
    if doc.get("dimensions").and_then(Value::as_array).is_some_and(|d| !d.is_empty()) {
        return invalid("zone data with extra dimensions has more than one value per sub-zone");
    }
    let zone = suid_to_cid(doc.get("zoneId").and_then(Value::as_str).ok_or_else(|| crate::Error::Invalid("missing zoneId".into()))?)?;
    let fields = doc.get("values").and_then(Value::as_object).ok_or_else(|| crate::Error::Invalid("missing values".into()))?;
    let entries = pick("field", fields.iter().map(|(k, v)| (k.clone(), v)).collect(), field.map(str::to_string))?;
    let entries = entries.as_array().ok_or_else(|| crate::Error::Invalid("a field's values must be an array".into()))?;
    let doc_depths: Vec<u64> = doc.get("depths").and_then(Value::as_array).map(|d| d.iter().filter_map(Value::as_u64).collect()).unwrap_or_default();
    let depth_of = |e: &Value| e.get("depth").and_then(Value::as_u64).or((doc_depths.len() == 1).then(|| doc_depths[0]));
    let mut labelled = Vec::new();
    for e in entries {
        let d = depth_of(e).ok_or_else(|| crate::Error::Invalid("a value entry names no depth".into()))?;
        labelled.push((d.to_string(), e));
    }
    let entry = pick("depth", labelled, depth.map(|d| d.to_string()))?;
    let k = depth_of(entry).expect("labelled above") as u32;
    let count = count_at(k)?;
    let shape = entry.get("shape");
    let shape_n = |key: &str| shape.and_then(|s| s.get(key)).and_then(Value::as_u64);
    if shape_n("subZones").is_some_and(|s| s != count) || shape_n("count").is_some_and(|c| c != count) {
        return invalid(format!("shape {shape:?} does not describe the {count} sub-zones at depth {k}"));
    }
    let data = entry.get("data").and_then(Value::as_array).ok_or_else(|| crate::Error::Invalid("missing data".into()))?;
    if data.len() as u64 != count {
        return invalid(format!("{} values for the {count} sub-zones at depth {k}", data.len()));
    }
    let mut covered = Vec::with_capacity(data.len());
    for v in data {
        covered.push(match (v, presence) {
            (Value::Null, _) => false,
            (Value::Number(_), Presence::NonNull) => true,
            (Value::Number(x), Presence::NonZero) => x.as_f64() != Some(0.0),
            _ => return invalid(format!("zone data value {v} is not a number or null")),
        });
    }
    let pz = position(&SPACE, zone)?;
    let d = pz.level + k;
    let ctx = Arc::new(profile.ctx(d)?);
    let n = side(&SPACE)?;
    let width = n.pow(k);
    let at = |_: u64, row: u64, col: u64| covered[((row - pz.row * width) * width + (col - pz.col * width)) as usize];
    let mut node = block_node(&at, n, pz.base as u64, pz.level, pz.row, pz.col, d, &ctx);
    let mut cid = zone;
    while let Some(parent) = SPACE.parent(cid) {
        let mut children = vec![Node::Empty; SPACE.a as usize];
        children[(cid % SPACE.a as u64) as usize] = node;
        node = make_branch(children, d - SPACE.level(parent), &ctx);
        cid = parent;
    }
    let mut tree = Tree::empty(SPACE, d, ctx)?;
    tree.bases[pz.base as usize] = node;
    Ok(tree)
}
