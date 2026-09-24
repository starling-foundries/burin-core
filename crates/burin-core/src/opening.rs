//! The opening: the one proof shape (coverage claims only).
//!
//! An opening is a tree of steps. A step is either a claim ("this node is constant: empty or
//! full") or `entries`: A child entries (B at the root level), each a bare hash or a nested
//! opening. The verifier recomputes hashes bottom-up and compares with the committed root. The
//! cell a proof speaks about is the sequence of opened positions, so it is structural: a proof
//! moved to another position verifies only where its claim is also true, as between the identical
//! children of a constant node.
//!
//! Wire form: `{"claim": "full" | "empty"}` or `{"entries": [hex | opening, ...]}`, wrapped in
//! `{"v": 1, "hash": id, "A": a, "B": b, "D": d, "root": hex, "opening": {...}}`.

use crate::error::{invalid, Result};
use crate::hash::{hex, unhex, Ctx, Digest};
use crate::hierarchy::Cid;
use crate::tree::{node_hash, Node, Tree};
use serde_json::{json, Value};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Claim {
    Empty,
    Full,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Entry {
    Hash(Digest),
    Open(Box<Opening>),
}

#[derive(Debug, Clone, PartialEq)]
pub enum Opening {
    Claim(Claim),
    Entries(Vec<Entry>),
}

impl Opening {
    /// The opened positions `(b, d_1, ...)` if exactly one entry is opened per level.
    pub fn path(&self) -> Option<Vec<u32>> {
        let mut out = Vec::new();
        let mut op = self;
        loop {
            match op {
                Opening::Claim(_) => return Some(out),
                Opening::Entries(es) => {
                    let opened: Vec<(usize, &Opening)> = es
                        .iter()
                        .enumerate()
                        .filter_map(|(i, e)| match e {
                            Entry::Open(o) => Some((i, o.as_ref())),
                            _ => None,
                        })
                        .collect();
                    if opened.len() != 1 {
                        return None;
                    }
                    out.push(opened[0].0 as u32);
                    op = opened[0].1;
                }
            }
        }
    }

    /// The claim at the end of a single-path opening.
    pub fn terminal(&self) -> Option<Claim> {
        let mut op = self;
        loop {
            match op {
                Opening::Claim(c) => return Some(*c),
                Opening::Entries(es) => {
                    let opened: Vec<&Opening> = es
                        .iter()
                        .filter_map(|e| match e {
                            Entry::Open(o) => Some(o.as_ref()),
                            _ => None,
                        })
                        .collect();
                    if opened.len() != 1 {
                        return None;
                    }
                    op = opened[0];
                }
            }
        }
    }

    pub fn to_json(&self) -> Value {
        match self {
            Opening::Claim(Claim::Empty) => json!({"claim": "empty"}),
            Opening::Claim(Claim::Full) => json!({"claim": "full"}),
            Opening::Entries(es) => Value::Array(
                es.iter()
                    .map(|e| match e {
                        Entry::Hash(h) => Value::String(hex(h)),
                        Entry::Open(o) => o.to_json(),
                    })
                    .collect(),
            )
            .pipe(|arr| json!({"entries": arr})),
        }
    }

    pub fn from_json(v: &Value) -> Result<Opening> {
        if let Some(c) = v.get("claim") {
            return match c.as_str() {
                Some("empty") => Ok(Opening::Claim(Claim::Empty)),
                Some("full") => Ok(Opening::Claim(Claim::Full)),
                _ => invalid(format!("unknown claim {c}")),
            };
        }
        if let Some(es) = v.get("entries").and_then(Value::as_array) {
            let mut out = Vec::with_capacity(es.len());
            for e in es {
                if let Some(s) = e.as_str() {
                    match unhex(s) {
                        Some(h) => out.push(Entry::Hash(h)),
                        None => return invalid(format!("bad hash {s}")),
                    }
                } else {
                    out.push(Entry::Open(Box::new(Opening::from_json(e)?)));
                }
            }
            return Ok(Opening::Entries(out));
        }
        invalid("an opening is a claim or entries")
    }
}

trait Pipe: Sized {
    fn pipe<T>(self, f: impl FnOnce(Self) -> T) -> T {
        f(self)
    }
}
impl<T> Pipe for T {}

/// Open `tree` along the path to `cid`. The terminal claim is the constant found at `cid`, which
/// may be inherited from a constant ancestor. `None` iff `cid` is a partial node.
pub fn open_path(tree: &Tree, cid: Cid) -> Result<Option<Opening>> {
    tree.h.check(cid, Some(tree.d))?;
    let path = tree.h.path(cid)?;
    let a = tree.h.a as usize;
    let ctx = &tree.ctx;

    fn go(node: &Node, digits: &[u32], d: u32, a: usize, ctx: &Ctx) -> Option<Opening> {
        if digits.is_empty() {
            return match node {
                Node::Branch(_) => None,
                Node::Empty => Some(Opening::Claim(Claim::Empty)),
                Node::Full => Some(Opening::Claim(Claim::Full)),
            };
        }
        let k = digits[0] as usize;
        // a constant node passes through: its children are A copies of itself and its siblings
        // hash to the ladder one level down, so the recomputation stays exact
        let child = node.child(k);
        let nested = go(&child, &digits[1..], d - 1, a, ctx)?;
        let entries = (0..a)
            .map(|i| if i == k { Entry::Open(Box::new(nested.clone())) } else { Entry::Hash(node_hash(&node.child(i), d - 1, ctx)) })
            .collect();
        Some(Opening::Entries(entries))
    }

    let b = path[0] as usize;
    let nested = match go(&tree.bases[b], &path[1..], tree.d, a, ctx) {
        Some(n) => n,
        None => return Ok(None),
    };
    let entries = (0..tree.h.b as usize)
        .map(|i| if i == b { Entry::Open(Box::new(nested.clone())) } else { Entry::Hash(tree.subtree_hash(i)) })
        .collect();
    Ok(Some(Opening::Entries(entries)))
}

/// The node hash an opening recomputes at `d` levels above the leaves, or `None` if malformed.
pub fn opening_hash(op: &Opening, d: u32, ctx: &Ctx) -> Option<Digest> {
    match op {
        Opening::Claim(Claim::Empty) => Some(ctx.ladders.empty(d)),
        Opening::Claim(Claim::Full) => Some(ctx.ladders.full(d)),
        Opening::Entries(es) => {
            if d < 1 || es.len() != ctx.a as usize {
                return None;
            }
            let mut kids = Vec::with_capacity(es.len());
            for e in es {
                kids.push(match e {
                    Entry::Hash(h) => *h,
                    Entry::Open(o) => opening_hash(o, d - 1, ctx)?,
                });
            }
            Some(ctx.hasher.node(&kids))
        }
    }
}

/// Stateless: recompute the root from the opening. Never panics on malformed input.
pub fn verify_opening(op: &Opening, root: &Digest, d: u32, ctx: &Ctx) -> bool {
    if d > ctx.ladders.depth() {
        return false;
    }
    let es = match op {
        Opening::Entries(es) if es.len() == ctx.b as usize => es,
        _ => return false,
    };
    let mut bases = Vec::with_capacity(es.len());
    for e in es {
        bases.push(match e {
            Entry::Hash(h) => *h,
            Entry::Open(o) => match opening_hash(o, d, ctx) {
                Some(h) => h,
                None => return false,
            },
        });
    }
    ctx.hasher.root(&bases) == *root
}

/// A self-describing opening record: what the proof is about and how to check it.
#[derive(Debug, Clone, PartialEq)]
pub struct OpeningRecord {
    pub hash: String,
    pub a: u32,
    pub b: u32,
    pub d: u32,
    pub root: Digest,
    pub opening: Opening,
}

pub const OPENING_WIRE_VERSION: u64 = 1;

/// A record's `A`, `B` and `D`: a hierarchy (SPEC §2) and a depth no deeper than its deepest level.
/// A number that does not fit is refused, never truncated, so each field has one spelling.
pub(crate) fn wire_shape(v: &Value) -> Result<(u32, u32, u32)> {
    let num = |k: &str| -> Result<u32> {
        let n = v.get(k).and_then(Value::as_u64).ok_or_else(|| crate::Error::Invalid(format!("missing {k}")))?;
        u32::try_from(n).map_err(|_| crate::Error::Invalid(format!("{k} = {n} is out of range")))
    };
    let (a, b, d) = (num("A")?, num("B")?, num("D")?);
    let h = crate::hierarchy::Hierarchy::new(a, b)?;
    if d > h.max_level() {
        return invalid(format!("D = {d} is deeper than the deepest level, {}", h.max_level()));
    }
    Ok((a, b, d))
}

/// Whether a record's shape is one `wire_shape` accepts.
pub(crate) fn valid_shape(a: u32, b: u32, d: u32) -> bool {
    crate::hierarchy::Hierarchy::new(a, b).is_ok_and(|h| d <= h.max_level())
}

impl OpeningRecord {
    pub fn new(tree: &Tree, opening: Opening) -> OpeningRecord {
        OpeningRecord { hash: tree.ctx.hasher.id().to_string(), a: tree.h.a, b: tree.h.b, d: tree.d, root: tree.root(), opening }
    }

    /// The cid this record opens, if it is a single path.
    pub fn cid(&self) -> Option<Cid> {
        let path = self.opening.path()?;
        crate::hierarchy::Hierarchy { a: self.a, b: self.b }.cid(&path).ok()
    }

    pub fn to_json(&self) -> Value {
        json!({
            "v": OPENING_WIRE_VERSION, "hash": self.hash, "A": self.a, "B": self.b, "D": self.d,
            "root": hex(&self.root), "opening": self.opening.to_json(),
        })
    }

    pub fn from_json(v: &Value) -> Result<OpeningRecord> {
        let ver = v.get("v").and_then(Value::as_u64).unwrap_or(0);
        if ver != OPENING_WIRE_VERSION {
            return invalid(format!("opening record is version v{ver}; this build reads v{OPENING_WIRE_VERSION}"));
        }
        let (a, b, d) = wire_shape(v)?;
        let root = v.get("root").and_then(Value::as_str).and_then(unhex).ok_or_else(|| crate::Error::Invalid("bad root".into()))?;
        Ok(OpeningRecord {
            hash: v.get("hash").and_then(Value::as_str).unwrap_or("").to_string(),
            a,
            b,
            d,
            root,
            opening: Opening::from_json(v.get("opening").ok_or_else(|| crate::Error::Invalid("missing opening".into()))?)?,
        })
    }

    /// Verify against a context built from the record's own parameters. `false` on any defect.
    pub fn verify(&self) -> bool {
        let h = match crate::hash::hasher_for(&self.hash) {
            Some(h) => h,
            None => return false,
        };
        if !valid_shape(self.a, self.b, self.d) {
            return false;
        }
        let ctx = Ctx::new(h, self.a, self.b, self.d);
        verify_opening(&self.opening, &self.root, self.d, &ctx)
    }
}
