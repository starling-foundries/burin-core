//! Set algebra over coverage trees.
//!
//! Node-local identities on the ladder values decide a combination without descending:
//! ```text
//! union:      FULL ∪ X = FULL     EMPTY ∪ X = X
//! intersect:  EMPTY ∩ X = EMPTY   FULL ∩ X = X
//! difference: EMPTY \ X = EMPTY   X \ FULL = EMPTY   X \ EMPTY = X
//! ```
//! and equal hashes mean equal sets. A walk descends only over the interface (both partial and
//! different). The transcript records `(h_a, h_b, h_c)` per step; a verifier replays the
//! identities with no access to either tree.

use crate::error::{invalid, Result};
use crate::hash::{hex, unhex, Ctx, Digest};
use crate::hierarchy::Cid;
use crate::tree::{make_branch, node_hash, Node, Tree};
use serde_json::{json, Value};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    Union,
    Intersect,
    Difference,
}

impl Op {
    pub fn name(self) -> &'static str {
        match self {
            Op::Union => "union",
            Op::Intersect => "intersect",
            Op::Difference => "difference",
        }
    }
    pub fn parse(s: &str) -> Option<Op> {
        match s {
            "union" => Some(Op::Union),
            "intersect" => Some(Op::Intersect),
            "difference" => Some(Op::Difference),
            _ => None,
        }
    }
}

/// The result hash at a node, or `None` if both operands are partial and differ.
pub fn rule(op: Op, ha: &Digest, hb: &Digest, d: u32, ctx: &Ctx) -> Option<Digest> {
    let (f, e) = (ctx.ladders.full(d), ctx.ladders.empty(d));
    if ha == hb {
        return Some(if op == Op::Difference { e } else { *ha });
    }
    match op {
        Op::Union => {
            if *ha == f || *hb == f {
                Some(f)
            } else if *ha == e {
                Some(*hb)
            } else if *hb == e {
                Some(*ha)
            } else {
                None
            }
        }
        Op::Intersect => {
            if *ha == e || *hb == e {
                Some(e)
            } else if *ha == f {
                Some(*hb)
            } else if *hb == f {
                Some(*ha)
            } else {
                None
            }
        }
        Op::Difference => {
            if *ha == e || *hb == f {
                Some(e)
            } else if *hb == e {
                Some(*ha)
            } else {
                None
            }
        }
    }
}

fn compatible(a: &Tree, b: &Tree) -> Result<()> {
    if a.h != b.h {
        return invalid(format!("topology mismatch: {:?} vs {:?}", a.h, b.h));
    }
    if a.d != b.d {
        return invalid(format!("depth mismatch: {} vs {}; depth is part of a root's identity", a.d, b.d));
    }
    if a.ctx.hasher.id() != b.ctx.hasher.id() {
        return invalid("hash mismatch");
    }
    Ok(())
}

fn combine_nodes(na: &Node, nb: &Node, op: Op, d: u32, ctx: &Ctx) -> Node {
    let (ha, hb) = (node_hash(na, d, ctx), node_hash(nb, d, ctx));
    if let Some(r) = rule(op, &ha, &hb, d, ctx) {
        if r == ctx.ladders.empty(d) {
            return Node::Empty;
        }
        if r == ctx.ladders.full(d) {
            return Node::Full;
        }
        return if r == ha { na.clone() } else { nb.clone() };
    }
    debug_assert!(d > 0, "leaves are always decided");
    let kids = (0..ctx.a as usize).map(|i| combine_nodes(&na.child(i), &nb.child(i), op, d - 1, ctx)).collect();
    make_branch(kids, d, ctx)
}

pub fn combine(a: &Tree, b: &Tree, op: Op) -> Result<Tree> {
    compatible(a, b)?;
    let ctx = &a.ctx;
    let bases = (0..a.h.b as usize).map(|i| combine_nodes(&a.bases[i], &b.bases[i], op, a.d, ctx)).collect();
    Ok(Tree { h: a.h, d: a.d, bases, ctx: a.ctx.clone() })
}

pub fn union(a: &Tree, b: &Tree) -> Result<Tree> {
    combine(a, b, Op::Union)
}
pub fn intersect(a: &Tree, b: &Tree) -> Result<Tree> {
    combine(a, b, Op::Intersect)
}
pub fn difference(a: &Tree, b: &Tree) -> Result<Tree> {
    combine(a, b, Op::Difference)
}

/// Union of one or more trees.
pub fn merge<'a>(trees: impl IntoIterator<Item = &'a Tree>) -> Result<Tree> {
    let mut it = trees.into_iter();
    let mut acc = match it.next() {
        Some(t) => t.clone(),
        None => return invalid("merge() needs at least one tree"),
    };
    for t in it {
        acc = union(&acc, t)?;
    }
    Ok(acc)
}

pub fn intersects(a: &Tree, b: &Tree) -> Result<bool> {
    Ok(!intersect(a, b)?.is_empty())
}
pub fn disjoint(a: &Tree, b: &Tree) -> Result<bool> {
    Ok(intersect(a, b)?.is_empty())
}
/// B ⊆ A.
pub fn contains(a: &Tree, b: &Tree) -> Result<bool> {
    Ok(difference(b, a)?.is_empty())
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Divergence {
    pub cid: Cid,
    pub hash_a: Digest,
    pub hash_b: Digest,
}

/// Maximal subtrees where two coverages differ; O(differences).
pub fn divergence(a: &Tree, b: &Tree) -> Result<Vec<Divergence>> {
    compatible(a, b)?;
    fn walk(na: &Node, nb: &Node, d: u32, cid: Cid, ctx: &Ctx, out: &mut Vec<Divergence>) {
        let (ha, hb) = (node_hash(na, d, ctx), node_hash(nb, d, ctx));
        if ha == hb {
            return;
        }
        if d == 0 || !(na.is_branch() && nb.is_branch()) {
            out.push(Divergence { cid, hash_a: ha, hash_b: hb });
            return;
        }
        for i in 0..ctx.a as usize {
            walk(&na.child(i), &nb.child(i), d - 1, ctx.a as u64 * cid + i as u64, ctx, out);
        }
    }
    let mut out = Vec::new();
    for i in 0..a.h.b as usize {
        walk(&a.bases[i], &b.bases[i], a.d, (a.h.a + i as u32) as u64, &a.ctx, &mut out);
    }
    Ok(out)
}

#[derive(Debug, Clone, PartialEq)]
pub struct Step {
    pub a: Digest,
    pub b: Digest,
    pub c: Digest,
    pub children: Option<Vec<Step>>,
}

impl Step {
    fn to_json(&self) -> Value {
        let mut v = json!({"a": hex(&self.a), "b": hex(&self.b), "c": hex(&self.c)});
        if let Some(k) = &self.children {
            v["children"] = Value::Array(k.iter().map(Step::to_json).collect());
        }
        v
    }
    fn from_json(v: &Value) -> Result<Step> {
        let h = |k: &str| v.get(k).and_then(Value::as_str).and_then(unhex).ok_or_else(|| crate::Error::Invalid(format!("bad step hash {k}")));
        let children = match v.get("children") {
            Some(Value::Array(ks)) => Some(ks.iter().map(Step::from_json).collect::<Result<Vec<_>>>()?),
            Some(_) => return invalid("children must be an array"),
            None => None,
        };
        Ok(Step { a: h("a")?, b: h("b")?, c: h("c")?, children })
    }
}

pub const SETOP_WIRE_VERSION: u64 = 1;

/// Evidence that `root_c = root_a ⊕ root_b` for `op`; size O(interface).
#[derive(Debug, Clone, PartialEq)]
pub struct SetOpProof {
    pub op: Op,
    pub hash: String,
    pub root_a: Digest,
    pub root_b: Digest,
    pub root_c: Digest,
    pub a: u32,
    pub b: u32,
    pub d: u32,
    pub steps: Vec<Step>,
}

impl SetOpProof {
    pub fn size(&self) -> usize {
        fn count(s: &Step) -> usize {
            1 + s.children.as_ref().map_or(0, |k| k.iter().map(count).sum())
        }
        self.steps.iter().map(count).sum()
    }

    /// Replay the identities; `false` on any defect. Never panics.
    pub fn verify(&self) -> bool {
        let h = match crate::hash::hasher_for(&self.hash) {
            Some(h) => h,
            None => return false,
        };
        if !crate::opening::valid_shape(self.a, self.b, self.d) || self.steps.len() != self.b as usize {
            return false;
        }
        let ctx = Ctx::new(h, self.a, self.b, self.d);
        fn check(s: &Step, d: u32, op: Op, ctx: &Ctx) -> bool {
            if let Some(r) = rule(op, &s.a, &s.b, d, ctx) {
                return s.children.is_none() && s.c == r; // a decided node must be a leaf step
            }
            let kids = match &s.children {
                Some(k) if d > 0 && k.len() == ctx.a as usize => k,
                _ => return false,
            };
            if !kids.iter().all(|c| check(c, d - 1, op, ctx)) {
                return false;
            }
            let col = |f: fn(&Step) -> Digest| -> Digest { ctx.hasher.node(&kids.iter().map(f).collect::<Vec<_>>()) };
            s.a == col(|c| c.a) && s.b == col(|c| c.b) && s.c == col(|c| c.c)
        }
        if !self.steps.iter().all(|s| check(s, self.d, self.op, &ctx)) {
            return false;
        }
        let col = |f: fn(&Step) -> Digest| -> Digest { ctx.hasher.root(&self.steps.iter().map(f).collect::<Vec<_>>()) };
        col(|s| s.a) == self.root_a && col(|s| s.b) == self.root_b && col(|s| s.c) == self.root_c
    }

    pub fn to_json(&self) -> Value {
        json!({
            "v": SETOP_WIRE_VERSION, "op": self.op.name(), "hash": self.hash,
            "root_a": hex(&self.root_a), "root_b": hex(&self.root_b), "root_c": hex(&self.root_c),
            "A": self.a, "B": self.b, "D": self.d,
            "steps": self.steps.iter().map(Step::to_json).collect::<Vec<_>>(),
        })
    }

    pub fn from_json(v: &Value) -> Result<SetOpProof> {
        let ver = v.get("v").and_then(Value::as_u64).unwrap_or(0);
        if ver != SETOP_WIRE_VERSION {
            return invalid(format!("set-operation proof is version v{ver}; this build reads v{SETOP_WIRE_VERSION}"));
        }
        let op = v.get("op").and_then(Value::as_str).and_then(Op::parse).ok_or_else(|| crate::Error::Invalid("bad op".into()))?;
        let h = |k: &str| v.get(k).and_then(Value::as_str).and_then(unhex).ok_or_else(|| crate::Error::Invalid(format!("bad {k}")));
        let (a, b, d) = crate::opening::wire_shape(v)?;
        let steps = match v.get("steps") {
            Some(Value::Array(ss)) => ss.iter().map(Step::from_json).collect::<Result<Vec<_>>>()?,
            _ => return invalid("missing steps"),
        };
        Ok(SetOpProof {
            op,
            hash: v.get("hash").and_then(Value::as_str).unwrap_or("").to_string(),
            root_a: h("root_a")?,
            root_b: h("root_b")?,
            root_c: h("root_c")?,
            a,
            b,
            d,
            steps,
        })
    }
}

pub fn prove(a: &Tree, b: &Tree, op: Op) -> Result<SetOpProof> {
    compatible(a, b)?;
    let ctx = &a.ctx;
    fn build(na: &Node, nb: &Node, d: u32, op: Op, ctx: &Ctx) -> Step {
        let (ha, hb) = (node_hash(na, d, ctx), node_hash(nb, d, ctx));
        if let Some(r) = rule(op, &ha, &hb, d, ctx) {
            return Step { a: ha, b: hb, c: r, children: None };
        }
        let kids: Vec<Step> = (0..ctx.a as usize).map(|i| build(&na.child(i), &nb.child(i), d - 1, op, ctx)).collect();
        let c = ctx.hasher.node(&kids.iter().map(|k| k.c).collect::<Vec<_>>());
        Step { a: ha, b: hb, c, children: Some(kids) }
    }
    let steps: Vec<Step> = (0..a.h.b as usize).map(|i| build(&a.bases[i], &b.bases[i], a.d, op, ctx)).collect();
    let root_c = ctx.hasher.root(&steps.iter().map(|s| s.c).collect::<Vec<_>>());
    Ok(SetOpProof { op, hash: ctx.hasher.id().to_string(), root_a: a.root(), root_b: b.root(), root_c, a: a.h.a, b: a.h.b, d: a.d, steps })
}
