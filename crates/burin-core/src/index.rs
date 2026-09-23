//! A coverage index: the catalog-side half of the primitive.
//!
//! A STAC item carries a root; the index keeps the canonical tree the root names, once per
//! distinct root (every scene of one tile shares a footprint and so shares a tree), and answers
//! spatial questions by set algebra rather than geometry:
//!
//! - `query(q)`: for every indexed key, how much of `q` it covers and how much of itself lies in
//!   `q`, both as exact equal-area fractions, plus the set relation.
//! - `attention(level)`: for every cell at `level`, how many indexed trees touch it. Rolled up
//!   to a coarse level this is a fair density of where the catalog's effort is.
//!
//! Every answer is a function of trees only; the index never sees a coordinate.

use crate::error::{invalid, Result};
use crate::hash::{hex, Digest};
use crate::hierarchy::Cid;
use crate::setops::{difference, intersect};
use crate::tree::{Node, Tree};
use std::collections::{BTreeMap, HashMap};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Relation {
    Disjoint,
    Equal,
    /// The item lies within the query.
    Within,
    /// The item contains the query.
    Contains,
    Overlaps,
}

impl Relation {
    pub fn name(self) -> &'static str {
        match self {
            Relation::Disjoint => "disjoint",
            Relation::Equal => "equal",
            Relation::Within => "within",
            Relation::Contains => "contains",
            Relation::Overlaps => "overlaps",
        }
    }
}

#[derive(Debug, Clone)]
pub struct Hit {
    pub key: String,
    pub root: Digest,
    pub relation: Relation,
    /// |item ∩ query| / |query|: how much of the asked-for ground this item covers.
    pub fraction_of_query: f64,
    /// |item ∩ query| / |item|: how much of the item lies in the asked-for ground.
    pub fraction_of_item: f64,
    pub shared_cells: u64,
}

#[derive(Default)]
pub struct Index {
    trees: HashMap<Digest, Tree>,
    keys: BTreeMap<String, Digest>,
}

impl Index {
    pub fn new() -> Index {
        Index::default()
    }

    /// Register `key` (an item id, say) under its tree. Trees with the same root are stored once.
    pub fn add(&mut self, key: &str, tree: &Tree) -> Result<()> {
        if let Some(t) = self.trees.values().next() {
            if t.h != tree.h || t.d != tree.d || t.ctx.hasher.id() != tree.ctx.hasher.id() {
                return invalid("every tree in an index must share the hierarchy, depth and hash");
            }
        }
        let root = tree.root();
        self.trees.entry(root).or_insert_with(|| tree.clone());
        self.keys.insert(key.to_string(), root);
        Ok(())
    }

    pub fn len(&self) -> usize {
        self.keys.len()
    }

    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    pub fn distinct_roots(&self) -> usize {
        self.trees.len()
    }

    pub fn root_of(&self, key: &str) -> Option<Digest> {
        self.keys.get(key).copied()
    }

    /// Every key whose coverage meets `q`, with exact fractions, best coverage of the query first.
    pub fn query(&self, q: &Tree) -> Result<Vec<Hit>> {
        let q_leaves = q.leaf_count();
        // one set-op pass per distinct root, then fan out to the keys that share it
        let mut per_root: HashMap<Digest, (Relation, u64, u64)> = HashMap::new();
        for (root, t) in &self.trees {
            let shared = intersect(t, q)?;
            let n = shared.leaf_count();
            if n == 0 {
                continue;
            }
            let t_leaves = t.leaf_count();
            let item_in_q = difference(t, q)?.is_empty();
            let q_in_item = difference(q, t)?.is_empty();
            let rel = match (item_in_q, q_in_item) {
                (true, true) => Relation::Equal,
                (true, false) => Relation::Within,
                (false, true) => Relation::Contains,
                (false, false) => Relation::Overlaps,
            };
            per_root.insert(*root, (rel, n, t_leaves));
        }
        let mut hits = Vec::new();
        for (key, root) in &self.keys {
            if let Some(&(rel, n, t_leaves)) = per_root.get(root) {
                hits.push(Hit {
                    key: key.clone(),
                    root: *root,
                    relation: rel,
                    fraction_of_query: if q_leaves > 0 { n as f64 / q_leaves as f64 } else { 0.0 },
                    fraction_of_item: if t_leaves > 0 { n as f64 / t_leaves as f64 } else { 0.0 },
                    shared_cells: n,
                });
            }
        }
        hits.sort_by(|a, b| b.fraction_of_query.partial_cmp(&a.fraction_of_query).unwrap().then(a.key.cmp(&b.key)));
        Ok(hits)
    }

    /// For every cell at `level`, the number of indexed keys whose coverage touches it.
    pub fn attention(&self, level: u32) -> Result<BTreeMap<Cid, u64>> {
        let mut counts: BTreeMap<Cid, u64> = BTreeMap::new();
        let mut per_root: HashMap<Digest, Vec<Cid>> = HashMap::new();
        for (root, t) in &self.trees {
            if level > t.d {
                return invalid(format!("attention level {level} is deeper than the trees ({})", t.d));
            }
            per_root.insert(*root, touched_cells(t, level));
        }
        for root in self.keys.values() {
            for &c in &per_root[root] {
                *counts.entry(c).or_insert(0) += 1;
            }
        }
        Ok(counts)
    }

    pub fn roots(&self) -> Vec<(String, String)> {
        self.keys.iter().map(|(k, r)| (k.clone(), hex(r))).collect()
    }
}

/// Cells at `level` under which the tree has any coverage (full or partial).
fn touched_cells(t: &Tree, level: u32) -> Vec<Cid> {
    fn walk(node: &Node, cid: Cid, depth: u32, level: u32, a: u64, out: &mut Vec<Cid>) {
        match node {
            Node::Empty => {}
            _ if depth == level => out.push(cid),
            Node::Full => {
                // every descendant at `level` is touched
                let k = level - depth;
                out.extend(cid * a.pow(k)..(cid + 1) * a.pow(k));
            }
            Node::Branch(b) => {
                for (k, c) in b.children.iter().enumerate() {
                    walk(c, a * cid + k as u64, depth + 1, level, a, out);
                }
            }
        }
    }
    let mut out = Vec::new();
    for (b, node) in t.bases.iter().enumerate() {
        walk(node, (t.h.a + b as u32) as u64, 0, level, t.h.a as u64, &mut out);
    }
    out
}
