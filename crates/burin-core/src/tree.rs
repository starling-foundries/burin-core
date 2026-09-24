//! The coverage tree: an immutable, canonical, sparse hierarchy tree.
//!
//! ```text
//! Empty          no leaf under this node is covered
//! Full           every leaf under this node is covered
//! Branch         A children, not all the same constant; hash stored on the node
//! ```
//! The structure is always in its coarsest form (coalesced on write) and never mutated: every
//! operation returns a new tree sharing untouched subtrees. The root is a function of the covered
//! leaf set alone, whatever order or tiling built it (the roll-up theorem).

use crate::error::{invalid, Result};
use crate::hash::{Ctx, Digest};
use crate::hierarchy::{Cid, Hierarchy};
use std::sync::Arc;

#[derive(Debug)]
pub struct Branch {
    pub children: Vec<Node>,
    pub hash: Digest,
}

#[derive(Debug, Clone)]
pub enum Node {
    Empty,
    Full,
    Branch(Arc<Branch>),
}

impl PartialEq for Node {
    fn eq(&self, other: &Node) -> bool {
        match (self, other) {
            (Node::Empty, Node::Empty) | (Node::Full, Node::Full) => true,
            (Node::Branch(a), Node::Branch(b)) => a.hash == b.hash,
            _ => false,
        }
    }
}

impl Node {
    pub fn is_branch(&self) -> bool {
        matches!(self, Node::Branch(_))
    }

    /// The k-th child: a constant's children are itself.
    pub fn child(&self, k: usize) -> Node {
        match self {
            Node::Branch(b) => b.children[k].clone(),
            c => c.clone(),
        }
    }
}

/// Hash of any node `d` levels above the leaves.
pub fn node_hash(node: &Node, d: u32, ctx: &Ctx) -> Digest {
    match node {
        Node::Empty => ctx.ladders.empty(d),
        Node::Full => ctx.ladders.full(d),
        Node::Branch(b) => b.hash,
    }
}

/// The canonical node for these children at `d` levels above the leaves.
pub fn make_branch(children: Vec<Node>, d: u32, ctx: &Ctx) -> Node {
    if children.iter().all(|c| matches!(c, Node::Empty)) {
        return Node::Empty;
    }
    if children.iter().all(|c| matches!(c, Node::Full)) {
        return Node::Full;
    }
    let hashes: Vec<Digest> = children.iter().map(|c| node_hash(c, d - 1, ctx)).collect();
    let hash = ctx.hasher.node(&hashes);
    Node::Branch(Arc::new(Branch { children, hash }))
}

#[derive(Clone)]
pub struct Tree {
    pub h: Hierarchy,
    pub d: u32,
    pub bases: Vec<Node>,
    pub ctx: Arc<Ctx>,
}

impl core::fmt::Debug for Tree {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Tree(H({},{}), D={}, root={})", self.h.a, self.h.b, self.d, crate::hash::hex(&self.root()))
    }
}

impl Tree {
    pub fn empty(h: Hierarchy, d: u32, ctx: Arc<Ctx>) -> Result<Tree> {
        if ctx.a != h.a || ctx.b != h.b {
            return invalid("hashing context does not match the hierarchy");
        }
        if ctx.ladders.depth() < d {
            return invalid(format!("context ladders reach depth {}, tree needs {d}", ctx.ladders.depth()));
        }
        if d > h.max_level() {
            return invalid(format!("depth {d} is deeper than the deepest level, {}", h.max_level()));
        }
        Ok(Tree { h, d, bases: vec![Node::Empty; h.b as usize], ctx })
    }

    /// Cover every cid (any levels up to `d`, any order, duplicates allowed). Built in one pass:
    /// each cell becomes its range of leaves, the ranges are merged, and the canonical tree is
    /// read off them top-down, so each node is hashed once.
    pub fn from_cells(h: Hierarchy, d: u32, cids: impl IntoIterator<Item = Cid>, ctx: Arc<Ctx>) -> Result<Tree> {
        let empty = Tree::empty(h, d, ctx)?;
        let a = h.a as u64;
        let mut spans = Vec::new();
        for c in cids {
            let r = h.check(c, Some(d))?;
            let width = a.pow(d - r);
            spans.push((c * width, (c + 1) * width));
        }
        spans.sort_unstable();
        let mut merged: Vec<(u64, u64)> = Vec::with_capacity(spans.len());
        for (lo, hi) in spans {
            match merged.last_mut() {
                Some(last) if lo <= last.1 => last.1 = last.1.max(hi),
                _ => merged.push((lo, hi)),
            }
        }
        fn build(lo: u64, hi: u64, above: u32, spans: &[(u64, u64)], a: u64, ctx: &Ctx) -> Node {
            let spans = &spans[spans.partition_point(|s| s.1 <= lo)..];
            let spans = &spans[..spans.partition_point(|s| s.0 < hi)];
            match spans.first() {
                None => return Node::Empty,
                Some(&(s, e)) if s <= lo && e >= hi => return Node::Full,
                _ => {}
            }
            let step = (hi - lo) / a;
            let children = (0..a).map(|k| build(lo + k * step, lo + (k + 1) * step, above - 1, spans, a, ctx)).collect();
            make_branch(children, above, ctx)
        }
        let span = a.pow(d);
        let bases = (0..h.b as u64).map(|b| build((a + b) * span, (a + b + 1) * span, d, &merged, a, &empty.ctx)).collect();
        Ok(Tree { bases, ..empty })
    }

    fn set_rec(&self, node: &Node, digits: &[u32], d: u32) -> Node {
        if digits.is_empty() {
            return Node::Full;
        }
        if let Node::Full = node {
            return Node::Full; // already covered by a constant ancestor
        }
        let a = self.h.a as usize;
        let k = digits[0] as usize;
        let mut children: Vec<Node> = match node {
            Node::Branch(b) => b.children.clone(),
            other => vec![other.clone(); a],
        };
        let new = self.set_rec(&children[k], &digits[1..], d - 1);
        if new == children[k] {
            return node.clone();
        }
        children[k] = new;
        make_branch(children, d, &self.ctx)
    }

    /// A new tree with the whole subtree at `cid` covered.
    pub fn set_full(&self, cid: Cid) -> Result<Tree> {
        self.h.check(cid, Some(self.d))?;
        let path = self.h.path(cid)?;
        let b = path[0] as usize;
        let new_base = self.set_rec(&self.bases[b], &path[1..], self.d);
        if new_base == self.bases[b] {
            return Ok(self.clone());
        }
        let mut bases = self.bases.clone();
        bases[b] = new_base;
        Ok(Tree { h: self.h, d: self.d, bases, ctx: self.ctx.clone() })
    }

    pub fn base_hashes(&self) -> Vec<Digest> {
        self.bases.iter().map(|n| node_hash(n, self.d, &self.ctx)).collect()
    }

    pub fn root(&self) -> Digest {
        self.ctx.hasher.root(&self.base_hashes())
    }

    pub fn root_hex(&self) -> String {
        crate::hash::hex(&self.root())
    }

    /// The hash of one base subtree.
    pub fn subtree_hash(&self, base: usize) -> Digest {
        node_hash(&self.bases[base], self.d, &self.ctx)
    }

    /// `(node, d, exact)`: the node at `cid` with its levels-below count, or the constant
    /// ancestor reached first (`exact` false).
    pub fn node_at(&self, cid: Cid) -> Result<(Node, u32, bool)> {
        self.h.check(cid, Some(self.d))?;
        let path = self.h.path(cid)?;
        let mut node = self.bases[path[0] as usize].clone();
        let mut d = self.d;
        for &k in &path[1..] {
            match node {
                Node::Branch(ref b) => {
                    let next = b.children[k as usize].clone();
                    node = next;
                    d -= 1;
                }
                _ => return Ok((node, d, false)),
            }
        }
        Ok((node, d, true))
    }

    /// Whether the leaf `cid` (at any level: the whole subtree) is covered.
    pub fn covers(&self, cid: Cid) -> Result<bool> {
        Ok(matches!(self.node_at(cid)?.0, Node::Full))
    }

    /// The canonical (coarsest) covered cids, sorted.
    pub fn cells(&self) -> Vec<Cid> {
        fn walk(node: &Node, cid: Cid, a: u64, out: &mut Vec<Cid>) {
            match node {
                Node::Empty => {}
                Node::Full => out.push(cid),
                Node::Branch(b) => {
                    for (k, c) in b.children.iter().enumerate() {
                        walk(c, a * cid + k as u64, a, out);
                    }
                }
            }
        }
        let mut out = Vec::new();
        for (b, node) in self.bases.iter().enumerate() {
            walk(node, (self.h.a + b as u32) as u64, self.h.a as u64, &mut out);
        }
        out.sort_unstable();
        out
    }

    /// Every covered leaf at depth `d`, sorted (expands the canonical cells).
    pub fn leaves(&self) -> Vec<Cid> {
        let mut out = Vec::new();
        for c in self.cells() {
            let k = self.d - self.h.level(c);
            out.extend(self.h.descendants(c, k));
        }
        out
    }

    pub fn is_empty(&self) -> bool {
        self.bases.iter().all(|b| matches!(b, Node::Empty))
    }

    pub fn is_full(&self) -> bool {
        self.bases.iter().all(|b| matches!(b, Node::Full))
    }

    /// Number of covered leaves at depth `d` (exact, from the canonical cells).
    pub fn leaf_count(&self) -> u64 {
        self.cells().iter().map(|&c| (self.h.a as u64).pow(self.d - self.h.level(c))).sum()
    }
}

impl PartialEq for Tree {
    fn eq(&self, other: &Tree) -> bool {
        self.h == other.h && self.d == other.d && self.root() == other.root()
    }
}
