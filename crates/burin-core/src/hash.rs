//! The node hash and the two constant ladders.
//!
//! A `Hasher` gives four domain-separated digests. With SHA-256 (the default, id `"sha256"`):
//! ```text
//! leaf_empty         = H(0x00)
//! leaf_full          = H(0x01)
//! node(c_0..c_A-1)   = H(0x02 ‖ c_0 ‖ … ‖ c_A-1)
//! root(b_0..b_B-1)   = H(0x03 ‖ b_0 ‖ … ‖ b_B-1)
//! EMPTY[0] = leaf_empty,  EMPTY[d] = node([EMPTY[d-1]] * A)      (likewise FULL)
//! ```
//! A constant subtree `d` levels above the leaves hashes to its ladder value, so a region held as
//! one cell and the same region spelled out to the leaves hash identically (the roll-up).

use sha2::{Digest as _, Sha256};

pub type Digest = [u8; 32];

pub trait Hasher: Send + Sync {
    fn id(&self) -> &'static str;
    fn leaf_empty(&self) -> Digest;
    fn leaf_full(&self) -> Digest;
    fn node(&self, children: &[Digest]) -> Digest;
    fn root(&self, bases: &[Digest]) -> Digest;
}

#[derive(Debug, Clone, Copy, Default)]
pub struct Sha256Tagged;

fn sha_tagged(tag: u8, parts: &[Digest]) -> Digest {
    let mut h = Sha256::new();
    h.update([tag]);
    for p in parts {
        h.update(p);
    }
    h.finalize().into()
}

impl Hasher for Sha256Tagged {
    fn id(&self) -> &'static str {
        "sha256"
    }
    fn leaf_empty(&self) -> Digest {
        sha_tagged(0x00, &[])
    }
    fn leaf_full(&self) -> Digest {
        sha_tagged(0x01, &[])
    }
    fn node(&self, children: &[Digest]) -> Digest {
        sha_tagged(0x02, children)
    }
    fn root(&self, bases: &[Digest]) -> Digest {
        sha_tagged(0x03, bases)
    }
}

/// The hasher a `hash` id names, if this crate implements it.
pub fn hasher_for(id: &str) -> Option<Box<dyn Hasher>> {
    match id {
        "sha256" => Some(Box::new(Sha256Tagged)),
        _ => None,
    }
}

/// `EMPTY[0..=d_max]` and `FULL[0..=d_max]` for aperture `a`.
#[derive(Debug, Clone)]
pub struct Ladders {
    empty: Vec<Digest>,
    full: Vec<Digest>,
}

impl Ladders {
    pub fn new(h: &dyn Hasher, a: u32, d_max: u32) -> Ladders {
        let mut empty = vec![h.leaf_empty()];
        let mut full = vec![h.leaf_full()];
        for d in 1..=d_max as usize {
            empty.push(h.node(&vec![empty[d - 1]; a as usize]));
            full.push(h.node(&vec![full[d - 1]; a as usize]));
        }
        Ladders { empty, full }
    }
    pub fn empty(&self, d: u32) -> Digest {
        self.empty[d as usize]
    }
    pub fn full(&self, d: u32) -> Digest {
        self.full[d as usize]
    }
    pub fn depth(&self) -> u32 {
        (self.empty.len() - 1) as u32
    }
}

/// The hashing context every tree, opening and transcript is evaluated against.
pub struct Ctx {
    pub hasher: Box<dyn Hasher>,
    pub a: u32,
    pub b: u32,
    pub ladders: Ladders,
}

impl Ctx {
    pub fn new(hasher: Box<dyn Hasher>, a: u32, b: u32, d_max: u32) -> Ctx {
        let ladders = Ladders::new(hasher.as_ref(), a, d_max);
        Ctx { hasher, a, b, ladders }
    }
    pub fn sha256(a: u32, b: u32, d_max: u32) -> Ctx {
        Ctx::new(Box::new(Sha256Tagged), a, b, d_max)
    }
}

pub fn hex(d: &Digest) -> String {
    ::hex::encode(d)
}

/// 64 lowercase hex digits, as `hex` writes them: a digest has one spelling on the wire.
pub fn unhex(s: &str) -> Option<Digest> {
    if s.bytes().any(|c| c.is_ascii_uppercase()) {
        return None;
    }
    let v = ::hex::decode(s).ok()?;
    v.try_into().ok()
}
