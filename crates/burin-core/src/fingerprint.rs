//! The fingerprint: a root with everything it is read against, and the one-call checks.
//!
//! ```text
//! burin:<profile>:space:<depth>:<root>      burin:ogc-rhealpix:space:11:f0a5…
//! burin:<profile>:time:<root>               burin:ogc-rhealpix:time:0e85…
//! ```
//! `<profile>` is a registered name or the profile's 64-hex id. A check answers a question in the
//! world's terms (a point, an instant, a cell): `Ok(true)` or `Ok(false)` when the proof settles
//! it, and `Err(Error::InvalidProof)` with the reason when it proves nothing.

use crate::error::{invalid, Error, Result};
use crate::hash::{hex, hasher_for, unhex, Ctx, Digest};
use crate::hierarchy::{Cid, Hierarchy};
use crate::opening::{Claim, OpeningRecord};
use crate::profile::{Profile, BURIN_1, OGC_RHEALPIX};
use crate::setops::{empty_root, Op, SetOpProof, SetRelation};
use crate::time::{self, TIME, TIME_DEPTH};
use crate::tree::Tree;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Axis {
    Space,
    Time,
}

impl Axis {
    pub fn name(self) -> &'static str {
        match self {
            Axis::Space => "space",
            Axis::Time => "time",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fingerprint {
    pub profile: Profile,
    pub axis: Axis,
    pub depth: u32,
    pub root: Digest,
}

fn refused<T>(why: impl Into<String>) -> Result<T> {
    Err(Error::InvalidProof(why.into()))
}

/// The registered profiles, by the names fingerprints use.
pub fn registered_profiles() -> [(&'static str, Profile); 2] {
    [("ogc-rhealpix", OGC_RHEALPIX.with_hash("sha256")), ("burin-1", BURIN_1.with_hash("sha256"))]
}

impl Fingerprint {
    /// The fingerprint of `tree`, read against `profile`.
    pub fn of(tree: &Tree, profile: &Profile) -> Result<Fingerprint> {
        profile.validate()?;
        if tree.ctx.hasher.id() != profile.hash {
            return invalid(format!("the tree is hashed with {:?}, the profile names {:?}", tree.ctx.hasher.id(), profile.hash));
        }
        let axis = if tree.h == TIME {
            Axis::Time
        } else if tree.h == profile.hierarchy() {
            Axis::Space
        } else {
            return invalid("the tree is neither the profile's grid nor time");
        };
        Ok(Fingerprint { profile: profile.clone(), axis, depth: tree.d, root: tree.root() })
    }

    pub fn hierarchy(&self) -> Hierarchy {
        match self.axis {
            Axis::Space => self.profile.hierarchy(),
            Axis::Time => TIME,
        }
    }

    /// The name of a registered profile, or its 64-hex id.
    pub fn profile_token(&self) -> String {
        registered_profiles().into_iter().find(|(_, p)| *p == self.profile).map_or_else(|| self.profile.id_hex(), |(n, _)| n.to_string())
    }

    /// The text form, which has one spelling.
    pub fn to_text(&self) -> String {
        match self.axis {
            Axis::Space => format!("burin:{}:space:{}:{}", self.profile_token(), self.depth, hex(&self.root)),
            Axis::Time => format!("burin:{}:time:{}", self.profile_token(), hex(&self.root)),
        }
    }

    /// Read the text form. A profile named by its id must be among `known`.
    pub fn parse(text: &str, known: &[Profile]) -> Result<Fingerprint> {
        let parts: Vec<&str> = text.split(':').collect();
        let (token, axis, depth, root) = match parts[..] {
            ["burin", p, "space", d, r] => (p, Axis::Space, d, r),
            ["burin", p, "time", r] => (p, Axis::Time, "61", r),
            _ => return invalid(format!("{text:?} is not a fingerprint (burin:<profile>:space:<depth>:<root> or burin:<profile>:time:<root>)")),
        };
        let profile = match registered_profiles().into_iter().find(|(n, _)| *n == token) {
            Some((_, p)) => p,
            None => match known.iter().find(|p| p.id_hex() == token) {
                Some(p) => p.clone(),
                None => return invalid(format!("profile {token:?} is neither registered nor among the profiles given")),
            },
        };
        let depth: u32 = depth.parse().map_err(|_| Error::Invalid(format!("depth {depth:?} is not a number")))?;
        let root = unhex(root).ok_or_else(|| Error::Invalid(format!("root {root:?} is not 64 lowercase hex digits")))?;
        let fp = Fingerprint { profile, axis, depth, root };
        let max = match axis {
            Axis::Space => fp.profile.hierarchy().max_level(),
            Axis::Time => TIME_DEPTH,
        };
        if fp.depth > max {
            return invalid(format!("depth {} is deeper than the deepest level, {max}", fp.depth));
        }
        if fp.to_text() != text {
            return invalid(format!("{text:?} is not spelled canonically; it is {:?}", fp.to_text()));
        }
        Ok(fp)
    }

    fn check_proof_shape(&self, proof: &OpeningRecord) -> Result<()> {
        let h = self.hierarchy();
        if proof.hash != self.profile.hash || (proof.a, proof.b, proof.d) != (h.a, h.b, self.depth) {
            return refused(format!(
                "the proof is for {} H({},{}) depth {}; this fingerprint is {} H({},{}) depth {}",
                proof.hash, proof.a, proof.b, proof.d, self.profile.hash, h.a, h.b, self.depth
            ));
        }
        if proof.root != self.root {
            return refused("the proof is against another fingerprint");
        }
        if !proof.verify() {
            return refused("the proof does not verify");
        }
        Ok(())
    }

    /// Whether `cell` is covered, as `proof` shows.
    pub fn check_cell(&self, proof: &OpeningRecord, cell: Cid) -> Result<bool> {
        self.hierarchy().check(cell, Some(self.depth))?;
        self.check_proof_shape(proof)?;
        match (proof.cid(), proof.opening.terminal()) {
            (Some(c), Some(claim)) if c == cell => Ok(claim == Claim::Full),
            (Some(c), Some(_)) => refused(format!("the proof is about cell {c}, not {cell}")),
            _ => refused("the proof does not open exactly one cell"),
        }
    }

    /// Whether the point `(lon, lat)`, in degrees, is in the set: its cell at this fingerprint's
    /// depth, on this fingerprint's grid, is covered.
    pub fn check_point(&self, proof: &OpeningRecord, lon: f64, lat: f64) -> Result<bool> {
        if self.axis != Axis::Space {
            return invalid("a time fingerprint has no points; check an instant");
        }
        self.check_cell(proof, crate::zone::cell_from_point(&self.profile, lon, lat, self.depth)?)
    }

    /// Whether the POSIX instant `posix_us` is in the set: its tick is covered.
    pub fn check_instant(&self, proof: &OpeningRecord, posix_us: i64) -> Result<bool> {
        if self.axis != Axis::Time {
            return invalid("a space fingerprint has no instants; check a point");
        }
        self.check_cell(proof, time::tick_cell(time::tick(&self.profile, posix_us)?)?)
    }

    /// Whether `relation` holds between this set and `other`'s, as the transcript shows. `Equals`
    /// takes no transcript.
    pub fn check_relation(&self, relation: SetRelation, other: &Fingerprint, transcript: Option<&SetOpProof>) -> Result<bool> {
        if (&self.profile, self.axis, self.depth) != (&other.profile, other.axis, other.depth) {
            return invalid("the two fingerprints are of different profiles, axes or depths; no relation between them is defined");
        }
        let p = match (relation, transcript) {
            (SetRelation::Equals, None) => return Ok(self.root == other.root),
            (SetRelation::Equals, Some(_)) => return refused("equals takes no transcript; compare the roots"),
            (_, None) => return refused(format!("{} needs a transcript", relation.name())),
            (_, Some(p)) => p,
        };
        let (op, a, b) = match relation {
            SetRelation::Disjoint | SetRelation::Intersects => (Op::Intersect, &self.root, &other.root),
            SetRelation::Within => (Op::Difference, &self.root, &other.root),
            SetRelation::Contains => (Op::Difference, &other.root, &self.root),
            SetRelation::Equals => unreachable!("handled above"),
        };
        let h = self.hierarchy();
        if p.hash != self.profile.hash || (p.a, p.b, p.d) != (h.a, h.b, self.depth) {
            return refused("the transcript is for another hierarchy or depth");
        }
        if p.op != op || p.root_a != *a || p.root_b != *b {
            return refused(format!("the transcript is not the {} of these two sets that {} needs", p.op.name(), relation.name()));
        }
        if !p.verify() {
            return refused("the transcript does not verify");
        }
        let empty = empty_root(&Ctx::new(hasher_for(&p.hash).expect("verified"), p.a, p.b, p.d), p.d);
        Ok(match relation {
            SetRelation::Intersects => p.root_c != empty,
            _ => p.root_c == empty,
        })
    }
}
