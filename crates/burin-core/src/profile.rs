//! The profile: every parameter a cell identifier and a root are read against, stored as
//! integers so that its identity is reproducible in any language with no float formatting.
//!
//! ```text
//! id = SHA-256( "BURIN-PROFILE-1" ‖ u8 len(hash) ‖ hash ‖ u32 A ‖ u32 B ‖ i64 lon_0_udeg
//!               ‖ u8 north_square ‖ u8 south_square ‖ u64 a_um ‖ u64 inv_f_nano
//!               ‖ u64 tick_us ‖ i64 epoch_us )          all big-endian
//! ```
//! The id is derived, never transmitted: a wire record carries the parameters and the receiver
//! recomputes it, so a record and its identity cannot drift apart.

use crate::error::{invalid, Result};
use crate::geo::{Ellipsoid, Grid};
use crate::hash::{hasher_for, Ctx, Digest};
use crate::hierarchy::Hierarchy;
use serde::{Deserialize, Serialize};
use sha2::{Digest as _, Sha256};

pub const LON_SCALE: i64 = 1_000_000; // degrees -> micro-degrees
pub const AXIS_SCALE: u64 = 1_000_000; // metres -> micro-metres
pub const INVF_SCALE: u64 = 1_000_000_000; // 1/f -> nano

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Profile {
    pub hash: String,
    pub aperture: u32,
    pub n_base: u32,
    pub lon_0_udeg: i64,
    pub north_square: u8,
    pub south_square: u8,
    pub a_um: u64,
    pub inv_f_nano: u64,
    pub tick_us: u64,
    pub epoch_us: i64,
}

/// The OGC API-DGGS Annex B rHEALPix DGGRS: WGS84, `lon_0 = 50`, `N_side = 3`. The default.
pub const OGC_RHEALPIX: Profile = Profile {
    hash: String::new(),
    aperture: 9,
    n_base: 6,
    lon_0_udeg: 50_000_000,
    north_square: 0,
    south_square: 0,
    a_um: 6_378_137_000_000,
    inv_f_nano: 298_257_223_563,
    tick_us: 1,
    epoch_us: 0,
};

/// rHEALPix with a 0° prime meridian: as above with `lon_0 = 0`.
pub const BURIN_1: Profile = Profile {
    hash: String::new(),
    aperture: 9,
    n_base: 6,
    lon_0_udeg: 0,
    north_square: 0,
    south_square: 0,
    a_um: 6_378_137_000_000,
    inv_f_nano: 298_257_223_563,
    tick_us: 1,
    epoch_us: 0,
};

impl Default for Profile {
    fn default() -> Profile {
        OGC_RHEALPIX.with_hash("sha256")
    }
}

impl Profile {
    /// The same parameters under a named hash. Constants above carry an empty hash id so they can
    /// be `const`; use this (or `Default`) to obtain a usable profile.
    pub fn with_hash(&self, hash: &str) -> Profile {
        Profile { hash: hash.to_string(), ..self.clone() }
    }

    pub fn ogc() -> Profile {
        Profile::default()
    }

    pub fn burin_1() -> Profile {
        BURIN_1.with_hash("sha256")
    }

    pub fn validate(&self) -> Result<()> {
        Hierarchy::new(self.aperture, self.n_base)?;
        if hasher_for(&self.hash).is_none() {
            return invalid(format!("unknown hash id {:?}", self.hash));
        }
        if self.lon_0_udeg.abs() > 180 * LON_SCALE {
            return invalid(format!("lon_0 {} micro-degrees outside [-180, 180] degrees", self.lon_0_udeg));
        }
        if self.north_square > 3 || self.south_square > 3 {
            return invalid("north_square/south_square must be in 0..3");
        }
        if self.a_um == 0 || self.inv_f_nano == 0 {
            return invalid("ellipsoid a and inv_f must be positive");
        }
        if self.tick_us == 0 {
            return invalid("tick_us must be >= 1");
        }
        let n = self.n_side();
        if n * n != self.aperture {
            return invalid(format!("aperture {} is not a square; no rHEALPix N_side", self.aperture));
        }
        Ok(())
    }

    pub fn n_side(&self) -> u32 {
        (self.aperture as f64).sqrt() as u32
    }

    pub fn hierarchy(&self) -> Hierarchy {
        Hierarchy { a: self.aperture, b: self.n_base }
    }

    pub fn lon_0(&self) -> f64 {
        self.lon_0_udeg as f64 / LON_SCALE as f64
    }

    pub fn a(&self) -> f64 {
        self.a_um as f64 / AXIS_SCALE as f64
    }

    pub fn inv_f(&self) -> f64 {
        self.inv_f_nano as f64 / INVF_SCALE as f64
    }

    pub fn ellipsoid(&self) -> Ellipsoid {
        Ellipsoid::new(self.a(), 1.0 / self.inv_f())
    }

    pub fn grid(&self) -> Grid {
        Grid::new(self.ellipsoid(), self.lon_0(), self.north_square as i32, self.south_square as i32, self.n_side())
    }

    /// Exact ellipsoidal cell area at `level`, m².
    pub fn area_m2(&self, level: u32) -> f64 {
        self.grid().cell_area(level)
    }

    /// A hashing context for trees of depth `d_max` under this profile.
    pub fn ctx(&self, d_max: u32) -> Result<Ctx> {
        self.validate()?;
        let h = hasher_for(&self.hash).expect("validated");
        Ok(Ctx::new(h, self.aperture, self.n_base, d_max))
    }

    /// The canonical bytes the id is the SHA-256 of.
    pub fn id_preimage(&self) -> Vec<u8> {
        let mut v = Vec::with_capacity(64);
        v.extend_from_slice(b"BURIN-PROFILE-1");
        v.push(self.hash.len() as u8);
        v.extend_from_slice(self.hash.as_bytes());
        v.extend_from_slice(&self.aperture.to_be_bytes());
        v.extend_from_slice(&self.n_base.to_be_bytes());
        v.extend_from_slice(&self.lon_0_udeg.to_be_bytes());
        v.push(self.north_square);
        v.push(self.south_square);
        v.extend_from_slice(&self.a_um.to_be_bytes());
        v.extend_from_slice(&self.inv_f_nano.to_be_bytes());
        v.extend_from_slice(&self.tick_us.to_be_bytes());
        v.extend_from_slice(&self.epoch_us.to_be_bytes());
        v
    }

    pub fn id(&self) -> Digest {
        Sha256::digest(self.id_preimage()).into()
    }

    pub fn id_hex(&self) -> String {
        crate::hash::hex(&self.id())
    }

    pub fn describe(&self) -> String {
        format!(
            "H({},{}) rhealpix lon_0={} ns={} ss={} a={} 1/f={} hash={} tick={}us epoch={}",
            self.aperture,
            self.n_base,
            self.lon_0(),
            self.north_square,
            self.south_square,
            self.a(),
            self.inv_f(),
            self.hash,
            self.tick_us,
            self.epoch_us
        )
    }
}
