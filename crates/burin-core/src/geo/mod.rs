//! The minimal rHEALPix geometry: authalic latitude, the HEALPix and rHEALPix projections
//! (both directions), and the planar cell layout. A port of the corresponding functions in
//! `rhealpixdggs-py` 0.8.6 (Raichev, Gibb, Law and contributors; MIT), kept line-for-line so that
//! the two agree to the last bit wherever the platform transcendental kernels do.
//!
//! Every transcendental call goes through the `libm` crate so the same bits come out on every
//! target, wasm included.

pub mod ellipsoid;
pub mod grid;
pub mod healpix;
pub mod rhealpix;
pub mod wrap;

pub use ellipsoid::Ellipsoid;
pub use grid::{Grid, Region, Shape};
