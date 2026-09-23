//! The HEALPix projection of the unit authalic sphere (rhealpixdggs `pj_healpix.py`).
#![allow(clippy::manual_clamp, clippy::neg_cmp_op_on_partial_ord)]

use super::ellipsoid::auth_lat;
use crate::error::{Error, Result};
use core::f64::consts::PI;
use libm::{asin, floor, sin, sqrt};

const IMAGE_EPS: f64 = 1e-10;

fn sign(x: f64) -> f64 {
    if x > 0.0 {
        1.0
    } else if x < 0.0 {
        -1.0
    } else {
        0.0
    }
}

/// The polar cap (0-3) a longitude falls in, clamped so rounding past ±π stays in the outer cap.
fn cap_number(lam: f64) -> f64 {
    let c = floor(2.0 * lam / PI + 2.0);
    c.max(0.0).min(3.0)
}

/// Forward, `(lam, phi)` in radians with `-π <= lam < π`, `-π/2 <= phi <= π/2`.
pub fn healpix_sphere(lam: f64, phi: f64) -> (f64, f64) {
    let phi0 = asin(2.0 / 3.0);
    if phi.abs() <= phi0 {
        (lam, 3.0 * PI / 8.0 * sin(phi))
    } else {
        let sigma = sqrt(3.0 * (1.0 - sin(phi).abs()));
        let c = cap_number(lam);
        let lamc = -3.0 * PI / 4.0 + (PI / 2.0) * c;
        let x = lamc + (lam - lamc) * sigma;
        let y = sign(phi) * PI / 4.0 * (2.0 - sigma);
        (x, y)
    }
}

/// Inverse of [`healpix_sphere`]. Errors outside the image (with a 1e-10 margin).
pub fn healpix_sphere_inverse(x: f64, y: f64) -> Result<(f64, f64)> {
    if !in_healpix_image(x, y) {
        return Err(Error::OutOfImage(x, y));
    }
    let y0 = PI / 4.0;
    if y.abs() <= y0 {
        Ok((x, asin(8.0 * y / (3.0 * PI))))
    } else if y.abs() < PI / 2.0 {
        let c = cap_number(x);
        let xc = -3.0 * PI / 4.0 + (PI / 2.0) * c;
        let tau = 2.0 - 4.0 * y.abs() / PI;
        let mut lam = xc + (x - xc) / tau;
        let phi = sign(y) * asin(1.0 - tau * tau / 3.0);
        if lam < -PI {
            lam = -PI;
        } else if lam > PI {
            lam = PI;
        }
        Ok((lam, phi))
    } else {
        // The poles: longitude is -π by convention.
        Ok((-PI, sign(y) * PI / 2.0))
    }
}

pub fn healpix_ellipsoid(lam: f64, phi: f64, e: f64) -> (f64, f64) {
    let beta = auth_lat(phi, e, false);
    healpix_sphere(lam, beta)
}

pub fn healpix_ellipsoid_inverse(x: f64, y: f64, e: f64) -> Result<(f64, f64)> {
    let (lam, beta) = healpix_sphere_inverse(x, y)?;
    Ok((lam, auth_lat(beta, e, true)))
}

/// True iff `(x, y)` lies in the image of the HEALPix projection, with a 1e-10 margin.
pub fn in_healpix_image(x: f64, y: f64) -> bool {
    let eps = IMAGE_EPS;
    let abs_y = y.abs();
    if !(x.abs() < PI + eps && abs_y < PI / 2.0 + eps) {
        return false;
    }
    if abs_y < PI / 4.0 + eps {
        return true;
    }
    // Truncation (Python `int()`), not floor, as in the reference.
    let c = ((2.0 * x / PI + 2.0).trunc() as i64).max(0).min(3) as f64;
    let x_c = -3.0 * PI / 4.0 + (PI / 2.0) * c;
    (x - x_c).abs() + abs_y < PI / 2.0 + eps
}
