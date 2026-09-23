//! `wrap_longitude` / `wrap_latitude` (rhealpixdggs `utils.py`). Degrees unless `radians`.

use core::f64::consts::PI;

/// Angle in `[-half, half)`; an input already in range is returned unchanged (no modulo).
pub fn wrap_longitude(lam: f64, radians: bool) -> f64 {
    let half = if radians { PI } else { 180.0 };
    if lam < -half || lam >= half {
        // Python's float `%` has the sign of the divisor, as `rem_euclid` does for a positive one.
        let mut r = lam.rem_euclid(2.0 * half);
        if r >= half {
            r -= 2.0 * half;
        }
        r
    } else {
        lam
    }
}

/// Reflect through the pole into `[-half/2, half/2]`.
pub fn wrap_latitude(phi: f64, radians: bool) -> f64 {
    let phi = wrap_longitude(phi, radians);
    let half = if radians { PI } else { 180.0 };
    if phi.abs() <= half / 2.0 {
        phi
    } else {
        phi - half.copysign(phi)
    }
}
