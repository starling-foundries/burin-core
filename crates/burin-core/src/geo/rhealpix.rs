//! The rHEALPix projection: HEALPix with the polar triangles combined into two squares
//! (rhealpixdggs `pj_rhealpix.py`), on the unit authalic sphere.
#![allow(clippy::manual_clamp, clippy::neg_cmp_op_on_partial_ord)]

use super::healpix::{healpix_ellipsoid, healpix_ellipsoid_inverse};
use crate::error::{Error, Result};
use core::f64::consts::PI;

const IMAGE_EPS: f64 = 1e-15;
const TRIANGLE_EPS: f64 = 1e-15;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Region {
    Equatorial,
    NorthPolar,
    SouthPolar,
}

/// Anticlockwise rotation by `k` quarter turns (`k` may be negative), exact.
fn rotate(k: i32, (x, y): (f64, f64)) -> (f64, f64) {
    match k.rem_euclid(4) {
        0 => (x, y),
        1 => (-y, x),
        2 => (-x, -y),
        _ => (y, -x),
    }
}

fn sign(x: f64) -> f64 {
    if x > 0.0 {
        1.0
    } else if x < 0.0 {
        -1.0
    } else {
        0.0
    }
}

/// The polar triangle number and region of `(x, y)`: in the HEALPix image (forward) or in the
/// `(north_square, south_square)`-rHEALPix image (inverse).
pub fn triangle(x: f64, y: f64, north_square: i32, south_square: i32, inverse: bool) -> (Option<i32>, Region) {
    let region = if y > PI / 4.0 {
        Region::NorthPolar
    } else if y < -PI / 4.0 {
        Region::SouthPolar
    } else {
        Region::Equatorial
    };
    if region == Region::Equatorial {
        return (None, region);
    }
    let t = if !inverse {
        if x < -PI / 2.0 {
            0
        } else if x < 0.0 {
            1
        } else if x < PI / 2.0 {
            2
        } else {
            3
        }
    } else {
        let eps = TRIANGLE_EPS;
        if region == Region::NorthPolar {
            let ns = north_square;
            let l1 = x - (-3.0 * PI / 4.0 + ((ns - 1) as f64) * PI / 2.0);
            let l2 = -x + (-3.0 * PI / 4.0 + ((ns + 1) as f64) * PI / 2.0);
            if y < l1 - eps && y >= l2 - eps {
                (ns + 1) % 4
            } else if y >= l1 - eps && y > l2 + eps {
                (ns + 2) % 4
            } else if y > l1 + eps && y <= l2 + eps {
                (ns + 3) % 4
            } else {
                ns
            }
        } else {
            let ss = south_square;
            let l1 = x - (-3.0 * PI / 4.0 + ((ss + 1) as f64) * PI / 2.0);
            let l2 = -x + (-3.0 * PI / 4.0 + ((ss - 1) as f64) * PI / 2.0);
            if y <= l1 + eps && y > l2 + eps {
                (ss + 1) % 4
            } else if y < l1 - eps && y <= l2 + eps {
                (ss + 2) % 4
            } else if y >= l1 - eps && y < l2 - eps {
                (ss + 3) % 4
            } else {
                ss
            }
        }
    };
    (Some(t), region)
}

/// Rearrange the polar triangles into squares (forward) or back (inverse).
pub fn combine_triangles(x: f64, y: f64, north_square: i32, south_square: i32, inverse: bool) -> (f64, f64) {
    let ns = north_square.rem_euclid(4);
    let ss = south_square.rem_euclid(4);
    let (c, region) = triangle(x, y, ns, ss, inverse);
    let c = match c {
        None => return (x, y),
        Some(c) => c,
    };
    let tc = (-3.0 * PI / 4.0 + (c as f64) * PI / 2.0, sign(y) * PI / 2.0);
    let (u, k_fwd) = if region == Region::NorthPolar {
        ((-3.0 * PI / 4.0 + (ns as f64) * PI / 2.0, PI / 2.0), c - ns)
    } else {
        ((-3.0 * PI / 4.0 + (ss as f64) * PI / 2.0, -PI / 2.0), -(c - ss))
    };
    if !inverse {
        let (rx, ry) = rotate(k_fwd, (x - tc.0, y - tc.1));
        (rx + u.0, ry + u.1)
    } else {
        let (rx, ry) = rotate(-k_fwd, (x - u.0, y - u.1));
        (rx + tc.0, ry + tc.1)
    }
}

/// Forward on the ellipsoid with eccentricity `e` whose authalic sphere is the unit sphere.
/// `region` is a hint: `Some(Equatorial)` skips the triangle rearrangement.
pub fn rhealpix_ellipsoid(lam: f64, phi: f64, e: f64, ns: i32, ss: i32, region: Option<Region>) -> (f64, f64) {
    let (x, y) = healpix_ellipsoid(lam, phi, e);
    if region != Some(Region::Equatorial) {
        combine_triangles(x, y, ns, ss, false)
    } else {
        (x, y)
    }
}

/// Inverse of [`rhealpix_ellipsoid`]; errors outside the image (1e-15 margin).
pub fn rhealpix_ellipsoid_inverse(x: f64, y: f64, e: f64, ns: i32, ss: i32, region: Option<Region>) -> Result<(f64, f64)> {
    if !in_rhealpix_image(x, y, ns, ss) {
        return Err(Error::OutOfImage(x, y));
    }
    let (x, y) = if region != Some(Region::Equatorial) {
        combine_triangles(x, y, ns, ss, true)
    } else {
        (x, y)
    };
    healpix_ellipsoid_inverse(x, y, e)
}

/// True iff `(x, y)` lies in the image of the rHEALPix projection, with a 1e-15 margin.
pub fn in_rhealpix_image(x: f64, y: f64, north_square: i32, south_square: i32) -> bool {
    let eps = IMAGE_EPS;
    if y.abs() < PI / 4.0 + eps {
        return -PI - eps < x && x < PI + eps;
    }
    if !(y.abs() < 3.0 * PI / 4.0 + eps) {
        return false;
    }
    let square = if y > 0.0 { north_square } else { south_square } as f64;
    -PI + square * PI / 2.0 - eps < x && x < -PI + (square + 1.0) * PI / 2.0 + eps
}
