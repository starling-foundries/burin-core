//! Authalic latitude and radius (rhealpixdggs `utils.py`, `ellipsoids.py`).

use libm::{asin, log, sin, sqrt};

/// Series coefficients of sin(2φ)..sin(12φ) for third flattening `n`: eq. A19 (forward) or A20
/// (inverse) of arXiv 2212.05818. The nesting and the evaluation order are those of the
/// reference, including the three terms of the form `n * K / D` which are `(n * K) / D`.
#[allow(clippy::excessive_precision)]
pub fn auth_lat_coefficients(n: f64, inverse: bool) -> [f64; 6] {
    if !inverse {
        [
            n * (-4.0 / 3.0
                + n * (-4.0 / 45.0
                    + n * (88.0 / 315.0
                        + n * (538.0 / 4725.0 + n * (20824.0 / 467775.0 + n * (-44732.0 / 2837835.0)))))),
            n * (n
                * (34.0 / 45.0
                    + n * (8.0 / 105.0
                        + n * (-2482.0 / 14175.0
                            + n * (-37192.0 / 467775.0 + n * (-12467764.0 / 212837625.0)))))),
            n * (n
                * (n
                    * (-1532.0 / 2835.0
                        + n * (-898.0 / 14175.0
                            + n * (54968.0 / 467775.0 + n * 100320856.0 / 1915538625.0))))),
            n * (n
                * (n
                    * (n
                        * (6007.0 / 14175.0
                            + n * (24496.0 / 467775.0 + n * (-5884124.0 / 70945875.0)))))),
            n * (n * (n * (n * (n * (-23356.0 / 66825.0 + n * (-839792.0 / 19348875.0)))))),
            n * (n * (n * (n * (n * (n * 570284222.0 / 1915538625.0))))),
        ]
    } else {
        [
            n * (4.0 / 3.0
                + n * (4.0 / 45.0
                    + n * (-16.0 / 35.0
                        + n * (-2582.0 / 14175.0
                            + n * (60136.0 / 467775.0 + n * 28112932.0 / 212837625.0))))),
            n * (n
                * (46.0 / 45.0
                    + n * (152.0 / 945.0
                        + n * (-11966.0 / 14175.0
                            + n * (-21016.0 / 51975.0 + n * 251310128.0 / 638512875.0))))),
            n * (n
                * (n
                    * (3044.0 / 2835.0
                        + n * (3802.0 / 14175.0
                            + n * (-94388.0 / 66825.0 + n * (-8797648.0 / 10945935.0)))))),
            n * (n
                * (n
                    * (n
                        * (6059.0 / 4725.0
                            + n * (41072.0 / 93555.0 + n * (-1472637812.0 / 638512875.0)))))),
            n * (n * (n * (n * (n * (768272.0 / 467775.0 + n * 455935736.0 / 638512875.0))))),
            n * (n * (n * (n * (n * (n * 4210684958.0 / 1915538625.0))))),
        ]
    }
}

/// Authalic latitude of geodetic `phi` (radians) on an ellipse of eccentricity `e`, or its
/// inverse. The forward direction uses the closed form for |f| > 1/150 and the series otherwise;
/// the inverse is always the series (so the two are not exact mutual inverses).
pub fn auth_lat(phi: f64, e: f64, inverse: bool) -> f64 {
    if e == 0.0 {
        return phi;
    }
    let s = sqrt(1.0 - e * e);
    let f = 1.0 - s;
    let n = (1.0 - s) / (1.0 + s);
    if !inverse && f.abs() > 1.0 / 150.0 {
        let sp = sin(phi);
        let q = ((1.0 - e * e) * sp) / (1.0 - (e * sp) * (e * sp))
            - (1.0 - e * e) / (2.0 * e) * log((1.0 - e * sp) / (1.0 + e * sp));
        let qp = 1.0 - (1.0 - e * e) / (2.0 * e) * log((1.0 - e) / (1.0 + e));
        let mut ratio = q / qp;
        if ratio.abs() > 1.0 {
            ratio = 1f64.copysign(ratio);
        }
        return asin(ratio);
    }
    let [c2, c4, c6, c8, c10, c12] = auth_lat_coefficients(n, inverse);
    phi + (c2 * sin(2.0 * phi)
        + c4 * sin(4.0 * phi)
        + c6 * sin(6.0 * phi)
        + c8 * sin(8.0 * phi)
        + c10 * sin(10.0 * phi)
        + c12 * sin(12.0 * phi))
}

/// Radius of the authalic sphere of the ellipsoid with major radius `a` and eccentricity `e`.
pub fn auth_rad(a: f64, e: f64) -> f64 {
    if e == 0.0 {
        return a;
    }
    let k = sqrt(0.5 * (1.0 - (1.0 - e * e) / (2.0 * e) * log((1.0 - e) / (1.0 + e))));
    a * k
}

/// An ellipsoid of revolution given by `a` and `f`, as `rhealpixdggs.ellipsoids.Ellipsoid(a=, f=)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipsoid {
    pub a: f64,
    pub f: f64,
    pub e: f64,
    /// Authalic radius; every planar coordinate scales by it.
    pub r_a: f64,
}

impl Ellipsoid {
    pub fn new(a: f64, f: f64) -> Ellipsoid {
        let e = sqrt(f * (2.0 - f));
        Ellipsoid { a, f, e, r_a: auth_rad(a, e) }
    }

    /// WGS84 with `f = 1/298.257223563`.
    pub fn wgs84() -> Ellipsoid {
        Ellipsoid::new(6378137.0, 1.0 / 298.257223563)
    }
}
