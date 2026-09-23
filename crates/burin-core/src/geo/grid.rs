//! The planar cell layout of an rHEALPix DGGS and the projection wrapper that carries `lon_0`
//! (rhealpixdggs `dggs.py`, `cell.py`, `projection_wrapper.py`). Angles are degrees.

use super::ellipsoid::Ellipsoid;
use super::rhealpix::{rhealpix_ellipsoid, rhealpix_ellipsoid_inverse};
use super::wrap::{wrap_latitude, wrap_longitude};
use crate::error::Result;
use core::f64::consts::PI;

pub use super::rhealpix::Region;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Shape {
    Quad,
    Cap,
    Dart,
    SkewQuad,
}

/// One rHEALPix grid: ellipsoid, central meridian, polar-square placement, `N_side`.
#[derive(Debug, Clone, PartialEq)]
pub struct Grid {
    pub ell: Ellipsoid,
    pub lon_0: f64,
    pub north_square: i32,
    pub south_square: i32,
    pub n_side: u32,
    /// Upper-left vertices of the six base cells N, O, P, Q, R, S, in metres.
    pub ul0: [(f64, f64); 6],
}

impl Grid {
    pub fn new(ell: Ellipsoid, lon_0: f64, north_square: i32, south_square: i32, n_side: u32) -> Grid {
        let ns = north_square.rem_euclid(4);
        let ss = south_square.rem_euclid(4);
        let unit = [
            (-PI + (ns as f64) * PI / 2.0, 3.0 * PI / 4.0),
            (-PI, PI / 4.0),
            (-PI / 2.0, PI / 4.0),
            (0.0, PI / 4.0),
            (PI / 2.0, PI / 4.0),
            (-PI + (ss as f64) * PI / 2.0, -PI / 4.0),
        ];
        let r = ell.r_a;
        let ul0 = unit.map(|(x, y)| (r * x, r * y));
        Grid { ell, lon_0, north_square: ns, south_square: ss, n_side, ul0 }
    }

    /// `N^(-r)` as the correctly rounded reciprocal of the exact integer `N^r`. (A libm `pow`
    /// is one ulp off this at r = 6; the reciprocal is what every platform's pow and CPython
    /// agree on, and it needs no transcendental function to reproduce.)
    pub fn n_pow_neg(&self, r: u32) -> f64 {
        1.0 / ((self.n_side as u64).pow(r) as f64)
    }

    /// `R_A · (π/2) · N^(-r)`, metres.
    pub fn cell_width(&self, resolution: u32) -> f64 {
        self.ell.r_a * (PI / 2.0) * self.n_pow_neg(resolution)
    }

    /// Exact ellipsoidal cell area at `resolution`, m².
    pub fn cell_area(&self, resolution: u32) -> f64 {
        let w = self.cell_width(resolution);
        8.0 / (3.0 * PI) * (w * w)
    }

    /// Forward projection of `(lon, lat)` degrees to planar metres.
    pub fn forward(&self, lon: f64, lat: f64, region: Option<Region>) -> (f64, f64) {
        let lam = wrap_longitude(lon - self.lon_0, false);
        let phi = wrap_latitude(lat, false);
        let (lam, phi) = (lam * (PI / 180.0), phi * (PI / 180.0));
        let (x, y) = rhealpix_ellipsoid(lam, phi, self.ell.e, self.north_square, self.south_square, region);
        (self.ell.r_a * x, self.ell.r_a * y)
    }

    /// Inverse projection of planar metres to `(lon, lat)` degrees.
    pub fn inverse(&self, x: f64, y: f64, region: Option<Region>) -> Result<(f64, f64)> {
        let (lam, phi) = rhealpix_ellipsoid_inverse(
            x / self.ell.r_a,
            y / self.ell.r_a,
            self.ell.e,
            self.north_square,
            self.south_square,
            region,
        )?;
        let (lam, phi) = (lam * (180.0 / PI), phi * (180.0 / PI));
        Ok((wrap_longitude(lam + self.lon_0, false), wrap_latitude(phi, false)))
    }

    pub fn region(suid: &[u32]) -> Region {
        match suid[0] {
            0 => Region::NorthPolar,
            5 => Region::SouthPolar,
            _ => Region::Equatorial,
        }
    }

    pub fn shape(&self, suid: &[u32]) -> Shape {
        let n = self.n_side;
        if (1..=4).contains(&suid[0]) {
            return Shape::Quad;
        }
        let digits = &suid[1..];
        if digits.is_empty() {
            return Shape::Cap;
        }
        if n % 2 == 1 && digits.iter().all(|&d| d == (n * n - 1) / 2) {
            return Shape::Cap;
        }
        if digits.iter().all(|&d| d % (n + 1) == 0 && d / (n + 1) < n) {
            return Shape::Dart;
        }
        if digits.iter().all(|&d| d >= n - 1 && d % (n - 1) == 0 && d / (n - 1) <= n) {
            return Shape::Dart;
        }
        Shape::SkewQuad
    }

    /// Planar upper-left vertex of a cell, metres.
    pub fn ul_vertex(&self, suid: &[u32]) -> (f64, f64) {
        let (x0, y0) = self.ul0[suid[0] as usize];
        let n = self.n_side;
        let res = (suid.len() - 1) as u32;
        let (mut sx, mut sy): (u64, u64) = (0, 0);
        for (i, &d) in suid[1..].iter().enumerate() {
            let (row, col) = ((d / n) as u64, (d % n) as u64);
            let scale = (n as u64).pow(res - (i as u32 + 1));
            sx += scale * col;
            sy += scale * row;
        }
        let k = self.n_pow_neg(res);
        let (dx, dy) = if res == 0 { (0.0, 0.0) } else { (sx as f64 * k, sy as f64 * k) };
        let w0 = self.cell_width(0);
        (x0 + w0 * dx, y0 - w0 * dy)
    }

    /// Planar nucleus (cell centre), metres.
    pub fn nucleus_planar(&self, suid: &[u32]) -> (f64, f64) {
        let (x, y) = self.ul_vertex(suid);
        let w = self.cell_width((suid.len() - 1) as u32);
        (x + w / 2.0, y - w / 2.0)
    }

    /// Ellipsoidal nucleus, `(lon, lat)` degrees: the cell's reference point.
    pub fn nucleus(&self, suid: &[u32]) -> Result<(f64, f64)> {
        let (x, y) = self.nucleus_planar(suid);
        self.inverse(x, y, None)
    }

    /// The `4n - 4 = 8` boundary points of `Cell.boundary(n=3, plane=False)`, `(lon, lat)`
    /// degrees, as a set (the reference's clockwise start point is not reproduced).
    pub fn boundary3(&self, suid: &[u32]) -> Result<Vec<(f64, f64)>> {
        let (ulx, uly) = self.ul_vertex(suid);
        let w = self.cell_width((suid.len() - 1) as u32);
        let n = 3usize;
        let delta = w / ((n - 1) as f64);
        let region = Grid::region(suid);
        if region == Region::Equatorial {
            // `_quad_boundary`: longitude depends only on x and latitude only on y here, so only
            // the north and west edges are projected and the rest is assembled from them.
            let (x_west, y_north) = (ulx, uly);
            let (x_east, y_south) = (ulx + w, uly - w);
            let hint = Some(Region::Equatorial);
            let mut xs = vec![x_west];
            for j in 0..(n - 1) {
                if j > 0 {
                    xs.push(x_west + delta * (j as f64));
                }
            }
            xs.push(x_east);
            let mut ys = vec![y_north];
            for j in 1..(n - 1) {
                ys.push(y_north - delta * (j as f64));
            }
            ys.push(y_south);
            let lons: Vec<f64> = xs.iter().map(|&x| self.inverse(x, y_north, hint).map(|p| p.0)).collect::<Result<_>>()?;
            let lats: Vec<f64> = ys.iter().map(|&y| self.inverse(x_west, y, hint).map(|p| p.1)).collect::<Result<_>>()?;
            let mut out = Vec::with_capacity(8);
            for &lon in &lons {
                out.push((lon, lats[0]));
                out.push((lon, lats[lats.len() - 1]));
            }
            for &lat in &lats[1..lats.len() - 1] {
                out.push((lons[0], lat));
                out.push((lons[lons.len() - 1], lat));
            }
            return Ok(out);
        }
        let (mut x, mut y) = (ulx, uly);
        let mut pts = vec![(x, y)];
        for (dx, dy) in [(1.0, 0.0), (0.0, -1.0), (-1.0, 0.0), (0.0, 1.0)] {
            for j in 1..n {
                pts.push((x + (j as f64) * delta * dx, y + (j as f64) * delta * dy));
            }
            let last = *pts.last().unwrap();
            x = last.0;
            y = last.1;
        }
        pts.pop();
        pts.iter().map(|&(px, py)| self.inverse(px, py, Some(region))).collect()
    }

    /// The cell boundary as a closed ring of `4(n-1)` points in `(lon, lat)` degrees, clockwise
    /// from the planar upper-left corner, `n` points per edge (`n >= 2`). For display and for
    /// GeoJSON export; the commitment never uses it.
    pub fn cell_ring(&self, suid: &[u32], n: usize) -> Result<Vec<(f64, f64)>> {
        let n = n.max(2);
        let (x0, y0) = self.ul_vertex(suid);
        let w = self.cell_width((suid.len() - 1) as u32);
        let delta = w / ((n - 1) as f64);
        let region = Grid::region(suid);
        let (mut x, mut y) = (x0, y0);
        let mut pts = vec![(x, y)];
        for (dx, dy) in [(1.0, 0.0), (0.0, -1.0), (-1.0, 0.0), (0.0, 1.0)] {
            for j in 1..n {
                pts.push((x + (j as f64) * delta * dx, y + (j as f64) * delta * dy));
            }
            let last = *pts.last().unwrap();
            x = last.0;
            y = last.1;
        }
        let mut ring: Vec<(f64, f64)> = pts.iter().map(|&(px, py)| self.inverse(px, py, Some(region))).collect::<Result<_>>()?;
        // the last point is the first point again (closed ring), projected once
        let first = ring[0];
        *ring.last_mut().unwrap() = first;
        Ok(ring)
    }
}
