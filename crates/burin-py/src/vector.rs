//! Array functions: one point or cell per element, numpy arrays in and out, the work spread over
//! threads with the interpreter released. The public, shape-preserving wrappers are in
//! `python/burin/_vector.py`; these take and return contiguous 1-D arrays.

use crate::{err, profile_or_default, Profile};
use burin_core::polyfill::MAX_RESOLUTION;
use burin_core::zone;
use numpy::prelude::*;
use numpy::{PyArray1, PyArray2, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;

/// Below this many elements a call stays on the calling thread.
const PARALLEL_MIN: usize = 1 << 14;

/// Two arrays returned together.
type Pair<'py, A, B> = PyResult<(Bound<'py, A>, Bound<'py, B>)>;

/// `f(0..n)` in order, on the rayon pool (`nthreads == 0`), on a pool of `nthreads` threads, or
/// serially for small inputs and `nthreads == 1`. The first error, by position, is returned.
fn run<T: Send>(n: usize, nthreads: usize, f: impl Fn(usize) -> Result<T, String> + Sync + Send) -> Result<Vec<T>, String> {
    if n < PARALLEL_MIN || nthreads == 1 {
        return (0..n).map(f).collect();
    }
    if nthreads == 0 {
        return (0..n).into_par_iter().map(f).collect();
    }
    let pool = rayon::ThreadPoolBuilder::new().num_threads(nthreads).build().map_err(|e| e.to_string())?;
    pool.install(|| (0..n).into_par_iter().map(&f).collect())
}

fn slice<'a, T: numpy::Element>(a: &'a PyReadonlyArray1<'_, T>) -> PyResult<&'a [T]> {
    a.as_slice().map_err(|e| PyValueError::new_err(format!("array is not contiguous: {e}")))
}

fn check_level(level: u32) -> PyResult<()> {
    if level > MAX_RESOLUTION {
        return Err(PyValueError::new_err(format!("level {level} is above the maximum {MAX_RESOLUTION}")));
    }
    Ok(())
}

/// The cell at `level` containing each `(lon[i], lat[i])`, degrees.
#[pyfunction]
#[pyo3(name = "_cells_from_lonlat", signature = (lon, lat, level, profile=None, nthreads=0))]
pub(crate) fn cells_from_lonlat<'py>(
    py: Python<'py>,
    lon: PyReadonlyArray1<'py, f64>,
    lat: PyReadonlyArray1<'py, f64>,
    level: u32,
    profile: Option<&Profile>,
    nthreads: usize,
) -> PyResult<Bound<'py, PyArray1<u64>>> {
    check_level(level)?;
    let p = profile_or_default(profile);
    p.validate().map_err(err)?;
    let (lon, lat) = (slice(&lon)?, slice(&lat)?);
    if lon.len() != lat.len() {
        return Err(PyValueError::new_err(format!("{} longitudes but {} latitudes", lon.len(), lat.len())));
    }
    let (grid, h) = (p.grid(), p.hierarchy());
    let cells = py.detach(|| {
        run(lon.len(), nthreads, |i| match grid.cell_from_lonlat(lon[i], lat[i], level) {
            Some((b, r, c)) => zone::cell_at(&h, b, level, r, c).map_err(|e| e.to_string()),
            None => Err(format!("point {i} ({}, {}) is not on the grid", lon[i], lat[i])),
        })
    });
    Ok(cells.map_err(PyValueError::new_err)?.into_pyarray(py))
}

/// The nucleus of each cell, `(lon, lat)` degrees.
#[pyfunction]
#[pyo3(name = "_cells_to_lonlat", signature = (cids, profile=None, nthreads=0))]
pub(crate) fn cells_to_lonlat<'py>(
    py: Python<'py>,
    cids: PyReadonlyArray1<'py, u64>,
    profile: Option<&Profile>,
    nthreads: usize,
) -> Pair<'py, PyArray1<f64>, PyArray1<f64>> {
    let p = profile_or_default(profile);
    p.validate().map_err(err)?;
    let cids = slice(&cids)?;
    let (grid, h) = (p.grid(), p.hierarchy());
    let points = py.detach(|| {
        run(cids.len(), nthreads, |i| {
            let path = h.path(cids[i]).map_err(|e| format!("cell {i}: {e}"))?;
            grid.nucleus(&path).map_err(|e| format!("cell {i}: {e}"))
        })
    });
    let (lon, lat): (Vec<f64>, Vec<f64>) = points.map_err(PyValueError::new_err)?.into_iter().unzip();
    Ok((lon.into_pyarray(py), lat.into_pyarray(py)))
}

/// Each cell's polygon (`Grid::cell_polygon`, `n` points per polar edge) as ragged arrays:
/// `coords` of shape `(M, 2)` holding every ring in turn, and `offsets` of length `N + 1` so that
/// cell `i` is `coords[offsets[i]:offsets[i + 1]]`.
#[pyfunction]
#[pyo3(name = "_cell_boundaries", signature = (cids, n=5, profile=None, nthreads=0))]
pub(crate) fn cell_boundaries<'py>(
    py: Python<'py>,
    cids: PyReadonlyArray1<'py, u64>,
    n: usize,
    profile: Option<&Profile>,
    nthreads: usize,
) -> Pair<'py, PyArray2<f64>, PyArray1<i64>> {
    let p = profile_or_default(profile);
    p.validate().map_err(err)?;
    let cids = slice(&cids)?;
    let (grid, h) = (p.grid(), p.hierarchy());
    let rings = py.detach(|| {
        run(cids.len(), nthreads, |i| {
            let path = h.path(cids[i]).map_err(|e| format!("cell {i}: {e}"))?;
            grid.cell_polygon(&path, n).map_err(|e| format!("cell {i}: {e}"))
        })
    });
    let rings = rings.map_err(PyValueError::new_err)?;
    let mut offsets = Vec::with_capacity(rings.len() + 1);
    offsets.push(0i64);
    let mut coords = Vec::with_capacity(rings.iter().map(|r| 2 * r.len()).sum());
    for ring in &rings {
        coords.extend(ring.iter().flat_map(|&(x, y)| [x, y]));
        offsets.push(offsets[offsets.len() - 1] + ring.len() as i64);
    }
    let points = coords.len() / 2;
    Ok((coords.into_pyarray(py).reshape([points, 2])?, offsets.into_pyarray(py)))
}

/// The four edge neighbours of each cell, columns up, right, down, left.
#[pyfunction]
#[pyo3(name = "_cell_neighbours", signature = (cids, profile=None, nthreads=0))]
pub(crate) fn cell_neighbours<'py>(
    py: Python<'py>,
    cids: PyReadonlyArray1<'py, u64>,
    profile: Option<&Profile>,
    nthreads: usize,
) -> PyResult<Bound<'py, PyArray2<u64>>> {
    let p = profile_or_default(profile);
    p.validate().map_err(err)?;
    let cids = slice(&cids)?;
    let rows = py.detach(|| run(cids.len(), nthreads, |i| zone::neighbours(&p, cids[i]).map_err(|e| format!("cell {i}: {e}"))));
    let flat: Vec<u64> = rows.map_err(PyValueError::new_err)?.into_iter().flatten().collect();
    let n = flat.len() / 4;
    flat.into_pyarray(py).reshape([n, 4])
}
