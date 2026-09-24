//! Time (SPEC §9–§12): ticks, cells, intervals and Allen's relations. The public wrappers are in
//! `python/burin/time.py`.

use crate::{err, profile_or_default, Profile};
use burin_core::time as core;
use numpy::prelude::*;
use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Two arrays returned together.
type Pair<'py, A, B> = PyResult<(Bound<'py, A>, Bound<'py, B>)>;

fn slice<'a, T: numpy::Element>(a: &'a PyReadonlyArray1<'_, T>) -> PyResult<&'a [T]> {
    a.as_slice().map_err(|e| PyValueError::new_err(format!("array is not contiguous: {e}")))
}

/// The tick containing each POSIX instant (µs).
#[pyfunction]
#[pyo3(name = "_ticks", signature = (posix_us, profile=None))]
pub(crate) fn ticks<'py>(py: Python<'py>, posix_us: PyReadonlyArray1<'py, i64>, profile: Option<&Profile>) -> PyResult<Bound<'py, PyArray1<u64>>> {
    let p = profile_or_default(profile);
    let t = slice(&posix_us)?;
    let out = py.detach(|| t.iter().enumerate().map(|(i, &x)| core::tick(&p, x).map_err(|e| format!("instant {i}: {e}"))).collect::<Result<Vec<_>, _>>());
    Ok(out.map_err(PyValueError::new_err)?.into_pyarray(py))
}

/// The POSIX instant (µs) at which each tick starts; tick 2^61 is the end of the line.
#[pyfunction]
#[pyo3(name = "_tick_starts", signature = (ticks, profile=None))]
pub(crate) fn tick_starts<'py>(py: Python<'py>, ticks: PyReadonlyArray1<'py, u64>, profile: Option<&Profile>) -> PyResult<Bound<'py, PyArray1<i64>>> {
    let p = profile_or_default(profile);
    let k = slice(&ticks)?;
    let out = py.detach(|| k.iter().enumerate().map(|(i, &x)| core::tick_start(&p, x).map_err(|e| format!("tick {i}: {e}"))).collect::<Result<Vec<_>, _>>());
    Ok(out.map_err(PyValueError::new_err)?.into_pyarray(py))
}

/// The compact cells of the tick interval `[lo, hi)`, in time order (the HINT boundary walk).
#[pyfunction]
#[pyo3(name = "_interval_cells")]
pub(crate) fn interval_cells(lo: u64, hi: u64) -> PyResult<Vec<u64>> {
    core::interval_cells(lo, hi).map_err(err)
}

/// The level-`level` cells whose midpoint tick lies in `[lo, hi)`.
#[pyfunction]
#[pyo3(name = "_coarse_cells")]
pub(crate) fn coarse_cells(lo: u64, hi: u64, level: u32) -> PyResult<Vec<u64>> {
    core::coarse_cells(lo, hi, level).map_err(err)
}

/// The ticks `[lo, hi)` of each time cell.
#[pyfunction]
#[pyo3(name = "_cell_ticks")]
pub(crate) fn cell_ticks<'py>(py: Python<'py>, cids: PyReadonlyArray1<'py, u64>) -> Pair<'py, PyArray1<u64>, PyArray1<u64>> {
    let c = slice(&cids)?;
    let spans = c.iter().enumerate().map(|(i, &x)| core::cell_ticks(x).map_err(|e| PyValueError::new_err(format!("cell {i}: {e}")))).collect::<PyResult<Vec<_>>>()?;
    let (lo, hi): (Vec<u64>, Vec<u64>) = spans.into_iter().unzip();
    Ok((lo.into_pyarray(py), hi.into_pyarray(py)))
}

/// Allen's relation of the tick interval `a` to `b`, by name.
#[pyfunction]
#[pyo3(name = "_allen")]
pub(crate) fn allen(a: (u64, u64), b: (u64, u64)) -> PyResult<&'static str> {
    Ok(core::allen(a, b).map_err(err)?.name())
}
