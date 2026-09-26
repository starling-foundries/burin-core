//! `Fingerprint`, `Proof`, `Transcript` and `InvalidProof`: verification that answers a question
//! (is this point, instant or cell in the set?) in one call.

use crate::{err, from_json, to_json, Profile};
use burin_core::fingerprint::Fingerprint as CoreFingerprint;
use burin_core::opening::{Claim, OpeningRecord};
use burin_core::setops::{SetOpProof, SetRelation};
use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

create_exception!(burin, InvalidProof, PyValueError, "The proof proves nothing about the question asked.");

/// Core errors from a check: a proof that proves nothing is `InvalidProof`, a bad question `ValueError`.
pub(crate) fn check_err(e: burin_core::Error) -> PyErr {
    match e {
        burin_core::Error::InvalidProof(m) => InvalidProof::new_err(m),
        other => err(other),
    }
}

/// POSIX microseconds of one instant, read as `burin.time.posix_us` reads it.
pub(crate) fn posix_us(instant: &Bound<'_, PyAny>) -> PyResult<i64> {
    let us = instant.py().import("burin.time")?.getattr("posix_us")?.call1((instant,))?;
    us.call_method0("__int__")?.extract()
}

/// A proof that one cell is, or is not, in a set: an opening record (SPEC §5).
#[pyclass(name = "Proof", module = "burin", frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct Proof {
    pub(crate) inner: OpeningRecord,
}

#[pymethods]
impl Proof {
    /// Read a proof from its JSON form (a dict or a string); `InvalidProof` if it is not one.
    #[staticmethod]
    fn from_json(obj: &Bound<'_, PyAny>) -> PyResult<Proof> {
        let inner = OpeningRecord::from_json(&to_json(obj)?).map_err(|e| InvalidProof::new_err(format!("not a proof: {e}")))?;
        Ok(Proof { inner })
    }
    fn to_json<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        from_json(py, &self.inner.to_json())
    }
    /// The cell this proof says it is about. Unverified: only a `Fingerprint.check_*` settles it.
    #[getter]
    fn claimed_cell(&self) -> Option<u64> {
        self.inner.cid()
    }
    /// Whether this proof says its cell is covered. Unverified, as `claimed_cell`.
    #[getter]
    fn claimed_covered(&self) -> Option<bool> {
        self.inner.opening.terminal().map(|c| c == Claim::Full)
    }
    fn __eq__(&self, other: &Proof) -> bool {
        self.inner == other.inner
    }
    fn __repr__(&self) -> String {
        let show = |v: Option<String>| v.unwrap_or_else(|| "None".into());
        let covered = self.claimed_covered().map(|c| if c { "True".to_string() } else { "False".to_string() });
        format!("Proof(claimed_cell={}, claimed_covered={})", show(self.claimed_cell().map(|c| c.to_string())), show(covered))
    }
}

/// A set-operation transcript (SPEC §5), the evidence for a relation between two sets.
#[pyclass(name = "Transcript", module = "burin", frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct Transcript {
    pub(crate) inner: SetOpProof,
}

#[pymethods]
impl Transcript {
    /// Read a transcript from its JSON form; `InvalidProof` if it is not one.
    #[staticmethod]
    fn from_json(obj: &Bound<'_, PyAny>) -> PyResult<Transcript> {
        let inner = SetOpProof::from_json(&to_json(obj)?).map_err(|e| InvalidProof::new_err(format!("not a transcript: {e}")))?;
        Ok(Transcript { inner })
    }
    fn to_json<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        from_json(py, &self.inner.to_json())
    }
    fn __eq__(&self, other: &Transcript) -> bool {
        self.inner == other.inner
    }
}

fn proof_of(obj: &Bound<'_, PyAny>) -> PyResult<OpeningRecord> {
    match obj.cast::<Proof>() {
        Ok(p) => Ok(p.get().inner.clone()),
        Err(_) => Ok(Proof::from_json(obj)?.inner),
    }
}

fn transcript_of(obj: Option<&Bound<'_, PyAny>>) -> PyResult<Option<SetOpProof>> {
    match obj {
        None => Ok(None),
        Some(o) => match o.cast::<Transcript>() {
            Ok(t) => Ok(Some(t.get().inner.clone())),
            Err(_) => Ok(Some(Transcript::from_json(o)?.inner)),
        },
    }
}

pub(crate) fn relation(name: &str) -> PyResult<SetRelation> {
    SetRelation::parse(name).ok_or_else(|| PyValueError::new_err(format!("unknown relation {name:?}; one of equals, disjoint, intersects, within, contains")))
}

/// What a published root is read against: the profile, the axis (space or time), the depth.
/// Its text form is what is published: `burin:ogc-rhealpix:space:11:<root>` or
/// `burin:ogc-rhealpix:time:<root>`. A check shows the proof is consistent with the set; it says
/// nothing of whether the point or instant asked about is true.
#[pyclass(name = "Fingerprint", module = "burin", frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct Fingerprint {
    pub(crate) inner: CoreFingerprint,
}

#[pymethods]
impl Fingerprint {
    /// Read the text form. A profile that is not registered is named by its id and must be in `profiles`.
    #[staticmethod]
    #[pyo3(signature = (text, profiles=None))]
    fn parse(text: &str, profiles: Option<Vec<PyRef<'_, Profile>>>) -> PyResult<Fingerprint> {
        let known: Vec<_> = profiles.unwrap_or_default().iter().map(|p| p.inner.clone()).collect();
        Ok(Fingerprint { inner: CoreFingerprint::parse(text, &known).map_err(err)? })
    }
    #[getter]
    fn profile(&self) -> Profile {
        Profile { inner: self.inner.profile.clone() }
    }
    #[getter]
    fn axis(&self) -> &'static str {
        self.inner.axis.name()
    }
    #[getter]
    fn depth(&self) -> u32 {
        self.inner.depth
    }
    #[getter]
    fn root_hex(&self) -> String {
        burin_core::hash::hex(&self.inner.root)
    }

    /// Whether `cell` is in the set. `InvalidProof` unless `proof` settles exactly this cell
    /// against exactly this fingerprint.
    fn check_cell(&self, proof: &Bound<'_, PyAny>, cell: u64) -> PyResult<bool> {
        self.inner.check_cell(&proof_of(proof)?, cell).map_err(check_err)
    }
    /// Whether the point `(lon, lat)` is in the set, its cell taken on this fingerprint's grid at
    /// its depth. `InvalidProof` unless `proof` is about that cell.
    fn check_point(&self, proof: &Bound<'_, PyAny>, lon: f64, lat: f64) -> PyResult<bool> {
        self.inner.check_point(&proof_of(proof)?, lon, lat).map_err(check_err)
    }
    /// Whether `instant` (datetime64, ISO 8601 text or POSIX µs) is in the set. `InvalidProof`
    /// unless `proof` is about its tick.
    fn check_instant(&self, proof: &Bound<'_, PyAny>, instant: &Bound<'_, PyAny>) -> PyResult<bool> {
        self.inner.check_instant(&proof_of(proof)?, posix_us(instant)?).map_err(check_err)
    }
    /// Whether `relation` holds between this set and `other`'s, as `transcript` shows (none for
    /// `equals`). `InvalidProof` unless the transcript is the one this relation needs.
    #[pyo3(signature = (relation, other, transcript=None))]
    fn check_relation(&self, relation: &str, other: &Fingerprint, transcript: Option<&Bound<'_, PyAny>>) -> PyResult<bool> {
        let t = transcript_of(transcript)?;
        self.inner.check_relation(self::relation(relation)?, &other.inner, t.as_ref()).map_err(check_err)
    }

    fn __str__(&self) -> String {
        self.inner.to_text()
    }
    fn __repr__(&self) -> String {
        format!("Fingerprint({:?})", self.inner.to_text())
    }
    fn __eq__(&self, other: &Fingerprint) -> bool {
        self.inner == other.inner
    }
    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut h = std::collections::hash_map::DefaultHasher::new();
        self.inner.to_text().hash(&mut h);
        h.finish()
    }
}

/// The fingerprint of a tree built under `profile`.
pub(crate) fn fingerprint_of(tree: &burin_core::tree::Tree, profile: &burin_core::profile::Profile) -> PyResult<Fingerprint> {
    Ok(Fingerprint { inner: CoreFingerprint::of(tree, profile).map_err(err)? })
}
