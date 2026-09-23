//! The extension module `burin._burin`; `python/burin` is the package around it.

use burin_core::hierarchy::{cid_to_suid as core_cid_to_suid, suid_to_cid as core_suid_to_cid, Cid};
use burin_core::index::Index as CoreIndex;
use burin_core::opening::{open_path, OpeningRecord};
use burin_core::polyfill as pf;
use burin_core::profile::Profile as CoreProfile;
use burin_core::setops::{self, Op, SetOpProof};
use burin_core::tree::Tree as CoreTree;
use burin_core::zone;
use burin_core::zone_data::{self, Presence};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyList};
use serde_json::Value;
use std::sync::Arc;

mod vector;

fn err(e: burin_core::Error) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// A Python object (dict/list/str) → JSON value, via the `json` module for anything but str.
fn to_json(obj: &Bound<'_, PyAny>) -> PyResult<Value> {
    let s: String = if let Ok(s) = obj.extract::<String>() {
        s
    } else {
        let json = obj.py().import("json")?;
        json.call_method1("dumps", (obj,))?.extract()?
    };
    serde_json::from_str(&s).map_err(|e| PyValueError::new_err(format!("not JSON: {e}")))
}

/// Cell ids from a uint64 numpy array (read directly) or any sequence of ints.
fn cid_list(obj: &Bound<'_, PyAny>) -> PyResult<Vec<u64>> {
    use numpy::PyArrayMethods;
    if let Ok(a) = obj.cast::<numpy::PyArray1<u64>>() {
        if let Ok(v) = a.to_vec() {
            return Ok(v);
        }
    }
    obj.extract()
}

fn from_json<'py>(py: Python<'py>, v: &Value) -> PyResult<Bound<'py, PyAny>> {
    let json = py.import("json")?;
    json.call_method1("loads", (v.to_string(),))
}

/// The parameters a cell identifier and a root are read against. Integers only; the identity
/// (`id_hex`) is derived from them and never transmitted.
#[pyclass(name = "Profile", module = "burin", frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct Profile {
    inner: CoreProfile,
}

#[pymethods]
impl Profile {
    #[new]
    #[pyo3(signature = (hash="sha256", aperture=9, n_base=6, lon_0_udeg=50_000_000, north_square=0, south_square=0, a_um=6_378_137_000_000, inv_f_nano=298_257_223_563, tick_us=1, epoch_us=0))]
    #[allow(clippy::too_many_arguments)]
    fn new(hash: &str, aperture: u32, n_base: u32, lon_0_udeg: i64, north_square: u8, south_square: u8, a_um: u64, inv_f_nano: u64, tick_us: u64, epoch_us: i64) -> PyResult<Profile> {
        let inner = CoreProfile { hash: hash.to_string(), aperture, n_base, lon_0_udeg, north_square, south_square, a_um, inv_f_nano, tick_us, epoch_us };
        inner.validate().map_err(err)?;
        Ok(Profile { inner })
    }

    /// The OGC API-DGGS Annex B rHEALPix grid (lon_0 = 50), SHA-256. The default.
    #[staticmethod]
    fn ogc() -> Profile {
        Profile { inner: CoreProfile::ogc() }
    }

    /// rHEALPix with a 0° prime meridian (lon_0 = 0), SHA-256.
    #[staticmethod]
    fn burin_1() -> Profile {
        Profile { inner: CoreProfile::burin_1() }
    }

    #[staticmethod]
    fn from_dict(d: &Bound<'_, PyAny>) -> PyResult<Profile> {
        let v = to_json(d)?;
        let inner: CoreProfile = serde_json::from_value(v).map_err(|e| PyValueError::new_err(format!("bad profile record: {e}")))?;
        inner.validate().map_err(err)?;
        Ok(Profile { inner })
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        from_json(py, &serde_json::to_value(&self.inner).unwrap())
    }

    fn id_hex(&self) -> String {
        self.inner.id_hex()
    }

    fn describe(&self) -> String {
        self.inner.describe()
    }

    /// Exact ellipsoidal cell area at `level`, m².
    fn area_m2(&self, level: u32) -> f64 {
        self.inner.area_m2(level)
    }

    #[getter]
    fn hash(&self) -> String {
        self.inner.hash.clone()
    }
    #[getter]
    fn lon_0(&self) -> f64 {
        self.inner.lon_0()
    }
    #[getter]
    fn aperture(&self) -> u32 {
        self.inner.aperture
    }
    #[getter]
    fn n_base(&self) -> u32 {
        self.inner.n_base
    }

    fn __repr__(&self) -> String {
        format!("Profile({})", self.inner.describe())
    }

    fn __eq__(&self, other: &Profile) -> bool {
        self.inner == other.inner
    }
}

fn profile_or_default(p: Option<&Profile>) -> CoreProfile {
    p.map(|p| p.inner.clone()).unwrap_or_default()
}

/// An immutable canonical coverage tree: a set of equal-area cells and its 32-byte root.
#[pyclass(name = "Tree", module = "burin", frozen)]
pub struct Tree {
    inner: CoreTree,
    profile: CoreProfile,
}

#[pymethods]
impl Tree {
    /// Cover every cid (any levels up to `depth`, any order, duplicates allowed).
    #[staticmethod]
    #[pyo3(signature = (cids, depth, profile=None))]
    fn from_cells(cids: &Bound<'_, PyAny>, depth: u32, profile: Option<&Profile>) -> PyResult<Tree> {
        let cids = cid_list(cids)?;
        let p = profile_or_default(profile);
        let ctx = Arc::new(p.ctx(depth).map_err(err)?);
        let inner = CoreTree::from_cells(p.hierarchy(), depth, cids, ctx).map_err(err)?;
        Ok(Tree { inner, profile: p })
    }

    /// The canonical coverage of a GeoJSON Polygon/MultiPolygon at `resolution` (nucleus rule).
    #[staticmethod]
    #[pyo3(signature = (geojson, resolution, profile=None))]
    fn from_geojson(geojson: &Bound<'_, PyAny>, resolution: u32, profile: Option<&Profile>) -> PyResult<Tree> {
        let p = profile_or_default(profile);
        let geom = pf::geometry_from_geojson(&to_json(geojson)?).map_err(err)?;
        let ctx = Arc::new(p.ctx(resolution).map_err(err)?);
        let inner = pf::cover(&geom, resolution, &p, ctx).map_err(err)?;
        Ok(Tree { inner, profile: p })
    }

    #[getter]
    fn root_hex(&self) -> String {
        self.inner.root_hex()
    }
    #[getter]
    fn depth(&self) -> u32 {
        self.inner.d
    }
    #[getter]
    fn profile(&self) -> Profile {
        Profile { inner: self.profile.clone() }
    }

    /// The canonical (coarsest) covered cids, sorted.
    fn cells(&self) -> Vec<u64> {
        self.inner.cells()
    }
    /// Every covered leaf at `depth`, sorted.
    fn leaves(&self) -> Vec<u64> {
        self.inner.leaves()
    }
    fn leaf_count(&self) -> u64 {
        self.inner.leaf_count()
    }
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    /// Whether the cell `cid` (its whole subtree) is covered.
    fn covers(&self, cid: u64) -> PyResult<bool> {
        self.inner.covers(cid).map_err(err)
    }

    /// The covered leaves as six base-cell rasters at `depth`, `uint8`, shape `(6, n, n)` with
    /// `n = 3**depth`, base-major then scanline: `numpy.frombuffer(t.raster(), numpy.uint8)`.
    fn raster<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        Ok(PyBytes::new(py, &zone::raster(&self.inner).map_err(err)?))
    }

    /// The canonical tree of a raster laid out as `raster` returns it (any nonzero byte is
    /// covered), e.g. `Tree.from_raster(mask.astype(numpy.uint8).tobytes(), depth)`.
    #[staticmethod]
    #[pyo3(signature = (data, depth, profile=None))]
    fn from_raster(data: &[u8], depth: u32, profile: Option<&Profile>) -> PyResult<Tree> {
        let p = profile_or_default(profile);
        let ctx = Arc::new(p.ctx(depth).map_err(err)?);
        let inner = zone::from_raster(p.hierarchy(), depth, data, ctx).map_err(err)?;
        Ok(Tree { inner, profile: p })
    }

    /// OGC API - DGGS zone data (DGGS-JSON) for the sub-zones of `zone` at the tree's depth, in
    /// the DGGRS's sub-zone order: 1 where covered, null elsewhere. The profile must be a
    /// registered DGGRS (the default is).
    #[pyo3(signature = (zone, field="coverage"))]
    fn to_dggs_json<'py>(&self, py: Python<'py>, zone: u64, field: &str) -> PyResult<Bound<'py, PyAny>> {
        from_json(py, &zone_data::to_dggs_json(&self.inner, &self.profile, zone, field).map_err(err)?)
    }

    /// The coverage a DGGS-JSON document describes: the sub-zones of its zone at one depth whose
    /// value in `field` is present (not null; also not zero if `nonzero`). `field` and `depth`
    /// may be omitted when the document has only one.
    #[staticmethod]
    #[pyo3(signature = (doc, field=None, depth=None, nonzero=false, profile=None))]
    fn from_dggs_json(doc: &Bound<'_, PyAny>, field: Option<&str>, depth: Option<u32>, nonzero: bool, profile: Option<&Profile>) -> PyResult<Tree> {
        let p = profile_or_default(profile);
        let presence = if nonzero { Presence::NonZero } else { Presence::NonNull };
        let inner = zone_data::from_dggs_json(&to_json(doc)?, &p, field, depth, presence).map_err(err)?;
        Ok(Tree { inner, profile: p })
    }

    fn union(&self, other: &Tree) -> PyResult<Tree> {
        Ok(Tree { inner: setops::union(&self.inner, &other.inner).map_err(err)?, profile: self.profile.clone() })
    }
    fn intersect(&self, other: &Tree) -> PyResult<Tree> {
        Ok(Tree { inner: setops::intersect(&self.inner, &other.inner).map_err(err)?, profile: self.profile.clone() })
    }
    fn difference(&self, other: &Tree) -> PyResult<Tree> {
        Ok(Tree { inner: setops::difference(&self.inner, &other.inner).map_err(err)?, profile: self.profile.clone() })
    }

    /// Maximal subtrees where two coverages differ: `[{"cid", "suid", "hash_a", "hash_b"}, ...]`.
    fn divergence<'py>(&self, py: Python<'py>, other: &Tree) -> PyResult<Bound<'py, PyList>> {
        let out = PyList::empty(py);
        for d in setops::divergence(&self.inner, &other.inner).map_err(err)? {
            let item = PyDict::new(py);
            item.set_item("cid", d.cid)?;
            item.set_item("suid", core_cid_to_suid(d.cid).map_err(err)?)?;
            item.set_item("hash_a", burin_core::hash::hex(&d.hash_a))?;
            item.set_item("hash_b", burin_core::hash::hex(&d.hash_b))?;
            out.append(item)?;
        }
        Ok(out)
    }

    /// A self-describing opening record for `cid` (membership or non-membership), or `None`
    /// if `cid` is a partial node.
    fn open<'py>(&self, py: Python<'py>, cid: u64) -> PyResult<Option<Bound<'py, PyAny>>> {
        match open_path(&self.inner, cid).map_err(err)? {
            None => Ok(None),
            Some(op) => Ok(Some(from_json(py, &OpeningRecord::new(&self.inner, op).to_json())?)),
        }
    }

    /// A transcript proving `self ⊕ other` for `op` in `"union" | "intersect" | "difference"`.
    fn prove<'py>(&self, py: Python<'py>, other: &Tree, op: &str) -> PyResult<Bound<'py, PyAny>> {
        let op = Op::parse(op).ok_or_else(|| PyValueError::new_err(format!("unknown op {op:?}")))?;
        from_json(py, &setops::prove(&self.inner, &other.inner, op).map_err(err)?.to_json())
    }

    /// GeoJSON of the canonical cells (see `cells_geojson`).
    #[pyo3(signature = (n=3))]
    fn geojson<'py>(&self, py: Python<'py>, n: usize) -> PyResult<Bound<'py, PyAny>> {
        let prof = Profile { inner: self.profile.clone() };
        cells_geojson(py, self.inner.cells(), Some(&prof), n)
    }

    fn __eq__(&self, other: &Tree) -> bool {
        self.inner == other.inner
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

/// A coverage index: register trees under keys, then ask by tree. The catalog-side half.
#[pyclass(name = "Index", module = "burin")]
pub struct Index {
    inner: CoreIndex,
}

#[pymethods]
impl Index {
    #[new]
    fn new() -> Index {
        Index { inner: CoreIndex::new() }
    }

    /// Register `key` under `tree`. Trees with the same root are stored once.
    fn add(&mut self, key: &str, tree: &Tree) -> PyResult<()> {
        self.inner.add(key, &tree.inner).map_err(err)
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    #[getter]
    fn distinct_roots(&self) -> usize {
        self.inner.distinct_roots()
    }

    /// Keys whose coverage meets `tree`: `[{"key", "root", "relation", "fraction_of_query",
    /// "fraction_of_item", "shared_cells"}, ...]`, best coverage of the query first.
    fn query<'py>(&self, py: Python<'py>, tree: &Tree) -> PyResult<Bound<'py, PyList>> {
        let out = PyList::empty(py);
        for h in self.inner.query(&tree.inner).map_err(err)? {
            let d = PyDict::new(py);
            d.set_item("key", h.key)?;
            d.set_item("root", burin_core::hash::hex(&h.root))?;
            d.set_item("relation", h.relation.name())?;
            d.set_item("fraction_of_query", h.fraction_of_query)?;
            d.set_item("fraction_of_item", h.fraction_of_item)?;
            d.set_item("shared_cells", h.shared_cells)?;
            out.append(d)?;
        }
        Ok(out)
    }

    /// `{cid: count}` for every cell at `level` touched by an indexed key's coverage.
    fn attention<'py>(&self, py: Python<'py>, level: u32) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        for (cid, n) in self.inner.attention(level).map_err(err)? {
            d.set_item(cid, n)?;
        }
        Ok(d)
    }

    /// `[(key, root_hex), ...]`.
    fn roots(&self) -> Vec<(String, String)> {
        self.inner.roots()
    }
}

/// Sorted depth-`resolution` cids whose nucleus the GeoJSON polygon contains.
#[pyfunction]
#[pyo3(signature = (geojson, resolution, profile=None))]
fn polyfill(geojson: &Bound<'_, PyAny>, resolution: u32, profile: Option<&Profile>) -> PyResult<Vec<u64>> {
    let geom = pf::geometry_from_geojson(&to_json(geojson)?).map_err(err)?;
    pf::polyfill(&geom, resolution, &profile_or_default(profile)).map_err(err)
}

/// Hex root of the canonical coverage of a GeoJSON polygon at `resolution`.
#[pyfunction]
#[pyo3(signature = (geojson, resolution, profile=None))]
fn fingerprint_polygon(geojson: &Bound<'_, PyAny>, resolution: u32, profile: Option<&Profile>) -> PyResult<String> {
    let d = pf::fingerprint_geojson(&to_json(geojson)?, resolution, &profile_or_default(profile)).map_err(err)?;
    Ok(burin_core::hash::hex(&d))
}

/// Verify an opening record (as returned by `Tree.open`). `False` on any defect; never raises.
#[pyfunction]
fn verify_opening(record: &Bound<'_, PyAny>) -> PyResult<bool> {
    let v = to_json(record)?;
    Ok(OpeningRecord::from_json(&v).map(|r| r.verify()).unwrap_or(false))
}

/// Verify a set-operation transcript (as returned by `Tree.prove`).
#[pyfunction]
fn verify_setop(proof: &Bound<'_, PyAny>) -> PyResult<bool> {
    let v = to_json(proof)?;
    Ok(SetOpProof::from_json(&v).map(|p| p.verify()).unwrap_or(false))
}

#[pyfunction]
#[pyo3(signature = (level, profile=None))]
fn cell_area_m2(level: u32, profile: Option<&Profile>) -> f64 {
    profile_or_default(profile).area_m2(level)
}

/// GeoJSON FeatureCollection of cell polygons for display: counterclockwise rings, `n` points
/// per edge on polar cells (equatorial cells are their four corners), rings that cross the
/// antimeridian shifted east of 180° so tile maps draw them whole, and each polar cap closed
/// through its pole. Properties: `cid`, `suid`, `level`; feature `id` is the cid as a string.
#[pyfunction]
#[pyo3(signature = (cids, profile=None, n=3))]
fn cells_geojson<'py>(py: Python<'py>, cids: Vec<u64>, profile: Option<&Profile>, n: usize) -> PyResult<Bound<'py, PyAny>> {
    let p = profile_or_default(profile);
    let grid = p.grid();
    let h = p.hierarchy();
    let mut features = Vec::with_capacity(cids.len());
    for cid in cids {
        let path = h.path(cid).map_err(err)?;
        let ring = grid.cell_polygon(&path, n).map_err(err)?;
        let coords: Vec<[f64; 2]> = ring.iter().map(|&(x, y)| [x, y]).collect();
        features.push(serde_json::json!({
            "type": "Feature", "id": cid.to_string(),
            "properties": {"cid": cid, "suid": core_cid_to_suid(cid).map_err(err)?, "level": h.level(cid)},
            "geometry": {"type": "Polygon", "coordinates": [coords]},
        }));
    }
    from_json(py, &serde_json::json!({"type": "FeatureCollection", "features": features}))
}

#[pyfunction]
fn suid_to_cid(suid: &str) -> PyResult<Cid> {
    core_suid_to_cid(suid).map_err(err)
}

#[pyfunction]
fn cid_to_suid(cid: u64) -> PyResult<String> {
    core_cid_to_suid(cid).map_err(err)
}

/// `(base, level, row, col)`: where `cid` sits in its base cell (row 0 at the top).
#[pyfunction]
#[pyo3(signature = (cid, profile=None))]
fn position(cid: u64, profile: Option<&Profile>) -> PyResult<(u32, u32, u64, u64)> {
    let p = zone::position(&profile_or_default(profile).hierarchy(), cid).map_err(err)?;
    Ok((p.base, p.level, p.row, p.col))
}

/// The cid at `(row, col)` of base cell `base` at `level`.
#[pyfunction]
#[pyo3(signature = (base, level, row, col, profile=None))]
fn cell_at(base: u32, level: u32, row: u64, col: u64, profile: Option<&Profile>) -> PyResult<Cid> {
    zone::cell_at(&profile_or_default(profile).hierarchy(), base, level, row, col).map_err(err)
}

/// The scanline index of `cid` among the sub-zones of `ancestor`.
#[pyfunction]
#[pyo3(signature = (ancestor, cid, profile=None))]
fn subzone_index(ancestor: u64, cid: u64, profile: Option<&Profile>) -> PyResult<u64> {
    zone::subzone_index(&profile_or_default(profile).hierarchy(), ancestor, cid).map_err(err)
}

/// The sub-zone of `ancestor` at relative `depth` with scanline index `index`.
#[pyfunction]
#[pyo3(signature = (ancestor, depth, index, profile=None))]
fn subzone_at(ancestor: u64, depth: u32, index: u64, profile: Option<&Profile>) -> PyResult<Cid> {
    zone::subzone_at(&profile_or_default(profile).hierarchy(), ancestor, depth, index).map_err(err)
}

/// Every sub-zone of `ancestor` at relative `depth`, in scanline order (the OGC API - DGGS
/// sub-zone order for rHEALPix).
#[pyfunction]
#[pyo3(signature = (ancestor, depth, profile=None))]
fn subzones(ancestor: u64, depth: u32, profile: Option<&Profile>) -> PyResult<Vec<Cid>> {
    zone::subzones(&profile_or_default(profile).hierarchy(), ancestor, depth).map_err(err)
}

/// The four edge neighbours of `cid`, keyed `up`, `right`, `down`, `left` in the cell's planar
/// frame; base-cell edges are crossed as the grid folds.
#[pyfunction]
#[pyo3(signature = (cid, profile=None))]
fn neighbours<'py>(py: Python<'py>, cid: u64, profile: Option<&Profile>) -> PyResult<Bound<'py, PyDict>> {
    let nb = zone::neighbours(&profile_or_default(profile), cid).map_err(err)?;
    let d = PyDict::new(py);
    for (dir, c) in zone::DIRECTIONS.iter().zip(nb) {
        d.set_item(dir.name(), c)?;
    }
    Ok(d)
}

/// The raster index of each cid (all at one level), `int64` little-endian:
/// `values_raster.ravel()[numpy.frombuffer(raster_index(cids), numpy.int64)] = values`.
#[pyfunction]
#[pyo3(signature = (cids, profile=None))]
fn raster_index<'py>(py: Python<'py>, cids: Vec<u64>, profile: Option<&Profile>) -> PyResult<Bound<'py, PyBytes>> {
    let idx = zone::raster_index(&profile_or_default(profile).hierarchy(), &cids).map_err(err)?;
    let bytes: Vec<u8> = idx.iter().flat_map(|x| (*x as i64).to_le_bytes()).collect();
    Ok(PyBytes::new(py, &bytes))
}

/// The cid at every position of the raster at `level`, `uint64` little-endian, shape `(6, n, n)`.
#[pyfunction]
#[pyo3(signature = (level, profile=None))]
fn raster_cells<'py>(py: Python<'py>, level: u32, profile: Option<&Profile>) -> PyResult<Bound<'py, PyBytes>> {
    let cells = zone::raster_cells(&profile_or_default(profile).hierarchy(), level).map_err(err)?;
    let bytes: Vec<u8> = cells.iter().flat_map(|x| x.to_le_bytes()).collect();
    Ok(PyBytes::new(py, &bytes))
}

/// A gather index padding each base-cell raster at `level` with `width` halo cells, `int64`
/// little-endian, shape `(6, n + 2w, n + 2w)`; each entry indexes the flat `(6, n, n)` raster, or
/// is -1 beyond a base cell's corners: `padded = values.ravel()[idx]` (mask where `idx < 0`).
#[pyfunction]
#[pyo3(signature = (level, width=1, profile=None))]
fn halo_index<'py>(py: Python<'py>, level: u32, width: u32, profile: Option<&Profile>) -> PyResult<Bound<'py, PyBytes>> {
    let idx = zone::halo_index(&profile_or_default(profile), level, width).map_err(err)?;
    let bytes: Vec<u8> = idx.iter().flat_map(|x| x.to_le_bytes()).collect();
    Ok(PyBytes::new(py, &bytes))
}

#[pymodule]
#[pyo3(name = "_burin")]
fn burin(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Profile>()?;
    m.add_class::<Tree>()?;
    m.add_class::<Index>()?;
    m.add_function(wrap_pyfunction!(polyfill, m)?)?;
    m.add_function(wrap_pyfunction!(fingerprint_polygon, m)?)?;
    m.add_function(wrap_pyfunction!(verify_opening, m)?)?;
    m.add_function(wrap_pyfunction!(verify_setop, m)?)?;
    m.add_function(wrap_pyfunction!(cell_area_m2, m)?)?;
    m.add_function(wrap_pyfunction!(cells_geojson, m)?)?;
    m.add_function(wrap_pyfunction!(suid_to_cid, m)?)?;
    m.add_function(wrap_pyfunction!(cid_to_suid, m)?)?;
    m.add_function(wrap_pyfunction!(position, m)?)?;
    m.add_function(wrap_pyfunction!(cell_at, m)?)?;
    m.add_function(wrap_pyfunction!(subzone_index, m)?)?;
    m.add_function(wrap_pyfunction!(subzone_at, m)?)?;
    m.add_function(wrap_pyfunction!(subzones, m)?)?;
    m.add_function(wrap_pyfunction!(neighbours, m)?)?;
    m.add_function(wrap_pyfunction!(halo_index, m)?)?;
    m.add_function(wrap_pyfunction!(vector::cells_from_lonlat, m)?)?;
    m.add_function(wrap_pyfunction!(vector::cells_to_lonlat, m)?)?;
    m.add_function(wrap_pyfunction!(vector::cell_boundaries, m)?)?;
    m.add_function(wrap_pyfunction!(vector::cell_neighbours, m)?)?;
    m.add_function(wrap_pyfunction!(raster_index, m)?)?;
    m.add_function(wrap_pyfunction!(raster_cells, m)?)?;
    m.add("OGC_RHEALPIX", Profile::ogc())?;
    m.add("BURIN_1", Profile::burin_1())?;
    m.add("MAX_RESOLUTION", pf::MAX_RESOLUTION)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
