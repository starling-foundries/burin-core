"""The Python/Rust boundary and the geometry as a whole: inputs as they arrive from files and
threads, errors that name the offending element, JSON strictness, and cell polygons that tile
the globe exactly once."""
import json
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

import burin as bc

Q453 = bc.suid_to_cid("Q453")


def _points(n, seed=0):
    rng = np.random.default_rng(seed)
    return rng.uniform(-180, 180, n), np.degrees(np.arcsin(rng.uniform(-1, 1, n)))


def test_byte_order_layout_and_strides_do_not_change_results():
    lon, lat = _points(4000)
    want = bc.cells_from_lonlat(lon, lat, 11)
    assert np.array_equal(bc.cells_from_lonlat(lon.astype(">f8"), lat.astype(">f8"), 11), want), "big-endian, as read from files"
    assert np.array_equal(bc.cells_from_lonlat(lon.astype("f4"), lat.astype("f4"), 5), bc.cells_from_lonlat(lon.astype("f4").astype("f8"), lat.astype("f4").astype("f8"), 5))
    grid_lon, grid_lat = np.asfortranarray(lon.reshape(40, 100)), np.asfortranarray(lat.reshape(40, 100))
    assert np.array_equal(bc.cells_from_lonlat(grid_lon, grid_lat, 11), want.reshape(40, 100)), "Fortran order"
    assert np.array_equal(bc.cells_from_lonlat(lon[::3], lat[::3], 11), want[::3]), "strided"
    ids = want.reshape(40, 100)
    lon_c, lat_c = bc.cells_to_lonlat(ids.astype(">u8"))
    assert lon_c.shape == (40, 100) and np.array_equal(lon_c, bc.cells_to_lonlat(ids)[0])
    assert np.array_equal(bc.cell_neighbours(np.asfortranarray(ids)), bc.cell_neighbours(ids))


def test_scalars_stay_scalars():
    lon, lat = bc.cells_to_lonlat(np.uint64(Q453))
    assert lon.shape == lat.shape == ()
    assert bc.cells_from_lonlat(float(lon), float(lat), 3).shape == ()
    assert int(bc.cells_from_lonlat(float(lon), float(lat), 3)) == Q453
    assert bc.cell_neighbours(Q453).shape == (4,)
    assert bc.zoom_to(Q453, 3, 1).shape == ()


def test_ids_are_integers_or_nothing():
    for bad in ([1.0], np.array([Q453], dtype="f8"), [True], [Q453, None], [2**70], ["Q453"]):
        with pytest.raises(TypeError):
            bc.as_cell_ids(bad)
    with pytest.raises(ValueError, match="non-negative"):
        bc.as_cell_ids(np.array([Q453, -1]))
    for good in ([Q453], np.array([Q453], "i8"), np.array([Q453], "u4"), np.array([Q453], ">u8")):
        assert bc.as_cell_ids(good).tolist() == [Q453]
    gap = 15 * 9**3  # after the last cell at level 3, before the first at level 4
    for not_a_cell in (0, 8, gap, 9**5 - 1, 15 * 9**15):
        with pytest.raises(ValueError, match="not a cell id"):
            bc.cell_levels([not_a_cell])


def test_errors_name_the_offending_element():
    lon, lat = np.zeros(10), np.zeros(10)
    lat[7] = np.nan
    with pytest.raises(ValueError, match="point 7"):
        bc.cells_from_lonlat(lon, lat, 4)
    lat[7], lat[3] = 0.0, 95.0
    with pytest.raises(ValueError, match="point 3"):
        bc.cells_from_lonlat(lon, lat, 4)
    with pytest.raises(ValueError, match=r"cell 2"):
        bc.cells_to_lonlat(np.array([Q453, Q453, 5], dtype=np.uint64))
    with pytest.raises(ValueError):
        bc.cells_from_lonlat(np.zeros(3), np.zeros(4), 4)
    with pytest.raises(OverflowError):
        bc.cells_from_lonlat([0.0], [0.0], 3, nthreads=-1)


def test_latitudes_are_bounded_and_longitudes_wrap():
    for lat in (90.0, -90.0):
        bc.cells_from_lonlat(0.0, lat, 6)
    for lat in (90.000001, -90.5, 180.0):
        with pytest.raises(ValueError, match="not on the grid"):
            bc.cells_from_lonlat(0.0, lat, 6)
    lon = np.arange(-180, 180, 0.25)
    lat = np.full_like(lon, 12.5)
    here = bc.cells_from_lonlat(lon, lat, 13)
    assert np.array_equal(bc.cells_from_lonlat(lon + 360, lat, 13), here)
    assert np.array_equal(bc.cells_from_lonlat(lon - 720, lat, 13), here)


def test_levels():
    assert bc.cells_from_lonlat(1.0, 1.0, np.int64(4)) == bc.cells_from_lonlat(1.0, 1.0, 4.0)
    for bad in (True, 4.5, -1, 16):
        with pytest.raises(ValueError, match="level"):
            bc.cells_from_lonlat(1.0, 1.0, bad)


def test_concurrent_callers_get_the_same_answers():
    lon, lat = _points(60_000, seed=3)
    want = bc.cells_from_lonlat(lon, lat, 12, nthreads=1)
    centres = bc.cells_to_lonlat(want, nthreads=1)
    def work(k):
        cells = bc.cells_from_lonlat(lon, lat, 12, nthreads=(0, 2, 5)[k % 3])
        return np.array_equal(cells, want) and all(np.array_equal(a, b) for a, b in zip(bc.cells_to_lonlat(cells), centres))
    with ThreadPoolExecutor(8) as pool:
        assert all(pool.map(work, range(16)))


def test_zoom_to_guards_its_size_and_nests():
    q = bc.suid_to_cid("Q")
    with pytest.raises(ValueError, match="MAX_ZOOM_CELLS"):
        bc.zoom_to([q], 0, 15)
    kids = bc.zoom_to([q], 0, 3)[0]
    assert np.array_equal(kids, np.sort(kids)), "nested (cid) order is sorted"
    assert bc.Tree.from_cells(kids, 3).cells() == [q], "the descendants are exactly the cell"
    assert np.array_equal(bc.zoom_to(kids, 3, 0), np.full(kids.size, q, dtype=np.uint64))


def test_trees_accept_every_integer_container_and_the_deepest_ids():
    deepest = 15 * 9**18 - 1  # beyond 2**53: only exact integer paths keep it
    as_list = bc.Tree.from_cells([deepest, Q453 * 9**15], 18).root_hex
    assert bc.Tree.from_cells(np.array([deepest, Q453 * 9**15], dtype=np.uint64), 18).root_hex == as_list
    assert bc.Tree.from_cells(np.array([deepest, Q453 * 9**15], dtype=np.uint64)[::-1], 18).root_hex == as_list
    assert bc.Tree.from_cells((c for c in [deepest, Q453 * 9**15]), 18).root_hex == as_list
    assert bc.Tree.from_cells({deepest, Q453 * 9**15, deepest}, 18).root_hex == as_list, "a set"
    with pytest.raises((TypeError, OverflowError)):
        bc.Tree.from_cells([Q453, -1], 3)
    with pytest.raises(TypeError):
        bc.Tree.from_cells([Q453, 1.5], 3)
    with pytest.raises(ValueError, match="deeper than the deepest level"):
        bc.Tree.from_cells([], 19)


def test_dggs_json_is_strict_json():
    t = bc.Tree.from_cells([bc.suid_to_cid("Q45")], 2)
    doc = t.to_dggs_json(bc.suid_to_cid("Q4"))
    assert bc.Tree.from_dggs_json(json.dumps(doc, allow_nan=False)).root_hex == t.root_hex, "a JSON string works too"
    for bad in (float("nan"), float("inf")):
        broken = json.loads(json.dumps(doc))
        broken["values"]["coverage"][0]["data"][0] = bad
        with pytest.raises(ValueError, match="null, not NaN"):
            bc.Tree.from_dggs_json(broken)
    boolean = json.loads(json.dumps(doc))
    boolean["values"]["coverage"][0]["data"] = [True] * 9
    with pytest.raises(ValueError, match="not a number"):
        bc.Tree.from_dggs_json(boolean)


@pytest.mark.parametrize("level", [0, 1, 2])
def test_cell_polygons_tile_the_globe_exactly_once(level):
    shapely = pytest.importorskip("shapely")
    coords, offsets = bc.cell_boundaries(bc.full_domain(level), n=7)
    polygons = shapely.polygons([shapely.linearrings(coords[a:b]) for a, b in zip(offsets[:-1], offsets[1:])])
    assert shapely.is_valid(polygons).all()
    total = shapely.area(polygons).sum()
    union = shapely.area(shapely.union_all(polygons))
    globe = 360.0 * 180.0
    assert abs(union - globe) < 1e-9 * globe, f"the cells cover {union} square degrees, not {globe}"
    assert abs(total - union) < 1e-9 * globe, f"the cells overlap by {total - union} square degrees"
