"""The array functions: agreement with the reference fixtures and the scalar functions, shapes,
dtypes, threads and errors."""
import json
from pathlib import Path

import numpy as np
import pytest

import burin as bc

FIX = Path(__file__).resolve().parents[2] / "crates" / "burin-core" / "tests" / "fixtures"
RANDOM = 160  # random points per level in points_*.json


def _f(h):
    return np.frombuffer(bytes.fromhex(h), dtype=">f8")[0]


@pytest.mark.parametrize("name,profile", [("ogc", bc.OGC_RHEALPIX), ("burin1", bc.BURIN_1)])
def test_points_match_the_reference(name, profile):
    rows = json.loads((FIX / f"points_{name}.json").read_text())["lonlat"]
    levels = sorted({r[2] for r in rows})
    per = len(rows) // len(levels)
    for k, level in enumerate(levels):
        block = rows[k * per:k * per + RANDOM]
        lon = np.array([_f(r[0]) for r in block]); lat = np.array([_f(r[1]) for r in block])
        got = bc.cells_from_lonlat(lon, lat, level, profile=profile)
        assert got.dtype == np.uint64
        assert [bc.cid_to_suid(int(c)) for c in got] == [r[3] for r in block], f"level {level}"


def test_centres_round_trip_and_agree_with_the_scalar_geometry():
    cells = bc.full_domain(3)
    lon, lat = bc.cells_to_lonlat(cells)
    assert np.array_equal(bc.cells_from_lonlat(lon, lat, 3), cells)
    geo = bc.cells_geojson([int(c) for c in cells[:50]])
    for c, f in zip(cells[:50], geo["features"]):
        assert f["properties"]["cid"] == int(c)


def test_shapes_broadcasting_and_contiguity():
    lon = np.linspace(-170, 170, 12).reshape(3, 4)
    got = bc.cells_from_lonlat(lon, 10.0, 5)
    assert got.shape == (3, 4)
    assert np.array_equal(bc.cells_from_lonlat(lon[:, ::2], 10.0, 5), got[:, ::2]), "strided input"
    clon, clat = bc.cells_to_lonlat(got)
    assert clon.shape == clat.shape == (3, 4)
    assert bc.cell_neighbours(got).shape == (3, 4, 4)
    assert bc.cells_from_lonlat([], [], 4).shape == (0,)


def test_threads_do_not_change_results():
    rng = np.random.default_rng(1)
    lon, lat = rng.uniform(-180, 180, 50_000), np.degrees(np.arcsin(rng.uniform(-1, 1, 50_000)))
    base = bc.cells_from_lonlat(lon, lat, 9, nthreads=1)
    for n in (0, 3):
        assert np.array_equal(bc.cells_from_lonlat(lon, lat, 9, nthreads=n), base)


def test_boundaries_are_ragged_closed_rings():
    shapely = pytest.importorskip("shapely")
    cells = bc.full_domain(2)
    coords, offsets = bc.cell_boundaries(cells, n=5)
    assert offsets.shape == (cells.size + 1,) and offsets[0] == 0 and offsets[-1] == len(coords)
    lon, lat = bc.cells_to_lonlat(cells)
    for i, c in enumerate(cells):
        ring = coords[offsets[i]:offsets[i + 1]]
        assert np.array_equal(ring[0], ring[-1])
        poly = shapely.Polygon(ring)
        assert poly.is_valid and poly.exterior.is_ccw, bc.cid_to_suid(int(c))
        x = lon[i] + 360 if ring[:, 0].max() > 180 and lon[i] < 0 else lon[i]
        is_cap = bc.cid_to_suid(int(c))[0] in "NS" and set(bc.cid_to_suid(int(c))[1:]) <= {"4"}
        if not is_cap:
            assert poly.contains(shapely.Point(x, lat[i])), bc.cid_to_suid(int(c))
    eq = bc.suid_to_cid("Q45")
    _, off = bc.cell_boundaries([eq])
    assert off[1] == 5, "an equatorial cell is its four corners"


def test_neighbours_match_the_scalar_function():
    cells = bc.full_domain(2)[::7]
    got = bc.cell_neighbours(cells)
    for c, row in zip(cells, got):
        nb = bc.neighbours(int(c))
        assert row.tolist() == [nb[k] for k in ("up", "right", "down", "left")]


def test_hierarchy_helpers():
    for r in range(4):
        cells = bc.full_domain(r)
        assert cells.size == 6 * 9**r and np.all(np.diff(cells.astype(np.int64)) == 1)
        assert np.all(bc.cell_levels(cells) == r)
    q453 = bc.suid_to_cid("Q453")
    assert bc.zoom_to([q453], 3, 1).tolist() == [bc.suid_to_cid("Q4")]
    kids = bc.zoom_to([bc.suid_to_cid("Q4")], 1, 3)
    assert kids.shape == (1, 81) and sorted(kids[0].tolist()) == sorted(bc.subzones(bc.suid_to_cid("Q4"), 2))
    assert bc.zoom_to(np.array([q453], dtype=np.int64), 3, 3).tolist() == [q453]


def test_errors():
    with pytest.raises(ValueError):
        bc.cells_from_lonlat([np.nan], [0.0], 3)
    with pytest.raises(ValueError):
        bc.cells_from_lonlat([0.0], [0.0], bc.MAX_LEVEL + 1)
    with pytest.raises(TypeError):
        bc.as_cell_ids([1.5])
    with pytest.raises(ValueError):
        bc.as_cell_ids([-1])
    with pytest.raises(ValueError):
        bc.as_cell_ids([bc.suid_to_cid("Q4")], level=2)
    with pytest.raises(ValueError):
        bc.cell_levels([5])
    with pytest.raises(ValueError):
        bc.cells_to_lonlat([5])


def test_tree_from_a_numpy_array():
    cells = bc.cells_from_lonlat(np.linspace(-0.2, 0.1, 200), np.linspace(51.4, 51.6, 200), 9)
    assert bc.Tree.from_cells(cells, 9).root_hex == bc.Tree.from_cells([int(c) for c in cells], 9).root_hex
