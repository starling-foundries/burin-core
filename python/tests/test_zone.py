"""Zone topology from Python: scanline order, neighbours, rasters and halos, against the fixtures."""
import json
from array import array
from pathlib import Path

import pytest

import burin as bc

FIX = Path(__file__).resolve().parents[2] / "crates" / "burin-core" / "tests" / "fixtures"


def test_subzones_follow_dggal_order():
    fx = json.loads((FIX / "dggal_rhealpix.json").read_text())
    for key, want in fx["subzones"].items():
        parent, depth = key.split("/")
        got = [bc.cid_to_suid(c) for c in bc.subzones(bc.suid_to_cid(parent), int(depth))]
        assert got == want, key
    q4 = bc.suid_to_cid("Q4")
    for i, c in enumerate(bc.subzones(q4, 2)):
        assert bc.subzone_index(q4, c) == i and bc.subzone_at(q4, 2, i) == c
    with pytest.raises(ValueError):
        bc.subzone_index(q4, bc.suid_to_cid("R0"))


def test_position_and_cell_at():
    c = bc.suid_to_cid("S0812")
    base, level, row, col = bc.position(c)
    assert (base, level) == (5, 4)
    assert bc.cell_at(base, level, row, col) == c
    with pytest.raises(ValueError):
        bc.cell_at(0, 1, 3, 0)


def test_neighbours_match_the_reference_and_depend_on_the_profile():
    fx = json.loads((FIX / "topology.json").read_text())
    rows = {r[0]: r[1:] for r in fx["neighbours"]["12"]}
    p = bc.Profile(north_square=1, south_square=2)
    for suid in ("N0", "N48", "S", "S26", "O30", "R85"):
        nb = bc.neighbours(bc.suid_to_cid(suid), p)
        assert [bc.cid_to_suid(nb[k]) for k in ("up", "right", "down", "left")] == rows[suid], suid
    assert bc.neighbours(bc.suid_to_cid("N0")) != bc.neighbours(bc.suid_to_cid("N0"), p)


def test_halo_index_shape_and_corners():
    level, width = 2, 1
    n, m = 3**level, 3**level + 2 * width
    idx = array("q")
    idx.frombytes(bc.halo_index(level, width))
    assert len(idx) == 6 * m * m
    assert sum(1 for x in idx if x < 0) == 24
    at = lambda b, i, j: idx[b * m * m + i * m + j]
    cell = lambda k: bc.cell_at(k // (n * n), level, (k % (n * n)) // n, k % n)
    for b in range(6):
        for r in range(n):
            for c in range(n):
                nb = bc.neighbours(cell(at(b, r + 1, c + 1)))
                assert cell(at(b, r, c + 1)) == nb["up"] and cell(at(b, r + 2, c + 1)) == nb["down"]
                assert cell(at(b, r + 1, c)) == nb["left"] and cell(at(b, r + 1, c + 2)) == nb["right"]
    with pytest.raises(ValueError):
        bc.halo_index(1, 4)


def test_raster_index_and_cells_move_values_both_ways():
    level = 2
    cells = array("Q")
    cells.frombytes(bc.raster_cells(level))
    assert len(cells) == 6 * 81 and cells[0] == bc.cell_at(0, level, 0, 0)
    some = [bc.suid_to_cid(s) for s in ("S88", "N04", "Q40")]
    idx = array("q")
    idx.frombytes(bc.raster_index(some))
    assert [cells[i] for i in idx] == some
    with pytest.raises(ValueError):
        bc.raster_index([bc.suid_to_cid("Q4"), bc.suid_to_cid("Q45")])


def test_raster_round_trip():
    aoi = {"type": "Polygon", "coordinates": [[[-0.13, 51.50], [-0.10, 51.50], [-0.10, 51.52], [-0.13, 51.52], [-0.13, 51.50]]]}
    t = bc.Tree.from_geojson(aoi, 6)
    r = t.raster()
    assert len(r) == 6 * 729 * 729 and sum(r) == t.leaf_count()
    back = bc.Tree.from_raster(r, 6)
    assert back.root_hex == t.root_hex and back.cells() == t.cells()
    with pytest.raises(ValueError):
        bc.Tree.from_raster(r[:-1], 6)


def test_dggs_json_round_trip_and_refusals():
    aoi = {"type": "Polygon", "coordinates": [[[-0.13, 51.50], [-0.10, 51.50], [-0.10, 51.52], [-0.13, 51.52], [-0.13, 51.50]]]}
    t = bc.Tree.from_geojson(aoi, 8)
    zone = t.cells()[0] // 9**2                        # an ancestor two levels above a covered cell
    doc = t.to_dggs_json(zone)
    assert json.loads(json.dumps(doc)) == doc          # plain JSON
    assert doc["dggrs"] == "https://www.opengis.net/def/dggrs/OGC/1.0/rHEALPix"
    assert doc["zoneId"] == bc.cid_to_suid(zone) and doc["depths"] == [8 - bc.position(zone)[1]]
    back = bc.Tree.from_dggs_json(doc)
    assert back.root_hex == t.intersect(bc.Tree.from_cells([zone], 8)).root_hex
    data = doc["values"]["coverage"][0]["data"]
    doc["values"]["coverage"][0]["data"] = [0 if v is None else v for v in data]
    assert bc.Tree.from_dggs_json(doc).leaf_count() == len(data)
    assert bc.Tree.from_dggs_json(doc, nonzero=True).root_hex == back.root_hex
    doc["dggrs"] = "https://www.opengis.net/def/dggrs/OGC/1.0/ISEA3H"
    with pytest.raises(ValueError):
        bc.Tree.from_dggs_json(doc)
    with pytest.raises(ValueError):
        bc.Tree.from_geojson(aoi, 8, bc.BURIN_1).to_dggs_json(zone)
