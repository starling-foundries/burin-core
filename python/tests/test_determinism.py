"""The array functions against burin's frozen answers, bit for bit (points_burin.json): every
level of every adversarial point, and the nuclei, under each profile. No tolerance; the release
workflow runs this on every wheel's own platform."""
import json
from pathlib import Path

import numpy as np
import pytest

import burin as bc

FIX = Path(__file__).resolve().parents[2] / "crates" / "burin-core" / "tests" / "fixtures"
GOLDENS = json.loads((FIX / "points_burin.json").read_text())
DEEPEST = GOLDENS["deepest"]


def _floats(hexes):
    return np.frombuffer(bytes.fromhex("".join(hexes)), dtype=">f8").astype("f8")


@pytest.mark.parametrize("name", sorted(GOLDENS["profiles"]))
def test_every_level_of_every_point_is_frozen(name):
    g = GOLDENS["profiles"][name]
    profile = bc.Profile(lon_0_udeg=g["lon_0_udeg"], north_square=g["north_square"], south_square=g["south_square"])
    rows = g["points"]
    lon, lat = _floats([r[0] for r in rows]), _floats([r[1] for r in rows])
    deepest = np.array([r[2] for r in rows], dtype=np.uint64)
    on_grid = deepest != 0
    for level in range(DEEPEST + 1):
        got = bc.cells_from_lonlat(lon[on_grid], lat[on_grid], level, profile=profile)
        want = deepest[on_grid] // np.uint64(9 ** (DEEPEST - level))
        bad = np.nonzero(got != want)[0]
        assert bad.size == 0, f"{name}: level {level}, first difference at point {np.nonzero(on_grid)[0][bad[0]]}"
    for i in np.nonzero(~on_grid)[0][:50]:
        with pytest.raises(ValueError):
            bc.cells_from_lonlat(lon[i], lat[i], 3, profile=profile)


@pytest.mark.parametrize("name", sorted(GOLDENS["profiles"]))
def test_nuclei_are_frozen(name):
    g = GOLDENS["profiles"][name]
    profile = bc.Profile(lon_0_udeg=g["lon_0_udeg"], north_square=g["north_square"], south_square=g["south_square"])
    rows = g["nuclei"]
    lon, lat = bc.cells_to_lonlat(np.array([r[0] for r in rows], dtype=np.uint64), profile=profile)
    assert np.array_equal(lon.view("u8"), _floats([r[1] for r in rows]).view("u8"))
    assert np.array_equal(lat.view("u8"), _floats([r[2] for r in rows]).view("u8"))
