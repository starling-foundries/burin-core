"""burin.time and the set relations, through the public Python API."""
import json
from pathlib import Path

import numpy as np
import pytest

import burin
from burin import time as T

FIXTURE = json.loads((Path(__file__).resolve().parents[2] / "crates/burin-core/tests/fixtures/time.json").read_text())


def test_every_datetime64_unit_floors_to_the_microsecond():
    for unit, value, us in [("ns", -1, -1), ("ns", -1000, -1), ("ns", 999, 0), ("ps", -1, -1), ("s", -3, -3_000_000),
                            ("D", -1, -86_400_000_000), ("W", 1, 604_800_000_000), ("M", 1, 31 * 86_400_000_000),
                            ("Y", 1, 365 * 86_400_000_000), ("us", 7, 7)]:
        assert int(T.posix_us(np.datetime64(value, unit))) == us, (unit, value)
    ten_seconds = np.array([1, -1], dtype="datetime64[10s]")
    assert T.posix_us(ten_seconds).tolist() == [10_000_000, -10_000_000]
    assert T.posix_us(np.array(["2013-07-15"]))[0] == 1_373_846_400_000_000


def test_what_is_not_an_instant_is_refused():
    with pytest.raises(TypeError):
        T.ticks(1.5)
    with pytest.raises(TypeError):
        T.ticks(np.array([True]))
    with pytest.raises(ValueError):
        T.ticks(np.array(["NaT"], dtype="datetime64[s]"))
    with pytest.raises(ValueError):
        T.ticks(np.uint64(2**63))
    with pytest.raises(ValueError):
        T.ticks(np.datetime64(2**62, "D"))  # beyond the 64-bit microsecond range
    with pytest.raises(ValueError):
        T.ticks(np.int64(-(2**60) - 1))  # off the default line


def test_the_public_api_reproduces_the_reference_roots():
    for case in FIXTURE["intervals"]:
        ivs = [tuple(i) for i in case["intervals"]]
        tree = burin.Tree.from_intervals(ivs)
        assert tree.root_hex == case["root"]
        assert sorted(tree.cells()) == sorted(case["cells"])
        assert burin.Tree.from_time_cells(np.array(case["cells"], dtype=np.uint64)).root_hex == case["root"]
    day = T.tree([("1970-01-01", "1970-01-02")])
    assert day.root_hex == FIXTURE["intervals"][5]["root"]
    assert burin.Tree.from_intervals([]).root_hex == FIXTURE["empty_root"]


def test_cells_contain_their_instants_at_every_level():
    times = np.array(["1969-12-31T23:59:59.999999", "2013-07-15T12:34:56.789", "1900-01-01", "2100-06-30T23:00"],
                     dtype="datetime64[us]")
    for level in (0, 1, 20, 40, 61):
        cells = T.cells_from_times(times, level)
        assert (T.cell_levels(cells) == level).all()
        start, end = T.times_of_cells(cells)
        assert ((start <= times) & (times < end)).all()
        assert (T.cells_from_times(start, level) == cells).all()
        if level:
            assert (T.zoom_to(cells, level, level - 1) == T.cells_from_times(times, level - 1)).all()
    kids = T.zoom_to(T.cells_from_times(times, 40), 40, 43)
    assert kids.shape == (4, 8) and (T.cell_levels(kids) == 43).all()


def test_an_interval_is_every_tick_it_touches():
    day = burin.Profile(tick_us=86_400_000_000)
    lo, hi = T.interval("2013-07-15T12:00", "2013-07-17T00:00", profile=day)
    assert hi - lo == 2, "the half day and the whole next day"
    lo, hi = T.interval("2013-07-15T12:00", "2013-07-17T00:00:00.000001", profile=day)
    assert hi - lo == 3, "a microsecond into a third day"
    lo, hi = T.interval(np.datetime64("2013-07-15"), np.datetime64("2013-07-16"))
    assert hi - lo == 86_400_000_000
    with pytest.raises(ValueError):
        T.interval("2013-07-16", "2013-07-15")
    cells = T.interval_cells(lo, hi)
    assert len(cells) <= 122 and burin.Tree.from_time_cells(cells).root_hex == burin.Tree.from_intervals([(lo, hi)]).root_hex


def test_allen_names_every_relation():
    names = {T.allen((a, b), (2, 5)) for a in range(0, 7) for b in range(a + 1, 8)}
    assert len(names) == 13
    assert T.allen((0, 2), (2, 4)) == "meets" and T.allen((2, 4), (0, 2)) == "met-by"
    with pytest.raises(ValueError):
        T.allen((3, 3), (0, 1))


def _space(*suids):
    return burin.Tree.from_cells([burin.suid_to_cid(s) for s in suids], 3)


@pytest.mark.parametrize("x,y", [
    (_space("Q4"), _space("Q41", "Q418")),
    (_space("Q41"), _space("Q4", "N")),
    (_space("Q4"), _space("S")),
    (T.tree([("2013-07-01", "2013-08-01")]), T.tree([("2013-07-15", "2013-07-16")])),
    (T.tree([("2013-07-01", "2013-08-01")]), T.tree([("2014-01-01", "2014-02-01")])),
])
def test_relation_evidence_verifies_exactly_when_the_relation_holds(x, y):
    for rel in ("equals", "disjoint", "intersects", "within", "contains"):
        evidence = x.relation_evidence(y, rel)
        assert burin.verify_relation(rel, x.root_hex, y.root_hex, evidence) == x.holds(y, rel), rel
        if evidence is not None:
            forged = json.loads(json.dumps(evidence))
            forged["root_c"] = forged["root_a"]
            assert not burin.verify_relation(rel, x.root_hex, y.root_hex, forged) or x.holds(y, rel)
            assert not burin.verify_relation(rel, "00" * 32, y.root_hex, evidence)
            assert not burin.verify_relation("equals", x.root_hex, y.root_hex, evidence)
    with pytest.raises(ValueError):
        x.holds(y, "overlaps")


def test_time_and_space_trees_stay_apart():
    t, s = T.tree([("2013-07-15", "2013-07-16")]), _space("Q4")
    assert (t.axis, s.axis) == ("time", "space")
    for call in (lambda: t.union(s), lambda: t.holds(s, "within"), t.geojson, t.raster,
                 lambda: t.to_dggs_json(burin.suid_to_cid("Q4")), burin.Tree.from_intervals([(0, T.TICKS)]).leaves):
        with pytest.raises(ValueError):
            call()
    record = t.open(int(T.cells_from_times(np.datetime64("2013-07-15T08:00"))))
    assert burin.verify_opening(record)
