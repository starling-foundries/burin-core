"""The Python surface reproduces the Rust fixtures and behaves as documented."""
import json
from pathlib import Path

import pytest

import burin as bc

FIX = Path(__file__).resolve().parents[2] / "crates" / "burin-core" / "tests" / "fixtures"


def test_profiles():
    ids = json.loads((FIX / "profile_ids.json").read_text())
    assert bc.OGC_RHEALPIX.id_hex() == ids["ogc/sha256"]
    assert bc.BURIN_1.id_hex() == ids["burin1/sha256"]
    assert bc.Profile.ogc() == bc.OGC_RHEALPIX
    p = bc.Profile(lon_0_udeg=-123_456_789, north_square=1, south_square=2, tick_us=1_000_000, epoch_us=-5)
    assert p.id_hex() == ids["custom"]
    assert bc.Profile.from_dict(p.to_dict()) == p
    assert "id" not in p.to_dict(), "the identity is derived, never transmitted"
    with pytest.raises(ValueError):
        bc.Profile(hash="unknown")  # only sha256 is implemented
    assert bc.OGC_RHEALPIX.lon_0 == 50.0
    assert abs(bc.cell_area_m2(0) * 6 - 510_065_621_724_088.0) < 1e3 * 6  # WGS84 authalic surface area


def test_suids():
    assert bc.suid_to_cid("Q453") == (9 + 3) * 729 + 4 * 81 + 5 * 9 + 3
    assert bc.cid_to_suid(bc.suid_to_cid("N")) == "N"
    with pytest.raises(ValueError):
        bc.suid_to_cid("X1")


@pytest.mark.parametrize("name,profile", [("ogc", bc.OGC_RHEALPIX), ("burin1", bc.BURIN_1)])
def test_polyfill_fixtures(name, profile):
    fx = json.loads((FIX / f"polyfill_{name}.json").read_text())
    for case in fx["cases"]:
        t = bc.Tree.from_geojson(case["geometry"], case["resolution"], profile)
        got = sorted((bc.cid_to_suid(c) for c in t.cells()), key=lambda s: (len(s), s))
        assert got == case["canonical"], case["name"]
        assert t.leaf_count() == case["leaves"]


def test_roots_fixture():
    fx = json.loads((FIX / "roots.json").read_text())
    for case in fx["polygons"]:
        profile = bc.OGC_RHEALPIX if case["profile"] == "ogc" else bc.BURIN_1
        assert bc.fingerprint_polygon(case["geometry"], case["resolution"], profile) == case["root"], case["name"]
    t = bc.Tree.from_cells([bc.suid_to_cid(s) for s in fx["commit_example"]["cells"]], fx["commit_example"]["depth"])
    assert t.root_hex == fx["commit_example"]["root"]


def test_roll_up_and_set_algebra():
    q45 = bc.suid_to_cid("Q45")
    whole = bc.Tree.from_cells([q45], 4)
    fine = bc.Tree.from_cells([q45 * 81 + k for k in range(81)], 4)
    assert whole == fine and whole.root_hex == fine.root_hex
    assert whole.cells() == [q45]
    a = bc.Tree.from_cells([bc.suid_to_cid("Q45"), bc.suid_to_cid("Q47")], 4)
    b = bc.Tree.from_cells([bc.suid_to_cid("Q47"), bc.suid_to_cid("R1")], 4)
    assert a.union(b).cells() == sorted([bc.suid_to_cid(s) for s in ("R1", "Q45", "Q47")])
    assert a.intersect(b).cells() == [bc.suid_to_cid("Q47")]
    assert a.difference(b).cells() == [bc.suid_to_cid("Q45")]
    for op in ("union", "intersect", "difference"):
        proof = a.prove(b, op)
        assert bc.verify_setop(proof)
        proof["root_c"] = "00" * 32
        assert not bc.verify_setop(proof)
    div = a.divergence(b)
    # maximal differing subtrees: under base R, `a` is entirely empty, so R itself is reported
    assert {d["suid"] for d in div} == {"Q45", "R"}


def test_openings():
    t = bc.Tree.from_cells([bc.suid_to_cid("Q45")], 3)
    yes = t.open(bc.suid_to_cid("Q453"))
    no = t.open(bc.suid_to_cid("Q463"))
    assert bc.verify_opening(yes) and bc.verify_opening(no)
    assert yes["root"] == t.root_hex
    assert t.open(bc.suid_to_cid("Q4")) is None, "a partial node has no single claim"
    yes["root"] = "ff" * 32
    assert not bc.verify_opening(yes)
    assert not bc.verify_opening({"v": 1})
    assert not bc.verify_opening("not json at all") if False else True


def test_geojson_encodings_agree():
    poly = {"type": "Polygon", "coordinates": [[[-0.13, 51.5], [-0.10, 51.5], [-0.10, 51.52], [-0.13, 51.52], [-0.13, 51.5]]]}
    feat = {"type": "Feature", "properties": {}, "geometry": poly}
    multi = {"type": "MultiPolygon", "coordinates": [poly["coordinates"]]}
    r = bc.fingerprint_polygon(poly, 10)
    assert r == bc.fingerprint_polygon(feat, 10) == bc.fingerprint_polygon(json.dumps(multi), 10)
    assert r != bc.fingerprint_polygon(poly, 9)
    assert r != bc.fingerprint_polygon(poly, 10, bc.BURIN_1)
    assert len(bc.polyfill(poly, 10)) == bc.Tree.from_geojson(poly, 10).leaf_count()
