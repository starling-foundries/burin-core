"""Fingerprint, Proof and InvalidProof: one call answers a question or refuses, never both."""
import json

import numpy as np
import pytest

import burin
from burin import time as bt

PARK = {"type": "Polygon", "coordinates": [[[-73.9819, 40.7681], [-73.9730, 40.7644], [-73.9493, 40.7967],
                                            [-73.9580, 40.8003], [-73.9819, 40.7681]]]}
BETHESDA, TIMES_SQUARE = (-73.9712, 40.7740), (-73.9855, 40.7580)


@pytest.fixture(scope="module")
def park():
    return burin.Tree.from_geojson(PARK, 11)


@pytest.fixture(scope="module")
def hours():
    return bt.tree([("2026-07-18T10:00", "2026-07-19T05:00"), ("2026-07-19T10:00", "2026-07-20T05:00")])


def test_a_fingerprint_is_published_as_text_and_read_back_exactly(park, hours):
    for tree, axis in ((park, "space"), (hours, "time")):
        fp = tree.fingerprint
        assert fp.axis == axis and fp.root_hex == tree.root_hex and fp.depth == tree.depth
        assert burin.Fingerprint.parse(str(fp)) == fp and hash(burin.Fingerprint.parse(str(fp))) == hash(fp)
    text = str(park.fingerprint)
    assert text.startswith("burin:ogc-rhealpix:space:11:")
    for bad in (text.upper(), text.replace(":11:", ":011:"), text + "0", text.replace("space", "time")):
        with pytest.raises(ValueError):
            burin.Fingerprint.parse(bad)
    odd = burin.Profile(lon_0_udeg=10_000_000)
    fp = burin.Tree.from_geojson(PARK, 9, odd).fingerprint
    with pytest.raises(ValueError):
        burin.Fingerprint.parse(str(fp))
    assert burin.Fingerprint.parse(str(fp), [odd]) == fp


def test_questions_are_answered_in_the_world_s_terms(park, hours):
    place, open_ = park.fingerprint, hours.fingerprint
    assert place.check_point(park.prove_point(*BETHESDA), *BETHESDA) is True
    assert place.check_point(park.prove_point(*TIMES_SQUARE), *TIMES_SQUARE) is False
    for moment, expected in (("2026-07-18T18:32", True), ("2026-07-18T06:10", False), ("2026-07-19T04:59:59.999999", True),
                             ("2026-07-19T05:00", False), (np.datetime64("2026-07-19T12:00", "ns"), True)):
        assert open_.check_instant(hours.prove_instant(moment), moment) is expected, moment
    cell = int(burin.cells_from_lonlat(*BETHESDA, 11))
    proof = park.prove_cell(cell)
    assert (proof.claimed_cell, proof.claimed_covered) == (cell, True)
    assert place.check_cell(proof, cell) is True
    assert place.check_cell(proof.to_json(), cell) is True, "a proof may arrive as JSON"
    assert place.check_cell(json.dumps(proof.to_json()), cell) is True


def test_a_proof_that_answers_another_question_is_refused(park, hours):
    place, open_ = park.fingerprint, hours.fingerprint
    bethesda = park.prove_point(*BETHESDA)
    with pytest.raises(burin.InvalidProof, match="is about cell"):
        place.check_point(bethesda, *TIMES_SQUARE)
    just_bethesda = burin.Tree.from_cells([int(burin.cells_from_lonlat(*BETHESDA, 11))], 11)
    with pytest.raises(burin.InvalidProof, match="another fingerprint"):
        just_bethesda.fingerprint.check_point(bethesda, *BETHESDA)
    with pytest.raises(burin.InvalidProof, match="depth"):
        burin.Tree.from_geojson(PARK, 10).fingerprint.check_point(bethesda, *BETHESDA)
    with pytest.raises(burin.InvalidProof, match="H\\(9,6\\)"):
        open_.check_instant(bethesda, "2026-07-18T18:32")
    forged = bethesda.to_json()
    entries = forged["opening"]["entries"]
    i = next(i for i, e in enumerate(entries) if isinstance(e, str))
    entries[i] = entries[i][:-1] + ("1" if entries[i][-1] == "0" else "0")
    with pytest.raises(burin.InvalidProof, match="does not verify"):
        place.check_point(forged, *BETHESDA)
    with pytest.raises(burin.InvalidProof, match="not a proof"):
        place.check_point({"v": 1}, *BETHESDA)
    assert issubclass(burin.InvalidProof, ValueError)


def test_asking_the_wrong_kind_of_question_is_an_error_not_an_answer(park, hours):
    with pytest.raises(ValueError):
        park.fingerprint.check_instant(park.prove_point(*BETHESDA), "2026-07-18T18:32")
    with pytest.raises(ValueError):
        hours.fingerprint.check_point(hours.prove_instant("2026-07-18T18:32"), *BETHESDA)
    with pytest.raises(ValueError):
        park.prove_instant("2026-07-18T18:32")
    with pytest.raises(ValueError):
        hours.prove_point(*BETHESDA)
    coarse = burin.suid_to_cid(burin.cid_to_suid(int(burin.cells_from_lonlat(*BETHESDA, 11)))[:4])
    with pytest.raises(ValueError, match="partly covered"):
        park.prove_cell(coarse)


def test_relations_are_checked_against_both_fingerprints(park):
    reservoir = burin.Tree.from_geojson({"type": "Polygon", "coordinates": [[[-73.9660, 40.7820], [-73.9610, 40.7820],
                                        [-73.9610, 40.7870], [-73.9660, 40.7870], [-73.9660, 40.7820]]]}, 11)
    assert reservoir.holds(park, "within") and not reservoir.is_empty()
    midtown = burin.Tree.from_geojson({"type": "Polygon", "coordinates": [[[-73.99, 40.75], [-73.97, 40.75], [-73.97, 40.76],
                                      [-73.99, 40.76], [-73.99, 40.75]]]}, 11)
    for x, y in ((reservoir, park), (park, reservoir), (midtown, park), (park, park)):
        for rel in ("equals", "disjoint", "intersects", "within", "contains"):
            t = x.prove_relation(y, rel)
            assert x.fingerprint.check_relation(rel, y.fingerprint, t) is x.holds(y, rel), rel
    t = reservoir.prove_relation(park, "within")
    with pytest.raises(burin.InvalidProof):
        park.fingerprint.check_relation("within", reservoir.fingerprint, t)
    with pytest.raises(burin.InvalidProof):
        reservoir.fingerprint.check_relation("equals", park.fingerprint, t)
    with pytest.raises(burin.InvalidProof):
        reservoir.fingerprint.check_relation("within", park.fingerprint)
    assert reservoir.fingerprint.check_relation("within", park.fingerprint, burin.Transcript.from_json(t.to_json()))
