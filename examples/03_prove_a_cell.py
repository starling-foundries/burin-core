"""Prove a cell: membership and non-membership openings, and tamper evidence.

An opening is a small self-describing record: a verifier replays it against
the root with no geometry, no tree, no trust. Flip one hex digit and the
proof dies.
"""

import json

import burin

RESOLUTION = 10

greenwich = {"type": "Polygon", "coordinates": [[
    [-0.010, 51.500], [0.000, 51.500], [0.000, 51.508], [-0.010, 51.508], [-0.010, 51.500],
]]}
elsewhere = {"type": "Polygon", "coordinates": [[
    [0.010, 51.500], [0.020, 51.500], [0.020, 51.508], [0.010, 51.508], [0.010, 51.500],
]]}

tree = burin.Tree.from_geojson(greenwich, RESOLUTION)
print(f"root: {tree.root_hex[:24]}…  ({tree.leaf_count()} cells at resolution {RESOLUTION})")


def claim_of(record):
    """The one-word verdict buried at the bottom of the opening ladder."""
    if isinstance(record, dict):
        if "claim" in record:
            return record["claim"]
        for v in record.values():
            found = claim_of(v)
            if found:
                return found
    elif isinstance(record, list):
        for v in record:
            found = claim_of(v)
            if found:
                return found
    return None


def show(label, cid):
    record = tree.open(cid)
    suid = burin.cid_to_suid(cid)
    ok = burin.verify_opening(record)
    print(f"\n{label}: cell {suid}")
    print(f"  claim: {claim_of(record)!r}   record: {len(json.dumps(record)):,} bytes of JSON")
    print(f"  verify_opening(record) -> {ok}")
    return record


# A cell deep inside the coverage: provably full.
inside = tree.leaves()[0]
good = show("membership", inside)

# A cell from a disjoint box across town: provably empty.
outside = burin.Tree.from_geojson(elsewhere, RESOLUTION).leaves()[0]
show("non-membership", outside)

# A coarser cell that straddles the boundary has no single claim to make.
straddler = tree.cells()[0] // 9  # parent of a coarsest covered cell is always partial
print(f"\nparent cell {burin.cid_to_suid(straddler)} straddles the edge: "
      f"open() -> {tree.open(straddler)}")

# The 10-line tamper-evidence demo: flip one hex digit deep in the record.
tampered = json.loads(json.dumps(good))  # deep copy
entries = tampered["opening"]["entries"]
k = next(i for i, v in enumerate(entries) if isinstance(v, str))  # a sibling hash
entries[k] = ("0" if entries[k][0] != "0" else "1") + entries[k][1:]
print(f"\nsame record, one hex digit flipped: verify_opening -> "
      f"{burin.verify_opening(tampered)}")
print("the root commits to every bit of the coverage; the opening carries that commitment down.")
