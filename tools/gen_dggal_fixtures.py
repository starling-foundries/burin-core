#!/usr/bin/env python3
"""Generate the OGC reference fixture from DGGAL, the OGC API - DGGS reference library.

DGGAL's registered rHEALPix DGGRS (lon_0 = 50, both polar squares at column 0) defines the
sub-zone order that DGGS-JSON zone data uses. This records that order for sample parents and the
edge neighbours of every zone at levels 0-2; `tests/zone.rs` checks the crate against it.

Run with a Python that can load the `dggal` wheel (on Apple silicon the macOS wheel is x86_64
only, so use an x86_64 interpreter under Rosetta):
    uv venv dggal-venv --python cpython-3.12-macos-x86_64-none && uv pip install --python dggal-venv/bin/python dggal
    arch -x86_64 dggal-venv/bin/python tools/gen_dggal_fixtures.py
Unlike tools/gen_fixtures.py this is not regenerated in CI; the output is committed.
"""
from __future__ import annotations

import json
from pathlib import Path

from dggal import Application, pydggal_setup, rHEALPix

OUT = Path(__file__).resolve().parent.parent / "crates" / "burin-core" / "tests" / "fixtures" / "dggal_rhealpix.json"
PARENTS = [("N", 1), ("N", 2), ("O", 2), ("P", 1), ("Q", 3), ("R", 2), ("S", 1), ("S", 2),
           ("Q4", 2), ("N81", 2), ("S0", 2), ("N0", 1), ("S8", 1), ("O2", 1), ("R26", 1)]


def zones(level: int) -> list[str]:
    out = [b for b in "NOPQRS"]
    for _ in range(level):
        out = [z + str(k) for z in out for k in range(9)]
    return out


def main() -> None:
    app = Application(appGlobals=globals())
    pydggal_setup(app)
    d = rHEALPix()
    tid = d.getZoneTextID
    subzones = {f"{p}/{k}": [tid(z) for z in d.getSubZones(d.getZoneFromTextID(p), k)] for p, k in PARENTS}
    neighbours = {z: sorted(tid(n) for n in d.getZoneNeighbors(d.getZoneFromTextID(z)))
                  for level in range(3) for z in zones(level)}
    OUT.write_text(json.dumps({"subzones": subzones, "neighbours": neighbours}, indent=1, sort_keys=True) + "\n")
    print(f"wrote {OUT.name}: {len(subzones)} parents, {len(neighbours)} zones")


if __name__ == "__main__":
    main()
