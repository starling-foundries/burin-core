# Changelog

What counts as a wire output, and what may change between releases, is set out in
[STABILITY.md](STABILITY.md).

## 0.3.0 (unreleased)

No wire output changed. Time is new.

### Added

- **The time line** (SPEC §9–§12), with the same program as space:
  - instants are POSIX microseconds, never floats, and fall in the profile's ticks (1 µs by
    default, centred on 1970, about 36,500 years either side);
  - time cells form the binary hierarchy `H(2, 1)` at a fixed depth of 61;
  - an interval's cells come from the HINT boundary walk, and a time set has one root with no
    resolution parameter;
  - Allen's thirteen relations between intervals.
- **In Rust:** `time::{tick, tick_start, interval_cells, coarse_cells, intervals_tree, allen, …}`.
- **In Python:** the `burin.time` module (`ticks`, `cells_from_times`, `times_of_cells`,
  `zoom_to`, `interval`, `tree`, `allen`, …) on datetime64 arrays of any unit, and
  `Tree.from_intervals`, `Tree.from_time_cells` and `Tree.axis`.
- **Five set relations** for space and time alike: `equals`, `disjoint`, `intersects`, `within`
  and `contains`. The evidence is a set-operation transcript, and it verifies against the two
  roots alone (`relation_transcript`, `verify_relation`, `Tree.holds`, `Tree.relation_evidence`).
- **Testing:** `time.json`, a conformance fixture computed from the specification's text by an
  independent reference (`tools/gen_time_fixtures.py`), plus a time fuzz target.
- `SPEC-time.md`, the draft of the space-time tree, zone data with time and the signed record.

### Changed

- `Tree.leaves()` refuses to list more than 2^28 leaves instead of trying to allocate them.
- `Tree.geojson()` refuses a time tree.

## 0.2.0 (2026-09-24)

No wire output changed. The frozen roots, polygon cell sets and topology fixtures of 0.1.0
reproduce exactly.

### Added

- `zone::cell_from_point` and `burin.cells_from_lonlat`: the cell containing a point. The rule
  is the reference's, bound for bound (SPEC §8).
- Array functions over points and cells in `burin`: `cells_from_lonlat`, `cells_to_lonlat`,
  `cell_boundaries`, `cell_neighbours`, `zoom_to`, `full_domain`, `level_range`, `cell_levels`
  and `as_cell_ids`. They take numpy arrays and run on all cores with the interpreter released.
- `Grid::cell_polygon` for drawing: equatorial cells as four corners, densified polar edges, an
  antimeridian shift and closed polar caps.
- `Tree::from_cells` builds in one pass. With the `parallel` feature, which the Python wheel
  enables, it builds on all cores. `Tree.from_cells` accepts a uint64 array.
- The answer is defined to the bit (SPEC §6): IEEE binary64 with no fused multiply-add, and
  every transcendental from the `libm` crate. `points_burin.json` freezes it, and
  `reference_disagreements.json` lists each 1–2 ulp difference from the reference.
- Python type stubs (`burin/_burin.pyi`, `py.typed`), checked against the compiled module.
- Fuzz targets for every parser and verifier (`fuzz/`), plus a stable test that applies every
  single edit to real records.
- Wheels for Linux, macOS and Windows on x86_64 and arm64 (abi3, Python 3.9+).
- STABILITY.md and this changelog.

### Changed

- The Python extension is now `burin._burin`, inside a `burin` package. `import burin` works as
  before.

### Refused

Input that was malformed, or had two spellings, and used to be accepted:

- Record fields `A`, `B` and `D` beyond a u32. They were read as their low 32 bits, so
  `"A": 4294967305` was read as 9.
- Records whose `A` and `B` are not a hierarchy with `2 ≤ A ≤ 4096`, or whose `D` is deeper
  than the deepest level (18 for rHEALPix). An unbounded `A` let a record make the verifier
  allocate and hash without limit. A 0.1.0 tree deeper than 18 levels could not hold a cell id
  at its leaves.
- Digests in upper-case hex. The crate writes lower case, and each digest now has one spelling.
- DGGS-JSON depths beyond a u32 (they were truncated), and depths that reach below the deepest
  level under their zone.
- Cell ids below the deepest level. Some level-19 ids fit in 64 bits, but the level is not
  whole. `Hierarchy::cid` already refused the matching paths.
- Latitudes outside [−90, 90], which used to fold into the other hemisphere, and non-finite
  coordinates.
- Levels above 18 in cell lookups, where the arithmetic used to wrap silently.

### Fixed

- `Hierarchy::max_level` no longer loops forever on a hand-built hierarchy with `A < 2`.

## 0.1.0

The first public release. It has:

- the coverage commitment (profile, cid and suid, polygon rule, canonical tree, SHA-256 root);
- openings and set-operation proofs;
- the zone topology (scanline sub-zone order, edge neighbours, rasters and halos);
- DGGS-JSON zone data;
- the conformance fixtures against rhealpixdggs-py 0.8.6 and DGGAL.
