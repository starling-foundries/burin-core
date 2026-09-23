# burin-core

**A deterministic, equal-area fingerprint of a geographic region.** Give it a polygon and a
resolution; it returns a 32-byte root that anyone, in any language, recomputes from the same
polygon and the same profile, down to the last bit. Two encodings of the same footprint (a COG,
a Zarr, a different tiling) produce one root. Two roots can be unioned, intersected, differenced
and compared without exchanging the geometry, and a cell can be proven inside or outside a root
with a short opening.

The ground is cut into the **rHEALPix** discrete global grid (the OGC API-DGGS Annex B
registered grid, `lon_0 = 50`, WGS84), which is exactly equal-area *and* exactly hierarchical:
nine children tile their parent, so a region held as one cell and the same region spelled out to
the leaves hash identically (the roll-up). Hashing is SHA-256 with domain-tagged nodes. The
hash, the grid parameters, and the time convention are bound into a **profile** whose identity
is derived from its integers, never transmitted.

Apache-2.0. Rust core with a Python wheel (not yet on PyPI: build it with `maturin build --release`, or point uv at this checkout with `burin-core = { path = "../burin-core" }` under `[tool.uv.sources]`); a JavaScript build is
planned. The vocabulary of OGC API - DGGS (zone, refinement level, compact zone list, DGGRS) maps onto the spec one to one; see SPEC.md §7. The projection is a port of the reference implementation
[rhealpixdggs-py](https://github.com/manaakiwhenua/rhealpixdggs-py) (MIT), kept in bit-level
parity by the fixtures under `crates/burin-core/tests/fixtures`.

## Python

```python
import burin as bc

aoi = {"type": "Polygon", "coordinates": [[[-0.13, 51.50], [-0.10, 51.50], [-0.10, 51.52], [-0.13, 51.52], [-0.13, 51.50]]]}

bc.fingerprint_polygon(aoi, resolution=10)          # '…' 64 hex chars; same polygon → same bytes, anywhere
bc.OGC_RHEALPIX.id_hex()                            # the profile identity the root is read against

t = bc.Tree.from_geojson(aoi, 10)                   # the canonical coverage
t.cells()                                           # coarsest covered cell ids; t.leaves() expands them
u = t.union(bc.Tree.from_geojson(other, 10))        # set algebra over coverages
proof = t.prove(u, "union"); bc.verify_setop(proof) # a transcript a verifier replays with no geometry
op = t.open(bc.suid_to_cid("Q453")); bc.verify_opening(op)   # membership / non-membership
bc.subzones(bc.suid_to_cid("Q4"), 2)              # sub-zones in the OGC scanline order (SPEC §8)
bc.neighbours(bc.suid_to_cid("Q453"))             # {'up': …, 'right': …, 'down': …, 'left': …}
bc.Tree.from_geojson(aoi, 6).raster()              # uint8 bytes, shape (6, 3**6, 3**6); halo_index(6) pads it
t.to_dggs_json(bc.suid_to_cid("Q"))                 # OGC DGGS-JSON zone data; bc.Tree.from_dggs_json reads it back
```

## Rust

```rust
use burin_core::{fingerprint_geojson, Profile};
let root = fingerprint_geojson(&geojson, 10, &Profile::ogc())?;
```

## The rule, in one paragraph

A polygon's cells at resolution R are the depth-R rHEALPix cells whose **nucleus** (the inverse
projection of the planar cell centre) the polygon contains, boundary excluded. The set is held
in its coarsest form and hashed bottom-up: `leaf_empty = H(0x00)`, `leaf_full = H(0x01)`,
`node = H(0x02 ‖ c₀‖…‖c₈)`, `root = H(0x03 ‖ b₀‖…‖b₅)`; a constant subtree hashes to its ladder
value, so tiling does not matter. The profile identity is
`SHA-256("BURIN-PROFILE-1" ‖ hash id ‖ A ‖ B ‖ lon_0 in µdeg ‖ polar squares ‖ a in µm ‖ 1/f in nano ‖ tick ‖ epoch)`.
`N^-r` is the reciprocal of the exact integer `N^r`. See `SPEC.md`.

## Development

```bash
cargo test -p burin-core                      # parity fixtures + structural properties
uv run --no-project --with rhealpixdggs==0.8.6 --with shapely python tools/gen_fixtures.py   # regenerate set 1
cargo run -p burin-core --example freeze     # regenerate set 2 (golden roots; a wire-format bump)
uv venv && uv pip install maturin pytest && uv run maturin develop && uv run pytest
```
