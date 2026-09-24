# Coverage commitment: specification (draft 1)

*Everything an implementer needs to produce the same 32 bytes from the same polygon. The Rust
crate is the reference implementation; the fixtures under `crates/burin-core/tests/fixtures`
are the conformance suite.*

## 1. Profile

A profile is the tuple of integers a cell identifier and a root are read against:

| field | type | default (`ogc-rhealpix`) | meaning |
|---|---|---|---|
| `hash` | ASCII id | `sha256` | the node hash (§4); any other id is refused |
| `aperture` A | u32 | 9 | children per cell (`N_side²`) |
| `n_base` B | u32 | 6 | base cells |
| `lon_0_udeg` | i64 | 50 000 000 | central meridian in micro-degrees |
| `north_square`, `south_square` | u8 | 0, 0 | polar square placement |
| `a_um` | u64 | 6 378 137 000 000 | semi-major axis in micro-metres (WGS84) |
| `inv_f_nano` | u64 | 298 257 223 563 | inverse flattening × 10⁹ (WGS84) |
| `tick_us`, `epoch_us` | u64, i64 | 1, 0 | time convention (reserved) |

Identity: `id = SHA-256("BURIN-PROFILE-1" ‖ u8 len(hash) ‖ hash ‖ u32 A ‖ u32 B ‖ i64 lon_0_udeg ‖ u8 ns ‖ u8 ss ‖ u64 a_um ‖ u64 inv_f_nano ‖ u64 tick_us ‖ i64 epoch_us)`, all big-endian. A wire record carries the fields; the receiver re-derives the id. `burin-1` is the same profile with `lon_0_udeg = 0`.

Derived floats: `lon_0 = lon_0_udeg / 10⁶`, `a = a_um / 10⁶`, `inv_f = inv_f_nano / 10⁹`, `f = 1/inv_f`, `e = sqrt(f(2−f))`.

## 2. Cells

The hierarchy H(A, B): a **path** `(b, d₁ … d_r)` with `b ∈ [0,B)`, `d_i ∈ [0,A)`; its **cid** is the integer `(A + b)·A^r + Σ d_i·A^(r−1−i)`. Parent = `cid div A`, child k = `A·cid + k`, descendants k levels down = `[cid·A^k, (cid+1)·A^k)`. The rHEALPix **suid** `Q453` is base letter `NOPQRS[b]` followed by the digits.

Cell geometry follows rhealpixdggs-py 0.8.6 exactly (§6): base upper-left vertices on the unit
authalic sphere `N (−π + ns·π/2, 3π/4)`, `O (−π, π/4)`, `P (−π/2, π/4)`, `Q (0, π/4)`, `R (π/2, π/4)`, `S (−π + ss·π/2, −π/4)`, scaled by `R_A`; child k of a cell occupies row `k div N`, column `k mod N` of its parent's N×N split; `width(r) = R_A·(π/2)·N^(−r)`; the **nucleus** is the inverse projection of the planar cell centre `(ul.x + w/2, ul.y − w/2)`.

`N^(−r)` **is the correctly rounded reciprocal of the exact integer `N^r`** (`1.0 / (N^r as f64)`), not a call to `pow`.

## 3. The polygon rule

`cells(P, R)` = the depth-R cells whose nucleus lies in the **interior** of P (boundary excluded, exact orientation arithmetic). P is a GeoJSON Polygon or MultiPolygon in WGS84 lon/lat with longitudes in `[−180, 180]`. A pruning descent from the six base cells is permitted only with a bounding box that contains the whole cell and a shortcut that emits a subtree only when its box lies in the **interior** of P (DE-9IM `T**FF*FF*`), so that pruning never changes the set (the reference does this; see `polyfill.rs`).

**Known tie surface.** A polygon split at the antimeridian has seam edges at `x = ±180`. A nucleus whose longitude is exactly `−180.0` lies on that seam and is excluded by the rule. Such nuclei exist (the fixtures list them under `seam_nuclei`), and an implementation must reproduce their longitude to the bit to agree with the reference; the conformance test checks exactly this.

## 4. The commitment

The set is held in canonical form: a node is `EMPTY`, `FULL`, or a branch of A children not all the same constant. With `H` the profile's hash:

```
leaf_empty = H(0x00)            EMPTY[0] = leaf_empty   EMPTY[d] = node(EMPTY[d−1] × A)
leaf_full  = H(0x01)            FULL[0]  = leaf_full    FULL[d]  = node(FULL[d−1]  × A)
node(c₀…c_{A−1}) = H(0x02 ‖ c₀ ‖ … ‖ c_{A−1})
root(b₀…b_{B−1})  = H(0x03 ‖ b₀ ‖ … ‖ b_{B−1})
```

A constant node d levels above the leaves hashes to `EMPTY[d]` or `FULL[d]`; a branch hashes to `node` over its children's hashes; the root is `root` over the B base hashes, each at depth D = R. **The root depends only on the covered leaf set, D, and the profile**: insertion order, duplicates and tiling do not matter. `fingerprint(P, R, profile) = root(cells(P, R))`.

## 5. Openings and set algebra

An **opening** is `{"claim": "empty"|"full"}` or `{"entries": [hex | opening] × A}` (× B at the top). The verifier recomputes bottom-up and compares with the root; the opened cell is the sequence of opened positions and cannot be relabelled. A constant node passes through as A copies of itself. Record: `{"v":1,"hash","A","B","D","root","opening"}`.

Set operations combine two trees of the same profile and depth by node-local rules on hashes (`FULL ∪ X = FULL`, `EMPTY ∪ X = X`; `EMPTY ∩ X = EMPTY`, `FULL ∩ X = X`; `EMPTY \ X = EMPTY`, `X \ FULL = EMPTY`, `X \ EMPTY = X`; equal hashes ⇒ equal sets, `X \ X = EMPTY`), descending only where both are partial and different. A **transcript** records `(h_a, h_b, h_c)` per step with children only at undecided steps; a decided step must have no children. Record: `{"v":1,"op","hash","root_a","root_b","root_c","A","B","D","steps"}`.

## 6. Numerics that must be reproduced

The projection is the rHEALPix projection of the WGS84 authalic sphere as in rhealpixdggs-py 0.8.6, including: the authalic-latitude series (arXiv 2212.05818 A19/A20) with its evaluation order; the polar cap clamp; the pole convention `lon = −π`; the polar-triangle tie-breaks with ε = 10⁻¹⁵ and their north/south comparison senses; `lon_0` applied as a degree shift outside the projection with wrap to `[−180, 180)`. The fixtures pin every one of these.

**The answer is defined to the bit.** Where a result is a float, the defining evaluation is this
crate's: IEEE-754 binary64 arithmetic in the order written, with no fused multiply-add, and every
transcendental function from the `libm` crate (a port of musl's libm), never the platform's. That
evaluation gives the same bits on every target, and `points_burin.json` freezes it: the cell of
every adversarial point at every level (points on cell corners and edges, the seams between base
cells, the band edge and the poles, each with its one-ulp neighbours), the forward projection, and
the nuclei. A conforming implementation reproduces that file exactly, and CI checks it on Linux
(x86_64 and aarch64), macOS (Intel and Apple silicon) and Windows. The reference implementation
evaluates with numpy and the platform's math library, so it is not itself bit-stable across
machines; where it differs from this crate on the committed fixtures, `reference_disagreements.json`
lists the difference exactly (no point lookup differs; two forward projections and about fifty of
2,200 nuclei per profile differ by one or two ulps in latitude), and the parity tests require every
other value to match bit for bit. `polyfill_*.json` is reproduced exactly.

## 7. Correspondence with OGC API - DGGS

The vocabulary of OGC API - DGGS Part 1 (OGC 21-038r1) maps onto this document as follows, so
that a fingerprint can be described in that standard's terms without translation.

| this document | OGC API - DGGS | note |
|---|---|---|
| cell | **zone** | Topic 21 prefers *zone*; this document keeps *cell* only in the code |
| resolution / depth R | **refinement level** (`zone-level`) | level 0 is the six base zones |
| suid `Q453` | **text zone identifier** | the `ogc-rhealpix` profile's identifiers are those of the OGC-registered rHEALPix DGGRS: N,O,P,Q,R,S then digits 0–8 in row-major order. Checked against DGGAL (the OGC reference library) on 66 probes at levels 0–8, polar zones included; the only differences are points lying exactly on a zone edge, where the two libraries break the tie differently. The nucleus rule never places a test point on an edge, so fingerprints are unaffected. |
| canonical cells (§4) | **compact zone list** | "children zones completely covering a parent are replaced by that parent, recursively"; the canonical form of §4 is exactly that list |
| scanline order (§8) | **sub-zone order** | the deterministic order DGGS-JSON zone data lists sub-zones in; for rHEALPix it is scanline, and §8's order is checked against DGGAL on 15 parents, both polar zones included |
| edge neighbours (§8) | zone neighbours | the four zones sharing an edge; checked against DGGAL for every zone at levels 0–2 |
| `ogc-rhealpix` profile | DGGRS `https://www.opengis.net/def/dggrs/OGC/1.0/rHEALPix` | `+proj=rhealpix +lon_0=50 +ellps=WGS84`, refinement ratio 9, `north_square = south_square = 0` |
| profile id | (no equivalent) | binds what the DGGRS leaves open and the root depends on: the hash, the ellipsoid integers, the nucleus rule |
| `burin-1` profile | not registered | rHEALPix with a 0° prime meridian |

**Zone list encoding.** A tree is exchanged as the standard's JSON zone list
(`https://schemas.opengis.net/ogcapi/dggs/part1/1.0/openapi/schemas/dggs-core/dggs-zones.json`):

```json
{"zones": ["Q174783", "Q1747806", …], "returnedAreaMetersSquare": 1.9e7,
 "links": [{"rel": "https://www.opengis.net/def/rel/ogc/1.0/dggrs-definition",
            "href": "https://www.opengis.net/def/dggrs/OGC/1.0/rHEALPix"}]}
```

`zones` is the compact list (mixed levels, each identifier a canonical FULL node), in the order
`Tree::cells` produces. Rebuilding the tree from the list at the declared refinement level and
hashing it reproduces the root; the list is therefore self-verifying against a root and needs
neither the polygon nor the raster it came from. `returnedAreaMetersSquare` is the zone count at
the refinement level times the zone area, which is exact and constant at a level.

**Zone data encoding.** A coverage over the sub-zones of one zone is exchanged as the standard's
DGGS-JSON zone data (`https://schemas.opengis.net/ogcapi/dggs/1.0/core/schemas/dggs-json/dggs-json.json`):

```json
{"dggrs": "https://www.opengis.net/def/dggrs/OGC/1.0/rHEALPix", "zoneId": "Q4", "depths": [2],
 "values": {"coverage": [{"depth": 2, "shape": {"count": 81, "subZones": 81}, "data": [1, null, …]}]}}
```

`data` lists the sub-zones at the stated relative depth in the scanline order of §8: `1` where
covered, `null` elsewhere. A reader takes one field and one depth and covers the sub-zones whose
value is present (not null; optionally, also not zero), which lets any DGGS-JSON field be read as
the extent of its data. It refuses a document whose `dggrs` is not the DGGRS the profile's grid
is (only `ogc-rhealpix` is registered), whose `shape` or `data` length is not `9^depth`, whose
values are not numbers or null, or that carries extra `dimensions`.

## 8. Zone topology

Nothing in this section changes a root. It fixes how cells lie in the plane and which cells
share an edge, so that computation over arrays of cells gives the same answer in every
implementation.

**Position.** Child digit `d = N·row + col` places a child in row `row`, column `col` of its
parent's N×N split (row 0 at the top), so a cell at level r has
`row = Σ ⌊d_i / N⌋·N^(r−i)` and `col = Σ (d_i mod N)·N^(r−i)` inside its base cell, whose side is
`n = N^r`. The cid of `(base, r, row, col)` interleaves them back.

**Scanline sub-zone order.** The sub-zones of a zone Z at relative depth k are ordered by
`index = row·N^k + col`, with row and column relative to Z: top to bottom, then left to right in
the planar layout. This is the sub-zone order of the OGC-registered rHEALPix DGGRS (§7), in which
DGGS-JSON zone data lists its values. Conformance: `dggal_rhealpix.json` (DGGAL's own order, 15
parents including N and S) and `topology.json` (the reference projection's planar cell centres
sorted top to bottom, then left to right, under every shipped profile).

**Edge neighbours.** The base cells lie in the plane as rHEALPix lays them: O, P, Q, R in band
columns 0–3, the band wrapping at both ends; N above column `north_square`, S below column
`south_square`. A cell's neighbours are the cells one step up, right, down and left in its own
planar frame. A step that stays inside the base cell lands on the adjacent row or column. A step
off the base cell, at the stepped point `p = (r, c)`:

| leaving | lands in | at |
|---|---|---|
| band column b, sideways | band column b ± 1 (mod 4) | `(r, c mod n)` |
| band column b, upwards | N | `(n + r, c)` turned counterclockwise `(b − north_square) mod 4` times |
| band column b, downwards | S | `(r − n, c)` turned clockwise `(b − south_square) mod 4` times |
| N through edge k (bottom 0, right 1, top 2, left 3) | band column `north_square + k` | `p` turned clockwise k times, then `(r − n, c)` |
| S through edge k (top 0, right 1, bottom 2, left 3) | band column `south_square + k` | `p` turned counterclockwise k times, then `(r + n, c)` |

A quarter turn of an n×n square is `(r, c) → (c, n−1−r)` clockwise and `(r, c) → (n−1−c, r)`
counterclockwise. Every cell has four neighbours. Conformance: rhealpixdggs-py
`Cell.neighbors(plane=True)`, direction for direction, for every cell at levels 0–2 under all
sixteen polar placements (`topology.json`), and DGGAL's zone neighbours for every zone at levels
0–2 (`dggal_rhealpix.json`).

**Points.** The cell at level r containing a point is found in the plane: project the point
(the forward projection of §6), take the base cell whose square contains it, and read the row and
column off the truncated distances from that square's upper-left corner,
`row = ⌊|y − y0| / w0 · N^r⌋`, `col = ⌊|x − x0| / w0 · N^r⌋`. The polar squares are open, the
equatorial band is closed in y and half-open in x (`x0 ≤ x < x0 + w0`), and a point on an edge
between two cells therefore belongs to the one to its south or east; a distance of exactly `w0`
is nudged in by half a cell width at the reference's finest resolution. This rule alone decides
which cell a point belongs to. Conformance:
`points_*.json`, exact on planar points placed on every kind of edge and vertex, and on
ellipsoidal points more than a micrometre from an edge.

**Polygons.** For drawing and for polygon libraries a cell is a closed counterclockwise ring in
`(lon, lat)`: the four corners of an equatorial cell, whose edges are meridians and parallels;
`n` points per edge of a polar cell; western longitudes moved past 180 on a ring that crosses the
antimeridian; and the cap around each pole closed through the pole along ±180. Polygons are
drawings, not definitions: their vertices come from the inverse projection and their polar edges
are chords of curves, so a point on or very near an edge may fall on the other side of a polygon
from the cell the point rule gives it. Membership is always the point rule's; no root depends on a
polygon.

**Rasters and halos.** A raster at level r is one n×n array per base cell, flattened base-major
then scanline: `index = base·n² + row·n + col`; any per-cell value array moves between cids and
this layout by that index alone. A halo of width `w ≤ n` pads each array to
`(n + 2w)²`; the rules above, applied to a point up to n outside one edge, give the cell drawn at
each padded position. Beyond a base cell's corner, outside two edges at once, there is no cell:
three base cells meet at every corner, so each base cell has four `w × w` blocks without one
(24 cells in all when `w = 1`).

