# Time and space-time: specification (draft)

*Status: **Draft**. Nothing here is implemented. Each section moves into SPEC.md, with fixtures,
when its code lands. Section numbers continue SPEC.md's.*

Time is given the same program as space: a hierarchy, a rule that turns a place in it into cells,
a canonical tree and its root, openings, set algebra and relations. Space-time is the spatial tree
with time sets at its leaves, defined so that every spatial root of SPEC §4 is already the
space-time root of "this region, at all times". Existing roots, records and profile ids are
unchanged.

## 9. The time line

**Scale.** Instants are POSIX time: microseconds since 1970-01-01T00:00:00Z, leap seconds not
counted. A leap second (`23:59:60`) has no POSIX instant and is refused, as is any calendar
without real instants (CF `noleap`, `360_day` and similar).

**Ticks.** Under a profile (§1), an instant `t` falls in the tick

```
tick(t) = floor((t − epoch_us) / tick_us) + 2^60        computed exactly (128-bit), floor toward −∞
```

and is refused unless `0 ≤ tick(t) < 2^61`. With the default `tick_us = 1`, `epoch_us = 0`, the
line is 1 µs ticks centred on 1970, about 36,500 years either side. The profile's time fields
were reserved at these values, so every profile id stays as it is. A profile with
`tick_us = 86 400 000 000` has day ticks; its id differs, so day-tick and µs-tick records cannot
be confused.

An instant given as text is RFC 3339 with a `Z` or numeric offset and at most 6 fractional digits.
A date alone (`2020-01-01`) is the start of that day in UTC. `datetime64[ns]` converts by
`floor(ns / 1000)`. A float is never accepted: time is integer arithmetic throughout, so every
result below is exact on every platform.

## 10. The time hierarchy

Time is the hierarchy **H(2, 1)** of §2 with the fixed depth **D_T = 61**. Cell `k` at level `r`
(`0 ≤ k < 2^r`) is the ticks `[k·2^(61−r), (k+1)·2^(61−r))`, and its cid is `2^(r+1) + k`. A
level-61 cid is `2^62 + tick`. Parent, children and descendants are those of §2: coarsening an
instant to level `r` is `cid >> (61 − r)`. Neighbours are `cid ± 1` within the level.

The resolution and span match the IVOA Time-MOC (MOC 2.0: 1 µs cells, order 61). The time scale
and origin differ, since the Time-MOC counts from Julian Day 0 in TCB, so the two agree cell for
cell only at the tick level, through a fixed offset.

## 11. The interval rule

An **interval** is `[lo, hi)` in ticks, `0 ≤ lo < hi ≤ 2^61`. Its cells are its ticks. Held
canonically (§4), those ticks are the compact cell list, at most `2·61` cells. The HINT boundary
walk produces that list directly:

```
walk(lo, hi, level = 61):   while lo < hi:
    if lo is odd:  emit (level, lo);  lo += 1
    if hi is odd:  hi -= 1;  emit (level, hi)
    lo, hi, level = lo/2, hi/2, level − 1
```

A **time set** is any union of intervals. Its **time root** is the root (§4) of its time tree at
depth 61, under the profile's hash, with `A = 2`, `B = 1`. The time root depends on the set
alone, not on how it was split into intervals. Because the depth is fixed, a set has exactly one
time root, and no resolution parameter enters it.

The empty set and the whole line have the constant roots `E_T` and `F_T`.

**Coarse rule.** Where a set is wanted at level `r < 61`, a level-`r` cell belongs to
`cells(I, r)` iff its midpoint tick lies in `I`. This is the time counterpart of the nucleus rule
of §3. It is optional: the exact set (level 61) is the default.

## 12. Relations

**Between two intervals.** Allen's thirteen relations are decided by comparing endpoints. For
non-empty `a = [a₁, a₂)` and `b = [b₁, b₂)`:

| relation | condition | inverse | condition |
|---|---|---|---|
| before | `a₂ < b₁` | after | `b₂ < a₁` |
| meets | `a₂ = b₁` | met-by | `b₂ = a₁` |
| overlaps | `a₁ < b₁ < a₂ < b₂` | overlapped-by | `b₁ < a₁ < b₂ < a₂` |
| starts | `a₁ = b₁ ∧ a₂ < b₂` | started-by | `a₁ = b₁ ∧ b₂ < a₂` |
| during | `b₁ < a₁ ∧ a₂ < b₂` | contains | `a₁ < b₁ ∧ b₂ < a₂` |
| finishes | `b₁ < a₁ ∧ a₂ = b₂` | finished-by | `a₁ < b₁ ∧ a₂ = b₂` |
| equals | `a₁ = b₁ ∧ a₂ = b₂` | | |

With half-open intervals, `meets` means adjacent and sharing no tick. A signed record states its
interval's endpoints (§15), so any reader can decide all thirteen against any other interval from
the record alone.

**Between two sets.** Space, time and space-time sets share five relations, each proven by one
set-operation transcript (§5) of the same hierarchy and depth. `∅` is the empty tree's root.

| relation | holds iff | evidence |
|---|---|---|
| equals | `root(X) = root(Y)` | none |
| disjoint | `X ∩ Y = ∅` | intersect transcript with `root_c = ∅` |
| intersects | `X ∩ Y ≠ ∅` | intersect transcript with `root_c ≠ ∅` |
| within | `X ∖ Y = ∅` | difference transcript `X ∖ Y` with `root_c = ∅` |
| contains | `Y ∖ X = ∅` | difference transcript `Y ∖ X` with `root_c = ∅` |

Only these lift from intervals to sets: whether two sets *meet* or one *starts* the other is not a
property of their cells alone. Where Allen and set relations share a name (`equals`, `contains`),
the Allen relation implies the set relation.

## 13. The space-time tree

A space-time set is a set of pairs (depth-R spatial cell, tick). It is held as the spatial tree of
§4 at depth R whose leaves carry time sets, so space is the outer structure, as in DGGS-JSON zone
data (§14) and in datasets indexed by zone with time as a dimension.

A spatial node whose descendant leaves all carry the same time set τ is **uniform in τ**. It is
held as one constant, with the hash

```
U_τ[d] = EMPTY[d]                    if τ = E_T        (§4 ladders, unchanged)
       = FULL[d]                     if τ = F_T
       = H(0x04 ‖ τ)                 if d = 0
       = node(U_τ[d−1] × A)          otherwise
```

and every other node is a branch, hashed by `node` over its children as in §4. The root is `root`
over the B base hashes. Canonical form is §4's, with "constant" read as "uniform in some τ".

Three things follow:
- **Spatial trees embed.** A spatial tree is a space-time tree whose leaves are all empty or
  all-time. Its root is unchanged, and it means that region at all times.
- **Records are refused if they spell τ two ways.** A record that claims `H(0x04 ‖ τ)` with
  `τ ∈ {E_T, F_T}` has a second spelling of `EMPTY[0]` or `FULL[0]`, and is refused.
- **Products are computable.** The product `S × τ` of a spatial tree and a time set is `S` with
  every `FULL[d]` node replaced by `U_τ[d]`. It costs one pass over S's canonical cells, so a
  record carrying a *where* and a *when* already commits to their product (§15).

**Projections.** The **footprint** of a space-time set is its spatial cells whose time set is not
empty. Its **span** is the union of its time sets. Both are ordinary spatial and time trees.

**Openings.** A space-time opening is a spatial opening (§5) whose terminal claim is `empty`,
`full`, or `{"claim": "timed", "time_root": τ, "time": …}`, where `time` is a time opening (§5,
`A = 2`, `B = 1`, `D = 61`) against τ. The spatial recomputation takes `U_τ[d]` at the claim. The
time opening then proves a tick, or the ticks of an interval, in or out of τ.

**Set algebra.** The node-local rules of §5 hold with EMPTY and FULL read as `U_{E_T}` and
`U_{F_T}`. Two uniform nodes combine into `U_{τ ⊕ σ}`, and a step there carries a nested time
transcript proving `τ ⊕ σ`. The transcript format is to be fixed with the implementation.

## 14. Zone data with time

DGGS-JSON (§7) carries time as one entry of `dimensions`:

```json
"dimensions": [{"name": "time", "interval": ["2020-01-01", "2020-04-30"],
                "grid": {"cellsCount": 4, "firstCoordinate": "2020-01-01", "resolution": "P1M"}}]
```

The data array is ordered sub-zone major (OGC API - DGGS, clause 16), so the value for sub-zone
`z` at step `i` is `data[z · T + i]`, with `T = cellsCount`. This is the space-major order of
§13.

**Step intervals.** Step `i` is an interval of §11:
- **Regular grid:** `[first + i·resolution, first + (i+1)·resolution)`, where durations are added
  to instants by the algorithm of XML Schema 1.1 Part 2, Appendix E. That algorithm is calendar
  arithmetic, in which a month is a calendar month.
- **Irregular grid:** `boundsCoordinates[i]` where given. Otherwise `coordinates[i]` is an
  instant, one tick.

The document's `interval` must agree with the grid: its start is the first step's start and its
end lies within the last step. Otherwise it is refused.

**Reading.** A reader takes one field and covers (sub-zone, step interval) wherever the value is
present, as in §7. The result is a space-time tree. A dimension other than time is refused,
because a sub-zone would then have more than one value per step.

## 15. The record

A record is one signed message stating where, when, and what. Each component can be verified on
its own, and together they are a space-time commitment.

| field | content |
|---|---|
| `profile` | the profile's fields (§1); the reader re-derives its id |
| `key` | Ed25519 public key (RFC 8032), 32 bytes |
| `seq` | u64; this key's records are numbered 0, 1, 2, … |
| `prev` | the id of this key's record `seq − 1`; 32 zero bytes when `seq = 0` |
| `issued` | tick (§9) at which the record was signed; never decreases along `seq` |
| `what` | 32 bytes chosen by the signer, typically a content digest; never interpreted |
| `where` | absent, a cell (cid), or a spatial root with its depth R |
| `when` | absent, an interval `[lo, hi)`, or a time root |
| `joint` | absent, or a space-time root (§13) when the set is not `where × when` |

When `joint` is present, `where` and `when` are its footprint and span. A reader who needs one axis
checks only that component. A reader of both recomputes the product, or checks `joint`.

**Two times.** `issued` is when the signer signed. `when` is the time the content concerns, which
may be long past. The chain constrains only `issued`.

**Encoding.** The signed bytes are the tag `BURIN-RECORD-1`, then the profile identity preimage
(§1), then each field in the order above. Integers are big-endian. Each optional field is
preceded by a kind byte: 0 means absent; 1 and 2 are the forms in the table, in order. The **id**
is SHA-256 of those bytes. The **signature** is Ed25519 over them. JSON records spell each field as
the other records of this specification do, with lowercase hex.

**Evidence of misbehaviour.** Two validly signed records under one key are self-contained proof:
- **equivocation:** equal `seq`, different ids;
- **backdating:** `seq_a < seq_b` but `issued_b < issued_a`.

Withheld records cannot be detected. That limit is inherent to any signed log.

**Compact forms.**
- **Words:** the id as 24 BIP-39 words, with a SHA-256 checksum byte.
- **Burst:** 70 bytes, `version(1) ‖ flags(1) ‖ seq(8) ‖ issued(8) ‖ keyid(20) ‖ id(32)`. Here
  `keyid` is the first 20 bytes of SHA-256(key), and `flags` has one presence bit each for
  where, when and joint.

Both identify a record; neither replaces it.

## Open questions

1. **Instants or intervals.** Should a time step without bounds (a CF time coordinate with no
   `bounds`) mean an instant (one tick, as above), or run until the next step?
2. **Coarse rule.** Is the midpoint rule for coarse time cells wanted in the standard, or only as a
   convenience?
3. **Space-time transcripts.** The nested transcript format of §13.
4. **Record version.** Whether records carry an explicit version (`v`) beside the tag.
5. **OGC Part 3.** OGC Topic 21 Part 3 (spatio-temporal DGGS, OGC 20-048r0, in draft) has not been
   read. Whether it constrains time hierarchies should be checked before §10 is final.
6. **Registered name.** A registered name for the time line, as §7 names the DGGRS.
