# Time and space-time: specification (draft)

*Status: **Draft**. §9–§12 (the time line, its hierarchy, the interval rule and relations) are
implemented and now live in SPEC.md. What remains here is not yet implemented; each section moves
into SPEC.md, with fixtures, when its code lands.*

Time is given the same program as space: a hierarchy, a rule that turns a place in it into cells,
a canonical tree and its root, openings, set algebra and relations. Space-time is the spatial tree
with time sets at its leaves, defined so that every spatial root of SPEC §4 is already the
space-time root of "this region, at all times". Existing roots, records and profile ids are
unchanged.

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

**Set algebra.** A space-time transcript is a §5 transcript over the spatial tree in which a step
may also carry `time`:

```json
{"a": hex, "b": hex, "c": hex,
 "time": {"ta": hex, "tb": hex, "tc": hex, "steps": [time step]}}
```

A step has at most one of `children` and `time`. A verifier checks a step at spatial depth `d` as
follows:
1. **Decided by §5.** If the rules of §5 decide it (an operand is `EMPTY[d]` or `FULL[d]` in a
   case the rules cover, or `a = b`), it must have neither `children` nor `time`, and `c` must be
   the rule's result.
2. **Time step.** Otherwise, if it has `time`, then:
   - `a = U_ta[d]` and `b = U_tb[d]`;
   - `{op, ta, tb, tc, steps}` verifies as a §5 transcript over time (`A = 2`, `B = 1`, `D = 61`);
   - `c = U_tc[d]`.
3. **Branch step.** Otherwise it has `children`, checked as in §5. A uniform operand passes
   through as A copies of `U_τ[d−1]`, as constants do in §5.

A prover takes a time step at the first node where both operands are uniform. `FULL ∖ U_τ`, for
example, is a time step with `ta = F_T`, whose result is the complement of τ. The statement a
transcript proves is its op and three roots; two transcripts of one statement may differ in
where they take their time steps.

## 14. Zone data with time

DGGS-JSON (§7) carries time as one entry of `dimensions`:

```json
"dimensions": [{"name": "time", "interval": ["2020-01-01", "2020-04-30"],
                "grid": {"cellsCount": 4, "firstCoordinate": "2020-01-01", "resolution": "P1M"}}]
```

A writer sets the dimension's `definition` to `https://www.opengis.net/def/crs/OGC/0/UnixTime`.
A reader refuses any other `definition`; when the field is absent, the instants are read as
POSIX time.

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

**Versions.** The tag is the version: a changed encoding gets a new tag
(`BURIN-RECORD-2`, …), and no separate version field is needed. `prev` is the id of the key's
previous record, whatever its tag. A chain therefore runs unbroken across format versions, and
each upgrade is itself a signed, chained step.

**The chain is per key.** It orders one signer's records and cannot place them against another
signer's. For that, a record names the other's id in `what`, or both sign into a shared public
log. Either is built from records, and neither needs anything further in this specification.

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
