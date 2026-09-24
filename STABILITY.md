# Stability

burin-core makes two kinds of promise. Its **wire outputs** are values other people store, sign,
compare and recompute elsewhere, so they are frozen. Its **API** is ordinary code and follows
semantic versioning.

## Wire outputs

These are defined by [SPEC.md](SPEC.md) and pinned by the fixtures under
`crates/burin-core/tests/fixtures`. CI checks them on every supported platform.

| output | defined in | pinned by |
|---|---|---|
| the profile identity | SPEC §1 | `profile_ids.json`, `roots.json` |
| the cid of a path, and the suid | SPEC §2 | `topology.json`, `dggal_rhealpix.json` |
| the cells of a polygon at a resolution | SPEC §3 | `polyfill_*.json` |
| the 32-byte root of a leaf set | SPEC §4 | `roots.json` |
| the opening record and the set-operation proof, `v: 1` | SPEC §5 | `roots.json` |
| the cell containing a point, the forward projection and the nuclei, to the bit | SPEC §6 | `points_burin.json` |
| the sub-zone order, the raster layout and DGGS-JSON zone data | SPEC §7–8 | `topology.json`, `dggal_rhealpix.json` |
| the order of the edge neighbours (up, right, down, left) | SPEC §8 | `topology.json`, `dggal_rhealpix.json` |

A wire output changes only under a new name: a record version (`v: 2`), a new hash id in the
profile, or a new profile. A release never changes what an existing name produces. If one ever
must, because an output is found to be wrong, the release says so in the changelog under
**Changed output**, and the old output gets a new name.

A reader may become **stricter**. It may refuse input it used to accept, when that input was
malformed or could be spelled two ways. A record this crate wrote is never refused by a later
release that reads the same version. Each refusal is listed in the changelog under **Refused**.

### Not wire outputs

- **Cell polygons** (`cell_polygon`, `cell_boundaries`, `cells_geojson`) are drawings. Their
  vertex count, densification and antimeridian handling may change. Membership is decided by
  the point rule.
- **Error messages** may change.
- **Performance and threading** may change. The parallel build hashes exactly what the serial
  one does.

## API

The version is `MAJOR.MINOR.PATCH`. Rust and Python share one version.

- **While the major version is 0**, a minor release may change the Rust or Python API and says
  how in the changelog. A patch release never does.
- **From 1.0**, the API changes incompatibly only in a major release.
- **Supported surface.** It is what `burin_core` re-exports at its crate root, and every name in
  `burin` that does not begin with an underscore. Private names (`burin._burin._cells_*`, and so
  on) may change in any release.
- **Minimum versions.** Raising the minimum Python (now 3.9) or Rust (stable) is a minor
  release.
