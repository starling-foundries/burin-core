"""A mini catalog: register scene footprints, query by study area.

Four satellite "scenes" over the lower Thames. The index ranks them by the
exact fraction of the study area each one covers — ground truth, not bbox
hope — and attention() shows where the catalog's effort already sits.
"""

import burin

RESOLUTION = 10
CELL = burin.cell_area_m2(RESOLUTION)


def box(west, south, east, north):
    return {"type": "Polygon", "coordinates": [[
        [west, south], [east, south], [east, north], [west, north], [west, south],
    ]]}


SCENES = {
    "s2-tile-west": box(-0.030, 51.490, -0.008, 51.510),
    "s2-tile-east": box(-0.008, 51.490, 0.014, 51.510),
    "s2-tile-north": box(-0.020, 51.508, 0.010, 51.522),
    "l9-adjacent-path": box(0.050, 51.490, 0.072, 51.510),  # misses the study area
}

index = burin.Index()
for key, footprint in SCENES.items():
    index.add(key, burin.Tree.from_geojson(footprint, RESOLUTION))
print(f"catalog: {len(index)} scenes registered, {index.distinct_roots} distinct coverages")

# The study area: a ~2 km box over Greenwich.
study = burin.Tree.from_geojson(box(-0.015, 51.495, 0.005, 51.515), RESOLUTION)
print(f"study area: {study.leaf_count()} cells, {study.leaf_count() * CELL:,.0f} m²\n")

print(f"{'scene':18s} {'relation':9s} {'of query':>9s} {'of scene':>9s} {'covers':>10s}")
for hit in index.query(study):
    shared_m2 = hit["shared_cells"] * CELL
    print(f"{hit['key']:18s} {hit['relation']:9s} "
          f"{hit['fraction_of_query']:8.1%} {hit['fraction_of_item']:8.1%} "
          f"{shared_m2:>9,.0f} m²")
print("\n(l9-adjacent-path is absent: zero shared cells, so it never enters the shortlist)")

# Where has the catalog already spent effort? Counts per coarse cell.
print("\nattention at resolution 7 (~4.2 km cells):")
for cid, count in sorted(index.attention(7).items()):
    print(f"  {burin.cid_to_suid(cid):10s} {count} scene(s)")
