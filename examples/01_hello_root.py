"""Hello, root: a polygon becomes 32 bytes anyone can recompute.

A ~1 km box in Greenwich, London, at resolution 10 (~156 m cells). No files,
no network: the fingerprint is a pure function of the polygon, the resolution,
and the profile.
"""

import burin

RESOLUTION = 10

# A ~1 km box over Greenwich Park, London (lon/lat, GeoJSON order).
greenwich = {
    "type": "Polygon",
    "coordinates": [[
        [-0.010, 51.500],
        [0.000, 51.500],
        [0.000, 51.508],
        [-0.010, 51.508],
        [-0.010, 51.500],
    ]],
}

root = burin.fingerprint_polygon(greenwich, RESOLUTION)
print(f"polygon: a ~1 km box over Greenwich, resolution {RESOLUTION}")
print(f"root:    {root}")
print("         64 hex chars. Anyone with the same polygon, resolution and")
print("         profile recomputes these exact bytes, in any language.")

tree = burin.Tree.from_geojson(greenwich, RESOLUTION)
area = tree.leaf_count() * burin.cell_area_m2(RESOLUTION)
print(f"\ncoverage: {tree.leaf_count()} cells of {burin.cell_area_m2(RESOLUTION):,.0f} m² each"
      f" (stored as {len(tree.cells())} coarse cells)")
print(f"exact area: {area:,.0f} m² — equal-area grid, so this is a count times a constant")

# The root does not care how the ring is spelled.
ring = greenwich["coordinates"][0]
reversed_ring = {"type": "Polygon", "coordinates": [list(reversed(ring))]}
rotated_ring = {"type": "Polygon", "coordinates": [ring[2:] + ring[1:3]]}
print(f"\nsame ring, vertices reversed: same root? "
      f"{burin.fingerprint_polygon(reversed_ring, RESOLUTION) == root}")
print(f"same ring, start vertex moved:  same root? "
      f"{burin.fingerprint_polygon(rotated_ring, RESOLUTION) == root}")

# But the resolution is part of the identity.
print(f"\nat resolution {RESOLUTION - 1} instead: "
      f"{burin.fingerprint_polygon(greenwich, RESOLUTION - 1)[:16]}…  (different — "
      "resolution is part of the identity)")
