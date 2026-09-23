"""Profiles: the same polygon under two grids is two different identities.

A root is meaningless without the profile it is read against. The profile's
identity is derived from its integer parameters — grid, hash, epoch — and is
never transmitted: both sides derive it and compare.
"""

import burin

RESOLUTION = 10

greenwich = {"type": "Polygon", "coordinates": [[
    [-0.010, 51.500], [0.000, 51.500], [0.000, 51.508], [-0.010, 51.508], [-0.010, 51.500],
]]}

ogc = burin.Profile.ogc()       # OGC API-DGGS Annex B rHEALPix, lon_0 = 50 (the default)
burin1 = burin.Profile.burin_1()  # rHEALPix with a 0° prime meridian, lon_0 = 0

print(f"Profile.ogc():     {ogc}")
print(f"  id: {ogc.id_hex()}")
print(f"Profile.burin_1(): {burin1}")
print(f"  id: {burin1.id_hex()}")

root_ogc = burin.fingerprint_polygon(greenwich, RESOLUTION, ogc)
root_b1 = burin.fingerprint_polygon(greenwich, RESOLUTION, burin1)
print(f"\nsame polygon, resolution {RESOLUTION}:")
print(f"  under ogc:     {root_ogc[:24]}…  ({burin.Tree.from_geojson(greenwich, RESOLUTION, ogc).leaf_count()} cells)")
print(f"  under burin_1: {root_b1[:24]}…  ({burin.Tree.from_geojson(greenwich, RESOLUTION, burin1).leaf_count()} cells)")
print(f"  roots equal? {root_ogc == root_b1}")

print(f"\nomitting the profile uses OGC: "
      f"{burin.fingerprint_polygon(greenwich, RESOLUTION) == root_ogc}")
print("the grid parameters are inside the identity, so a root can never be "
      "silently read against the wrong grid.")
