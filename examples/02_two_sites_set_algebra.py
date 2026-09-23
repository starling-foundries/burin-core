"""Two sites, set algebra: union, intersect, difference over 32-byte roots.

Two overlapping survey plots. Their coverages combine without ever
re-exchanging the geometry — the algebra runs on the trees themselves.
"""

import burin

RESOLUTION = 10
CELL = burin.cell_area_m2(RESOLUTION)


def box(west, south, east, north):
    return {"type": "Polygon", "coordinates": [[
        [west, south], [east, south], [east, north], [west, north], [west, south],
    ]]}


# Two overlapping plots near Greenwich: plot_b straddles plot_a's NE corner.
plot_a = box(-0.010, 51.500, 0.000, 51.508)
plot_b = box(-0.005, 51.504, 0.005, 51.512)

a = burin.Tree.from_geojson(plot_a, RESOLUTION)
b = burin.Tree.from_geojson(plot_b, RESOLUTION)

u = a.union(b)
i = a.intersect(b)
d = a.difference(b)

print(f"plot A: {a.leaf_count()} cells, {a.leaf_count() * CELL:,.0f} m²")
print(f"plot B: {b.leaf_count()} cells, {b.leaf_count() * CELL:,.0f} m²")
print(f"\nA ∪ B: {u.leaf_count()} cells, {u.leaf_count() * CELL:,.0f} m²")
print(f"A ∩ B: {i.leaf_count()} cells, {i.leaf_count() * CELL:,.0f} m²")
print(f"A ∖ B: {d.leaf_count()} cells, {d.leaf_count() * CELL:,.0f} m²")
print(f"\nconservation: |A| + |B| = |A ∪ B| + |A ∩ B|? "
      f"{a.leaf_count() + b.leaf_count() == u.leaf_count() + i.leaf_count()}")

# Roots are the identity: the intersection is the same object from either side.
print(f"commutativity: root(A ∩ B) == root(B ∩ A)? "
      f"{a.intersect(b).root_hex == b.intersect(a).root_hex}")

# Where do the two plots actually disagree? Maximal subtrees, not every leaf.
div = a.divergence(b)
print(f"\ndivergence: {len(div)} maximal cells where A and B differ "
      f"(vs {u.leaf_count() - i.leaf_count()} differing leaves):")
for entry in div[:5]:
    print(f"  {entry['suid']:12s} A:{entry['hash_a'][:12]}…  B:{entry['hash_b'][:12]}…")
