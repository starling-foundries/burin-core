#!/usr/bin/env python3
"""Generate the parity fixtures (set 1) from the reference implementation, rhealpixdggs-py.

Run:  uv run --with /path/to/rhealpixdggs-py --with shapely python tools/gen_fixtures.py
The output under crates/burin-core/tests/fixtures/ is committed; re-running must produce an
empty diff. Every value here is what the reference computes; the Rust crate must agree.
"""
from __future__ import annotations

import hashlib
import json
import random
import struct
import sys
from pathlib import Path

from numpy import pi
from shapely.geometry import MultiPolygon, Point, Polygon, box
from shapely.prepared import prep

from rhealpixdggs import pj_healpix as H
from rhealpixdggs import pj_rhealpix as R
from rhealpixdggs.dggs import RHEALPixDGGS
from rhealpixdggs.ellipsoids import WGS84_E, Ellipsoid
from rhealpixdggs.utils import auth_lat, auth_rad

OUT = Path(__file__).resolve().parent.parent / "crates" / "burin-core" / "tests" / "fixtures"
OUT.mkdir(parents=True, exist_ok=True)
A_WGS84, INV_F = 6378137.0, 298.257223563
LETTERS = "NOPQRS"


def dump(name: str, obj) -> None:
    (OUT / name).write_text(json.dumps(obj, indent=1, sort_keys=True) + "\n")
    print(f"wrote {name}")


def f(x) -> float:
    return float(x)


# ----------------------------------------------------------------------------- profiles
PROFILES = {
    "ogc": dict(lon_0=50.0, ns=0, ss=0),
    "burin1": dict(lon_0=0.0, ns=0, ss=0),
    "ns1ss2": dict(lon_0=0.0, ns=1, ss=2),
}


def dggs_for(p) -> RHEALPixDGGS:
    ell = Ellipsoid(a=A_WGS84, f=1.0 / INV_F, lon_0=p["lon_0"])
    return RHEALPixDGGS(ellipsoid=ell, north_square=p["ns"], south_square=p["ss"], N_side=3)


def profile_id(hash_id: str, A: int, B: int, lon_0_udeg: int, ns: int, ss: int,
               a_um: int, inv_f_nano: int, tick_us: int, epoch_us: int) -> str:
    pre = b"BURIN-PROFILE-1" + bytes([len(hash_id)]) + hash_id.encode()
    pre += struct.pack(">IIqBBQQQq", A, B, lon_0_udeg, ns, ss, a_um, inv_f_nano, tick_us, epoch_us)
    return hashlib.sha256(pre).hexdigest()


# ----------------------------------------------------------------------------- 1. auth_lat
def gen_auth_lat():
    phis = [0.0, pi / 6, -pi / 6, pi / 4, -pi / 4, pi / 3, -pi / 3, pi / 2, -pi / 2]
    phis += [-pi / 2 + k * pi / 40 for k in range(1, 40)]
    rows = []
    for e in (f(WGS84_E), 0.8, 0.0):
        for phi in phis:
            for inv in (False, True):
                rows.append({"phi": phi, "e": e, "inverse": inv,
                             "value": f(auth_lat(phi, e, inverse=inv, radians=True))})
    dump("auth_lat.json", {"rows": rows, "auth_rad_wgs84": f(auth_rad(A_WGS84, f(WGS84_E))),
                           "e_wgs84": f(WGS84_E)})


# ----------------------------------------------------------------------------- 2. healpix
def gen_healpix():
    lams = [-pi, -3 * pi / 4, -pi / 2, -pi / 4, 0.0, pi / 4, pi / 2, 3 * pi / 4, pi - 1e-12,
            -pi + 1e-12, -2.9, -1.3, 0.7, 2.2, 3.0]
    phis = [-pi / 2, -pi / 3, -pi / 4, -0.5, 0.0, 0.3, pi / 4, pi / 3, pi / 2, 0.7297276562, -0.7297276562,
            0.9, -1.2, 1.5]
    fwd, inv = [], []
    for e in (0.0, f(WGS84_E), 0.8):
        for lam in lams:
            for phi in phis:
                x, y = H.healpix_ellipsoid(lam, phi, e=e) if e else H.healpix_sphere(lam, phi)
                fwd.append({"lam": lam, "phi": phi, "e": e, "x": f(x), "y": f(y)})
                lam2, phi2 = H.healpix_ellipsoid_inverse(f(x), f(y), e=e) if e else H.healpix_sphere_inverse(f(x), f(y))
                inv.append({"x": f(x), "y": f(y), "e": e, "lam": f(lam2), "phi": f(phi2)})
    # pole rounding cases pinned exactly (lon = -pi)
    poles = []
    for (lam, phi) in ((-7 * pi / 8, 3 * pi / 8), (-5 * pi / 6, 5 * pi / 12)):
        x, y = H.healpix_sphere(lam, phi)
        lam2, phi2 = H.healpix_sphere_inverse(f(x), f(y))
        poles.append({"x": f(x), "y": f(y), "lam": f(lam2), "phi": f(phi2)})
    for (x, y) in ((-3 * pi / 4, pi / 2), (pi / 4, -pi / 2), (-3 * pi / 4, pi / 2 + 5e-11), (3 * pi / 4, -pi / 2 - 5e-11)):
        lam2, phi2 = H.healpix_sphere_inverse(f(x), f(y))
        poles.append({"x": f(x), "y": f(y), "lam": f(lam2), "phi": f(phi2)})
    image = []
    corners = [(-pi, pi / 4), (-3 * pi / 4, pi / 2), (-pi / 2, pi / 4), (-pi / 4, pi / 2), (0, pi / 4),
               (pi / 4, pi / 2), (pi / 2, pi / 4), (3 * pi / 4, pi / 2), (pi, pi / 4), (pi, -pi / 4),
               (3 * pi / 4, -pi / 2), (pi / 2, -pi / 4), (pi / 4, -pi / 2), (0, -pi / 4),
               (-pi / 4, -pi / 2), (-pi / 2, -pi / 4), (-3 * pi / 4, -pi / 2), (-pi, -pi / 4), (0, 0)]
    for (x, y) in corners:
        for eps in (0.0, 5e-11, -5e-11, 2e-10, -2e-10, 0.1):
            for (dx, dy) in ((eps, 0.0), (0.0, eps), (eps, eps)):
                px, py = f(x + dx), f(y + dy)
                image.append({"x": px, "y": py, "inside": bool(H.in_healpix_image(px, py))})
    dump("healpix.json", {"forward": fwd, "inverse": inv, "poles": poles, "image": image})


# ----------------------------------------------------------------------------- 3. rhealpix
def gen_rhealpix():
    lams = [-pi, -3 * pi / 4 + 1e-3, -pi / 2, -pi / 4, 0.0, pi / 4, pi / 2, 3 * pi / 4, pi - 1e-9, -2.0, 1.1, 2.9]
    phis = [-pi / 2, -1.3, -pi / 3, -0.75, 0.0, 0.75, pi / 3, 1.3, pi / 2, 0.7297276562 + 1e-9]
    combine, tri, ell = [], [], []
    for ns in range(4):
        for ss in range(4):
            for lam in lams:
                for phi in phis:
                    hx, hy = H.healpix_sphere(lam, phi)
                    hx, hy = f(hx), f(hy)
                    t, region = R.triangle(hx, hy, ns, ss, inverse=False)
                    tri.append({"x": hx, "y": hy, "ns": ns, "ss": ss, "inverse": False, "t": t, "region": region})
                    rx, ry = R.combine_triangles(hx, hy, ns, ss)
                    rx, ry = f(rx), f(ry)
                    combine.append({"x": hx, "y": hy, "ns": ns, "ss": ss, "inverse": False, "rx": rx, "ry": ry})
                    t2, region2 = R.triangle(rx, ry, ns, ss, inverse=True)
                    tri.append({"x": rx, "y": ry, "ns": ns, "ss": ss, "inverse": True, "t": t2, "region": region2})
                    ux, uy = R.combine_triangles(rx, ry, ns, ss, inverse=True)
                    combine.append({"x": rx, "y": ry, "ns": ns, "ss": ss, "inverse": True, "rx": f(ux), "ry": f(uy)})
            # points on and just off the polar-square diagonals (the inverse tie-breaks)
            for sq, sgn in ((ns, 1.0), (ss, -1.0)):
                cx = -pi + sq * pi / 2 + pi / 4          # square centre x
                cy = sgn * pi / 2                        # square centre y
                for u in (-0.3, -0.1, 0.1, 0.3):
                    for (dx, dy) in ((u, u), (u, -u)):
                        for eps in (0.0, 1e-16, -1e-16, 1e-14, -1e-14, 3e-15, -3e-15):
                            px, py = f(cx + dx + eps), f(cy + dy)
                            t2, region2 = R.triangle(px, py, ns, ss, inverse=True)
                            tri.append({"x": px, "y": py, "ns": ns, "ss": ss, "inverse": True, "t": t2, "region": region2})
                            if R.in_rhealpix_image(px, py, ns, ss):
                                ux, uy = R.combine_triangles(px, py, ns, ss, inverse=True)
                                combine.append({"x": px, "y": py, "ns": ns, "ss": ss, "inverse": True, "rx": f(ux), "ry": f(uy)})
    image = []
    for ns in range(4):
        for ss in range(4):
            pts = [(-pi, pi / 4), (-pi + ns * pi / 2, pi / 4), (-pi + ns * pi / 2, 3 * pi / 4),
                   (-pi + (ns + 1) * pi / 2, 3 * pi / 4), (-pi + (ns + 1) * pi / 2, pi / 4), (pi, pi / 4),
                   (pi, -pi / 4), (-pi + (ss + 1) * pi / 2, -pi / 4), (-pi + (ss + 1) * pi / 2, -3 * pi / 4),
                   (-pi + ss * pi / 2, -3 * pi / 4), (-pi + ss * pi / 2, -pi / 4), (-pi, -pi / 4), (0, 0)]
            for (x, y) in pts:
                for eps in (0.0, 5e-16, -5e-16, 2e-15, -2e-15, 0.1):
                    for (dx, dy) in ((eps, 0.0), (0.0, eps), (eps, eps)):
                        px, py = f(x + dx), f(y + dy)
                        image.append({"x": px, "y": py, "ns": ns, "ss": ss,
                                      "inside": bool(R.in_rhealpix_image(px, py, ns, ss))})
    e = f(WGS84_E)
    for ns, ss in ((0, 0), (1, 2), (3, 3)):
        for lam in lams:
            for phi in phis:
                x, y = R.rhealpix_ellipsoid(lam, phi, e=e, north_square=ns, south_square=ss)
                x, y = f(x), f(y)
                lam2, phi2 = R.rhealpix_ellipsoid_inverse(x, y, e=e, north_square=ns, south_square=ss)
                ell.append({"lam": lam, "phi": phi, "ns": ns, "ss": ss, "x": x, "y": y,
                            "lam_back": f(lam2), "phi_back": f(phi2)})
    dump("rhealpix.json", {"triangle": tri, "combine": combine, "image": image, "ellipsoid": ell, "e": e})


# ----------------------------------------------------------------------------- 4. cells
def all_suids(res: int):
    out = []
    for b in LETTERS:
        stack = [(b,)]
        while stack:
            s = stack.pop()
            if len(s) - 1 == res:
                out.append(s)
            else:
                stack.extend(s + (k,) for k in range(9))
    return out


def suid_str(s) -> str:
    return s[0] + "".join(str(d) for d in s[1:])


def gen_cells(name: str, p):
    d = dggs_for(p)
    rng = random.Random(20260916)
    suids = []
    for r in (0, 1, 2):
        suids += all_suids(r)
    suids += [s for s in all_suids(3) if s[0] in "NS"]
    for _ in range(200):
        r = rng.randint(4, 10)
        suids.append((rng.choice(LETTERS), *[rng.randint(0, 8) for _ in range(r)]))
    rows = []
    for s in suids:
        c = d.cell(list(s))
        ul = c.ul_vertex(plane=True)
        nuc_p = c.nucleus(plane=True)
        nuc = c.nucleus(plane=False)
        bnd = c.boundary(n=3, plane=False)
        rows.append({
            "suid": suid_str(s),
            "ul": [f(ul[0]), f(ul[1])],
            "nucleus_planar": [f(nuc_p[0]), f(nuc_p[1])],
            "nucleus": [f(nuc[0]), f(nuc[1])],
            "boundary3": sorted([[f(q[0]), f(q[1])] for q in bnd]),
            "region": c.region(),
            "shape": c.ellipsoidal_shape,
        })
    seam = []
    for r in range(5):
        for s in all_suids(r):
            lon, lat = d.cell(list(s)).nucleus(plane=False)
            if float(lon) == -180.0:
                seam.append({"suid": suid_str(s), "nucleus": [f(lon), f(lat)]})
    dump(f"cells_{name}.json", {
        "profile": p,
        "seam_nuclei": seam,
        "cell_width": [f(d.cell_width(r)) for r in range(16)],
        "cell_area": [f(d.cell_area(r, plane=False)) for r in range(16)],
        "cells": rows,
    })


# ----------------------------------------------------------------------------- 5. polyfill (nucleus rule)
def _lon_range(lons):
    lo, hi = min(lons), max(lons)
    if hi - lo <= 180.0:
        return lo, hi
    shifted = [x + 360.0 if x < 0.0 else x for x in lons]
    slo, shi = min(shifted), max(shifted)
    if shi - slo <= 180.0:
        return slo, shi
    return -180.0, 180.0


def _bbox_geom(lon_lo, lat_lo, lon_hi, lat_hi):
    if lon_hi <= 180.0:
        return box(lon_lo, lat_lo, lon_hi, lat_hi)
    return MultiPolygon([box(lon_lo, lat_lo, 180.0, lat_hi), box(-180.0, lat_lo, lon_hi - 360.0, lat_hi)])


class Grid:
    """The conservative lon/lat bounding box of a cell (mirrors polyfill.rs `cell_bbox`).

    Equatorial cells: longitude depends only on x and latitude only on y, so the 8 boundary
    samples (corners included) bound the cell exactly. Polar cells: latitude is a function of the
    L-infinity distance from the polar square's centre, so the minimum latitude is at the farthest
    corner and the maximum at the nearest point of the planar rectangle (analytic); longitude
    is monotone along the perimeter as seen from the centre, so its extremes are at corners; the
    cap cell (centre inside) spans every longitude."""

    def __init__(self, d: RHEALPixDGGS):
        self.d = d
        self.cap = min(float(pt[1]) for pt in d.cell(["N"]).boundary(n=3, plane=False))
        self.centre = {b: float(d.cell([b]).nucleus(plane=False)[0]) for b in "OPQR"}
        w0 = d.cell_width(0)
        self.pole = {b: (d.ul_vertex[b][0] + w0 / 2, d.ul_vertex[b][1] - w0 / 2) for b in "NS"}

    def bbox(self, suid):
        base = suid[0]
        if len(suid) == 1:
            if base == "N":
                return (-180.0, self.cap, 180.0, 90.0)
            if base == "S":
                return (-180.0, -90.0, 180.0, -self.cap)
            c = self.centre[base]
            lo, hi = c - 45.0, c + 45.0
            if lo < -180.0:
                lo, hi = lo + 360.0, hi + 360.0
            return (lo, -self.cap, hi, self.cap)
        cell = self.d.cell(list(suid))
        if base in "OPQR":
            pts = cell.boundary(n=3, plane=False)
            lats = [float(q[1]) for q in pts]
            lon_min, lon_max = _lon_range([float(q[0]) for q in pts])
            return lon_min, min(lats), lon_max, max(lats)
        x0, y0 = cell.ul_vertex(plane=True)
        w = cell.width()
        x1, y1 = x0 + w, y0 - w
        cx, cy = self.pole[base]
        if x0 < cx < x1 and y1 < cy < y0:
            return (-180.0, self.cap, 180.0, 90.0) if base == "N" else (-180.0, -90.0, 180.0, -self.cap)
        corners = [self.d.rhealpix(x, y, inverse=True) for (x, y) in ((x0, y0), (x1, y0), (x1, y1), (x0, y1))]
        lon_min, lon_max = _lon_range([float(q[0]) for q in corners])
        lat_far = min(float(q[1]) for q in corners) if base == "N" else max(float(q[1]) for q in corners)
        nx, ny = min(max(cx, x0), x1), min(max(cy, y1), y0)
        lat_near = float(self.d.rhealpix(nx, ny, inverse=True)[1])
        if base == "N":
            lat_min, lat_max = lat_far, lat_near
            if lat_max > 60.0:
                lat_max = 90.0
        else:
            lat_min, lat_max = lat_near, lat_far
            if lat_min < -60.0:
                lat_min = -90.0
        return lon_min, lat_min, lon_max, lat_max

    def nucleus(self, suid):
        x, y = self.d.cell(list(suid)).nucleus(plane=False)
        return float(x), float(y)


def cover_nodes(polygon, resolution: int, grid: Grid, margin: list):
    prepared = prep(polygon)
    boundary = seam_free_boundary(polygon)
    out = []

    def visit(suid, level):
        bb = _bbox_geom(*grid.bbox(suid))
        if not prepared.intersects(bb):
            return
        if level == resolution:
            pt = Point(*grid.nucleus(suid))
            margin.append(boundary.distance(pt))
            if polygon.contains(pt):
                out.append(suid)
            return
        if prepared.contains_properly(bb):      # interior only: consistent with the leaf rule
            out.append(suid)
            return
        for k in range(9):
            visit(suid + (k,), level + 1)

    for b in LETTERS:
        visit((b,), 0)
    return out


def seam_free_boundary(geom):
    """The polygon boundary without the segments lying on the antimeridian seam (x = ±180), which
    are artifacts of the GeoJSON split rather than edges of the region."""
    from shapely.geometry import LineString, MultiLineString
    segs = []
    polys = geom.geoms if hasattr(geom, "geoms") else [geom]
    for poly in polys:
        for ring in [poly.exterior, *poly.interiors]:
            cs = list(ring.coords)
            for a, b in zip(cs, cs[1:]):
                if abs(a[0]) == 180.0 and abs(b[0]) == 180.0:
                    continue
                segs.append(LineString([a, b]))
    return MultiLineString(segs)


def expand(nodes, resolution):
    leaves = set()
    for s in nodes:
        stack = [s]
        while stack:
            t = stack.pop()
            if len(t) - 1 == resolution:
                leaves.add(t)
            else:
                stack.extend(t + (k,) for k in range(9))
    return leaves


def canonical(leaves, resolution):
    cur = set(leaves)
    for _ in range(resolution):
        by_parent = {}
        for s in cur:
            if len(s) > 1:
                by_parent.setdefault(s[:-1], set()).add(s[-1])
        nxt = set()
        merged = set()
        for parent, kids in by_parent.items():
            if len(kids) == 9:
                nxt.add(parent)
                merged.add(parent)
        for s in cur:
            if not (len(s) > 1 and s[:-1] in merged):
                nxt.add(s)
        if nxt == cur:
            break
        cur = nxt
    return sorted(cur, key=lambda s: (len(s), suid_str(s)))


def brute(polygon, resolution, grid):
    out = set()
    for s in all_suids(resolution):
        if polygon.contains(Point(*grid.nucleus(s))):
            out.add(s)
    return out


def polyfill_cases(name: str):
    seam = 95.0 if name == "ogc" else 45.0   # lon_0 + 45: an equatorial face edge
    cases = [
        ("london", 10, box(-0.13, 51.50, -0.10, 51.52)),
        ("manhattan", 11, box(-74.02, 40.70, -73.93, 40.80)),
        ("central_park", 12, box(-73.981, 40.768, -73.949, 40.800)),
        ("london_hole", 10, Polygon(box(-0.13, 51.50, -0.10, 51.52).exterior.coords,
                                    [list(box(-0.122, 51.505, -0.108, 51.515).exterior.coords)])),
        ("face_edge", 8, box(seam - 5.0, 10.0, seam + 5.0, 20.0)),
        ("face_edge_r4", 4, box(seam - 5.0, 10.0, seam + 5.0, 20.0)),
        ("north_pole", 6, box(-181.0, 85.0, 181.0, 95.0)),
        ("north_pole_r3", 3, box(-181.0, 85.0, 181.0, 95.0)),
        ("south_cap", 6, box(-181.0, -95.0, 181.0, -87.0)),
        ("south_cap_r4", 4, box(-181.0, -95.0, 181.0, -87.0)),
        ("subcell_empty", 5, box(10.0001, 10.0001, 10.0002, 10.0002)),
        # edges carry small offsets so no nucleus lands exactly on a boundary (175.0 is one)
        ("antimeridian_multi", 6, MultiPolygon([box(175.0123, -5.0123, 180.0, 5.0321), box(-180.0, -5.0123, -174.9877, 5.0321)])),
        ("antimeridian_multi_r4", 4, MultiPolygon([box(175.0123, -5.0123, 180.0, 5.0321), box(-180.0, -5.0123, -174.9877, 5.0321)])),
    ]
    for i, (lat_lo, lat_hi) in enumerate(((-10.0, 10.0), (40.0, 60.0), (-60.0, -40.0), (70.0, 85.0), (-85.0, -70.0))):
        lat_lo, lat_hi = lat_lo + 0.0123, lat_hi - 0.0321
        cases.append((f"antimeridian_{i}", 6, MultiPolygon([box(170.0123, lat_lo, 180.0, lat_hi), box(-180.0, lat_lo, -169.9877, lat_hi)])))
        cases.append((f"antimeridian_{i}_r4", 4, MultiPolygon([box(170.0123, lat_lo, 180.0, lat_hi), box(-180.0, lat_lo, -169.9877, lat_hi)])))
    return cases


def gen_polyfill(name: str, p):
    grid = Grid(dggs_for(p))
    rows = []
    for cname, res, geom in polyfill_cases(name):
        margin = []
        nodes = cover_nodes(geom, res, grid, margin)
        leaves = expand(nodes, res)
        sfb = seam_free_boundary(geom)
        margin += [sfb.distance(Point(*grid.nucleus(s))) for s in leaves]
        if res <= 4:
            b = brute(geom, res, grid)
            assert b == leaves, f"{name}/{cname}: pruned cover != brute force ({len(b)} vs {len(leaves)})"
        m = min(margin) if margin else 1.0
        assert m >= 1e-9, f"{name}/{cname}: a nucleus lies {m} deg from the boundary; nudge the fixture"
        canon = canonical(leaves, res)
        rows.append({"name": cname, "resolution": res, "geometry": geom.__geo_interface__,
                     "canonical": [suid_str(s) for s in canon], "leaves": len(leaves),
                     "min_margin_deg": m})
        print(f"  {name}/{cname}: {len(leaves)} leaves, {len(canon)} canonical, margin {m:.3g}")
    dump(f"polyfill_{name}.json", {"profile": p, "cases": rows})


# ----------------------------------------------------------------------------- 6. profile ids
def gen_profile_ids():
    ids = {}
    for name, lon_0_udeg in (("ogc", 50_000_000), ("burin1", 0)):
        for h in ("sha256",):
            ids[f"{name}/{h}"] = profile_id(h, 9, 6, lon_0_udeg, 0, 0, 6_378_137_000_000, 298_257_223_563, 1, 0)
    ids["custom"] = profile_id("sha256", 9, 6, -123_456_789, 1, 2, 6_378_137_000_000, 298_257_223_563, 1_000_000, -5)
    dump("profile_ids.json", ids)


# ----------------------------------------------------------------------------- 7. zone topology
def gen_topology():
    """Edge neighbours (planar up/right/down/left) of every cell at levels 0-2 under all sixteen
    polar placements, and the sub-zones of sample parents sorted by their planar centres (top to
    bottom, then left to right): the reference for the scanline order."""
    neighbours = {}
    for ns in range(4):
        for ss in range(4):
            d = dggs_for(dict(lon_0=0.0, ns=ns, ss=ss))
            rows = []
            for res in range(3):
                for s in sorted(all_suids(res), key=suid_str):
                    nb = d.cell(list(s)).neighbors(plane=True)
                    rows.append([suid_str(s)] + [str(nb[k]) for k in ("up", "right", "down", "left")])
            neighbours[f"{ns}{ss}"] = rows
    scanline = {}
    for name, prof in PROFILES.items():
        d = dggs_for(prof)
        for parent, depth in (("N", 2), ("O", 2), ("Q", 3), ("S", 2), ("Q4", 2), ("N81", 2), ("S0", 2), ("R26", 1)):
            base = (parent[0], *(int(c) for c in parent[1:]))
            subs = [s for s in all_suids(len(base) - 1 + depth) if s[:len(base)] == base]
            key = lambda s: (-round(float(d.cell(list(s)).nucleus(plane=True)[1]), 6),
                             round(float(d.cell(list(s)).nucleus(plane=True)[0]), 6))
            scanline[f"{name}/{parent}/{depth}"] = [suid_str(s) for s in sorted(subs, key=key)]
    dump("topology.json", {"neighbours": neighbours, "scanline": scanline})


if __name__ == "__main__":
    which = set(sys.argv[1:]) or {"all"}
    if which & {"all", "geo"}:
        gen_auth_lat()
        gen_healpix()
        gen_rhealpix()
    if which & {"all", "cells"}:
        for name, p in PROFILES.items():
            gen_cells(name, p)
    if which & {"all", "polyfill"}:
        for name in ("ogc", "burin1"):
            gen_polyfill(name, PROFILES[name])
    if which & {"all", "ids"}:
        gen_profile_ids()
    if which & {"all", "topology"}:
        gen_topology()
