"""Was this photo taken in Central Park, while the park was open?

A publisher fingerprints a place and its opening hours, 32 bytes each. For one photo, the
photographer gets two small proofs from the publisher. A verifier who holds only the two
fingerprints checks them, without ever seeing the park's outline or its hours.

    python examples/central_park.py
"""
import json
from datetime import date, datetime, timedelta, timezone

import numpy as np

import burin
from burin import time as bt

LEVEL = 11  # cells about 50 m across


def instant(text):
    """An ISO 8601 time with its UTC offset, as datetime64[us] in UTC."""
    utc = datetime.fromisoformat(text).astimezone(timezone.utc).replace(tzinfo=None)
    return np.datetime64(utc, "us")


# The publisher ----------------------------------------------------------------------------------
# Central Park's corners: 59th St at Central Park West, at Fifth Avenue; 110th St at Fifth Avenue,
# at Central Park West.
outline = {"type": "Polygon", "coordinates": [[[-73.9819, 40.7681], [-73.9730, 40.7644],
                                               [-73.9493, 40.7967], [-73.9580, 40.8003],
                                               [-73.9819, 40.7681]]]}
park = burin.Tree.from_geojson(outline, LEVEL)

# Open 6 am to 1 am, New York time (UTC-4 all July), every day of July 2026.
days = [date(2026, 7, 1) + timedelta(days=i) for i in range(31)]
hours = bt.tree([(instant(f"{d}T06:00-04:00"), instant(f"{d + timedelta(days=1)}T01:00-04:00")) for d in days])

published = {"place": park.root_hex, "hours": hours.root_hex}
print("PUBLISHER")
print(f"  the park is {len(park.leaves()):,} cells of {burin.cell_area_m2(LEVEL):,.0f} m², "
      f"held as {len(park.cells())} canonical cells")
print(f"  its hours are {len(days)} daily windows, held as {len(hours.cells())} time cells")
print(f"  publishes:  place {published['place']}")
print(f"              hours {published['hours']}")


# The photographer -------------------------------------------------------------------------------
photos = [
    ("Bethesda Terrace, 2:32 pm", -73.9712, 40.7740, "2026-07-18T14:32:07-04:00"),
    ("Bethesda Terrace, 2:10 am", -73.9712, 40.7740, "2026-07-18T02:10:00-04:00"),
    ("Times Square, 2:32 pm", -73.9855, 40.7580, "2026-07-18T14:32:07-04:00"),
]


def where_and_when(lon, lat, taken):
    """The photo's cell and its tick: what everyone derives from its metadata."""
    return int(burin.cells_from_lonlat(lon, lat, LEVEL)), int(bt.cells_from_times(instant(taken)))


def ask_publisher(cell, moment):
    """The publisher opens both trees at the photo's cell and tick; a proof of absence is a proof too."""
    return park.open(cell), hours.open(moment)


# The verifier -----------------------------------------------------------------------------------
def check(proof, fingerprint, cell):
    """True or False if the proof settles `cell` against `fingerprint`; None if it proves nothing."""
    if not burin.verify_opening(proof) or proof["root"] != fingerprint:
        return None
    stated, covered = burin.opening_statement(proof)
    return covered if stated == cell else None


def verdict(result, yes, no):
    return {True: yes, False: no, None: "REJECTED: the proof does not hold"}[result]


print("\nPHOTOGRAPHER and VERIFIER")
proofs = {}
for name, lon, lat, taken in photos:
    cell, moment = where_and_when(lon, lat, taken)
    place_proof, hours_proof = proofs[name] = ask_publisher(cell, moment)
    size = len(json.dumps(place_proof)) + len(json.dumps(hours_proof))
    in_park = check(place_proof, published["place"], cell)
    is_open = check(hours_proof, published["hours"], moment)
    print(f"  {name}  ({size:,} bytes of proof)")
    print(f"    place: {verdict(in_park, 'in the park', 'not in the park')}")
    print(f"    time:  {verdict(is_open, 'during opening hours', 'outside opening hours')}")
    print(f"    →      {'yes' if in_park and is_open else 'no'}")


print("\nTRYING TO CHEAT")
name, lon, lat, taken = photos[0]
cell, moment = where_and_when(lon, lat, taken)
place_proof, hours_proof = proofs[name]

forged = json.loads(json.dumps(place_proof))
entries = forged["opening"]["entries"]
i = next(i for i, e in enumerate(entries) if isinstance(e, str))  # a sibling's hash
entries[i] = entries[i][:-1] + ("1" if entries[i][-1] == "0" else "0")
print(f"  one hex digit of the place proof changed:     {verdict(check(forged, published['place'], cell), 'in', 'out')}")

times_square = where_and_when(*photos[2][1:])
print(f"  Bethesda's proof shown for the Times Square photo: "
      f"{verdict(check(place_proof, published['place'], times_square[0]), 'in', 'out')}")

other_park = burin.Tree.from_geojson(outline, LEVEL - 1)
print(f"  a proof checked against another fingerprint:  "
      f"{verdict(check(place_proof, other_park.root_hex, cell), 'in', 'out')}")
