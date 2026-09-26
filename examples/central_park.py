"""Was this photo taken in Central Park, while the park was open?

A publisher publishes two fingerprints: the park, and its opening hours. For one photo, the
photographer gets a proof for its place and one for its time. A verifier holding only the two
fingerprints asks each one its question and gets yes or no, or an InvalidProof if the proof
does not answer it.

This shows the photo's claimed place and time are consistent with what the park published. It
does not show the photo was taken there: the metadata is whatever the photographer supplied.

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
    return np.datetime64(datetime.fromisoformat(text).astimezone(timezone.utc).replace(tzinfo=None), "us")


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

published = [str(park.fingerprint), str(hours.fingerprint)]
print("PUBLISHER")
print(f"  the park is {len(park.leaves()):,} cells of {burin.cell_area_m2(LEVEL):,.0f} m²; "
      f"its hours are {len(days)} daily windows")
for line in published:
    print(f"  publishes {line}")


# The photographer, asking the publisher ----------------------------------------------------------
photos = [
    ("Bethesda Terrace, 2:32 pm", (-73.9712, 40.7740), instant("2026-07-18T14:32:07-04:00")),
    ("Bethesda Terrace, 2:10 am", (-73.9712, 40.7740), instant("2026-07-18T02:10:00-04:00")),
    ("Times Square, 2:32 pm", (-73.9855, 40.7580), instant("2026-07-18T14:32:07-04:00")),
]
proofs = {name: (park.prove_point(*where), hours.prove_instant(when)) for name, where, when in photos}


# The verifier: holds two strings, asks two questions --------------------------------------------
place = burin.Fingerprint.parse(published[0])
open_ = burin.Fingerprint.parse(published[1])

print("\nVERIFIER")
for name, where, when in photos:
    place_proof, time_proof = proofs[name]
    in_park = place.check_point(place_proof, *where)
    is_open = open_.check_instant(time_proof, when)
    size = len(json.dumps(place_proof.to_json())) + len(json.dumps(time_proof.to_json()))
    print(f"  {name:28} in the park: {'yes' if in_park else 'no ':3}   open: {'yes' if is_open else 'no ':3}"
          f"   so: {'yes' if in_park and is_open else 'no'}   ({size:,} bytes of proof)")


print("\nTRYING TO CHEAT")
name, where, when = photos[0]
place_proof, time_proof = proofs[name]


def attempt(label, check):
    try:
        print(f"  {label:44} answered {check()}  (this should never happen)")
    except burin.InvalidProof as refused:
        print(f"  {label:44} refused: {refused}")


forged = place_proof.to_json()
entries = forged["opening"]["entries"]
i = next(i for i, e in enumerate(entries) if isinstance(e, str))  # a sibling's hash
entries[i] = entries[i][:-1] + ("1" if entries[i][-1] == "0" else "0")
attempt("one hex digit of the place proof changed", lambda: place.check_point(forged, *where))
attempt("Bethesda's proof for the Times Square photo", lambda: place.check_point(place_proof, *photos[2][1]))
attempt("the place proof offered as the time proof", lambda: open_.check_instant(place_proof, when))
other = burin.Tree.from_geojson(outline, LEVEL - 1).fingerprint
attempt("checked against a coarser park's fingerprint", lambda: other.check_point(place_proof, *where))
