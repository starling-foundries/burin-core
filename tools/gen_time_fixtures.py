"""Conformance fixture for time (SPEC §9–§12), computed from the specification's text alone.

An independent reference: ticks by exact integer floor division, an interval's cells by
recursive halving (not the boundary walk), roots from the hash definitions of §4, Allen's
relations from the endpoint conditions of §12. Standard library only; integers only.

    python tools/gen_time_fixtures.py     # writes crates/burin-core/tests/fixtures/time.json
"""
import hashlib
import json
from functools import lru_cache
from pathlib import Path

DEPTH, TICKS, EPOCH_TICK = 61, 1 << 61, 1 << 60
OUT = Path(__file__).resolve().parents[1] / "crates/burin-core/tests/fixtures/time.json"


def tick(t, tick_us, epoch_us):
    k = (t - epoch_us) // tick_us + EPOCH_TICK  # Python's // floors toward -inf
    return k if 0 <= k < TICKS else None


def tick_start(k, tick_us, epoch_us):
    t = (k - EPOCH_TICK) * tick_us + epoch_us
    return t if -(1 << 63) <= t < (1 << 63) else None


def merge(intervals):
    out = []
    for a, b in sorted(intervals):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def clip(intervals, lo, hi):
    parts = [(max(a, lo), min(b, hi)) for a, b in intervals if a < hi and b > lo]
    return parts, sum(b - a for a, b in merge(parts))


def cells_of(intervals, level=0, k=0):
    """The compact cells of a union of intervals: a cell is kept whole if covered, else split."""
    s = DEPTH - level
    parts, covered = clip(intervals, k << s, (k + 1) << s)
    if covered == 0:
        return []
    if covered == 1 << s:
        return [(2 << level) + k]
    return cells_of(parts, level + 1, 2 * k) + cells_of(parts, level + 1, 2 * k + 1)


def sha(tag, *parts):
    return hashlib.sha256(bytes([tag]) + b"".join(parts)).digest()


@lru_cache(maxsize=None)
def ladder(full, d):
    return sha(1 if full else 0) if d == 0 else sha(2, ladder(full, d - 1), ladder(full, d - 1))


def node_hash(intervals, level=0, k=0):
    s = DEPTH - level
    parts, covered = clip(intervals, k << s, (k + 1) << s)
    if covered in (0, 1 << s):
        return ladder(covered != 0, s)
    return sha(2, node_hash(parts, level + 1, 2 * k), node_hash(parts, level + 1, 2 * k + 1))


def root(intervals):
    return sha(3, node_hash(intervals)).hex()


def allen(a, b):
    (a1, a2), (b1, b2) = a, b
    table = [
        ("before", a2 < b1), ("meets", a2 == b1), ("overlaps", a1 < b1 < a2 < b2),
        ("starts", a1 == b1 and a2 < b2), ("during", b1 < a1 and a2 < b2),
        ("finishes", b1 < a1 and a2 == b2), ("equals", a1 == b1 and a2 == b2),
        ("finished-by", a1 < b1 and a2 == b2), ("contains", a1 < b1 and b2 < a2),
        ("started-by", a1 == b1 and b2 < a2), ("overlapped-by", b1 < a1 < b2 < a2),
        ("met-by", b2 == a1), ("after", b2 < a1),
    ]
    names = [n for n, holds in table if holds]
    assert len(names) == 1, (a, b, names)
    return names[0]


def main():
    e = 1 << 60
    profiles = {
        "microseconds": (1, 0),
        "days": (86_400_000_000, 0),
        "millis_from_year_1": (1_000, -62_135_596_800_000_000),
    }
    instants = [0, 1, -1, -999, 999, 1_000, -1_001, 1_373_846_400_000_000, -2_208_988_800_000_000,
                e - 1, e, -e, -e - 1, (1 << 63) - 1, -(1 << 63)]
    ticks = [0, 1, EPOCH_TICK - 1, EPOCH_TICK, EPOCH_TICK + 1, TICKS - 1, TICKS]
    lines = {
        name: {
            "tick_us": tick_us, "epoch_us": epoch_us,
            "ticks": [[t, tick(t, tick_us, epoch_us)] for t in instants],
            "starts": [[k, tick_start(k, tick_us, epoch_us) if k <= TICKS else None] for k in ticks + [TICKS + 1]],
        }
        for name, (tick_us, epoch_us) in profiles.items()
    }
    sets = [
        [(0, 1)], [(0, TICKS)], [(TICKS - 1, TICKS)], [(5, 9)], [(1, TICKS - 1)],
        [(EPOCH_TICK, EPOCH_TICK + 86_400_000_000)],
        [(EPOCH_TICK + 1_373_846_400_000_000, EPOCH_TICK + 1_373_932_800_000_000)],
        [(3, 17), (40, 41), (17, 30)],
        [(12_345_678_901_234, 98_765_432_109_876), (2**59, 2**59 + 7), (2**60 - 3, 2**60 + 3)],
    ]
    small = [(a, b) for a in range(6) for b in range(a + 1, 7)]
    fx = {
        "empty_root": sha(3, ladder(False, DEPTH)).hex(),
        "full_root": sha(3, ladder(True, DEPTH)).hex(),
        "lines": lines,
        "intervals": [{"intervals": [list(i) for i in s], "cells": cells_of(s), "root": root(s)} for s in sets],
        "allen": [[list(a), list(b), allen(a, b)] for a in small for b in small],
    }
    OUT.write_text(json.dumps(fx, indent=1) + "\n")
    print(f"wrote {OUT.name}: {len(sets)} sets, {len(fx['allen'])} Allen pairs")


if __name__ == "__main__":
    main()
