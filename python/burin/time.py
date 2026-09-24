"""Time (SPEC §9–§12), shaped like the spatial functions: numpy in, numpy out.

Instants are datetime64 of any unit, ISO 8601 strings, or integers taken as POSIX microseconds; a
float is never an instant. Under a profile, each falls in a tick (1 µs by default), and ticks are
the leaves of the dyadic time hierarchy: a level-``r`` cell spans ``2**(61 - r)`` ticks and its
cid is ``2**(r + 1) + k``. Everything is integer arithmetic, exact on every platform.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Any, Iterable

import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import _burin

DEPTH = 61
TICKS = 1 << DEPTH
EPOCH_TICK = 1 << (DEPTH - 1)
#: The OGC-registered temporal CRS the line counts in.
UNIX_TIME_CRS = "https://www.opengis.net/def/crs/OGC/0/UnixTime"
#: The most cells ``zoom_to`` returns in one call (2 GiB of ids).
MAX_ZOOM_CELLS = 2**28

__all__ = [
    "DEPTH",
    "TICKS",
    "EPOCH_TICK",
    "UNIX_TIME_CRS",
    "MAX_ZOOM_CELLS",
    "posix_us",
    "ticks",
    "level_range",
    "cell_levels",
    "cells_from_times",
    "times_of_cells",
    "zoom_to",
    "interval",
    "interval_cells",
    "coarse_cells",
    "tree",
    "allen",
]

_US_PER: dict[str, int | Fraction] = {"W": 604_800_000_000, "D": 86_400_000_000, "h": 3_600_000_000, "m": 60_000_000,
           "s": 1_000_000, "ms": 1_000, "us": 1, "ns": Fraction(1, 10**3), "ps": Fraction(1, 10**6),
           "fs": Fraction(1, 10**9), "as": Fraction(1, 10**12)}
_I64_MAX = 2**63 - 1


def posix_us(times: ArrayLike) -> NDArray[np.int64]:
    """POSIX microseconds of each instant, same shape.

    datetime64 of any unit is floored to the microsecond; years and months go through days, which
    is exact. Strings are read by numpy's ISO 8601 parser. Integers are microseconds already.

    Raises
    ------
    TypeError
        For floats or other non-time values.
    ValueError
        For NaT, or an instant outside the 64-bit microsecond range.
    """
    a = np.asarray(times)
    if a.dtype.kind in "iu":
        if a.dtype.kind == "u" and (a > _I64_MAX).any():
            raise ValueError("an instant is outside the 64-bit microsecond range")
        return a.astype(np.int64)
    if a.dtype.kind == "U":
        a = a.astype("datetime64")
    if a.dtype.kind != "M":
        raise TypeError(f"instants are datetime64, ISO 8601 strings or integer microseconds, got dtype {a.dtype}; "
                        "a float is never an instant")
    if np.isnat(a).any():
        raise ValueError("NaT is not an instant")
    unit, count = np.datetime_data(a.dtype)
    if unit in ("Y", "M"):
        a, unit, count = a.astype("datetime64[D]"), "D", 1
    ratio = Fraction(_US_PER[unit]) * count
    v = a.view(np.int64)
    if ratio.denominator == 1:
        if (np.abs(v) > _I64_MAX // ratio.numerator).any():
            raise ValueError("an instant is outside the 64-bit microsecond range")
        return v * np.int64(ratio.numerator)
    if ratio.numerator == 1:
        return np.floor_divide(v, np.int64(ratio.denominator))
    flat = [(int(x) * ratio.numerator) // ratio.denominator for x in v.ravel()]
    return np.array(flat, dtype=np.int64).reshape(v.shape)


def ticks(times: ArrayLike, *, profile: _burin.Profile | None = None) -> NDArray[np.uint64]:
    """The tick containing each instant, same shape; ``ValueError`` off the line."""
    us = posix_us(times)
    return _burin._ticks(np.ascontiguousarray(us).ravel(), profile).reshape(us.shape)


def _check_level(level: Any) -> int:
    if isinstance(level, bool) or int(level) != level or not 0 <= int(level) <= DEPTH:
        raise ValueError(f"time level must be an integer between 0 and {DEPTH}, got {level!r}")
    return int(level)


def level_range(level: int) -> tuple[int, int]:
    """The half-open range of time cids at ``level``: ``(2**(level + 1), 2**(level + 2))``."""
    level = _check_level(level)
    return 2 ** (level + 1), 2 ** (level + 2)


def _as_cids(cids: ArrayLike) -> NDArray[np.uint64]:
    ids = np.asarray(cids)
    if ids.dtype.kind not in "iu":
        raise TypeError(f"cell ids must be integers, got dtype {ids.dtype}")
    if ids.dtype.kind == "i" and (ids < 0).any():
        raise ValueError("cell ids are non-negative")
    return np.asarray(ids, dtype=np.uint64, order="C")


def cell_levels(cids: ArrayLike) -> NDArray[np.int64]:
    """The level of each time cid; ``ValueError`` if any is not one."""
    ids = _as_cids(cids)
    levels = np.array([int(c).bit_length() - 2 for c in ids.ravel()], dtype=np.int64).reshape(ids.shape)
    bad = (ids < 2) | (levels > DEPTH)
    if bad.any():
        raise ValueError(f"{int(ids[bad].ravel()[0])} is not a time cell at levels 0-{DEPTH}")
    return levels


def cells_from_times(times: ArrayLike, level: int = DEPTH, *, profile: _burin.Profile | None = None) -> NDArray[np.uint64]:
    """The time cell at ``level`` containing each instant; at level 61 it is the instant's tick."""
    shift = np.uint64(DEPTH - _check_level(level))
    return (np.uint64(1 << (DEPTH + 1)) + ticks(times, profile=profile)) >> shift


def times_of_cells(cids: ArrayLike, *, profile: _burin.Profile | None = None) -> tuple[NDArray[np.datetime64], NDArray[np.datetime64]]:
    """The ``[start, end)`` of each time cell as ``datetime64[us]``, each of the shape of ``cids``."""
    ids = _as_cids(cids)
    lo, hi = _burin._cell_ticks(ids.ravel())
    start = _burin._tick_starts(lo, profile).view("datetime64[us]")
    end = _burin._tick_starts(hi, profile).view("datetime64[us]")
    return start.reshape(ids.shape), end.reshape(ids.shape)


def zoom_to(cids: ArrayLike, level: int, new_level: int) -> NDArray[np.uint64]:
    """Ancestors (same shape) or descendants (shape ``cids.shape + (2**(new_level - level),)``,
    in time order) of time cells at ``level``."""
    ids = _as_cids(cids)
    level, new_level = _check_level(level), _check_level(new_level)
    lo, hi = level_range(level)
    if ((ids < lo) | (ids >= hi)).any():
        raise ValueError(f"not every id is a time cell at level {level}")
    step = abs(new_level - level)
    if new_level <= level:
        return ids >> np.uint64(step)
    if ids.size * 2**step > MAX_ZOOM_CELLS:
        raise ValueError(f"{ids.size * 2**step} descendants is more than MAX_ZOOM_CELLS ({MAX_ZOOM_CELLS})")
    return (ids[..., None] << np.uint64(step)) + np.arange(2**step, dtype=np.uint64)


def interval(start: Any, end: Any, *, profile: _burin.Profile | None = None) -> tuple[int, int]:
    """The tick interval ``[lo, hi)`` of the instants ``[start, end)``: every tick they touch."""
    s, e = (int(posix_us(x)) for x in (start, end))
    if e <= s:
        raise ValueError(f"[{start}, {end}) is empty")
    lo, last = (int(t) for t in ticks(np.array([s, e - 1], dtype=np.int64), profile=profile))
    return lo, last + 1


def interval_cells(lo: int, hi: int) -> NDArray[np.uint64]:
    """The compact cells of the tick interval ``[lo, hi)``, in time order (at most 122)."""
    return np.array(_burin._interval_cells(int(lo), int(hi)), dtype=np.uint64)


def coarse_cells(lo: int, hi: int, level: int) -> NDArray[np.uint64]:
    """The level-``level`` cells whose midpoint tick lies in ``[lo, hi)``, in time order."""
    return np.array(_burin._coarse_cells(int(lo), int(hi), _check_level(level)), dtype=np.uint64)


def tree(spans: Iterable[tuple[Any, Any]], *, profile: _burin.Profile | None = None) -> _burin.Tree:
    """The canonical time tree of a union of ``(start, end)`` instant intervals (see ``interval``)."""
    return _burin.Tree.from_intervals([interval(s, e, profile=profile) for s, e in spans], profile)


def allen(a: tuple[int, int], b: tuple[int, int]) -> str:
    """Allen's relation of the tick interval ``a`` to ``b``: ``"before"``, ``"meets"``, … ``"after"``."""
    return _burin._allen((int(a[0]), int(a[1])), (int(b[0]), int(b[1])))
