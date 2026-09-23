"""Array functions over points and cell ids: numpy in, numpy out, spread over threads.

Cell ids here are the integer cids of one refinement level. At level ``r`` they fill the
half-open range ``[9**(r + 1), 15 * 9**r)``, so the whole grid at a level, a cell's ancestors and
its descendants are integer arithmetic; only the geometry goes to the compiled core.
"""
from __future__ import annotations

import numpy as np

from . import _burin

MAX_LEVEL = int(_burin.MAX_RESOLUTION)
_A = 9

__all__ = [
    "MAX_LEVEL",
    "level_range",
    "full_domain",
    "cell_levels",
    "as_cell_ids",
    "zoom_to",
    "cells_from_lonlat",
    "cells_to_lonlat",
    "cell_boundaries",
    "cell_neighbours",
]


def _check_level(level) -> int:
    if isinstance(level, bool) or int(level) != level or not 0 <= int(level) <= MAX_LEVEL:
        raise ValueError(f"level must be an integer between 0 and {MAX_LEVEL}, got {level!r}")
    return int(level)


def level_range(level: int) -> tuple[int, int]:
    """The half-open range of cids at ``level``.

    Parameters
    ----------
    level : int
        Refinement level, 0 to ``MAX_LEVEL``.

    Returns
    -------
    tuple of int
        ``(9**(level + 1), 15 * 9**level)``.
    """
    level = _check_level(level)
    return _A ** (level + 1), 15 * _A**level


def full_domain(level: int) -> np.ndarray:
    """Every cell at ``level``, sorted, as uint64 (``6 * 9**level`` of them)."""
    lo, hi = level_range(level)
    return np.arange(lo, hi, dtype=np.uint64)


def cell_levels(cell_ids) -> np.ndarray:
    """The level of each cid, as int64; raises ``ValueError`` if any id is not a cid."""
    ids = as_cell_ids(cell_ids)
    starts = np.array([_A ** (r + 1) for r in range(MAX_LEVEL + 1)], dtype=np.uint64)
    ends = np.array([15 * _A**r for r in range(MAX_LEVEL + 1)], dtype=np.uint64)
    levels = np.searchsorted(starts, ids, side="right").astype(np.int64) - 1
    bad = (levels < 0) | (ids >= ends[np.clip(levels, 0, MAX_LEVEL)])
    if bad.any():
        raise ValueError(f"{int(ids[bad].ravel()[0])} is not a cell id at levels 0-{MAX_LEVEL}")
    return levels


def as_cell_ids(cell_ids, level: int | None = None) -> np.ndarray:
    """Cell ids as a contiguous uint64 array of the same shape.

    Integer arrays of any width are accepted if non-negative. With ``level``, every id must be a
    cell at that level.

    Raises
    ------
    TypeError
        If the ids are not integers.
    ValueError
        If an id is negative, or not at ``level`` when one is given.
    """
    ids = np.asarray(cell_ids)
    if ids.dtype.kind not in "iu":
        raise TypeError(f"cell ids must be integers, got dtype {ids.dtype}")
    if ids.dtype.kind == "i" and (ids < 0).any():
        raise ValueError("cell ids are non-negative")
    ids = np.ascontiguousarray(ids, dtype=np.uint64)
    if level is not None:
        lo, hi = level_range(level)
        outside = (ids < lo) | (ids >= hi)
        if outside.any():
            raise ValueError(f"{int(ids[outside].ravel()[0])} is not a cell at level {level}")
    return ids


def zoom_to(cell_ids, level: int, new_level: int) -> np.ndarray:
    """Ancestors or descendants of cells at ``level``.

    Parameters
    ----------
    cell_ids : array_like of int
        Cells at ``level``.
    level, new_level : int
        The current level and the level to move to.

    Returns
    -------
    numpy.ndarray of uint64
        For ``new_level <= level``, the ancestors, same shape as ``cell_ids``. For
        ``new_level > level``, shape ``cell_ids.shape + (9**(new_level - level),)``: each cell's
        descendants in nested (cid) order.
    """
    ids = as_cell_ids(cell_ids, level)
    new_level = _check_level(new_level)
    step = abs(new_level - level)
    if new_level <= level:
        return ids // np.uint64(_A**step)
    children = np.arange(_A**step, dtype=np.uint64)
    return ids[..., None] * np.uint64(_A**step) + children


def cells_from_lonlat(lon, lat, level: int, *, profile=None, nthreads: int = 0) -> np.ndarray:
    """The cell at ``level`` containing each point.

    Parameters
    ----------
    lon, lat : array_like of float
        Degrees; broadcast against each other.
    level : int
        Refinement level, 0 to ``MAX_LEVEL``.
    profile : burin.Profile, optional
        The grid; the OGC rHEALPix profile by default.
    nthreads : int, default 0
        Worker threads; 0 uses the shared pool, 1 stays on the calling thread.

    Returns
    -------
    numpy.ndarray of uint64
        The broadcast shape of ``lon`` and ``lat``.

    Raises
    ------
    ValueError
        If a coordinate is not finite.
    """
    lon, lat = np.broadcast_arrays(np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64))
    cells = _burin._cells_from_lonlat(np.ascontiguousarray(lon).ravel(), np.ascontiguousarray(lat).ravel(),
                                      _check_level(level), profile, nthreads)
    return cells.reshape(lon.shape)


def cells_to_lonlat(cell_ids, *, profile=None, nthreads: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """The nucleus of each cell as ``(lon, lat)`` degrees, each of the shape of ``cell_ids``."""
    ids = as_cell_ids(cell_ids)
    lon, lat = _burin._cells_to_lonlat(ids.ravel(), profile, nthreads)
    return lon.reshape(ids.shape), lat.reshape(ids.shape)


def cell_boundaries(cell_ids, *, n: int = 5, profile=None, nthreads: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Each cell's polygon, as ragged arrays.

    Rings are closed and counterclockwise in ``(lon, lat)`` degrees. Equatorial cells are their
    four corners (their edges are meridians and parallels); polar cells have ``n`` points per
    edge; a ring crossing the antimeridian has its western longitudes moved past 180; the cap
    around each pole is closed through the pole along ±180.

    Returns
    -------
    coords : numpy.ndarray of float64, shape (M, 2)
        Every ring in turn, for the flattened ``cell_ids``.
    offsets : numpy.ndarray of int64, shape (N + 1,)
        Cell ``i`` is ``coords[offsets[i]:offsets[i + 1]]``.
    """
    if n < 2:
        raise ValueError("n must be at least 2")
    return _burin._cell_boundaries(as_cell_ids(cell_ids).ravel(), int(n), profile, nthreads)


def cell_neighbours(cell_ids, *, profile=None, nthreads: int = 0) -> np.ndarray:
    """The four edge neighbours of each cell, shape ``cell_ids.shape + (4,)``: up, right, down,
    left in the cell's planar frame."""
    ids = as_cell_ids(cell_ids)
    return _burin._cell_neighbours(ids.ravel(), profile, nthreads).reshape(ids.shape + (4,))
