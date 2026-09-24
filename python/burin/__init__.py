"""Deterministic, equal-area coverage commitments over the rHEALPix DGGS.

The compiled core is ``burin._burin``; everything public is re-exported here.
"""
from ._burin import *  # noqa: F401,F403
from ._burin import __version__  # noqa: F401
from ._vector import *  # noqa: F401,F403
from . import time  # noqa: F401,E402  (burin.time: the same program for time)
