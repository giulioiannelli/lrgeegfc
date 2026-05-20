"""Graph-level surrogate utilities (matched-strength, etc.).

Distinct from `lrg_eegfc.utils.fc.msc.surrogates` which operates on raw
timeseries (Welch / phase randomization). This subpackage operates on
already-computed FC adjacency matrices.
"""
from .matched_strength import (
    SURROGATE_LRG_CACHE,
    load_or_compute_surrogate_eigs,
    strength_preserving_shuffle,
    surrogate_cache_path,
    verify_strengths,
)

__all__ = [
    "SURROGATE_LRG_CACHE",
    "load_or_compute_surrogate_eigs",
    "strength_preserving_shuffle",
    "surrogate_cache_path",
    "verify_strengths",
]
