"""DEPRECATED — use ``lrg_eegfc.utils.metrics.hypothesis`` instead.

Thin re-export shim kept for backcompat. Will be removed after the
failed-session archive in phase 3.7 of the reorg. New code must import
from the library directly.
"""
from __future__ import annotations

from lrg_eegfc.utils.metrics.hypothesis import (
    bh_fdr,
    boot_ci_mean,
    cluster_stats,
    rank_biserial,
    wilcoxon_z,
)

__all__ = ["wilcoxon_z", "rank_biserial", "boot_ci_mean", "bh_fdr", "cluster_stats"]
