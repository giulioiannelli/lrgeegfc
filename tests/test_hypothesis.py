"""Smoke tests for lrg_eegfc.utils.metrics.hypothesis."""
from __future__ import annotations

import numpy as np
import pytest

from lrg_eegfc.utils.metrics.hypothesis import (
    bh_fdr,
    boot_ci_mean,
    cluster_stats,
    rank_biserial,
    wilcoxon_z,
)


def test_wilcoxon_z_all_positive():
    z, p = wilcoxon_z(np.array([1.0, 2.0, 3.0, 4.0, 5.0]))
    assert np.isfinite(z) and z > 0
    assert p < 0.05


def test_wilcoxon_z_short_returns_nan():
    z, p = wilcoxon_z(np.array([1.0, 2.0]))
    assert np.isnan(z) and np.isnan(p)


def test_rank_biserial_bounds_and_sign():
    pos = rank_biserial(np.array([1.0, 2.0, 3.0, 4.0, 5.0]))
    neg = rank_biserial(np.array([-1.0, -2.0, -3.0, -4.0, -5.0]))
    mix = rank_biserial(np.array([1.0, -1.0, 2.0, -2.0, 3.0]))
    assert pos == pytest.approx(1.0, abs=1e-9)
    assert neg == pytest.approx(-1.0, abs=1e-9)
    assert -1.0 <= mix <= 1.0


def test_boot_ci_mean_brackets_mean():
    rng = np.random.default_rng(42)
    x = rng.normal(loc=1.5, scale=0.5, size=200)
    mean, lo, hi = boot_ci_mean(x, B=2000, rng=rng)
    assert lo < mean < hi
    assert abs(mean - 1.5) < 0.1


def test_bh_fdr_monotone_and_clipped():
    qs = bh_fdr([0.001, 0.01, 0.02, 0.05, 0.1, 0.5])
    assert all(0.0 <= q <= 1.0 for q in qs)
    assert qs[0] <= qs[-1]  # ordering preserved after BH


def test_bh_fdr_preserves_input_order():
    pvals = [0.5, 0.001, 0.02, 0.1, 0.01, 0.05]
    qs = bh_fdr(pvals)
    assert len(qs) == len(pvals)
    assert qs[1] == min(qs)  # smallest p was at index 1


def test_cluster_stats_detects_run():
    z = np.array([0.0, 0.0, 2.0, 3.0, 2.5, 0.0, 1.5, 0.0])
    out = cluster_stats(z, thresh=1.0)
    assert len(out) == 2
    (s0, e0, m0), (s1, e1, m1) = out
    assert (s0, e0) == (2, 4)
    assert m0 == pytest.approx(7.5)
    assert (s1, e1) == (6, 6)
    assert m1 == pytest.approx(1.5)


def test_cluster_stats_no_supra():
    out = cluster_stats(np.array([0.0, 0.5, 0.9]), thresh=1.0)
    assert out == []
