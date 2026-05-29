"""Smoke tests for lrg_eegfc.utils.metrics.hypothesis."""
from __future__ import annotations

import numpy as np
import pytest

from lrg_eegfc.utils.metrics.hypothesis import (
    bh_fdr,
    boot_ci_mean,
    cluster_stats,
    loo_sensitivity,
    rank_biserial,
    surrogate_p_value,
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


def test_surrogate_p_value_upper_unbiased():
    surr = np.arange(100, dtype=float)  # 0..99
    p = surrogate_p_value(obs=95.0, surr_array=surr, tail="upper")
    # 5 surrogate values (95..99) >= 95 → (5 + 1) / (100 + 1) = 6/101
    assert p == pytest.approx(6 / 101, abs=1e-9)


def test_surrogate_p_value_lower():
    surr = np.arange(100, dtype=float)
    p = surrogate_p_value(obs=4.0, surr_array=surr, tail="lower")
    # 5 surrogate values (0..4) <= 4 → (5 + 1) / (100 + 1)
    assert p == pytest.approx(6 / 101, abs=1e-9)


def test_surrogate_p_value_naive_no_addone():
    surr = np.arange(100, dtype=float)
    p = surrogate_p_value(obs=95.0, surr_array=surr, tail="upper", add_one=False)
    assert p == pytest.approx(5 / 100, abs=1e-9)


def test_surrogate_p_value_empty_returns_nan():
    p = surrogate_p_value(obs=1.0, surr_array=np.array([]))
    assert np.isnan(p)


def test_surrogate_p_value_nans_ignored():
    surr = np.array([0.0, np.nan, 5.0, np.nan, 10.0])
    p = surrogate_p_value(obs=4.0, surr_array=surr, tail="upper")
    # finite surr: [0,5,10], 2 ≥ 4 → (2 + 1) / (3 + 1) = 3/4
    assert p == pytest.approx(3 / 4, abs=1e-9)


def test_surrogate_p_value_bad_tail_raises():
    with pytest.raises(ValueError, match="tail"):
        surrogate_p_value(obs=1.0, surr_array=np.array([0.0, 1.0]), tail="both")


def test_loo_sensitivity_basic_shape_and_full_value():
    rng = np.random.default_rng(0)
    vals = rng.normal(size=10)
    res = loo_sensitivity(vals, test_fn=lambda a: float(a.mean()))
    assert res["full"] == pytest.approx(vals.mean(), abs=1e-12)
    assert res["loo"].shape == (10,)
    # LOO of mean: mean of n-1 elements ≠ overall mean
    assert not np.allclose(res["loo"], res["full"])
    # worst is the max LOO mean → drops the smallest value
    expected_worst_idx = int(np.argmin(vals))
    assert res["worst_patient"] == str(expected_worst_idx)


def test_loo_sensitivity_with_labels():
    vals = np.array([10.0, 1.0, 1.0, 1.0, 1.0])
    labels = ["A", "B", "C", "D", "E"]
    res = loo_sensitivity(vals, test_fn=lambda a: float(a.mean()), labels=labels)
    # Dropping "A" (10) leaves mean=1.0; dropping any other leaves mean=3.25.
    # "Worst" by np.nanargmax is the highest LOO → keeping the 10 → dropping a non-A.
    assert res["worst_patient"] in {"B", "C", "D", "E"}
    assert res["worst"] == pytest.approx(3.25, abs=1e-9)


def test_loo_sensitivity_labels_length_mismatch_raises():
    vals = np.array([1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="labels length"):
        loo_sensitivity(vals, test_fn=lambda a: 0.0, labels=["A", "B"])


def test_loo_sensitivity_handles_test_fn_failure():
    vals = np.array([1.0, 2.0, 3.0])

    def flaky(a):
        if a.shape[0] == 2 and a[0] == 2.0:  # fails on a specific subset
            raise RuntimeError("boom")
        return float(a.mean())

    res = loo_sensitivity(vals, test_fn=flaky)
    assert np.isnan(res["loo"][0])  # dropping idx 0 leaves [2,3], triggers flaky
    assert np.isfinite(res["loo"][1])
    assert np.isfinite(res["loo"][2])
