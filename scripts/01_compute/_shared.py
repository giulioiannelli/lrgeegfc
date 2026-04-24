"""Shared statistics helpers for the H2 family of cross-patient tests.

Canonical implementations of the four helpers that every cross-patient
hypothesis test in `scripts/01_compute/` needs. Extracted 2026-04-24 from
three scripts that each had a local copy:

  * h2c_ultrametric_drift.py
  * h2d_coactivation_persistence.py
  * rigorous_hypothesis_test.py

The versions below match rigorous_hypothesis_test.py's signatures
(cleanest defaults, includes the `rng` parameter on the bootstrap for
reproducibility). If a caller wants different defaults it can override
them at the call site.

See `.agents/guides/02_methods/H2_METRICS.md` for how these feed into
the overall hypothesis-testing framework.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from scipy import stats


__all__ = ["wilcoxon_z", "rank_biserial", "boot_ci_mean", "bh_fdr"]


N_BOOT_DEFAULT = 10_000


def wilcoxon_z(x: np.ndarray) -> tuple[float, float]:
    """One-sample Wilcoxon signed-rank test against 0, one-sided 'greater'.

    Returns ``(z, p)``. The z-score is computed from the normal
    approximation (no ties correction — adequate for continuous
    contrasts where exact zeros are vanishingly rare). Exact zeros are
    dropped; n<3 returns ``(nan, nan)``.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    x = x[x != 0.0]
    n = len(x)
    if n < 3:
        return np.nan, np.nan
    res = stats.wilcoxon(x, alternative="greater", method="approx")
    mu = n * (n + 1) / 4.0
    sigma = np.sqrt(n * (n + 1) * (2 * n + 1) / 24.0)
    z = (float(res.statistic) - mu) / sigma if sigma > 0 else np.nan
    return float(z), float(res.pvalue)


def rank_biserial(x: np.ndarray) -> float:
    """Rank-biserial correlation for one-sample signed-rank.

    Unbiased effect-size estimator, range ``[-1, +1]``. ``+1`` means every
    sample is positive and larger in absolute value than any negative
    sample (if any). Drops exact zeros; n<3 returns ``nan``.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    x = x[x != 0.0]
    n = len(x)
    if n < 3:
        return np.nan
    ranks = stats.rankdata(np.abs(x))
    r_plus = ranks[x > 0].sum()
    r_minus = ranks[x < 0].sum()
    return float((r_plus - r_minus) / (n * (n + 1) / 2.0))


def boot_ci_mean(
    x: np.ndarray,
    B: int = N_BOOT_DEFAULT,
    rng: Optional[np.random.Generator] = None,
) -> tuple[float, float, float]:
    """Return ``(mean, lo, hi)`` 95 % bootstrap CI of the sample mean.

    B = number of bootstrap resamples (default 10 000). ``rng`` is an
    optional ``numpy.random.Generator`` for reproducibility; if omitted
    uses ``default_rng(0)``.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return float(np.nan), float(np.nan), float(np.nan)
    rng = rng or np.random.default_rng(0)
    idx = rng.integers(0, len(x), size=(B, len(x)))
    means = x[idx].mean(axis=1)
    lo, hi = np.quantile(means, [0.025, 0.975])
    return float(x.mean()), float(lo), float(hi)


def bh_fdr(pvals: list[float]) -> list[float]:
    """Benjamini–Hochberg FDR-adjusted q-values.

    Input: list of raw p-values (unsorted, any order).
    Output: list of q-values in the SAME order, clipped to ``[0, 1]``.
    """
    p = np.array(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = q
    return np.clip(out, 0.0, 1.0).tolist()
