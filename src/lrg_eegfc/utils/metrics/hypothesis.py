"""Cross-patient hypothesis-testing primitives.

Five helpers used by every H2-family cross-patient test. Elevated to the
library in 2026-04 after being duplicated across 3+ scripts under
`scripts/01_compute/`. See `.agents/guides/02_methods/h2-metrics.md` for
the scientific role.

Signatures follow the `rigorous_hypothesis_test.py` defaults (includes an
optional `rng` on the bootstrap for reproducibility).
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from scipy import stats


__all__ = [
    "wilcoxon_z",
    "rank_biserial",
    "boot_ci_mean",
    "bh_fdr",
    "cluster_stats",
]


N_BOOT_DEFAULT = 10_000


def wilcoxon_z(x: np.ndarray) -> tuple[float, float]:
    """One-sample Wilcoxon signed-rank test against 0, one-sided 'greater'.

    Returns ``(z, p)``. z is the normal-approximation z-score (no ties
    correction — adequate for continuous contrasts where exact zeros are
    vanishingly rare). Exact zeros are dropped; n<3 returns ``(nan, nan)``.
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

    Unbiased effect-size estimator in ``[-1, +1]``. ``+1`` means every
    sample is positive and larger in absolute value than any negative
    sample. Drops exact zeros; n<3 returns ``nan``.
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
    """Return ``(mean, lo, hi)`` 95% bootstrap CI of the sample mean.

    ``B`` = bootstrap resamples (default 10 000). ``rng`` is an optional
    ``numpy.random.Generator`` for reproducibility; if omitted uses
    ``default_rng(0)``.
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

    Input: raw p-values (any order). Output: q-values in the SAME order,
    clipped to ``[0, 1]``.
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


def cluster_stats(z: np.ndarray, thresh: float) -> list[tuple[int, int, float]]:
    """Supra-threshold run identification for cluster-based permutation tests.

    Returns ``[(start_idx, end_idx_inclusive, cluster_mass), ...]`` for each
    contiguous run where ``z > thresh``. Used by the Maris-Oostenveld
    sign-flip null on 1D axes (k, scale, λ).
    """
    sup = np.isfinite(z) & (z > thresh)
    out: list[tuple[int, int, float]] = []
    in_run = False
    start = 0
    mass = 0.0
    for i, s in enumerate(sup):
        if s and not in_run:
            in_run = True
            start = i
            mass = float(z[i])
        elif s:
            mass += float(z[i])
        elif in_run:
            out.append((start, i - 1, mass))
            in_run = False
    if in_run:
        out.append((start, len(sup) - 1, mass))
    return out
