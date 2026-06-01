"""Cross-patient hypothesis-testing primitives.

Helpers used by every cross-patient test. The first five (wilcoxon_z,
rank_biserial, boot_ci_mean, bh_fdr, cluster_stats) were elevated to the
library in 2026-04 after being duplicated across 3+ scripts under
`scripts/01_compute/`. `surrogate_p_value` and `loo_sensitivity` were
elevated 2026-05-28 after being re-rolled in audit_62/65/66/67/70 (the
matched-strength surrogate family). See
`.agents/guides/02_methods/h2-metrics.md` for the scientific role and
`feedback_matched_strength_mandatory` for why every FC-derived cohort
claim must run a surrogate p-value.

Signatures follow the `rigorous_hypothesis_test.py` defaults (includes an
optional `rng` on the bootstrap for reproducibility).
"""
from __future__ import annotations

from typing import Callable, Optional, Sequence

import numpy as np
from scipy import stats


__all__ = [
    "wilcoxon_z",
    "rank_biserial",
    "boot_ci_mean",
    "bh_fdr",
    "cluster_stats",
    "surrogate_p_value",
    "loo_sensitivity",
    "regression_slope_through_origin",
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


def surrogate_p_value(
    obs: float,
    surr_array: np.ndarray,
    tail: str = "upper",
    *,
    add_one: bool = True,
) -> float:
    """One-tailed empirical p-value from a surrogate / permutation array.

    Canonical formula across the matched-strength surrogate family
    (audit_62/65/66/67/70). Computes the fraction of surrogate values
    that are at least as extreme as the observed value in the requested
    tail.

    Parameters
    ----------
    obs
        Observed statistic (scalar).
    surr_array
        Surrogate / permutation distribution of the statistic. NaNs are
        silently dropped. If empty after NaN removal, returns ``nan``.
    tail
        ``"upper"`` (default): p = P(surr ≥ obs). Use when "larger =
        more extreme" (trace mass under the locked T_d > 0 = TRACE
        convention). ``"lower"``: p = P(surr ≤ obs).
    add_one
        Phipson–Smyth (2010) unbiased estimator: numerator and
        denominator each receive +1 (counts the observed value as one
        permutation). Default True. Set False for naive
        ``mean(surr ≥ obs)`` — only appropriate if the surrogate
        ensemble already includes the observed configuration.

    Returns
    -------
    float
        p-value in ``[0, 1]``.
    """
    surr = np.asarray(surr_array, dtype=float).ravel()
    surr = surr[np.isfinite(surr)]
    n = surr.size
    if n == 0 or not np.isfinite(obs):
        return float("nan")
    if tail == "upper":
        n_extreme = int(np.sum(surr >= obs))
    elif tail == "lower":
        n_extreme = int(np.sum(surr <= obs))
    else:
        raise ValueError(f"tail must be 'upper' or 'lower', got {tail!r}")
    if add_one:
        return (n_extreme + 1) / (n + 1)
    return n_extreme / n


def loo_sensitivity(
    cohort_values: np.ndarray,
    test_fn: Callable[[np.ndarray], float],
    *,
    labels: Optional[Sequence[str]] = None,
) -> dict:
    """Leave-one-out sensitivity for any cohort-level scalar test.

    Mandated since 2026-05-19 (see ``feedback_no_single_patient_p_driven``):
    every cohort p-value (Wilcoxon, cluster-permutation, anatomy A1) must
    be paired with a leave-one-out report. n=10 Wilcoxon is vulnerable to
    direction-outliers; LOO surfaces single-patient-leveraged verdicts
    (e.g., "p=0.005 overall, p=0.07 dropping Pat_XX"). LOO is descriptive,
    never a hardcoded gate.

    Parameters
    ----------
    cohort_values
        1D array of per-patient values, length ``n`` (rows = patients).
        2D arrays are passed to ``test_fn`` row-sliced (i.e. ``arr[mask]``
        where ``mask`` removes one row). Higher-dim arrays index along
        axis 0.
    test_fn
        Callable that takes a (n-1)-row subset of ``cohort_values`` and
        returns a scalar (typically a p-value, sometimes an effect size).
    labels
        Optional patient labels matching ``cohort_values`` axis 0.
        Default ``range(n)``.

    Returns
    -------
    dict with keys
        ``full`` — ``test_fn(cohort_values)`` (the n-patient value).
        ``loo`` — ``np.ndarray`` of n LOO values.
        ``worst`` — max of the LOO values (worst-case p when the metric
            IS a p-value, smaller-is-better; caller swaps sign if needed).
        ``worst_patient`` — label of the patient whose removal produced
            ``worst``. ``None`` if all LOO values are NaN.
    """
    arr = np.asarray(cohort_values)
    n = arr.shape[0]
    if labels is None:
        labels = [str(i) for i in range(n)]
    else:
        labels = list(labels)
        if len(labels) != n:
            raise ValueError(
                f"labels length {len(labels)} ≠ cohort_values axis 0 ({n})"
            )

    full = float(test_fn(arr))
    loo = np.full(n, np.nan)
    for i in range(n):
        keep = np.ones(n, dtype=bool)
        keep[i] = False
        try:
            loo[i] = float(test_fn(arr[keep]))
        except Exception:
            loo[i] = np.nan

    finite = loo[np.isfinite(loo)]
    if finite.size == 0:
        return {"full": full, "loo": loo, "worst": float("nan"), "worst_patient": None}
    worst_idx = int(np.nanargmax(loo))
    return {
        "full": full,
        "loo": loo,
        "worst": float(loo[worst_idx]),
        "worst_patient": labels[worst_idx],
    }


def regression_slope_through_origin(
    x: np.ndarray, y: np.ndarray
) -> tuple[float, float]:
    """Through-origin regression slope and r-squared of ``y ≈ slope · x``.

    The slope is the ordinary-least-squares coefficient of the no-intercept
    model ``y = slope · x + ε``, i.e. the projection of ``y`` onto ``x``::

        slope     = <x, y> / <x, x>
        r_squared = pearson(x, y) ** 2

    ``slope`` is **asymmetric** in ``x ↔ y`` (``slope(y, x) = <x, y> / <y, y>``);
    by convention ``x`` is the cause / predictor and ``y`` the response. When
    ``(x, y) = (Δ_task, Δ_rest)`` on a per-pair LRG distance, ``slope`` reads as
    the fraction of the task-induced per-pair shift recovered in the post-task
    rest. ``r_squared`` is the mean-centred Pearson squared (fraction of the
    response variance linearly explained); it is rotation-invariant and so the
    cross-primitive-comparable companion to the magnitude-bearing ``slope``.

    NOTE — naming: this is the ``s_TR`` measure. The bare letter ``β`` is
    reserved for the 13–30 Hz band project-wide; never name a regression slope
    ``beta`` (see ``feedback_no_beta_for_regression_slope``).

    Parameters
    ----------
    x, y : np.ndarray
        Paired 1-D vectors of equal length. Non-finite entries in either
        vector are dropped pairwise before the fit.

    Returns
    -------
    tuple[float, float]
        ``(slope, r_squared)``. Returns ``(nan, nan)`` when the inputs are
        empty, length-mismatched, reduce to fewer than two finite pairs, or
        when ``<x, x> == 0`` (no predictor variation, slope undefined).
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    if x.shape != y.shape or x.size == 0:
        return np.nan, np.nan

    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    if x.size < 2:
        return np.nan, np.nan

    xx = float(np.dot(x, x))
    if xx == 0.0:
        return np.nan, np.nan
    slope = float(np.dot(x, y) / xx)

    # r_squared = mean-centred Pearson², guarded against zero-variance y.
    if np.ptp(x) == 0.0 or np.ptp(y) == 0.0:
        r_squared = np.nan
    else:
        r = float(np.corrcoef(x, y)[0, 1])
        r_squared = r * r
    return slope, r_squared
