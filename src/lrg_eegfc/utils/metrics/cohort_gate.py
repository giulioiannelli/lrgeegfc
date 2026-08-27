"""Canonical cross-patient gate for surrogate-referenced cohort claims.

A *cohort cell* is one (statistic, condition) pair -- e.g. the cross-phase trace
``rho_sym`` at one (band, scale) -- for which every patient contributes

* one observed value ``o_k``, and
* a length-``R`` surrogate ensemble ``S_k = {s_{k,1}, ..., s_{k,R}}`` produced by
  running the *identical* pipeline on a null-transformed input.

This module turns ``K`` such (obs, ensemble) pairs into one reportable verdict.
It is deliberately the **only** place that decision is made, so that every lane
of a multi-lane analysis inherits the same contract instead of re-rolling a
Wilcoxon with its own conventions.

The contract, and why each clause exists
----------------------------------------
1. **Margin, never raw value.** The tested quantity is
   ``m_k = o_k - median_r s_{k,r}`` (:func:`patient_margin`), never ``o_k``.
   Per-patient observed values co-vary with the *height of that patient's own
   null* -- graphs differ in size, density and strength distribution, so the
   null floor is patient-specific and a raw-value spread across patients mixes
   baseline non-identifiability into what looks like effect size.
   :func:`null_height_covariation` measures that covariation for both the raw
   value and the margin, so the confound is quantified rather than assumed away.
2. **Per condition, never best-of.** :func:`gate_grid` evaluates every cell of
   the (condition x condition) grid and returns all of them. Nothing in this
   module selects a maximum, and no scalar collapse of a swept axis is offered.
3. **Whole-grid multiplicity.** :func:`gate_grid` applies Benjamini-Hochberg
   across the entire submitted grid by default, because the grid *is* the
   coordinated family: the scientific claim is the pattern of which cells hold
   and which do not, so every cell was interrogated in service of one claim and
   every cell must pay for it. Per-group BH is available but must be justified
   by an argument that the groups are separate claims.
4. **Leave-one-patient-out of the verdict.** :func:`cohort_margin_gate` reports
   the p-value under each single-patient drop, and :func:`loo_grid_flips`
   recomputes the *whole grid including its BH step* under each drop and counts
   how many cells change verdict. A LOO number that is constant by construction
   is worthless; what is reported here is the swing in the conclusion.
5. **Bootstrap CI on the cohort effect**, from
   :func:`~lrg_eegfc.utils.metrics.hypothesis.boot_ci_mean`, alongside the
   rank-biserial effect size.
6. **Counts are descriptive.** ``n_above`` / ``frac_pos`` are reported because
   readers want them, and are never consulted by the decision. The signed-rank
   test is the gate.

Statistic, substrate and null are all injected, not chosen here
--------------------------------------------------------------
This module never loads data, never builds a graph, never draws a surrogate and
never computes the statistic. It consumes ``(obs, surr)`` arrays, so it gates a
cross-phase trace, a partial correlation, a detector AUC or a localization score
identically. Swapping the FC transform, the backbone, the surrogate family or
the functional changes what is handed in and changes nothing here -- the
decoupling that lets a locked gate outlive an unlocked pipeline.

Multiplicity
------------
The default is whole-grid BH, because the grid is the coordinated family. But a
*swept axis* -- diffusion scale, frequency, lag -- is a smooth curve whose
positions are near-duplicates, and BH charges full price for every one of them.
:func:`effective_tests` measures how many independent tests the sweep is really
worth, and :func:`axis_cluster_gate` offers the alternative that respects the
smoothness: one sign-flip cluster-mass test for the whole axis, so a family of
six bands is six tests rather than six times the axis length. Report the grid in
full (its shape is a result) and choose the family deliberately.

Calibration
-----------
:func:`calibrate_from_surrogates` runs the gate on data where the null is known
to be true *by construction*: it promotes one surrogate realization to the role
of "observed" and tests it against the remaining ensemble. Under the null the
observed value is exchangeable with its own surrogates, so the resulting
p-values must be uniform; the empirical rejection rate at a nominal level is the
gate's false-positive rate on the real graphs, not on a Gaussian toy.
:func:`calibrate_synthetic` provides the toy as a cross-check.

**Conditional statistics must be calibrated before they are gated.** A partial
correlation entangles the conditioning variable with the estimator, and can
return systematically positive values from input containing no signal at all --
this is not hypothetical here: a sham arc constructed entirely inside pre-task
rest, where no consolidation can exist, produced a significantly positive
conditional trace (p = 0.007) that exceeded the real value. So
``gate_grid(..., calibrate=True)`` measures each cell's false-positive rate and,
by default, **withholds** ``p`` and ``q`` for any cell whose statistic/null pair
fails. A p-value from a miscalibrated pair is not conservative, not liberal and
not interpretable, and reporting it with a caveat is worse than not reporting
it. Note what this does and does not cover: the held-out-realization test
calibrates the statistic against *the null actually being used*. It cannot
detect a bias that the null shares -- for that the statistic needs a
**data-based placebo** (a no-signal arc built from the recording itself), which
is a separate deliverable and a separate dependency.

General statistics primitive -- no dataset-, band- or manuscript-local scope.
"""
from __future__ import annotations

from typing import Iterable, Mapping, Optional, Sequence

import numpy as np
from scipy import stats
from scipy.stats import kstest, spearmanr

from lrg_eegfc.utils.metrics.hypothesis import (
    bh_fdr,
    boot_ci_mean,
    cluster_stats,
    rank_biserial,
    wilcoxon_z,
)

__all__ = [
    "patient_margin",
    "null_height_covariation",
    "cohort_margin_gate",
    "gate_grid",
    "loo_grid_flips",
    "effective_tests",
    "axis_cluster_gate",
    "calibrate_from_surrogates",
    "calibrate_synthetic",
    "DESCRIPTIVE_ALPHA",
]

#: Level at which the *descriptive* per-patient counts are tabulated. It is not
#: an acceptance threshold: nothing in this module branches on it.
DESCRIPTIVE_ALPHA = 0.05

_EPS = 1e-12


def _as_cell(obs, surr) -> tuple[np.ndarray, np.ndarray]:
    """Coerce a cell to ``(obs: (K,), surr: (K, R))`` float arrays.

    Rows are patients; columns are surrogate realizations. NaNs are permitted
    in either (failed realizations, missing cells) and are handled by the
    ``nan``-aware reductions downstream.
    """
    obs = np.asarray(obs, dtype=float).ravel()
    surr = np.asarray(surr, dtype=float)
    if surr.ndim == 1:
        surr = surr[None, :]
    if surr.shape[0] != obs.shape[0]:
        raise ValueError(
            f"surr must be (K, R) matching obs (K,); got {surr.shape} vs {obs.shape}"
        )
    return obs, surr


def patient_margin(obs, surr) -> np.ndarray:
    """Per-patient margin ``m_k = o_k - median_r s_{k,r}``.

    The quantity every cohort test in this project is run on. Subtracting the
    patient's *own* surrogate median removes the patient-specific null floor,
    which is what makes margins comparable across patients whose graphs differ
    in size, density and strength distribution.
    """
    obs, surr = _as_cell(obs, surr)
    return obs - np.nanmedian(surr, axis=1)


def null_height_covariation(obs, surr) -> dict:
    """Spearman covariation of raw value and of margin with the null height.

    Returns ``rho_obs`` / ``p_obs`` for ``Spearman(o_k, median_r s_{k,r})`` and
    ``rho_margin`` / ``p_margin`` for ``Spearman(m_k, median_r s_{k,r})``.

    A large positive ``rho_obs`` is the diagnostic that raw observed values are
    partly a readout of each patient's own baseline rather than of effect size;
    ``rho_margin`` shows how much of that survives the margin subtraction. Both
    are reported so the correction is auditable instead of assumed.
    """
    obs, surr = _as_cell(obs, surr)
    p50 = np.nanmedian(surr, axis=1)
    m = obs - p50
    ok = np.isfinite(obs) & np.isfinite(p50)
    out = dict(rho_obs=np.nan, p_obs=np.nan, rho_margin=np.nan, p_margin=np.nan)
    if ok.sum() >= 3:
        r1, q1 = spearmanr(obs[ok], p50[ok])
        r2, q2 = spearmanr(m[ok], p50[ok])
        out.update(rho_obs=float(r1), p_obs=float(q1),
                   rho_margin=float(r2), p_margin=float(q2))
    return out


def _wilcoxon_margin_p(margin: np.ndarray) -> float:
    """One-sided (``greater``) signed-rank p on a margin vector; ``nan`` if degenerate.

    Uses SciPy's default method, which is the **exact** permutation distribution
    at the cohort sizes this project works with (``n <= 50``). At ``n = 10`` the
    normal approximation is materially wrong in the tail that matters -- and the
    exact test is what every incumbent script in the repository ran, so the
    locked gate stays numerically comparable to them.
    :func:`~lrg_eegfc.utils.metrics.hypothesis.wilcoxon_z` is still used for the
    reported ``z``, which is a normal-approximation effect summary, not the gate.
    """
    x = np.asarray(margin, dtype=float)
    x = x[np.isfinite(x)]
    x = x[x != 0.0]
    if x.size < 3:
        return float("nan")
    try:
        return float(stats.wilcoxon(x, alternative="greater").pvalue)
    except ValueError:
        return float("nan")


def cohort_margin_gate(
    obs,
    surr,
    *,
    labels: Optional[Sequence[str]] = None,
    n_boot: int = 10_000,
    rng: Optional[np.random.Generator] = None,
    full: bool = True,
    expect_n: Optional[int] = None,
    exclusion_reason: Optional[str] = None,
) -> dict:
    """The canonical single-cell cohort gate.

    Parameters
    ----------
    obs
        ``(K,)`` observed statistic, one per patient.
    surr
        ``(K, R)`` surrogate ensemble, one row per patient. Rows need not be
        aligned across patients (this gate never compares realization ``r`` of
        one patient with realization ``r`` of another).
    labels
        Optional patient identifiers, length ``K``. Used only to name the
        leave-one-out worst case.
    n_boot
        Bootstrap resamples for the CI on the mean margin.
    rng
        Optional generator for the bootstrap.
    full
        When ``False``, skips the bootstrap and the leave-one-out sweep and
        returns only ``p`` plus the summaries. Used by the calibration loops,
        which call the gate thousands of times.
    expect_n
        The cohort the analysis is entitled to. Pass it and a shortfall raises
        unless ``exclusion_reason`` is also given. The standing policy is that
        patients are not dropped; the two narrow exceptions (a single
        cohort-level anti-aligned patient at the probe under test, or genuinely
        corrupt data) must be argued at the point of use, not discovered later
        in a row count. This is the mechanism that makes a silent exclusion
        impossible rather than merely discouraged.
    exclusion_reason
        Free text recorded in the returned dict as ``exclusion_reason``. Its
        presence is what authorises ``n_pat < expect_n``.

    Returns
    -------
    dict
        ``p``             -- the gate. One-sided Wilcoxon signed-rank
                             (``greater``) on the per-patient margins.
        ``z``             -- normal-approximation z of that test.
        ``rank_biserial`` -- effect size of the same test, in ``[-1, 1]``.
        ``n_pat``         -- patients contributing a finite margin.
        ``obs_med``, ``surr_med``, ``margin_med`` -- cohort medians.
        ``margin_mean``, ``ci_lo``, ``ci_hi``     -- bootstrap 95% CI of the
                             mean margin (``full=True`` only).
        ``rho_obs_null``, ``rho_margin_null`` -- the covariation diagnostic of
                             :func:`null_height_covariation`.
        ``loo_p``         -- ``(K,)`` p under each single-patient drop.
        ``loo_p_max``, ``loo_p_min``, ``loo_worst_patient``,
        ``loo_log10_swing`` -- the swing in the verdict under LOO
                             (``full=True`` only).
        ``n_above``, ``frac_pos`` -- **descriptive only**: number of patients
                             whose own upper-tail surrogate p is below
                             :data:`DESCRIPTIVE_ALPHA`, and fraction of patients
                             whose observed value exceeds their own surrogate
                             median. Never consulted by the decision.
    """
    obs, surr = _as_cell(obs, surr)
    if labels is None:
        labels = [str(i) for i in range(obs.shape[0])]
    labels = list(labels)
    if expect_n is not None and obs.shape[0] < expect_n and not exclusion_reason:
        raise ValueError(
            f"cohort gate received {obs.shape[0]} patients but the analysis is "
            f"entitled to {expect_n}. Pass exclusion_reason=... to record why "
            f"they were dropped; silent exclusions are not permitted."
        )

    margin = patient_margin(obs, surr)
    finite = np.isfinite(margin)
    m = margin[finite]
    lab = [labels[i] for i in np.flatnonzero(finite)]

    p = _wilcoxon_margin_p(m)
    z, _ = wilcoxon_z(m)

    # descriptive per-patient counts (never a gate)
    tail_p = np.full(obs.shape[0], np.nan)
    for k in range(obs.shape[0]):
        col = surr[k][np.isfinite(surr[k])]
        if col.size and np.isfinite(obs[k]):
            tail_p[k] = (1 + int(np.sum(col >= obs[k]))) / (col.size + 1)

    out = dict(
        p=p,
        z=float(z) if np.isfinite(z) else np.nan,
        rank_biserial=float(rank_biserial(m)),
        n_pat=int(m.size),
        obs_med=float(np.nanmedian(obs)) if np.isfinite(obs).any() else np.nan,
        surr_med=float(np.nanmedian(np.nanmedian(surr, axis=1))),
        margin_med=float(np.median(m)) if m.size else np.nan,
        n_above=int(np.nansum(tail_p < DESCRIPTIVE_ALPHA)),
        frac_pos=float(np.mean(m > 0)) if m.size else np.nan,
        exclusion_reason=exclusion_reason,
    )
    if not full:
        # the covariation diagnostic and the CI/LOO sweep are the expensive part;
        # the calibration loops call this thousands of times and only need `p`.
        out.update(rho_obs_null=np.nan, rho_margin_null=np.nan,
                   margin_mean=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                   loo_p=np.array([]), loo_p_max=np.nan, loo_p_min=np.nan,
                   loo_worst_patient=None, loo_log10_swing=np.nan)
        return out

    cov = null_height_covariation(obs, surr)
    out["rho_obs_null"] = cov["rho_obs"]
    out["rho_margin_null"] = cov["rho_margin"]

    mean, lo, hi = boot_ci_mean(m, B=n_boot, rng=rng)
    loo = np.full(m.size, np.nan)
    for i in range(m.size):
        keep = np.ones(m.size, bool)
        keep[i] = False
        loo[i] = _wilcoxon_margin_p(m[keep])
    fin = np.isfinite(loo)
    if fin.any():
        i_worst = int(np.nanargmax(loo))
        swing = float(np.log10(max(np.nanmax(loo), _EPS))
                      - np.log10(max(np.nanmin(loo), _EPS)))
        out.update(loo_p=loo, loo_p_max=float(np.nanmax(loo)),
                   loo_p_min=float(np.nanmin(loo)),
                   loo_worst_patient=lab[i_worst], loo_log10_swing=swing)
    else:
        out.update(loo_p=loo, loo_p_max=np.nan, loo_p_min=np.nan,
                   loo_worst_patient=None, loo_log10_swing=np.nan)
    out.update(margin_mean=mean, ci_lo=lo, ci_hi=hi)
    return out


def gate_grid(
    cells: Iterable[tuple[Mapping, np.ndarray, np.ndarray]],
    *,
    labels: Optional[Sequence[str]] = None,
    q_group: Optional[Sequence[str]] = None,
    n_boot: int = 10_000,
    rng: Optional[np.random.Generator] = None,
    full: bool = True,
    expect_n: Optional[int] = None,
    exclusion_reason: Optional[str] = None,
    calibrate: bool = False,
    calibration_draws: int = 200,
    fpr_tolerance: float = 0.10,
    refuse_uncalibrated: bool = True,
):
    """Run :func:`cohort_margin_gate` over a grid of cells and add BH q-values.

    Parameters
    ----------
    cells
        Iterable of ``(key_mapping, obs, surr)``. ``key_mapping`` supplies the
        identifying columns (e.g. ``{"band": "beta", "s": 5.6}``).
    labels
        Patient identifiers shared by all cells.
    q_group
        ``None`` (default) -> **whole-grid** BH across every submitted cell,
        because the grid is the coordinated family being interrogated. Passing
        column names restricts BH to within each group; do that only with an
        explicit argument that the groups are independent claims.
    full
        Forwarded to :func:`cohort_margin_gate`.
    calibrate
        Run :func:`calibrate_from_surrogates` on every cell and add ``fpr_05``,
        ``fpr_01`` and ``calibrated``. **Required for any statistic whose
        estimator can manufacture signal from null input** -- conditional and
        partial correlations above all, which entangle the conditioning variable
        with the estimator and are known in this project to have returned a
        significantly positive value on a no-signal placebo. An ordinary
        difference-of-correlations statistic is much less exposed, but the check
        is cheap enough to run on everything.
    calibration_draws
        Held-out realizations per cell.
    fpr_tolerance
        A cell is marked ``calibrated=False`` when its measured false-positive
        rate at the 0.05 level exceeds this. It is a *calibration* tolerance --
        a statement about whether the test is a test at all -- and never an
        acceptance threshold on the science.
    refuse_uncalibrated
        When ``True`` (default) an uncalibrated cell's ``p`` and ``q`` are set to
        NaN rather than reported. A p-value from a miscalibrated statistic/null
        pair is not conservative, not liberal, and not interpretable; it is
        withheld. The measured FPR is still reported so the failure is visible.

    Returns
    -------
    pandas.DataFrame
        One row per cell: the key columns, every scalar the gate returns, and
        ``q``. The per-patient ``loo_p`` vector is dropped (use
        :func:`loo_grid_flips` for the grid-level LOO).
    """
    import pandas as pd

    rows = []
    for keys, obs, surr in cells:
        res = cohort_margin_gate(obs, surr, labels=labels, n_boot=n_boot,
                                 rng=rng, full=full, expect_n=expect_n,
                                 exclusion_reason=exclusion_reason)
        res.pop("loo_p", None)
        if calibrate:
            cal = calibrate_from_surrogates(np.asarray(surr, float),
                                            n_draws=calibration_draws, rng=rng)
            res["fpr_05"] = cal["fpr_0.05"]
            res["fpr_01"] = cal["fpr_0.01"]
            res["calibrated"] = bool(np.isfinite(cal["fpr_0.05"])
                                     and cal["fpr_0.05"] <= fpr_tolerance)
            if refuse_uncalibrated and not res["calibrated"]:
                res["p"] = np.nan
        rows.append({**dict(keys), **res})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["q"] = np.nan
    if q_group is None:
        ok = df["p"].notna()
        if ok.any():
            df.loc[ok, "q"] = bh_fdr(df.loc[ok, "p"].to_numpy())
    else:
        for _, idx in df.groupby(list(q_group)).groups.items():
            sub = df.loc[idx]
            ok = sub["p"].notna()
            if ok.any():
                df.loc[sub.index[ok], "q"] = bh_fdr(sub.loc[ok, "p"].to_numpy())
    return df


def loo_grid_flips(
    cells: Sequence[tuple[Mapping, np.ndarray, np.ndarray]],
    *,
    labels: Sequence[str],
    q_level: float = DESCRIPTIVE_ALPHA,
    q_group: Optional[Sequence[str]] = None,
):
    """Leave-one-patient-out of the **verdict**: recompute the whole grid, count flips.

    For every patient in turn the entire grid -- p-values *and* the BH step,
    since dropping a patient changes the multiplicity correction as well as the
    individual tests -- is recomputed without that patient, and each cell's
    ``q < q_level`` verdict is compared with the full-cohort verdict.

    This is the LOO the project needs. Leaving a patient out of the *statistic*
    and reporting how the statistic moved is uninformative (and in one earlier
    incarnation was constant by construction); leaving them out of the
    *conclusion* is the thing a reader wants to know.

    ``q_level`` is a reporting level for tabulating flips, not an acceptance
    threshold: no function in this module branches on the returned verdicts.

    Returns
    -------
    (per_patient, per_cell) : tuple[pandas.DataFrame, pandas.DataFrame]
        ``per_patient`` -- one row per dropped patient: ``n_clear_full``,
        ``n_clear_loo``, ``n_lost`` (cleared with everyone, not without them),
        ``n_gained``, ``n_flip``.
        ``per_cell`` -- the key columns plus ``clear_full``, ``n_drop_lost``
        (how many single-patient drops cost this cell its verdict) and
        ``loo_p_max``.
    """
    import pandas as pd

    cells = list(cells)
    labels = list(labels)
    full_df = gate_grid(cells, labels=labels, q_group=q_group, full=False)
    if full_df.empty:
        return full_df, full_df
    key_cols = [c for c in dict(cells[0][0]).keys()]
    clear_full = (full_df["q"] < q_level).to_numpy()

    lost = np.zeros(len(cells), dtype=int)
    loo_p_max = np.full(len(cells), np.nan)
    rows = []
    for i, pat in enumerate(labels):
        keep = np.ones(len(labels), bool)
        keep[i] = False
        sub_cells = [(k, np.asarray(o, float)[keep], np.asarray(s, float)[keep])
                     for k, o, s in cells]
        sub_labels = [labels[j] for j in np.flatnonzero(keep)]
        d = gate_grid(sub_cells, labels=sub_labels, q_group=q_group, full=False)
        clear_loo = (d["q"] < q_level).to_numpy()
        loo_p_max = np.fmax(loo_p_max, d["p"].to_numpy())
        lost += (clear_full & ~clear_loo).astype(int)
        rows.append(dict(dropped=pat,
                         n_clear_full=int(clear_full.sum()),
                         n_clear_loo=int(clear_loo.sum()),
                         n_lost=int((clear_full & ~clear_loo).sum()),
                         n_gained=int((~clear_full & clear_loo).sum()),
                         n_flip=int((clear_full != clear_loo).sum())))
    per_cell = full_df[key_cols].copy()
    per_cell["clear_full"] = clear_full
    per_cell["n_drop_lost"] = lost
    per_cell["loo_p_max"] = loo_p_max
    return pd.DataFrame(rows), per_cell


# --------------------------------------------------------------------------- #
# multiplicity: how many tests were really run
# --------------------------------------------------------------------------- #
def effective_tests(M) -> dict:
    """How many independent tests a correlated sweep is really worth.

    ``M`` is ``(K, n_axis)`` -- one row per patient, one column per position on a
    swept axis (diffusion scale, frequency, lag). Benjamini-Hochberg treats the
    ``n_axis`` columns as ``n_axis`` separate tests, but a swept axis is a smooth
    curve: neighbouring positions are near-duplicates, and correcting as if they
    were independent can cost an order of magnitude of power for nothing.

    Two standard readouts of the correlation matrix ``C`` of the columns:

    ``n_eff_pr``
        participation ratio of the eigenvalues, ``(sum l)^2 / sum l^2`` -- the
        effective number of independent directions.
    ``n_eff_cn``
        the Cheverud-Nyholt style estimate ``1 + (n - 1)(1 - var(l)/n)``.

    Returns both plus ``n_axis`` and the mean off-diagonal correlation. This is a
    *diagnostic*, not a correction: report it alongside a whole-grid BH so the
    reader knows how conservative that correction is. Substituting ``n_eff`` for
    ``n`` in a BH denominator is not a validated procedure and this module does
    not do it.
    """
    M = np.asarray(M, dtype=float)
    ok = np.isfinite(M).all(axis=0)
    X = M[:, ok]
    n = X.shape[1]
    out = dict(n_axis=int(M.shape[1]), n_finite=int(n),
               n_eff_pr=np.nan, n_eff_cn=np.nan, mean_offdiag=np.nan)
    if n < 2 or X.shape[0] < 3:
        return out
    C = np.corrcoef(X, rowvar=False)
    C = np.where(np.isfinite(C), C, 0.0)
    lam = np.linalg.eigvalsh(C)
    lam = np.clip(lam, 0.0, None)
    ssum = lam.sum()
    if ssum > 0:
        out["n_eff_pr"] = float(ssum ** 2 / np.sum(lam ** 2))
    out["n_eff_cn"] = float(1.0 + (n - 1) * (1.0 - np.var(lam) / n))
    iu = np.triu_indices(n, 1)
    out["mean_offdiag"] = float(np.mean(C[iu]))
    return out


def axis_cluster_gate(
    M,
    *,
    n_perm: int = 10_000,
    rng: Optional[np.random.Generator] = None,
    z_thresh: float = 1.0,
) -> dict:
    """One p-value for a whole swept axis, via sign-flip cluster mass.

    ``M`` is ``(K, n_axis)`` of per-patient margins. Instead of asking "does the
    margin clear at position ``j``?" ``n_axis`` times, this asks the single
    question the sweep was designed to answer -- "is there a contiguous stretch
    of this axis where the cohort margin is positive?" -- and answers it against
    a sign-flip null that respects the axis's smoothness (Maris-Oostenveld).

    Procedure: form the per-position one-sample z of the margins, find
    supra-threshold runs with :func:`~lrg_eegfc.utils.metrics.hypothesis.cluster_stats`,
    take the largest cluster mass, and compare it with the same statistic under
    random sign flips of each patient's entire profile. Flipping whole profiles,
    not individual positions, is what preserves the correlation along the axis.

    This collapses a 28-position sweep to ONE test, so a family of six bands is
    six tests rather than 168 -- the honest multiplicity when the scientific
    question is per-band. It is a companion to, not a replacement for, the
    per-position grid: the grid is still reported in full because the shape of
    the curve is itself a result.

    ``z_thresh`` is the cluster-forming threshold. It selects which clusters are
    *formed*, not which are significant -- the permutation distribution is built
    with the same threshold, so the test remains exact whatever it is set to.

    Returns ``p``, ``mass``, ``cluster`` (start, end indices), ``n_clusters``,
    and ``z`` (the per-position z profile).
    """
    M = np.asarray(M, dtype=float)
    rng = rng or np.random.default_rng(0)
    ok = np.isfinite(M).all(axis=0)
    X = M[:, ok]
    K, n = X.shape
    out = dict(p=np.nan, mass=np.nan, cluster=None, n_clusters=0,
               z=np.full(M.shape[1], np.nan))

    def _mass(Y):
        mu = Y.mean(axis=0)
        sd = Y.std(axis=0, ddof=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            z = np.where(sd > _EPS, mu / (sd / np.sqrt(Y.shape[0])), 0.0)
        cl = cluster_stats(z, z_thresh)
        return (max((c[2] for c in cl), default=0.0), cl, z)

    if K < 3 or n < 1:
        return out
    obs_mass, cl, z = _mass(X)
    out["z"][ok] = z
    out["mass"] = float(obs_mass)
    out["n_clusters"] = len(cl)
    if cl:
        best = max(cl, key=lambda c: c[2])
        idx = np.flatnonzero(ok)
        out["cluster"] = (int(idx[best[0]]), int(idx[best[1]]))
    null = np.empty(n_perm)
    for i in range(n_perm):
        sgn = rng.choice((-1.0, 1.0), size=(K, 1))
        null[i] = _mass(X * sgn)[0]
    out["p"] = float((1 + int(np.sum(null >= obs_mass))) / (n_perm + 1))
    return out


# --------------------------------------------------------------------------- #
# calibration
# --------------------------------------------------------------------------- #
def calibrate_from_surrogates(
    surr,
    *,
    n_draws: int = 200,
    rng: Optional[np.random.Generator] = None,
    levels: Sequence[float] = (0.01, 0.05, 0.10),
) -> dict:
    """False-positive rate of the gate on data where the null is true by construction.

    Promotes surrogate realization ``r`` to the role of "observed" and tests it
    against the remaining ``R - 1`` realizations of the *same* patients. Under
    the null the true observation is exchangeable with its own surrogates, so
    the p-values this produces are draws from the gate's null distribution on
    the real graphs -- with the real ``K``, the real ``R``, the real
    heteroscedasticity across patients, and the real correlation structure of
    the surrogate ensemble. A Gaussian toy tests none of that.

    ``surr`` is ``(K, R)`` for one cell. Realization columns must be aligned
    across patients only in the weak sense that column ``r`` of every patient
    is drawn independently -- which is what the shuffle families here do.

    Returns ``p`` (the draws), ``ks_D`` / ``ks_p`` (Kolmogorov-Smirnov against
    ``Uniform(0, 1)``), and ``fpr_<level>`` for each nominal level.

    **Read the FPRs, not the KS p.** The exact signed-rank p at ``K = 10`` takes
    only ``2^10`` distinct values, so its null distribution is discrete and a KS
    test against the *continuous* uniform rejects by construction, no matter how
    well calibrated the gate is. ``ks_D`` is retained as a descriptive
    calibration error; the statement that matters is that ``P(p <= a) <= a`` at
    each nominal ``a``, which is what ``fpr_<level>`` measures.
    """
    surr = np.asarray(surr, dtype=float)
    if surr.ndim != 2:
        raise ValueError(f"surr must be (K, R); got shape {surr.shape}")
    K, R = surr.shape
    rng = rng or np.random.default_rng(0)
    draws = rng.choice(R, size=min(n_draws, R), replace=False)
    ps = []
    for r in draws:
        keep = np.ones(R, bool)
        keep[r] = False
        res = cohort_margin_gate(surr[:, r], surr[:, keep], full=False)
        if np.isfinite(res["p"]):
            ps.append(res["p"])
    ps = np.asarray(ps, float)
    out = dict(p=ps, n_draws=int(ps.size), K=int(K), R=int(R),
               ks_D=np.nan, ks_p=np.nan)
    if ps.size >= 5:
        k = kstest(ps, "uniform")
        out["ks_D"], out["ks_p"] = float(k.statistic), float(k.pvalue)
    for a in levels:
        out[f"fpr_{a:g}"] = float(np.mean(ps < a)) if ps.size else np.nan
    return out


def calibrate_synthetic(
    *,
    n_patients: int = 10,
    n_surrogates: int = 200,
    n_draws: int = 2000,
    rng: Optional[np.random.Generator] = None,
    levels: Sequence[float] = (0.01, 0.05, 0.10),
    heteroscedastic: bool = True,
) -> dict:
    """Toy-null calibration: independent Gaussian obs and surrogates.

    ``heteroscedastic=True`` gives each patient its own null mean and scale, the
    situation the margin subtraction exists to handle; ``False`` gives the
    homoscedastic ideal. Returns the same keys as
    :func:`calibrate_from_surrogates`.
    """
    rng = rng or np.random.default_rng(0)
    ps = np.full(n_draws, np.nan)
    for d in range(n_draws):
        if heteroscedastic:
            mu = rng.normal(0.0, 1.0, size=n_patients)
            sd = rng.uniform(0.3, 3.0, size=n_patients)
        else:
            mu = np.zeros(n_patients)
            sd = np.ones(n_patients)
        obs = mu + sd * rng.normal(size=n_patients)
        srr = mu[:, None] + sd[:, None] * rng.normal(size=(n_patients, n_surrogates))
        ps[d] = cohort_margin_gate(obs, srr, full=False)["p"]
    ps = ps[np.isfinite(ps)]
    out = dict(p=ps, n_draws=int(ps.size), K=int(n_patients),
               R=int(n_surrogates), ks_D=np.nan, ks_p=np.nan)
    if ps.size >= 5:
        k = kstest(ps, "uniform")
        out["ks_D"], out["ks_p"] = float(k.statistic), float(k.pvalue)
    for a in levels:
        out[f"fpr_{a:g}"] = float(np.mean(ps < a)) if ps.size else np.nan
    return out
