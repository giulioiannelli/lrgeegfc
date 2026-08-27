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

Substrate and null are injected, not chosen here
------------------------------------------------
This module never loads data, never builds a graph and never draws a surrogate.
It consumes ``(obs, surr)`` arrays. Swapping the FC transform, the backbone or
the surrogate family changes what is handed in and changes nothing here, which
is exactly the decoupling that lets a locked gate outlive an unlocked pipeline.

Calibration
-----------
:func:`calibrate_from_surrogates` runs the gate on data where the null is known
to be true *by construction*: it promotes one surrogate realization to the role
of "observed" and tests it against the remaining ensemble. Under the null the
observed value is exchangeable with its own surrogates, so the resulting
p-values must be uniform; the empirical rejection rate at a nominal level is the
gate's false-positive rate on the real graphs, not on a Gaussian toy.
:func:`calibrate_synthetic` provides the toy as a cross-check.

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
    rank_biserial,
    wilcoxon_z,
)

__all__ = [
    "patient_margin",
    "null_height_covariation",
    "cohort_margin_gate",
    "gate_grid",
    "loo_grid_flips",
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
