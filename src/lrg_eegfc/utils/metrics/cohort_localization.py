"""Cohort aggregation statistics for per-unit localization scores.

A *localization cell* is a set of patients, each contributing one observed
per-unit score ``m_u(p)`` (e.g. an endpoint-incidence demeaned per-system mean of
a per-pair concordance) and a matched-strength surrogate ensemble
``M_u(p) = {m_u(p; r)}_{r=1..R}``. The scientific question is whether the cohort's
observed score is enriched above the strength-preserving null.

There is no single "right" way to turn ``K`` per-patient (obs, surrogate-ensemble)
pairs into one cohort p-value, and the choice is **not innocent** — it can flip a
verdict on identical data. This module makes the choice explicit by computing four
cohort statistics side by side, each with a known failure mode:

- :func:`pooled_median_p` — **(a) POOLED** cohort-median-of-obs vs ``R`` pooled
  surrogate-medians. Floor ``1/(R+1)``, *decoupled from patient count* → testable at
  ``K=5``, but LIBERAL: one or two strong patients drag the cohort median; no
  cross-patient consistency is required. (The whole-graph β→OFC statistic:
  audit_83/151/158/160/171.)
- :func:`wilcoxon_p` — **(b) WILCOXON** signed-rank of ``m_u(p) − median_r M_u(p;r)``
  across patients. Floor ``1/2^K`` → structurally DOA at ``K=5`` (``p_min=0.031``,
  BH-dead over ≥2 units). CONSERVATIVE to the point of being uninformative at sparse
  coverage. (The mst@0.20 backbone statistic: scripts 16/17.)
- :func:`stouffer_p` — **(c) STOUFFER** fixed-effect combine of per-patient
  matched-strength z-scores. Powered (R-resolution via the per-patient z), but an
  OMNIBUS test ("evidence somewhere in the cohort") that a single patient can drive.
- :func:`sign_consistency_p` — **(d) SIGN-CONSISTENCY permutation.** Cohort statistic
  ``Θ = median_p z_u(p)`` (median, not mean → robust to one outlier patient); its null
  is the matched-strength ensemble itself (``Θ_null[r] = median_p z_u(p;r)``), giving an
  ``R``-resolution permutation p. Reported together with ``frac_pos`` = fraction of
  patients with obs above their own surrogate median. Powered AND consistency-favouring.

None of these is privileged here — :func:`all_statistics` returns all four so the
reader adjudicates by convergence. BH-FDR (whole-grid vs per-cell) and
leave-one-patient-out are applied by the caller via :func:`bh_columns` and
:func:`loo_worst` because they need the full cell grid / the cell's patient set.

General graph/statistics primitive — no manuscript-local scope. Reused by the
localization statistic bake-off (script 21) and the raw-vs-cophenetic value-add
(script "C4").
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm, wilcoxon

from lrg_eegfc.utils.metrics.hypothesis import bh_fdr

__all__ = [
    "per_patient_z",
    "per_patient_tail_p",
    "pooled_median_p",
    "wilcoxon_p",
    "stouffer_p",
    "sign_consistency_p",
    "frac_positive",
    "all_statistics",
    "loo_worst",
    "bh_columns",
]

_EPS = 1e-12


def _as_cell(obs, surr):
    """Coerce a cell to ``(obs: (K,), surr: (K, R))`` float arrays.

    ``obs`` is one score per patient; ``surr`` is that patient's length-``R``
    matched-strength surrogate ensemble (NaN allowed for failed realizations).
    Rows are aligned by realization index ``r`` across patients (the pooled and
    sign-consistency nulls compare like realizations), so every row must share the
    same length ``R`` — pad short rows with NaN before calling.
    """
    obs = np.asarray(obs, dtype=float)
    surr = np.asarray(surr, dtype=float)
    if surr.ndim != 2 or surr.shape[0] != obs.shape[0]:
        raise ValueError(f"surr must be (K, R) matching obs (K,); got {surr.shape} vs {obs.shape}")
    return obs, surr


def per_patient_z(obs, surr):
    """Per-patient matched-strength z ``(m − mean_r M) / sd_r M``; NaN if ``sd=0``."""
    obs, surr = _as_cell(obs, surr)
    mu = np.nanmean(surr, axis=1)
    sd = np.nanstd(surr, axis=1, ddof=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(sd > _EPS, (obs - mu) / sd, np.nan)
    return z


def per_patient_tail_p(obs, surr):
    """Per-patient upper-tail surrogate p ``(1 + #{M ≥ m}) / (n_finite + 1)``."""
    obs, surr = _as_cell(obs, surr)
    out = np.full(obs.shape[0], np.nan)
    for i in range(obs.shape[0]):
        col = surr[i][np.isfinite(surr[i])]
        if col.size:
            out[i] = (1 + int(np.sum(col >= obs[i]))) / (col.size + 1)
    return out


def frac_positive(obs, surr):
    """Fraction of patients whose obs exceeds their own surrogate median."""
    obs, surr = _as_cell(obs, surr)
    med = np.nanmedian(surr, axis=1)
    ok = np.isfinite(med)
    return float(np.mean(obs[ok] > med[ok])) if ok.any() else np.nan


def pooled_median_p(obs, surr, upper=True):
    """(a) POOLED cohort-median-of-obs vs ``R`` cross-patient surrogate-medians.

    ``M_surr[r] = nanmedian_p M_u(p; r)`` over realizations gives an ``R``-length null
    of the cohort-median statistic; ``p = (1 + #{M_surr ⋛ M_obs}) / (n_finite + 1)``.
    Floor ``1/(R+1)`` independent of ``K`` (LIBERAL — see module docstring).
    """
    obs, surr = _as_cell(obs, surr)
    m_obs = float(np.nanmedian(obs))
    m_surr = np.nanmedian(surr, axis=0)
    m_surr = m_surr[np.isfinite(m_surr)]
    if m_surr.size == 0:
        return np.nan
    hits = np.sum(m_surr >= m_obs) if upper else np.sum(m_surr <= m_obs)
    return float((1 + int(hits)) / (m_surr.size + 1))


def wilcoxon_p(obs, surr, alternative="greater"):
    """(b) WILCOXON signed-rank of ``obs − median_r surr`` across patients.

    Floor ``1/2^K`` → DOA at sparse coverage (``K=5`` → ``p_min=0.031``). Returns 1.0
    when all differences are zero / degenerate.
    """
    obs, surr = _as_cell(obs, surr)
    diff = obs - np.nanmedian(surr, axis=1)
    diff = diff[np.isfinite(diff)]
    if diff.size < 1 or np.allclose(diff, 0.0):
        return 1.0
    try:
        _, p = wilcoxon(diff, alternative=alternative)
    except ValueError:
        return 1.0
    return float(p)


def stouffer_p(obs, surr, upper=True):
    """(c) STOUFFER fixed-effect combine of per-patient matched-strength z.

    ``Z_c = Σ_p z_p / √K_finite``; ``p = Φ(−Z_c)`` (upper). Powered but OMNIBUS —
    one strong patient can carry ``Z_c`` (see module docstring).
    """
    z = per_patient_z(obs, surr)
    z = z[np.isfinite(z)]
    if z.size == 0:
        return np.nan
    zc = float(np.sum(z) / np.sqrt(z.size))
    return float(norm.sf(zc)) if upper else float(norm.cdf(zc))


def sign_consistency_p(obs, surr, upper=True):
    """(d) SIGN-CONSISTENCY permutation p of the cohort median-of-z statistic.

    ``Θ = median_p z_p(obs)`` with ``z`` the per-patient matched-strength z; the null
    is the matched-strength ensemble itself, ``Θ_null[r] = median_p z_p(surr r)`` where
    ``z_p(surr r) = (M_u(p;r) − mean_r M_u(p)) / sd_r M_u(p)`` (each surrogate realization
    standardised against its own patient's ensemble). ``p = (1 + #{Θ_null ⋛ Θ}) /
    (n_finite + 1)``. The MEDIAN (not mean) over patients makes ``Θ`` robust to a single
    outlier patient → powered AND consistency-favouring. Pair with :func:`frac_positive`.
    """
    obs, surr = _as_cell(obs, surr)
    mu = np.nanmean(surr, axis=1, keepdims=True)
    sd = np.nanstd(surr, axis=1, ddof=1, keepdims=True)
    good = (sd[:, 0] > _EPS) & np.isfinite(obs)
    if good.sum() == 0:
        return np.nan
    obs, surr, mu, sd = obs[good], surr[good], mu[good], sd[good]
    z_obs = (obs - mu[:, 0]) / sd[:, 0]
    theta = float(np.median(z_obs))
    with np.errstate(invalid="ignore", divide="ignore"):
        z_surr = (surr - mu) / sd                     # (K, R) standardised realizations
    theta_null = np.nanmedian(z_surr, axis=0)         # (R,) cohort median per realization
    theta_null = theta_null[np.isfinite(theta_null)]
    if theta_null.size == 0:
        return np.nan
    hits = np.sum(theta_null >= theta) if upper else np.sum(theta_null <= theta)
    return float((1 + int(hits)) / (theta_null.size + 1))


def all_statistics(obs, surr, upper=True):
    """All four cohort statistics + diagnostics for one localization cell.

    Returns a dict: ``pooled_p, wilcoxon_p, stouffer_p, signconsist_p, frac_pos,
    n_patients, obs_med, z_med`` (``z_med`` = median per-patient z). ``upper=True``
    tests enrichment (carrier); ``upper=False`` tests depletion.
    """
    obs, surr = _as_cell(obs, surr)
    z = per_patient_z(obs, surr)
    return {
        "pooled_p": pooled_median_p(obs, surr, upper=upper),
        "wilcoxon_p": wilcoxon_p(obs, surr, alternative="greater" if upper else "less"),
        "stouffer_p": stouffer_p(obs, surr, upper=upper),
        "signconsist_p": sign_consistency_p(obs, surr, upper=upper),
        "frac_pos": frac_positive(obs, surr),
        "n_patients": int(np.isfinite(obs).sum()),
        "obs_med": float(np.nanmedian(obs)),
        "z_med": float(np.nanmedian(z[np.isfinite(z)])) if np.isfinite(z).any() else np.nan,
    }


def loo_worst(obs, surr, statistic="signconsist_p", upper=True):
    """Worst (max) cohort p over all single-patient drops — leave-one-patient-out.

    ``statistic`` ∈ the keys of :func:`all_statistics`. A localization that survives
    LOO has ``loo_worst < threshold``; a single-patient-driven one crosses it.
    """
    obs, surr = _as_cell(obs, surr)
    K = obs.shape[0]
    if K <= 2:
        return np.nan
    func = {"pooled_p": pooled_median_p, "wilcoxon_p": wilcoxon_p,
            "stouffer_p": stouffer_p, "signconsist_p": sign_consistency_p}[statistic]
    worst = 0.0
    for drop in range(K):
        keep = np.ones(K, bool); keep[drop] = False
        if statistic == "wilcoxon_p":
            p = func(obs[keep], surr[keep], alternative="greater" if upper else "less")
        else:
            p = func(obs[keep], surr[keep], upper=upper)
        if np.isfinite(p):
            worst = max(worst, p)
    return float(worst)


def bh_columns(df, p_cols, group_cols=None, suffix="_q"):
    """Add BH-FDR q-value columns for each p-column in ``df``.

    ``group_cols=None`` → **whole-grid** BH over the entire frame (the honest
    multiplicity family). ``group_cols=[...]`` → BH **within** each group (e.g. the
    per-cell reference over systems within a (target, band, scale)). Returns a new
    frame; NaN p-values pass through as NaN q.
    """
    import pandas as pd

    out = df.copy()
    for pc in p_cols:
        qc = pc.replace("_p", suffix) if pc.endswith("_p") else pc + suffix
        out[qc] = np.nan
        if group_cols is None:
            m = out[pc].notna()
            if m.any():
                out.loc[m, qc] = bh_fdr(out.loc[m, pc].to_numpy())
        else:
            for _, idx in out.groupby(group_cols).groups.items():
                sub = out.loc[idx]
                m = sub[pc].notna()
                if m.any():
                    out.loc[sub.index[m], qc] = bh_fdr(sub.loc[m, pc].to_numpy())
    return out
