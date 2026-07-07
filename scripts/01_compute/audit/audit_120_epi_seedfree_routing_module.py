#!/usr/bin/env python3
"""Audit 120 — SEED-FREE routing-module marker: the premise (audit_119) holds — SOZ-SOZ paths
are privately routed beyond endpoint strength. Can we exploit that with NO labels? Build the
strength-residualised routing affinity between every pair, find the most privately-routed
MODULE spectrally, and score each node by it — then ask whether that label-free score recovers
SOZ. This is the actual fold-into-a-seed-free-marker test.

Why this can work where node-intrinsic failed (audit_118): the signal is RELATIONAL (which
nodes are privately routed TO EACH OTHER), carried by the leading community of the routing-
residual matrix, not by any per-node scalar. Why it might still fail (the audit_94 wall):
multiple physiological modules may be privately routed; the SOZ module may not be the dominant
eigenvector. We MEASURE, label-free, per patient and cohort.

Marker (label-free, per patient, per band):
  M = routeff_ij residualised on endpoint strength (s_i+s_j, s_i*s_j) -> symmetric pair matrix.
  Score variants:
    eig1      = |leading eigenvector of M| (the dominant routing module)
    eigK_part = sum_{m=1..K} eigvec_m^2 weighted by eigenvalue (participation in top routing modules)
    rcore     = mean of node's top-q routing-residual affinities (local routing coreness)
  Rank all nodes by score -> AUC vs SOZ (NO labels used to build the score). Baselines: strength,
  and a phase-randomised / label-shuffle module null.

Critical preamble
=================
(1) Claim: a label-free routing-module score ranks SOZ above healthy (AUC>0.5) cohort-wide,
    beating node strength, in the bands where the premise holds (delta/beta/low_gamma).
(2) Null: label-free routing AUC = 0.5 (the SOZ module is not the recoverable routing module);
    = strength (routing score is just hubness re-expressed).
(3) Strongest alternative: hubness -> M is strength-residualised by construction + explicit
    strength baseline; the dominant module is the highest-strength block -> the residual removes
    the strength gradient, and we test eig1 AND deeper modes.
(4) Reach: NO labels enter the score; per-patient unsupervised; cohort sign test; per band;
    hub patients (Pat_10/15) kept. Multiple score variants reported (multiple-comparison noted).
(5) Falsification: if no label-free routing score beats 0.5 / beats strength cohort-wide, the
    premise does NOT fold into a seed-free marker (audit_94 wall confirmed) and we report it.

Outputs (data/audit/epi_seedfree_routing/)
    seedfree_routing_per_patient.csv ; seedfree_routing_cohort.csv ; README.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
from numpy.linalg import eigh, pinv
from scipy.linalg import expm
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.utils.io.patient import build_epi_masks

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc, strength_preserving_shuffle  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI  # type: ignore
from audit_101_epi_marker_library import _conc  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_seedfree_routing"
SCORES = ["eig1", "eigK_part", "rcore", "strength"]
KMODES = 5
QCORE = 0.10
NULL_SEED = 20260619


def _eigK_part(M):
    w, V = eigh(M)
    order = np.argsort(-w)
    V = V[:, order]; w = w[order]
    return (V[:, :KMODES] ** 2 * np.clip(w[:KMODES], 0, None)).sum(1)


def _auc(score, epi):
    s = np.asarray(score, float)
    c, t = _conc(s[epi], s[~epi])
    return c / t if t else np.nan


def routing_residual_matrix(W):
    N = W.shape[0]
    d = W.sum(1)
    L = np.diag(d) - W
    Lp = pinv(L)
    diagLp = np.diag(Lp)
    R = diagLp[:, None] + diagLp[None, :] - 2 * Lp
    rhoW = max(abs(np.linalg.eigvalsh(W)).max(), 1e-9)
    G = expm(W / rhoW)
    routeff = G / (R + 1e-9)
    np.fill_diagonal(routeff, 0.0)
    iu = np.triu_indices(N, 1)
    si, sj = d[iu[0]], d[iu[1]]
    X = np.column_stack([np.ones_like(si), si + sj, si * sj,
                         np.log1p(si) + np.log1p(sj)])
    f = routeff[iu]
    beta, *_ = np.linalg.lstsq(X, f, rcond=None)
    res = f - X @ beta                                   # routing beyond endpoint strength
    M = np.zeros((N, N))
    M[iu] = res; M = M + M.T
    return M, d


def per_patient_band(pat, band):
    try:
        W = load_phase_fc(pat, "rest_post", band)
    except Exception:
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    if epi.sum() < MIN_EPI or (~epi).sum() < 5:
        return None
    M, d = routing_residual_matrix(W)
    # symmetric eighowever: M may be indefinite; use eigh
    w, V = eigh(M)
    order = np.argsort(-w)                                # descending eigenvalue
    V = V[:, order]; w = w[order]
    v1 = V[:, 0]
    if np.sum(v1) < 0:
        v1 = -v1
    eig1 = np.abs(v1)
    eigK_part = (V[:, :KMODES] ** 2 * np.clip(w[:KMODES], 0, None)).sum(1)
    # routing coreness: mean of each node's top-q routing-residual affinities
    q = max(2, int(QCORE * N))
    Ms = np.sort(M, axis=1)[:, ::-1]
    rcore = Ms[:, :q].mean(1)
    scores = {"eig1": eig1, "eigK_part": eigK_part, "rcore": rcore, "strength": d}
    rows = []
    for name, sc in scores.items():
        rows.append(dict(patient=pat, band=band, score=name, auc=float(_auc(sc, epi)),
                         n_soz=int(epi.sum())))
    return pd.DataFrame(rows)


def cohort(df):
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(df.band)]:
        for sc in SCORES:
            d = df[(df.band == band) & (df.score == sc)]
            v = d.auc.to_numpy(float); v = v[np.isfinite(v)]
            if v.size < 4:
                continue
            p = np.nan
            try:
                p = float(wilcoxon(v - 0.5, alternative="greater")[1])
            except ValueError:
                pass
            out.append(dict(band=band, score=sc, n=int(v.size), median_auc=float(np.median(v)),
                            mean_auc=float(np.mean(v)), n_above_half=int((v > 0.5).sum()), p=p))
    return pd.DataFrame(out)


def matched_strength_null(pat, band, n_null):
    """observed eigK_part AUC vs strength-preserving-surrogate eigK_part AUC (same SOZ labels)."""
    try:
        W = load_phase_fc(pat, "rest_post", band)
    except Exception:
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    if epi.sum() < MIN_EPI or (~epi).sum() < 5:
        return None
    M, _ = routing_residual_matrix(W)
    obs = _auc(_eigK_part(M), epi)
    rng = np.random.default_rng(NULL_SEED + hash(pat + band) % 9999)
    surr = []
    for _ in range(n_null):
        Ws = strength_preserving_shuffle(W, 10 * N, rng, w_max=1.0)
        Ms, _ = routing_residual_matrix(Ws)
        surr.append(_auc(_eigK_part(Ms), epi))
    surr = np.array([s for s in surr if np.isfinite(s)])
    return dict(patient=pat, band=band, auc_obs=float(obs),
                surr_mean=float(np.mean(surr)), surr_p95=float(np.percentile(surr, 95)),
                p_null=float(np.mean(surr >= obs)), n_null=int(surr.size))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=["delta", "beta", "low_gamma"])
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--null", type=int, default=0, help="matched-strength surrogates per cell")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    frames = [per_patient_band(p, b) for b in args.bands for p in args.patients]
    df = pd.concat([f for f in frames if f is not None], ignore_index=True)
    df.to_csv(OUT / "seedfree_routing_per_patient.csv", index=False)
    coh = cohort(df)
    coh.to_csv(OUT / "seedfree_routing_cohort.csv", index=False)

    if args.null > 0:
        nrows = [matched_strength_null(p, b, args.null) for b in args.bands for p in args.patients]
        nd = pd.DataFrame([r for r in nrows if r is not None])
        nd.to_csv(OUT / "seedfree_routing_matched_strength_null.csv", index=False)
        print(f"\n[audit_120] MATCHED-STRENGTH NULL on eigK_part ({args.null} surrogates/cell):")
        print("  band        median AUC_obs  median AUC_surr  cohort p_null  n>obs_excess")
        for band in args.bands:
            d = nd[nd.band == band]
            if d.empty:
                continue
            excess = (d.auc_obs - d.surr_mean)
            try:
                pc = float(wilcoxon(excess, alternative="greater")[1]) if (excess != 0).any() else 1.0
            except ValueError:
                pc = np.nan
            print(f"  {band:11s} {d.auc_obs.median():.3f}           {d.surr_mean.median():.3f}"
                  f"            {pc:.3f}         {int((excess>0).sum())}/{len(d)}")

    print(f"[audit_120] SEED-FREE routing-module marker, LABEL-FREE ({time.time()-t0:.1f}s)")
    print("  Q: does the routing module recover SOZ with NO labels, beating strength?\n")
    print("  band        score        median AUC   n>0.5   p(>0.5)")
    for band in [b for b in es.ALL_BANDS if b in set(coh.band)]:
        for sc in SCORES:
            r = coh[(coh.band == band) & (coh.score == sc)]
            if r.empty:
                continue
            r = r.iloc[0]
            star = "***" if r.p < 0.01 else ("**" if r.p < 0.05 else "")
            print(f"  {band:11s} {sc:11s}  {r.median_auc:.3f}        {int(r.n_above_half)}/{int(r.n)}    {r.p:.3f}{star}")
        print()
    print("  per-patient eig1 / rcore AUC (delta):")
    d = df[(df.band == 'delta') & (df.score.isin(['eig1', 'rcore', 'strength']))]
    piv = d.pivot(index='patient', columns='score', values='auc')
    print(piv.round(3).to_string())
    print(f"\n[audit_120] -> {OUT}")


if __name__ == "__main__":
    main()
