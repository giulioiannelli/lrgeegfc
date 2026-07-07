#!/usr/bin/env python3
"""Audit 119 — PREMISE TEST for a path/routing SOZ marker: do the PATHS connecting SOZ
contacts have UNIQUE ROUTING properties (beyond raw connection magnitude / endpoint strength)?
This is the i-j / pair view of the propagator, not the node view. If SOZ-SOZ pairs are
distinctively routed AFTER controlling for endpoint strength, a seed-free routing module-hunt
is worth building; if not, the routing path is dead before we invest.

Why pairs, not nodes. Every per-NODE intrinsic marker reduced to strength + electrode depth
(audit_80/81/85/118). But e^{-tau L}_{ij} is a PATH integral i->j; routing properties live on
PAIRS, not nodes. Routing features (resistance distance = path redundancy; communicability =
weighted walk count; heat coupling) are normalized/topological and can be orthogonal to the
endpoint strengths. The hypothesis (literature: SOZ = inward/privately-routed module): SOZ-SOZ
pairs are MORE redundantly / privately routed than their endpoint strengths predict.

Design. Per patient/band compute pair routing features (resistance R_ij, communicability G_ij,
slow-heat K_ij), RESIDUALISE each on endpoint strength (regress on s_i+s_j and s_i*s_j) so the
test is routing-BEYOND-magnitude, then compare the residual over pair classes EE (SOZ-SOZ), EN
(SOZ-healthy), NN (healthy-healthy). Cohort: is EE distinctive vs NN and vs EN, paired across
patients (Wilcoxon)? Per band.

Critical preamble
=================
(1) Claim: SOZ-SOZ connecting paths have distinctive routing (e.g. lower resistance / higher
    communicability = more redundant private routing) BEYOND endpoint strength, cohort-wide.
(2) Null: EE routing residual = NN routing residual (SOZ pairs route like any pair of equal
    endpoint strength) -> no path signature.
(3) Strongest alternative: SOZ pairs just have strong endpoints (magnitude, hubness) ->
    killed by residualising every pair feature on s_i+s_j AND s_i*s_j; proximity -> |ImCoh|
    substrate is zero-lag immune (Nolte 2004), and resistance/communicability are global path
    properties, not adjacency.
(4) Reach: pair classes from the labels (this is a PREMISE test, labels allowed); residual is
    per-patient (own strength fit); cohort Wilcoxon paired over patients; per band; Pat_15
    has no NN-vs-EE issue (kept). Reported EE-NN and EE-EN separately.
(5) Falsification: if EE residual approx NN residual cohort-wide, SOZ paths are NOT uniquely
    routed and the seed-free routing-module idea is reported dead. No pre-registered gate.

Outputs (data/audit/epi_path_routing/)
    routing_pairclass_per_patient.csv ; routing_cohort.csv ; README.md
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
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI, N_TAU  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_path_routing"
ROUTING = ["resist", "comm", "heat_slow", "routeff"]   # routeff = comm/resistance (routing efficiency)
MIN_CLASS_PAIRS = 10


def pair_routing(W):
    """pairwise routing features (condensed triu order) + endpoint strengths."""
    N = W.shape[0]
    d = W.sum(1)
    L = np.diag(d) - W
    lam, U = eigh(L)
    lam = np.clip(lam, 0, None)
    Lp = pinv(L)
    diagLp = np.diag(Lp)
    R = diagLp[:, None] + diagLp[None, :] - 2 * Lp           # resistance distance (path redundancy)
    rhoW = max(abs(np.linalg.eigvalsh(W)).max(), 1e-9)
    G = expm(W / rhoW)                                       # communicability (weighted walk count)
    taus = np.geomspace(1.0 / lam[-1], 10.0 / lam[-1], N_TAU)
    Kslow = U @ (np.exp(-taus[-1] * lam)[:, None] * U.T)
    iu = np.triu_indices(N, 1)
    feats = {
        "resist": -R[iu],                                   # higher = MORE redundantly routed
        "comm": G[iu],
        "heat_slow": Kslow[iu],
        "routeff": G[iu] / (R[iu] + 1e-9),                  # walks per resistance (private routing)
    }
    si, sj = d[iu[0]], d[iu[1]]
    return feats, si, sj, iu


def _resid_on_strength(f, si, sj):
    """residualise a pair feature on endpoint strengths (sum and product) -> routing beyond magnitude."""
    X = np.column_stack([np.ones_like(si), si + sj, si * sj,
                         np.log1p(si) + np.log1p(sj)])
    f = np.asarray(f, float)
    ok = np.isfinite(f)
    if ok.sum() < 5:
        return f - np.nanmean(f)
    beta, *_ = np.linalg.lstsq(X[ok], f[ok], rcond=None)
    return f - X @ beta


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
    feats, si, sj, iu = pair_routing(W)
    ei, ej = epi[iu[0]], epi[iu[1]]
    cls_EE = ei & ej
    cls_NN = (~ei) & (~ej)
    cls_EN = ei ^ ej
    if cls_EE.sum() < MIN_CLASS_PAIRS or cls_NN.sum() < MIN_CLASS_PAIRS:
        return None
    rows = []
    for name in ROUTING:
        r = _resid_on_strength(feats[name], si, sj)
        # z-score residual within patient so cross-patient pooling is scale-free
        r = (r - np.nanmean(r)) / (np.nanstd(r) + 1e-12)
        mEE, mEN, mNN = np.nanmean(r[cls_EE]), np.nanmean(r[cls_EN]), np.nanmean(r[cls_NN])
        rows.append(dict(patient=pat, band=band, feature=name,
                         mean_EE=float(mEE), mean_EN=float(mEN), mean_NN=float(mNN),
                         EE_minus_NN=float(mEE - mNN), EE_minus_EN=float(mEE - mEN),
                         n_EE=int(cls_EE.sum()), n_NN=int(cls_NN.sum())))
    return pd.DataFrame(rows)


def cohort(df):
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(df.band)]:
        for feat in ROUTING:
            d = df[(df.band == band) & (df.feature == feat)]
            if len(d) < 4:
                continue
            for contrast in ("EE_minus_NN", "EE_minus_EN"):
                v = d[contrast].to_numpy(float); v = v[np.isfinite(v)]
                p = np.nan
                try:
                    p = float(wilcoxon(v, alternative="two-sided")[1]) if (v != 0).any() else 1.0
                except ValueError:
                    pass
                out.append(dict(band=band, feature=feat, contrast=contrast,
                                n=int(v.size), median=float(np.median(v)),
                                n_pos=int((v > 0).sum()), p=p))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+",
                    default=["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"])
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    frames = [per_patient_band(p, b) for b in args.bands for p in args.patients]
    df = pd.concat([f for f in frames if f is not None], ignore_index=True)
    df.to_csv(OUT / "routing_pairclass_per_patient.csv", index=False)
    coh = cohort(df)
    coh.to_csv(OUT / "routing_cohort.csv", index=False)

    print(f"[audit_119] SOZ path-routing premise test ({time.time()-t0:.1f}s)")
    print("  Q: do SOZ-SOZ paths route distinctively BEYOND endpoint strength?\n")
    print("  cohort EE−NN (routing residual, z; >0 = SOZ pairs more redundantly/privately routed):")
    print("  band        feature      median(EE−NN)  n_pos/n   p        | median(EE−EN)")
    for band in [b for b in es.ALL_BANDS if b in set(coh.band)]:
        for feat in ROUTING:
            a = coh[(coh.band == band) & (coh.feature == feat) & (coh.contrast == "EE_minus_NN")]
            b2 = coh[(coh.band == band) & (coh.feature == feat) & (coh.contrast == "EE_minus_EN")]
            if a.empty:
                continue
            r = a.iloc[0]; rb = b2.iloc[0] if not b2.empty else None
            star = "***" if r.p < 0.01 else ("**" if r.p < 0.05 else "")
            print(f"  {band:11s} {feat:11s}  {r['median']:+.3f}        {int(r.n_pos)}/{int(r.n)}    "
                  f"{r.p:.3f}{star}   | {rb['median']:+.3f}" if rb is not None else "")
    print(f"\n[audit_119] -> {OUT}")


if __name__ == "__main__":
    main()
