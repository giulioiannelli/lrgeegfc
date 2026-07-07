#!/usr/bin/env python3
"""Audit 118 — SEED-FREE intrinsic SOZ marker: can multiscale/multiband Laplacian-propagator
properties flag SOZ with NO seeds, beating the strength + along-shaft-depth baselines that
killed every prior intrinsic marker (audit_80/81/85)? The go/no-go for a seedless detector.

Motivation. The seed-based detector (audit_117) works because of a node's RELATION to known
SOZ; remove the seeds and the anchor is gone. Two prior seed-free attempts FAILED: label-free
within-patient selection at chance (audit_94), cross-patient intrinsic node classifier at
chance (audit_80/81; everything reduced to hubness + electrode depth). So a naive "drop the
seed features" is expected at chance. The only way seed-free can work is a NEW discriminative
property. We test the literature's strongest candidate — the epileptogenic zone as an
ISOLATED / inward-connected module (high internal, low external connectivity) — read in
propagator language as SEED-FREE intrinsic features, multiband-fused:

  ISOLATION (Path 1): row-entropy of e^{-tau L} (low = concentrated diffusion = isolated);
      heat self-return rho_ii(tau); heat-row IPR sum_j p_ij^2 (high = concentrated).
  EIGENMODE LOCALIZATION (Path 2): per-node participation in localized slow modes
      sum_{k=1..K} U[i,k]^4 (IPR-weighted).
  MULTISCALE PERSISTENCE (Path 3-lite): mean cophenetic distance (late merge = isolated).

scored against the MANDATORY confound baselines: node strength (hubness) and normalized
along-shaft contact depth (the audit_85/§10 killer). Within-patient z-score, LOPO logistic.

Critical preamble
=================
(1) Claim: an intrinsic (seed-free) multiscale/multiband Laplacian feature set ranks SOZ
    above chance cross-patient (LOPO) AND beats strength-only, depth-only, strength+depth.
(2) Null: seed-free intrinsic AUC = 0.5; OR = the strength+depth baseline (i.e. the new
    features add nothing beyond the two known confounds); OR label-shuffle ~ 0.5.
(3) Strongest alternative: hubness (strength is an explicit feature + ablation rung) and
    electrode depth (along-shaft index is an explicit feature + ablation rung). EVERY prior
    intrinsic marker reduced to these two; the new features must beat strength+depth TOGETHER.
(4) Reach: features computed per band, z-scored WITHIN patient (kills FC-scale drift), LOPO
    (held-out patient never trains); ablations isolate intrinsic-new vs strength vs depth vs
    strength+depth; label-shuffle null permutes SOZ within patient; hub patients kept.
(5) Falsification: if intrinsic-new does NOT beat strength+depth (or is at chance / = shuffle),
    seed-free intrinsic marking is NOT achievable on this data and we report it dead, and the
    seed-based detector (audit_117) remains the ceiling. No pre-registered acceptance gate.

Outputs (data/audit/epi_seedfree_marker/)
    seedfree_lopo_per_patient.csv ; seedfree_ablation.csv ; feature_importance.csv ; README.md
"""
from __future__ import annotations

import argparse
import re
import sys
import time

import numpy as np
import pandas as pd
from numpy.linalg import eigh
from scipy.spatial.distance import squareform
from sklearn.linear_model import LogisticRegression

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.utils.io.patient import build_epi_masks

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc, lrg_ultrametric_condensed  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI, N_TAU  # type: ignore
from audit_101_epi_marker_library import _conc  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_seedfree_marker"
ISO_FEATS = ["ent_t0", "ent_t2", "ent_t5", "heatret_t5", "iprrow_t2",
             "iprslow_K5", "iprslow_K10", "cophdepth"]
SEED = 20260619


def _depth_along_shaft(channels, probes):
    """normalized along-shaft contact index in [0,1] per probe (the audit_85 depth baseline)."""
    chans = np.asarray(channels, object); pr = np.asarray(probes, object)
    num = np.array([int(m.group()) if (m := re.search(r"\d+", str(c))) else -1 for c in chans])
    depth = np.full(len(chans), np.nan)
    for p in set(pr.tolist()):
        idx = np.where(pr == p)[0]
        v = num[idx].astype(float); v[v < 0] = np.nan
        if np.isfinite(v).sum() >= 2 and np.nanmax(v) > np.nanmin(v):
            depth[idx] = (v - np.nanmin(v)) / (np.nanmax(v) - np.nanmin(v))
        else:
            depth[idx] = 0.5
    return depth


def _intrinsic_features(W):
    """SEED-FREE per-node intrinsic Laplacian-propagator features."""
    N = W.shape[0]
    d = W.sum(1)
    L = np.diag(d) - W
    lam, U = eigh(L)
    lam = np.clip(lam, 0, None)
    taus = np.geomspace(1.0 / lam[-1], 10.0 / lam[-1], N_TAU)
    feat = {}
    for tag, ti in (("t0", 0), ("t2", 2), ("t5", 5)):
        K = U @ (np.exp(-taus[ti] * lam)[:, None] * U.T)
        K = np.clip(K, 0, None)
        rs = K.sum(1, keepdims=True); rs[rs < 1e-300] = 1e-300
        P = K / rs                                              # row-normalized diffusion
        with np.errstate(divide="ignore", invalid="ignore"):
            ent = -(P * np.log(np.clip(P, 1e-300, None))).sum(1)  # diffusion spread (low=isolated)
        feat[f"ent_{tag}"] = ent
        if tag == "t2":
            feat["iprrow_t2"] = (P ** 2).sum(1)                 # row concentration (high=isolated)
        if tag == "t5":
            feat["heatret_t5"] = np.diag(K).copy()              # slow self-return
    for Kk in (5, 10):
        feat[f"iprslow_K{Kk}"] = (U[:, 1:1 + Kk] ** 4).sum(1)   # localized slow-mode participation
    Dc = squareform(lrg_ultrametric_condensed(W))
    feat["cophdepth"] = Dc.mean(1)                              # mean cophenetic distance (late merge)
    return feat, d


def patient_matrix(pat, bands):
    pm = build_epi_masks(pat)
    probes = np.asarray(pm.probes, object)
    blocks, names = [], []
    N0 = None
    strength_b0 = None
    for band in bands:
        try:
            W = load_phase_fc(pat, "rest_post", band)
        except Exception:
            return None
        N = W.shape[0]
        if len(pm.channels) != N:
            return None
        if N0 is None:
            N0 = N
        elif N != N0:
            return None
        feat, d = _intrinsic_features(W)
        for f in ISO_FEATS:
            blocks.append(feat[f]); names.append(f"{band}:{f}")
        blocks.append(d); names.append(f"{band}:strength")
        if strength_b0 is None:
            strength_b0 = d
    depth = _depth_along_shaft(pm.channels, probes)
    blocks.append(depth); names.append("geom:depth")
    X = np.column_stack(blocks).astype(float)
    # within-patient z-score
    mu = np.nanmean(X, 0); sd = np.nanstd(X, 0); sd[sd < 1e-12] = 1.0
    Xz = (X - mu) / sd
    Xz[~np.isfinite(Xz)] = 0.0
    epi = np.asarray(pm.epi_mask, bool)
    y = epi.astype(int)
    if epi.sum() < MIN_EPI or (~epi).sum() < 5:
        return None
    return Xz, y, np.array([pat] * N0), names


def _auc(p, yv):
    c, t = _conc(np.asarray(p, float)[yv == 1], np.asarray(p, float)[yv == 0])
    return c / t if t else np.nan


def _prec_at(p, yv, k):
    o = np.argsort(-p)[:k]
    return float(yv[o].sum() / k)


def lopo(X, y, pid, pats, cols=None, C=1.0):
    Xu = X[:, cols] if cols is not None else X
    out = {}
    for pt in pats:
        tr, te = pid != pt, pid == pt
        if te.sum() == 0 or y[tr].sum() == 0 or y[te].sum() == 0:
            continue
        clf = LogisticRegression(max_iter=3000, C=C).fit(Xu[tr], y[tr])
        out[pt] = (clf.predict_proba(Xu[te])[:, 1], y[te])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+",
                    default=["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"])
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--C", type=float, default=1.0)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    t0 = time.time()

    mats = {}
    names = None
    for pat in args.patients:
        r = patient_matrix(pat, args.bands)
        if r is None:
            continue
        Xz, y, pid, names = r
        mats[pat] = (Xz, y, pid)
    pats = list(mats.keys())
    X = np.vstack([mats[p][0] for p in pats])
    y = np.concatenate([mats[p][1] for p in pats])
    pid = np.concatenate([mats[p][2] for p in pats])
    idx = {n: i for i, n in enumerate(names)}

    iso_cols = [i for n, i in idx.items() if any(n.endswith(":" + f) for f in ISO_FEATS)]
    str_cols = [i for n, i in idx.items() if n.endswith(":strength")]
    depth_col = [idx["geom:depth"]]
    groups = {
        "intrinsic_new (iso+ipr+coph)": iso_cols,
        "strength_only": str_cols,
        "depth_only": depth_col,
        "strength+depth": str_cols + depth_col,
        "intrinsic+strength+depth (all)": None,
    }

    rows = []
    preds_full = lopo(X, y, pid, pats, cols=None, C=args.C)
    for g, cols in groups.items():
        pr = lopo(X, y, pid, pats, cols=cols, C=args.C)
        aucs = [_auc(p, yt) for p, yt in pr.values()]
        p5 = [_prec_at(p, yt, 5) for p, yt in pr.values()]
        rows.append(dict(group=g, n_feat=(len(cols) if cols else X.shape[1]),
                         mean_auc=float(np.nanmean(aucs)), median_auc=float(np.nanmedian(aucs)),
                         n_above_half=int(np.sum(np.array(aucs) > 0.5)), n=len(aucs),
                         mean_prec5=float(np.nanmean(p5))))
    dab = pd.DataFrame(rows)
    dab.to_csv(OUT / "seedfree_ablation.csv", index=False)

    # label-shuffle null on intrinsic_new
    yshuf = y.copy()
    for pt in pats:
        ii = np.where(pid == pt)[0]
        yshuf[ii] = rng.permutation(yshuf[ii])
    pr_sh = lopo(X, yshuf, pid, pats, cols=iso_cols, C=args.C)
    auc_sh = float(np.nanmean([_auc(p, yt) for p, yt in pr_sh.values()]))

    perpat = [dict(patient=pt, auc=_auc(p, yt), prec5=_prec_at(p, yt, 5),
                   n_soz=int(yt.sum())) for pt, (p, yt) in preds_full.items()]
    dpp = pd.DataFrame(perpat)
    dpp.to_csv(OUT / "seedfree_lopo_per_patient.csv", index=False)

    clf = LogisticRegression(max_iter=3000, C=args.C).fit(X[:, iso_cols], y)
    fi = pd.DataFrame({"feature": [names[i] for i in iso_cols], "coef": clf.coef_[0],
                       "abs_coef": np.abs(clf.coef_[0])}).sort_values("abs_coef", ascending=False)
    fi.to_csv(OUT / "feature_importance.csv", index=False)

    print(f"[audit_118] SEED-FREE intrinsic marker, LOPO ({len(pats)} pts, {time.time()-t0:.1f}s)")
    print(f"  bands={args.bands}\n")
    print("  ABLATION (mean LOPO AUC, the go/no-go is intrinsic_new vs strength+depth):")
    for _, r in dab.iterrows():
        print(f"    {r.group:32s} AUC {r.mean_auc:.3f} (med {r.median_auc:.3f}, "
              f"{r.n_above_half}/{r.n}>0.5, prec@5 {r.mean_prec5*100:3.0f}%, {int(r.n_feat)}f)")
    print(f"\n  label-shuffle null (intrinsic_new): {auc_sh:.3f}")
    print("\n  per-patient (intrinsic+strength+depth full model):")
    for _, r in dpp.sort_values("auc", ascending=False).iterrows():
        print(f"    {r.patient}: AUC {r.auc:.3f}  prec@5 {r.prec5*100:3.0f}%  (#SOZ {int(r.n_soz)})")
    print("\n  top intrinsic features (|coef|):")
    for _, r in fi.head(8).iterrows():
        print(f"    {r.feature:22s} {r.coef:+.3f}")
    print(f"\n[audit_118] -> {OUT}")


if __name__ == "__main__":
    main()
