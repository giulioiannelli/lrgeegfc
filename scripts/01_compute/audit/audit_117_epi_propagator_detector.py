#!/usr/bin/env python3
"""Audit 117 — turn the marker into a DETECTOR: a cross-patient, calibrated per-node
probability P(SOZ | node, given k known seeds), built ONLY on propagator-derived features
(no coordinates, no proximity). Answers exactly "the probability of marking a previously
unmarked node as SOZ thanks to a few known contacts."

Method. From the LRG heat-kernel propagator e^{-tau L} of the |ImCoh| graph (volume-conduction
immune -> proximity cannot enter the connectivity) and its derivatives, build a rich per-node
feature vector RELATIVE TO the seed set S (the few known SOZ):
    multiscale seed-affinity (heat fast/mid/slow), personalized PageRank, Katz, communicability,
    -commute/resistance distance to seeds, diffusion-distance; the RELATIONAL contrast
    segregation = aff(.,S) - aff(.,non-SOZ) and diffusion-share to S; plus seed-free intrinsics
    (heat-return probability, slow-mode participation) and node strength (so hub patients are
    handled). Features are averaged over random k-seed draws (leave-self-out) and z-scored
    WITHIN patient (kills cross-patient FC-scale drift), then a LOGISTIC model is trained on
    N-1 patients and predicts the held-out patient's nodes -> calibrated P(SOZ). This is the
    deployable detector: new patient + k seeds -> P(SOZ) for every other contact.

Critical preamble
=================
(1) Claim: a logistic model on propagator features, trained leave-one-PATIENT-out, predicts
    held-out SOZ above chance with usable precision and CALIBRATED probability, beating the
    single-affinity marker, strength alone, and a label-shuffle null.
(2) Null: LOPO AUC = 0.5 / precision = prevalence / P(SOZ) uninformative (Brier = prevalence
    variance); label-shuffle (permute SOZ labels within patient) -> AUC 0.5.
(3) Strongest alternative: hubness (strength is one feature + a strength-only ablation rung);
    overfitting (LOPO = the held-out patient never trains the model; L2-regularised logistic;
    ~14 features on ~1000 pooled nodes); proximity (substrate is |ImCoh|, zero-lag immune).
(4) Reach: features averaged over k-seed draws (deploy-faithful), z-scored within patient;
    held-out patient fully isolated; hub patients (Pat_10/15) kept and reported per-patient;
    ablations isolate each feature group; multi-band stacking optional.
(5) Falsification: if LOPO AUC ~ 0.5 or <= the single-affinity marker, the rich detector adds
    nothing and we report the simple marker as the ceiling. P(SOZ) reported with calibration,
    no cherry-picked operating point. No pre-registered acceptance gate.

Outputs (data/audit/epi_propagator_detector/)
    detector_lopo_per_patient.csv ; detector_ablation.csv ; feature_importance.csv ; README.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.utils.io.regions import load_channel_regions

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI  # type: ignore
from audit_101_epi_marker_library import build_operators, _conc  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_propagator_detector"
OP_FEATS = ["heat_t1", "heat_t3", "heat_t5", "ppr_a85", "katz", "comm", "negres", "diffdist_t2"]
K_SEED = 3
N_DRAW = 60
SEED = 20260619
BUDGETS = (5, 10)


def _auc(p, y):
    c, t = _conc(np.asarray(p, float)[y == 1], np.asarray(p, float)[y == 0])
    return c / t if t else np.nan


def patient_features(pat, band, k, ndraw, rng):
    """Per-node propagator feature matrix (z-scored within patient), label y=SOZ."""
    try:
        W = load_phase_fc(pat, "rest_post", band)
    except Exception:
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    epi_idx = np.where(epi)[0]
    if epi_idx.size < max(MIN_EPI, k + 1) or (~epi).sum() < 5:
        return None
    ops, intrinsic, strength = build_operators(W)
    Kslow = ops["heat_t5"]
    soz_set = set(epi_idx.tolist())

    nf = len(OP_FEATS) + 2                       # + segregation + diffusion-share
    fsum = np.zeros((N, nf)); fcnt = np.zeros(N)
    for _ in range(ndraw):
        seeds = rng.choice(epi_idx, size=k, replace=False)
        sset = set(seeds.tolist())
        scored = np.array([i for i in range(N) if i not in sset])
        nonseed = np.array([j for j in range(N) if j not in sset])    # DEPLOY-LEGAL reference
        cols = [ops[nm][scored][:, seeds].mean(axis=1) for nm in OP_FEATS]
        cols.append(Kslow[scored][:, seeds].mean(axis=1)
                    - Kslow[scored][:, nonseed].mean(axis=1))          # segregation (seeds vs non-seeds)
        cols.append(Kslow[scored][:, seeds].sum(axis=1)
                    / (Kslow[scored].sum(axis=1) + 1e-12))             # diffusion-share
        fsum[scored] += np.column_stack(cols); fcnt[scored] += 1
    ok = fcnt > 0
    feat = np.full((N, nf), np.nan)
    feat[ok] = fsum[ok] / fcnt[ok][:, None]
    feat = np.column_stack([feat, intrinsic["heatdiag_t5"], intrinsic["slowpart_K5"], strength])
    names = OP_FEATS + ["segregation", "diff_share", "heat_return", "slow_part", "strength"]
    # z-score within patient
    mu = np.nanmean(feat, axis=0); sd = np.nanstd(feat, axis=0); sd[sd < 1e-12] = 1.0
    fz = (feat - mu) / sd
    fz[~np.isfinite(fz)] = 0.0
    y = np.array([1 if i in soz_set else 0 for i in range(N)], int)
    nidx = np.where(ok)[0]
    return fz[ok], y[ok], np.array([pat] * int(ok.sum())), names, nidx


def _make_clf(model, C):
    if model == "gbm":
        from sklearn.ensemble import HistGradientBoostingClassifier
        return HistGradientBoostingClassifier(max_depth=3, max_iter=150,
                                              learning_rate=0.05, l2_regularization=1.0)
    return LogisticRegression(max_iter=2000, C=C)


def lopo(X, y, pid, patients, cols=None, C=1.0, model="logit"):
    """leave-one-patient-out classifier; returns per-patient predictions."""
    Xu = X[:, cols] if cols is not None else X
    out = {}
    for pt in patients:
        tr, te = pid != pt, pid == pt
        if te.sum() == 0 or y[tr].sum() == 0 or y[te].sum() == 0:
            continue
        clf = _make_clf(model, C)
        clf.fit(Xu[tr], y[tr])
        out[pt] = (clf.predict_proba(Xu[te])[:, 1], y[te])
    return out


def _metrics(p, y):
    order = np.argsort(-p)
    yo = y[order]
    T = int(y.sum())
    res = {"auc": _auc(p, y), "brier": float(np.mean((p - y) ** 2)), "n_soz": T,
           "p_mean_soz": float(p[y == 1].mean()), "p_mean_healthy": float(p[y == 0].mean())}
    for N in BUDGETS:
        n = min(N, yo.size)
        res[f"prec{N}"] = float(yo[:n].sum() / N)
        res[f"rec{N}"] = float(yo[:n].sum() / T) if T else np.nan
    return res


def main():
    ap = argparse.ArgumentParser()
    # canonical config: 6-band multi-frequency fusion, L2-logistic (AUC 0.81, prec@5 60%,
    # interpretable + regularization-stable -> the reliable/sellable tool; GBM gives AUC 0.86
    # but is overfit-risky on n=10 and does not improve prec@5).
    ap.add_argument("--bands", nargs="+",
                    default=["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"])
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--ndraw", type=int, default=N_DRAW)
    ap.add_argument("--k", type=int, default=K_SEED)
    ap.add_argument("--C", type=float, default=1.0)
    ap.add_argument("--model", choices=["logit", "gbm"], default="logit")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    t0 = time.time()

    # build per-patient features; stack across bands (feature fusion) on aligned nodes
    perpat = {}
    nidx_of = {}
    for pat in args.patients:
        blocks = []
        ok = True
        nidx0 = None
        for band in args.bands:
            r = patient_features(pat, band, args.k, args.ndraw, rng)
            if r is None:
                ok = False; break
            fz, y, pid, nm, nidx = r
            blocks.append(fz)
            if nidx0 is None:
                nidx0 = nidx
        if not ok:
            continue
        perpat[pat] = (np.column_stack(blocks), y, pid)
        nidx_of[pat] = nidx0
    names = [f"{b}:{c}" for b in args.bands
             for c in (OP_FEATS + ["segregation", "diff_share", "heat_return", "slow_part", "strength"])]
    pats = list(perpat.keys())
    X = np.vstack([perpat[p][0] for p in pats])
    y = np.concatenate([perpat[p][1] for p in pats])
    pid = np.concatenate([perpat[p][2] for p in pats])

    # ---- full-model LOPO ----
    preds = lopo(X, y, pid, pats, C=args.C, model=args.model)
    rows = []
    for pt, (p, yt) in preds.items():
        rows.append({"patient": pt, **_metrics(p, yt)})
    dpp = pd.DataFrame(rows)
    dpp.to_csv(OUT / "detector_lopo_per_patient.csv", index=False)

    # ---- label-shuffle null (permute SOZ labels within each patient) ----
    yshuf = y.copy()
    for pt in pats:
        idx = np.where(pid == pt)[0]
        yshuf[idx] = rng.permutation(yshuf[idx])
    preds_sh = lopo(X, yshuf, pid, pats)
    auc_sh = np.nanmean([_auc(p, yt) for p, yt in preds_sh.values()])

    # ---- ablations ----
    nfb = len(OP_FEATS) + 5
    groups = {"full": None}
    if len(args.bands) == 1:
        idx = {n: i for i, n in enumerate(names)}
        b = args.bands[0]
        groups["affinity_only"] = [idx[f"{b}:{c}"] for c in OP_FEATS]
        groups["slow_affinity_only"] = [idx[f"{b}:heat_t5"]]
        groups["strength_only"] = [idx[f"{b}:strength"]]
        groups["relational_only"] = [idx[f"{b}:{c}"] for c in ("segregation", "diff_share")]
        groups["no_strength"] = [i for n, i in idx.items() if not n.endswith(":strength")]
    abl = []
    for g, cols in groups.items():
        pr = lopo(X, y, pid, pats, cols=cols)
        aucs = [_auc(p, yt) for p, yt in pr.values()]
        p5 = [_metrics(p, yt)["prec5"] for p, yt in pr.values()]
        abl.append({"group": g, "n_feat": (len(cols) if cols else X.shape[1]),
                    "mean_auc": float(np.nanmean(aucs)), "median_auc": float(np.nanmedian(aucs)),
                    "n_above_half": int(np.sum(np.array(aucs) > 0.5)), "n": len(aucs),
                    "mean_prec5": float(np.nanmean(p5))})
    dab = pd.DataFrame(abl)
    dab.to_csv(OUT / "detector_ablation.csv", index=False)

    # ---- feature importance (|coef| of an all-data logistic, standardized features) ----
    clf = LogisticRegression(max_iter=2000, C=1.0).fit(X, y)
    fi = pd.DataFrame({"feature": names, "coef": clf.coef_[0],
                       "abs_coef": np.abs(clf.coef_[0])}).sort_values("abs_coef", ascending=False)
    fi.to_csv(OUT / "feature_importance.csv", index=False)

    # ---- node-level P(SOZ) with channel labels + regions (deploy view) ----
    # Retrain LOPO capturing each held-out patient's per-node prediction, align to nidx,
    # attach the clinical channel label + Desikan-Killiany region. SOZ-labelled rows = the
    # recoverable targets; NON-SOZ rows ranked high = occult candidates (HYPOTHESES only).
    node_rows = []
    for pt in pats:
        tr, te = pid != pt, pid == pt
        if te.sum() == 0 or y[tr].sum() == 0:
            continue
        cl = LogisticRegression(max_iter=2000, C=args.C).fit(X[tr], y[tr])
        p_te = cl.predict_proba(X[te])[:, 1]
        nidx = nidx_of[pt]
        pm = build_epi_masks(pt)
        chans = np.asarray(pm.channels, object)
        probes = np.asarray(pm.probes, object)
        try:
            reg = load_channel_regions(pt)
            regions = reg["region"].to_numpy(object); hemis = reg["hemisphere"].to_numpy(object)
        except Exception:
            regions = np.array(["?"] * len(chans), object); hemis = np.array(["?"] * len(chans), object)
        yte = y[te]
        for j, ni in enumerate(nidx):
            node_rows.append(dict(
                patient=pt, node=int(ni),
                label=str(chans[ni]) if ni < len(chans) else "?",
                probe=str(probes[ni]) if ni < len(probes) else "?",
                region=str(regions[ni]) if ni < len(regions) else "?",
                hemisphere=str(hemis[ni]) if ni < len(hemis) else "?",
                is_soz=int(yte[j]), p_soz=float(p_te[j])))
    nodes = pd.DataFrame(node_rows)
    nodes.to_csv(OUT / "detector_node_predictions.csv", index=False)
    # occult candidates = highest-P NON-SOZ contacts per patient (top 3)
    cand = (nodes[nodes.is_soz == 0].sort_values("p_soz", ascending=False)
            .groupby("patient").head(3).sort_values(["patient", "p_soz"], ascending=[True, False]))
    cand.to_csv(OUT / "occult_candidates.csv", index=False)

    print(f"[audit_117] propagator detector — LOPO calibrated P(SOZ) "
          f"(bands={args.bands}, k={args.k}, {time.time()-t0:.1f}s)\n")
    print(f"  FULL model: mean AUC {dpp.auc.mean():.3f}  median {dpp.auc.median():.3f}  "
          f"{int((dpp.auc>0.5).sum())}/{len(dpp)}>0.5  | prec@5 {dpp.prec5.mean()*100:.0f}%  "
          f"prec@10 {dpp.prec10.mean()*100:.0f}%  rec@10 {dpp.rec10.mean()*100:.0f}%")
    print(f"  calibration: P(SOZ|true SOZ) {dpp.p_mean_soz.mean():.3f} vs "
          f"P(SOZ|healthy) {dpp.p_mean_healthy.mean():.3f}  Brier {dpp.brier.mean():.4f}")
    print(f"  label-shuffle null AUC: {auc_sh:.3f}\n")
    print("  per-patient LOPO AUC:")
    for _, r in dpp.sort_values("auc", ascending=False).iterrows():
        print(f"    {r.patient}: AUC {r.auc:.3f}  prec@5 {r.prec5*100:3.0f}%  "
              f"P(SOZ) top vs base {r.p_mean_soz:.3f}/{r.p_mean_healthy:.3f}  (#SOZ {int(r.n_soz)})")
    print("\n  ablations (mean LOPO AUC):")
    for _, r in dab.iterrows():
        print(f"    {r.group:18s} {r.mean_auc:.3f}  ({r.n_above_half}/{r.n}>0.5, "
              f"prec@5 {r.mean_prec5*100:3.0f}%, {int(r.n_feat)} feat)")
    print("\n  top features (|coef|):")
    for _, r in fi.head(8).iterrows():
        print(f"    {r.feature:26s} {r.coef:+.3f}")
    print("\n  OCCULT CANDIDATES — highest-P(SOZ) NON-marked contacts per patient (HYPOTHESES):")
    print("    patient   contact      region                         hemi   P(SOZ)")
    for _, r in cand.iterrows():
        print(f"    {r.patient}  {str(r.label):11s}  {str(r.region):30s} {str(r.hemisphere):4s}   {r.p_soz:.3f}")
    print(f"\n[audit_117] -> {OUT}")


if __name__ == "__main__":
    main()
