#!/usr/bin/env python3
"""Audit 123 — SEED-FREE SOZ marker from the SIGNED (anti-symmetric) ImCoh FLOW via its
Helmholtz-Hodge source/sink decomposition, inside the Laplacian propagator framework.

The user's idea: stop throwing away the SIGN. `imcoh_abs` = <|signed|>_f is blind to anti-phase
coupling, which is exactly why module-isolation failed (audit_122): the SOZ looks like a hub in
magnitude. The raw signed ImCoh A (range [-1,1]) is provably ANTI-SYMMETRIC (A = -A^T, Im(C_ij) =
-Im(C_ji)) — a DIRECTED phase-lead/lag FLOW, not a symmetric signed graph. So the literal real
Kunegis L_s = |D| - A does not apply (A is not symmetric); its complex/magnetic analogue was tried
(audit_88) and came back negative.

The CORRECT in-frame object for an anti-symmetric flow is the graph Helmholtz-Hodge decomposition:
    divergence  b_i = sum_j A_ij                          net source(>0)/sink(<0) at node i
    potential   s   = L^+ b ,  L = D_mag - |A|             source-sink RANK (HodgeRank potential)
    residual    R   = A - |A| * (s_i - s_j)                rotational (cyclic / re-entrant) flow
where L is the Laplacian of the |ImCoh| magnitude graph. Crucially s = L^+ b SOLVES L s = b — the
steady state of DIFFUSION WITH SOURCES/SINKS b. So this is simultaneously (a) the seedless way to
leverage the negative links and (b) the concrete realization of the "diffusion with sinks"
generalization of the LRG propagator: b (from the signs) is the sink field, L^+ (the propagator)
solves it. SOZ hypothesis (drives/leads seizure propagation -> consistent phase lead): SOZ nodes
are EXTREME in the source-sink potential (high |s|), and/or carry anomalous cyclic (re-entrant)
flow, and are RIGID across phases. b is NOT node strength: a strongly-but-symmetrically coupled
node (equal lead and lag) has b ~ 0; only DIRECTIONAL consistency makes |b| large.

Critical preamble
=================
(1) Claim: a label-free signed-flow source-sink score (Hodge potential |s| and/or its cross-phase
    rigidity, and/or cyclic-flow) ranks SOZ above healthy CROSS-PATIENT (LOPO), beating |ImCoh|
    node strength and a label-shuffle null -- using sign information magnitude markers cannot.
(2) Null: LOPO AUC = 0.5; = strength-only (|ImCoh| abs degree); label-shuffle -> 0.5. Per-feature:
    SOZ |s|/|b| = matched-strength random nodes' (strength-stratified label permutation, the
    coupling-safe null; per-phase-independent graph-rewiring is the fake-positive trap, audit_121).
(3) Strongest alternatives: (a) hubness -> b cancels for symmetric coupling, NOT strength; explicit
    strength baseline + strat-null; (b) the sign is noise / band-average cancels it -> tested
    directly (if so, |b|~0 and AUC~0.5, reported); (c) proximity -> ImCoh zero-lag immune by
    construction (the sign is pure lagged phase). audit_88 prior: signed/magnetic at NODE level
    within-patient = no-beyond-hubness; THIS is the seedless source-sink Hodge potential, cross-
    phase, cross-patient -- a different object and framing, reported honestly either way.
(4) Reach: s/b/R use NO labels; potential is per-patient (own graph); held-out patient isolated;
    hub patients kept, read per-patient; cross-phase rigidity from the 4 phases; per band + fused.
(5) Falsification: if source-sink/cyclic LOPO ~ 0.5 / <= strength and the strat-null is at the
    median null per band, the signed flow does NOT yield a seedless marker on this data and we
    report it FLAT -- closing the last propagator-frame seedless angle. No acceptance gate.

Outputs (data/audit/epi_signed_flow/)
    signed_flow_per_patient.csv ; signed_flow_stratnull.csv ; signed_flow_lopo.csv ;
    signed_flow_node_predictions.csv ; README.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
from numpy.linalg import pinv
from scipy.stats import wilcoxon
from sklearn.linear_model import LogisticRegression

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.workflow.fc import load_fc_matrix

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_85_epi_propagator_recovery import MIN_EPI  # type: ignore
from audit_101_epi_marker_library import _conc  # type: ignore
import _epi_stratify as es  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_signed_flow"
CROSS_PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
SCORES = ["absdiv_mean", "abshodge_mean", "hodge_mean", "curl_mean", "curl_frac_mean", "strength"]
HEADLINE_FEATS = ["curl_mean", "curl_frac_mean", "abshodge_mean", "absdiv_mean"]
SRCSINK_FEATS = ["hodge_mean", "abshodge_mean", "hodge_rig", "absdiv_mean", "absdiv_rig"]
CYCLIC_FEATS = ["curl_mean", "curl_rig"]
CYCLIC_NORM_FEATS = ["curl_frac_mean", "curl_frac_rig"]   # strength-normalized rotational fraction
LOPO_FEATS = SRCSINK_FEATS + CYCLIC_FEATS + CYCLIC_NORM_FEATS + ["strength"]
N_PERM_STRAT = 500
N_STRENGTH_BINS = 4
EPS = 1e-12
SEED = 20260622


def _auc(score, epi):
    s = np.asarray(score, float)
    c, t = _conc(s[epi], s[~epi])
    return c / t if t else np.nan


def load_signed(pat, phase, band):
    A = load_fc_matrix(pat, phase, band, fc_method="imcoh")   # signed band-avg, anti-symmetric
    if A is None:
        raise FileNotFoundError(f"{pat}/{phase}/{band} signed imcoh not cached")
    A = np.asarray(A, float)
    np.fill_diagonal(A, 0.0)
    return 0.5 * (A - A.T)                                     # enforce exact anti-symmetry


def hodge_features(A):
    """Helmholtz-Hodge source/sink decomposition of an anti-symmetric flow A."""
    Mag = np.abs(A)
    deg = Mag.sum(1)
    L = np.diag(deg) - Mag                                     # |ImCoh| graph Laplacian
    b = A.sum(1)                                               # divergence = net source(>0)/sink(<0)
    s = pinv(L) @ b                                            # Hodge potential = L^+ b (diffusion-with-sources)
    s = s - np.median(s)
    grad = s[:, None] - s[None, :]
    R = A - Mag * grad                                        # rotational (cyclic / re-entrant) residual flow
    curl = np.sqrt((R ** 2).sum(1))
    return dict(div=b, absdiv=np.abs(b), hodge=s, abshodge=np.abs(s), curl=curl), deg


def per_patient_band(pat, band):
    pm = build_epi_masks(pat)
    chans = np.asarray(pm.channels, object)
    probes = np.asarray(pm.probes, object)
    N0 = len(chans)
    feats_per_phase, degs, avail = [], [], []
    for ph in CROSS_PHASES:
        try:
            A = load_signed(pat, ph, band)
        except Exception:
            continue
        if A.shape[0] != N0:
            continue
        f, deg = hodge_features(A)
        feats_per_phase.append(f); degs.append(deg); avail.append(ph)
    if len(avail) < 3:
        return None
    N = N0
    epi = np.asarray(pm.epi_mask, bool)
    if epi.sum() < MIN_EPI or (~epi).sum() < 5:
        return None
    strength = np.mean(degs, axis=0)

    def stack(key):
        return np.array([f[key] for f in feats_per_phase])      # (n_phase, N)
    feat = {"strength": strength}
    for key in ("div", "absdiv", "hodge", "abshodge", "curl"):
        M = stack(key)
        feat[f"{key}_mean"] = M.mean(0)
        feat[f"{key}_rig"] = -M.std(0)                          # high = state-invariant
    # strength-normalized rotational fraction (curl per unit |ImCoh| degree) — removes the trivial
    # "more flow -> more residual" scaling so we test rotationality BEYOND magnitude.
    cf = np.array([f["curl"] / (d + EPS) for f, d in zip(feats_per_phase, degs)])
    feat["curl_frac_mean"] = cf.mean(0)
    feat["curl_frac_rig"] = -cf.std(0)
    rows = []
    for sc in SCORES:
        rows.append(dict(patient=pat, band=band, score=sc, auc=float(_auc(feat[sc], epi)),
                         n_soz=int(epi.sum()), n_nodes=int(N), n_phases=len(avail)))
    return dict(rows=pd.DataFrame(rows), feat=feat, epi=epi, strength=strength,
                chans=chans, probes=probes, avail=avail)


def strat_null(score, epi, strength, rng, n_perm=N_PERM_STRAT, nbins=N_STRENGTH_BINS):
    s = np.asarray(score, float); ok = np.isfinite(s)
    s = s[ok]; epi = epi[ok]; st = np.asarray(strength, float)[ok]
    if epi.sum() < MIN_EPI or (~epi).sum() < 5:
        return np.nan, np.nan, np.nan
    obs = _auc(s, epi)
    edges = np.quantile(st, np.linspace(0, 1, nbins + 1)); edges[0] -= 1e-9; edges[-1] += 1e-9
    binid = np.clip(np.digitize(st, edges[1:-1]), 0, nbins - 1)
    idx_by_bin = [np.where(binid == bb)[0] for bb in range(nbins)]
    soz_per_bin = [int(epi[ib].sum()) for ib in idx_by_bin]
    perm = np.empty(n_perm)
    for k in range(n_perm):
        fake = np.zeros(epi.size, bool)
        for ib, ns in zip(idx_by_bin, soz_per_bin):
            if ns > 0 and ib.size > 0:
                fake[rng.choice(ib, size=min(ns, ib.size), replace=False)] = True
        perm[k] = _auc(s, fake)
    return obs, float(np.mean(perm >= obs)), float(obs - np.median(perm))


def _zscore_within(feat_dict, feats):
    cols = []
    for nm in feats:
        v = np.asarray(feat_dict[nm], float)
        mu, sd = np.nanmean(v), np.nanstd(v); sd = sd if sd > 1e-12 else 1.0
        z = (v - mu) / sd; z[~np.isfinite(z)] = 0.0
        cols.append(z)
    return np.column_stack(cols)


def lopo_auc(X, y, pid, pats, cols=None):
    Xu = X[:, cols] if cols is not None else X
    aucs = {}
    for pt in pats:
        tr, te = pid != pt, pid == pt
        if te.sum() == 0 or y[tr].sum() == 0 or y[te].sum() == 0:
            continue
        clf = LogisticRegression(max_iter=2000, C=1.0).fit(Xu[tr], y[tr])
        aucs[pt] = _auc(clf.predict_proba(Xu[te])[:, 1], y[te] == 1)
    return aucs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--n-perm", type=int, default=N_PERM_STRAT)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    t0 = time.time()

    cells, auc_rows, strat_rows = {}, [], []
    for band in args.bands:
        for pat in args.patients:
            res = per_patient_band(pat, band)
            if res is None:
                continue
            cells[(pat, band)] = res
            auc_rows.append(res["rows"])
            for sc in HEADLINE_FEATS:
                obs, p, exc = strat_null(res["feat"][sc], res["epi"], res["strength"],
                                         rng, n_perm=args.n_perm)
                strat_rows.append(dict(patient=pat, band=band, score=sc,
                                       auc_obs=obs, p_strat=p, excess_vs_matched=exc))
        print(f"[audit_123] {band}: {sum(1 for k in cells if k[1]==band)} cells ({time.time()-t0:.0f}s)")

    df = pd.concat(auc_rows, ignore_index=True)
    df.to_csv(OUT / "signed_flow_per_patient.csv", index=False)
    strat = pd.DataFrame(strat_rows)
    strat.to_csv(OUT / "signed_flow_stratnull.csv", index=False)

    coh_rows = []
    for band in [b for b in es.ALL_BANDS if b in set(df.band)]:
        for sc in SCORES:
            d = df[(df.band == band) & (df.score == sc)]
            v = d.auc.to_numpy(float); v = v[np.isfinite(v)]
            if v.size < 4:
                continue
            try:
                p = float(wilcoxon(v - 0.5, alternative="greater")[1])
            except ValueError:
                p = np.nan
            sd = strat[(strat.band == band) & (strat.score == sc)]
            coh_rows.append(dict(band=band, score=sc, n=int(v.size), median_auc=float(np.median(v)),
                                 n_above_half=int((v > 0.5).sum()), p_sign=p,
                                 median_excess_matched=float(sd.excess_vs_matched.median()) if not sd.empty else np.nan,
                                 n_strat_pass=int((sd.p_strat < 0.05).sum()) if not sd.empty else 0))
    coh = pd.DataFrame(coh_rows)
    coh.to_csv(OUT / "signed_flow_cohort.csv", index=False)

    lopo_rows, node_rows = [], []
    bands_present = [b for b in args.bands if any(k[1] == b for k in cells)]
    loop_bands = bands_present + (["__multiband__"] if len(bands_present) > 1 else [])
    for band in loop_bands:
        if band == "__multiband__":
            pats = [p for p in args.patients if all((p, b) in cells for b in bands_present)]
            if len(pats) < 4:
                continue
            X = np.vstack([np.column_stack([_zscore_within(cells[(p, b)]["feat"], LOPO_FEATS)
                                            for b in bands_present]) for p in pats])
            y = np.concatenate([cells[(p, bands_present[0])]["epi"].astype(int) for p in pats])
            pid = np.concatenate([np.array([p] * cells[(p, bands_present[0])]["epi"].size) for p in pats])
            names = [f"{b}:{f}" for b in bands_present for f in LOPO_FEATS]
        else:
            pats = [p for p in args.patients if (p, band) in cells]
            if len(pats) < 4:
                continue
            X = np.vstack([_zscore_within(cells[(p, band)]["feat"], LOPO_FEATS) for p in pats])
            y = np.concatenate([cells[(p, band)]["epi"].astype(int) for p in pats])
            pid = np.concatenate([np.array([p] * cells[(p, band)]["epi"].size) for p in pats])
            names = [f"{band}:{f}" for f in LOPO_FEATS]
        idx = {n: i for i, n in enumerate(names)}

        def cols_for(suff):
            return [i for n, i in idx.items() if n.split(":", 1)[1] in suff]
        groups = {
            "signed_full": cols_for(set(SRCSINK_FEATS) | set(CYCLIC_FEATS) | set(CYCLIC_NORM_FEATS)),
            "sourcesink": cols_for(set(SRCSINK_FEATS)),
            "cyclic": cols_for(set(CYCLIC_FEATS)),
            "cyclic_norm": cols_for(set(CYCLIC_NORM_FEATS)),
            "strength_only": cols_for({"strength"}),
            "signed_plus_strength": cols_for(set(SRCSINK_FEATS) | set(CYCLIC_FEATS) | set(CYCLIC_NORM_FEATS) | {"strength"}),
        }
        for g, cols in groups.items():
            av = np.array(list(lopo_auc(X, y, pid, pats, cols=cols).values()))
            lopo_rows.append(dict(band=band, group=g, n_feat=len(cols),
                                  mean_auc=float(np.nanmean(av)), median_auc=float(np.nanmedian(av)),
                                  n_above_half=int((av > 0.5).sum()), n=int(av.size)))
        yshuf = y.copy()
        for pt in pats:
            ii = np.where(pid == pt)[0]; yshuf[ii] = rng.permutation(yshuf[ii])
        av = np.array(list(lopo_auc(X, yshuf, pid, pats, cols=groups["signed_full"]).values()))
        lopo_rows.append(dict(band=band, group="label_shuffle_null", n_feat=len(groups["signed_full"]),
                              mean_auc=float(np.nanmean(av)), median_auc=float(np.nanmedian(av)),
                              n_above_half=int((av > 0.5).sum()), n=int(av.size)))
        for pt, a in lopo_auc(X, y, pid, pats, cols=groups["signed_full"]).items():
            lopo_rows.append(dict(band=band, group=f"  pp:{pt}", n_feat=np.nan, mean_auc=a,
                                  median_auc=a, n_above_half=np.nan, n=1))
        if band in bands_present:
            ci = groups["signed_full"]
            for pt in pats:
                tr, te = pid != pt, pid == pt
                if te.sum() == 0 or y[tr].sum() == 0:
                    continue
                clf = LogisticRegression(max_iter=2000, C=1.0).fit(X[tr][:, ci], y[tr])
                pte = clf.predict_proba(X[te][:, ci])[:, 1]
                cell = cells[(pt, band)]; chans, probes, yte = cell["chans"], cell["probes"], y[te]
                for j in range(te.sum()):
                    node_rows.append(dict(band=band, patient=pt, node=int(j),
                                          label=str(chans[j]) if j < len(chans) else "?",
                                          probe=str(probes[j]) if j < len(probes) else "?",
                                          is_soz=int(yte[j]), p_soz=float(pte[j])))
    lopo = pd.DataFrame(lopo_rows)
    lopo.to_csv(OUT / "signed_flow_lopo.csv", index=False)
    if node_rows:
        pd.DataFrame(node_rows).to_csv(OUT / "signed_flow_node_predictions.csv", index=False)

    rt = time.time() - t0
    print(f"\n[audit_123] SIGNED-FLOW source/sink (Hodge) cross-phase, LABEL-FREE ({rt:.1f}s)")
    print("  Per-score within-patient AUC + matched-strength (strat) null:")
    print("  band        score          med AUC  n>0.5  p_sign | med excess  n strat-pass")
    for band in [b for b in es.ALL_BANDS if (not coh.empty and b in set(coh.band))]:
        for sc in SCORES:
            r = coh[(coh.band == band) & (coh.score == sc)]
            if r.empty:
                continue
            r = r.iloc[0]
            print(f"  {band:11s} {sc:14s} {r.median_auc:.3f}   {int(r.n_above_half)}/{int(r.n)}"
                  f"   {r.p_sign:.3f} | {r.median_excess_matched:+.3f}     {int(r.n_strat_pass)}/{int(r.n)}")
        print()
    print("  ===== CROSS-PATIENT LOPO (HEADLINE — per band; SOZ vs healthy) =====")
    print("  band        group                   mean AUC  median  n>0.5")
    for band in loop_bands:
        sub = lopo[(lopo.band == band) & (~lopo.group.str.startswith("  pp:"))]
        for _, r in sub.iterrows():
            print(f"  {str(band):11s} {r.group:23s} {r.mean_auc:.3f}     {r.median_auc:.3f}   "
                  f"{int(r.n_above_half) if np.isfinite(r.n_above_half) else '?'}/{int(r.n)}")
        print()
    print(f"[audit_123] -> {OUT}")


if __name__ == "__main__":
    main()
