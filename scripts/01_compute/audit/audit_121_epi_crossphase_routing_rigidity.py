#!/usr/bin/env python3
"""Audit 121 — SEED-FREE SOZ marker from the CROSS-PHASE relation of i->j routing.

The dimension every prior seed-free attempt threw away. audit_118 (intrinsic), audit_120
(routing module), audit_94 (label-free community) all scored a SINGLE phase (rest_post) and
asked "which module is anomalous now?" — and all failed (the SOZ module is not separable from
physiological modules in one snapshot). This audit uses the PHASE axis: the propagator's i->j
routing is measured in all four phases (rest_pre, task_learn, task_test, rest_post) and the
marker is the CROSS-PHASE BEHAVIOUR of each node's routing profile, not its single-phase value.

Hypothesis (anchor taxonomy + epilepsy literature). The seizure zone is a STATE-INVARIANT
pathological module: its internal routing is RIGID — it does not reorganise when the brain
switches rest<->task — whereas healthy, functionally-engaged tissue reorganises its routing
with task demand. So SOZ nodes should keep the same routing partners across phases (high
cross-phase routing-profile persistence) MORE than healthy nodes of equal strength do. This is
the per-node, routing-view reading of the project's anchor/trace/reset cross-phase taxonomy.

The marker (LABEL-FREE, per node i, per band). Per phase p, build the strength-residualised
routing matrix M_p (routeff_ij = communicability_ij / resistance_ij, regressed off endpoint
strengths -> topology BEYOND magnitude; the audit_119 definition). Node i's routing profile in
phase p = row M_p[i, :]. Then:
    floor       = Spearman(profile_preA, profile_preB)              within-state reproducibility
    rig_prepost = Spearman(profile_pre,  profile_post)              rest->rest across the task
    rig_mean    = mean Spearman over all cross-phase pairs          state-invariance
    rig_min     = min  Spearman over all cross-phase pairs          survives the WORST transition
    reorg       = floor - rig_mean   (>=0; high = reorganiser)      anchor<->reset axis
    rig_excess  = resid(rig_mean ~ floor)                           rigid BEYOND own noise floor
NO labels enter any score. `floor` (split-half of the SAME state) is the per-node, strength- and
SNR-matched noise floor: a node is genuinely RIGID only if it survives a STATE change as well as
it survives mere RESAMPLING. That built-in baseline is what makes this different from a single-
phase score and is the honest control for "is it just a clean (high-SNR / hub) channel".

Critical preamble
=================
(1) Claim: a label-free CROSS-PHASE routing-rigidity score ranks SOZ above healthy, cohort-wide
    and CROSS-PATIENT (LOPO), beating node strength and a label-shuffle null — i.e. a deployable
    seed-free marker, which single-phase routing (audit_120) was NOT.
(2) Null: cross-patient LOPO AUC = 0.5; = strength-only LOPO; = floor-only LOPO (rigidity is just
    within-state reproducibility); label-shuffle -> 0.5. Per-feature: SOZ rigidity = matched-
    strength random nodes' rigidity (strength-stratified label permutation).
(3) Strongest alternatives: (a) hubness — every M_p is strength-residualised by construction,
    PLUS an explicit strength baseline + strength-only LOPO + a strength-stratified permutation
    null; (b) high-SNR/clean channels look rigid — controlled by `floor` (a clean channel is
    reproducible WITHIN a state too, so reorg/rig_excess remove it); (c) proximity — the |ImCoh|
    substrate is zero-lag immune (Nolte 2004), routing features are global path properties, no
    coordinates enter.
(4) WHY NOT a graph-rewiring matched-strength null here (the key design choice): a per-phase-
    INDEPENDENT strength-preserving rewiring DESTROYS cross-phase coupling, so surrogate rigidity
    collapses to noise and the observed score "beats" it TRIVIALLY for every node — a guaranteed
    fake positive (audit_63 README documents exactly this; there is no well-defined coordinated
    cross-phase strength surrogate). The fair, coupling-safe matched-strength null is therefore
    a STRENGTH-STRATIFIED LABEL PERMUTATION on the real coupled data: is the SOZ more rigid than
    random nodes of the SAME strength? Read PER-PATIENT (excess + p_strat + who carries it), never
    a cohort sign-test (the audit_120 mirage lesson).
(5) Falsification: if LOPO AUC ~ 0.5 or <= strength-only / floor-only, and the per-feature
    strat-null is at chance per-patient, the cross-phase routing relation does NOT yield a
    seed-free marker on this data and we report it FLAT — a clean negative that closes the
    phase-axis idea, not a dressed-up wobble. No pre-registered acceptance gate.

Outputs (data/audit/epi_crossphase_routing/)
    crossphase_rigidity_per_patient.csv ; crossphase_rigidity_cohort.csv ;
    crossphase_rigidity_stratnull.csv ; crossphase_lopo.csv ;
    crossphase_node_predictions.csv ; README.md
"""
from __future__ import annotations

import argparse
import sys
import time
from itertools import combinations

import numpy as np
import pandas as pd
from numpy.linalg import pinv
from scipy.linalg import expm
from scipy.stats import spearmanr, wilcoxon
from sklearn.linear_model import LogisticRegression

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.utils.io.patient import build_epi_masks

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc, ensure_half_fcs  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI  # type: ignore
from audit_101_epi_marker_library import _conc, _resid  # type: ignore
import _epi_stratify as es  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_crossphase_routing"
CROSS_PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
# label-free per-node rigidity scores (SOZ hypothesis: rigid => HIGH, except reorg => LOW)
SCORES = ["rig_prepost", "rig_mean", "rig_min", "rig_excess", "neg_reorg", "strength"]
LOPO_FEATS = ["floor", "rig_prepost", "rig_mean", "rig_min", "reorg", "rig_excess", "strength"]
HEADLINE_FEATS = ["rig_mean", "rig_excess"]   # pre-specified for the strat-null diagnostic
N_PERM_STRAT = 500
N_STRENGTH_BINS = 4
SEED = 20260622


# ---------------------------------------------------------------------------
# Routing matrix (audit_119/120 definition: routing residual beyond endpoint strength)
# ---------------------------------------------------------------------------
def _routing_matrix(W):
    """Strength-residualised routing-efficiency pair matrix M (symmetric, zero diag) + strength."""
    N = W.shape[0]
    d = W.sum(1)
    L = np.diag(d) - W
    Lp = pinv(L)
    dLp = np.diag(Lp)
    R = dLp[:, None] + dLp[None, :] - 2 * Lp                 # effective-resistance distance
    rho = max(abs(np.linalg.eigvalsh(W)).max(), 1e-9)
    G = expm(W / rho)                                        # communicability
    routeff = G / (R + 1e-9)                                 # walks per resistance = private routing
    np.fill_diagonal(routeff, 0.0)
    iu = np.triu_indices(N, 1)
    si, sj = d[iu[0]], d[iu[1]]
    X = np.column_stack([np.ones_like(si), si + sj, si * sj,
                         np.log1p(si) + np.log1p(sj)])
    f = routeff[iu]
    beta, *_ = np.linalg.lstsq(X, f, rcond=None)
    res = f - X @ beta                                       # routing BEYOND endpoint strength
    M = np.zeros((N, N)); M[iu] = res; M = M + M.T
    return M, d


def _profile_corr_vec(Mp, Mq):
    """Per-node Spearman of routing profile (row) between two phases, self excluded."""
    N = Mp.shape[0]
    out = np.full(N, np.nan)
    allidx = np.arange(N)
    for i in range(N):
        m = allidx != i
        a, b = Mp[i, m], Mq[i, m]
        if np.std(a) < 1e-12 or np.std(b) < 1e-12:
            continue
        out[i] = spearmanr(a, b).correlation
    return out


def _auc(score, epi):
    s = np.asarray(score, float)
    c, t = _conc(s[epi], s[~epi])
    return c / t if t else np.nan


# ---------------------------------------------------------------------------
# Per-patient/band cross-phase rigidity features (all label-free)
# ---------------------------------------------------------------------------
def per_patient_band(pat, band):
    pm = build_epi_masks(pat)
    chans = np.asarray(pm.channels, object)
    N0 = len(chans)
    # within-state split-half floor (rest_pre_A / rest_pre_B)
    try:
        MA = _routing_matrix(load_phase_fc(pat, "rest_pre_A", band))[0]
        MB = _routing_matrix(load_phase_fc(pat, "rest_pre_B", band))[0]
    except Exception:
        return None
    N = MA.shape[0]
    if N != N0:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    if epi.sum() < MIN_EPI or (~epi).sum() < 5:
        return None
    # cross-phase routing matrices (whatever is available; need >=3 incl pre & post)
    Mp, dvecs, avail = {}, [], []
    for ph in CROSS_PHASES:
        try:
            W = load_phase_fc(pat, ph, band)
            M, d = _routing_matrix(W)
        except Exception:
            continue
        if M.shape[0] != N:
            continue
        Mp[ph] = M; dvecs.append(d); avail.append(ph)
    if "rest_pre" not in avail or "rest_post" not in avail or len(avail) < 3:
        return None
    strength = np.mean(np.array(dvecs), axis=0)

    floor = _profile_corr_vec(MA, MB)
    pair_corrs = [_profile_corr_vec(Mp[p], Mp[q]) for p, q in combinations(avail, 2)]
    PC = np.array(pair_corrs)                                # (n_pairs, N)
    with np.errstate(invalid="ignore"):
        rig_mean = np.nanmean(PC, axis=0)
        rig_min = np.nanmin(PC, axis=0)
    rig_prepost = _profile_corr_vec(Mp["rest_pre"], Mp["rest_post"])
    reorg = floor - rig_mean                                  # >0 reorganiser, <=0 anchor/rigid
    rig_excess = _resid(rig_mean, floor)                     # rigid beyond own within-state floor
    neg_reorg = -reorg

    feat = dict(floor=floor, rig_prepost=rig_prepost, rig_mean=rig_mean,
                rig_min=rig_min, reorg=reorg, rig_excess=rig_excess,
                neg_reorg=neg_reorg, strength=strength)
    rows = []
    for sc in SCORES:
        rows.append(dict(patient=pat, band=band, score=sc,
                         auc=float(_auc(feat[sc], epi)),
                         n_soz=int(epi.sum()), n_nodes=int(N),
                         n_phases=len(avail)))
    return dict(rows=pd.DataFrame(rows), feat=feat, epi=epi, strength=strength,
                chans=chans, probes=np.asarray(pm.probes, object), avail=avail)


# ---------------------------------------------------------------------------
# Strength-stratified label-permutation null (coupling-safe matched-strength)
# ---------------------------------------------------------------------------
def strat_null(score, epi, strength, rng, n_perm=N_PERM_STRAT, nbins=N_STRENGTH_BINS):
    """Is SOZ rigidity > random nodes of the SAME strength? Permute SOZ labels within
    strength bins, keeping per-bin SOZ counts -> strength-matched fake-SOZ sets."""
    s = np.asarray(score, float)
    ok = np.isfinite(s)
    s = s[ok]; epi = epi[ok]; st = np.asarray(strength, float)[ok]
    if epi.sum() < MIN_EPI or (~epi).sum() < 5:
        return np.nan, np.nan, np.nan
    obs = _auc(s, epi)
    # strength bins (quantile edges; fall back to fewer bins if degenerate)
    edges = np.quantile(st, np.linspace(0, 1, nbins + 1))
    edges[0] -= 1e-9; edges[-1] += 1e-9
    binid = np.clip(np.digitize(st, edges[1:-1]), 0, nbins - 1)
    idx_by_bin = [np.where(binid == b)[0] for b in range(nbins)]
    soz_per_bin = [int(epi[ib].sum()) for ib in idx_by_bin]
    perm_aucs = np.empty(n_perm)
    for k in range(n_perm):
        fake = np.zeros(epi.size, bool)
        for ib, nsoz in zip(idx_by_bin, soz_per_bin):
            if nsoz > 0 and ib.size > 0:
                fake[rng.choice(ib, size=min(nsoz, ib.size), replace=False)] = True
        perm_aucs[k] = _auc(s, fake)
    p_upper = float(np.mean(perm_aucs >= obs))               # one-sided: SOZ more rigid
    excess = float(obs - np.median(perm_aucs))
    return obs, p_upper, excess


# ---------------------------------------------------------------------------
# Cross-patient LOPO (the headline bar)
# ---------------------------------------------------------------------------
def _zscore_within(feat_dict, feats):
    cols = []
    for nm in feats:
        v = np.asarray(feat_dict[nm], float)
        mu, sd = np.nanmean(v), np.nanstd(v)
        sd = sd if sd > 1e-12 else 1.0
        z = (v - mu) / sd
        z[~np.isfinite(z)] = 0.0
        cols.append(z)
    return np.column_stack(cols)


def lopo_auc(X, y, pid, pats, cols=None):
    Xu = X[:, cols] if cols is not None else X
    aucs = {}
    preds = {}
    for pt in pats:
        tr, te = pid != pt, pid == pt
        if te.sum() == 0 or y[tr].sum() == 0 or y[te].sum() == 0:
            continue
        clf = LogisticRegression(max_iter=2000, C=1.0).fit(Xu[tr], y[tr])
        p = clf.predict_proba(Xu[te])[:, 1]
        aucs[pt] = _auc(p, y[te] == 1)
        preds[pt] = (p, y[te])
    return aucs, preds


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--n-perm", type=int, default=N_PERM_STRAT)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    t0 = time.time()

    print("[audit_121] pre-flight: ensure rest_pre half FCs cached")
    for pat in args.patients:
        try:
            ensure_half_fcs(pat, args.bands)
        except Exception as e:
            print(f"  {pat}: half-FC cache failed ({e})")

    # ---- per (patient, band) features ----
    cells = {}
    auc_rows, strat_rows = [], []
    for band in args.bands:
        for pat in args.patients:
            res = per_patient_band(pat, band)
            if res is None:
                continue
            cells[(pat, band)] = res
            auc_rows.append(res["rows"])
            for sc in HEADLINE_FEATS:
                obs, p, exc = strat_null(res["feat"][sc], res["epi"],
                                         res["strength"], rng, n_perm=args.n_perm)
                strat_rows.append(dict(patient=pat, band=band, score=sc,
                                       auc_obs=obs, p_strat=p, excess_vs_matched=exc))
        print(f"[audit_121] {band}: {sum(1 for k in cells if k[1]==band)} cells")

    df = pd.concat(auc_rows, ignore_index=True)
    df.to_csv(OUT / "crossphase_rigidity_per_patient.csv", index=False)
    strat = pd.DataFrame(strat_rows)
    strat.to_csv(OUT / "crossphase_rigidity_stratnull.csv", index=False)

    # ---- per-band per-score cohort (median AUC + sign vs 0.5, DIAGNOSTIC ONLY) ----
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
            # strat-null cohort view (median per-patient excess over matched-strength nodes)
            sd = strat[(strat.band == band) & (strat.score == sc)]
            med_exc = float(sd.excess_vs_matched.median()) if not sd.empty else np.nan
            n_pass = int((sd.p_strat < 0.05).sum()) if not sd.empty else 0
            coh_rows.append(dict(band=band, score=sc, n=int(v.size),
                                 median_auc=float(np.median(v)),
                                 n_above_half=int((v > 0.5).sum()), p_sign=p,
                                 median_excess_matched=med_exc, n_strat_pass=n_pass))
    coh = pd.DataFrame(coh_rows)
    coh.to_csv(OUT / "crossphase_rigidity_cohort.csv", index=False)

    # ---- CROSS-PATIENT LOPO (headline): multiband rigidity vs strength-only vs floor vs shuffle ----
    # patients with ALL requested bands present (aligned-node multiband fusion)
    full_pats = [p for p in args.patients
                 if all((p, b) in cells for b in args.bands)]
    lopo_summary = []
    node_pred_rows = []
    if len(full_pats) >= 4:
        Xblocks, ylist, pidlist = [], [], []
        feat_names = []
        for p in full_pats:
            perband = [_zscore_within(cells[(p, b)]["feat"], LOPO_FEATS) for b in args.bands]
            Xblocks.append(np.column_stack(perband))
            ylist.append(cells[(p, args.bands[0])]["epi"].astype(int))
            pidlist.append(np.array([p] * cells[(p, args.bands[0])]["epi"].size))
        feat_names = [f"{b}:{f}" for b in args.bands for f in LOPO_FEATS]
        X = np.vstack(Xblocks); y = np.concatenate(ylist); pid = np.concatenate(pidlist)

        idx = {n: i for i, n in enumerate(feat_names)}
        groups = {
            "full_rigidity": None,
            "strength_only": [idx[f"{b}:strength"] for b in args.bands],
            "floor_only": [idx[f"{b}:floor"] for b in args.bands],
            "rigidity_no_strength": [i for n, i in idx.items() if not n.endswith(":strength")],
        }
        for g, cols in groups.items():
            aucs, preds = lopo_auc(X, y, pid, full_pats, cols=cols)
            av = np.array(list(aucs.values()))
            lopo_summary.append(dict(group=g, n_feat=(len(cols) if cols else X.shape[1]),
                                     mean_auc=float(np.nanmean(av)),
                                     median_auc=float(np.nanmedian(av)),
                                     n_above_half=int((av > 0.5).sum()), n=int(av.size)))
            if g == "full_rigidity":
                for pt, a in aucs.items():
                    lopo_summary.append(dict(group=f"  per-patient:{pt}", n_feat=np.nan,
                                             mean_auc=a, median_auc=a, n_above_half=np.nan, n=1))
                # node-level P(SOZ) for candidate inspection
                for pt in full_pats:
                    tr, te = pid != pt, pid == pt
                    if te.sum() == 0 or y[tr].sum() == 0:
                        continue
                    clf = LogisticRegression(max_iter=2000, C=1.0).fit(X[tr], y[tr])
                    pte = clf.predict_proba(X[te])[:, 1]
                    cell0 = cells[(pt, args.bands[0])]
                    chans, probes = cell0["chans"], cell0["probes"]
                    yte = y[te]
                    for j in range(te.sum()):
                        node_pred_rows.append(dict(patient=pt, node=int(j),
                                                   label=str(chans[j]) if j < len(chans) else "?",
                                                   probe=str(probes[j]) if j < len(probes) else "?",
                                                   is_soz=int(yte[j]), p_soz=float(pte[j])))
        # label-shuffle null on the full model
        yshuf = y.copy()
        for pt in full_pats:
            ii = np.where(pid == pt)[0]
            yshuf[ii] = rng.permutation(yshuf[ii])
        aucs_sh, _ = lopo_auc(X, yshuf, pid, full_pats)
        lopo_summary.append(dict(group="label_shuffle_null", n_feat=X.shape[1],
                                 mean_auc=float(np.nanmean(list(aucs_sh.values()))),
                                 median_auc=float(np.nanmedian(list(aucs_sh.values()))),
                                 n_above_half=int(np.sum(np.array(list(aucs_sh.values())) > 0.5)),
                                 n=len(aucs_sh)))
    lopo_df = pd.DataFrame(lopo_summary)
    lopo_df.to_csv(OUT / "crossphase_lopo.csv", index=False)
    if node_pred_rows:
        pd.DataFrame(node_pred_rows).to_csv(OUT / "crossphase_node_predictions.csv", index=False)

    # ---- report ----
    rt = time.time() - t0
    print(f"\n[audit_121] CROSS-PHASE routing rigidity, LABEL-FREE ({rt:.1f}s)")
    print("  Per-score cohort AUC (SOZ vs healthy; DIAGNOSTIC, not the headline):")
    print("  band        score        med AUC  n>0.5   p_sign  | matched-strength: med excess  n strat-pass")
    coh_bands = [b for b in es.ALL_BANDS if (not coh.empty and b in set(coh.band))]
    for band in coh_bands:
        for sc in SCORES:
            r = coh[(coh.band == band) & (coh.score == sc)]
            if r.empty:
                continue
            r = r.iloc[0]
            print(f"  {band:11s} {sc:11s}  {r.median_auc:.3f}    {int(r.n_above_half)}/{int(r.n)}"
                  f"     {r.p_sign:.3f}   |  {r.median_excess_matched:+.3f}        {int(r.n_strat_pass)}/{int(r.n)}")
        print()

    if not lopo_df.empty:
        print("  ===== CROSS-PATIENT LOPO (THE HEADLINE — multiband fusion) =====")
        for _, r in lopo_df[~lopo_df.group.str.startswith("  per-patient")].iterrows():
            print(f"    {r.group:22s} mean AUC {r.mean_auc:.3f}  median {r.median_auc:.3f}  "
                  f"{r.n_above_half if np.isfinite(r.n_above_half) else '?'}/{int(r.n)}>0.5")
        print("    per-patient (full rigidity model):")
        for _, r in lopo_df[lopo_df.group.str.startswith("  per-patient")].sort_values(
                "mean_auc", ascending=False).iterrows():
            print(f"      {r.group.strip():18s} AUC {r.mean_auc:.3f}")
    print(f"\n[audit_121] -> {OUT}")


if __name__ == "__main__":
    main()
