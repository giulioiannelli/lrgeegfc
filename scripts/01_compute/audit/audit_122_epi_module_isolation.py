#!/usr/bin/env python3
"""Audit 122 — SEED-FREE SOZ marker from MODULE ISOLATION (graph conductance) and its
CROSS-PHASE rigidity, entirely inside the LRG/Laplacian propagator framework.

The last structurally-different seedless angle. Every prior seedless attempt scored a NODE and
asked "is it anomalous"; all reduced to node strength or hit the label-free selection wall
(audit_80/81/94/118/120/121). This one scores a MODULE property that is provably NOT node
strength: ISOLATION (conductance / boundary-to-volume). The "epileptic network" hypothesis says
the SOZ is a near-CLOSED diffusion TRAP — low outflow conductance: heat starting inside is slow
to leave. A uniformly high-strength (hub) module does NOT have low conductance (its boundary
edges are strong too), so isolation can survive the matched-strength wall that sank node-rigidity.
Cross-phase: a pathological trap should be STATE-INVARIANT (rigid) — isolated in rest AND task.

Connection to a future "diffusion-with-sinks" generalization. Conductance phi(M) = (heat leaving
M per unit time)/(volume) is the absorbing-boundary (Dirichlet) special case of the killed-
diffusion propagator e^{-tau(L+Gamma)}, Gamma=diag(sink rates): heat-retention in M = survival
probability of a walk absorbed at the boundary of M. So this seedless isolation marker is the
undirected precursor of a node-resolved sinks model; this audit measures the trap, the future
generalization learns where Gamma>0 (which nodes absorb). Stays in the propagator frame.

Method (no dendrogram cut -> avoids the "one cluster dominates at any K" degeneracy; memory
feedback_no_partition_metrics_use_rho_coph). Local graph clustering, Andersen-Chung-Lang:
  - Reference graph A_ref = mean over the 4 phases (rest_pre, task_learn, task_test, rest_post).
  - For each node i and teleport alpha in {0.85, 0.95} (two diffusion scales): personalized
    PageRank vector pi_i = row i of (1-a)(I - a P)^{-1}, P = row-stochastic walk on A_ref. Sweep:
    sort nodes by pi_i(v)/deg(v), take the prefix of MINIMUM conductance -> S_i^alpha, node i's
    own best-isolated community (seedless: i seeds its OWN sweep; NO labels).
  - For each phase p, compute conductance of the FIXED set S_i^alpha on that phase's graph A_p.
    Per-node features: neg_cond = -mean_phase(phi) (high = isolated trap); cond_stab =
    -std_phase(phi) (high = state-invariantly isolated = rigid); internal_strength and size of
    S_i (the "dense/large module" BASELINE isolation must beat); node strength.
  - SOZ hypothesis: SOZ nodes sit in low-conductance (high neg_cond), state-invariant (high
    cond_stab) communities. Test cross-patient (LOPO), per band and fused, vs strength-only, vs
    the internal-strength+size baseline, vs label-shuffle, with a strength-stratified null.

Critical preamble
=================
(1) Claim: a label-free MODULE-ISOLATION score (local conductance + its cross-phase rigidity)
    ranks SOZ above healthy CROSS-PATIENT (LOPO), beating node strength, an internal-density
    baseline, and a label-shuffle null -- the seedless marker single-phase/node scores were not.
(2) Null: LOPO AUC = 0.5; = strength-only; = internal_strength+size (isolated == merely dense/big);
    label-shuffle -> 0.5. Per-feature: SOZ conductance = matched-strength random nodes'
    conductance (strength-stratified label permutation, coupling-safe -- see audit_121 preamble
    for why a per-phase-independent graph-rewiring null is the fake-positive trap here too).
(3) Strongest alternatives: (a) hubness -> conductance is volume-normalised (a hub module is NOT
    low-conductance) + explicit strength baseline + strat-null; (b) isolated == just dense ->
    internal_strength+size baseline rung; (c) small-community artifact -> sweep min size >= 2,
    size included as a control feature; (d) proximity -> |ImCoh| zero-lag immune (Nolte 2004),
    conductance is a global path property.
(4) Reach: communities defined with NO labels (each node seeds its own sweep); membership fixed
    on A_ref then conductance recomputed per phase (clean isolation vs rigidity split); held-out
    patient never trains; hub patients kept and read per-patient; two diffusion scales.
(5) Falsification: if isolation LOPO ~ 0.5 / <= strength / <= internal-density, and the strat-null
    is at the median null per band, module isolation does NOT yield a seedless marker on this data
    and we report it FLAT -- the principled close of the seedless question. No acceptance gate.

Outputs (data/audit/epi_module_isolation/)
    module_isolation_per_patient.csv ; module_isolation_stratnull.csv ;
    module_isolation_lopo.csv ; module_isolation_node_predictions.csv ; README.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.linear_model import LogisticRegression

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.utils.io.patient import build_epi_masks

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI  # type: ignore
from audit_101_epi_marker_library import _conc, _resid  # type: ignore
import _epi_stratify as es  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_module_isolation"
CROSS_PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
ALPHAS = [0.85, 0.95]                  # two diffusion / locality scales for the PPR sweep
SCORES = ["neg_cond_a85", "neg_cond_a95", "cond_stab_a85", "internal_strength_a85",
          "size_a85", "strength"]
HEADLINE_FEATS = ["neg_cond_a85", "neg_cond_a95"]
ISO_FEATS = ["neg_cond_a85", "neg_cond_a95", "cond_stab_a85", "cond_stab_a95"]
DENS_FEATS = ["internal_strength_a85", "size_a85", "internal_strength_a95", "size_a95"]
LOPO_FEATS = ISO_FEATS + DENS_FEATS + ["strength"]
N_PERM_STRAT = 500
N_STRENGTH_BINS = 4
EPS = 1e-12
SEED = 20260622


def _auc(score, epi):
    s = np.asarray(score, float)
    c, t = _conc(s[epi], s[~epi])
    return c / t if t else np.nan


def _ppr_rows(W, alpha):
    """Row i = personalized PageRank seeded at node i (row-stochastic walk)."""
    N = W.shape[0]
    deg = W.sum(1)
    P = W / np.maximum(deg, EPS)[:, None]
    return (1.0 - alpha) * np.linalg.inv(np.eye(N) - alpha * P)


def _sweep_community(ppr_vec, W, deg, min_size=2):
    """Andersen-Chung-Lang sweep: prefix of minimum conductance along pi/deg order -> node set."""
    N = W.shape[0]
    order = np.argsort(-(ppr_vec / np.maximum(deg, EPS)))
    A = W[np.ix_(order, order)]
    dego = deg[order]
    voltot = deg.sum()
    cumvol = np.cumsum(dego)                      # vol of first t nodes = cumvol[t-1]
    rowto = np.tril(A, -1).sum(axis=1)            # rowto[m] = edges from node m to earlier nodes
    twiceint = 2.0 * np.cumsum(rowto)             # A[:t,:t].sum() = twiceint[t-1]
    vol = cumvol[:-1]                             # t = 1 .. N-1
    cut = np.maximum(vol - twiceint[:-1], 0.0)
    denom = np.minimum(vol, voltot - vol)
    phi = np.where(denom > 0, cut / np.maximum(denom, EPS), np.inf)
    phi[: max(0, min_size - 1)] = np.inf          # require community size >= min_size
    phi[-(min_size - 1):] = np.inf if min_size > 1 else phi[-1]
    t = int(np.argmin(phi)) + 1
    S = np.zeros(N, bool); S[order[:t]] = True
    return S


def _conductance(S, W, deg):
    volS = deg[S].sum()
    volC = deg.sum() - volS
    denom = min(volS, volC)
    if denom <= 0:
        return np.nan
    cut = W[np.ix_(S, ~S)].sum()
    return float(cut / denom)


def per_patient_band(pat, band):
    pm = build_epi_masks(pat)
    chans = np.asarray(pm.channels, object)
    probes = np.asarray(pm.probes, object)
    N0 = len(chans)
    Ws, degs, avail = [], [], []
    for ph in CROSS_PHASES:
        try:
            W = load_phase_fc(pat, ph, band)
        except Exception:
            continue
        if W.shape[0] != N0:
            continue
        Ws.append(W); degs.append(W.sum(1)); avail.append(ph)
    if len(avail) < 3:
        return None
    N = N0
    epi = np.asarray(pm.epi_mask, bool)
    if epi.sum() < MIN_EPI or (~epi).sum() < 5:
        return None
    Wref = np.mean(Ws, axis=0)
    degref = Wref.sum(1)
    strength = np.mean(degs, axis=0)

    feat = {"strength": strength}
    for a in ALPHAS:
        tag = f"a{int(a * 100)}"
        M = _ppr_rows(Wref, a)
        cond_per_phase = np.full((len(Ws), N), np.nan)
        internal = np.full(N, np.nan)
        size = np.full(N, np.nan)
        for i in range(N):
            S = _sweep_community(M[i], Wref, degref)
            size[i] = S.sum()
            internal[i] = Wref[np.ix_(S, S)].sum() / max(1, S.sum())   # mean internal degree
            for pi, (W, dg) in enumerate(zip(Ws, degs)):
                cond_per_phase[pi, i] = _conductance(S, W, dg)
        with np.errstate(invalid="ignore"):
            cmean = np.nanmean(cond_per_phase, axis=0)
            cstd = np.nanstd(cond_per_phase, axis=0)
        feat[f"neg_cond_{tag}"] = -cmean              # high = isolated trap
        feat[f"cond_stab_{tag}"] = -cstd              # high = state-invariantly isolated (rigid)
        feat[f"internal_strength_{tag}"] = internal   # baseline: isolated vs merely dense
        feat[f"size_{tag}"] = size.astype(float)

    rows = []
    for sc in SCORES:
        rows.append(dict(patient=pat, band=band, score=sc,
                         auc=float(_auc(feat[sc], epi)),
                         n_soz=int(epi.sum()), n_nodes=int(N), n_phases=len(avail)))
    return dict(rows=pd.DataFrame(rows), feat=feat, epi=epi, strength=strength,
                chans=chans, probes=probes, avail=avail)


# --- coupling-safe matched-strength null: strength-stratified label permutation (audit_121) ---
def strat_null(score, epi, strength, rng, n_perm=N_PERM_STRAT, nbins=N_STRENGTH_BINS):
    s = np.asarray(score, float)
    ok = np.isfinite(s)
    s = s[ok]; epi = epi[ok]; st = np.asarray(strength, float)[ok]
    if epi.sum() < MIN_EPI or (~epi).sum() < 5:
        return np.nan, np.nan, np.nan
    obs = _auc(s, epi)
    edges = np.quantile(st, np.linspace(0, 1, nbins + 1)); edges[0] -= 1e-9; edges[-1] += 1e-9
    binid = np.clip(np.digitize(st, edges[1:-1]), 0, nbins - 1)
    idx_by_bin = [np.where(binid == b)[0] for b in range(nbins)]
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
        mu, sd = np.nanmean(v), np.nanstd(v)
        sd = sd if sd > 1e-12 else 1.0
        z = (v - mu) / sd
        z[~np.isfinite(z)] = 0.0
        cols.append(z)
    return np.column_stack(cols)


def lopo_auc(X, y, pid, pats, cols=None):
    Xu = X[:, cols] if cols is not None else X
    aucs, preds = {}, {}
    for pt in pats:
        tr, te = pid != pt, pid == pt
        if te.sum() == 0 or y[tr].sum() == 0 or y[te].sum() == 0:
            continue
        clf = LogisticRegression(max_iter=2000, C=1.0).fit(Xu[tr], y[tr])
        p = clf.predict_proba(Xu[te])[:, 1]
        aucs[pt] = _auc(p, y[te] == 1); preds[pt] = (p, y[te])
    return aucs, preds


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
        print(f"[audit_122] {band}: {sum(1 for k in cells if k[1]==band)} cells ({time.time()-t0:.0f}s)")

    df = pd.concat(auc_rows, ignore_index=True)
    df.to_csv(OUT / "module_isolation_per_patient.csv", index=False)
    strat = pd.DataFrame(strat_rows)
    strat.to_csv(OUT / "module_isolation_stratnull.csv", index=False)

    # per-band per-score cohort (diagnostic)
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
            med_exc = float(sd.excess_vs_matched.median()) if not sd.empty else np.nan
            n_pass = int((sd.p_strat < 0.05).sum()) if not sd.empty else 0
            coh_rows.append(dict(band=band, score=sc, n=int(v.size),
                                 median_auc=float(np.median(v)),
                                 n_above_half=int((v > 0.5).sum()), p_sign=p,
                                 median_excess_matched=med_exc, n_strat_pass=n_pass))
    coh = pd.DataFrame(coh_rows)
    coh.to_csv(OUT / "module_isolation_cohort.csv", index=False)

    # ---- per-band LOPO (HEADLINE; multiband as secondary) ----
    lopo_rows, node_rows = [], []
    bands_present = [b for b in args.bands if any(k[1] == b for k in cells)]
    for band in bands_present + (["__multiband__"] if len(bands_present) > 1 else []):
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

        def cols_for(suffixes):
            return [i for n, i in idx.items() if n.split(":", 1)[1] in suffixes]
        groups = {
            "isolation_full": cols_for(set(ISO_FEATS)),
            "strength_only": cols_for({"strength"}),
            "internal_density": cols_for(set(DENS_FEATS)),
            "isolation_plus_strength": cols_for(set(ISO_FEATS) | {"strength"}),
        }
        for g, cols in groups.items():
            aucs, _ = lopo_auc(X, y, pid, pats, cols=cols)
            av = np.array(list(aucs.values()))
            lopo_rows.append(dict(band=band, group=g, n_feat=len(cols),
                                  mean_auc=float(np.nanmean(av)), median_auc=float(np.nanmedian(av)),
                                  n_above_half=int((av > 0.5).sum()), n=int(av.size)))
        # shuffle null + per-patient + node preds on isolation_full
        yshuf = y.copy()
        for pt in pats:
            ii = np.where(pid == pt)[0]; yshuf[ii] = rng.permutation(yshuf[ii])
        aucs_sh, _ = lopo_auc(X, yshuf, pid, pats, cols=groups["isolation_full"])
        av = np.array(list(aucs_sh.values()))
        lopo_rows.append(dict(band=band, group="label_shuffle_null", n_feat=len(groups["isolation_full"]),
                              mean_auc=float(np.nanmean(av)), median_auc=float(np.nanmedian(av)),
                              n_above_half=int((av > 0.5).sum()), n=int(av.size)))
        aucs_if, _ = lopo_auc(X, y, pid, pats, cols=groups["isolation_full"])
        for pt, a in aucs_if.items():
            lopo_rows.append(dict(band=band, group=f"  pp:{pt}", n_feat=np.nan,
                                  mean_auc=a, median_auc=a, n_above_half=np.nan, n=1))
        if band in bands_present:                 # node-level P(SOZ) per band for inspection
            ci = groups["isolation_full"]
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
    lopo.to_csv(OUT / "module_isolation_lopo.csv", index=False)
    if node_rows:
        pd.DataFrame(node_rows).to_csv(OUT / "module_isolation_node_predictions.csv", index=False)

    # ---- report ----
    rt = time.time() - t0
    print(f"\n[audit_122] MODULE ISOLATION (local conductance) cross-phase, LABEL-FREE ({rt:.1f}s)")
    print("  Per-score within-patient AUC + matched-strength (strat) null:")
    print("  band        score                 med AUC  n>0.5  p_sign | med excess  n strat-pass")
    for band in [b for b in es.ALL_BANDS if (not coh.empty and b in set(coh.band))]:
        for sc in SCORES:
            r = coh[(coh.band == band) & (coh.score == sc)]
            if r.empty:
                continue
            r = r.iloc[0]
            print(f"  {band:11s} {sc:21s} {r.median_auc:.3f}   {int(r.n_above_half)}/{int(r.n)}"
                  f"   {r.p_sign:.3f} | {r.median_excess_matched:+.3f}     {int(r.n_strat_pass)}/{int(r.n)}")
        print()
    print("  ===== CROSS-PATIENT LOPO (HEADLINE — per band; SOZ vs healthy) =====")
    print("  band        group                    mean AUC  median  n>0.5")
    for band in bands_present + (["__multiband__"] if "__multiband__" in set(lopo.band) else []):
        sub = lopo[(lopo.band == band) & (~lopo.group.str.startswith("  pp:"))]
        for _, r in sub.iterrows():
            print(f"  {str(band):11s} {r.group:24s} {r.mean_auc:.3f}     {r.median_auc:.3f}   "
                  f"{int(r.n_above_half) if np.isfinite(r.n_above_half) else '?'}/{int(r.n)}")
        print()
    print(f"[audit_122] -> {OUT}")


if __name__ == "__main__":
    main()
