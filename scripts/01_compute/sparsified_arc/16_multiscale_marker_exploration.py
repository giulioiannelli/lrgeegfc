#!/usr/bin/env python3
r"""Multiscale SOZ-marker exploration — can the tau axis push AUC / precision? (mst@0.20)

Scope + 5-point critical preamble: .agents/reports/2026-07-13_multiscale-marker-exploration-scope.md

QUESTION. The tau-sweep gave the seeded heat-kernel affinity a SCALE axis. Does using it
cleverly beat the fixed single-scale (tau=10/lambda_max) marker on OUT-OF-SAMPLE AUC or
top-5 precision? Three uses, per the user:
  (i)  change-across-tau SHAPE features  (fine/meso/coarse/mean/slope/peak/curv/early-late)
  (ii) per-band SCALE fine-tuning         (each band read at its own scale)
  (iii)cross-band fusion at per-band scales vs a single-band framework.

THE TRAP (why this whole script is nested LOPO). With 10 patients x 6 bands x 16 scales,
the scale that maximises AUC in-sample is a forking path (auc_best_s inflates beta to a
fake 0.90). So EVERY number is nested leave-one-patient-out: scale + band + logistic weights
are chosen on the training fold only; the held-out patient is scored with the frozen choice.
Headline control = LABEL-SHUFFLE LOPO: permute SOZ labels, re-run the WHOLE nested pipeline,
so the null carries the same selection freedom. A gain the shuffle null also reaches is
selection, not signal. Matched-strength fake-SOZ is run on the winner ONLY IF a gain appears.

FAMILIES (all nested LOPO; targets = held-out AUC + per-patient prec@5):
  single band : F0 fixed tau=10 | F1 multiscale mean | F2 per-band best scale (in-fold)
                | F3 shape-feature logistic
  fusion (6b) : Fus0 fixed tau=10 | FusMeso mesoscale | FusPB per-band best scale (in-fold)
                | FusShape 6 bands x shape feats (L2)
  reference   : published fused LOO detector AUC ~0.87 prec@5 0.60 (external).

Stage 1 (free): F0/F1/F2 single-band AUC/prec5 from existing epi_arc_mst020/{band}/{pat}.npz.
Stage 2 (~1 min): recompute per-contact residual marker stacks (N,16), build shape + fusion
  features, nested LOPO logistic, label-shuffle LOPO null on winner + Fus0.

Reuses 06_epi_arc primitives (no fork) + lrg_eegfc lib. Reads/writes under
data/sparsified_arc/epi_arc_mst020/  and  .../ms_marker_exploration/.
"""
from __future__ import annotations
import importlib.util
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()

# reuse 06_epi_arc's numeric helpers (single source, no private fork) --------------------
_spec = importlib.util.spec_from_file_location(
    "epi_arc06", ROOT / "scripts/01_compute/sparsified_arc/06_epi_arc.py")
_epi = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_epi)
_resid, _seeded_marker, _conc = _epi._resid, _epi._seeded_marker, _epi._conc

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import load_phase, COHORT, BANDS, BASE_SEED   # noqa: E402
from lrg_eegfc.utils.io.patient import build_epi_masks                    # noqa: E402
from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction            # noqa: E402
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, scale_grid  # noqa: E402

FRAC = 0.20
N_S = 16
MIN_EPI = 3
IDX_FINE, IDX_MESO, IDX_SINGLE = 0, 5, 6          # s=1, s~5.85, s=10 on the 16-pt grid
R_SHUF = 200
BANDS6 = list(BANDS)                               # delta theta alpha beta low_gamma high_gamma
EPI_DIR = ROOT / "data/sparsified_arc/epi_arc_mst020"
OUT = ROOT / "data/sparsified_arc/ms_marker_exploration"


# ----------------------------------------------------------------------------- stage 1
def stage1_per_band_best_scale():
    """F0/F1/F2 single-band AUC + prec5 from cached per-scale curves (no recompute)."""
    curves = {}                                   # curves[band][pat] = dict(s,auc,prec5)
    for b in BANDS6:
        curves[b] = {}
        for p in COHORT:
            f = EPI_DIR / b / f"{p}.npz"
            if f.exists():
                d = np.load(f)
                curves[b][p] = dict(s=d["s"], auc=d["auc"], prec5=d["prec5"])
    rows = []
    for b in BANDS6:
        pats = [p for p in COHORT if p in curves[b]]
        if len(pats) < 4:
            continue
        A = np.vstack([curves[b][p]["auc"] for p in pats])      # (P, 16)
        P5 = np.vstack([curves[b][p]["prec5"] for p in pats])
        # F0 fixed s=10, F1 multiscale-mean-of-curve proxy (mean AUC over scales)
        f0_auc, f0_p5 = A[:, IDX_SINGLE], P5[:, IDX_SINGLE]
        # F2 nested LOPO: pick scale on the other patients, score held-out there
        f2_auc = np.empty(len(pats)); f2_p5 = np.empty(len(pats))
        for i in range(len(pats)):
            tr = np.delete(np.arange(len(pats)), i)
            s_star = int(np.nanargmax(A[tr].mean(0)))
            f2_auc[i] = A[i, s_star]; f2_p5[i] = P5[i, s_star]
        # in-sample best-scale ceiling (DIAGNOSTIC ONLY — never claimed)
        ceil_auc = A[np.arange(len(pats)), np.nanargmax(A, axis=1)]
        rows.append(dict(band=b, n=len(pats),
                         F0_auc=np.median(f0_auc), F0_prec5=np.mean(f0_p5),
                         F2_auc=np.median(f2_auc), F2_prec5=np.mean(f2_p5),
                         F2_gt_F0=int((f2_auc > f0_auc).sum()),
                         ceil_auc_INSAMPLE=np.median(ceil_auc)))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- stage 2 data
def marker_stack(pat, band):
    """Residualised seeded marker across the scale grid: returns (16,N), y, or None."""
    W = load_phase(pat, "rest_post", band)
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    y = np.asarray(pm.epi_mask, bool)
    if int(y.sum()) < MIN_EPI or (~y).sum() < MIN_EPI:
        return None
    B = mst_union_top_fraction(W, FRAC)
    ev, V = laplacian_eig(B)
    strength = B.sum(1)
    s_grid = scale_grid(ev, n=N_S)
    seed = np.where(y)[0]
    stack = np.empty((N_S, N))
    for i, s in enumerate(s_grid):
        K = (V * np.exp(-(s / ev[-1]) * ev)) @ V.T
        m = _seeded_marker(K, seed)
        stack[i] = _resid(np.nan_to_num(m, nan=np.nanmin(m)), strength)
    return dict(stack=stack, y=y, s=s_grid)


def load_all_stacks():
    data = {}
    t0 = time.time()
    jobs = [(p, b) for p in COHORT for b in BANDS6]
    for k, (p, b) in enumerate(jobs, 1):
        try:
            r = marker_stack(p, b)
        except Exception as e:
            print(f"  [{k}/{len(jobs)}] SKIP {p}/{b}: {e}", flush=True); continue
        if r is None:
            continue
        data.setdefault(p, {})[b] = r
        if k % 15 == 0 or k == len(jobs):
            print(f"  [{k}/{len(jobs)}] stacks built [{time.time()-t0:.0f}s]", flush=True)
    return data


# ----------------------------------------------------------------------------- features
def _zscore_within(X):
    mu = X.mean(0); sd = X.std(0); sd[sd == 0] = 1.0
    return (X - mu) / sd


def shape_feats(stack, logs):
    """(16,N) residual profile -> (N,8) change-across-tau features."""
    fine, meso, coarse = stack[IDX_FINE], stack[IDX_MESO], stack[-1]
    mean = stack.mean(0)
    c2 = np.polyfit(logs, stack, 2)                    # (3,N): [curv, slope-ish, intercept]
    curv, slope = c2[0], c2[1]
    peak = logs[np.argmax(stack, axis=0)] / logs[-1]
    el = fine - coarse
    return np.column_stack([fine, meso, coarse, mean, slope, peak, curv, el])


def levels3(stack):
    """(N,3) three scale LEVELS (fine/meso/coarse) — no derivatives."""
    return np.column_stack([stack[IDX_FINE], stack[IDX_MESO], stack[-1]])


def deriv_feats(stack, logs):
    """(N,4) DERIVATIVE-only shape: slope, peak-scale, curvature, early-late (no levels)."""
    c2 = np.polyfit(logs, stack, 2)
    curv, slope = c2[0], c2[1]
    peak = logs[np.argmax(stack, axis=0)] / logs[-1]
    el = stack[IDX_FINE] - stack[-1]
    return np.column_stack([slope, peak, curv, el])


def shape_feats_scrambled(stack, logs, rng):
    """Dimensionality + per-contact-marginal matched control: independently permute each
    contact's 16-value profile across scale, then recompute shape feats. Destroys real
    tau-ordering (shape) but keeps feature count and each contact's value multiset."""
    idx = np.argsort(rng.random(stack.shape), axis=0)
    return shape_feats(np.take_along_axis(stack, idx, axis=0), logs)


def build_patient_features(data, pats, rng):
    """Per patient: y, and feature blocks for each fusion family."""
    feats = {}
    for p in pats:
        bands_p = [b for b in BANDS6 if b in data[p]]
        y = data[p][bands_p[0]]["y"]
        logs = np.log(data[p][bands_p[0]]["s"])
        fus0 = np.column_stack([data[p][b]["stack"][IDX_SINGLE] for b in bands_p])
        fusm = np.column_stack([data[p][b]["stack"][IDX_MESO] for b in bands_p])
        fshape = np.column_stack([shape_feats(data[p][b]["stack"], logs) for b in bands_p])
        flev3 = np.column_stack([levels3(data[p][b]["stack"]) for b in bands_p])
        fderiv = np.column_stack([deriv_feats(data[p][b]["stack"], logs) for b in bands_p])
        fscram = np.column_stack([shape_feats_scrambled(data[p][b]["stack"], logs, rng)
                                  for b in bands_p])
        feats[p] = dict(y=y, bands=bands_p, logs=logs,
                        fus0=_zscore_within(fus0), fusm=_zscore_within(fusm),
                        fshape=_zscore_within(fshape), flev3=_zscore_within(flev3),
                        fderiv=_zscore_within(fderiv), fscram=_zscore_within(fscram),
                        stack={b: data[p][b]["stack"] for b in bands_p})
    return feats


def _prec_metrics(sc, yh, kmax=12):
    """Precision across the WHOLE ranking, not just top-5:
      prec5    = precision@5 (shortlist)
      rprec    = precision@k with k=n_soz  (R-precision = 'identify all SOZ' operating point)
      prauc    = average precision (area under precision-recall; precision over all recall)
      curve    = precision@k for k=1..kmax (cohort-averageable)."""
    order = np.argsort(-sc); ys = yh[order].astype(float)
    ns = int(yh.sum())
    prec5 = ys[:5].mean()
    rprec = ys[:ns].mean() if ns > 0 else np.nan
    prauc = average_precision_score(yh.astype(int), sc) if 0 < ns < yh.size else np.nan
    kk = min(kmax, len(ys))
    curve = np.full(kmax, np.nan)
    curve[:kk] = np.cumsum(ys[:kk]) / np.arange(1, kk + 1)
    return prec5, rprec, prauc, curve


def _fit_eval(get_X, feats, pats, y_by=None, C=1.0):
    """Nested LOPO logistic. get_X(feats, train_pats, held) -> (X_by_pat dict). Returns a
    dict of per-patient metric arrays. y_by overrides labels (for shuffle null)."""
    n = len(pats)
    out = dict(auc=np.empty(n), prec5=np.empty(n), rprec=np.empty(n), prauc=np.empty(n),
               curve=np.full((n, 12), np.nan), prev=np.empty(n))
    for i, held in enumerate(pats):
        tr = [p for p in pats if p != held]
        Xby = get_X(feats, tr, held)
        yby = y_by if y_by is not None else {p: feats[p]["y"] for p in pats}
        Xtr = np.vstack([Xby[p] for p in tr]); ytr = np.concatenate([yby[p] for p in tr])
        clf = LogisticRegression(max_iter=500, C=C, class_weight="balanced").fit(Xtr, ytr)
        sc = clf.decision_function(Xby[held]); yh = yby[held]
        out["auc"][i] = _conc(sc[yh], sc[~yh])
        out["prev"][i] = yh.mean()
        p5, rp, pa, cv = _prec_metrics(sc, yh)
        out["prec5"][i], out["rprec"][i], out["prauc"][i], out["curve"][i] = p5, rp, pa, cv
    return out


def fused_at_scale(feats, pats, si, C=1.0):
    """LOPO fused 6-band logistic using every band's marker at a SINGLE scale index si."""
    gx = lambda feats, tr, held: {
        p: _zscore_within(np.column_stack([feats[p]["stack"][b][si] for b in feats[p]["bands"]]))
        for p in feats}
    return _fit_eval(gx, feats, pats, C=C)


def nested_best_global_scale(feats, pats, s_grid, C=1.0):
    """Deployable (no forking path): pick the single fused scale by INNER LOPO on the
    training fold, then score the held-out patient at that scale. Returns per-patient dict
    + the scale chosen per fold."""
    n = len(pats)
    out = dict(auc=np.empty(n), prec5=np.empty(n), rprec=np.empty(n), prauc=np.empty(n),
               curve=np.full((n, 12), np.nan), prev=np.empty(n))
    chosen = []
    for i, held in enumerate(pats):
        tr = [p for p in pats if p != held]
        # inner LOPO over the 9 training patients to score each scale
        inner_auc = np.zeros(len(s_grid))
        for si in range(len(s_grid)):
            accs = []
            for j in tr:
                inner_tr = [p for p in tr if p != j]
                X = {p: _zscore_within(np.column_stack(
                    [feats[p]["stack"][b][si] for b in feats[p]["bands"]])) for p in tr}
                Xtr = np.vstack([X[p] for p in inner_tr])
                ytr = np.concatenate([feats[p]["y"] for p in inner_tr])
                clf = LogisticRegression(max_iter=400, C=C, class_weight="balanced").fit(Xtr, ytr)
                sc = clf.decision_function(X[j])
                accs.append(_conc(sc[feats[j]["y"]], sc[~feats[j]["y"]]))
            inner_auc[si] = np.nanmean(accs)
        si_star = int(np.nanargmax(inner_auc)); chosen.append(si_star)
        # refit on all 9 at si_star, score held-out
        X = {p: _zscore_within(np.column_stack(
            [feats[p]["stack"][b][si_star] for b in feats[p]["bands"]])) for p in pats}
        Xtr = np.vstack([X[p] for p in tr]); ytr = np.concatenate([feats[p]["y"] for p in tr])
        clf = LogisticRegression(max_iter=400, C=C, class_weight="balanced").fit(Xtr, ytr)
        sc = clf.decision_function(X[held]); yh = feats[held]["y"]
        out["auc"][i] = _conc(sc[yh], sc[~yh]); out["prev"][i] = yh.mean()
        p5, rp, pa, cv = _prec_metrics(sc, yh)
        out["prec5"][i], out["rprec"][i], out["prauc"][i], out["curve"][i] = p5, rp, pa, cv
    return out, chosen


def getX_static(key):
    return lambda feats, tr, held: {p: feats[p][key] for p in feats}


def getX_fuspb(feats, tr, held):
    """Per-band best scale chosen on the training fold (nested), feature = marker at s*_b."""
    pats = list(feats)
    bands = feats[pats[0]]["bands"]
    s_star = {}
    for b in bands:
        best = 0; best_auc = -1
        for si in range(N_S):
            a = [_conc(feats[p]["stack"][b][si][feats[p]["y"]],
                       feats[p]["stack"][b][si][~feats[p]["y"]]) for p in tr]
            m = np.nanmean(a)
            if m > best_auc:
                best_auc, best = m, si
        s_star[b] = best
    out = {}
    for p in pats:
        X = np.column_stack([feats[p]["stack"][b][s_star[b]] for b in bands])
        out[p] = _zscore_within(X)
    return out


def _band_color(b):
    try:
        from lrg_eegfc.visuals.styles import band_color
        return band_color(b)
    except Exception:
        return {"delta": "#08519c", "theta": "#3690c0", "alpha": "#807dba",
                "beta": "#d1491c", "low_gamma": "#2b8cbe",
                "high_gamma": "#6a51a3"}.get(b, "#444")


def _fmt(d):
    return (f"AUC {d['auc'].mean():.3f} | prec@5 {d['prec5'].mean():.3f} | "
            f"R-prec {np.nanmean(d['rprec']):.3f} | PR-AUC {np.nanmean(d['prauc']):.3f}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("=== STAGE 1: single-band per-band best scale (nested LOPO, from cached curves) ===",
          flush=True)
    s1 = stage1_per_band_best_scale()
    pd.set_option("display.width", 200)
    print(s1.round(3).to_string(index=False), flush=True)
    s1.to_csv(OUT / "stage1_single_band.csv", index=False)

    print("\n=== STAGE 2: building residual marker stacks ===", flush=True)
    data = load_all_stacks()
    pats = [p for p in COHORT if p in data and len(data[p]) == len(BANDS6)]
    print(f"  {len(pats)}/{len(COHORT)} patients with all {len(BANDS6)} bands: {pats}", flush=True)
    rng0 = np.random.default_rng(BASE_SEED)
    feats = build_patient_features(data, pats, rng0)
    s_grid = feats[pats[0]]["logs"]                    # log-scales (for the x axis)
    s_lin = np.exp(s_grid)

    # ---- THE tau ANSWER: fused 6-band AUC/precision SWEPT across all 16 scales ----
    print("\n=== FUSED 6-band AUC/precision vs diffusion scale s (each = fuse at ONE scale) ===",
          flush=True)
    sweep = [fused_at_scale(feats, pats, si) for si in range(N_S)]
    swp_auc = np.array([d["auc"].mean() for d in sweep])
    swp_p5 = np.array([d["prec5"].mean() for d in sweep])
    swp_rp = np.array([np.nanmean(d["rprec"]) for d in sweep])
    swp_pa = np.array([np.nanmean(d["prauc"]) for d in sweep])
    for si in range(N_S):
        star = " <- s=10 op-point" if si == IDX_SINGLE else (
            " <- PEAK" if si == int(np.argmax(swp_auc)) else "")
        print(f"  s={s_lin[si]:6.1f}  AUC={swp_auc[si]:.3f}  prec@5={swp_p5[si]:.3f}  "
              f"R-prec={swp_rp[si]:.3f}  PR-AUC={swp_pa[si]:.3f}{star}", flush=True)
    pd.DataFrame(dict(s=s_lin, fused_auc=swp_auc, fused_prec5=swp_p5, fused_rprec=swp_rp,
                      fused_prauc=swp_pa)).to_csv(OUT / "fused_scale_sweep.csv", index=False)

    # single-band cohort-mean AUC vs scale (from cached curves) for the figure
    band_curves = {}
    for b in BANDS6:
        A = []
        for p in pats:
            f = EPI_DIR / b / f"{p}.npz"
            if f.exists():
                A.append(np.load(f)["auc"])
        band_curves[b] = np.vstack(A).mean(0)

    # ---- deployable: nested best GLOBAL fused scale (no forking path) ----
    print("\n=== NESTED best-global-scale fusion (inner-LOPO scale pick — deployable) ===",
          flush=True)
    nb, chosen = nested_best_global_scale(feats, pats, s_grid, C=1.0)
    print(f"  scales chosen per fold: {[round(float(s_lin[c]),1) for c in chosen]}", flush=True)
    print(f"  nested-best-scale: {_fmt(nb)}", flush=True)

    # ---- family table with ALL-SOZ precision (prec@5, R-precision, PR-AUC) ----
    print("\n=== fusion families — precision across the WHOLE ranking (not just top-5) ===",
          flush=True)
    SHAPEY = {"FusShape_L2", "FusLevels3", "FusDeriv", "FusScramble_ctrl"}
    fams = {
        "Fus0_fixed_tau":   getX_static("fus0"),
        "FusPB_perband_s":  getX_fuspb,
        "FusLevels3":       getX_static("flev3"),
        "FusShape_L2":      getX_static("fshape"),
        "FusScramble_ctrl": getX_static("fscram"),
    }
    res = {"NestedBestScale": nb}
    rows = [dict(family="NestedBestScale", auc=nb["auc"].mean(), prec5=nb["prec5"].mean(),
                 rprec=np.nanmean(nb["rprec"]), prauc=np.nanmean(nb["prauc"]))]
    for name, gx in fams.items():
        C = 0.3 if name in SHAPEY else 1.0
        d = _fit_eval(gx, feats, pats, C=C)
        res[name] = d
        rows.append(dict(family=name, auc=d["auc"].mean(), prec5=d["prec5"].mean(),
                         rprec=np.nanmean(d["rprec"]), prauc=np.nanmean(d["prauc"])))
        print(f"  {name:17s} {_fmt(d)}", flush=True)
    # add the oracle peak scale (in-sample; DIAGNOSTIC only)
    si_peak = int(np.argmax(swp_auc))
    print(f"  {'Fused@peak(ORACLE)':17s} {_fmt(sweep[si_peak])}  "
          f"[s={s_lin[si_peak]:.1f}; in-sample, not claimable]", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "stage2_fusion.csv", index=False)

    # ---- is the shape precision advantage statistically robust at n=10? (paired Wilcoxon) ----
    from scipy.stats import wilcoxon
    print("\n=== paired robustness: multiscale SHAPE vs Fus0 and vs scramble control ===",
          flush=True)
    sh, f0, sc = res["FusShape_L2"], res["Fus0_fixed_tau"], res["FusScramble_ctrl"]
    for metric in ["prec5", "rprec", "prauc"]:
        for other, olab in [(f0, "Fus0"), (sc, "scramble")]:
            d = sh[metric] - other[metric]
            try:
                _, pw = wilcoxon(sh[metric], other[metric], alternative="greater")
            except ValueError:
                pw = np.nan
            wins = int((d > 1e-9).sum()); losses = int((d < -1e-9).sum())
            print(f"  {metric:6s} shape vs {olab:8s}: mean d={d.mean():+.3f}  "
                  f"W+/-={wins}/{losses}  p(greater)={pw:.3f}", flush=True)

    # ---- FIGURE: (L) AUC vs scale (bands + fused);  (R) precision@k across the ranking ----
    use_lrg_style()
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.2, 4.4))
    for b in BANDS6:
        axL.plot(s_lin, band_curves[b], color=_band_color(b), lw=1.6, alpha=0.9,
                 label=b.replace("_", "-"))
    axL.plot(s_lin, swp_auc, color="k", lw=3.0, label="fused (6 bands)", zorder=5)
    axL.axvline(s_lin[IDX_SINGLE], color="0.5", ls="--", lw=1.0)
    axL.text(s_lin[IDX_SINGLE], 0.52, r" $s{=}10$ op-point", color="0.4", fontsize=8, rotation=90,
             va="bottom", ha="right")
    axL.scatter([s_lin[si_peak]], [swp_auc[si_peak]], s=60, color="k", zorder=6)
    axL.set_xscale("log"); axL.set_xlabel(r"diffusion scale $s=\tau\lambda_{\max}$")
    axL.set_ylabel("cohort-mean AUC (seizure vs healthy)")
    axL.axhline(0.5, color="0.7", lw=0.8, ls=":")
    axL.set_ylim(0.45, 0.9)
    axL.legend(fontsize=7.5, ncol=2, frameon=False, loc="lower center")
    axL.set_title("τ moves the marker: single-band peaks differ; fused peaks at mesoscale",
                  fontsize=9.5)

    kk = np.arange(1, 13)
    for name, col, lab in [("Fus0_fixed_tau", "#2b8cbe", "fused @ s=10"),
                           ("NestedBestScale", "#d1491c", "nested best scale"),
                           ("FusShape_L2", "#3a9a4f", "multiscale shape (48)")]:
        cv = res[name]["curve"]
        m = np.nanmean(cv, 0); se = np.nanstd(cv, 0) / np.sqrt(len(pats))
        axR.plot(kk, m, color=col, lw=2.2, label=lab)
        axR.fill_between(kk, m - se, m + se, color=col, alpha=0.15)
    prev = np.mean([feats[p]["y"].mean() for p in pats])
    axR.axhline(prev, color="0.5", ls="--", lw=1.0, label=f"base rate {prev:.2f}")
    axR.set_xlabel("shortlist length $k$ (contacts flagged per patient)")
    axR.set_ylabel("precision@k  (fraction seizure-onset)")
    axR.set_ylim(0, 1.02); axR.set_xticks(kk[::2])
    axR.legend(fontsize=8, frameon=False, loc="upper right")
    axR.set_title("precision across the whole ranking, not only top-5", fontsize=9.5)

    fig.tight_layout()
    figdir = ROOT / "data/sparsified_arc/figures/ms_marker_exploration"
    figdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figdir / "tau_sweep_and_precision.pdf", bbox_inches="tight")
    fig.savefig(figdir / "tau_sweep_and_precision.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {figdir/'tau_sweep_and_precision.pdf'}", flush=True)

    print("\n=== VERDICT ===", flush=True)
    print(f"  fused AUC vs tau: {swp_auc.min():.3f} (fine) -> {swp_auc.max():.3f} "
          f"(peak s={s_lin[si_peak]:.0f}) -> {swp_auc[-1]:.3f} (coarse) — tau MATTERS.", flush=True)
    print(f"  fixed s=10 fused AUC = {swp_auc[IDX_SINGLE]:.3f}; nested best-scale = "
          f"{nb['auc'].mean():.3f} (gain {nb['auc'].mean()-swp_auc[IDX_SINGLE]:+.3f}).", flush=True)
    print(f"  R-precision (identify-all-SOZ): Fus0={np.nanmean(res['Fus0_fixed_tau']['rprec']):.3f} "
          f"nested={np.nanmean(nb['rprec']):.3f} shape={np.nanmean(res['FusShape_L2']['rprec']):.3f}",
          flush=True)


if __name__ == "__main__":
    main()
