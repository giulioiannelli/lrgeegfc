#!/usr/bin/env python3
"""Audit 127 — N4 replay states via FOUR-PHASE rho^coph contrasts (in-framework).

Stays entirely in the LRG cophenetic (rho^coph) on |ImCoh| time windows — the N1
substrate (windowed recipe == N1, anchored audit_124). The refinement (PI
directive 2026-06-22): don't template-match windows to a single task target;
leverage the FOUR-PHASE design — differences between task_learn (encoding) and
task_test (inference), and between task_test and rest_pre/post (trace). Reuses the
per-window cophenetic cache from audit_125 (rest_pre, rest_post windows) + the
full-phase references c_pre, c_learn, c_test (cached imcoh_abs -> cophenetic).

Per rest window w (cophenetic c_w), with rho = Spearman over pairs:
  a_infer(w) = rho(c_w - c_learn, c_test - c_learn)   # projection on encoding->inference axis
  a_trace(w) = rho(c_w - c_pre,   c_test - c_pre)      # projection on pre->test (trace) axis
  g_IE(w)    = rho(c_w, c_test) - rho(c_w, c_learn)    # inference-vs-encoding similarity
  g_TP(w)    = rho(c_w, c_test) - rho(c_w, c_pre)      # trace similarity
Multiband = mean over {alpha, beta} of per-band z-scored (pooled pre+post) detector.

The headline four-phase question: does rest_post reinstate the INFERENCE (test)
configuration PREFERENTIALLY over encoding (learn), more than rest_pre — and is it
carried by discrete STATES (burstiness) or a stationary SHIFT?

Battery per detector (each the right test for its claim):
  - rest_post vs rest_pre: mean + occupancy (Wilcoxon cohort, LOO)              [asymmetry]
  - burstiness: AC1 vs time-shuffle null (per-patient p; cohort post>0, post>pre) [states]
  - shift-vs-states: mean_shift, spread_excess, occupancy beyond uniform shift   [states vs N1]

Outputs
-------
    data/audit/replay_states/fourphase_per_patient.csv
    data/audit/replay_states/fourphase_cohort.csv
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.utils.surrogate.matched_strength import cophenetic_condensed_from_adjacency
from lrg_eegfc.workflow.fc import load_fc_matrix

OUT = ROOT / "data" / "audit" / "replay_states"
CACHE = OUT / "cache"
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["alpha", "beta"]
DETECTORS = ["a_infer", "a_trace", "g_IE", "g_TP"]


def coph(W):
    return cophenetic_condensed_from_adjacency(np.clip(np.asarray(W, float), 0.0, 1.0))


def load_cached(pat, band, phase):
    f = CACHE / f"{pat}_{band}_{phase}_coph.npz"
    return np.load(f)["coph"].astype(np.float64) if f.exists() else None


def ac1(x):
    x = np.asarray(x, float); x = x - x.mean()
    d = np.dot(x, x)
    return float(np.dot(x[:-1], x[1:]) / d) if d > 0 else 0.0


def perm_p(x, fn, n, rng):
    obs = fn(x)
    null = np.array([fn(rng.permutation(x)) for _ in range(n)])
    return obs, (1 + np.sum(null >= obs)) / (n + 1)


def detector_values(arr, refs, name):
    """Compute one detector's per-window values for a cophenetic array (n_win, P)."""
    cpre, clearn, ctest = refs
    out = []
    for c in arr:
        if name == "a_infer":
            out.append(spearmanr(c - clearn, ctest - clearn)[0])
        elif name == "a_trace":
            out.append(spearmanr(c - cpre, ctest - cpre)[0])
        elif name == "g_IE":
            out.append(spearmanr(c, ctest)[0] - spearmanr(c, clearn)[0])
        elif name == "g_TP":
            out.append(spearmanr(c, ctest)[0] - spearmanr(c, cpre)[0])
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=20260622)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    t0 = time.time()
    rows = []
    for pat in args.patients:
        refs = {}
        ok = True
        for b in BANDS:
            Wp = load_fc_matrix(pat, "rest_pre", b, fc_method="imcoh_abs")
            Wl = load_fc_matrix(pat, "task_learn", b, fc_method="imcoh_abs")
            Wt = load_fc_matrix(pat, "task_test", b, fc_method="imcoh_abs")
            if any(W is None for W in (Wp, Wl, Wt)):
                ok = False; break
            refs[b] = (coph(Wp), coph(Wl), coph(Wt))
        A = {b: {ph: load_cached(pat, b, ph) for ph in ("rest_pre", "rest_post")} for b in BANDS}
        if not ok or any(A[b][ph] is None for b in BANDS for ph in A[b]):
            print(f"[skip] {pat}"); continue

        for det in DETECTORS:
            # per-band raw, then z-score on pooled pre+post, average bands
            z = {"rest_pre": [], "rest_post": []}
            for b in BANDS:
                vpre = detector_values(A[b]["rest_pre"], refs[b], det)
                vpost = detector_values(A[b]["rest_post"], refs[b], det)
                pooled = np.concatenate([vpre, vpost])
                mu, sd = pooled.mean(), pooled.std()
                z["rest_pre"].append((vpre - mu) / sd if sd > 0 else vpre * 0)
                z["rest_post"].append((vpost - mu) / sd if sd > 0 else vpost * 0)
            gpre = np.mean(z["rest_pre"], axis=0)
            gpost = np.mean(z["rest_post"], axis=0)

            theta = np.median(np.concatenate([gpre, gpost]))
            occ_pre = float(np.mean(gpre > theta)); occ_post = float(np.mean(gpost > theta))
            shift = float(gpost.mean() - gpre.mean())
            occ_pre_shift = float(np.mean((gpre + shift) > theta))
            ac_post, p_post = perm_p(gpost, ac1, args.n_perm, rng)
            ac_pre, _ = perm_p(gpre, ac1, args.n_perm, rng)
            rows.append(dict(
                patient=pat, detector=det,
                mean_pre=float(gpre.mean()), mean_post=float(gpost.mean()), mean_shift=shift,
                occ_pre=occ_pre, occ_post=occ_post, docc=occ_post - occ_pre,
                occ_beyond_shift=occ_post - occ_pre_shift,
                spread_excess=float(gpost.std() - gpre.std()),
                AC1_pre=ac_pre, AC1_post=ac_post, p_time_post=p_post,
            ))
        rd = {r["detector"]: r for r in rows if r["patient"] == pat}
        print(f"  {pat}: " + "  ".join(
            f"{d}:Δμ{rd[d]['mean_shift']:+.2f}/Δocc{rd[d]['docc']:+.2f}/AC1p{rd[d]['p_time_post']:.2f}"
            for d in DETECTORS))

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "fourphase_per_patient.csv", index=False)

    # ---- cohort verdict ----
    def wtest(sub, col):
        d = sub[col].to_numpy(); d = d[np.isfinite(d)]
        if len(d) < 5:
            return np.nan, 0, len(d), np.nan, np.nan
        try:
            p = wilcoxon(d, alternative="greater")[1]
        except Exception:
            p = np.nan
        loo = []
        for i in range(len(d)):
            try:
                loo.append(wilcoxon(np.delete(d, i), alternative="greater")[1])
            except Exception:
                loo.append(np.nan)
        return float(np.median(d)), int((d > 0).sum()), len(d), float(p), float(np.nanmax(loo)) if loo else np.nan

    coh = []
    for det in DETECTORS:
        sub = df[df.detector == det]
        sub = sub.assign(dAC1=sub.AC1_post - sub.AC1_pre)
        for metric, col in [("shift_post_gt_pre(=N1)", "mean_shift"),
                            ("occupancy_post_gt_pre", "docc"),
                            ("occ_BEYOND_shift(states)", "occ_beyond_shift"),
                            ("spread_excess(states)", "spread_excess"),
                            ("burstiness_AC1_post_gt_pre(states)", "dAC1")]:
            med, npos, n, p, loo = wtest(sub, col)
            coh.append(dict(detector=det, metric=metric, median=med,
                            n_pos=f"{npos}/{n}", wilcoxon_p=p, loo_max_p=loo,
                            n_sig_burst=int((sub.p_time_post < 0.05).sum()) if metric.startswith("burstiness") else ""))
    cohdf = pd.DataFrame(coh)
    cohdf.to_csv(OUT / "fourphase_cohort.csv", index=False)
    pd.set_option("display.width", 160)
    print("\n=== FOUR-PHASE COHORT VERDICT ===")
    print(cohdf.to_string(index=False))
    print(f"\n[audit_127] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
