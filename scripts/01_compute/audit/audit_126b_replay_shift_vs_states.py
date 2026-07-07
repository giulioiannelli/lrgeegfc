#!/usr/bin/env python3
"""Audit 126b — N4: is the occupancy excess STATES or just N1's stationary SHIFT?

audit_126 (nearest-centroid) found task-state occupancy higher in rest_post than
rest_pre (p=0.005, 8/10, LOO-robust). But a UNIFORM mean shift toward task (= N1's
stationary trace) mechanically raises occupancy too — no discrete states required.
This script separates the two, per patient, on the multiband NC discriminant
g_w = rho^coph(c_w, c_task) - rho^coph(c_w, c_pre):

  (1) STATIONARY SHIFT (= N1):   mean(g_post) - mean(g_pre)            > 0 ?
  (2) EXCESS SPREAD (states):    std(g_post) - std(g_pre)             > 0 ?
      (discrete states make some windows snap to task -> heavier tail/variance)
  (3) OCCUPANCY BEYOND SHIFT:    occ_post - occ_pre_SHIFTED            > 0 ?
      occ_pre_SHIFTED = occupancy of rest_pre AFTER adding the (mean_post-mean_pre)
      shift. If the occupancy excess is a pure uniform shift, this is ~0; if there
      is extra state structure, it is > 0.

Verdict logic: shift>0 + spread~0 + beyond-shift~0  =>  stationary trace (N1),
NOT discrete replay states. Reuses audit_125 cophenetic cache.
"""
from __future__ import annotations

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


def coph(W):
    return cophenetic_condensed_from_adjacency(np.clip(np.asarray(W, float), 0.0, 1.0))


def load_cached(pat, band, phase):
    f = CACHE / f"{pat}_{band}_{phase}_coph.npz"
    return np.load(f)["coph"].astype(np.float64) if f.exists() else None


def main():
    t0 = time.time()
    rows = []
    for pat in COHORT:
        refs, ok = {}, True
        for b in BANDS:
            Wp = load_fc_matrix(pat, "rest_pre", b, fc_method="imcoh_abs")
            Wt = load_fc_matrix(pat, "task_test", b, fc_method="imcoh_abs")
            if Wp is None or Wt is None:
                ok = False; break
            refs[b] = (coph(Wp), coph(Wt))
        A = {b: {ph: load_cached(pat, b, ph) for ph in ("rest_pre", "rest_post")} for b in BANDS}
        if not ok or any(A[b][ph] is None for b in BANDS for ph in A[b]):
            print(f"[skip] {pat}"); continue

        # multiband NC discriminant, z-scored per band on pooled pre+post
        gm = {}
        for ph in ("rest_pre", "rest_post"):
            zb = []
            for b in BANDS:
                cpre, ctask = refs[b]
                gb = np.array([spearmanr(c, ctask)[0] - spearmanr(c, cpre)[0] for c in A[b][ph]])
                pooled = np.concatenate([
                    [spearmanr(c, ctask)[0] - spearmanr(c, cpre)[0] for c in A[b]["rest_pre"]],
                    [spearmanr(c, ctask)[0] - spearmanr(c, cpre)[0] for c in A[b]["rest_post"]],
                ])
                mu, sd = pooled.mean(), pooled.std()
                zb.append((gb - mu) / sd if sd > 0 else gb * 0)
            gm[ph] = np.mean(zb, axis=0)

        gpre, gpost = gm["rest_pre"], gm["rest_post"]
        theta = np.median(np.concatenate([gpre, gpost]))
        occ_pre = float(np.mean(gpre > theta))
        occ_post = float(np.mean(gpost > theta))
        shift = float(gpost.mean() - gpre.mean())
        occ_pre_shifted = float(np.mean((gpre + shift) > theta))
        rows.append(dict(patient=pat,
                         mean_shift=shift,
                         spread_excess=float(gpost.std() - gpre.std()),
                         std_ratio=float(gpost.std() / gpre.std()) if gpre.std() > 0 else np.nan,
                         occ_pre=occ_pre, occ_post=occ_post,
                         docc=occ_post - occ_pre,
                         occ_beyond_shift=occ_post - occ_pre_shifted))
        print(f"  {pat}: shift={shift:+.2f} spread_excess={gpost.std()-gpre.std():+.2f} "
              f"(ratio {gpost.std()/gpre.std():.2f})  Δocc={occ_post-occ_pre:+.2f} "
              f"beyond_shift={occ_post-occ_pre_shifted:+.2f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "shift_vs_states.csv", index=False)

    def w(col, alt="greater"):
        d = df[col].to_numpy(); d = d[np.isfinite(d)]
        try:
            p = wilcoxon(d, alternative=alt)[1]
        except Exception:
            p = np.nan
        return float(np.median(d)), int((d > 0).sum()), len(d), float(p)

    print("\n=== SHIFT vs STATES cohort verdict (n={}) ===".format(len(df)))
    for label, col, alt in [
        ("(1) stationary mean shift toward task (= N1)", "mean_shift", "greater"),
        ("(2) excess spread in rest_post (states)", "spread_excess", "greater"),
        ("(3) occupancy excess BEYOND uniform shift (states)", "occ_beyond_shift", "greater"),
        ("    raw occupancy excess docc", "docc", "greater"),
    ]:
        med, npos, n, p = w(col, alt)
        print(f"  {label:52s}: median={med:+.3f}  n_pos={npos}/{n}  Wilcoxon_p={p:.4f}")
    print(f"\n[audit_126b] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
