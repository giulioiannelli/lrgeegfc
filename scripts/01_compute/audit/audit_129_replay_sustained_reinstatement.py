#!/usr/bin/env python3
"""Audit 129 — N4 v2 / Direction A: SUSTAINED, STABILIZED reinstatement of the trace.

Reframe (PI 2026-06-22): the replay/reinstatement is not transient flashes (ruled
out, audit_125-127) but a SUSTAINED STATE — rest_post dwells in the task connectivity
configuration AND more tightly (lower window-to-window dispersion = a tighter
attractor) than rest_pre. Tested in-framework on the 20 s ρ^coph window cache.

Per patient, per band (alpha, beta), comparing rest_pre vs rest_post windows c_w
(cophenetic vectors) against the full-phase task centroid c_test:
  A1 LEVEL      : mean_w rho^coph(c_w, c_test)              post > pre   (= N1, window-level)
  A2 TIGHTNESS  : within-phase dispersion = median pairwise (1 - Spearman) among the
                  phase's windows                            post < pre   (tighter STATE)
  A2b SPREAD    : SD across windows of rho^coph(c_w, c_test) post < pre   (settles on task)

Cohort one-sided Wilcoxon (LOO-max). A1 ≈ N1 (the level); A2/A2b are the new
"it is a STATE, not just a shift" content: rest_post is a tighter task-shaped
attractor. Multiband = average the per-band per-patient metric over {alpha, beta}.

Outputs: data/audit/replay_states/sustained_reinstatement_{per_patient,cohort}.csv
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, wilcoxon

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


def internal_dispersion(arr):
    """median pairwise (1 - Spearman) among rows of arr (n_win, P), via rank+corrcoef."""
    if arr.shape[0] < 3:
        return np.nan
    R = np.vstack([rankdata(row) for row in arr])
    C = np.corrcoef(R)                                  # Spearman matrix
    iu = np.triu_indices(C.shape[0], 1)
    return float(np.median(1.0 - C[iu]))


def main():
    t0 = time.time()
    rows = []
    for pat in COHORT:
        per_band = {}
        ok = True
        for b in BANDS:
            Wt = load_fc_matrix(pat, "task_test", b, fc_method="imcoh_abs")
            Wp = load_fc_matrix(pat, "rest_pre", b, fc_method="imcoh_abs")
            pre = load_cached(pat, b, "rest_pre")
            post = load_cached(pat, b, "rest_post")
            if Wt is None or Wp is None or pre is None or post is None:
                ok = False; break
            ctest = coph(Wt); cpre = coph(Wp)
            # A1 LEVEL = contrast similarity (toward task RELATIVE to pre) — the
            # N1-aligned measure; raw spearman(c, ctest) is swamped by generic
            # window typicality and does NOT isolate the trace.
            sim_pre = np.array([spearmanr(c, ctest)[0] - spearmanr(c, cpre)[0] for c in pre])
            sim_post = np.array([spearmanr(c, ctest)[0] - spearmanr(c, cpre)[0] for c in post])
            per_band[b] = dict(
                A1_pre=sim_pre.mean(), A1_post=sim_post.mean(),
                A2_pre=internal_dispersion(pre), A2_post=internal_dispersion(post),
                A2b_pre=sim_pre.std(), A2b_post=sim_post.std(),
            )
        if not ok:
            print(f"[skip] {pat}"); continue
        # multiband average
        def avg(key):
            return float(np.mean([per_band[b][key] for b in BANDS]))
        A1_pre, A1_post = avg("A1_pre"), avg("A1_post")
        A2_pre, A2_post = avg("A2_pre"), avg("A2_post")
        A2b_pre, A2b_post = avg("A2b_pre"), avg("A2b_post")
        rows.append(dict(patient=pat,
                         A1_level_post_minus_pre=A1_post - A1_pre,
                         A2_tightness_pre_minus_post=A2_pre - A2_post,
                         A2b_settle_pre_minus_post=A2b_pre - A2b_post,
                         A1_pre=A1_pre, A1_post=A1_post,
                         A2_pre=A2_pre, A2_post=A2_post))
        print(f"  {pat}: A1 level Δ={A1_post-A1_pre:+.3f}  A2 tightness(pre-post)={A2_pre-A2_post:+.3f}  "
              f"A2b settle(pre-post)={A2b_pre-A2b_post:+.4f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "sustained_reinstatement_per_patient.csv", index=False)

    def wtest(col):
        d = df[col].to_numpy(); d = d[np.isfinite(d)]
        try:
            p = wilcoxon(d, alternative="greater")[1]
        except Exception:
            p = np.nan
        loo = [wilcoxon(np.delete(d, i), alternative="greater")[1] for i in range(len(d))]
        return float(np.median(d)), int((d > 0).sum()), len(d), float(p), float(np.nanmax(loo))

    coh = []
    for label, col in [("A1 LEVEL post>pre (=N1)", "A1_level_post_minus_pre"),
                       ("A2 TIGHTER post<pre (state)", "A2_tightness_pre_minus_post"),
                       ("A2b SETTLES post<pre (state)", "A2b_settle_pre_minus_post")]:
        med, npos, n, p, loo = wtest(col)
        coh.append(dict(test=label, median=med, n_pos=f"{npos}/{n}", wilcoxon_p=p, loo_max_p=loo))
    cohdf = pd.DataFrame(coh)
    cohdf.to_csv(OUT / "sustained_reinstatement_cohort.csv", index=False)
    print("\n=== DIRECTION A — SUSTAINED, STABILIZED REINSTATEMENT (cohort) ===")
    print(cohdf.to_string(index=False))
    print(f"\n[audit_129] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
