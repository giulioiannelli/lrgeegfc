#!/usr/bin/env python3
"""audit_134 — N4 v3 / R2: dual-target (learn vs test) replay TRAJECTORY over rest.

Scope: .agents/guides/task-persistence-investigation/2026-06-22_replay-states-v3-scale-target-pair.md

audit_127 nulled the POOLED sustained inference-preference (g_IE: rest_post is no more
test-like than learn-like on average). R2's NEW content is the TEMPORAL TRAJECTORY: even
if the average preference is flat, the replay TARGET may DRIFT across the rest period —
e.g. early rest_post replays the just-finished inference (task_test) config, then
consolidates toward the underlying learned hierarchy (task_learn), or vice versa. That is
a dynamic, in-framework signature the pooled mean cannot see.

Per patient, band (τ₀=1/λmax, from the audit_125 cophenetic cache):
  refs: c_L=coph(task_learn), c_T=coph(task_test), c_P=coph(rest_pre)  [load_fc_matrix]
  per window w:  d(w) = [ρcoph(c_w,c_T) − ρcoph(c_w,c_P)] − [ρcoph(c_w,c_L) − ρcoph(c_w,c_P)]
                      = ρcoph(c_w,c_T) − ρcoph(c_w,c_L)   (>0 ⇒ more test- than learn-like)
  TRAJECTORY : trend = Spearman(window_index, d(w)) over rest_post  (drift of target)
               control: same trend over rest_pre (should be flat)
  T-SPEC BURST: AC1(d(w)) over rest_post vs time-shuffle  (transient snaps to test-specific)
  level (sanity vs audit_127): mean d_post − mean d_pre

Cohort: trend post-vs-pre (two-sided, drift can go either way), LOO; burst dAC1 vs shuffle.
Output: data/audit/replay_states/dual_target_{per_patient,cohort}.csv
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


def coph(W):
    return cophenetic_condensed_from_adjacency(np.clip(np.asarray(W, float), 0.0, 1.0))


def load_windows(pat, band, phase):
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260622)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    t0 = time.time()
    rows = []
    for pat in COHORT:
        for band in BANDS:
            WL = load_fc_matrix(pat, "task_learn", band, fc_method="imcoh_abs")
            WT = load_fc_matrix(pat, "task_test", band, fc_method="imcoh_abs")
            pre = load_windows(pat, band, "rest_pre")
            post = load_windows(pat, band, "rest_post")
            if WL is None or WT is None or pre is None or post is None:
                print(f"  [skip] {pat} {band}: missing inputs"); continue
            cL, cT = coph(WL), coph(WT)

            def d_series(arr):
                return np.array([spearmanr(c, cT)[0] - spearmanr(c, cL)[0] for c in arr])
            d_pre, d_post = d_series(pre), d_series(post)
            tr_post = spearmanr(np.arange(len(d_post)), d_post)[0]
            tr_pre = spearmanr(np.arange(len(d_pre)), d_pre)[0]
            ac_post, p_post = perm_p(d_post, ac1, args.n_perm, rng)
            ac_pre, _ = perm_p(d_pre, ac1, args.n_perm, rng)
            rows.append(dict(patient=pat, band=band,
                             trend_post=float(tr_post), trend_pre=float(tr_pre),
                             trend_diff=float(tr_post - tr_pre),
                             level_post=float(d_post.mean()), level_pre=float(d_pre.mean()),
                             level_diff=float(d_post.mean() - d_pre.mean()),
                             AC1_post=ac_post, AC1_pre=ac_pre, dAC1=ac_post - ac_pre,
                             p_time_post=p_post, n_post=len(d_post)))
            print(f"  {pat} {band}: trend_post={tr_post:+.3f} (pre {tr_pre:+.3f}) "
                  f"level_Δ={d_post.mean()-d_pre.mean():+.3f} AC1_post={ac_post:+.3f}(p={p_post:.3f})")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "dual_target_per_patient.csv", index=False)

    def w2(d, alt):
        d = np.asarray(d, float); d = d[np.isfinite(d)]
        if d.size < 5 or np.allclose(d, 0):
            return np.nan, 0, len(d)
        try:
            return float(wilcoxon(d, alternative=alt).pvalue), int((d > 0).sum()), len(d)
        except ValueError:
            return np.nan, int((d > 0).sum()), len(d)

    coh = []
    for band in BANDS:
        sub = df[df.band == band]
        if sub.empty:
            continue
        # trend: two-sided (drift either way) — report |trend| consistency + both tails
        td = sub.trend_diff.to_numpy()
        p_up, n_up, K = w2(td, "greater")
        p_dn, _, _ = w2(-td, "greater")
        p_b, n_b, _ = w2(sub.dAC1, "greater")
        p_lvl, n_lvl, _ = w2(sub.level_diff, "greater")
        coh.append(dict(band=band, K=K,
                        trend_post_med=float(sub.trend_post.median()),
                        trend_diff_med=float(sub.trend_diff.median()),
                        trend_npos=f"{n_up}/{K}", trend_p_up=p_up, trend_p_down=p_dn,
                        burst_dAC1_med=float(sub.dAC1.median()), burst_npos=f"{n_b}/{K}", burst_p=p_b,
                        n_indiv_bursty=int((sub.p_time_post < 0.05).sum()),
                        level_diff_med=float(sub.level_diff.median()), level_npos=f"{n_lvl}/{K}", level_p=p_lvl))
    cohdf = pd.DataFrame(coh)
    cohdf.to_csv(OUT / "dual_target_cohort.csv", index=False)
    print("\n=== R2 dual-target (test−learn) TRAJECTORY + burst — cohort ===")
    print("  trend>0 ⇒ drifts toward TEST over rest; <0 ⇒ toward LEARN. burst = transient.")
    for _, r in cohdf.iterrows():
        print(f"  [{r.band}] K={r.K}  trend_diff med={r.trend_diff_med:+.3f} {r.trend_npos} "
              f"(p_up={r.trend_p_up:.3f}/p_down={r.trend_p_down:.3f})  "
              f"burst dAC1={r.burst_dAC1_med:+.3f} p={r.burst_p:.3f} indiv={r.n_indiv_bursty}/{r.K}  "
              f"level_Δ={r.level_diff_med:+.3f} p={r.level_p:.3f}")
    print(f"\n[audit_134] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
