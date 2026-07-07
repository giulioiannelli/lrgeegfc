#!/usr/bin/env python3
"""audit_141 — N4 reopened: SUB-SECOND high-γ replay (the ripple regime).

Scope: .agents/guides/task-persistence-investigation/2026-06-22_replay-states-v3-scale-target-pair.md

PI pushback (2026-06-23, CORRECT): the ≥10 s floor is α/β-ONLY. At fs=2048 Hz high γ is
feasible down to 0.2 s (audit_128: n_bins=3 @0.2 s, 14 @0.5 s) — and classic replay is
high-freq ripples (~80–150 Hz) lasting ~50–150 ms, i.e. EXACTLY sub-second high γ. audit_130/
131 stopped gamma at 2 s out of habit, not necessity. This pushes high+low γ to L∈{0.1..2}s
and runs the make-or-break burstiness battery at the ripple timescale.

⚠ Caveat: high γ is NOT the β trace → a hit is a NEW fast-replay phenomenon (the actual
neural-replay regime), not replay OF the trace. Lead to sharpen/kill: 2 s high-γ excess-flash
was 8/10 p=0.065 (audit_131).

Per patient, γ band, window L (nperseg floored at 64 for high-γ feasibility, ~12+ segments):
  c_w = cophenetic of windowed |ImCoh|_γ LRG; canonical full-phase refs c_T, c_{R-}
  s_w = ρcoph(c_w, c_T) − ρcoph(c_w, c_{R-})   (z-scored pooled)
  shift = mean_post − mean_pre                                 (level / N1-analog)
  AC1(s_post) vs TIME-SHUFFLE                                  (clustered replay STATE)
  excess_flash = extra high-task-likeness post windows beyond the uniform shift (lone events)
  PLACEBO = node-permuted task target (specificity: real burst must beat placebo burst)
Cohort one-sided Wilcoxon per (band, L). A (band, L) cell is a replay STATE only if
burst > shuffle AND target-specific (real > placebo).

Output: data/audit/replay_states/subsecond_gamma_{per_patient,cohort}.csv
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.fc.coherence import segment_ffts, imcoh_abs_cube, band_abs_average
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.utils.surrogate.matched_strength import cophenetic_condensed_from_adjacency

OUT = ROOT / "data" / "audit" / "replay_states"
OUT.mkdir(parents=True, exist_ok=True)
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
GAMMA = ["high_gamma", "low_gamma"]
L_GRID = [0.1, 0.2, 0.5, 1.0, 2.0]


def coph(W):
    return cophenetic_condensed_from_adjacency(np.clip(np.asarray(W, float), 0.0, 1.0))


def pick_nperseg(L, fs, target_seg=14, floor=64, cap=1024):
    nps = int(2 * L * fs / (target_seg + 1))
    nps = 2 ** int(np.floor(np.log2(max(floor, nps))))
    return int(min(max(nps, floor), cap))


def ac1(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]; x = x - x.mean()
    d = np.dot(x, x)
    return float(np.dot(x[:-1], x[1:]) / d) if d > 0 and x.size > 3 else 0.0


def perm_p(x, fn, n, rng):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    obs = fn(x)
    null = np.array([fn(rng.permutation(x)) for _ in range(n)])
    return obs, (1 + np.sum(null >= obs)) / (n + 1)


def windowed_cophs(X, fs, L, nperseg, band):
    wlen = int(round(L * fs))
    lo, hi = BRAIN_BANDS[band]
    out = []
    for s in range(0, X.shape[1] - wlen + 1, wlen):
        freqs, ff = segment_ffts(X[:, s:s + wlen], fs, nperseg, nperseg // 2, hi + 5.0)
        if ff.shape[1] < 4:
            continue
        cube = imcoh_abs_cube(ff)
        if cube is None:
            continue
        W = band_abs_average(cube, freqs, lo, hi)
        if W is None:
            continue
        out.append(coph(W))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--bands", nargs="+", default=GAMMA)
    ap.add_argument("--L-grid", nargs="+", type=float, default=L_GRID)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260623)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    t0 = time.time()
    rows = []
    for pat in args.patients:
        fs = FS_OVERRIDES.get(pat, 2048.0)
        try:
            def _load(ph):
                Xx = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
                return Xx if Xx.shape[0] <= Xx.shape[1] else Xx.T
            Xpre, Xpost = _load("rest_pre"), _load("rest_post")
        except Exception as e:
            print(f"[skip] {pat}: {e}"); continue
        for band in args.bands:
            Wt = load_fc_matrix(pat, "task_test", band, fc_method="imcoh_abs")
            Wp = load_fc_matrix(pat, "rest_pre", band, fc_method="imcoh_abs")
            if Wt is None or Wp is None:
                print(f"  [skip] {pat} {band}: refs missing"); continue
            n = Wt.shape[0]
            cT, cP = coph(Wt), coph(Wp)
            perm = rng.permutation(n)
            cTp = coph(np.asarray(Wt, float)[np.ix_(perm, perm)])
            for L in args.L_grid:
                nperseg = pick_nperseg(L, fs)
                cw_pre = windowed_cophs(Xpre, fs, L, nperseg, band)
                cw_post = windowed_cophs(Xpost, fs, L, nperseg, band)
                if len(cw_pre) < 8 or len(cw_post) < 8:
                    continue
                s_pre = np.array([spearmanr(c, cT)[0] - spearmanr(c, cP)[0] for c in cw_pre])
                s_post = np.array([spearmanr(c, cT)[0] - spearmanr(c, cP)[0] for c in cw_post])
                sp_post = np.array([spearmanr(c, cTp)[0] - spearmanr(c, cP)[0] for c in cw_post])
                pooled = np.concatenate([s_pre, s_post]); mu, sd = pooled.mean(), pooled.std()
                if sd > 0:
                    zpre, zpost = (s_pre - mu) / sd, (s_post - mu) / sd
                else:
                    zpre, zpost = s_pre * 0, s_post * 0
                shift = float(zpost.mean() - zpre.mean())
                q95 = np.quantile(zpre, 0.95)
                excess_flash = float(np.mean(zpost > q95) - np.mean((zpre + shift) > q95))
                ac_post, p_post = perm_p(zpost, ac1, args.n_perm, rng)
                ac_pre, _ = perm_p(zpre, ac1, args.n_perm, rng)
                spp = (sp_post - sp_post.mean()) / sp_post.std() if sp_post.std() > 0 else sp_post * 0
                ac_plac = ac1(spp)
                rows.append(dict(patient=pat, band=band, L_sec=L, nperseg=nperseg,
                                 n_win_post=len(cw_post), shift=shift,
                                 AC1_post=ac_post, AC1_pre=ac_pre, AC1_placebo=ac_plac,
                                 p_time_post=p_post, excess_flash=excess_flash,
                                 spread_excess=float(zpost.std() - zpre.std())))
            print(f"  {pat} {band}: done ({time.time()-t0:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "subsecond_gamma_per_patient.csv", index=False)

    def wt(d):
        d = np.asarray(d, float); d = d[np.isfinite(d)]
        if d.size < 5 or np.allclose(d, 0):
            return np.nan, 0, len(d)
        try:
            return float(wilcoxon(d, alternative="greater").pvalue), int((d > 0).sum()), len(d)
        except ValueError:
            return np.nan, int((d > 0).sum()), len(d)

    coh = []
    for band in args.bands:
        for L in args.L_grid:
            sub = df[(df.band == band) & (df.L_sec == L)]
            if sub.empty:
                continue
            p_sh, n_sh, K = wt(sub["shift"])
            dac = (sub.AC1_post - sub.AC1_pre).to_numpy()
            p_b, n_b, _ = wt(dac)
            p_sp, n_sp, _ = wt((sub.AC1_post - sub.AC1_placebo).to_numpy())
            p_ef, n_ef, _ = wt(sub.excess_flash)
            coh.append(dict(band=band, L_sec=L, K=K,
                            med_nwin=int(sub.n_win_post.median()), nperseg=int(sub.nperseg.median()),
                            shift_med=float(sub["shift"].median()), shift_npos=f"{n_sh}/{K}", shift_p=p_sh,
                            burst_dAC1=float(np.nanmedian(dac)), burst_npos=f"{n_b}/{K}", burst_p=p_b,
                            specific_p=p_sp, flash_npos=f"{n_ef}/{K}", flash_p=p_ef,
                            n_indiv_bursty=int((sub.p_time_post < 0.05).sum())))
    cohdf = pd.DataFrame(coh)
    cohdf.to_csv(OUT / "subsecond_gamma_cohort.csv", index=False)
    print("\n=== SUB-SECOND γ replay — cohort (the ripple regime) ===")
    print("  burst = AC1 post>pre vs time-shuffle; specific = real>placebo; flash = lone events")
    print("  * = burst_p<0.05 AND specific_p<0.05  (a genuine replay STATE)")
    for band in args.bands:
        blk = cohdf[cohdf.band == band].sort_values("L_sec")
        if blk.empty:
            continue
        print(f"\n  [{band}]")
        print("   L_sec nps  nwin  shift(npos,p)      burst(npos,p)       spec_p flash(npos,p)  indiv")
        for _, r in blk.iterrows():
            hit = "*" if (r.burst_p < 0.05 and r.specific_p < 0.05) else (
                  "+" if r.flash_p < 0.05 else " ")
            print(f" {hit}{r.L_sec:5.1f} {r.nperseg:4d} {r.med_nwin:4d}  "
                  f"{r.shift_med:+.3f} {r.shift_npos:>5} p={r.shift_p:.3f}  "
                  f"{r.burst_dAC1:+.3f} {r.burst_npos:>5} p={r.burst_p:.3f}  "
                  f"{r.specific_p:.3f}  {r.flash_npos:>5} p={r.flash_p:.3f}  {r.n_indiv_bursty:2d}/{r.K}")
    print(f"\n[audit_141] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
