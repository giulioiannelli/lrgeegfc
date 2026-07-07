#!/usr/bin/env python3
"""Audit 130 — N4 v2 / Direction B2: SHORT-WINDOW (fast) replay test in GAMMA.

The short-window precondition (audit_128) shows alpha/beta are unmeasurable below
~10 s (band too narrow vs coarse frequency resolution); only GAMMA survives at
short windows. So the only feasible "fast replay" test (the 2 s / sub-second
timescale the PI asked for) is in gamma — which is NOT a trace band, so this is an
honest exploration of whether a FAST, ripple-timescale reinstatement event exists
at all, even outside the slow trace bands.

Per short window L (default 2 s), bands {low_gamma, high_gamma}, nperseg chosen for
~12-16 segments at that L: windowed |ImCoh|_gamma -> LRG cophenetic c_w. References
c_test, c_pre = cophenetic of FULL-PHASE task/pre gamma |ImCoh| computed at the SAME
nperseg (consistent spectral resolution). Task-likeness contrast
  g(w) = rho^coph(c_w, c_test) - rho^coph(c_w, c_pre)
multiband z-scored mean over the two gamma bands. Then the make-or-break battery:
  - burstiness: AC1(g) vs time-shuffle null (per-patient p)   [fast STATES]
  - rest_post vs rest_pre: mean shift + occupancy             [consolidation]
  - shift-vs-states: spread excess, occupancy beyond shift    [states vs uniform]

Many short windows per phase (e.g. 2 s -> ~300 windows in 10 min rest) gives strong
temporal statistics; per-window reliability is low (precondition) but the
time-shuffle null accounts for it.

Outputs: data/audit/replay_states/shortwindow_gamma_{per_patient,cohort}.csv
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
from lrg_eegfc.utils.fc.coherence import (
    compute_imcoh, segment_ffts, imcoh_abs_cube, band_abs_average,
)
from lrg_eegfc.utils.surrogate.matched_strength import cophenetic_condensed_from_adjacency

OUT = ROOT / "data" / "audit" / "replay_states"
OUT.mkdir(parents=True, exist_ok=True)
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
GAMMA = ["low_gamma", "high_gamma"]


def coph(W):
    return cophenetic_condensed_from_adjacency(np.clip(np.asarray(W, float), 0.0, 1.0))


def pick_nperseg(L, fs, target_seg=12, cap=4096):
    nps = int(2 * L * fs / (target_seg + 1))
    nps = 2 ** int(np.floor(np.log2(max(16, nps))))
    return min(nps, cap)


def full_phase_coph(X, fs, lo, hi, nperseg):
    freqs, signed = compute_imcoh(X, fs, nperseg=nperseg, noverlap=nperseg // 2)
    mask = (freqs >= lo) & (freqs <= hi)
    W = np.abs(signed[:, :, mask]).mean(-1)
    W = 0.5 * (W + W.T); np.fill_diagonal(W, 0.0)
    return coph(W)


def ac1(x):
    x = np.asarray(x, float); x = x - x.mean()
    d = np.dot(x, x)
    return float(np.dot(x[:-1], x[1:]) / d) if d > 0 else 0.0


def perm_p(x, fn, n, rng):
    obs = fn(x)
    null = np.array([fn(rng.permutation(x)) for _ in range(n)])
    return obs, (1 + np.sum(null >= obs)) / (n + 1)


def windowed_cophs(X, fs, L, nperseg, bands):
    wlen = int(round(L * fs))
    starts = range(0, X.shape[1] - wlen + 1, wlen)
    fmax = max(BRAIN_BANDS[b][1] for b in bands) + 5.0
    out = {b: [] for b in bands}
    for s in starts:
        freqs, ff = segment_ffts(X[:, s:s + wlen], fs, nperseg, nperseg // 2, fmax)
        if ff.shape[1] < 4:
            continue
        cube = imcoh_abs_cube(ff)
        for b in bands:
            W = band_abs_average(cube, freqs, *BRAIN_BANDS[b])
            out[b].append(coph(W))
    return {b: np.array(v) for b, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=["Pat_05", "Pat_02", "Pat_13"])
    ap.add_argument("--L-sec", type=float, default=2.0)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260622)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    t0 = time.time()
    rows = []
    for pat in args.patients:
        fs = FS_OVERRIDES.get(pat, 2048.0)
        nperseg = pick_nperseg(args.L_sec, fs)
        print(f"\n=== {pat} (fs={fs:.0f}, L={args.L_sec}s, nperseg={nperseg}) ===")
        # load each phase ONCE
        def _load(ph):
            X = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
            return X if X.shape[0] <= X.shape[1] else X.T
        try:
            Xpre = _load("rest_pre"); Xtask = _load("task_test"); Xpost = _load("rest_post")
        except Exception as e:
            print(f"  [load fail] {e}"); continue
        # references at matching nperseg
        refs = {}
        for b in GAMMA:
            lo, hi = BRAIN_BANDS[b]
            refs[b] = (full_phase_coph(Xpre, fs, lo, hi, nperseg),
                       full_phase_coph(Xtask, fs, lo, hi, nperseg))
        # windowed cophenetics for rest phases
        cw = {}
        for ph, X in (("rest_pre", Xpre), ("rest_post", Xpost)):
            cw[ph] = windowed_cophs(X, fs, args.L_sec, nperseg, GAMMA)
            print(f"  {ph}: {len(next(iter(cw[ph].values())))} windows")

        # multiband contrast task-likeness g, z-scored pooled
        g = {}
        for ph in ("rest_pre", "rest_post"):
            zb = []
            for b in GAMMA:
                cpre, ctest = refs[b]
                vals = np.array([spearmanr(c, ctest)[0] - spearmanr(c, cpre)[0] for c in cw[ph][b]])
                pooled = np.concatenate([
                    [spearmanr(c, ctest)[0] - spearmanr(c, cpre)[0] for c in cw["rest_pre"][b]],
                    [spearmanr(c, ctest)[0] - spearmanr(c, cpre)[0] for c in cw["rest_post"][b]],
                ])
                mu, sd = pooled.mean(), pooled.std()
                zb.append((vals - mu) / sd if sd > 0 else vals * 0)
            g[ph] = np.mean(zb, axis=0)

        gpre, gpost = g["rest_pre"], g["rest_post"]
        theta = np.median(np.concatenate([gpre, gpost]))
        shift = float(gpost.mean() - gpre.mean())
        occ_pre = float(np.mean(gpre > theta)); occ_post = float(np.mean(gpost > theta))
        occ_pre_shift = float(np.mean((gpre + shift) > theta))
        ac_post, p_post = perm_p(gpost, ac1, args.n_perm, rng)
        ac_pre, _ = perm_p(gpre, ac1, args.n_perm, rng)
        rows.append(dict(patient=pat, L_sec=args.L_sec, nperseg=nperseg,
                         n_win_pre=len(gpre), n_win_post=len(gpost),
                         mean_shift=shift, docc=occ_post - occ_pre,
                         occ_beyond_shift=occ_post - occ_pre_shift,
                         spread_excess=float(gpost.std() - gpre.std()),
                         AC1_pre=ac_pre, AC1_post=ac_post, p_time_post=p_post))
        print(f"  shift={shift:+.3f} Δocc={occ_post-occ_pre:+.3f} beyond={occ_post-occ_pre_shift:+.3f} "
              f"spread_exc={gpost.std()-gpre.std():+.3f}  AC1 post={ac_post:+.3f}(p={p_post:.3f}) pre={ac_pre:+.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "shortwindow_gamma_per_patient.csv", index=False)

    if len(df) >= 5:
        def wtest(col):
            d = df[col].to_numpy(); d = d[np.isfinite(d)]
            try:
                p = wilcoxon(d, alternative="greater")[1]
            except Exception:
                p = np.nan
            return float(np.median(d)), int((d > 0).sum()), len(d), float(p)
        df["dAC1"] = df.AC1_post - df.AC1_pre
        coh = []
        for label, col in [("shift_post_gt_pre", "mean_shift"),
                           ("occupancy_post_gt_pre", "docc"),
                           ("occ_beyond_shift(states)", "occ_beyond_shift"),
                           ("spread_excess(states)", "spread_excess"),
                           ("burstiness_AC1_post_gt_pre(states)", "dAC1")]:
            med, npos, n, p = wtest(col)
            coh.append(dict(metric=label, median=med, n_pos=f"{npos}/{n}", wilcoxon_p=p))
        cohdf = pd.DataFrame(coh)
        cohdf.to_csv(OUT / "shortwindow_gamma_cohort.csv", index=False)
        print("\n=== SHORT-WINDOW GAMMA cohort verdict ===")
        print(cohdf.to_string(index=False))
        print(f"  n_patients individually bursty (p_time<0.05): {int((df.p_time_post<0.05).sum())}/{len(df)}")
    print(f"\n[audit_130] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
