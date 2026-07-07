#!/usr/bin/env python3
"""Audit 131 — N4 v2 / Direction B2b: fast replay in GAMMA, LOW vs HIGH, two detectors.

Answers "are there fast replay STATES in gamma, and high or low?". Extends audit_130
by (a) keeping low_gamma and high_gamma SEPARATE, and (b) adding an ISOLATED-FLASH
detector, because at 2 s windows the AC1 (clustering) test is blind to lone
single-window events (one 2 s flash = one high window, not autocorrelated).

Per gamma band, L=2 s, references at matching nperseg (full-phase task/pre cophenetic):
  g(w) = rho^coph(c_w, c_test) - rho^coph(c_w, c_pre)   (z-scored, pooled pre+post)
Two state detectors, both vs the rest_pre within-subject control:
  - CLUSTERING:  AC1(g_post) vs time-shuffle null -> p_time      (runs of states)
  - EXCESS FLASHES: extra high-task-likeness windows in rest_post BEYOND the uniform
      shift:  mean(g_post > q95_pre) - mean((g_pre + shift) > q95_pre)
      (>0 ⇒ heavier upper tail in post = isolated reinstatement events; the
      shift-correction removes the trivial "post is uniformly higher" effect)
Also report the stationary shift + spread for context.

Outputs: data/audit/replay_states/gamma_perband_{per_patient,cohort}.csv
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
    return min(2 ** int(np.floor(np.log2(max(16, nps)))), cap)


def _band_abs(signed, freqs, lo, hi):
    mask = (freqs >= lo) & (freqs <= hi)
    W = np.abs(signed[:, :, mask]).mean(-1)
    W = 0.5 * (W + W.T); np.fill_diagonal(W, 0.0)
    return W


def full_phase_signed(X, fs, nperseg, fmax):
    """Full-phase signed ImCoh once (band-independent); band-average per band after."""
    freqs, signed = compute_imcoh(X, fs, nperseg=nperseg, noverlap=nperseg // 2)
    keep = freqs <= fmax
    return freqs[keep], signed[:, :, keep]


def ac1(x):
    x = np.asarray(x, float); x = x - x.mean()
    d = np.dot(x, x)
    return float(np.dot(x[:-1], x[1:]) / d) if d > 0 else 0.0


def perm_p(x, fn, n, rng):
    obs = fn(x)
    null = np.array([fn(rng.permutation(x)) for _ in range(n)])
    return obs, (1 + np.sum(null >= obs)) / (n + 1)


def windowed_cophs_multiband(X, fs, L, nperseg, bands):
    """Per-window cophenetic for several bands; cube computed once per window."""
    wlen = int(round(L * fs))
    fmax = max(BRAIN_BANDS[b][1] for b in bands) + 5.0
    out = {b: [] for b in bands}
    for s in range(0, X.shape[1] - wlen + 1, wlen):
        freqs, ff = segment_ffts(X[:, s:s + wlen], fs, nperseg, nperseg // 2, fmax)
        if ff.shape[1] < 4:
            continue
        cube = imcoh_abs_cube(ff)
        for b in bands:
            out[b].append(coph(band_abs_average(cube, freqs, *BRAIN_BANDS[b])))
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
        print(f"\n=== {pat} (L={args.L_sec}s nperseg={nperseg}) ===")
        try:
            def _load(ph):
                X = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
                return X if X.shape[0] <= X.shape[1] else X.T
            Xpre, Xtask, Xpost = _load("rest_pre"), _load("task_test"), _load("rest_post")
        except Exception as e:
            print(f"  [load fail] {e}"); continue

        fmax = max(BRAIN_BANDS[b][1] for b in GAMMA) + 5.0
        fr_p, sg_p = full_phase_signed(Xpre, fs, nperseg, fmax)
        fr_t, sg_t = full_phase_signed(Xtask, fs, nperseg, fmax)
        cwp = windowed_cophs_multiband(Xpre, fs, args.L_sec, nperseg, GAMMA)
        cwo = windowed_cophs_multiband(Xpost, fs, args.L_sec, nperseg, GAMMA)
        for b in GAMMA:
            lo, hi = BRAIN_BANDS[b]
            cpre_ref = coph(_band_abs(sg_p, fr_p, lo, hi))
            ctest_ref = coph(_band_abs(sg_t, fr_t, lo, hi))
            gp = np.array([spearmanr(c, ctest_ref)[0] - spearmanr(c, cpre_ref)[0] for c in cwp[b]])
            go = np.array([spearmanr(c, ctest_ref)[0] - spearmanr(c, cpre_ref)[0] for c in cwo[b]])
            pooled = np.concatenate([gp, go]); mu, sd = pooled.mean(), pooled.std()
            if sd > 0:
                gp = (gp - mu) / sd; go = (go - mu) / sd
            shift = float(go.mean() - gp.mean())
            q95 = np.quantile(gp, 0.95)
            excess_flash = float(np.mean(go > q95) - np.mean((gp + shift) > q95))
            ac_post, p_post = perm_p(go, ac1, args.n_perm, rng)
            ac_pre, _ = perm_p(gp, ac1, args.n_perm, rng)
            rows.append(dict(patient=pat, band=b, n_win_pre=len(gp), n_win_post=len(go),
                             shift=shift, spread_excess=float(go.std() - gp.std()),
                             AC1_post=ac_post, AC1_pre=ac_pre, p_time_post=p_post,
                             excess_flash=excess_flash))
            print(f"  {b:11s}: shift={shift:+.3f} spread_exc={go.std()-gp.std():+.3f}  "
                  f"AC1_post={ac_post:+.3f}(p={p_post:.3f})  excess_flash={excess_flash:+.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "gamma_perband_per_patient.csv", index=False)

    if df.patient.nunique() >= 5:
        def wtest(sub, col):
            d = sub[col].to_numpy(); d = d[np.isfinite(d)]
            try:
                p = wilcoxon(d, alternative="greater")[1]
            except Exception:
                p = np.nan
            return float(np.median(d)), int((d > 0).sum()), len(d), float(p)
        coh = []
        for b in GAMMA:
            sub = df[df.band == b].copy()
            sub["dAC1"] = sub.AC1_post - sub.AC1_pre
            for metric, col in [("shift_post_gt_pre", "shift"),
                                ("burstiness_AC1_post_gt_pre(states)", "dAC1"),
                                ("excess_flash(states)", "excess_flash"),
                                ("spread_excess(states)", "spread_excess")]:
                med, npos, n, p = wtest(sub, col)
                coh.append(dict(band=b, metric=metric, median=med, n_pos=f"{npos}/{n}",
                                wilcoxon_p=p,
                                n_sig=int((sub.p_time_post < 0.05).sum()) if "burst" in metric else ""))
        cohdf = pd.DataFrame(coh)
        cohdf.to_csv(OUT / "gamma_perband_cohort.csv", index=False)
        print("\n=== GAMMA per-band cohort verdict (low vs high) ===")
        print(cohdf.to_string(index=False))
    print(f"\n[audit_131] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
