#!/usr/bin/env python3
"""audit_138 — N4 v3 / P4: SUBSPACE (Grassmann) replay — does the leading diffusion
eigenmode transiently rotate toward the task mode?

Scope: .agents/guides/task-persistence-investigation/2026-06-22_replay-states-v3-scale-target-pair.md

P1 (directional) and P5 (raw propagator) both add clean transient-NEGATIVES. The last
genuinely-different propagator-based REPRESENTATION is the spectral SUBSPACE: instead of a
pairwise distance (cophenetic / raw propagator), use the k leading diffusion eigenmodes
(the low-λ Laplacian subspace the propagator e^{−τL} weights most). A replay state could be
a transient ROTATION of the window's dominant community subspace toward the task's, even if
the full pairwise geometry is stationary.

Per patient, band, subspace dim k (the scale analog):
  U_w^(k) = k lowest non-trivial Laplacian eigenvectors of windowed |ImCoh| (20 s).
  refs U_T^(k), U_{R-}^(k) full-phase. Projection alignment
    a(w, X) = ||U_w^(k)ᵀ U_X^(k)||_F² / k ∈ [0,1]   (mean cos² principal angle)
  s_w(k) = a(w, T) − a(w, R-)
  shift(k) = mean_post − mean_pre   (subspace sustained = dynamic-Grassmann N1; NB static
             Grassmann was NULL for the trace → informative either way)
  AC1/excess/spread vs TIME-SHUFFLE = transient subspace STATE
  PLACEBO: node-permuted task subspace (anatomy scrambled) — coarse-k drift control.

Verdict: a k* clears only if burstiness > time-shuffle AND target-specific (real>placebo).
Output: data/audit/replay_states/subspace_{per_patient,cohort}.csv
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.fc.coherence import segment_ffts, imcoh_abs_cube, band_abs_average
from lrg_eegfc.workflow.fc import load_fc_matrix

OUT = ROOT / "data" / "audit" / "replay_states"
OUT.mkdir(parents=True, exist_ok=True)
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
KGRID = [2, 4, 8, 16, 32]


def eig_of_W(W):
    W = 0.5 * (W + W.T); np.fill_diagonal(W, 0.0)
    L = np.diag(W.sum(1)) - W
    return np.linalg.eigh(L)            # ascending; col 0 ≈ trivial constant mode


def subspace(U, k):
    """k lowest NON-trivial eigenvectors (skip the λ≈0 constant mode)."""
    return U[:, 1:k + 1]


def align(Uw, Ut):
    """Projection metric ||Uwᵀ Ut||_F² / k ∈ [0,1] (mean squared principal-angle cos)."""
    k = Uw.shape[1]
    M = Uw.T @ Ut
    return float(np.sum(M * M) / k)


def ac1(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]; x = x - x.mean()
    d = np.dot(x, x)
    return float(np.dot(x[:-1], x[1:]) / d) if d > 0 and x.size > 3 else 0.0


def perm_p(x, fn, n, rng):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    obs = fn(x)
    null = np.array([fn(rng.permutation(x)) for _ in range(n)])
    return obs, (1 + np.sum(null >= obs)) / (n + 1)


def windowed_eigs(X, fs, L_sec, nperseg, band):
    wlen = int(round(L_sec * fs))
    fmax = BRAIN_BANDS[band][1] + 5.0
    out = []
    for s in range(0, X.shape[1] - wlen + 1, wlen):
        freqs, ff = segment_ffts(X[:, s:s + wlen], fs, nperseg, nperseg // 2, fmax)
        if ff.shape[1] < 4:
            continue
        W = band_abs_average(imcoh_abs_cube(ff), freqs, *BRAIN_BANDS[band])
        ev, U = eig_of_W(W)
        out.append(U)
    return out


def zscore_pair(a, b):
    pooled = np.concatenate([a, b])
    mu, sd = np.nanmean(pooled), np.nanstd(pooled)
    return (a - mu) / sd if sd > 0 else a * 0, (b - mu) / sd if sd > 0 else b * 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--bands", nargs="+", default=["alpha", "beta"])
    ap.add_argument("--L-sec", type=float, default=20.0)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260622)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    t0 = time.time()
    rows = []
    for pat in args.patients:
        fs = FS_OVERRIDES.get(pat, 2048.0)
        nperseg = nperseg_for_fs(fs)
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
            _, Ut = eig_of_W(np.asarray(Wt, float))
            _, Up = eig_of_W(np.asarray(Wp, float))
            perm = rng.permutation(n)
            _, Utp = eig_of_W(np.asarray(Wt, float)[np.ix_(perm, perm)])
            Uw_pre = windowed_eigs(Xpre, fs, args.L_sec, nperseg, band)
            Uw_post = windowed_eigs(Xpost, fs, args.L_sec, nperseg, band)
            for k in KGRID:
                if k + 1 >= n:
                    continue
                Tk, Pk, TPk = subspace(Ut, k), subspace(Up, k), subspace(Utp, k)

                def contrast(Ulist, Xk):
                    return np.array([align(subspace(U, k), Xk) - align(subspace(U, k), Pk)
                                     for U in Ulist])
                s_pre, s_post = contrast(Uw_pre, Tk), contrast(Uw_post, Tk)
                sp_post = contrast(Uw_post, TPk)
                if len(s_pre) < 5 or len(s_post) < 5:
                    continue
                zpre, zpost = zscore_pair(s_pre, s_post)
                shift = float(np.nanmean(zpost) - np.nanmean(zpre))
                ac_post, p_post = perm_p(zpost, ac1, args.n_perm, rng)
                ac_pre, _ = perm_p(zpre, ac1, args.n_perm, rng)
                spp = sp_post[np.isfinite(sp_post)]
                spp = (spp - spp.mean()) / spp.std() if spp.std() > 0 else spp * 0
                ac_plac = ac1(spp)
                gp = zpost[np.isfinite(zpost)]; q95 = np.nanquantile(zpre, 0.95)
                excess_flash = float(np.mean(gp > q95) - np.mean((zpre + shift) > q95))
                rows.append(dict(patient=pat, band=band, k=k, shift=shift,
                                 AC1_post=ac_post, AC1_pre=ac_pre, AC1_placebo=ac_plac,
                                 p_time_post=p_post, excess_flash=excess_flash,
                                 spread_excess=float(np.nanstd(zpost) - np.nanstd(zpre)),
                                 n_post=len(s_post)))
            print(f"  {pat} {band}: {len(Uw_post)} post-win × {len(KGRID)} k ({time.time()-t0:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "subspace_per_patient.csv", index=False)

    def wtest(d):
        d = np.asarray(d, float); d = d[np.isfinite(d)]
        if d.size < 5 or np.allclose(d, 0):
            return np.nan, 0, len(d)
        try:
            return float(wilcoxon(d, alternative="greater").pvalue), int((d > 0).sum()), len(d)
        except ValueError:
            return np.nan, int((d > 0).sum()), len(d)

    coh = []
    for band in args.bands:
        for k in KGRID:
            sub = df[(df.band == band) & (df.k == k)]
            if sub.empty:
                continue
            p_sh, n_sh, K = wtest(sub["shift"])
            dac = (sub.AC1_post - sub.AC1_pre).to_numpy()
            p_b, n_b, _ = wtest(dac)
            dspec = (sub.AC1_post - sub.AC1_placebo).to_numpy()
            p_sp, n_sp, _ = wtest(dspec)
            p_ef, n_ef, _ = wtest(sub.excess_flash)
            coh.append(dict(band=band, k=k, K=K,
                            shift_med=float(sub["shift"].median()), shift_npos=f"{n_sh}/{K}", shift_p=p_sh,
                            burst_dAC1_med=float(np.nanmedian(dac)), burst_npos=f"{n_b}/{K}", burst_p=p_b,
                            specific_p=p_sp, flash_npos=f"{n_ef}/{K}", flash_p=p_ef,
                            n_indiv_bursty=int((sub.p_time_post < 0.05).sum())))
    cohdf = pd.DataFrame(coh)
    cohdf.to_csv(OUT / "subspace_cohort.csv", index=False)

    print("\n=== P4 SUBSPACE (Grassmann) replay — cohort k-profile ===")
    print("  shift = subspace reinstatement LEVEL; burst = transient subspace STATE;")
    print("  specific = real>placebo; * = burst_p<0.05 AND specific_p<0.05")
    for band in args.bands:
        blk = cohdf[cohdf.band == band].sort_values("k")
        if blk.empty:
            continue
        print(f"\n  [{band}]")
        print("    k   shift(npos,p)        burst dAC1(npos,p)    specific  flash_p  indiv")
        for _, r in blk.iterrows():
            hit = "*" if (r.burst_p < 0.05 and r.specific_p < 0.05) else " "
            print(f" {hit}{r.k:4d}  {r.shift_med:+.3f} {r.shift_npos:>5} p={r.shift_p:.3f}  "
                  f"{r.burst_dAC1_med:+.3f} {r.burst_npos:>5} p={r.burst_p:.3f}  "
                  f"p={r.specific_p:.3f}   {r.flash_p:.3f}   {r.n_indiv_bursty:2d}/{r.K}")
    print(f"\n[audit_138] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
