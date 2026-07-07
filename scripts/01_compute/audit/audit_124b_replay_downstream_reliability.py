#!/usr/bin/env python3
"""Audit 124b — N4 PRECONDITION (decisive): downstream reliability of the
per-window cophenetic config and task-likeness `s`.

audit_124 showed the RAW per-window |ImCoh| edge ranking is noisy for the trace
bands (split-half Spearman ~0.1-0.2). But |ImCoh| has many near-zero edges
(volume-conduction immune), so ranking the full edge vector is pessimistic — what
N4 actually consumes is the LRG COPHENETIC vector (driven by the strong-edge
backbone) and, on top of it, the TASK-LIKENESS scalar

    s[w] = Spearman( c[w] - cbar_pre ,  cbar_T - cbar_pre )      (T = task_test)

If `s[w]` carries reproducible between-window signal, replay states are
detectable (the burstiness test has power). If `s[w]` is pure estimator noise,
the S time-course is white noise and N4 cannot be established (a null result would
be unpowered, not informative).

Method (per patient, band, window length L, phase=rest_post unless --phases):
  For each non-overlapping window, split its Welch segments even/odd, build
  |ImCoh|_band on each half, run the LRG cophenetic pipeline -> c_even, c_odd, and
  also on ALL segments -> c_full. Then:
    - cophenetic split-half reliability  = Spearman(c_even, c_odd)            [per window]
    - task-likeness s_even, s_odd, s_full                                     [per window]
    - r_half(s)  = Pearson(s_even[:], s_odd[:]) across windows
    - Spearman-Brown full-window reliability  r_full = 2 r_half / (1 + r_half)
    - variance decomposition (independent of SB):
         noise_var_half = Var(s_even - s_odd) / 2
         noise_var_full ~ noise_var_half / 2     (full window ~2x segments)
         signal_var     = Var(s_full) - noise_var_full
         reliability_vd = max(0, signal_var) / Var(s_full)

DECISION NUMBER: r_full(s) and reliability_vd. >~0.5 => detectable; ~0 => either
no transient states OR noise-swamped (report Var(s_full) to disambiguate: large
spread + low reliability = noise; small spread = little to detect).

NO acceptance gate pre-registered (feedback_no_pre_registered_acceptance) — we
report the numbers and decide window length / band set with the PI.

Outputs
-------
    data/audit/replay_windowability/downstream_reliability.csv
    data/audit/replay_windowability/figures/downstream_reliability.pdf
"""
from __future__ import annotations

import argparse
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.fc.coherence import segment_ffts, imcoh_abs_cube, band_abs_average
from lrg_eegfc.utils.surrogate.matched_strength import cophenetic_condensed_from_adjacency
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.visuals.styles import use_lrg_style

use_lrg_style()

OUT = ROOT / "data" / "audit" / "replay_windowability"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)


def coph(W):
    return cophenetic_condensed_from_adjacency(np.clip(np.asarray(W, float), 0.0, 1.0))


def task_likeness(c, cbar_pre, cbar_T):
    """Spearman over pairs of (c - cbar_pre) against (cbar_T - cbar_pre)."""
    r, _ = spearmanr(c - cbar_pre, cbar_T - cbar_pre)
    return float(r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=["Pat_05", "Pat_02", "Pat_13"])
    ap.add_argument("--bands", nargs="+", default=["theta", "alpha", "beta", "low_gamma"])
    ap.add_argument("--phases", nargs="+", default=["rest_post"])
    ap.add_argument("--lengths", nargs="+", type=float, default=[20.0, 30.0, 45.0])
    ap.add_argument("--max-windows", type=int, default=40)
    ap.add_argument("--target-phase", default="task_test")
    args = ap.parse_args()

    t0 = time.time()
    rows = []
    for pat in args.patients:
        fs = FS_OVERRIDES.get(pat, 2048.0)
        nperseg = nperseg_for_fs(fs)
        print(f"\n=== {pat} (fs={fs:.0f}) ===")
        # references (full-phase cached imcoh_abs -> cophenetic)
        refs = {}
        for band in args.bands:
            Wpre = load_fc_matrix(pat, "rest_pre", band, fc_method="imcoh_abs")
            WT = load_fc_matrix(pat, args.target_phase, band, fc_method="imcoh_abs")
            if Wpre is None or WT is None:
                print(f"  [ref miss] {band}")
                continue
            refs[band] = (coph(Wpre), coph(WT))

        for ph in args.phases:
            try:
                X = load_timeseries(pat, ph, SEEG_DATAPATH)
            except Exception as e:
                print(f"  [load fail] {ph}: {e}")
                continue
            X = np.asarray(X, float)
            if X.shape[0] > X.shape[1]:
                X = X.T
            N, T = X.shape

            for L in args.lengths:
                wlen = int(round(L * fs))
                if wlen > T:
                    continue
                starts = list(range(0, T - wlen + 1, wlen))   # NON-overlapping
                if len(starts) > args.max_windows:
                    starts = list(np.linspace(0, T - wlen, args.max_windows, dtype=int))
                fmax = max(BRAIN_BANDS[b][1] for b in args.bands) + 5.0

                acc = {b: dict(rcoph=[], se=[], so=[], sf=[]) for b in args.bands}
                for s in starts:
                    win = X[:, s:s + wlen]
                    freqs, ff = segment_ffts(win, fs, nperseg, fmax_keep=fmax)
                    nseg = ff.shape[1]
                    if nseg < 4:
                        continue
                    ev = np.arange(0, nseg, 2)
                    od = np.arange(1, nseg, 2)
                    cube_e = imcoh_abs_cube(ff, ev)
                    cube_o = imcoh_abs_cube(ff, od)
                    cube_f = imcoh_abs_cube(ff)
                    for band in args.bands:
                        if band not in refs:
                            continue
                        lo, hi = BRAIN_BANDS[band]
                        We = band_abs_average(cube_e, freqs, lo, hi)
                        Wo = band_abs_average(cube_o, freqs, lo, hi)
                        Wf = band_abs_average(cube_f, freqs, lo, hi)
                        if We is None or Wo is None or Wf is None:
                            continue
                        ce, co, cf = coph(We), coph(Wo), coph(Wf)
                        rc, _ = spearmanr(ce, co)
                        cbar_pre, cbar_T = refs[band]
                        acc[band]["rcoph"].append(rc)
                        acc[band]["se"].append(task_likeness(ce, cbar_pre, cbar_T))
                        acc[band]["so"].append(task_likeness(co, cbar_pre, cbar_T))
                        acc[band]["sf"].append(task_likeness(cf, cbar_pre, cbar_T))

                for band in args.bands:
                    d = acc[band]
                    if len(d["sf"]) < 5:
                        continue
                    se = np.array(d["se"]); so = np.array(d["so"]); sf = np.array(d["sf"])
                    rcoph = float(np.mean(d["rcoph"]))
                    # reliability of per-window task-likeness
                    if np.std(se) > 0 and np.std(so) > 0:
                        r_half = float(np.corrcoef(se, so)[0, 1])
                    else:
                        r_half = np.nan
                    r_full_sb = (2 * r_half / (1 + r_half)) if np.isfinite(r_half) and r_half > -1 else np.nan
                    var_sf = float(np.var(sf, ddof=1))
                    noise_half = float(np.var(se - so, ddof=1) / 2.0)
                    noise_full = noise_half / 2.0
                    signal_var = var_sf - noise_full
                    rel_vd = max(0.0, signal_var) / var_sf if var_sf > 0 else np.nan
                    rows.append(dict(
                        patient=pat, phase=ph, band=band, L_sec=L, n_win=len(sf),
                        coph_splithalf=rcoph,
                        s_full_mean=float(np.mean(sf)), s_full_sd=float(np.std(sf, ddof=1)),
                        r_half_s=r_half, r_full_s_SB=r_full_sb, reliability_vd=rel_vd,
                    ))
                    print(f"  [{ph} {band:9s} L={L:4.0f}s n={len(sf):2d}] "
                          f"coph_r={rcoph:.2f}  s_mean={np.mean(sf):+.3f} s_sd={np.std(sf):.3f}  "
                          f"r_half(s)={r_half:+.2f} r_full(s)={r_full_sb:+.2f} rel_vd={rel_vd:.2f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "downstream_reliability.csv", index=False)

    print("\n=== SUMMARY: full-window task-likeness reliability r_full(s) [Spearman-Brown] ===")
    if not df.empty:
        piv = df.pivot_table(index="band", columns="L_sec", values="r_full_s_SB", aggfunc="mean")
        print(piv.reindex([b for b in args.bands if b in piv.index]).to_string(float_format=lambda v: f"{v:.2f}"))
        print("\n=== variance-decomposition reliability ===")
        piv2 = df.pivot_table(index="band", columns="L_sec", values="reliability_vd", aggfunc="mean")
        print(piv2.reindex([b for b in args.bands if b in piv2.index]).to_string(float_format=lambda v: f"{v:.2f}"))

        # figure
        fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2))
        for band in args.bands:
            sub = df[df.band == band]
            if sub.empty:
                continue
            g1 = sub.groupby("L_sec")["r_full_s_SB"].mean()
            g2 = sub.groupby("L_sec")["coph_splithalf"].mean()
            axes[0].plot(g1.index, g1.values, "-o", ms=4, label=BRAIN_BAND_TEX_DICT.get(band, band))
            axes[1].plot(g2.index, g2.values, "-o", ms=4, label=BRAIN_BAND_TEX_DICT.get(band, band))
        for ax, ttl in zip(axes, [r"task-likeness $r_{\mathrm{full}}(s)$ (Spearman-Brown)",
                                   r"cophenetic split-half $r$"]):
            ax.set_xlabel("window length L (s)")
            ax.set_ylabel(ttl)
            ax.axhline(0.5, color="0.6", ls="--", lw=0.8)
            ax.axhline(0.0, color="0.8", ls="-", lw=0.6)
            ax.legend(frameon=False, fontsize=8)
            ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(FIG / "downstream_reliability.pdf", transparent=True)
        plt.close(fig)

    print(f"\n[audit_124b] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
