#!/usr/bin/env python3
"""Audit 128 — N4 v2 / Direction B0: SHORT-WINDOW |ImCoh| precondition (the gate).

20 s windows are long; fast replay would be brief and in fast bands. But |ImCoh|
needs several Welch segments each spanning enough cycles, so short windows force a
small nperseg -> coarse frequency resolution -> the narrow alpha/beta bands become
unrepresentable. This script measures, empirically, the timescale FLOOR per band:
for L in {0.2..20} s it picks nperseg to keep ~12-16 segments, then reports
  - df = fs/nperseg (frequency resolution) and n_bins in each band (feasible if >=2)
  - within-window split-half reliability r of |ImCoh| (even/odd segments)
so we know which (band, window-length) cells can carry any replay test at all.

EXPECTED: alpha/beta die below ~10 s; only gamma survives at <=2 s. Confirm.

Outputs: data/audit/replay_windowability/short_window.csv (+ console table)
"""
from __future__ import annotations

import argparse
import time

import numpy as np
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.fc.coherence import segment_ffts, imcoh_abs_cube, band_abs_average

OUT = ROOT / "data" / "audit" / "replay_windowability"
OUT.mkdir(parents=True, exist_ok=True)
BANDS = list(BRAIN_BANDS.keys())
L_GRID = [0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]


def pick_nperseg(L, fs, target_seg=12, cap=4096):
    nps = int(2 * L * fs / (target_seg + 1))           # 50% overlap -> n_seg ~ target
    nps = 2 ** int(np.floor(np.log2(max(16, nps))))    # power of two
    return min(nps, cap)


def triu(W):
    iu = np.triu_indices(W.shape[0], 1)
    return W[iu]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=["Pat_05", "Pat_02"])
    ap.add_argument("--phase", default="rest_post")
    ap.add_argument("--n-windows", type=int, default=15)
    args = ap.parse_args()

    import pandas as pd
    rows = []
    t0 = time.time()
    for pat in args.patients:
        fs = FS_OVERRIDES.get(pat, 2048.0)
        print(f"\n=== {pat} (fs={fs:.0f}) ===")
        X = load_timeseries(pat, args.phase, SEEG_DATAPATH)
        X = np.asarray(X, float)
        if X.shape[0] > X.shape[1]:
            X = X.T
        N, T = X.shape
        fmax = max(BRAIN_BANDS[b][1] for b in BANDS) + 5.0
        for L in L_GRID:
            wlen = int(round(L * fs))
            if wlen > T or wlen < 16:
                continue
            nperseg = pick_nperseg(L, fs)
            if nperseg > wlen:
                nperseg = 2 ** int(np.floor(np.log2(wlen)))
            df = fs / nperseg
            starts = np.linspace(0, T - wlen, args.n_windows, dtype=int)
            rs = {b: [] for b in BANDS}
            nbins = {b: int(((np.fft.rfftfreq(nperseg, 1/fs) >= BRAIN_BANDS[b][0]) &
                             (np.fft.rfftfreq(nperseg, 1/fs) <= BRAIN_BANDS[b][1])).sum())
                     for b in BANDS}
            nseg_used = None
            for s in starts:
                win = X[:, s:s + wlen]
                freqs, ff = segment_ffts(win, fs, nperseg, nperseg // 2, fmax)
                nseg = ff.shape[1]
                nseg_used = nseg
                if nseg < 4:
                    continue
                ev = np.arange(0, nseg, 2); od = np.arange(1, nseg, 2)
                ce = imcoh_abs_cube(ff, ev); co = imcoh_abs_cube(ff, od)
                if ce is None or co is None:
                    continue
                for b in BANDS:
                    if nbins[b] < 1:
                        continue
                    We = band_abs_average(ce, freqs, *BRAIN_BANDS[b])
                    Wo = band_abs_average(co, freqs, *BRAIN_BANDS[b])
                    if We is None or Wo is None:
                        continue
                    r, _ = spearmanr(triu(We), triu(Wo))
                    if np.isfinite(r):
                        rs[b].append(r)
            for b in BANDS:
                rows.append(dict(patient=pat, L_sec=L, nperseg=nperseg, df_Hz=df,
                                 band=b, n_bins=nbins[b], n_seg=nseg_used,
                                 feasible=nbins[b] >= 2,
                                 r_mean=float(np.mean(rs[b])) if rs[b] else np.nan))
            cells = "  ".join(
                f"{b[:5]}={'--' if nbins[b] < 1 else (f'{np.mean(rs[b]):.2f}' if rs[b] else 'x')}"
                f"({nbins[b]})" for b in BANDS)
            print(f"  L={L:5.1f}s nps={nperseg:4d} df={df:5.1f}Hz nseg={nseg_used}: {cells}")

    df_out = pd.DataFrame(rows)
    df_out.to_csv(OUT / "short_window.csv", index=False)
    print("\n=== reliability r by (band, L), feasible cells only (n_bins>=2), patient-avg ===")
    feas = df_out[df_out.feasible]
    piv = feas.pivot_table(index="band", columns="L_sec", values="r_mean", aggfunc="mean")
    piv = piv.reindex([b for b in BANDS if b in piv.index])
    print(piv.to_string(float_format=lambda v: f"{v:.2f}"))
    print(f"\n[audit_128] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
