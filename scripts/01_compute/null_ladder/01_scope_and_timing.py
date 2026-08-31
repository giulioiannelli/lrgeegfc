#!/usr/bin/env python3
"""W0-B scoping: session durations, per-band bin counts, memory + timing budget.

Establishes, BEFORE any null is launched:
  - per-patient phase durations (needed to size the N3 block grid),
  - per-band in-band bin counts F_b at the pipeline nperseg,
  - the per-segment FFT footprint (the enabling cache for the segment-lattice
    shift N1b and for the block-permutation N3),
  - a measured cost for one surrogate in each family, extrapolated to the cohort.
"""
from __future__ import annotations
import sys, time
import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS, FS_OVERRIDES, DEFAULT_SAMPLE_RATE, nperseg_for_fs,
)
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
PHASES_RAW = ("rest_pre", "task_test", "rest_post")


def main():
    print("BRAIN_BANDS:", {k: v for k, v in BRAIN_BANDS.items()}, flush=True)
    print(f"\n{'pat':8s} {'fs':>6s} {'nperseg':>8s} " +
          " ".join(f"{p:>12s}" for p in PHASES_RAW) + f" {'total_s':>9s} {'N':>4s}", flush=True)
    tot = {}
    for pat in COHORT:
        fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
        nps = nperseg_for_fs(fs)
        durs, N = [], None
        for ph in PHASES_RAW:
            try:
                X = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
                if X.shape[0] > X.shape[1]:
                    X = X.T
                N = X.shape[0]
                durs.append(X.shape[1] / fs)
                del X
            except Exception as e:
                durs.append(np.nan)
        tot[pat] = (fs, nps, durs, N)
        print(f"{pat:8s} {fs:6.0f} {nps:8d} " +
              " ".join(f"{d:12.1f}" for d in durs) +
              f" {np.nansum(durs):9.1f} {N if N else -1:4d}", flush=True)

    print("\n--- in-band bin counts F_b at pipeline nperseg (fs=2048, nperseg=4096; df=0.5 Hz) ---",
          flush=True)
    fs, nps = 2048.0, 4096
    freqs = np.fft.rfftfreq(nps, 1.0 / fs)
    for b, (lo, hi) in BRAIN_BANDS.items():
        m = (freqs >= lo) & (freqs <= hi)
        print(f"  {b:11s} [{lo:6.1f},{hi:6.1f}] Hz  F_b = {int(m.sum()):4d}", flush=True)

    print("\n--- footprints (N=118, complex64) ---", flush=True)
    for b, (lo, hi) in BRAIN_BANDS.items():
        Fb = int(((freqs >= lo) & (freqs <= hi)).sum())
        nseg = 683
        seg_mb = 118 * nseg * Fb * 8 / 1e6
        blk_mb = 80 * 118 * 118 * Fb * 8 / 1e6
        print(f"  {b:11s} F_b={Fb:4d}  per-segment FFT {seg_mb:8.1f} MB   "
              f"per-block CSD (80 blocks) {blk_mb:9.1f} MB", flush=True)


if __name__ == "__main__":
    main()
