#!/usr/bin/env python3
"""Cache rest_pre split halves under BOTH imcoh transforms (prep for the A2 sweep).

The incumbent halves cache (``data/cache/imcoh_halves_fc``) stores ``imcoh_abs``
only. The A2 stability sweep is run under both transforms (pre-registration
addendum R1-bis), so the ``imcoh_sq`` halves have to exist too. Both are written
here through the SAME library path
(:func:`lrg_eegfc.utils.fc.split_half.imcoh_split_half_adjacencies`) so the two
arms cannot diverge by an implementation difference.

Welch settings are identical to the incumbent: ``nperseg = nperseg_for_fs(fs) // 2``
on each contiguous half. The ``abs`` arm is written alongside and diffed against
the existing cache as a reproduction check (expected max|diff| ~ 1e-8, float32
round-trip of the incumbent files).

Output: IMCOH_HALVES_CACHE/{patient}/{band}_rest_pre_{A|B}_imcoh_{abs|sq}.npy
"""
from __future__ import annotations

import os
import time
from multiprocessing import Pool

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    DEFAULT_SAMPLE_RATE,
    FS_OVERRIDES,
    PATIENTS_4PHASE,
    nperseg_for_fs,
)
from lrg_eegfc.config.paths import CACHE_ROOT, IMCOH_HALVES_CACHE, SEEG_DATAPATH
from lrg_eegfc.utils.fc.split_half import imcoh_split_half_adjacencies
from lrg_eegfc.utils.io.patient import load_timeseries

OUT = IMCOH_HALVES_CACHE                       # canonical cache, config-driven
REF = CACHE_ROOT / "imcoh_halves_fc"           # legacy abs-only cache (reproduction check)
TRANSFORMS = ("abs", "sq")


def per_patient(pat: str) -> str:
    fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
    nper = max(256, nperseg_for_fs(fs) // 2)
    X = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH), float)
    if X.shape[0] > X.shape[1]:
        X = X.T
    d = OUT / pat
    d.mkdir(parents=True, exist_ok=True)
    worst = 0.0
    for tk in TRANSFORMS:
        halves = imcoh_split_half_adjacencies(X, fs, nper, BRAIN_BANDS, transform=tk)
        for (band, tag), A in halves.items():
            np.save(d / f"{band}_rest_pre_{tag}_imcoh_{tk}.npy", np.asarray(A, np.float64))
            if tk == "abs":
                ref = REF / pat / f"{band}_rest_pre_{tag}_imcoh_abs.npy"
                if ref.exists():
                    worst = max(worst, float(np.abs(np.asarray(np.load(ref), float) - A).max()))
    del X
    return f"  [{pat}] N={halves[('beta', 'A')].shape[0]} written; abs vs incumbent cache max|diff|={worst:.3e}"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ncpu = int(os.environ.get("SA_WORKERS", 5))
    pats = list(PATIENTS_4PHASE)
    print(f"[halves] {len(pats)} patients x {len(BRAIN_BANDS)} bands x 2 halves x "
          f"{len(TRANSFORMS)} transforms, {ncpu} workers -> {OUT}", flush=True)
    t0 = time.time()
    with Pool(ncpu) as pool:
        for i, msg in enumerate(pool.imap_unordered(per_patient, pats), 1):
            el = time.time() - t0
            print(f"[{i}/{len(pats)}] {el:.0f}s ETA {el / i * (len(pats) - i):.0f}s\n{msg}",
                  flush=True)
    print(f"\n[halves] done in {time.time() - t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
