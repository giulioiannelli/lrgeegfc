#!/usr/bin/env python3
"""Compute frequency-resolved signed ImCoh per (patient, phase, band).

Canonical Nolte 2004: Im(S_ij) / sqrt(S_ii * S_jj), range [-1, 1].

**Storage scheme:** per-band frequency-resolved signed ImCoh with shape
``(N, N, F_band)`` where F_band is the number of Welch frequency bins
inside the band. This preserves all information needed to derive any of
the standard band-averaged quantities at load time:

- ``imcoh``     → ``<signed>_f`` — sign-preserving average (Nolte 2004)
- ``imcoh_abs`` → ``<|ImCoh|>_f`` — connectivity strength (Ewald 2012)
- ``imcoh_sq`` → ``<|ImCoh|²>_f`` — squared imaginary coherence (Ewald 2012)

Order of operations matters: band-averaging signed values allows
positive/negative phase-lag bins to cancel (Jensen's inequality), so
``abs(<signed>_f) ≠ <|signed|>_f`` and ``(<signed>_f)² ≠ <signed²>_f``.
Hence we store the freq-resolved intermediate and band-average at load.

Storage footprint per patient: ~66 MB (dominated by high_gamma band).
Cohort total: ~330 MB.

Sequential over (patient, phase) with explicit memory freeing; a single
Welch call per (patient, phase) amortises CSD computation across all
bands.
"""
from __future__ import annotations

import argparse
import gc
import resource
import time
from pathlib import Path

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, PATIENTS_4PHASE, PHASE_LABELS, nperseg_for_fs
from lrg_eegfc.config.paths import IMCOH_CACHE, SEEG_DATAPATH
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch
from lrg_eegfc.utils.io.patient import (
    load_mat_pat_data,
    load_timeseries as _canonical_load_timeseries,
)

from lrg_eegfc.config.const import FS_OVERRIDES as FS_MAP  # canonical


def set_memory_limit_gb(gb: float) -> None:
    soft = int(gb * 1024 ** 3)
    resource.setrlimit(resource.RLIMIT_AS, (soft, soft))


def rss_mb() -> int:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024


def load_timeseries(patient: str, phase: str):
    try:
        X = _canonical_load_timeseries(patient, phase, SEEG_DATAPATH)
    except (FileNotFoundError, OSError):
        return None, None
    if X is None:
        return None, None
    X = np.asarray(X, dtype=np.float64)
    if X.shape[0] > X.shape[1]:
        X = X.T
    fs = FS_MAP.get(patient, 2048.0)
    return X, fs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--patients", nargs="+", default=PATIENTS_4PHASE)
    p.add_argument("--phases", nargs="+", default=list(PHASE_LABELS))
    p.add_argument("--bands", nargs="+", default=list(BRAIN_BANDS.keys()))
    p.add_argument("--mem-gb", type=float, default=12.0)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    set_memory_limit_gb(args.mem_gb)
    t_start = time.time()
    n_done = n_skip = n_fail = 0

    IMCOH_CACHE.mkdir(parents=True, exist_ok=True)

    for pat in args.patients:
        fs = FS_MAP.get(pat, 2048.0)
        nperseg = nperseg_for_fs(fs)
        pat_dir = IMCOH_CACHE / pat
        pat_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n=== {pat}  (fs={fs:.0f} Hz, nperseg={nperseg}) ===")

        for phase in args.phases:
            # Collect (band, out_path) for bands not yet cached
            to_compute = []
            for band in args.bands:
                out = pat_dir / (
                    f"{band}_{phase}_imcoh_freqresolved_"
                    f"nperseg-{nperseg}.npy"
                )
                if out.exists() and not args.overwrite:
                    n_skip += 1
                else:
                    to_compute.append((band, out))
            if not to_compute:
                continue

            t0 = time.time()
            X, fs_check = load_timeseries(pat, phase)
            if X is None:
                print(f"  {phase}: timeseries missing, skip")
                n_fail += len(to_compute)
                continue

            if args.verbose:
                print(f"  {phase}: X shape {X.shape}, computing signed Welch ImCoh...")

            # Single Welch call per (patient, phase) — CSD amortised across bands
            freqs, Coh = compute_msc_welch(
                X, fs_check, nperseg=nperseg, metric="imcoh",
            )  # shape (N, N, F), signed in [-1, 1]

            for band, out_path in to_compute:
                flo, fhi = BRAIN_BANDS[band]
                mask = (freqs >= flo) & (freqs <= fhi)
                if not mask.any():
                    print(f"  {band}: no freq bins in [{flo},{fhi}], skip")
                    n_fail += 1
                    continue
                sub = Coh[:, :, mask].astype(np.float32, copy=True)
                np.save(out_path, sub)
                if args.verbose:
                    print(
                        f"    {band:<10} -> {out_path.name} "
                        f"shape={sub.shape} "
                        f"range=[{sub.min():+.3f},{sub.max():+.3f}]"
                    )
                n_done += 1
                del sub

            del X, Coh
            gc.collect()
            if args.verbose:
                print(f"  {phase}: {time.time()-t0:.1f}s  RSS={rss_mb()} MB")

    dt = time.time() - t_start
    print(
        f"\nDone: {n_done} computed, {n_skip} skipped, {n_fail} failed "
        f"in {dt:.0f}s.\nCache root: {IMCOH_CACHE}"
    )


if __name__ == "__main__":
    main()
