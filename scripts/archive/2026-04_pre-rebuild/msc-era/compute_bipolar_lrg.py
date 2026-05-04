#!/usr/bin/env python3
"""Compute LRG analysis on alternative MSC matrices (bipolar or rescaled).

Reads from specified MSC cache, writes to specified LRG cache.

Run:
  python scripts/py/compute_bipolar_lrg.py --mode rescaled [--patients Pat_02] [-v]
  python scripts/py/compute_bipolar_lrg.py --mode bipolar [--patients Pat_02] [-v]
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, PATIENTS_4PHASE, PHASE_LABELS, nperseg_for_fs
from lrg_eegfc.config.paths import BIPOLAR_CACHE
from lrg_eegfc.workflow.lrg import compute_lrg_analysis

BIPOLAR_MSC = BIPOLAR_CACHE / "msc_cache"
BIPOLAR_LRG = BIPOLAR_CACHE / "lrg_cache"
from lrg_eegfc.config.const import FS_OVERRIDES as FS_MAP  # canonical


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--patients", nargs="+", default=PATIENTS_4PHASE)
    parser.add_argument("--phases", nargs="+", default=list(PHASE_LABELS))
    parser.add_argument("--bands", nargs="+", default=list(BRAIN_BANDS.keys()))
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    total_t = time.time()
    n_done = 0
    n_skip = 0
    n_fail = 0

    for pat in args.patients:
        fs = FS_MAP.get(pat, 2048.0)
        nperseg = nperseg_for_fs(fs)

        print(f"\n{'='*60}")
        print(f"{pat}")
        print(f"{'='*60}")

        for phase in args.phases:
            for band in args.bands:
                # Load bipolar MSC
                msc_path = (
                    BIPOLAR_MSC / pat
                    / f"{band}_{phase}_msc_sparsify-none_nperseg-{nperseg}.npy"
                )
                if not msc_path.exists():
                    if args.verbose:
                        print(f"  {phase}/{band}: no MSC, skipping")
                    n_skip += 1
                    continue

                # Check if LRG already cached
                out_dir = BIPOLAR_LRG / pat
                out_path = out_dir / f"{band}_{phase}_lrg_msc.npz"
                if out_path.exists():
                    if args.verbose:
                        print(f"  {phase}/{band}: LRG cached, skipping")
                    n_skip += 1
                    continue

                t0 = time.time()
                A = np.load(msc_path)
                np.fill_diagonal(A, 0)

                try:
                    result = compute_lrg_analysis(
                        A, patient=pat, phase=phase, band=band,
                        fc_method="msc",
                    )
                except Exception as e:
                    print(f"  {phase}/{band}: FAILED — {e}")
                    n_fail += 1
                    continue

                # Save manually (compute_lrg_analysis caches to default dir)
                out_dir.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(
                    out_path,
                    ultrametric_matrix=result.ultrametric_matrix,
                    linkage_matrix=result.linkage_matrix,
                    entropy_tau=result.entropy_tau,
                    entropy_1_minus_S=result.entropy_1_minus_S,
                    entropy_C=result.entropy_C,
                    optimal_threshold=result.optimal_threshold,
                    patient=pat,
                    phase=phase,
                    band=band,
                    fc_method="msc",
                    n_nodes=result.n_nodes,
                )

                dt = time.time() - t0
                if args.verbose:
                    print(f"  {phase}/{band}: n={result.n_nodes}, "
                          f"threshold={result.optimal_threshold:.4f}, "
                          f"{dt:.1f}s")
                n_done += 1

    print(f"\nDone: {n_done} computed, {n_skip} skipped, {n_fail} failed "
          f"in {time.time()-total_t:.0f}s")
    print(f"Output: {BIPOLAR_LRG}")


if __name__ == "__main__":
    main()
