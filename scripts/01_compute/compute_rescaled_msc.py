#!/usr/bin/env python3
"""Compute probe-rescaled MSC matrices for all patients, phases, bands.

Loads existing dense MSC caches and applies percentile rescaling to
same-probe edges (Mode C2), which eliminates the same-probe inflation
while preserving all channels and relative within-probe ordering.

Run: python scripts/py/compute_rescaled_msc.py [--patients Pat_02] [--verbose]
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, PATIENTS_4PHASE, PHASE_LABELS, nperseg_for_fs
from lrg_eegfc.config.paths import MSC_CACHE, RESCALED_CACHE, SEEG_DATAPATH
from lrg_eegfc.utils.fc.msc import rescale_same_probe_edges
from lrg_eegfc.workflow.msc import load_msc_matrix

RESCALED_MSC = RESCALED_CACHE / "msc_cache"
from lrg_eegfc.config.const import FS_OVERRIDES as FS_MAP  # canonical


def load_ch(pat):
    p = SEEG_DATAPATH / pat / "channel_labels.csv"
    with open(p) as f:
        first = f.readline().strip()
    skip = 1 if first.lower() == "label" else 0
    df = pd.read_csv(p, header=None, skiprows=skip)
    return [
        str(l).strip('"').split(",")[0].strip().replace(" ", "")
        for l in df.iloc[:, 0]
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--patients", nargs="+", default=PATIENTS_4PHASE)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    total_t = time.time()
    n_done = 0

    for pat in args.patients:
        fs = FS_MAP.get(pat, 2048.0)
        nperseg = nperseg_for_fs(fs)
        ch = load_ch(pat)

        print(f"\n{pat} ({len(ch)} channels)")

        for phase in PHASE_LABELS:
            for band in BRAIN_BANDS:
                out_dir = RESCALED_MSC / pat
                out_dir.mkdir(parents=True, exist_ok=True)
                out_path = out_dir / f"{band}_{phase}_msc_sparsify-none_nperseg-{nperseg}.npy"

                if out_path.exists():
                    continue

                A = load_msc_matrix(
                    pat, phase, band,
                    cache_root=MSC_CACHE,
                    sparsify="none",
                    n_surrogates=0,
                    nperseg=nperseg,
                )
                if A is None:
                    if args.verbose:
                        print(f"  {phase}/{band}: no source MSC")
                    continue

                A_rescaled = rescale_same_probe_edges(A, ch)
                np.save(out_path, A_rescaled)
                n_done += 1

                if args.verbose:
                    print(f"  {phase}/{band}: rescaled → {out_path.name}")

    print(f"\nDone: {n_done} matrices in {time.time()-total_t:.0f}s")
    print(f"Output: {RESCALED_MSC}")


if __name__ == "__main__":
    main()
