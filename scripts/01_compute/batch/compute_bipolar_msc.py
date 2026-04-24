#!/usr/bin/env python3
"""Compute bipolar-referenced MSC matrices for all patients and phases.

Applies bipolar re-referencing BEFORE computing MSC, then saves all
band-specific adjacency matrices to data/bipolar/msc_cache/.

Run: python scripts/py/compute_bipolar_msc.py [--patients Pat_02 Pat_05] [--verbose]
"""
from __future__ import annotations

import argparse
import gc
import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    PATIENTS_4PHASE,
    PHASE_LABELS,
    nperseg_for_fs,
)
from lrg_eegfc.config.paths import BIPOLAR_CACHE, SEEG_DATAPATH
from lrg_eegfc.utils.fc.msc import coherence_fc_pipeline
from lrg_eegfc.utils.io.patient import (
    bipolar_rereference,
    load_timeseries,
)

# ── Config ─────────────────────────────────────────────────────

DATA_ROOT = SEEG_DATAPATH
BIPOLAR_MSC_CACHE = BIPOLAR_CACHE / "msc_cache"
BIPOLAR_CHANNELS = BIPOLAR_CACHE / "channels"

# Sampling rates per patient
from lrg_eegfc.config.const import FS_OVERRIDES as FS_MAP  # canonical
DEFAULT_FS = 2048.0


def get_fs(patient: str) -> float:
    return FS_MAP.get(patient, DEFAULT_FS)


def load_channel_labels_robust(patient: str) -> list[str]:
    path = DATA_ROOT / patient / "channel_labels.csv"
    with open(path) as f:
        first = f.readline().strip()
    skip = 1 if first.lower() == "label" else 0
    df = pd.read_csv(path, header=None, skiprows=skip)
    return [
        str(l).strip('"').split(",")[0].strip().replace(" ", "")
        for l in df.iloc[:, 0]
    ]


def compute_patient_phase(
    patient: str,
    phase: str,
    verbose: bool = False,
) -> None:
    """Compute bipolar MSC for one (patient, phase) and save all bands."""
    fs = get_fs(patient)
    nperseg = nperseg_for_fs(fs)

    # Check if all bands already cached
    all_cached = True
    for band in BRAIN_BANDS:
        out = BIPOLAR_MSC_CACHE / patient / f"{band}_{phase}_msc_sparsify-none_nperseg-{nperseg}.npy"
        if not out.exists():
            all_cached = False
            break
    if all_cached:
        if verbose:
            print(f"  {patient}/{phase}: all bands cached, skipping")
        return

    # Load raw data
    t0 = time.time()
    labels = load_channel_labels_robust(patient)
    ts = load_timeseries(patient, phase, DATA_ROOT)
    if verbose:
        print(f"  {patient}/{phase}: loaded {ts.shape} in {time.time()-t0:.1f}s")

    # Bipolar re-reference
    t1 = time.time()
    bip_ts, bip_labels, _ = bipolar_rereference(ts, labels)
    if verbose:
        print(f"  {patient}/{phase}: bipolar {ts.shape[0]}→{bip_ts.shape[0]} in {time.time()-t1:.1f}s")

    # Save bipolar channel info (once per patient)
    ch_dir = BIPOLAR_CHANNELS / patient
    ch_dir.mkdir(parents=True, exist_ok=True)
    labels_path = ch_dir / "bipolar_labels.json"
    if not labels_path.exists():
        with open(labels_path, "w") as f:
            json.dump(bip_labels, f, indent=2)

    # Compute MSC
    t2 = time.time()
    band_matrices = coherence_fc_pipeline(
        bip_ts,
        fs,
        nperseg=nperseg,
        noverlap=nperseg // 2,
        bands=BRAIN_BANDS,
        sparsify="none",
        n_surrogates=0,
    )
    if verbose:
        print(f"  {patient}/{phase}: MSC computed in {time.time()-t2:.1f}s")

    # Save per-band matrices
    out_dir = BIPOLAR_MSC_CACHE / patient
    out_dir.mkdir(parents=True, exist_ok=True)
    for band, W in band_matrices.items():
        out = out_dir / f"{band}_{phase}_msc_sparsify-none_nperseg-{nperseg}.npy"
        np.save(out, W)

    if verbose:
        print(f"  {patient}/{phase}: saved {len(band_matrices)} bands, "
              f"total {time.time()-t0:.1f}s")

    # Free memory
    del ts, bip_ts, band_matrices
    gc.collect()


def main():
    parser = argparse.ArgumentParser(description="Compute bipolar MSC")
    parser.add_argument("--patients", nargs="+", default=PATIENTS_4PHASE)
    parser.add_argument("--phases", nargs="+", default=list(PHASE_LABELS))
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING)

    total_t = time.time()
    n_done = 0

    for patient in args.patients:
        print(f"\n{'='*60}")
        print(f"{patient} (fs={get_fs(patient)} Hz)")
        print(f"{'='*60}")

        for phase in args.phases:
            try:
                compute_patient_phase(patient, phase, verbose=args.verbose)
                n_done += 1
            except Exception as e:
                print(f"  ERROR {patient}/{phase}: {e}")

    print(f"\nDone: {n_done} patient×phase combos in {time.time()-total_t:.0f}s")
    print(f"Output: {BIPOLAR_MSC_CACHE}")


if __name__ == "__main__":
    main()
