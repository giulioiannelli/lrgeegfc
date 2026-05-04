#!/usr/bin/env python3
"""Audit Step 6 — restrict cross-patient H2a to the original N=5 core.

Core cohort: Pat_02, Pat_03, Pat_05, Pat_07, Pat_08.

Produces the same four 6 x K_max panels as Step 5, plus a small stats
file that compares the restricted picture to the N=9 case at the (band,
k) cells the user cares about (beta mid-k).

Outputs:
  data/audit/task_trace_old5_only.png
  data/audit/task_trace_old5_only.npz
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT

import sys
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_05_cross_patient import aggregate, plot  # noqa: E402

IN_NPZ = ROOT / "data" / "audit" / "task_trace_per_patient" / "task_trace_arrays.npz"
OUT_PNG = ROOT / "data" / "audit" / "task_trace_old5_only.png"
OUT_NPZ = ROOT / "data" / "audit" / "task_trace_old5_only.npz"

CORE_N5 = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
BAND_ORDER = list(BRAIN_BANDS.keys())


def main() -> None:
    d = np.load(IN_NPZ, allow_pickle=True)
    H2a_all = d["H2a"]
    pats = list(d["patient_ids"])
    bands = list(d["band_names"])
    k_values = np.asarray(d["k_values"], dtype=int)

    # Restrict to the original 5 core patients
    core_idx = [pats.index(p) for p in CORE_N5]
    H2a = H2a_all[core_idx]
    print(f"Restricted to {len(CORE_N5)} core patients: {CORE_N5}")
    print(f"Restricted H2a shape: {H2a.shape}")

    agg = aggregate(H2a)
    plot(agg, k_values, n_patients_total=len(CORE_N5), out_path=OUT_PNG,
          title_suffix="  (core N=5)")

    # Compare to the unrestricted aggregate at beta mid-k for the report
    ib_beta = bands.index("beta")
    beta_mean = agg["mean"][ib_beta]
    beta_npos = agg["n_pos"][ib_beta]
    total_contrib = agg["n_contrib"][ib_beta]

    # Report a small stats block: fraction of k with N5 unanimity in each band,
    # plus mean H2a across k and #patients>0.
    per_band_stats = []
    for ib, band in enumerate(bands):
        mat_contrib = agg["n_contrib"][ib]
        mat_mean = agg["mean"][ib]
        mat_npos = agg["n_pos"][ib]
        mat_nneg = agg["n_neg"][ib]
        # A cell is "unanimous +" if n_contrib == 5 AND n_pos == 5
        unanimous_pos = int(((mat_contrib == 5) & (mat_npos == 5)).sum())
        unanimous_neg = int(((mat_contrib == 5) & (mat_nneg == 5)).sum())
        valid = int((mat_contrib >= 2).sum())
        with np.errstate(invalid="ignore"):
            mean_mean = float(np.nanmean(mat_mean)) if valid else float("nan")
        per_band_stats.append({
            "band": band, "valid_cells": valid,
            "unanimous_pos": unanimous_pos,
            "unanimous_neg": unanimous_neg,
            "mean_of_mean": mean_mean,
        })

    print("\nPer-band core-N=5 summary:")
    for s in per_band_stats:
        print(f"  {s['band']:>10s}: mean={s['mean_of_mean']:+.4f}  "
              f"unan+={s['unanimous_pos']:3d}  "
              f"unan-={s['unanimous_neg']:3d}  "
              f"valid_k={s['valid_cells']}")

    np.savez_compressed(
        OUT_NPZ,
        mean=agg["mean"], median=agg["median"],
        n_pos=agg["n_pos"], n_neg=agg["n_neg"],
        n_contrib=agg["n_contrib"], sign_frac=agg["sign_frac"],
        band_names=np.array(bands), k_values=k_values,
        patient_ids=np.array(CORE_N5),
    )
    print(f"\nWrote {OUT_PNG.relative_to(ROOT)}")
    print(f"Wrote {OUT_NPZ.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
