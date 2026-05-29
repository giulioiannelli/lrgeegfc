#!/usr/bin/env python3
"""WP0 Task 2 — Compute all metrics on the full dataset.

Computes 11 LRG-based metrics + 4 raw FC baseline metrics for every
(patient, band, phase-pair) combination. Outputs a single CSV.

Produces:
  data/wp0_metric_exploration/task2_full_results.csv

Run: python scripts/wp0/task2_compute.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "wp0"))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd

from _common import (
    ALL_PAIRS,
    BANDS,
    FC_METHOD,
    OUT_ROOT,
    PATIENTS,
    PHASES,
    classify_pair,
    get_fc_metric_specs,
    get_lrg_metric_specs,
    load_all_fc,
    load_all_lrg,
)

CSV_OUT = OUT_ROOT / "task2_full_results.csv"


def main():
    CSV_OUT.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # ── Load data ─────────────────────────────────────────────────────────
    print("Loading LRG results...")
    lrg_data = load_all_lrg()
    print("Loading raw FC matrices...")
    fc_data = load_all_fc()

    # ── Build metric specs ────────────────────────────────────────────────
    lrg_specs = get_lrg_metric_specs()
    fc_specs = get_fc_metric_specs()

    # ── Compute ───────────────────────────────────────────────────────────
    rows = []
    n_computed = 0
    n_skipped = 0

    for pat in PATIENTS:
        for band in BANDS:
            for pa, pb in ALL_PAIRS:
                pair_type = classify_pair(pa, pb)

                # --- LRG metrics ---
                key_a = (pat, pa, band)
                key_b = (pat, pb, band)
                lrg_a = lrg_data.get(key_a)
                lrg_b = lrg_data.get(key_b)

                if lrg_a is not None and lrg_b is not None:
                    for mname, mspec in lrg_specs.items():
                        try:
                            val = mspec["fn"](lrg_a, lrg_b)
                        except Exception as e:
                            val = np.nan
                            print(f"  WARN: {mname} failed for {pat} {band} {pa}-{pb}: {e}")
                        rows.append({
                            "patient": pat,
                            "band": band,
                            "phase_a": pa,
                            "phase_b": pb,
                            "pair_type": pair_type,
                            "level": "lrg",
                            "metric_name": mname,
                            "value": val,
                        })
                        n_computed += 1
                else:
                    n_skipped += len(lrg_specs)

                # --- Raw FC metrics ---
                fc_a = fc_data.get(key_a)
                fc_b = fc_data.get(key_b)

                if fc_a is not None and fc_b is not None:
                    for mname, mspec in fc_specs.items():
                        try:
                            val = mspec["fn"](fc_a, fc_b)
                        except Exception as e:
                            val = np.nan
                            print(f"  WARN: {mname} failed for {pat} {band} {pa}-{pb}: {e}")
                        rows.append({
                            "patient": pat,
                            "band": band,
                            "phase_a": pa,
                            "phase_b": pb,
                            "pair_type": pair_type,
                            "level": "raw_fc",
                            "metric_name": mname,
                            "value": val,
                        })
                        n_computed += 1
                else:
                    n_skipped += len(fc_specs)

    df = pd.DataFrame(rows)

    # ── Sanity checks ─────────────────────────────────────────────────────
    print(f"\n── Sanity checks ──")
    print(f"  Rows computed: {n_computed}")
    print(f"  Rows skipped (missing data): {n_skipped}")
    print(f"  DataFrame shape: {df.shape}")
    print(f"  NaN values: {df['value'].isna().sum()}")
    print(f"  Unique patients: {sorted(df['patient'].unique())}")
    print(f"  Unique bands: {sorted(df['band'].unique())}")
    print(f"  Unique metrics: {sorted(df['metric_name'].unique())}")
    print(f"  Metrics per level:")
    for level, grp in df.groupby("level"):
        print(f"    {level}: {sorted(grp['metric_name'].unique())}")

    # Check value ranges for bounded metrics
    bounded = {
        "scaled_distance": (0, 1),
        "tree_robinson_foulds": (0, 1),
        "tree_cophenetic_corr": (-1, 1),
        "tree_baker_gamma": (-1, 1),
        "tree_fowlkes_mallows": (0, 1),
        "ari": (-1, 1),
        "cluster_swap": (0, 1),
        "fc_scaled_frobenius": (0, 1),
        "rank_distance": (0, 2),
        "fc_rank_distance": (0, 2),
    }
    for mname, (lo, hi) in bounded.items():
        sub = df[df["metric_name"] == mname]["value"].dropna()
        if len(sub) == 0:
            continue
        oob = ((sub < lo - 1e-9) | (sub > hi + 1e-9)).sum()
        if oob > 0:
            print(f"  WARNING: {mname} has {oob} out-of-range values!")

    # ── Save ──────────────────────────────────────────────────────────────
    df.to_csv(CSV_OUT, index=False)
    elapsed = time.time() - t0
    print(f"\nTask 2 complete: {CSV_OUT}")
    print(f"  {len(df)} rows, {elapsed:.1f}s")


if __name__ == "__main__":
    main()
