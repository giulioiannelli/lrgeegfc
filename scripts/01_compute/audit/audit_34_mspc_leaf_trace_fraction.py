#!/usr/bin/env python3
"""Audit 34 — MSPC f_trace per-leaf re-mining.

Re-mines ``data/audit/per_patient_hierarchy_mspc/multiscale_assignment.csv``
(already on disk; ~407k rows). Per-leaf trace fraction across the integer-k
grid:

    f_trace(p, b, l) = |{k : dominant_at_k(p, b, l, k) = "trace"}| / |k_grid|

Per-(patient, band): the distribution of f_trace across leaves; the count of
TRACE-leaves at threshold f_trace >= 0.5; the cohort distribution per band.
This complements the K-averaged leaf classification already at
``leaf_assignment.csv`` by exposing the scale-explicit leaf-level trace
fraction, which is what the localization layer (measure 09) consumes.

Outputs
-------
``data/audit/mspc_leaf_trace_fraction/leaf_trace_fraction.csv``
``data/audit/mspc_leaf_trace_fraction/cohort_summary.csv``
``data/outputs/figures/section_5_lrg_trace/mspc_f_trace/cohort_distribution.pdf``
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT

MSPC_CELLS = ROOT / "data" / "audit" / "per_patient_hierarchy_mspc" / "multiscale_assignment.csv"
OUT_DIR = ROOT / "data" / "audit" / "mspc_leaf_trace_fraction"
FIG_DIR = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "mspc_f_trace"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(MSPC_CELLS)
    df["is_trace"] = (df["dominant_at_k"] == "trace").astype(int)
    df["is_persist"] = (df["dominant_at_k"] == "persist").astype(int)
    df["is_reset"] = (df["dominant_at_k"] == "reset").astype(int)
    df["is_rearrange"] = (df["dominant_at_k"] == "rearrange").astype(int)

    leaf = (
        df.groupby(["patient", "band", "leaf_id"])
        .agg(
            n_k=("k", "count"),
            f_trace=("is_trace", "mean"),
            f_persist=("is_persist", "mean"),
            f_reset=("is_reset", "mean"),
            f_rearrange=("is_rearrange", "mean"),
        )
        .reset_index()
    )
    leaf["is_trace_leaf_05"] = (leaf["f_trace"] >= 0.5).astype(int)
    leaf["is_trace_leaf_03"] = (leaf["f_trace"] >= 0.3).astype(int)
    leaf.to_csv(OUT_DIR / "leaf_trace_fraction.csv", index=False)

    # Cohort summary per band
    summary_rows = []
    for b in BRAIN_BANDS_NAMES:
        sub = leaf[leaf["band"] == b]
        n_total_leaves = len(sub)
        per_pat = sub.groupby("patient").agg(
            n_leaves=("leaf_id", "count"),
            n_trace_05=("is_trace_leaf_05", "sum"),
            n_trace_03=("is_trace_leaf_03", "sum"),
            mean_f_trace=("f_trace", "mean"),
            median_f_trace=("f_trace", "median"),
        ).reset_index()
        per_pat["frac_trace_leaves_05"] = per_pat["n_trace_05"] / per_pat["n_leaves"]
        per_pat["frac_trace_leaves_03"] = per_pat["n_trace_03"] / per_pat["n_leaves"]
        per_pat.to_csv(OUT_DIR / f"per_patient_{b}.csv", index=False)
        # Cohort scalar: how many patients have frac_trace_leaves_03 above the cohort median
        # of the corresponding rest_pre half-baseline? Without that null we report descriptive
        # cohort distribution.
        summary_rows.append({
            "band": b,
            "n_total_leaves": n_total_leaves,
            "frac_trace_leaves_05_cohort_mean": float(per_pat["frac_trace_leaves_05"].mean()),
            "frac_trace_leaves_05_cohort_median": float(per_pat["frac_trace_leaves_05"].median()),
            "frac_trace_leaves_03_cohort_mean": float(per_pat["frac_trace_leaves_03"].mean()),
            "frac_trace_leaves_03_cohort_median": float(per_pat["frac_trace_leaves_03"].median()),
            "mean_f_trace_cohort_mean": float(per_pat["mean_f_trace"].mean()),
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_DIR / "cohort_summary.csv", index=False)
    print(summary.to_string(index=False))

    # Cohort distribution figure: per-band violin of leaf f_trace distribution
    fig, ax = plt.subplots(figsize=(9, 4.2))
    data = [leaf[leaf["band"] == b]["f_trace"].dropna().values for b in BRAIN_BANDS_NAMES]
    parts = ax.violinplot(data, positions=range(1, len(BRAIN_BANDS_NAMES) + 1),
                          showmeans=True, showmedians=False)
    for pc in parts["bodies"]:
        pc.set_facecolor("#9ec6e5")
        pc.set_edgecolor("#1a4f73")
        pc.set_alpha(0.8)
    ax.axhline(0.5, color="#d62728", lw=0.8, ls="--", label=r"$f_\mathrm{trace} \geq 0.5$ threshold")
    ax.axhline(0.3, color="#7f7f7f", lw=0.8, ls=":", label=r"$f_\mathrm{trace} \geq 0.3$ threshold")
    ax.set_xticks(range(1, len(BRAIN_BANDS_NAMES) + 1))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES])
    ax.set_ylabel(r"$f_\mathrm{trace}(\ell)$ -- fraction of $k$-grid where leaf is dominant-trace")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "cohort_distribution.pdf")
    plt.close(fig)
    print(f"[audit_34] outputs at {OUT_DIR}; figure at {FIG_DIR}")


if __name__ == "__main__":
    main()
