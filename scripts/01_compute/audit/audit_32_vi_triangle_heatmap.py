#!/usr/bin/env python3
"""Audit 32 — VI(k) phase-pair triangle as a multiscale band × k heatmap.

Re-mines ``data/reports/imcoh_vi/vi_raw_profiles.csv`` (already on disk).
Computes, per (patient, band, k):

    T_VI(p, b, k) = VI(c^TT, c^RPost; k) − VI(c^RPre, c^TT; k)

Negative T_VI = trace direction (test-tree closer to post than to pre at scale k).
No collapse over k. The result is reported as a per-patient heatmap and a cohort
``n_trace(b, k) = #{p : T_VI(p, b, k) < 0}`` heatmap. Same-k cross-patient
comparison is dirty (different |L_p|) — flagged in the body of the deliverable.

Outputs
-------
``data/audit/vi_triangle_heatmap/T_VI_per_patient_per_band_per_k.csv``
``data/audit/vi_triangle_heatmap/cohort_n_trace_band_k.csv``
``data/outputs/figures/section_5_lrg_trace/vi_triangle_heatmap_cohort.pdf``
``data/outputs/figures/section_5_lrg_trace/vi_triangle_heatmap_per_patient/<Pat_XX>.pdf``
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

VI_PROFILES = ROOT / "data" / "reports" / "imcoh_vi" / "vi_raw_profiles.csv"
OUT_DIR = ROOT / "data" / "audit" / "vi_triangle_heatmap"
FIG_DIR = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    (FIG_DIR / "vi_triangle_heatmap_per_patient").mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(VI_PROFILES)
    # Pivot to (patient, band, k) → {Pre-TT: vi, TT-Post: vi}
    pivot = (
        df[df["pair"].isin(["Pre-TT", "TT-Post"])]
        .pivot_table(index=["patient", "band", "k"], columns="pair", values="vi")
        .reset_index()
    )
    pivot["T_VI"] = pivot["TT-Post"] - pivot["Pre-TT"]
    pivot.to_csv(OUT_DIR / "T_VI_per_patient_per_band_per_k.csv", index=False)

    # Cohort n_trace(b, k)
    cohort = (
        pivot.assign(is_trace=(pivot["T_VI"] < 0).astype(int))
        .groupby(["band", "k"])["is_trace"]
        .sum()
        .reset_index(name="n_trace")
    )
    cohort["n_eligible"] = pivot.groupby(["band", "k"])["T_VI"].count().values
    cohort.to_csv(OUT_DIR / "cohort_n_trace_band_k.csv", index=False)

    # Cohort heatmap: 6 bands × k axis, color = n_trace / 10
    bands_ord = BRAIN_BANDS_NAMES
    k_max = int(cohort["k"].max())
    k_grid = np.arange(2, k_max + 1)
    M = np.full((len(bands_ord), len(k_grid)), np.nan)
    for i, b in enumerate(bands_ord):
        sub = cohort[cohort["band"] == b].set_index("k")
        for j, k in enumerate(k_grid):
            if k in sub.index:
                M[i, j] = sub.loc[k, "n_trace"]
    fig, ax = plt.subplots(figsize=(10, 3.5))
    im = ax.imshow(M, aspect="auto", origin="lower", cmap="viridis", vmin=0, vmax=10,
                   extent=(k_grid[0] - 0.5, k_grid[-1] + 0.5, -0.5, len(bands_ord) - 0.5))
    ax.set_yticks(range(len(bands_ord)))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in bands_ord])
    ax.set_xlabel(r"partition cut $k$")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label(r"$n_{\mathrm{trace}}(b, k)$ / 10")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "vi_triangle_heatmap_cohort.pdf")
    plt.close(fig)

    # Per-patient signed heatmaps
    for pat, sub in pivot.groupby("patient"):
        M_p = np.full((len(bands_ord), len(k_grid)), np.nan)
        for i, b in enumerate(bands_ord):
            sb = sub[sub["band"] == b].set_index("k")
            for j, k in enumerate(k_grid):
                if k in sb.index:
                    M_p[i, j] = sb.loc[k, "T_VI"]
        vmax = float(np.nanmax(np.abs(M_p))) if np.isfinite(np.nanmax(np.abs(M_p))) else 1.0
        fig, ax = plt.subplots(figsize=(10, 3.5))
        im = ax.imshow(M_p, aspect="auto", origin="lower", cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax,
                       extent=(k_grid[0] - 0.5, k_grid[-1] + 0.5, -0.5, len(bands_ord) - 0.5))
        ax.set_yticks(range(len(bands_ord)))
        ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in bands_ord])
        ax.set_xlabel(r"partition cut $k$")
        cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
        cbar.set_label(r"$T_{\mathrm{VI}}(p, b, k)$ — negative = trace")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "vi_triangle_heatmap_per_patient" / f"{pat}.pdf")
        plt.close(fig)

    # Compact cohort summary by band: integrated n_trace + frac of k-grid where ≥7/10
    summary_rows = []
    for b in bands_ord:
        sub = cohort[cohort["band"] == b]
        n_total_cells = len(sub)
        n_cells_ge7 = int((sub["n_trace"] >= 7).sum())
        n_cells_ge8 = int((sub["n_trace"] >= 8).sum())
        median_n_trace = float(sub["n_trace"].median())
        max_n_trace = int(sub["n_trace"].max())
        argmax_k = int(sub.loc[sub["n_trace"].idxmax(), "k"]) if len(sub) else -1
        summary_rows.append({
            "band": b,
            "k_grid_size": n_total_cells,
            "median_n_trace_over_k": median_n_trace,
            "max_n_trace": max_n_trace,
            "argmax_k": argmax_k,
            "n_cells_n_trace_ge_7": n_cells_ge7,
            "n_cells_n_trace_ge_8": n_cells_ge8,
            "frac_k_ge_7": n_cells_ge7 / max(n_total_cells, 1),
        })
    pd.DataFrame(summary_rows).to_csv(OUT_DIR / "cohort_band_summary.csv", index=False)

    print(f"[audit_32] wrote {OUT_DIR}")
    print(f"[audit_32] wrote cohort heatmap → {FIG_DIR / 'vi_triangle_heatmap_cohort.pdf'}")
    print(f"[audit_32] wrote {len(pivot.patient.unique())} per-patient PDFs")


if __name__ == "__main__":
    main()
