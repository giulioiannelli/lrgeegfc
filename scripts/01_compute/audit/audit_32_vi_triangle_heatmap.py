#!/usr/bin/env python3
"""Audit 32 — VI(k) phase-pair triangle as a multiscale band x k heatmap.

Re-mines ``data/reports/imcoh_vi/vi_raw_profiles.csv`` (already on disk).
Computes, per (patient, band, k):

    T_VI(p, b, k) = VI(c^TT, c^RPost; k) - VI(c^RPre, c^TT; k)

Negative T_VI = trace direction (test-tree closer to post than to pre at scale k).
No collapse over k. The result is reported as a per-patient heatmap and a cohort
``n_trace(b, k) = #{p : T_VI(p, b, k) < 0}`` heatmap. Same-k cross-patient
comparison is dirty (different |L_p|) -- flagged in the body of the deliverable.

Critical addition: singleton-share filter. At k close to N_p (the number of
contacts in patient p), most clusters are singletons and VI(c^a, c^b, k) becomes
a Hamming-like statistic on individual cluster ID assignments rather than a
module-structure comparison. We compute, per (p, b, k, phi), the fraction of
leaves that are in singleton clusters in phase phi at scale k. Cells where the
average across (RPre, TT, RPost) exceeds 0.5 are flagged ``singleton_dominant``;
they are NOT removed from the data but are visually annotated and reported
separately so the reader knows where the regime turns.

Outputs
-------
``data/audit/vi_triangle_heatmap/T_VI_per_patient_per_band_per_k.csv``
``data/audit/vi_triangle_heatmap/cohort_n_trace_band_k.csv``
``data/audit/vi_triangle_heatmap/singleton_share_band_k.csv``
``data/audit/vi_triangle_heatmap/cohort_band_summary.csv``
``data/outputs/figures/section_5_lrg_trace/vi_triangle_heatmap_cohort.pdf``
``data/outputs/figures/section_5_lrg_trace/vi_triangle_heatmap_smallmult.pdf``
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PATIENTS_LIST,
    PHASE_LABELS,
)
from lrg_eegfc.workflow.lrg import load_lrg_result

VI_PROFILES = ROOT / "data" / "reports" / "imcoh_vi" / "vi_raw_profiles.csv"
OUT_DIR = ROOT / "data" / "audit" / "vi_triangle_heatmap"
FIG_DIR = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace"


def singleton_share(Z, k):
    """Fraction of leaves that are in singleton clusters at the k-cluster cut."""
    labels = fcluster(Z, k, criterion="maxclust")
    _, counts = np.unique(labels, return_counts=True)
    n_leaves = labels.size
    return float((counts[counts == 1].sum()) / n_leaves)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(VI_PROFILES)
    pivot = (
        df[df["pair"].isin(["Pre-TT", "TT-Post"])]
        .pivot_table(index=["patient", "band", "k"], columns="pair", values="vi")
        .reset_index()
    )
    pivot["T_VI"] = pivot["TT-Post"] - pivot["Pre-TT"]
    pivot.to_csv(OUT_DIR / "T_VI_per_patient_per_band_per_k.csv", index=False)

    # ---- Singleton share per (p, b, k) averaged across (RPre, TT, RPost) phases.
    # This is the diagnostic that flags the fine-k regime where VI degenerates.
    sing_rows = []
    phases = ["rest_pre", "task_test", "rest_post"]
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            try:
                Z_phase = {}
                for phi in phases:
                    res = load_lrg_result(pat, phi, band, fc_method="imcoh_abs")
                    if res is None:
                        raise FileNotFoundError(f"missing LRG cache: {pat} {phi} {band}")
                    Z_phase[phi] = res.linkage_matrix
                N = Z_phase["rest_pre"].shape[0] + 1
                k_max = min(N - 1, int(pivot[(pivot.patient == pat) & (pivot.band == band)]["k"].max()))
                for k in range(2, k_max + 1):
                    shares = [singleton_share(Z_phase[phi], k) for phi in phases]
                    sing_rows.append({
                        "patient": pat,
                        "band": band,
                        "k": k,
                        "N_p": N,
                        "singleton_share_mean": float(np.mean(shares)),
                        "singleton_share_max": float(np.max(shares)),
                    })
            except Exception as e:
                print(f"[audit_32] WARN: {pat} {band}: {e}")
    if not sing_rows:
        raise RuntimeError("singleton-share computation produced zero rows; check load_lrg_result calls")
    sing_df = pd.DataFrame(sing_rows)
    sing_df.to_csv(OUT_DIR / "singleton_share_per_patient.csv", index=False)

    # Cohort-aggregate singleton share per (band, k): mean over patients.
    sing_band_k = (
        sing_df.groupby(["band", "k"])["singleton_share_mean"].mean().reset_index()
    )
    sing_band_k.to_csv(OUT_DIR / "singleton_share_band_k.csv", index=False)
    sing_lookup = sing_band_k.set_index(["band", "k"])["singleton_share_mean"].to_dict()

    # ---- Cohort n_trace(b, k)
    cohort = (
        pivot.assign(is_trace=(pivot["T_VI"] < 0).astype(int))
        .groupby(["band", "k"])["is_trace"]
        .sum()
        .reset_index(name="n_trace")
    )
    cohort["n_eligible"] = pivot.groupby(["band", "k"])["T_VI"].count().values
    cohort["singleton_share"] = [sing_lookup.get((b, k), np.nan) for b, k in zip(cohort.band, cohort.k)]
    cohort.to_csv(OUT_DIR / "cohort_n_trace_band_k.csv", index=False)

    # ---- Per-band summary, with singleton-mask split
    bands_ord = BRAIN_BANDS_NAMES
    summary_rows = []
    for b in bands_ord:
        sub = cohort[cohort["band"] == b].copy()
        n_total = len(sub)
        clean = sub[sub["singleton_share"] < 0.5]
        sing = sub[sub["singleton_share"] >= 0.5]
        summary_rows.append({
            "band": b,
            "k_grid_size": n_total,
            "n_clean_cells": len(clean),
            "n_singleton_dominant_cells": len(sing),
            "max_n_trace_clean": int(clean["n_trace"].max()) if len(clean) else 0,
            "argmax_k_clean": int(clean.loc[clean["n_trace"].idxmax(), "k"]) if len(clean) else -1,
            "max_n_trace_singleton": int(sing["n_trace"].max()) if len(sing) else 0,
            "argmax_k_singleton": int(sing.loc[sing["n_trace"].idxmax(), "k"]) if len(sing) else -1,
            "median_n_trace_clean": float(clean["n_trace"].median()) if len(clean) else float("nan"),
            "n_cells_clean_n_trace_ge_7": int((clean["n_trace"] >= 7).sum()),
            "n_cells_clean_n_trace_ge_8": int((clean["n_trace"] >= 8).sum()),
            "n_cells_clean_n_trace_ge_9": int((clean["n_trace"] >= 9).sum()),
        })
    pd.DataFrame(summary_rows).to_csv(OUT_DIR / "cohort_band_summary.csv", index=False)

    # ---- Cohort heatmap (with singleton overlay)
    k_max = int(cohort["k"].max())
    k_grid = np.arange(2, k_max + 1)

    def build_M(field):
        M = np.full((len(bands_ord), len(k_grid)), np.nan)
        for i, b in enumerate(bands_ord):
            sub = cohort[cohort["band"] == b].set_index("k")
            for j, k in enumerate(k_grid):
                if k in sub.index:
                    M[i, j] = sub.loc[k, field]
        return M

    M_cohort = build_M("n_trace")
    M_sing = build_M("singleton_share")

    fig = plt.figure(figsize=(11, 4.5))
    ax = fig.add_subplot(111)
    im = ax.imshow(
        M_cohort, aspect="auto", origin="lower", cmap="viridis", vmin=0, vmax=10,
        extent=(k_grid[0] - 0.5, k_grid[-1] + 0.5, -0.5, len(bands_ord) - 0.5),
    )
    # Singleton-dominant overlay: hatched where >=0.5
    Y, X = np.meshgrid(np.arange(len(bands_ord)), k_grid, indexing="ij")
    mask = M_sing >= 0.5
    ax.contourf(
        X, Y, mask.astype(float), levels=[0.5, 1.5], hatches=["///"], colors="none",
    )
    ax.set_yticks(range(len(bands_ord)))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in bands_ord])
    ax.set_xlabel(r"partition cut $k$")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label(r"$n_{\mathrm{trace}}(b, k)$ / 10")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "vi_triangle_heatmap_cohort.pdf")
    plt.close(fig)

    # ---- Small-multiples figure: cohort + 10 per-patient panels in one PDF
    pats_ord = sorted(pivot["patient"].unique())
    n_pat = len(pats_ord)
    fig = plt.figure(figsize=(13, 9))
    gs = gridspec.GridSpec(4, 5, figure=fig, height_ratios=[1.6, 1, 1, 0.05])
    # Top row: cohort heatmap spanning 5 cols
    ax_cohort = fig.add_subplot(gs[0, :])
    im_c = ax_cohort.imshow(
        M_cohort, aspect="auto", origin="lower", cmap="viridis", vmin=0, vmax=10,
        extent=(k_grid[0] - 0.5, k_grid[-1] + 0.5, -0.5, len(bands_ord) - 0.5),
    )
    ax_cohort.contourf(
        X, Y, mask.astype(float), levels=[0.5, 1.5], hatches=["///"], colors="none",
    )
    ax_cohort.set_yticks(range(len(bands_ord)))
    ax_cohort.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in bands_ord])
    ax_cohort.set_title(r"Cohort $n_{\mathrm{trace}}(b, k)$ -- hatch = singleton-dominant ($\geq 0.5$)")

    # Symmetric vmax for per-patient signed heatmaps (computed once across cohort)
    vmax_signed = float(np.nanpercentile(np.abs(pivot["T_VI"].values), 99))

    for idx, pat in enumerate(pats_ord):
        row = 1 + idx // 5
        col = idx % 5
        ax = fig.add_subplot(gs[row, col])
        sub = pivot[pivot["patient"] == pat]
        Mp = np.full((len(bands_ord), len(k_grid)), np.nan)
        Sp = np.full((len(bands_ord), len(k_grid)), np.nan)
        for i, b in enumerate(bands_ord):
            sb = sub[sub["band"] == b].set_index("k")
            ssh = sing_df[(sing_df.patient == pat) & (sing_df.band == b)].set_index("k")
            for j, k in enumerate(k_grid):
                if k in sb.index:
                    Mp[i, j] = sb.loc[k, "T_VI"]
                if k in ssh.index:
                    Sp[i, j] = ssh.loc[k, "singleton_share_mean"]
        ax.imshow(
            Mp, aspect="auto", origin="lower", cmap="RdBu_r",
            vmin=-vmax_signed, vmax=vmax_signed,
            extent=(k_grid[0] - 0.5, k_grid[-1] + 0.5, -0.5, len(bands_ord) - 0.5),
        )
        # Per-patient singleton overlay
        Yp, Xp = np.meshgrid(np.arange(len(bands_ord)), k_grid, indexing="ij")
        ax.contourf(
            Xp, Yp, (Sp >= 0.5).astype(float), levels=[0.5, 1.5],
            hatches=["///"], colors="none",
        )
        ax.set_title(pat, fontsize=9)
        if col == 0:
            ax.set_yticks(range(len(bands_ord)))
            ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in bands_ord], fontsize=7)
        else:
            ax.set_yticks([])
        if row == 2:
            ax.set_xlabel(r"$k$", fontsize=8)
        ax.tick_params(axis="x", labelsize=7)

    # Cohort colorbar (top axis)
    cbar_ax = fig.add_axes([0.86, 0.95, 0.10, 0.012])
    fig.colorbar(im_c, cax=cbar_ax, orientation="horizontal").set_label(
        r"$n_{\mathrm{trace}}/10$", fontsize=8,
    )

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIG_DIR / "vi_triangle_heatmap_smallmult.pdf")
    plt.close(fig)

    print(f"[audit_32] outputs at {OUT_DIR}")
    print(f"[audit_32] figures at {FIG_DIR}")
    print()
    print("Per-band summary (clean = singleton_share < 0.5):")
    print(pd.read_csv(OUT_DIR / "cohort_band_summary.csv").to_string(index=False))


if __name__ == "__main__":
    main()
