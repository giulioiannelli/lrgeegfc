#!/usr/bin/env python3
"""Surrogate validation: verify whether circular-shift surrogate filtering
is necessary for the MSC-based functional connectivity pipeline.

Produces:
  data/surrogate_validation/fig_surrogate_null_vs_empirical_Pat_02_rsPre.pdf
  data/surrogate_validation/fig_adjacency_raw_vs_hard_Pat_02_rsPre.pdf
  data/surrogate_validation/fig_surrogate_summary_heatmap.pdf
  data/surrogate_validation/surrogate_summary_table.csv
  data/surrogate_validation/surrogate_summary_report.txt
"""

from __future__ import annotations

import gc
import time
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    BRAIN_BANDS_NAMES,
    BRAIN_BANDS_TEX_NAMES,
    BRAIN_BAND_TEX_DICT,
    nperseg_for_fs,
    sEEG_DATAPATH,
)
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch, band_average_msc
from lrg_eegfc.utils.fc.msc.surrogates import surrogate_msc_null
from lrg_eegfc.utils.io.patient_robust import load_patient_dataset_robust

# ── constants ────────────────────────────────────────────────────────────────
from lrg_eegfc.config.paths import SURROGATE_VALIDATION_CACHE
OUTPUT_DIR = SURROGATE_VALIDATION_CACHE
N_SURROGATES = 100
PHASE = "rest_pre"
REFERENCE_PATIENT = "Pat_02"
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
BAND_NAMES = BRAIN_BANDS_NAMES  # ordered list

# ── helpers ──────────────────────────────────────────────────────────────────

def _upper_tri(W: np.ndarray) -> np.ndarray:
    """Return flattened upper triangle (excluding diagonal)."""
    idx = np.triu_indices_from(W, k=1)
    return W[idx]


def _load_timeseries(patient: str, phase: str):
    """Load timeseries + sampling rate for a patient/phase."""
    recordings = load_patient_dataset_robust(patient, phases=[phase])
    rec = recordings[phase]
    fs = float(rec.parameters.get("fs", 2048.0))
    return rec.timeseries, fs


def _compute_empirical_msc(X: np.ndarray, fs: float, nperseg: int):
    """Compute dense MSC matrices for all 6 bands."""
    freqs, Coh = compute_msc_welch(X, fs, nperseg=nperseg, noverlap=nperseg // 2)
    W_bands = band_average_msc(Coh, freqs, BRAIN_BANDS)
    del Coh
    gc.collect()
    return W_bands


def _compute_null(X, fs, nperseg, n_surrogates=N_SURROGATES, n_workers=None):
    """Compute surrogate null distribution for all bands."""
    return surrogate_msc_null(
        X, fs, BRAIN_BANDS, n_surrogates,
        nperseg=nperseg, noverlap=nperseg // 2,
        n_workers=n_workers,
    )


# ═════════════════════════════════════════════════════════════════════════════
#  TASK 1: Null distribution vs empirical MSC — one patient, all bands
# ═════════════════════════════════════════════════════════════════════════════

def task1(W_bands, W_null, report_lines):
    """Analyse and plot empirical vs surrogate null for Pat_02, rest_pre."""
    print("\n══════ TASK 1: Null distribution vs empirical MSC ══════")
    report_lines.append("=" * 70)
    report_lines.append("TASK 1: Null distribution vs empirical MSC (Pat_02, rest_pre)")
    report_lines.append("=" * 70)

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    axes = axes.ravel()

    for idx, band in enumerate(BAND_NAMES):
        W = W_bands[band].copy()
        np.fill_diagonal(W, 0.0)
        Wn = W_null[band]  # (n_surr, N, N)

        emp_vals = _upper_tri(W)
        N = W.shape[0]

        # Pooled surrogate values (all edges, all surrogates)
        triu_i, triu_j = np.triu_indices(N, k=1)
        surr_pooled = Wn[:, triu_i, triu_j].ravel()

        # Per-band statistics
        surr_mean = surr_pooled.mean()
        surr_std = surr_pooled.std()
        emp_mean = emp_vals.mean()

        # Edge-specific thresholds (percentiles across surrogates)
        surr_per_edge = Wn[:, triu_i, triu_j]  # (n_surr, n_edges)
        p95 = np.percentile(surr_per_edge, 95, axis=0)
        p99 = np.percentile(surr_per_edge, 99, axis=0)
        frac_above_95 = np.mean(emp_vals > p95)
        frac_above_99 = np.mean(emp_vals > p99)

        line = (
            f"  {band:>12s}:  surr mean={surr_mean:.4f} ± {surr_std:.4f}  |  "
            f"emp mean={emp_mean:.4f}  |  >p95={frac_above_95:.3f}  >p99={frac_above_99:.3f}"
        )
        print(line)
        report_lines.append(line)

        # Subplot
        ax = axes[idx]
        bins = np.linspace(0, max(emp_vals.max(), surr_pooled.max()) * 1.05, 80)
        ax.hist(surr_pooled, bins=bins, density=True, color="0.6", alpha=0.6, label="Surrogate null")
        ax.hist(emp_vals, bins=bins, density=True, color="steelblue", alpha=0.7, label="Empirical")
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=13)
        ax.set_xlabel("MSC")
        if idx % 3 == 0:
            ax.set_ylabel("Density")
        ax.legend(fontsize=8)

    fig.suptitle("Empirical MSC vs Surrogate Null — Pat_02, rest_pre", fontsize=14, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = OUTPUT_DIR / "fig_surrogate_null_vs_empirical_Pat_02_rsPre.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  → saved {out}")
    report_lines.append(f"  Figure saved: {out}")
    report_lines.append("")


# ═════════════════════════════════════════════════════════════════════════════
#  TASK 2: Effect on adjacency matrix structure — one patient
# ═════════════════════════════════════════════════════════════════════════════

def task2(W_bands, W_null, report_lines):
    """Hard/soft thresholding comparison for Pat_02, rest_pre."""
    print("\n══════ TASK 2: Effect on adjacency matrix structure ══════")
    report_lines.append("=" * 70)
    report_lines.append("TASK 2: Effect on adjacency matrix structure (Pat_02, rest_pre)")
    report_lines.append("=" * 70)

    results = {}

    for band in BAND_NAMES:
        W = W_bands[band].copy()
        np.fill_diagonal(W, 0.0)
        Wn = W_null[band]
        n_surr = Wn.shape[0]

        # Empirical p-values: p_ij = (1 + #{surr >= emp}) / (n_surr + 1)
        p_values = (1 + np.sum(Wn >= W[None, :, :], axis=0)) / (n_surr + 1)

        # Hard threshold
        W_hard = np.where(p_values < 0.05, W, 0.0)
        np.fill_diagonal(W_hard, 0.0)

        # Soft mask
        W_soft = W * (1.0 - p_values)
        np.fill_diagonal(W_soft, 0.0)

        # Metrics
        ut = _upper_tri
        frac_retained = np.mean(ut(W_hard) > 0)
        norm_W = np.linalg.norm(ut(W))
        frob_hard = np.linalg.norm(ut(W) - ut(W_hard)) / norm_W if norm_W > 0 else 0.0
        frob_soft = np.linalg.norm(ut(W) - ut(W_soft)) / norm_W if norm_W > 0 else 0.0
        rho, _ = spearmanr(ut(W), ut(W_soft))

        results[band] = {
            "frac_retained": frac_retained,
            "frob_hard": frob_hard,
            "frob_soft": frob_soft,
            "spearman": rho,
            "W": W,
            "W_hard": W_hard,
            "W_soft": W_soft,
        }

        line = (
            f"  {band:>12s}:  retained={frac_retained:.3f}  "
            f"‖W-W_hard‖/‖W‖={frob_hard:.4f}  "
            f"‖W-W_soft‖/‖W‖={frob_soft:.4f}  "
            f"ρ(W,W_soft)={rho:.4f}"
        )
        print(line)
        report_lines.append(line)

    # Figure: 2 rows × 6 cols (top=W, bottom=W_hard)
    fig, axes = plt.subplots(2, 6, figsize=(24, 8))
    vmin, vmax = 0, 0.6

    for col, band in enumerate(BAND_NAMES):
        r = results[band]
        im0 = axes[0, col].imshow(r["W"], vmin=vmin, vmax=vmax, cmap="viridis", aspect="equal")
        axes[0, col].set_title(f"{BRAIN_BAND_TEX_DICT[band]}\nW (raw)", fontsize=11)
        axes[0, col].set_xticks([]); axes[0, col].set_yticks([])

        im1 = axes[1, col].imshow(r["W_hard"], vmin=vmin, vmax=vmax, cmap="viridis", aspect="equal")
        axes[1, col].set_title(f"W_hard (p<0.05)\nret={r['frac_retained']:.2f}", fontsize=10)
        axes[1, col].set_xticks([]); axes[1, col].set_yticks([])

    fig.colorbar(im0, ax=axes[0, :].tolist(), shrink=0.7, label="MSC")
    fig.colorbar(im1, ax=axes[1, :].tolist(), shrink=0.7, label="MSC")
    fig.suptitle("Raw vs Hard-Thresholded Adjacency — Pat_02, rest_pre", fontsize=14, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = OUTPUT_DIR / "fig_adjacency_raw_vs_hard_Pat_02_rsPre.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  → saved {out}")
    report_lines.append(f"  Figure saved: {out}")
    report_lines.append("")

    return results


# ═════════════════════════════════════════════════════════════════════════════
#  TASK 3: Effect on LRG dendrogram — quick check
# ═════════════════════════════════════════════════════════════════════════════

def task3(task2_results, report_lines):
    """LRG dendrogram comparison for the two extreme bands."""
    import networkx as nx
    from lrgsglib.core import (
        compute_laplacian_properties,
        compute_normalized_linkage,
        compute_optimal_threshold,
        get_giant_component,
    )
    from scipy.cluster.hierarchy import cophenet

    print("\n══════ TASK 3: Effect on LRG dendrogram ══════")
    report_lines.append("=" * 70)
    report_lines.append("TASK 3: Effect on LRG dendrogram — extreme bands (Pat_02, rest_pre)")
    report_lines.append("=" * 70)

    # Identify bands with largest / smallest Spearman(W, W_soft)
    spearman_vals = {band: task2_results[band]["spearman"] for band in BAND_NAMES}
    band_min = min(spearman_vals, key=spearman_vals.get)  # largest effect
    band_max = max(spearman_vals, key=spearman_vals.get)  # smallest effect

    report_lines.append(
        f"  Largest surrogate effect:  {band_min} (ρ = {spearman_vals[band_min]:.4f})"
    )
    report_lines.append(
        f"  Smallest surrogate effect: {band_max} (ρ = {spearman_vals[band_max]:.4f})"
    )
    print(f"  Largest effect:  {band_min} (ρ = {spearman_vals[band_min]:.4f})")
    print(f"  Smallest effect: {band_max} (ρ = {spearman_vals[band_max]:.4f})")

    def _run_lrg(W):
        """Run LRG on adjacency matrix W, return linkage and n*."""
        np.fill_diagonal(W, 0.0)
        G = nx.from_numpy_array(W)
        giant = get_giant_component(G)
        _, _, _, Trho, _ = compute_laplacian_properties(giant)
        dists = squareform(Trho)
        linkage, labels, _ = compute_normalized_linkage(dists, giant)
        threshold, *_ = compute_optimal_threshold(linkage)
        return linkage, threshold, giant.number_of_nodes()

    for band in [band_min, band_max]:
        r = task2_results[band]
        print(f"\n  Band: {band}")
        report_lines.append(f"\n  Band: {band}")

        link_W, thresh_W, n_W = _run_lrg(r["W"].copy())
        link_Ws, thresh_Ws, n_Ws = _run_lrg(r["W_soft"].copy())

        print(f"    W:      n*={n_W}, threshold={thresh_W:.4f}")
        print(f"    W_soft: n*={n_Ws}, threshold={thresh_Ws:.4f}")
        report_lines.append(f"    W:      n_nodes={n_W}, ψ-optimal threshold={thresh_W:.4f}")
        report_lines.append(f"    W_soft: n_nodes={n_Ws}, ψ-optimal threshold={thresh_Ws:.4f}")

        # Cophenetic correlation between the two dendrograms
        # Both linkage matrices must describe the same number of nodes
        if n_W == n_Ws:
            d_W = squareform(cophenet(link_W))
            d_Ws = squareform(cophenet(link_Ws))
            coph_rho, _ = spearmanr(_upper_tri(d_W), _upper_tri(d_Ws))
            print(f"    Cophenetic correlation (W vs W_soft): {coph_rho:.4f}")
            report_lines.append(f"    Cophenetic correlation (W vs W_soft): {coph_rho:.4f}")
        else:
            # Giant components differ in size — compare on intersection
            print(f"    Giant components differ in size ({n_W} vs {n_Ws}), "
                  "cophenetic correlation not directly comparable.")
            report_lines.append(
                f"    Giant components differ ({n_W} vs {n_Ws}), "
                "cophenetic correlation skipped."
            )

    report_lines.append("")


# ═════════════════════════════════════════════════════════════════════════════
#  TASK 4: Cross-patient, cross-band summary
# ═════════════════════════════════════════════════════════════════════════════

def task4(report_lines, n_workers=None):
    """Repeat key metrics for all 5 patients × 6 bands, rest_pre only."""
    print("\n══════ TASK 4: Cross-patient summary ══════")
    report_lines.append("=" * 70)
    report_lines.append("TASK 4: Cross-patient, cross-band summary (rest_pre)")
    report_lines.append("=" * 70)

    rows = []

    for patient in ALL_PATIENTS:
        print(f"\n  Processing {patient}...")
        X, fs = _load_timeseries(patient, PHASE)
        nperseg = nperseg_for_fs(fs)
        print(f"    fs={fs}, nperseg={nperseg}, shape={X.shape}")

        W_bands = _compute_empirical_msc(X, fs, nperseg)
        W_null = _compute_null(X, fs, nperseg, n_workers=n_workers)

        for band in BAND_NAMES:
            W = W_bands[band].copy()
            np.fill_diagonal(W, 0.0)
            Wn = W_null[band]
            n_surr = Wn.shape[0]

            # Empirical p-values
            p_values = (1 + np.sum(Wn >= W[None, :, :], axis=0)) / (n_surr + 1)

            # Fraction significant
            frac_p05 = np.mean(_upper_tri(p_values) < 0.05)
            frac_p01 = np.mean(_upper_tri(p_values) < 0.01)

            # Spearman(W, W_soft)
            W_soft = W * (1.0 - p_values)
            np.fill_diagonal(W_soft, 0.0)
            rho, _ = spearmanr(_upper_tri(W), _upper_tri(W_soft))

            rows.append({
                "patient": patient,
                "band": band,
                "frac_significant_p05": frac_p05,
                "frac_significant_p01": frac_p01,
                "spearman_W_vs_Wsoft": rho,
            })

            line = f"    {band:>12s}: sig_p05={frac_p05:.3f}  sig_p01={frac_p01:.3f}  ρ={rho:.4f}"
            print(line)
            report_lines.append(line)

        del X, W_bands, W_null
        gc.collect()

    # Save CSV
    df = pd.DataFrame(rows)
    csv_path = OUTPUT_DIR / "surrogate_summary_table.csv"
    df.to_csv(csv_path, index=False)
    print(f"\n  → saved {csv_path}")
    report_lines.append(f"\n  Table saved: {csv_path}")

    # Heatmap figure
    pivot = df.pivot(index="patient", columns="band", values="frac_significant_p05")
    pivot = pivot[BAND_NAMES]  # enforce column order

    fig, ax = plt.subplots(figsize=(9, 5))
    im = ax.imshow(pivot.values, cmap="YlOrRd", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(BAND_NAMES)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_NAMES], fontsize=12)
    ax.set_yticks(range(len(ALL_PATIENTS)))
    ax.set_yticklabels(ALL_PATIENTS, fontsize=11)
    ax.set_xlabel("Frequency band")
    ax.set_ylabel("Patient")
    ax.set_title("Fraction of edges significant at p < 0.05\n(surrogate test, rest_pre)", fontsize=13)

    # Annotate cells
    for i in range(len(ALL_PATIENTS)):
        for j in range(len(BAND_NAMES)):
            val = pivot.values[i, j]
            color = "white" if val > 0.7 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=10, color=color)

    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    out = OUTPUT_DIR / "fig_surrogate_summary_heatmap.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  → saved {out}")
    report_lines.append(f"  Figure saved: {out}")

    return df


# ═════════════════════════════════════════════════════════════════════════════
#  DECISION CRITERION
# ═════════════════════════════════════════════════════════════════════════════

def decision_criterion(df, report_lines):
    """Classify each (patient, band) and print final decision."""
    print("\n══════ DECISION CRITERION ══════")
    report_lines.append("")
    report_lines.append("=" * 70)
    report_lines.append("DECISION CRITERION")
    report_lines.append("=" * 70)

    cat_A, cat_B, cat_C = 0, 0, 0

    for _, row in df.iterrows():
        rho = row["spearman_W_vs_Wsoft"]
        frac = row["frac_significant_p05"]
        patient = row["patient"]
        band = row["band"]

        if rho > 0.95:
            cat = "A"
            cat_A += 1
        elif frac > 0.90:
            cat = "B"
            cat_B += 1
        else:
            cat = "C"
            cat_C += 1

        line = f"  {patient} {band:>12s}:  ρ={rho:.4f}  frac_sig={frac:.3f}  → ({cat})"
        report_lines.append(line)

    total = cat_A + cat_B + cat_C
    summary = [
        "",
        f"  Category (A) ρ > 0.95 — filtering negligible:           {cat_A}/{total}",
        f"  Category (B) frac_sig > 0.90 — nearly all edges above null: {cat_B}/{total}",
        f"  Category (C) neither — structural impact:                  {cat_C}/{total}",
        "",
    ]

    if cat_C == 0:
        summary.append(
            "  CONCLUSION: All (patient, band) cells fall into category (A) or (B)."
        )
        summary.append(
            "  Surrogates are omitted. Justification: circular-shift surrogate filtering"
        )
        summary.append(
            "  has negligible effect on MSC adjacency structure (Spearman > 0.95 or"
        )
        summary.append(
            "  >90% edges significant) across all patients and frequency bands."
        )
    elif cat_C <= 3:
        summary.append(
            f"  CONCLUSION: {cat_C} cells show structural impact — mostly negligible."
        )
        summary.append(
            "  Consider adding a brief note to the manuscript acknowledging minor"
        )
        summary.append(
            "  differences in the affected bands."
        )
    else:
        summary.append(
            f"  CONCLUSION: {cat_C} cells show structural impact."
        )
        summary.append(
            "  Surrogate filtering should be added to the pipeline for affected bands."
        )

    for s in summary:
        print(s)
        report_lines.append(s)


# ═════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Surrogate validation for MSC pipeline")
    parser.add_argument("--workers", type=int, default=None,
                        help="Number of parallel workers for surrogate computation")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    report_lines: list[str] = [
        "SURROGATE VALIDATION REPORT",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"N_surrogates: {N_SURROGATES}",
        f"Phase: {PHASE}",
        "",
    ]

    t0 = time.time()

    # ── Tasks 1–3: Pat_02 analysis ──────────────────────────────────────────
    print("Loading Pat_02 rest_pre data...")
    X, fs = _load_timeseries(REFERENCE_PATIENT, PHASE)
    nperseg = nperseg_for_fs(fs)
    print(f"  fs={fs}, nperseg={nperseg}, X.shape={X.shape}")
    report_lines.append(f"Pat_02: fs={fs} Hz, nperseg={nperseg}, shape={X.shape}")
    report_lines.append(f"Welch: hann window, noverlap={nperseg // 2}, Δf={fs / nperseg:.2f} Hz")
    report_lines.append("")

    print("Computing empirical MSC for Pat_02...")
    W_bands = _compute_empirical_msc(X, fs, nperseg)

    print(f"Computing {N_SURROGATES} surrogate null distributions for Pat_02...")
    W_null = _compute_null(X, fs, nperseg, n_workers=args.workers)
    del X
    gc.collect()

    task1(W_bands, W_null, report_lines)
    task2_results = task2(W_bands, W_null, report_lines)
    task3(task2_results, report_lines)

    del W_bands, W_null
    gc.collect()

    # ── Task 4: Cross-patient summary ───────────────────────────────────────
    df = task4(report_lines, n_workers=args.workers)
    decision_criterion(df, report_lines)

    elapsed = time.time() - t0
    report_lines.append(f"\nTotal runtime: {elapsed / 60:.1f} min")
    print(f"\nTotal runtime: {elapsed / 60:.1f} min")

    # Save report
    report_path = OUTPUT_DIR / "surrogate_summary_report.txt"
    report_path.write_text("\n".join(report_lines))
    print(f"Report saved: {report_path}")


if __name__ == "__main__":
    main()
