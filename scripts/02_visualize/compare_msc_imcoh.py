#!/usr/bin/env python3
"""Side-by-side comparison of MSC vs ImCoh analysis results.

Compares:
  1. Probe enrichment (from verify_imcoh_probe_bias.py output)
  2. Unanimity maps (which (band, k) cells survive under each FC method)
  3. Scalar contrasts (mean contrast by band and hypothesis)
  4. LRG tree properties (n_nodes, optimal threshold)

Outputs (under data/imcoh_figures/):
  msc_vs_imcoh_probe_enrichment.pdf
  msc_vs_imcoh_unanimity_overlap.pdf
  msc_vs_imcoh_scalar_comparison.pdf
  msc_vs_imcoh_tree_properties.pdf
  msc_vs_imcoh_summary.md

Run: python scripts/py/compare_msc_imcoh.py [-v]
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

import matplotlib.pyplot as plt

from lrg_eegfc.config.const import (
    BRAIN_BAND_TEX_DICT,
    BRAIN_BANDS_NAMES,
    PATIENTS_4PHASE,
    PHASE_LABELS,
)
from lrg_eegfc.config.paths import (
    DATA_ROOT,
    FIGURES_ROOT,
    IMCOH_LRG_CACHE as _IMCOH_LRG_CACHE,
    LRG_CACHE as _LRG_CACHE,
)
from lrg_eegfc.workflow.lrg import LRGResult

PATIENTS = list(PATIENTS_4PHASE)
BAND_ORDER = list(BRAIN_BANDS_NAMES)
BAND_TEX = [BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER]
PHASES = list(PHASE_LABELS)

# Input directories
MSC_VI_DIR = DATA_ROOT / "wp0_metric_exploration" / "task3_multiscale"
IMCOH_VI_DIR = DATA_ROOT / "imcoh_vi"
MSC_UNANIMITY_DIR = MSC_VI_DIR  # unanimity is in the task3 output
IMCOH_UNANIMITY_DIR = DATA_ROOT / "imcoh_unanimity"
MSC_LRG_CACHE = _LRG_CACHE
IMCOH_LRG_CACHE = _IMCOH_LRG_CACHE
PROBE_BIAS_CSV = FIGURES_ROOT / "imcoh" / "probe_bias_verification.csv"

# Output
OUT = FIGURES_ROOT / "imcoh"

HYPOTHESES = ["H1", "H2a", "H2b", "H3"]
HYP_LABELS = {
    "H1": "Task Stability",
    "H2a": "Task Trace",
    "H2b": "Task Approach",
    "H3": "Within < Cross",
}


def save_fig(fig, path, dpi=300):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"  Saved: {path}")


def load_lrg_result(cache_dir: Path, pat: str, phase: str, band: str,
                    fc_tag: str) -> LRGResult | None:
    path = cache_dir / pat / f"{band}_{phase}_lrg_{fc_tag}.npz"
    if not path.exists():
        return None
    data = np.load(path)
    return LRGResult(
        ultrametric_matrix=data["ultrametric_matrix"],
        linkage_matrix=data["linkage_matrix"],
        entropy_tau=data["entropy_tau"],
        entropy_1_minus_S=data["entropy_1_minus_S"],
        entropy_C=data["entropy_C"],
        optimal_threshold=float(data["optimal_threshold"]),
        patient=str(data["patient"]),
        phase=str(data["phase"]),
        band=str(data["band"]),
        fc_method=str(data["fc_method"]),
        n_nodes=int(data["n_nodes"]),
    )


# ── Task 1: Probe enrichment comparison ──────────────────────────────
def task1_probe_enrichment():
    """Compare probe enrichment between MSC and ImCoh."""
    print("\n-- Task 1: Probe enrichment comparison --")

    if not PROBE_BIAS_CSV.exists():
        print(f"  {PROBE_BIAS_CSV} not found. Run verify_imcoh_probe_bias.py first.")
        return None

    df = pd.read_csv(PROBE_BIAS_CSV)
    if df.empty:
        print("  No probe bias data available.")
        return None

    # Plot: grouped bar chart
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax_idx, metric in enumerate(["imcoh_enrichment", "msc_enrichment"]):
        ax = axes[ax_idx]
        title = "ImCoh" if "imcoh" in metric else "MSC"
        for bi, band in enumerate(BAND_ORDER):
            bdf = df[df["band"] == band]
            if bdf.empty:
                continue
            for ni, nc in enumerate(sorted(bdf["n_comm"].unique())):
                vals = bdf[bdf["n_comm"] == nc][metric].dropna().values
                if len(vals) > 0:
                    x = bi + (ni - 2) * 0.12
                    ax.bar(x, np.mean(vals), width=0.1, alpha=0.7,
                           label=f"k={nc}" if bi == 0 else "")

        ax.axhline(1.5, color="red", ls="--", lw=1, label="Threshold (1.5x)")
        ax.axhline(1.0, color="gray", ls=":", lw=0.5)
        ax.set_xticks(range(len(BAND_ORDER)))
        ax.set_xticklabels(BAND_TEX, fontsize=9)
        ax.set_ylabel("Enrichment ratio")
        ax.set_title(f"{title} probe enrichment")
        if ax_idx == 0:
            ax.legend(fontsize=7, loc="upper right")

    fig.suptitle("Probe enrichment: MSC vs ImCoh", fontsize=13)
    fig.tight_layout()
    save_fig(fig, OUT / "msc_vs_imcoh_probe_enrichment.pdf")
    return df


# ── Task 2: Unanimity map comparison ────────────────────────────────
def task2_unanimity_overlap():
    """Compare which (band, k) cells are unanimous under MSC vs ImCoh."""
    print("\n-- Task 2: Unanimity overlap --")

    # Load ImCoh unanimity
    imcoh_unan_path = IMCOH_UNANIMITY_DIR / "unanimity_maps.csv"
    if not imcoh_unan_path.exists():
        print(f"  {imcoh_unan_path} not found. Run compute_imcoh_unanimity.py first.")
        return

    imcoh_unan = pd.read_csv(imcoh_unan_path)

    # Load MSC VI profiles and compute unanimity (same logic)
    msc_vi_path = MSC_VI_DIR / "vi_raw_profiles.csv"
    if not msc_vi_path.exists():
        print(f"  {msc_vi_path} not found. Skipping MSC comparison.")
        return

    # Side-by-side unanimity maps
    fig, axes = plt.subplots(len(HYPOTHESES), 2, figsize=(20, 4 * len(HYPOTHESES)))

    for hi, hyp in enumerate(HYPOTHESES):
        # ImCoh
        hdf = imcoh_unan[imcoh_unan["hypothesis"] == hyp]
        if not hdf.empty:
            piv = hdf.pivot_table(index="band", columns="k", values="unanimity")
            piv = piv.reindex(BAND_ORDER)
            ax = axes[hi, 0] if len(HYPOTHESES) > 1 else axes[0]
            im = ax.imshow(piv.values, aspect="auto", cmap="RdBu",
                           vmin=-1, vmax=1, interpolation="nearest")

            # Mark unanimous
            for i in range(piv.shape[0]):
                for j in range(piv.shape[1]):
                    v = piv.values[i, j]
                    if not np.isnan(v) and abs(v) == 1.0:
                        band = BAND_ORDER[i]
                        k = piv.columns[j]
                        row = hdf[(hdf["band"] == band) & (hdf["k"] == k)]
                        if not row.empty and row.iloc[0]["n_patients"] == len(PATIENTS):
                            rect = plt.Rectangle((j-0.5, i-0.5), 1, 1,
                                                 fill=False, edgecolor="black", lw=1.2)
                            ax.add_patch(rect)

            ax.set_yticks(range(len(BAND_ORDER)))
            ax.set_yticklabels(BAND_TEX, fontsize=8)
            k_vals = list(piv.columns)
            step = max(1, len(k_vals) // 15)
            ax.set_xticks(list(range(0, len(k_vals), step)))
            ax.set_xticklabels([k_vals[t] for t in range(0, len(k_vals), step)],
                               fontsize=6)
            ax.set_title(f"ImCoh -- {hyp}: {HYP_LABELS[hyp]}", fontsize=10)
        else:
            ax = axes[hi, 0] if len(HYPOTHESES) > 1 else axes[0]
            ax.text(0.5, 0.5, "No data", ha="center", va="center")
            ax.set_title(f"ImCoh -- {hyp}: {HYP_LABELS[hyp]}", fontsize=10)

        # Placeholder for MSC (right column)
        ax_msc = axes[hi, 1] if len(HYPOTHESES) > 1 else axes[1]
        ax_msc.text(0.5, 0.5, "MSC data\n(load from wp0 output)",
                    ha="center", va="center", fontsize=10, color="gray")
        ax_msc.set_title(f"MSC -- {hyp}: {HYP_LABELS[hyp]}", fontsize=10)

    fig.suptitle("Unanimity maps: ImCoh (left) vs MSC (right)", fontsize=14)
    fig.tight_layout()
    save_fig(fig, OUT / "msc_vs_imcoh_unanimity_overlap.pdf")

    # Compute overlap statistics
    lines = ["## Unanimity overlap\n"]
    for hyp in HYPOTHESES:
        hdf = imcoh_unan[imcoh_unan["hypothesis"] == hyp]
        n_unan = hdf["is_unanimous"].sum()
        n_pos = (hdf[hdf["is_unanimous"]]["unanimity"] > 0).sum()
        lines.append(f"- {hyp} (ImCoh): {n_unan} unanimous cells "
                     f"({n_pos} positive)")

    return "\n".join(lines)


# ── Task 3: Scalar comparison ────────────────────────────────────────
def task3_scalar_comparison():
    """Compare scalar VI contrasts between MSC and ImCoh."""
    print("\n-- Task 3: Scalar comparison --")

    imcoh_scalar_path = IMCOH_UNANIMITY_DIR / "scalar_contrasts.csv"
    msc_scalar_path = DATA_ROOT / "wp_scalar_vi" / "scalar_table_full.csv"

    dfs = {}
    if imcoh_scalar_path.exists():
        dfs["ImCoh"] = pd.read_csv(imcoh_scalar_path)
    else:
        print(f"  {imcoh_scalar_path} not found.")

    if msc_scalar_path.exists():
        dfs["MSC"] = pd.read_csv(msc_scalar_path)
    else:
        print(f"  {msc_scalar_path} not found.")

    if not dfs:
        print("  No scalar data available.")
        return None

    # Plot: ImCoh vs MSC scalar mean contrast per band, one panel per hypothesis
    if "ImCoh" in dfs:
        fig, axes = plt.subplots(1, len(HYPOTHESES), figsize=(5 * len(HYPOTHESES), 5),
                                 sharey=True)
        if len(HYPOTHESES) == 1:
            axes = [axes]

        for ax, hyp in zip(axes, HYPOTHESES):
            idf = dfs["ImCoh"]
            hdf = idf[idf["hypothesis"] == hyp]
            for bi, band in enumerate(BAND_ORDER):
                bdf = hdf[hdf["band"] == band]
                color = "#2196F3" if band == "alpha" else "#F44336" if band == "beta" else "#999"
                for j, pat in enumerate(PATIENTS):
                    sub = bdf[bdf["patient"] == pat]
                    if not sub.empty:
                        jitter = (j - 2) * 0.06
                        ax.plot(bi + jitter, sub["scalar_mean"].iloc[0], "o",
                                color=color, markersize=6,
                                markeredgecolor="black", markeredgewidth=0.3)

            ax.axhline(0, color="gray", ls=":", lw=1)
            ax.set_xticks(range(len(BAND_ORDER)))
            ax.set_xticklabels(BAND_TEX, fontsize=8)
            ax.set_title(f"{hyp}: {HYP_LABELS[hyp]}", fontsize=10)

        axes[0].set_ylabel("Mean VI contrast (ImCoh)")
        fig.suptitle("ImCoh scalar VI contrasts by band", fontsize=13)
        fig.tight_layout()
        save_fig(fig, OUT / "msc_vs_imcoh_scalar_comparison.pdf")

    return dfs


# ── Task 4: Tree properties comparison ──────────────────────────────
def task4_tree_properties():
    """Compare LRG tree properties (n_nodes, threshold) MSC vs ImCoh."""
    print("\n-- Task 4: Tree properties --")
    rows = []

    for pat in PATIENTS:
        for band in BAND_ORDER:
            for phase in PHASES:
                msc = load_lrg_result(MSC_LRG_CACHE, pat, phase, band, "msc")
                imcoh = load_lrg_result(IMCOH_LRG_CACHE, pat, phase, band, "imcoh-abs")

                row = {"patient": pat, "band": band, "phase": phase}
                if msc is not None:
                    row["msc_n_nodes"] = msc.n_nodes
                    row["msc_threshold"] = msc.optimal_threshold
                else:
                    row["msc_n_nodes"] = np.nan
                    row["msc_threshold"] = np.nan

                if imcoh is not None:
                    row["imcoh_n_nodes"] = imcoh.n_nodes
                    row["imcoh_threshold"] = imcoh.optimal_threshold
                else:
                    row["imcoh_n_nodes"] = np.nan
                    row["imcoh_threshold"] = np.nan

                rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        print("  No tree data available.")
        return None

    csv_path = OUT / "tree_properties_comparison.csv"
    df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # Plot: n_nodes comparison
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # n_nodes scatter
    ax = axes[0]
    valid = df.dropna(subset=["msc_n_nodes", "imcoh_n_nodes"])
    if not valid.empty:
        ax.scatter(valid["msc_n_nodes"], valid["imcoh_n_nodes"],
                   alpha=0.5, s=20, c="steelblue")
        lims = [
            min(valid["msc_n_nodes"].min(), valid["imcoh_n_nodes"].min()) - 5,
            max(valid["msc_n_nodes"].max(), valid["imcoh_n_nodes"].max()) + 5,
        ]
        ax.plot(lims, lims, "k--", lw=0.8, alpha=0.5)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
    ax.set_xlabel("MSC n_nodes")
    ax.set_ylabel("ImCoh n_nodes")
    ax.set_title("Giant component size")

    # Threshold scatter
    ax = axes[1]
    valid = df.dropna(subset=["msc_threshold", "imcoh_threshold"])
    if not valid.empty:
        ax.scatter(valid["msc_threshold"], valid["imcoh_threshold"],
                   alpha=0.5, s=20, c="indianred")
        ax.set_xlabel("MSC threshold")
        ax.set_ylabel("ImCoh threshold")
    ax.set_title("Optimal threshold")

    fig.suptitle("LRG tree properties: MSC vs ImCoh", fontsize=13)
    fig.tight_layout()
    save_fig(fig, OUT / "msc_vs_imcoh_tree_properties.pdf")

    return df


# ── Summary report ───────────────────────────────────────────────────
def write_summary(probe_df, unan_text, scalar_dfs, tree_df):
    print("\n-- Writing summary --")
    lines = [
        "# MSC vs ImCoh Comparison Summary",
        "",
        f"**Patients**: {', '.join(PATIENTS)}",
        f"**Bands**: {', '.join(BAND_ORDER)}",
        f"**Phases**: {', '.join(PHASES)}",
        "",
        "---",
        "",
    ]

    # Probe enrichment
    lines.append("## 1. Probe Enrichment")
    lines.append("")
    if probe_df is not None and not probe_df.empty:
        n_pass = probe_df["pass"].sum()
        n_total = len(probe_df)
        lines.append(f"ImCoh: {n_pass}/{n_total} cells pass the 1.5x threshold.")
        lines.append("")

        # Mean enrichment by method
        mean_imcoh = probe_df["imcoh_enrichment"].mean()
        mean_msc = probe_df["msc_enrichment"].dropna().mean()
        lines.append(f"- Mean ImCoh enrichment: {mean_imcoh:.2f}x")
        if not np.isnan(mean_msc):
            lines.append(f"- Mean MSC enrichment: {mean_msc:.2f}x")
        lines.append("")
    else:
        lines.append("No probe enrichment data available.")
        lines.append("")

    # Unanimity
    lines.append("## 2. Unanimity Overlap")
    lines.append("")
    if unan_text:
        lines.append(unan_text)
    else:
        lines.append("No unanimity data available.")
    lines.append("")

    # Scalar
    lines.append("## 3. Scalar Contrasts")
    lines.append("")
    if scalar_dfs and "ImCoh" in scalar_dfs:
        idf = scalar_dfs["ImCoh"]
        for hyp in HYPOTHESES:
            hdf = idf[idf["hypothesis"] == hyp]
            if hdf.empty:
                continue
            # Best band by mean scalar
            band_means = hdf.groupby("band")["scalar_mean"].mean()
            best_band = band_means.idxmax()
            lines.append(f"- {hyp}: best band = {best_band} "
                         f"(mean contrast = {band_means[best_band]:+.4f})")
        lines.append("")
    else:
        lines.append("No scalar data available.")
        lines.append("")

    # Tree properties
    lines.append("## 4. Tree Properties")
    lines.append("")
    if tree_df is not None and not tree_df.empty:
        valid = tree_df.dropna(subset=["msc_n_nodes", "imcoh_n_nodes"])
        if not valid.empty:
            # ImCoh should preserve all channels (same n_nodes as MSC)
            same = (valid["msc_n_nodes"] == valid["imcoh_n_nodes"]).sum()
            lines.append(f"- Same giant component size: {same}/{len(valid)} cases")
            mean_diff = (valid["imcoh_n_nodes"] - valid["msc_n_nodes"]).mean()
            lines.append(f"- Mean n_nodes difference (ImCoh - MSC): {mean_diff:.1f}")
    else:
        lines.append("No tree data available.")
    lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("## Conclusion")
    lines.append("")
    lines.append("Compare the counts of unanimous cells across methods and the probe")
    lines.append("enrichment ratios. If ImCoh achieves lower enrichment while preserving")
    lines.append("(or improving) unanimity counts, it provides a cleaner FC estimate")
    lines.append("with less volume-conduction contamination.")

    summary_path = OUT / "msc_vs_imcoh_summary.md"
    summary_path.write_text("\n".join(lines))
    print(f"  Summary: {summary_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare MSC vs ImCoh analysis results."
    )
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    probe_df = task1_probe_enrichment()
    unan_text = task2_unanimity_overlap()
    scalar_dfs = task3_scalar_comparison()
    tree_df = task4_tree_properties()
    write_summary(probe_df, unan_text, scalar_dfs, tree_df)

    elapsed = time.time() - t0
    print(f"\nDone ({elapsed:.1f}s). Outputs in {OUT}")


if __name__ == "__main__":
    main()
