#!/usr/bin/env python3
"""Section 5.2 — Metric Concordance Analysis (all patients).

Computes all 9 phase-reorganization metrics between every pair of cognitive
phases, across all 6 frequency bands, for each patient with LRG cache.
Then computes the inter-metric Spearman correlation matrix to identify
redundant vs complementary metrics.

Outputs per patient:
  - fig_metric_concordance_{pat}.pdf   (9x9 Spearman heatmap + dendrogram)
  - fig_metric_values_grid_{pat}.pdf   (3x3 strip plot, one panel per metric)
  - metric_concordance_{pat}.csv       (all metric values)
  - metric_correlation_matrix_{pat}.csv
  - metric_clusters_{pat}.csv

Run: python scripts/gen_metric_concordance.py
"""
import sys
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram, leaves_list
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PHASE_LABELS,
)
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.utils.metrics.reorganization import build_metric_specs

# ── Config ────────────────────────────────────────────────────────────
FC_METHOD = "msc"
LRG_CACHE = ROOT / "data" / "lrg_cache"

FIG_ROOT = ROOT / "data" / "figures" / "metric_concordance"
DATA_DIR = ROOT / "data" / "metric_concordance"
DATA_DIR.mkdir(parents=True, exist_ok=True)

PHASES = list(PHASE_LABELS)
BANDS = list(BRAIN_BANDS_NAMES)
BAND_TEX = [BRAIN_BAND_TEX_DICT[b] for b in BANDS]

# Display names for metrics
METRIC_DISPLAY = {
    "matrix_distance":      "Frobenius",
    "scaled_distance":      "Scaled (log)",
    "rank_distance":        "Rank (1-Spearman)",
    "quantile_rmse":        "Quantile RMSE",
    "permutation_robust":   "Permutation-robust",
    "tree_robinson_foulds": "Robinson-Foulds",
    "tree_cophenetic_corr": "Cophenetic corr.",
    "tree_baker_gamma":     "Baker gamma",
    "tree_fowlkes_mallows": "Fowlkes-Mallows",
}

# Metrics that are similarities (higher = more similar) — convert to distance
SIMILARITY_METRICS = {
    "tree_cophenetic_corr",
    "tree_baker_gamma",
    "tree_fowlkes_mallows",
}

# Discover patients with LRG cache
PATIENTS = sorted(
    p.name for p in LRG_CACHE.iterdir()
    if p.is_dir() and p.name.startswith("Pat_")
)

metric_specs = build_metric_specs()
metric_names = list(metric_specs.keys())
n_metrics = len(metric_names)


# ── Per-patient analysis ──────────────────────────────────────────────
def run_patient(patient):
    """Run full metric concordance analysis for one patient."""
    output_dir = FIG_ROOT / patient
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Load LRG results
    lrg_results = {}
    available_phases = set()
    for phase in PHASES:
        for band in BANDS:
            res = load_lrg_result(patient, phase, band, FC_METHOD, cache_root=LRG_CACHE)
            if res is not None:
                lrg_results[(phase, band)] = res
                available_phases.add(phase)

    available_phases = sorted(available_phases, key=PHASES.index)
    phase_pairs = list(combinations(available_phases, 2))

    if len(available_phases) < 2:
        print(f"  SKIP: only {len(available_phases)} phase(s) available")
        return

    print(f"  {len(lrg_results)} LRG results, "
          f"{len(available_phases)} phases ({', '.join(available_phases)}), "
          f"{len(phase_pairs)} pairs")

    # Step 2: Compute all metrics
    rows = []
    for band in BANDS:
        for p1, p2 in phase_pairs:
            res1 = lrg_results.get((p1, band))
            res2 = lrg_results.get((p2, band))
            if res1 is None or res2 is None:
                continue

            U1 = squareform(res1.ultrametric_matrix)
            U2 = squareform(res2.ultrametric_matrix)
            Z1 = res1.linkage_matrix
            Z2 = res2.linkage_matrix

            for mname, mspec in metric_specs.items():
                try:
                    val = mspec["fn"](U1, U2, Z1, Z2)
                except Exception:
                    val = np.nan

                if mname in SIMILARITY_METRICS and np.isfinite(val):
                    val = 1.0 - val

                rows.append({
                    "patient": patient,
                    "phase_A": p1,
                    "phase_B": p2,
                    "phase_pair": f"{p1}_vs_{p2}",
                    "band": band,
                    "metric": mname,
                    "value": val,
                })

    df = pd.DataFrame(rows)
    csv1_path = DATA_DIR / f"metric_concordance_{patient}.csv"
    df.to_csv(csv1_path, index=False)

    # Step 3: Inter-metric Spearman correlations
    pivot = df.pivot_table(index=["phase_pair", "band"], columns="metric", values="value")
    pivot = pivot[metric_names]

    corr_matrix = np.full((n_metrics, n_metrics), np.nan)
    for i in range(n_metrics):
        for j in range(n_metrics):
            vi = pivot.iloc[:, i].values
            vj = pivot.iloc[:, j].values
            valid = np.isfinite(vi) & np.isfinite(vj)
            if valid.sum() >= 3:
                corr_matrix[i, j], _ = spearmanr(vi[valid], vj[valid])

    corr_df = pd.DataFrame(corr_matrix, index=metric_names, columns=metric_names)
    csv2_path = DATA_DIR / f"metric_correlation_matrix_{patient}.csv"
    corr_df.to_csv(csv2_path)

    # Step 4: Hierarchical clustering of metrics
    dist_for_clust = 1.0 - np.abs(corr_matrix)
    dist_for_clust = (dist_for_clust + dist_for_clust.T) / 2
    np.fill_diagonal(dist_for_clust, 0)
    dist_for_clust = np.clip(dist_for_clust, 0, None)
    # Handle NaN: replace with max distance
    nan_mask = np.isnan(dist_for_clust)
    if nan_mask.any():
        dist_for_clust[nan_mask] = 1.0
    dist_condensed = squareform(dist_for_clust)

    Z_metrics = linkage(dist_condensed, method="average")

    for n_target in [3, 2, 4]:
        cluster_labels = fcluster(Z_metrics, n_target, criterion="maxclust")
        if len(np.unique(cluster_labels)) >= 2:
            break

    cluster_desc = {}
    for cid in np.unique(cluster_labels):
        member_idx = [i for i in range(n_metrics) if cluster_labels[i] == cid]
        members = [metric_names[i] for i in member_idx]
        n_matrix = sum(1 for m in members if not m.startswith("tree_"))
        n_tree = sum(1 for m in members if m.startswith("tree_") and
                     m not in ("tree_fowlkes_mallows",))
        n_part = sum(1 for m in members if m in ("tree_fowlkes_mallows",))
        total = len(members)
        if n_part > 0 and total <= 2:
            cluster_desc[cid] = "partition-sensitive"
        elif n_matrix > 0 and n_tree == 0 and n_part == 0:
            cluster_desc[cid] = "magnitude-sensitive"
        elif n_tree > 0 and n_matrix == 0:
            cluster_desc[cid] = "topology-sensitive"
        elif n_matrix >= n_tree:
            cluster_desc[cid] = "magnitude-topology"
        else:
            cluster_desc[cid] = "topology-magnitude"

    cluster_rows = []
    for i, mname in enumerate(metric_names):
        cid = int(cluster_labels[i])
        cluster_rows.append({
            "metric": mname,
            "cluster_id": cid,
            "cluster_label": cluster_desc[cid],
        })
    cluster_df = pd.DataFrame(cluster_rows)
    csv3_path = DATA_DIR / f"metric_clusters_{patient}.csv"
    cluster_df.to_csv(csv3_path, index=False)

    # ── Figure 1: Concordance heatmap ─────────────────────────────────
    leaf_order = leaves_list(Z_metrics)
    corr_reordered = corr_matrix[leaf_order][:, leaf_order]
    display_reordered = [METRIC_DISPLAY.get(metric_names[i], metric_names[i])
                         for i in leaf_order]
    clusters_reordered = cluster_labels[leaf_order]

    fig, ax_heat = plt.subplots(figsize=(8, 8))

    im = ax_heat.imshow(corr_reordered, cmap="RdBu_r", vmin=-1, vmax=1,
                         aspect="equal", interpolation="nearest")

    for i in range(n_metrics):
        for j in range(n_metrics):
            val = corr_reordered[i, j]
            if np.isfinite(val):
                color = "white" if abs(val) > 0.6 else "black"
                ax_heat.text(j, i, f"{val:.2f}", ha="center", va="center",
                            fontsize=7, color=color, fontweight="medium")

    ax_heat.set_xticks(range(n_metrics))
    ax_heat.set_xticklabels(display_reordered, rotation=45, ha="right", fontsize=8)
    ax_heat.set_yticks(range(n_metrics))
    ax_heat.set_yticklabels(display_reordered, fontsize=8)

    # Cluster brackets on the right margin
    unique_clusters_ordered = []
    cluster_ranges = {}
    for idx, cid in enumerate(clusters_reordered):
        if cid not in cluster_ranges:
            cluster_ranges[cid] = [idx, idx]
            unique_clusters_ordered.append(cid)
        else:
            cluster_ranges[cid][1] = idx

    bracket_colors = ["#e6194b", "#3cb44b", "#4363d8", "#f58231"]
    for ci, cid in enumerate(unique_clusters_ordered):
        r0, r1 = cluster_ranges[cid]
        color = bracket_colors[ci % len(bracket_colors)]
        ax_heat.add_patch(plt.Rectangle(
            (n_metrics - 0.5 + 0.15, r0 - 0.5), 0.3, r1 - r0 + 1,
            facecolor=color, edgecolor="none", clip_on=False, zorder=5,
        ))

    ax_heat.set_xlim(-0.5, n_metrics - 0.5)

    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.08, shrink=0.8)
    cbar.set_label("Spearman $r_S$", fontsize=10)

    # Draw dendrogram in heatmap coordinates (above the matrix)
    dn = dendrogram(Z_metrics, no_plot=True)
    y_top = -0.5
    dend_max_h = max(val for seg in dn["dcoord"] for val in seg)
    dend_scale = 2.5 / max(dend_max_h, 1e-10)

    for xcoords, ycoords, _ in zip(dn["icoord"], dn["dcoord"], dn["color_list"]):
        xs = [(x - 5) / 10 for x in xcoords]
        ys = [y_top - d * dend_scale for d in ycoords]
        ax_heat.plot(xs, ys, color="0.3", linewidth=1.2, clip_on=False, zorder=5)

    ax_heat.set_ylim(n_metrics - 0.5, y_top - 2.5 - 0.3)

    fig1_path = output_dir / f"fig_metric_concordance_{patient}.pdf"
    fig.savefig(fig1_path, bbox_inches="tight", dpi=300)
    plt.close(fig)

    # ── Figure 2: Metric values grid (3x3 strip plot) ────────────────
    pair_cmap = plt.cm.tab10
    all_pairs = sorted(df["phase_pair"].unique())
    pair_colors = {pp: pair_cmap(i / max(len(all_pairs) - 1, 1))
                   for i, pp in enumerate(all_pairs)}

    n_rows, n_cols = 3, 3
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(11, 9), sharey=False)

    for idx, mname in enumerate(metric_names):
        ax = axes.flat[idx]
        sub = df[df["metric"] == mname].copy()

        band_pos = {b: i for i, b in enumerate(BANDS)}
        sub["x"] = sub["band"].map(band_pos)

        for pp_str, grp in sub.groupby("phase_pair"):
            grp_sorted = grp.sort_values("x")
            ax.plot(grp_sorted["x"], grp_sorted["value"],
                    "-o", color=pair_colors[pp_str], markersize=4,
                    linewidth=1.2, alpha=0.85,
                    label=pp_str.replace("_vs_", " vs "))

        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels(BAND_TEX, fontsize=8)
        ax.set_title(METRIC_DISPLAY.get(mname, mname), fontsize=9, fontweight="bold")
        ax.tick_params(labelsize=7)
        ax.grid(axis="y", alpha=0.3)

    for idx in range(len(metric_names), n_rows * n_cols):
        axes.flat[idx].set_visible(False)

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=7,
               bbox_to_anchor=(0.5, 1.02), frameon=False)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig2_path = output_dir / f"fig_metric_values_grid_{patient}.pdf"
    fig.savefig(fig2_path, bbox_inches="tight", dpi=300)
    plt.close(fig)

    # Summary
    print(f"  -> {len(df)} values, {len(all_pairs)} pairs, "
          f"clusters: {dict(zip(cluster_df['metric'], cluster_df['cluster_label']))}")
    return df


# ── Main loop ─────────────────────────────────────────────────────────
print(f"Metric concordance analysis — {len(PATIENTS)} patients, "
      f"{n_metrics} metrics, {FC_METHOD}\n")

all_dfs = []
for patient in PATIENTS:
    print(f"\n{'─'*60}\n{patient}")
    patient_df = run_patient(patient)
    if patient_df is not None:
        all_dfs.append(patient_df)

# Save combined CSV
if all_dfs:
    combined = pd.concat(all_dfs, ignore_index=True)
    combined_path = DATA_DIR / "metric_concordance_all_patients.csv"
    combined.to_csv(combined_path, index=False)
    print(f"\nCombined CSV: {combined_path} ({len(combined)} rows)")

print(f"\nDone! Figures in {FIG_ROOT}, data in {DATA_DIR}")
