#!/usr/bin/env python3
"""WP1 — Updated figures for Section 5 of the report.

Task A: Patient-averaged 15x15 metric concordance matrix
Task B: Single-scale reorganization figures for Pat_05, Pat_07, Pat_08
        + n* summary table for all patients

Outputs in data/wp1_report_figures/

Run: python scripts/wp1/wp1_figures.py
"""
import sys
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec

# ── Global font scaling (matching gen_phase_reorg_figures.py) ──────────
plt.rcParams.update({
    "font.size":        14,
    "axes.titlesize":   15,
    "axes.labelsize":   14,
    "xtick.labelsize":  12,
    "ytick.labelsize":  12,
    "legend.fontsize":  12,
})

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import (
    linkage, fcluster, dendrogram, leaves_list,
)
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PHASE_LABELS,
)
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.utils.metrics.reorganization import build_metric_specs
from lrg_eegfc.visuals.lrg import (
    compute_partition_stability_index,
    _find_psi_optimal_partition,
)

# ── Config ────────────────────────────────────────────────────────────
FC_METHOD = "msc"
LRG_CACHE = ROOT / "data" / "lrg_cache"
WP0_CSV = ROOT / "data" / "wp0_metric_exploration" / "task2_full_results.csv"
METRIC_CONC_DIR = ROOT / "data" / "metric_concordance"
OUT = ROOT / "data" / "wp1_report_figures"
OUT.mkdir(parents=True, exist_ok=True)

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
TASK_B_PATIENTS = ["Pat_05", "Pat_07", "Pat_08"]

PHASES = list(PHASE_LABELS)
BANDS = list(BRAIN_BANDS_NAMES)
BAND_TEX = [BRAIN_BAND_TEX_DICT[b] for b in BANDS]
PHASE_SHORT = {
    "rest_pre": "rest_pre", "task_learn": "tLearn",
    "task_test": "tTest", "rest_post": "rest_post",
}

# 15 metrics in canonical order
ALL_METRICS = [
    "matrix_distance", "scaled_distance", "rank_distance", "quantile_rmse",
    "tree_robinson_foulds", "tree_cophenetic_corr", "tree_baker_gamma",
    "tree_fowlkes_mallows",
    "ari", "cluster_swap",
    "fc_frobenius", "fc_scaled_frobenius", "fc_rank_distance", "fc_mean_abs_diff",
]

METRIC_DISPLAY = {
    "matrix_distance":      "Frobenius",
    "scaled_distance":      "Scaled (log)",
    "rank_distance":        "Rank (1−ρ)",
    "quantile_rmse":        "Quantile RMSE",
    "tree_robinson_foulds": "Robinson–Foulds",
    "tree_cophenetic_corr": "Cophenetic corr.",
    "tree_baker_gamma":     "Baker γ",
    "tree_fowlkes_mallows": "Fowlkes–Mallows",
    "ari":                  "ARI",
    "cluster_swap":         "Cluster swap",
    "fc_frobenius":         "FC Frobenius",
    "fc_scaled_frobenius":  "FC scaled Frob.",
    "fc_rank_distance":     "FC rank (1−ρ)",
    "fc_mean_abs_diff":     "FC mean |Δ|",
}

# Family labels for grouping
METRIC_FAMILY = {}
for m in ["matrix_distance", "scaled_distance", "rank_distance",
          "quantile_rmse"]:
    METRIC_FAMILY[m] = "Ultrametric"
for m in ["tree_robinson_foulds", "tree_cophenetic_corr",
          "tree_baker_gamma", "tree_fowlkes_mallows"]:
    METRIC_FAMILY[m] = "Tree"
for m in ["ari", "cluster_swap"]:
    METRIC_FAMILY[m] = "Partition"
for m in ["fc_frobenius", "fc_scaled_frobenius",
          "fc_rank_distance", "fc_mean_abs_diff"]:
    METRIC_FAMILY[m] = "Raw FC"

FAMILY_COLORS = {
    "Ultrametric": "#1f77b4",
    "Tree": "#2ca02c",
    "Partition": "#ff7f0e",
    "Raw FC": "#d62728",
}

# Similarity metrics (higher = more similar) — need sign flip for distance convention
SIMILARITY_METRICS = {
    "tree_cophenetic_corr", "tree_baker_gamma", "tree_fowlkes_mallows",
    "ari", "cluster_swap",
}

# Metrics for reorganization matrices (Task B Fig 1)
SELECTED_METRICS = ["scaled_distance", "tree_robinson_foulds", "tree_fowlkes_mallows"]
METRIC_LABELS_B = {
    "scaled_distance": "Scaled distance",
    "tree_robinson_foulds": "Robinson\u2013Foulds",
    "tree_fowlkes_mallows": "Fowlkes\u2013Mallows",
}

# Bands for dendrograms (Task B Fig 3) — per-patient selection:
# reorganizes / stable / irregular
DENDRO_BANDS_PER_PATIENT = {
    "Pat_05": ["beta", "alpha", "high_gamma"],       # reorganizes (7→9→9→5) / stable (2→2→2→2) / irregular (7→2→2→2)
    "Pat_07": ["theta", "delta", "alpha"],            # reorganizes (8→8→7→2) / stable (5→6→5→4) / irregular (4→11→5→2)
    "Pat_08": ["delta", "beta", "high_gamma"],        # mild variation (6→5→4→6) / task unfolds (2→4→2→2) / stable (2→2→2→2)
}
DENDRO_BANDS_DEFAULT = ["theta", "beta", "high_gamma"]
COMM_COLORS = plt.cm.tab10.colors


# ── Helpers ───────────────────────────────────────────────────────────
def save_fig(fig, path, *, what, proves, how_to_read, dpi=300):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"  Saved: {path}")
    md = path.with_suffix(".md")
    md.write_text(
        f"# {path.stem}\n\n"
        f"## What the figure shows\n{what}\n\n"
        f"## What it proves\n{proves}\n\n"
        f"## How to read it\n{how_to_read}\n"
    )


def align_labels_hungarian(ref_labels, target_labels):
    ref = np.asarray(ref_labels)
    tgt = np.asarray(target_labels)
    u_ref, u_tgt = np.unique(ref), np.unique(tgt)
    dim = max(len(u_ref), len(u_tgt))
    cost = np.zeros((dim, dim))
    for i, r in enumerate(u_ref):
        for j, t in enumerate(u_tgt):
            cost[i, j] = -np.sum((ref == r) & (tgt == t))
    ri, ci = linear_sum_assignment(cost)
    mapping = {}
    for r, c in zip(ri, ci):
        if r < len(u_ref) and c < len(u_tgt):
            mapping[u_tgt[c]] = u_ref[r]
    used = set(mapping.values())
    nxt = max(max(used) + 1, max(u_ref) + 1) if used else 1
    for t in u_tgt:
        if t not in mapping:
            mapping[t] = nxt
            nxt += 1
    return np.array([mapping[l] for l in tgt])


def get_optimal_partition(linkage_matrix):
    psi, n_comm = compute_partition_stability_index(linkage_matrix)
    if len(psi) == 0:
        return 2
    return _find_psi_optimal_partition(psi, n_comm)


def build_node_to_comm(Z, labels):
    n_nodes = len(labels)
    node_to_comm = {}
    for i in range(n_nodes):
        node_to_comm[i] = int(labels[i])
    for i in range(len(Z)):
        left_id = int(Z[i, 0])
        right_id = int(Z[i, 1])
        lc = node_to_comm.get(left_id, -1)
        rc = node_to_comm.get(right_id, -1)
        node_to_comm[n_nodes + i] = lc if lc == rc else -1
    return node_to_comm


def make_link_color_func(node_to_comm):
    def fn(k):
        c = node_to_comm.get(k, -1)
        if c >= 0:
            return mcolors.to_hex(COMM_COLORS[c % len(COMM_COLORS)])
        return "#aaaaaa"
    return fn


def load_patient_lrg(patient):
    results = {}
    for phase in PHASES:
        for band in BANDS:
            res = load_lrg_result(patient, phase, band, FC_METHOD,
                                  cache_root=LRG_CACHE)
            if res is not None:
                results[(phase, band)] = res
    return results


# ======================================================================
# TASK A: Patient-averaged 15×15 metric concordance matrix
# ======================================================================
def task_a():
    print("\n" + "=" * 70)
    print("TASK A: Patient-averaged 15×15 metric concordance")
    print("=" * 70)

    df = pd.read_csv(WP0_CSV)

    # Convert similarity metrics to distance convention for consistent correlation
    def to_dist(row):
        if row["metric_name"] in SIMILARITY_METRICS:
            return -row["value"]
        return row["value"]

    df["dist_value"] = df.apply(to_dist, axis=1)

    # Per-patient Spearman correlation matrices
    patient_corrs = []
    for pat in PATIENTS:
        pdf = df[df["patient"] == pat]
        pivot = pdf.pivot_table(
            index=["band", "phase_a", "phase_b"],
            columns="metric_name",
            values="dist_value",
        )
        metrics_present = [m for m in ALL_METRICS if m in pivot.columns]
        pivot = pivot[metrics_present]

        n = len(metrics_present)
        corr = np.full((n, n), np.nan)
        for i in range(n):
            for j in range(n):
                valid = pivot.iloc[:, i].notna() & pivot.iloc[:, j].notna()
                if valid.sum() >= 3:
                    corr[i, j] = spearmanr(
                        pivot.iloc[:, i][valid], pivot.iloc[:, j][valid]
                    )[0]

        patient_corrs.append(corr)
        print(f"  {pat}: {valid.sum()} data points, {n} metrics")

    # Average across patients
    avg_corr = np.nanmean(patient_corrs, axis=0)
    n = avg_corr.shape[0]
    metrics_present = [m for m in ALL_METRICS if m in pivot.columns]

    # Hierarchical clustering for reordering
    dist = 1.0 - np.abs(avg_corr)
    dist = (dist + dist.T) / 2
    np.fill_diagonal(dist, 0)
    dist = np.clip(np.nan_to_num(dist, nan=1.0), 0, None)
    dist_condensed = squareform(dist)
    Z_metrics = linkage(dist_condensed, method="average")
    leaf_order = leaves_list(Z_metrics)

    corr_reordered = avg_corr[leaf_order][:, leaf_order]
    labels_reordered = [metrics_present[i] for i in leaf_order]
    display_reordered = [METRIC_DISPLAY.get(m, m) for m in labels_reordered]
    family_reordered = [METRIC_FAMILY.get(m, "?") for m in labels_reordered]

    # ── Figure ────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 10))

    im = ax.imshow(corr_reordered, cmap="RdBu_r", vmin=-1, vmax=1,
                   aspect="equal", interpolation="nearest")

    # Annotate cells
    for i in range(n):
        for j in range(n):
            val = corr_reordered[i, j]
            if np.isfinite(val):
                color = "white" if abs(val) > 0.6 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=5.5, color=color)

    # Tick labels colored by family
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    xlabels = ax.set_xticklabels(display_reordered, rotation=45, ha="right", fontsize=8)
    ylabels = ax.set_yticklabels(display_reordered, fontsize=8)
    for i, (xl, yl) in enumerate(zip(xlabels, ylabels)):
        color = FAMILY_COLORS.get(family_reordered[i], "black")
        xl.set_color(color)
        yl.set_color(color)

    # Dendrogram above matrix
    dn = dendrogram(Z_metrics, no_plot=True)
    y_top = -0.5
    dend_max_h = max(val for seg in dn["dcoord"] for val in seg)
    dend_scale = 2.5 / max(dend_max_h, 1e-10)
    for xcoords, ycoords in zip(dn["icoord"], dn["dcoord"]):
        xs = [(x - 5) / 10 for x in xcoords]
        ys = [y_top - d * dend_scale for d in ycoords]
        ax.plot(xs, ys, color="0.3", linewidth=1.0, clip_on=False, zorder=5)
    ax.set_ylim(n - 0.5, y_top - 2.5 - 0.3)

    # Colorbar — tight to the axes
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="Spearman $r_S$")

    # Family legend
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=c, label=f) for f, c in FAMILY_COLORS.items()]
    ax.legend(handles=handles, loc="lower left", fontsize=7,
              title="Metric family", title_fontsize=8)

    ax.set_title(f"Inter-metric concordance (patient-averaged, N={len(PATIENTS)})")

    save_fig(fig, OUT / "metric_concordance_15x15_avg.pdf",
        what=f"15×15 Spearman rank correlation matrix between all metrics (11 LRG + 4 raw FC), "
             f"averaged across {len(PATIENTS)} patients. Rows/columns reordered by hierarchical "
             f"clustering. Dendrogram shown above. Metric labels colored by family: "
             f"blue=Ultrametric, green=Tree, orange=Partition, red=Raw FC. "
             f"Black separator lines between families. "
             f"Each patient contributes 36 data points (6 bands × 6 phase pairs).",
        proves="Inter-metric redundancy structure at the population level. "
               "Highly correlated blocks (ρ > 0.9) are redundant; weakly correlated metrics "
               "carry independent information. The family clustering shows whether LRG metrics "
               "and raw FC metrics capture the same or different aspects of reorganization.",
        how_to_read="Red = positive correlation (redundant), Blue = negative correlation. "
                    "Family-colored labels help identify which metric types cluster together. "
                    "Separator lines show family boundaries after clustering reordering.",
    )

    # Save companion markdown with detailed findings
    md_lines = [
        "# Patient-averaged 15×15 metric concordance",
        f"",
        f"Averaged across {len(PATIENTS)} patients, 36 data points each.",
        "",
        "## Highly redundant pairs (avg |ρ| > 0.9):",
        "",
    ]
    for i in range(n):
        for j in range(i + 1, n):
            if abs(avg_corr[metrics_present.index(labels_reordered[i]),
                            metrics_present.index(labels_reordered[j])]) > 0.9:
                val = avg_corr[metrics_present.index(labels_reordered[i]),
                               metrics_present.index(labels_reordered[j])]
                md_lines.append(
                    f"- `{labels_reordered[i]}` ↔ `{labels_reordered[j]}`: ρ = {val:.3f}"
                )

    md_lines += ["", "## Weakly correlated (max |ρ| < 0.5):", ""]
    for i in range(n):
        off_diag = [abs(avg_corr[metrics_present.index(labels_reordered[i]),
                                  metrics_present.index(labels_reordered[j])])
                    for j in range(n) if j != i]
        if off_diag and max(off_diag) < 0.5:
            md_lines.append(f"- `{labels_reordered[i]}`: max |ρ| = {max(off_diag):.3f}")

    md_lines += ["", "## Cross-family correlations:", ""]
    for fam1 in ["Ultrametric", "Tree", "Partition", "Raw FC"]:
        for fam2 in ["Ultrametric", "Tree", "Partition", "Raw FC"]:
            if fam1 >= fam2:
                continue
            vals = []
            for i in range(n):
                for j in range(n):
                    mi, mj = labels_reordered[i], labels_reordered[j]
                    fi, fj = METRIC_FAMILY.get(mi), METRIC_FAMILY.get(mj)
                    if fi == fam1 and fj == fam2:
                        v = avg_corr[metrics_present.index(mi), metrics_present.index(mj)]
                        if np.isfinite(v):
                            vals.append(v)
            if vals:
                md_lines.append(
                    f"- {fam1} ↔ {fam2}: mean ρ = {np.mean(vals):.3f} "
                    f"(range [{min(vals):.3f}, {max(vals):.3f}])"
                )

    (OUT / "metric_concordance_15x15_avg_detailed.md").write_text("\n".join(md_lines))


# ======================================================================
# TASK B: Single-scale reorganization figures for multiple patients
# ======================================================================
def task_b():
    print("\n" + "=" * 70)
    print("TASK B: Single-scale reorganization figures")
    print("=" * 70)

    # First compute metric concordance CSVs for patients that don't have them yet
    metric_specs = build_metric_specs()
    metric_names = list(metric_specs.keys())

    for patient in TASK_B_PATIENTS:
        csv_path = METRIC_CONC_DIR / f"metric_concordance_{patient}.csv"
        if not csv_path.exists():
            print(f"\n  Computing metric concordance for {patient}...")
            compute_metric_concordance(patient, metric_specs, metric_names)
        else:
            print(f"  {patient} concordance CSV exists: {csv_path}")

    # Generate figures for each patient
    for patient in TASK_B_PATIENTS:
        print(f"\n── {patient} ──")
        lrg = load_patient_lrg(patient)
        print(f"  {len(lrg)} LRG results loaded")

        csv_path = METRIC_CONC_DIR / f"metric_concordance_{patient}.csv"
        df_metrics = pd.read_csv(csv_path)

        fig1_reorg_matrices(patient, lrg, df_metrics)
        fig2_ultrametric_sidebyside(patient, lrg)
        fig3_dendrograms(patient, lrg)

    # n* summary table for ALL patients
    fig4_nstar_table()


def compute_metric_concordance(patient, metric_specs, metric_names):
    """Compute metric concordance CSV (same as gen_metric_concordance.py)."""
    METRIC_CONC_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for band in BANDS:
        for p1, p2 in combinations(PHASES, 2):
            res1 = load_lrg_result(patient, p1, band, FC_METHOD, cache_root=LRG_CACHE)
            res2 = load_lrg_result(patient, p2, band, FC_METHOD, cache_root=LRG_CACHE)
            if res1 is None or res2 is None:
                continue
            U1 = squareform(res1.ultrametric_matrix)
            U2 = squareform(res2.ultrametric_matrix)
            Z1, Z2 = res1.linkage_matrix, res2.linkage_matrix

            for mname, mspec in metric_specs.items():
                try:
                    val = mspec["fn"](U1, U2, Z1, Z2)
                except Exception:
                    val = np.nan
                if mname in {"tree_cophenetic_corr", "tree_baker_gamma",
                             "tree_fowlkes_mallows"} and np.isfinite(val):
                    val = 1.0 - val
                rows.append({
                    "patient": patient, "phase_A": p1, "phase_B": p2,
                    "phase_pair": f"{p1}_vs_{p2}", "band": band,
                    "metric": mname, "value": val,
                })
    df = pd.DataFrame(rows)
    df.to_csv(METRIC_CONC_DIR / f"metric_concordance_{patient}.csv", index=False)
    print(f"    Saved {len(df)} rows")


def fig1_reorg_matrices(patient, lrg, df_metrics):
    """3×6 grid of 4×4 reorganization matrices (3 metrics × 6 bands).

    Per-row normalization with colorbars. Diagonal masked in gray.
    """
    fig = plt.figure(figsize=(16, 8))
    gs_outer = GridSpec(3, 1, figure=fig, hspace=0.40)

    for row_idx, metric in enumerate(SELECTED_METRICS):
        sub = df_metrics[df_metrics["metric"] == metric]

        # Row-wide normalization: use full range across all bands for this metric
        all_vals = sub["value"].dropna()
        vmin_row = 0  # all metrics have 0 as self-distance baseline
        vmax_row = all_vals.max() if len(all_vals) > 0 else 1.0

        gs_row = gs_outer[row_idx].subgridspec(
            1, 7, width_ratios=[1, 1, 1, 1, 1, 1, 0.06], wspace=0.12
        )

        last_im = None
        for col_idx, band in enumerate(BANDS):
            ax = fig.add_subplot(gs_row[0, col_idx])

            mat = np.full((4, 4), np.nan)
            bsub = sub[sub["band"] == band]
            for _, r in bsub.iterrows():
                i = PHASES.index(r["phase_A"])
                j = PHASES.index(r["phase_B"])
                mat[i, j] = r["value"]
                mat[j, i] = r["value"]

            # Mask diagonal for display
            mat_display = np.copy(mat)
            np.fill_diagonal(mat_display, np.nan)

            # Plot with masked diagonal (NaN → gray background)
            ax.set_facecolor("#e0e0e0")
            cmap = plt.cm.plasma.copy()
            cmap.set_bad("#e0e0e0")
            im = ax.imshow(mat_display, cmap=cmap, vmin=vmin_row, vmax=vmax_row,
                           aspect="equal", interpolation="nearest")
            last_im = im

            # Annotate off-diagonal with values
            for i in range(4):
                for j in range(4):
                    if i == j:
                        continue
                    val = mat[i, j]
                    if np.isnan(val):
                        continue
                    norm_val = (val - vmin_row) / max(vmax_row - vmin_row, 1e-10)
                    color = "white" if norm_val > 0.55 else "black"
                    ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                            fontsize=9, fontweight="bold", color=color)

            if row_idx == 2:
                ax.set_xticks(range(4))
                ax.set_xticklabels([PHASE_SHORT[p] for p in PHASES],
                                   rotation=35, ha="right")
            else:
                ax.set_xticks([])
            if col_idx == 0:
                ax.set_yticks(range(4))
                ax.set_yticklabels([PHASE_SHORT[p] for p in PHASES])
                ax.set_ylabel(METRIC_LABELS_B[metric], fontweight="bold")
            else:
                ax.set_yticks([])
            if row_idx == 0:
                ax.set_title(BAND_TEX[col_idx])

        # Colorbar for this row — matched to row height
        cax = fig.add_subplot(gs_row[0, 6])
        fig.colorbar(last_im, cax=cax)

    save_fig(fig, OUT / f"fig_reorg_matrices_{patient}.pdf",
        what=f"3×6 grid of 4×4 phase-reorganization distance matrices for {patient}. "
             f"Rows = metrics (Scaled distance, Robinson–Foulds, Fowlkes–Mallows). "
             f"Columns = frequency bands (δ through γh). Each cell = distance between "
             f"two phases (YlOrRd colormap, brighter = more reorganization).",
        proves=f"Extends Section 5.3 to {patient}, showing whether the reorganization "
               f"pattern is consistent with Pat_02.",
        how_to_read="Brighter off-diagonal cells = more reorganization between those phases. "
                    "The task_learn-task_test cell should be darkest (most similar). "
                    "Compare structure across bands and with Pat_02.",
    )


def fig2_ultrametric_sidebyside(patient, lrg):
    """Ultrametric D(τ') matrices for ALL 6 bands × 4 phases.

    Shows all bands so that trace (rest_post ≈ task) vs reset (rest_post ≈ rest_pre)
    is visually apparent across the full frequency spectrum.
    """
    # Use rest_pre of first available band for leaf ordering
    leaf_order = None
    for band in BANDS:
        ref_key = ("rest_pre", band)
        if ref_key in lrg:
            leaf_order = leaves_list(lrg[ref_key].linkage_matrix)
            break
    if leaf_order is None:
        print(f"  SKIP ultrametric sidebyside: no rest_pre for {patient}")
        return

    n_bands = len(BANDS)
    fig = plt.figure(figsize=(13, 2.8 * n_bands))
    gs = GridSpec(n_bands, 5, figure=fig, width_ratios=[1, 1, 1, 1, 0.05],
                  wspace=0.08, hspace=0.20)

    for row_idx, band in enumerate(BANDS):
        matrices = []
        for phase in PHASES:
            res = lrg.get((phase, band))
            if res is None:
                matrices.append(None)
                continue
            U = squareform(res.ultrametric_matrix)
            # Reorder to match leaf_order (handle size mismatch gracefully)
            if U.shape[0] == len(leaf_order):
                U_ordered = U[np.ix_(leaf_order, leaf_order)]
            else:
                U_ordered = U
            matrices.append(U_ordered)

        valid = [m for m in matrices if m is not None]
        if not valid:
            continue
        all_vals = np.concatenate([m[m > 0].ravel() for m in valid])
        if len(all_vals) == 0:
            continue
        vmin_log = np.log10(all_vals.min())
        vmax_log = np.log10(all_vals.max())

        last_im = None
        for col_idx, (phase, mat) in enumerate(zip(PHASES, matrices)):
            ax = fig.add_subplot(gs[row_idx, col_idx])
            if mat is None:
                ax.set_visible(False)
                continue
            mat_log = np.log10(np.where(mat > 0, mat, np.nan))
            im = ax.imshow(mat_log, cmap="magma", vmin=vmin_log, vmax=vmax_log,
                           aspect="equal", interpolation="nearest")
            last_im = im
            ax.set_xticks([])
            ax.set_yticks([])
            if row_idx == 0:
                ax.set_title(phase)
            if col_idx == 0:
                ax.set_ylabel(BAND_TEX[row_idx] + " band", fontweight="bold")

        if last_im is not None:
            cax = fig.add_subplot(gs[row_idx, 4])
            cbar = fig.colorbar(last_im, cax=cax)
            if row_idx == n_bands - 1:
                cbar.set_label(r"$\log_{10}\, D$")

    save_fig(fig, OUT / f"fig_ultrametric_sidebyside_{patient}.pdf",
        what=f"Ultrametric distance matrices D(τ') for {patient}, 4 phases × 6 bands. "
             f"Log₁₀ color scale (magma). Node order fixed by rest_pre leaf ordering. "
             f"Each row = one frequency band, each column = one cognitive phase.",
        proves=f"Visual comparison of hierarchical structure across ALL bands for {patient}. "
               f"Bands where rest_post resembles the task phases (task_learn/task_test) show TRACE. "
               f"Bands where rest_post resembles rest_pre show RESET. "
               f"Alpha should show trace (rest_post ≈ task), beta/high_gamma should show reset.",
        how_to_read="Each panel = one ultrametric matrix. Compare the PATTERN (not just brightness) "
                    "across columns. If rest_post (col 4) looks like task_learn/task_test (cols 2-3), "
                    "the task left a trace. If rest_post looks like rest_pre (col 1), the hierarchy reset. "
                    "Scroll down the rows to see which bands trace vs reset.",
    )


def _get_nstar_and_cut(linkage_matrix, n_nodes):
    """Get n* via PSI peak detection, capped to reasonable range.

    Uses _find_psi_optimal_partition (same as original Pat_02 script),
    with a cap at n_nodes/2 to avoid degenerate partitions.
    """
    n_star = get_optimal_partition(linkage_matrix)
    # Cap: if PSI gives a degenerate result (n* > N/3), fall back to
    # the largest-gap partition (n*=2 equivalent)
    if n_star > n_nodes // 3:
        n_star = 2
    labels = fcluster(linkage_matrix, n_star, criterion="maxclust")
    # Compute cut height: midpoint between the two merge heights around the cut
    merge_heights = linkage_matrix[:, 2]
    cut_idx = n_nodes - n_star
    if 0 < cut_idx < len(merge_heights):
        cut_height = (merge_heights[cut_idx - 1] + merge_heights[cut_idx]) / 2
    else:
        cut_height = merge_heights[0] * 0.5
    return n_star, labels, cut_height


def fig3_dendrograms(patient, lrg):
    """Phase-resolved dendrograms for 3 patient-specific bands.

    Band selection: one that reorganizes, one stable, one irregular.
    Uses PSI peak detection with cap for n*.
    """
    dendro_bands = DENDRO_BANDS_PER_PATIENT.get(patient, DENDRO_BANDS_DEFAULT)
    dendro_band_labels = [BRAIN_BAND_TEX_DICT[b] + " band" for b in dendro_bands]

    fig, axes = plt.subplots(len(dendro_bands), 4,
                              figsize=(16, 4 * len(dendro_bands)))

    for row_idx, band in enumerate(dendro_bands):
        partitions = {}
        n_stars = {}
        cut_heights = {}
        for phase in PHASES:
            res = lrg.get((phase, band))
            if res is None:
                continue
            n_star, labels, cut_h = _get_nstar_and_cut(
                res.linkage_matrix, res.n_nodes
            )
            partitions[phase] = labels
            n_stars[phase] = n_star
            cut_heights[phase] = cut_h

        # Align community labels to rest_pre reference
        if "rest_pre" in partitions:
            ref = partitions["rest_pre"]
            for phase in PHASES[1:]:
                if phase in partitions:
                    partitions[phase] = align_labels_hungarian(ref, partitions[phase])

        for col_idx, phase in enumerate(PHASES):
            ax = axes[row_idx, col_idx]
            res = lrg.get((phase, band))
            if res is None or phase not in partitions:
                ax.set_visible(False)
                continue

            labels = partitions[phase]
            n_star = n_stars[phase]
            cut_height = cut_heights[phase]
            Z = res.linkage_matrix
            n_nodes = len(labels)
            merge_heights = Z[:, 2]
            tmin = merge_heights[0] * 0.8
            tmax = merge_heights[-1] * 1.05

            node_to_comm = build_node_to_comm(Z, labels)
            link_color_fn = make_link_color_func(node_to_comm)

            dendrogram(
                Z, ax=ax, orientation="right",
                truncate_mode="lastp", p=min(40, n_nodes - 1),
                link_color_func=link_color_fn,
                above_threshold_color="#aaaaaa",
                no_labels=True, color_threshold=0,
            )
            ax.set_xscale("log")
            ax.set_xlim(tmin, tmax)
            ax.axvline(cut_height, color="k", linestyle="--", linewidth=1.0, alpha=0.7)

            if row_idx == 0:
                ax.set_title(f"{phase}\n$n^*={n_star}$")
            else:
                ax.set_title(f"$n^*={n_star}$")
            if col_idx == 0:
                ax.set_ylabel(dendro_band_labels[row_idx], fontweight="bold")
            for spine in ["top", "right"]:
                ax.spines[spine].set_visible(False)

    fig.tight_layout()
    save_fig(fig, OUT / f"fig_dendrograms_phase_comparison_{patient}.pdf",
        what=f"Phase-resolved dendrograms for {patient}: 3 bands (θ, β, γh) × 4 phases. "
             f"Horizontal dendrograms with log x-axis (merge distance). "
             f"Colors = community assignment at Ψ-optimal n*. "
             f"Dashed line = cut height for n* clusters.",
        proves=f"Visual comparison of hierarchical community structure across phases "
               f"for {patient}. n* values shown per panel.",
        how_to_read="Colors show communities (aligned to rest_pre via Hungarian matching). "
                    "If colors are consistent across phases, community structure is stable. "
                    "n* indicates the optimal number of communities at each phase.",
    )


def fig4_nstar_table():
    """n* summary table for all patients."""
    print("\n── n* summary table ──")
    all_patients = ["Pat_02"] + TASK_B_PATIENTS

    rows = []
    for patient in all_patients:
        lrg = load_patient_lrg(patient)
        for phase in PHASES:
            for band in BANDS:
                res = lrg.get((phase, band))
                if res is None:
                    continue
                n_star, _, cut_h = _get_nstar_and_cut(
                    res.linkage_matrix, res.n_nodes
                )
                rows.append({
                    "patient": patient, "phase": phase, "band": band,
                    "n_star": n_star, "n_nodes": res.n_nodes,
                    "threshold": cut_h,
                })

    df = pd.DataFrame(rows)

    # Format as markdown table
    lines = [
        "# Ψ-optimal community count n* across patients",
        "",
        f"Patients: {', '.join(all_patients)}",
        "",
    ]

    for patient in all_patients:
        pdf = df[df["patient"] == patient]
        lines.append(f"## {patient} (N = {pdf['n_nodes'].iloc[0] if len(pdf) > 0 else '?'} nodes)")
        lines.append("")
        lines.append("| Band | " + " | ".join(PHASES) + " |")
        lines.append("|------|" + "|".join(["------"] * 4) + "|")

        for band in BANDS:
            vals = []
            for phase in PHASES:
                sub = pdf[(pdf["band"] == band) & (pdf["phase"] == phase)]
                if not sub.empty:
                    vals.append(str(sub["n_star"].iloc[0]))
                else:
                    vals.append("—")
            lines.append(f"| {BRAIN_BAND_TEX_DICT[band]} | " + " | ".join(vals) + " |")
        lines.append("")

    # Cross-patient summary
    lines += [
        "## Cross-patient summary: mean n* per (band, phase)",
        "",
        "| Band | " + " | ".join(PHASES) + " |",
        "|------|" + "|".join(["------"] * 4) + "|",
    ]
    for band in BANDS:
        vals = []
        for phase in PHASES:
            sub = df[(df["band"] == band) & (df["phase"] == phase)]
            if not sub.empty:
                vals.append(f"{sub['n_star'].mean():.1f} ± {sub['n_star'].std():.1f}")
            else:
                vals.append("—")
        lines.append(f"| {BRAIN_BAND_TEX_DICT[band]} | " + " | ".join(vals) + " |")

    out_path = OUT / "nstar_summary_table.md"
    out_path.write_text("\n".join(lines))
    print(f"  Saved: {out_path}")


# ======================================================================
# TASK C: Task-task similarity investigation
# ======================================================================
def task_c():
    print("\n" + "=" * 70)
    print("TASK C: Task-task similarity tally")
    print("=" * 70)

    FIXES_DIR = ROOT / "data" / "wp1_fixes"
    FIXES_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(WP0_CSV)
    all_pats = ["Pat_02", "Pat_05", "Pat_07", "Pat_08"]  # exclude Pat_03

    # For similarity metrics, lower value = more similar (distance convention)
    # For similarity metrics stored as raw similarity, HIGHER = more similar
    # In task2_full_results.csv: all LRG metrics from build_metric_specs are
    # already in their native convention. Similarity metrics need to be handled.
    # tree_cophenetic_corr, tree_baker_gamma, tree_fowlkes_mallows, ari, cluster_swap
    # are similarities (higher = more similar). For finding minimum distance,
    # we want the MAXIMUM of these.

    REORG_METRICS = ["scaled_distance", "tree_robinson_foulds", "tree_fowlkes_mallows"]
    # All 11 LRG metrics
    LRG_METRICS = [m for m in ALL_METRICS if not m.startswith("fc_")]

    rows = []
    for pat in all_pats:
        for band in BANDS:
            for metric in LRG_METRICS:
                is_sim = metric in SIMILARITY_METRICS
                sub = df[(df["patient"] == pat) & (df["band"] == band) &
                         (df["metric_name"] == metric)]
                if sub.empty:
                    continue

                # Get all 6 pair values
                pair_vals = {}
                for _, r in sub.iterrows():
                    pair = (r["phase_a"], r["phase_b"])
                    pair_vals[pair] = r["value"]

                if len(pair_vals) < 6:
                    continue

                # taskL-taskT value
                tt_key = ("task_learn", "task_test")
                if tt_key not in pair_vals:
                    continue
                tt_val = pair_vals[tt_key]

                # For distance metrics: minimum = most similar
                # For similarity metrics: maximum = most similar
                if is_sim:
                    sorted_pairs = sorted(pair_vals.items(), key=lambda x: -x[1])
                else:
                    sorted_pairs = sorted(pair_vals.items(), key=lambda x: x[1])

                rank = [p for p, v in sorted_pairs].index(tt_key) + 1
                is_min = rank == 1

                rows.append({
                    "patient": pat,
                    "band": band,
                    "metric": metric,
                    "taskL_taskT_value": tt_val,
                    "is_most_similar": is_min,
                    "rank_among_6": rank,
                })

    result_df = pd.DataFrame(rows)
    csv_path = FIXES_DIR / "task_task_similarity_full_table.csv"
    result_df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # Tally per metric
    total_cells = len(all_pats) * len(BANDS)
    md = [
        "# Task-task similarity tally",
        "",
        f"For each (patient, band, metric): is task_learn–task_test the most similar pair?",
        f"Patients: {', '.join(all_pats)} (N={len(all_pats)}), 6 bands, total cells = {total_cells}",
        "",
        "## Tally per metric",
        "",
        "| Metric | Most similar | Fraction | Mean rank |",
        "|--------|-------------|----------|-----------|",
    ]

    for metric in LRG_METRICS:
        sub = result_df[result_df["metric"] == metric]
        n_min = sub["is_most_similar"].sum()
        frac = n_min / len(sub) if len(sub) > 0 else 0
        mean_rank = sub["rank_among_6"].mean() if len(sub) > 0 else 0
        bold = "**" if frac >= 0.5 else ""
        md.append(f"| {bold}`{metric}`{bold} | {bold}{n_min}/{len(sub)}{bold} | "
                  f"{bold}{frac:.2f}{bold} | {mean_rank:.1f} |")

    md += [
        "",
        "Chance level: 1/6 = 0.167",
        "",
        "## Focus: three reorg-matrix metrics",
        "",
        "| Patient | Band | scaled_dist rank | RF rank | FM rank | All agree? |",
        "|---------|------|-----------------|---------|---------|------------|",
    ]

    disagree_cells = []
    for pat in all_pats:
        for band in BANDS:
            ranks = {}
            for metric in REORG_METRICS:
                sub = result_df[(result_df["patient"] == pat) &
                                (result_df["band"] == band) &
                                (result_df["metric"] == metric)]
                if not sub.empty:
                    ranks[metric] = sub["rank_among_6"].iloc[0]

            if len(ranks) == 3:
                all_1 = all(r == 1 for r in ranks.values())
                any_1 = any(r == 1 for r in ranks.values())
                sd_r = int(ranks["scaled_distance"])
                rf_r = int(ranks["tree_robinson_foulds"])
                fm_r = int(ranks["tree_fowlkes_mallows"])
                agree_str = "YES" if all_1 else ("partial" if any_1 else "NO")
                bold = "**" if not all_1 else ""
                md.append(
                    f"| {bold}{pat}{bold} | {bold}{band}{bold} | "
                    f"{sd_r} | {rf_r} | {fm_r} | {bold}{agree_str}{bold} |"
                )
                if not all_1:
                    disagree_cells.append((pat, band, sd_r, rf_r, fm_r))

    md += [
        "",
        f"## Disagreement cells: {len(disagree_cells)}/{total_cells}",
        "",
    ]
    if disagree_cells:
        for pat, band, sd, rf, fm in disagree_cells:
            md.append(f"- {pat} {band}: scaled_dist rank={sd}, RF rank={rf}, FM rank={fm}")
    else:
        md.append("None — all three metrics agree on every cell.")

    md += [
        "",
        "## Per-patient summary",
        "",
    ]
    for pat in all_pats:
        sub = result_df[(result_df["patient"] == pat) &
                        (result_df["metric"].isin(REORG_METRICS))]
        n_min = sub["is_most_similar"].sum()
        total = len(sub)
        md.append(f"- **{pat}**: taskL-taskT is most similar in {n_min}/{total} "
                  f"(metric × band) cells ({n_min/total:.0%})")

    md_path = FIXES_DIR / "task_task_similarity_tally.md"
    md_path.write_text("\n".join(md))
    print(f"  Saved: {md_path}")


# ── Main ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    task_a()
    task_b()
    task_c()
    print(f"\nAll done. Outputs in {OUT}")
