#!/usr/bin/env python3
"""WP0 Task 3 — Metric behavior analysis (claims A-D, redundancy, normalization).

Reads task2_full_results.csv and produces figures + companion .md files
for each sub-analysis (3a through 3f).

Produces (under data/wp0_metric_exploration/task3_analysis/):
  metric_correlation_matrix.pdf + .md
  claim_A_band_dependence/  (2 figures + .md)
  claim_B_task_signature/   (1 figure + .md)
  claim_C_task_imprint/     (1 figure + .md)
  normalization_check.md
  cross_patient_consistency/ (1 figure + .md, signal_to_noise.md)

Run: python scripts/wp0/task3_analysis.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "wp0"))
sys.path.insert(0, str(ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from _common import (
    ALL_PAIRS,
    BANDS,
    PATIENTS,
    PHASES,
    OUT_ROOT,
    _METRIC_META,
    classify_pair,
    save_fig,
    to_distance,
)

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT

CSV_IN = OUT_ROOT / "task2_full_results.csv"
OUT = OUT_ROOT / "task3_analysis"

# Ordered metric names
LRG_METRICS = [
    "matrix_distance", "scaled_distance", "rank_distance", "quantile_rmse",
    "permutation_robust", "tree_robinson_foulds", "tree_cophenetic_corr",
    "tree_baker_gamma", "tree_fowlkes_mallows", "ari", "cluster_swap",
]
FC_METRICS = [
    "fc_frobenius", "fc_scaled_frobenius", "fc_rank_distance", "fc_mean_abs_diff",
]
ALL_METRICS = LRG_METRICS + FC_METRICS

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = [BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER]

PAIR_LABELS = [f"{a[:2]}{a[-4:]}-{b[:2]}{b[-4:]}" for a, b in ALL_PAIRS]

# Map metric to short display label
METRIC_SHORT = {
    "matrix_distance": "mat_dist",
    "scaled_distance": "scl_dist",
    "rank_distance": "rnk_dist",
    "quantile_rmse": "qnt_rmse",
    "permutation_robust": "perm_rob",
    "tree_robinson_foulds": "RF",
    "tree_cophenetic_corr": "coph_cor",
    "tree_baker_gamma": "baker_γ",
    "tree_fowlkes_mallows": "FM",
    "ari": "ARI",
    "cluster_swap": "swap",
    "fc_frobenius": "fc_frob",
    "fc_scaled_frobenius": "fc_s_frob",
    "fc_rank_distance": "fc_rnk",
    "fc_mean_abs_diff": "fc_mad",
}


def load_data() -> pd.DataFrame:
    df = pd.read_csv(CSV_IN)
    # Add distance-convention column
    df["dist_value"] = df.apply(
        lambda r: to_distance(r["value"], _METRIC_META[r["metric_name"]][1]),
        axis=1,
    )
    return df


# ══════════════════════════════════════════════════════════════════════════
# 3a — Inter-metric redundancy (Spearman correlation matrix)
# ══════════════════════════════════════════════════════════════════════════
def analysis_3a(df: pd.DataFrame):
    print("\n3a — Inter-metric redundancy...")
    out = OUT / "metric_correlation_matrix"

    # Pivot: each row = (patient, band, phase_a, phase_b), columns = metrics
    pivot = df.pivot_table(
        index=["patient", "band", "phase_a", "phase_b"],
        columns="metric_name",
        values="dist_value",
    ).dropna(axis=1, how="all")

    metrics_present = [m for m in ALL_METRICS if m in pivot.columns]
    pivot = pivot[metrics_present]

    n = len(metrics_present)
    corr_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            mask = pivot.iloc[:, i].notna() & pivot.iloc[:, j].notna()
            if mask.sum() > 3:
                corr_matrix[i, j] = spearmanr(
                    pivot.iloc[:, i][mask], pivot.iloc[:, j][mask]
                )[0]
            else:
                corr_matrix[i, j] = np.nan

    labels = [METRIC_SHORT.get(m, m) for m in metrics_present]

    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(corr_matrix, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=8)

    # Annotate cells
    for i in range(n):
        for j in range(n):
            v = corr_matrix[i, j]
            if not np.isnan(v):
                color = "white" if abs(v) > 0.7 else "black"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=5.5, color=color)

    # Draw separator between LRG and FC metrics
    n_lrg = sum(1 for m in metrics_present if m in LRG_METRICS)
    ax.axhline(n_lrg - 0.5, color="black", linewidth=2)
    ax.axvline(n_lrg - 0.5, color="black", linewidth=2)

    fig.colorbar(im, ax=ax, label="Spearman ρ", shrink=0.8)
    ax.set_title("Inter-metric Spearman correlation (distance convention)")
    fig.tight_layout()

    save_fig(fig, out.with_suffix(".pdf"),
        what="Spearman rank correlation matrix between all 15 metrics (11 LRG + 4 raw FC). "
             "All similarity metrics converted to distance convention (negated). "
             "Computed across all (patient, band, phase-pair) = 150 data points.",
        proves="Identifies redundant metric pairs (|ρ| > 0.9) and independent information axes. "
               "Black lines separate LRG-based from raw FC metrics.",
        how_to_read="Red = positively correlated (redundant), Blue = negatively correlated. "
                    "Look for blocks of high correlation (clusters of redundant metrics). "
                    "Metrics that are weakly correlated with all others carry unique information.",
    )

    # Also save the correlation values as markdown table
    md_lines = ["# Inter-metric Spearman correlation", "",
                "Computed in distance convention across 150 (patient, band, phase-pair) cells.", "",
                "## Highly redundant pairs (|ρ| > 0.9):", ""]
    for i in range(n):
        for j in range(i + 1, n):
            if abs(corr_matrix[i, j]) > 0.9:
                md_lines.append(
                    f"- `{metrics_present[i]}` ↔ `{metrics_present[j]}`: "
                    f"ρ = {corr_matrix[i, j]:.3f}"
                )
    md_lines += ["", "## Weakly correlated with others (max |ρ| < 0.5):", ""]
    for i in range(n):
        off_diag = [abs(corr_matrix[i, j]) for j in range(n) if j != i and not np.isnan(corr_matrix[i, j])]
        if off_diag and max(off_diag) < 0.5:
            md_lines.append(f"- `{metrics_present[i]}`: max |ρ| = {max(off_diag):.3f}")

    (OUT / "metric_correlation_matrix.md").write_text("\n".join(md_lines))


# ══════════════════════════════════════════════════════════════════════════
# 3b — Claim A: Band dependence
# ══════════════════════════════════════════════════════════════════════════
def analysis_3b(df: pd.DataFrame):
    print("\n3b — Claim A: Band dependence...")
    out_dir = OUT / "claim_A_band_dependence"
    out_dir.mkdir(parents=True, exist_ok=True)

    # For each metric & patient: mean reorganization across phase pairs, per band
    agg = df.groupby(["metric_name", "patient", "band"])["dist_value"].mean().reset_index()

    # -- Figure 1: heatmap band × patient for each metric --
    n_metrics = len(ALL_METRICS)
    n_cols = 5
    n_rows = (n_metrics + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.5 * n_rows))
    axes = axes.flatten()

    for idx, mname in enumerate(ALL_METRICS):
        ax = axes[idx]
        sub = agg[agg["metric_name"] == mname]
        piv = sub.pivot(index="band", columns="patient", values="dist_value")
        piv = piv.reindex(BAND_ORDER)

        if piv.empty:
            ax.set_visible(False)
            continue

        im = ax.imshow(piv.values, aspect="auto", cmap="viridis")
        ax.set_yticks(range(len(BAND_ORDER)))
        ax.set_yticklabels(BAND_TEX, fontsize=8)
        ax.set_xticks(range(len(piv.columns)))
        ax.set_xticklabels(piv.columns, fontsize=7, rotation=45, ha="right")
        ax.set_title(METRIC_SHORT.get(mname, mname), fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.7)

    for idx in range(len(ALL_METRICS), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle("Claim A: Mean reorganization per band × patient", fontsize=13, y=1.01)
    fig.tight_layout()
    save_fig(fig, out_dir / "band_ranking_per_patient.pdf",
        what="Heatmap grid (one panel per metric). Each panel shows bands (rows) vs patients (columns). "
             "Color = mean reorganization (distance) averaged over all 6 phase pairs.",
        proves="Claim A (band-dependence): if reorganization varies across bands, rows show different colors. "
               "If consistent across patients, column patterns are similar.",
        how_to_read="Brighter = more reorganization. Look for consistent row ordering across patients. "
                    "If delta is always brightest and high_gamma always darkest (or vice versa), "
                    "band dependence is confirmed.",
    )

    # -- Figure 2: raw FC vs LRG band ranking comparison --
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, (level, metric_list) in zip(axes, [("LRG", LRG_METRICS), ("Raw FC", FC_METRICS)]):
        # Average across metrics in level, then show band ranking per patient
        sub = agg[agg["metric_name"].isin(metric_list)]
        piv = sub.groupby(["patient", "band"])["dist_value"].mean().reset_index()
        piv2 = piv.pivot(index="band", columns="patient", values="dist_value")
        piv2 = piv2.reindex(BAND_ORDER)

        # Rank within each patient (1 = highest reorganization)
        ranks = piv2.rank(ascending=False)

        im = ax.imshow(ranks.values, aspect="auto", cmap="RdYlGn_r", vmin=1, vmax=6)
        ax.set_yticks(range(len(BAND_ORDER)))
        ax.set_yticklabels(BAND_TEX, fontsize=10)
        ax.set_xticks(range(len(ranks.columns)))
        ax.set_xticklabels(ranks.columns, fontsize=9)
        ax.set_title(f"{level} metrics (averaged)", fontsize=11)

        for i in range(ranks.shape[0]):
            for j in range(ranks.shape[1]):
                ax.text(j, i, f"{ranks.values[i, j]:.0f}",
                        ha="center", va="center", fontsize=9, fontweight="bold")

        fig.colorbar(im, ax=ax, label="Rank (1=most reorganization)", shrink=0.8)

    fig.suptitle("Claim A: Band ranking — LRG vs Raw FC", fontsize=13)
    fig.tight_layout()
    save_fig(fig, out_dir / "raw_fc_vs_lrg_comparison.pdf",
        what="Two panels: LRG metrics (left) and Raw FC metrics (right). "
             "Each shows band rank (1=most reorganization) per patient. "
             "Ranks computed from mean distance value averaged across all metrics within the level.",
        proves="Whether LRG and raw FC produce the same band ordering. "
               "If they differ, LRG reveals structure invisible to raw FC.",
        how_to_read="Red = rank 1 (most reorganization), Green = rank 6 (least). "
                    "Consistent coloring across patients within a panel = reproducible band dependence. "
                    "Different patterns between panels = LRG adds value.",
    )

    # Band ranking summary table
    md_lines = ["# Claim A — Band Dependence Summary", ""]
    for level, mlist in [("LRG", LRG_METRICS), ("Raw FC", FC_METRICS)]:
        sub = agg[agg["metric_name"].isin(mlist)]
        overall = sub.groupby("band")["dist_value"].mean().reindex(BAND_ORDER)
        ranking = overall.rank(ascending=False).astype(int)
        md_lines.append(f"## {level} metrics — overall band ranking (1=most reorganization)")
        md_lines.append("")
        md_lines.append("| Band | Mean distance | Rank |")
        md_lines.append("|------|--------------|------|")
        for b in BAND_ORDER:
            md_lines.append(f"| {b} | {overall[b]:.4f} | {ranking[b]} |")
        md_lines.append("")

    (out_dir / "band_ranking_summary.md").write_text("\n".join(md_lines))


# ══════════════════════════════════════════════════════════════════════════
# 3c — Claim B: Task-state signature (taskL-taskT = minimum distance)
# ══════════════════════════════════════════════════════════════════════════
def analysis_3c(df: pd.DataFrame):
    print("\n3c — Claim B: Task-state signature...")
    out_dir = OUT / "claim_B_task_signature"
    out_dir.mkdir(parents=True, exist_ok=True)

    # For each (metric, patient, band): which phase pair has minimum distance?
    # For similarity metrics (is_distance=False), min distance = max similarity → use dist_value (negated)
    # Since dist_value = -similarity, min dist_value = max similarity → consistent: argmin dist_value

    results = []
    for mname in ALL_METRICS:
        is_dist = _METRIC_META[mname][1]
        sub = df[df["metric_name"] == mname]
        for pat in PATIENTS:
            for band in BAND_ORDER:
                cell = sub[(sub["patient"] == pat) & (sub["band"] == band)]
                if cell.empty:
                    continue
                if is_dist:
                    # Distance: min value = most similar pair
                    min_idx = cell["value"].idxmin()
                else:
                    # Similarity: max value = most similar pair
                    min_idx = cell["value"].idxmax()
                min_row = cell.loc[min_idx]
                pair = (min_row["phase_a"], min_row["phase_b"])
                is_task_pair = pair == ("task_learn", "task_test") or pair == ("task_test", "task_learn")
                results.append({
                    "metric": mname,
                    "patient": pat,
                    "band": band,
                    "most_similar_pair": f"{pair[0]}-{pair[1]}",
                    "is_task_pair": is_task_pair,
                    "level": "lrg" if mname in LRG_METRICS else "raw_fc",
                })

    res_df = pd.DataFrame(results)
    total_cells = len(PATIENTS) * len(BAND_ORDER)

    # Tally per metric: how many cells have taskL-taskT as most similar?
    tally = res_df.groupby(["metric", "level"])["is_task_pair"].sum().reset_index()
    tally.columns = ["metric", "level", "count"]
    tally = tally.sort_values("count", ascending=False)

    # Bar chart
    fig, ax = plt.subplots(figsize=(14, 5))
    metrics_ordered = [m for m in ALL_METRICS if m in tally["metric"].values]
    counts = [tally[tally["metric"] == m]["count"].iloc[0] for m in metrics_ordered]
    colors = ["#2196F3" if m in LRG_METRICS else "#FF9800" for m in metrics_ordered]
    labels_short = [METRIC_SHORT.get(m, m) for m in metrics_ordered]

    bars = ax.bar(range(len(metrics_ordered)), counts, color=colors)
    ax.set_xticks(range(len(metrics_ordered)))
    ax.set_xticklabels(labels_short, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel(f"# cells where taskL–taskT is most similar (out of {total_cells})")
    ax.set_title("Claim B: Task-state signature — taskL–taskT as minimum distance pair")
    ax.axhline(total_cells, color="gray", ls="--", alpha=0.5, label=f"Maximum ({total_cells})")
    ax.axhline(total_cells / 6, color="red", ls=":", alpha=0.5, label=f"Chance (1/6 = {total_cells/6:.0f})")
    ax.legend(fontsize=8)

    # Color legend
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor="#2196F3", label="LRG metrics"),
        Patch(facecolor="#FF9800", label="Raw FC metrics"),
        plt.Line2D([0], [0], color="gray", ls="--", label=f"Maximum ({total_cells})"),
        plt.Line2D([0], [0], color="red", ls=":", label=f"Chance ({total_cells/6:.0f})"),
    ], fontsize=8, loc="upper right")

    fig.tight_layout()
    save_fig(fig, out_dir / "min_distance_tally.pdf",
        what=f"Bar chart: for each metric, how many of the {total_cells} (patient × band) cells "
             "have task_learn–task_test as the most similar phase pair (minimum distance or maximum similarity). "
             "Blue = LRG metrics, Orange = raw FC metrics. Gray dashed = maximum possible, Red dotted = chance level (1/6).",
        proves="Claim B (task-state signature): if taskL and taskT share a common FC architecture, "
               "they should be the most similar pair more often than chance. "
               "Higher bars = stronger task signature. Compare LRG vs FC to test claim E.",
        how_to_read="Bars well above the red chance line confirm claim B. "
                    "If LRG bars are higher than FC bars, LRG detects the task signature better.",
    )

    # Summary table
    md_lines = ["# Claim B — Task-State Signature Tally", "",
                f"Total cells: {total_cells} (5 patients × 6 bands)", "",
                "| Metric | Level | Count | Fraction |",
                "|--------|-------|-------|----------|"]
    for _, row in tally.iterrows():
        frac = row["count"] / total_cells
        md_lines.append(f"| `{row['metric']}` | {row['level']} | {int(row['count'])} | {frac:.2f} |")
    md_lines.append(f"\nChance level: {1/6:.3f}")
    (out_dir / "min_distance_tally.md").write_text("\n".join(md_lines))


# ══════════════════════════════════════════════════════════════════════════
# 3d — Claim C: Task imprint on rest_post
# ══════════════════════════════════════════════════════════════════════════
def analysis_3d(df: pd.DataFrame):
    print("\n3d — Claim C: Task imprint on rest_post...")
    out_dir = OUT / "claim_C_task_imprint"
    out_dir.mkdir(parents=True, exist_ok=True)

    # For each (metric, patient, band): compare d(rest_pre, rest_post) vs d(taskL, taskT)
    results = []
    for mname in ALL_METRICS:
        is_dist = _METRIC_META[mname][1]
        sub = df[df["metric_name"] == mname]

        for pat in PATIENTS:
            for band in BAND_ORDER:
                cell = sub[(sub["patient"] == pat) & (sub["band"] == band)]
                if cell.empty:
                    continue

                def _get_val(pa, pb):
                    row = cell[(cell["phase_a"] == pa) & (cell["phase_b"] == pb)]
                    if row.empty:
                        row = cell[(cell["phase_a"] == pb) & (cell["phase_b"] == pa)]
                    if row.empty:
                        return np.nan
                    return row["dist_value"].iloc[0]

                d_rest = _get_val("rest_pre", "rest_post")
                d_task = _get_val("task_learn", "task_test")
                d_taskT_rsPost = _get_val("task_test", "rest_post")

                results.append({
                    "metric": mname,
                    "patient": pat,
                    "band": band,
                    "d_rsPre_rsPost": d_rest,
                    "d_taskL_taskT": d_task,
                    "d_taskT_rsPost": d_taskT_rsPost,
                    "task_imprint": d_rest > d_task,  # rest_post changed more than task phases differ
                    "level": "lrg" if mname in LRG_METRICS else "raw_fc",
                })

    res_df = pd.DataFrame(results)
    total_cells = len(PATIENTS) * len(BAND_ORDER)

    # Scatter: d(rest_pre, rest_post) vs d(taskL, taskT) for selected metrics
    # Select a diverse subset: 2 best LRG + 2 best FC
    tally = res_df.groupby("metric")["task_imprint"].sum().sort_values(ascending=False)
    top_lrg = [m for m in tally.index if m in LRG_METRICS][:3]
    top_fc = [m for m in tally.index if m in FC_METRICS][:2]
    selected = top_lrg + top_fc

    n_sel = len(selected)
    fig, axes = plt.subplots(1, n_sel, figsize=(4.5 * n_sel, 4))
    if n_sel == 1:
        axes = [axes]
    band_colors = dict(zip(BAND_ORDER, plt.cm.tab10(np.linspace(0, 1, 6))))

    for ax, mname in zip(axes, selected):
        sub = res_df[res_df["metric"] == mname]
        for band in BAND_ORDER:
            bsub = sub[sub["band"] == band]
            ax.scatter(bsub["d_taskL_taskT"], bsub["d_rsPre_rsPost"],
                       c=[band_colors[band]], label=band, s=50, alpha=0.8,
                       edgecolors="black", linewidth=0.5)

        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1]),
        ]
        ax.plot(lims, lims, "k--", alpha=0.4, lw=1)
        ax.set_xlabel("d(taskL, taskT)", fontsize=9)
        ax.set_ylabel("d(rest_pre, rest_post)", fontsize=9)
        n_imprint = sub["task_imprint"].sum()
        ax.set_title(f"{METRIC_SHORT.get(mname, mname)} ({n_imprint}/{total_cells})", fontsize=10)

    axes[0].legend(fontsize=7, loc="lower right", ncol=2)
    fig.suptitle("Claim C: Task imprint — d(rest_pre,rest_post) vs d(taskL,taskT)", fontsize=12, y=1.02)
    fig.tight_layout()
    save_fig(fig, out_dir / "rsPre_rsPost_vs_taskL_taskT.pdf",
        what="Scatter plots for top-performing metrics. X-axis = d(taskL, taskT), "
             "Y-axis = d(rest_pre, rest_post). Points colored by frequency band. "
             "Diagonal dashed line = equality. Title shows count of cells above the line.",
        proves="Claim C (task imprint): points above the diagonal mean rest_post has changed more "
               "from rest_pre than the two task phases differ from each other — the task left a trace. "
               "If most points are above the line, the task experience altered resting-state organization.",
        how_to_read="Points above diagonal = task imprint detected. "
                    "Higher count in title = metric detects the imprint more reliably. "
                    "Band coloring shows whether the imprint is band-dependent.",
    )

    # Summary table
    md_lines = ["# Claim C — Task Imprint Summary", "",
                f"Total cells: {total_cells} (5 patients × 6 bands)", "",
                "Task imprint detected when d(rest_pre, rest_post) > d(taskL, taskT) [distance convention].", "",
                "| Metric | Level | Imprint count | Fraction |",
                "|--------|-------|---------------|----------|"]
    for mname in ALL_METRICS:
        sub = res_df[res_df["metric"] == mname]
        n_imp = sub["task_imprint"].sum()
        frac = n_imp / total_cells
        level = "lrg" if mname in LRG_METRICS else "raw_fc"
        md_lines.append(f"| `{mname}` | {level} | {n_imp} | {frac:.2f} |")

    md_lines += ["", "## Per-band breakdown (all metrics pooled)", ""]
    for band in BAND_ORDER:
        sub = res_df[res_df["band"] == band]
        n_imp = sub["task_imprint"].sum()
        total = len(sub)
        md_lines.append(f"- **{band}**: {n_imp}/{total} = {n_imp/max(total,1):.2f}")

    (out_dir / "rsPre_rsPost_vs_taskL_taskT.md").write_text("\n".join(md_lines))


# ══════════════════════════════════════════════════════════════════════════
# 3e — Normalization check
# ══════════════════════════════════════════════════════════════════════════
def analysis_3e(df: pd.DataFrame):
    print("\n3e — Normalization check...")

    md_lines = [
        "# Task 3e — Normalization and Interpretability Check", "",
        "| Metric | Range | Bounded? | Reference value | Cross-patient comparable? | Recommendation |",
        "|--------|-------|----------|-----------------|--------------------------|----------------|",
    ]

    for mname in ALL_METRICS:
        meta = _METRIC_META[mname]
        label, is_dist, rng, ref = meta
        sub = df[df["metric_name"] == mname]["value"]
        actual_min = sub.min()
        actual_max = sub.max()

        bounded = "∞" not in rng
        if bounded:
            comparable = "Yes"
            recommendation = "Use directly"
        elif "1 −" in label or "rank" in mname.lower():
            comparable = "Yes (rank-based)"
            recommendation = "Use directly"
        else:
            comparable = "No — depends on N and scale"
            recommendation = "Normalize per patient or use rank-based alternative"

        md_lines.append(
            f"| `{mname}` | {rng} | {'Yes' if bounded else 'No'} | "
            f"{ref} | {comparable} | {recommendation} |"
        )

    md_lines += [
        "", "## Actual value ranges in dataset", "",
        "| Metric | Min | Max | Mean | Std |",
        "|--------|-----|-----|------|-----|",
    ]
    for mname in ALL_METRICS:
        sub = df[df["metric_name"] == mname]["value"]
        md_lines.append(
            f"| `{mname}` | {sub.min():.4f} | {sub.max():.4f} | "
            f"{sub.mean():.4f} | {sub.std():.4f} |"
        )

    md_lines += [
        "", "## Metrics suitable for cross-patient comparison:", "",
        "Bounded and interpretable: `scaled_distance`, `rank_distance`, "
        "`tree_robinson_foulds`, `tree_cophenetic_corr`, `tree_baker_gamma`, "
        "`tree_fowlkes_mallows`, `ari`, `cluster_swap`, `fc_scaled_frobenius`, `fc_rank_distance`", "",
        "Unbounded (need normalization): `matrix_distance`, `quantile_rmse`, "
        "`permutation_robust`, `fc_frobenius`, `fc_mean_abs_diff`",
    ]

    (OUT / "normalization_check.md").write_text("\n".join(md_lines))


# ══════════════════════════════════════════════════════════════════════════
# 3f — Cross-patient consistency
# ══════════════════════════════════════════════════════════════════════════
def analysis_3f(df: pd.DataFrame):
    print("\n3f — Cross-patient consistency...")
    out_dir = OUT / "cross_patient_consistency"
    out_dir.mkdir(parents=True, exist_ok=True)

    # CoV = std/|mean| across patients for each (metric, band, phase_pair)
    cov_rows = []
    for mname in ALL_METRICS:
        sub = df[df["metric_name"] == mname]
        for band in BAND_ORDER:
            for pa, pb in ALL_PAIRS:
                cell = sub[(sub["band"] == band) &
                           (sub["phase_a"] == pa) & (sub["phase_b"] == pb)]
                vals = cell["dist_value"].values
                if len(vals) > 1:
                    mean_val = np.mean(vals)
                    std_val = np.std(vals)
                    cov = std_val / max(abs(mean_val), 1e-15)
                    cov_rows.append({
                        "metric": mname,
                        "band": band,
                        "pair": f"{pa[:2]}{pa[-4:]}-{pb[:2]}{pb[-4:]}",
                        "cov": cov,
                        "mean": mean_val,
                        "std": std_val,
                        "n_patients": len(vals),
                    })

    cov_df = pd.DataFrame(cov_rows)

    # Heatmap: metric × (band × pair)
    # Average CoV across pairs for each (metric, band)
    avg_cov = cov_df.groupby(["metric", "band"])["cov"].mean().reset_index()
    piv = avg_cov.pivot(index="metric", columns="band", values="cov")
    piv = piv.reindex(index=ALL_METRICS, columns=BAND_ORDER)

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(piv.values, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(BAND_ORDER)))
    ax.set_xticklabels(BAND_TEX, fontsize=10)
    ax.set_yticks(range(len(ALL_METRICS)))
    ax.set_yticklabels([METRIC_SHORT.get(m, m) for m in ALL_METRICS], fontsize=8)

    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if not np.isnan(v):
                color = "white" if v > 0.8 else "black"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=7, color=color)

    # Separator between LRG and FC
    n_lrg = len(LRG_METRICS)
    ax.axhline(n_lrg - 0.5, color="black", linewidth=2)

    fig.colorbar(im, ax=ax, label="Coefficient of Variation (across patients)", shrink=0.8)
    ax.set_title("Cross-patient consistency: CoV per metric × band")
    fig.tight_layout()
    save_fig(fig, out_dir / "coefficient_of_variation.pdf",
        what="Heatmap of coefficient of variation (std/|mean| across patients) for each metric × band. "
             "Lower CoV = more consistent across patients. "
             "Values averaged across the 6 phase pairs.",
        proves="A good metric should have LOW CoV within each condition (band × pair), "
               "meaning patients agree on the relative reorganization pattern. "
               "High CoV = metric is noisy or patient-dependent.",
        how_to_read="Yellow = low CoV (good consistency), Red = high CoV (poor consistency). "
                    "Metrics with uniformly low CoV are reliable for cross-patient reporting.",
    )

    # Signal-to-noise ratio: between-condition variance / within-condition variance
    md_lines = [
        "# Cross-patient consistency — Signal-to-Noise Ratio", "",
        "SNR = var(condition means) / mean(within-condition variance across patients)", "",
        "Higher SNR = metric differentiates conditions (band × pair) while patients agree.", "",
        "| Metric | Level | SNR_band | SNR_pair | SNR_combined | Mean CoV |",
        "|--------|-------|----------|----------|--------------|----------|",
    ]

    for mname in ALL_METRICS:
        sub = df[df["metric_name"] == mname]
        level = "lrg" if mname in LRG_METRICS else "raw_fc"

        # SNR for bands: var of band means / mean of within-band patient variance
        band_means = sub.groupby("band")["dist_value"].mean()
        band_vars = sub.groupby("band")["dist_value"].var()
        snr_band = band_means.var() / max(band_vars.mean(), 1e-15)

        # SNR for pairs
        sub_with_pair = sub.copy()
        sub_with_pair["pair"] = sub_with_pair["phase_a"] + "-" + sub_with_pair["phase_b"]
        pair_means = sub_with_pair.groupby("pair")["dist_value"].mean()
        pair_vars = sub_with_pair.groupby("pair")["dist_value"].var()
        snr_pair = pair_means.var() / max(pair_vars.mean(), 1e-15)

        snr_combined = (snr_band + snr_pair) / 2

        mean_cov = cov_df[cov_df["metric"] == mname]["cov"].mean()

        md_lines.append(
            f"| `{mname}` | {level} | {snr_band:.3f} | {snr_pair:.3f} | "
            f"{snr_combined:.3f} | {mean_cov:.3f} |"
        )

    md_lines += [
        "", "## Interpretation", "",
        "- **SNR_band**: how well the metric differentiates frequency bands",
        "- **SNR_pair**: how well the metric differentiates phase pairs",
        "- **SNR_combined**: average of both",
        "- **Mean CoV**: average coefficient of variation across all conditions",
        "", "Higher SNR + lower CoV = best metric for cross-patient reporting.",
    ]

    (out_dir / "signal_to_noise.md").write_text("\n".join(md_lines))


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("Loading data...")
    df = load_data()
    print(f"  {len(df)} rows loaded")

    analysis_3a(df)
    analysis_3b(df)
    analysis_3c(df)
    analysis_3d(df)
    analysis_3e(df)
    analysis_3f(df)

    print(f"\nTask 3 complete. Outputs in {OUT}")


if __name__ == "__main__":
    main()
