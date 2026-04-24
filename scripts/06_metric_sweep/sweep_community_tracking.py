#!/usr/bin/env python3
"""Community tracking across phases via Hungarian-matched fcluster partitions.

For each pair of dendrograms (same patient, same band, different phases):
1. Cut both dendrograms at the same k using fcluster(Z, k, criterion='maxclust')
2. Use Hungarian matching to align communities (maximize overlap)
3. Compute node-level stability metrics

Metrics (at fixed k):
  - fraction_stable_k: fraction of nodes that stay in matched community
  - jaccard_mean_k: mean Jaccard similarity between matched communities
  - largest_community_stability_k: stability of the largest community only
  - community_size_correlation_k: Spearman corr of sorted community size vectors

Multi-k aggregate:
  - mean_fraction_stable: average fraction_stable across k=2..10
  - min_fraction_stable: worst-case across k values

Natural k (from largest gap in merge heights):
  - fraction_stable_natural: fraction_stable at natural k

Output -> data/figures/metric_exploration/community_tracking/
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from scipy.cluster.hierarchy import fcluster
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
PATIENTS = list(PATIENT_PHASES.keys())
BANDS = BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import LRG_CACHE
OUT_DIR = FIGURES_ROOT / "metric_exploration" / "community_tracking"

K_VALUES = [2, 3, 4, 5, 6, 8, 10]

REST = {"rest_pre", "rest_post"}
TASK = {"task_learn", "task_test"}

# Phase pair display labels
PAIR_SHORT = {
    ("rest_pre", "task_learn"): "rest_pre-tL",
    ("rest_pre", "task_test"): "rest_pre-tT",
    ("rest_pre", "rest_post"): "rest_pre-rsP",
    ("task_learn", "task_test"): "tL-tT",
    ("task_learn", "rest_post"): "tL-rsP",
    ("task_test", "rest_post"): "tT-rsP",
}


def classify_pair(p1: str, p2: str) -> str:
    s = {p1, p2}
    if s <= REST or s <= TASK:
        return "within"
    return "cross"


# ---------------------------------------------------------------------------
# Core algorithms
# ---------------------------------------------------------------------------
def match_communities(labels1: np.ndarray, labels2: np.ndarray):
    """Match communities between two partitions using the Hungarian algorithm.

    Returns
    -------
    matched_labels2 : np.ndarray
        labels2 relabeled to maximally align with labels1
    fraction_same : float
        fraction of nodes where labels1 == matched_labels2
    jaccard_per_comm : list[float]
        Jaccard similarity for each matched community pair
    """
    unique1 = np.unique(labels1)
    unique2 = np.unique(labels2)
    k1, k2 = len(unique1), len(unique2)
    k_max = max(k1, k2)

    # Cost matrix: -overlap (we maximize overlap)
    cost = np.zeros((k_max, k_max))
    for i, c1 in enumerate(unique1):
        for j, c2 in enumerate(unique2):
            cost[i, j] = -np.sum((labels1 == c1) & (labels2 == c2))

    row_ind, col_ind = linear_sum_assignment(cost)

    # Build relabeling map
    label_map = {}
    for i, j in zip(row_ind, col_ind):
        if i < k1 and j < k2:
            label_map[unique2[j]] = unique1[i]

    matched = np.array([label_map.get(l, -1) for l in labels2])
    fraction_same = np.mean(labels1 == matched)

    # Per-community Jaccard
    jaccard_per_comm = []
    for c in unique1:
        s1 = set(np.where(labels1 == c)[0])
        s2 = set(np.where(matched == c)[0])
        if len(s1 | s2) == 0:
            jaccard_per_comm.append(0.0)
        else:
            jaccard_per_comm.append(len(s1 & s2) / len(s1 | s2))

    return matched, fraction_same, jaccard_per_comm


def find_natural_k(Z: np.ndarray, k_min: int = 2, k_max: int = 15) -> int:
    """Find natural number of clusters from largest gap in merge heights."""
    heights = Z[:, 2]
    if len(heights) < 2:
        return 2
    # Look at gaps between consecutive merge heights (top of dendrogram)
    # We want to find the largest gap, which suggests natural cluster count
    n = len(heights)
    # The number of clusters at merge step i (from top) is n+1-i
    # So we want the gap between height[i] and height[i-1]
    gaps = np.diff(heights)
    # Only consider the upper portion (last merges correspond to small k)
    # merge index i corresponds to k = n+1 - i clusters
    # We want k in [k_min, k_max]
    best_k = 2
    best_gap = -1.0
    for i in range(len(gaps)):
        k_at_merge = n - i  # number of clusters just before this merge
        if k_min <= k_at_merge <= k_max:
            if gaps[i] > best_gap:
                best_gap = gaps[i]
                best_k = k_at_merge
    return best_k


def community_size_correlation(labels1: np.ndarray, labels2: np.ndarray) -> float:
    """Spearman correlation of sorted community size vectors."""
    _, counts1 = np.unique(labels1, return_counts=True)
    _, counts2 = np.unique(labels2, return_counts=True)
    sizes1 = np.sort(counts1)[::-1]
    sizes2 = np.sort(counts2)[::-1]
    # Pad shorter vector with zeros
    max_len = max(len(sizes1), len(sizes2))
    s1 = np.zeros(max_len)
    s2 = np.zeros(max_len)
    s1[:len(sizes1)] = sizes1
    s2[:len(sizes2)] = sizes2
    if np.std(s1) == 0 or np.std(s2) == 0:
        return 0.0
    rho, _ = spearmanr(s1, s2)
    return rho


def largest_community_stability(labels1: np.ndarray, matched2: np.ndarray) -> float:
    """Stability of the largest community (by size in labels1)."""
    unique, counts = np.unique(labels1, return_counts=True)
    largest = unique[np.argmax(counts)]
    s1 = set(np.where(labels1 == largest)[0])
    s2 = set(np.where(matched2 == largest)[0])
    if len(s1 | s2) == 0:
        return 0.0
    return len(s1 & s2) / len(s1 | s2)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_all_linkages():
    """Load linkage matrices for all (patient, phase, band)."""
    data = {}
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for ph in phases:
                lrg = load_lrg_result(pat, ph, band, "msc", cache_root=LRG_CACHE)
                if lrg is None:
                    print(f"  WARNING: missing {pat} {ph} {band}")
                    continue
                data[(pat, ph, band)] = lrg.linkage_matrix
    return data


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------
def compute_all_metrics(data: dict) -> pd.DataFrame:
    """Compute all community tracking metrics for all pairs."""
    rows = []

    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            # Check all phases are available
            available_phases = [ph for ph in phases if (pat, ph, band) in data]
            if len(available_phases) < 2:
                continue

            for ph1, ph2 in combinations(available_phases, 2):
                Z1 = data[(pat, ph1, band)]
                Z2 = data[(pat, ph2, band)]

                n1 = Z1.shape[0] + 1
                n2 = Z2.shape[0] + 1

                if n1 != n2:
                    print(f"  WARNING: n_nodes mismatch {pat} {band}: "
                          f"{ph1}={n1}, {ph2}={n2} -- skipping")
                    continue

                n_nodes = n1
                pair_key = (ph1, ph2) if (ph1, ph2) in PAIR_SHORT else (ph2, ph1)
                pair_label = PAIR_SHORT.get(pair_key, f"{ph1}-{ph2}")
                pair_type = classify_pair(ph1, ph2)

                # Natural k
                nat_k1 = find_natural_k(Z1)
                nat_k2 = find_natural_k(Z2)
                nat_k = max(nat_k1, nat_k2)  # use the larger one

                # --- Metrics at each k ---
                k_metrics = {}
                for k in K_VALUES:
                    if k >= n_nodes:
                        continue
                    lab1 = fcluster(Z1, k, criterion="maxclust")
                    lab2 = fcluster(Z2, k, criterion="maxclust")
                    matched2, frac_stable, jaccards = match_communities(lab1, lab2)

                    k_metrics[k] = {
                        "fraction_stable": frac_stable,
                        "jaccard_mean": np.mean(jaccards),
                        "largest_stability": largest_community_stability(lab1, matched2),
                        "size_correlation": community_size_correlation(lab1, lab2),
                    }

                # Multi-k aggregates
                frac_vals = [v["fraction_stable"] for v in k_metrics.values()]
                mean_frac = np.mean(frac_vals) if frac_vals else np.nan
                min_frac = np.min(frac_vals) if frac_vals else np.nan

                # Natural k metric
                nat_k_clamped = min(nat_k, n_nodes - 1)
                nat_k_clamped = max(nat_k_clamped, 2)
                lab1_nat = fcluster(Z1, nat_k_clamped, criterion="maxclust")
                lab2_nat = fcluster(Z2, nat_k_clamped, criterion="maxclust")
                _, frac_nat, jaccards_nat = match_communities(lab1_nat, lab2_nat)

                # Store row with per-k values and aggregates
                row = {
                    "patient": pat,
                    "band": band,
                    "phase1": ph1,
                    "phase2": ph2,
                    "pair_label": pair_label,
                    "pair_type": pair_type,
                    "n_nodes": n_nodes,
                    "natural_k": nat_k_clamped,
                }

                # Per-k metrics
                for k, km in k_metrics.items():
                    for mname, mval in km.items():
                        row[f"{mname}_k{k}"] = mval

                # Aggregates
                row["mean_fraction_stable"] = mean_frac
                row["min_fraction_stable"] = min_frac
                row["fraction_stable_natural"] = frac_nat
                row["jaccard_mean_natural"] = np.mean(jaccards_nat)

                # Best k (highest fraction_stable)
                if k_metrics:
                    best_k = max(k_metrics, key=lambda k: k_metrics[k]["fraction_stable"])
                    row["best_k"] = best_k
                    row["fraction_stable_best_k"] = k_metrics[best_k]["fraction_stable"]
                    row["jaccard_mean_best_k"] = k_metrics[best_k]["jaccard_mean"]
                else:
                    row["best_k"] = np.nan
                    row["fraction_stable_best_k"] = np.nan
                    row["jaccard_mean_best_k"] = np.nan

                rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Hypothesis testing
# ---------------------------------------------------------------------------
def test_hypotheses(df: pd.DataFrame) -> pd.DataFrame:
    """Test the five hypotheses across patients and bands."""
    results = []

    # Only 4-phase patients for most tests
    df4 = df[df["patient"].isin(PATIENTS_4PHASE)].copy()

    for metric_col in ["mean_fraction_stable", "fraction_stable_best_k",
                       "fraction_stable_natural", "min_fraction_stable"]:

        for band in BANDS:
            dfb = df4[df4["band"] == band].copy()
            if dfb.empty:
                continue

            # H1: task_learn-task_test has highest fraction_stable
            h1_pass = 0
            h1_total = 0
            for pat in dfb["patient"].unique():
                dfp = dfb[dfb["patient"] == pat]
                tl_tt = dfp[dfp["pair_label"] == "tL-tT"][metric_col].values
                if len(tl_tt) == 0:
                    continue
                h1_total += 1
                if tl_tt[0] == dfp[metric_col].max():
                    h1_pass += 1

            # H2: task_test-rest_post > rest_pre-rest_post
            h2_pass = 0
            h2_total = 0
            for pat in dfb["patient"].unique():
                dfp = dfb[dfb["patient"] == pat]
                tt_rsp = dfp[dfp["pair_label"] == "tT-rsP"][metric_col].values
                rsp_rsp = dfp[dfp["pair_label"] == "rest_pre-rsP"][metric_col].values
                if len(tt_rsp) == 0 or len(rsp_rsp) == 0:
                    continue
                h2_total += 1
                if tt_rsp[0] > rsp_rsp[0]:
                    h2_pass += 1

            # H3: within > cross (mean)
            h3_pass = 0
            h3_total = 0
            for pat in dfb["patient"].unique():
                dfp = dfb[dfb["patient"] == pat]
                within_mean = dfp[dfp["pair_type"] == "within"][metric_col].mean()
                cross_mean = dfp[dfp["pair_type"] == "cross"][metric_col].mean()
                h3_total += 1
                if within_mean > cross_mean:
                    h3_pass += 1

            # H4: ALL pairs have high stability (fraction > 0.8)
            h4_pass = (dfb[metric_col] > 0.8).all()

            # H5: rest-rest pair has lowest stability
            h5_pass = 0
            h5_total = 0
            for pat in dfb["patient"].unique():
                dfp = dfb[dfb["patient"] == pat]
                rest_rest = dfp[dfp["pair_label"] == "rest_pre-rsP"][metric_col].values
                if len(rest_rest) == 0:
                    continue
                h5_total += 1
                if rest_rest[0] == dfp[metric_col].min():
                    h5_pass += 1

            results.append({
                "metric": metric_col,
                "band": band,
                "H1_taskLearn_taskTest_highest": f"{h1_pass}/{h1_total}",
                "H2_taskTest_rsPost_gt_rsPre_rsPost": f"{h2_pass}/{h2_total}",
                "H3_within_gt_cross": f"{h3_pass}/{h3_total}",
                "H4_all_high_stability": h4_pass,
                "H5_rest_rest_lowest": f"{h5_pass}/{h5_total}",
            })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_summary_heatmap(df: pd.DataFrame, out_dir: Path):
    """Band x pair heatmap of fraction_stable at best k, averaged over patients."""
    # All 6 canonical pair labels (for 4-phase patients)
    pair_order = ["rest_pre-tL", "rest_pre-tT", "rest_pre-rsP", "tL-tT", "tL-rsP", "tT-rsP"]
    band_tex = [BRAIN_BAND_TEX_DICT[b] for b in BANDS]

    # 4-phase patients only for the main heatmap
    df4 = df[df["patient"].isin(PATIENTS_4PHASE)].copy()

    fig, axes = plt.subplots(1, 2, figsize=(16, 5))

    for ax_idx, (metric, title) in enumerate([
        ("fraction_stable_best_k", "Fraction Stable (best k)"),
        ("mean_fraction_stable", "Mean Fraction Stable (k=2..10)")
    ]):
        ax = axes[ax_idx]
        matrix = np.full((len(BANDS), len(pair_order)), np.nan)

        for i, band in enumerate(BANDS):
            for j, pl in enumerate(pair_order):
                vals = df4[(df4["band"] == band) & (df4["pair_label"] == pl)][metric]
                if len(vals) > 0:
                    matrix[i, j] = vals.mean()

        im = ax.imshow(matrix, aspect="auto", cmap="RdYlGn", vmin=0.3, vmax=1.0)
        ax.set_xticks(range(len(pair_order)))
        ax.set_xticklabels(pair_order, rotation=45, ha="right", fontsize=9)
        ax.set_yticks(range(len(BANDS)))
        ax.set_yticklabels(band_tex, fontsize=11)
        ax.set_title(title, fontsize=12)

        # Annotate cells
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                val = matrix[i, j]
                if not np.isnan(val):
                    color = "white" if val < 0.5 else "black"
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                            fontsize=8, color=color)

        plt.colorbar(im, ax=ax, shrink=0.8)

    fig.suptitle("Community Tracking: Stability by Band and Phase Pair\n"
                 "(4-phase patients, averaged)", fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(out_dir / "summary_heatmap.pdf", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved summary_heatmap.pdf")


def plot_k_sweep(df: pd.DataFrame, out_dir: Path):
    """How fraction_stable changes with k, grouped by pair type."""
    df4 = df[df["patient"].isin(PATIENTS_4PHASE)].copy()

    fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharey=True)
    axes = axes.ravel()

    for idx, band in enumerate(BANDS):
        ax = axes[idx]
        dfb = df4[df4["band"] == band]

        # Group by pair type
        for ptype, color, marker in [("within", "C0", "o"), ("cross", "C1", "s")]:
            dfp = dfb[dfb["pair_type"] == ptype]
            if dfp.empty:
                continue

            means = []
            stds = []
            valid_ks = []
            for k in K_VALUES:
                col = f"fraction_stable_k{k}"
                if col in dfp.columns:
                    vals = dfp[col].dropna()
                    if len(vals) > 0:
                        means.append(vals.mean())
                        stds.append(vals.std())
                        valid_ks.append(k)

            if valid_ks:
                means = np.array(means)
                stds = np.array(stds)
                ax.plot(valid_ks, means, f"-{marker}", color=color,
                        label=ptype, markersize=5)
                ax.fill_between(valid_ks, means - stds, means + stds,
                                alpha=0.15, color=color)

        # Also plot individual pair labels for more detail
        pair_colors = {
            "tL-tT": "tab:green",
            "rest_pre-rsP": "tab:blue",
            "rest_pre-tL": "tab:orange",
            "tT-rsP": "tab:red",
        }
        for pl, pc in pair_colors.items():
            dfpl = dfb[dfb["pair_label"] == pl]
            if dfpl.empty:
                continue
            means = []
            valid_ks = []
            for k in K_VALUES:
                col = f"fraction_stable_k{k}"
                if col in dfpl.columns:
                    vals = dfpl[col].dropna()
                    if len(vals) > 0:
                        means.append(vals.mean())
                        valid_ks.append(k)
            if valid_ks:
                ax.plot(valid_ks, means, "--", color=pc, alpha=0.5,
                        label=pl, linewidth=1)

        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=12)
        ax.set_xlabel("k (number of communities)")
        if idx % 3 == 0:
            ax.set_ylabel("Fraction Stable")
        ax.set_ylim(0, 1.05)
        ax.axhline(0.5, color="gray", linestyle=":", alpha=0.5)
        if idx == 0:
            ax.legend(fontsize=7, loc="lower left", ncol=2)

    fig.suptitle("Community Stability vs k\n(4-phase patients, mean +/- std)",
                 fontsize=13, y=1.01)
    fig.tight_layout()
    fig.savefig(out_dir / "k_sweep.pdf", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved k_sweep.pdf")


def plot_node_tracking(df: pd.DataFrame, data: dict, out_dir: Path):
    """For the most interesting (patient, band), show node-level movements.

    Shows which nodes changed community between two phases.
    """
    # Find the case with the largest spread in fraction_stable across pairs
    df4 = df[df["patient"].isin(PATIENTS_4PHASE)].copy()

    best_spread = -1
    best_pat = None
    best_band = None

    for pat in df4["patient"].unique():
        for band in BANDS:
            vals = df4[(df4["patient"] == pat) & (df4["band"] == band)][
                "mean_fraction_stable"
            ]
            if len(vals) >= 4:
                spread = vals.max() - vals.min()
                if spread > best_spread:
                    best_spread = spread
                    best_pat = pat
                    best_band = band

    if best_pat is None:
        print("  No suitable case found for node_tracking")
        return

    print(f"  Node tracking for {best_pat} {best_band} (spread={best_spread:.3f})")

    phases = PATIENT_PHASES[best_pat]
    n_pairs = len(list(combinations(phases, 2)))

    # Choose a good k: use the best_k from the most stable pair
    dfpb = df4[(df4["patient"] == best_pat) & (df4["band"] == best_band)]
    best_k = int(dfpb["best_k"].mode().iloc[0]) if not dfpb.empty else 4

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.ravel()

    pair_idx = 0
    for ph1, ph2 in combinations(phases, 2):
        key1 = (best_pat, ph1, best_band)
        key2 = (best_pat, ph2, best_band)
        if key1 not in data or key2 not in data:
            continue

        Z1, Z2 = data[key1], data[key2]
        n = Z1.shape[0] + 1

        k_use = min(best_k, n - 1)
        lab1 = fcluster(Z1, k_use, criterion="maxclust")
        lab2 = fcluster(Z2, k_use, criterion="maxclust")
        matched2, frac, _ = match_communities(lab1, lab2)

        ax = axes[pair_idx]

        # Plot: each node colored by community in phase1, x=node index
        # Show which nodes stayed (circle) vs moved (X)
        stayed = lab1 == matched2
        moved = ~stayed

        cmap = plt.cm.Set2
        colors1 = [cmap(c % 8) for c in lab1]

        # Background: community assignment in phase1
        for node_idx in range(n):
            marker = "o" if stayed[node_idx] else "X"
            ax.scatter(node_idx, 0.5, c=[colors1[node_idx]],
                       marker=marker, s=30, edgecolors="k" if stayed[node_idx] else "red",
                       linewidths=0.5, zorder=2)
            # Draw an arrow for moved nodes showing new community
            if moved[node_idx]:
                new_color = cmap(matched2[node_idx] % 8)
                ax.scatter(node_idx, -0.5, c=[new_color],
                           marker="o", s=20, edgecolors="gray",
                           linewidths=0.3, alpha=0.7, zorder=2)
                ax.plot([node_idx, node_idx], [0.4, -0.4],
                        color="red", alpha=0.3, linewidth=0.5)

        pair_key = (ph1, ph2)
        pair_label = PAIR_SHORT.get(pair_key, f"{ph1}-{ph2}")

        ax.set_title(f"{pair_label}  (stable={frac:.2f}, k={k_use})", fontsize=10)
        ax.set_xlim(-1, n)
        ax.set_ylim(-1.2, 1.2)
        ax.set_yticks([0.5, -0.5])
        ax.set_yticklabels([ph1, ph2], fontsize=8)
        ax.set_xlabel("Node index", fontsize=8)
        ax.axhline(0, color="gray", linestyle=":", alpha=0.3)

        # Count moved
        n_moved = moved.sum()
        ax.text(0.98, 0.98, f"{n_moved}/{n} moved",
                transform=ax.transAxes, ha="right", va="top", fontsize=8,
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

        pair_idx += 1

    # Hide unused axes
    for i in range(pair_idx, len(axes)):
        axes[i].set_visible(False)

    fig.suptitle(f"Node Community Tracking: {best_pat} {BRAIN_BAND_TEX_DICT[best_band]}\n"
                 f"Top row = phase1 assignment (o=stayed, X=moved), "
                 f"Bottom row = phase2 reassignment",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(out_dir / "node_tracking.pdf", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved node_tracking.pdf")


def plot_consistency_table(hyp_df: pd.DataFrame, out_dir: Path):
    """Render hypothesis test results as a table figure."""
    # Pivot: for each metric, count how many bands pass each hypothesis
    metrics = hyp_df["metric"].unique()

    fig, ax = plt.subplots(figsize=(14, len(metrics) * len(BANDS) * 0.22 + 2))
    ax.axis("off")

    # Prepare table data
    col_labels = ["Metric", "Band", "H1: tL-tT highest", "H2: tT-rsP > rest_pre-rsP",
                  "H3: within > cross", "H4: all > 0.8", "H5: rest_pre-rsP lowest"]
    cell_data = []
    cell_colors = []

    for _, row in hyp_df.iterrows():
        cells = [
            row["metric"].replace("_", " "),
            BRAIN_BAND_TEX_DICT.get(row["band"], row["band"]),
            str(row["H1_taskLearn_taskTest_highest"]),
            str(row["H2_taskTest_rsPost_gt_rsPre_rsPost"]),
            str(row["H3_within_gt_cross"]),
            str(row["H4_all_high_stability"]),
            str(row["H5_rest_rest_lowest"]),
        ]
        cell_data.append(cells)

        # Color coding
        colors = ["white", "white"]
        for val_str in cells[2:]:
            if val_str == "True":
                colors.append("#90EE90")  # light green
            elif val_str == "False":
                colors.append("#FFB6C1")  # light pink
            elif "/" in val_str:
                parts = val_str.split("/")
                if int(parts[1]) > 0:
                    ratio = int(parts[0]) / int(parts[1])
                    if ratio >= 0.75:
                        colors.append("#90EE90")
                    elif ratio >= 0.5:
                        colors.append("#FFFACD")
                    else:
                        colors.append("#FFB6C1")
                else:
                    colors.append("white")
            else:
                colors.append("white")
        cell_colors.append(colors)

    if cell_data:
        table = ax.table(cellText=cell_data, colLabels=col_labels,
                         cellColours=cell_colors, loc="center",
                         cellLoc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(7)
        table.scale(1.0, 1.3)

        # Bold header
        for j in range(len(col_labels)):
            table[0, j].set_text_props(fontweight="bold")

    fig.suptitle("Hypothesis Consistency: Community Tracking Metrics", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_dir / "consistency_table.pdf", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved consistency_table.pdf")


def plot_per_patient_heatmaps(df: pd.DataFrame, out_dir: Path):
    """Per-patient heatmaps showing fraction_stable for each band x pair."""
    for pat in PATIENTS:
        dfp = df[df["patient"] == pat]
        if dfp.empty:
            continue

        phases = PATIENT_PHASES[pat]
        pair_labels = []
        for ph1, ph2 in combinations(phases, 2):
            key = (ph1, ph2)
            pair_labels.append(PAIR_SHORT.get(key, f"{ph1}-{ph2}"))

        band_tex = [BRAIN_BAND_TEX_DICT[b] for b in BANDS]
        matrix = np.full((len(BANDS), len(pair_labels)), np.nan)

        for i, band in enumerate(BANDS):
            for j, pl in enumerate(pair_labels):
                vals = dfp[(dfp["band"] == band) & (dfp["pair_label"] == pl)][
                    "mean_fraction_stable"
                ]
                if len(vals) > 0:
                    matrix[i, j] = vals.iloc[0]

        fig, ax = plt.subplots(figsize=(max(6, len(pair_labels) * 1.2), 5))
        im = ax.imshow(matrix, aspect="auto", cmap="RdYlGn", vmin=0.3, vmax=1.0)
        ax.set_xticks(range(len(pair_labels)))
        ax.set_xticklabels(pair_labels, rotation=45, ha="right", fontsize=9)
        ax.set_yticks(range(len(BANDS)))
        ax.set_yticklabels(band_tex, fontsize=11)

        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                val = matrix[i, j]
                if not np.isnan(val):
                    color = "white" if val < 0.5 else "black"
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                            fontsize=9, color=color)

        plt.colorbar(im, ax=ax, shrink=0.8)
        ax.set_title(f"{pat}: Mean Fraction Stable (k=2..10)", fontsize=12)
        fig.tight_layout()
        fig.savefig(out_dir / f"heatmap_{pat}.pdf", bbox_inches="tight", dpi=150)
        plt.close(fig)
        print(f"  Saved heatmap_{pat}.pdf")


def plot_all_patients_summary(df: pd.DataFrame, out_dir: Path):
    """Grand summary: band-level stability averaged over all patients and pairs."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel 1: Mean fraction stable by band (box plot across all patients/pairs)
    ax = axes[0]
    band_data = []
    for band in BANDS:
        vals = df[df["band"] == band]["mean_fraction_stable"].dropna().values
        band_data.append(vals)

    bp = ax.boxplot(band_data, labels=[BRAIN_BAND_TEX_DICT[b] for b in BANDS],
                    patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("lightblue")
    ax.set_ylabel("Mean Fraction Stable")
    ax.set_title("Stability by Band (all patients)")
    ax.axhline(0.5, color="gray", linestyle=":", alpha=0.5)

    # Panel 2: By pair type
    ax = axes[1]
    within = df[df["pair_type"] == "within"]["mean_fraction_stable"].dropna()
    cross = df[df["pair_type"] == "cross"]["mean_fraction_stable"].dropna()
    bp2 = ax.boxplot([within, cross], labels=["Within-type", "Cross-type"],
                     patch_artist=True)
    bp2["boxes"][0].set_facecolor("lightgreen")
    bp2["boxes"][1].set_facecolor("lightsalmon")
    ax.set_ylabel("Mean Fraction Stable")
    ax.set_title("Within vs Cross (all patients)")

    # Panel 3: By specific pair (4-phase patients)
    ax = axes[2]
    df4 = df[df["patient"].isin(PATIENTS_4PHASE)]
    pair_order = ["rest_pre-tL", "rest_pre-tT", "rest_pre-rsP", "tL-tT", "tL-rsP", "tT-rsP"]
    pair_data = []
    for pl in pair_order:
        vals = df4[df4["pair_label"] == pl]["mean_fraction_stable"].dropna().values
        pair_data.append(vals)
    bp3 = ax.boxplot(pair_data, labels=pair_order, patch_artist=True)
    pair_colors_map = {
        "rest_pre-tL": "lightsalmon", "rest_pre-tT": "lightsalmon",
        "rest_pre-rsP": "lightblue", "tL-tT": "lightgreen",
        "tL-rsP": "lightsalmon", "tT-rsP": "lightsalmon",
    }
    for patch, pl in zip(bp3["boxes"], pair_order):
        patch.set_facecolor(pair_colors_map.get(pl, "lightgray"))
    ax.set_xticklabels(pair_order, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Mean Fraction Stable")
    ax.set_title("Stability by Phase Pair")
    ax.axhline(0.5, color="gray", linestyle=":", alpha=0.5)

    fig.suptitle("Community Tracking: Grand Summary", fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(out_dir / "grand_summary.pdf", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved grand_summary.pdf")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=" * 70)
    print("Community Tracking Sweep")
    print("=" * 70)

    print("\n[1/6] Loading all linkage matrices...")
    data = load_all_linkages()
    print(f"  Loaded {len(data)} (patient, phase, band) entries")

    print("\n[2/6] Computing all metrics...")
    df = compute_all_metrics(data)
    print(f"  Computed {len(df)} pairwise comparisons")

    # Save raw results
    csv_path = OUT_DIR / "results.csv"
    df.to_csv(csv_path, index=False, float_format="%.4f")
    print(f"  Saved results.csv ({len(df)} rows)")

    # Quick summary to stdout
    print("\n--- Quick summary ---")
    print(f"  Patients: {sorted(df['patient'].unique())}")
    print(f"  Bands: {sorted(df['band'].unique())}")
    for metric in ["mean_fraction_stable", "fraction_stable_best_k",
                    "fraction_stable_natural"]:
        if metric in df.columns:
            print(f"  {metric}: mean={df[metric].mean():.3f}, "
                  f"std={df[metric].std():.3f}, "
                  f"range=[{df[metric].min():.3f}, {df[metric].max():.3f}]")

    # Within vs cross
    df4 = df[df["patient"].isin(PATIENTS_4PHASE)]
    within = df4[df4["pair_type"] == "within"]["mean_fraction_stable"]
    cross = df4[df4["pair_type"] == "cross"]["mean_fraction_stable"]
    print(f"\n  Within-type mean: {within.mean():.3f} +/- {within.std():.3f}")
    print(f"  Cross-type mean:  {cross.mean():.3f} +/- {cross.std():.3f}")

    # task_learn-task_test specifically
    tl_tt = df4[df4["pair_label"] == "tL-tT"]["mean_fraction_stable"]
    print(f"  task_learn-task_test mean: {tl_tt.mean():.3f} +/- {tl_tt.std():.3f}")

    # Check for bands where ALL pairs > 0.8
    print("\n  Bands where ALL pairs have mean_fraction_stable > 0.8:")
    for band in BANDS:
        vals = df4[df4["band"] == band]["mean_fraction_stable"]
        if len(vals) > 0 and vals.min() > 0.8:
            print(f"    {band}: min={vals.min():.3f}")
    print("  (none)" if not any(
        df4[df4["band"] == b]["mean_fraction_stable"].min() > 0.8
        for b in BANDS if len(df4[df4["band"] == b]) > 0
    ) else "")

    # Hypothesis testing
    print("\n[3/6] Testing hypotheses...")
    hyp_df = test_hypotheses(df)
    hyp_csv = OUT_DIR / "hypothesis_results.csv"
    hyp_df.to_csv(hyp_csv, index=False)
    print(f"  Saved hypothesis_results.csv")

    # Print summary
    for metric in hyp_df["metric"].unique():
        mdf = hyp_df[hyp_df["metric"] == metric]
        print(f"\n  {metric}:")
        # Count how many bands have H3 (within>cross) for majority of patients
        h3_vals = mdf["H3_within_gt_cross"].values
        n_strong = sum(1 for v in h3_vals if "/" in v and
                       int(v.split("/")[0]) >= int(v.split("/")[1]) * 0.75)
        print(f"    H3 (within>cross) strong in {n_strong}/{len(h3_vals)} bands")

    print("\n[4/6] Generating plots...")
    plot_summary_heatmap(df, OUT_DIR)
    plot_k_sweep(df, OUT_DIR)
    plot_node_tracking(df, data, OUT_DIR)
    plot_consistency_table(hyp_df, OUT_DIR)

    print("\n[5/6] Per-patient heatmaps...")
    plot_per_patient_heatmaps(df, OUT_DIR)

    print("\n[6/6] Grand summary...")
    plot_all_patients_summary(df, OUT_DIR)

    print("\n" + "=" * 70)
    print("DONE. Output in:", OUT_DIR)
    print("=" * 70)


if __name__ == "__main__":
    main()
