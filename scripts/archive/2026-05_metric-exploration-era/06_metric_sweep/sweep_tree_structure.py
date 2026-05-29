#!/usr/bin/env python3
"""Tree structure comparison metrics for LRG dendrograms.

Metrics based on tree TOPOLOGY (bipartitions, merge order) rather than
cophenetic distances or flat partitions.

Output: data/figures/metric_exploration/tree_structure/
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import to_tree, fcluster
from scipy.stats import spearmanr

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import LRG_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
from lrg_eegfc.config.const import PATIENTS_4PHASE
PATIENTS_4PH = PATIENTS_4PHASE
BANDS = BRAIN_BANDS_NAMES
OUT_DIR = FIGURES_ROOT / "metric_exploration" / "tree_structure"

ALL_PAIRS = [
    ("rest_pre", "rest_post"), ("task_learn", "task_test"),
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"),
    ("task_learn", "rest_post"), ("task_test", "rest_post"),
]

PAIR_SHORT = {
    ("rest_pre", "rest_post"): "Pre↔Post", ("task_learn", "task_test"): "TL↔TT",
    ("rest_pre", "task_learn"): "Pre↔TL", ("rest_pre", "task_test"): "Pre↔TT",
    ("task_learn", "rest_post"): "TL↔Post", ("task_test", "rest_post"): "TT↔Post",
}

PAIR_CAT = {
    ("rest_pre", "rest_post"): "within", ("task_learn", "task_test"): "within",
    ("rest_pre", "task_learn"): "cross", ("rest_pre", "task_test"): "cross",
    ("task_learn", "rest_post"): "cross", ("task_test", "rest_post"): "cross",
}

def bl(b):
    return BRAIN_BAND_TEX_DICT.get(b, b)


# ===================================================================
# Extract bipartitions from scipy linkage
# ===================================================================
def get_bipartitions(Z, n):
    """Extract all bipartitions (as frozensets of leaf ids) from linkage."""
    root = to_tree(Z)
    bips = {}  # frozenset → merge_height

    def traverse(node):
        if node.is_leaf():
            return frozenset([node.id])
        left_leaves = traverse(node.get_left())
        right_leaves = traverse(node.get_right())
        smaller = min(left_leaves, right_leaves, key=len)
        bips[smaller] = node.dist
        return left_leaves | right_leaves

    traverse(root)
    return bips


def get_merge_order(Z, n):
    """For each pair of leaves (i, j), get the merge step at which they first
    become co-clustered. Returns condensed vector of length n*(n-1)/2."""
    # Build cluster memberships at each step
    from scipy.spatial.distance import squareform

    # Use the linkage merge heights as cophenetic distances
    # The merge order is equivalent to ranking cophenetic distances
    from scipy.cluster.hierarchy import cophenet
    coph = cophenet(Z)
    # coph returns (condensed_coph_dists, coph_matrix) — we want condensed
    if isinstance(coph, tuple):
        coph_dists = coph[0]
    else:
        coph_dists = coph
    # Rank them (merge order)
    from scipy.stats import rankdata
    return rankdata(coph_dists)


# ===================================================================
# Metrics
# ===================================================================
def robinson_foulds_sim(Z1, Z2, n):
    """Normalized Robinson-Foulds similarity: fraction of shared bipartitions."""
    bips1 = set(get_bipartitions(Z1, n).keys())
    bips2 = set(get_bipartitions(Z2, n).keys())
    if not bips1 and not bips2:
        return 1.0
    shared = len(bips1 & bips2)
    total = max(len(bips1), len(bips2))
    return shared / total if total > 0 else 1.0


def weighted_rf_sim(Z1, Z2, n):
    """Weighted Robinson-Foulds: shared bipartitions weighted by merge height."""
    bips1 = get_bipartitions(Z1, n)
    bips2 = get_bipartitions(Z2, n)
    shared_keys = set(bips1.keys()) & set(bips2.keys())
    if not shared_keys:
        return 0.0
    # Weight = mean merge height of the bipartition in both trees
    weights = [0.5 * (bips1[k] + bips2[k]) for k in shared_keys]
    total_weight = sum(0.5 * (bips1[k] + bips2[k]) for k in set(bips1) | set(bips2))
    return sum(weights) / total_weight if total_weight > 0 else 0.0


def branch_score_sim(Z1, Z2, n):
    """Branch score: for shared bipartitions, compare their merge heights.
    Similarity = 1 - mean|h1-h2|/max_h."""
    bips1 = get_bipartitions(Z1, n)
    bips2 = get_bipartitions(Z2, n)
    shared = set(bips1.keys()) & set(bips2.keys())
    if not shared:
        return 0.0
    max_h = max(max(bips1.values()), max(bips2.values()))
    if max_h < 1e-12:
        return 1.0
    diffs = [abs(bips1[k] - bips2[k]) / max_h for k in shared]
    return 1.0 - np.mean(diffs)


def merge_order_corr(Z1, Z2, n):
    """Spearman correlation of merge orders (when do leaf pairs first co-cluster)."""
    ord1 = get_merge_order(Z1, n)
    ord2 = get_merge_order(Z2, n)
    if len(ord1) != len(ord2):
        return np.nan
    return spearmanr(ord1, ord2).correlation


def height_distribution_sim(Z1, Z2, n):
    """Compare distributions of merge heights using 1 - KS statistic."""
    from scipy.stats import ks_2samp
    h1 = np.sort(Z1[:, 2])
    h2 = np.sort(Z2[:, 2])
    ks_stat = ks_2samp(h1, h2).statistic
    return 1.0 - ks_stat


METRICS = {
    "RF_sim": robinson_foulds_sim,
    "weighted_RF": weighted_rf_sim,
    "branch_score": branch_score_sim,
    "merge_order": merge_order_corr,
    "height_dist": height_distribution_sim,
}


# ===================================================================
# Load data
# ===================================================================
def load_data():
    data = {}
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for ph in phases:
                lrg = load_lrg_result(pat, ph, band, "msc", cache_root=LRG_CACHE)
                if lrg is None:
                    continue
                data[(pat, ph, band)] = {
                    "Z": lrg.linkage_matrix,
                    "n": lrg.n_nodes,
                }
    return data


# ===================================================================
# Compute all metrics
# ===================================================================
def compute_all(data):
    rows = []
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for p1, p2 in combinations(phases, 2):
                k1, k2 = (pat, p1, band), (pat, p2, band)
                if k1 not in data or k2 not in data:
                    continue
                d1, d2 = data[k1], data[k2]
                if d1["n"] != d2["n"]:
                    continue
                n = d1["n"]
                row = {"patient": pat, "band": band, "phase1": p1, "phase2": p2,
                       "pair_type": PAIR_CAT.get((p1, p2), "?")}
                for mname, mfunc in METRICS.items():
                    try:
                        row[mname] = mfunc(d1["Z"], d2["Z"], n)
                    except Exception as e:
                        row[mname] = np.nan
                rows.append(row)
    return pd.DataFrame(rows)


# ===================================================================
# Analysis
# ===================================================================
def analyze_hypotheses(df):
    """Test H1-H4 for each metric."""
    results = {}
    for mname in METRICS:
        h1_unan, h1_relax = 0, 0  # task persistence
        h2_unan, h2_relax = 0, 0  # task trace
        h3_unan, h3_relax = 0, 0  # within > cross
        for band in BANDS:
            bd = df[(df["band"] == band) & (df["patient"].isin(PATIENTS_4PH))]
            if bd.empty:
                continue
            # Z-score within each patient
            h1_votes, h2_votes, h3_votes = [], [], []
            for pat in PATIENTS_4PH:
                pd_pat = bd[bd["patient"] == pat]
                if len(pd_pat) < 6:
                    continue
                vals = pd_pat[mname].values
                mu, sigma = vals.mean(), vals.std()
                if sigma < 1e-10:
                    continue
                zs = {(r["phase1"], r["phase2"]): (r[mname] - mu) / sigma
                      for _, r in pd_pat.iterrows()}

                # H1: task_learn-task_test has highest z
                tt_z = zs.get(("task_learn", "task_test"), None)
                if tt_z is not None:
                    other_zs = [v for k, v in zs.items() if k != ("task_learn", "task_test")]
                    h1_votes.append(tt_z > max(other_zs) if other_zs else True)

                # H2: task_test-rest_post z > rest_pre-rest_post z
                tt_post = zs.get(("task_test", "rest_post"), None)
                pre_post = zs.get(("rest_pre", "rest_post"), None)
                if tt_post is not None and pre_post is not None:
                    h2_votes.append(tt_post > pre_post)

                # H3: within z > cross z
                within = [v for k, v in zs.items()
                          if k in [("rest_pre", "rest_post"), ("task_learn", "task_test")]]
                cross = [v for k, v in zs.items()
                         if k not in [("rest_pre", "rest_post"), ("task_learn", "task_test")]]
                if within and cross:
                    h3_votes.append(np.mean(within) > np.mean(cross))

            if h1_votes:
                if all(h1_votes): h1_unan += 1
                if sum(h1_votes) >= len(h1_votes) - 1: h1_relax += 1
            if h2_votes:
                if all(h2_votes): h2_unan += 1
                if sum(h2_votes) >= len(h2_votes) - 1: h2_relax += 1
            if h3_votes:
                if all(h3_votes): h3_unan += 1
                if sum(h3_votes) >= len(h3_votes) - 1: h3_relax += 1

        results[mname] = {
            "H1_unan": h1_unan, "H1_relax": h1_relax,
            "H2_unan": h2_unan, "H2_relax": h2_relax,
            "H3_unan": h3_unan, "H3_relax": h3_relax,
            "score": h1_unan + h2_unan + h3_unan,
        }
    return results


# ===================================================================
# Figures
# ===================================================================
def plot_summary(hyp_results, df, out_dir):
    """Summary heatmap of hypothesis results."""
    metrics = list(hyp_results.keys())
    fig, ax = plt.subplots(figsize=(10, 5))
    mat = np.zeros((len(metrics), 6))
    cols = ["H1_unan", "H1_relax", "H2_unan", "H2_relax", "H3_unan", "H3_relax"]
    for i, m in enumerate(metrics):
        for j, c in enumerate(cols):
            mat[i, j] = hyp_results[m][c]

    im = ax.imshow(mat, cmap='YlGn', aspect='auto', vmin=0, vmax=6)
    for i in range(len(metrics)):
        for j in range(6):
            ax.text(j, i, f"{int(mat[i,j])}/6", ha='center', va='center', fontsize=9)
    ax.set_xticks(range(6))
    ax.set_xticklabels(["H1\nunan", "H1\nrelax", "H2\nunan", "H2\nrelax",
                         "H3\nunan", "H3\nrelax"], fontsize=9)
    ax.set_yticks(range(len(metrics)))
    ax.set_yticklabels(metrics, fontsize=10)
    plt.colorbar(im, ax=ax, shrink=0.8, label="# bands passing")
    ax.set_title("Tree structure metrics: hypothesis testing\n"
                 "H1=task persist, H2=task trace, H3=within>cross")
    fig.tight_layout()
    fig.savefig(out_dir / "summary_heatmap.pdf", bbox_inches='tight')
    plt.close(fig)

    # Per-pair profiles for each metric
    fig2, axes = plt.subplots(2, 3, figsize=(18, 10))
    for m_idx, mname in enumerate(metrics):
        if m_idx >= 6:
            break
        ax = axes[m_idx // 3, m_idx % 3]
        for pat in PATIENTS_4PH:
            for pair in ALL_PAIRS:
                p1, p2 = pair
                row = df[(df["patient"] == pat) & (df["phase1"] == p1) & (df["phase2"] == p2)]
                vals = [row[row["band"] == b][mname].values[0] if len(row[row["band"] == b]) > 0
                        else np.nan for b in BANDS]
                cat = PAIR_CAT[pair]
                ls = '-' if cat == "within" else '--'
                lw = 1.5 if cat == "within" else 0.8
                ax.plot(range(len(BANDS)), vals, f'o{ls}', markersize=3,
                        linewidth=lw, alpha=0.5)

        # Mean per pair
        for pair in ALL_PAIRS:
            p1, p2 = pair
            means = []
            for b in BANDS:
                rows = df[(df["patient"].isin(PATIENTS_4PH)) &
                          (df["phase1"] == p1) & (df["phase2"] == p2) &
                          (df["band"] == b)]
                means.append(rows[mname].mean() if len(rows) > 0 else np.nan)
            cat = PAIR_CAT[pair]
            ls = '-' if cat == "within" else '--'
            lw = 2.5 if cat == "within" else 1.5
            ax.plot(range(len(BANDS)), means, f's{ls}', markersize=6,
                    linewidth=lw, label=PAIR_SHORT[pair])

        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([bl(b) for b in BANDS], fontsize=7, rotation=30)
        ax.set_title(mname, fontsize=10, fontweight='bold')
        ax.legend(fontsize=6, ncol=2)

    for idx in range(len(metrics), 6):
        axes[idx // 3, idx % 3].set_visible(False)

    fig2.suptitle("Tree structure metrics: per-pair profiles\n"
                  "(solid=within, dashed=cross, squares=mean)", fontsize=12)
    fig2.tight_layout(rect=[0, 0, 1, 0.93])
    fig2.savefig(out_dir / "pair_profiles.pdf", bbox_inches='tight')
    plt.close(fig2)

    # Within-cross gap heatmap
    fig3, axes3 = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 5))
    if len(metrics) == 1:
        axes3 = [axes3]
    for m_idx, mname in enumerate(metrics):
        ax = axes3[m_idx]
        gap_mat = np.full((len(BANDS), len(PATIENTS_4PH)), np.nan)
        for b_idx, band in enumerate(BANDS):
            for p_idx, pat in enumerate(PATIENTS_4PH):
                pd_pat = df[(df["patient"] == pat) & (df["band"] == band)]
                if len(pd_pat) < 6:
                    continue
                within_vals = pd_pat[pd_pat["pair_type"] == "within"][mname].values
                cross_vals = pd_pat[pd_pat["pair_type"] == "cross"][mname].values
                if len(within_vals) > 0 and len(cross_vals) > 0:
                    gap_mat[b_idx, p_idx] = np.mean(within_vals) - np.mean(cross_vals)

        vmax = max(0.05, np.nanmax(np.abs(gap_mat)))
        im = ax.imshow(gap_mat, cmap='RdBu', vmin=-vmax, vmax=vmax, aspect='auto')
        for i in range(len(BANDS)):
            for j in range(len(PATIENTS_4PH)):
                v = gap_mat[i, j]
                if not np.isnan(v):
                    ax.text(j, i, f"{v:+.3f}", ha='center', va='center', fontsize=7)
        ax.set_xticks(range(len(PATIENTS_4PH)))
        ax.set_xticklabels([p.replace("Pat_0", "P") for p in PATIENTS_4PH], fontsize=8)
        ax.set_yticks(range(len(BANDS)))
        ax.set_yticklabels([bl(b) for b in BANDS] if m_idx == 0 else [], fontsize=9)
        ax.set_title(mname, fontsize=9)

    fig3.suptitle("Within−Cross gap per (patient, band)", fontsize=12)
    fig3.tight_layout(rect=[0, 0, 1, 0.93])
    fig3.savefig(out_dir / "within_cross_gap.pdf", bbox_inches='tight')
    plt.close(fig3)


# ===================================================================
# Main
# ===================================================================
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Tree Structure Comparison Sweep")
    print("=" * 70)

    print("\n[1/4] Loading data...")
    data = load_data()
    print(f"  Loaded {len(data)} entries")

    print("\n[2/4] Computing metrics...")
    df = compute_all(data)
    print(f"  Computed {len(df)} comparisons")
    df.to_csv(OUT_DIR / "results.csv", index=False)

    print("\n[3/4] Testing hypotheses...")
    hyp = analyze_hypotheses(df)

    print("\n  Metric         H1u H1r H2u H2r H3u H3r  Score")
    print("  " + "-" * 55)
    for m in sorted(hyp, key=lambda x: -hyp[x]["score"]):
        h = hyp[m]
        print(f"  {m:<16} {h['H1_unan']:>3} {h['H1_relax']:>3} "
              f"{h['H2_unan']:>3} {h['H2_relax']:>3} "
              f"{h['H3_unan']:>3} {h['H3_relax']:>3}  {h['score']:>5}")

    # Within-cross gap detail for best metric
    best = max(hyp, key=lambda x: hyp[x]["H3_unan"])
    print(f"\n  Best H3 metric: {best}")
    print(f"  Within-cross gap per (patient, band):")
    for band in BANDS:
        print(f"    {bl(band):<14}", end="")
        all_pos = True
        for pat in PATIENTS_4PH:
            pd_pat = df[(df["patient"] == pat) & (df["band"] == band)]
            if len(pd_pat) < 6:
                print(f"  {'N/A':>8}", end="")
                continue
            wv = pd_pat[pd_pat["pair_type"] == "within"][best].values
            cv = pd_pat[pd_pat["pair_type"] == "cross"][best].values
            gap = np.mean(wv) - np.mean(cv)
            if gap <= 0:
                all_pos = False
            print(f"  {gap:>+8.4f}", end="")
        print(f"  {'ALL+' if all_pos else 'MIXED'}")

    print("\n[4/4] Generating figures...")
    plot_summary(hyp, df, OUT_DIR)

    print("\n" + "=" * 70)
    print("DONE. Output in:", OUT_DIR)
    print("=" * 70)


if __name__ == "__main__":
    main()
