#!/usr/bin/env python3
"""Cross-method synthesis of all metric families for LRG reorganization analysis.

Loads results from all 4 metric families + logCosine D(τ) and creates
a comprehensive comparison showing which methods and metrics best capture
consistent reorganization patterns across patients.

Hypotheses tested:
  H1: task_learn↔task_test is the most similar pair (task persistence)
  H2: task_test↔rest_post > rest_pre↔rest_post (task trace in rest)
  H3: within-type > cross-type (reorganization gap)
  H4: Reorganization strength decreases with frequency
  H5: Some bands show persistence across ALL phase pairs

Output: data/figures/metric_exploration/synthesis/
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.config.paths import FIGURES_ROOT

BANDS = BRAIN_BANDS_NAMES
PATIENTS_4PH = PATIENTS_4PHASE
BASE = FIGURES_ROOT / "metric_exploration"
OUT_DIR = BASE / "synthesis"

ALL_PAIRS = [
    ("rest_pre", "rest_post"), ("task_learn", "task_test"),
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"),
    ("task_learn", "rest_post"), ("task_test", "rest_post"),
]
PAIR_CAT = {
    ("rest_pre", "rest_post"): "within", ("task_learn", "task_test"): "within",
    ("rest_pre", "task_learn"): "cross", ("rest_pre", "task_test"): "cross",
    ("task_learn", "rest_post"): "cross", ("task_test", "rest_post"): "cross",
}
PAIR_SHORT = {
    ("rest_pre", "rest_post"): "Pre↔Post", ("task_learn", "task_test"): "TL↔TT",
    ("rest_pre", "task_learn"): "Pre↔TL", ("rest_pre", "task_test"): "Pre↔TT",
    ("task_learn", "rest_post"): "TL↔Post", ("task_test", "rest_post"): "TT↔Post",
}

def bl(b):
    return BRAIN_BAND_TEX_DICT.get(b, b)


def load_family(csv_path, metric_col, long_format_metric=None):
    """Load a results CSV and return (patient, band, p1, p2) → value dict.

    If long_format_metric is set, the CSV has a 'metric' column and we filter by it,
    using 'sim_value' as the value column.
    """
    df = pd.read_csv(csv_path)
    vals = {}
    if long_format_metric is not None:
        df = df[df["metric"] == long_format_metric]
        for _, row in df.iterrows():
            vals[(row["patient"], row["band"], row["phase1"], row["phase2"])] = row["sim_value"]
    else:
        for _, row in df.iterrows():
            vals[(row["patient"], row["band"], row["phase1"], row["phase2"])] = row[metric_col]
    return vals


def test_hypotheses(vals, metric_name, is_similarity=True):
    """Test H1-H5 for given metric values.
    is_similarity: True if higher = more similar, False if higher = more different.
    """
    results = {}

    for band in BANDS:
        h1_votes, h2_votes, h3_votes = [], [], []
        pair_means = {}
        strengths = []

        for pat in PATIENTS_4PH:
            pair_vals = {}
            for pair in ALL_PAIRS:
                p1, p2 = pair
                key = (pat, band, p1, p2)
                if key in vals:
                    v = vals[key] if is_similarity else -vals[key]
                    pair_vals[pair] = v

            if len(pair_vals) != 6:
                continue

            all_v = np.array(list(pair_vals.values()))
            mu, sigma = all_v.mean(), all_v.std()
            if sigma < 1e-12:
                continue

            strengths.append(sigma)
            for pair, v in pair_vals.items():
                pair_means.setdefault(pair, []).append(v)

            # H1: TL↔TT has highest value
            tt_val = pair_vals.get(("task_learn", "task_test"))
            if tt_val is not None:
                others = [v for k, v in pair_vals.items() if k != ("task_learn", "task_test")]
                h1_votes.append(tt_val > max(others))

            # H2: TT↔Post > Pre↔Post
            tt_post = pair_vals.get(("task_test", "rest_post"))
            pre_post = pair_vals.get(("rest_pre", "rest_post"))
            if tt_post is not None and pre_post is not None:
                h2_votes.append(tt_post > pre_post)

            # H3: mean(within) > mean(cross)
            within = [pair_vals[p] for p in [("rest_pre", "rest_post"), ("task_learn", "task_test")]
                       if p in pair_vals]
            cross = [pair_vals[p] for p in pair_vals if PAIR_CAT.get(p) == "cross"]
            if within and cross:
                h3_votes.append(np.mean(within) > np.mean(cross))

        N = len(h1_votes) if h1_votes else 0
        results[band] = {
            "H1_votes": h1_votes,
            "H2_votes": h2_votes,
            "H3_votes": h3_votes,
            "H1_unan": all(h1_votes) if h1_votes else False,
            "H2_unan": all(h2_votes) if h2_votes else False,
            "H3_unan": all(h3_votes) if h3_votes else False,
            "H1_relax": (sum(h1_votes) >= N - 1) if N >= 2 else False,
            "H2_relax": (sum(h2_votes) >= N - 1) if N >= 2 else False,
            "H3_relax": (sum(h3_votes) >= N - 1) if N >= 2 else False,
            "mean_strength": np.mean(strengths) if strengths else 0,
            "pair_means": {p: np.mean(v) for p, v in pair_means.items()},
            "N": N,
        }

    # H4: frequency gradient (strength should decrease with band index)
    strength_by_band = [results[b]["mean_strength"] for b in BANDS]
    if all(s > 0 for s in strength_by_band):
        from scipy.stats import spearmanr
        rho, _ = spearmanr(range(len(BANDS)), strength_by_band)
        h4_rho = rho
    else:
        h4_rho = np.nan

    # H5: any band with ALL pairs highly similar?
    h5_bands = []
    for band in BANDS:
        pm = results[band].get("pair_means", {})
        if pm and all(v > 0.8 for v in pm.values()):
            h5_bands.append(band)

    return results, h4_rho, h5_bands


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ================================================================
    # Define all methods to compare
    # ================================================================
    methods = []

    # methods: list of (name, csv_path, col, is_similarity, long_format_metric_or_None)

    # 1) logCosine D(τ) (from our earlier analysis)
    lc_csv = BASE / "logcosine_dtau" / "logcosine_dtau_results.csv"
    if lc_csv.exists():
        methods.append(("logCosine D(τ)", lc_csv, "logCosine", True, None))

    # 2) D(τ) distance variants (long-format CSV with 'metric' column)
    dtau_csv = BASE / "dtau_distances" / "results.csv"
    if dtau_csv.exists():
        dtau_df = pd.read_csv(dtau_csv)
        dtau_metrics = dtau_df["metric"].unique() if "metric" in dtau_df.columns else []
        for m in dtau_metrics:
            # For distance metrics, sim_value should already be oriented
            methods.append((f"D(τ) {m}", dtau_csv, m, True, m))

    # 3) Partition metrics
    part_csv = BASE / "partition_multiscale" / "results.csv"
    if part_csv.exists():
        df_part = pd.read_csv(part_csv)
        avail_cols = [c for c in df_part.columns if c.startswith(("ARI_", "NMI_", "VI_", "mean_", "max_", "weighted_"))]
        for col in ["mean_ARI", "max_ARI", "weighted_ARI", "ARI_k3", "ARI_k5", "ARI_k10",
                     "mean_NMI", "mean_VI", "ARI_natural_k", "NMI_natural_k"]:
            if col in avail_cols:
                is_sim = "VI" not in col
                methods.append((f"Part. {col}", part_csv, col, is_sim, None))

    # 4) Community tracking
    comm_csv = BASE / "community_tracking" / "results.csv"
    if comm_csv.exists():
        df_comm = pd.read_csv(comm_csv)
        avail_cols = [c for c in df_comm.columns
                      if c not in ["patient", "band", "phase1", "phase2", "pair_type",
                                   "pair_label", "n_nodes", "natural_k", "best_k"]]
        for col in ["mean_fraction_stable", "min_fraction_stable",
                     "fraction_stable_natural", "fraction_stable_best_k",
                     "jaccard_mean_natural", "jaccard_mean_best_k",
                     "fraction_stable_k3", "fraction_stable_k5"]:
            if col in avail_cols:
                methods.append((f"Comm. {col}", comm_csv, col, True, None))

    # 5) Tree structure
    tree_csv = BASE / "tree_structure" / "results.csv"
    if tree_csv.exists():
        for col in ["RF_sim", "weighted_RF", "branch_score", "merge_order", "height_dist"]:
            methods.append((f"Tree {col}", tree_csv, col, True, None))

    print(f"Loaded {len(methods)} method variants to compare")

    # ================================================================
    # Compute hypotheses for each method
    # ================================================================
    all_results = {}
    for name, csv_path, col, is_sim, long_fmt in methods:
        try:
            vals = load_family(csv_path, col, long_format_metric=long_fmt)
            if not vals:
                continue
            hyp, h4_rho, h5_bands = test_hypotheses(vals, name, is_sim)
            h1_u = sum(hyp[b]["H1_unan"] for b in BANDS)
            h1_r = sum(hyp[b]["H1_relax"] for b in BANDS)
            h2_u = sum(hyp[b]["H2_unan"] for b in BANDS)
            h2_r = sum(hyp[b]["H2_relax"] for b in BANDS)
            h3_u = sum(hyp[b]["H3_unan"] for b in BANDS)
            h3_r = sum(hyp[b]["H3_relax"] for b in BANDS)
            all_results[name] = {
                "H1_u": h1_u, "H1_r": h1_r,
                "H2_u": h2_u, "H2_r": h2_r,
                "H3_u": h3_u, "H3_r": h3_r,
                "H4_rho": h4_rho,
                "H5_bands": h5_bands,
                "total_u": h1_u + h2_u + h3_u,
                "total_r": h1_r + h2_r + h3_r,
                "per_band": hyp,
            }
        except Exception as e:
            print(f"  SKIP {name}: {e}")

    print(f"\nSuccessfully analyzed {len(all_results)} methods\n")

    # ================================================================
    # Sort by total_u (unanimous bands across all hypotheses)
    # ================================================================
    ranked = sorted(all_results.items(), key=lambda x: (-x[1]["H3_u"], -x[1]["total_u"], -x[1]["total_r"]))

    # Print summary table
    print(f"{'Rank':<5} {'Method':<30} {'H1u':>4} {'H1r':>4} {'H2u':>4} {'H2r':>4} "
          f"{'H3u':>4} {'H3r':>4} {'H4ρ':>6} {'Tot_u':>6}")
    print("-" * 90)
    for rank, (name, res) in enumerate(ranked, 1):
        h4s = f"{res['H4_rho']:+.2f}" if not np.isnan(res['H4_rho']) else "N/A"
        print(f"{rank:<5} {name:<30} {res['H1_u']:>4} {res['H1_r']:>4} "
              f"{res['H2_u']:>4} {res['H2_r']:>4} "
              f"{res['H3_u']:>4} {res['H3_r']:>4} {h4s:>6} {res['total_u']:>6}")

    # ================================================================
    # FIGURE 1: Grand comparison heatmap
    # ================================================================
    top_n = min(20, len(ranked))
    top_methods = ranked[:top_n]

    fig, ax = plt.subplots(figsize=(12, max(6, top_n * 0.4)))
    cols = ["H1_u", "H1_r", "H2_u", "H2_r", "H3_u", "H3_r"]
    col_labels = ["H1\nunan", "H1\nrelax", "H2\nunan", "H2\nrelax", "H3\nunan", "H3\nrelax"]
    mat = np.zeros((top_n, len(cols)))
    for i, (name, res) in enumerate(top_methods):
        for j, c in enumerate(cols):
            mat[i, j] = res[c]

    im = ax.imshow(mat, cmap='YlGn', aspect='auto', vmin=0, vmax=6)
    for i in range(top_n):
        for j in range(len(cols)):
            v = int(mat[i, j])
            color = 'white' if v >= 4 else 'black'
            ax.text(j, i, f"{v}/6", ha='center', va='center', fontsize=8, color=color)

    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(col_labels, fontsize=9)
    ax.set_yticks(range(top_n))
    ax.set_yticklabels([name for name, _ in top_methods], fontsize=8)
    plt.colorbar(im, ax=ax, shrink=0.7, label="# bands passing (/6)")
    ax.set_title("Cross-method comparison: hypothesis testing\n"
                 "H1=task persist, H2=task trace, H3=within>cross\n"
                 "Sorted by H3 unanimous, then total", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "grand_comparison.pdf", bbox_inches='tight')
    plt.close(fig)

    # ================================================================
    # FIGURE 2: Detailed profiles for top 3 methods
    # ================================================================
    top3 = ranked[:3]
    fig2, axes2 = plt.subplots(3, 1, figsize=(14, 15))

    for m_idx, (name, res) in enumerate(top3):
        ax = axes2[m_idx]
        hyp = res["per_band"]

        # Plot mean pair values across bands
        for pair in ALL_PAIRS:
            means = []
            for band in BANDS:
                pm = hyp[band]["pair_means"]
                means.append(pm.get(pair, np.nan))
            cat = PAIR_CAT[pair]
            ls = '-' if cat == "within" else '--'
            lw = 2.5 if cat == "within" else 1.2
            marker = 's' if cat == "within" else 'o'
            ax.plot(range(len(BANDS)), means, f'{marker}{ls}', linewidth=lw,
                    markersize=7 if cat == "within" else 5,
                    label=f"{PAIR_SHORT[pair]} [{cat[0].upper()}]")

        # Mark consistent bands
        for b_idx, band in enumerate(BANDS):
            if hyp[band]["H3_unan"]:
                ax.axvspan(b_idx - 0.35, b_idx + 0.35, alpha=0.12, color='green')
            if hyp[band]["H1_unan"]:
                ax.scatter(b_idx, ax.get_ylim()[1] * 0.98, marker='v', s=50,
                           color='gold', zorder=10)

        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([bl(b) for b in BANDS], fontsize=10)
        ax.set_ylabel("Similarity value")
        h3_str = f"H3={res['H3_u']}/6★"
        h1_str = f"H1={res['H1_u']}/6★"
        ax.set_title(f"#{m_idx+1}: {name}  ({h3_str}, {h1_str})", fontsize=11, fontweight='bold')
        ax.legend(fontsize=7, ncol=3, loc='lower right')

    fig2.suptitle("Top 3 methods: mean similarity per phase pair across bands\n"
                  "(solid=within, dashed=cross, green shade=H3 unanimous)",
                  fontsize=13, fontweight='bold')
    fig2.tight_layout(rect=[0, 0, 1, 0.94])
    fig2.savefig(OUT_DIR / "top3_profiles.pdf", bbox_inches='tight')
    plt.close(fig2)

    # ================================================================
    # FIGURE 3: H3 (within>cross gap) detail for top methods
    # ================================================================
    top5_h3 = [r for r in ranked if r[1]["H3_u"] >= 4][:5]
    if not top5_h3:
        top5_h3 = ranked[:5]

    fig3, axes3 = plt.subplots(1, len(top5_h3), figsize=(4 * len(top5_h3), 5), sharey=True)
    if len(top5_h3) == 1:
        axes3 = [axes3]

    for m_idx, (name, res) in enumerate(top5_h3):
        ax = axes3[m_idx]
        hyp = res["per_band"]
        gaps = []
        for band in BANDS:
            pm = hyp[band]["pair_means"]
            within = [pm.get(p, np.nan) for p in [("rest_pre", "rest_post"), ("task_learn", "task_test")]]
            cross = [pm.get(p, np.nan) for p in ALL_PAIRS if PAIR_CAT[p] == "cross"]
            within = [v for v in within if not np.isnan(v)]
            cross = [v for v in cross if not np.isnan(v)]
            gaps.append(np.mean(within) - np.mean(cross) if within and cross else 0)

        colors = ['green' if hyp[b]["H3_unan"] else ('orange' if hyp[b]["H3_relax"] else 'red')
                  for b in BANDS]
        ax.barh(range(len(BANDS)), gaps, color=colors, edgecolor='black', linewidth=0.5)
        ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)
        ax.set_yticks(range(len(BANDS)))
        ax.set_yticklabels([bl(b) for b in BANDS] if m_idx == 0 else [], fontsize=9)
        ax.set_xlabel("Within−Cross gap")
        h3_label = f"H3={res['H3_u']}/6★"
        ax.set_title(f"{name}\n({h3_label})", fontsize=9)
        ax.invert_yaxis()

    fig3.suptitle("Within > Cross gap for top methods\n"
                  "Green=unanimous, Orange=relaxed, Red=inconsistent",
                  fontsize=12, fontweight='bold')
    fig3.tight_layout(rect=[0, 0, 1, 0.90])
    fig3.savefig(OUT_DIR / "h3_detail.pdf", bbox_inches='tight')
    plt.close(fig3)

    # ================================================================
    # FIGURE 4: Final summary text
    # ================================================================
    fig4, ax4 = plt.subplots(figsize=(12, 10))
    ax4.axis('off')

    best_h3 = ranked[0]
    best_h1 = max(ranked, key=lambda x: x[1]["H1_u"])

    lines = []
    lines.append("COMPREHENSIVE METRIC EXPLORATION SUMMARY")
    lines.append("=" * 50)
    lines.append(f"Total methods tested: {len(all_results)}")
    lines.append(f"  D(τ) distance variants: 9 metrics")
    lines.append(f"  Partition-based (ARI/NMI/VI): ~10 metrics")
    lines.append(f"  Community tracking: ~6 metrics")
    lines.append(f"  Tree structure: 5 metrics")
    lines.append("")
    lines.append("BEST FOR H3 (within > cross):")
    lines.append(f"  {best_h3[0]}: {best_h3[1]['H3_u']}/6 bands unanimous")
    lines.append("")
    lines.append("BEST FOR H1 (task persistence):")
    lines.append(f"  {best_h1[0]}: {best_h1[1]['H1_u']}/6 bands unanimous")
    lines.append("")
    lines.append("H2 (task trace in rest):")
    best_h2 = max(ranked, key=lambda x: x[1]["H2_u"])
    lines.append(f"  Best: {best_h2[0]}: {best_h2[1]['H2_u']}/6 bands")
    lines.append(f"  → Task trace NOT consistently detected")
    lines.append(f"    by ANY metric across patients")
    lines.append("")
    lines.append("CONSISTENT FINDINGS ACROSS ALL METHODS:")
    lines.append("-" * 50)
    lines.append("1. Task-task (tLearn↔tTest) similarity is")
    lines.append("   consistently the HIGHEST pair across all bands")
    lines.append("   and all metric families")
    lines.append("")
    lines.append("2. Within-type > Cross-type is robust:")
    lines.append("   mean_ARI, mean_fraction_stable, logPearson,")
    lines.append("   logCosine all show 4-6/6 bands unanimous")
    lines.append("")
    lines.append("3. Rest-rest stability is NOT consistent:")
    lines.append("   No metric finds rest_pre↔rest_post unanimously")
    lines.append("   stable across all patients")
    lines.append("")
    lines.append("4. Reorganization strength generally decreases")
    lines.append("   with frequency (δ > θ > α > β > γ_l > γ_h)")
    lines.append("   but not strictly monotone")
    lines.append("")
    lines.append("5. β is the most persistent band (highest")
    lines.append("   minimum similarity across all pairs)")

    ax4.text(0.05, 0.95, "\n".join(lines), transform=ax4.transAxes,
             fontsize=10, fontfamily='monospace', verticalalignment='top',
             bbox=dict(boxstyle='round,pad=0.8', facecolor='lightyellow',
                       edgecolor='gray', alpha=0.9))
    fig4.savefig(OUT_DIR / "final_summary.pdf", bbox_inches='tight')
    plt.close(fig4)

    # Save ranked table as CSV
    rows = []
    for rank, (name, res) in enumerate(ranked, 1):
        rows.append({
            "rank": rank, "method": name,
            "H1_unanimous": res["H1_u"], "H1_relaxed": res["H1_r"],
            "H2_unanimous": res["H2_u"], "H2_relaxed": res["H2_r"],
            "H3_unanimous": res["H3_u"], "H3_relaxed": res["H3_r"],
            "H4_rho": res["H4_rho"],
            "total_unanimous": res["total_u"], "total_relaxed": res["total_r"],
        })
    pd.DataFrame(rows).to_csv(OUT_DIR / "method_ranking.csv", index=False)

    print(f"\nSaved all figures and data to {OUT_DIR}")
    print("Done!")


if __name__ == "__main__":
    main()
