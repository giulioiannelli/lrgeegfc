#!/usr/bin/env python3
"""Sweep D(tau) distance/similarity metrics across patients, phases, bands.

For each of 8 metrics, compute pairwise (dis)similarity of the condensed
ultrametric vector between every pair of phases, for every (patient, band).
Then test four neuroscience hypotheses:

  H1  sim(task_learn, task_test) is highest among all phase pairs  (task persistence)
  H2  sim(task_test, rest_post) > sim(rest_pre, rest_post)                (task trace in rest)
  H3  within-type similarity > cross-type similarity             (reorganization gap)
  H4  Reorganization strength (std of pairwise values) decreases
      with frequency: delta > theta > ... > high_gamma

Outputs (in data/figures/metric_exploration/dtau_distances/):
  summary_table.pdf         -- heatmap: metric x hypothesis, bands passing unanimously
  best_metrics_detail.pdf   -- detailed panels for top 2-3 metrics
  per_metric_profiles.pdf   -- pair profiles across bands for each metric
  results.csv               -- all raw pairwise values
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import kendalltau, ks_2samp, pearsonr, spearmanr, wasserstein_distance

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

# -----------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------
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
ALL_PATIENTS = list(PATIENT_PHASES.keys())
BANDS = BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import LRG_CACHE
OUT_DIR = FIGURES_ROOT / "metric_exploration" / "dtau_distances"

REST = {"rest_pre", "rest_post"}
TASK = {"task_learn", "task_test"}

PHASE_ORDER = ["rest_pre", "task_learn", "task_test", "rest_post"]

# Ordered pair labels for 4-phase patients (6 pairs)
PAIR_LABELS_4PH = [
    ("rest_pre", "task_learn"),
    ("rest_pre", "task_test"),
    ("rest_pre", "rest_post"),
    ("task_learn", "task_test"),
    ("task_learn", "rest_post"),
    ("task_test", "rest_post"),
]
PAIR_SHORT = [f"{a[:2]}{a[-4:]}-{b[:2]}{b[-4:]}" for a, b in PAIR_LABELS_4PH]


def classify_pair(p1, p2):
    s = {p1, p2}
    if s <= REST or s <= TASK:
        return "within"
    return "cross"


# -----------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------
def load_all():
    """Load ultrametric data for all patients/phases/bands.

    Returns dict keyed by (patient, phase, band) with values containing
    the raw condensed ultrametric and its log transform.
    """
    data = {}
    n_loaded = 0
    n_failed = 0
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for ph in phases:
                lrg = load_lrg_result(pat, ph, band, "msc", cache_root=LRG_CACHE)
                if lrg is None:
                    n_failed += 1
                    continue
                U = lrg.ultrametric_matrix
                # Ensure condensed
                if U.ndim == 2:
                    from scipy.spatial.distance import squareform
                    U = squareform(U)
                U_log = np.log(np.clip(U, 1e-15, None))
                data[(pat, ph, band)] = {
                    "U": U,
                    "U_log": U_log,
                    "n_nodes": lrg.n_nodes,
                }
                n_loaded += 1
    print(f"  Loaded {n_loaded} entries ({n_failed} missing)")
    return data


# -----------------------------------------------------------------------
# Metrics
# -----------------------------------------------------------------------
# Convention: similarity metrics return higher = more similar.
# Distance metrics are negated at analysis time (see _to_similarity).

def m_log_spearman(d1, d2):
    """Spearman rank correlation of log(D(tau)).
    Rank-invariant so identical to Spearman(D), but using log for clarity."""
    return spearmanr(d1["U_log"], d2["U_log"])[0]


def m_log_pearson(d1, d2):
    """Pearson correlation of log(D(tau))."""
    return pearsonr(d1["U_log"], d2["U_log"])[0]


def m_log_frobenius(d1, d2):
    """Normalized L2 distance of log(D): ||log U1 - log U2||_2 / ||log U1||_2.
    DISTANCE metric (lower = more similar)."""
    a, b = d1["U_log"], d2["U_log"]
    diff = np.linalg.norm(a - b)
    norm = np.linalg.norm(a)
    return diff / max(norm, 1e-15)


def m_log_l1(d1, d2):
    """Normalized L1 distance of log(D): ||log U1 - log U2||_1 / ||log U1||_1.
    DISTANCE metric (lower = more similar)."""
    a, b = d1["U_log"], d2["U_log"]
    diff = np.sum(np.abs(a - b))
    norm = np.sum(np.abs(a))
    return diff / max(norm, 1e-15)


def m_ks_statistic(d1, d2):
    """Kolmogorov-Smirnov statistic between D(tau) distributions.
    DISTANCE metric (lower = more similar)."""
    return ks_2samp(d1["U"], d2["U"]).statistic


def m_emd(d1, d2):
    """Earth Mover's Distance (Wasserstein-1) between D(tau) distributions.
    DISTANCE metric (lower = more similar)."""
    return wasserstein_distance(d1["U"], d2["U"])


def m_rank_corr(d1, d2):
    """Kendall tau of cophenetic distances (no log)."""
    return kendalltau(d1["U"], d2["U"])[0]


def m_raw_spearman(d1, d2):
    """Spearman correlation of raw D(tau) (no log)."""
    return spearmanr(d1["U"], d2["U"])[0]


# Registry: (function, is_distance)
# is_distance=True means lower = more similar, will be negated for analysis
METRICS = {
    "logSpearman":   (m_log_spearman,   False),
    "logPearson":    (m_log_pearson,     False),
    "logFrobenius":  (m_log_frobenius,   True),
    "logL1":         (m_log_l1,          True),
    "KS_statistic":  (m_ks_statistic,    True),
    "EMD":           (m_emd,             True),
    "rankCorr":      (m_rank_corr,       False),
    "rawSpearman":   (m_raw_spearman,    False),
}
METRIC_NAMES = list(METRICS.keys())


# -----------------------------------------------------------------------
# Computation
# -----------------------------------------------------------------------
def compute_all_pairwise(data):
    """Compute all pairwise metric values.

    Returns a DataFrame with columns:
      metric, patient, band, phase1, phase2, raw_value, sim_value
    where sim_value has distance metrics negated so higher = more similar.
    """
    rows = []
    for metric_name, (metric_fn, is_distance) in METRICS.items():
        for pat, phases in PATIENT_PHASES.items():
            for band in BANDS:
                for p1, p2 in combinations(phases, 2):
                    key1 = (pat, p1, band)
                    key2 = (pat, p2, band)
                    if key1 not in data or key2 not in data:
                        continue
                    try:
                        val = metric_fn(data[key1], data[key2])
                        if val is None or np.isnan(val):
                            continue
                    except Exception as e:
                        continue
                    sim = -val if is_distance else val
                    rows.append({
                        "metric": metric_name,
                        "patient": pat,
                        "band": band,
                        "phase1": p1,
                        "phase2": p2,
                        "raw_value": val,
                        "sim_value": sim,
                        "pair_type": classify_pair(p1, p2),
                        "pair_label": f"{p1}-{p2}",
                    })
    df = pd.DataFrame(rows)
    print(f"  Computed {len(df)} pairwise values across {len(METRICS)} metrics")
    return df


# -----------------------------------------------------------------------
# Hypothesis testing
# -----------------------------------------------------------------------

def test_h1(df):
    """H1: sim(task_learn, task_test) is highest among all pairs.

    For each (metric, patient, band), check if the task-task pair has max similarity.
    Returns dict[metric][band] = (n_agree, n_total, patients_agree).
    """
    results = {}
    for metric in METRIC_NAMES:
        results[metric] = {}
        for band in BANDS:
            agree = []
            for pat in PATIENTS_4PH:
                sub = df[(df["metric"] == metric) & (df["patient"] == pat) & (df["band"] == band)]
                if len(sub) < 6:
                    continue
                tt_row = sub[(sub["phase1"] == "task_learn") & (sub["phase2"] == "task_test")]
                if tt_row.empty:
                    continue
                tt_sim = tt_row["sim_value"].values[0]
                max_sim = sub["sim_value"].max()
                agree.append(tt_sim >= max_sim - 1e-12)
            results[metric][band] = (sum(agree), len(agree), agree)
    return results


def test_h2(df):
    """H2: sim(task_test, rest_post) > sim(rest_pre, rest_post).

    Task leaves a trace in post-task rest.
    Returns dict[metric][band] = (n_agree, n_total, agree_list).
    """
    results = {}
    for metric in METRIC_NAMES:
        results[metric] = {}
        for band in BANDS:
            agree = []
            for pat in PATIENTS_4PH:
                sub = df[(df["metric"] == metric) & (df["patient"] == pat) & (df["band"] == band)]
                tt_post = sub[(sub["phase1"] == "task_test") & (sub["phase2"] == "rest_post")]
                pre_post = sub[(sub["phase1"] == "rest_pre") & (sub["phase2"] == "rest_post")]
                if tt_post.empty or pre_post.empty:
                    continue
                agree.append(tt_post["sim_value"].values[0] > pre_post["sim_value"].values[0])
            results[metric][band] = (sum(agree), len(agree), agree)
    return results


def test_h3(df):
    """H3: within-type similarity > cross-type similarity.

    For 4-phase patients, compare mean(within) vs mean(cross).
    Returns dict[metric][band] = (n_agree, n_total, agree_list).
    """
    results = {}
    for metric in METRIC_NAMES:
        results[metric] = {}
        for band in BANDS:
            agree = []
            for pat in PATIENTS_4PH:
                sub = df[(df["metric"] == metric) & (df["patient"] == pat) & (df["band"] == band)]
                within_mean = sub[sub["pair_type"] == "within"]["sim_value"].mean()
                cross_mean = sub[sub["pair_type"] == "cross"]["sim_value"].mean()
                if np.isnan(within_mean) or np.isnan(cross_mean):
                    continue
                agree.append(within_mean > cross_mean)
            results[metric][band] = (sum(agree), len(agree), agree)
    return results


def test_h4(df):
    """H4: Reorganization strength (std of pairwise sim) decreases with frequency.

    For each (metric, patient), compute std of sim_value per band, then check
    if the sequence delta>theta>...>high_gamma is monotonically decreasing.
    Returns dict[metric] = (n_patients_monotone, n_patients_total, rank_corr_with_band_idx).
    """
    results = {}
    for metric in METRIC_NAMES:
        band_stds_per_patient = []
        for pat in PATIENTS_4PH:
            stds = []
            for band in BANDS:
                sub = df[(df["metric"] == metric) & (df["patient"] == pat) & (df["band"] == band)]
                if len(sub) >= 3:
                    stds.append(sub["sim_value"].std())
                else:
                    stds.append(np.nan)
            band_stds_per_patient.append(stds)

        # Check monotone decrease
        n_monotone = 0
        all_corrs = []
        for stds in band_stds_per_patient:
            if any(np.isnan(stds)):
                continue
            # Spearman correlation of std with band index (should be negative)
            rho, _ = spearmanr(range(len(BANDS)), stds)
            all_corrs.append(rho)
            # Check strict monotone decrease
            diffs = np.diff(stds)
            if all(d < 0 for d in diffs):
                n_monotone += 1

        mean_corr = np.mean(all_corrs) if all_corrs else np.nan
        results[metric] = (n_monotone, len(band_stds_per_patient), mean_corr)
    return results


# -----------------------------------------------------------------------
# Scoring
# -----------------------------------------------------------------------
def compute_scores(df):
    """Run all hypothesis tests and produce summary."""
    h1 = test_h1(df)
    h2 = test_h2(df)
    h3 = test_h3(df)
    h4 = test_h4(df)

    # Build summary: for each (metric, hypothesis), count bands with unanimous agreement
    summary = {}
    for metric in METRIC_NAMES:
        row = {}

        # H1: unanimous = all 4 patients agree
        h1_unanimous_bands = sum(
            1 for band in BANDS
            if h1[metric][band][1] >= 4 and h1[metric][band][0] == h1[metric][band][1]
        )
        h1_relaxed_bands = sum(
            1 for band in BANDS
            if h1[metric][band][1] >= 4 and h1[metric][band][0] >= 3
        )
        row["H1_unan"] = h1_unanimous_bands
        row["H1_relax"] = h1_relaxed_bands

        # H2
        h2_unanimous_bands = sum(
            1 for band in BANDS
            if h2[metric][band][1] >= 4 and h2[metric][band][0] == h2[metric][band][1]
        )
        h2_relaxed_bands = sum(
            1 for band in BANDS
            if h2[metric][band][1] >= 4 and h2[metric][band][0] >= 3
        )
        row["H2_unan"] = h2_unanimous_bands
        row["H2_relax"] = h2_relaxed_bands

        # H3
        h3_unanimous_bands = sum(
            1 for band in BANDS
            if h3[metric][band][1] >= 4 and h3[metric][band][0] == h3[metric][band][1]
        )
        h3_relaxed_bands = sum(
            1 for band in BANDS
            if h3[metric][band][1] >= 4 and h3[metric][band][0] >= 3
        )
        row["H3_unan"] = h3_unanimous_bands
        row["H3_relax"] = h3_relaxed_bands

        # H4
        row["H4_monotone"] = h4[metric][0]
        row["H4_mean_rho"] = h4[metric][2]

        # Composite score: sum of unanimous bands across H1-H3 + bonus for H4
        row["composite"] = (
            row["H1_unan"] + row["H2_unan"] + row["H3_unan"]
            + (1 if h4[metric][0] >= 2 else 0)
        )

        summary[metric] = row

    return summary, h1, h2, h3, h4


def compute_zscore_consistency(df):
    """For each metric, compute z-scored similarity within each (patient, band).

    Returns DataFrame with z-scores and sign consistency info.
    """
    rows = []
    for metric in METRIC_NAMES:
        for pat in PATIENTS_4PH:
            for band in BANDS:
                sub = df[(df["metric"] == metric) & (df["patient"] == pat) & (df["band"] == band)].copy()
                if len(sub) < 3:
                    continue
                mu = sub["sim_value"].mean()
                sd = sub["sim_value"].std()
                if sd < 1e-15:
                    continue
                for _, r in sub.iterrows():
                    z = (r["sim_value"] - mu) / sd
                    rows.append({
                        "metric": metric,
                        "patient": pat,
                        "band": band,
                        "phase1": r["phase1"],
                        "phase2": r["phase2"],
                        "sim_value": r["sim_value"],
                        "z_score": z,
                        "pair_type": r["pair_type"],
                    })
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------
def plot_summary_table(summary, out_dir):
    """Heatmap: metric x hypothesis, showing how many bands pass unanimously."""
    metrics = METRIC_NAMES
    hyp_cols = ["H1_unan", "H1_relax", "H2_unan", "H2_relax", "H3_unan", "H3_relax"]
    hyp_labels = [
        "H1 unan\n(task-task\nhighest)",
        "H1 relax\n(3/4)",
        "H2 unan\n(task_test-rest_post\n> rest_pre-rest_post)",
        "H2 relax\n(3/4)",
        "H3 unan\n(within >\ncross)",
        "H3 relax\n(3/4)",
    ]

    mat = np.zeros((len(metrics), len(hyp_cols)))
    for i, m in enumerate(metrics):
        for j, h in enumerate(hyp_cols):
            mat[i, j] = summary[m][h]

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=6)
    ax.set_xticks(range(len(hyp_cols)))
    ax.set_xticklabels(hyp_labels, fontsize=8, ha="center")
    ax.set_yticks(range(len(metrics)))
    ax.set_yticklabels(metrics, fontsize=9)
    ax.set_xlabel("Hypothesis test")
    ax.set_ylabel("Metric")
    ax.set_title("D(tau) metric comparison: bands passing each hypothesis (out of 6)")

    # Annotate cells
    for i in range(len(metrics)):
        for j in range(len(hyp_cols)):
            v = int(mat[i, j])
            color = "white" if v >= 4 else "black"
            ax.text(j, i, str(v), ha="center", va="center", fontsize=10,
                    fontweight="bold" if v >= 4 else "normal", color=color)

    # Add H4 info as right-side annotations
    for i, m in enumerate(metrics):
        h4_mono = summary[m]["H4_monotone"]
        h4_rho = summary[m]["H4_mean_rho"]
        ax.text(len(hyp_cols) + 0.3, i,
                f"H4: {h4_mono}/4 mono, rho={h4_rho:+.2f}",
                va="center", fontsize=7, color="darkblue")

    # Add composite score
    for i, m in enumerate(metrics):
        ax.text(-0.7, i, f"[{summary[m]['composite']}]",
                va="center", ha="right", fontsize=8, fontweight="bold", color="darkred")

    plt.colorbar(im, ax=ax, label="Bands passing (out of 6)", shrink=0.8)
    fig.tight_layout()
    fig.savefig(out_dir / "summary_table.pdf", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved summary_table.pdf")


def plot_per_metric_profiles(df, out_dir):
    """For each metric, show per-pair similarity profiles across bands for 4-phase patients."""
    n_metrics = len(METRIC_NAMES)
    fig, axes = plt.subplots(2, 4, figsize=(20, 10), sharey=False)
    axes = axes.flat

    band_tex = [BRAIN_BAND_TEX_DICT[b] for b in BANDS]

    # Color map for pairs
    pair_colors = {
        "rest_pre-task_learn": "#1f77b4",
        "rest_pre-task_test": "#ff7f0e",
        "rest_pre-rest_post": "#2ca02c",
        "task_learn-task_test": "#d62728",
        "task_learn-rest_post": "#9467bd",
        "task_test-rest_post": "#8c564b",
    }

    for idx, metric in enumerate(METRIC_NAMES):
        ax = axes[idx]

        for p1, p2 in PAIR_LABELS_4PH:
            pair_key = f"{p1}-{p2}"
            vals = []
            for band in BANDS:
                patient_vals = []
                for pat in PATIENTS_4PH:
                    sub = df[(df["metric"] == metric) & (df["patient"] == pat)
                             & (df["band"] == band)
                             & (df["phase1"] == p1) & (df["phase2"] == p2)]
                    if not sub.empty:
                        patient_vals.append(sub["sim_value"].values[0])
                vals.append(np.mean(patient_vals) if patient_vals else np.nan)

            ax.plot(range(len(BANDS)), vals, "o-", label=pair_key,
                    color=pair_colors.get(pair_key, "gray"), markersize=4, linewidth=1.5)

        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels(band_tex, fontsize=9)
        ax.set_title(metric, fontsize=11, fontweight="bold")
        ax.set_ylabel("Similarity" if idx % 4 == 0 else "")
        ax.grid(True, alpha=0.3)

    # Legend on last panel
    axes[n_metrics - 1].legend(fontsize=7, loc="best")

    fig.suptitle("Per-metric pair profiles across bands (mean over 4-phase patients)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_dir / "per_metric_profiles.pdf", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved per_metric_profiles.pdf")


def plot_best_metrics_detail(df, summary, h1, h2, h3, out_dir):
    """Detailed analysis panels for the top 2-3 metrics by composite score."""
    # Sort by composite score
    ranking = sorted(METRIC_NAMES, key=lambda m: -summary[m]["composite"])
    top_metrics = ranking[:3]

    band_tex = [BRAIN_BAND_TEX_DICT[b] for b in BANDS]

    fig, axes = plt.subplots(3, 4, figsize=(22, 15))

    for row_idx, metric in enumerate(top_metrics):
        is_distance = METRICS[metric][1]

        # --- Column 0: Per-patient pairwise heatmap (averaged over bands) ---
        ax = axes[row_idx, 0]
        # Build 4x4 phase matrix averaged over bands for each 4-phase patient
        mat = np.full((4, 4), np.nan)
        for i, p1 in enumerate(PHASE_ORDER):
            for j, p2 in enumerate(PHASE_ORDER):
                if i == j:
                    continue
                pi, pj = (p1, p2) if PHASE_ORDER.index(p1) < PHASE_ORDER.index(p2) else (p2, p1)
                vals = []
                for pat in PATIENTS_4PH:
                    for band in BANDS:
                        sub = df[(df["metric"] == metric) & (df["patient"] == pat)
                                 & (df["band"] == band)
                                 & (df["phase1"] == pi) & (df["phase2"] == pj)]
                        if not sub.empty:
                            vals.append(sub["sim_value"].values[0])
                if vals:
                    mat[i, j] = np.mean(vals)

        vmin, vmax = np.nanmin(mat), np.nanmax(mat)
        if vmin < 0 < vmax:
            norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
            cmap = "RdBu_r"
        else:
            norm = None
            cmap = "viridis"
        im = ax.imshow(mat, cmap=cmap, norm=norm)
        ax.set_xticks(range(4))
        ax.set_xticklabels(["rest_pre", "tLearn", "tTest", "rest_post"], fontsize=7, rotation=45)
        ax.set_yticks(range(4))
        ax.set_yticklabels(["rest_pre", "tLearn", "tTest", "rest_post"], fontsize=7)
        for i in range(4):
            for j in range(4):
                if not np.isnan(mat[i, j]):
                    ax.text(j, i, f"{mat[i,j]:.3f}", ha="center", va="center", fontsize=6)
        ax.set_title(f"{metric}\nMean sim matrix", fontsize=9, fontweight="bold")
        plt.colorbar(im, ax=ax, shrink=0.7)

        # --- Column 1: H1 + H2 detail per band ---
        ax = axes[row_idx, 1]
        x = np.arange(len(BANDS))
        h1_fracs = []
        h2_fracs = []
        for band in BANDS:
            n_ok, n_tot, _ = h1[metric][band]
            h1_fracs.append(n_ok / max(n_tot, 1))
            n_ok2, n_tot2, _ = h2[metric][band]
            h2_fracs.append(n_ok2 / max(n_tot2, 1))
        w = 0.35
        ax.bar(x - w/2, h1_fracs, w, label="H1 (task-task highest)", color="#d62728", alpha=0.8)
        ax.bar(x + w/2, h2_fracs, w, label="H2 (task trace)", color="#2ca02c", alpha=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(band_tex, fontsize=9)
        ax.set_ylabel("Fraction of patients agreeing")
        ax.set_ylim(0, 1.1)
        ax.axhline(0.75, ls="--", color="gray", alpha=0.5, label="3/4 threshold")
        ax.axhline(1.0, ls=":", color="gray", alpha=0.3)
        ax.legend(fontsize=7)
        ax.set_title(f"H1 & H2 per band", fontsize=9)

        # --- Column 2: H3 gap (within - cross) per patient per band ---
        ax = axes[row_idx, 2]
        for pat in PATIENTS_4PH:
            gaps = []
            for band in BANDS:
                sub = df[(df["metric"] == metric) & (df["patient"] == pat) & (df["band"] == band)]
                w_mean = sub[sub["pair_type"] == "within"]["sim_value"].mean()
                c_mean = sub[sub["pair_type"] == "cross"]["sim_value"].mean()
                gaps.append(w_mean - c_mean if not (np.isnan(w_mean) or np.isnan(c_mean)) else np.nan)
            ax.plot(range(len(BANDS)), gaps, "o-", label=pat.replace("Pat_0", "P"),
                    markersize=4, linewidth=1.5)
        ax.axhline(0, ls="-", color="black", alpha=0.3)
        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels(band_tex, fontsize=9)
        ax.set_ylabel("Gap (within - cross)")
        ax.legend(fontsize=7)
        ax.set_title("H3: reorganization gap", fontsize=9)
        ax.grid(True, alpha=0.3)

        # --- Column 3: H4 reorganization strength (std) across bands ---
        ax = axes[row_idx, 3]
        for pat in PATIENTS_4PH:
            stds = []
            for band in BANDS:
                sub = df[(df["metric"] == metric) & (df["patient"] == pat) & (df["band"] == band)]
                stds.append(sub["sim_value"].std() if len(sub) >= 3 else np.nan)
            ax.plot(range(len(BANDS)), stds, "s-", label=pat.replace("Pat_0", "P"),
                    markersize=4, linewidth=1.5)
        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels(band_tex, fontsize=9)
        ax.set_ylabel("Std of pairwise sim")
        ax.legend(fontsize=7)
        ax.set_title(f"H4: reorg strength vs freq (rho={summary[metric]['H4_mean_rho']:+.2f})", fontsize=9)
        ax.grid(True, alpha=0.3)

    fig.suptitle("Top 3 metrics -- detailed hypothesis analysis",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_dir / "best_metrics_detail.pdf", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved best_metrics_detail.pdf")


# -----------------------------------------------------------------------
# Print summary
# -----------------------------------------------------------------------
def print_summary(summary, h1, h2, h3, h4):
    """Print a clear console summary."""
    print("\n" + "=" * 90)
    print("D(tau) DISTANCE METRIC SWEEP -- RESULTS SUMMARY")
    print("=" * 90)

    # Overall ranking
    ranking = sorted(METRIC_NAMES, key=lambda m: -summary[m]["composite"])

    print(f"\n{'Rank':<5} {'Metric':<15} {'H1u':>4} {'H1r':>4} {'H2u':>4} {'H2r':>4} "
          f"{'H3u':>4} {'H3r':>4} {'H4m':>4} {'H4rho':>7} {'Score':>6}")
    print("-" * 80)
    for rank, m in enumerate(ranking, 1):
        s = summary[m]
        print(f"{rank:<5} {m:<15} {s['H1_unan']:>4} {s['H1_relax']:>4} "
              f"{s['H2_unan']:>4} {s['H2_relax']:>4} "
              f"{s['H3_unan']:>4} {s['H3_relax']:>4} "
              f"{s['H4_monotone']:>4} {s['H4_mean_rho']:>+7.3f} {s['composite']:>6}")

    # Per-band detail for top 3
    print("\n--- Per-band detail for top 3 metrics ---")
    for m in ranking[:3]:
        print(f"\n  {m} (composite={summary[m]['composite']}):")
        print(f"    {'Band':<12} {'H1':>10} {'H2':>10} {'H3':>10}")
        print(f"    {'-'*44}")
        for band in BANDS:
            h1_n, h1_t, _ = h1[m][band]
            h2_n, h2_t, _ = h2[m][band]
            h3_n, h3_t, _ = h3[m][band]
            h1_str = f"{h1_n}/{h1_t}" + (" ***" if h1_n == h1_t >= 4 else " *" if h1_n >= 3 else "")
            h2_str = f"{h2_n}/{h2_t}" + (" ***" if h2_n == h2_t >= 4 else " *" if h2_n >= 3 else "")
            h3_str = f"{h3_n}/{h3_t}" + (" ***" if h3_n == h3_t >= 4 else " *" if h3_n >= 3 else "")
            print(f"    {band:<12} {h1_str:>10} {h2_str:>10} {h3_str:>10}")

    # H4 detail
    print("\n--- H4 detail (reorganization strength vs frequency) ---")
    print(f"  {'Metric':<15} {'Monotone':>10} {'Mean rho':>10}")
    for m in ranking:
        mono, total, rho = h4[m]
        print(f"  {m:<15} {mono}/{total:>9} {rho:>+10.3f}")

    # Best metric recommendation
    best = ranking[0]
    print(f"\n{'='*60}")
    print(f"BEST METRIC: {best}")
    print(f"  H1 (task persistence): {summary[best]['H1_unan']}/6 bands unanimous, "
          f"{summary[best]['H1_relax']}/6 relaxed")
    print(f"  H2 (task trace):       {summary[best]['H2_unan']}/6 bands unanimous, "
          f"{summary[best]['H2_relax']}/6 relaxed")
    print(f"  H3 (reorg gap):        {summary[best]['H3_unan']}/6 bands unanimous, "
          f"{summary[best]['H3_relax']}/6 relaxed")
    print(f"  H4 (freq gradient):    {h4[best][0]}/4 patients monotone, "
          f"mean rho={h4[best][2]:+.3f}")
    print(f"{'='*60}")


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 90)
    print("Sweep D(tau) distance metrics")
    print("=" * 90)

    print("\n[1/5] Loading data...")
    data = load_all()

    print("\n[2/5] Computing all pairwise metrics...")
    df = compute_all_pairwise(data)

    print("\n[3/5] Testing hypotheses...")
    summary, h1, h2, h3, h4 = compute_scores(df)

    print("\n[4/5] Generating figures...")
    plot_summary_table(summary, OUT_DIR)
    plot_per_metric_profiles(df, OUT_DIR)
    plot_best_metrics_detail(df, summary, h1, h2, h3, OUT_DIR)

    print("\n[5/5] Saving results...")
    df.to_csv(OUT_DIR / "results.csv", index=False)
    print(f"  Saved results.csv ({len(df)} rows)")

    # Summary table as CSV too
    summary_df = pd.DataFrame(summary).T
    summary_df.index.name = "metric"
    summary_df.to_csv(OUT_DIR / "summary_scores.csv")
    print(f"  Saved summary_scores.csv")

    print_summary(summary, h1, h2, h3, h4)


if __name__ == "__main__":
    main()
