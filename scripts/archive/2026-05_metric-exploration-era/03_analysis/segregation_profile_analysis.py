#!/usr/bin/env python3
"""Scale-Dependent Segregation-Integration Analysis.

For each hierarchical level k (number of clusters from the LRG dendrogram),
compute within- and between-community FC and a segregation index:

    S(k) = (within_FC - between_FC) / (within_FC + between_FC)

This links the LRG hierarchical partition directly to the MSC functional
connectivity, giving a segregation-vs-scale profile.

Output: data/figures/multiscale_investigation/segregation_profiles/
"""
from __future__ import annotations

import warnings
from pathlib import Path
from collections import defaultdict

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    BRAIN_BANDS_NAMES,
    PHASE_LABELS,
    BRAIN_BAND_TEX_DICT,
)
from lrg_eegfc.config.paths import FIGURES_ROOT

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
ALL_PHASES = list(PHASE_LABELS)
BANDS = list(BRAIN_BANDS_NAMES)
FC_METHOD = "msc"
K_RANGE = range(2, 31)  # k=2..30

OUTPUT_DIR = FIGURES_ROOT / "multiscale_investigation" / "segregation_profiles"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Phase display config
PHASE_COLORS = {
    "rest_pre": "#1f77b4",
    "task_learn": "#ff7f0e",
    "task_test": "#2ca02c",
    "rest_post": "#d62728",
}
PHASE_DISPLAY = {
    "rest_pre": "Rest Pre",
    "task_learn": "Task Learn",
    "task_test": "Task Test",
    "rest_post": "Rest Post",
}

CONDITION_MAP = {
    "rest_pre": "rest",
    "rest_post": "rest",
    "task_learn": "task",
    "task_test": "task",
}


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------
def compute_segregation_profile(
    linkage_matrix: np.ndarray,
    fc_matrix: np.ndarray,
    k_range: range,
) -> dict:
    """Compute segregation index S(k) for each k in k_range.

    Returns dict with keys: k_values, within_fc, between_fc, segregation.
    """
    n = fc_matrix.shape[0]
    # Pre-compute upper-triangle mask
    triu_i, triu_j = np.triu_indices(n, k=1)
    fc_upper = fc_matrix[triu_i, triu_j]

    k_values = []
    within_fc_list = []
    between_fc_list = []
    segregation_list = []

    for k in k_range:
        if k > n:
            break
        labels = fcluster(linkage_matrix, k, criterion="maxclust")

        # Classify edges
        same_cluster = labels[triu_i] == labels[triu_j]

        within_vals = fc_upper[same_cluster]
        between_vals = fc_upper[~same_cluster]

        if len(within_vals) == 0 or len(between_vals) == 0:
            within_mean = np.nan
            between_mean = np.nan
            seg = np.nan
        else:
            within_mean = np.mean(within_vals)
            between_mean = np.mean(between_vals)
            denom = within_mean + between_mean
            seg = (within_mean - between_mean) / denom if denom != 0 else np.nan

        k_values.append(k)
        within_fc_list.append(within_mean)
        between_fc_list.append(between_mean)
        segregation_list.append(seg)

    return {
        "k_values": np.array(k_values),
        "within_fc": np.array(within_fc_list),
        "between_fc": np.array(between_fc_list),
        "segregation": np.array(segregation_list),
    }


def extract_giant_component_fc(fc_matrix: np.ndarray, n_target: int) -> np.ndarray:
    """Extract giant-component submatrix from full MSC matrix to match LRG size."""
    if fc_matrix.shape[0] == n_target:
        return fc_matrix
    G = nx.from_numpy_array(fc_matrix)
    gcc_nodes = max(nx.connected_components(G), key=len)
    gcc_indices = sorted(gcc_nodes)
    fc_gcc = fc_matrix[np.ix_(gcc_indices, gcc_indices)]
    if fc_gcc.shape[0] != n_target:
        raise ValueError(
            f"Giant component has {fc_gcc.shape[0]} nodes but LRG expects {n_target}"
        )
    return fc_gcc


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_all_profiles() -> dict:
    """Load LRG + MSC for all patient/band/phase combos and compute profiles.

    Returns nested dict: profiles[patient][band][phase] = profile_dict
    Also returns a flat list of records for CSV export.
    """
    profiles = defaultdict(lambda: defaultdict(dict))
    records = []
    n_loaded = 0
    n_skipped = 0

    for patient in PATIENTS:
        for band in BANDS:
            for phase in ALL_PHASES:
                # Load LRG result
                lrg = load_lrg_result(patient, phase, band, FC_METHOD)
                if lrg is None:
                    n_skipped += 1
                    continue

                # Load MSC matrix (dense, nperseg=4096)
                msc = load_msc_matrix(patient, phase, band)
                if msc is None:
                    n_skipped += 1
                    continue

                # Extract giant component to match LRG
                try:
                    fc_gcc = extract_giant_component_fc(msc, lrg.n_nodes)
                except ValueError as e:
                    print(f"  SKIP {patient}/{band}/{phase}: {e}")
                    n_skipped += 1
                    continue

                # Compute segregation profile
                profile = compute_segregation_profile(
                    lrg.linkage_matrix, fc_gcc, K_RANGE
                )
                profiles[patient][band][phase] = profile
                n_loaded += 1

                # Store records
                for i, k in enumerate(profile["k_values"]):
                    records.append({
                        "patient": patient,
                        "band": band,
                        "phase": phase,
                        "k": int(k),
                        "within_fc": profile["within_fc"][i],
                        "between_fc": profile["between_fc"][i],
                        "segregation": profile["segregation"][i],
                    })

    print(f"Loaded {n_loaded} profiles, skipped {n_skipped}")
    return profiles, records


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------
def _style_ax(ax, xlabel="Number of clusters (k)", ylabel="Segregation index S(k)"):
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _mark_peak(ax, k_values, seg_values, color, marker="v", ms=8):
    """Mark the k of maximum segregation."""
    valid = ~np.isnan(seg_values)
    if not np.any(valid):
        return None
    idx = np.nanargmax(seg_values)
    ax.plot(k_values[idx], seg_values[idx], marker=marker, color=color,
            markersize=ms, zorder=5, markeredgecolor="black", markeredgewidth=0.5)
    return k_values[idx]


# ---------------------------------------------------------------------------
# Figure 1: Pat_02 S(k) profiles, one subplot per band
# ---------------------------------------------------------------------------
def plot_fig1(profiles):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey=True)
    fig.suptitle("Segregation Profile S(k) -- Pat_02", fontsize=14, fontweight="bold")

    for idx, band in enumerate(BANDS):
        ax = axes.flat[idx]
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)

        for phase in ALL_PHASES:
            if phase not in profiles.get("Pat_02", {}).get(band, {}):
                continue
            p = profiles["Pat_02"][band][phase]
            ax.plot(p["k_values"], p["segregation"],
                    color=PHASE_COLORS[phase], label=PHASE_DISPLAY[phase],
                    linewidth=1.5, alpha=0.85)
            _mark_peak(ax, p["k_values"], p["segregation"], PHASE_COLORS[phase])

        ax.axhline(0, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)
        _style_ax(ax)

    # Single legend
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])

    path = OUTPUT_DIR / "fig1_pat02_segregation_profiles.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")


# ---------------------------------------------------------------------------
# Figure 2: Scale of max segregation (bar chart)
# ---------------------------------------------------------------------------
def plot_fig2(profiles):
    """Bar chart: k* (scale of max segregation) per band, averaged over patients."""
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=False, sharey=True)
    fig.suptitle("Scale of Maximum Segregation (k*)", fontsize=14, fontweight="bold")

    for idx, band in enumerate(BANDS):
        ax = axes.flat[idx]
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)

        phase_k_stars = {}
        for phase in ALL_PHASES:
            k_stars = []
            for patient in PATIENTS:
                if phase not in profiles.get(patient, {}).get(band, {}):
                    continue
                p = profiles[patient][band][phase]
                valid = ~np.isnan(p["segregation"])
                if np.any(valid):
                    k_stars.append(p["k_values"][np.nanargmax(p["segregation"])])
            phase_k_stars[phase] = k_stars

        phases_present = [ph for ph in ALL_PHASES if phase_k_stars.get(ph)]
        x = np.arange(len(phases_present))
        means = [np.mean(phase_k_stars[ph]) for ph in phases_present]
        sems = [np.std(phase_k_stars[ph]) / np.sqrt(len(phase_k_stars[ph]))
                for ph in phases_present]
        colors = [PHASE_COLORS[ph] for ph in phases_present]

        ax.bar(x, means, yerr=sems, color=colors, alpha=0.8, capsize=3,
               edgecolor="black", linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels([PHASE_DISPLAY[ph] for ph in phases_present],
                           fontsize=7, rotation=25, ha="right")
        _style_ax(ax, xlabel="", ylabel="k* (max segregation scale)")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    path = OUTPUT_DIR / "fig2_max_segregation_scale.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")


# ---------------------------------------------------------------------------
# Figure 3: Maximum segregation value heatmap (phases x bands)
# ---------------------------------------------------------------------------
def plot_fig3(profiles):
    """Heatmap: max S(k) value for each phase x band (averaged over patients)."""
    data = np.full((len(ALL_PHASES), len(BANDS)), np.nan)
    counts = np.zeros_like(data)

    for pi, phase in enumerate(ALL_PHASES):
        for bi, band in enumerate(BANDS):
            vals = []
            for patient in PATIENTS:
                if phase not in profiles.get(patient, {}).get(band, {}):
                    continue
                p = profiles[patient][band][phase]
                if np.any(~np.isnan(p["segregation"])):
                    vals.append(np.nanmax(p["segregation"]))
            if vals:
                data[pi, bi] = np.mean(vals)
                counts[pi, bi] = len(vals)

    fig, ax = plt.subplots(figsize=(10, 4))
    im = ax.imshow(data, aspect="auto", cmap="YlOrRd", interpolation="nearest")
    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=10)
    ax.set_yticks(range(len(ALL_PHASES)))
    ax.set_yticklabels([PHASE_DISPLAY[ph] for ph in ALL_PHASES], fontsize=10)
    ax.set_title("Maximum Segregation S(k*) -- Patient Average", fontsize=13,
                 fontweight="bold")

    # Annotate cells
    for pi in range(len(ALL_PHASES)):
        for bi in range(len(BANDS)):
            val = data[pi, bi]
            if not np.isnan(val):
                ax.text(bi, pi, f"{val:.3f}", ha="center", va="center",
                        fontsize=8, color="black" if val < 0.15 else "white")

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("max S(k)", fontsize=10)

    fig.tight_layout()
    path = OUTPUT_DIR / "fig3_max_segregation_heatmap.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")


# ---------------------------------------------------------------------------
# Figure 4: Grand average rest vs task S(k)
# ---------------------------------------------------------------------------
def plot_fig4(profiles):
    """Grand-average S(k): rest vs task, averaged across patients and bands."""
    rest_profiles = defaultdict(list)  # k -> list of S values
    task_profiles = defaultdict(list)

    for patient in PATIENTS:
        for band in BANDS:
            for phase in ALL_PHASES:
                if phase not in profiles.get(patient, {}).get(band, {}):
                    continue
                p = profiles[patient][band][phase]
                target = rest_profiles if CONDITION_MAP[phase] == "rest" else task_profiles
                for i, k in enumerate(p["k_values"]):
                    target[int(k)].append(p["segregation"][i])

    k_vals = sorted(rest_profiles.keys() & task_profiles.keys())
    rest_mean = np.array([np.nanmean(rest_profiles[k]) for k in k_vals])
    rest_sem = np.array([np.nanstd(rest_profiles[k]) / np.sqrt(
        np.sum(~np.isnan(rest_profiles[k]))) for k in k_vals])
    task_mean = np.array([np.nanmean(task_profiles[k]) for k in k_vals])
    task_sem = np.array([np.nanstd(task_profiles[k]) / np.sqrt(
        np.sum(~np.isnan(task_profiles[k]))) for k in k_vals])

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(k_vals, rest_mean, color="#1f77b4", linewidth=2, label="Rest")
    ax.fill_between(k_vals, rest_mean - rest_sem, rest_mean + rest_sem,
                    color="#1f77b4", alpha=0.2)
    ax.plot(k_vals, task_mean, color="#ff7f0e", linewidth=2, label="Task")
    ax.fill_between(k_vals, task_mean - task_sem, task_mean + task_sem,
                    color="#ff7f0e", alpha=0.2)

    _mark_peak(ax, np.array(k_vals), rest_mean, "#1f77b4", ms=10)
    _mark_peak(ax, np.array(k_vals), task_mean, "#ff7f0e", ms=10)

    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)
    _style_ax(ax)
    ax.set_title("Grand Average Segregation Profile: Rest vs Task",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=10, frameon=False)

    fig.tight_layout()
    path = OUTPUT_DIR / "fig4_rest_vs_task_grand_average.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")


# ---------------------------------------------------------------------------
# Figure 5: Per-patient profiles for alpha band
# ---------------------------------------------------------------------------
def plot_fig5(profiles):
    """Per-patient S(k) profiles for the alpha band."""
    patients_with_data = [p for p in PATIENTS
                          if "alpha" in profiles.get(p, {})]

    n = len(patients_with_data)
    cols = min(3, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows),
                             sharex=True, sharey=True)
    if n == 1:
        axes = np.array([axes])
    axes_flat = axes.flat if hasattr(axes, 'flat') else [axes]

    fig.suptitle(r"Segregation Profile -- $\alpha$ band (per patient)",
                 fontsize=14, fontweight="bold")

    for idx, patient in enumerate(patients_with_data):
        ax = axes_flat[idx]
        ax.set_title(patient, fontsize=11)
        for phase in ALL_PHASES:
            if phase not in profiles[patient].get("alpha", {}):
                continue
            p = profiles[patient]["alpha"][phase]
            ax.plot(p["k_values"], p["segregation"],
                    color=PHASE_COLORS[phase], label=PHASE_DISPLAY[phase],
                    linewidth=1.5, alpha=0.85)
            _mark_peak(ax, p["k_values"], p["segregation"], PHASE_COLORS[phase])
        ax.axhline(0, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)
        _style_ax(ax)

    # Hide unused axes
    for idx in range(n, len(list(axes_flat))):
        axes_flat[idx].set_visible(False)

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])

    path = OUTPUT_DIR / "fig5_alpha_per_patient.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")


# ---------------------------------------------------------------------------
# Figure 6: Delta-segregation (task - rest) profiles
# ---------------------------------------------------------------------------
def plot_fig6(profiles):
    """Delta-segregation: task minus rest S(k) profiles per band.

    Average (task_learn+task_test)/2 - (rest_pre+rest_post)/2 per patient,
    then average across patients.
    """
    rest_phases = ["rest_pre", "rest_post"]
    task_phases = ["task_learn", "task_test"]

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey=True)
    fig.suptitle(r"$\Delta$Segregation: Task $-$ Rest",
                 fontsize=14, fontweight="bold")

    for idx, band in enumerate(BANDS):
        ax = axes.flat[idx]
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)

        # Compute per-patient delta profiles
        delta_per_patient = {}
        for patient in PATIENTS:
            band_data = profiles.get(patient, {}).get(band, {})
            rest_available = [ph for ph in rest_phases if ph in band_data]
            task_available = [ph for ph in task_phases if ph in band_data]

            if not rest_available or not task_available:
                continue

            # Get common k range
            all_k = [band_data[ph]["k_values"] for ph in rest_available + task_available]
            common_k = all_k[0]
            for k_arr in all_k[1:]:
                common_k = np.intersect1d(common_k, k_arr)

            if len(common_k) == 0:
                continue

            # Average rest and task profiles over available phases
            rest_seg = np.zeros(len(common_k))
            for ph in rest_available:
                p = band_data[ph]
                mask = np.isin(p["k_values"], common_k)
                rest_seg += p["segregation"][mask]
            rest_seg /= len(rest_available)

            task_seg = np.zeros(len(common_k))
            for ph in task_available:
                p = band_data[ph]
                mask = np.isin(p["k_values"], common_k)
                task_seg += p["segregation"][mask]
            task_seg /= len(task_available)

            delta_per_patient[patient] = (common_k, task_seg - rest_seg)

        if not delta_per_patient:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
            _style_ax(ax, ylabel=r"$\Delta$S(k)")
            continue

        # Average across patients using common k
        all_common_k = list(delta_per_patient.values())[0][0]
        for _, (ck, _) in delta_per_patient.items():
            all_common_k = np.intersect1d(all_common_k, ck)

        delta_matrix = []
        for patient, (ck, delta) in delta_per_patient.items():
            mask = np.isin(ck, all_common_k)
            delta_matrix.append(delta[mask])

        delta_matrix = np.array(delta_matrix)
        mean_delta = np.nanmean(delta_matrix, axis=0)
        sem_delta = np.nanstd(delta_matrix, axis=0) / np.sqrt(delta_matrix.shape[0])

        ax.plot(all_common_k, mean_delta, color="#9467bd", linewidth=2)
        ax.fill_between(all_common_k, mean_delta - sem_delta, mean_delta + sem_delta,
                        color="#9467bd", alpha=0.25)
        ax.axhline(0, color="gray", linewidth=0.8, linestyle="--", alpha=0.6)

        # Shade regions where task > rest
        ax.fill_between(all_common_k, 0, mean_delta,
                        where=mean_delta > 0, color="#2ca02c", alpha=0.1)
        ax.fill_between(all_common_k, 0, mean_delta,
                        where=mean_delta < 0, color="#d62728", alpha=0.1)

        _style_ax(ax, ylabel=r"$\Delta$S(k) [task $-$ rest]")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    path = OUTPUT_DIR / "fig6_delta_segregation_task_minus_rest.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")


# ---------------------------------------------------------------------------
# CSV export + FINDINGS.md
# ---------------------------------------------------------------------------
def export_csv(records):
    df = pd.DataFrame(records)
    csv_path = OUTPUT_DIR / "segregation_profiles.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved {csv_path} ({len(df)} rows)")
    return df


def write_findings(profiles, df):
    """Write FINDINGS.md summarizing the key results."""
    lines = [
        "# Scale-Dependent Segregation-Integration Analysis",
        "",
        "## Method",
        "For each hierarchical level k (number of clusters from the LRG dendrogram),",
        "we compute within- and between-community MSC and a segregation index:",
        "",
        "    S(k) = (within_FC - between_FC) / (within_FC + between_FC)",
        "",
        f"- Patients: {', '.join(PATIENTS)}",
        f"- Bands: {', '.join(BANDS)}",
        f"- Phases: {', '.join(ALL_PHASES)}",
        f"- k range: {K_RANGE.start}--{K_RANGE.stop - 1}",
        f"- FC method: MSC (dense, nperseg=4096)",
        "",
        "## Key Findings",
        "",
    ]

    # 1. Scale of max segregation per band/phase
    lines.append("### Scale of Maximum Segregation (k*)")
    lines.append("")
    lines.append("| Band | Phase | k* (mean +/- SEM) | max S(k*) |")
    lines.append("|------|-------|-------------------|-----------|")

    for band in BANDS:
        for phase in ALL_PHASES:
            k_stars = []
            max_segs = []
            for patient in PATIENTS:
                if phase not in profiles.get(patient, {}).get(band, {}):
                    continue
                p = profiles[patient][band][phase]
                valid = ~np.isnan(p["segregation"])
                if np.any(valid):
                    best_idx = np.nanargmax(p["segregation"])
                    k_stars.append(p["k_values"][best_idx])
                    max_segs.append(p["segregation"][best_idx])
            if k_stars:
                km = np.mean(k_stars)
                ks = np.std(k_stars) / np.sqrt(len(k_stars))
                sm = np.mean(max_segs)
                tex = BRAIN_BAND_TEX_DICT[band].replace("$", "")
                lines.append(
                    f"| {band} | {PHASE_DISPLAY[phase]} | "
                    f"{km:.1f} +/- {ks:.1f} | {sm:.4f} |"
                )
    lines.append("")

    # 2. Rest vs task comparison
    lines.append("### Rest vs Task Comparison")
    lines.append("")

    rest_segs = []
    task_segs = []
    rest_k_stars = []
    task_k_stars = []
    for patient in PATIENTS:
        for band in BANDS:
            for phase in ALL_PHASES:
                if phase not in profiles.get(patient, {}).get(band, {}):
                    continue
                p = profiles[patient][band][phase]
                valid = ~np.isnan(p["segregation"])
                if np.any(valid):
                    max_s = np.nanmax(p["segregation"])
                    k_star = p["k_values"][np.nanargmax(p["segregation"])]
                    if CONDITION_MAP[phase] == "rest":
                        rest_segs.append(max_s)
                        rest_k_stars.append(k_star)
                    else:
                        task_segs.append(max_s)
                        task_k_stars.append(k_star)

    if rest_segs and task_segs:
        lines.append(
            f"- **Rest**: max S = {np.mean(rest_segs):.4f} +/- {np.std(rest_segs)/np.sqrt(len(rest_segs)):.4f}, "
            f"k* = {np.mean(rest_k_stars):.1f} +/- {np.std(rest_k_stars)/np.sqrt(len(rest_k_stars)):.1f}"
        )
        lines.append(
            f"- **Task**: max S = {np.mean(task_segs):.4f} +/- {np.std(task_segs)/np.sqrt(len(task_segs)):.4f}, "
            f"k* = {np.mean(task_k_stars):.1f} +/- {np.std(task_k_stars)/np.sqrt(len(task_k_stars)):.1f}"
        )
        diff = np.mean(task_segs) - np.mean(rest_segs)
        direction = "higher" if diff > 0 else "lower"
        lines.append(
            f"- Task segregation is **{direction}** than rest by {abs(diff):.4f}"
        )
    lines.append("")

    # 3. Band ranking by max segregation
    lines.append("### Band Ranking by Maximum Segregation")
    lines.append("")
    band_means = {}
    for band in BANDS:
        vals = []
        for patient in PATIENTS:
            for phase in ALL_PHASES:
                if phase not in profiles.get(patient, {}).get(band, {}):
                    continue
                p = profiles[patient][band][phase]
                if np.any(~np.isnan(p["segregation"])):
                    vals.append(np.nanmax(p["segregation"]))
        if vals:
            band_means[band] = np.mean(vals)
    sorted_bands = sorted(band_means.items(), key=lambda x: x[1], reverse=True)
    for rank, (band, mean_s) in enumerate(sorted_bands, 1):
        lines.append(f"{rank}. **{band}**: mean max S = {mean_s:.4f}")
    lines.append("")

    # 4. Cross-condition similarity
    lines.append("### Profile Similarity Within vs Between Conditions")
    lines.append("")
    lines.append("To assess whether S(k) profiles are more similar within-condition")
    lines.append("(rest-rest, task-task) vs cross-condition (rest-task), we compute")
    lines.append("profile correlations:")
    lines.append("")

    within_corrs = []
    between_corrs = []
    for patient in PATIENTS:
        for band in BANDS:
            bd = profiles.get(patient, {}).get(band, {})
            available = [ph for ph in ALL_PHASES if ph in bd]
            for i, ph1 in enumerate(available):
                for ph2 in available[i + 1:]:
                    p1 = bd[ph1]["segregation"]
                    p2 = bd[ph2]["segregation"]
                    # Align lengths
                    min_len = min(len(p1), len(p2))
                    v1 = p1[:min_len]
                    v2 = p2[:min_len]
                    valid = ~(np.isnan(v1) | np.isnan(v2))
                    if np.sum(valid) > 3:
                        r = np.corrcoef(v1[valid], v2[valid])[0, 1]
                        same_cond = CONDITION_MAP[ph1] == CONDITION_MAP[ph2]
                        if same_cond:
                            within_corrs.append(r)
                        else:
                            between_corrs.append(r)

    if within_corrs and between_corrs:
        lines.append(
            f"- **Within-condition** correlation: r = {np.mean(within_corrs):.3f} +/- {np.std(within_corrs)/np.sqrt(len(within_corrs)):.3f} (n={len(within_corrs)})"
        )
        lines.append(
            f"- **Between-condition** correlation: r = {np.mean(between_corrs):.3f} +/- {np.std(between_corrs)/np.sqrt(len(between_corrs)):.3f} (n={len(between_corrs)})"
        )
    lines.append("")

    findings_path = OUTPUT_DIR / "FINDINGS.md"
    findings_path.write_text("\n".join(lines))
    print(f"Saved {findings_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 60)
    print("Scale-Dependent Segregation-Integration Analysis")
    print("=" * 60)

    print("\n[1/3] Loading data and computing segregation profiles...")
    profiles, records = load_all_profiles()

    print(f"\n[2/3] Exporting CSV...")
    df = export_csv(records)

    print(f"\n[3/3] Generating figures...")
    plot_fig1(profiles)
    plot_fig2(profiles)
    plot_fig3(profiles)
    plot_fig4(profiles)
    plot_fig5(profiles)
    plot_fig6(profiles)

    print(f"\nWriting findings...")
    write_findings(profiles, df)

    print("\nDone!")


if __name__ == "__main__":
    main()
