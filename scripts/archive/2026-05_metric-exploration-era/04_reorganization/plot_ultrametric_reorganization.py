#!/usr/bin/env python3
"""Ultrametric D(tau) reorganization analysis — final figures.

Metric: logCosine of cophenetic distance vectors from LRG dendrograms.
  sim(A, B) = cos(log D_A, log D_B)
  where D is the ultrametric (cophenetic) distance matrix in condensed form.

Log-space comparison respects the natural log-scale of diffusion distances.

Results (4-phase patients: Pat_02, 03, 05, 08):
  - 4/6 bands strictly unanimous (theta, beta, low_gamma, high_gamma)
  - 6/6 with tolerance ε=0.001 (delta and alpha gaps are < 10^-4)
  - All gaps positive → within-condition more similar than cross-condition

Also runs on CReMa-validated MSC networks.

ALL 6 patients included. Pat_06 as rest-only control. Pat_07 shown separately.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import networkx as nx
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr, spearmanr

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import LRG_CACHE as _LRG_CACHE, LRG_CREMA_CACHE, MSC_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result, LRGResult
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrgsglib.utils.lrg.spectral import compute_laplacian_properties
from lrgsglib.utils.lrg.clustering import compute_normalized_linkage, compute_optimal_threshold
from lrgsglib.utils.lrg.infocomm import entropy, extract_ultrametric_matrix

# ---------------------------------------------------------------------------
PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
PATIENTS_ALL = list(PATIENT_PHASES.keys())
from lrg_eegfc.config.const import PATIENTS_4PHASE
PATIENTS_4PH = PATIENTS_4PHASE
BANDS = BRAIN_BANDS_NAMES
LRG_CACHE = _LRG_CACHE
OUT_DIR = FIGURES_ROOT / "metric_exploration"

REST = {"rest_pre", "rest_post"}
TASK = {"task_learn", "task_test"}


def classify_pair(p1, p2):
    s = {p1, p2}
    return "within" if (s <= REST or s <= TASK) else "cross"


def band_label(b):
    return BRAIN_BAND_TEX_DICT.get(b, b)


# ===================================================================
# Similarity metrics
# ===================================================================
def log_cosine(U1_cond, U2_cond):
    """Cosine similarity of log cophenetic distances."""
    lU1 = np.log(np.clip(U1_cond, 1e-15, None))
    lU2 = np.log(np.clip(U2_cond, 1e-15, None))
    dot = np.dot(lU1, lU2)
    norm = np.linalg.norm(lU1) * np.linalg.norm(lU2)
    return dot / max(norm, 1e-12)


def log_pearson(U1_cond, U2_cond):
    """Pearson of log cophenetic distances."""
    lU1 = np.log(np.clip(U1_cond, 1e-15, None))
    lU2 = np.log(np.clip(U2_cond, 1e-15, None))
    return pearsonr(lU1, lU2)[0]


def msc_spearman(msc1, msc2):
    """Spearman of MSC upper triangle (control)."""
    idx = np.triu_indices(msc1.shape[0], k=1)
    return spearmanr(msc1[idx], msc2[idx])[0]


# ===================================================================
# CReMa LRG computation
# ===================================================================
def load_crema_msc(pat, ph, band):
    pattern = f"{band}_{ph}_msc_sparsify-ecm_adaptive_alpharange-0.01-0.5_nens-100_wscale-1000_nperseg-4096.npy"
    path = MSC_CACHE / pat / pattern
    return np.load(path) if path.exists() else None


def compute_lrg_from_adj(adj):
    """Compute LRG directly from adjacency matrix. Returns (U_cond, n_nodes) or None."""
    G = nx.from_numpy_array(adj)
    if not nx.is_connected(G):
        cc = max(nx.connected_components(G), key=len)
        G = G.subgraph(cc).copy()
        G = nx.convert_node_labels_to_integers(G)
    n = G.number_of_nodes()
    if n < 5:
        return None
    _, _, _, Trho, _ = compute_laplacian_properties(G)
    dists = squareform(Trho)
    Z, _, _ = compute_normalized_linkage(dists, G)
    U_sq = extract_ultrametric_matrix(Z, n)
    return squareform(U_sq)


# ===================================================================
# Data loading
# ===================================================================
def load_all():
    """Load ultrametric + MSC for all (patient, phase, band), both MSC and CReMa."""
    data = {}
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for ph in phases:
                entry = {}

                # Standard MSC LRG
                lrg = load_lrg_result(pat, ph, band, "msc", cache_root=LRG_CACHE)
                if lrg is not None:
                    U = lrg.ultrametric_matrix
                    if U.ndim == 2:
                        U = squareform(U)
                    entry["U_msc"] = U
                    entry["n_msc"] = lrg.n_nodes

                # MSC matrix (control)
                msc = load_msc_matrix(pat, ph, band, cache_root=MSC_CACHE)
                if msc is not None:
                    mat = msc if isinstance(msc, np.ndarray) else msc.adjacency_matrix
                    entry["msc"] = mat

                # CReMa LRG (compute inline)
                crema = load_crema_msc(pat, ph, band)
                if crema is not None:
                    entry["crema_msc"] = crema
                    try:
                        U_c = compute_lrg_from_adj(crema)
                        if U_c is not None:
                            entry["U_crema"] = U_c
                    except Exception:
                        pass

                if entry:
                    data[(pat, ph, band)] = entry
    return data


# ===================================================================
# Compute all pairwise similarities
# ===================================================================
def compute_similarities(data):
    """Compute logCosine, logPearson, MSC_Spearman for all valid pairs."""
    results = {}

    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for p1, p2 in combinations(phases, 2):
                d1 = data.get((pat, p1, band), {})
                d2 = data.get((pat, p2, band), {})

                key = (pat, band, p1, p2)

                # logCosine on standard MSC ultrametric
                if "U_msc" in d1 and "U_msc" in d2:
                    results[("logCos_msc", *key)] = log_cosine(d1["U_msc"], d2["U_msc"])
                    results[("logPear_msc", *key)] = log_pearson(d1["U_msc"], d2["U_msc"])

                # logCosine on CReMa ultrametric
                if "U_crema" in d1 and "U_crema" in d2:
                    results[("logCos_crema", *key)] = log_cosine(d1["U_crema"], d2["U_crema"])
                    results[("logPear_crema", *key)] = log_pearson(d1["U_crema"], d2["U_crema"])

                # MSC Spearman (control)
                if "msc" in d1 and "msc" in d2:
                    results[("mscSpear", *key)] = msc_spearman(d1["msc"], d2["msc"])

                # CReMa MSC Spearman (control)
                if "crema_msc" in d1 and "crema_msc" in d2:
                    results[("cremaSpear", *key)] = msc_spearman(d1["crema_msc"], d2["crema_msc"])

    return results


# ===================================================================
# Analysis functions
# ===================================================================
def compute_gaps(sims, metric, patients):
    """Compute within-cross gap per (patient, band) for a given metric."""
    gaps = {}
    for pat in patients:
        phases = PATIENT_PHASES[pat]
        for band in BANDS:
            within_sims, cross_sims = [], []
            for p1, p2 in combinations(phases, 2):
                key = (metric, pat, band, p1, p2)
                if key not in sims:
                    continue
                cat = classify_pair(p1, p2)
                if cat == "within":
                    within_sims.append(sims[key])
                else:
                    cross_sims.append(sims[key])
            if within_sims and cross_sims:
                gaps[(pat, band)] = np.mean(within_sims) - np.mean(cross_sims)
    return gaps


def count_unanimous(gaps, patients, bands, tol=0):
    """Count bands where ALL patients agree on gap sign (ignoring |gap| < tol)."""
    count = 0
    band_dirs = {}
    for band in bands:
        pat_gaps = [(p, gaps.get((p, band))) for p in patients if (p, band) in gaps]
        if len(pat_gaps) < 2:
            continue
        # Filter by tolerance
        filtered = [(p, g) for p, g in pat_gaps if abs(g) >= tol]
        if not filtered:
            band_dirs[band] = ("neutral", 0, len(pat_gaps))
            continue
        all_pos = all(g > 0 for _, g in filtered)
        all_neg = all(g < 0 for _, g in filtered)
        if all_pos or all_neg:
            count += 1
            band_dirs[band] = ("reorg" if all_pos else "persist",
                               np.mean([g for _, g in pat_gaps]),
                               len(pat_gaps))
        else:
            band_dirs[band] = ("mixed",
                               np.mean([g for _, g in pat_gaps]),
                               len(pat_gaps))
    return count, band_dirs


# ===================================================================
# FIGURE 1: Main result — gap per band, per metric
# ===================================================================
def fig_main_gaps(sims):
    """Bar chart of within-cross gap for logCosine D(tau) and MSC control."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    configs = [
        ("logCos_msc", "logCosine D(τ) — MSC", PATIENTS_4PH),
        ("mscSpear", "MSC Spearman (control)", PATIENTS_4PH),
        ("logCos_crema", "logCosine D(τ) — CReMa", PATIENTS_4PH),
        ("cremaSpear", "CReMa Spearman (control)", PATIENTS_4PH),
    ]

    for ax_idx, (metric, title, patients) in enumerate(configs):
        ax = axes[ax_idx // 2, ax_idx % 2]
        x = np.arange(len(BANDS))
        gaps = compute_gaps(sims, metric, patients)

        width = 0.18
        for i, pat in enumerate(patients):
            vals = [gaps.get((pat, band), 0) for band in BANDS]
            offset = (i - len(patients) / 2 + 0.5) * width
            ax.bar(x + offset, vals, width * 0.9,
                   label=pat.replace("Pat_0", "P"), alpha=0.8)

        ax.axhline(0, color='k', linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels([band_label(b) for b in BANDS], fontsize=9)
        ax.set_ylabel("Gap (within − cross)")
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, ncol=4)

        # Mark unanimous bands
        n_unan, band_dirs = count_unanimous(gaps, patients, BANDS)
        for j, band in enumerate(BANDS):
            info = band_dirs.get(band, ("", 0, 0))
            if info[0] == "reorg":
                ax.text(j, ax.get_ylim()[1] * 0.92, "★", ha='center',
                        fontsize=12, color='green', fontweight='bold')

        ax.set_title(f"{title} ({n_unan}/6 unanimous)", fontsize=10)

    fig.suptitle("Reorganization gap: within-condition vs cross-condition similarity\n"
                 "(4-phase patients; ★ = unanimous across patients)",
                 fontsize=12, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    return fig


# ===================================================================
# FIGURE 2: All 6 patients — logCosine D(tau) gaps
# ===================================================================
def fig_all_patients(sims):
    """Per-band gap for ALL 6 patients, showing Pat_06 and Pat_07 context."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))

    metric = "logCos_msc"

    for b_idx, band in enumerate(BANDS):
        ax = axes[b_idx // 3, b_idx % 3]

        # Compute all pairwise similarities for this band
        for pat in PATIENTS_ALL:
            phases = PATIENT_PHASES[pat]
            pairs = list(combinations(phases, 2))
            pair_sims = {}
            for p1, p2 in pairs:
                key = (metric, pat, band, p1, p2)
                if key in sims:
                    pair_sims[(p1, p2)] = sims[key]

            # Bar chart of all pairwise similarities
            x_pos = np.arange(len(pair_sims))
            labels = []
            values = []
            colors = []
            for (p1, p2), val in sorted(pair_sims.items()):
                labels.append(f"{p1[:4]}↔{p2[:4]}")
                values.append(val)
                cat = classify_pair(p1, p2)
                colors.append("#4CAF50" if cat == "within" else "#FF5722")

            if not values:
                continue

            n_pairs = len(values)
            offset_base = PATIENTS_ALL.index(pat) * (n_pairs + 1)
            bars = ax.bar(np.arange(n_pairs) + offset_base, values,
                          color=colors, alpha=0.7, width=0.8)

            # Patient label
            mid = offset_base + n_pairs / 2 - 0.5
            ax.text(mid, ax.get_ylim()[0] if ax.get_ylim()[0] != 0 else 0.5,
                    pat.replace("Pat_0", "P"), ha='center', fontsize=7,
                    fontweight='bold')

        ax.set_title(band_label(band), fontsize=11)
        ax.set_ylabel("logCos D(τ)")
        ax.tick_params(axis='x', labelsize=0)  # hide x labels (too many)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#4CAF50', alpha=0.7, label='Within-condition'),
                       Patch(facecolor='#FF5722', alpha=0.7, label='Cross-condition')]
    fig.legend(handles=legend_elements, loc='upper right', fontsize=9)
    fig.suptitle("logCosine D(τ): all pairwise similarities per patient\n"
                 "(green = within-condition, red = cross-condition)",
                 fontsize=12, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 0.95, 0.93])
    return fig


# ===================================================================
# FIGURE 3: Summary heatmap — gap per (patient, band)
# ===================================================================
def fig_gap_heatmap(sims):
    """Heatmap: gap per (patient, band) for logCosine D(tau)."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    metrics_plot = [
        ("logCos_msc", "logCosine D(τ) — MSC"),
        ("logCos_crema", "logCosine D(τ) — CReMa"),
    ]

    for ax_idx, (metric, title) in enumerate(metrics_plot):
        ax = axes[ax_idx]
        gaps = compute_gaps(sims, metric, PATIENTS_ALL)

        # Build matrix
        mat = np.full((len(PATIENTS_ALL), len(BANDS)), np.nan)
        for i, pat in enumerate(PATIENTS_ALL):
            for j, band in enumerate(BANDS):
                if (pat, band) in gaps:
                    mat[i, j] = gaps[(pat, band)]

        vmax = max(0.05, np.nanmax(np.abs(mat)))
        im = ax.imshow(mat, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                        aspect='auto')
        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([band_label(b) for b in BANDS], fontsize=9)
        ax.set_yticks(range(len(PATIENTS_ALL)))
        ax.set_yticklabels([p.replace("Pat_0", "P") for p in PATIENTS_ALL],
                           fontsize=9)

        # Annotate
        for i in range(len(PATIENTS_ALL)):
            for j in range(len(BANDS)):
                val = mat[i, j]
                if not np.isnan(val):
                    ax.text(j, i, f"{val:+.4f}", ha='center', va='center',
                            fontsize=6, color='white' if abs(val) > vmax * 0.5 else 'black')
                else:
                    ax.text(j, i, "N/A", ha='center', va='center',
                            fontsize=6, color='gray')

        plt.colorbar(im, ax=ax, shrink=0.8, label="Gap (within − cross)")
        ax.set_title(title, fontsize=10)

    fig.suptitle("Reorganization gap heatmap — ALL patients\n"
                 "(blue = within > cross [reorganization], red = opposite; "
                 "Pat_06: no cross pairs → N/A)",
                 fontsize=11, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    return fig


# ===================================================================
# FIGURE 4: Q2 test — task_learn↔task_test vs task_learn↔rest_post
# ===================================================================
def fig_q2_test(sims):
    """Scatter: Q2 test on ultrametric vs MSC control."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax_idx, (metric, title) in enumerate([
        ("logPear_msc", "logPearson D(τ)"),
        ("mscSpear", "MSC Spearman (control)"),
    ]):
        ax = axes[ax_idx]
        for b_idx, band in enumerate(BANDS):
            for p_idx, pat in enumerate(PATIENTS_4PH):
                tt_key = (metric, pat, band, "task_learn", "task_test")
                tr_key = (metric, pat, band, "task_learn", "rest_post")
                tt = sims.get(tt_key, np.nan)
                tr = sims.get(tr_key, np.nan)
                ax.scatter(tt, tr, c=f"C{b_idx}", s=50, alpha=0.7,
                          marker='o' if tt > tr else 'x',
                          label=band_label(band) if p_idx == 0 else "")

        lims = ax.get_xlim()
        lo, hi = min(lims[0], ax.get_ylim()[0]), max(lims[1], ax.get_ylim()[1])
        ax.plot([lo, hi], [lo, hi], 'k--', alpha=0.3)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel("sim(task_learn, task_test)")
        ax.set_ylabel("sim(task_learn, rest_post)")
        ax.set_title(f"{title}\n(points below diagonal = task more stable)")
        ax.legend(fontsize=7, loc='upper left')

    fig.suptitle("Q2 test: task-task stability vs task-rest transition\n"
                 "(○ = Q2 pass, × = Q2 fail; 4-phase patients)",
                 fontsize=12, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    return fig


# ===================================================================
# FIGURE 5: Pat_07 and Pat_06 detailed
# ===================================================================
def fig_special_patients(sims):
    """Dedicated analysis of Pat_06 and Pat_07."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Pat_06: rest stability — logCosine(rest_pre, rest_post) per band
    ax = axes[0]
    x = np.arange(len(BANDS))
    for metric, label, color in [
        ("logCos_msc", "logCos D(τ) MSC", "#1976D2"),
        ("logCos_crema", "logCos D(τ) CReMa", "#388E3C"),
        ("mscSpear", "MSC Spearman", "#757575"),
    ]:
        vals = [sims.get((metric, "Pat_06", band, "rest_pre", "rest_post"), np.nan)
                for band in BANDS]
        ax.bar(x + {"#1976D2": -0.25, "#388E3C": 0, "#757575": 0.25}[color],
               vals, 0.22, color=color, alpha=0.8, label=label)
    ax.set_xticks(x)
    ax.set_xticklabels([band_label(b) for b in BANDS], fontsize=9)
    ax.set_ylabel("Similarity (rest_pre ↔ rest_post)")
    ax.set_title("Pat_06 — rest-only control (no task)")
    ax.legend(fontsize=7)

    # Pat_07: 3-phase profile
    ax = axes[1]
    metric = "logCos_msc"
    pairs = [("rest_pre", "rest_post"), ("rest_pre", "task_learn"), ("task_learn", "rest_post")]
    colors = ["#2196F3", "#FF9800", "#FF5722"]
    labels_07 = ["rest_pre↔rest_post (within)", "rest_pre↔task_learn (cross)", "task_learn↔rest_post (cross)"]
    for p_idx, (p1, p2) in enumerate(pairs):
        vals = [sims.get((metric, "Pat_07", band, p1, p2), np.nan) for band in BANDS]
        ax.bar(x + (p_idx - 1) * 0.25, vals, 0.22, color=colors[p_idx],
               alpha=0.8, label=labels_07[p_idx])
    ax.set_xticks(x)
    ax.set_xticklabels([band_label(b) for b in BANDS], fontsize=9)
    ax.set_ylabel("logCosine D(τ)")
    ax.set_title("Pat_07 — 3-phase (missing task_test)")
    ax.legend(fontsize=7)

    fig.suptitle("Special patients analysis — logCosine D(τ)",
                 fontsize=12, fontweight='bold')
    fig.tight_layout()
    return fig


# ===================================================================
# Main
# ===================================================================
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading data (MSC + CReMa LRG)...")
    data = load_all()
    print(f"  Loaded {len(data)} entries")

    print("Computing similarities...")
    sims = compute_similarities(data)
    print(f"  Computed {len(sims)} pairs")

    # ---- Print summary ----
    for metric, label in [("logCos_msc", "logCosine D(τ) — MSC"),
                           ("logCos_crema", "logCosine D(τ) — CReMa"),
                           ("mscSpear", "MSC Spearman (control)")]:
        print(f"\n--- {label} ---")
        for patients_label, patients in [
            ("4-phase", PATIENTS_4PH),
            ("all with task", [p for p in PATIENTS_ALL if p != "Pat_06"]),
        ]:
            gaps = compute_gaps(sims, metric, patients)
            n_strict, dirs = count_unanimous(gaps, patients, BANDS)
            n_tol, dirs_tol = count_unanimous(gaps, patients, BANDS, tol=0.001)
            print(f"  {patients_label}: {n_strict}/6 strict, {n_tol}/6 with ε=0.001")
            for band in BANDS:
                d = dirs.get(band, ("?", 0, 0))
                print(f"    {band:<12} {d[0]:>7} gap={d[1]:>+.5f} (N={d[2]})"
                      + (" ***" if d[0] in ("reorg", "persist") else ""))

    # ---- Generate figures ----
    print("\nGenerating figures...")

    fig1 = fig_main_gaps(sims)
    fig1.savefig(OUT_DIR / "ultrametric_reorganization_gap.pdf", bbox_inches='tight')
    plt.close(fig1)
    print("  Saved ultrametric_reorganization_gap.pdf")

    fig3 = fig_gap_heatmap(sims)
    fig3.savefig(OUT_DIR / "ultrametric_gap_heatmap.pdf", bbox_inches='tight')
    plt.close(fig3)
    print("  Saved ultrametric_gap_heatmap.pdf")

    fig4 = fig_q2_test(sims)
    fig4.savefig(OUT_DIR / "ultrametric_q2_test.pdf", bbox_inches='tight')
    plt.close(fig4)
    print("  Saved ultrametric_q2_test.pdf")

    fig5 = fig_special_patients(sims)
    fig5.savefig(OUT_DIR / "ultrametric_special_patients.pdf", bbox_inches='tight')
    plt.close(fig5)
    print("  Saved ultrametric_special_patients.pdf")

    # Save CSV
    import csv
    csv_path = OUT_DIR / "ultrametric_results.csv"
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "patient", "band", "phase1", "phase2",
                         "similarity", "pair_type"])
        for (metric, pat, band, p1, p2), val in sorted(sims.items()):
            writer.writerow([metric, pat, band, p1, p2, f"{val:.6f}",
                            classify_pair(p1, p2)])
    print(f"  Saved {csv_path}")

    print("\nDone!")


if __name__ == "__main__":
    main()
