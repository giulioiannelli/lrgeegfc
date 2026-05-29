#!/usr/bin/env python3
"""Scale-Resolved Modularity Profiles via LRG Dendrogram Partitions.

For each patient/band/phase, cut the LRG dendrogram at many thresholds,
compute modularity Q of the resulting partition using the original MSC
adjacency matrix, and analyze how Q varies across hierarchical scales.
"""

import warnings
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms.community.quality import modularity
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    BRAIN_BANDS_NAMES,
    PHASE_LABELS,
    BRAIN_BAND_TEX_DICT,
    PATIENTS_4PHASE,
)
from lrg_eegfc.config.paths import FIGURES_ROOT

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Configuration ──────────────────────────────────────────────────────
FC_METHOD = "msc"
N_THRESHOLDS = 40
ALL_PATIENTS = PATIENTS_4PHASE  # all 4 phases
BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)

PHASE_COLORS = {
    "rest_pre": "#1f77b4",
    "task_learn": "#ff7f0e",
    "task_test": "#2ca02c",
    "rest_post": "#d62728",
}
PHASE_LABELS_NICE = {
    "rest_pre": "rest_pre",
    "task_learn": "task_learn",
    "task_test": "task_test",
    "rest_post": "rest_post",
}

OUT_DIR = FIGURES_ROOT / "multiscale_investigation" / "modularity_profiles"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ── Helper: compute modularity profile ────────────────────────────────
def compute_modularity_profile(linkage_matrix, adj_matrix, n_nodes, n_thresholds=40):
    """Compute Q(scale) by cutting dendrogram at evenly-spaced distance thresholds.

    Returns
    -------
    n_clusters_arr : ndarray, shape (n_points,)
        Number of clusters at each threshold.
    q_arr : ndarray, shape (n_points,)
        Modularity Q at each threshold.
    thresholds : ndarray, shape (n_thresholds,)
        The distance thresholds used.
    """
    # Extract the giant component from the MSC matrix to build the graph
    G_full = nx.from_numpy_array(adj_matrix)
    giant_nodes = max(nx.connected_components(G_full), key=len)
    G = G_full.subgraph(sorted(giant_nodes)).copy()

    # Map giant component nodes to 0..n_gc-1
    gc_nodes = sorted(G.nodes())
    n_gc = len(gc_nodes)

    # If n_gc != n_nodes from LRG, we need to be careful
    # The LRG was computed on the giant component, so n_nodes should match n_gc
    if n_gc != n_nodes:
        print(f"  WARNING: giant component size {n_gc} != LRG n_nodes {n_nodes}")
        # Use the LRG n_nodes - the dendrogram was built on that
        # We need to extract matching subgraph
        if n_gc > n_nodes:
            # LRG used a smaller giant component - just use first n_nodes of gc_nodes
            gc_nodes = gc_nodes[:n_nodes]
            G = G_full.subgraph(gc_nodes).copy()
            n_gc = n_nodes

    # Relabel to 0..n_gc-1 for consistency with fcluster output
    mapping = {old: new for new, old in enumerate(gc_nodes)}
    G = nx.relabel_nodes(G, mapping)

    # Get linkage distance range
    distances = linkage_matrix[:, 2]
    d_min, d_max = distances.min(), distances.max()
    # Avoid endpoints where we get 1 cluster or n_nodes clusters
    eps = (d_max - d_min) * 0.01
    thresholds = np.linspace(d_min + eps, d_max - eps, n_thresholds)

    n_clusters_list = []
    q_list = []

    for t in thresholds:
        labels = fcluster(linkage_matrix, t=t, criterion="distance")
        unique_labels = np.unique(labels)
        n_clust = len(unique_labels)

        # Build community sets
        communities = []
        for lab in unique_labels:
            members = set(np.where(labels == lab)[0])
            # Only include nodes that exist in the graph
            members = members & set(G.nodes())
            if members:
                communities.append(frozenset(members))

        # Ensure all graph nodes are covered
        covered = set()
        for c in communities:
            covered |= set(c)
        uncovered = set(G.nodes()) - covered
        if uncovered:
            # Add uncovered nodes as singletons
            for node in uncovered:
                communities.append(frozenset([node]))

        if len(communities) <= 1:
            # Modularity is 0 for a single community
            q_list.append(0.0)
        else:
            try:
                q = modularity(G, communities)
                q_list.append(q)
            except Exception:
                q_list.append(np.nan)

        n_clusters_list.append(n_clust)

    return (
        np.array(n_clusters_list),
        np.array(q_list),
        thresholds,
    )


# ── Step 1: Compute all profiles ──────────────────────────────────────
print("=" * 70)
print("SCALE-RESOLVED MODULARITY PROFILES")
print("=" * 70)

# Storage: results[patient][band][phase] = (n_clusters, q_values, thresholds, optimal_threshold)
results = {}
all_rows = []  # for CSV

for patient in ALL_PATIENTS:
    results[patient] = {}
    for band in BANDS:
        results[patient][band] = {}
        for phase in PHASES:
            lrg = load_lrg_result(patient, phase, band, FC_METHOD)
            msc = load_msc_matrix(patient, phase, band)

            if lrg is None or msc is None:
                print(f"  SKIP {patient} {band} {phase}: missing data")
                results[patient][band][phase] = None
                continue

            n_clust, q_vals, thresholds = compute_modularity_profile(
                lrg.linkage_matrix, msc, lrg.n_nodes, N_THRESHOLDS
            )

            # Also compute Q at the LRG optimal threshold
            labels_opt = fcluster(
                lrg.linkage_matrix, t=lrg.optimal_threshold, criterion="distance"
            )
            n_clust_opt = len(np.unique(labels_opt))

            results[patient][band][phase] = {
                "n_clusters": n_clust,
                "q_values": q_vals,
                "thresholds": thresholds,
                "optimal_threshold": lrg.optimal_threshold,
                "n_clust_opt": n_clust_opt,
                "n_nodes": lrg.n_nodes,
            }

            # Find max Q
            valid = ~np.isnan(q_vals)
            if valid.any():
                idx_max = np.nanargmax(q_vals)
                q_max = q_vals[idx_max]
                n_clust_max = n_clust[idx_max]
                t_max = thresholds[idx_max]
            else:
                q_max = np.nan
                n_clust_max = np.nan
                t_max = np.nan

            all_rows.append({
                "patient": patient,
                "band": band,
                "phase": phase,
                "q_max": q_max,
                "n_clusters_at_q_max": n_clust_max,
                "threshold_at_q_max": t_max,
                "optimal_threshold": lrg.optimal_threshold,
                "n_clusters_at_optimal": n_clust_opt,
                "n_nodes": lrg.n_nodes,
            })

            print(f"  {patient} {band:>12s} {phase:>10s}: "
                  f"Q_max={q_max:.4f} @ k={n_clust_max:.0f}, "
                  f"LRG_opt @ k={n_clust_opt}")

df = pd.DataFrame(all_rows)
df.to_csv(OUT_DIR / "modularity_profiles_summary.csv", index=False)
print(f"\nSaved summary CSV: {OUT_DIR / 'modularity_profiles_summary.csv'}")


# ── Figure 1: Q(n_clusters) for Pat_02, one subplot per band ─────────
print("\n--- Figure 1: Pat_02 Q(n_clusters) profiles ---")
fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharex=False, sharey=True)
axes = axes.flatten()

for i, band in enumerate(BANDS):
    ax = axes[i]
    for phase in PHASES:
        data = results["Pat_02"][band][phase]
        if data is None:
            continue
        ax.plot(
            data["n_clusters"],
            data["q_values"],
            color=PHASE_COLORS[phase],
            label=PHASE_LABELS_NICE[phase],
            linewidth=1.8,
            alpha=0.85,
        )
        # Mark max Q
        idx_max = np.nanargmax(data["q_values"])
        ax.plot(
            data["n_clusters"][idx_max],
            data["q_values"][idx_max],
            "o",
            color=PHASE_COLORS[phase],
            markersize=7,
            zorder=5,
        )

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=13)
    ax.set_xlabel("Number of clusters", fontsize=10)
    if i % 3 == 0:
        ax.set_ylabel("Modularity Q", fontsize=10)
    ax.grid(True, alpha=0.3)
    if i == 0:
        ax.legend(fontsize=9, loc="best")

fig.suptitle("Scale-Resolved Modularity Profiles - Pat_02", fontsize=15, y=0.98)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUT_DIR / "fig1_pat02_modularity_profiles.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: fig1_pat02_modularity_profiles.png")


# ── Figure 2: Scale of max Q across phases, averaged over patients ───
print("\n--- Figure 2: Scale of max Q (averaged over patients) ---")
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
axes = axes.flatten()

for i, band in enumerate(BANDS):
    ax = axes[i]
    phase_means = []
    phase_stds = []
    for phase in PHASES:
        vals = []
        for patient in ALL_PATIENTS:
            data = results[patient][band][phase]
            if data is not None:
                idx_max = np.nanargmax(data["q_values"])
                vals.append(data["n_clusters"][idx_max])
        phase_means.append(np.mean(vals) if vals else np.nan)
        phase_stds.append(np.std(vals) if vals else np.nan)

    x = np.arange(len(PHASES))
    bars = ax.bar(
        x,
        phase_means,
        yerr=phase_stds,
        color=[PHASE_COLORS[p] for p in PHASES],
        edgecolor="black",
        linewidth=0.5,
        capsize=5,
        alpha=0.8,
    )
    ax.set_xticks(x)
    ax.set_xticklabels([PHASE_LABELS_NICE[p] for p in PHASES], fontsize=9, rotation=30)
    ax.set_ylabel("k at max Q", fontsize=10)
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=13)
    ax.grid(True, alpha=0.3, axis="y")

fig.suptitle(
    "Scale of Maximum Modularity (n_clusters) - Mean over patients",
    fontsize=14, y=0.98,
)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUT_DIR / "fig2_scale_max_q_by_phase.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: fig2_scale_max_q_by_phase.png")


# ── Figure 3: Max Q heatmap (phases x bands), averaged over patients ─
print("\n--- Figure 3: Max Q heatmap ---")
q_max_matrix = np.full((len(PHASES), len(BANDS)), np.nan)
for j, band in enumerate(BANDS):
    for i, phase in enumerate(PHASES):
        vals = []
        for patient in ALL_PATIENTS:
            data = results[patient][band][phase]
            if data is not None:
                vals.append(np.nanmax(data["q_values"]))
        if vals:
            q_max_matrix[i, j] = np.mean(vals)

fig, ax = plt.subplots(figsize=(10, 5))
im = ax.imshow(q_max_matrix, aspect="auto", cmap="YlOrRd", interpolation="nearest")
ax.set_xticks(range(len(BANDS)))
ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=12)
ax.set_yticks(range(len(PHASES)))
ax.set_yticklabels([PHASE_LABELS_NICE[p] for p in PHASES], fontsize=12)

# Annotate cells
for i in range(len(PHASES)):
    for j in range(len(BANDS)):
        val = q_max_matrix[i, j]
        if not np.isnan(val):
            ax.text(j, i, f"{val:.3f}", ha="center", va="center", fontsize=10,
                    color="white" if val > np.nanmedian(q_max_matrix) else "black")

cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label("Max Modularity Q", fontsize=11)
ax.set_title("Maximum Modularity Q (mean over patients)", fontsize=14)
fig.tight_layout()
fig.savefig(OUT_DIR / "fig3_max_q_heatmap.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: fig3_max_q_heatmap.png")


# ── Figure 4: Grand average Q(scale) - rest vs task ─────────────────
print("\n--- Figure 4: Grand average rest vs task ---")
REST_PHASES = ["rest_pre", "rest_post"]
TASK_PHASES = ["task_learn", "task_test"]

fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharey=True)
axes = axes.flatten()

for i, band in enumerate(BANDS):
    ax = axes[i]

    # Collect profiles for rest and task, interpolated to common n_clusters grid
    # Use n_clusters as x-axis; interpolate all profiles to a common grid
    all_n_clust = []
    for patient in ALL_PATIENTS:
        for phase in PHASES:
            data = results[patient][band][phase]
            if data is not None:
                all_n_clust.extend(data["n_clusters"])

    if not all_n_clust:
        continue

    k_min = max(2, min(all_n_clust))
    k_max = max(all_n_clust)
    k_grid = np.linspace(k_min, k_max, 50)

    for condition, cond_phases, color, label in [
        ("rest", REST_PHASES, "#1f77b4", "Rest (rest_pre + rest_post)"),
        ("task", TASK_PHASES, "#ff7f0e", "Task (task_learn + task_test)"),
    ]:
        interp_profiles = []
        for patient in ALL_PATIENTS:
            for phase in cond_phases:
                data = results[patient][band][phase]
                if data is not None:
                    # Sort by n_clusters for interpolation
                    sort_idx = np.argsort(data["n_clusters"])
                    nc = data["n_clusters"][sort_idx]
                    qv = data["q_values"][sort_idx]
                    # Interpolate
                    q_interp = np.interp(k_grid, nc, qv)
                    interp_profiles.append(q_interp)

        if interp_profiles:
            arr = np.array(interp_profiles)
            mean = np.nanmean(arr, axis=0)
            sem = np.nanstd(arr, axis=0) / np.sqrt(arr.shape[0])
            ax.plot(k_grid, mean, color=color, linewidth=2, label=label)
            ax.fill_between(k_grid, mean - sem, mean + sem, color=color, alpha=0.2)

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=13)
    ax.set_xlabel("Number of clusters", fontsize=10)
    if i % 3 == 0:
        ax.set_ylabel("Modularity Q", fontsize=10)
    ax.grid(True, alpha=0.3)
    if i == 0:
        ax.legend(fontsize=9, loc="best")

fig.suptitle("Grand Average Modularity Profiles: Rest vs Task", fontsize=15, y=0.98)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUT_DIR / "fig4_rest_vs_task_profiles.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: fig4_rest_vs_task_profiles.png")


# ── Figure 5: Per-patient Q profiles for alpha band ──────────────────
print("\n--- Figure 5: Per-patient alpha profiles ---")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

for idx, patient in enumerate(ALL_PATIENTS):
    ax = axes[idx]
    for phase in PHASES:
        data = results[patient]["alpha"][phase]
        if data is None:
            continue
        ax.plot(
            data["n_clusters"],
            data["q_values"],
            color=PHASE_COLORS[phase],
            label=PHASE_LABELS_NICE[phase],
            linewidth=1.8,
            alpha=0.85,
        )
        # Mark max Q
        idx_max = np.nanargmax(data["q_values"])
        ax.plot(
            data["n_clusters"][idx_max],
            data["q_values"][idx_max],
            "o",
            color=PHASE_COLORS[phase],
            markersize=7,
            zorder=5,
        )
        # Mark LRG optimal threshold
        # Find closest n_clusters to the optimal threshold
        opt_t = data["optimal_threshold"]
        dists = np.abs(data["thresholds"] - opt_t)
        idx_opt = np.argmin(dists)
        ax.axvline(
            data["n_clusters"][idx_opt],
            color=PHASE_COLORS[phase],
            linestyle="--",
            alpha=0.4,
            linewidth=1,
        )

    ax.set_title(f"{patient} - {BRAIN_BAND_TEX_DICT['alpha']}", fontsize=13)
    ax.set_xlabel("Number of clusters", fontsize=10)
    ax.set_ylabel("Modularity Q", fontsize=10)
    ax.grid(True, alpha=0.3)
    if idx == 0:
        ax.legend(fontsize=9, loc="best")

fig.suptitle(
    "Per-Patient Modularity Profiles (Alpha Band)\n"
    "Circles = max Q, dashed lines = LRG optimal threshold",
    fontsize=14, y=1.0,
)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT_DIR / "fig5_per_patient_alpha.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: fig5_per_patient_alpha.png")


# ── Analysis & Findings ──────────────────────────────────────────────
print("\n" + "=" * 70)
print("KEY FINDINGS")
print("=" * 70)

findings = []

# 1. Where is Q maximized relative to LRG optimal threshold?
print("\n1. Q_max scale vs LRG optimal threshold scale:")
q_max_at_opt = []
q_max_anywhere = []
for _, row in df.iterrows():
    q_max_anywhere.append(row["q_max"])

# Compare n_clusters at Q_max vs at optimal threshold
k_at_qmax = df["n_clusters_at_q_max"].values
k_at_opt = df["n_clusters_at_optimal"].values
ratio = k_at_qmax / k_at_opt
print(f"   Ratio k(Q_max) / k(optimal): mean={np.nanmean(ratio):.2f}, "
      f"std={np.nanstd(ratio):.2f}")
print(f"   k(Q_max) mean={np.nanmean(k_at_qmax):.1f}, "
      f"k(optimal) mean={np.nanmean(k_at_opt):.1f}")
finding1 = (f"The scale of maximum modularity (mean k={np.nanmean(k_at_qmax):.1f}) "
            f"differs from the LRG optimal threshold scale (mean k={np.nanmean(k_at_opt):.1f}). "
            f"Ratio k(Q_max)/k(optimal) = {np.nanmean(ratio):.2f} +/- {np.nanstd(ratio):.2f}.")
findings.append(("Q_max vs LRG optimal", finding1))
print(f"   -> {finding1}")

# 2. Does scale of max Q shift between rest and task?
print("\n2. Rest vs Task scale of max Q:")
for band in BANDS:
    rest_k = []
    task_k = []
    for patient in ALL_PATIENTS:
        for phase in REST_PHASES:
            data = results[patient][band][phase]
            if data is not None:
                rest_k.append(data["n_clusters"][np.nanargmax(data["q_values"])])
        for phase in TASK_PHASES:
            data = results[patient][band][phase]
            if data is not None:
                task_k.append(data["n_clusters"][np.nanargmax(data["q_values"])])
    if rest_k and task_k:
        print(f"   {band:>12s}: rest k={np.mean(rest_k):.1f}+/-{np.std(rest_k):.1f}, "
              f"task k={np.mean(task_k):.1f}+/-{np.std(task_k):.1f}")

rest_df = df[df["phase"].isin(REST_PHASES)]
task_df = df[df["phase"].isin(TASK_PHASES)]
rest_k_mean = rest_df["n_clusters_at_q_max"].mean()
task_k_mean = task_df["n_clusters_at_q_max"].mean()
finding2 = (f"Rest phases have max Q at k={rest_k_mean:.1f} clusters on average, "
            f"task phases at k={task_k_mean:.1f}. "
            f"{'Task shows finer-grained' if task_k_mean > rest_k_mean else 'Rest shows finer-grained'} "
            f"optimal modular structure.")
findings.append(("Rest vs Task scale shift", finding2))
print(f"   -> {finding2}")

# 3. Profile shape consistency
print("\n3. Profile shape consistency across patients:")
# Compute pairwise correlation of Q profiles between patients (for alpha band)
from itertools import combinations
for band in ["alpha", "beta"]:
    corrs = []
    for p1, p2 in combinations(ALL_PATIENTS, 2):
        for phase in PHASES:
            d1 = results[p1][band][phase]
            d2 = results[p2][band][phase]
            if d1 is not None and d2 is not None:
                # Interpolate to common grid
                k_min = max(d1["n_clusters"].min(), d2["n_clusters"].min())
                k_max = min(d1["n_clusters"].max(), d2["n_clusters"].max())
                if k_min < k_max:
                    k_grid = np.linspace(k_min, k_max, 30)
                    s1 = np.argsort(d1["n_clusters"])
                    s2 = np.argsort(d2["n_clusters"])
                    q1 = np.interp(k_grid, d1["n_clusters"][s1], d1["q_values"][s1])
                    q2 = np.interp(k_grid, d2["n_clusters"][s2], d2["q_values"][s2])
                    r = np.corrcoef(q1, q2)[0, 1]
                    corrs.append(r)
    if corrs:
        print(f"   {band:>12s}: mean inter-patient r={np.mean(corrs):.3f}, "
              f"min={np.min(corrs):.3f}, max={np.max(corrs):.3f}")

finding3 = "Profile shape consistency evaluated via pairwise correlation of Q(k) curves between patients."
findings.append(("Profile consistency", finding3))

# 4. Max Q differences between conditions
print("\n4. Max Q values by condition:")
for band in BANDS:
    rest_q = []
    task_q = []
    for patient in ALL_PATIENTS:
        for phase in REST_PHASES:
            data = results[patient][band][phase]
            if data is not None:
                rest_q.append(np.nanmax(data["q_values"]))
        for phase in TASK_PHASES:
            data = results[patient][band][phase]
            if data is not None:
                task_q.append(np.nanmax(data["q_values"]))
    if rest_q and task_q:
        print(f"   {band:>12s}: rest Q_max={np.mean(rest_q):.4f}+/-{np.std(rest_q):.4f}, "
              f"task Q_max={np.mean(task_q):.4f}+/-{np.std(task_q):.4f}")

rest_qmax = rest_df["q_max"].mean()
task_qmax = task_df["q_max"].mean()
finding4 = (f"Overall max Q: rest={rest_qmax:.4f}, task={task_qmax:.4f}. "
            f"{'Task' if task_qmax > rest_qmax else 'Rest'} condition shows higher peak modularity.")
findings.append(("Max Q rest vs task", finding4))
print(f"   -> {finding4}")

# 5. Band-specific patterns
print("\n5. Band-specific patterns:")
band_qmax = {}
band_kmax = {}
for band in BANDS:
    band_q = df[df["band"] == band]["q_max"].values
    band_k = df[df["band"] == band]["n_clusters_at_q_max"].values
    band_qmax[band] = np.nanmean(band_q)
    band_kmax[band] = np.nanmean(band_k)
    print(f"   {band:>12s}: mean Q_max={np.nanmean(band_q):.4f}, mean k={np.nanmean(band_k):.1f}")

best_band = max(band_qmax, key=band_qmax.get)
worst_band = min(band_qmax, key=band_qmax.get)
finding5 = (f"Highest mean Q_max in {best_band} ({band_qmax[best_band]:.4f}), "
            f"lowest in {worst_band} ({band_qmax[worst_band]:.4f}). "
            f"Band with most clusters at Q_max: "
            f"{max(band_kmax, key=band_kmax.get)} (k={max(band_kmax.values()):.1f}).")
findings.append(("Band-specific patterns", finding5))
print(f"   -> {finding5}")


# ── Write FINDINGS.md ────────────────────────────────────────────────
findings_md = OUT_DIR / "FINDINGS.md"
with open(findings_md, "w") as f:
    f.write("# Scale-Resolved Modularity Profiles: Key Findings\n\n")
    f.write("## Method\n\n")
    f.write("For each patient/band/phase combination, the LRG dendrogram was cut at "
            f"{N_THRESHOLDS} evenly-spaced distance thresholds. At each cut, "
            "modularity Q was computed using the original MSC adjacency matrix as "
            "the underlying graph and the dendrogram-derived partition as communities. "
            "This produces a Q(scale) profile revealing at which hierarchical level "
            "the network is most modular.\n\n")
    f.write(f"**Patients:** {', '.join(ALL_PATIENTS)} (all with 4 recording phases)\n\n")
    f.write(f"**Bands:** {', '.join(BANDS)}\n\n")
    f.write(f"**Phases:** {', '.join(PHASES)}\n\n")
    f.write(f"**FC method:** {FC_METHOD}\n\n")

    f.write("## Key Findings\n\n")
    for title, text in findings:
        f.write(f"### {title}\n\n{text}\n\n")

    f.write("## Figures\n\n")
    f.write("1. **fig1_pat02_modularity_profiles.png** - Q(n_clusters) profiles for Pat_02 across all bands/phases\n")
    f.write("2. **fig2_scale_max_q_by_phase.png** - Bar chart of n_clusters at max Q, averaged over patients\n")
    f.write("3. **fig3_max_q_heatmap.png** - Heatmap of max Q values (phases x bands)\n")
    f.write("4. **fig4_rest_vs_task_profiles.png** - Grand average Q profiles: rest vs task with SEM shading\n")
    f.write("5. **fig5_per_patient_alpha.png** - Per-patient alpha Q profiles with LRG optimal threshold markers\n")

print(f"\nSaved: {findings_md}")
print(f"\nAll outputs in: {OUT_DIR}")
print("Done!")
