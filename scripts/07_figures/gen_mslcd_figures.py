#!/usr/bin/env python3
"""Generate figures for Section: Multiscale community detection (sec:mslcd).

All fast — uses cached MSC matrix + eigendecomposition. ~10 seconds.

Figures (PDF, for LaTeX):
  1. fig_dendrogram_example.pdf      — Dendrogram with optimal cut highlighted
  2. fig_psi_profile.pdf             — Ψ(n;τ) for several τ values
  3. fig_entropy_susceptibility.pdf   — S(τ) and C(τ) with peaks marked
  4. fig_nstar_tau.pdf                — n*(τ) step function
  5. fig_metastable_nodes.pdf         — Brain-space μ_i visualization
  6. fig_community_alluvial.pdf       — Alluvial diagram across τ*

Interactive (HTML, Plotly/nilearn):
  7. fig_sankey_community.html        — Plotly Sankey of cluster evolution
  8. fig_brain_connectome.html        — nilearn interactive 3D brain connectome

Run: python scripts/gen_mslcd_figures.py
"""
from pathlib import Path
from collections import Counter

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import networkx as nx
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, leaves_list
from scipy.spatial.distance import squareform
from scipy.signal import find_peaks
from scipy.optimize import linear_sum_assignment

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.visuals.lrg import compute_partition_stability_index, _load_channel_labels

# ── Config ────────────────────────────────────────────────────────────
PATIENT = "Pat_02"
PHASE = "rest_pre"
BAND = "beta"
FC_METHOD = "msc"
from lrg_eegfc.config.paths import MSC_CACHE, LRG_CACHE, SEEG_DATAPATH, FIGURES_ROOT
DATASET_ROOT = SEEG_DATAPATH
OUTPUT_DIR = FIGURES_ROOT / "report_mslcd_section" / PATIENT
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Shared tau/n_clusters pairs for Sankey + alluvial (fine → coarse)
TAU_NCLUST_PAIRS = [
    (0.1, 15),
    (0.3, 12),
    (0.5, 8),
    (0.8, 6),
    (1.0, 5),
    (2.0, 3),
    (5.0, 2),
]

# ── Load MSC and compute Laplacian eigendecomposition ─────────────────
print(f"Loading MSC matrix: {PATIENT}, {PHASE}, {BAND}...")
A = load_msc_matrix(PATIENT, PHASE, BAND,
                    cache_root=MSC_CACHE, sparsify="none", n_surrogates=0, nperseg=4096)
np.fill_diagonal(A, 0)
N = A.shape[0]
print(f"  N = {N} channels")

D_diag = A.sum(axis=1)
L = np.diag(D_diag) - A
eigenvalues, eigenvectors = np.linalg.eigh(L)
eigenvalues = np.maximum(eigenvalues, 0.0)

lambda_gap = eigenvalues[1]
lambda_max = eigenvalues[-1]
tau_min = 1.0 / lambda_max
tau_max = 1.0 / lambda_gap
print(f"  λ_2 = {lambda_gap:.4f}, λ_max = {lambda_max:.2f}")
print(f"  τ window: [{tau_min:.4f}, {tau_max:.4f}]")

# ── Load cached LRG result (for dendrogram + brain connectome) ────────
print("Loading cached LRG result...")
lrg = load_lrg_result(PATIENT, PHASE, BAND, FC_METHOD, LRG_CACHE)
if lrg is None:
    raise RuntimeError(f"No cached LRG result for {PATIENT}/{PHASE}/{BAND}/{FC_METHOD}")
print(f"  LRG: {lrg.n_nodes} nodes, optimal_threshold = {lrg.optimal_threshold:.4f}")

# ── Build NetworkX graph (shared by Sankey + alluvial) ────────────────
G = nx.from_numpy_array(A)
Gcc_nodes = max(nx.connected_components(G), key=len)
Gcc = G.subgraph(Gcc_nodes).copy()
print(f"  Giant component: {Gcc.number_of_nodes()} nodes")

# ── Load channel labels ──────────────────────────────────────────────
print("Loading channel labels...")
channel_labels = _load_channel_labels(PATIENT, DATASET_ROOT)
assert len(channel_labels) == N, f"Label count {len(channel_labels)} != N={N}"
print(f"  Loaded {len(channel_labels)} labels: {channel_labels[:5]}...")


# ── Helper functions ──────────────────────────────────────────────────
# Curated maximally-distinct colors (Kelly / Boynton, reordered for contrast)
# Curated maximally-distinct colors — all dark/saturated, no pastels
_DISTINCT_COLORS = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#d63b00", "#f032e6", "#469990", "#9A6324", "#800000",
    "#808000", "#000075", "#e03080", "#1a7820", "#7030c0",
    "#b05010", "#206060", "#c04080", "#305090", "#704020",
]


def get_distinct_colors(n):
    """Generate n maximally perceptually distinct colors."""
    if n <= len(_DISTINCT_COLORS):
        return [matplotlib.colors.to_rgba(c) for c in _DISTINCT_COLORS[:n]]
    # Fall back to golden-angle HSV for very large n
    colors = []
    for i in range(n):
        hue = (i * 0.618033988749895) % 1.0
        colors.append(matplotlib.colors.hsv_to_rgb([hue, 0.75, 0.85]))
    return colors


def get_subtree_cluster(Z, clusters, n_nodes):
    """For each node (leaves + internal), determine its fcluster community.

    Returns dict: node_id → cluster_label (or -1 if mixed / above cut).
    """
    node_cluster = {}
    for i in range(n_nodes):
        node_cluster[i] = clusters[i]
    for i in range(len(Z)):
        left, right = int(Z[i, 0]), int(Z[i, 1])
        cl = node_cluster.get(left, -1)
        cr = node_cluster.get(right, -1)
        node_cluster[n_nodes + i] = cl if (cl == cr and cl != -1) else -1
    return node_cluster


def compute_propagator(eigenvalues, eigenvectors, tau):
    """K(τ) = exp(-τ L) via eigendecomposition."""
    exp_vals = np.exp(-tau * eigenvalues)
    return (eigenvectors * exp_vals[None, :]) @ eigenvectors.T


def compute_ultrametric_distance(K):
    """D_ij(τ) = (1 - δ_ij) / K_ij(τ)."""
    D = np.zeros_like(K)
    mask = ~np.eye(K.shape[0], dtype=bool)
    D[mask] = 1.0 / np.where(K[mask] > 1e-30, K[mask], 1e-30)
    return D


def entropy_at_tau(eigenvalues, tau):
    """S(τ) = -Σ p_ℓ ln p_ℓ, numerically stable."""
    log_boltz = -tau * eigenvalues
    log_Z = np.logaddexp.reduce(log_boltz)
    log_p = log_boltz - log_Z
    p = np.exp(log_p)
    return -np.sum(p * np.where(p > 1e-300, log_p, 0))


def linkage_and_psi(eigenvalues, eigenvectors, tau):
    """Compute D(τ) → linkage → PSI at a single τ."""
    K = compute_propagator(eigenvalues, eigenvectors, tau)
    D = compute_ultrametric_distance(K)
    D_condensed = squareform(D, checks=False)
    Z = linkage(D_condensed, method="average")
    psi_vals, n_comms = compute_partition_stability_index(Z)
    return Z, psi_vals, n_comms


def sensible_psi_peaks(psi_vals, n_comms, log_drop_max=0.5, min_n_spacing=5):
    """Find PSI peaks using mean-floor + log-drop rule.

    1. Find ALL local maxima of PSI(n).
    2. Discard peaks below mean(PSI) — noise floor from the curve shape.
    3. Sort survivors by PSI descending.
    4. Keep the top peak, then keep each subsequent peak if:
       - log10(PSI_prev) - log10(PSI_k) <= log_drop_max  (0.5 decades)
       - |n_k - n_kept| >= min_n_spacing for all already-kept peaks
    5. Stop at first peak that fails the log-drop.

    Returns (kept_indices, count).
    """
    if len(psi_vals) == 0:
        return [], 0
    peak_idx, _ = find_peaks(psi_vals)
    if len(peak_idx) == 0:
        return [np.argmax(psi_vals)], 1
    # Noise floor: only peaks above mean PSI are candidates
    psi_mean = psi_vals.mean()
    peak_idx = peak_idx[psi_vals[peak_idx] > psi_mean]
    if len(peak_idx) == 0:
        return [np.argmax(psi_vals)], 1
    # Sort by PSI descending
    order = np.argsort(psi_vals[peak_idx])[::-1]
    peak_idx = peak_idx[order]
    kept = [peak_idx[0]]
    kept_n = {int(n_comms[peak_idx[0]])}
    for ci in peak_idx[1:]:
        n_val = int(n_comms[ci])
        if any(abs(n_val - kn) < min_n_spacing for kn in kept_n):
            continue
        log_drop = (np.log10(psi_vals[kept[-1]])
                    - np.log10(max(psi_vals[ci], 1e-30)))
        if log_drop <= log_drop_max:
            kept.append(ci)
            kept_n.add(n_val)
        else:
            break
    return kept, len(kept)


def align_cluster_labels(labels_list):
    """Align cluster labels across scales by maximum overlap (Hungarian)."""
    aligned = [labels_list[0].copy()]
    for s in range(1, len(labels_list)):
        prev, curr = aligned[-1], labels_list[s]
        u_prev, u_curr = np.unique(prev), np.unique(curr)
        overlap = np.zeros((len(u_prev), len(u_curr)))
        for i, lp in enumerate(u_prev):
            for j, lc in enumerate(u_curr):
                overlap[i, j] = np.sum((prev == lp) & (curr == lc))
        dim = max(len(u_prev), len(u_curr))
        cost = np.zeros((dim, dim))
        cost[:len(u_prev), :len(u_curr)] = -overlap
        row_ind, col_ind = linear_sum_assignment(cost)
        label_map = {}
        for r, c in zip(row_ind, col_ind):
            if r < len(u_prev) and c < len(u_curr):
                label_map[u_curr[c]] = u_prev[r]
        used = set(label_map.values())
        next_lbl = max(max(used) + 1, max(u_prev) + 1) if used else 1
        for lc in u_curr:
            if lc not in label_map:
                label_map[lc] = next_lbl
                next_lbl += 1
        aligned.append(np.array([label_map[l] for l in curr]))
    return aligned


# ── Compute entropy S(τ) and susceptibility C(τ) ─────────────────────
print("Computing entropy curve...")
# Extend tau far enough that S → 0 (1 − S/ln N → 1)
tau_entropy = np.logspace(np.log10(tau_min * 0.3), np.log10(tau_max * 50), 700)
S_values = np.array([entropy_at_tau(eigenvalues, t) for t in tau_entropy])
lnN = np.log(N)

log10_tau = np.log10(tau_entropy)
C_values = -np.gradient(S_values, log10_tau)

# Find peaks of C(τ) within resolution window
in_window = (tau_entropy >= tau_min) & (tau_entropy <= tau_max)
C_masked = C_values.copy()
C_masked[~in_window] = 0
peaks_idx, peak_props = find_peaks(C_masked, prominence=0.05 * C_masked.max())
if len(peaks_idx) > 0:
    proms = peak_props["prominences"]
    peaks_idx = peaks_idx[proms > 0.1 * proms.max()]

tau_stars = tau_entropy[peaks_idx]
print(f"  S range: [{S_values.min():.4f}, {S_values.max():.2f}] (ln N = {lnN:.2f})")
print(f"  1−S/lnN range: [{(1 - S_values.max()/lnN):.4f}, {(1 - S_values.min()/lnN):.4f}]")
print(f"  Found {len(tau_stars)} characteristic scales:")
for k, ts in enumerate(tau_stars):
    print(f"    τ*_{k+1} = {ts:.4f} (C = {C_values[peaks_idx[k]]:.2f})")

# ── Compute multiscale community structure ────────────────────────────
print("Computing D(τ) → linkage → PSI → n* across fine τ grid...")
tau_fine = np.logspace(np.log10(tau_min), np.log10(tau_max), 150)
nstar_tau = np.zeros(len(tau_fine), dtype=int)
n_sensible_tau = np.zeros(len(tau_fine), dtype=int)

for i, tau in enumerate(tau_fine):
    Z, psi_vals, n_comms = linkage_and_psi(eigenvalues, eigenvectors, tau)
    if len(psi_vals) > 0:
        nstar_tau[i] = int(n_comms[np.argmax(psi_vals)])
        _, n_sensible_tau[i] = sensible_psi_peaks(psi_vals, n_comms)
    else:
        nstar_tau[i] = 1
        n_sensible_tau[i] = 1

# PSI at τ_min = 1/λ_max specifically — used for peak selection + red curve
print("Computing PSI at τ_min = 1/λ_max for peak selection...")
Z_tmin, psi_tmin, n_comms_tmin = linkage_and_psi(eigenvalues, eigenvectors, tau_min)
psi_peaks_tmin_idx, n_sensible_tmin = sensible_psi_peaks(psi_tmin, n_comms_tmin)
psi_peak_n_tmin = [int(n_comms_tmin[i]) for i in psi_peaks_tmin_idx]
psi_peak_psi_tmin = [psi_tmin[i] for i in psi_peaks_tmin_idx]
print(f"  PSI peaks at τ_min ({n_sensible_tmin} sensible): n* = {psi_peak_n_tmin}")

# Store PSI at ~6 representative τ values for fig 2 background curves
psi_collection = {}
psi_tau_indices = np.linspace(0, len(tau_fine) - 1, 6, dtype=int)
for idx in psi_tau_indices:
    tau = tau_fine[idx]
    if abs(tau - tau_min) / tau_min > 0.05:
        Z, psi_vals, n_comms = linkage_and_psi(eigenvalues, eigenvectors, tau)
        psi_collection[tau] = (psi_vals, n_comms)

# ── Compute communities for alluvial (fixed n, for visualization) ──────
print("Computing communities for alluvial/Sankey (fixed n)...")

from lrg_eegfc.visuals.metastable import (
    compute_clustering_across_tau,
    create_sankey_diagram,
)

alluvial_taus = [t for t, _ in TAU_NCLUST_PAIRS]
alluvial_nclust = [n for _, n in TAU_NCLUST_PAIRS]

partitions, n_clusters_dict = compute_clustering_across_tau(
    None, np.array(alluvial_taus), Gcc,
    n_clusters_list=alluvial_nclust,
)

# Extract ordered label arrays (one per τ scale)
communities_raw = [partitions[t] for t in alluvial_taus]
communities_aligned_alluvial = align_cluster_labels(communities_raw)

# ── Data-driven metastability: n_max(τ) = max S(τ), capped ───────────
# At each τ, S(τ) is the set of sensible PSI peaks. n_max(τ) = max S(τ)
# is the finest partition the network still supports at that diffusion scale.
# Cap at N/2: peaks above this are trivial single-node splits, not
# genuine mesoscopic structure (ultrametric distance becomes near-uniform
# at moderate τ, generating a forest of high-n PSI artifacts).
N_TAU_META = 20
NMAX_CAP = N // 2
tau_meta = np.logspace(np.log10(tau_min), np.log10(tau_max), N_TAU_META)
nmax_meta = np.zeros(N_TAU_META, dtype=int)

print(f"Computing data-driven metastability ({N_TAU_META} τ in "
      f"[{tau_min:.4f}, {tau_max:.4f}], n_max cap={NMAX_CAP})...")

meta_communities = []
for i, tau in enumerate(tau_meta):
    Z_m, psi_m, n_comms_m = linkage_and_psi(eigenvalues, eigenvectors, tau)
    peaks_idx_m, n_sensible_m = sensible_psi_peaks(psi_m, n_comms_m)
    if n_sensible_m > 0:
        sensible_n = [int(n_comms_m[pi]) for pi in peaks_idx_m]
        # Cap: discard trivial high-n peaks (single-node splitting noise)
        sensible_n = [n for n in sensible_n if n <= NMAX_CAP]
        nmax_meta[i] = max(sensible_n) if sensible_n else 1
    else:
        nmax_meta[i] = 1
    labels_m = fcluster(Z_m, max(nmax_meta[i], 1), criterion="maxclust")
    meta_communities.append(labels_m)

print(f"  n_max(τ) trajectory: {nmax_meta.tolist()}")

meta_communities_aligned = align_cluster_labels(meta_communities)

n_transitions = N_TAU_META - 1
mu = np.zeros(N)
if n_transitions > 0:
    for s in range(n_transitions):
        mu += (meta_communities_aligned[s] != meta_communities_aligned[s + 1]).astype(float)
    mu /= n_transitions

n_metastable = np.sum(mu > 0)
print(f"  Metastable nodes (μ > 0): {n_metastable}/{N}")
print(f"  Max μ = {mu.max():.3f}, Mean μ = {mu.mean():.3f}")
print("  All objects computed.\n")

# Dendrogram: use cached LRG linkage (consistent with plot_lrg_full_panel)
Z_dendro = lrg.linkage_matrix
merge_heights = Z_dendro[:, 2]
tmin_ax = merge_heights[0] * 0.8
tmax_ax = merge_heights[-1] * 1.05

# ── Find prominent PSI peaks for multi-cut dendrogram ────────────────
# PSI from the LRG dendrogram linkage
psi_vals_dendro, n_comms_dendro = compute_partition_stability_index(Z_dendro)

# Candidates: global argmax (boundary) + local maxima from find_peaks
global_max_idx = np.argmax(psi_vals_dendro)
local_peaks_idx, _ = find_peaks(psi_vals_dendro, prominence=0.001)
candidate_idx = np.unique(np.concatenate([[global_max_idx], local_peaks_idx]))

# Sort candidates by PSI descending
candidate_idx = candidate_idx[np.argsort(psi_vals_dendro[candidate_idx])[::-1]]

# Half-rule in LOG scale: drop between consecutive kept peaks <= 0.5 decades
# i.e. PSI_i >= PSI_prev * 10^(-0.5) ≈ PSI_prev / 3.16
# + minimum n-spacing of 5 to avoid redundant nearby peaks
LOG_DROP_MAX = 0.5
MIN_N_SPACING = 5
MAX_PANELS = 4

prominent_idx = [candidate_idx[0]]
prominent_n_set = {int(n_comms_dendro[candidate_idx[0]])}
for ci in candidate_idx[1:]:
    if len(prominent_idx) >= MAX_PANELS:
        break
    n_val = int(n_comms_dendro[ci])
    # Skip if too close in n to any already-kept peak
    if any(abs(n_val - kn) < MIN_N_SPACING for kn in prominent_n_set):
        continue
    # Check log-scale drop from previous kept peak
    log_drop = np.log10(psi_vals_dendro[prominent_idx[-1]]) - np.log10(psi_vals_dendro[ci])
    if log_drop <= LOG_DROP_MAX:
        prominent_idx.append(ci)
        prominent_n_set.add(n_val)
    else:
        break

prominent_n = [int(n_comms_dendro[i]) for i in prominent_idx]
prominent_psi = [psi_vals_dendro[i] for i in prominent_idx]
print(f"  Prominent PSI peaks ({len(prominent_n)}):")
for n_val, psi_val in zip(prominent_n, prominent_psi):
    print(f"    n* = {n_val}, PSI = {psi_val:.4f}")

# ======================================================================
# FIG 1: Separate dendrogram per prominent PSI peak
# ======================================================================
print("Fig 1: dendrogram_example (one PDF per cut)...")

# Sort prominent peaks by n ascending for consistent ordering
prom_order = np.argsort(prominent_n)
prominent_n_sorted = [prominent_n[i] for i in prom_order]
prominent_psi_sorted = [prominent_psi[i] for i in prom_order]

# Pre-compute leaf order (same ordering dendrogram will use)
leaf_order = leaves_list(Z_dendro)

# Pre-compute network layout for inset (shared across all dendrogram cuts)
print("  Computing network layout for inset...")
G_net = nx.from_numpy_array(A)
pos_spectral = nx.spectral_layout(G_net)
pos_net = nx.spring_layout(G_net, pos=pos_spectral, k=3.0 / np.sqrt(N),
                           iterations=100, seed=42)

# Edge data: only show top 10% strongest edges for visual clarity
triu_i, triu_j = np.triu_indices(N, k=1)
all_weights = A[triu_i, triu_j]
weight_threshold = np.percentile(all_weights, 90)
strong_mask = all_weights >= weight_threshold
inset_edges = list(zip(triu_i[strong_mask].tolist(), triu_j[strong_mask].tolist()))
inset_weights = all_weights[strong_mask]

# Normalize edge widths: thinnest → 0.15, thickest → 4.0
w_min_e, w_max_e = inset_weights.min(), inset_weights.max()
inset_widths = 0.15 + 3.85 * (inset_weights - w_min_e) / (w_max_e - w_min_e + 1e-10)

# Per-edge RGBA (gray with alpha ∝ weight)
inset_edge_rgba = []
for w in inset_weights:
    frac = (w - w_min_e) / (w_max_e - w_min_e + 1e-10)
    inset_edge_rgba.append((0.4, 0.4, 0.4, 0.12 + 0.58 * frac))

# Node strengths for sizing
node_strengths = A.sum(axis=1)
s_min, s_max = node_strengths.min(), node_strengths.max()
inset_node_sizes = 12 + 140 * (node_strengths - s_min) / (s_max - s_min + 1e-10)

print(f"  Inset: {len(inset_edges)} edges (top 10%), layout ready.")

for p, (n_cut, psi_val) in enumerate(zip(prominent_n_sorted, prominent_psi_sorted)):
    fig_w = max(N * 0.18, 16)
    fig, ax = plt.subplots(figsize=(fig_w, 7))

    # Cut height for n_cut clusters: midpoint between adjacent merge steps
    if n_cut < N:
        cut_height = (merge_heights[N - n_cut - 1] + merge_heights[N - n_cut]) / 2
    else:
        cut_height = merge_heights[0] * 0.5

    # Get cluster assignments and build color mapping
    clusters = fcluster(Z_dendro, n_cut, criterion="maxclust")
    palette = get_distinct_colors(n_cut)

    # Determine left-to-right cluster appearance order
    first_seen = {}
    for leaf in leaf_order:
        c = clusters[leaf]
        if c not in first_seen:
            first_seen[c] = len(first_seen)
    cluster_lr_order = sorted(first_seen.keys(), key=lambda c: first_seen[c])
    cluster_color_map = {c: palette[i] for i, c in enumerate(cluster_lr_order)}

    # Robust per-link coloring via link_color_func (handles singletons correctly)
    node_cluster = get_subtree_cluster(Z_dendro, clusters, N)

    def _link_color(cluster_id):
        c = node_cluster.get(cluster_id, -1)
        if c != -1:
            return matplotlib.colors.to_hex(cluster_color_map[c])
        return "0.5"

    dn = dendrogram(
        Z_dendro, ax=ax, no_labels=True,
        link_color_func=_link_color,
    )

    # Cut line
    ax.axhline(cut_height, color="0.3", ls="--", lw=2)

    ax.set_yscale("log")
    ax.set_ylim(tmin_ax, tmax_ax)

    # Add staggered channel name labels colored by community
    for i, leaf_idx in enumerate(leaf_order):
        x = 5 + i * 10  # default dendrogram leaf spacing
        y_frac = -0.02 if i % 2 == 0 else -0.08
        label = channel_labels[leaf_idx]
        color = cluster_color_map[clusters[leaf_idx]]
        ax.text(x, y_frac, label, transform=ax.get_xaxis_transform(),
                rotation=90, ha="center", va="top", fontsize=11,
                fontweight="bold", color=color, clip_on=False)

    ax.set_xticks([])
    ax.set_ylabel(r"$\mathcal{D}\,/\,\mathcal{D}_{\max}$", fontsize=18)
    ax.tick_params(axis="y", labelsize=13)
    ax.grid(axis="y", alpha=0.3)

    # ── Network inset (top-right corner) ─────────────────────────────
    ax_ins = ax.inset_axes([0.70, 0.45, 0.28, 0.52])

    # Edges
    nx.draw_networkx_edges(G_net, pos_net, edgelist=inset_edges,
                           width=inset_widths.tolist(),
                           edge_color=inset_edge_rgba,
                           ax=ax_ins)

    # Nodes (color by this cut's clustering, size by strength)
    node_colors = [cluster_color_map[clusters[i]] for i in range(N)]
    nx.draw_networkx_nodes(G_net, pos_net, node_color=node_colors,
                           node_size=inset_node_sizes,
                           edgecolors="0.3", linewidths=0.4,
                           ax=ax_ins)

    ax_ins.set_axis_off()

    # Legend (upper left to avoid inset)
    ax.legend(
        handles=[Line2D([0], [0], color="0.3", ls="--", lw=2.5,
                        label=rf"$n^* = {n_cut}$  ($\Psi = {psi_val:.3f}$)")],
        fontsize=15, loc="upper left",
    )

    out = OUTPUT_DIR / f"fig_dendrogram_n{n_cut}_{PHASE}_{BAND}.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved: {out}")

# ======================================================================
# FIG 2: PSI profile Ψ(n;τ) at several τ
# ======================================================================
print("Fig 2: psi_profile...")

fig, ax = plt.subplots(figsize=(12, 6))
cmap = plt.cm.viridis
tau_vals_sorted = sorted(psi_collection.keys())

# Background curves at representative τ values
for k, tau in enumerate(tau_vals_sorted):
    psi_vals, n_comms = psi_collection[tau]
    if len(psi_vals) == 0:
        continue
    color = cmap(k / max(len(tau_vals_sorted) - 1, 1))
    nstar = int(n_comms[np.argmax(psi_vals)])
    ax.plot(n_comms, psi_vals, "o-", color=color, markersize=3, lw=1.2, alpha=0.7,
            label=rf"$\tau = {tau:.3f}$ → $n^* = {nstar}$")
    ax.plot(nstar, psi_vals[np.argmax(psi_vals)], "*", color=color,
            markersize=10, zorder=4)

# Highlighted RED curve: PSI at τ_min = 1/λ_max
nstar_tmin = int(n_comms_tmin[np.argmax(psi_tmin)])
ax.plot(n_comms_tmin, psi_tmin, "o-", color="crimson", markersize=4, lw=2.5, zorder=5,
        label=rf"$\tau_\min = 1/\lambda_{{\max}} = {tau_min:.4f}$ → $n^* = {nstar_tmin}$")
ax.plot(nstar_tmin, psi_tmin[np.argmax(psi_tmin)], "*", color="crimson",
        markersize=16, zorder=6)

# Vertical dashed lines at PSI peaks from the τ_min curve
for n_pk, psi_pk in zip(psi_peak_n_tmin, psi_peak_psi_tmin):
    ax.axvline(n_pk, color="crimson", ls="--", lw=1.5, alpha=0.6)
    ax.text(n_pk, 0.97, rf"$n^*={n_pk}$", color="crimson", fontsize=11,
            fontweight="bold", rotation=90, ha="right", va="top",
            transform=ax.get_xaxis_transform())

ax.set_xlabel(r"Number of communities $n$", fontsize=16)
ax.set_ylabel(r"Partition stability index $\Psi(n;\tau)$", fontsize=16)
ax.tick_params(labelsize=13)
ax.legend(fontsize=12, loc="upper right")
ax.grid(alpha=0.3)
ax.set_xlim(1, None)

out = OUTPUT_DIR / f"fig_psi_profile_{PHASE}_{BAND}.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")

# ======================================================================
# FIG 3: Entropy (1 − S̃) and susceptibility C̃(τ) — single twinx plot
# ======================================================================
print("Fig 3: entropy_susceptibility...")

COLOR_ENT = "#4363d8"   # blue for entropy
COLOR_C   = "#e6194b"   # red for susceptibility

def _make_entropy_susceptibility_fig(band_name, eig_vals, out_path):
    """Generate entropy/susceptibility figure for a given band."""
    N_b = len(eig_vals)
    lnN_b = np.log(N_b)
    tau_min_b = 1.0 / eig_vals[-1]
    tau_max_b = 1.0 / eig_vals[1]

    tau_grid = np.logspace(np.log10(tau_min_b * 0.3), np.log10(tau_max_b * 50), 700)
    S_vals = np.array([entropy_at_tau(eig_vals, t) for t in tau_grid])
    log10_t = np.log10(tau_grid)

    # Normalized: S̃ = S/ln N,  C̃ = -(1/ln N) dS/d(log10 τ)
    S_norm = S_vals / lnN_b
    C_norm = -np.gradient(S_norm, log10_t)

    # Peaks of C̃ within resolution window
    in_win = (tau_grid >= tau_min_b) & (tau_grid <= tau_max_b)
    C_masked = C_norm.copy()
    C_masked[~in_win] = 0
    pk_idx, pk_props = find_peaks(C_masked, prominence=0.05 * C_masked.max())
    if len(pk_idx) > 0:
        proms = pk_props["prominences"]
        pk_idx = pk_idx[proms > 0.1 * proms.max()]
    tau_star = tau_grid[pk_idx]

    fig, ax_ent = plt.subplots(figsize=(12, 6))
    ax_C = ax_ent.twinx()

    one_minus_S = 1.0 - S_norm

    # Left axis: 1 − S̃(τ)
    ax_ent.plot(tau_grid, one_minus_S, color=COLOR_ENT, lw=2.5,
                label=r"$1 - \tilde{S}(\tau)$")
    ax_ent.set_ylabel(r"$1 - \tilde{S}(\tau)$", fontsize=16, color=COLOR_ENT)
    ax_ent.tick_params(axis="y", labelcolor=COLOR_ENT, labelsize=13)
    ax_ent.set_ylim(-0.05, 1.05)

    # Right axis: C̃(τ) = -(1/ln N) dS/d(log10 τ)
    ax_C.plot(tau_grid, C_norm, color=COLOR_C, lw=2.5,
              label=r"$\tilde{C}(\tau) = -d\tilde{S}/d\log_{10}\tau$")
    for k_pk, ts in enumerate(tau_star):
        ax_C.plot(ts, C_norm[pk_idx[k_pk]], "v", color=COLOR_C,
                  markersize=12, zorder=5,
                  label=rf"$\tau^*_{k_pk+1} = {ts:.3f}$" if k_pk < 4 else None)
    ax_C.set_ylabel(r"$\tilde{C}(\tau)$", fontsize=16, color=COLOR_C)
    ax_C.tick_params(axis="y", labelcolor=COLOR_C, labelsize=13)
    ax_C.set_ylim(bottom=0)

    # Shared x-axis
    ax_ent.set_xscale("log")
    ax_ent.set_xlabel(r"Diffusion time $\tau$", fontsize=16)
    ax_ent.tick_params(axis="x", labelsize=13)

    # Combined legend
    h1, l1 = ax_ent.get_legend_handles_labels()
    h2, l2 = ax_C.get_legend_handles_labels()
    ax_ent.legend(h1 + h2, l1 + l2, fontsize=12, loc="center right")
    ax_ent.grid(alpha=0.3)

    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")

# Generate for beta (already have eigenvalues)
_make_entropy_susceptibility_fig(
    BAND, eigenvalues,
    OUTPUT_DIR / f"fig_entropy_susceptibility_{PHASE}_{BAND}.pdf",
)

# Generate for high_gamma
print("  Computing high_gamma eigenvalues...")
A_hg = load_msc_matrix(PATIENT, PHASE, "high_gamma",
                        cache_root=MSC_CACHE, sparsify="none",
                        n_surrogates=0, nperseg=4096)
np.fill_diagonal(A_hg, 0)
D_hg = A_hg.sum(axis=1)
L_hg = np.diag(D_hg) - A_hg
eig_hg, _ = np.linalg.eigh(L_hg)
eig_hg = np.maximum(eig_hg, 0.0)
_make_entropy_susceptibility_fig(
    "high_gamma", eig_hg,
    OUTPUT_DIR / f"fig_entropy_susceptibility_{PHASE}_high_gamma.pdf",
)

# ======================================================================
# FIG 4: Number of sensible partitions vs τ
# ======================================================================
print("Fig 4: sensible partitions...")

fig, ax = plt.subplots(figsize=(12, 5))

ax.plot(tau_fine, n_sensible_tau, color="0.3", lw=2.5, zorder=3)
ax.fill_between(tau_fine, 0, n_sensible_tau, alpha=0.10, color="steelblue")

# Mark characteristic C(τ) peaks
for k, ts in enumerate(tau_stars):
    idx = np.argmin(np.abs(tau_fine - ts))
    ns = n_sensible_tau[idx]
    ax.plot(ts, ns, "v", color="crimson", markersize=12, zorder=5,
            label=rf"$\tau^* = {ts:.3f}$, sensible peaks$={ns}$" if k < 5 else None)

# τ window
ax.axvline(tau_min, color="0.4", ls="--", lw=1, alpha=0.6,
           label=rf"$1/\lambda_{{\max}}$")
ax.axvline(tau_max, color="0.4", ls="--", lw=1, alpha=0.6,
           label=rf"$1/\lambda_2$")
ax.axvspan(tau_min, tau_max, alpha=0.06, color="steelblue",
           label="Resolution window")

ax.set_xscale("log")
ax.set_xlabel(r"Diffusion time $\tau$", fontsize=16)
ax.set_ylabel(r"Sensible partitions $|\{n : \Psi_n \geq \Psi_{\mathrm{prev}}/2\}|$",
              fontsize=15)
ax.tick_params(labelsize=13)
ax.legend(fontsize=12, loc="upper right")
ax.grid(alpha=0.3, which="both")
ax.set_ylim(bottom=0)

out = OUTPUT_DIR / f"fig_sensible_partitions_{PHASE}_{BAND}.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")

# ======================================================================
# FIG 5: Metastable nodes on brain
# ======================================================================
print("Fig 5: metastable_nodes...")

from lrg_eegfc.visuals.spatial import load_spatial_metadata, prepare_spatial_coordinates

metadata = load_spatial_metadata(PATIENT, DATASET_ROOT)
coords_mni = prepare_spatial_coordinates(metadata, scale="mm", to_mni=True)
assert coords_mni.shape[0] == N, f"Coordinate count {coords_mni.shape[0]} != N={N}"

try:
    from nilearn import plotting as ni_plot

    display = ni_plot.plot_markers(
        node_values=mu,
        node_coords=coords_mni,
        node_cmap="RdYlGn_r",
        node_size=mu * 80 + 15,
        display_mode="lzry",
        title=(
            rf"{PATIENT}, {PHASE}, {BRAIN_BAND_TEX_DICT[BAND]} — "
            rf"Metastability $\mu_i$ "
            rf"({N_TAU_META} $\tau$ in $[1/\lambda_{{\max}},\,1/\lambda_2]$, "
            rf"$n_{{\max}}$: {nmax_meta[0]}$\to${nmax_meta[-1]})"
        ),
        colorbar=True,
    )
    out = OUTPUT_DIR / f"fig_metastable_nodes_{PHASE}_{BAND}.pdf"
    display.savefig(str(out), dpi=150)
    display.close()
    print(f"  Saved (nilearn): {out}")

except Exception as e:
    print(f"  nilearn plot_markers failed ({e}), using matplotlib fallback")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    projections = [
        ("Sagittal (L)", 1, 2, "y (mm)", "z (mm)"),
        ("Coronal", 0, 2, "x (mm)", "z (mm)"),
        ("Axial", 0, 1, "x (mm)", "y (mm)"),
    ]
    for ax_p, (title, xi, yi, xlabel, ylabel) in zip(axes, projections):
        sc = ax_p.scatter(coords_mni[:, xi], coords_mni[:, yi], c=mu, cmap="RdYlGn_r",
                          s=40, edgecolors="0.3", lw=0.5, vmin=0, vmax=max(mu.max(), 0.01))
        ax_p.set_title(title, fontsize=11, fontweight="bold")
        ax_p.set_xlabel(xlabel); ax_p.set_ylabel(ylabel)
        ax_p.set_aspect("equal"); ax_p.grid(alpha=0.2)
    fig.colorbar(sc, ax=axes, shrink=0.7, label=r"Metastability $\mu_i$")
    fig.suptitle(
        rf"{PATIENT}, {PHASE}, {BRAIN_BAND_TEX_DICT[BAND]} — "
        rf"Metastability $\mu_i$ "
        rf"({N_TAU_META} $\tau$ in $[1/\lambda_{{\max}},\,1/\lambda_2]$, "
        rf"$n_{{\max}}$: {nmax_meta[0]}$\to${nmax_meta[-1]})",
        fontsize=13, fontweight="bold", y=1.02,
    )
    out = OUTPUT_DIR / f"fig_metastable_nodes_{PHASE}_{BAND}.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved (matplotlib): {out}")

# ======================================================================
# FIG 6: Alluvial diagram with pin labels packed side-by-side
# ======================================================================
print("Fig 6: community_alluvial (with pin labels)...")

n_scales = len(alluvial_taus)
all_labels = sorted(set(l for arr in communities_aligned_alluvial for l in arr))
comm_colors = {l: get_distinct_colors(len(all_labels))[i]
               for i, l in enumerate(all_labels)}

# Sort nodes by trajectory for clean visual grouping
trajectories = [tuple(communities_aligned_alluvial[s][i] for s in range(n_scales)) for i in range(N)]
sort_order = sorted(range(N), key=lambda i: trajectories[i])
labels_sorted = [communities_aligned_alluvial[s][sort_order] for s in range(n_scales)]

# ── Pre-compute text rows for each block ──────────────────────────────
# Pack pin labels side-by-side; wrap only when the line exceeds the box.
# Estimate ~5 chars fit per unit of bar_width using a "chars_per_line"
# heuristic calibrated to the figure. Each text row occupies ROW_HEIGHT
# y-units; block height = max(n_members_in_block, n_text_rows * ROW_HEIGHT).

PIN_FONTSIZE = 11
ROW_HEIGHT_MIN = 1.8   # minimum y-units per text row (prevents overlap)
CHARS_PER_LINE = 22    # fewer chars → horizontal padding inside box
BOX_PAD_Y = 1.0        # vertical padding (y-units) top+bottom inside box
bar_width = 0.75
gap = N * 0.015


def _pack_pins_into_rows(pin_names, chars_per_line=CHARS_PER_LINE):
    """Pack pin names side-by-side, wrapping when the line is full."""
    rows = []
    current = ""
    for pin in pin_names:
        candidate = f"{current}  {pin}" if current else pin
        if len(candidate) > chars_per_line and current:
            rows.append(current)
            current = pin
        else:
            current = candidate
    if current:
        rows.append(current)
    return rows


# First pass: compute block heights (may be taller than n_members)
block_info = []  # per scale: list of (comm, member_pins, text_rows, height)
for s in range(n_scales):
    labels = labels_sorted[s]
    seen, unique = set(), []
    for l in labels:
        if l not in seen:
            unique.append(l); seen.add(l)
    info = []
    for comm in unique:
        member_mask = labels == comm
        member_sorted_idx = np.where(member_mask)[0]
        pins = [channel_labels[sort_order[si]] for si in member_sorted_idx]
        rows = _pack_pins_into_rows(pins)
        height = max(len(pins), len(rows) * ROW_HEIGHT_MIN + 2 * BOX_PAD_Y)
        info.append((comm, pins, rows, height))
    block_info.append(info)

# Second pass: build block positions using computed heights
block_pos = []
for s in range(n_scales):
    positions = {}
    y = 0
    for comm, pins, rows, height in block_info[s]:
        positions[comm] = (y, y + height, rows)
        y += height + gap
    block_pos.append(positions)

fig, ax = plt.subplots(figsize=(max(4.5 * n_scales, 16), 14))

# Draw blocks + text rows
for s in range(n_scales):
    for comm, (y0, y1, rows) in block_pos[s].items():
        height = y1 - y0
        rect = plt.Rectangle(
            (s - bar_width / 2, y0), bar_width, height,
            facecolor=comm_colors[comm], edgecolor="white", lw=1.5, zorder=3,
        )
        ax.add_patch(rect)

        # Distribute text rows evenly inside the block (with padding)
        n_rows = len(rows)
        usable_h = height - 2 * BOX_PAD_Y
        row_h = usable_h / max(n_rows, 1)
        y_start = y0 + BOX_PAD_Y + row_h / 2
        for r, row_text in enumerate(rows):
            ax.text(s, y_start + r * row_h, row_text,
                    ha="center", va="center",
                    fontsize=PIN_FONTSIZE, fontweight="bold",
                    color="white", zorder=4, family="monospace")

# Bezier flows — use proportional offsets based on member counts
for s in range(n_scales - 1):
    left, right = labels_sorted[s], labels_sorted[s + 1]
    trans = Counter(zip(left, right))
    # Track y-offset within each block for stacking flows
    off_left = {}
    for comm, (y0, y1, _) in block_pos[s].items():
        n_mem = int(np.sum(labels_sorted[s] == comm))
        off_left[comm] = (y0, (y1 - y0) / n_mem)  # (current_y, scale)
    off_right = {}
    for comm, (y0, y1, _) in block_pos[s + 1].items():
        n_mem = int(np.sum(labels_sorted[s + 1] == comm))
        off_right[comm] = (y0, (y1 - y0) / n_mem)

    for (c_l, c_r), count in sorted(trans.items(), key=lambda x: (-x[0][0], -x[1])):
        cy_l, sc_l = off_left[c_l]
        y_src_bot = cy_l
        y_src_top = cy_l + count * sc_l
        off_left[c_l] = (y_src_top, sc_l)

        cy_r, sc_r = off_right[c_r]
        y_tgt_bot = cy_r
        y_tgt_top = cy_r + count * sc_r
        off_right[c_r] = (y_tgt_top, sc_r)

        x_left, x_right = s + bar_width / 2, (s + 1) - bar_width / 2
        t = np.linspace(0, 1, 80)
        x_curve = x_left + t * (x_right - x_left)
        h = 3 * t**2 - 2 * t**3
        ax.fill_between(x_curve,
                        y_src_bot + h * (y_tgt_bot - y_src_bot),
                        y_src_top + h * (y_tgt_top - y_src_top),
                        color=comm_colors[c_l], alpha=0.3, edgecolor="none", zorder=1)

ax.set_xticks(range(n_scales))
ax.set_xticklabels(
    [rf"$\tau={t:.1f}$" + "\n" + rf"$n={n}$" for t, n in TAU_NCLUST_PAIRS],
    fontsize=14,
)
y_top = max(bp[1] for bps in block_pos for _, bp in bps.items() if len(bp) >= 2)
ax.set_xlim(-0.8, n_scales - 0.2)
ax.set_ylim(-gap, y_top + gap * 2)
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

out = OUTPUT_DIR / f"fig_community_alluvial_{PHASE}_{BAND}.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")

# ======================================================================
# FIG 7: Interactive Plotly Sankey HTML (standard, no pin labels)
# ======================================================================
print("Fig 7: sankey_community (Plotly HTML)...")

import plotly.graph_objects as go
from plotly.express import colors as px_colors

sankey_labels = []
sankey_hovers = []
sankey_colors = []
color_palette = (px_colors.qualitative.Set3 + px_colors.qualitative.Dark2
                 + px_colors.qualitative.Pastel)
tau_to_offset = {}
node_count = 0

for i, tau in enumerate(alluvial_taus):
    tau_to_offset[tau] = node_count
    clusters_arr = partitions[tau]
    unique_c = np.unique(clusters_arr)
    for c in unique_c:
        members = sorted([channel_labels[j] for j in range(N) if clusters_arr[j] == c])
        n_mem = len(members)
        sankey_labels.append(f"C{c} ({n_mem})")
        hover_pins = "<br>".join(members) if n_mem <= 20 else (
            "<br>".join(members[:20]) + f"<br>... +{n_mem - 20} more")
        sankey_hovers.append(f"τ={tau:.2f} C{c}: {n_mem} pins<br>{hover_pins}")
        sankey_colors.append(color_palette[(c - 1) % len(color_palette)])
    node_count += len(unique_c)

sources, targets, values, link_colors = [], [], [], []
for i in range(len(alluvial_taus) - 1):
    tau1, tau2 = alluvial_taus[i], alluvial_taus[i + 1]
    c1_arr, c2_arr = partitions[tau1], partitions[tau2]
    trans_sk = Counter(zip(c1_arr, c2_arr))
    uc1 = list(np.unique(c1_arr))
    uc2 = list(np.unique(c2_arr))
    for (cl, cr), cnt in trans_sk.items():
        sources.append(tau_to_offset[tau1] + uc1.index(cl))
        targets.append(tau_to_offset[tau2] + uc2.index(cr))
        values.append(cnt)
        tc = sankey_colors[tau_to_offset[tau2] + uc2.index(cr)]
        if tc.startswith("rgb"):
            link_colors.append(tc.replace("rgb", "rgba").replace(")", ",0.4)"))
        elif tc.startswith("#"):
            hx = tc.lstrip("#")
            r, g, b = (int(hx[j:j+2], 16) for j in (0, 2, 4))
            link_colors.append(f"rgba({r},{g},{b},0.4)")
        else:
            link_colors.append(f"rgba(180,180,180,0.4)")

sankey_fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=20, thickness=30,
        line=dict(color="black", width=0.5),
        label=sankey_labels, color=sankey_colors,
        hovertemplate="%{customdata}<extra></extra>",
        customdata=sankey_hovers,
    ),
    link=dict(source=sources, target=targets, value=values, color=link_colors),
)])
sankey_fig.update_layout(
    title_text=f"{PATIENT} {PHASE} {BAND.upper()} — Cluster Evolution",
    font_size=10, width=1400, height=700,
)

out_html = OUTPUT_DIR / f"fig_sankey_community_{PHASE}_{BAND}.html"
sankey_fig.write_html(str(out_html))
print(f"  Saved (plotly): {out_html}")

out_pdf_sk = OUTPUT_DIR / f"fig_sankey_community_{PHASE}_{BAND}.pdf"
try:
    sankey_fig.write_image(str(out_pdf_sk), format="pdf", width=1400, height=700, scale=2)
    print(f"  Saved: {out_pdf_sk}")
except Exception as e:
    print(f"  PDF export failed ({e}), trying png fallback...")
    try:
        out_png_sk = out_pdf_sk.with_suffix(".png")
        sankey_fig.write_image(str(out_png_sk), format="png", width=1400, height=700, scale=2)
        print(f"  Saved: {out_png_sk}")
    except Exception as e2:
        print(f"  Image export failed ({e2}). Install kaleido: pip install kaleido")

# ======================================================================
# FIG 8: 3D brain connectomes — one per PSI partition (n* values)
# ======================================================================
print("Fig 8: brain_connectome (one per PSI n*)...")

import re
import plotly.graph_objects as go
from scipy.spatial import ConvexHull

# ── Brain surface (fsaverage pial) ───────────────────────────────────
from nilearn.datasets import fetch_surf_fsaverage
from nilearn.surface import load_surf_mesh

fsaverage = fetch_surf_fsaverage()
brain_meshes = []
for hemi_name, mesh_key in [("L", "pial_left"), ("R", "pial_right")]:
    mesh = load_surf_mesh(fsaverage[mesh_key])
    brain_meshes.append((hemi_name, mesh.coordinates, mesh.faces))

# ── Electrode shafts (shared across all n*) ──────────────────────────
electrode_groups = {}
for i, lbl in enumerate(channel_labels):
    m = re.match(r"([A-Za-z]+)", lbl)
    if m:
        electrode_groups.setdefault(m.group(1), []).append(i)

def _make_shaft_traces():
    traces = []
    for prefix, indices in electrode_groups.items():
        if len(indices) < 2:
            continue
        indices_sorted = sorted(indices,
            key=lambda idx: int(re.search(r"(\d+)", channel_labels[idx]).group(1)))
        x = [coords_mni[i, 0] for i in indices_sorted]
        y = [coords_mni[i, 1] for i in indices_sorted]
        z = [coords_mni[i, 2] for i in indices_sorted]
        traces.append(go.Scatter3d(
            x=x, y=y, z=z, mode="lines",
            line=dict(width=3.5, color="rgba(60,60,60,0.35)"),
            hoverinfo="none", showlegend=False,
        ))
    return traces

# ── Edges: top 20% strongest, curved arcs ────────────────────────────
triu_idx = np.triu_indices(N, k=1)
edge_i_all, edge_j_all = triu_idx[0], triu_idx[1]
edge_w_all = A[triu_idx]

EDGE_PERCENTILE = 85  # show top 15% — fewer but more visible
edge_threshold = np.percentile(edge_w_all, EDGE_PERCENTILE)
strong_mask = edge_w_all >= edge_threshold
edge_i = edge_i_all[strong_mask]
edge_j = edge_j_all[strong_mask]
edge_w = edge_w_all[strong_mask]
print(f"  Edges: {len(edge_w)}/{len(edge_w_all)} shown (top {100 - EDGE_PERCENTILE}%)")

w_min_e, w_max_e = edge_w.min(), edge_w.max()
w_norm = (edge_w - w_min_e) / (w_max_e - w_min_e + 1e-10)

# Centroid of all nodes (for arc offset direction)
centroid = coords_mni.mean(axis=0)

# Arc interpolation: push midpoint outward from centroid
ARC_NPTS = 14   # points per arc
ARC_BULGE = 6.0  # mm offset at midpoint — larger bulge for visibility

def _make_arc(pi, pj):
    """Create a curved arc from pi to pj, bulging outward from centroid."""
    mid = (pi + pj) / 2
    # Direction: midpoint → away from centroid
    outward = mid - centroid
    norm = np.linalg.norm(outward)
    if norm < 1e-6:
        outward = np.array([0, 0, 1.0])
    else:
        outward = outward / norm
    # Also add perpendicular component to the edge direction
    edge_dir = pj - pi
    edge_len = np.linalg.norm(edge_dir)
    if edge_len > 1e-6:
        perp = np.cross(edge_dir / edge_len, outward)
        perp_norm = np.linalg.norm(perp)
        if perp_norm > 1e-6:
            perp = perp / perp_norm
            outward = 0.7 * outward + 0.3 * perp
            outward = outward / np.linalg.norm(outward)
    # Quadratic Bezier: P(t) = (1-t)^2 * pi + 2t(1-t) * control + t^2 * pj
    control = mid + ARC_BULGE * outward
    t = np.linspace(0, 1, ARC_NPTS)
    pts = (np.outer((1 - t)**2, pi)
           + np.outer(2 * t * (1 - t), control)
           + np.outer(t**2, pj))
    return pts

EDGE_POWER = 1.8
WIDTH_MIN, WIDTH_MAX = 0.5, 16.0

def _edge_rgba(w_n):
    gray = int(80 * (1 - w_n))  # darker: 80→0
    alpha = 0.15 + 0.75 * np.power(w_n, 1.0)  # 0.15 → 0.90
    return f"rgba({gray},{gray},{gray},{alpha:.3f})"

N_BUCKETS = 15
bucket_bounds = np.linspace(0, 1.001, N_BUCKETS + 1)

def _make_edge_traces():
    traces = []
    for b in range(N_BUCKETS):
        bmask = (w_norm >= bucket_bounds[b]) & (w_norm < bucket_bounds[b + 1])
        if not np.any(bmask):
            continue
        bw = w_norm[bmask].mean()
        bwidth = WIDTH_MIN + np.power(bw, EDGE_POWER) * (WIDTH_MAX - WIDTH_MIN)
        bcolor = _edge_rgba(bw)
        x, y, z = [], [], []
        for ei, ej in zip(edge_i[bmask], edge_j[bmask]):
            arc = _make_arc(coords_mni[ei], coords_mni[ej])
            x.extend(arc[:, 0].tolist() + [None])
            y.extend(arc[:, 1].tolist() + [None])
            z.extend(arc[:, 2].tolist() + [None])
        traces.append(go.Scatter3d(
            x=x, y=y, z=z, mode="lines",
            line=dict(width=bwidth, color=bcolor),
            hoverinfo="none", showlegend=False,
        ))
    return traces

edge_traces_shared = _make_edge_traces()
shaft_traces_shared = _make_shaft_traces()
print(f"  {len(shaft_traces_shared)} electrode shafts, "
      f"{len(edge_traces_shared)} edge buckets")

# ── Unit sphere mesh (reused per node) ───────────────────────────────
NODE_RADIUS = 0.9  # mm — lives in data coordinates, scales with zoom

def _unit_sphere(n_lat=8, n_lon=12):
    """Create unit sphere vertices and triangle faces."""
    verts = [[0, 0, 1]]  # top pole
    for i in range(1, n_lat):
        lat = np.pi * i / n_lat
        for j in range(n_lon):
            lon = 2 * np.pi * j / n_lon
            verts.append([np.sin(lat)*np.cos(lon),
                          np.sin(lat)*np.sin(lon),
                          np.cos(lat)])
    verts.append([0, 0, -1])  # bottom pole
    verts = np.array(verts)

    faces = []
    # Top cap
    for j in range(n_lon):
        faces.append([0, 1 + j, 1 + (j + 1) % n_lon])
    # Middle strips
    for i in range(n_lat - 2):
        ring = 1 + i * n_lon
        nxt = ring + n_lon
        for j in range(n_lon):
            j1 = (j + 1) % n_lon
            faces.append([ring + j, nxt + j, ring + j1])
            faces.append([ring + j1, nxt + j, nxt + j1])
    # Bottom cap
    bottom = len(verts) - 1
    last_ring = 1 + (n_lat - 2) * n_lon
    for j in range(n_lon):
        faces.append([bottom, last_ring + (j + 1) % n_lon, last_ring + j])
    return verts, np.array(faces)

_sph_v, _sph_f = _unit_sphere()
_n_sph_v = len(_sph_v)

# ── Generate one figure per PSI n* ───────────────────────────────────
# Use n* values from PSI profile at τ_min
brain_n_values = sorted(set(psi_peak_n_tmin))
print(f"  Generating {len(brain_n_values)} brain figures: n* = {brain_n_values}")

for n_cut in brain_n_values:
    print(f"  n*={n_cut}...", end=" ")

    # Cluster assignment
    cl = fcluster(lrg.linkage_matrix, n_cut, criterion="maxclust")
    palette = get_distinct_colors(n_cut)
    unique_cl = sorted(np.unique(cl))
    cl_cmap = {c: palette[i % len(palette)] for i, c in enumerate(unique_cl)}
    node_hex = [matplotlib.colors.to_hex(cl_cmap[c]) for c in cl]

    # ── Brain surface traces ──
    btrs = []
    for hemi_name, vtx, tri in brain_meshes:
        btrs.append(go.Mesh3d(
            x=vtx[:, 0], y=vtx[:, 1], z=vtx[:, 2],
            i=tri[:, 0], j=tri[:, 1], k=tri[:, 2],
            color="lightgray", opacity=0.18,
            hoverinfo="none", showlegend=False,
            lighting=dict(ambient=0.9, diffuse=0.1),
        ))

    # ── Cluster surfaces: per-electrode tubes + cross-electrode bridges ──
    # Two layers:
    #   1. Per-electrode inflated hulls (tubes) — higher opacity, show local grouping
    #   2. Whole-cluster hull for multi-electrode clusters — lower opacity, show
    #      the spatial extent connecting far-apart electrode segments
    HULL_INFLATE = 2.0  # mm radius for per-electrode tubes
    BRIDGE_INFLATE = 3.0  # mm radius for cross-electrode bridges (bigger to connect)
    # Dense sphere of offsets for smooth, curvy hulls (Fibonacci lattice)
    def _sphere_offsets(radius, n_pts=26):
        pts = []
        golden = (1 + np.sqrt(5)) / 2
        for i in range(n_pts):
            theta = np.arccos(1 - 2 * (i + 0.5) / n_pts)
            phi = 2 * np.pi * i / golden
            pts.append([radius * np.sin(theta) * np.cos(phi),
                        radius * np.sin(theta) * np.sin(phi),
                        radius * np.cos(theta)])
        return np.array(pts)
    _offsets_tube = _sphere_offsets(HULL_INFLATE, 26)
    _offsets_bridge = _sphere_offsets(BRIDGE_INFLATE, 26)

    # Map each node → electrode prefix
    node_electrode = {}
    for i, lbl in enumerate(channel_labels):
        m = re.match(r"([A-Za-z]+)", lbl)
        node_electrode[i] = m.group(1) if m else "?"

    hull_trs = []
    for c in unique_cl:
        c_indices = np.where(cl == c)[0]
        if len(c_indices) < 2:
            continue
        # Group by electrode
        elec_groups = {}
        for idx in c_indices:
            elec_groups.setdefault(node_electrode[idx], []).append(idx)

        hex_c = matplotlib.colors.to_hex(cl_cmap[c])

        # Layer 1: per-electrode tubes (solid, higher opacity)
        for elec, indices in elec_groups.items():
            pts = coords_mni[indices]
            if pts.shape[0] < 2:
                continue
            inflated = [pts]
            for off in _offsets_tube:
                inflated.append(pts + off[None, :])
            inflated = np.vstack(inflated)
            try:
                hull = ConvexHull(inflated)
            except Exception:
                continue
            hull_trs.append(go.Mesh3d(
                x=inflated[hull.vertices, 0],
                y=inflated[hull.vertices, 1],
                z=inflated[hull.vertices, 2],
                alphahull=0,
                color=hex_c, opacity=0.30,
                hoverinfo="none", showlegend=False,
                lighting=dict(ambient=0.7, diffuse=0.3, specular=0.15,
                              fresnel=0.1),
                flatshading=False,
            ))

        # Layer 2: cross-electrode bridge (only if cluster spans 2+ electrodes)
        if len(elec_groups) >= 2:
            all_pts = coords_mni[c_indices]
            inflated = [all_pts]
            for off in _offsets_bridge:
                inflated.append(all_pts + off[None, :])
            inflated = np.vstack(inflated)
            try:
                hull = ConvexHull(inflated)
            except Exception:
                pass
            else:
                hull_trs.append(go.Mesh3d(
                    x=inflated[hull.vertices, 0],
                    y=inflated[hull.vertices, 1],
                    z=inflated[hull.vertices, 2],
                    alphahull=0,
                    color=hex_c, opacity=0.14,
                    hoverinfo="none", showlegend=False,
                    lighting=dict(ambient=0.8, diffuse=0.2, specular=0.05,
                                  fresnel=0.1),
                    flatshading=False,
                ))

    # ── Node spheres (Mesh3d — scale with zoom) ──
    node_sphere_trs = []
    for c in unique_cl:
        c_mask = cl == c
        c_idxs = np.where(c_mask)[0]
        all_v, all_fi, all_fj, all_fk = [], [], [], []
        for k_node, idx in enumerate(c_idxs):
            sv = _sph_v * NODE_RADIUS + coords_mni[idx]
            offset = k_node * _n_sph_v
            all_v.append(sv)
            all_fi.append(_sph_f[:, 0] + offset)
            all_fj.append(_sph_f[:, 1] + offset)
            all_fk.append(_sph_f[:, 2] + offset)
        all_v = np.vstack(all_v)
        all_fi = np.concatenate(all_fi)
        all_fj = np.concatenate(all_fj)
        all_fk = np.concatenate(all_fk)
        hex_c = matplotlib.colors.to_hex(cl_cmap[c])
        node_sphere_trs.append(go.Mesh3d(
            x=all_v[:, 0], y=all_v[:, 1], z=all_v[:, 2],
            i=all_fi, j=all_fj, k=all_fk,
            color=hex_c, opacity=1.0,
            hoverinfo="none", showlegend=False,
            lighting=dict(ambient=0.6, diffuse=0.4, specular=0.3,
                          fresnel=0.2),
        ))

    # ── Node labels + hover (text-only Scatter3d) ──
    # Offset label positions in data coords so they scale with zoom
    LABEL_OFFSET_Z = NODE_RADIUS + 1.5  # mm above sphere top
    label_coords = coords_mni.copy()
    label_coords[:, 2] += LABEL_OFFSET_Z
    node_label_tr = go.Scatter3d(
        x=label_coords[:, 0], y=label_coords[:, 1], z=label_coords[:, 2],
        mode="text",
        text=[channel_labels[i] for i in range(N)],
        textposition="middle center",
        textfont=dict(size=11, color="rgba(30,30,30,0.85)"),
        hoverinfo="text",
        hovertext=[f"{channel_labels[i]}<br>C{cl[i]}" for i in range(N)],
        showlegend=False,
    )

    # ── Colorbar ──
    cb_tr = go.Scatter3d(
        x=[None], y=[None], z=[None], mode="markers",
        marker=dict(
            size=0.1,
            colorscale=[[0, "rgba(100,100,100,0.08)"], [1, "rgba(0,0,0,0.85)"]],
            cmin=w_min_e, cmax=w_max_e,
            colorbar=dict(title="MSC", thickness=20, len=0.6, x=0.98,
                          titlefont=dict(size=16),
                          tickfont=dict(size=13)),
            color=[w_min_e],
        ),
        hoverinfo="none", showlegend=False,
    )

    # ── Assemble ──
    fig_b = go.Figure(
        data=(btrs + hull_trs + edge_traces_shared + shaft_traces_shared
              + [cb_tr] + node_sphere_trs + [node_label_tr])
    )
    fig_b.update_layout(
        scene=dict(
            xaxis=dict(visible=False), yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            aspectmode="data", bgcolor="white",
            camera=dict(
                eye=dict(x=-1.8, y=0.0, z=0.3),  # left lateral
                up=dict(x=0, y=0, z=1),
            ),
        ),
        showlegend=False,
        width=1400, height=900,
        paper_bgcolor="white",
        margin=dict(l=0, r=0, t=0, b=0),
    )

    out_html = OUTPUT_DIR / f"fig_brain_connectome_n{n_cut}_{PHASE}_{BAND}.html"
    fig_b.write_html(str(out_html))
    print(f"html saved")

print(f"\nDone! Figures saved to {OUTPUT_DIR}")
