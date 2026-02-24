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
import sys
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from scipy.signal import find_peaks
from scipy.optimize import linear_sum_assignment

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.visuals.lrg import compute_partition_stability_index

# ── Config ────────────────────────────────────────────────────────────
PATIENT = "Pat_02"
PHASE = "rsPre"
BAND = "beta"
FC_METHOD = "msc"
MSC_CACHE = ROOT / "data" / "msc_cache"
LRG_CACHE = ROOT / "data" / "lrg_cache"
DATASET_ROOT = ROOT / "data" / "stereoeeg_patients"
OUTPUT_DIR = ROOT / "data" / "figures" / "report_mslcd_section" / PATIENT
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


# ── Helper functions ──────────────────────────────────────────────────
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
tau_entropy = np.logspace(np.log10(tau_min * 0.3), np.log10(tau_max * 3), 500)
S_values = np.array([entropy_at_tau(eigenvalues, t) for t in tau_entropy])

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
print(f"  S range: [{S_values.min():.2f}, {S_values.max():.2f}] (ln N = {np.log(N):.2f})")
print(f"  Found {len(tau_stars)} characteristic scales:")
for k, ts in enumerate(tau_stars):
    print(f"    τ*_{k+1} = {ts:.4f} (C = {C_values[peaks_idx[k]]:.2f})")

# ── Compute multiscale community structure ────────────────────────────
print("Computing D(τ) → linkage → PSI → n* across τ grid...")
tau_coarse = np.logspace(np.log10(tau_min), np.log10(tau_max), 60)
nstar_tau = np.zeros(len(tau_coarse), dtype=int)
psi_collection = {}

for i, tau in enumerate(tau_coarse):
    Z, psi_vals, n_comms = linkage_and_psi(eigenvalues, eigenvectors, tau)
    if len(psi_vals) > 0:
        nstar_tau[i] = int(n_comms[np.argmax(psi_vals)])
    else:
        nstar_tau[i] = 1

# Store PSI at ~5 representative τ values for fig 2
psi_tau_indices = np.linspace(0, len(tau_coarse) - 1, 5, dtype=int)
for idx in psi_tau_indices:
    tau = tau_coarse[idx]
    Z, psi_vals, n_comms = linkage_and_psi(eigenvalues, eigenvectors, tau)
    psi_collection[tau] = (psi_vals, n_comms)

# ── Compute communities for alluvial and metastability ────────────────
# Use same tau/n_clusters pairs as the Sankey (via compute_clustering_across_tau)
# so that alluvial PDF and Sankey HTML show identical data.
print("Computing communities (shared Sankey + alluvial)...")

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
communities_aligned = align_cluster_labels(communities_raw)

# Compute metastability μ_i
n_transitions = len(alluvial_taus) - 1
mu = np.zeros(N)
if n_transitions > 0:
    for s in range(n_transitions):
        mu += (communities_aligned[s] != communities_aligned[s + 1]).astype(float)
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

for p, (n_cut, psi_val) in enumerate(zip(prominent_n_sorted, prominent_psi_sorted)):
    fig, ax = plt.subplots(figsize=(14, 7))

    # Cut height for n_cut clusters: midpoint between adjacent merge steps
    if n_cut < N:
        cut_height = (merge_heights[N - n_cut - 1] + merge_heights[N - n_cut]) / 2
    else:
        cut_height = merge_heights[0] * 0.5

    dendrogram(
        Z_dendro, ax=ax, no_labels=True,
        color_threshold=cut_height,
        above_threshold_color="black",
        leaf_rotation=0,
    )

    ax.axhline(cut_height, color="blue", ls="--", lw=2,
               label=rf"$n^* = {n_cut}$  ($\Psi = {psi_val:.3f}$)")

    ax.set_yscale("log")
    ax.set_ylim(tmin_ax, tmax_ax)
    ax.set_xlabel("Channel index (reordered by dendrogram)", fontsize=12)
    ax.set_ylabel(r"Merge height $\Delta_n$ (ultrametric distance)", fontsize=12)
    ax.set_title(
        rf"{PATIENT}, {PHASE}, {BRAIN_BAND_TEX_DICT[BAND]} — "
        rf"LRG dendrogram, cut at $n^* = {n_cut}$",
        fontsize=13, fontweight="bold",
    )
    ax.legend(fontsize=11, loc="upper right")
    ax.grid(axis="y", alpha=0.3)

    out = OUTPUT_DIR / f"fig_dendrogram_n{n_cut}.pdf"
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

for k, tau in enumerate(tau_vals_sorted):
    psi_vals, n_comms = psi_collection[tau]
    if len(psi_vals) == 0:
        continue
    color = cmap(k / max(len(tau_vals_sorted) - 1, 1))
    nstar = int(n_comms[np.argmax(psi_vals)])
    ax.plot(n_comms, psi_vals, "o-", color=color, markersize=3, lw=1.5,
            label=rf"$\tau = {tau:.3f}$ → $n^* = {nstar}$")
    ax.plot(nstar, psi_vals[np.argmax(psi_vals)], "*", color=color,
            markersize=12, zorder=5)

# Vertical red dashed lines at prominent dendrogram cuts
for n_cut, psi_val in zip(prominent_n_sorted, prominent_psi_sorted):
    ax.axvline(n_cut, color="red", ls="--", lw=1.5, alpha=0.7)
    ax.text(n_cut + 0.3, ax.get_ylim()[1] * 0.95 if ax.get_ylim()[1] > 0 else psi_val,
            rf"$n^*={n_cut}$", color="red", fontsize=9, va="top")

ax.set_xlabel(r"Number of communities $n$", fontsize=12)
ax.set_ylabel(r"Partition stability index $\Psi(n;\tau)$", fontsize=12)
ax.set_title(
    rf"{PATIENT}, {PHASE}, {BRAIN_BAND_TEX_DICT[BAND]} — "
    r"$\Psi(n;\tau)$ across diffusion scales",
    fontsize=13, fontweight="bold",
)
ax.legend(fontsize=9, loc="upper right")
ax.grid(alpha=0.3)
ax.set_xlim(1, None)

out = OUTPUT_DIR / "fig_psi_profile.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")

# ======================================================================
# FIG 3: Entropy S(τ) and susceptibility C(τ)
# ======================================================================
print("Fig 3: entropy_susceptibility...")

fig, (ax_S, ax_C) = plt.subplots(2, 1, figsize=(12, 8), sharex=True,
                                  gridspec_kw={"height_ratios": [1, 1], "hspace": 0.08})

# --- Top: S(τ) ---
ax_S.plot(log10_tau, S_values, color="0.3", lw=2)
ax_S.axhline(np.log(N), color="steelblue", ls=":", lw=1, alpha=0.7,
             label=rf"$\ln N = {np.log(N):.2f}$")
ax_S.axvline(np.log10(tau_min), color="C3", ls="--", lw=1, alpha=0.5)
ax_S.axvline(np.log10(tau_max), color="C0", ls="--", lw=1, alpha=0.5)
ax_S.axvspan(np.log10(tau_min), np.log10(tau_max), alpha=0.06, color="steelblue")
for ts in tau_stars:
    ax_S.axvline(np.log10(ts), color="crimson", ls=":", lw=1, alpha=0.5)
ax_S.set_ylabel(r"Von Neumann entropy $S(\tau)$", fontsize=12)
ax_S.set_title(
    rf"{PATIENT}, {PHASE}, {BRAIN_BAND_TEX_DICT[BAND]} — "
    r"Entropy and entropic susceptibility",
    fontsize=13, fontweight="bold",
)
ax_S.legend(fontsize=9)
ax_S.grid(alpha=0.3)
ax_S.set_ylim(bottom=0)

# --- Bottom: C(τ) ---
ax_C.plot(log10_tau, C_values, color="0.3", lw=2)
ax_C.axvline(np.log10(tau_min), color="C3", ls="--", lw=1, alpha=0.5,
             label=rf"$\tau' = 1/\lambda_{{\max}}$")
ax_C.axvline(np.log10(tau_max), color="C0", ls="--", lw=1, alpha=0.5,
             label=rf"$1/\lambda_2$")
ax_C.axvspan(np.log10(tau_min), np.log10(tau_max), alpha=0.06, color="steelblue",
             label="Resolution window")
for k, ts in enumerate(tau_stars):
    ax_C.axvline(np.log10(ts), color="crimson", ls=":", lw=1.5)
    ax_C.plot(np.log10(ts), C_values[peaks_idx[k]], "v", color="crimson",
              markersize=10, zorder=5,
              label=rf"$\tau^*_{k+1} = {ts:.3f}$")
ax_C.set_xlabel(r"$\log_{10}\,\tau$", fontsize=12)
ax_C.set_ylabel(r"Entropic susceptibility $C(\tau)$", fontsize=12)
ax_C.legend(fontsize=8, ncol=2, loc="upper right")
ax_C.grid(alpha=0.3)
ax_C.set_ylim(bottom=0)

out = OUTPUT_DIR / "fig_entropy_susceptibility.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")

# ======================================================================
# FIG 4: n*(τ) step function  (clean version)
# ======================================================================
print("Fig 4: nstar_tau...")

fig, ax = plt.subplots(figsize=(12, 5))

ax.step(np.log10(tau_coarse), nstar_tau, where="mid", color="0.3", lw=2.5,
        label=r"$n^*(\tau) = \arg\max_n\,\Psi(n;\tau)$")
ax.fill_between(np.log10(tau_coarse), 0, nstar_tau, step="mid",
                alpha=0.08, color="steelblue")

# Mark characteristic scales
for k, ts in enumerate(tau_stars):
    idx = np.argmin(np.abs(tau_coarse - ts))
    ns = nstar_tau[idx]
    ax.plot(np.log10(ts), ns, "v", color="crimson", markersize=12, zorder=5,
            label=rf"$\tau^*_{k+1} = {ts:.3f}$, $n^* = {ns}$" if k < 5 else None)

# τ window
ax.axvline(np.log10(tau_min), color="C3", ls="--", lw=1, alpha=0.5,
           label=rf"$\tau' = 1/\lambda_{{\max}}$")
ax.axvline(np.log10(tau_max), color="C0", ls="--", lw=1, alpha=0.5,
           label=rf"$1/\lambda_2$")
ax.axvspan(np.log10(tau_min), np.log10(tau_max), alpha=0.06, color="steelblue")

ax.set_xlabel(r"$\log_{10}\,\tau$", fontsize=12)
ax.set_ylabel(r"Optimal communities $n^*(\tau)$", fontsize=12)
ax.set_title(
    rf"{PATIENT}, {PHASE}, {BRAIN_BAND_TEX_DICT[BAND]} — "
    r"Optimal community count $n^*(\tau)$",
    fontsize=13, fontweight="bold",
)
ax.legend(fontsize=9, loc="upper right")
ax.grid(alpha=0.3)
ax.set_ylim(bottom=0)

out = OUTPUT_DIR / "fig_nstar_tau.pdf"
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
            rf"Metastability $\mu_i$ ({len(alluvial_taus)} scales, n={alluvial_nclust[0]}→{alluvial_nclust[-1]})"
        ),
        colorbar=True,
    )
    out = OUTPUT_DIR / "fig_metastable_nodes.pdf"
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
        rf"Metastability $\mu_i$ ({len(alluvial_taus)} scales, n={alluvial_nclust[0]}→{alluvial_nclust[-1]})",
        fontsize=13, fontweight="bold", y=1.02,
    )
    out = OUTPUT_DIR / "fig_metastable_nodes.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved (matplotlib): {out}")

# ======================================================================
# FIG 6: Alluvial diagram (community flow across scales)
# ======================================================================
print("Fig 6: community_alluvial...")

n_scales = len(alluvial_taus)
all_labels = sorted(set(l for arr in communities_aligned for l in arr))
cmap_comm = matplotlib.colormaps["tab10"]
comm_colors = {l: cmap_comm(i % 10) for i, l in enumerate(all_labels)}

# Sort nodes by trajectory for clean visual grouping
trajectories = [tuple(communities_aligned[s][i] for s in range(n_scales)) for i in range(N)]
sort_order = sorted(range(N), key=lambda i: trajectories[i])
labels_sorted = [communities_aligned[s][sort_order] for s in range(n_scales)]

fig, ax = plt.subplots(figsize=(max(4 * n_scales, 10), 8))
bar_width = 0.25
gap = N * 0.015

block_pos = []
for s in range(n_scales):
    labels = labels_sorted[s]
    seen, unique = set(), []
    for l in labels:
        if l not in seen:
            unique.append(l); seen.add(l)
    positions = {}
    y = 0
    for comm in unique:
        count = int(np.sum(labels == comm))
        positions[comm] = (y, y + count)
        y += count + gap
    block_pos.append(positions)

for s in range(n_scales):
    for comm, (y0, y1) in block_pos[s].items():
        count = int(y1 - y0)
        rect = plt.Rectangle(
            (s - bar_width / 2, y0), bar_width, count,
            facecolor=comm_colors[comm], edgecolor="white", lw=1.5, zorder=3,
        )
        ax.add_patch(rect)
        if count >= 5:
            ax.text(s, (y0 + y1) / 2, str(count),
                    ha="center", va="center", fontsize=8,
                    fontweight="bold", color="white", zorder=4)

for s in range(n_scales - 1):
    left, right = labels_sorted[s], labels_sorted[s + 1]
    trans = Counter(zip(left, right))
    off_left = {c: block_pos[s][c][0] for c in block_pos[s]}
    off_right = {c: block_pos[s + 1][c][0] for c in block_pos[s + 1]}
    for (c_l, c_r), count in sorted(trans.items(), key=lambda x: (-x[0][0], -x[1])):
        y_src_bot = off_left[c_l]; y_src_top = y_src_bot + count; off_left[c_l] = y_src_top
        y_tgt_bot = off_right[c_r]; y_tgt_top = y_tgt_bot + count; off_right[c_r] = y_tgt_top
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
    [rf"$\tau={t:.1f}$, $n={n}$" for t, n in TAU_NCLUST_PAIRS], fontsize=10,
)
ax.set_ylabel("Nodes (sorted by trajectory)", fontsize=11)
ax.set_xlim(-0.6, n_scales - 0.4)
ax.set_ylim(-gap, max(bp[1] for bp in block_pos[-1].values()) + gap)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
handles = [plt.Rectangle((0, 0), 1, 1, facecolor=comm_colors[l]) for l in all_labels]
ax.legend(handles, [f"Community {l}" for l in all_labels], fontsize=8,
          loc="upper right", ncol=2)
ax.set_title(
    rf"{PATIENT}, {PHASE}, {BRAIN_BAND_TEX_DICT[BAND]} — "
    r"Community evolution across $\tau$ (hierarchical merging)",
    fontsize=13, fontweight="bold",
)

out = OUTPUT_DIR / "fig_community_alluvial.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")

# ======================================================================
# FIG 7: Interactive Plotly Sankey (like fig5_sankey_*.html)
# ======================================================================
print("Fig 7: sankey_community (Plotly)...")

# Reuse partitions already computed for the alluvial (identical data)
tau_values_sankey = np.array(alluvial_taus)

sankey_fig = create_sankey_diagram(
    partitions, tau_values_sankey,
    title=f"{PATIENT} {PHASE} {BAND.upper()} — Cluster Evolution Across τ Scales",
    width=1400, height=700,
)

out = OUTPUT_DIR / f"fig_sankey_community_{BAND}.html"
sankey_fig.write_html(str(out))
print(f"  Saved: {out}")

# ======================================================================
# FIG 8: Interactive brain connectome (Plotly 3D + fsaverage brain surface)
# ======================================================================
print("Fig 8: brain_connectome (Plotly + brain mesh)...")

import plotly.graph_objects as go
from nilearn.datasets import fetch_surf_fsaverage
from nilearn.surface import load_surf_mesh

# ── Brain surface mesh (semi-transparent fsaverage pial) ──
fsaverage = fetch_surf_fsaverage()
brain_traces = []
for hemi_name, mesh_key in [("L", "pial_left"), ("R", "pial_right")]:
    mesh = load_surf_mesh(fsaverage[mesh_key])
    vtx, tri = mesh.coordinates, mesh.faces
    brain_traces.append(go.Mesh3d(
        x=vtx[:, 0], y=vtx[:, 1], z=vtx[:, 2],
        i=tri[:, 0], j=tri[:, 1], k=tri[:, 2],
        color="lightgray", opacity=0.08,
        hoverinfo="none", showlegend=False, name=f"brain_{hemi_name}",
    ))

# ── LRG cluster coloring (use finest prominent PSI cut) ──
n_clusters_brain = prominent_n_sorted[-1]  # finest meaningful n from PSI criterion
cluster_labels = fcluster(lrg.linkage_matrix, n_clusters_brain, criterion="maxclust")
print(f"  LRG clusters: {n_clusters_brain} (finest PSI peak)")

cmap_brain = matplotlib.colormaps["tab20"]
node_colors_hex = [matplotlib.colors.to_hex(cmap_brain(c % 20)) for c in cluster_labels]

# ── ALL edges with variable width + alpha (no threshold) ──
triu_idx = np.triu_indices(N, k=1)
edge_i_all, edge_j_all = triu_idx[0], triu_idx[1]
edge_w_all = A[triu_idx]
print(f"  Edges: {len(edge_w_all)} total (all displayed, alpha~0 for weakest)")

# Normalize to [0, 1]
w_min, w_max = edge_w_all.min(), edge_w_all.max()
w_norm = (edge_w_all - w_min) / (w_max - w_min + 1e-10)

# Width: power-law for contrast (ultra-thin → very thick)
EDGE_POWER = 2.5
WIDTH_MIN, WIDTH_MAX = 0.3, 20.0

# Color: light gray → black, alpha: ~0 (invisible) → 0.9 (opaque)
def _w_to_rgba(w_n):
    gray = int(200 * (1 - w_n))
    alpha = 0.02 + 0.88 * np.power(w_n, 1.5)  # ~0 at w_min, ~0.9 at w_max
    return f"rgba({gray},{gray},{gray},{alpha:.3f})"

# Bucket edges by normalized weight for efficient rendering (~25 traces)
N_BUCKETS = 25
bucket_bounds = np.linspace(0, 1.001, N_BUCKETS + 1)
edge_traces = []

for b in range(N_BUCKETS):
    bmask = (w_norm >= bucket_bounds[b]) & (w_norm < bucket_bounds[b + 1])
    if not np.any(bmask):
        continue
    bw = w_norm[bmask].mean()
    bwidth = WIDTH_MIN + np.power(bw, EDGE_POWER) * (WIDTH_MAX - WIDTH_MIN)
    bcolor = _w_to_rgba(bw)

    x, y, z = [], [], []
    for ei, ej in zip(edge_i_all[bmask], edge_j_all[bmask]):
        x.extend([coords_mni[ei, 0], coords_mni[ej, 0], None])
        y.extend([coords_mni[ei, 1], coords_mni[ej, 1], None])
        z.extend([coords_mni[ei, 2], coords_mni[ej, 2], None])

    edge_traces.append(go.Scatter3d(
        x=x, y=y, z=z, mode="lines",
        line=dict(width=bwidth, color=bcolor),
        hoverinfo="none", showlegend=False,
    ))

# ── Colorbar: invisible scatter trace mapping weight → color ──
colorbar_trace = go.Scatter3d(
    x=[None], y=[None], z=[None], mode="markers",
    marker=dict(
        size=0.1,
        colorscale=[[0, "rgba(200,200,200,0.02)"], [1, "rgba(0,0,0,0.9)"]],
        cmin=w_min, cmax=w_max,
        colorbar=dict(title="MSC weight", thickness=15, len=0.5),
        color=[w_min],
    ),
    hoverinfo="none", showlegend=False,
)

# ── Node markers with LRG cluster colors (larger, stable size) ──
node_trace = go.Scatter3d(
    x=coords_mni[:, 0], y=coords_mni[:, 1], z=coords_mni[:, 2],
    mode="markers",
    marker=dict(
        size=8, color=node_colors_hex,
        line=dict(width=1, color="gray"),
        sizemode="diameter",
    ),
    hoverinfo="text",
    text=[f"Ch {i} (cluster {cluster_labels[i]})" for i in range(N)],
    showlegend=False,
)

# ── Assemble figure: brain surface + edges + colorbar + nodes ──
fig_brain = go.Figure(data=brain_traces + edge_traces + [colorbar_trace, node_trace])
fig_brain.update_layout(
    title=dict(
        text=f"{PATIENT} {BAND.upper()} {PHASE} — Brain Connectome "
             f"({n_clusters_brain} LRG Clusters)",
        font=dict(size=16),
    ),
    scene=dict(
        xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False),
        aspectmode="data",
        bgcolor="white",
    ),
    width=1200, height=800,
    paper_bgcolor="white",
)

out = OUTPUT_DIR / f"fig_brain_connectome_{BAND}.html"
fig_brain.write_html(str(out))
print(f"  Saved: {out}")

print(f"\nDone! 8 figures saved to {OUTPUT_DIR}")
