#!/usr/bin/env python3
"""Generate figures for Section 3: Laplacian diffusion states and communication ultrametrics.

Two bands shown per figure (high_gamma = strong modular, theta = weak/diffuse)
to demonstrate the framework across network regimes.

Figures:
  1. fig_lrg_overview.pdf          — Schematic pipeline (band-independent)
  2. fig_spectrum_tau_window.pdf    — Laplacian spectrum + tau window (2 rows)
  3. fig_rho_tau.pdf                — Density matrix rho(tau) at 3 tau values (2 rows)

Run: python scripts/gen_ultrametric_figures.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.workflow.msc import load_msc_matrix

# ── Config ────────────────────────────────────────────────────────────
PATIENT = "Pat_02"
PHASE = "rsPre"
BANDS = ["high_gamma", "theta"]
MSC_CACHE = ROOT / "data" / "msc_cache"
NPERSEG = 4096
OUTPUT_DIR = ROOT / "data" / "figures" / "report_ultrametric_section" / PATIENT
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CMAP_RHO = "inferno"
CMAP_DIST = "cividis"

BAND_ROW_COLORS = {"high_gamma": "#64B5CD", "theta": "#55A868"}


# ── Helper functions ──────────────────────────────────────────────────
def compute_propagator(eigenvalues, eigenvectors, tau):
    """K(tau) = exp(-tau L) via eigendecomposition."""
    exp_vals = np.exp(-tau * eigenvalues)
    return (eigenvectors * exp_vals[None, :]) @ eigenvectors.T


def compute_rho(K):
    """rho(tau) = K(tau) / Tr[K(tau)]."""
    return K / np.trace(K)


def compute_ultrametric_distance(K):
    """D_ij(tau) = (1 - delta_ij) / K_ij(tau)."""
    D = np.zeros_like(K)
    mask = ~np.eye(K.shape[0], dtype=bool)
    D[mask] = 1.0 / np.where(K[mask] > 1e-30, K[mask], 1e-30)
    return D


# ── Load MSC and compute Laplacian objects for each band ──────────────
band_data = {}
for band in BANDS:
    print(f"Loading MSC: {PATIENT}, {PHASE}, {band}...")
    A = load_msc_matrix(PATIENT, PHASE, band,
                        cache_root=MSC_CACHE, sparsify="none",
                        n_surrogates=0, nperseg=NPERSEG)
    np.fill_diagonal(A, 0)
    N = A.shape[0]

    # Combinatorial Laplacian
    D_diag = A.sum(axis=1)
    L = np.diag(D_diag) - A

    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(L)
    eigenvalues = np.maximum(eigenvalues, 0.0)

    lambda_gap = eigenvalues[1]
    lambda_max = eigenvalues[-1]
    tau_min = 1.0 / lambda_max
    tau_max = 1.0 / lambda_gap
    tau_mid = np.sqrt(tau_min * tau_max)

    print(f"  N={N}, lambda_2={lambda_gap:.4f}, lambda_max={lambda_max:.2f}")
    print(f"  tau'={tau_min:.4f}, tau_mid={tau_mid:.4f}, 1/lambda_2={tau_max:.4f}")

    tau_values = [tau_min, tau_mid, tau_max]
    tau_labels = [
        rf"$\tau' = 1/\lambda_{{\max}}$ = {tau_min:.3f}",
        rf"$\tau_{{\mathrm{{mid}}}}$ = {tau_mid:.3f}",
        rf"$1/\lambda_2$ = {tau_max:.2f}",
    ]

    Ks = [compute_propagator(eigenvalues, eigenvectors, t) for t in tau_values]
    rhos = [compute_rho(K) for K in Ks]
    Ds = [compute_ultrametric_distance(K) for K in Ks]

    # Dendrogram leaf order from mid-tau
    D_mid_condensed = squareform(Ds[1], checks=False)
    Z = linkage(D_mid_condensed, method="average")
    leaf_order = leaves_list(Z)

    band_data[band] = dict(
        N=N, eigenvalues=eigenvalues, eigenvectors=eigenvectors,
        lambda_gap=lambda_gap, lambda_max=lambda_max,
        tau_min=tau_min, tau_max=tau_max, tau_mid=tau_mid,
        tau_values=tau_values, tau_labels=tau_labels,
        Ks=Ks, rhos=rhos, Ds=Ds, leaf_order=leaf_order,
    )

print("\nAll spectral objects computed.\n")


# ======================================================================
# FIG 1: LRG overview schematic (band-independent, unchanged)
# ======================================================================
print("Fig 1: lrg_overview (schematic)...")
from scipy.cluster.hierarchy import dendrogram as scipy_dendrogram

# 6 boxes equally spaced with generous gaps
# Centers at 1.5, 4.5, 7.5, 10.5, 13.5, 16.5  (spacing = 3.0)
box_w, box_h = 1.8, 1.6
box_y = 0.5
box_xs = [1.5, 4.5, 7.5, 10.5, 13.5, 16.5]

fig, ax = plt.subplots(figsize=(18, 3.5))
ax.set_xlim(0, 18)
ax.set_ylim(-1.6, 1.8)
ax.axis("off")

boxes = [
    (box_xs[0], r"$A$", "MSC\nadjacency", "#E8D5B7"),
    (box_xs[1], r"$L$", "Combinatorial\nLaplacian", "#D5E8D4"),
    (box_xs[2], r"$K(\tau)$", "Diffusion\npropagator", "#DAE8FC"),
    (box_xs[3], r"$\rho(\tau)$", "Density\nmatrix", "#E1D5E7"),
    (box_xs[4], r"$D(\tau)$", "Ultrametric\ndistance", "#FFE6CC"),
]

for x, math_label, desc, color in boxes:
    rect = FancyBboxPatch((x - box_w/2, box_y - box_h/2), box_w, box_h,
                           boxstyle="round,pad=0.1", facecolor=color,
                           edgecolor="0.3", lw=1.5)
    ax.add_patch(rect)
    ax.text(x, box_y + 0.25, math_label, ha="center", va="center",
            fontsize=18, fontweight="bold")
    ax.text(x, box_y - 0.32, desc, ha="center", va="center",
            fontsize=9.5, color="0.3", linespacing=1.3)

# Last box: dendrogram — draw box, then place inset inside
dend_x = box_xs[5]
rect = FancyBboxPatch((dend_x - box_w/2, box_y - box_h/2), box_w, box_h,
                       boxstyle="round,pad=0.1", facecolor="#F8CECC",
                       edgecolor="0.3", lw=1.5)
ax.add_patch(rect)

# Force draw so coordinate transforms are valid
fig.canvas.draw()

# Inset axes for mini dendrogram
trans = ax.transData + fig.transFigure.inverted()
pad_x, pad_bot, pad_top = 0.08, 0.08, 0.05
bl = trans.transform((dend_x - box_w/2 + pad_x, box_y - box_h/2 + pad_bot))
tr = trans.transform((dend_x + box_w/2 - pad_x, box_y + box_h/2 - pad_top))
inset = fig.add_axes([bl[0], bl[1], tr[0] - bl[0], tr[1] - bl[1]])

# Use cached LRG linkage (normalized heights 0–1, proper dendrogram)
from lrg_eegfc.workflow.lrg import load_lrg_result
lrg_res = load_lrg_result("Pat_02", "rsPre", "high_gamma", "msc")
Z_lrg = lrg_res.linkage_matrix

# Extract a sub-branch with ~15 leaves for a clean schematic tree
from scipy.cluster.hierarchy import fcluster
clusters = fcluster(Z_lrg, t=lrg_res.optimal_threshold, criterion="distance")
# Pick the largest cluster
from collections import Counter
cluster_counts = Counter(clusters)
biggest_cluster_id = cluster_counts.most_common(1)[0][0]
sub_indices = np.where(clusters == biggest_cluster_id)[0]
# If too large, take first 15
if len(sub_indices) > 20:
    sub_indices = sub_indices[:15]
# Build sub-linkage from the ultrametric condensed distances
from scipy.spatial.distance import squareform as sqf
um_sq = sqf(lrg_res.ultrametric_matrix)
um_sub = um_sq[np.ix_(sub_indices, sub_indices)]
um_sub_cond = sqf(um_sub)
Z_sub = linkage(um_sub_cond, method="average")

scipy_dendrogram(Z_sub, ax=inset, no_labels=True,
                 color_threshold=0)
inset.set_xticks([])
inset.set_yscale("log")
inset.set_ylim(Z_sub[:, 2].min() * 0.8, Z_sub[:, 2].max() * 1.1)
inset.yaxis.set_major_locator(plt.NullLocator())
inset.yaxis.set_minor_locator(plt.NullLocator())
inset.xaxis.set_major_locator(plt.NullLocator())
inset.xaxis.set_minor_locator(plt.NullLocator())
for spine in inset.spines.values():
    spine.set_visible(False)
inset.set_facecolor("#F8CECC")

# Arrows between boxes — labels below, thicker arrows
arrow_kw = dict(arrowstyle="-|>,head_length=0.6,head_width=0.3",
                color="0.4", lw=2.5, connectionstyle="arc3,rad=0")
arrow_labels = [
    r"$L = D - A$",
    r"$e^{-\tau L}$",
    r"$\div\, Z(\tau)$",
    r"$\propto 1/\rho(\tau)$",
    "cluster",
]
for i in range(5):
    x0 = box_xs[i] + box_w / 2
    x1 = box_xs[i + 1] - box_w / 2
    ax.annotate("", xy=(x1, box_y), xytext=(x0, box_y),
                arrowprops=arrow_kw)
    ax.text((x0 + x1) / 2, box_y - box_h / 2 - 0.2, arrow_labels[i],
            ha="center", va="top", fontsize=11, color="0.4", style="italic")

# Tau scale bar at bottom
ax.annotate("", xy=(16.5, -1.1), xytext=(1.5, -1.1),
            arrowprops=dict(arrowstyle="<->", color="steelblue", lw=1.5))
ax.text(9.0, -1.35, r"diffusion time $\tau$ (resolution parameter)",
        ha="center", fontsize=12, color="steelblue")

out = OUTPUT_DIR / "fig_lrg_overview.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# FIG 2: Laplacian spectrum + tau window (2 rows: high_gamma, theta)
# ======================================================================
print("Fig 2: spectrum_tau_window...")

fig, axes = plt.subplots(2, 2, figsize=(14, 9),
                         gridspec_kw={"width_ratios": [1, 1]})

for row, band in enumerate(BANDS):
    bd = band_data[band]
    N = bd["N"]
    eigenvalues = bd["eigenvalues"]
    lambda_gap = bd["lambda_gap"]
    lambda_max = bd["lambda_max"]
    tau_min = bd["tau_min"]
    tau_max = bd["tau_max"]

    # --- Left: eigenvalue spectrum ---
    ax_spec = axes[row, 0]
    ax_spec.plot(range(1, N+1), eigenvalues, "o-", color="0.3",
                 markersize=2.5, lw=0.8)
    ax_spec.axhline(lambda_gap, color="C0", ls="--", lw=1.5, alpha=0.7,
                    label=rf"$\lambda_2 = {lambda_gap:.4f}$")
    ax_spec.axhline(lambda_max, color="C3", ls="--", lw=1.5, alpha=0.7,
                    label=rf"$\lambda_{{\max}} = {lambda_max:.2f}$")
    ax_spec.set_ylabel(r"$\lambda_\ell$", fontsize=13)
    ax_spec.set_yscale("log")
    ax_spec.set_xlim(1, N)
    ax_spec.legend(fontsize=10)
    ax_spec.grid(alpha=0.3)

    if row == 0:
        ax_spec.set_title("(a)  Laplacian spectrum", fontsize=13, fontweight="bold")
    if row == 1:
        ax_spec.set_xlabel(r"Eigenvalue index $\ell$", fontsize=12)

    # Band label on the left
    ax_spec.annotate(BRAIN_BAND_TEX_DICT[band], xy=(-0.22, 0.5),
                     xycoords="axes fraction", fontsize=16, fontweight="bold",
                     ha="center", va="center", rotation=90)

    # --- Right: partition function Z(tau) ---
    ax_tau = axes[row, 1]
    tau_scan = np.logspace(np.log10(tau_min * 0.3), np.log10(tau_max * 3), 300)
    Z_vals = np.array([np.sum(np.exp(-t * eigenvalues)) for t in tau_scan])
    ax_tau.plot(tau_scan, Z_vals, color="0.3", lw=1.5, label=r"$Z(\tau)$")

    ax_tau.axvline(tau_min, color="C3", ls="--", lw=1.5,
                   label=rf"$\tau' = 1/\lambda_{{\max}} = {tau_min:.4f}$")
    ax_tau.axvline(tau_max, color="C0", ls="--", lw=1.5,
                   label=rf"$1/\lambda_2 = {tau_max:.2f}$")
    ax_tau.axvspan(tau_min, tau_max, alpha=0.1, color="steelblue",
                   label="Relevant window")

    ax_tau.set_ylabel(r"$Z(\tau)$", fontsize=13)
    ax_tau.set_xscale("log")
    ax_tau.set_yscale("log")
    ax_tau.legend(fontsize=9, loc="upper right")
    ax_tau.grid(alpha=0.3)

    if row == 0:
        ax_tau.set_title("(b)  Partition function and scale window",
                         fontsize=13, fontweight="bold")
    if row == 1:
        ax_tau.set_xlabel(r"Diffusion time $\tau$", fontsize=12)

plt.tight_layout(rect=[0.03, 0, 1, 1])
out = OUTPUT_DIR / "fig_spectrum_tau_window.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# FIG 3: rho(tau) heatmaps (2 rows x 3 cols)
# ======================================================================
print("Fig 3: rho_tau...")
from matplotlib.colors import LogNorm

fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)

for row, band in enumerate(BANDS):
    bd = band_data[band]
    rhos = bd["rhos"]
    leaf_order = bd["leaf_order"]
    tau_labels = bd["tau_labels"]

    for k in range(3):
        ax = axes[row, k]
        rho_ord = rhos[k][np.ix_(leaf_order, leaf_order)]
        # Clamp zeros for log scale
        rho_pos = np.where(rho_ord > 0, rho_ord, np.nan)
        im = ax.imshow(rho_pos, cmap=CMAP_RHO, aspect="equal",
                       interpolation="none", norm=LogNorm())
        ax.set_xticks([])
        ax.set_yticks([])

        # Tau label below each column (same for both rows)
        ax.set_xlabel(tau_labels[k], fontsize=14)

        if k == 0:
            ax.set_ylabel(BRAIN_BAND_TEX_DICT[band], fontsize=18,
                          fontweight="bold")

    # One colorbar per row
    cbar = fig.colorbar(im, ax=axes[row, :], shrink=0.85, pad=0.02, aspect=25)
    cbar.set_label(r"$\rho_{ij}(\tau)$  (log scale)", fontsize=13)
    cbar.ax.tick_params(labelsize=11)

out = OUTPUT_DIR / "fig_rho_tau.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


print(f"\nDone! 3 figures saved to {OUTPUT_DIR}")
