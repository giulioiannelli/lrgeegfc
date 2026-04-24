#!/usr/bin/env python3
"""Generate cross-comparison figures for Section 4: Multiscale Community Detection.

Figures 7–11: comparisons across patients, phases, and bands.

  Fig 7:  fig_C_tau_across_patients.pdf     — C(τ) for all patients (fixed band/phase)
  Fig 8:  fig_nstar_across_patients.pdf     — n*(τ) for all patients (fixed band/phase)
  Fig 9:  fig_metastability_across_patients.pdf — μ_i distributions per patient
  Fig 10: fig_C_tau_across_phases.pdf       — C(τ) across 4 phases (Pat_02, fixed band)
  Fig 11: fig_C_tau_across_bands.pdf        — C(τ) across 6 bands (Pat_02, rest_pre)

Run: python scripts/gen_mslcd_cross_figures.py
"""
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.signal import find_peaks
from scipy.optimize import linear_sum_assignment

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PHASE_LABELS,
)
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.visuals.lrg import compute_partition_stability_index

# ── Config ────────────────────────────────────────────────────────────
from lrg_eegfc.config.paths import MSC_CACHE, LRG_CACHE, FIGURES_ROOT
NPERSEG = 4096
FC_METHOD = "msc"
OUTPUT_DIR = FIGURES_ROOT / "report_mslcd_section" / "cross"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Patients with full MSC+LRG at nperseg=4096
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
PHASE_COLORS = {"rest_pre": "#4C72B0", "task_learn": "#55A868",
                "task_test": "#C44E52", "rest_post": "#8172B3"}
PATIENT_COLORS = {
    "Pat_02": "#4C72B0", "Pat_03": "#55A868", "Pat_05": "#C44E52",
    "Pat_07": "#8172B3", "Pat_08": "#CCB974",
}
BAND_COLORS = {
    "delta": "#4C72B0", "theta": "#55A868", "alpha": "#C44E52",
    "beta": "#8172B3", "low_gamma": "#CCB974", "high_gamma": "#64B5CD",
}

# Fixed choices for cross-comparisons
CROSS_BAND = "beta"        # for cross-patient and cross-phase
CROSS_PHASE = "rest_pre"      # for cross-patient and cross-band
CROSS_PATIENT = "Pat_02"   # for cross-phase and cross-band


# ── Helper functions ──────────────────────────────────────────────────
def entropy_at_tau(eigenvalues, tau):
    """S(τ) = -Σ p_ℓ ln p_ℓ, numerically stable."""
    log_boltz = -tau * eigenvalues
    log_Z = np.logaddexp.reduce(log_boltz)
    log_p = log_boltz - log_Z
    p = np.exp(log_p)
    return -np.sum(p * np.where(p > 1e-300, log_p, 0))


def compute_entropy_curve(eigenvalues, n_points=500, pad_factor=0.3):
    """Compute S(τ), C(τ), and resolution window for a set of eigenvalues."""
    eigenvalues = np.maximum(eigenvalues, 0.0)
    lambda_gap = eigenvalues[1]
    lambda_max = eigenvalues[-1]
    tau_min = 1.0 / lambda_max
    tau_max = 1.0 / lambda_gap

    tau_grid = np.logspace(
        np.log10(tau_min * pad_factor),
        np.log10(tau_max * (1.0 / pad_factor)),
        n_points,
    )
    S = np.array([entropy_at_tau(eigenvalues, t) for t in tau_grid])
    log10_tau = np.log10(tau_grid)
    C = -np.gradient(S, log10_tau)

    return dict(
        tau=tau_grid, log10_tau=log10_tau, S=S, C=C,
        tau_min=tau_min, tau_max=tau_max,
        lambda_gap=lambda_gap, lambda_max=lambda_max,
    )


def load_eigenvalues(patient, phase, band):
    """Load MSC matrix and return Laplacian eigenvalues."""
    A = load_msc_matrix(patient, phase, band,
                        cache_root=MSC_CACHE, sparsify="none",
                        n_surrogates=0, nperseg=NPERSEG)
    if A is None:
        return None
    np.fill_diagonal(A, 0)
    D_diag = A.sum(axis=1)
    L = np.diag(D_diag) - A
    eigenvalues, _ = np.linalg.eigh(L)
    return np.maximum(eigenvalues, 0.0)


def compute_nstar_curve(eigenvalues, n_tau=60):
    """Compute n*(τ) across the resolution window."""
    eigenvalues = np.maximum(eigenvalues, 0.0)
    lambda_gap = eigenvalues[1]
    lambda_max = eigenvalues[-1]
    tau_min = 1.0 / lambda_max
    tau_max = 1.0 / lambda_gap
    N = len(eigenvalues)

    tau_grid = np.logspace(np.log10(tau_min), np.log10(tau_max), n_tau)
    nstar = np.zeros(n_tau, dtype=int)

    eigenvectors = None  # lazy compute
    for i, tau in enumerate(tau_grid):
        if eigenvectors is None:
            # Need eigenvectors for propagator
            D_diag_recon = np.zeros(N)  # dummy — recompute from scratch
            # Actually just recompute eigh to get eigenvectors
            pass
        # Use cached LRG linkage if available, else compute from eigenvalues
        # For speed, use the entropy-based approach: PSI from LRG linkage
        pass

    return tau_grid, nstar


def compute_nstar_from_lrg(patient, phase, band, n_tau=60):
    """Compute n*(τ) using eigendecomposition + PSI."""
    A = load_msc_matrix(patient, phase, band,
                        cache_root=MSC_CACHE, sparsify="none",
                        n_surrogates=0, nperseg=NPERSEG)
    if A is None:
        return None, None, None
    np.fill_diagonal(A, 0)
    N = A.shape[0]
    D_diag = A.sum(axis=1)
    L_mat = np.diag(D_diag) - A
    eigenvalues, eigenvectors = np.linalg.eigh(L_mat)
    eigenvalues = np.maximum(eigenvalues, 0.0)

    tau_min = 1.0 / eigenvalues[-1]
    tau_max = 1.0 / eigenvalues[1]
    tau_grid = np.logspace(np.log10(tau_min), np.log10(tau_max), n_tau)
    nstar = np.zeros(n_tau, dtype=int)

    for i, tau in enumerate(tau_grid):
        exp_vals = np.exp(-tau * eigenvalues)
        K = (eigenvectors * exp_vals[None, :]) @ eigenvectors.T
        D_ultra = np.zeros_like(K)
        mask = ~np.eye(N, dtype=bool)
        D_ultra[mask] = 1.0 / np.where(K[mask] > 1e-30, K[mask], 1e-30)
        D_cond = squareform(D_ultra, checks=False)
        Z = linkage(D_cond, method="average")
        psi_vals, n_comms = compute_partition_stability_index(Z)
        if len(psi_vals) > 0:
            nstar[i] = int(n_comms[np.argmax(psi_vals)])
        else:
            nstar[i] = 1

    return tau_grid, nstar, eigenvalues


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


def compute_metastability(patient, phase, band, tau_nclust_pairs):
    """Compute metastability scores for a patient/phase/band."""
    A = load_msc_matrix(patient, phase, band,
                        cache_root=MSC_CACHE, sparsify="none",
                        n_surrogates=0, nperseg=NPERSEG)
    if A is None:
        return None
    np.fill_diagonal(A, 0)
    N = A.shape[0]

    G = nx.from_numpy_array(A)
    Gcc_nodes = max(nx.connected_components(G), key=len)
    Gcc = G.subgraph(Gcc_nodes).copy()

    from lrg_eegfc.visuals.metastable import compute_clustering_across_tau
    taus = [t for t, _ in tau_nclust_pairs]
    nclust = [n for _, n in tau_nclust_pairs]

    partitions, _ = compute_clustering_across_tau(
        None, np.array(taus), Gcc, n_clusters_list=nclust,
    )

    communities_raw = [partitions[t] for t in taus]
    communities_aligned = align_cluster_labels(communities_raw)

    n_transitions = len(taus) - 1
    mu = np.zeros(N)
    if n_transitions > 0:
        for s in range(n_transitions):
            mu += (communities_aligned[s] != communities_aligned[s + 1]).astype(float)
        mu /= n_transitions

    return mu


# ======================================================================
# PRECOMPUTE: Eigenvalues for all needed triplets
# ======================================================================
print("=" * 60)
print("Precomputing eigenvalues...")
print("=" * 60)

# Cache eigenvalues to avoid recomputation
eigen_cache = {}


def get_eigenvalues(patient, phase, band):
    key = (patient, phase, band)
    if key not in eigen_cache:
        eig = load_eigenvalues(patient, phase, band)
        if eig is not None:
            eigen_cache[key] = eig
            print(f"  {patient}/{phase}/{band}: N={len(eig)}")
        else:
            print(f"  {patient}/{phase}/{band}: MISSING")
            return None
    return eigen_cache[key]


# Preload for cross-patient (fixed band + phase)
print(f"\nCross-patient: {CROSS_PHASE}, {CROSS_BAND}")
for patient in PATIENTS:
    get_eigenvalues(patient, CROSS_PHASE, CROSS_BAND)

# Preload for cross-phase (fixed patient + band)
print(f"\nCross-phase: {CROSS_PATIENT}, {CROSS_BAND}")
for phase in PHASES:
    get_eigenvalues(CROSS_PATIENT, phase, CROSS_BAND)

# Preload for cross-band (fixed patient + phase)
print(f"\nCross-band: {CROSS_PATIENT}, {CROSS_PHASE}")
for band in BRAIN_BANDS_NAMES:
    get_eigenvalues(CROSS_PATIENT, CROSS_PHASE, band)

print(f"\n{len(eigen_cache)} triplets loaded.\n")


# ======================================================================
# FIG 7: C(τ) across patients (fixed band, fixed phase)
# ======================================================================
print("Fig 7: C(tau) across patients...")

fig, ax = plt.subplots(figsize=(12, 6))

for patient in PATIENTS:
    eig = get_eigenvalues(patient, CROSS_PHASE, CROSS_BAND)
    if eig is None:
        continue
    ec = compute_entropy_curve(eig)
    ax.plot(ec["log10_tau"], ec["C"], lw=2, color=PATIENT_COLORS[patient],
            label=f"{patient} (N={len(eig)})")

    # Mark resolution window boundaries
    ax.axvline(np.log10(ec["tau_min"]), color=PATIENT_COLORS[patient],
               ls=":", lw=0.8, alpha=0.4)
    ax.axvline(np.log10(ec["tau_max"]), color=PATIENT_COLORS[patient],
               ls=":", lw=0.8, alpha=0.4)

ax.set_xlabel(r"$\log_{10}\,\tau$", fontsize=14)
ax.set_ylabel(r"Entropic susceptibility $C(\tau)$", fontsize=14)
ax.legend(fontsize=11)
ax.grid(alpha=0.3)
ax.set_ylim(bottom=0)
ax.tick_params(labelsize=12)

out = OUTPUT_DIR / "fig_C_tau_across_patients.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# FIG 8: n*(τ) across patients (fixed band, fixed phase)
# ======================================================================
print("Fig 8: n*(tau) across patients...")

fig, ax = plt.subplots(figsize=(12, 6))

for patient in PATIENTS:
    eig = get_eigenvalues(patient, CROSS_PHASE, CROSS_BAND)
    if eig is None:
        continue
    print(f"  Computing n*(tau) for {patient}...")
    tau_grid, nstar, _ = compute_nstar_from_lrg(patient, CROSS_PHASE, CROSS_BAND,
                                                 n_tau=40)
    if tau_grid is None:
        continue
    ax.step(np.log10(tau_grid), nstar, where="mid", lw=2,
            color=PATIENT_COLORS[patient],
            label=f"{patient} (N={len(eig)})")

ax.set_xlabel(r"$\log_{10}\,\tau$", fontsize=14)
ax.set_ylabel(r"Optimal communities $n^*(\tau)$", fontsize=14)
ax.legend(fontsize=11)
ax.grid(alpha=0.3)
ax.set_ylim(bottom=0)
ax.tick_params(labelsize=12)

out = OUTPUT_DIR / "fig_nstar_across_patients.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# FIG 9: Metastability distributions across patients
# ======================================================================
print("Fig 9: metastability across patients...")

TAU_NCLUST_PAIRS = [
    (0.1, 15), (0.3, 12), (0.5, 8), (0.8, 6),
    (1.0, 5), (2.0, 3), (5.0, 2),
]

fig, ax = plt.subplots(figsize=(10, 6))
mu_data = []
mu_labels = []

for patient in PATIENTS:
    print(f"  Computing metastability for {patient}...")
    mu = compute_metastability(patient, CROSS_PHASE, CROSS_BAND, TAU_NCLUST_PAIRS)
    if mu is not None:
        mu_data.append(mu)
        mu_labels.append(f"{patient}\n(N={len(mu)})")

bp = ax.boxplot(mu_data, labels=mu_labels, patch_artist=True,
                showfliers=True, widths=0.6,
                flierprops=dict(marker=".", markersize=3, alpha=0.5),
                medianprops=dict(color="black", lw=2))
for patch, patient in zip(bp["boxes"], PATIENTS[:len(mu_data)]):
    patch.set_facecolor(PATIENT_COLORS[patient])
    patch.set_alpha(0.7)

ax.set_ylabel(r"Metastability score $\mu_i$", fontsize=14)
ax.set_xlabel("Patient", fontsize=14)
ax.grid(axis="y", alpha=0.3)
ax.tick_params(labelsize=12)

# Annotate fraction of metastable nodes
for i, mu in enumerate(mu_data):
    frac = np.mean(mu > 0)
    ax.text(i + 1, ax.get_ylim()[1] * 0.95, f"{frac:.0%}",
            ha="center", fontsize=10, color="0.3")

out = OUTPUT_DIR / "fig_metastability_across_patients.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# FIG 10: C(τ) across phases (Pat_02, fixed band)
# ======================================================================
print("Fig 10: C(tau) across phases...")

fig, ax = plt.subplots(figsize=(12, 6))

for phase in PHASES:
    eig = get_eigenvalues(CROSS_PATIENT, phase, CROSS_BAND)
    if eig is None:
        continue
    ec = compute_entropy_curve(eig)
    ax.plot(ec["log10_tau"], ec["C"], lw=2, color=PHASE_COLORS[phase],
            label=phase)

    ax.axvline(np.log10(ec["tau_min"]), color=PHASE_COLORS[phase],
               ls=":", lw=0.8, alpha=0.4)
    ax.axvline(np.log10(ec["tau_max"]), color=PHASE_COLORS[phase],
               ls=":", lw=0.8, alpha=0.4)

ax.set_xlabel(r"$\log_{10}\,\tau$", fontsize=14)
ax.set_ylabel(r"Entropic susceptibility $C(\tau)$", fontsize=14)
ax.legend(fontsize=12)
ax.grid(alpha=0.3)
ax.set_ylim(bottom=0)
ax.tick_params(labelsize=12)

out = OUTPUT_DIR / "fig_C_tau_across_phases.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# FIG 11: C(τ) across bands (Pat_02, rest_pre)
# ======================================================================
print("Fig 11: C(tau) across bands...")

fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)

for k, band in enumerate(BRAIN_BANDS_NAMES):
    ax = axes.ravel()[k]
    eig = get_eigenvalues(CROSS_PATIENT, CROSS_PHASE, band)
    if eig is None:
        continue
    ec = compute_entropy_curve(eig)

    ax.plot(ec["log10_tau"], ec["C"], lw=2, color=BAND_COLORS[band])

    # Mark resolution window
    ax.axvline(np.log10(ec["tau_min"]), color="C3", ls="--", lw=1, alpha=0.5)
    ax.axvline(np.log10(ec["tau_max"]), color="C0", ls="--", lw=1, alpha=0.5)
    ax.axvspan(np.log10(ec["tau_min"]), np.log10(ec["tau_max"]),
               alpha=0.06, color="steelblue")

    # Find and mark peaks
    in_window = (ec["tau"] >= ec["tau_min"]) & (ec["tau"] <= ec["tau_max"])
    C_masked = ec["C"].copy()
    C_masked[~in_window] = 0
    peaks_idx, _ = find_peaks(C_masked, prominence=0.05 * C_masked.max())
    for pi in peaks_idx:
        ax.plot(ec["log10_tau"][pi], ec["C"][pi], "v", color="crimson",
                markersize=8, zorder=5)

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=16, fontweight="bold")
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.3)
    ax.tick_params(labelsize=11)

    if k >= 3:
        ax.set_xlabel(r"$\log_{10}\,\tau$", fontsize=13)
    if k % 3 == 0:
        ax.set_ylabel(r"$C(\tau)$", fontsize=13)

out = OUTPUT_DIR / "fig_C_tau_across_bands.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


print(f"\nDone! All cross-comparison figures saved to {OUTPUT_DIR}")
