#!/usr/bin/env python3
"""Generate metastable nodes brain figures for all phase/band combos.

Uses data-driven metastability: at each τ in [1/λ_max, 1/λ_2], compute
the set of sensible PSI peaks S(τ), take n_max(τ) = max S(τ) (capped at
N/2 to exclude single-node splitting noise), cut the dendrogram there,
align labels via Hungarian matching, and compute μ_i across the sequence.

Run: python scripts/gen_metastable_grid.py
"""
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.signal import find_peaks
from scipy.optimize import linear_sum_assignment

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PHASE_LABELS
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.visuals.lrg import compute_partition_stability_index
from lrg_eegfc.visuals.spatial import load_spatial_metadata, prepare_spatial_coordinates

# ── Config ────────────────────────────────────────────────────────────
from lrg_eegfc.config.paths import MSC_CACHE, SEEG_DATAPATH, FIGURES_ROOT
PATIENT = "Pat_02"
DATASET_ROOT = SEEG_DATAPATH
OUTPUT_DIR = FIGURES_ROOT / "report_mslcd_section" / PATIENT / "metastable_nodes"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# All phases x all bands
COMBOS = [(phase, band) for phase in PHASE_LABELS for band in BRAIN_BANDS_NAMES]

N_TAU_META = 20  # log-spaced τ values in the resolution window


# ── Helper functions ──────────────────────────────────────────────────
def compute_propagator(eigenvalues, eigenvectors, tau):
    exp_vals = np.exp(-tau * eigenvalues)
    return (eigenvectors * exp_vals[None, :]) @ eigenvectors.T


def compute_ultrametric_distance(K):
    D = np.zeros_like(K)
    mask = ~np.eye(K.shape[0], dtype=bool)
    D[mask] = 1.0 / np.where(K[mask] > 1e-30, K[mask], 1e-30)
    return D


def sensible_psi_peaks(psi_vals, n_comms, log_drop_max=0.5, min_n_spacing=5):
    """Find PSI peaks using mean-floor + log-drop chain rule."""
    if len(psi_vals) == 0:
        return [], 0
    peak_idx, _ = find_peaks(psi_vals)
    if len(peak_idx) == 0:
        return [np.argmax(psi_vals)], 1
    psi_mean = psi_vals.mean()
    peak_idx = peak_idx[psi_vals[peak_idx] > psi_mean]
    if len(peak_idx) == 0:
        return [np.argmax(psi_vals)], 1
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


def compute_metastability_datadriven(A, N):
    """Data-driven metastability: n_max(τ) = max S(τ), capped at N/2.

    At each of N_TAU_META log-spaced τ in [1/λ_max, 1/λ_2], compute the
    sensible PSI peaks, take the largest (capped at N/2), cut the
    dendrogram there, align labels, and compute μ_i.

    Returns (mu, nmax_trajectory, tau_grid) or (None, None, None) on failure.
    """
    D_diag = A.sum(axis=1)
    L = np.diag(D_diag) - A
    eigenvalues, eigenvectors = np.linalg.eigh(L)
    eigenvalues = np.maximum(eigenvalues, 0.0)

    if eigenvalues[1] < 1e-10:
        return None, None, None  # disconnected graph

    tau_min = 1.0 / eigenvalues[-1]
    tau_max = 1.0 / eigenvalues[1]
    nmax_cap = N // 2

    tau_grid = np.logspace(np.log10(tau_min), np.log10(tau_max), N_TAU_META)
    nmax_traj = np.zeros(N_TAU_META, dtype=int)
    communities = []

    for i, tau in enumerate(tau_grid):
        K = compute_propagator(eigenvalues, eigenvectors, tau)
        D = compute_ultrametric_distance(K)
        D_cond = squareform(D, checks=False)
        Z = linkage(D_cond, method="average")
        psi_vals, n_comms = compute_partition_stability_index(Z)

        peaks_idx, n_sensible = sensible_psi_peaks(psi_vals, n_comms)
        if n_sensible > 0:
            sensible_n = [int(n_comms[pi]) for pi in peaks_idx
                          if int(n_comms[pi]) <= nmax_cap]
            nmax_traj[i] = max(sensible_n) if sensible_n else 1
        else:
            nmax_traj[i] = 1

        labels = fcluster(Z, max(nmax_traj[i], 1), criterion="maxclust")
        communities.append(labels)

    communities_aligned = align_cluster_labels(communities)

    n_transitions = N_TAU_META - 1
    mu = np.zeros(N)
    if n_transitions > 0:
        for s in range(n_transitions):
            mu += (communities_aligned[s] != communities_aligned[s + 1]).astype(float)
        mu /= n_transitions

    return mu, nmax_traj, tau_grid


# ── Spatial coordinates (shared) ──────────────────────────────────────
metadata = load_spatial_metadata(PATIENT, DATASET_ROOT)
coords_mni = prepare_spatial_coordinates(metadata, scale="mm", to_mni=True)
N = coords_mni.shape[0]
print(f"{PATIENT}: {N} channels, {len(COMBOS)} combos to generate\n")


# ── Generate figures ──────────────────────────────────────────────────
try:
    from nilearn import plotting as ni_plot
    USE_NILEARN = True
except ImportError:
    USE_NILEARN = False
    print("  nilearn not available, skipping")

for phase, band in COMBOS:
    print(f"  {phase}/{band}...", end=" ")

    A = load_msc_matrix(PATIENT, phase, band,
                        cache_root=MSC_CACHE, sparsify="none",
                        n_surrogates=0, nperseg=4096)
    if A is None:
        print("SKIP (no MSC cache)")
        continue
    np.fill_diagonal(A, 0)
    assert A.shape[0] == N

    mu, nmax_traj, tau_grid = compute_metastability_datadriven(A, N)
    if mu is None:
        print("SKIP (disconnected)")
        continue

    n_meta = int(np.sum(mu > 0))
    band_tex = BRAIN_BAND_TEX_DICT[band]

    if USE_NILEARN:
        display = ni_plot.plot_markers(
            node_values=mu,
            node_coords=coords_mni,
            node_cmap="RdYlGn_r",
            node_size=mu * 80 + 15,
            display_mode="lzry",
            title=(
                rf"{PATIENT}, {phase}, {band_tex} — "
                rf"Metastability $\mu_i$ "
                rf"({N_TAU_META} $\tau$ in "
                rf"$[1/\lambda_{{\max}},\,1/\lambda_2]$, "
                rf"$n_{{\max}}$: {nmax_traj[0]}$\to${nmax_traj[-1]})"
            ),
            colorbar=True,
        )
        out = OUTPUT_DIR / f"fig_metastable_nodes_{phase}_{band}.pdf"
        display.savefig(str(out), dpi=150)
        display.close()
        print(f"μ>0: {n_meta}/{N}, n_max: {nmax_traj[0]}→{nmax_traj[-1]}")
    else:
        print("SKIP (nilearn)")

print(f"\nDone! Figures in {OUTPUT_DIR}")
