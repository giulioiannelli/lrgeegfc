#!/usr/bin/env python3
"""C1 PSI-overlay validation figure.

Overlays the sensible PSI peak positions (from the full unthresholded network)
onto the C(tau) curves from the thresholded networks.

For each band, the sensible n* values are obtained from the PSI profile at
tau_min of the full MSC matrix.  The corresponding tau positions are read from
the precomputed coarsening trajectories: for each n*, the first tau where
n_max(tau) drops to that level marks the scale at which that community
structure becomes dominant.

Vertical dashed red lines at these tau positions are overlaid on the threshold
C(tau) curves.  If threshold-induced peaks in C(tau) align with the red lines,
the PSI community levels from the full network explain the emerging peaks.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from scipy.signal import find_peaks

# -- project imports -----------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
# Ensure we can import from the project
sys.path.insert(0, "src")

from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.workflow.diagnostics import (
    _compute_propagator,
    _ultrametric_from_propagator,
    sensible_psi_peaks,
)
from lrg_eegfc.visuals.lrg import compute_partition_stability_index
from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT

# -- constants -----------------------------------------------------------------
BANDS = list(BRAIN_BANDS_NAMES)
PATIENT = "Pat_02"
PHASE = "rsPre"

# -- load precomputed data -----------------------------------------------------
threshold_data = dict(
    np.load("data/tables/threshold_analysis_Pat02.npz", allow_pickle=True)
)
percentiles = threshold_data["percentiles"]

traj_data = dict(
    np.load("data/tables/mslcd_coarsening_trajectories.npz", allow_pickle=True)
)

# -- colour maps ---------------------------------------------------------------
thresh_cmap = cm.get_cmap("plasma", len(percentiles))

_RC = {
    "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 14,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 10,
    "text.usetex": False, "mathtext.fontset": "cm",
}

# ==============================================================================
# Main
# ==============================================================================
with plt.rc_context(_RC):
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)

    for k, band in enumerate(BANDS):
        ax = axes.ravel()[k]

        # ------------------------------------------------------------------
        # 1. Load full MSC matrix -> Laplacian -> eigendecomposition
        # ------------------------------------------------------------------
        A = load_msc_matrix(
            PATIENT, PHASE, band,
            cache_root=Path("data/msc_cache"),
            sparsify="none", n_surrogates=0, nperseg=4096,
        )
        if A is None:
            ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                    transform=ax.transAxes)
            print(f"\n=== {band} ===  SKIP (no MSC cache)")
            continue

        np.fill_diagonal(A, 0)
        N = A.shape[0]
        D_diag = A.sum(axis=1)
        L = np.diag(D_diag) - A
        eigenvalues, eigenvectors = np.linalg.eigh(L)
        eigenvalues = np.maximum(eigenvalues, 0.0)

        tau_min = 1.0 / eigenvalues[-1]
        tau_max_val = 1.0 / eigenvalues[1]

        # ------------------------------------------------------------------
        # 2. PSI profile at tau_min -> sensible n* values
        # ------------------------------------------------------------------
        K = _compute_propagator(eigenvalues, eigenvectors, tau_min)
        D_ultra = _ultrametric_from_propagator(K)
        D_cond = squareform(D_ultra, checks=False)
        Z = linkage(D_cond, method="average")
        psi_vals, n_comms = compute_partition_stability_index(Z)
        kept_idx, n_sens = sensible_psi_peaks(psi_vals, n_comms)
        sensible_n_vals = sorted([int(n_comms[ki]) for ki in kept_idx])

        # ------------------------------------------------------------------
        # 3. Map each sensible n* to its tau position via coarsening trajectory
        # ------------------------------------------------------------------
        tau_key = f"{PATIENT}__{PHASE}__{band}__tau"
        nmax_key = f"{PATIENT}__{PHASE}__{band}__nmax"

        psi_tau_positions = []
        if tau_key in traj_data:
            traj_tau = traj_data[tau_key]
            traj_nmax = traj_data[nmax_key]

            for n_star in sensible_n_vals:
                # Find the first tau where n_max drops to <= n_star
                idx = np.where(traj_nmax <= n_star)[0]
                if len(idx) > 0:
                    psi_tau_positions.append((n_star, float(traj_tau[idx[0]])))
                else:
                    # n_star is larger than all n_max values -- use tau_min
                    psi_tau_positions.append((n_star, float(traj_tau[0])))
        else:
            print(f"  WARNING: no coarsening trajectory for {band}")

        # ------------------------------------------------------------------
        # 4. Print diagnostics
        # ------------------------------------------------------------------
        print(f"\n=== {band} ===")
        print(f"  N={N}, lambda_2={eigenvalues[1]:.4f}, "
              f"lambda_max={eigenvalues[-1]:.4f}")
        print(f"  tau_min={tau_min:.6f}, tau_max={tau_max_val:.6f}")
        print(f"  Sensible n* = {sensible_n_vals}")
        print(f"  PSI-predicted tau positions:")
        for n_star, tau_pos in psi_tau_positions:
            print(f"    n*={n_star:3d} -> tau={tau_pos:.6f} "
                  f"(log10={np.log10(tau_pos):.3f})")

        # ------------------------------------------------------------------
        # 5. Check for C(tau) peak alignment at each threshold
        # ------------------------------------------------------------------
        tau_key_t = f"{band}__tau"
        if tau_key_t not in threshold_data:
            ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                    transform=ax.transAxes)
            continue
        tau = threshold_data[tau_key_t]
        log10_tau = np.log10(tau)

        print(f"  C(tau) peak alignment check:")
        for ti, pct in enumerate(percentiles):
            c_key = f"{band}__C__{pct}"
            if c_key not in threshold_data:
                continue
            C = threshold_data[c_key]
            # Find peaks in C(tau)
            prom = max(0.05 * C.max(), 1e-6) if C.max() > 0 else 1e-6
            c_peaks, props = find_peaks(C, prominence=prom)

            if len(c_peaks) > 0 and len(psi_tau_positions) > 0:
                c_peak_taus = tau[c_peaks]
                alignments = []
                for n_star, tau_pos in psi_tau_positions:
                    # Check if any C(tau) peak is within 0.3 decades of
                    # the PSI-predicted tau
                    log_dists = np.abs(np.log10(c_peak_taus) - np.log10(tau_pos))
                    min_dist = log_dists.min()
                    closest_peak_tau = c_peak_taus[log_dists.argmin()]
                    aligned = min_dist < 0.3
                    alignments.append(
                        f"n*={n_star}:{'YES' if aligned else 'no '} "
                        f"(d={min_dist:.2f} dec, "
                        f"nearest_peak_tau={closest_peak_tau:.4f})"
                    )
                print(f"    pct={pct:2d}%: "
                      f"{len(c_peaks)} peaks; {'; '.join(alignments)}")

        # ------------------------------------------------------------------
        # 6. Plot C(tau) curves from threshold data
        # ------------------------------------------------------------------
        for ti, pct in enumerate(percentiles):
            c_key = f"{band}__C__{pct}"
            if c_key not in threshold_data:
                continue
            C = threshold_data[c_key]
            ax.plot(
                tau, C, lw=1.5,
                color=thresh_cmap(ti / max(len(percentiles) - 1, 1)),
                label=f"{pct}%" if k == 0 else None,
                alpha=0.8,
            )

        # ------------------------------------------------------------------
        # 7. Overlay PSI-predicted tau positions as vertical lines
        # ------------------------------------------------------------------
        ymax = ax.get_ylim()[1]
        for n_star, tau_pos in psi_tau_positions:
            ax.axvline(tau_pos, color="red", ls="--", lw=1.2, alpha=0.7)
            ax.text(
                tau_pos, ymax * 0.97,
                f"  $n^*$={n_star}",
                fontsize=7, color="red", rotation=90,
                va="top", ha="left",
            )

        # ------------------------------------------------------------------
        # 8. Axis formatting
        # ------------------------------------------------------------------
        ax.set_xscale("log")
        band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
        ax.text(0.04, 0.94, band_tex, transform=ax.transAxes,
                fontsize=13, fontweight="bold", va="top")
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.3, which="both")
        if k >= 3:
            ax.set_xlabel(r"$\tau$")
        if k % 3 == 0:
            ax.set_ylabel(r"$C(\tau)$")

    # -- Shared colorbar -------------------------------------------------------
    sm = cm.ScalarMappable(
        cmap=thresh_cmap,
        norm=plt.Normalize(vmin=percentiles[0], vmax=percentiles[-1]),
    )
    sm.set_array([])
    cb = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.02, pad=0.02)
    cb.set_label("Edge-weight percentile threshold (%)", fontsize=11)

    fig.suptitle(
        f"{PATIENT} / {PHASE} -- C1 with PSI-predicted community levels",
        fontsize=14, fontweight="bold",
    )

    # -- Save -------------------------------------------------------------------
    out = Path(
        "data/figures/report_mslcd_section/threshold_analysis/"
        "fig_C1_psi_overlay_validation.png"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", dpi=200)
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {out}")
    print(f"Saved: {out.with_suffix('.pdf')}")
