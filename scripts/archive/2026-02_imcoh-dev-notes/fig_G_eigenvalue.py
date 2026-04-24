#!/usr/bin/env python3
"""Section 2 Figure G: Laplacian eigenvalue spectrum comparison.

Pat_02, all 6 bands, rest_pre. Two rows (MSC top, ImCoh bottom) x 6 columns.
All 4 phases overlaid as colored curves per panel.
Marks lambda_2 (red triangle) and lambda_max (grey triangle).

Run:
  python scripts/10_notes_imcoh/fig_G_eigenvalue.py [-v]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import sys
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib.pyplot as plt

from lrg_eegfc.utils.metrics.hypothesis import (
    BANDS, BAND_TEX, PHASE_LABELS,
    PHASE_COLORS, PHASE_MARKERS,
    CLR_MSC, CLR_IMCOH,
    REPR_PATIENTS,
    load_fc_for_patient, compute_laplacian_eigenvalues,
    apply_pub_style, save_fig, SECTION2_ROOT,
    SECTION2_METHODS as METHODS,
    SECTION2_METHOD_LABELS as METHOD_LABELS
)



def fig_g_eigenvalue_spectrum(
    patient: str = "Pat_02",
    output_dir: Path = SECTION2_ROOT / "fig_G",
    verbose: bool = False,
):
    """2-row x 6-col eigenvalue spectrum: MSC (top) vs |ImCoh| (bottom)."""

    fig, axes = plt.subplots(2, 6, figsize=(28, 8), sharex=True)

    for row, method in enumerate(METHODS):
        data = load_fc_for_patient(patient, method)

        for col, band in enumerate(BANDS):
            ax = axes[row, col]

            for phase in PHASE_LABELS:
                mat = data[band].get(phase)
                if mat is None:
                    continue

                eigvals = compute_laplacian_eigenvalues(mat)
                N = len(eigvals)
                x = np.arange(1, N + 1) / N  # normalized index

                ax.plot(x, eigvals, color=PHASE_COLORS[phase],
                        lw=1.2, alpha=0.8, label=phase if col == 0 else "")

                # Mark lambda_2 and lambda_max for rest_pre only
                if phase == "rest_pre":
                    ax.plot(2 / N, eigvals[1], "^", color="red",
                            markersize=6, zorder=10)
                    ax.plot(1.0, eigvals[-1], "v", color="gray",
                            markersize=6, zorder=10)

                    if verbose and row == 0:
                        print(f"    {method} {band} {phase}: "
                              f"λ₂={eigvals[1]:.4f}, λ_max={eigvals[-1]:.1f}")

            if row == 0:
                ax.set_title(BAND_TEX[col], fontsize=14, fontweight="bold")
            if col == 0:
                clr = CLR_MSC if method == "msc" else CLR_IMCOH
                ax.set_ylabel(f"{METHOD_LABELS[method]}\n$\\lambda$",
                              fontsize=12, color=clr, fontweight="bold")
            if row == 1:
                ax.set_xlabel("$\\ell / N$", fontsize=10)

            ax.grid(True, alpha=0.2)
            ax.tick_params(labelsize=8)

    # Add legend to first panel
    axes[0, 0].legend(fontsize=8, loc="upper left")

    # Legend for markers
    axes[0, 5].plot([], [], "^", color="red", markersize=6, label="$\\lambda_2$")
    axes[0, 5].plot([], [], "v", color="gray", markersize=6, label="$\\lambda_{\\max}$")
    axes[0, 5].legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    save_fig(fig, output_dir / f"fig_G1_eigenvalue_spectrum_{patient}_MSC_vs_ImCoh")


def main():
    parser = argparse.ArgumentParser(description="Section 2 eigenvalue spectrum figure.")
    parser.add_argument("--output-dir", type=Path, default=SECTION2_ROOT)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    apply_pub_style()

    for patient in REPR_PATIENTS:
        print(f"--- Figure G: Eigenvalue spectrum ({patient}) ---")
        fig_g_eigenvalue_spectrum(patient, args.output_dir / "fig_G", args.verbose)
    print("\nDone.")


if __name__ == "__main__":
    main()
