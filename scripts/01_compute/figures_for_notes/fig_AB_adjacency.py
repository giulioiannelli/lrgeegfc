#!/usr/bin/env python3
"""Section 2 Figures A & B: adjacency heatmaps (MSC vs ImCoh).

Figure A — One figure PER fc_method: phases x bands grid. Per-band color normalization.
Figure B — All patients, one phase, one band, MSC vs ImCoh rows.
           Produced for multiple (band, phase) combinations.

Same-probe blocks outlined in yellow on every heatmap.

Run:
  python scripts/10_notes_imcoh/fig_AB_adjacency.py [-v]
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
import networkx as nx
from _shared import (
    ALL_PATIENTS, BANDS, BAND_TEX, CMAP_MATRIX,
    PHASE_LABELS, BRAIN_BAND_TEX_DICT,
    REPR_PATIENTS, REPR_BANDS, REPR_PHASES,
    load_channel_labels, load_fc_for_patient, probe_sort_indices,
    draw_probe_outlines, apply_pub_style, save_fig, SECTION2_ROOT,
    extract_probe_labels, draw_network_edges, spring_k_for,
    compute_network_layout, NODE_SIZE_DEFAULT,
    SECTION2_METHODS as METHODS,
    SECTION2_METHOD_LABELS as METHOD_LABELS
)
from lrg_eegfc.workflow.fc import load_fc_matrix



# ---------------------------------------------------------------------------
# Figure A — One patient, separate figure per fc_method, per-band normalization
# ---------------------------------------------------------------------------

def fig_a_adjacency_single_method(
    patient: str,
    fc_method: str,
    output_dir: Path,
    verbose: bool = False,
):
    """Produce one figure: phases x bands grid for a single FC method.
    Each column (band) has its own color normalization based on the max
    across phases for that band.
    """
    ch = load_channel_labels(patient)
    sort_idx = probe_sort_indices(ch)

    # Detect available phases
    phases = [p for p in PHASE_LABELS
              if load_fc_matrix(patient, p, "alpha", fc_method) is not None]
    if not phases:
        print(f"  SKIP {patient} {fc_method}: no data")
        return

    data = load_fc_for_patient(patient, fc_method, phases=phases)
    n_phases = len(phases)
    n_bands = len(BANDS)

    # Per-band vmax: max absolute value across all phases for that band
    band_vmax = {}
    for band in BANDS:
        mx = 0
        for phase in phases:
            mat = data[band].get(phase)
            if mat is not None:
                vals = np.abs(mat[np.triu_indices_from(mat, k=1)])
                if len(vals) > 0:
                    mx = max(mx, np.percentile(vals[vals > 0], 99) if (vals > 0).any() else 0)
        band_vmax[band] = mx * 1.05 if mx > 0 else 1.0

    fig, axes = plt.subplots(
        n_phases, n_bands, figsize=(3.2 * n_bands, 3.0 * n_phases),
        constrained_layout=True,
    )
    if n_phases == 1:
        axes = axes.reshape(1, -1)

    for j, band in enumerate(BANDS):
        for i, phase in enumerate(phases):
            ax = axes[i, j]
            mat = data[band].get(phase)
            if mat is None:
                ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                        transform=ax.transAxes, fontsize=14, color="0.5")
                ax.set_xticks([]); ax.set_yticks([])
                continue

            mat_sorted = np.abs(mat[np.ix_(sort_idx, sort_idx)])
            # PowerNorm(gamma=0.5) stretches the low end so that phases
            # with 10x-smaller dynamic range than rest_pre (Pat_05 high_gamma
            # task phases, e.g.) still show visible structure at the
            # same per-band vmax.
            from matplotlib.colors import PowerNorm
            norm = PowerNorm(gamma=0.5, vmin=0, vmax=band_vmax[band])
            im = ax.imshow(mat_sorted, cmap=CMAP_MATRIX, norm=norm,
                           aspect="equal", interpolation="none")
            draw_probe_outlines(ax, ch, sort_idx)
            ax.set_xticks([]); ax.set_yticks([])

            if j == 0:
                ax.set_ylabel(phase, fontsize=16, fontweight="bold")
            if i == 0:
                ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=18,
                             fontweight="bold")

    label = METHOD_LABELS[fc_method]
    save_fig(fig, output_dir / f"fig_A_{fc_method}_{patient}")
    if verbose:
        for b in BANDS:
            print(f"    {b} vmax={band_vmax[b]:.3f}")


# ---------------------------------------------------------------------------
# Figure A (network view) — same data as A, but as network drawings
# ---------------------------------------------------------------------------

def fig_a_network_grid(
    patient: str,
    fc_method: str,
    output_dir: Path = SECTION2_ROOT / "fig_A",
    verbose: bool = False,
):
    """Network-drawing counterpart to the adjacency grid.

    Same phases x bands layout, but each cell is a spring-layout network
    of the top `edge_pct` edges. Nodes colored by electrode shaft.
    """
    ch = load_channel_labels(patient)
    probes = extract_probe_labels(ch)
    N = len(ch)

    phases = [p for p in PHASE_LABELS
              if load_fc_matrix(patient, p, "alpha", fc_method) is not None]
    if not phases:
        return
    n_phases = len(phases)
    n_bands = len(BANDS)

    # Shaft colors
    unique_shafts = sorted(set(probes))
    shaft_cmap = plt.get_cmap("tab20", len(unique_shafts))
    node_colors = [shaft_cmap(unique_shafts.index(p)) for p in probes]

    fig, axes = plt.subplots(
        n_phases, n_bands, figsize=(3.2 * n_bands, 3.0 * n_phases),
        constrained_layout=True,
    )
    if n_phases == 1:
        axes = axes.reshape(1, -1)

    for j, band in enumerate(BANDS):
        for i, phase in enumerate(phases):
            ax = axes[i, j]
            mat = load_fc_matrix(patient, phase, band, fc_method)
            if mat is None:
                ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                        transform=ax.transAxes, fontsize=12, color="0.5")
                ax.axis("off")
                continue

            A = np.abs(mat.copy())
            np.fill_diagonal(A, 0)
            pos_arr = compute_network_layout(
                A, fc_method=fc_method,
                patient=patient, phase=phase, band=band,
            )
            draw_network_edges(ax, pos_arr, A, fc_method=fc_method,
                               probe_labels=probes, highlight_same_probe=True)
            ax.scatter(pos_arr[:, 0], pos_arr[:, 1],
                       c=node_colors, s=NODE_SIZE_DEFAULT, edgecolors="white",
                       linewidths=0.3, zorder=5)
            ax.axis("off")

            if j == 0:
                ax.set_ylabel(phase, fontsize=14, fontweight="bold")
                # ylabel doesn't show on axis off, add text instead
                ax.text(-0.08, 0.5, phase, transform=ax.transAxes,
                        fontsize=14, fontweight="bold", rotation=90,
                        ha="right", va="center")
            if i == 0:
                ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=16,
                             fontweight="bold")

    label = METHOD_LABELS[fc_method]
    save_fig(fig, output_dir / f"fig_A_network_{fc_method}_{patient}")


# ---------------------------------------------------------------------------
# Figure B — All patients, one (band, phase), MSC vs ImCoh rows
# ---------------------------------------------------------------------------

def fig_b_all_patients(
    patients: list[str] | None = None,
    phase: str = "rest_pre",
    band: str = "beta",
    output_dir: Path = SECTION2_ROOT / "fig_B",
    verbose: bool = False,
):
    """2-row (MSC top, ImCoh bottom) x N-patient columns.
    Per-method normalization on the max across patients.
    """
    if patients is None:
        patients = ALL_PATIENTS

    n_pats = len(patients)

    fc_data = {}
    ch_data = {}
    sizes = []
    method_max = {m: 0.0 for m in METHODS}
    for pat in patients:
        ch = load_channel_labels(pat)
        ch_data[pat] = ch
        sizes.append(len(ch))
        for method in METHODS:
            mat = load_fc_matrix(pat, phase, band, method)
            fc_data[(method, pat)] = mat
            if mat is not None:
                vals = np.abs(mat[np.triu_indices_from(mat, k=1)])
                p99 = np.percentile(vals[vals > 0], 99) if (vals > 0).any() else 0
                method_max[method] = max(method_max[method], p99)

    vmax = {m: max(v * 1.05, 0.01) for m, v in method_max.items()}

    width_ratios = sizes
    fig, axes = plt.subplots(
        2, n_pats,
        figsize=(sum(sizes) * 0.035 + 2, max(sizes) * 0.035 * 2 + 1),
        gridspec_kw={"width_ratios": width_ratios, "height_ratios": [1, 1]},
        constrained_layout=True,
    )
    if n_pats == 1:
        axes = axes.reshape(2, 1)

    for row, method in enumerate(METHODS):
        vm = vmax[method]
        for col, pat in enumerate(patients):
            ax = axes[row, col]
            mat = fc_data[(method, pat)]
            ch = ch_data[pat]
            sort_idx = probe_sort_indices(ch)

            if mat is None:
                ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                        transform=ax.transAxes, fontsize=10, color="0.5")
                ax.set_xticks([]); ax.set_yticks([])
                continue

            mat_sorted = np.abs(mat[np.ix_(sort_idx, sort_idx)])
            im = ax.imshow(mat_sorted, cmap=CMAP_MATRIX, vmin=0, vmax=vm,
                           aspect="auto", interpolation="none")
            draw_probe_outlines(ax, ch, sort_idx)
            ax.set_xticks([]); ax.set_yticks([])

            if row == 0:
                ax.set_title(f"{pat}\nN={len(ch)}", fontsize=11, fontweight="bold")
            if col == 0:
                ax.set_ylabel(METHOD_LABELS[method], fontsize=14, fontweight="bold")

    band_tex = BRAIN_BAND_TEX_DICT[band]
    save_fig(fig, output_dir / f"fig_B_{band}_{phase}_MSC_vs_ImCoh")


# ---------------------------------------------------------------------------
# Main — produces multiple (patient, band, phase) combinations
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Section 2 adjacency heatmap figures.")
    parser.add_argument("--output-dir", type=Path, default=SECTION2_ROOT)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    apply_pub_style()

    # Figure A: separate MSC/ImCoh for representative patients (adjacency grid)
    for patient in REPR_PATIENTS:
        for method in METHODS:
            print(f"--- Figure A (adjacency): {method} {patient} ---")
            fig_a_adjacency_single_method(
                patient, method,
                args.output_dir / "fig_A", args.verbose,
            )

    # Figure A: same data but as network drawings
    for patient in REPR_PATIENTS:
        for method in METHODS:
            print(f"--- Figure A (network): {method} {patient} ---")
            fig_a_network_grid(
                patient, method,
                output_dir=args.output_dir / "fig_A", verbose=args.verbose,
            )

    # Figure B: all patients, multiple (band, phase)
    for band in REPR_BANDS:
        for phase in REPR_PHASES:
            print(f"\n--- Figure B: all patients, {band}, {phase} ---")
            fig_b_all_patients(
                phase=phase, band=band,
                output_dir=args.output_dir / "fig_B", verbose=args.verbose,
            )

    print("\nDone.")


if __name__ == "__main__":
    main()
