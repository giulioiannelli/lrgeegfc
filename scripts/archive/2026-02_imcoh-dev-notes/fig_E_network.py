#!/usr/bin/env python3
"""Section 2 Figure E: brain connectome and spring layout.

E1 — Brain connectome (nilearn) for patients with valid coords, MSC vs ImCoh.
E2 — Spring layout (patients x methods), nodes colored by shaft.

Uses FC-method-aware spring constants and edge scaling from config.const.

Run:
  python scripts/10_notes_imcoh/fig_E_network.py [-v]
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

from lrg_eegfc.utils.metrics.hypothesis import (
    BRAIN_BAND_TEX_DICT,
    REPR_BANDS, REPR_PHASES, PATIENTS_WITH_COORDS,
    load_channel_labels, extract_probe_labels,
    apply_pub_style, save_fig, SECTION2_ROOT,
    draw_network_edges, scale_edge_weights,
    spring_k_for, compute_network_layout, NODE_SIZE_DEFAULT,
    render_sbm_panel, render_lrg_panel, compute_percolation_threshold,
    SECTION2_METHODS as METHODS,
    SECTION2_METHOD_LABELS as METHOD_LABELS
)
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.visuals.spatial import load_spatial_metadata, prepare_spatial_coordinates
from lrg_eegfc.config.paths import SEEG_DATAPATH



def _get_shaft_colors(probes: list[str]) -> tuple[list, dict]:
    """Assign a distinct color to each shaft."""
    unique = sorted(set(probes))
    cmap = plt.get_cmap("tab20", len(unique))
    color_map = {s: cmap(i) for i, s in enumerate(unique)}
    return [color_map[p] for p in probes], color_map


# ---------------------------------------------------------------------------
# Figure E1 — Brain connectome (nilearn glass brain)
# ---------------------------------------------------------------------------

def fig_e1_brain_connectome(
    patients: list[str] | None = None,
    band: str = "beta",
    phase: str = "rest_pre",
    output_dir: Path = SECTION2_ROOT / "fig_E",
    verbose: bool = False,
    *,
    gamma: float | None = 3.0,
    wmin: float | None = 0.1,
    wmax: float | None = 1.0,
    amin: float | None = 0.01,
    amax: float | None = 0.7,
    scaling: str = "value",
    display_mode: str = "z",
):
    """Nilearn glass brain + rank-based edge overlay (same recipe as fig_E2).

    Nilearn is used ONLY for the brain outline + node markers (we pass an
    all-zero adjacency so it draws no edges). Edges are overlaid manually
    onto the projected axes using :func:`draw_network_edges`, which honours
    the γ / wmin / wmax / αmin / αmax knobs (defaults from config.const).
    """
    if patients is None:
        patients = PATIENTS_WITH_COORDS

    try:
        from nilearn.plotting import plot_connectome as nilearn_plot_connectome
    except ImportError:
        print("  SKIP E1: nilearn not installed")
        return

    # 3D → 2D projection per display_mode (for overlay coords)
    proj_ij = {"x": (1, 2), "y": (0, 2), "z": (0, 1)}[display_mode]

    n_pats = len(patients)
    fig, axes = plt.subplots(2, n_pats, figsize=(3.0 * n_pats, 5.2))
    fig.subplots_adjust(left=0.03, right=0.99, top=0.94, bottom=0.02,
                        wspace=0.0, hspace=0.0)
    if n_pats == 1:
        axes = axes.reshape(2, 1)


    for col, pat in enumerate(patients):
        metadata = load_spatial_metadata(pat, SEEG_DATAPATH)
        if metadata is None:
            for r in range(2):
                axes[r, col].axis("off")
            continue
        try:
            coords = prepare_spatial_coordinates(metadata, scale="mm",
                                                 center=False, to_mni=True)
        except ValueError:
            for r in range(2):
                axes[r, col].axis("off")
            continue

        ch = load_channel_labels(pat)
        probes = extract_probe_labels(ch)
        node_colors, _ = _get_shaft_colors(probes)
        N = len(ch)
        coords_arr = np.asarray(coords[:N] if len(coords) >= N else coords)
        M = len(coords_arr)
        pos2d = coords_arr[:, list(proj_ij)]

        for row, method in enumerate(METHODS):
            ax = axes[row, col]
            mat = load_fc_matrix(pat, phase, band, method)
            if mat is None:
                ax.axis("off")
                continue

            A = np.abs(mat.copy())
            np.fill_diagonal(A, 0)
            A = A[:M, :M]

            # Glass brain + nodes only (zero adjacency → nilearn draws no edges)
            disp = nilearn_plot_connectome(
                np.zeros_like(A), coords_arr,
                node_color=node_colors[:M],
                node_size=30,
                display_mode=display_mode,
                axes=ax,
                colorbar=False,
                annotate=False,
            )

            # Overlay our edges on the projected matplotlib axes
            overlay_ax = list(disp.axes.values())[0].ax
            draw_network_edges(
                overlay_ax, pos2d, A, fc_method=method,
                probe_labels=probes[:M], highlight_same_probe=True,
                gamma=gamma, wmin=wmin, wmax=wmax, amin=amin, amax=amax,
                scaling=scaling,
            )

            # Zoom to electrode bounding box (+ margin)
            xmin, ymin = pos2d.min(axis=0)
            xmax, ymax = pos2d.max(axis=0)
            xspan = max(xmax - xmin, 1.0)
            yspan = max(ymax - ymin, 1.0)
            pad = 0.10 * max(xspan, yspan)
            overlay_ax.set_xlim(xmin - pad, xmax + pad)
            overlay_ax.set_ylim(ymin - pad, ymax + pad)

            if row == 0:
                ax.set_title(pat, fontsize=12, fontweight="bold")
            if col == 0:
                ax.text(-0.02, 0.5, METHOD_LABELS[method],
                        transform=ax.transAxes, rotation=90,
                        ha="right", va="center",
                        fontsize=12, fontweight="bold")

    save_fig(fig, output_dir / f"fig_E1_brain_connectome_{band}_{phase}")


# ---------------------------------------------------------------------------
# Figure E2 — Spring network layout
# ---------------------------------------------------------------------------

def fig_e2_sbm_grid(
    patients: list[str],
    band: str,
    phase: str,
    output_dir: Path,
    verbose: bool = False,
    *,
    threshold_scale: float = 1.0,
    panel_size: tuple = (700, 700),
) -> None:
    """Render fig_E2 via graph-tool's native state.draw() — one PDF per
    (patient, method) panel, then compose into a grid with pdfjam.

    Memory-guarded: each panel is rendered sequentially with explicit gc
    between iterations; RSS stays below ~1 GB.
    """
    import subprocess
    import tempfile
    import resource
    import gc

    # Memory cap: 6 GB virtual address space
    try:
        resource.setrlimit(resource.RLIMIT_AS, (6 * 1024**3, 6 * 1024**3))
    except (ValueError, OSError):
        pass

    output_dir.mkdir(parents=True, exist_ok=True)
    n_pats = len(patients)

    # METHODS defines row order: row 0 = MSC, row 1 = |ImCoh|.
    # pdfjam --nup COLS x ROWS reads input files in row-major order
    # (left-to-right, top-to-bottom). So we iterate methods (rows) as
    # the outer loop and patients (cols) as the inner loop.
    tmpdir = Path(tempfile.mkdtemp(prefix="sbm_"))
    panel_files = []
    try:
        for row, method in enumerate(METHODS):
            for col, pat in enumerate(patients):
                mat = load_fc_matrix(pat, phase, band, method)
                if mat is None:
                    if verbose:
                        print(f"  SKIP {pat} {method}: no data")
                    continue
                A = np.abs(mat.copy())
                np.fill_diagonal(A, 0)

                ch = load_channel_labels(pat)
                probes = extract_probe_labels(ch)
                node_colors, _ = _get_shaft_colors(probes)

                panel_path = tmpdir / f"p_{row}_{col:02d}_{pat}_{method}.pdf"
                if verbose:
                    theta = (compute_percolation_threshold(A)
                             * threshold_scale)
                    n_draw = int((A >= theta).sum() // 2)
                    print(f"    [{row},{col}] {pat} {method}: "
                          f"θ={theta:.4f}, edges={n_draw}", flush=True)

                render_sbm_panel(
                    panel_path, A,
                    probe_labels=probes,
                    shaft_colors=node_colors,
                    threshold_scale=threshold_scale,
                    output_size=panel_size,
                )
                panel_files.append(str(panel_path))
                del mat, A, ch, probes, node_colors
                gc.collect()

        if not panel_files:
            print("  no panels produced; skipping compose step")
            return

        combined = output_dir / (
            f"fig_E2_sbm_{'_'.join(patients)}_{band}_{phase}.pdf"
        )
        cmd = [
            "pdfjam", "--nup", f"{n_pats}x{len(METHODS)}",
            "--papersize", f"{{{panel_size[0] * n_pats}pt,"
                           f"{panel_size[1] * len(METHODS)}pt}}",
            "--outfile", str(combined),
        ] + panel_files
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  pdfjam failed: {result.stderr}", flush=True)
            # Fallback: use pdfunite (just stacks pages vertically)
            subprocess.run(
                ["pdfunite", *panel_files, str(combined)],
                check=False,
            )
        print(f"  Saved: {combined}")
    finally:
        # Clean up temp panels
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)


def fig_e2_lrg_grid(
    patients: list[str],
    band: str,
    phase: str,
    output_dir: Path,
    verbose: bool = False,
    *,
    threshold_scale: float = 1.0,
    panel_size: tuple = (700, 700),
    beta: float = 0.8,
    k_clusters: int = 30,
) -> None:
    """Render fig_E2 using the LRG dendrogram as the hierarchy tree for
    graph-tool's bezier edges — no SBM fit.

    Nodes sit on the outer ring in dendrogram order (same-cluster nodes
    angularly adjacent); edges curve through the LRG hierarchy. Each
    (patient, method) panel is one PDF, composed with pdfjam.
    """
    import subprocess
    import tempfile
    import resource
    import gc

    try:
        resource.setrlimit(resource.RLIMIT_AS, (6 * 1024**3, 6 * 1024**3))
    except (ValueError, OSError):
        pass

    output_dir.mkdir(parents=True, exist_ok=True)
    n_pats = len(patients)

    tmpdir = Path(tempfile.mkdtemp(prefix="lrg_"))
    panel_files = []
    try:
        for row, method in enumerate(METHODS):
            for col, pat in enumerate(patients):
                mat = load_fc_matrix(pat, phase, band, method)
                if mat is None:
                    if verbose:
                        print(f"  SKIP {pat} {method}: no FC cache")
                    continue
                lrg = load_lrg_result(pat, phase, band, method)
                if lrg is None:
                    if verbose:
                        print(f"  SKIP {pat} {method}: no LRG cache")
                    continue

                A = np.abs(mat.copy())
                np.fill_diagonal(A, 0)
                ch = load_channel_labels(pat)
                probes = extract_probe_labels(ch)
                node_colors, _ = _get_shaft_colors(probes)

                panel_path = tmpdir / f"p_{row}_{col:02d}_{pat}_{method}.pdf"
                if verbose:
                    # Match renderer: threshold is cross-probe-only.
                    Nv = A.shape[0]
                    same_mask = np.array([
                        [probes[i] == probes[j] for j in range(Nv)]
                        for i in range(Nv)
                    ])
                    A_cross = np.where(same_mask, 0.0, A)
                    theta = (compute_percolation_threshold(A_cross)
                             * threshold_scale)
                    n_cross_kept = int(((A_cross >= theta)
                                        & (A_cross > 0)).sum() // 2)
                    n_same = int((np.triu(same_mask, k=1)
                                  & (A > 0)).sum())
                    print(f"    [{row},{col}] {pat} {method}: "
                          f"θ_cross={theta:.4f}, "
                          f"cross-probe kept={n_cross_kept}, "
                          f"same-probe={n_same}", flush=True)

                render_lrg_panel(
                    panel_path, A, lrg,
                    probe_labels=probes,
                    shaft_colors=node_colors,
                    threshold_scale=threshold_scale,
                    output_size=panel_size,
                    beta=beta,
                    k_clusters=k_clusters,
                )
                panel_files.append(str(panel_path))
                del mat, A, lrg, ch, probes, node_colors
                gc.collect()

        if not panel_files:
            print("  no panels produced; skipping compose step")
            return

        combined = output_dir / (
            f"fig_E2_lrg_{'_'.join(patients)}_{band}_{phase}.pdf"
        )
        cmd = [
            "pdfjam", "--nup", f"{n_pats}x{len(METHODS)}",
            "--papersize", f"{{{panel_size[0] * n_pats}pt,"
                           f"{panel_size[1] * len(METHODS)}pt}}",
            "--outfile", str(combined),
        ] + panel_files
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  pdfjam failed: {result.stderr}", flush=True)
            subprocess.run(
                ["pdfunite", *panel_files, str(combined)],
                check=False,
            )
        print(f"  Saved: {combined}")
    finally:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)


def fig_e2_spring_layout(
    patients: list[str] | None = None,
    band: str = "beta",
    phase: str = "rest_pre",
    edge_pct: float = 0.15,
    output_dir: Path = SECTION2_ROOT / "fig_E",
    verbose: bool = False,
    *,
    layout_override: str | None = None,
    spring_k: float | None = None,
    spring_iter: int | None = None,
    sbm_threshold_scale: float = 1.0,
    lrg_k_clusters: int = 30,
):
    """patients x methods grid, force-directed layout, shaft coloring.
    Uses FC-method-aware spring constants and edge scaling from config.

    Parameters
    ----------
    layout_override : str, optional
        Force the same layout for ALL methods (e.g. ``"spring"``,
        ``"mds_lrg_continuous"``).  ``None`` = per-method default from config.
    spring_k : float, optional
        Override the spring repulsion constant (only for spring layouts).
    spring_iter : int, optional
        Override spring iterations (only for spring layouts).
    """
    if patients is None:
        patients = ["Pat_02", "Pat_05", "Pat_08"]

    # Nested-SBM layout uses a completely different rendering path:
    # graph-tool's native state.draw() per panel + pdfjam grid compose.
    # Skip the matplotlib path entirely.
    if layout_override == "sbm_nested":
        fig_e2_sbm_grid(
            patients=patients, band=band, phase=phase,
            output_dir=output_dir, verbose=verbose,
            threshold_scale=sbm_threshold_scale,
        )
        return

    if layout_override == "lrg_dendrogram":
        fig_e2_lrg_grid(
            patients=patients, band=band, phase=phase,
            output_dir=output_dir, verbose=verbose,
            threshold_scale=sbm_threshold_scale,
            k_clusters=lrg_k_clusters,
        )
        return

    # Horizontal layout: 2 rows (MSC top, |ImCoh| bottom) x n_pats columns.
    n_pats = len(patients)
    fig, axes = plt.subplots(2, n_pats, figsize=(5.5 * n_pats, 11))
    if n_pats == 1:
        axes = axes.reshape(2, 1)

    shaft_cmap = None  # will be set from last patient for legend

    for col, pat in enumerate(patients):
        ch = load_channel_labels(pat)
        probes = extract_probe_labels(ch)
        node_colors, shaft_cmap = _get_shaft_colors(probes)
        N = len(ch)

        for row, method in enumerate(METHODS):
            ax = axes[row, col]
            mat = load_fc_matrix(pat, phase, band, method)
            if mat is None:
                ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                        transform=ax.transAxes)
                ax.axis("off")
                continue

            A = np.abs(mat.copy())
            np.fill_diagonal(A, 0)

            layout_kw = {}
            if layout_override is not None:
                layout_kw["layout_method"] = layout_override
            if spring_k is not None:
                layout_kw["spring_k"] = spring_k
            if spring_iter is not None:
                layout_kw["spring_iter"] = spring_iter
            pos_arr = compute_network_layout(
                A, fc_method=method,
                patient=pat, phase=phase, band=band,
                **layout_kw,
            )

            # Draw all edges with gamma-scaled width/alpha (weak edges fade)
            draw_network_edges(ax, pos_arr, A, fc_method=method,
                               probe_labels=probes, highlight_same_probe=True)

            # Draw nodes
            ax.scatter(pos_arr[:, 0], pos_arr[:, 1],
                       c=node_colors, s=NODE_SIZE_DEFAULT, edgecolors="white",
                       linewidths=0.5, zorder=5)

            ax.set_title(f"{pat} — {METHOD_LABELS[method]}",
                         fontsize=12, fontweight="bold")
            ax.axis("off")

            margin = 0.1
            xmin, xmax = pos_arr[:, 0].min(), pos_arr[:, 0].max()
            ymin, ymax = pos_arr[:, 1].min(), pos_arr[:, 1].max()
            span = max(xmax - xmin, ymax - ymin)
            ax.set_xlim(xmin - margin * span, xmax + margin * span)
            ax.set_ylim(ymin - margin * span, ymax + margin * span)

            if verbose:
                n_edges = (A > 0).sum() // 2
                print(f"    {pat} {method}: edges={n_edges}")

    # Legend
    from matplotlib.patches import Patch
    fig.legend(handles=[
        Patch(facecolor=(0.8, 0.2, 0.2, 0.5), label="same-probe edge"),
        Patch(facecolor=(0.5, 0.5, 0.5, 0.4), label="cross-probe edge"),
    ], loc="lower center", fontsize=9, ncol=2, bbox_to_anchor=(0.5, -0.01))

    band_tex = BRAIN_BAND_TEX_DICT[band]
    fig.tight_layout()
    save_fig(fig, output_dir / f"fig_E2_spring_{'_'.join(patients)}_{band}_{phase}")


# ---------------------------------------------------------------------------
# Main — variety over patients, bands, phases
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Section 2 network figures.")
    parser.add_argument("--output-dir", type=Path, default=SECTION2_ROOT)
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--e1-only", action="store_true",
                        help="Only regenerate fig_E1 (brain connectome).")
    parser.add_argument("--e2-only", action="store_true",
                        help="Only regenerate fig_E2 (network layout).")
    parser.add_argument("--bands", nargs="+", default=None,
                        help="Override bands (default REPR_BANDS).")
    parser.add_argument("--phases", nargs="+", default=None,
                        help="Override phases (default REPR_PHASES).")
    parser.add_argument("--gamma", type=float, default=3.0,
                        help="Rank exponent γ for edges (fig_E1 default 3.0).")
    parser.add_argument("--wmin", type=float, default=0.1)
    parser.add_argument("--wmax", type=float, default=1.0)
    parser.add_argument("--amin", type=float, default=0.01)
    parser.add_argument("--amax", type=float, default=0.7)
    parser.add_argument("--scaling", choices=["rank", "value"], default="value",
                        help="Edge scaling: 'rank' (per-patient rank) or "
                             "'value' (per-patient min-max normalised weight, "
                             "default — makes γ/w/α knobs comparable across "
                             "subjects).")
    parser.add_argument("--layout", default=None,
                        choices=["spring", "mds_lrg_continuous",
                                 "sbm_nested", "lrg_dendrogram",
                                 "mds_ultrametric", "kk_ultrametric"],
                        help="Force same layout for ALL methods in fig_E2 "
                             "(default: per-method from config).")
    parser.add_argument("--spring-k", type=float, default=None,
                        help="Override spring repulsion constant (spring layouts only).")
    parser.add_argument("--spring-iter", type=int, default=None,
                        help="Override spring iterations (spring layouts only).")
    parser.add_argument("--sbm-threshold-scale", type=float, default=1.0,
                        help="Scale the auto percolation threshold for "
                             "sbm_nested / lrg_dendrogram layouts: 1.0 = "
                             "just-connected, >1 sparser, <1 denser.")
    parser.add_argument("--lrg-k-clusters", type=int, default=30,
                        help="Number of communities to cut the LRG "
                             "dendrogram at (lrg_dendrogram layout only). "
                             "Low k -> coarser arcs, high k -> same-shaft "
                             "mixing becomes visible for |ImCoh|.")
    args = parser.parse_args()

    apply_pub_style()
    out = args.output_dir / "fig_E"
    bands = args.bands or REPR_BANDS
    phases = args.phases or REPR_PHASES

    # E1: brain connectome for patients with coords
    if not args.e2_only:
        for band in bands:
            for phase in phases:
                print(f"--- Figure E1: brain connectome ({band}, {phase}) ---")
                fig_e1_brain_connectome(
                    patients=PATIENTS_WITH_COORDS,
                    band=band, phase=phase, output_dir=out, verbose=args.verbose,
                    gamma=args.gamma, wmin=args.wmin, wmax=args.wmax,
                    amin=args.amin, amax=args.amax, scaling=args.scaling,
                )
    if args.e1_only:
        return

    # E2: spring layout for representative patients
    for band in bands:
        for phase in phases:
            print(f"\n--- Figure E2: spring layout ({band}, {phase}) ---")
            fig_e2_spring_layout(
                patients=["Pat_02", "Pat_05", "Pat_08"],
                band=band, phase=phase, output_dir=out, verbose=args.verbose,
                layout_override=args.layout,
                spring_k=args.spring_k,
                spring_iter=args.spring_iter,
                sbm_threshold_scale=args.sbm_threshold_scale,
                lrg_k_clusters=args.lrg_k_clusters,
            )

    print("\nDone.")


if __name__ == "__main__":
    main()
