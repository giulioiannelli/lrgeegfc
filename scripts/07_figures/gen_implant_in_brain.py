"""Per-patient implant figure: shafts in a 2D glass brain + 3D HTML.

For each patient in the n=10 cohort, draws every sEEG shaft as a
coloured stick (markers at each contact + line connecting same-shaft
contacts in numeric order). Two outputs per patient:

- ``Pat_NN_implant_glass.pdf`` — nilearn 4-view glass brain, vector PDF.
- ``Pat_NN_implant_3d.html``  — nilearn ``view_connectome`` HTML, interactive.

No connectivity data is loaded — purely electrode geometry. Intended
as the foundation for later figures that overlay LRG measures
(dendrogram cuts, ρ̂, eigenmodes, ...) on the same brain template.

Usage::

    conda activate lapbrain
    python scripts/07_figures/gen_implant_in_brain.py
    # or one patient at a time:
    python scripts/07_figures/gen_implant_in_brain.py --patient Pat_05
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors as mcolors
from nilearn.plotting import plot_glass_brain, view_connectome

from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import parse_seeg_label
from lrg_eegfc.visuals.spatial import (
    load_spatial_metadata,
    prepare_spatial_coordinates,
)
from lrg_eegfc.visuals.styles import use_lrg_style

COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]


def _shaft_order(labels):
    """Return (shaft_id, sort_key) for each label so contacts plot along the shaft."""
    out = []
    for lab in labels:
        probe, contact = parse_seeg_label(str(lab))
        if probe is None:
            # Unparseable label — treat as its own singleton shaft.
            out.append((str(lab), 0))
        else:
            out.append((probe, int(contact) if contact is not None else 0))
    return out


def _shaft_color_map(shafts):
    """Map shaft id → hex colour using tab20 (enough for ≤14 shafts/patient)."""
    cmap = plt.get_cmap("tab20")
    uniq = sorted(set(shafts))
    return {s: mcolors.to_hex(cmap(i % cmap.N)) for i, s in enumerate(uniq)}


def _filter_valid(metadata):
    """Drop channels with missing coordinates; return cleaned metadata."""
    mask = metadata[["x", "y", "z"]].notna().all(axis=1)
    return metadata.loc[mask].reset_index(drop=True)


def _patient_geometry(patient):
    metadata = _filter_valid(load_spatial_metadata(patient, SEEG_DATAPATH))
    coords = prepare_spatial_coordinates(
        metadata, scale="mm", center=False, to_mni=True
    )
    labels = list(metadata["label"])
    shafts = [s for s, _ in _shaft_order(labels)]
    return metadata, coords, labels, np.asarray(shafts)


def _build_shaft_adjacency(coords, shafts, labels):
    """Binary adjacency: 1 between consecutive contacts on the same shaft."""
    n = len(labels)
    A = np.zeros((n, n), dtype=float)
    by_shaft = {}
    for i, (lab, sh) in enumerate(zip(labels, shafts)):
        _, contact = parse_seeg_label(str(lab))
        by_shaft.setdefault(sh, []).append((contact if contact is not None else 0, i))
    for nodes in by_shaft.values():
        nodes.sort()
        for (_, a), (_, b) in zip(nodes[:-1], nodes[1:]):
            A[a, b] = A[b, a] = 1.0
    return A


def _plot_glass(patient, coords, shafts, labels, out_pdf):
    color_map = _shaft_color_map(shafts)
    node_colors = [color_map[s] for s in shafts]

    fig = plt.figure(figsize=(11, 3.2))
    display = plot_glass_brain(
        None,
        display_mode="lyrz",
        figure=fig,
        annotate=True,
        alpha=0.25,
    )

    # Shaft lines first (so contacts sit on top).
    by_shaft = {}
    for i, (lab, sh) in enumerate(zip(labels, shafts)):
        _, contact = parse_seeg_label(str(lab))
        by_shaft.setdefault(sh, []).append((contact if contact is not None else 0, i))
    for sh, nodes in by_shaft.items():
        if len(nodes) < 2:
            continue
        nodes.sort()
        idx = [i for _, i in nodes]
        # Build a sub-adjacency among just these contacts (chain).
        sub_A = np.zeros((len(idx), len(idx)))
        for k in range(len(idx) - 1):
            sub_A[k, k + 1] = sub_A[k + 1, k] = 1.0
        try:
            display.add_graph(
                sub_A,
                coords[idx],
                node_color=color_map[sh],
                node_size=0,           # markers drawn separately below
                edge_kwargs={"color": color_map[sh], "linewidth": 1.6, "alpha": 0.9},
            )
        except TypeError:
            # Older nilearn — fall back to add_markers + manual lines via add_edges
            pass

    # Markers per shaft (consistent colour with line).
    for sh in sorted(set(shafts)):
        m = shafts == sh
        display.add_markers(
            coords[m],
            marker_color=color_map[sh],
            marker_size=20,
            edgecolors="black",
            linewidths=0.4,
        )

    # Small legend below the glass-brain row.
    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="-", color=color_map[s],
                   markerfacecolor=color_map[s], markeredgecolor="black",
                   markersize=5, label=s)
        for s in sorted(set(shafts))
    ]
    n = len(handles)
    ncol = min(n, 8)
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.10),
        ncol=ncol,
        frameon=False,
        fontsize=7,
        title=f"{patient} shafts (n={n}, contacts={len(shafts)})",
        title_fontsize=8,
    )

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    return out_pdf


def _plot_3d_html(patient, coords, shafts, labels, out_html):
    color_map = _shaft_color_map(shafts)
    node_colors = [color_map[s] for s in shafts]
    A = _build_shaft_adjacency(coords, shafts, labels)
    title = f"{patient} — implant ({len(shafts)} contacts, {len(set(shafts))} shafts)"
    view = view_connectome(
        adjacency_matrix=A,
        node_coords=coords,
        edge_threshold=None,
        edge_cmap="Greys",
        symmetric_cmap=False,
        linewidth=4.0,
        node_color=node_colors,
        node_size=4.0,
        colorbar=False,
        title=title,
    )
    out_html.parent.mkdir(parents=True, exist_ok=True)
    view.save_as_html(str(out_html))
    return out_html


def _plot_compound(patients_data, out_pdf, *, display_mode="z", ncols=5):
    """Single PDF with one glass-brain panel per patient.

    ``display_mode`` is a *single-view* nilearn mode (``"z"`` axial,
    ``"l"``/``"r"`` sagittal, ``"y"`` coronal). Shaft colouring is
    re-mapped per patient so each panel is internally consistent.
    """
    # Per-view aspect: sagittal brain is taller than wide; axial/coronal
    # are wider than tall — pick panel size so the brain roughly fills it.
    panel_wh = {
        "z": (2.6, 2.2),  # axial — wider
        "y": (2.6, 2.2),  # coronal — wider
        "l": (2.4, 2.6),  # sagittal — taller
        "r": (2.4, 2.6),
        "x": (2.4, 2.6),
    }.get(display_mode, (2.6, 2.4))
    pw, ph = panel_wh

    n = len(patients_data)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(pw * ncols, ph * nrows))
    axes = np.atleast_1d(axes).flatten()
    fig.subplots_adjust(wspace=0.05, hspace=0.18)

    for ax, (patient, coords, labels, shafts) in zip(axes, patients_data):
        color_map = _shaft_color_map(shafts)
        display = plot_glass_brain(
            None,
            display_mode=display_mode,
            axes=ax,
            annotate=False,
            alpha=0.25,
        )

        by_shaft = {}
        for i, (lab, sh) in enumerate(zip(labels, shafts)):
            _, contact = parse_seeg_label(str(lab))
            by_shaft.setdefault(sh, []).append((contact if contact is not None else 0, i))
        for sh, nodes in by_shaft.items():
            if len(nodes) < 2:
                continue
            nodes.sort()
            idx = [i for _, i in nodes]
            sub_A = np.zeros((len(idx), len(idx)))
            for k in range(len(idx) - 1):
                sub_A[k, k + 1] = sub_A[k + 1, k] = 1.0
            display.add_graph(
                sub_A,
                coords[idx],
                node_color=color_map[sh],
                node_size=0,
                edge_kwargs={"color": color_map[sh], "linewidth": 0.9, "alpha": 0.9},
            )
        for sh in sorted(set(shafts)):
            m = shafts == sh
            display.add_markers(
                coords[m],
                marker_color=color_map[sh],
                marker_size=6,
                edgecolors="black",
                linewidths=0.25,
            )

        ax.set_title(
            f"{patient}\n{len(set(shafts))} shafts · {len(shafts)} contacts",
            fontsize=8,
            pad=2,
        )

    for ax in axes[n:]:
        ax.axis("off")

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    return out_pdf


def _plot_cohort_overlay_2d(patients_data, out_pdf):
    """All patients in one 4-view glass brain, contacts coloured by patient."""
    cmap = plt.get_cmap("tab10")
    pat_color = {p: mcolors.to_hex(cmap(i % cmap.N))
                 for i, (p, *_) in enumerate(patients_data)}

    fig = plt.figure(figsize=(11.5, 3.4))
    display = plot_glass_brain(
        None,
        display_mode="lyrz",
        figure=fig,
        annotate=True,
        alpha=0.22,
    )

    for patient, coords, labels, shafts in patients_data:
        c = pat_color[patient]
        # Shaft chains in patient colour.
        by_shaft = {}
        for k, (lab, sh) in enumerate(zip(labels, shafts)):
            _, contact = parse_seeg_label(str(lab))
            by_shaft.setdefault(sh, []).append((contact if contact is not None else 0, k))
        for sh, nodes in by_shaft.items():
            if len(nodes) < 2:
                continue
            nodes.sort()
            idx = [j for _, j in nodes]
            sub_A = np.zeros((len(idx), len(idx)))
            for k in range(len(idx) - 1):
                sub_A[k, k + 1] = sub_A[k + 1, k] = 1.0
            display.add_graph(
                sub_A,
                coords[idx],
                node_color=c,
                node_size=0,
                edge_kwargs={"color": c, "linewidth": 0.7, "alpha": 0.55},
            )
        # Markers in patient colour.
        display.add_markers(
            coords,
            marker_color=c,
            marker_size=6,
            edgecolors="black",
            linewidths=0.15,
        )

    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="-",
                   color=pat_color[p], markerfacecolor=pat_color[p],
                   markeredgecolor="black", markersize=5,
                   label=f"{p} ({len(c)})")
        for p, c, *_ in patients_data
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=min(len(handles), 5),
        frameon=False,
        fontsize=7,
        title=f"Cohort n={len(patients_data)} — contacts per patient",
        title_fontsize=8,
    )
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    return out_pdf


def _plot_cohort_overlay_3d(patients_data, out_html):
    """All patients in one interactive 3D HTML, contacts coloured by patient."""
    cmap = plt.get_cmap("tab10")
    coord_blocks = []
    color_blocks = []
    shaft_pairs = []  # (global_i, global_j) for each shaft edge
    offset = 0
    for i, (patient, coords, labels, shafts) in enumerate(patients_data):
        c = mcolors.to_hex(cmap(i % cmap.N))
        coord_blocks.append(coords)
        color_blocks.extend([c] * coords.shape[0])
        by_shaft = {}
        for k, (lab, sh) in enumerate(zip(labels, shafts)):
            _, contact = parse_seeg_label(str(lab))
            by_shaft.setdefault(sh, []).append((contact if contact is not None else 0, k))
        for sh, nodes in by_shaft.items():
            if len(nodes) < 2:
                continue
            nodes.sort()
            idx = [j for _, j in nodes]
            for a, b in zip(idx[:-1], idx[1:]):
                shaft_pairs.append((offset + a, offset + b))
        offset += coords.shape[0]

    coords_all = np.vstack(coord_blocks)
    n = coords_all.shape[0]
    A = np.zeros((n, n), dtype=float)
    for a, b in shaft_pairs:
        A[a, b] = A[b, a] = 1.0

    title = f"Cohort n={len(patients_data)} implants — coloured by patient"
    view = view_connectome(
        adjacency_matrix=A,
        node_coords=coords_all,
        edge_threshold=None,
        edge_cmap="Greys",
        symmetric_cmap=False,
        linewidth=2.5,
        node_color=color_blocks,
        node_size=3.5,
        colorbar=False,
        title=title,
    )
    out_html.parent.mkdir(parents=True, exist_ok=True)
    view.save_as_html(str(out_html))
    return out_html


def run(patients, *, overlay: bool = True, mosaic: bool = False, per_patient: bool = True):
    use_lrg_style()
    out_root = FIGURES_ROOT / "implant_in_brain"
    out_root.mkdir(parents=True, exist_ok=True)

    cache = []
    for p in patients:
        metadata, coords, labels, shafts = _patient_geometry(p)
        cache.append((p, coords, labels, shafts))
        if per_patient:
            pdf = out_root / f"{p}_implant_glass.pdf"
            html = out_root / f"{p}_implant_3d.html"
            _plot_glass(p, coords, shafts, labels, pdf)
            _plot_3d_html(p, coords, shafts, labels, html)
            print(f"[{p}] {len(shafts)} contacts, {len(set(shafts))} shafts → {pdf.name}, {html.name}")

    if overlay:
        pdf = out_root / "cohort_implants_overlay.pdf"
        html = out_root / "cohort_implants_overlay.html"
        _plot_cohort_overlay_2d(cache, pdf)
        _plot_cohort_overlay_3d(cache, html)
        print(f"[overlay] {len(cache)} patients → {pdf.name}, {html.name}")

    if mosaic:
        for mode, tag in [("z", "axial"), ("l", "sagittal_l"), ("y", "coronal")]:
            out = out_root / f"cohort_implants_mosaic_{tag}.pdf"
            _plot_compound(cache, out, display_mode=mode)
            print(f"[mosaic:{tag}] {len(cache)} patients → {out.name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", action="append",
                    help="Patient id (repeatable). Default: all 10 cohort patients.")
    ap.add_argument("--no-per-patient", action="store_true",
                    help="Skip per-patient figures.")
    ap.add_argument("--no-overlay", action="store_true",
                    help="Skip the single cohort-overlay figure.")
    ap.add_argument("--mosaic", action="store_true",
                    help="Also emit per-patient mosaic PDFs (one panel per patient).")
    args = ap.parse_args()
    patients = args.patient or COHORT
    run(patients,
        per_patient=not args.no_per_patient,
        overlay=not args.no_overlay,
        mosaic=args.mosaic)


if __name__ == "__main__":
    main()
