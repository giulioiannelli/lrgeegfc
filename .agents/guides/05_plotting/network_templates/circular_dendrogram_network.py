"""circular_dendrogram_network — chord network + LRG dendrogram overlay.

Template 1 in the "hierarchy-bundled chord" family.  Crystallised from
panel (c) of ``preprint_07_beta_rho_split_figure_test2.py`` —
specifically the "backbone" pass before any per-pair highlight is
applied: leaves equiangularly placed on the outer ring in
``leaves_list(Z)`` order, every cross-probe and same-shaft pair drawn
as a curvy graph-tool bezier bundled through the depth-2 LRG cut, and
the full circular dendrogram (every merge in ``Z``) overlaid as
matplotlib polylines on top.

Three independent layers (mirrors ``single_network``):

  1.  **Chord layout**  — graph-tool ``get_hierarchy_control_points``
      with a depth-2 tree at ``k_clusters`` (default 7).  Higher
      ``beta`` → tighter bundling through the cluster centroids.
  2.  **Edge style**    — same-probe pairs painted in the shaft's
      ``tab20`` colour; cross-probe in mid-gray.  Width and alpha scale
      as ``t = (|w|/|w|_max)**edge_gamma`` so strong edges read above
      the dendrogram overlay without overwhelming it.
  3.  **Dendrogram overlay** — full LRG dendrogram (every merge in
      ``lrg.linkage_matrix``) drawn as radial + arc polylines from the
      leaf ring inward to the root.  Toggleable via ``--no-dendrogram``.

Companion ``circular_dendrogram_network.md`` carries the style sheet.

Usage
-----
    conda activate lapbrain
    python .agents/guides/05_plotting/network_templates/circular_dendrogram_network.py \\
        --patient Pat_05 --band beta --phase rest_pre

    # No dendrogram overlay — chord-only baseline
    python .agents/guides/05_plotting/network_templates/circular_dendrogram_network.py \\
        --patient Pat_05 --band beta --phase rest_pre --no-dendrogram

    # Tighter / looser bundling
    python .agents/guides/05_plotting/network_templates/circular_dendrogram_network.py \\
        --patient Pat_05 --band beta --phase rest_pre --beta 0.75 --k-clusters 10
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT, PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_patient_metadata
from lrg_eegfc.visuals.network_templates import (
    HIERARCHY_DEFAULT_BETA,
    HIERARCHY_DEFAULT_FIT_VIEW,
    HIERARCHY_DEFAULT_K_CLUSTERS,
    HIERARCHY_DEFAULT_RENDER_PX,
    _load_inputs,
    plot_chord_with_dendrogram,
)
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result


PHASE_TEX = {
    "rest_pre": r"rest$_{\mathrm{pre}}$",
    "task_learn": r"task$_{\mathrm{learn}}$",
    "task_test": r"task$_{\mathrm{test}}$",
    "rest_post": r"rest$_{\mathrm{post}}$",
}


def _parse_contact_label(raw: str) -> str:
    """``F 1,G2`` → ``F1``; ``G' 4,X1`` → ``G4``; mojibake-tolerant."""
    head = str(raw).split(",")[0]
    return "".join(c for c in head if c.isascii() and c.isalnum())


def _load_contact_labels(patient: str) -> list[str]:
    meta = load_patient_metadata(patient, SEEG_DATAPATH)
    if meta is None or "label" not in meta.columns:
        return []
    return [_parse_contact_label(s) for s in meta["label"].tolist()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--patient", default="Pat_05")
    parser.add_argument("--band", default="beta")
    parser.add_argument("--phase", default="rest_pre", choices=PHASE_LABELS)
    parser.add_argument(
        "--fc-method", default="imcoh_abs",
        choices=["imcoh_abs", "imcoh_sq", "msc"],
        help="Magnitude FC only (chord edge recipe assumes ``|w|``).",
    )
    parser.add_argument(
        "--k-clusters", type=int, default=HIERARCHY_DEFAULT_K_CLUSTERS,
        help="Depth-2 cut of the LRG dendrogram for chord bundling "
             "(default 7).  More clusters → thinner bundles.",
    )
    parser.add_argument(
        "--beta", type=float, default=HIERARCHY_DEFAULT_BETA,
        help="Bezier bundle tightness (0 chord, 1 fully bundled through "
             "cluster centres).  Default 0.92.",
    )
    parser.add_argument(
        "--render-px", type=int, default=HIERARCHY_DEFAULT_RENDER_PX,
        help="PNG raster the graph-tool chord renders into.  Default "
             "2200 — crisp embed in a PDF, ~250 kB.",
    )
    parser.add_argument(
        "--fit-view", type=float, default=HIERARCHY_DEFAULT_FIT_VIEW,
        help="graph-tool fit_view (data bbox → canvas fraction).  "
             "Default 0.92.",
    )
    parser.add_argument(
        "--edge-gamma", type=float, default=2.0,
        help="γ in ``t = (|w|/|w|_max)**γ`` for the edge width/alpha "
             "scaling.  Default 2.0 with zero-floor width/alpha — "
             "suppresses the bottom ~half of edges into invisibility, "
             "leaves the top ~30% as visible context against the "
             "dendrogram polyline overlay.  Push to 3+ to make the "
             "dendrogram fully dominant; drop to 1 if you want a "
             "denser FC backdrop.",
    )
    parser.add_argument(
        "--no-dendrogram", action="store_true",
        help="Suppress the circular-dendrogram overlay (chord only).",
    )
    parser.add_argument(
        "--no-labels", action="store_true",
        help="Suppress the radial contact-label ring.",
    )
    parser.add_argument(
        "--vertex-size", type=float, default=24.0,
        help="Leaf marker size in graph-tool device pixels.  Default 24.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(
            FIGURES_ROOT / "network_templates" / "circular_dendrogram_network"
        ),
    )
    args = parser.parse_args()

    use_lrg_style()

    A, probes = _load_inputs(
        args.patient, args.band, args.phase, args.fc_method,
    )
    lrg = load_lrg_result(args.patient, args.phase, args.band, args.fc_method)
    if lrg is None or lrg.linkage_matrix is None:
        raise SystemExit(
            f"No LRG cache for {args.patient}/{args.band}/{args.phase}/"
            f"{args.fc_method}; cannot build chord backbone."
        )

    labels = None if args.no_labels else _load_contact_labels(args.patient)

    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    plot_chord_with_dendrogram(
        ax, A, probes, lrg,
        k_clusters=args.k_clusters, beta=args.beta,
        render_px=args.render_px, fit_view=args.fit_view,
        show_dendrogram=not args.no_dendrogram,
        labels=labels,
        edge_gamma=args.edge_gamma,
        vertex_size=args.vertex_size,
    )

    band_tex = BRAIN_BAND_TEX_DICT.get(args.band, args.band)
    phase_tex = PHASE_TEX.get(args.phase, args.phase)
    overlay = "+ dendro" if not args.no_dendrogram else "no dendro"
    ax.set_title(
        f"{args.patient}  ·  {band_tex}  ·  {phase_tex}  ·  "
        f"{args.fc_method}    chord (K={args.k_clusters}, "
        rf"$\beta$={args.beta:g}, {overlay})",
        fontsize=10, loc="left", pad=6,
    )

    tag = "with_dendro" if not args.no_dendrogram else "no_dendro"
    out = (
        Path(args.out_dir)
        / f"{args.patient}_{args.band}_{args.phase}_{args.fc_method}_"
          f"K{args.k_clusters}_b{args.beta:g}_g{args.edge_gamma:g}_{tag}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    print(f"size: {out.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()
