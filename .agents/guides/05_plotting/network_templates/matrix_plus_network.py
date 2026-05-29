"""matrix_plus_network — FC matrix imshow + cmap-linked network drawing.

Crystallised from the MSC-era figure
``data/outputs/figures/presentation_figures/Pat_02/fig1_msc_alpha.pdf``
— matrix on the left with an attached network on the right whose
edges are coloured by the SAME imshow colormap.  The reader's eye
reads both panels in the same colour language: dark cells map to dark
edges, bright cells to bright edges.

Three independent layers (mirrors the network_templates family):

  1.  **Matrix**   — :func:`lrg_eegfc.visuals.fc_templates.plot_fc_adjacency`
      with the canonical colorbar label, tick mode (generic / index /
      chnames), optional ``LogNorm``.
  2.  **Layout**   — any name in :data:`LAYOUT_REGISTRY` (``spring``,
      ``kk``, ``spectral``, ``laplacian_pca``, ``lrg_kk``,
      ``lrg_sfdp``, ``arf``, …).  ``spring`` uses the canonical
      ``k_base/√N`` auto-tuning by default; pass ``--k`` to override.
  3.  **Edge style** — :func:`draw_gamma_edges` with γ-power on
      width/alpha.  Default ``coloring="cmap"`` paints each edge as
      ``cmap(norm(|w|))``, locked to the matrix.  Other modes
      (``probe``, ``shaft``, ``signed``, ``weight``) are available for
      diagnostic comparison.

Phases are an ``nargs='+'`` list: pass one for a single panel
(``--phases rest_pre``), pass multiple for the row-per-phase mosaic
that matches the reference figure (``--phases rest_pre rest_post``).

Companion ``matrix_plus_network.md`` carries the style sheet.

Usage
-----
    conda activate lapbrain

    # Single panel — rsPre, default magma cmap, spring k=k_base/√N
    python .agents/guides/05_plotting/network_templates/matrix_plus_network.py \\
        --patient Pat_02 --band alpha --phases rest_pre

    # Reference figure layout — two-row rsPre / rsPost, shared scale
    python .agents/guides/05_plotting/network_templates/matrix_plus_network.py \\
        --patient Pat_02 --band alpha --phases rest_pre rest_post

    # Switch layout to Kamada-Kawai, keep cmap-linked edges
    python .agents/guides/05_plotting/network_templates/matrix_plus_network.py \\
        --patient Pat_02 --band alpha --phases rest_pre --layout kk

    # Probe-coloured edges instead of cmap (diagnostic)
    python .agents/guides/05_plotting/network_templates/matrix_plus_network.py \\
        --patient Pat_02 --band alpha --phases rest_pre \\
        --coloring probe --node-color shaft

    # LogNorm cmap so weak edges read above the floor
    python .agents/guides/05_plotting/network_templates/matrix_plus_network.py \\
        --patient Pat_02 --band alpha --phases rest_pre --log-scale

    # Spring k sweep — fixed k for diagnostic comparison
    python .agents/guides/05_plotting/network_templates/matrix_plus_network.py \\
        --patient Pat_02 --band alpha --phases rest_pre \\
        --spring-k-scale fixed --spring-k-fixed 0.05
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT, PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.visuals.network_templates import (
    DEFAULT_ALPHA_RANGE,
    DEFAULT_GAMMA,
    DEFAULT_NODE_SIZE,
    DEFAULT_WIDTH_RANGE,
    LAYOUT_REGISTRY,
    plot_fc_matrix_and_network_rows,
    spring_auto_k,
)
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.fc import load_fc_matrix


PHASE_TEX = {
    "rest_pre": r"rest$_{\mathrm{pre}}$",
    "task_learn": r"task$_{\mathrm{learn}}$",
    "task_test": r"task$_{\mathrm{test}}$",
    "rest_post": r"rest$_{\mathrm{post}}$",
}


def _load_channel_labels(patient: str) -> list[str]:
    csv = SEEG_DATAPATH / patient / "channel_labels.csv"
    return pd.read_csv(csv)["label"].astype(str).tolist()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--patient", default="Pat_02")
    parser.add_argument("--band", default="alpha")
    parser.add_argument(
        "--phases", nargs="+", default=["rest_pre", "rest_post"],
        choices=list(PHASE_LABELS),
        help="One or more phases.  N phases → N-row mosaic.",
    )
    parser.add_argument(
        "--fc-method", default="imcoh_abs",
        choices=["corr", "msc", "imcoh", "imcoh_abs", "imcoh_sq"],
    )
    # ----- matrix -----
    parser.add_argument(
        "--cmap", default=None,
        help="Colormap shared by matrix and network edges.  Default = "
             "magma for magnitude FC, RdBu_r for signed.",
    )
    parser.add_argument(
        "--vmin", type=float, default=0.0,
        help="Matrix ``imshow`` vmin and edge norm lower bound.  Default 0.",
    )
    parser.add_argument(
        "--vmax", type=float, default=None,
        help="Matrix ``imshow`` vmax.  Default = cohort-wide max across "
             "selected phases (shared scale).",
    )
    parser.add_argument(
        "--log-scale", action="store_true",
        help="LogNorm on the colour scale (matrix + network) — useful "
             "when most |w| cluster near zero.",
    )
    parser.add_argument(
        "--tick-labels", default="index",
        choices=["generic", "index", "chnames"],
        help="Axis-tick mode for the matrix panel.  ``index`` (default) "
             "matches the reference figure.",
    )
    parser.add_argument(
        "--no-shared-scale", dest="shared_scale", action="store_false",
        help="Give each phase its own vmax (no shared colorbar scale).",
    )
    parser.set_defaults(shared_scale=True)
    # ----- layout -----
    parser.add_argument(
        "--layout", default="spring",
        choices=list(LAYOUT_REGISTRY.keys()),
        help="Network layout.  Default: spring.",
    )
    parser.add_argument(
        "--spring-k-scale", default="msc_era",
        choices=["msc_era", "kbase_sqrtN", "inv_sqrtN", "fixed"],
        help="Spring auto-k rule (only used when layout=spring and "
             "--k is not set).  Default: msc_era (k=0.1) — the trick "
             "that produced the structured MSC presentation figures.  "
             "Pair with --iterations 50 (the default here).  Switch to "
             "kbase_sqrtN + --iterations 600 for a fully-converged "
             "layout when you actually want it.",
    )
    parser.add_argument(
        "--spring-k-base", type=float, default=5.0,
        help="``k_base`` in ``k = k_base / √N``.  Default 5.0.",
    )
    parser.add_argument(
        "--spring-k-fixed", type=float, default=0.1,
        help="``k`` when --spring-k-scale=fixed.  Default 0.1.",
    )
    parser.add_argument(
        "--k", type=float, default=None,
        help="Raw spring k, bypassing --spring-k-scale.  Default None.",
    )
    parser.add_argument(
        "--iterations", type=int, default=100,
        help="Spring solver iterations.  Default 100 — paired with "
             "msc_era k=0.1, the strong heavy-tail edges form clusters "
             "before relaxation equilibrates them.  Push to 600+ for a "
             "fully converged layout (and brace for a ball on "
             "heavy-tailed FC).",
    )
    parser.add_argument(
        "--no-shared-layout", dest="shared_layout", action="store_false",
        help="Compute layout per phase instead of once on the cross-phase "
             "average FC.  Default: shared layout — load-bearing for "
             "phase-to-phase comparison.",
    )
    parser.set_defaults(shared_layout=True)
    parser.add_argument(
        "--n-communities", type=int, default=10,
        help="LRG dendrogram cut for LRG-seeded layouts and "
             "node-color=community.  Default 10.",
    )
    # ----- edge style -----
    parser.add_argument(
        "--coloring", default="cmap",
        choices=["cmap", "probe", "shaft", "signed", "weight"],
        help="Edge coloring.  Default: cmap (locked to matrix).",
    )
    parser.add_argument(
        "--node-color", default="shaft",
        choices=["shaft", "community", "uniform"],
        help="Node coloring.  Default: shaft.",
    )
    parser.add_argument(
        "--gamma", type=float, default=2.0,
        help="γ in ``t = (|w|/|w|_max)**γ`` for width/alpha scaling.  "
             "Default 2 — matches the MSC-presentation recipe.  γ=6 "
             "(single_network's default) over-suppresses bulk edges "
             "into a dense floor on dense FC and looks like a ball.",
    )
    parser.add_argument(
        "--width-min", type=float, default=0.0,
        help="Default 0 — weak edges become invisible so only the "
             "strong-edge structural skeleton survives.",
    )
    parser.add_argument(
        "--width-max", type=float, default=4.0,
    )
    parser.add_argument(
        "--alpha-min", type=float, default=0.0,
        help="Default 0 (same reason as width-min).",
    )
    parser.add_argument(
        "--alpha-max", type=float, default=1.0,
    )
    parser.add_argument(
        "--min-alpha", type=float, default=0.0,
        help="Drop edges below this rendered alpha for file weight.  0 = no filter.",
    )
    parser.add_argument(
        "--node-size", type=int, default=DEFAULT_NODE_SIZE,
    )
    # ----- output -----
    parser.add_argument(
        "--out-dir",
        default=str(FIGURES_ROOT / "network_templates" / "matrix_plus_network"),
    )
    args = parser.parse_args()

    use_lrg_style()

    # Build spring layout_kwargs with the auto-tuned k unless --k is given.
    layout_kwargs = {}
    if args.layout == "spring":
        if args.k is not None:
            layout_kwargs["k"] = args.k
        else:
            # Probe N from the first phase to feed into the auto-k.
            A_probe = load_fc_matrix(
                patient=args.patient, phase=args.phases[0],
                band=args.band, fc_method=args.fc_method,
            )
            if A_probe is None:
                raise SystemExit(
                    f"No cached FC for {args.patient}/{args.band}/"
                    f"{args.phases[0]}/{args.fc_method}."
                )
            N = A_probe.shape[0]
            layout_kwargs["k"] = spring_auto_k(
                N, scale=args.spring_k_scale,
                k_base=args.spring_k_base,
                k_fixed=args.spring_k_fixed,
            )
        layout_kwargs["iterations"] = args.iterations

    fig, axes = plot_fc_matrix_and_network_rows(
        args.patient, args.band,
        phases=args.phases,
        fc_method=args.fc_method,
        cmap=args.cmap,
        vmin=args.vmin, vmax=args.vmax,
        log_scale=args.log_scale,
        tick_labels=args.tick_labels,
        layout=args.layout,
        layout_kwargs=layout_kwargs,
        shared_layout=args.shared_layout,
        coloring=args.coloring,
        node_color=args.node_color,
        node_size=args.node_size,
        n_communities=args.n_communities,
        gamma=args.gamma,
        width_range=(args.width_min, args.width_max),
        alpha_range=(args.alpha_min, args.alpha_max),
        min_alpha=args.min_alpha,
        shared_scale=args.shared_scale,
    )

    band_tex = BRAIN_BAND_TEX_DICT.get(args.band, args.band)
    for row, phase in enumerate(args.phases):
        phase_tex = PHASE_TEX.get(phase, phase)
        axes[row, 0].set_title(
            f"{args.patient}  ·  {band_tex}  ·  {phase_tex}  ·  "
            f"{args.fc_method}",
            fontsize=10, loc="left", pad=4,
        )
        axes[row, 1].set_title(
            f"network    layout: {args.layout}    edge: {args.coloring}",
            fontsize=9.5, loc="left", pad=4,
        )

    fig.tight_layout()

    k_tag = ""
    if args.layout == "spring":
        k_val = layout_kwargs["k"]
        k_tag = f"_spring-k{k_val:.4g}_it{args.iterations}"
    phases_tag = "-".join(p.replace("rest_", "r").replace("task_", "t")
                          for p in args.phases)
    out = (
        Path(args.out_dir)
        / f"{args.patient}_{args.band}_{phases_tag}_{args.fc_method}_"
          f"{args.layout}{k_tag}_{args.coloring}_g{args.gamma:g}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    print(f"size: {out.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()
