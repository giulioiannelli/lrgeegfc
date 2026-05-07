"""single_network — canonical FC network drawing (one panel).

Thin wrapper around ``lrg_eegfc.visuals.plot_fc_network``.  Loads one
(patient, band, phase, fc_method) FC matrix from cache, computes the
chosen layout, renders edges with parameterised width and alpha
scaling functions, and saves a vector PDF.

Edge scaling — the formulas exposed for iteration
-------------------------------------------------

Let ``W = max_{ij} |w_{ij}|`` and ``t_{ij} = (|w_{ij}| / W)^γ``.

    width(w)  = w_min + (w_max − w_min) · t
    alpha(w)  = α_min + (α_max − α_min) · t

``γ = 1`` is the **linear** baseline (default).  ``γ = 6`` is the
power-law variant inherited from the archive notes; we'll evaluate
both side-by-side and decide.

The chosen ``(γ, w_min, w_max, α_min, α_max)`` are stamped into the
panel title and the file name so iterations are traceable at a glance.

Usage
-----
    conda activate lapbrain
    python .agents/guides/05_plotting/network_templates/single_network.py \
        --patient Pat_05 --band beta --phase rest_pre

See sibling ``single_network.md`` for the style sheet.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.visuals.layout import add_provenance_footer
from lrg_eegfc.visuals.network_templates import (
    LAYOUT_REGISTRY,
    is_lrg_layout,
    plot_fc_network,
    _load_inputs,
)
from lrg_eegfc.workflow.lrg import load_lrg_result


PHASE_TEX = {
    "rest_pre": r"rest$_{\mathrm{pre}}$",
    "task_learn": r"task$_{\mathrm{learn}}$",
    "task_test": r"task$_{\mathrm{test}}$",
    "rest_post": r"rest$_{\mathrm{post}}$",
}


def _scaling_title(
    *,
    patient: str,
    band: str,
    phase: str,
    fc_method: str,
    layout: str,
    coloring: str,
    gamma: float,
    width_range: tuple[float, float],
    alpha_range: tuple[float, float],
    k: float | None = None,
    iterations: int | None = None,
    t_mode: str = "magnitude",
) -> str:
    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    phase_tex = PHASE_TEX.get(phase, phase)
    layout_tag = layout
    if layout == "spring":
        bits = []
        if k is not None:
            bits.append(f"k={k:g}")
        if iterations is not None:
            bits.append(f"iter={iterations}")
        if bits:
            layout_tag = f"spring ({', '.join(bits)})"
    head = (
        f"{patient}  ·  {band_tex}  ·  {phase_tex}  ·  {fc_method}  ·  "
        f"layout: {layout_tag}  ·  edge-color: {coloring}"
    )
    w0, w1 = width_range
    a0, a1 = alpha_range
    if t_mode == "rank":
        t_pre = r"$t_0 = \mathrm{rank}(|w|)/E$"
    else:
        t_pre = r"$t_0 = |w|/|w|_{\max}$"
    if gamma == 1.0:
        kind = "linear"
        t_str = rf"{t_pre},  $t = t_0$"
    else:
        kind = f"γ = {gamma:g}"
        t_str = rf"{t_pre},  $t = t_0^{{{gamma:g}}}$"
    formula = (
        rf"width = {w0:g} + {w1 - w0:g}$\cdot t$,   "
        rf"$\alpha$ = {a0:g} + {a1 - a0:g}$\cdot t$,   "
        rf"{t_str}   ({kind})"
    )
    return f"{head}\n{formula}"


def _scaling_filename_tag(
    *,
    gamma: float,
    width_range: tuple[float, float],
    alpha_range: tuple[float, float],
    t_mode: str = "magnitude",
    min_alpha: float = 0.0,
) -> str:
    mode = "rank" if t_mode == "rank" else "mag"
    g = "lin" if gamma == 1.0 else f"g{gamma:g}"
    w0, w1 = width_range
    a0, a1 = alpha_range
    tag = f"{mode}_{g}_w{w0:g}-{w1:g}_a{a0:g}-{a1:g}"
    if min_alpha > 0.0:
        tag += f"_minA{min_alpha:g}"
    return tag


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--patient", default="Pat_05")
    parser.add_argument("--band", default="beta")
    parser.add_argument("--phase", default="rest_pre")
    parser.add_argument(
        "--fc-method", default="imcoh_abs",
        choices=["corr", "msc", "imcoh", "imcoh_abs", "imcoh_sq"],
    )
    parser.add_argument(
        "--layout", default="spring",
        choices=list(LAYOUT_REGISTRY.keys()),
        help="Layout algorithm.  Default: spring.",
    )
    parser.add_argument(
        "--node-color", default="shaft",
        choices=["shaft", "community", "uniform"],
        help="Node coloring.  Default: shaft.",
    )
    parser.add_argument(
        "--coloring", default=None,
        choices=["probe", "shaft", "signed", "weight"],
        help=(
            "Edge coloring.  'probe' = same-probe red / cross-probe gray "
            "(default for magnitude FC).  'shaft' = same-probe painted "
            "in the shaft's tab20 colour (matches node colour); "
            "cross-probe gray.  'signed' = RdBu_r diverging on signed "
            "weight (default for corr / raw imcoh).  'weight' = "
            "sequential gray scale by |w|."
        ),
    )
    # ----- spring layout parameters (raw, for iteration) -----
    parser.add_argument(
        "--k", type=float, default=None,
        help="Raw spring k (fed to nx.spring_layout verbatim).  Default "
             "None → k_base/√N (canonical).",
    )
    parser.add_argument(
        "--iterations", type=int, default=600,
        help="Spring solver iterations.  Default 600.",
    )
    parser.add_argument(
        "--min-alpha", type=float, default=0.0,
        help="Drop edges below this rendered alpha for file weight "
             "(layout still uses full graph).  Default 0 (no filter).",
    )
    # ----- edge scaling parameters (exposed for iteration) -----
    parser.add_argument(
        "--t-mode", default="magnitude", choices=["magnitude", "rank"],
        help="Pre-γ transform.  'magnitude' (default): t = |w|/|w|_max. "
             "'rank': t = rank(|w|)/E (uniform on (0,1] before γ; γ then "
             "acts on a well-conditioned distribution and concentrates "
             "ink on top edges without driving the bulk to zero).",
    )
    parser.add_argument(
        "--gamma", type=float, default=1.0,
        help="Power exponent in t' = t^γ.  γ=1 → linear (default); γ=3-4 "
             "with --t-mode rank concentrates ink on top edges while "
             "keeping the bulk visible as faint hairlines.",
    )
    parser.add_argument(
        "--width-min", type=float, default=0.15,
        help="Minimum edge width.  Default 0.15.",
    )
    parser.add_argument(
        "--width-max", type=float, default=4.0,
        help="Maximum edge width.  Default 4.0.",
    )
    parser.add_argument(
        "--alpha-min", type=float, default=0.03,
        help="Minimum edge alpha.  Default 0.03.",
    )
    parser.add_argument(
        "--alpha-max", type=float, default=0.9,
        help="Maximum edge alpha.  Default 0.9.",
    )
    # -----------------------------------------------------------
    parser.add_argument(
        "--n-communities", type=int, default=10,
        help="LRG dendrogram cut for LRG-seeded layouts and "
             "node-color=community.  Default 10.",
    )
    parser.add_argument(
        "--node-size", type=int, default=22,
        help="Node marker size (default 22 for single-panel).",
    )
    parser.add_argument(
        "--watermark", action="store_true",
        help="Add a small grey provenance string bottom-right "
             "(off by default).",
    )
    parser.add_argument(
        "--out-dir",
        default=str(FIGURES_ROOT / "network_templates" / "single_network"),
    )
    args = parser.parse_args()

    A, probes = _load_inputs(args.patient, args.band, args.phase, args.fc_method)

    community_labels = None
    ultrametric_matrix = None
    if is_lrg_layout(args.layout) or args.node_color == "community":
        lrg = load_lrg_result(args.patient, args.phase, args.band, args.fc_method)
        if lrg is None or lrg.linkage_matrix is None:
            raise SystemExit(
                f"No LRG cache for {args.patient}/{args.band}/{args.phase}/"
                f"{args.fc_method}; cannot run layout={args.layout} or "
                f"node-color=community."
            )
        N = A.shape[0]
        community_labels = fcluster(
            lrg.linkage_matrix, t=args.n_communities, criterion="maxclust",
        )[:N]
        ultrametric_matrix = getattr(lrg, "ultrametric_matrix", None)

    width_range = (args.width_min, args.width_max)
    alpha_range = (args.alpha_min, args.alpha_max)

    layout_kwargs = {}
    if args.layout == "spring":
        if args.k is not None:
            layout_kwargs["k"] = args.k
        layout_kwargs["iterations"] = args.iterations

    fig, ax = plt.subplots(figsize=(6.5, 6.8))
    plot_fc_network(
        A, probes,
        ax=ax,
        layout=args.layout,
        layout_kwargs=layout_kwargs,
        coloring=args.coloring,
        fc_method=args.fc_method,
        community_labels=community_labels,
        ultrametric_matrix=ultrametric_matrix,
        node_color=args.node_color,
        node_size=args.node_size,
        gamma=args.gamma,
        width_range=width_range,
        alpha_range=alpha_range,
        t_mode=args.t_mode,
        min_alpha=args.min_alpha,
    )
    coloring_resolved = (
        args.coloring
        or ("signed" if args.fc_method in ("corr", "imcoh") else "probe")
    )
    ax.set_title(
        _scaling_title(
            patient=args.patient, band=args.band, phase=args.phase,
            fc_method=args.fc_method, layout=args.layout,
            coloring=coloring_resolved,
            gamma=args.gamma,
            width_range=width_range, alpha_range=alpha_range,
            k=args.k if args.layout == "spring" else None,
            iterations=args.iterations if args.layout == "spring" else None,
            t_mode=args.t_mode,
        ),
        fontsize=9, fontweight="normal",
    )

    if args.watermark:
        add_provenance_footer(
            fig,
            f"single_network · {args.patient} · "
            f"{BRAIN_BAND_TEX_DICT.get(args.band, args.band)} · "
            f"{PHASE_TEX.get(args.phase, args.phase)} · {args.fc_method} · "
            f"{args.layout}",
        )

    scaling_tag = _scaling_filename_tag(
        gamma=args.gamma,
        width_range=width_range,
        alpha_range=alpha_range,
        t_mode=args.t_mode,
        min_alpha=args.min_alpha,
    )
    layout_tag = args.layout
    if args.layout == "spring" and args.k is not None:
        layout_tag = f"spring_k{args.k:g}_it{args.iterations}"
    elif args.layout == "spring":
        layout_tag = f"spring_it{args.iterations}"
    out = (
        Path(args.out_dir)
        / f"{args.patient}_{args.band}_{args.phase}_"
          f"{args.fc_method}_{layout_tag}_{args.node_color}_"
          f"edge-{coloring_resolved}_{scaling_tag}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
