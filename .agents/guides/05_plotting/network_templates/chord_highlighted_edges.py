"""chord_highlighted_edges — chord network with user-defined edge highlights.

Template 2 in the "hierarchy-bundled chord" family.  Crystallised from
the "trace + anti" overlay pass in panel (c) of
``preprint_07_beta_rho_split_figure_test2.py``, generalised so any
caller can name one or more highlight classes (each = a colour + a
list of pairs + optional per-pair magnitudes for width/alpha scaling)
and have every other pair drawn as a faint grey "null backbone"
underneath.

The default demo highlights the **top-K |w| cross-probe pairs** in
forest green and the **top-K same-probe pairs** in brick red, which is
a useful sanity check for any new patient/band combination — it tells
you whether the strongest cross-probe communication aligns with the
LRG hierarchy, or whether same-shaft bias dominates the top of the
distribution.

For the trace/anti use case, call ``plot_chord_with_highlight`` from
your own script with ``highlights=[{"pairs": ..., "rgb": (...)}, ...]``
and your loaded per-pair signed magnitudes.  See the README of this
folder for the trace-direction palette anchors.

Three independent layers (mirrors ``single_network``):

  1.  **Chord layout**      — graph-tool hierarchy-bundled bezier
      bundles through a depth-2 LRG cut (``k_clusters``, ``beta``).
  2.  **Background**        — every non-highlighted pair drawn as a
      faint grey hairline (``HIERARCHY_BG_RGB`` at
      ``HIERARCHY_BG_ALPHA``).  Optional cross-probe-only filter so
      same-shaft pairs don't drown the grey wash.
  3.  **Highlights**        — one or more colour classes, each drawn
      on top of the background with per-class ``width_range`` and
      ``alpha_range``.  If per-pair magnitudes are passed they scale
      width/alpha within the class via
      ``t = (m/m_max)**class_gamma``.

Companion ``chord_highlighted_edges.md`` carries the style sheet.

Usage
-----
    conda activate lapbrain
    # Default demo — top-K cross-probe green + top-K same-probe red
    python .agents/guides/05_plotting/network_templates/chord_highlighted_edges.py \\
        --patient Pat_05 --band beta --phase rest_pre

    # Highlight top-200 strongest |w| in a single colour (no red/green split)
    python .agents/guides/05_plotting/network_templates/chord_highlighted_edges.py \\
        --patient Pat_05 --band beta --phase rest_pre \\
        --highlight-mode topk --k-trace 200 --k-anti 0

    # Cross-probe-only background (kills the same-shaft grey clump)
    python .agents/guides/05_plotting/network_templates/chord_highlighted_edges.py \\
        --patient Pat_05 --band beta --phase rest_pre --cross-probe-bg
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT, PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_patient_metadata
from lrg_eegfc.visuals.network_templates import (
    HIERARCHY_DEFAULT_BETA,
    HIERARCHY_DEFAULT_FIT_VIEW,
    HIERARCHY_DEFAULT_K_CLUSTERS,
    HIERARCHY_DEFAULT_RENDER_PX,
    _load_inputs,
    plot_chord_with_highlight,
)
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result


# Anchored to PALETTE_COPH navy/forest from preprint_11_beta_anatomy_brain
# so paired figures share a hue family.  Override at the API level by
# passing your own ``highlights[k]["rgb"]``.
CLR_HL_PRIMARY = "#1a7c3e"      # forest green
CLR_HL_SECONDARY = "#c0392b"    # brick red

PHASE_TEX = {
    "rest_pre": r"rest$_{\mathrm{pre}}$",
    "task_learn": r"task$_{\mathrm{learn}}$",
    "task_test": r"task$_{\mathrm{test}}$",
    "rest_post": r"rest$_{\mathrm{post}}$",
}


def _parse_contact_label(raw: str) -> str:
    head = str(raw).split(",")[0]
    return "".join(c for c in head if c.isascii() and c.isalnum())


def _load_contact_labels(patient: str) -> list[str]:
    meta = load_patient_metadata(patient, SEEG_DATAPATH)
    if meta is None or "label" not in meta.columns:
        return []
    return [_parse_contact_label(s) for s in meta["label"].tolist()]


def _topk_demo_highlights(A, probes, k_trace: int, k_anti: int):
    """Demo: top-k_trace cross-probe + top-k_anti same-probe in |w|.

    Returns a list of highlight dicts matching
    :func:`plot_chord_with_highlight`'s ``highlights`` API.  Each
    highlight carries per-pair magnitudes so width/alpha scale within
    the class.
    """
    import matplotlib.colors as mcolors

    N = min(A.shape[0], len(probes))
    r, c = np.triu_indices(N, k=1)
    w = np.abs(A[r, c])

    probes_arr = np.asarray(probes[:N])
    cross_mask = probes_arr[r] != probes_arr[c]
    same_mask = ~cross_mask

    highlights = []

    if k_trace > 0 and cross_mask.any():
        cross_idx = np.where(cross_mask)[0]
        top = cross_idx[np.argsort(w[cross_idx])[-k_trace:][::-1]]
        rgb_primary = mcolors.to_rgb(CLR_HL_PRIMARY)
        highlights.append(dict(
            pairs=top,
            rgb=rgb_primary,
            magnitudes=w[top],
            width_range=(1.4, 6.4),
            alpha_range=(0.45, 0.95),
            gamma=0.55,
            draw_order_offset=2000.0,
            label=f"top-{k_trace} cross-probe |w|",
        ))

    if k_anti > 0 and same_mask.any():
        same_idx = np.where(same_mask)[0]
        top = same_idx[np.argsort(w[same_idx])[-k_anti:][::-1]]
        rgb_secondary = mcolors.to_rgb(CLR_HL_SECONDARY)
        highlights.append(dict(
            pairs=top,
            rgb=rgb_secondary,
            magnitudes=w[top],
            width_range=(1.4, 6.4),
            alpha_range=(0.45, 0.95),
            gamma=0.55,
            draw_order_offset=1000.0,
            label=f"top-{k_anti} same-probe |w|",
        ))

    return highlights


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--patient", default="Pat_05")
    parser.add_argument("--band", default="beta")
    parser.add_argument("--phase", default="rest_pre", choices=PHASE_LABELS)
    parser.add_argument(
        "--fc-method", default="imcoh_abs",
        choices=["imcoh_abs", "imcoh_sq", "msc"],
    )
    parser.add_argument(
        "--highlight-mode", default="topk",
        choices=["topk"],
        help="Currently only 'topk' is wired up — top-k cross-probe "
             "(primary) + top-k same-probe (secondary).  Custom "
             "highlight sets (e.g. trace/anti from a per-pair NPZ) "
             "are built by your own script calling "
             "``plot_chord_with_highlight`` directly.",
    )
    parser.add_argument(
        "--k-trace", type=int, default=120,
        help="Number of primary (default cross-probe) highlight edges.",
    )
    parser.add_argument(
        "--k-anti", type=int, default=40,
        help="Number of secondary (default same-probe) highlight edges. "
             "Set to 0 for a single-class highlight.",
    )
    parser.add_argument(
        "--cross-probe-bg", action="store_true",
        help="Restrict the BACKGROUND grey wash to cross-probe pairs "
             "(default: all pairs).  Highlight pairs are drawn either way.",
    )
    parser.add_argument(
        "--k-clusters", type=int, default=HIERARCHY_DEFAULT_K_CLUSTERS,
    )
    parser.add_argument(
        "--beta", type=float, default=HIERARCHY_DEFAULT_BETA,
    )
    parser.add_argument(
        "--render-px", type=int, default=HIERARCHY_DEFAULT_RENDER_PX,
    )
    parser.add_argument(
        "--fit-view", type=float, default=HIERARCHY_DEFAULT_FIT_VIEW,
    )
    parser.add_argument(
        "--show-dendrogram", action="store_true",
        help="Also overlay the full circular-dendrogram (Template 1+2 hybrid).",
    )
    parser.add_argument(
        "--no-labels", action="store_true",
        help="Suppress the radial contact-label ring.",
    )
    parser.add_argument(
        "--vertex-size", type=float, default=24.0,
    )
    parser.add_argument(
        "--out-dir",
        default=str(
            FIGURES_ROOT / "network_templates" / "chord_highlighted_edges"
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

    if args.highlight_mode == "topk":
        highlights = _topk_demo_highlights(
            A, probes, args.k_trace, args.k_anti,
        )
    else:
        raise ValueError(f"unknown highlight_mode: {args.highlight_mode}")

    labels = None if args.no_labels else _load_contact_labels(args.patient)

    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    plot_chord_with_highlight(
        ax, A, probes, lrg,
        highlights=highlights,
        k_clusters=args.k_clusters, beta=args.beta,
        render_px=args.render_px, fit_view=args.fit_view,
        cross_probe_only=args.cross_probe_bg,
        show_dendrogram=args.show_dendrogram,
        labels=labels,
        vertex_size=args.vertex_size,
    )

    band_tex = BRAIN_BAND_TEX_DICT.get(args.band, args.band)
    phase_tex = PHASE_TEX.get(args.phase, args.phase)
    class_str = " + ".join(h.get("label", f"class{k}")
                            for k, h in enumerate(highlights))
    bg_str = "cross-probe bg" if args.cross_probe_bg else "all-pair bg"
    extra = " + dendro" if args.show_dendrogram else ""
    ax.set_title(
        f"{args.patient}  ·  {band_tex}  ·  {phase_tex}  ·  "
        f"{args.fc_method}    chord (K={args.k_clusters}, "
        rf"$\beta$={args.beta:g})    {class_str}    {bg_str}{extra}",
        fontsize=9, loc="left", pad=6,
    )

    tag = f"topk{args.k_trace}-{args.k_anti}"
    if args.cross_probe_bg:
        tag += "_xbg"
    if args.show_dendrogram:
        tag += "_dendro"
    out = (
        Path(args.out_dir)
        / f"{args.patient}_{args.band}_{args.phase}_{args.fc_method}_"
          f"K{args.k_clusters}_b{args.beta:g}_{tag}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    print(f"size: {out.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()
