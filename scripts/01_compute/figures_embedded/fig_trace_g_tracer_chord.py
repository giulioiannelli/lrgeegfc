#!/usr/bin/env python3
r"""fig:trace_g (hero) — the held multiscale trace, drawn on the network.

IDEA: show what a "held multiscale trace" looks like ON the network. We draw one exemplar
patient's post-task-rest beta network as a hierarchy-bundled chord (graph-tool bezier edges
routed through the LRG depth-2 hierarchy): every coupling is a faint grey backbone, and the
GREEN links are the "tracer" edges — the contact pairs whose multiscale (cophenetic)
relationship the task reshaped and post-task rest KEPT.

A pair is a tracer if its cophenetic distance moved the same way at task and at rest,
relative to two independent rest_pre halves (the per-pair reading of rho_sym):

    score = 0.5 * [ rankz(D_task - D_preA) . rankz(D_post - D_preB)
                  + rankz(D_task - D_preB) . rankz(D_post - D_preA) ]

score > 0  = the pair's task-induced cophenetic shift persists into rest_post (a tracer).
Green = the top-K tracer pairs; edge width/alpha scale with the score.

ILLUSTRATIVE, NOT A PER-EDGE CLAIM. Individual per-pair tracer links are noisy; the
load-bearing evidence for the trace is node/system-level (beta -> OFC), never a single arc.
This panel is an exemplar (Pat_08, beta) that makes the measure legible on the network
itself — the cohort statistics live in the other panels. Nodes are drawn neutral grey (NOT
colored by probe) so the eye follows the traced edges.

NB: the chord backbone is rendered by graph-tool to a high-res PNG and embedded (this genre
is raster-by-construction, unlike the vector panels).

Reads : imcoh_abs FC + LRG cophenetic caches (task_test/rest_post + rest_pre halves), Pat_08 beta.
Writes: data/preprint/figures/results_section1/fig_trace_g_tracer_chord.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import squareform
from scipy.stats import rankdata

from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.network_templates import _load_inputs, plot_chord_with_highlight
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result

ROOT = setup_script_env()
use_lrg_style()

OUT = ROOT / "data" / "preprint" / "figures" / "results_section1" / "fig_trace_g_tracer_chord.pdf"
LRG_HALVES = CACHE_ROOT / "imcoh_lrg_halves"

PATIENT, BAND, FC = "Pat_08", "beta", "imcoh_abs"
PHASE_BG = "rest_post"                 # backbone network = the held state
K_GREEN = 120                          # number of tracer links highlighted
CLR_TRACE = (0.10, 0.49, 0.24)         # forest green
NODE_RGBA = (0.74, 0.74, 0.74, 1.0)    # neutral grey nodes (not probe-colored)


def _Dsq(phase, cache_root=None):
    r = (load_lrg_result(PATIENT, phase, BAND, fc_method=FC, cache_root=cache_root)
         if cache_root else load_lrg_result(PATIENT, phase, BAND, fc_method=FC))
    um = np.asarray(r.ultrametric_matrix)
    return squareform(um) if um.ndim == 1 else um


def _rankz(x):
    """Rank-transform to [-1, 1] — matches the Spearman the trace is measured with."""
    return 2.0 * (rankdata(x) - 1.0) / (x.size - 1) - 1.0


def tracer_scores(N):
    """Per-pair persistent-mover score (per-pair reading of rho_sym); > 0 = tracer."""
    Dt = _Dsq("task_test")[:N, :N]
    Dp = _Dsq("rest_post")[:N, :N]
    Da = _Dsq("rest_pre_A", LRG_HALVES)[:N, :N]
    Db = _Dsq("rest_pre_B", LRG_HALVES)[:N, :N]
    iu = np.triu_indices(N, k=1)
    c1 = _rankz((Dt - Da)[iu]) * _rankz((Dp - Db)[iu])
    c2 = _rankz((Dt - Db)[iu]) * _rankz((Dp - Da)[iu])
    return 0.5 * (c1 + c2)


def build_panel(ax):
    """Render the tracer chord onto ``ax`` (importable by the Section-1 mosaic)."""
    A, probes = _load_inputs(PATIENT, BAND, PHASE_BG, FC)
    lrg = load_lrg_result(PATIENT, PHASE_BG, BAND, fc_method=FC)
    N = min(A.shape[0], len(probes))

    score = tracer_scores(N)
    pos = np.where(score > 0)[0]
    top = pos[np.argsort(score[pos])[-K_GREEN:][::-1]]

    highlights = [dict(
        pairs=top, rgb=CLR_TRACE, magnitudes=score[top],
        width_range=(1.2, 6.2), alpha_range=(0.50, 0.95), gamma=0.55,
        draw_order_offset=2000.0, label="tracer links",
    )]
    plot_chord_with_highlight(
        ax, A, probes, lrg, highlights=highlights,
        background_alpha=0.075, vertex_size=18.0,
        vertex_rgba=np.tile(NODE_RGBA, (N, 1)), labels=None,
    )
    ax.text(0.0, 1.0, r"$\mathbf{a}$", transform=ax.transAxes, fontsize=17,
            va="top", ha="left", fontweight="bold")
    return N, int(top.size)


def main():
    fig, ax = plt.subplots(figsize=(7.4, 7.4))
    N, n_green = build_panel(ax)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"fig:trace_g tracer chord — {PATIENT} {BAND} ({PHASE_BG} backbone), "
          f"N={N} contacts, {n_green} green tracer links")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
