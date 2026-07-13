#!/usr/bin/env python3
r"""fig:ribbon-tree (mst@0.20 backbone) — the retained rest_post hierarchy
(branch-origin dendrogram) + per-leaf encoding/inference ribbon, rebuilt on the
**sparsified multiscale backbone** at a single diffusion scale s = tau*lambda_max ~ 5.65.

FORK of ``figures_embedded/fig_ribbon_tree.py``. The ONLY change is the tree/cophenetic
substrate: the dense single-scale LRG tree (``tau = 1/lambda_max`` on the fully-connected
FC, via ``audit_65.lrg_linkage`` / ``fig_branch_origin_tree.coph_square``) is replaced by
the ``mst_union_top_fraction(W, 0.20)`` backbone addressed at the reporting scale
``s = S_REPORT`` (``_common.tree_at_scale`` / ``coph_square_at_scale``). The
encoding/inference branch-origin attribution, the diverging e/f ribbon, the log
cophenetic-height axis and every style choice are inherited unchanged from the machinery
module (DRY: ``attribute`` / ``link_color`` / ``leaf_origin_weights`` are imported, not
copied).

Does NOT modify ``_common.py``, ``fig_branch_origin_tree.py`` or the original figure script.
Reads : phase FC via ``_common.load_phase`` (imcoh_abs; identical normalisation to
        ``audit_63.load_phase_fc`` for the four attribution phases).
Writes: data/preprint/figures/new_results_sec1/_drafts/
        fig_ribbon_tree_{patient}_{band}_mst020_DRAFT.pdf
"""
from __future__ import annotations

import argparse
import sys

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import squareform

# --- sparsified multiscale backbone layer (single source; do not edit) ---
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded" / "new_results_sec1"))
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded"))
from _common import (                                                  # type: ignore
    load_phase, S_REPORT, FRAC, eig_backbone, use_lrg_style, DRAFTDIR,
)
from lrg_eegfc.utils.fc.heat_multiscale import linkage_at_scale, cophenetic_at_scale

# --- shared branch-origin machinery (pure helpers; DRY, unchanged) ---
from fig_branch_origin_tree import (                                   # type: ignore
    attribute, link_color, leaf_origin_weights,
    AMBER, TEAL, GREY, BAND_SYM, PHASES,           # PHASES = 4-phase attribution list
)


def coph_tree_ms(W):
    """``(N x N cophenetic matrix, linkage Z)`` of the mst@FRAC tree at scale ``S_REPORT``.

    Drop-in for ``fig_branch_origin_tree.coph_square(W)`` (which is fixed at the dense
    ``tau = 1/lambda_max``). One eigendecomposition drives both the linkage skeleton and
    the cophenetic square; byte-identical to calling ``_common.tree_at_scale`` and
    ``coph_square_at_scale`` separately (``cophenetic_at_scale``/``linkage_at_scale`` share
    one body), just without the redundant second eig.
    """
    ev, V = eig_backbone(W, FRAC)
    Z = linkage_at_scale(ev, V, S_REPORT)
    Csq = squareform(cophenetic_at_scale(ev, V, S_REPORT))
    return Csq, Z


def build_tree_ms(pat, band):
    """rest_post mst@0.20 linkage + per-phase cophenetic matrices (aligned leaves)."""
    Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    ns = {ph: W.shape[0] for ph, W in Ws.items()}
    if len(set(ns.values())) != 1:
        print(f"  !! {pat} {band}: phase N mismatch {ns} — skipping (cannot align pairs)")
        return None
    coph = {}
    Z_rp = None
    for ph, W in Ws.items():
        C, Z = coph_tree_ms(W)
        coph[ph] = C
        if ph == "rest_post":
            Z_rp = Z
    return Z_rp, coph


def draw(pat, band):
    got = build_tree_ms(pat, band)
    if got is None:
        return
    Z, coph = got
    N = Z.shape[0] + 1
    efp, members = attribute(Z, coph)
    pscale = np.percentile([abs(p) for _, _, p in efp.values()], 75) or 1e-9
    enc, inf, tot = leaf_origin_weights(Z, N, efp)

    dend = dendrogram(Z, no_plot=True)
    order = dend["leaves"]
    xpos = {leaf: i for i, leaf in enumerate(order)}
    heights = np.sort(Z[:, 2])
    base = heights[0] * 0.8

    def cx(c):
        return float(np.mean([xpos[l] for l in members[c]]))

    def ch(c):
        return base if c < N else Z[c - N, 2]

    use_lrg_style()
    fig = plt.figure(figsize=(9.2, 6.4))
    gs = GridSpec(2, 1, height_ratios=[3.0, 1.25], hspace=0.06, figure=fig)
    ax_t = fig.add_subplot(gs[0])
    ax_r = fig.add_subplot(gs[1], sharex=ax_t)

    # --- top: branch-origin dendrogram (vivid/persistent branches on top) ---
    links = []
    for k in range(N - 1):
        c1, c2 = int(Z[k, 0]), int(Z[k, 1])
        col, lw, pw = link_color(*efp[N + k], pscale)
        links.append((pw, cx(c1), cx(c2), ch(c1), ch(c2), Z[k, 2], col, lw))
    for pw, x1, x2, h1, h2, h, col, lw in sorted(links, key=lambda t: t[0]):
        ax_t.plot([x1, x1], [h1, h], color=col, lw=lw, solid_capstyle="round", zorder=2 + pw)
        ax_t.plot([x2, x2], [h2, h], color=col, lw=lw, solid_capstyle="round", zorder=2 + pw)
        ax_t.plot([x1, x2], [h, h], color=col, lw=lw, solid_capstyle="round", zorder=2 + pw)
    ax_t.set_yscale("log")
    ax_t.set_ylim(base, heights[-1] * 1.05)
    ax_t.set_ylabel(r"cophenetic height $1/\hat\rho$", fontsize=10)
    for s in ("top", "right", "bottom"):
        ax_t.spines[s].set_visible(False)
    ax_t.tick_params(axis="x", length=0, labelbottom=False)

    # --- bottom: diverging per-leaf e/f ribbon (inference up, encoding down) ---
    x = np.arange(N)
    inf_o, enc_o = inf[order], enc[order]
    s = 0.96 / max(inf_o.max(), enc_o.max(), 1e-12)        # per-figure relative scale
    ax_r.fill_between(x, 0, inf_o * s, color=TEAL, alpha=0.9, lw=0, zorder=3)
    ax_r.fill_between(x, 0, -enc_o * s, color=AMBER, alpha=0.9, lw=0, zorder=3)
    ax_r.plot(x, inf_o * s, color=np.array(TEAL) * 0.7, lw=0.6, zorder=4)
    ax_r.plot(x, -enc_o * s, color=np.array(AMBER) * 0.8, lw=0.6, zorder=4)
    ax_r.axhline(0, color="#3a3f47", lw=0.8, zorder=5)
    ax_r.set_ylim(-1.05, 1.05)
    ax_r.set_xlim(-1.0, N)
    ax_r.set_yticks([-0.5, 0.5])
    ax_r.set_yticklabels(["encoding", "inference"], fontsize=8.5)
    ax_r.set_xticks([])
    ax_r.set_xlabel("contacts  (retained-hierarchy / dendrogram order)", fontsize=9.5)
    for s_ in ("top", "right"):
        ax_r.spines[s_].set_visible(False)
    ax_r.spines["left"].set_visible(False)
    ax_r.tick_params(axis="y", length=0)

    handles = [Patch(facecolor=AMBER, label="encoding (premises)"),
               Patch(facecolor=TEAL, label="inference (reasoning)"),
               Patch(facecolor=GREY, label="baseline / reverted"),
               Line2D([], [], color=(0.4, 0.42, 0.46), lw=3.4,
                      label="branch thickness ∝ persistence into rest")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.03),
               ncol=4, frameon=False, fontsize=9, handlelength=1.4, columnspacing=1.8)
    ax_t.text(0.004, 0.98, f"{pat} · {BAND_SYM.get(band, band)}  ·  mst@0.20  s≈5.6",
              transform=ax_t.transAxes, ha="left", va="top", fontsize=11.5,
              color="#2b2f36", fontweight="bold")

    DRAFTDIR.mkdir(parents=True, exist_ok=True)
    out = DRAFTDIR / f"fig_ribbon_tree_{pat}_{band}_mst020_DRAFT.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True, pad_inches=0.03)
    plt.close(fig)
    fi = float(np.mean(inf_o > enc_o))
    print(f"wrote {out.name}   inference-leaning leaves: {fi:.2f}   N={N}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--patients", nargs="+", default=["Pat_06", "Pat_05"],
                    help="patient id(s); default the theta exemplars")
    ap.add_argument("--bands", nargs="+", default=["theta"],
                    help="band(s); default theta. any of "
                         "delta theta alpha beta low_gamma high_gamma")
    a = ap.parse_args()
    for pat in a.patients:
        for band in a.bands:
            draw(pat, band)


if __name__ == "__main__":
    main()
