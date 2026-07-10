#!/usr/bin/env python3
r"""fig:ribbon-tree — the retained rest_post hierarchy (branch-origin dendrogram) with a
per-leaf encoding/inference RIBBON aligned beneath it (§2, single patient).

CORE MESSAGE. The branch-origin tree colours each *merge* by the stage that wrote it; this
figure pushes the same decomposition down to each *contact*. Under the tree, a diverging
ribbon shows — for every leaf, in tree order so clades stay together — how much of the
structure it sits under was written at INFERENCE (teal, up) vs ENCODING (amber, down):

    for leaf l:  inf[l] = Σ_{merges ν on root-path} |f_ν|·|p_ν|   (inference-written, held)
                 enc[l] = Σ_{merges ν on root-path} |e_ν|·|p_ν|   (encoding-written, held)
    e_ν = h_learn−h_pre (encoding)   f_ν = h_test−h_learn (inference)   p_ν = h_post−h_pre
    (persistence);  h_X(ν) = mean cophenetic height of ν's leaf-pair block in phase X.

So a clade that reorganised while REASONING and held it shows a fat teal ribbon; a clade
carried from ENCODING shows amber. The ribbon is the concrete, per-contact version of the
single node tint in the chord-network figure (same leaf_origin_weights).

HONEST REGISTER. Illustrative single-patient decomposition (full-phase cophenetic distances,
descriptive attribution), NOT a gated test — the cohort claim stays the matched-strength
inference gate T_infspec_pe (β p=0.010, 6/10). Ribbon magnitude is scaled per-figure
(relative), so read the BALANCE (teal vs amber) and the clade pattern, not absolute heights.
Runs for ANY band / patient; off-trace bands show a mostly-grey tree + flat ribbon.

Reads : phase FC via audit_63.load_phase_fc (imcoh_abs); LRG tree via audit_65.lrg_linkage.
Writes: data/preprint/figures/_drafts/fig_ribbon_tree_{patient}_{band}_DRAFT.pdf
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

from lrg_eegfc.visuals.styles import use_lrg_style

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded"))
from fig_branch_origin_tree import (                                 # type: ignore  (DRY)
    build_tree, attribute, link_color, leaf_origin_weights,
    AMBER, TEAL, GREY, BAND_SYM,
)
from scipy.cluster.hierarchy import dendrogram

OUTDIR = ROOT / "data/preprint/figures/_drafts"


def draw(pat, band):
    got = build_tree(pat, band)
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
    ax_t.text(0.004, 0.98, f"{pat} · {BAND_SYM.get(band, band)}", transform=ax_t.transAxes,
              ha="left", va="top", fontsize=11.5, color="#2b2f36", fontweight="bold")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / f"fig_ribbon_tree_{pat}_{band}_DRAFT.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True, pad_inches=0.03)
    qa = ("/tmp/claude-1000/-home-giulio-Documents-research-neural-networks-lrgeegfc/"
          "c2197e94-c496-4e4e-97cc-89e43564e7da/scratchpad")
    fig.savefig(f"{qa}/ribbontree_{pat}_{band}.png", dpi=150, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    fi = float(np.mean(inf_o > enc_o))
    print(f"wrote {out.name}   inference-leaning leaves: {fi:.2f}   N={N}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--patients", nargs="+", default=["Pat_06"],
                    help="patient id(s); default a beta tracer (Pat_06)")
    ap.add_argument("--bands", nargs="+", default=["beta"],
                    help="band(s); default beta. any of "
                         "delta theta alpha beta low_gamma high_gamma")
    a = ap.parse_args()
    for pat in a.patients:
        for band in a.bands:
            draw(pat, band)


if __name__ == "__main__":
    main()
