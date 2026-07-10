#!/usr/bin/env python3
r"""fig:trace_c (companion) — the task's cophenetic fingerprint, and its echo (Results §1, para c).

IDEA (no correlation on screen — the eye does the correlating): for one patient, lay the
change the task wrote into the LRG cophenetic hierarchy next to the change still present at
post-task rest. Each map is the per-pair change in cophenetic distance, contacts ordered by
the rest_pre dendrogram so the stable backbone is the pale bulk and the movers are the
saturated blocks. When the task's fingerprint (top) reappears in the held map (bottom) —
same red/blue blocks — that IS the trace, shown as a pattern match instead of a small ρ.
It reappears for β and α (green labels, clear the cohort gate); it does not for θ and δ
(grey labels). The pooled Spearman is only ~0.2 because the backbone is unchanged; the point
of this figure is that among the pairs that move, the pattern genuinely persists.

HONESTY (why split-half baselines): the naive maps D(task)−D(pre) and D(post)−D(pre) share
the same D(pre), which manufactures a matching pattern from the backbone alone — full-phase
δ/θ then look as "matched" as β (probe: δ +0.51, θ +0.44 vs β +0.72). So each map is built
against an INDEPENDENT half of rest_pre (task vs pre_A, post vs pre_B) — the ρ_sym
construction — which removes the shared-baseline artifact and leaves only the real, persistent
reorganization (β +0.50, α +0.34; θ −0.05, δ −0.11 on this patient). Exemplar Pat_08 (the
cleanest tracer); the cohort verdict is the gate/forest in (a).

Reads : LRG cophenetic caches (full task/post/pre + rest_pre halves), rho_sym gate summary.
Writes: data/preprint/figures/results_section1/fig_trace_c_diffmaps.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list
from scipy.spatial.distance import squareform
from scipy.stats import rankdata, spearmanr

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result

ROOT = setup_script_env()
use_lrg_style()

LRG_HALVES = CACHE_ROOT / "imcoh_lrg_halves"
GATE_SUMMARY = ROOT / "data" / "audit" / "rho_sym_gate" / "cohort_summary.csv"
OUT = ROOT / "data" / "preprint" / "figures" / "results_section1" / "fig_trace_c_diffmaps.pdf"

PATIENT = "Pat_08"
FC = "imcoh_abs"
BANDS = ["beta", "alpha", "low_gamma", "delta", "theta", "high_gamma"]  # all six bands
GATE_PASS, GATE_FAIL = "#1f7a34", "#565656"
# diverging: pair pulled together (blue) ↔ pulled apart (red); pale centre = backbone (no change)
DIVCMAP = "coolwarm"


def _Dsq(phase, band, cache_root=None):
    r = (load_lrg_result(PATIENT, phase, band, fc_method=FC, cache_root=cache_root)
         if cache_root else load_lrg_result(PATIENT, phase, band, fc_method=FC))
    um = np.asarray(r.ultrametric_matrix)
    return squareform(um) if um.ndim == 1 else um


def _rankmap(delta_full):
    """Signed-rank image of the off-diagonal Δ, in [-1, 1]. Displaying RANKS makes the
    picture identical to what Spearman ρ measures: a global offset (θ: everything pulled
    together at task, apart at post) becomes a uniform rank field and disappears, leaving
    only the co-varying block structure — the actual trace."""
    n = delta_full.shape[0]
    iu = np.triu_indices(n, k=1)
    scaled = 2.0 * (rankdata(delta_full[iu]) - 1.0) / (iu[0].size - 1) - 1.0
    R = np.zeros((n, n))
    R[iu] = scaled
    return R + R.T


def band_maps(band):
    """(Δtask, Δrest) reordered by the rest_pre dendrogram + the arm-1 ρ_split (diagnostic)."""
    Z_pre = load_lrg_result(PATIENT, "rest_pre", band, fc_method=FC).linkage_matrix
    order = leaves_list(Z_pre)
    Dtask = _Dsq("task_test", band)
    Dpost = _Dsq("rest_post", band)
    DpreA = _Dsq("rest_pre_A", band, LRG_HALVES)
    DpreB = _Dsq("rest_pre_B", band, LRG_HALVES)
    dtask = Dtask - DpreA                              # task vs an independent half of rest_pre
    drest = Dpost - DpreB                              # held change vs the OTHER half
    iu = np.triu_indices(dtask.shape[0], k=1)
    rho = float(spearmanr(dtask[iu], drest[iu]).correlation)
    ix = np.ix_(order, order)
    return _rankmap(dtask)[ix], _rankmap(drest)[ix], rho


def draw(target):
    ncol = len(BANDS)
    gate = pd.read_csv(GATE_SUMMARY).set_index("band")["gate_p_sym"].to_dict()
    maps = {b: band_maps(b) for b in BANDS}

    print(f"fig:trace_c diffmaps — {PATIENT}, split-half baselines (ρ_sym arm-1 diagnostic):")
    for b in BANDS:
        print(f"  {b:11s} pattern-match ρ = {maps[b][2]:+.3f}   "
              f"gate_p_sym={gate.get(b, float('nan')):.3f}")

    # Two rows of square (aspect="equal") imshow maps with a small, clean inter-row gap.
    # The imshow images fill their axes edge-to-edge (unlike the former pcolormesh, whose
    # margin left a hairline), so hspace=0 makes the rows read as touching/overlapping; a
    # modest hspace separates them. The two-row block height is sized from the tile w/h so
    # the cells stay ~square, with the gap folded in so nothing overflows into panel c.
    left, right, wspace = 0.052, 0.906, 0.06
    hspace = 0.10                                        # inter-row gap (frac of cell height)
    aw = (right - left) / (ncol + (ncol - 1) * wspace)   # one square's side as fig-fraction of width
    tile_ar = target.bbox.width / target.bbox.height     # tile aspect (w/h)
    top = 0.90
    span = min((2.0 + hspace) * aw * tile_ar, top - 0.045)  # two rows + gap; never overflow bottom
    bottom = top - span
    gs = target.add_gridspec(2, ncol, left=left, right=right, top=top, bottom=bottom,
                             wspace=wspace, hspace=hspace)

    row_tex = [r"$\mathrm{task_{test}}$ vs $\mathrm{rest_{pre}}$",
               r"$\mathrm{rest_{post}}$ vs $\mathrm{rest_{pre}}$"]
    last = None
    for j, b in enumerate(BANDS):
        dtask, drest, _ = maps[b]
        lab_col = GATE_PASS if gate.get(b, 1.0) < 0.05 else GATE_FAIL
        for i, M in enumerate((dtask, drest)):
            ax = target.add_subplot(gs[i, j])
            # Native-resolution raster (imshow) rather than vector pcolormesh: an N^2
            # fingerprint embeds as an ~N x N image (tens of KB) instead of ~N^2 vector
            # quads, nearest-neighbour upscaled crisp and lossless at the data resolution.
            # Deliberate, authorized exception to the "FC matrices stay vector" rule, for
            # file size (twelve matrices in one compound). origin='upper' puts M[0,0] at
            # top-left, matching the former pcolormesh(M)+invert_yaxis, so no invert here.
            last = ax.imshow(M, cmap=DIVCMAP, vmin=-1.0, vmax=1.0,
                             interpolation="nearest", aspect="equal")
            ax.set_xticks([])
            ax.set_yticks([])
            for s in ax.spines.values():
                s.set_edgecolor(lab_col if i == 0 else "0.6")
                s.set_linewidth(1.6 if i == 0 else 0.8)
            if i == 0:
                ax.set_title(BRAIN_BAND_TEX_DICT.get(b, b), fontsize=19,
                             fontweight="bold", color=lab_col, pad=5)
            if j == 0:
                ax.set_ylabel(row_tex[i], fontsize=12.5, color="0.20")

    target.text(0.012, 0.925, r"$\mathbf{c}$", fontsize=17, va="bottom", ha="left",
                fontweight="bold")

    # colorbar carries ONLY the three symbols (-, 0, +); their meaning is in the
    # caption, not on the figure.
    cax = target.add_axes([0.918, bottom + 0.06 * span, 0.013, 0.88 * span])
    cb = target.colorbar(last, cax=cax, extend="both")
    cb.set_ticks([-1, 0, 1])
    cb.set_ticklabels([r"$-$", "0", r"$+$"])
    cb.ax.tick_params(labelsize=16)


def main():
    ncol = len(BANDS)
    # figure height chosen so the square (box_aspect=1) maps fill their grid
    # cells with no vertical slack — otherwise the mismatch reads as whitespace.
    fig = plt.figure(figsize=(3.02 * ncol, 6.05))
    draw(fig)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
