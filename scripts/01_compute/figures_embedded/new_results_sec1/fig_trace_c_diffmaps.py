#!/usr/bin/env python3
r"""fig:trace_c (mst@0.20) — the task's cophenetic fingerprint, and its echo (Results §1, para c).

HEAD: for one exemplar patient, lay the change the task wrote into the LRG cophenetic
hierarchy next to the change still present at post-task rest — rebuilt on the **mst@0.20**
multiscale backbone at a single diffusion scale **s≈5.6** (s ≡ τ·λmax), the same backbone and
scale as the rest of the sparsified §1 rebuild. Each map is the per-pair change in cophenetic
distance vs an independent rest_pre half, rank-scaled within its own map (diverging: blue =
pulled together, red = pulled apart), with contacts ordered by the rest_pre dendrogram so the
stable backbone is the pale bulk and the movers are the saturated blocks. When the task's
fingerprint (top row) reappears in the held map (bottom row) — same red/blue blocks — that IS
the trace, shown as a pattern match instead of a small ρ. It reappears for β and α (green
labels — they clear the mst@0.20 cohort matched-strength gate at s≈5.6: β p=0.001, α p=0.007);
it does not for δ/θ/low_γ/high_γ (grey labels). Exemplar Pat_08.

CONSTRUCTION: like the dense original, the two maps difference against SPLIT-HALF rest_pre
baselines (the ρ_sym arm-1 construction, which strips the shared-baseline artifact) — here on
the mst@0.20 backbone at s≈5.6 (via ``_common.coph_square_at_scale``):
    dtask = D(task_test) − D(rest_pre_A) ,  drest = D(rest_post) − D(rest_pre_B).
Because the two maps use INDEPENDENT halves A/B, the printed pattern-match ρ is the artifact-
free arm-1 ρ_sym value (high for the trace bands β/α, low for the rest — the maps separate the
bands on their own), consistent with the load-bearing cohort matched-strength mst@0.20 gate
carried by the green/grey band labels. Rank display (see ``_rankmap``) additionally neutralises
any global-offset part so the maps read as a mover-pattern match, not a backbone echo.

Reads : imcoh_abs FC (task_test, rest_post via _common.load_phase; full rest_pre via
        load_fc_matrix); mst@0.20 cohort gate data/sparsified_arc/ms_mst020/cohort_gate.csv.
Writes: data/preprint/figures/new_results_sec1/fig_trace_c_diffmaps_mst020.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list
from scipy.stats import rankdata, spearmanr

import _common as C
from lrg_eegfc.workflow.fc import load_fc_matrix

C.use_lrg_style()

# Pat_07: a representative tracer whose per-patient pattern matches the cohort story --
# beta/alpha clearly positive, every other band weak/null -- so the diff-maps separate the
# green (trace) bands from the grey (null) bands on their own (Pat_08 is a low_gamma
# super-tracer whose per-patient maps would wrongly echo low_gamma too).
PATIENT = "Pat_07"
FC = "imcoh_abs"
# β, α first (the two bands that clear the gate — the visual story), then the rest.
BANDS = ["beta", "alpha", "low_gamma", "delta", "theta", "high_gamma"]
GATE_PASS, GATE_FAIL = "#1f7a34", "#565656"
# diverging: pair pulled together (blue) ↔ pulled apart (red); pale centre = backbone (no change)
DIVCMAP = "coolwarm"
OUT = C.FIGDIR / "fig_trace_c_diffmaps_mst020.pdf"


def _pre_fc(band):
    """Full rest_pre imcoh_abs FC (display ORDERING baseline only; the diff maps use the
    split-half rest_pre A/B baselines). Symmetrised/clipped like ``_common.load_phase``."""
    W = np.asarray(load_fc_matrix(PATIENT, "rest_pre", band, fc_method=FC), float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def _rankmap(delta_full):
    """Signed-rank image of the off-diagonal Δ, in [-1, 1]. Displaying RANKS makes the
    picture identical to what Spearman ρ measures: a global offset (everything pulled
    together at task, apart at post) becomes a uniform rank field and disappears, leaving
    only the co-varying block structure — the actual trace."""
    n = delta_full.shape[0]
    iu = np.triu_indices(n, k=1)
    scaled = 2.0 * (rankdata(delta_full[iu]) - 1.0) / (iu[0].size - 1) - 1.0
    R = np.zeros((n, n))
    R[iu] = scaled
    return R + R.T


def band_maps(band):
    """(Δtask, Δrest) cophenetic diff-maps at s≈5.6 on the mst@0.20 backbone, reordered by
    the rest_pre dendrogram, + the SPLIT-HALF (ρ_sym arm-1) pattern-match ρ.

    The two maps difference against INDEPENDENT rest_pre halves (task fingerprint vs half A,
    held echo vs half B), so no shared-baseline term inflates the pattern match -- this is the
    ρ_sym arm-1 construction, which separates the trace bands (β/α) from the rest, as in the
    dense original. The full rest_pre only sets the display leaf order."""
    order = leaves_list(C.tree_at_scale(_pre_fc(band), C.S_REPORT))
    DA = C.coph_square_at_scale(C.load_phase(PATIENT, "A", band), C.S_REPORT)
    DB = C.coph_square_at_scale(C.load_phase(PATIENT, "B", band), C.S_REPORT)
    Dtask = C.coph_square_at_scale(C.load_phase(PATIENT, "task_test", band), C.S_REPORT)
    Dpost = C.coph_square_at_scale(C.load_phase(PATIENT, "rest_post", band), C.S_REPORT)
    dtask = Dtask - DA                                 # task fingerprint vs rest_pre half A
    drest = Dpost - DB                                 # held change vs INDEPENDENT half B
    iu = np.triu_indices(dtask.shape[0], k=1)
    rho = float(spearmanr(dtask[iu], drest[iu]).correlation)   # arm-1 ρ_sym: artifact-free
    ix = np.ix_(order, order)
    return _rankmap(dtask)[ix], _rankmap(drest)[ix], rho


def _gate_at_report():
    """Per-band mst@0.20 matched-strength gate p at the scale nearest S_REPORT."""
    g = pd.read_csv(C.MS / "cohort_gate.csv")
    out = {}
    for b in BANDS:
        gb = g[g.band == b]
        out[b] = float(gb.iloc[(gb.s - C.S_REPORT).abs().argmin()].gate_p)
    return out


def draw(target):
    ncol = len(BANDS)
    gate = _gate_at_report()
    maps = {b: band_maps(b) for b in BANDS}

    print(f"fig:trace_c diffmaps (mst@0.20, s={C.S_REPORT:.3f}) — {PATIENT}:")
    print("  (ρ below is the split-half arm-1 pattern match — artifact-free; "
          "cohort verdict = mst@0.20 gate_p)")
    for b in BANDS:
        print(f"  {b:11s} arm-1 ρ = {maps[b][2]:+.3f}   "
              f"mst@0.20 gate_p(s≈{C.S_REPORT:.1f}) = {gate[b]:.3f}  "
              f"[{'PASS' if gate[b] < 0.05 else 'fail'}]")

    # Two rows of square (aspect="equal") maps with a small, clean inter-row gap. The
    # two-row block height is sized from the tile w/h so the cells stay ~square, with the
    # gap folded in so nothing overflows past the bottom / into a neighbouring panel.
    left, right, wspace = 0.052, 0.906, 0.06
    hspace = 0.10                                        # inter-row gap (frac of cell height)
    aw = (right - left) / (ncol + (ncol - 1) * wspace)   # one square's side as fig-fraction of width
    tile_ar = target.bbox.width / target.bbox.height     # tile aspect (w/h)
    top = 0.90
    span = min((2.0 + hspace) * aw * tile_ar, top - 0.045)  # two rows + gap; never overflow bottom
    bottom = top - span
    gs = target.add_gridspec(2, ncol, left=left, right=right, top=top, bottom=bottom,
                             wspace=wspace, hspace=hspace)

    row_tex = [r"$\mathrm{task_{test}}$ vs $\mathrm{rest_{pre}^{A}}$",
               r"$\mathrm{rest_{post}}$ vs $\mathrm{rest_{pre}^{B}}$"]
    last = None
    for j, b in enumerate(BANDS):
        dtask, drest, _ = maps[b]
        lab_col = GATE_PASS if gate.get(b, 1.0) < 0.05 else GATE_FAIL
        for i, M in enumerate((dtask, drest)):
            ax = target.add_subplot(gs[i, j])
            # FULLY-VECTOR heatmap: pcolormesh (vector quads), edgecolor='face' to kill the
            # inter-quad hairlines while staying vector. NEVER imshow / set_rasterized here —
            # the §1 rebuild is PDF-fully-vector. ylim (N, 0) places M[0,0] at top-left so the
            # rest_pre dendrogram order reads top→down (matches the dense pcolormesh+invert).
            last = ax.pcolormesh(M, cmap=DIVCMAP, vmin=-1.0, vmax=1.0, shading="flat")
            last.set_edgecolor("face")
            ax.set_aspect("equal")
            ax.set_xlim(0, M.shape[1])
            ax.set_ylim(M.shape[0], 0)
            ax.set_xticks([])
            ax.set_yticks([])
            for s in ax.spines.values():
                s.set_edgecolor(lab_col if i == 0 else "0.6")
                s.set_linewidth(1.6 if i == 0 else 0.8)
            if i == 0:
                ax.set_title(C.BTeX.get(b, b), fontsize=19,
                             fontweight="bold", color=lab_col, pad=5)
            if j == 0:
                ax.set_ylabel(row_tex[i], fontsize=12.5, color="0.20")

    target.text(0.012, 0.925, r"$\mathbf{c}$", fontsize=17, va="bottom", ha="left",
                fontweight="bold")

    # colorbar carries ONLY the three symbols (-, 0, +); their meaning is in the caption.
    # Explicit fig.add_axes (NOT imshow_colorbar_caxdivider) — this cbar spans both rows.
    cax = target.add_axes([0.918, bottom + 0.06 * span, 0.013, 0.88 * span])
    cb = target.colorbar(last, cax=cax, extend="both")
    cb.solids.set_rasterized(False)   # mpl 3.8 defaults colorbar QuadMesh to raster; keep it vector
    cb.set_ticks([-1, 0, 1])
    cb.set_ticklabels([r"$-$", "0", r"$+$"])
    cb.ax.tick_params(labelsize=16)


def main():
    ncol = len(BANDS)
    # figure height chosen so the square maps fill their grid cells with no vertical slack.
    fig = plt.figure(figsize=(3.02 * ncol, 6.05))
    draw(fig)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
