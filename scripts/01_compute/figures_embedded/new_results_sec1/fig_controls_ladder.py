#!/usr/bin/env python3
r"""fig_controls_ladder -- the multiscale read-out is necessary, shown WITHOUT a ground truth.

The band-selective held trace is unique to the multiscale cophenetic read-out. We
demonstrate that without ever asserting another descriptor picked the "wrong" band
-- the discredited dense/full-graph result is NOT a ground truth. Every descriptor is
read off the SAME mst@0.20 backbone (apples-to-apples; controls_ladder_apples) against
the SAME matched-strength null, and a cell counts as a REPRESENTATIVE clear only if it
passes a gate with no target band in it (cohort p<0.05 AND leave-one-out AND >=3/10
patients above their own null; see _common.representativeness_gate).

Reading after the gate (T_test / trace), all on the backbone:
  - node strength / clustering / resistance : no representative clear anywhere.
  - raw FC edges     : representative in beta only -> recovers the beta carrier but not alpha.
  - graph geodesic   : one representative clear, in low_gamma -- a NON-carrier band.
  - cophenetic (ours): the ONLY read-out representative in BOTH carrier bands (alpha, beta)
                       and silent elsewhere; alpha is cophenetic-EXCLUSIVE, beta a
                       convergence band (raw backbone edges also carry it).

The band identity of the trace lives in the scale (tau) dependence (fig_tau_role); the
competitors are single-scale/edge/spectral and have no axis on which that signature can
exist. Which bands the hierarchy flags is anchored by external anatomy (beta->OFC).

Encoding: dot AREA = per-patient support; FILL = representative (band colour, solid),
RING = cohort-significant but not representative (band colour, hollow), grey open = null.

Reads (CSV only): controls_ladder_apples/per_cell.csv (functional == "T_test").
Writes: data/preprint/figures/new_results_sec1/fig_controls_ladder.pdf
"""
from __future__ import annotations
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import _common as C

F = "T_test"


def main():
    C.use_lrg_style()
    df = pd.read_csv(C.LADDER / "per_cell.csv")

    descs, bands = C.LADDER_DESCS, list(C.BANDS)
    n_desc = len(descs)
    fig, ax = plt.subplots(figsize=(9.0, 5.6))

    # faint separator setting off the nested-hierarchy (cophenetic) row at the bottom
    ax.axhline(len(C.LADDER_HIER) - 0.5, color="0.85", lw=0.9, ls=(0, (4, 3)), zorder=0)

    for i, (desc, _) in enumerate(descs):
        y = n_desc - 1 - i                       # first descriptor on top
        for j, band in enumerate(bands):
            C.plot_gate_cell(ax, j, y, band, C.representativeness_gate(df, desc, band, F))

    ax.set_xlim(-0.7, len(bands) - 0.3)
    ax.set_ylim(-0.7, n_desc - 0.3)
    ax.set_xticks(range(len(bands)))
    ax.set_xticklabels([C.BTeX[b] for b in bands], fontsize=13)
    ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
    ax.set_yticks([n_desc - 1 - i for i in range(n_desc)])
    ax.set_yticklabels([lab for _, lab in descs], fontsize=11)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
    ax.set_title("Held trace across the descriptor ladder "
                 r"(T$_\mathrm{test}$, mst@0.20; cohort gate $+$ LOO $+$ consistency)",
                 fontsize=11, fontweight="bold", loc="left", pad=26)

    size_h = [plt.scatter([], [], s=C.gate_area(k), facecolor="0.45", edgecolor="white",
                          lw=0.6, label=f"{k}/10") for k in (2, 5, 8)]
    lvl_h = [
        Line2D([], [], marker="o", ls="none", ms=10, mfc="0.45", mec="white", label="representative"),
        Line2D([], [], marker="o", ls="none", ms=10, mfc="none", mec="0.35", mew=1.9,
               label="sig. but not representative"),
        Line2D([], [], marker="o", ls="none", ms=10, mfc="none", mec="0.6", label="null"),
    ]
    leg1 = fig.legend(handles=size_h, title="patients above own null",
                      loc="lower left", bbox_to_anchor=(0.09, -0.04), ncol=3,
                      frameon=False, fontsize=9, title_fontsize=9, handletextpad=0.2,
                      columnspacing=1.1)
    fig.add_artist(leg1)
    fig.legend(handles=lvl_h, loc="lower right", bbox_to_anchor=(0.98, -0.04),
               ncol=3, frameon=False, fontsize=9, handletextpad=0.3, columnspacing=1.1)

    C.FIGDIR.mkdir(parents=True, exist_ok=True)
    out = C.FIGDIR / "fig_controls_ladder.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


if __name__ == "__main__":
    main()
