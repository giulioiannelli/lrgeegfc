#!/usr/bin/env python3
r"""talk_fig_taxonomy_quad — the four cross-phase fates as ONE side-by-side panel for the talk
slide "A 4-fold task-induced taxonomy of neuronal populations".

Composites the four exemplar tanglegrams — reset (Pat_10 β) · reorganised (Pat_08 δ) ·
anchor (Pat_08 low-γ) · trace (Pat_06 δ), the exact clades chosen by
fig_cophenetic_tanglegram_class and shipped as the current talk PNGs — into a single figure of
four EQUAL-HEIGHT columns. Text is stripped to the minimum the user asked for:
    1. in-bead contact (node) labels
    2. per-gap ρ + crossing counts
    3. the phase label above each dendrogram
    4. the category title
Removed: the italic Pat/N/descriptor subtitle, the grey ρ_pre,post note, the sub-clade colour
legend, and both honesty footnotes.

DRY: reuses fig_cophenetic_tanglegram.build_row verbatim (same untangling, colours, annotations);
only the multi-column composition + minimal-text framing are new.

REGISTER (inherits the explainer's preamble). Illustrative, single-clade, subset ρ^coph on the
dense per-phase cophenetics — NOT a gated test; the cohort claim stays the matched-strength gates.
Exemplars/clades are disclosed-selected. This is the taxonomy VOCABULARY slide, not evidence.

Writes: data/outputs/figures/talk/taxonomy_quad_tanglegram.{pdf,png}
"""
from __future__ import annotations

import sys
import time

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.config.paths import FIGURES_ROOT

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded"))
from fig_cophenetic_tanglegram import (                    # type: ignore  (DRY)
    load_cophs, best_clades, build_row, contact_names, ROW_META,
)

# fate -> (patient, band): the exact exemplars the class figure picked (== the current talk PNGs)
EXEMPLARS = [
    ("reset",       "Pat_10", "beta"),
    ("reorganized", "Pat_08", "delta"),
    ("anchor",      "Pat_08", "low_gamma"),
    ("trace",       "Pat_06", "delta"),
]
TASK = "task_test"
LO, HI = 20, 30                       # same clade-size window as fig_cophenetic_tanglegram_class
OUTDIR = FIGURES_ROOT / "talk"
TITLE_INK = "#2b2f36"

# Tighter rails + wider dendrograms than the per-fate default: give the tree structure more width
# and the cross-phase ribbons less (rail gap 1.0 → 0.74; ext 0.34 → 0.52). show_outer_rho is off,
# so we also trim the wide right pad the default reserved for that note (XPAD_R small).
RAIL_Q = (0.0, 0.74, 1.48)
EXT_Q = 0.52
EXT_MID_Q = 0.55 * EXT_Q
XPAD_L, XPAD_R = 0.10, 0.12


def main():
    t0 = time.perf_counter()
    use_lrg_style()

    # resolve each exemplar's clade exactly as the class figure does (deterministic: no RNG here —
    # the fleet-scan RNG only chose the patient/band, which we hard-code from the shipped PNGs)
    cols = []
    for fate, pat, band in EXEMPLARS:
        got = load_cophs(pat, band, TASK)
        if got is None:
            raise SystemExit(f"{pat} {band}: phase-N mismatch — cannot draw {fate}")
        C, Zpost = got
        S = best_clades(C, Zpost, TASK, LO, HI)[fate][1]
        cols.append((fate, pat, band, C, S))
        print(f"  {fate:11s} {pat} {band}  N={len(S)}", flush=True)

    fig = plt.figure(figsize=(21.6, 10.8))                     # 2:1 landscape — wide columns so the
    gs = GridSpec(1, len(cols), top=0.915, bottom=0.02, left=0.008, right=0.992,   # dendrograms
                  wspace=0.025, figure=fig)                                         # don't squash

    for j, (fate, pat, band, C, S) in enumerate(cols):
        ax = fig.add_subplot(gs[j])
        labels = [contact_names(pat, C["rest_post"].shape[0])[l] for l in S]
        phase_order = ("rest_pre", "rest_post", TASK) if fate == "reset" else None
        build_row(ax, C, S, TASK, is_top=True, labels=labels, phase_order=phase_order,
                  show_outer_rho=False, rail=RAIL_Q, ext=EXT_Q, ext_mid=EXT_MID_Q)
        ax.set_xlim(RAIL_Q[0] - EXT_Q - XPAD_L, RAIL_Q[2] + EXT_Q + XPAD_R)   # trim dead margin
        # category title, horizontal, centred above the phase labels of its column
        pos = ax.get_position()
        fig.text(0.5 * (pos.x0 + pos.x1), 0.965, ROW_META[fate][0], ha="center", va="top",
                 fontsize=16, fontweight="bold", color=TITLE_INK)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    pdf = OUTDIR / "taxonomy_quad_tanglegram.pdf"
    fig.savefig(pdf, bbox_inches="tight", transparent=True, pad_inches=0.05)   # vector primary
    png = OUTDIR / "taxonomy_quad_tanglegram.png"                              # Canva import (opt-in PNG)
    fig.savefig(png, dpi=300, bbox_inches="tight", transparent=True, pad_inches=0.05)
    plt.close(fig)
    print(f"wrote {pdf}\n      {png}   [{time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
