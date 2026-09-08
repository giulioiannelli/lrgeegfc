#!/usr/bin/env python3
r"""talk_fig_grassmann_contrast -- the trace is a MULTISCALE reorganization: the standard toolkit
can't resolve it, spectral clustering is MISLED, only the multiscale hierarchy names the carriers.

A dot-matrix of the whole read-out ladder (rows) x band (cols), each cell gated on the SAME
mst@0.20 backbone (Grassmann on the whole graph) against the SAME matched-strength null:

  STANDARD NETWORK DESCRIPTORS (the toolkit we tested):
    raw FC edges  -> beta (a convergence carrier); alpha/high_g fragile -- pairwise, non-selective
    clustering coef. / effective resistance -> blind (nothing representative)
    graph geodesic -> low_gamma ONLY -- a GLOBAL path metric picking up a global structure
    (node strength is NOT shown: the matched-strength null fixes it by construction, so it is
     the control variable -- testing it against that null is circular.)
  FIELD-STANDARD SPECTRAL CLUSTERING:
    Grassmann chordal (whole-graph Laplacian eigen-subspaces, one fixed dim k)
      -> beta (agrees) + low_gamma; MISSES alpha
  MULTISCALE (ours):
    cophenetic (diffusion across ALL scales) -> alpha + beta, and NOTHING else

The point (read the low_gamma column): low_gamma is flagged ONLY by the two GLOBAL methods --
graph geodesic and Grassmann -- and NEVER by the multiscale read. So low_gamma is a whole-graph /
global spectral shift, NOT a multiscale reorganization: single-scale spectral clustering conflates
the two. It also MISSES alpha, which lives at a scale a fixed subspace cannot isolate. Only the
multiscale hierarchy resolves the true carriers (alpha + beta) and is not fooled by the global
low_gamma. The Grassmann formula lives in the slide text.

Verdicts computed FRESH:
  standard + cophenetic : data/sparsified_arc/controls_ladder_apples/per_cell.csv (representativeness)
  Grassmann             : data/audit/grassmann_cluster_extent/cohort_summary.csv (cluster-extent + LOO)
Writes: data/outputs/figures/talk/fig_grassmann_contrast.{png,pdf}  (transparent, dark, Canva-ready)
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[1] /
                       "01_compute" / "figures_embedded" / "new_results_sec1"))
import _common as C  # noqa: E402

C.use_lrg_style()
INK = "white"
plt.rcParams.update({
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
    "axes.titlecolor": INK, "xtick.color": INK, "ytick.color": INK,
})
SEP = "0.6"
NONCLEAR = "white"      # non-clearing cell -> white ring (was grey); clearing -> band colour
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
OUT = C.ROOT / "data" / "outputs" / "figures" / "talk" / "fig_grassmann_contrast"


def _grass_marks():
    g = pd.read_csv(C.ROOT / "data" / "audit" / "grassmann_cluster_extent" / "cohort_summary.csv")
    out = {}
    for b in BANDS:
        r = g[g.band == b].iloc[0]
        strong = str(r.verdict_cluster_extent) == "strong"
        loo_ok = float(r.cluster_p_mass_loo_max) < 0.05
        out[b] = "fill" if (strong and loo_ok) else ("ring" if strong else "open")
    return out


def _backbone_marks(desc):
    """Representativeness verdict (rep/fragile/null) per band for a backbone descriptor."""
    df = pd.read_csv(C.LADDER / "per_cell.csv")
    out = {}
    for b in BANDS:
        lvl = C.representativeness_gate(df, desc, b, "T_test")["level"]
        out[b] = "fill" if lvl == "rep" else ("ring" if lvl == "fragile" else "open")
    return out


def _cell(ax, j, y, band, mark):
    col = C.band_color(band)
    if mark == "fill":                                     # clears -> filled band colour
        ax.scatter(j, y, s=560, facecolor=col, edgecolor="white", lw=1.2, zorder=4)
    elif mark == "ring":                                   # sig-but-LOO-fragile -> band ring
        ax.scatter(j, y, s=560, facecolor="none", edgecolor=col, lw=3.0, zorder=3)
    else:                                                  # non-clearing -> white ring
        ax.scatter(j, y, s=470, facecolor="none", edgecolor=NONCLEAR, lw=1.7, zorder=2)


def main():
    # the full read-out ladder: standard toolkit -> field-standard spectral clustering -> ours
    # NB: node strength is intentionally NOT a row -- the matched-strength null
    # preserves the node-strength sequence by construction, so testing it against
    # that null is circular (it can never clear). It is the control, not a probe.
    std = [("raw_fc",     "raw FC edges\n(pairwise)"),
           ("clustering", "clustering coef.\n(local)"),
           ("geodesic",   "graph geodesic\n(paths · global)"),
           ("resistance", "effective resistance\n(spectral · global)")]
    grass = _grass_marks()
    coph = _backbone_marks("coph_meso")

    n = len(std) + 2                       # + Grassmann + cophenetic
    rows, y = [], n - 1                     # raw at top (y=n-1), cophenetic at bottom (y=0)
    for desc, lab in std:
        rows.append((lab, _backbone_marks(desc), y)); y -= 1
    rows.append(("Grassmann chordal\n(spectral clustering · global)", grass, y)); y -= 1
    rows.append(("cophenetic\n(multiscale · ours)", coph, y))

    fig = plt.figure(figsize=(13.4, 8.8))
    ax = fig.add_axes([0.30, 0.11, 0.66, 0.72])
    jb, ja, jg = (BANDS.index(x) for x in ("beta", "alpha", "low_gamma"))
    ax.axvspan(jb - 0.45, jb + 0.45, color="0.85", alpha=0.09, zorder=0)                  # β convergence
    ax.axvspan(jg - 0.45, jg + 0.45, color=C.band_color("low_gamma"), alpha=0.07, zorder=0)  # γ_low global

    for _, marks, yy in rows:
        for j, b in enumerate(BANDS):
            _cell(ax, j, yy, b, marks[b])

    # tier separators: standard toolkit | spectral clustering | multiscale (ours)
    ax.axhline(1.5, color="0.5", lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.axhline(0.5, color="0.5", lw=1.0, ls=(0, (4, 3)), zorder=1)

    ax.set_xlim(-0.6, len(BANDS) - 0.4)
    ax.set_ylim(-0.7, n - 0.2)                              # tight top: nothing above the dots
    ax.set_xticks(range(len(BANDS)))
    xtl = ax.set_xticklabels([C.BTeX[b] for b in BANDS], fontsize=33)
    for t, b in zip(xtl, BANDS):                            # colour each band name
        t.set_color(C.band_color(b))

    # interpretive callouts relocated BELOW the band names (top stays clean)
    tr = ax.get_xaxis_transform()                          # x in data, y in axes frac
    ax.text(ja, -0.14, "multiscale-\nonly", transform=tr, ha="center", va="top",
            fontsize=16, color=C.band_color("alpha"), clip_on=False)
    ax.text(jg, -0.14, "global —\nnot multiscale", transform=tr, ha="center", va="top",
            fontsize=16, color=C.band_color("low_gamma"), clip_on=False)
    ax.set_yticks([r[2] for r in rows])
    ax.set_yticklabels([r[0] for r in rows], fontsize=17)
    for t in ax.get_yticklabels():
        if "cophenetic" in t.get_text():
            t.set_fontweight("bold")                    # ours
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)

    leg = [Line2D([], [], marker="o", ls="none", ms=17, mfc="0.8", mec="white",
                  label="robust trace"),
           Line2D([], [], marker="o", ls="none", ms=17, mfc="none", mec="0.8", mew=3.0,
                  label="sig. but LOO-fragile"),
           Line2D([], [], marker="o", ls="none", ms=16, mfc="none", mec=NONCLEAR,
                  label="no trace")]
    ax.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.5, -0.30), ncol=3,
              frameon=False, fontsize=16, handletextpad=0.3, columnspacing=1.8)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}.{ext}", transparent=True, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"wrote {OUT}.png + .pdf")
    for lab, marks, _ in rows:
        print(f"  {lab.splitlines()[0]:26s}", marks)


if __name__ == "__main__":
    main()
