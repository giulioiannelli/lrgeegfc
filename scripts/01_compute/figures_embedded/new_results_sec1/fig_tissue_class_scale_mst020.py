#!/usr/bin/env python3
r"""fig_tissue_class_scale_mst020 -- which tissue carries the beta trace? ALL of it, bar the
diseased core (mst@0.20).

HEAD: labelling every contact pair by the tissue it links and reading the beta cophenetic trace
across all 16 diffusion scales on the honest mst@0.20 backbone (matched-strength null, BH), the
held trace is a GRAY + WHITE + MIXED phenomenon that thins only over the seizure-onset core:

  (a) TISSUE. gray<->gray, gray<->white and white<->white coupling ALL clear the strength null
      at essentially every scale (14-16 / 16 scales, BH q~=0.015). The beta trace is tissue-
      NON-SPECIFIC -- it does not prefer gray over white; it is broadly distributed, the tissue
      face of the delocalization seen anatomically (fig_localization_scale_mst020).
  (b) SEIZURE-ONSET ZONE. non-SOZ<->non-SOZ pairs clear at 15/16 scales, but SOZ<->SOZ pairs
      clear at only 1/16 (q~=0.10) -- the beta trace THINS over the diseased core, echoing the
      delta SOZ-avoidance (fig_soz_scale_mst020) at the carrier band.

Band-specificity (console + caption): only beta clears; alpha hovers just under the null
everywhere (best BH q=0.06, a near-miss at every class), delta/low_gamma do not clear -- so this
is the beta trace's tissue signature, not a generic effect. Filled marker = clears BH (q<0.05) at
that scale; open = does not.

Reads : data/sparsified_arc/tissue_pairclass_mst020/cohort.csv
        (band, pair_class, s, obs_med, surr_med, n_pos, p_upper, q_upper; matched-strength R=200,
         BH; script 19_tissue_pairclass_mst020).
Writes: data/preprint/figures/new_results_sec1/fig_tissue_class_scale_mst020.pdf
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import _common as C

C.use_lrg_style()

SRC = C.SA / "tissue_pairclass_mst020" / "cohort.csv"
CARRIER = "beta"
# pair-class categories are TISSUE types (not frequency bands) -> a dedicated categorical
# palette here is correct; the band-palette rule governs band colours only.
TISSUE = [("gray_gray", r"gray $\leftrightarrow$ gray", "#3a6fb0"),
          ("cross_wm", r"gray $\leftrightarrow$ white", "#5fa06b"),
          ("wm_wm", r"white $\leftrightarrow$ white", "#b0763a")]
SOZ = [("nonepi_nonepi", r"non-SOZ $\leftrightarrow$ non-SOZ", "#2f9e52"),
       ("cross_epi", r"SOZ $\leftrightarrow$ cortex", "#d99a1f"),
       ("epi_epi", r"SOZ $\leftrightarrow$ SOZ", "#d1352b")]


def cls_curve(df, band, cls):
    d = df[(df.band == band) & (df.pair_class == cls)].sort_values("s")
    return (d.s.to_numpy(float), d.obs_med.to_numpy(float),
            d.surr_med.to_numpy(float), d.q_upper.to_numpy(float))


def draw_panel(ax, df, classes, title, tag):
    for cls, lab, col in classes:
        s, obs, surr, q = cls_curve(df, CARRIER, cls)
        if s.size == 0:
            continue
        ax.plot(s, obs, "-", color=col, lw=1.9, alpha=0.9, label=lab, zorder=3)
        clr = q < 0.05
        ax.scatter(s[clr], obs[clr], s=42, color=col, edgecolor="black", linewidth=1.0,
                   zorder=5)                                    # clears BH
        ax.scatter(s[~clr], obs[~clr], s=30, facecolor="white", edgecolor=col,
                   linewidth=1.3, zorder=4)                     # within null
        nclr = int(clr.sum())
        ax.annotate(f"{nclr}/16", (s[-1], obs[-1]), xytext=(5, 0),
                    textcoords="offset points", fontsize=8.5, color=col, va="center",
                    fontweight="bold")
    # matched-strength null band (median surrogate across classes, for context)
    s0, _, surr0, _ = cls_curve(df, CARRIER, classes[0][0])
    ax.plot(s0, surr0, color="0.55", lw=1.0, ls="--", zorder=2)
    ax.text(s0[0], surr0[0], "MS null", fontsize=8, color="0.5", va="bottom", ha="left")
    ax.set_xscale("log")
    ax.set_xlabel(C.XLAB_S)
    ax.set_ylabel(r"held trace  $\rho_{\mathrm{sym}}^{\mathrm{coph}}$  (cohort)")
    ax.set_ylim(-0.03, 0.40)
    ax.set_title(title, fontsize=13, fontweight="bold", loc="left")
    ax.text(0.012, 1.015, rf"$\mathbf{{{tag}}}$", transform=ax.transAxes, fontsize=16,
            va="bottom", ha="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=9.5, loc="lower right", handletextpad=0.5)
    ax.spines[["top", "right"]].set_visible(False)


def main():
    df = pd.read_csv(SRC)
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.2), sharey=True)
    draw_panel(axes[0], df, TISSUE, r"$\beta$ trace is tissue-non-specific", "a")
    draw_panel(axes[1], df, SOZ, r"$\beta$ trace thins over the seizure-onset core", "b")

    handles = [Line2D([0], [0], marker="o", ls="", mfc="0.35", mec="black", ms=8,
                      label="clears matched-strength BH ($q<0.05$)"),
               Line2D([0], [0], marker="o", ls="", mfc="white", mec="0.4", mew=1.3, ms=8,
                      label="within null")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.04),
               ncol=2, frameon=False, fontsize=9.5)
    fig.text(0.5, 0.955, r"only $\beta$ clears (right margin = scales cleared / 16); "
             r"$\alpha$ is a near-miss everywhere (BH $q=0.06$), $\delta/\gamma_{\rm low}$ null",
             ha="center", va="bottom", fontsize=10, color="0.25", fontstyle="italic")
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))

    C.FIGDIR.mkdir(parents=True, exist_ok=True)
    out = C.FIGDIR / "fig_tissue_class_scale_mst020.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)

    print("fig_tissue_class_scale_mst020 — beta trace by tissue/SOZ class (mst@0.20):")
    for cls, lab, _ in TISSUE + SOZ:
        s, obs, surr, q = cls_curve(df, CARRIER, cls)
        print(f"  {cls:14s} clears {int((q < 0.05).sum()):2d}/16  best q={np.nanmin(q):.3f}  "
              f"obs@s5.6={obs[5]:+.3f}")
    print("  band-selectivity (best q over all classes/scales):")
    for band in ["beta", "alpha", "delta", "low_gamma"]:
        qb = df[df.band == band].q_upper
        print(f"     {band:10s} best q={qb.min():.3f}  clears {(qb < 0.05).sum()} cells")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
