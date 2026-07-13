#!/usr/bin/env python3
"""fig_tau_dendrogram_morph -- the ONE tau-morph (multiscale made concrete).

One exemplar network's diffusion cophenetic dendrogram read at four scales
``s = tau*lambda_max in {1, 5.6, 30, 90}``. As ``s`` grows the heat kernel integrates
longer paths, so the hierarchy coarsens: many tight fine clusters (s=1, tau_min) fold
into a few broad communities (s=90). This is what "multiscale" means -- one nested
dendrogram whose cophenetic distance folds every scale into one number; sweeping ``s``
just reads it at different depths.

Single exemplar, mst@0.20 backbone. Each panel is drawn in its own average-linkage
order; the coarsening is read from the merge heights, the community colouring and the
community count at a fixed relative cut.
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

import _common as C
from lrg_eegfc.utils.fc.heat_multiscale import linkage_at_scale

EXEMPLAR = ("Pat_06", "rest_post", "beta")   # super-responder, flagship band
SCALES = [1.0, C.S_REPORT, 30.0, 90.0]


def main():
    C.use_lrg_style()
    pat, phase, band = EXEMPLAR
    ev, V = C.eig_backbone(C.load_phase(pat, phase, band))
    col = C.band_color(band)

    n = V.shape[0]
    fig, axes = plt.subplots(1, len(SCALES), figsize=(16, 4.4))
    for ax, s in zip(axes, SCALES):
        Z = linkage_at_scale(ev, V, s)
        heights = np.sort(Z[:, 2])
        # cophenetic heights span many decades -> cut at the geometric (log) midpoint,
        # so the community count is meaningful and tracks the coarsening across scales.
        lo, hi = np.log10(heights[0]), np.log10(heights[-1])
        thr = 10.0 ** (lo + 0.5 * (hi - lo))
        # consistent look: fine (below-cut) links in the band colour, coarse backbone grey
        dendrogram(Z, ax=ax, no_labels=True,
                   link_color_func=lambda i: col if Z[i - n, 2] <= thr else "0.7")
        ax.set_yscale("log")
        ax.set_ylim(heights[0] * 0.8, heights[-1] * 1.05)
        ax.set_xticks([])
        # height dynamic range (decades) is the clean, monotonic depth signal: a rich
        # multiscale hierarchy at tau_min flattens toward a single-scale near-star at s=90.
        tag = "  collapsing" if (hi - lo) < 2.0 else ""
        ax.set_title(rf"$s={s:.1f}$   ({hi - lo:.1f} height decades){tag}",
                     fontsize=12, fontweight="bold", loc="left",
                     color=col if abs(s - C.S_REPORT) < 1e-6 else "0.15")
        ax.set_xlabel("contacts")
    axes[0].set_ylabel(r"cophenetic height $1/\hat\rho$")
    # mark the reporting-scale panel
    idx = int(np.argmin([abs(s - C.S_REPORT) for s in SCALES]))
    for spine in axes[idx].spines.values():
        spine.set_edgecolor(col)
        spine.set_linewidth(2.0)

    C.FIGDIR.mkdir(parents=True, exist_ok=True)
    out = C.FIGDIR / "fig_tau_dendrogram_morph.pdf"
    fig.tight_layout()
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}  ({pat} {band} {phase}, N={V.shape[0]}, "
          f"scales {SCALES})")


if __name__ == "__main__":
    main()
