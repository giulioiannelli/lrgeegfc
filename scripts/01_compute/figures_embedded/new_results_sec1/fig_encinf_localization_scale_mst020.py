#!/usr/bin/env python3
r"""fig_encinf_localization_scale_mst020 -- the encoding/inference ANATOMY dissolves (mst@0.20).

HEAD: the whole-graph double dissociation -- encoding -> orbitofrontal cortex (q=.010-.040),
inference-specific -> cingulate (q=.030-.040) -- does NOT survive on the honest mst@0.20
backbone. Swept across all 16 diffusion scales under the matched-strength null (BH per band x
scale, the localization ARC, script 17), the two a-priori hotspots never clear:

  * ENCODING -> OFC: the beta OFC enrichment is flat and non-significant (best p=0.156, obs
    even positive but never separating from the strength null); NOTHING clears BH for encoding
    in any band/system/scale. Encoding delocalizes, exactly like the whole-task trace.
  * INFERENCE -> cingulate: the strongest surviving lean is a sub-threshold ALPHA trend at a
    coarse-mesoscale s~=16 (p=0.008, BH q=0.06, 7/8 patients) -- it MISSES BH, and it is alpha,
    not beta (beta cingulate p=0.23, n.s.). A directional hint, not an established home.

So R2's anatomy is a degenerate-dense-graph artifact (audit_174): removing the degeneracy
dissolves both hotspots, the same fate as the beta TRACE (fig_localization_scale_mst020). The
load-bearing encoding/inference result on mst@0.20 is the SCALE SIGNATURE (fig_encinf_tau_gate:
beta encoding continuous-multiscale, alpha not) and the band DISSOCIATION -- not a place. Shown
against the retired whole-graph q's (struck through) so the reader sees exactly what changed.

Reads : data/sparsified_arc/localization_arc_mst020/system_enrichment.csv
        (target in {encoding, inference}, band, s, unit, obs_med, n_pos, p, q; matched-strength
         R=200, BH per band x scale; script 17_localization_arc_mst020).
Writes: data/preprint/figures/new_results_sec1/fig_encinf_localization_scale_mst020.pdf
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import _common as C

C.use_lrg_style()

SRC = C.SA / "localization_arc_mst020" / "system_enrichment.csv"
BANDS = ["delta", "alpha", "beta", "low_gamma"]
# (target, a-priori hotspot, retired whole-graph q, panel title)
PANELS = [
    ("encoding", "OFC", "0.010", r"encoding $\rightarrow$ OFC"),
    ("inference", "cingulate", "0.030", r"inference-specific $\rightarrow$ cingulate"),
]
PCAP = 3.0


def hotspot_curve(df, target, unit, band):
    """(-log10 p, q) vs scale for one system/band/target, scale-sorted."""
    d = df[(df.target == target) & (df.band == band) & (df.unit == unit)].sort_values("s")
    s = d.s.to_numpy(float)
    lp = np.minimum(-np.log10(np.clip(d.p.to_numpy(float), 1e-3, 1.0)), PCAP)
    q = d.q.to_numpy(float)
    return s, lp, q


def main():
    df = pd.read_csv(SRC)
    p05 = -np.log10(0.05)

    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.2), sharey=True)
    for ax, (target, unit, wg_q, title) in zip(axes, PANELS):
        for band in BANDS:
            col = C.band_color(band)
            s, lp, q = hotspot_curve(df, target, unit, band)
            if s.size == 0:
                continue
            ax.plot(s, lp, "-o", color=col, lw=2.0, ms=4, alpha=0.9, label=C.BTeX[band])
            # mark BH-clearing points (none expected) + the strongest sub-threshold lean
            clr = q < 0.05
            if clr.any():
                ax.scatter(s[clr], lp[clr], s=150, facecolor=col, edgecolor="black",
                           linewidth=1.6, zorder=5)
            # callout the best near-miss on the a-priori hotspot
            ib = int(np.argmin(q))
            if q[ib] < 0.10:                                 # only annotate a genuine near-miss
                ax.annotate(rf"{C.BTeX[band]}: $p={10**(-lp[ib]):.3f}$, BH $q={q[ib]:.2f}$"
                            + ("  (misses BH)" if q[ib] >= 0.05 else ""),
                            (s[ib], lp[ib]), xytext=(6, 12), textcoords="offset points",
                            fontsize=8.5, color=col, fontweight="bold",
                            arrowprops=dict(arrowstyle="->", color=col, lw=1.0))
        ax.axhline(p05, color="0.5", lw=1.0, ls="--", zorder=1)
        ax.text(1.05, p05 + 0.03, r"nominal $p=0.05$", fontsize=8, color="0.5", va="bottom")
        ax.set_xscale("log")
        ax.set_xlabel(C.XLAB_S)
        ax.set_title(title, fontsize=13, fontweight="bold", loc="left")
        # retired whole-graph reference, struck through
        ax.text(0.97, 0.955, f"whole-graph: q={wg_q}", transform=ax.transAxes, ha="right",
                va="top", fontsize=9.5, color="0.5",
                bbox=dict(boxstyle="round,pad=0.2", fc="0.93", ec="none"))
        ax.plot([0.80, 0.965], [0.945, 0.945], transform=ax.transAxes, color="#c23",
                lw=1.4, zorder=6, clip_on=False)                  # strike-through
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(r"$-\log_{10}p_{\mathrm{MS}}$  (a-priori hotspot)")
    axes[0].set_ylim(0, PCAP + 0.2)
    axes[0].legend(frameon=False, fontsize=10, ncol=4, loc="upper left",
                   handletextpad=0.4, columnspacing=1.0)

    fig.text(0.5, -0.02, "no band/scale clears BH on either hotspot — the enc→OFC / inf→"
             "cingulate dissociation is a degenerate-graph artifact; only a sub-threshold "
             r"$\alpha$→cingulate lean ($q=0.06$) survives", ha="center", va="top",
             fontsize=10, color="0.25", fontstyle="italic")
    fig.tight_layout()

    C.FIGDIR.mkdir(parents=True, exist_ok=True)
    out = C.FIGDIR / "fig_encinf_localization_scale_mst020.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)

    print("fig_encinf_localization_scale_mst020 — enc/inf hotspots vs scale (mst@0.20):")
    for target, unit, _, _ in PANELS:
        print(f"  {target} -> {unit}:")
        for band in BANDS:
            s, lp, q = hotspot_curve(df, target, unit, band)
            if s.size == 0:
                continue
            ib = int(np.argmin(q))
            nclr = int((q < 0.05).sum())
            print(f"     {band:10s} best p={10**(-lp[ib]):.3f} q={q[ib]:.2f} @s={s[ib]:.1f}"
                  f"  clears {nclr}/16 scales")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
