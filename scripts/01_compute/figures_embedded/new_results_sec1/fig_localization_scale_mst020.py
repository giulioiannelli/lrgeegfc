#!/usr/bin/env python3
r"""fig_localization_scale_mst020 -- the held trace is strong but PLACELESS (mst@0.20).

HEAD: on the honest **mst@0.20** backbone, swept across diffusion scale, the cophenetic
task->rest trace has NO anatomical home. For every band that carries a trace (delta, alpha,
beta, low_gamma) and every a-priori cortical system, the per-system enrichment of the
per-pair rho_sym concordance is tested against the matched-strength null and BH-corrected
within band. Read across the three scales s in {1.0, 2.83, 5.65}:

  * NOTHING clears BH at any scale, in any band. The strongest single hits are only TRENDS
    that fail correction: low_gamma->PFC (p=0.023, q=0.19 @ s=1), alpha->insula (p=0.031,
    q=0.25 @ s=2.83), delta->cingulate (p=0.055 @ s=5.65).
  * The beta trace -- the flagship -- is the clearest non-result: its orbitofrontal
    enrichment is p=0.50 -> 0.78 -> 0.59 across scale and the sign FLIPS +91k -> +130k ->
    -322k. There is no OFC concentration on the recovered scheme.

This SUPERSEDES the whole-graph "beta -> OFC" localization (q=0.009-0.013): that was read on
the degenerate dense graph, which is not a baseline (audit_174). Removing the degeneracy
dissolves the anatomical home -- the trace delocalizes, consistent with the ~88% node
co-movement of a delocalised reorganization. The localization axis is a WEAK strength on
mst@0.20; the load-bearing results are the band-selective cohort gate (fig_tau_role) and its
scale shape, not a place.

The single BH-significant spatial result on this backbone is a DISEASE contrast, not a
cognitive home -- the delta trace AVOIDS the seizure-onset zone (companion
fig_soz_scale_mst020).

Reads : data/sparsified_arc/localization_mst020_s{01.0,02.8,05.7}/system_enrichment.csv
        (band, unit, n, obs_med, n_pos, p, q; matched-strength R=200, BH per band).
Writes: data/preprint/figures/new_results_sec1/fig_localization_scale_mst020.pdf
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import _common as C

C.use_lrg_style()

SA = C.SA
SCALES = [("s01.0", 1.0), ("s02.8", 2.83), ("s05.7", 5.65)]      # dir suffix, s value
BANDS = ["delta", "alpha", "beta", "low_gamma"]                   # bands with a trace + on disk
# limbic/paralimbic core first (OFC = the retired whole-graph "home")
SYS_ORDER = ["OFC", "cingulate", "MTL", "insula", "lateral_temporal",
             "PFC", "parietal", "sensorimotor", "occipital", "subcortical_other"]
# the honest trend cells (fail BH) worth naming on the figure
TRENDS = {("low_gamma", "PFC"): "s=1", ("alpha", "insula"): "s=2.83",
          ("delta", "cingulate"): "s=5.65"}
PCAP = 2.2                                                        # -log10 p display cap
# fade fine->meso so the scale sweep reads as one connected track
SCALE_ALPHA = {1.0: 0.32, 2.83: 0.60, 5.65: 1.0}
SCALE_MS = {1.0: 34, 2.83: 52, 5.65: 82}


def load_scale(suffix):
    df = pd.read_csv(SA / f"localization_mst020_{suffix}" / "system_enrichment.csv")
    return {(r.band, r.unit): (float(r.p), float(r.q), float(r.obs_med), int(r.n_pos), int(r.n))
            for r in df.itertuples()}


def signed_logp(p, obs_med):
    """signed -log10 p: enrichment (obs>0) to the right, avoidance to the left, capped."""
    return np.sign(obs_med if obs_med != 0 else 1.0) * min(-np.log10(max(p, 1e-3)), PCAP)


def main():
    data = {s: load_scale(suf) for suf, s in SCALES}
    systems = [s for s in SYS_ORDER
               if any((b, s) in data[5.65] for b in BANDS)]
    y0 = np.arange(len(systems))[::-1]
    ypos = dict(zip(systems, y0))

    fig, axes = plt.subplots(1, len(BANDS), figsize=(15.0, 5.6), sharey=True)
    p_bh = -np.log10(0.05)                                        # nominal .05 reference

    for ax, band in zip(axes, BANDS):
        col = C.band_color(band)
        ax.axvspan(-p_bh, p_bh, color="0.5", alpha=0.07, zorder=0)   # sub-.05 "no-signal" band
        ax.axvline(0, color="0.55", lw=0.9, zorder=1)
        for sys in systems:
            y = ypos[sys]
            track = []
            for suf, s in SCALES:
                rec = data[s].get((band, sys))
                if rec is None:
                    continue
                p, q, obs, npos, n = rec
                x = signed_logp(p, obs)
                track.append((x, y, s, q))
            if not track:
                continue
            xs = [t[0] for t in track]
            ax.plot(xs, [y] * len(xs), color=col, lw=1.0, alpha=0.35, zorder=2)  # scale track
            for x, yy, s, q in track:
                clears = q < 0.05
                ax.scatter(x, yy, s=SCALE_MS[s], color=col if x >= 0 else "0.5",
                           edgecolor="black" if clears else "white",
                           linewidth=1.3 if clears else 0.5,
                           alpha=SCALE_ALPHA[s], zorder=4 if s == 5.65 else 3)
            # name the honest trend cell
            if (band, sys) in TRENDS:
                xm = max(xs, key=abs)
                ax.annotate(f"trend ({TRENDS[(band, sys)]})\nfails BH", (xm, y),
                            xytext=(6 if xm >= 0 else -6, 10), textcoords="offset points",
                            fontsize=7.0, color="0.30", ha="left" if xm >= 0 else "right",
                            va="bottom")
        # beta OFC sign-flip callout (the flagship non-result)
        if band == "beta":
            ax.annotate(r"OFC sign flips $+\!\to\!-$" + "\n(no home)", (0, ypos["OFC"]),
                        xytext=(0.30, 0.90), textcoords="axes fraction", fontsize=8.0,
                        color="0.20", ha="center", fontweight="bold",
                        arrowprops=dict(arrowstyle="->", color="0.4", lw=1.0))
        ax.set_title(rf"{C.BTeX[band]}", fontsize=16, fontweight="bold", color=col,
                     loc="center", pad=6)
        ax.set_xlim(-PCAP - 0.35, PCAP + 0.35)
        ax.set_xlabel(r"$\leftarrow$ avoids   signed $-\log_{10}p_{\mathrm{MS}}$   enriched $\rightarrow$",
                      fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)

    axes[0].set_yticks(y0)
    axes[0].set_yticklabels(systems, fontsize=10.5)
    for lab in axes[0].get_yticklabels():
        if lab.get_text() in ("OFC", "cingulate"):
            lab.set_fontweight("bold")
    axes[0].set_ylim(y0.min() - 0.7, y0.max() + 0.7)

    # figure-level: scale legend + the empty-sky verdict
    handles = [Line2D([0], [0], marker="o", ls="", mfc="0.35", mec="white",
                      ms=np.sqrt(SCALE_MS[s]) * 1.05, alpha=SCALE_ALPHA[s],
                      label=rf"$s={s:g}$" + ("  (reporting)" if s == 5.65 else ""))
               for _, s in SCALES]
    handles.append(Line2D([0], [0], marker="o", ls="", mfc="0.35", mec="black", mew=1.3,
                          ms=8, label=r"BH $q<0.05$ (none occur)"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.03),
               ncol=4, frameon=False, fontsize=9.5, handletextpad=0.4, columnspacing=1.6)
    fig.text(0.5, 0.965, "no system clears BH at any scale — the mst@0.20 trace is "
             "delocalized (whole-graph β→OFC was a degenerate-graph artifact)",
             ha="center", va="top", fontsize=11, color="0.15", fontstyle="italic")
    fig.tight_layout(rect=(0, 0.03, 1, 0.94))

    C.FIGDIR.mkdir(parents=True, exist_ok=True)
    out = C.FIGDIR / "fig_localization_scale_mst020.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)

    print("fig_localization_scale_mst020 — system localization vs scale (mst@0.20):")
    for band in BANDS:
        print(f"  {band}:")
        for _, s in SCALES:
            recs = [(sys, *data[s][(band, sys)]) for sys in systems if (band, sys) in data[s]]
            best = min(recs, key=lambda r: r[1])            # min p
            nclr = sum(1 for r in recs if r[2] < 0.05)      # q<.05 count
            print(f"     s={s:<5g} best={best[0]:16s} p={best[1]:.3f} q={best[2]:.2f}  "
                  f"| BH q<.05 systems: {nclr}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
