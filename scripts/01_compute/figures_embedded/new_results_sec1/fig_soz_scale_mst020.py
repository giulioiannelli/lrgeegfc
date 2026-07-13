#!/usr/bin/env python3
r"""fig_soz_scale_mst020 -- does the seizure-onset zone carry its own scale? (mst@0.20).

HEAD: the user's question -- does the SOZ have a band/scale of its own in the trace, and is
it alpha? Swept across diffusion scale on the honest mst@0.20 backbone, the answer is a
DISSOCIATION, and it is not alpha:

  (a) TRACE vs SOZ (per-pair rho_sym concordance, matched-strength null, BH per band, three
      localization scales s in {1, 2.83, 5.65}). The cognitive carrier bands alpha and beta
      lean weakly INTO the SOZ (SOZ-SOZ pairs slightly enriched) but never clear -- alpha
      p=0.14->0.22, beta p=0.35->0.39. The one band with a SIGNIFICANT, scale-growing SOZ
      relationship is delta, and it points the OTHER way: the delta trace AVOIDS the SOZ,
      concentrating in non-SOZ tissue (p 0.39 @ s=1 -> BH q=0.006 @ meso, 9/10 patients).

  (b) DISEASE marker (seeded heat-kernel SOZ-detection AUC vs scale, epi_arc). The same
      slow bands that dissociate in (a) are the strongest structural SOZ markers -- delta
      (peak AUC 0.85 @ s~=3) and low_gamma (0.84 @ s=1) -- while alpha is the WEAKEST marker
      (0.67). The SOZ's "own scale" is a mesoscale delta phenomenon on both axes: delta reads
      the diseased core best AND its cognitive trace steers around it.

So the SOZ does have a privileged scale, but it is a DELTA / mesoscale disease signature, not
an alpha cognitive trace. This replaces the whole-graph "alpha recruits the SOZ (rho=+0.41,
p=0.005)" reading, which does not reach significance once the degenerate dense graph is
removed (audit_174); the honest mst@0.20 result is delta SOZ-avoidance.

Reads : data/sparsified_arc/localization_mst020_s{01.0,02.8,05.7}/soz_enrichment.csv
        data/sparsified_arc/epi_arc_mst020/<band>/<Pat>.npz  (s, auc; cohort median)
Writes: data/preprint/figures/new_results_sec1/fig_soz_scale_mst020.pdf
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import _common as C

C.use_lrg_style()

SCALES = [("s01.0", 1.0), ("s02.8", 2.83), ("s05.7", 5.65)]
BANDS = ["delta", "alpha", "beta", "low_gamma"]
PCAP = 3.0


def soz_pref(band, suffix):
    """One SOZ-preference stat per (band, scale): sign +1 = trace leans INTO SOZ,
    -1 = AVOIDS SOZ (concentrates in non-SOZ). |value| = -log10 p of the leaning
    direction; also returns BH q of that direction."""
    df = pd.read_csv(C.SA / f"localization_mst020_{suffix}" / "soz_enrichment.csv")
    d = df[df.band == band]
    soz = d[d.unit == "soz"].iloc[0]
    non = d[d.unit == "non_soz"].iloc[0]
    if float(soz.p) <= float(non.p):                 # leans into SOZ
        sign, p, q = +1.0, float(soz.p), float(soz.q)
    else:                                            # avoids SOZ
        sign, p, q = -1.0, float(non.p), float(non.q)
    return sign * min(-np.log10(max(p, 1e-3)), PCAP), q


def auc_vs_scale(band):
    dd = C.SA / "epi_arc_mst020" / band
    aucs, svals = [], None
    for p in C.COHORT:
        fp = dd / f"{p}.npz"
        if fp.exists():
            z = np.load(fp)
            aucs.append(z["auc"])
            svals = z["s"]
    A = np.vstack(aucs)
    return svals, np.nanmedian(A, axis=0), A.shape[0]


def panel_a(ax):
    sx = [s for _, s in SCALES]
    ax.axhspan(-np.log10(0.05), np.log10(0.05), color="0.5", alpha=0.08, zorder=0)  # ns band
    ax.axhline(0, color="0.55", lw=0.9, zorder=1)
    for band in BANDS:
        col = C.band_color(band)
        ys, qs = zip(*[soz_pref(band, suf) for suf, _ in SCALES])
        ax.plot(sx, ys, "-", color=col, lw=1.8, alpha=0.9, zorder=3)
        for x, y, q in zip(sx, ys, qs):
            clears = q < 0.05
            ax.scatter(x, y, s=150 if clears else 70, color=col,
                       edgecolor="black" if clears else "white",
                       linewidth=1.6 if clears else 0.6, zorder=4)
            if clears:
                ax.annotate(rf"$q={q:.3f}$" + "\n" + r"$\delta$ avoids SOZ", (x, y),
                            xytext=(-8, -4), textcoords="offset points", ha="right",
                            va="center", fontsize=8.5, color=col, fontweight="bold")
    ax.set_xscale("log")
    ax.set_xticks(sx)
    ax.set_xticklabels([f"{s:g}" for s in sx])
    ax.set_xlabel(C.XLAB_S)
    ax.set_ylabel(r"$\leftarrow$ avoids SOZ   signed $-\log_{10}p$   into SOZ $\rightarrow$")
    ax.set_ylim(-PCAP - 0.3, PCAP * 0.55)
    ax.set_title(r"$\mathbf{a}$   trace vs seizure-onset zone", fontsize=13,
                 fontweight="bold", loc="left")
    ax.spines[["top", "right"]].set_visible(False)


def panel_b(ax):
    for band in BANDS:
        col = C.band_color(band)
        s, med, n = auc_vs_scale(band)
        ax.plot(s, med, "-", color=col, lw=2.0, alpha=0.9,
                label=rf"{C.BTeX[band]}  (peak {med.max():.2f})")
        ax.scatter(s[np.nanargmax(med)], med.max(), s=60, color=col,
                   edgecolor="white", lw=0.7, zorder=4)
    ax.axhline(0.5, color="0.55", lw=0.9, ls="--", zorder=1)
    ax.text(1.05, 0.505, "chance", fontsize=8, color="0.5", va="bottom")
    ax.set_xscale("log")
    ax.set_xlabel(C.XLAB_S)
    ax.set_ylabel("SOZ-detection AUC (cohort median)")
    ax.set_ylim(0.45, 0.92)
    ax.set_title(r"$\mathbf{b}$   SOZ as a structural disease marker", fontsize=13,
                 fontweight="bold", loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower right", frameon=False, fontsize=9.5, ncol=2,
              handletextpad=0.5, columnspacing=1.2)


def main():
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.2))
    panel_a(axes[0])
    panel_b(axes[1])
    fig.text(0.5, -0.01, r"the SOZ's privileged scale is a mesoscale $\delta$ DISEASE "
             r"signature (avoidance + strongest marker), not an $\alpha$ cognitive trace",
             ha="center", va="top", fontsize=10, color="0.25", fontstyle="italic")
    fig.tight_layout()

    C.FIGDIR.mkdir(parents=True, exist_ok=True)
    out = C.FIGDIR / "fig_soz_scale_mst020.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)

    print("fig_soz_scale_mst020 — SOZ vs scale (mst@0.20):")
    print("  (a) trace SOZ-preference (sign: + into / - avoids), q of leaning direction:")
    for band in BANDS:
        cells = [(s, *soz_pref(band, suf)) for suf, s in SCALES]
        for s, v, q in cells:
            print(f"     {band:10s} s={s:<5g} signed_logp={v:+.2f} q={q:.3f}"
                  f"{'  ***BH***' if q < 0.05 else ''}")
    print("  (b) disease-marker AUC peaks:")
    for band in BANDS:
        s, med, n = auc_vs_scale(band)
        print(f"     {band:10s} peak AUC={med.max():.2f} @ s={s[np.nanargmax(med)]:.1f} (n={n})")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
