#!/usr/bin/env python3
"""fig_encinf_tau_gate -- encoding vs inference, resolved across the diffusion scale.

The four-phase decomposition (rest -> learn -> test -> rest) read through the mst@0.20
cophenetic hierarchy, tau-resolved. One panel per functional; each plots the cohort
matched-strength gate -log10(p) against s = tau*lambda_max for all six bands, with the
p<0.05 line marked. The dissociations the sweep exposes:

  - encoding  T_learn : BETA is continuous-multiscale (clears at every scale); ALPHA is
    NOT multiscale (a narrow fine-scale flicker only). Encoding is a real trace, band-
    dissociated by scale -- not the strength artifact an earlier read claimed.
  - held trace T_test : beta scale-broad, alpha mesoscale (as in fig_tau_role).
  - inference T_infspec (raw) : beta-only, fine scales.
  - inference T_infspec_pe (encoding partialled out) : delta/alpha/beta.

Reads: data/sparsified_arc/enc_inf_arc_mst020/{cohort_gate_vs_s.csv, band_nature.csv}.
"""
from __future__ import annotations
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import _common as C

FUNCS = [("T_learn",        "encoding  ($D_{\\mathrm{learn}}-D_{\\mathrm{pre}}$)",
          r"$\beta$ multiscale, $\alpha$ NOT"),
         ("T_test",         "held trace  ($D_{\\mathrm{test}}-D_{\\mathrm{pre}}$)",
          r"$\beta$ scale-broad, $\alpha$ mesoscale"),
         ("T_infspec",      "inference, raw  ($D_{\\mathrm{test}}-D_{\\mathrm{learn}}$)",
          r"$\beta$-only, fine"),
         ("T_infspec_pe",   "inference, enc-partialled",
          r"$\delta/\alpha/\beta$")]
P05 = -np.log10(0.05)


def _emph(band):
    """beta/alpha are the claimed bands -> heavier; the rest muted."""
    return (2.6, 1.0) if band in ("alpha", "beta") else (1.4, 0.5)


def main():
    C.use_lrg_style()
    g = pd.read_csv(C.ENCINF / "cohort_gate_vs_s.csv")

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    for ax, (fn, title, tag) in zip(axes.ravel(), FUNCS):
        gf = g[g.functional == fn]
        for band in C.BANDS:
            x = gf[gf.band == band].sort_values("s")
            if x.empty:
                continue
            y = -np.log10(np.clip(x.gate_p.values, 1e-4, 1.0))
            lw, al = _emph(band)
            ax.plot(x.s, y, "-o", color=C.band_color(band), lw=lw, ms=3.2,
                    alpha=al, label=C.BTeX[band])
        ax.axhline(P05, color="0.5", lw=1.0, ls="--")
        ax.text(1.02, P05, r"$p{=}0.05$", color="0.4", fontsize=8, va="bottom")
        ax.axvline(1.0, color="0.85", lw=0.7, ls=":")
        ax.set_xscale("log")
        ax.set_ylim(0, 3.2)
        ax.set_xlabel(C.XLAB_S)
        ax.set_ylabel(r"cohort gate  $-\log_{10}p$")
        ax.set_title(f"{title}\n[{tag}]", fontsize=11, fontweight="bold", loc="left")

    handles = [Line2D([0], [0], color=C.band_color(b), lw=2.4, marker="o", ms=4,
                      label=C.BTeX[b]) for b in C.BANDS]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=len(handles), frameon=False, fontsize=11)
    fig.tight_layout(rect=(0, 0.02, 1, 1))

    C.FIGDIR.mkdir(parents=True, exist_ok=True)
    out = C.FIGDIR / "fig_encinf_tau_gate.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


if __name__ == "__main__":
    main()
