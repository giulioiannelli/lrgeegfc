#!/usr/bin/env python3
"""fig_trace_a_band_violin -- per-band held-trace distribution + the tau landscape.

Two panels on the mst@0.20 backbone:
  (a) at the reporting scale s~=5.6, a violin per band of the 10 per-patient
      rho_sym^coph, with the patients overplotted and the matched-strength null level
      (cohort-median p95) as a tick; cohort gate p annotated. beta and alpha sit above
      their null, the rest straddle it.
  (b) the band x scale gate landscape: -log10(cohort gate p) over all 16 scales, so the
      alpha-mesoscale vs beta-scale-broad shape is visible at a glance (delta/high_gamma
      only light up patchily / into the collapse tail).

Reads: data/sparsified_arc/ms_mst020/{per_patient_scale.csv, cohort_gate.csv}.
"""
from __future__ import annotations
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize

import _common as C

RNG = np.random.default_rng(0)


def violin_panel(ax, pp, gate):
    s0 = C.SGRID[C.I_REPORT]
    for k, band in enumerate(C.BANDS):
        x = pp[(pp.band == band) & np.isclose(pp.s, s0)]
        vals = x.obs_rho.dropna().values
        if len(vals) == 0:
            continue
        col = C.band_color(band)
        vp = ax.violinplot(vals, positions=[k], widths=0.72, showextrema=False)
        for b in vp["bodies"]:
            b.set_facecolor(col); b.set_alpha(0.28); b.set_edgecolor(col); b.set_lw(1.2)
        jit = (RNG.random(len(vals)) - 0.5) * 0.22
        ax.scatter(k + jit, vals, s=16, color=col, edgecolor="white", lw=0.4, zorder=3)
        ax.plot([k], [np.median(vals)], "_", color=col, ms=22, mew=2.4, zorder=4)
        # matched-strength null level (cohort-median p95 at this scale)
        p95 = np.nanmedian(x.surr_p95.values)
        ax.plot([k - 0.36, k + 0.36], [p95, p95], color="0.45", lw=1.3, ls="--", zorder=2)
        # cohort gate p annotation
        gp = gate[(gate.band == band) & np.isclose(gate.s, s0)]
        if not gp.empty:
            p = float(gp.gate_p.iloc[0])
            ax.text(k, ax.get_ylim()[1], f"p={p:.3f}",
                    ha="center", va="bottom", fontsize=8,
                    color=col if p < 0.05 else "0.5",
                    fontweight="bold" if p < 0.05 else "normal")
    ax.axhline(0, color="0.7", lw=0.7)
    ax.set_xticks(range(len(C.BANDS)))
    ax.set_xticklabels([C.BTeX[b] for b in C.BANDS], fontsize=13)
    ax.set_ylabel(C.YLAB_RHO)
    ax.set_title(rf"held trace by band  @ $s={s0:.1f}$   "
                 r"(dashed = matched-strength null $p_{95}$)",
                 fontsize=11, fontweight="bold", loc="left")


def gate_heatmap(ax, gate):
    s_vals = np.sort(gate.s.unique())
    M = np.full((len(C.BANDS), len(s_vals)), np.nan)
    for i, b in enumerate(C.BANDS):
        gb = gate[gate.band == b]
        for _, r in gb.iterrows():
            j = int(np.argmin(np.abs(s_vals - r.s)))
            M[i, j] = -np.log10(max(r.gate_p, 1e-4))
    im = ax.imshow(M, aspect="auto", cmap="magma", norm=Normalize(0, 3),
                   extent=[np.log10(s_vals.min()), np.log10(s_vals.max()),
                           len(C.BANDS) - 0.5, -0.5])
    ax.set_yticks(range(len(C.BANDS)))
    ax.set_yticklabels([C.BTeX[b] for b in C.BANDS], fontsize=12)
    ax.set_xlabel(r"$\log_{10} s = \log_{10}(\tau\lambda_{\max})$")
    ax.axvline(0.0, color="w", lw=1.0, ls=":")
    ax.axvline(np.log10(C.S_REPORT), color="c", lw=1.0, ls="--")
    ax.set_title(r"$\tau$ landscape: $-\log_{10}$ cohort gate $p$  "
                 r"(dotted $=\tau_{\min}$, dashed $=s_{\mathrm{report}}$)",
                 fontsize=11, fontweight="bold", loc="left")
    cb = ax.figure.colorbar(im, ax=ax, pad=0.015, fraction=0.045)
    cb.set_label(r"$-\log_{10}p$")
    cb.ax.axhline(-np.log10(0.05), color="c", lw=1.5)


def main():
    C.use_lrg_style()
    pp = pd.read_csv(C.MS / "per_patient_scale.csv")
    gate = pd.read_csv(C.MS / "cohort_gate.csv")

    fig = plt.figure(figsize=(12, 9))
    gs = GridSpec(2, 1, height_ratios=[1.35, 1.0], hspace=0.32, figure=fig)
    violin_panel(fig.add_subplot(gs[0]), pp, gate)
    gate_heatmap(fig.add_subplot(gs[1]), gate)

    C.FIGDIR.mkdir(parents=True, exist_ok=True)
    out = C.FIGDIR / "fig_trace_a_band_violin.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


if __name__ == "__main__":
    main()
