#!/usr/bin/env python3
r"""fig_perpatient_bands_forest -- per-patient own-null trace forest, tau-aware (mst@0.20).

One panel per band. Each of the 10 patients is a row read at the reporting scale
s~=5.6: a stem to its rho_sym^coph, a dot (filled band-colour = clears its OWN matched-
strength null p95; open = within; red = anti), the null p95 as a tick, and -- the
tau-aware addition -- a right-margin count "k/16" of how many diffusion scales that
patient clears, so a patient that traces broadly across scales is distinguished from one
that only flickers at a single scale. The per-band corner glyph marks the cohort verdict
at s~=5.6.

Reads: data/sparsified_arc/ms_mst020/{per_patient_scale.csv, cohort_gate.csv}.
"""
from __future__ import annotations
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import _common as C
from lrg_eegfc.visuals.styles import band_marker

C_ANTI = "#d1352b"


def main():
    C.use_lrg_style()
    pp = pd.read_csv(C.MS / "per_patient_scale.csv")
    gate = pd.read_csv(C.MS / "cohort_gate.csv")
    s0 = C.SGRID[C.I_REPORT]
    ymap = {p: len(C.COHORT) - 1 - i for i, p in enumerate(C.COHORT)}
    # scales cleared per (patient, band): per-patient own-null p<0.05 over the 16 scales
    cleared = (pp.assign(hit=pp.p < 0.05)
                 .groupby(["patient", "band"]).hit.sum().astype(int))

    fig, axes = plt.subplots(2, 3, figsize=(14, 8.2), sharex=True, sharey=True)
    for ax, band in zip(axes.ravel(), C.BANDS):
        col = C.band_color(band)
        mk = band_marker(band)
        x = pp[(pp.band == band) & np.isclose(pp.s, s0)].set_index("patient")
        ax.axvline(0.0, color="0.55", lw=0.7, ls="--", zorder=1)
        for pat in C.COHORT:
            if pat not in x.index:
                continue
            r = x.loc[pat]; y = ymap[pat]
            clears = r.p < 0.05
            sign_col = col if r.obs_rho >= 0 else C_ANTI
            ax.plot([r.surr_p95], [y], "|", color="0.5", ms=11, mew=1.8, zorder=2)
            ax.plot([0.0, r.obs_rho], [y, y], color=sign_col, lw=1.0, zorder=3)
            if clears:
                ax.scatter(r.obs_rho, y, s=26, marker=mk, color=col,
                           edgecolor="white", lw=0.5, zorder=4)
            else:
                ax.scatter(r.obs_rho, y, s=26, marker=mk, facecolor="white",
                           edgecolor=sign_col, lw=1.0, zorder=4)
            # tau-aware count: scales cleared /16
            k = int(cleared.get((pat, band), 0))
            ax.text(0.99, y, f"{k}/16", transform=ax.get_yaxis_transform(),
                    ha="right", va="center", fontsize=6.5,
                    color=col if k >= 8 else "0.55")
        gp = gate[(gate.band == band) & np.isclose(gate.s, s0)]
        p = float(gp.gate_p.iloc[0]) if not gp.empty else np.nan
        glyph = "◉" if p < 0.05 else ("◐" if p < 0.15 else "○")
        gcol = col if p < 0.05 else ("#d99a1f" if p < 0.15 else "0.6")
        ax.set_title(rf"{C.BTeX[band]}  {glyph}  cohort $p={p:.3f}$",
                     color=gcol, fontsize=12, fontweight="bold", loc="left")

    for ax in axes[:, 0]:
        ax.set_yticks(list(ymap.values()))
        ax.set_yticklabels([p.replace("Pat_", "") for p in C.COHORT], fontsize=8)
    for ax in axes[-1, :]:
        ax.set_xlabel(C.YLAB_RHO)
    for ax in axes.ravel():
        ax.tick_params(axis="y", length=0)
        ax.set_ylim(-0.6, len(C.COHORT) - 0.4)

    handles = [
        Line2D([0], [0], marker="|", ls="", color="0.5", mew=1.8, ms=9,
               label=r"own null $p_{95}$"),
        Line2D([0], [0], marker="o", ls="", mfc="0.3", mec="white", ms=6, label="clears own null"),
        Line2D([0], [0], marker="o", ls="", mfc="white", mec="0.3", mew=1.0, ms=6, label="within null"),
        Line2D([0], [0], marker="o", ls="", mfc=C_ANTI, mec="white", ms=6, label="anti"),
        Line2D([0], [0], ls="", marker="", label=r"right margin  $k/16$ = # scales cleared"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.03),
               ncol=5, frameon=False, fontsize=9, handletextpad=0.5, columnspacing=1.4)
    fig.tight_layout(rect=(0, 0.02, 1, 1))

    C.FIGDIR.mkdir(parents=True, exist_ok=True)
    out = C.FIGDIR / "fig_perpatient_bands_forest.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


if __name__ == "__main__":
    main()
