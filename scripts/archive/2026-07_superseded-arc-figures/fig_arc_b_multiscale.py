#!/usr/bin/env python3
r"""fig:arc_b — the inference-specific component is MULTISCALE (Results §2, R2.2).

CORE MESSAGE: the persistent beta inference-specific component is not tied to the
hierarchy's finest scale. Re-tested against the same strength-matched surrogate at the
coarsest STILL-RESOLVED diffusion time (tau ~ 2.6/lambda_max, the scale of multi-step
paths), it clears the null just as at the fine scale:
  tau = 1.00 : median +0.09, surrogate +0.01, p = 0.010 (6/10 patients)
  tau = 2.61 : median +0.13, surrogate +0.01, p = 0.005 (7/10 patients)
No other band clears at either scale (strongest control p ~ 0.25). The surrogate median
stays flat (+0.005 -> +0.009) across scales, so the coarse-scale signal is a real margin
above the null, not a generic coarse-graining rise.

The scale bound is set by the operator: with Z_eff the number of diffusion modes the
propagator retains (of N ~ 118 contacts), the tree is finely resolved at tau = 1
(Z_eff ~ 72) and still resolved at tau ~ 2.6 (Z_eff ~ 35), but by tau ~ 7 the propagator
has relaxed toward the uniform mode (Z_eff ~ 6) with no hierarchy left to read -- so that
coarsest point is in a COLLAPSE ZONE and is not part of the test.

Reads : data/audit/consolidation_arc_rhosym/arc_scale_cohort_verdict.csv
        data/audit/consolidation_arc_rhosym/arc_scale_null_R200.csv
Writes: data/reports/results_section2/fig_arc_b_multiscale.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

BASE = ROOT / "data/audit/consolidation_arc_rhosym"
OUT = ROOT / "data/reports/results_section2/fig_arc_b_multiscale.pdf"

FUNC = "T_infspec_pe"
C_BETA = "#d1491c"           # the inference-carrying band
C_SURR = "#6f6f6f"           # matched-strength surrogate median
C_CTRL = "#c2c2c2"           # control bands (null at every scale)
# diffusion modes retained (paragraph), by sorted-tau position: fine, mesoscale, collapse
Z_EFF_BY_IDX = [72, 35, 6]


def main():
    cv = pd.read_csv(BASE / "arc_scale_cohort_verdict.csv")
    cv = cv[cv.functional == FUNC].copy()
    pp = pd.read_csv(BASE / "arc_scale_null_R200.csv")

    taus = sorted(cv.tau_mult.unique())               # [1.0, 2.61, 6.81] (exact floats)
    Z_EFF = dict(zip(taus, Z_EFF_BY_IDX))
    tested = taus[:2]                                 # fine + coarsest still-resolved
    collapse = taus[2:]                               # near-uniform propagator (dropped)

    fig, ax = plt.subplots(figsize=(8.6, 6.3))

    # collapse zone: coarser than the resolved range -> propagator near uniform
    # (annotation anchored with get_xaxis_transform: x in data, y in axes fraction)
    if collapse:
        x0 = np.sqrt(taus[1] * collapse[0])           # geometric midpoint on log-x
        ax.axvspan(x0, taus[-1] * 1.35, color="0.87", alpha=0.55, zorder=0)
        ax.text(collapse[0], 0.985,
                "collapse zone\n(near-uniform propagator,\nno hierarchy to read)",
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=9.0, color="0.35", style="italic")

    # control bands (everything but beta): thin grey, null at every scale
    for b in [x for x in cv.band.unique() if x != "beta"]:
        d = cv[cv.band == b].sort_values("tau_mult")
        ax.plot(d.tau_mult, d.med_obs, color=C_CTRL, lw=1.3, marker="o", ms=4,
                zorder=2, alpha=0.9)

    # beta IQR ribbon (per-patient spread -- honest about the large IQR)
    bpp = pp[pp.band == "beta"]
    q25 = [bpp[np.isclose(bpp.tau_mult, t)].T_infspec_pe_obs.quantile(0.25) for t in taus]
    q75 = [bpp[np.isclose(bpp.tau_mult, t)].T_infspec_pe_obs.quantile(0.75) for t in taus]
    ax.fill_between(taus, q25, q75, color=C_BETA, alpha=0.14, zorder=1,
                    label="_nolegend_")

    # beta surrogate median (the flat null) and beta observed median (the signal)
    bcv = cv[cv.band == "beta"].sort_values("tau_mult")
    ax.plot(bcv.tau_mult, bcv.med_surr, color=C_SURR, lw=1.8, ls="--", marker="s",
            ms=5, zorder=3)
    ax.plot(bcv.tau_mult, bcv.med_obs, color=C_BETA, lw=3.0, zorder=4)

    # markers: filled star at tested+cleared scales, open ring in the collapse zone
    for _, r in bcv.iterrows():
        t, m, p = r.tau_mult, r.med_obs, r.wilcoxon_p
        if t in collapse:
            ax.scatter(t, m, s=150, marker="o", facecolors="white",
                       edgecolors=C_BETA, linewidths=2.0, zorder=6)
        else:
            filled = p < 0.05
            ax.scatter(t, m, s=230, marker="*",
                       facecolors=C_BETA if filled else "white",
                       edgecolors=C_BETA, linewidths=1.6, zorder=6)
            ax.annotate(f"$p={p:.3f}$", (t, m), textcoords="offset points",
                        xytext=(9, 11), fontsize=12, color=C_BETA, fontweight="bold")

    ax.axhline(0.0, color="0.55", lw=1.0, ls=":", zorder=1)

    # Z_eff annotation strip (x in data, y in axes fraction)
    ymax = max(q75) * 1.28
    for t in taus:
        ax.text(t, 0.86, rf"$Z_{{\mathrm{{eff}}}}\!\approx\!{Z_EFF[t]}$",
                transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=10.5, color="0.3")

    ax.set_xscale("log")
    ax.set_xticks(taus)
    ax.set_xticklabels([f"{t:.1f}".rstrip("0").rstrip(".") for t in taus], fontsize=13)
    ax.set_xlim(taus[0] * 0.82, taus[-1] * 1.4)
    ax.set_ylim(min(q25) - 0.03, ymax + 0.10)
    ax.set_xlabel(r"diffusion time $\tau$  (units of $1/\lambda_{\max}$)"
                  r"       coarser $\rightarrow$", fontsize=14)
    ax.set_ylabel(r"inference-specific persistence  $T_{\mathrm{inf}\cdot e}$",
                  fontsize=14)
    ax.spines[["top", "right"]].set_visible(False)

    handles = [
        Line2D([0], [0], color=C_BETA, lw=3.0, marker="*", ms=14, mfc=C_BETA,
               mec=C_BETA, label=rf"{BRAIN_BAND_TEX_DICT['beta']}  (inference-specific)"),
        Line2D([0], [0], color=C_SURR, lw=1.8, ls="--", marker="s", ms=6,
               label="matched-strength surrogate median"),
        Line2D([0], [0], color=C_CTRL, lw=1.3, marker="o", ms=5,
               label=r"control bands (null at every scale, min $p\approx0.25$)"),
        Line2D([0], [0], marker="*", ls="", mfc=C_BETA, mec=C_BETA, ms=15,
               label=r"clears matched-strength ($p<0.05$)"),
        Patch(facecolor=C_BETA, alpha=0.14, label=r"$\beta$ inter-quartile range (patients)"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.06),
               ncol=2, frameon=False, fontsize=10.5, handletextpad=0.6,
               columnspacing=1.6)

    fig.subplots_adjust(left=0.13, right=0.965, top=0.93, bottom=0.24)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("fig:arc_b — multiscale robustness of the inference-specific beta component\n")
    for _, r in bcv.iterrows():
        tag = "COLLAPSE" if r.tau_mult in collapse else "resolved"
        print(f"  tau={r.tau_mult:5.2f} ({tag}, Z_eff~{Z_EFF[r.tau_mult]:2d}) "
              f"med_obs={r.med_obs:+.3f} med_surr={r.med_surr:+.3f} "
              f"p={r.wilcoxon_p:.3f}")
    ctrl_min = cv[(cv.band != "beta") & (cv.tau_mult.isin(tested))]
    print(f"  strongest control (tested scales): p = {ctrl_min.wilcoxon_p.min():.3f}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
