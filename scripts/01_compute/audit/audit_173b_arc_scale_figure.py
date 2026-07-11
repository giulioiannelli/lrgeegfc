#!/usr/bin/env python3
r"""audit_173b — Part B verdict figure: encoding vs inference over the scale sweep.

Two panels, -log10 cohort matched-strength gate p vs s = tau*lam_max (log x):
  (a) T_infspec_pe (inference-specific | encoding): beta clears at EVERY scale
      (beta solid; all controls dashed, null) -> inference is beta-only and
      scale-broad, mildly mesoscale-favouring.
  (b) T_learn (encoding echo): alpha/beta solid (candidates), rest dashed. alpha
      is FINE-scale-favouring (weakens monotonically toward coarse); delta
      STRENGTHENS toward collapse (anatomy channel, rejected).

Together: encoding lives fine and dissolves at coarse; inference persists broadly
-> a component-level echo of the Part A alpha-mesoscale / beta-broad trace split.

Reads : data/audit/arc_scale_sweep_rhosym/arc_scale_cohort_verdict_s.csv
        data/audit/tau_sweep_trace/cohort.csv (median Fiedler s_F)
Writes: data/audit/arc_scale_sweep_rhosym/audit_173_arc_scale_verdict.pdf (+ .md)
"""
from __future__ import annotations
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

ROOT = setup_script_env()
use_lrg_style()

V = ROOT / "data" / "audit" / "arc_scale_sweep_rhosym" / "arc_scale_cohort_verdict_s.csv"
FIED = ROOT / "data" / "audit" / "tau_sweep_trace" / "cohort.csv"
OUT = ROOT / "data" / "audit" / "arc_scale_sweep_rhosym" / "audit_173_arc_scale_verdict.pdf"

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
# candidate band(s) per functional (solid); others dashed
CAND = {"T_infspec_pe": {"beta"}, "T_learn": {"alpha", "beta"}}
P_SIG = 0.05


def _panel(ax, v, k, sF_med):
    cand = CAND[k]
    for band in BANDS:
        g = v[(v.functional == k) & (v.band == band)].sort_values("s")
        if g.empty:
            continue
        s = g.s.values
        c = band in cand
        st = dict(color=band_color(band), lw=2.6 if c else 1.5,
                  ls="-" if c else (0, (5, 2)), alpha=1.0 if c else 0.6,
                  zorder=5 if c else 3)
        y = -np.log10(np.clip(g.wilcoxon_p.values, 1e-6, 1.0))
        ax.plot(s, y, **st)
        sig = g.wilcoxon_p.values < P_SIG
        ax.scatter(s[sig], y[sig], s=32, facecolor=st["color"], edgecolor="white",
                   linewidth=0.5, zorder=st["zorder"] + 1)
        ax.scatter(s[~sig], y[~sig], s=24, facecolor="white", edgecolor=st["color"],
                   linewidth=1.0, zorder=st["zorder"] + 1)
    ax.set_xscale("log")
    ax.set_xlabel(r"diffusion scale $s \equiv \tau\,\lambda_{\max}$")
    ax.axvline(1.0, color="0.35", lw=1.0, ls=":", zorder=1)
    ax.axvline(sF_med, color="0.6", lw=1.0, ls=(0, (1, 2)), zorder=1)
    ax.axhline(-np.log10(P_SIG), color="0.7", lw=0.8, zorder=0)
    ax.set_xticks([0.5, 1, 2, 3, 5, 8]); ax.set_xticklabels(["0.5", "1", "2", "3", "5", "8"])
    ax.set_ylabel(r"$-\log_{10}$ cohort gate $p(s)$")


def main():
    v = pd.read_csv(V)
    sF_med = float(np.median(pd.read_csv(FIED).groupby("band")["median_alpha_fiedler"].first()))
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.4, 4.5), sharey=True)
    _panel(axA, v, "T_infspec_pe", sF_med)
    _panel(axB, v, "T_learn", sF_med)
    axA.set_title(r"inference-specific $\,\mathrm{(T_{inf\,|\,enc})}$", fontsize=13)
    axB.set_title(r"encoding echo $\,\mathrm{(T_{learn})}$", fontsize=13)
    axA.text(0.02, 1.06, r"$\mathbf{a}$", transform=axA.transAxes, fontsize=15,
             fontweight="bold", va="bottom")
    axB.text(0.02, 1.06, r"$\mathbf{b}$", transform=axB.transAxes, fontsize=15,
             fontweight="bold", va="bottom")
    handles = [Line2D([0], [0], color=band_color(b), lw=2.0, label=BRAIN_BAND_TEX_DICT[b])
               for b in BANDS]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.03),
               ncol=6, frameon=False, handlelength=2.4)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True)
    plt.close(fig)
    print(f"wrote {OUT}")
    OUT.with_suffix(".md").write_text(
        "# audit_173 — encoding vs inference over the scale sweep\n\n"
        "`-log10` cohort matched-strength gate p vs `s = tau*lam_max` (log x). Dotted "
        "`s=1` = tau_min; fine dotted = median Fiedler. Filled = clears (p<0.05).\n\n"
        "**(a) inference-specific (T_inf|enc):** beta (solid) clears at EVERY scale "
        "(best 0.003 @ s~2.1, also 0.010 at tau_min); every other band stays null at "
        "every s -> inference is beta-only and scale-broad, mildly mesoscale-favouring "
        "(confirms audit_157).\n\n"
        "**(b) encoding echo (T_learn):** alpha/beta (solid) are FINE-scale-favouring — "
        "alpha strongest at s~0.85 (0.010), weakening to null past Fiedler. delta "
        "(dashed) STRENGTHENS toward collapse (0.014->0.003): the anatomy channel, "
        "rejected.\n\nTogether: encoding lives fine and dissolves at coarse; inference "
        "persists broadly — a component echo of the Part A alpha-meso / beta-broad split. "
        "Anchors: arm1 == audit_103d, sym(s=1) == audit_152 (both 1e-16). No paper "
        "change: all components clear at tau_min. Build: audit_173b_arc_scale_figure.py.\n")
    print(f"wrote {OUT.with_suffix('.md')}")


if __name__ == "__main__":
    main()
