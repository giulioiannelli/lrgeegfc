#!/usr/bin/env python3
r"""fig:trace_c (companion) — the task's cophenetic fingerprint, and its echo (Results §1, para c).

IDEA (no correlation on screen — the eye does the correlating): for one patient, lay the
change the task wrote into the LRG cophenetic hierarchy next to the change still present at
post-task rest. Each map is the per-pair change in cophenetic distance, contacts ordered by
the rest_pre dendrogram so the stable backbone is the pale bulk and the movers are the
saturated blocks. When the task's fingerprint (top) reappears in the held map (bottom) —
same red/blue blocks — that IS the trace, shown as a pattern match instead of a small ρ.
It reappears for β and α (green labels, clear the cohort gate); it does not for θ and δ
(grey labels). The pooled Spearman is only ~0.2 because the backbone is unchanged; the point
of this figure is that among the pairs that move, the pattern genuinely persists.

HONESTY (why split-half baselines): the naive maps D(task)−D(pre) and D(post)−D(pre) share
the same D(pre), which manufactures a matching pattern from the backbone alone — full-phase
δ/θ then look as "matched" as β (probe: δ +0.51, θ +0.44 vs β +0.72). So each map is built
against an INDEPENDENT half of rest_pre (task vs pre_A, post vs pre_B) — the ρ_sym
construction — which removes the shared-baseline artifact and leaves only the real, persistent
reorganization (β +0.50, α +0.34; θ −0.05, δ −0.11 on this patient). Exemplar Pat_08 (the
cleanest tracer); the cohort verdict is the gate/forest in (a).

Reads : LRG cophenetic caches (full task/post/pre + rest_pre halves), rho_sym gate summary.
Writes: data/reports/results_section1/fig_trace_c_diffmaps.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list
from scipy.spatial.distance import squareform
from scipy.stats import rankdata, spearmanr

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result

ROOT = setup_script_env()
use_lrg_style()

LRG_HALVES = CACHE_ROOT / "imcoh_lrg_halves"
GATE_SUMMARY = ROOT / "data" / "audit" / "rho_sym_gate" / "cohort_summary.csv"
OUT = ROOT / "data" / "reports" / "results_section1" / "fig_trace_c_diffmaps.pdf"

PATIENT = "Pat_08"
FC = "imcoh_abs"
BANDS = ["beta", "alpha", "theta", "delta"]      # 2 traces (green) then 2 nulls (grey)
GATE_PASS, GATE_FAIL = "#1f7a34", "#565656"
# diverging: pair pulled together (blue) ↔ pulled apart (red); pale centre = backbone (no change)
DIVCMAP = "coolwarm"


def _Dsq(phase, band, cache_root=None):
    r = (load_lrg_result(PATIENT, phase, band, fc_method=FC, cache_root=cache_root)
         if cache_root else load_lrg_result(PATIENT, phase, band, fc_method=FC))
    um = np.asarray(r.ultrametric_matrix)
    return squareform(um) if um.ndim == 1 else um


def _rankmap(delta_full):
    """Signed-rank image of the off-diagonal Δ, in [-1, 1]. Displaying RANKS makes the
    picture identical to what Spearman ρ measures: a global offset (θ: everything pulled
    together at task, apart at post) becomes a uniform rank field and disappears, leaving
    only the co-varying block structure — the actual trace."""
    n = delta_full.shape[0]
    iu = np.triu_indices(n, k=1)
    scaled = 2.0 * (rankdata(delta_full[iu]) - 1.0) / (iu[0].size - 1) - 1.0
    R = np.zeros((n, n))
    R[iu] = scaled
    return R + R.T


def band_maps(band):
    """(Δtask, Δrest) reordered by the rest_pre dendrogram + the arm-1 ρ_split (diagnostic)."""
    Z_pre = load_lrg_result(PATIENT, "rest_pre", band, fc_method=FC).linkage_matrix
    order = leaves_list(Z_pre)
    Dtask = _Dsq("task_test", band)
    Dpost = _Dsq("rest_post", band)
    DpreA = _Dsq("rest_pre_A", band, LRG_HALVES)
    DpreB = _Dsq("rest_pre_B", band, LRG_HALVES)
    dtask = Dtask - DpreA                              # task vs an independent half of rest_pre
    drest = Dpost - DpreB                              # held change vs the OTHER half
    iu = np.triu_indices(dtask.shape[0], k=1)
    rho = float(spearmanr(dtask[iu], drest[iu]).correlation)
    ix = np.ix_(order, order)
    return _rankmap(dtask)[ix], _rankmap(drest)[ix], rho


def main():
    gate = pd.read_csv(GATE_SUMMARY).set_index("band")["gate_p_sym"].to_dict()
    maps = {b: band_maps(b) for b in BANDS}

    print(f"fig:trace_c diffmaps — {PATIENT}, split-half baselines (ρ_sym arm-1 diagnostic):")
    for b in BANDS:
        print(f"  {b:11s} pattern-match ρ = {maps[b][2]:+.3f}   "
              f"gate_p_sym={gate.get(b, float('nan')):.3f}")

    ncol = len(BANDS)
    fig = plt.figure(figsize=(3.15 * ncol, 6.9))
    gs = fig.add_gridspec(2, ncol, left=0.085, right=0.90, top=0.86, bottom=0.06,
                          wspace=0.10, hspace=0.10)

    row_tex = [r"$\mathrm{task_{test}}$ vs $\mathrm{rest_{pre}}$",
               r"$\mathrm{rest_{post}}$ vs $\mathrm{rest_{pre}}$"]
    last = None
    for j, b in enumerate(BANDS):
        dtask, drest, _ = maps[b]
        lab_col = GATE_PASS if gate.get(b, 1.0) < 0.05 else GATE_FAIL
        for i, M in enumerate((dtask, drest)):
            ax = fig.add_subplot(gs[i, j])
            last = ax.pcolormesh(M, cmap=DIVCMAP, vmin=-1.0, vmax=1.0,
                                 rasterized=False)
            ax.set_box_aspect(1.0)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.invert_yaxis()
            for s in ax.spines.values():
                s.set_edgecolor(lab_col if i == 0 else "0.6")
                s.set_linewidth(1.6 if i == 0 else 0.8)
            if i == 0:
                ax.set_title(BRAIN_BAND_TEX_DICT.get(b, b), fontsize=22,
                             fontweight="bold", color=lab_col, pad=8)
            if j == 0:
                ax.set_ylabel(row_tex[i], fontsize=12.5, color="0.20")

    fig.text(0.02, 0.95, r"$\mathbf{c}$", fontsize=17, va="bottom", ha="left",
             fontweight="bold")
    fig.text(0.492, 0.965,
             r"the task's cophenetic fingerprint $\;\longrightarrow\;$ and its echo at rest",
             fontsize=13.5, ha="center", va="center", color="0.25")

    cax = fig.add_axes([0.915, 0.14, 0.014, 0.60])
    cb = fig.colorbar(last, cax=cax, extend="both")
    cb.set_label("per-pair cophenetic change, ranked within each map\n"
                 "blue = most pulled together   ·   red = most pulled apart",
                 rotation=270, labelpad=30, fontsize=10.5)
    cb.set_ticks([-1, 0, 1])
    cb.set_ticklabels(["−", "0", "+"])

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
