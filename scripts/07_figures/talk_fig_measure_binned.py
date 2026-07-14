#!/usr/bin/env python3
r"""talk_fig_measure_binned -- rho^coph binned + honestly calibrated (S2 companion).

TOP: the correlation, binned (not a point cloud). Pairs are grouped into quantile bins by
Delta^task; each bin shows mean rank of Delta^rest +/- 95% CI. A rising line = the task's
reshaping predicts rest's. References: flat matched-strength null, perfect rho=1 diagonal.

BOTTOM: how to read rho=0.48 (a number line, calibrated per the adversarial-verification verdict
2026-07-13). The LOAD-BEARING reference is the NULL: 0.48 is ~8x the matched-strength surrogate
median and outside its 95th percentile. The honest reason 1.0 is unreachable is the model-free
split-half reliability 0.53 (two halves of the SAME rest recording only reproduce the hierarchy
at 0.53). The perfect-persistence "ceiling" 0.62 is shown ONLY as a lighter generous upper bound
-- it shares the D_task term (task-noise treated as signal), which is why it sits ABOVE 0.53; the
true trace ceiling is lower. NO "% of ceiling" progress bar (that ratio is spin). Cohort strip
shows Pat_08 (0.48) is a top-of-range exemplar (cohort median 0.20); the cohort trace claim rests
on the null-based 10/10 panel (Wilcoxon p=0.001, 16/16 scales), not this number.

Exemplar Pat_08 beta, s=5.6.  Writes: data/outputs/figures/talk/fig_measure_binned.png
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from scipy.stats import spearmanr, rankdata

sys.path.insert(0, str(Path(__file__).resolve().parents[1] /
                       "01_compute" / "figures_embedded" / "new_results_sec1"))
import _common as C  # noqa: E402

C.use_lrg_style()

PAT, BAND = "Pat_08", "beta"
NBIN = 12
C_OBS, C_REL, C_CEIL, C_NULL = "#136b32", "#444444", "#c2c2c2", "#9a9a9a"
CAP = "0.45"
OUT = C.ROOT / "data" / "outputs" / "figures" / "talk" / "fig_measure_binned.png"


def _triu(M):
    return M[np.triu_indices(M.shape[0], 1)]


def main():
    D = {ph: C.coph_square_at_scale(C.load_phase(PAT, ph, BAND))
         for ph in ("A", "B", "task_test", "rest_post")}
    DA, DB, DT, DP = (_triu(D[k]) for k in ("A", "B", "task_test", "rest_post"))
    dtask = np.concatenate([DT - DA, DT - DB])
    drest = np.concatenate([DP - DB, DP - DA])
    rho = float(spearmanr(dtask, drest).correlation)
    rel = float(spearmanr(DA, DB).correlation)               # honest model-free split-half reliability
    ceil = float(spearmanr(DT - DA, DT - DB).correlation)    # generous upper bound (shares D_task)

    # matched-strength null + cohort context from the canonical per-patient sweep
    pp = pd.read_csv(C.MS / "per_patient_scale.csv"); s0 = C.SGRID[C.I_REPORT]
    r8 = pp[(pp.patient == PAT) & (pp.band == BAND) & np.isclose(pp.s, s0)].iloc[0]
    null_med, null_p95, obs = float(r8.surr_p50), float(r8.surr_p95), float(r8.obs_rho)
    coh = pp[(pp.band == BAND) & np.isclose(pp.s, s0)].obs_rho
    coh_med, coh_lo, coh_hi = float(coh.median()), float(coh.min()), float(coh.max())
    x_null = obs / null_med

    # jittered ranks (break tied-merge-height stripes so the density reads as a diagonal ridge)
    jr = np.random.default_rng(7)
    rx = np.clip(rankdata(dtask) / dtask.size + jr.uniform(-0.012, 0.012, dtask.size), 0, 1)
    ry = np.clip(rankdata(drest) / drest.size + jr.uniform(-0.012, 0.012, drest.size), 0, 1)

    from matplotlib.colors import LinearSegmentedColormap
    GREENS = LinearSegmentedColormap.from_list("g", ["#e9f4ea", "#8fca90", "#136b32"])

    from scipy.ndimage import gaussian_filter
    from matplotlib.colors import LogNorm
    fig = plt.figure(figsize=(8.1, 6.7))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.4, 3.3], wspace=0.30,
                          left=0.085, right=0.92, top=0.96, bottom=0.11)

    # ---- LEFT: vertical filled bar, cut to the reachable max ----
    g = fig.add_subplot(gs[0])
    g.set_xlim(0, 1); g.set_ylim(-0.008, rel * 1.015)
    W = 0.62; x0 = (1 - W) / 2
    g.add_patch(Rectangle((x0, 0), W, rel, facecolor="0.9", edgecolor="0.5", lw=1.0, zorder=1))
    g.add_patch(Rectangle((x0, 0), W, obs, facecolor=C_OBS, edgecolor="none", zorder=2))
    g.text(0.5, obs * 0.5, rf"$\mathbf{{{obs:.2f}}}$", rotation=90, ha="center", va="center",
           color="white", fontweight="bold", fontsize=13, zorder=3)
    # dashed white null mark: matched-strength surrogate 95th pct = the level obs must clear
    g.plot([x0, x0 + W], [null_p95, null_p95], color="white", lw=1.4,
           ls=(0, (4, 3)), zorder=4, solid_capstyle="butt")
    g.text(x0 + W + 0.05, null_p95, "null", ha="left", va="center", fontsize=8, color="0.5")
    g.text(x0 - 0.14, 0, "0", ha="right", va="center", fontsize=9, color="0.45")
    g.text(x0 - 0.14, rel, f"{rel:.2f}", ha="right", va="center", fontsize=9, color="0.45")
    for sp in g.spines.values():
        sp.set_visible(False)
    g.set_xticks([]); g.set_yticks([])

    # ---- RIGHT: normalized, log-scaled density map of the correlation ----
    ax = fig.add_subplot(gs[1])
    H2, _, _ = np.histogram2d(rx, ry, bins=46, range=[[0, 1], [0, 1]])
    Hs = gaussian_filter(H2.T, sigma=1.7)
    Hs = Hs / Hs.max()                            # relative density, peak = 1
    vmin = 0.1                                     # one clean decade [10^-1, 10^0]
    im = ax.imshow(np.clip(Hs, vmin, None), origin="lower", extent=[0, 1, 0, 1], cmap=GREENS,
                   aspect="equal", norm=LogNorm(vmin=vmin, vmax=1.0), zorder=2)
    ax.plot([0, 1], [0, 1], color="0.6", lw=1.0, ls="--", alpha=0.55, zorder=3)
    ax.text(0.035, 0.965, rf"$\rho^{{\mathrm{{coph}}}} = {rho:.2f}$", transform=ax.transAxes,
            ha="left", va="top", fontsize=18, color="#0f5427", fontweight="bold")
    ax.set_xlabel(r"pairs ranked by $\Delta^{\mathrm{task}}$", fontsize=11)
    ax.set_ylabel(r"pairs ranked by $\Delta^{\mathrm{rest}}$", fontsize=11)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    cax = ax.inset_axes([1.03, 0.0, 0.03, 1.0])
    cb = fig.colorbar(im, cax=cax); cb.set_label("pair density", fontsize=9)
    cb.ax.minorticks_off()
    cb.set_ticks([0.1, 1.0]); cb.set_ticklabels([r"$10^{-1}$", r"$10^{0}$"])
    cb.ax.tick_params(labelsize=8); cb.outline.set_visible(False)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"{PAT} {BAND}: rho={rho:.3f} obs={obs:.3f} null_med={null_med:.3f} p95={null_p95:.3f} "
          f"x_null={x_null:.1f} rel={rel:.3f} ceil={ceil:.3f} cohort_med={coh_med:.3f}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
