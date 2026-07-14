#!/usr/bin/env python3
r"""talk_fig_measure_scatter -- rho^coph is just the correlation in this cloud (S2 companion).

The colored dendrograms show WHERE the hierarchy is preserved (per-contact); this shows the
measure ITSELF. Each dot is one contact-pair: x = how the TASK moved its cophenetic distance
(Delta^task = D_task - D_rsPreA), y = how REST-post moved it (Delta^rest = D_post - D_rsPreB).
The rank-trend of the cloud IS rho^coph -- a positive tilt = the task's reshaping persists.
One exemplar (Pat_08 beta, s=5.6), disclosed. Purely conceptual: the literal picture of
    rho^coph(tau) = 1/2 [ rho_S(Delta^task_A, Delta^rest_B) + (A<->B) ].

Writes: data/outputs/figures/talk/fig_measure_scatter.png  (transparent, Canva-ready)
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, rankdata

sys.path.insert(0, str(Path(__file__).resolve().parents[1] /
                       "01_compute" / "figures_embedded" / "new_results_sec1"))
import _common as C  # noqa: E402

C.use_lrg_style()

PAT, BAND = "Pat_08", "beta"
C_PT, C_TREND = "#2f8f4f", "#136b32"
CAP = "0.5"
OUT = C.ROOT / "data" / "outputs" / "figures" / "talk" / "fig_measure_scatter.png"


def _triu(D):
    return D[np.triu_indices(D.shape[0], 1)]


def main():
    D = {ph: C.coph_square_at_scale(C.load_phase(PAT, ph, BAND))
         for ph in ("A", "B", "task_test", "rest_post")}
    # pool the two symmetric baseline-half arms so the shown cloud IS rho^coph (not one arm)
    dtask = np.concatenate([_triu(D["task_test"] - D["A"]), _triu(D["task_test"] - D["B"])])
    drest = np.concatenate([_triu(D["rest_post"] - D["B"]), _triu(D["rest_post"] - D["A"])])
    rho = float(spearmanr(dtask, drest).correlation)   # ~ the symmetric rho^coph (0.48)

    # Only pairs the TASK actually moved (upper half by |Delta^task|). Unbiased -- we condition
    # on the task, NOT on rest, so reverters ("task moved it, rest let it go") are KEPT, they
    # just land near the rest-midline. Removing the immobile pairs de-clutters without inflating.
    thr = np.percentile(np.abs(dtask), 50.0)
    m = np.abs(dtask) >= thr
    dt, dr = dtask[m], drest[m]
    n = dt.size
    rho_sub = float(spearmanr(dt, dr).correlation)
    agree = np.sign(dt) == np.sign(dr)
    frac_agree = float(np.mean(agree))

    # signed-rank cloud: rank of signed Delta in [0,1] (0.5 ~ no move). Agree = same side of 0.5.
    # Cophenetic merge-heights are heavily tied (pairs under the same top merge share a height),
    # so raw ranks stack into stripes; a small jitter breaks the ties into a readable cloud.
    jrng = np.random.default_rng(7)
    rx = np.clip(rankdata(dt) / n + jrng.uniform(-0.014, 0.014, n), 0, 1)
    ry = np.clip(rankdata(dr) / n + jrng.uniform(-0.014, 0.014, n), 0, 1)

    fig, ax = plt.subplots(figsize=(6.6, 6.2))
    ax.plot([0, 1], [0, 1], color="0.72", lw=1.1, ls="--", zorder=1)
    ax.text(0.985, 0.985, r"$\rho=1$", color="0.6", fontsize=10, ha="right", va="top",
            rotation=45, rotation_mode="anchor")
    ax.scatter(rx[~agree], ry[~agree], s=5, color="0.72", alpha=0.07, lw=0, zorder=2)
    ax.scatter(rx[agree], ry[agree], s=5, color=C_PT, alpha=0.13, lw=0, zorder=3)
    xx = np.array([0.0, 1.0])
    ax.plot(xx, 0.5 + rho_sub * (xx - 0.5), color=C_TREND, lw=3.0, zorder=4)

    ax.set_xlabel(r"pairs ranked by $\Delta^{\mathrm{task}}$" "\n(how the task moved them)",
                  fontsize=12)
    ax.set_ylabel(r"pairs ranked by $\Delta^{\mathrm{rest}}$" "\n(how rest$_{\mathrm{post}}$ moved them)",
                  fontsize=12)
    ax.text(0.04, 0.965, rf"$\rho^{{\mathrm{{coph}}}} = {rho:.2f}$",
            transform=ax.transAxes, ha="left", va="top", fontsize=18, color=C_TREND,
            fontweight="bold")
    ax.text(0.04, 0.895, "the pairs the task moved:\n"
            rf"${frac_agree*100:.0f}\%$ keep their direction at rest",
            transform=ax.transAxes, ha="left", va="top", fontsize=10.5, color=CAP)
    ax.text(0.97, 0.05, "tilt toward the diagonal\n= task's reshaping persists",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=10.5, color=CAP)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
    ax.set_aspect("equal")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"{PAT} {BAND} s={C.S_REPORT:.2f}  rho^coph(pooled)={rho:+.2f}  "
          f"task-movers(50%)={n}  subset rho={rho_sub:+.2f}  sign-agree={frac_agree:.2f}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
