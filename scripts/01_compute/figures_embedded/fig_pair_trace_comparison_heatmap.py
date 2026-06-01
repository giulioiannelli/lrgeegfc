"""Primary verdict figure for audit_73 — the 4-measure × 3-primitive comparison.

A verdict heatmap: 12 rows (3 primitives × 4 measures, grouped) × 6 band
columns, cell shade = cohort-paired one-sided Wilcoxon ``−log10(p)`` for
``obs > surrogate-median`` under the matched-strength null. Because the test is
one-sided "greater", anti-trace bands (negative obs − surrogate) fall to high p
→ dark cells automatically, so brightness reads unambiguously as "trace-
positive AND matched-strength-significant". Cells with p < 0.05 carry a bold
outline; each cell is annotated with the p-value and per-patient ``n_above/10``.

Reading the three D_coph rows (bottom group) against each other is the whole
point: ``ρ_Spearman`` / ``ρ_Pearson`` light up at α, β (the current §5.3
signature); the asymmetric ``s_TR`` slope does NOT; ``R²`` lights up everywhere
(including anti-trace θ), so it cannot discriminate trace from drift.

No n/10 reference lines on a per-patient axis (the n_above text is a count, not
a gate). PDF only, vector. No suptitle.

Output: data/audit/pair_trace_measures_comparison/figures/comparison_heatmap.pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.notebook import move_to_rootf
from lrg_eegfc.visuals.styles import use_lrg_style

move_to_rootf(pathname="lrgeegfc")
use_lrg_style()

AUDIT_DIR = Path("data/audit/pair_trace_measures_comparison")
FIG_DIR = AUDIT_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

PRIMITIVES = ["rhohat", "D", "Dcoph"]
PRIMITIVE_TEX = {
    "rhohat": r"$\hat{\rho}(\tau_{\max})$",
    "D": r"$D(\tau_{\max})$",
    "Dcoph": r"$D_{\mathrm{coph}}$",
}
MEASURES = ["spearman", "pearson", "slope", "rsq"]
MEASURE_TEX = {
    "spearman": r"$\rho_{\mathrm{Spear}}$",
    "pearson": r"$\rho_{\mathrm{Pear}}$",
    "slope": r"$s_{\mathrm{TR}}$",
    "rsq": r"$R^2$",
}
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]

P_FLOOR = 1e-3  # -log10 clip so the colour scale is not blown out by ~0 p


def main() -> None:
    co = pd.read_csv(AUDIT_DIR / "cohort_summary.csv")
    key = co.set_index(["primitive", "measure", "band"])

    # row order: primitive-major, measure-minor
    row_specs = [(p, m) for p in PRIMITIVES for m in MEASURES]
    nrows, ncols = len(row_specs), len(BANDS)

    neglogp = np.full((nrows, ncols), np.nan)
    anti = np.zeros((nrows, ncols), dtype=bool)
    pvals = np.full((nrows, ncols), np.nan)
    nabove = np.full((nrows, ncols), np.nan)
    for i, (prim, meas) in enumerate(row_specs):
        for j, band in enumerate(BANDS):
            try:
                r = key.loc[(prim, meas, band)]
            except KeyError:
                continue
            p = float(r["wilcoxon_p_greater"])
            pvals[i, j] = p
            nabove[i, j] = float(r["n_above_p05"])
            anti[i, j] = float(r["obs_median"]) < float(r["surr_median"])
            if np.isfinite(p):
                neglogp[i, j] = -np.log10(max(p, P_FLOOR))

    fig, ax = plt.subplots(figsize=(8.6, 9.4))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("0.92")
    im = ax.imshow(neglogp, cmap=cmap, aspect="auto",
                   vmin=0.0, vmax=-np.log10(P_FLOOR))

    # cell annotations
    for i in range(nrows):
        for j in range(ncols):
            p = pvals[i, j]
            if not np.isfinite(p):
                continue
            sig = p < 0.05
            shade = neglogp[i, j] / (-np.log10(P_FLOOR))
            tcol = "white" if shade > 0.55 else "black"
            ptxt = f"{p:.3f}" if p >= P_FLOOR else f"<{P_FLOOR:g}"
            mark = "▲" if not anti[i, j] else "▽"
            ax.text(j, i - 0.16, f"{mark} {ptxt}", ha="center", va="center",
                    fontsize=6.6, color=tcol,
                    fontweight=("bold" if sig else "normal"))
            ax.text(j, i + 0.22, f"{int(nabove[i, j])}/10", ha="center",
                    va="center", fontsize=6.0, color=tcol, alpha=0.85)
            if sig:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       edgecolor="#d62728", lw=2.0, zorder=5))

    # group separators between primitives
    for g in (len(MEASURES), 2 * len(MEASURES)):
        ax.axhline(g - 0.5, color="white", lw=2.5)

    ax.set_xticks(range(ncols))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=11)
    ax.set_yticks(range(nrows))
    ax.set_yticklabels([MEASURE_TEX[m] for _, m in row_specs], fontsize=10)
    ax.tick_params(length=0)

    # primitive group labels on the right
    for gi, prim in enumerate(PRIMITIVES):
        yc = gi * len(MEASURES) + (len(MEASURES) - 1) / 2
        ax.text(ncols - 0.35, yc, PRIMITIVE_TEX[prim], rotation=270,
                ha="left", va="center", fontsize=12,
                transform=ax.transData, clip_on=False)

    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.10)
    cb.set_label(r"$-\log_{10}\,p$ (matched-strength Wilcoxon, one-sided)",
                 fontsize=9)
    # 0.05 reference on the colorbar
    cb.ax.axhline(-np.log10(0.05), color="#d62728", lw=1.5)
    cb.ax.text(1.6, -np.log10(0.05), "  p=.05", color="#d62728",
               fontsize=7, va="center", transform=cb.ax.get_yaxis_transform())

    ax.set_title("▲ trace-positive   ▽ anti   red box: p < 0.05",
                 fontsize=9, pad=6)
    fig.tight_layout()
    out = FIG_DIR / "comparison_heatmap.pdf"
    fig.savefig(out, transparent=True)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
