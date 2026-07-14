#!/usr/bin/env python3
r"""talk_fig_epi_band_dissociation -- one operator, two readouts (slide 17, dark).

The SAME diffusion operator sorts the frequency bands into two jobs. Read every band twice:

  COGNITION  = the cophenetic trace (does the task reorganization persist into rest_post,
               representative under matched-strength + LOO + >=3/10?).  -> alpha, beta ONLY.
  EPILEPSY   = the relational seizure-zone marker (strength-residual seed-affinity AUC, how
               well the seizure contacts co-diffuse, above their own matched-strength null).
               -> strongest + most consistent in delta, low_gamma, beta.

The dissociation reads off the two columns:
  * cognition is SHARP -- only alpha + beta clear the trace gate.
  * epilepsy is broader but peaks at delta / low_gamma / beta; theta is the weakest in BOTH
    (the selectivity control: a global artefact would light every band).
  * beta is the one band lit in BOTH columns -> the BRIDGE.
  * delta is the cleanest face of the split: it marks the seizure zone (AUC 0.83) yet carries
    NO cognitive trace.

Honest note: high_gamma also beats its own null in 7/10 (AUC 0.75) -- shown at full value; it is
treated as patient-specific in the band taxonomy, not part of the clean cohort epilepsy family.
Nothing is hidden; the reader sees every band's real number.

Reads : data/sparsified_arc/epi_arc_mst020/per_cell.csv      (relational marker AUC + beats-null)
        data/sparsified_arc/controls_ladder_apples/per_cell.csv (cophenetic trace, rep gate)
Writes: data/outputs/figures/talk/fig_epi_band_dissociation.{png,pdf}  (transparent, dark)
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] /
                       "01_compute" / "figures_embedded" / "new_results_sec1"))
import _common as C  # noqa: E402

C.use_lrg_style()
INK = "white"
plt.rcParams.update({
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
    "axes.titlecolor": INK, "xtick.color": INK, "ytick.color": INK,
})
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BEAT_GATE = 7                                  # beats own null in >= 7/10 -> filled epilepsy dot
OUT = C.ROOT / "data" / "outputs" / "figures" / "talk" / "fig_epi_band_dissociation"


def _cog_rep():
    lad = pd.read_csv(C.LADDER / "per_cell.csv")
    return {b: C.representativeness_gate(lad, "coph_meso", b, "T_test")["level"] for b in BANDS}


def _epi_auc():
    pc = pd.read_csv(C.ROOT / "data" / "sparsified_arc" / "epi_arc_mst020" / "per_cell.csv")
    out = {}
    for b in BANDS:
        d = pc[pc.band == b]
        out[b] = (float(d.auc_multi.median()), int((d.p_null_multi < 0.05).sum()))
    return out


def main():
    cog, epi = _cog_rep(), _epi_auc()
    ys = {b: len(BANDS) - 1 - i for i, b in enumerate(BANDS)}   # delta on top

    fig = plt.figure(figsize=(11.4, 6.4))
    axc = fig.add_axes([0.235, 0.14, 0.165, 0.72])             # cognition column (dots)
    axe = fig.add_axes([0.470, 0.14, 0.470, 0.72])             # epilepsy column (lollipops)

    # beta = the bridge: a soft full-width highlight behind its row in both panels
    yb = ys["beta"]
    for ax, x0, x1 in ((axc, -0.6, 0.6), (axe, 0.47, 1.02)):
        ax.axhspan(yb - 0.45, yb + 0.45, xmin=0, xmax=1, color="0.82", alpha=0.10, zorder=0)

    # ---- COGNITION column: filled band dot if representative (alpha, beta), else white ring
    for b in BANDS:
        y, col = ys[b], C.band_color(b)
        rep = cog[b] == "rep"
        if rep:
            axc.scatter(0, y, s=560, facecolor=col, edgecolor="white", lw=1.2, zorder=4)
        else:
            axc.scatter(0, y, s=470, facecolor="none", edgecolor="white", lw=1.6, zorder=3)
    axc.set_xlim(-0.6, 0.6)
    axc.set_ylim(-0.7, len(BANDS) - 0.3)
    axc.set_xticks([])
    axc.set_yticks(list(ys.values()))
    axc.set_yticklabels([C.BTeX[b] for b in BANDS], fontsize=26)
    for t, b in zip(axc.get_yticklabels(), BANDS):
        t.set_color(C.band_color(b))
    axc.set_title("cognition\n(persisting trace)", fontsize=15, pad=12)
    for s in ("top", "right", "left", "bottom"):
        axc.spines[s].set_visible(False)
    axc.tick_params(length=0)

    # ---- EPILEPSY column: lollipop 0.5 -> AUC, filled if beats own null in >= BEAT_GATE/10
    axe.axvline(0.5, color="0.55", lw=1.1, ls="--", zorder=1)       # strength alone / chance
    for b in BANDS:
        y, col = ys[b], C.band_color(b)
        auc, nb = epi[b]
        filled = nb >= BEAT_GATE
        axe.plot([0.5, auc], [y, y], color=col, lw=4.0, zorder=3, solid_capstyle="round",
                 alpha=0.95 if filled else 0.5)
        axe.scatter(auc, y, s=300, facecolor=col if filled else "none",
                    edgecolor="white" if filled else col, lw=1.6 if filled else 2.4, zorder=5)
        axe.text(auc + 0.02, y, f"{auc:.2f}", va="center", ha="left", fontsize=15,
                 color=col, fontweight="bold" if filled else "normal")
        axe.text(1.04, y, f"{nb}/10", va="center", ha="right", fontsize=12.5, color="0.8")
    axe.text(0.5, len(BANDS) - 0.35, "chance", fontsize=11, color="0.6", ha="center", va="bottom")
    axe.text(1.04, len(BANDS) - 0.35, "beats\nown null", fontsize=10.5, color="0.8",
             ha="right", va="bottom")
    axe.set_xlim(0.47, 1.08)
    axe.set_ylim(-0.7, len(BANDS) - 0.3)
    axe.set_yticks([])
    axe.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9])
    axe.set_xticklabels(["0.5", "0.6", "0.7", "0.8", "0.9"], fontsize=13)
    axe.set_xlabel("seizure-zone marker  —  strength-residual seed-affinity AUC", fontsize=13)
    axe.set_title("epilepsy\n(seizure-zone marker)", fontsize=15, pad=12)
    for s in ("top", "right", "left"):
        axe.spines[s].set_visible(False)
    axe.tick_params(axis="y", length=0)

    # bridge / control callouts
    axe.annotate(r"$\beta$ = the BRIDGE", xy=(epi["beta"][0], yb), xytext=(0.86, yb + 0.42),
                 fontsize=14, color=C.band_color("beta"), fontweight="bold",
                 ha="center", va="bottom",
                 arrowprops=dict(arrowstyle="-", color=C.band_color("beta"), lw=1.4))
    yt = ys["theta"]
    axe.text(0.52, yt + 0.34, "clears neither gate — the control", fontsize=11, color="0.72",
             va="bottom", ha="left", style="italic")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}.{ext}", transparent=True, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"wrote {OUT}.png + .pdf")
    print("  band        cognition      epilepsy AUC (beats/10)")
    for b in BANDS:
        print(f"  {b:11s} {cog[b]:6s}        {epi[b][0]:.3f} ({epi[b][1]}/10)")


if __name__ == "__main__":
    main()
