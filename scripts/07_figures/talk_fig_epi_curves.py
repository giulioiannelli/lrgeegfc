#!/usr/bin/env python3
r"""talk_fig_epi_curves -- the epilepsy marker as CLEAN curves (slide 16, dark).

THREE text-free data panels (all explanatory words live on the slide, not in the figure), all on
the mst@0.20 backbone -- the SAME single backbone that carries the cognitive trace (TMFG dropped
2026-07-16, one backbone only) -- against the same matched-strength null. Each shows per-patient
fluctuations (thin) under a bold cohort average:

  A  auc_by_band      -- strength-residual seed-affinity AUC across the six bands. One faint
                         polyline per patient, a bold cohort-median profile, IQR shaded, chance
                         at 0.5. Reads: the marker peaks at delta / low-gamma; alpha/theta weak.
  B  detector_roc     -- leave-one-patient-out node ROC. One thin ROC per patient, bold mean ROC
                         (AUC ~0.85); one right-hemi hub (Pat_15) hugs the diagonal, shown not
                         stated.
  C  precision_at_k   -- precision of the top-k shortlist. Thin per-patient, bold cohort mean,
                         dashed base rate. Top of the list ~7x enriched (prec@5 = 0.60).

Reads : data/sparsified_arc/epi_arc_mst020/per_cell.csv                              (per-band AUC)
        data/sparsified_arc/marker_detector_mst020/detector_node_predictions.csv     (p_soz / node)
Writes: data/outputs/figures/talk/fig_epi_auc_by_band.{png,pdf}
        data/outputs/figures/talk/fig_epi_detector_roc.{png,pdf}
        data/outputs/figures/talk/fig_epi_precision_at_k.{png,pdf}   (transparent, dark)
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
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
FIGDIR = C.ROOT / "data" / "outputs" / "figures" / "talk"
DET_DIR = C.ROOT / "data" / "sparsified_arc" / "marker_detector_mst020"     # LOPO detector (script 27)
EPI_SRC = C.ROOT / "data" / "sparsified_arc" / "epi_arc_mst020" / "per_cell.csv"

C_PAT = "#8a929e"        # per-patient fluctuation (neutral slate)
C_MEAN = "#f2f4f7"       # cohort average (bright)
C_REF = "#6a7079"        # chance / base-rate reference


def _finish(fig, ax, out):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", transparent=True, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"wrote {out}.png")


# ------------------------------------------------------------------ A: AUC by band
def auc_by_band():
    pc = pd.read_csv(EPI_SRC)
    piv = pc.pivot_table(index="patient", columns="band", values="auc_multi")[BANDS]
    x = np.arange(len(BANDS))

    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    ax.axhline(0.5, color=C_REF, lw=1.2, ls=(0, (5, 4)), zorder=1)          # chance
    for _, row in piv.iterrows():
        ax.plot(x, row.values, color=C_PAT, lw=1.1, alpha=0.38, zorder=2)   # fluctuation
    med = piv.median(0).values
    q1, q3 = piv.quantile(0.25).values, piv.quantile(0.75).values
    ax.fill_between(x, q1, q3, color=C_MEAN, alpha=0.12, zorder=3)          # IQR
    ax.plot(x, med, color=C_MEAN, lw=3.4, zorder=5, solid_capstyle="round")  # cohort median
    for j, b in enumerate(BANDS):                                           # band-coloured nodes
        ax.scatter(j, med[j], s=150, color=C.band_color(b), edgecolor="white",
                   lw=1.3, zorder=6)

    ax.set_xticks(x)
    ax.set_xticklabels([C.BTeX[b] for b in BANDS], fontsize=22)
    for t, b in zip(ax.get_xticklabels(), BANDS):
        t.set_color(C.band_color(b))
    ax.set_ylim(0.30, 1.0)
    ax.set_xlim(-0.4, len(BANDS) - 0.6)
    ax.set_ylabel("seizure-zone AUC", fontsize=15)
    ax.tick_params(axis="x", length=0)
    _finish(fig, ax, FIGDIR / "fig_epi_auc_by_band")
    print(f"  band AUC (mst@0.20, cohort-median): "
          + "  ".join(f"{b} {med[j]:.2f}" for j, b in enumerate(BANDS)))


# ------------------------------------------------------------------ B: detector ROC
def _roc(y, s):
    order = np.argsort(-s)
    y = y[order]
    P, N = y.sum(), (~y.astype(bool)).sum()
    tpr = np.concatenate([[0], np.cumsum(y) / max(P, 1)])
    fpr = np.concatenate([[0], np.cumsum(1 - y) / max(N, 1)])
    return fpr, tpr


def _conc(y, s):
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    return float((pos[:, None] > neg[None, :]).mean() + 0.5 * (pos[:, None] == neg[None, :]).mean())


def detector_roc():
    nd = pd.read_csv(DET_DIR / "detector_node_predictions.csv")
    grid = np.linspace(0, 1, 101)

    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    ax.plot([0, 1], [0, 1], color=C_REF, lw=1.2, ls=(0, (5, 4)), zorder=1)   # chance
    interp, aucs = [], []
    for pat, g in nd.groupby("patient"):
        yy, ss = g.is_soz.to_numpy(int), g.p_soz.to_numpy(float)
        fpr, tpr = _roc(yy, ss)
        ax.plot(fpr, tpr, color=C_PAT, lw=1.3, alpha=0.5, zorder=3)          # per-patient (one pop.)
        interp.append(np.interp(grid, fpr, tpr))
        aucs.append(_conc(yy, ss))
    ax.plot(grid, np.mean(interp, 0), color=C_MEAN, lw=3.6, zorder=5,
            solid_capstyle="round")                                          # mean ROC

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("false-positive rate", fontsize=14)
    ax.set_ylabel("true-positive rate", fontsize=14)
    ax.set_box_aspect(1.0)
    _finish(fig, ax, FIGDIR / "fig_epi_detector_roc")
    print(f"  mean per-patient AUC (node-level) = {np.nanmean(aucs):.2f}  "
          f"(worst {np.nanmin(aucs):.2f} = the right-hemi hub near chance)")


# ------------------------------------------------------------------ D: precision @ k
def precision_at_k(kmax=20, n_rand=2000):
    nd = pd.read_csv(DET_DIR / "detector_node_predictions.csv")
    ks = np.arange(1, kmax + 1)
    base, curves, ys = [], {}, {}
    for pat, g in nd.groupby("patient"):
        s, y = g.p_soz.to_numpy(float), g.is_soz.to_numpy(int)
        yy = y[np.argsort(-s)]
        curves[pat] = np.array([yy[:k].sum() / k for k in ks])
        base.append(y.mean()); ys[pat] = y

    # RANDOM-RANKING BAND: what cohort-mean precision@k looks like if the ranker were random.
    # Per trial, draw a random ranking per patient (label permutation), take precision@k, average
    # over patients -> one random cohort-mean curve; the 5-95% envelope over trials is the band.
    rng = np.random.default_rng(0)
    pats = list(ys)
    rand = np.empty((n_rand, kmax))
    for t in range(n_rand):
        per = np.empty((len(pats), kmax))
        for i, p in enumerate(pats):
            yr = rng.permutation(ys[p])
            per[i] = np.cumsum(yr)[:kmax] / ks
        rand[t] = per.mean(0)
    r_lo, r_md, r_hi = np.percentile(rand, [5, 50, 95], axis=0)

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.fill_between(ks, r_lo, r_hi, color=C_REF, alpha=0.22, lw=0, zorder=1)    # random 5-95% band
    ax.plot(ks, r_md, color=C_REF, lw=1.3, ls=(0, (5, 4)), zorder=2)           # random median
    ax.text(kmax * 0.62, r_hi[int(kmax * 0.62) - 1] + 0.015, "random",
            color=C_REF, fontsize=12, style="italic", zorder=2)
    for pat, prec in curves.items():
        ax.plot(ks, prec, color=C_PAT, lw=1.2, alpha=0.4, zorder=3)             # per-patient detector
    mean = np.mean(list(curves.values()), 0)
    ax.plot(ks, mean, color=C_MEAN, lw=3.4, zorder=5, solid_capstyle="round")   # detector cohort mean

    ax.set_xlim(1, kmax)
    ax.set_ylim(0, 1.0)
    ax.set_xticks([1, 5, 10, 15, 20])
    ax.set_xlabel("shortlist size  $k$", fontsize=14)
    ax.set_ylabel("precision @ $k$", fontsize=14)
    _finish(fig, ax, FIGDIR / "fig_epi_precision_at_k")
    print(f"  prec@5 detector = {mean[4]:.2f}   random 5-95% @5 = [{r_lo[4]:.2f}, {r_hi[4]:.2f}]"
          f"   (~{mean[4]/r_md[4]:.0f}x median)")


def main():
    auc_by_band()
    detector_roc()
    precision_at_k()


if __name__ == "__main__":
    main()
