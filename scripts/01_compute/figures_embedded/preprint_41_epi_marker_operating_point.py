"""
preprint_41 — Epi |ImCoh| connectivity marker: the marker and where it lives (Figure 1).

ALL-CONTACTS framing (corrected 2026-06-23). |ImCoh| is zero-lag-immune (Nolte 2004), so
proximity cannot enter the connectivity: NO shaft holdout, NO spatial control. The marker is
judged only against the MATCHED-STRENGTH null, on all contacts (audit_132, leave-one-SOZ-
contact-out, strength-residualised). Titles are neutral; the result is read off the geometry.

  A  band specificity   — heat tau5 all-contacts AUC beyond strength, per band, each band's
                          own matched-strength null band: delta/beta/low-gamma clear (MULTIBAND),
                          alpha marginal
  B  diffusion-time tau  — delta per-patient curves: slow-tau plateau; cohort mean peaks at tau5
  C  operator comparison — slow heat-kernel wins; PPR/Katz/comm/strength trail
  D  per-patient         — delta heat tau5 all-contacts AUC, responders vs the 2 hub-patients

Source caches (data/audit/epi_marker_allcontacts/): allcontacts_per_patient.csv,
allcontacts_cohort.csv, allcontacts_nulls.csv.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

AC = ROOT / "data/audit/epi_marker_allcontacts"
OUT = ROOT / "data/preprint/figures/all_bands"
OUT.mkdir(parents=True, exist_ok=True)

HEAD = "heat_t5"
BAND = "delta"
BANDS = ["delta", "beta", "low_gamma", "alpha"]
RESPONDERS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
              "Pat_08", "Pat_13", "Pat_14"]
C_BAND = {"delta": "#08519c", "beta": "#6a51a3", "low_gamma": "#238b45",
          "alpha": "#969696"}
C_RESP = "#1b7837"
C_HUB = "#c0392b"
C_GREY = "#9e9e9e"
NULLBAND = "0.86"
STAR = "#d4a017"


def _nullband(nulls, band):
    r = nulls[(nulls.band == band) & (nulls.marker == HEAD)]
    if not len(r):
        return np.nan, np.nan
    return float(r.null_median_mean.iloc[0]), float(r.null_median_p95.iloc[0])


# ===========================================================================
def panel_band(ax, pp, mlc, nulls):
    rng = np.random.default_rng(0)
    ax.axhline(0.5, color="black", ls="--", lw=0.8, zorder=1)
    for x, b in enumerate(BANDS):
        lo, hi = _nullband(nulls, b)
        if np.isfinite(lo):
            ax.add_patch(Rectangle((x - 0.34, lo), 0.68, hi - lo, facecolor=NULLBAND,
                                   edgecolor="0.6", linewidth=0.4, zorder=0))
        ys = pp[(pp.band == b) & (pp.marker == HEAD)].auc_resid.dropna().to_numpy()
        xs = x + rng.uniform(-0.13, 0.13, size=len(ys))
        ax.scatter(xs, ys, s=24, c=C_BAND[b], edgecolor="k", linewidth=0.3,
                   alpha=0.9, zorder=3)
        med = mlc[(mlc.band == b) & (mlc.marker == HEAD)
                  & (mlc.metric == "auc_resid")].median_auc
        med = float(med.iloc[0]) if len(med) else float(np.median(ys))
        ax.plot([x - 0.28, x + 0.28], [med, med], color="k", lw=1.7, zorder=4)
        if b == BAND:
            ax.scatter([x], [med], marker="*", s=180, c=STAR, edgecolor="k",
                       linewidth=0.5, zorder=5)
    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS])
    ax.set_ylabel("all-contacts AUC beyond strength")
    ax.set_title("A   band specificity (heat $\\tau_5$)", loc="left", fontweight="bold")
    ax.set_ylim(0.30, 1.0)


def panel_tau(ax, pp, nulls):
    lo, hi = _nullband(nulls, BAND)
    taus = list(range(6))
    cols = [f"heat_t{t}" for t in taus]
    if np.isfinite(lo):
        ax.axhspan(lo, hi, color=NULLBAND, zorder=0)
    ax.axhline(0.5, color="black", ls="--", lw=0.8, zorder=1)
    piv = (pp[(pp.band == BAND) & (pp.marker.isin(cols))]
           .pivot_table(index="patient", columns="marker", values="auc_resid")[cols])
    for pat in piv.index:
        c = C_RESP if pat in RESPONDERS else C_HUB
        ax.plot(taus, piv.loc[pat].to_numpy(), "-", color=c, lw=1.0, alpha=0.5, zorder=2)
    mean = piv.mean().to_numpy()
    ax.plot(taus, mean, "-", color="k", lw=2.6, zorder=4)
    ax.scatter(taus, mean, s=26, color="k", zorder=4)
    ax.scatter([5], [mean[5]], marker="*", s=210, c=STAR, edgecolor="k",
               linewidth=0.5, zorder=5)
    ax.set_xticks(taus)
    ax.set_xticklabels([f"$\\tau_{t}$" for t in taus])
    ax.set_xlabel("diffusion time  (fast $\\rightarrow$ slow)")
    ax.set_ylabel("all-contacts AUC beyond strength")
    ax.set_title("B   diffusion-time dependence ($\\delta$)", loc="left", fontweight="bold")
    ax.set_ylim(0.30, 1.0)
    handles = [Line2D([], [], color="k", lw=2.6, label="cohort mean"),
               Line2D([], [], color=C_RESP, lw=1.2, alpha=0.7, label="responder"),
               Line2D([], [], color=C_HUB, lw=1.2, alpha=0.7, label="hub-patient")]
    ax.legend(handles=handles, frameon=False, fontsize=7, loc="lower right",
              handlelength=1.4)


OPS = [("heat_t5", "heat ($\\tau_5$)"), ("heatN_t3", "heat$_N$ ($\\tau_3$)"),
       ("ppr_a95", "PPR ($\\alpha$.95)"), ("comm", "communicab."),
       ("diffdist_t2", "diffusion dist."), ("negres", "neg-resist."),
       ("katz", "Katz"), ("ppr_a50", "PPR ($\\alpha$.50)"),
       ("strength_baseline", "node strength")]


def panel_ops(ax, mlc, nulls):
    lo, hi = _nullband(nulls, BAND)
    rows = []
    for mk, lab in OPS:
        m = mlc[(mlc.band == BAND) & (mlc.marker == mk)
                & (mlc.metric == "auc_resid")].median_auc
        if len(m) and np.isfinite(m.iloc[0]):
            rows.append((float(m.iloc[0]), lab, mk))
    rows.sort(key=lambda r: r[0])
    y = np.arange(len(rows))
    if np.isfinite(lo):
        ax.axvspan(lo, hi, color=NULLBAND, zorder=0)
    ax.axvline(0.5, color="black", ls="--", lw=0.8, zorder=1)
    for yi, (v, lab, mk) in zip(y, rows):
        col = C_RESP if mk == HEAD else C_GREY
        ax.barh(yi, v, height=0.66, color=col, edgecolor="k", linewidth=0.4,
                zorder=3, alpha=0.9)
        if mk == HEAD:
            ax.scatter([v + 0.012], [yi], marker="*", s=130, c=STAR, edgecolor="k",
                       linewidth=0.5, zorder=5)
    ax.set_yticks(y)
    ax.set_yticklabels([r[1] for r in rows], fontsize=7.5)
    ax.set_xlabel("all-contacts AUC beyond strength  (strength $=0.5$)")
    ax.set_title("C   operator comparison ($\\delta$)", loc="left", fontweight="bold")
    ax.set_xlim(0.40, 0.85)


def panel_perpatient(ax, pp, nulls):
    lo, hi = _nullband(nulls, BAND)
    d = (pp[(pp.band == BAND) & (pp.marker == HEAD)]
         .dropna(subset=["auc_resid"]).sort_values("auc_resid"))
    med = float(d.auc_resid.median())
    y = np.arange(len(d))
    if np.isfinite(lo):
        ax.axvspan(lo, hi, color=NULLBAND, zorder=0)
    ax.axvline(0.5, color="black", ls="--", lw=0.9, zorder=2)
    ax.axvline(med, color=STAR, lw=1.6, zorder=2)
    for yi, (_, r) in zip(y, d.iterrows()):
        c = C_RESP if r.patient in RESPONDERS else C_HUB
        ax.barh(yi, r.auc_resid - 0.5, left=0.5, height=0.70, color=c,
                edgecolor="k", linewidth=0.4, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([p.replace("Pat_", "") for p in d.patient])
    ax.set_xlabel("all-contacts AUC beyond strength")
    ax.set_ylabel("patient")
    ax.set_title("D   per-patient ($\\delta$, heat $\\tau_5$)", loc="left", fontweight="bold")
    ax.set_xlim(0.15, 1.0)
    ax.set_ylim(-0.7, len(d) - 0.3)


def main():
    pp = pd.read_csv(AC / "allcontacts_per_patient.csv")
    mlc = pd.read_csv(AC / "allcontacts_cohort.csv")
    nulls = pd.read_csv(AC / "allcontacts_nulls.csv")

    fig = plt.figure(figsize=(11.0, 7.8))
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.30)
    panel_band(fig.add_subplot(gs[0, 0]), pp, mlc, nulls)
    panel_tau(fig.add_subplot(gs[0, 1]), pp, nulls)
    panel_ops(fig.add_subplot(gs[1, 0]), mlc, nulls)
    panel_perpatient(fig.add_subplot(gs[1, 1]), pp, nulls)
    handles = [Line2D([], [], marker="s", ls="", color=C_RESP, markeredgecolor="k",
                      markersize=8, label="responder (8/10)"),
               Line2D([], [], marker="s", ls="", color=C_HUB, markeredgecolor="k",
                      markersize=8, label="hub-patient (2/10)")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.01),
               ncol=2, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    out = OUT / "fig_epi_marker_connectivity.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
