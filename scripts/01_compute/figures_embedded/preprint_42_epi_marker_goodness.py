"""
preprint_42 — Epi |ImCoh| connectivity marker: beyond strength + per-patient (Figure 2).

ALL-CONTACTS framing (corrected 2026-06-23): matched-strength null only, no shaft holdout.
Shows the marker is not hubness and is consistent per patient. Neutral titles.

  A  strength decomposition — heat tau5 all-contacts AUC: raw (with strength) vs residual
                              (beyond strength); the gap is the strength share. delta/low-gamma
                              keep almost all of it; beta/alpha lose more
  B  propagator vs strength  — per-patient raw AUC: responders above identity; Pat_15 sits far
                              below (strength AUC 0.93, propagator 0.61 -> its SOZ is the hub)
  C  matched-strength verdict — cohort-median residual AUC (star) vs each band's matched-strength
                              null box: delta/beta/low-gamma clear, alpha marginal
  D  per-patient excess       — residual AUC minus each patient's OWN strength-matched null:
                              positive for responders, negative for the 2 hub-patients

Source caches (data/audit/epi_marker_allcontacts/).
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
C_RAW = "#bdbdbd"
NULLBAND = "0.82"
STAR = "#d4a017"


def _col(pat):
    return C_RESP if pat in RESPONDERS else C_HUB


def panel_decomp(ax, mlc):
    ax.axvline(0.5, color="black", ls="--", lw=0.8, zorder=1)
    for y, b in enumerate(BANDS):
        raw = mlc[(mlc.band == b) & (mlc.marker == HEAD) & (mlc.metric == "auc_raw")].median_auc
        res = mlc[(mlc.band == b) & (mlc.marker == HEAD) & (mlc.metric == "auc_resid")].median_auc
        if not (len(raw) and len(res)):
            continue
        raw, res = float(raw.iloc[0]), float(res.iloc[0])
        ax.plot([res, raw], [y, y], color="0.6", lw=1.4, zorder=2)
        ax.scatter([raw], [y], s=42, facecolor="white", edgecolor=C_BAND[b],
                   linewidth=1.4, zorder=3)
        ax.scatter([res], [y], s=58, color=C_BAND[b], edgecolor="k", linewidth=0.4,
                   zorder=4)
    ax.set_yticks(range(len(BANDS)))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS])
    ax.set_xlabel("all-contacts AUC")
    ax.set_title("A   strength decomposition (heat $\\tau_5$)", loc="left", fontweight="bold")
    ax.set_xlim(0.45, 0.95)
    ax.set_ylim(-0.6, len(BANDS) - 0.4)
    handles = [Line2D([], [], marker="o", ls="", markerfacecolor="white",
                      markeredgecolor="0.3", markersize=7, label="raw (with strength)"),
               Line2D([], [], marker="o", ls="", color="0.3", markeredgecolor="k",
                      markersize=7, label="beyond strength")]
    ax.legend(handles=handles, frameon=False, fontsize=7, loc="lower right")


def panel_twopop(ax, pp):
    piv = pp[pp.band == BAND].pivot_table(index="patient", columns="marker",
                                          values="auc_raw")
    x = piv["strength_baseline"]; y = piv[HEAD]
    lo = float(min(x.min(), y.min())) - 0.04
    hi = float(max(x.max(), y.max())) + 0.04
    ax.fill_between([lo, hi], [lo, hi], hi, color=C_RESP, alpha=0.07, zorder=0)
    ax.plot([lo, hi], [lo, hi], color="k", lw=0.7, ls="--", alpha=0.6, zorder=1)
    for pat in piv.index:
        ax.scatter([x[pat]], [y[pat]], s=54, c=_col(pat), edgecolor="k",
                   linewidth=0.4, zorder=3)
        ax.text(x[pat], y[pat], pat.replace("Pat_", ""), fontsize=5.2, va="bottom",
                ha="left", alpha=0.85)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel("node-strength AUC (raw)")
    ax.set_ylabel("propagator AUC (raw)")
    ax.set_title("B   propagator vs node strength ($\\delta$)", loc="left", fontweight="bold")


def panel_verdict(ax, mlc, nulls):
    ax.axvline(0.5, color="black", ls="--", lw=0.8, zorder=1)
    for y, b in enumerate(BANDS):
        nr = nulls[(nulls.band == b) & (nulls.marker == HEAD)]
        if len(nr):
            lo = float(nr.null_median_mean.iloc[0]); hi = float(nr.null_median_p95.iloc[0])
            ax.add_patch(Rectangle((lo, y - 0.26), hi - lo, 0.52, facecolor=NULLBAND,
                                   edgecolor="0.55", linewidth=0.5, zorder=2))
        obs = mlc[(mlc.band == b) & (mlc.marker == HEAD)
                  & (mlc.metric == "auc_resid")].median_auc
        if len(obs):
            ax.scatter([float(obs.iloc[0])], [y], marker="*", s=150, c=C_BAND[b],
                       edgecolor="k", linewidth=0.5, zorder=4)
    ax.set_yticks(range(len(BANDS)))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS])
    ax.set_xlabel("all-contacts AUC beyond strength")
    ax.set_title("C   matched-strength verdict", loc="left", fontweight="bold")
    ax.set_xlim(0.40, 0.95)
    ax.set_ylim(-0.6, len(BANDS) - 0.4)
    handles = [Line2D([], [], marker="*", ls="", color="0.3", markeredgecolor="k",
                      markersize=11, label="cohort median"),
               Line2D([], [], marker="s", ls="", color=NULLBAND, markeredgecolor="0.55",
                      markersize=9, label="matched-strength null")]
    ax.legend(handles=handles, frameon=False, fontsize=7, loc="lower right")


def panel_excess(ax, sn):
    rng = np.random.default_rng(3)
    ax.axhline(0.0, color="black", lw=0.8, zorder=1)
    for x, b in enumerate(BANDS):
        d = sn[sn.band == b].dropna(subset=["obs_resid", "null_p50"])
        exc = (d.obs_resid - d.null_p50).to_numpy()
        cols = [_col(p) for p in d.patient]
        xs = x + rng.uniform(-0.14, 0.14, size=len(exc))
        ax.scatter(xs, exc, s=34, c=cols, edgecolor="k", linewidth=0.3, zorder=3)
        ax.plot([x - 0.26, x + 0.26], [np.median(exc)] * 2, color="k", lw=1.5, zorder=4)
    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS])
    ax.set_ylabel("AUC $-$ own matched-strength null")
    ax.set_title("D   per-patient excess over own null", loc="left", fontweight="bold")


def main():
    pp = pd.read_csv(AC / "allcontacts_per_patient.csv")
    mlc = pd.read_csv(AC / "allcontacts_cohort.csv")
    nulls = pd.read_csv(AC / "allcontacts_nulls.csv")
    sn = pd.read_csv(AC / "allcontacts_stratnull_per_patient.csv")

    fig = plt.figure(figsize=(11.0, 7.9))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.30)
    panel_decomp(fig.add_subplot(gs[0, 0]), mlc)
    panel_twopop(fig.add_subplot(gs[0, 1]), pp)
    panel_verdict(fig.add_subplot(gs[1, 0]), mlc, nulls)
    panel_excess(fig.add_subplot(gs[1, 1]), sn)
    handles = [Line2D([], [], marker="s", ls="", color=C_RESP, markeredgecolor="k",
                      markersize=8, label="responder (8/10)"),
               Line2D([], [], marker="s", ls="", color=C_HUB, markeredgecolor="k",
                      markersize=8, label="hub-patient (2/10)")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.01),
               ncol=2, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    out = OUT / "fig_epi_marker_beyond_strength.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
