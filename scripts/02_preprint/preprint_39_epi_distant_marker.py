#!/usr/bin/env python3
"""preprint_39 — the verified delta distant-SOZ marker (LRG heat-kernel propagator).

The LRG heat-kernel propagator rho(tau)=e^{-tau L}/Z at slow diffusion time, in the delta
band, given known SOZ on some electrode shafts, ranks SOZ on OTHER shafts above healthy
contacts -- the non-trivial *distant* discovery where spatial proximity is useless by
construction. All four panels read the audit_101/102 result CSVs (read-only):

(a) Per-patient leave-one-shaft-out distant-discovery AUC (delta, slow heat kernel heat_t5),
    strength-residualised; the label-shuffle null band + chance line. 8/10 clear the null;
    the two that fail are the hub-patients (Pat_10, Pat_15) whose SOZ ARE the network hubs.
(b) Slow-diffusion mechanism: cohort-median AUC vs diffusion time tau for the heat kernel
    e^{-tau L} and its normalized-Laplacian variant -- distant discovery needs SLOW tau
    (fast tau sits at/below chance); the best PageRank operator is shown for reference.
(c) Operator comparison: cohort-median delta AUC, best member of each Laplacian operator
    family -- the heat kernel (the LRG propagator) beats PageRank / communicability / Katz /
    effective-resistance / diffusion-distance, and the node-strength baseline sits at chance.
(d) Band specificity: the SAME fixed marker (heat_t5) across delta/alpha/beta/low-gamma --
    only delta clears the null band.

Verification (label-shuffle null p=0.000; nested leave-one-patient-out held-out median 0.716)
is printed to stdout and shown only as reference lines / a shaded null band (house style:
no in-axes stat text). PDF only, full vector, transparent.

Inputs (read-only) data/audit/epi_marker_library/:
    marker_library_per_patient.csv, marker_library_cohort.csv, verify_nulls.csv, verify_lopo.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

ML = ROOT / "data/audit/epi_marker_library"
OUT = ROOT / "data/preprint/figures/all_bands"
OUT.mkdir(parents=True, exist_ok=True)

BAND = "delta"
HEAD = "heat_t5"                       # headline marker: slowest tau, selected 9/10 LOPO folds
BANDS = ["delta", "alpha", "beta", "low_gamma"]
BAND_C = {"delta": "#1b7837", "beta": "#2166ac", "low_gamma": "#762a83", "alpha": "#b35806"}
C_HEAT = "#1b7837"      # heat kernel (the propagator)
C_HEATN = "#7fbf7b"     # normalized-Laplacian heat
C_OTHER = "#9e9e9e"     # non-heat operators
C_BASE = "#b2182b"      # strength baseline
C_FAIL = "#d6604d"      # hub-patient (marker fails)

# operator families: (display name, member markers, colour)
FAMILIES = [
    ("heat kernel", [f"heat_t{i}" for i in range(6)], C_HEAT),
    ("norm.-Lap. heat", [f"heatN_t{i}" for i in range(6)], C_HEATN),
    ("PageRank", ["ppr_a50", "ppr_a85", "ppr_a95"], C_OTHER),
    ("communicability", ["comm"], C_OTHER),
    ("diffusion-dist.", ["diffdist_t2"], C_OTHER),
    ("neg. resistance", ["negres"], C_OTHER),
    ("Katz", ["katz"], C_OTHER),
    ("strength (baseline)", ["strength_baseline"], C_BASE),
]
TAU = np.geomspace(1.0, 10.0, 6)       # diffusion time in units of 1/lambda_max (t0..t5)


def _tex(b):
    return BRAIN_BAND_TEX_DICT.get(b, b)


def _null(nulls, marker=HEAD):
    r = nulls[nulls.marker == marker].iloc[0]
    return float(r.null_median_mean), float(r.null_median_p95)


# ----------------------------------------------------------------------------- panels
def panel_a(ax, pp, nulls):
    d = pp[(pp.band == BAND) & (pp.marker == HEAD)].dropna(subset=["auc_resid"])
    d = d.sort_values("auc_resid").reset_index(drop=True)
    lo, hi = _null(nulls)
    ax.axhspan(lo, hi, color="0.86", zorder=0)
    ax.axhline(hi, color="0.55", ls=":", lw=0.8, zorder=1)
    ax.axhline(0.5, color="black", ls="--", lw=0.8, zorder=1)
    x = np.arange(len(d))
    col = [C_HEAT if v >= hi else C_FAIL for v in d.auc_resid]
    ax.scatter(x, d.auc_resid, c=col, s=72, edgecolor="black", linewidth=0.5, zorder=3)
    med = float(d.auc_resid.median())
    ax.axhline(med, color=C_HEAT, ls="-", lw=1.1, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels([p.replace("Pat_", "") for p in d.patient], fontsize=6.6)
    ax.set_xlabel("patient")
    ax.set_ylabel("distant-discovery AUC\n(off-shaft, strength-residual)")
    ax.set_ylim(0.12, 1.0)
    ax.text(0.02, 1.03, f"a  per-patient ({_tex(BAND)}, slow heat kernel)",
            transform=ax.transAxes, va="bottom", fontsize=8.3)
    return med, int((d.auc_resid >= hi).sum()), len(d)


def panel_b(ax, coh, nulls):
    cd = coh[(coh.band == BAND) & (coh.metric == "auc_resid")]
    heat = [float(cd[cd.marker == f"heat_t{i}"].median_auc.iloc[0]) for i in range(6)]
    heatN = [float(cd[cd.marker == f"heatN_t{i}"].median_auc.iloc[0]) for i in range(6)]
    ppr = float(cd[cd.marker == "ppr_a95"].median_auc.iloc[0])
    lo, hi = _null(nulls)
    ax.axhspan(lo, hi, color="0.86", zorder=0)
    ax.axhline(0.5, color="black", ls="--", lw=0.8, zorder=1)
    ax.plot(TAU, heat, "-o", color=C_HEAT, lw=2.1, ms=5.5,
            label=r"heat kernel $e^{-\tau L}$", zorder=3)
    ax.plot(TAU, heatN, "-s", color=C_HEATN, lw=1.7, ms=4.5,
            label="norm.-Laplacian heat", zorder=3)
    ax.axhline(ppr, color="0.45", ls="-.", lw=1.0, label="best PageRank (ref.)", zorder=2)
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 5, 10])
    ax.set_xticklabels(["1", "2", "5", "10"])
    ax.set_xlabel(r"diffusion time $\tau$  ($1/\lambda_{\max}$ units; fast $\to$ slow)")
    ax.set_ylabel("cohort-median AUC")
    ax.set_ylim(0.40, 0.80)
    ax.legend(loc="lower right", frameon=False, fontsize=6.6)
    ax.text(0.02, 1.03, "b  slow diffusion reaches far", transform=ax.transAxes,
            va="bottom", fontsize=8.3)


def panel_c(ax, coh, nulls):
    cd = coh[(coh.band == BAND) & (coh.metric == "auc_resid")]
    rows = []
    for name, members, color in FAMILIES:
        if name.startswith("strength"):
            v = 0.5                                      # residual on itself = chance
        else:
            vals = cd[cd.marker.isin(members)].median_auc.dropna()
            v = float(vals.max()) if len(vals) else np.nan
        rows.append((name, v, color))
    rows.sort(key=lambda r: r[1])
    lo, hi = _null(nulls)
    ax.axvspan(lo, hi, color="0.86", zorder=0)
    ax.axvline(0.5, color="black", ls="--", lw=0.8, zorder=1)
    y = np.arange(len(rows))
    ax.barh(y, [r[1] for r in rows], color=[r[2] for r in rows],
            edgecolor="black", linewidth=0.4, height=0.62, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.8)
    ax.set_xlim(0.38, 0.80)
    ax.set_xlabel("cohort-median AUC (best of family)")
    ax.text(0.02, 1.03, f"c  operator comparison ({_tex(BAND)})", transform=ax.transAxes,
            va="bottom", fontsize=8.3)


def panel_d(ax, pp, nulls):
    lo, hi = _null(nulls)
    ax.axhspan(lo, hi, color="0.86", zorder=0)
    ax.axhline(hi, color="0.55", ls=":", lw=0.8, zorder=1)
    ax.axhline(0.5, color="black", ls="--", lw=0.8, zorder=1)
    meds = []
    for i, b in enumerate(BANDS):
        d = pp[(pp.band == b) & (pp.marker == HEAD)].auc_resid.dropna().to_numpy()
        jit = np.linspace(-0.14, 0.14, len(d)) if len(d) > 1 else np.zeros(len(d))
        ax.scatter(np.full(len(d), i) + jit, d, s=24, color=BAND_C[b],
                   alpha=0.55, edgecolor="none", zorder=2)
        m = float(np.median(d))
        meds.append((b, m))
        ax.scatter(i, m, s=95, color=BAND_C[b], edgecolor="black", linewidth=0.6,
                   marker="D", zorder=3)
    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([_tex(b) for b in BANDS])
    ax.set_xlabel("band (same fixed marker)")
    ax.set_ylabel("distant-discovery AUC")
    ax.set_ylim(0.12, 1.0)
    ax.text(0.02, 1.03, "d  band specificity", transform=ax.transAxes,
            va="bottom", fontsize=8.3)
    return meds


def main():
    pp = pd.read_csv(ML / "marker_library_per_patient.csv")
    coh = pd.read_csv(ML / "marker_library_cohort.csv")
    nulls = pd.read_csv(ML / "verify_nulls.csv")
    lopo = pd.read_csv(ML / "verify_lopo.csv")

    fig = plt.figure(figsize=(11.0, 7.6))
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.30)
    med, nwin, ntot = panel_a(fig.add_subplot(gs[0, 0]), pp, nulls)
    panel_b(fig.add_subplot(gs[0, 1]), coh, nulls)
    panel_c(fig.add_subplot(gs[1, 0]), coh, nulls)
    band_meds = panel_d(fig.add_subplot(gs[1, 1]), pp, nulls)

    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=C_HEAT,
               markeredgecolor="black", markersize=8, label=r"recovered SOZ ($\geq$ null)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=C_FAIL,
               markeredgecolor="black", markersize=8, label="hub-patient (marker fails)"),
        Line2D([0], [0], color="black", ls="--", lw=0.9, label="chance (AUC 0.5)"),
        Patch(facecolor="0.86", edgecolor="0.55", label="label-shuffle null (mean -> 95th pct)"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, -0.015),
               ncol=4, frameon=False, fontsize=7.6)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    out = OUT / "fig_epi_distant_soz_marker.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)

    # ---- stdout: every number behind the figure (house style: stats to stdout) ----
    print(f"[preprint_39] -> {out.name}")
    print(f"\n(a) per-patient {BAND} {HEAD} strength-residual LOSO distant-discovery AUC:")
    d = pp[(pp.band == BAND) & (pp.marker == HEAD)].sort_values("auc_resid", ascending=False)
    for _, r in d.iterrows():
        print(f"      {r.patient:7s} {r.auc_resid:.3f}")
    print(f"    median={med:.3f}  recovered(>=null p95)={nwin}/{ntot}")
    print("\n(b/c) verification (audit_102):")
    for _, r in nulls.iterrows():
        print(f"      label-shuffle null {r.marker:9s}: real={r.real_median_auc:.3f} "
              f"null~{r.null_median_mean:.3f} (p95={r.null_median_p95:.3f}) p_emp={r.p_empirical:.3f}")
    print(f"      nested-LOPO held-out: median={lopo.held_auc.median():.3f} "
          f"mean={lopo.held_auc.mean():.3f} {int((lopo.held_auc>0.5).sum())}/{len(lopo)}>0.5 "
          f"(selected {lopo.selected_marker.value_counts().to_dict()})")
    print(f"\n(d) {HEAD} cohort-median AUC by band:")
    for b, m in band_meds:
        print(f"      {b:10s} {m:.3f}")


if __name__ == "__main__":
    main()
