"""Build the figures that are missing from the raw-FC audit folder.

1. ``z_inflation_diagnostic.pdf``
   Per band × distance × phase pair: scatter of `d_obs` (y) vs
   `null_median` (y′) per patient, with the within-patient null IQR as
   error bars. Visualises why Z = 12-685 in audit_25 is a tiny-MAD
   artefact rather than a huge-distance finding.

2. ``Td_dS_vs_dF_scatter.pdf``
   Per band: per-patient (T_d^(d_S), T_d^(d_F)) scatter with identity
   line. Shows the structural-shift verdict (points cluster on the
   identity line) versus an amplitude-only confound (points would be on
   the y-axis).

3. ``Td_swarm_per_band.pdf``
   Per band: strip plot of per-patient T_d^(d_S) and T_d^(d_F) with
   median markers, sign-count annotation, and Pat_03 highlighted.
   Shows the cohort dispersion behind the headline median T_d.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE
from lrg_eegfc.config.paths import DATA_ROOT
from lrg_eegfc.visuals.layout import figure_legend, add_provenance_footer

OUT_BASE = DATA_ROOT / "audit" / "raw_fc_phase_distance"
N_COHORT = len(PATIENTS_4PHASE)  # single source of truth — never hardcode

# Reusable colour map per phase pair
PAIR_COLOR = {
    ("rest_pre", "task_test"): "tab:red",
    ("rest_pre", "rest_post"): "tab:blue",
    ("task_test", "rest_post"): "tab:green",
}
PAIR_LABEL = {
    ("rest_pre", "task_test"): "RPre→TT",
    ("rest_pre", "rest_post"): "RPre→RPost",
    ("task_test", "rest_post"): "TT→RPost",
}


# ---------------------------------------------------------------------
# Figure 1 — Z-inflation diagnostic
# ---------------------------------------------------------------------
scale = pd.read_csv(OUT_BASE / "per_patient_with_scale.csv")
out_pdf = OUT_BASE / "z_inflation_diagnostic.pdf"
with PdfPages(out_pdf) as pdf:
    nb = len(BRAIN_BANDS_NAMES)
    nrows = 3  # one row per distance
    fig, axes = plt.subplots(
        nrows, nb, figsize=(2.4 * nb + 0.6, 2.0 * nrows + 0.4),
        squeeze=False,
    )
    distances = ("P", "S", "F")
    for r, d_label in enumerate(distances):
        for c, band in enumerate(BRAIN_BANDS_NAMES):
            ax = axes[r, c]
            sub = scale[(scale.band == band) & (scale.distance == d_label)]
            if sub.empty:
                ax.set_xticks([]); ax.set_yticks([]); continue
            patients = sorted(sub.patient.unique())
            x_pos = np.arange(len(patients))
            # Plot null median ± IQR as a grey ribbon per patient
            for i, p in enumerate(patients):
                rows = sub[sub.patient == p]
                if rows.empty:
                    continue
                nm = rows.iloc[0].null_median
                q1 = rows.iloc[0].null_q1
                q3 = rows.iloc[0].null_q3
                ax.add_patch(plt.Rectangle((i - 0.4, q1), 0.8, q3 - q1,
                                            facecolor="lightgray", edgecolor="k",
                                            linewidth=0.3, zorder=1))
                ax.plot([i - 0.4, i + 0.4], [nm, nm],
                        color="k", lw=0.7, zorder=2)
                # observed cross-phase distances per pair
                for pair_key, color in PAIR_COLOR.items():
                    phi_A, phi_B = pair_key
                    row = rows[(rows.phase_A == phi_A) & (rows.phase_B == phi_B)]
                    if not row.empty:
                        ax.scatter(
                            i, row.iloc[0].d_obs, color=color,
                            edgecolor="k", linewidth=0.3, s=18, zorder=3,
                        )
            ax.set_xticks(x_pos)
            ax.set_xticklabels([p.replace("Pat_", "") for p in patients],
                               fontsize=5, rotation=45)
            ax.tick_params(axis="y", labelsize=6)
            if r == 0:
                ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
            if c == 0:
                ax.set_ylabel(f"d_{d_label}\nnull IQR + obs",
                              fontsize=8)
    # Figure-level shared legend (NEVER per-axis for shared content) —
    # see .agents/guides/05_plotting/legends.md
    handles = [
        plt.Line2D([], [], marker="o", linestyle="",
                   color=PAIR_COLOR[k],
                   markeredgecolor="k", markersize=6,
                   label=PAIR_LABEL[k])
        for k in PAIR_COLOR
    ]
    handles.append(plt.Rectangle((0, 0), 1, 1,
                                  facecolor="lightgray",
                                  edgecolor="k", label="null IQR"))
    figure_legend(fig, handles, where="bottom", fontsize=8)
    add_provenance_footer(fig, "z_inflation_diagnostic — null IQR vs observed")
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    pdf.savefig(fig, dpi=200); plt.close(fig)
print(f"wrote {out_pdf}")


# ---------------------------------------------------------------------
# Figure 2 — d_S vs d_F triangle scatter per band
# ---------------------------------------------------------------------
td = pd.read_csv(OUT_BASE / "Td_per_patient_per_band.csv")
out_pdf = OUT_BASE / "Td_dS_vs_dF_scatter.pdf"
with PdfPages(out_pdf) as pdf:
    nb = len(BRAIN_BANDS_NAMES)
    fig, axes = plt.subplots(1, nb, figsize=(2.2 * nb + 0.6, 2.6),
                              squeeze=False)
    for c, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes[0, c]
        sub = td[td.band == band].dropna(subset=["S", "F"])
        if sub.empty:
            continue
        in_pool = sub  # n=10; Pat_03 included
        in_pool_excl = sub[sub.patient != "Pat_03"]  # for blue circles
        out_p = sub[sub.patient == "Pat_03"]
        # Identity
        lim_lo = float(min(sub.S.min(), sub.F.min())) * 1.1
        lim_hi = float(max(sub.S.max(), sub.F.max())) * 1.1
        ax.axhline(0, color="gray", lw=0.4, alpha=0.6)
        ax.axvline(0, color="gray", lw=0.4, alpha=0.6)
        ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi],
                color="k", lw=0.6, ls="--", alpha=0.6, label="identity")
        ax.scatter(in_pool_excl.S, in_pool_excl.F, c="tab:blue",
                   edgecolor="k", linewidth=0.3, s=32, zorder=3)
        if not out_p.empty:
            ax.scatter(out_p.S, out_p.F, c="tab:orange",
                       marker="^", edgecolor="k", linewidth=0.3, s=42,
                       zorder=4, label="Pat_03")
        # Patient labels (n=10 total)
        for _, r in in_pool_excl.iterrows():
            ax.text(r.S, r.F, r.patient.replace("Pat_", ""),
                    fontsize=4.5, va="bottom", ha="left", alpha=0.7)
        if not out_p.empty:
            for _, r in out_p.iterrows():
                ax.text(r.S, r.F, r.patient.replace("Pat_", ""),
                        fontsize=4.5, va="bottom", ha="left", alpha=0.7,
                        color="tab:orange")
        ax.set_xlim(lim_lo, lim_hi); ax.set_ylim(lim_lo, lim_hi)
        ax.set_aspect("equal")
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
        ax.set_xlabel(r"$T_d^{(d_S)}$", fontsize=8)
        if c == 0:
            ax.set_ylabel(r"$T_d^{(d_F)}$", fontsize=9)
        ax.tick_params(labelsize=6)
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", color="tab:blue",
                   markeredgecolor="k", markersize=6, label="in-pool patients"),
        plt.Line2D([], [], marker="^", linestyle="", color="tab:orange",
                   markeredgecolor="k", markersize=7, label="Pat_03 (1024 Hz)"),
        plt.Line2D([], [], color="k", lw=0.6, ls="--", label="identity"),
    ]
    figure_legend(fig, handles, where="bottom", fontsize=8)
    add_provenance_footer(fig, "Td_dS_vs_dF_scatter — structural (on identity) vs amplitude (off)")
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    pdf.savefig(fig, dpi=200); plt.close(fig)
print(f"wrote {out_pdf}")


# ---------------------------------------------------------------------
# Figure 2b — T_d^(d_S) vs T_d^(d_P) per band  (orthogonal interpretable pair)
# ---------------------------------------------------------------------
out_pdf = OUT_BASE / "Td_dS_vs_dP_scatter.pdf"
with PdfPages(out_pdf) as pdf:
    nb = len(BRAIN_BANDS_NAMES)
    fig, axes = plt.subplots(1, nb, figsize=(2.2 * nb + 0.6, 2.6),
                              squeeze=False)
    for c, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes[0, c]
        sub = td[td.band == band].dropna(subset=["S", "P"])
        if sub.empty:
            continue
        in_pool = sub  # n=10; Pat_03 included
        in_pool_excl = sub[sub.patient != "Pat_03"]  # for blue circles
        out_p = sub[sub.patient == "Pat_03"]
        lim_lo = float(min(sub.S.min(), sub.P.min())) * 1.1
        lim_hi = float(max(sub.S.max(), sub.P.max())) * 1.1
        ax.axhline(0, color="gray", lw=0.4, alpha=0.6)
        ax.axvline(0, color="gray", lw=0.4, alpha=0.6)
        ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi],
                color="k", lw=0.6, ls="--", alpha=0.6, label="identity")
        ax.scatter(in_pool_excl.S, in_pool_excl.P, c="tab:blue",
                   edgecolor="k", linewidth=0.3, s=32, zorder=3)
        if not out_p.empty:
            ax.scatter(out_p.S, out_p.P, c="tab:orange",
                       marker="^", edgecolor="k", linewidth=0.3, s=42,
                       zorder=4, label="Pat_03")
        for _, r in in_pool_excl.iterrows():
            ax.text(r.S, r.P, r.patient.replace("Pat_", ""),
                    fontsize=4.5, va="bottom", ha="left", alpha=0.7)
        if not out_p.empty:
            for _, r in out_p.iterrows():
                ax.text(r.S, r.P, r.patient.replace("Pat_", ""),
                        fontsize=4.5, va="bottom", ha="left", alpha=0.7,
                        color="tab:orange")
        ax.set_xlim(lim_lo, lim_hi); ax.set_ylim(lim_lo, lim_hi)
        ax.set_aspect("equal")
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
        ax.set_xlabel(r"$T_d^{(d_S)}$  (topology)", fontsize=8)
        if c == 0:
            ax.set_ylabel(r"$T_d^{(d_P)}$  (volume + topology)", fontsize=9)
        ax.tick_params(labelsize=6)
    # Figure-level legend (shared across panels) — see 05_plotting/legends.md
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", color="tab:blue",
                   markeredgecolor="k", markersize=6, label="in-pool patients"),
        plt.Line2D([], [], marker="^", linestyle="", color="tab:orange",
                   markeredgecolor="k", markersize=7, label="Pat_03 (1024 Hz)"),
        plt.Line2D([], [], color="k", lw=0.6, ls="--", label="identity"),
    ]
    figure_legend(fig, handles, where="bottom", fontsize=8)
    add_provenance_footer(fig, "Td_dS_vs_dP_scatter — orthogonal interpretable pair")
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    pdf.savefig(fig, dpi=200); plt.close(fig)
print(f"wrote {out_pdf}")


# ---------------------------------------------------------------------
# Figure 3 — T_d swarm per band
# ---------------------------------------------------------------------
out_pdf = OUT_BASE / "Td_swarm_per_band.pdf"
with PdfPages(out_pdf) as pdf:
    nb = len(BRAIN_BANDS_NAMES)
    fig, axes = plt.subplots(1, nb, figsize=(2.4 * nb + 0.6, 2.7),
                              squeeze=False, sharey=True)
    for c, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes[0, c]
        sub = td[td.band == band].dropna(subset=["S", "P", "F"])
        if sub.empty:
            continue
        in_pool = sub  # n=10; Pat_03 included
        in_pool_excl = sub[sub.patient != "Pat_03"]  # for blue circles
        out_p = sub[sub.patient == "Pat_03"]
        # x positions: 0 = d_S, 1 = d_P, 2 = d_F
        rng = np.random.default_rng(42)
        for x, col in zip([0, 1, 2], ["S", "P", "F"]):
            ys_blue = in_pool_excl[col].to_numpy()
            ys_full = in_pool[col].to_numpy()  # n=10 for stats
            xs = x + rng.uniform(-0.12, 0.12, size=len(ys_blue))
            ax.scatter(xs, ys_blue, c="tab:blue", edgecolor="k",
                       linewidth=0.3, s=24, zorder=2)
            if not out_p.empty:
                ax.scatter([x], [out_p[col].iloc[0]],
                           c="tab:orange", marker="^",
                           edgecolor="k", linewidth=0.3, s=36, zorder=4)
            # Median bar over n=10
            med = float(np.median(ys_full))
            ax.plot([x - 0.25, x + 0.25], [med, med],
                    color="k", lw=1.2, zorder=3)
            # n_neg / n_tot annotation (n=10)
            n_neg = int((ys_full < 0).sum()); n_tot = len(ys_full)
            ax.text(x, ax.get_ylim()[1] if c == 0 else 0,
                    f"{n_neg}/{n_tot}", fontsize=6, ha="center", va="top",
                    bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                              edgecolor="gray", linewidth=0.3))
        ax.axhline(0, color="gray", lw=0.4, alpha=0.6)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["$d_S$\ntopo.", "$d_P$\nvol.+topo.", "$d_F$\namp."], fontsize=7)
        ax.set_xlim(-0.5, 2.5)
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
        if c == 0:
            ax.set_ylabel(r"$T_d = d(\mathrm{TT,RPost}) - d(\mathrm{RPre,TT})$",
                          fontsize=8)
        ax.tick_params(labelsize=6)
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", color="tab:blue",
                   markeredgecolor="k", markersize=6, label="in-pool patients"),
        plt.Line2D([], [], marker="^", linestyle="", color="tab:orange",
                   markeredgecolor="k", markersize=7, label="Pat_03 (1024 Hz)"),
        plt.Line2D([], [], color="k", lw=1.2, label=f"median (n={N_COHORT})"),
    ]
    figure_legend(fig, handles, where="bottom", fontsize=8)
    add_provenance_footer(fig, "Td_swarm — per-band cohort dispersion, T_d<0 sign count")
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    pdf.savefig(fig, dpi=200); plt.close(fig)
print(f"wrote {out_pdf}")
