#!/usr/bin/env python3
"""Audit round-3 Section-5 redos — CTM headline.

Split out of audit_round3_section5_redo.py on 2026-05-29 (Phase 4-B
split 2/7). Imports shared constants + helpers from
_audit_round3_shared and runs: redo4.
"""
from __future__ import annotations

import matplotlib  # noqa: F401
import matplotlib.lines as mlines  # noqa: F401
import matplotlib.patches as mpatches  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram  # noqa: F401

from _audit_round3_shared import (
    BAND_ORDER, BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, DEND_LABEL_DROP_HEIGHTS,
    DEND_LABEL_FONTSIZE, DEND_LABEL_Y_FRAC, FIG, OUT, PATIENTS_LIST,
    PHASES, R2, R2_TBL, ROOT, S5, TBL, asterisks, kc_vectors,
    load_lrg_result,
)


# ===================================================================
# Redo 4 — CTM headline (per-patient lines + asterisks above split)
# ===================================================================
def redo4_ctm_headline() -> None:
    CTM = pd.read_csv(S5 / "05_ctm_sigma_aggregate/tables/cohort_summary.csv")
    CTM_PT = pd.read_csv(S5 / "05_ctm_sigma_aggregate/tables/Td_per_patient_per_band.csv")

    fig, ax = plt.subplots(figsize=(13, 5.0))
    width = 0.25
    x_band = np.arange(len(BAND_ORDER))
    colors = {"split": "#1f77b4", "drift": "#888888", "xprobe": "#e07b00"}
    box_keys = (("split", "rho_split"), ("drift", "rho_null_drift"),
                ("xprobe", "rho_split_cross_probe"))

    # Boxes
    for off, (name, key) in zip((-width, 0.0, width), box_keys):
        data = [CTM_PT[CTM_PT["band"] == b][key].dropna().values for b in BAND_ORDER]
        ax.boxplot(
            data, positions=x_band + off, widths=width * 0.85,
            patch_artist=True, showmeans=False, showfliers=False,
            boxprops=dict(facecolor=colors[name], alpha=0.50, edgecolor="#444"),
            medianprops=dict(color="#111", lw=1.4),
            whiskerprops=dict(color="#666"), capprops=dict(color="#666"),
        )

    # Per-patient connected dots: split → drift → xprobe
    pats = sorted(CTM_PT["patient"].unique())
    for xi, b in enumerate(BAND_ORDER):
        sub = CTM_PT[CTM_PT["band"] == b].set_index("patient")
        for p in pats:
            if p not in sub.index:
                continue
            ys = [
                float(sub.loc[p, "rho_split"]),
                float(sub.loc[p, "rho_null_drift"]),
                float(sub.loc[p, "rho_split_cross_probe"]),
            ]
            xs = [xi - width, xi, xi + width]
            is_pro = ys[0] > 0
            line_c = "#1a7c3e" if is_pro else "#c0392b"
            ax.plot(xs, ys, color=line_c, alpha=0.45, lw=0.8, zorder=3)
            ax.scatter(xs, ys, s=10, color=line_c, alpha=0.65,
                       edgecolor="white", linewidth=0.3, zorder=4)

    # Significance asterisks ABOVE the split (blue) box at top of each box
    for xi, b in enumerate(BAND_ORDER):
        rs = CTM[CTM["band"] == b].iloc[0]
        p = float(rs["wilcoxon_split_gt_drift_p"])
        # Anchor above the IQR of the split box
        v = CTM_PT[CTM_PT["band"] == b]["rho_split"].dropna().values
        if not len(v):
            continue
        q3 = float(np.nanpercentile(v, 75))
        whisk_top = float(np.nanmax(v))
        anchor = max(q3, whisk_top) + 0.03
        stars = asterisks(p)
        txt = f"p={p:.3f}{stars}"
        col = "#c0392b" if p <= 0.05 else "#444"
        ax.text(xi - width, anchor, txt, ha="center", va="bottom",
                fontsize=8, color=col, fontweight="bold" if stars else "normal")

    ax.axhline(0, color="0.3", lw=0.7, ls="--")
    ax.set_xticks(x_band)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER], fontsize=11)
    ax.set_ylabel(r"$\rho$ (CTM split / drift / x-probe)", fontsize=10)

    # y-limits: tighten using p1/p99 of split + xprobe (drift is small)
    all_v = np.concatenate([
        CTM_PT[k].dropna().values for _, k in box_keys
    ])
    qmin = float(np.nanpercentile(all_v, 1))
    qmax = float(np.nanpercentile(all_v, 99))
    ax.set_ylim(qmin - 0.15, qmax + 0.18)

    handles = [
        mpatches.Patch(facecolor=colors["split"], alpha=0.50, edgecolor="#444",
                       label=r"$\rho_{\mathrm{split}}$ (controlled)"),
        mpatches.Patch(facecolor=colors["drift"], alpha=0.50, edgecolor="#444",
                       label=r"$\rho_{\mathrm{null\,drift}}$"),
        mpatches.Patch(facecolor=colors["xprobe"], alpha=0.50, edgecolor="#444",
                       label=r"$\rho_{\mathrm{xprobe}}$"),
        mlines.Line2D([], [], color="#1a7c3e", marker="o", lw=0.8, alpha=0.7,
                      markersize=4, label="patient: pro (ρ_split > 0)"),
        mlines.Line2D([], [], color="#c0392b", marker="o", lw=0.8, alpha=0.7,
                      markersize=4, label="patient: anti (ρ_split ≤ 0)"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.03),
               ncol=5, frameon=False, fontsize=8.5)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(FIG / "ctm_headline.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[redo4] wrote ctm_headline.pdf")


def main() -> None:
    print("[redo] 4 — ctm_headline")
    redo4_ctm_headline()
    print(f"[ctm] done · outputs at {OUT}")


if __name__ == "__main__":
    main()
