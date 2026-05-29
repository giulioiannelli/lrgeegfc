#!/usr/bin/env python3
"""Audit round-3 Section-5 redos — cross-probe signs.

Split out of audit_round3_section5_redo.py on 2026-05-29 (Phase 4-B
split 2/7). Imports shared constants + helpers from
_audit_round3_shared and runs: redo5.
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
# Redo 5 — Cross-probe signs, sorted by trace-agreement count
# ===================================================================
def redo5_cross_probe_signs() -> None:
    df = pd.read_csv(R2_TBL / "cross_probe_signs.csv")
    probes = ["d_rank_dS", "kc_lam0", "grassmann_k13", "ctm", "vi_clean"]
    sign_cols = [f"sign_{x}" for x in probes]
    probe_labels = ["d_rank d_S", "KC λ=0", "Grass k=13", "CTM", "VI(k)"]

    # trace direction = sign == -1 ; agreement count per (patient, band)
    df["agree_trace"] = (df[sign_cols] == -1).sum(axis=1)

    # Patient ordering: sort by mean agreement across all 6 bands
    pat_avg = df.groupby("patient")["agree_trace"].mean().sort_values(ascending=False)
    cohort_pat_order = list(pat_avg.index)

    fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.0), sharex=False, sharey=False)
    cmap = plt.cm.RdYlGn_r

    for ax, b in zip(axes.flat, BAND_ORDER):
        sub = df[df["band"] == b].copy()
        # In-panel sort: by agree_trace desc; ties broken by cohort order
        sub["pat_rank"] = sub["patient"].map(
            {p: i for i, p in enumerate(cohort_pat_order)}
        )
        sub = sub.sort_values(["agree_trace", "pat_rank"],
                              ascending=[False, True]).reset_index(drop=True)
        M = sub[sign_cols].values.astype(float)
        ax.imshow(M, aspect="auto", cmap=cmap, vmin=-1, vmax=1,
                  interpolation="nearest")
        ax.set_xticks(range(len(probes)))
        ax.set_xticklabels(probe_labels, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(
            [f"{p}  ({a}/5)" for p, a in zip(sub["patient"], sub["agree_trace"])],
            fontsize=7,
        )
        ax.set_title(BRAIN_BAND_TEX_DICT[b], fontsize=10)
        # Cell glyphs
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                v = M[i, j]
                txt = "−" if v == -1 else ("+" if v == 1 else "·")
                ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                        color="white" if v != 0 else "black")
        # Block dividers: between agree>=4, agree∈{2,3}, agree<=1
        agree = sub["agree_trace"].values
        for thresh in (4, 2):
            # Find first row index where agree < thresh
            below = np.where(agree < thresh)[0]
            if len(below):
                y = float(below[0]) - 0.5
                ax.axhline(y, color="#222", lw=0.9)

    fig.text(
        0.5, 0.005,
        "rows sorted by # probes in trace direction (≥ 4 / 5 = consistently pro · "
        "2–3 = mixed · ≤ 1 = consistently anti) · "
        "green/− = trace direction · red/+ = anti-trace · gray/· = ~zero",
        ha="center", fontsize=8, color="0.4",
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(FIG / "cross_probe_signs.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[redo5] wrote cross_probe_signs.pdf")


def main() -> None:
    print("[redo] 5 — cross_probe_signs")
    redo5_cross_probe_signs()
    print(f"[cross-probe] done · outputs at {OUT}")


if __name__ == "__main__":
    main()
