#!/usr/bin/env python3
"""Forward vs backward trace comparison.

Forward  = J*(post_subtree, T_task_test)  − J*(post_subtree, T_rest_pre)
Backward = J*(pre_subtree,  T_task_learn) − J*(pre_subtree,  T_rest_post)

If the effect is task-induced (rest_post inherits task structure that
rest_pre is untouched by), forward >> backward. If the effect is
session-order or anatomy (each rest phase is more similar to its
temporally-adjacent task phase, regardless of task content), forward
≈ backward.

Reads:  ``trace_sweep.csv``, ``trace_null.csv``,
        ``trace_backward_sweep.csv``, ``trace_backward_null.csv``
Writes: ``data/reports/imcoh_mrl/figures/trace_forward_vs_backward.pdf``
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT


IN_DIR = REPORTS_ROOT / "imcoh_mrl"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT

C_FW = "#d62728"
C_BW = "#1f77b4"


def main() -> None:
    fw = pd.read_csv(IN_DIR / "trace_sweep.csv")
    fw_null = pd.read_csv(IN_DIR / "trace_null.csv")
    bw = pd.read_csv(IN_DIR / "trace_backward_sweep.csv")
    bw_null = pd.read_csv(IN_DIR / "trace_backward_null.csv")

    fw_full = fw.merge(fw_null[["patient", "band", "p_empirical", "z_observed"]],
                       on=["patient", "band"], how="left")
    bw_full = bw.merge(bw_null[["patient", "band", "p_empirical", "z_observed"]],
                       on=["patient", "band"], how="left")

    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.0), dpi=150,
                             sharey=True)
    rng = np.random.default_rng(0)

    for ax, df, color, title in [
        (axes[0], fw_full, C_FW,
         r"Forward — $J^*(\mathrm{post}{\to}\mathrm{test}) - J^*(\mathrm{post}{\to}\mathrm{pre})$"),
        (axes[1], bw_full, C_BW,
         r"Backward — $J^*(\mathrm{pre}{\to}\mathrm{learn}) - J^*(\mathrm{pre}{\to}\mathrm{post})$"),
    ]:
        for j, band in enumerate(BANDS):
            sub = df[df["band"] == band]
            if sub.empty:
                continue
            scores = sub["score"].to_numpy()
            xs = j + rng.uniform(-0.18, 0.18, size=len(scores))
            ax.scatter(xs, scores, s=24, color=color, alpha=0.55,
                       edgecolor="black", linewidth=0.3, zorder=3)
            med = float(np.median(scores))
            ax.hlines(med, j - 0.22, j + 0.22, color="#404040",
                      lw=2.2, zorder=5)
            n_01 = int((sub["p_empirical"] < 0.01).sum())
            ax.text(j, 1.04, f"{n_01}/{len(sub)}",
                    ha="center", va="bottom", fontsize=10, color=color,
                    transform=ax.get_xaxis_transform(), fontweight="bold")

        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([TEX[b] for b in BANDS], fontsize=12)
        ax.set_xlim(-0.5, len(BANDS) - 0.5)
        ax.axhline(0, color="black", lw=0.6, ls="-", zorder=1)
        ax.set_ylim(-0.05, 1.08)
        ax.set_xlabel("band", fontsize=11)
        ax.set_title(title + "  (numbers = patients with p<0.01)",
                     fontsize=10.5)

    axes[0].set_ylabel("score", fontsize=11)

    fig.tight_layout()
    out = OUT_DIR / "trace_forward_vs_backward.pdf"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"saved {out}")
    plt.close(fig)

    # Quick numerical compare
    print()
    print("[fwd vs bwd] per-band median score and #patients p<0.01:")
    print(f"{'band':<10s}  {'fw_med':>7s}  {'bw_med':>7s}  "
          f"{'Δ(fw-bw)':>9s}  {'fw n01':>6s}  {'bw n01':>6s}")
    for band in BANDS:
        f = fw_full[fw_full["band"] == band]
        b = bw_full[bw_full["band"] == band]
        if f.empty or b.empty:
            continue
        fm = f["score"].median()
        bm = b["score"].median()
        f01 = int((f["p_empirical"] < 0.01).sum())
        b01 = int((b["p_empirical"] < 0.01).sum())
        print(f"{band:<10s}  {fm:>+.3f}  {bm:>+.3f}  {fm-bm:>+.3f}  "
              f"{f01:>3d}/{len(f):<3d}  {b01:>3d}/{len(b):<3d}")


if __name__ == "__main__":
    main()
