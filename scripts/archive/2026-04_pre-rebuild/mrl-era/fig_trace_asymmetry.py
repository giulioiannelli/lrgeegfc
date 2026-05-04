#!/usr/bin/env python3
"""Forward−backward asymmetry per band: the genuinely task-specific signal.

For each (patient, band):
    Δ(p, b) := score_forward(p, b) − score_backward(p, b)
where
    forward  = J*(post_subtree, T_test) − J*(post_subtree, T_pre)
    backward = J*(pre_subtree,  T_learn) − J*(pre_subtree,  T_post).

If task induces a genuinely asymmetric trace, Δ > 0 cohort-wide; if
the effect is symmetric session-order/anatomy, Δ ≈ 0.

Per-band test: one-sided Wilcoxon signed-rank on per-patient Δ
(H₁: Δ > 0), BH-FDR over 6 bands, rank-biserial effect size,
95% bootstrap CI.

Reads:  ``trace_sweep.csv``, ``trace_backward_sweep.csv``
Writes: ``data/reports/imcoh_mrl/figures/trace_asymmetry.pdf``,
        ``data/reports/imcoh_mrl/trace_asymmetry.csv``
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT
from lrg_eegfc.utils.metrics.hypothesis import (
    wilcoxon_z, rank_biserial, boot_ci_mean, bh_fdr,
)


IN_DIR = REPORTS_ROOT / "imcoh_mrl"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT

C_POS = "#2ca02c"   # Δ > 0 (task-asymmetric)
C_NEG = "#d62728"   # Δ < 0
C_NULL = "#cccccc"


def main() -> None:
    fw = pd.read_csv(IN_DIR / "trace_sweep.csv")[["patient", "band", "score"]] \
            .rename(columns={"score": "fw"})
    bw = pd.read_csv(IN_DIR / "trace_backward_sweep.csv")[
            ["patient", "band", "score"]].rename(columns={"score": "bw"})
    df = fw.merge(bw, on=["patient", "band"], how="inner")
    df["delta"] = df["fw"] - df["bw"]

    # Per-band Wilcoxon (one-sided greater).
    band_stats = []
    for band in BANDS:
        vals = df[df["band"] == band]["delta"].to_numpy()
        z, p = wilcoxon_z(vals)
        r_rb = rank_biserial(vals)
        m, lo, hi = boot_ci_mean(vals)
        band_stats.append({
            "band": band, "n": len(vals),
            "median": float(np.median(vals)),
            "mean": m, "ci_lo": lo, "ci_hi": hi,
            "z": z, "p_raw": p, "r_rb": r_rb,
            "n_pos": int((vals > 0).sum()),
        })
    p_raw = [r["p_raw"] for r in band_stats]
    p_bh  = bh_fdr(p_raw)
    for r, q in zip(band_stats, p_bh):
        r["q_bh"] = q

    out_csv = pd.DataFrame(band_stats)
    out_csv.to_csv(IN_DIR / "trace_asymmetry.csv", index=False)
    print("[asym] per-band paired Wilcoxon (H₁: Δ > 0):")
    print(f"{'band':<10s}  {'med':>7s}  {'mean':>7s}  "
          f"{'95% CI':>17s}  {'r_rb':>5s}  {'p_raw':>6s}  {'q_BH':>6s}  "
          f"{'pos':>5s}")
    for r in band_stats:
        ci = f"[{r['ci_lo']:+.3f},{r['ci_hi']:+.3f}]"
        star = ""
        if r["q_bh"] < 0.001: star = "★★★"
        elif r["q_bh"] < 0.01: star = "★★"
        elif r["q_bh"] < 0.05: star = "★"
        print(f"{r['band']:<10s}  {r['median']:>+6.3f}  {r['mean']:>+6.3f}  "
              f"{ci:>15s}  {r['r_rb']:>+.2f}  {r['p_raw']:>6.3f}  "
              f"{r['q_bh']:>6.3f} {star}  "
              f"{r['n_pos']:>2d}/{r['n']}")

    # ───── Figure ─────
    fig, ax = plt.subplots(figsize=(10.0, 5.4), dpi=150)
    rng = np.random.default_rng(0)

    for j, band in enumerate(BANDS):
        vals = df[df["band"] == band]["delta"].to_numpy()
        if len(vals) == 0:
            continue
        cols = [C_POS if v > 0 else C_NEG for v in vals]
        xs = j + rng.uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(xs, vals, s=44, color=cols, alpha=0.85,
                   edgecolor="black", linewidth=0.3, zorder=4)

        med = float(np.median(vals))
        q1, q3 = float(np.quantile(vals, 0.25)), float(np.quantile(vals, 0.75))
        ax.vlines(j, q1, q3, color="#404040", lw=2.4, zorder=5, alpha=0.85)
        ax.hlines(med, j - 0.22, j + 0.22, color="#404040", lw=2.6, zorder=6)

        bs = next(r for r in band_stats if r["band"] == band)
        star = ("★★★" if bs["q_bh"] < 0.001
                else "★★" if bs["q_bh"] < 0.01
                else "★"  if bs["q_bh"] < 0.05
                else "")
        n_pos = bs["n_pos"]
        ax.text(j, 1.04, f"{n_pos}/{bs['n']}", ha="center", va="bottom",
                fontsize=10, color="#404040", fontweight="bold",
                transform=ax.get_xaxis_transform())
        ax.text(j, 1.10,
                f"q = {bs['q_bh']:.3f} {star}" if bs["q_bh"] < 0.5
                else f"q = {bs['q_bh']:.2f}",
                ha="center", va="bottom", fontsize=9,
                color=(C_POS if bs["q_bh"] < 0.05 else "#888888"),
                transform=ax.get_xaxis_transform())

    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([TEX[b] for b in BANDS], fontsize=12)
    ax.set_xlim(-0.5, len(BANDS) - 0.5)
    ax.axhline(0, color="black", lw=0.8, ls="-", zorder=1)
    ymin = float(df["delta"].min()) - 0.05
    ymax = float(df["delta"].max()) + 0.05
    ax.set_ylim(ymin, max(ymax, 0.25))
    ax.set_ylabel(r"$\Delta = $ score$_{\mathrm{forward}}$ "
                  r"$- $ score$_{\mathrm{backward}}$", fontsize=11)
    ax.set_xlabel("band", fontsize=11)

    handles = [
        plt.Line2D([0], [0], marker='o', color=C_POS, lw=0,
                   markersize=7, markeredgecolor="black",
                   markeredgewidth=0.3, label="Δ > 0 (task-asymmetric)"),
        plt.Line2D([0], [0], marker='o', color=C_NEG, lw=0,
                   markersize=7, markeredgecolor="black",
                   markeredgewidth=0.3, label="Δ < 0"),
        plt.Line2D([0], [0], color="#404040", lw=2.4,
                   label="cohort median + IQR"),
    ]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=9)

    ax.set_title("Forward − backward asymmetry per band — "
                 "paired Wilcoxon (H₁: Δ > 0), BH-FDR over 6 bands.  "
                 "Numbers above: #patients with Δ > 0",
                 fontsize=10)

    fig.tight_layout()
    out = OUT_DIR / "trace_asymmetry.pdf"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print()
    print(f"saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
