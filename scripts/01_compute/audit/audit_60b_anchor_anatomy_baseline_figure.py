#!/usr/bin/env python3
"""Audit 60b — anchor anatomy baseline first-pass figure.

Two-row figure consuming
``data/audit/anchor_anatomy_baseline/per_module_same_probe_fraction.csv``:

  Top row    boxplot per band of (anchor / trace / reset / rearrange)
             same-probe-pair-fraction with the per-patient cohort
             baseline overlaid as a horizontal reference at each band.
             6 bands x 4 classes = 24 boxes; the cohort baseline is the
             dotted reference.
  Bottom row per-patient scatter at beta only, x = anchor same-probe
             fraction, y = trace same-probe fraction, with the diagonal
             and the cohort-baseline cross drawn. One marker per
             (patient, anchor module) where a sibling trace module
             exists in the same cell; if multiple, take the size-weighted
             mean per patient.

Output
------
``data/reports/notes_verification_2026-05-08/figures/anchor_anatomy_baseline.pdf``

Convention: PDF only, full vector (CLAUDE.md never-list).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT


BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
CLASSES = ["anchor", "trace", "reset", "rearrange"]
CLASS_COLORS = {
    "anchor": "#807dba",
    "trace": "#2c7fb8",
    "reset": "#41ab5d",
    "rearrange": "#fc9272",
}

CSV_IN = ROOT / "data/audit/anchor_anatomy_baseline/per_module_same_probe_fraction.csv"
PDF_OUT = (
    ROOT / "data/reports/notes_verification_2026-05-08"
    / "figures/anchor_anatomy_baseline.pdf"
)


def _draw_top_panel(ax, df: pd.DataFrame) -> None:
    """Boxplot per (band, class) of same_probe_frac. Cohort baseline drawn
    as a black dot per band (median across patients) plus a thin gray band
    spanning the inter-patient min-max baseline range."""
    n_bands = len(BANDS)
    n_classes = len(CLASSES)
    width = 0.18
    pad = 0.025
    group_width = (width + pad) * n_classes - pad

    band_to_x = {b: i for i, b in enumerate(BANDS)}
    legend_handles = []

    # Boxplots per (band, class)
    for ci, cls in enumerate(CLASSES):
        positions = []
        data = []
        for b in BANDS:
            sub = df[(df["band"] == b) & (df["class"] == cls)]
            x_center = band_to_x[b]
            offset = -group_width / 2 + ci * (width + pad) + width / 2
            positions.append(x_center + offset)
            data.append(sub["same_probe_frac"].dropna().to_numpy())

        bp = ax.boxplot(
            data,
            positions=positions,
            widths=width,
            patch_artist=True,
            showfliers=False,
            medianprops=dict(color="black", linewidth=1.4),
            boxprops=dict(facecolor=CLASS_COLORS[cls], edgecolor="black",
                          linewidth=0.7, alpha=0.85),
            whiskerprops=dict(color="black", linewidth=0.6),
            capprops=dict(color="black", linewidth=0.6),
        )
        legend_handles.append(
            plt.Rectangle(
                (0, 0), 1, 1, fc=CLASS_COLORS[cls], ec="black",
                linewidth=0.7, label=cls,
            )
        )

    # Cohort baseline: per-band, draw a dashed horizontal segment spanning
    # the band's group, at the cohort-MEDIAN baseline. Different patients have
    # slightly different baselines; thin gray band for the inter-patient
    # min-max range.
    for b in BANDS:
        sub = df[df["band"] == b].drop_duplicates("patient")
        if sub.empty:
            continue
        baselines = sub["cohort_baseline_same_probe_pair_frac"].to_numpy()
        med = float(np.median(baselines))
        lo, hi = float(np.min(baselines)), float(np.max(baselines))
        x_center = band_to_x[b]
        x_lo = x_center - group_width / 2 - pad
        x_hi = x_center + group_width / 2 + pad
        ax.fill_between(
            [x_lo, x_hi], [lo, lo], [hi, hi],
            color="0.3", alpha=0.18, zorder=0, linewidth=0,
        )
        ax.plot(
            [x_lo, x_hi], [med, med],
            color="0.15", linestyle=":", linewidth=1.2, zorder=1,
        )
    legend_handles.append(
        plt.Line2D(
            [0], [0], color="0.15", linestyle=":", linewidth=1.2,
            label=r"cohort same-probe baseline (per-patient median $\pm$ range)",
        )
    )

    ax.set_xticks(list(range(n_bands)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=10)
    ax.set_ylim(-0.02, 1.02)
    ax.set_ylabel(r"Within-module same-probe pair fraction", fontsize=10)
    ax.set_title(
        "Per-module same-probe enrichment vs cohort baseline by band and class",
        fontsize=11, pad=6,
    )
    ax.tick_params(axis="y", labelsize=9)
    for sn in ("top", "right"):
        ax.spines[sn].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.grid(axis="y", linewidth=0.4, alpha=0.35, zorder=-1)
    ax.legend(
        handles=legend_handles,
        loc="upper left", bbox_to_anchor=(0.0, 1.0),
        ncol=5, frameon=False, fontsize=9, columnspacing=1.2,
        handletextpad=0.5,
    )


def _draw_bottom_panel(ax, df: pd.DataFrame) -> None:
    """Per-patient scatter at beta only.

    x = patient mean anchor same-probe frac (size-weighted by n_pairs)
    y = patient mean trace same-probe frac (size-weighted by n_pairs)

    Diagonal y = x drawn as a thin reference; cohort same-probe baseline
    drawn as a thin dashed cross (vertical + horizontal).
    """
    df_b = df[(df["band"] == "beta")]
    points = []
    for patient, sub in df_b.groupby("patient"):
        anchor_rows = sub[(sub["class"] == "anchor") & (sub["n_pairs"] > 0)]
        trace_rows = sub[(sub["class"] == "trace") & (sub["n_pairs"] > 0)]
        if anchor_rows.empty or trace_rows.empty:
            continue

        def _wmean(rows):
            w = rows["n_pairs"].astype(float).to_numpy()
            v = rows["same_probe_frac"].astype(float).to_numpy()
            ok = np.isfinite(v) & np.isfinite(w) & (w > 0)
            if not ok.any():
                return float("nan")
            return float(np.average(v[ok], weights=w[ok]))

        x = _wmean(anchor_rows)
        y = _wmean(trace_rows)
        baseline = float(sub["cohort_baseline_same_probe_pair_frac"].iloc[0])
        if not (np.isfinite(x) and np.isfinite(y)):
            continue
        points.append({
            "patient": patient,
            "anchor_same_probe": x,
            "trace_same_probe": y,
            "cohort_baseline": baseline,
        })
    P = pd.DataFrame(points)

    # Thin cohort-baseline cross at median baseline
    if not P.empty:
        med_baseline = float(np.median(P["cohort_baseline"]))
        ax.axvline(
            med_baseline, color="0.5", linestyle=":", linewidth=1.0,
            zorder=0, label=r"cohort baseline (median)",
        )
        ax.axhline(med_baseline, color="0.5", linestyle=":", linewidth=1.0, zorder=0)

    ax.plot(
        [0, 1], [0, 1], color="0.3", linestyle="--", linewidth=0.8, zorder=1,
        label=r"y = x",
    )

    pat_color = plt.get_cmap("tab10")
    P = P.reset_index(drop=True).sort_values("patient").reset_index(drop=True)
    # offset close labels via a per-point bias (small jitter only on annotation, not on marker)
    seen_xy = []
    for i, row in P.iterrows():
        ax.scatter(
            row["anchor_same_probe"], row["trace_same_probe"],
            s=72, color=pat_color(i % 10), edgecolors="black", linewidths=0.6,
            zorder=3, label=row["patient"],
        )
        # label offset: +5,+4 by default; if a previous point sits within 0.04
        # of this one, shift the offset by sign of the index to avoid stacking.
        dx, dy = 6, 5
        for px, py in seen_xy:
            if abs(px - row["anchor_same_probe"]) < 0.04 and abs(py - row["trace_same_probe"]) < 0.04:
                dy = -10
                break
        ax.annotate(
            row["patient"].replace("Pat_", ""),
            (row["anchor_same_probe"], row["trace_same_probe"]),
            xytext=(dx, dy), textcoords="offset points",
            fontsize=8, color="0.15", zorder=4,
        )
        seen_xy.append((row["anchor_same_probe"], row["trace_same_probe"]))

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    ax.set_xlabel(
        r"$\beta$ anchor same-probe pair fraction (size-weighted mean per patient)",
        fontsize=10,
    )
    ax.set_ylabel(
        r"$\beta$ trace same-probe pair fraction (size-weighted mean per patient)",
        fontsize=10,
    )
    ax.set_title(
        r"Per-patient $\beta$ anchor vs trace same-probe enrichment",
        fontsize=11, pad=6,
    )
    ax.tick_params(axis="both", labelsize=9)
    for sn in ("top", "right"):
        ax.spines[sn].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.grid(linewidth=0.4, alpha=0.3, zorder=-2)
    ax.legend(
        loc="upper left", bbox_to_anchor=(1.02, 1.0),
        frameon=False, fontsize=8.5,
    )


def main():
    df = pd.read_csv(CSV_IN)

    fig = plt.figure(figsize=(12.0, 11.5))
    gs = fig.add_gridspec(
        2, 1, height_ratios=[1.0, 1.05], hspace=0.32,
        left=0.08, right=0.85, top=0.96, bottom=0.06,
    )
    ax_top = fig.add_subplot(gs[0])
    ax_bot = fig.add_subplot(gs[1])

    _draw_top_panel(ax_top, df)
    _draw_bottom_panel(ax_bot, df)

    PDF_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF_OUT)
    plt.close(fig)
    print(f"[audit_60b] wrote {PDF_OUT}")


if __name__ == "__main__":
    main()
