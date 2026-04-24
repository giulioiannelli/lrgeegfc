#!/usr/bin/env python3
"""Compute unanimity maps and scalar VI contrasts for ImCoh LRG results.

Reads VI profiles from data/imcoh_vi/vi_raw_profiles.csv, then:
  1. For each hypothesis (H1, H2a, H2b, H3) at each (band, k):
     - Compute contrast for each patient
     - Check if all 5 patients agree on sign (unanimity)
  2. Compute scalar summaries (mean contrast, integrated, fraction positive)

Outputs (under data/imcoh_unanimity/):
  unanimity_maps.csv           -- (band, k, hypothesis, unanimity, ...)
  scalar_contrasts.csv         -- per-patient scalar summaries
  unanimity_summary.md         -- formatted report
  unanimity_map_{hyp}.pdf      -- heatmap figures

Run: python scripts/py/compute_imcoh_unanimity.py [-v]
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    ALL_PHASE_PAIRS,
    BRAIN_BAND_TEX_DICT,
    BRAIN_BANDS_NAMES,
    PATIENTS_4PHASE,
    PHASE_LABELS,
)
from lrg_eegfc.config.paths import REPORTS_ROOT

PATIENTS = PATIENTS_4PHASE
BAND_ORDER = BRAIN_BANDS_NAMES
BAND_TEX = [BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER]

PHASES = list(PHASE_LABELS)
ALL_PAIRS = list(ALL_PHASE_PAIRS)
WITHIN_PAIRS = [("rest_pre", "rest_post"), ("task_learn", "task_test")]
CROSS_PAIRS = [
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"),
    ("task_learn", "rest_post"), ("task_test", "rest_post"),
]

VI_DIR = REPORTS_ROOT / "imcoh_vi"
OUT = REPORTS_ROOT / "imcoh_unanimity"


# ── VI helpers ──────────────────────────────────────────────────────────
def _get_vi(pdf: pd.DataFrame, pa: str, pb: str) -> float:
    row = pdf[
        ((pdf["phase_a"] == pa) & (pdf["phase_b"] == pb)) |
        ((pdf["phase_a"] == pb) & (pdf["phase_b"] == pa))
    ]
    if row.empty:
        return np.nan
    return row["vi"].iloc[0]


# ── Hypothesis contrasts ───────────────────────────────────────────────
def _contrast_H1(pdf: pd.DataFrame) -> float:
    """H1: task stability. POSITIVE = supported.
    mean(other pairs VI) - VI(taskL, taskT)."""
    vi_tt = _get_vi(pdf, "task_learn", "task_test")
    other_vis = []
    for pa, pb in ALL_PAIRS:
        if (pa, pb) == ("task_learn", "task_test"):
            continue
        v = _get_vi(pdf, pa, pb)
        if not np.isnan(v):
            other_vis.append(v)
    if np.isnan(vi_tt) or len(other_vis) == 0:
        return np.nan
    return np.mean(other_vis) - vi_tt


def _contrast_H2a(pdf: pd.DataFrame) -> float:
    """H2a: task trace. POSITIVE = trace detected.
    VI(rest_pre,rest_post) - VI(taskT,rest_post)."""
    return _get_vi(pdf, "rest_pre", "rest_post") - _get_vi(pdf, "task_test", "rest_post")


def _contrast_H2b(pdf: pd.DataFrame) -> float:
    """H2b: task approach. POSITIVE = approach detected.
    VI(rest_pre,taskT) - VI(taskT,rest_post)."""
    return _get_vi(pdf, "rest_pre", "task_test") - _get_vi(pdf, "task_test", "rest_post")


def _contrast_H3(pdf: pd.DataFrame) -> float:
    """H3: within < cross. POSITIVE = supported.
    mean(cross VI) - mean(within VI)."""
    within_vis = []
    for pa, pb in WITHIN_PAIRS:
        v = _get_vi(pdf, pa, pb)
        if not np.isnan(v):
            within_vis.append(v)
    cross_vis = []
    for pa, pb in CROSS_PAIRS:
        v = _get_vi(pdf, pa, pb)
        if not np.isnan(v):
            cross_vis.append(v)
    if not within_vis or not cross_vis:
        return np.nan
    return np.mean(cross_vis) - np.mean(within_vis)


HYPOTHESES = [
    ("H1", "Task Stability", _contrast_H1),
    ("H2a", "Task Trace", _contrast_H2a),
    ("H2b", "Task Approach", _contrast_H2b),
    ("H3", "Within < Cross", _contrast_H3),
]


# ── Step 1: Unanimity maps ────────────────────────────────────────────
def compute_unanimity_maps(vi_df: pd.DataFrame) -> pd.DataFrame:
    """Compute signed unanimity at every (hypothesis, band, k)."""
    rows = []
    for hyp_name, hyp_label, contrast_fn in HYPOTHESES:
        for band in BAND_ORDER:
            bdf = vi_df[vi_df["band"] == band]
            k_values = sorted(bdf["k"].unique())

            for k in k_values:
                kdf = bdf[bdf["k"] == k]
                signs = []

                for pat in PATIENTS:
                    pdf = kdf[kdf["patient"] == pat]
                    if pdf.empty:
                        continue
                    contrast = contrast_fn(pdf)
                    if contrast is not None and not np.isnan(contrast):
                        signs.append(np.sign(contrast))

                n_patients = len(signs)
                if n_patients == 0:
                    continue

                unanimity = np.mean(signs)
                n_positive = sum(1 for s in signs if s > 0)
                n_negative = sum(1 for s in signs if s < 0)

                rows.append({
                    "hypothesis": hyp_name,
                    "band": band,
                    "k": k,
                    "unanimity": unanimity,
                    "n_positive": n_positive,
                    "n_negative": n_negative,
                    "n_patients": n_patients,
                    "is_unanimous": abs(unanimity) == 1.0 and n_patients == len(PATIENTS),
                })

    return pd.DataFrame(rows)


# ── Step 2: Scalar summaries ──────────────────────────────────────────
def compute_scalar_contrasts(vi_df: pd.DataFrame) -> pd.DataFrame:
    """For each (patient, band, hypothesis), compute scalar summaries
    over all k values: mean contrast, integrated, fraction positive."""
    rows = []

    for hyp_name, hyp_label, contrast_fn in HYPOTHESES:
        for pat in PATIENTS:
            for band in BAND_ORDER:
                bdf = vi_df[(vi_df["band"] == band) & (vi_df["patient"] == pat)]
                k_values = sorted(bdf["k"].unique())
                if not k_values:
                    continue

                contrasts = []
                for k in k_values:
                    kdf = bdf[bdf["k"] == k]
                    c = contrast_fn(kdf)
                    if c is not None and not np.isnan(c):
                        contrasts.append(c)

                if not contrasts:
                    continue

                contrasts = np.array(contrasts)
                rows.append({
                    "hypothesis": hyp_name,
                    "patient": pat,
                    "band": band,
                    "scalar_mean": float(np.mean(contrasts)),
                    "scalar_integrated": float(np.sum(contrasts)),
                    "scalar_fraction_pos": float(np.mean(contrasts > 0)),
                    "n_k": len(contrasts),
                })

    return pd.DataFrame(rows)


# ── Step 3: Figures ────────────────────────────────────────────────────
def plot_unanimity_maps(unan_df: pd.DataFrame, out_dir: Path):
    """Generate one heatmap per hypothesis: band x k colored by unanimity."""
    import matplotlib.pyplot as plt

    for hyp_name, hyp_label, _ in HYPOTHESES:
        hdf = unan_df[unan_df["hypothesis"] == hyp_name]
        if hdf.empty:
            print(f"  {hyp_name}: no data, skipping figure")
            continue

        fig, ax = plt.subplots(figsize=(16, 4.5))

        piv = hdf.pivot_table(index="band", columns="k", values="unanimity")
        piv = piv.reindex(BAND_ORDER)

        im = ax.imshow(
            piv.values, aspect="auto", cmap="RdBu", vmin=-1, vmax=1,
            interpolation="nearest",
        )

        # Mark unanimous cells with black border
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                v = piv.values[i, j]
                if not np.isnan(v) and abs(v) == 1.0:
                    band = BAND_ORDER[i]
                    k = piv.columns[j]
                    row = hdf[(hdf["band"] == band) & (hdf["k"] == k)]
                    if not row.empty and row.iloc[0]["n_patients"] == len(PATIENTS):
                        rect = plt.Rectangle(
                            (j - 0.5, i - 0.5), 1, 1,
                            fill=False, edgecolor="black", lw=1.5,
                        )
                        ax.add_patch(rect)

        ax.set_yticks(range(len(BAND_ORDER)))
        ax.set_yticklabels(BAND_TEX, fontsize=10)

        k_vals = list(piv.columns)
        step = max(1, len(k_vals) // 20)
        tick_idx = list(range(0, len(k_vals), step))
        ax.set_xticks(tick_idx)
        ax.set_xticklabels([k_vals[t] for t in tick_idx], fontsize=7)
        ax.set_xlabel("k (number of clusters)", fontsize=10)

        fig.colorbar(im, ax=ax, label="Signed unanimity (5 patients)", shrink=0.8)
        ax.set_title(f"{hyp_name}: {hyp_label} -- ImCoh unanimity map", fontsize=12)
        fig.tight_layout()

        fig_path = out_dir / f"unanimity_map_{hyp_name}.pdf"
        fig.savefig(fig_path, bbox_inches="tight", dpi=300)
        plt.close(fig)
        print(f"  Saved: {fig_path}")


# ── Step 4: Summary report ────────────────────────────────────────────
def write_summary(unan_df: pd.DataFrame, scalar_df: pd.DataFrame,
                  out_dir: Path):
    """Write formatted markdown summary."""
    lines = [
        "# ImCoh -- Multiscale Hypothesis Testing Summary",
        "",
        f"**N = {len(PATIENTS)} patients**: {', '.join(PATIENTS)}",
        "**FC method**: Imaginary Coherence (ImCoh)",
        "**Criterion**: Strict unanimity -- 5/5 patients must agree on the sign",
        "",
        "With N=5, 5/5 agreement gives Wilcoxon one-sided p = 0.031.",
        "",
        "---",
        "",
    ]

    for hyp_name, hyp_label, _ in HYPOTHESES:
        hdf = unan_df[unan_df["hypothesis"] == hyp_name]
        lines.append(f"## {hyp_name}: {hyp_label}")
        lines.append("")

        unan_cells = hdf[hdf["is_unanimous"] == True]
        pos_cells = unan_cells[unan_cells["unanimity"] > 0]
        neg_cells = unan_cells[unan_cells["unanimity"] < 0]

        if pos_cells.empty and neg_cells.empty:
            lines.append("**No unanimous (band, k) cells.** Not supported at any scale.")
            lines.append("")
        else:
            if not pos_cells.empty:
                lines.append(f"### Supported (+1): {len(pos_cells)} cells")
                lines.append("")
                lines.append("| Band | k range | # cells |")
                lines.append("|------|---------|---------|")
                for band in BAND_ORDER:
                    bc = pos_cells[pos_cells["band"] == band]
                    if bc.empty:
                        continue
                    k_min, k_max = bc["k"].min(), bc["k"].max()
                    lines.append(f"| {band} | {k_min}--{k_max} | {len(bc)} |")
                lines.append("")

            if not neg_cells.empty:
                lines.append(f"### Reversed (-1): {len(neg_cells)} cells")
                lines.append("")
                lines.append("| Band | k range | # cells |")
                lines.append("|------|---------|---------|")
                for band in BAND_ORDER:
                    bc = neg_cells[neg_cells["band"] == band]
                    if bc.empty:
                        continue
                    k_min, k_max = bc["k"].min(), bc["k"].max()
                    lines.append(f"| {band} | {k_min}--{k_max} | {len(bc)} |")
                lines.append("")

        # Scalar summary
        sdf = scalar_df[scalar_df["hypothesis"] == hyp_name]
        if not sdf.empty:
            lines.append("### Scalar mean contrast by band")
            lines.append("")
            lines.append("| Band | " + " | ".join(PATIENTS) + " |")
            lines.append("|------|" + "|".join(["------"] * len(PATIENTS)) + "|")
            for band in BAND_ORDER:
                vals = []
                for pat in PATIENTS:
                    sub = sdf[(sdf["patient"] == pat) & (sdf["band"] == band)]
                    if not sub.empty:
                        vals.append(f"{sub['scalar_mean'].iloc[0]:+.4f}")
                    else:
                        vals.append("--")
                lines.append(f"| {band} | " + " | ".join(vals) + " |")
            lines.append("")

        lines.append("---")
        lines.append("")

    (out_dir / "unanimity_summary.md").write_text("\n".join(lines))
    print(f"  Summary: {out_dir / 'unanimity_summary.md'}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute unanimity maps and scalar contrasts for ImCoh."
    )
    parser.add_argument("--no-figures", action="store_true",
                        help="Skip figure generation")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # Load VI profiles
    vi_path = VI_DIR / "vi_raw_profiles.csv"
    if not vi_path.exists():
        print(f"ERROR: {vi_path} not found.")
        print("Run compute_imcoh_vi.py first.")
        sys.exit(1)

    vi_df = pd.read_csv(vi_path)
    print(f"Loaded {len(vi_df)} VI profile points from {vi_path}")

    # Step 1: Unanimity maps
    print("\nComputing unanimity maps...")
    unan_df = compute_unanimity_maps(vi_df)
    unan_path = OUT / "unanimity_maps.csv"
    unan_df.to_csv(unan_path, index=False)
    print(f"  Saved: {unan_path}")

    for hyp_name, _, _ in HYPOTHESES:
        hdf = unan_df[unan_df["hypothesis"] == hyp_name]
        n_unan = hdf["is_unanimous"].sum()
        n_pos = (hdf[hdf["is_unanimous"]]["unanimity"] > 0).sum()
        n_neg = (hdf[hdf["is_unanimous"]]["unanimity"] < 0).sum()
        print(f"  {hyp_name}: {n_unan} unanimous cells "
              f"({n_pos} positive, {n_neg} negative)")

    # Step 2: Scalar contrasts
    print("\nComputing scalar contrasts...")
    scalar_df = compute_scalar_contrasts(vi_df)
    scalar_path = OUT / "scalar_contrasts.csv"
    scalar_df.to_csv(scalar_path, index=False)
    print(f"  Saved: {scalar_path}")

    # Step 3: Figures
    if not args.no_figures:
        print("\nGenerating figures...")
        plot_unanimity_maps(unan_df, OUT)

    # Step 4: Summary
    print("\nWriting summary...")
    write_summary(unan_df, scalar_df, OUT)

    elapsed = time.time() - t0
    print(f"\nDone ({elapsed:.1f}s). Outputs in {OUT}")


if __name__ == "__main__":
    main()
