#!/usr/bin/env python3
"""Scale-specific analysis of H2 (task trace in rest_post).

Mean VI averaged across all k gives mixed H2 results. This script tests
whether the task trace is visible only at specific hierarchical scales.

H2 contrasts (VI is a distance -- lower = more similar):
  H2a: VI(Pre,Post) - VI(TT,Post)
        Positive => rest_post closer to task_test than to rest_pre (task persists)
  H2b: VI(Pre,TT)  - VI(TT,Post)
        Positive => rest_post closer to task_test than rest_pre was (task persists)

We examine each k individually and also macro/meso/micro scale ranges.
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

from lrg_eegfc.config.paths import FIGURES_ROOT

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
DATA_CSV = FIGURES_ROOT / "metric_exploration" / "partition_multiscale" / "results.csv"
OUT_DIR = FIGURES_ROOT / "metric_exploration" / "scale_specific_h2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

from lrg_eegfc.config.const import PATIENTS_4PHASE
PATIENTS = PATIENTS_4PHASE
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_LABELS = {
    "delta": r"$\delta$",
    "theta": r"$\theta$",
    "alpha": r"$\alpha$",
    "beta": r"$\beta$",
    "low_gamma": r"low $\gamma$",
    "high_gamma": r"high $\gamma$",
}
K_VALUES = [2, 3, 4, 5, 6, 8, 10, 15, 20]
SCALE_RANGES = {
    "macro\n(k=2-4)": [2, 3, 4],
    "meso\n(k=5-8)": [5, 6, 8],
    "micro\n(k=10-20)": [10, 15, 20],
}

N_PATIENTS = len(PATIENTS)


def load_data() -> pd.DataFrame:
    """Load results CSV and filter to 4-phase patients."""
    df = pd.read_csv(DATA_CSV)
    df = df[df["patient"].isin(PATIENTS)].copy()
    return df


def get_vi(df: pd.DataFrame, patient: str, band: str, pair: str, k: int) -> float:
    """Extract VI at a given k for a specific (patient, band, pair)."""
    row = df[(df["patient"] == patient) & (df["band"] == band) & (df["pair"] == pair)]
    if row.empty:
        return np.nan
    return row.iloc[0][f"VI_k{k}"]


def compute_contrasts_single_k(df: pd.DataFrame):
    """Compute H2a and H2b for each (patient, band, k)."""
    records = []
    for pat in PATIENTS:
        for band in BANDS:
            for k in K_VALUES:
                vi_pre_post = get_vi(df, pat, band, "rest_pre-rest_post", k)
                vi_tt_post = get_vi(df, pat, band, "task_test-rest_post", k)
                vi_pre_tt = get_vi(df, pat, band, "rest_pre-task_test", k)

                h2a = vi_pre_post - vi_tt_post  # positive = task persists
                h2b = vi_pre_tt - vi_tt_post     # positive = task persists

                records.append({
                    "patient": pat,
                    "band": band,
                    "k": k,
                    "VI_PrePost": vi_pre_post,
                    "VI_TTPost": vi_tt_post,
                    "VI_PreTT": vi_pre_tt,
                    "H2a": h2a,
                    "H2b": h2b,
                })
    return pd.DataFrame(records)


def compute_contrasts_scale_ranges(df: pd.DataFrame):
    """Average VI within macro/meso/micro ranges, then compute contrasts."""
    records = []
    for pat in PATIENTS:
        for band in BANDS:
            for range_name, ks in SCALE_RANGES.items():
                vi_pre_post = np.mean([get_vi(df, pat, band, "rest_pre-rest_post", k) for k in ks])
                vi_tt_post = np.mean([get_vi(df, pat, band, "task_test-rest_post", k) for k in ks])
                vi_pre_tt = np.mean([get_vi(df, pat, band, "rest_pre-task_test", k) for k in ks])

                h2a = vi_pre_post - vi_tt_post
                h2b = vi_pre_tt - vi_tt_post

                records.append({
                    "patient": pat,
                    "band": band,
                    "scale_range": range_name,
                    "VI_PrePost_avg": vi_pre_post,
                    "VI_TTPost_avg": vi_tt_post,
                    "VI_PreTT_avg": vi_pre_tt,
                    "H2a": h2a,
                    "H2b": h2b,
                })
    return pd.DataFrame(records)


def compute_unanimity(contrasts_df: pd.DataFrame, group_col: str, contrast_col: str):
    """For each (band, group), count how many patients have positive contrast.

    Returns a DataFrame with unanimity score: +N (all persist) to -N (all recover).
    Score = 2 * n_positive - N_PATIENTS  (range: -N to +N)
    """
    rows = []
    groups = sorted(contrasts_df[group_col].unique(),
                    key=lambda x: K_VALUES.index(x) if x in K_VALUES else 0)
    for band in BANDS:
        for grp in groups:
            sub = contrasts_df[(contrasts_df["band"] == band) & (contrasts_df[group_col] == grp)]
            n_pos = (sub[contrast_col] > 0).sum()
            n_neg = (sub[contrast_col] < 0).sum()
            n_zero = (sub[contrast_col] == 0).sum()
            unanimity = n_pos - n_neg  # range: -4 to +4
            mean_contrast = sub[contrast_col].mean()
            rows.append({
                "band": band,
                group_col: grp,
                "n_positive": int(n_pos),
                "n_negative": int(n_neg),
                "unanimity": unanimity,
                "mean_contrast": mean_contrast,
            })
    return pd.DataFrame(rows)


def plot_heatmap(unan_df: pd.DataFrame, group_col: str, contrast_name: str,
                 title: str, filename: str, group_labels=None):
    """Plot a heatmap: rows=bands, columns=group (k or scale range)."""
    groups = sorted(unan_df[group_col].unique(),
                    key=lambda x: K_VALUES.index(x) if x in K_VALUES else
                    list(SCALE_RANGES.keys()).index(x))

    # Build matrix
    mat = np.full((len(BANDS), len(groups)), np.nan)
    for i, band in enumerate(BANDS):
        for j, grp in enumerate(groups):
            row = unan_df[(unan_df["band"] == band) & (unan_df[group_col] == grp)]
            if not row.empty:
                mat[i, j] = row.iloc[0]["unanimity"]

    fig, ax = plt.subplots(figsize=(max(6, len(groups) * 0.8 + 2), 4.5))

    # Diverging colormap centered at 0
    vmax = N_PATIENTS
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = ax.imshow(mat, cmap="RdBu_r", norm=norm, aspect="auto")

    # Labels
    if group_labels is None:
        group_labels = [str(g) for g in groups]
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(group_labels, fontsize=10)
    ax.set_yticks(range(len(BANDS)))
    ax.set_yticklabels([BAND_LABELS[b] for b in BANDS], fontsize=11)
    ax.set_xlabel("Number of clusters (k)", fontsize=11)

    # Annotate cells with unanimity value and mean contrast
    for i in range(len(BANDS)):
        for j in range(len(groups)):
            val = mat[i, j]
            if not np.isnan(val):
                band = BANDS[i]
                grp = groups[j]
                row = unan_df[(unan_df["band"] == band) & (unan_df[group_col] == grp)]
                mean_c = row.iloc[0]["mean_contrast"]
                sign_str = f"{int(val):+d}"
                # Color text for readability
                text_color = "white" if abs(val) >= 3 else "black"
                ax.text(j, i, f"{sign_str}\n({mean_c:.3f})",
                        ha="center", va="center", fontsize=8,
                        color=text_color, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, label="Unanimity score", shrink=0.85)
    cbar.set_ticks(range(-N_PATIENTS, N_PATIENTS + 1))

    ax.set_title(title, fontsize=12, fontweight="bold", pad=10)
    fig.tight_layout()
    fig.savefig(OUT_DIR / filename, bbox_inches="tight")
    print(f"  Saved: {OUT_DIR / filename}")
    plt.close(fig)


def plot_patient_detail(contrasts_df: pd.DataFrame, contrast_col: str,
                        title: str, filename: str):
    """Heatmap showing each patient's contrast at each k, paneled by band."""
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey=True)
    axes = axes.flatten()

    for idx, band in enumerate(BANDS):
        ax = axes[idx]
        sub = contrasts_df[contrasts_df["band"] == band]

        # Build matrix: rows=patients, cols=k
        mat = np.full((len(PATIENTS), len(K_VALUES)), np.nan)
        for i, pat in enumerate(PATIENTS):
            for j, k in enumerate(K_VALUES):
                row = sub[(sub["patient"] == pat) & (sub["k"] == k)]
                if not row.empty:
                    mat[i, j] = row.iloc[0][contrast_col]

        vmax = max(0.3, np.nanmax(np.abs(mat)))
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
        im = ax.imshow(mat, cmap="RdBu_r", norm=norm, aspect="auto")

        ax.set_title(BAND_LABELS[band], fontsize=12)
        ax.set_xticks(range(len(K_VALUES)))
        ax.set_xticklabels([str(k) for k in K_VALUES], fontsize=8)
        ax.set_yticks(range(len(PATIENTS)))
        ax.set_yticklabels([p.replace("Pat_0", "P") for p in PATIENTS], fontsize=9)

        # Annotate
        for i in range(len(PATIENTS)):
            for j in range(len(K_VALUES)):
                val = mat[i, j]
                if not np.isnan(val):
                    text_color = "white" if abs(val) > vmax * 0.6 else "black"
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                            fontsize=6.5, color=text_color)

        fig.colorbar(im, ax=ax, shrink=0.7)

    # Common labels
    fig.text(0.5, 0.02, "Number of clusters (k)", ha="center", fontsize=11)
    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    fig.savefig(OUT_DIR / filename, bbox_inches="tight")
    print(f"  Saved: {OUT_DIR / filename}")
    plt.close(fig)


def plot_scale_range_comparison(range_df: pd.DataFrame, contrast_col: str,
                                title: str, filename: str):
    """Grouped bar chart: for each band, show unanimity at each scale range."""
    unan = compute_unanimity(range_df, "scale_range", contrast_col)
    scale_names = list(SCALE_RANGES.keys())

    fig, ax = plt.subplots(figsize=(9, 5))

    x = np.arange(len(BANDS))
    width = 0.25
    colors = ["#4393c3", "#f4a582", "#d6604d"]

    for s_idx, scale_name in enumerate(scale_names):
        vals = []
        for band in BANDS:
            row = unan[(unan["band"] == band) & (unan["scale_range"] == scale_name)]
            vals.append(row.iloc[0]["unanimity"] if not row.empty else 0)
        offset = (s_idx - 1) * width
        bars = ax.bar(x + offset, vals, width, label=scale_name.replace("\n", " "),
                       color=colors[s_idx], edgecolor="black", linewidth=0.5)
        # Annotate
        for b, v in zip(bars, vals):
            if v != 0:
                ax.text(b.get_x() + b.get_width() / 2, v + 0.1 * np.sign(v),
                        f"{int(v):+d}", ha="center", va="bottom" if v > 0 else "top",
                        fontsize=8, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([BAND_LABELS[b] for b in BANDS], fontsize=11)
    ax.set_ylabel("Unanimity score", fontsize=11)
    ax.set_ylim(-N_PATIENTS - 0.8, N_PATIENTS + 0.8)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.axhline(N_PATIENTS, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.axhline(-N_PATIENTS, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.legend(fontsize=9, loc="upper right")
    ax.set_title(title, fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / filename, bbox_inches="tight")
    print(f"  Saved: {OUT_DIR / filename}")
    plt.close(fig)


def find_best_scales(unan_single: pd.DataFrame, group_col: str = "k"):
    """Find k values with the most unanimous bands."""
    groups = sorted(unan_single[group_col].unique(),
                    key=lambda x: K_VALUES.index(x) if x in K_VALUES else 0)

    print("\n  Per-k summary (number of bands with |unanimity| >= 3):")
    best_k = None
    best_count = -1
    for grp in groups:
        sub = unan_single[unan_single[group_col] == grp]
        n_strong = (sub["unanimity"].abs() >= 3).sum()
        n_perfect = (sub["unanimity"].abs() == N_PATIENTS).sum()
        direction_summary = []
        for _, row in sub.iterrows():
            u = row["unanimity"]
            if u >= 3:
                direction_summary.append(f"{BAND_LABELS[row['band']]}+")
            elif u <= -3:
                direction_summary.append(f"{BAND_LABELS[row['band']]}-")
        dirs = ", ".join(direction_summary) if direction_summary else "none"
        print(f"    k={grp:>2d}: {n_strong} strong bands ({n_perfect} perfect) -- {dirs}")
        if n_strong > best_count:
            best_count = n_strong
            best_k = grp
    return best_k


def print_contrast_table(contrasts_df: pd.DataFrame, contrast_col: str, k: int):
    """Print a detailed table for a given k."""
    sub = contrasts_df[contrasts_df["k"] == k]
    print(f"\n  Detailed {contrast_col} values at k={k}:")
    print(f"  {'Band':<12s}", end="")
    for pat in PATIENTS:
        print(f"  {pat:>8s}", end="")
    print(f"  {'Mean':>8s}  {'Sign':>4s}")
    print("  " + "-" * 68)
    for band in BANDS:
        print(f"  {band:<12s}", end="")
        vals = []
        for pat in PATIENTS:
            row = sub[(sub["patient"] == pat) & (sub["band"] == band)]
            v = row.iloc[0][contrast_col] if not row.empty else np.nan
            vals.append(v)
            print(f"  {v:>+8.4f}", end="")
        m = np.nanmean(vals)
        n_pos = sum(1 for v in vals if v > 0)
        print(f"  {m:>+8.4f}  {n_pos}/{len(vals)}")


def main():
    print("=" * 70)
    print("Scale-Specific H2 Analysis")
    print("=" * 70)

    # Load data
    df = load_data()
    print(f"\nLoaded {len(df)} rows for {len(PATIENTS)} patients.")

    # --- Single-k contrasts ---
    print("\n--- Single-k contrasts ---")
    contrasts = compute_contrasts_single_k(df)

    # Unanimity for H2a and H2b
    unan_h2a = compute_unanimity(contrasts, "k", "H2a")
    unan_h2b = compute_unanimity(contrasts, "k", "H2b")

    # --- Heatmaps: unanimity ---
    print("\nPlotting unanimity heatmaps...")
    plot_heatmap(
        unan_h2a, "k", "H2a",
        "H2a: VI(Pre,Post) - VI(TT,Post)\nPositive = task persists in rest_post",
        "h2a_unanimity_heatmap.pdf",
        group_labels=[f"k={k}" for k in K_VALUES],
    )
    plot_heatmap(
        unan_h2b, "k", "H2b",
        "H2b: VI(Pre,TT) - VI(TT,Post)\nPositive = rest_post closer to TT than rest_pre was",
        "h2b_unanimity_heatmap.pdf",
        group_labels=[f"k={k}" for k in K_VALUES],
    )

    # --- Patient detail heatmaps ---
    print("\nPlotting patient-level detail...")
    plot_patient_detail(
        contrasts, "H2a",
        "H2a per patient: VI(Pre,Post) - VI(TT,Post)",
        "h2a_patient_detail.pdf",
    )
    plot_patient_detail(
        contrasts, "H2b",
        "H2b per patient: VI(Pre,TT) - VI(TT,Post)",
        "h2b_patient_detail.pdf",
    )

    # --- Scale ranges ---
    print("\n--- Scale-range contrasts (macro / meso / micro) ---")
    range_contrasts = compute_contrasts_scale_ranges(df)

    plot_scale_range_comparison(
        range_contrasts, "H2a",
        "H2a unanimity by scale range",
        "h2a_scale_range_bars.pdf",
    )
    plot_scale_range_comparison(
        range_contrasts, "H2b",
        "H2b unanimity by scale range",
        "h2b_scale_range_bars.pdf",
    )

    # Scale-range heatmaps
    unan_range_h2a = compute_unanimity(range_contrasts, "scale_range", "H2a")
    unan_range_h2b = compute_unanimity(range_contrasts, "scale_range", "H2b")

    plot_heatmap(
        unan_range_h2a, "scale_range", "H2a",
        "H2a: VI(Pre,Post) - VI(TT,Post) by scale range\nPositive = task persists",
        "h2a_scale_range_heatmap.pdf",
        group_labels=[n.replace("\n", " ") for n in SCALE_RANGES.keys()],
    )
    plot_heatmap(
        unan_range_h2b, "scale_range", "H2b",
        "H2b: VI(Pre,TT) - VI(TT,Post) by scale range\nPositive = rest_post closer to TT than rest_pre was",
        "h2b_scale_range_heatmap.pdf",
        group_labels=[n.replace("\n", " ") for n in SCALE_RANGES.keys()],
    )

    # --- Find best scales ---
    print("\n" + "=" * 70)
    print("Best scales for H2a:")
    best_k_h2a = find_best_scales(unan_h2a)
    print(f"\n  => Best single k for H2a: k={best_k_h2a}")
    print_contrast_table(contrasts, "H2a", best_k_h2a)

    print("\n" + "=" * 70)
    print("Best scales for H2b:")
    best_k_h2b = find_best_scales(unan_h2b)
    print(f"\n  => Best single k for H2b: k={best_k_h2b}")
    print_contrast_table(contrasts, "H2b", best_k_h2b)

    # --- Print all k tables for reference ---
    print("\n" + "=" * 70)
    print("Full H2a contrast tables:")
    for k in K_VALUES:
        print_contrast_table(contrasts, "H2a", k)

    print("\n" + "=" * 70)
    print("Full H2b contrast tables:")
    for k in K_VALUES:
        print_contrast_table(contrasts, "H2b", k)

    # --- Scale-range detail tables ---
    print("\n" + "=" * 70)
    print("Scale-range H2a detail:")
    for rng_name in SCALE_RANGES:
        sub = range_contrasts[range_contrasts["scale_range"] == rng_name]
        rng_label = rng_name.replace("\n", " ")
        print(f"\n  {rng_label}:")
        print(f"  {'Band':<12s}", end="")
        for pat in PATIENTS:
            print(f"  {pat:>8s}", end="")
        print(f"  {'Mean':>8s}  {'Sign':>4s}")
        print("  " + "-" * 68)
        for band in BANDS:
            print(f"  {band:<12s}", end="")
            vals = []
            for pat in PATIENTS:
                row = sub[(sub["patient"] == pat) & (sub["band"] == band)]
                v = row.iloc[0]["H2a"] if not row.empty else np.nan
                vals.append(v)
                print(f"  {v:>+8.4f}", end="")
            m = np.nanmean(vals)
            n_pos = sum(1 for v in vals if v > 0)
            print(f"  {m:>+8.4f}  {n_pos}/{len(vals)}")

    print("\n" + "=" * 70)
    print("Scale-range H2b detail:")
    for rng_name in SCALE_RANGES:
        sub = range_contrasts[range_contrasts["scale_range"] == rng_name]
        rng_label = rng_name.replace("\n", " ")
        print(f"\n  {rng_label}:")
        print(f"  {'Band':<12s}", end="")
        for pat in PATIENTS:
            print(f"  {pat:>8s}", end="")
        print(f"  {'Mean':>8s}  {'Sign':>4s}")
        print("  " + "-" * 68)
        for band in BANDS:
            print(f"  {band:<12s}", end="")
            vals = []
            for pat in PATIENTS:
                row = sub[(sub["patient"] == pat) & (sub["band"] == band)]
                v = row.iloc[0]["H2b"] if not row.empty else np.nan
                vals.append(v)
                print(f"  {v:>+8.4f}", end="")
            m = np.nanmean(vals)
            n_pos = sum(1 for v in vals if v > 0)
            print(f"  {m:>+8.4f}  {n_pos}/{len(vals)}")

    print("\n" + "=" * 70)
    print("Done. All figures saved to:", OUT_DIR)


if __name__ == "__main__":
    main()
