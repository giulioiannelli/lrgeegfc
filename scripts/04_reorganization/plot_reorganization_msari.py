#!/usr/bin/env python3
"""Reorganization analysis using multiscale ARI with ALL patients.

Uses the multiscale structural distance: D = integral(1 - ARI(h)) dh
across normalized dendrogram heights.

Patients and their available phases:
  Pat_02, Pat_03, Pat_05, Pat_08: rest_pre, task_learn, task_test, rest_post (4 phases)
  Pat_07: rest_pre, task_learn, rest_post (3 phases, missing task_test)
  Pat_06: rest_pre, rest_post (2 rest phases only — no task data)

For the RI computation:
  REST = {rest_pre, rest_post},  TASK = available task phases
  D_within = mean of within-group distances
  D_cross  = mean of cross-group distances
  RI = D_cross / D_within

Pat_06 has no task phases, so only D(rest_pre, rest_post) is reported.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import adjusted_rand_score

from lrg_eegfc.config import BRAIN_BANDS_NAMES, PHASE_LABELS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

# ---------------------------------------------------------------------------
# Patient definitions with available phases
# ---------------------------------------------------------------------------
PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
PATIENTS = list(PATIENT_PHASES.keys())
BANDS = BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import LRG_CACHE
OUTPUT_DIR = FIGURES_ROOT / "metric_exploration"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

BAND_TEX = [BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS]

REST_PHASES = {"rest_pre", "rest_post"}
TASK_PHASES = {"task_learn", "task_test"}


# ---------------------------------------------------------------------------
# Multiscale ARI distance (the user's metric)
# ---------------------------------------------------------------------------
def multiscale_ari_distance(Z1, Z2, n_nodes, n_thresholds=100):
    """Compute D = integral(1 - ARI(h)) dh across normalized heights.

    At each normalized threshold t in [0.01, 1.0], cut both dendrograms
    at t * max_height and compute ARI. Then integrate (1-ARI).
    """
    max_h1 = Z1[:, 2].max()
    max_h2 = Z2[:, 2].max()
    thresholds = np.linspace(0.01, 1.0, n_thresholds)
    ari_vals = []
    for t in thresholds:
        l1 = fcluster(Z1, t=t * max_h1, criterion="distance")
        l2 = fcluster(Z2, t=t * max_h2, criterion="distance")
        ari_vals.append(adjusted_rand_score(l1, l2))
    ari_vals = np.array(ari_vals)
    return float(np.trapz(1 - ari_vals, thresholds) / (thresholds[-1] - thresholds[0]))


# ---------------------------------------------------------------------------
# Load and compute
# ---------------------------------------------------------------------------
def load_all_lrg():
    """Load all available LRG results."""
    data = {}
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for phase in phases:
                lrg = load_lrg_result(pat, phase, band, "msc", cache_root=LRG_CACHE)
                if lrg is not None:
                    data[(pat, phase, band)] = lrg
    return data


def compute_all_distances(lrg_data):
    """Compute pairwise MS-ARI distances for each patient/band."""
    dists = {}
    for pat, phases in PATIENT_PHASES.items():
        pairs = list(combinations(phases, 2))
        for band in BANDS:
            for p1, p2 in pairs:
                k1 = (pat, p1, band)
                k2 = (pat, p2, band)
                if k1 in lrg_data and k2 in lrg_data:
                    d = multiscale_ari_distance(
                        lrg_data[k1].linkage_matrix,
                        lrg_data[k2].linkage_matrix,
                        lrg_data[k1].n_nodes,
                    )
                    dists[(pat, band, (p1, p2))] = d
                    print(f"  {pat}/{band}/{p1}↔{p2}: {d:.4f}")
    return dists


def normalize_per_patient(dists):
    """Min-max normalize distances within each patient to [0,1]."""
    normed = {}
    for pat in PATIENTS:
        vals = [v for (p, b, pair), v in dists.items() if p == pat]
        if not vals:
            continue
        vmin, vmax = min(vals), max(vals)
        span = vmax - vmin if vmax > vmin else 1e-12
        for (p, b, pair), v in dists.items():
            if p == pat:
                normed[(p, b, pair)] = (v - vmin) / span
    return normed


def classify_pair(pair):
    """Classify a phase pair as within-rest, within-task, or cross."""
    s = set(pair)
    if s <= REST_PHASES:
        return "within-rest"
    elif s <= TASK_PHASES:
        return "within-task"
    else:
        return "cross"


def compute_ri(dists, pat, band):
    """Compute RI = D_cross / D_within for a patient/band.

    Adapts to available phases:
    - 4 phases: D_within = mean(D_rest, D_task), D_cross = mean of 4 cross pairs
    - 3 phases (no task_test): D_within = D_rest only, D_cross = mean of 2 cross pairs
    - 2 rest phases only: cannot compute RI (returns NaN)
    """
    phases = PATIENT_PHASES[pat]
    pairs = list(combinations(phases, 2))

    within_vals = []
    cross_vals = []
    for pair in pairs:
        key = (pat, band, pair)
        if key not in dists:
            continue
        cat = classify_pair(pair)
        if cat in ("within-rest", "within-task"):
            within_vals.append(dists[key])
        elif cat == "cross":
            cross_vals.append(dists[key])

    if not within_vals or not cross_vals:
        return np.nan

    d_within = np.mean(within_vals)
    d_cross = np.mean(cross_vals)
    if d_within == 0:
        return np.nan
    return d_cross / d_within


def build_dataframe(raw, normed):
    """Build results DataFrame with raw and normalized RI + distances."""
    rows = []
    for pat in PATIENTS:
        phases = PATIENT_PHASES[pat]
        pairs = list(combinations(phases, 2))
        for band in BANDS:
            ri_raw = compute_ri(raw, pat, band)
            ri_norm = compute_ri(normed, pat, band)

            # Aggregate distances by category
            d = {"within-rest": [], "within-task": [], "cross": []}
            d_norm = {"within-rest": [], "within-task": [], "cross": []}
            for pair in pairs:
                cat = classify_pair(pair)
                if (pat, band, pair) in raw:
                    d[cat].append(raw[(pat, band, pair)])
                    d_norm[cat].append(normed[(pat, band, pair)])

            rows.append({
                "patient": pat,
                "band": band,
                "n_phases": len(phases),
                "RI_raw": ri_raw,
                "RI_norm": ri_norm,
                "D_rest_raw": np.mean(d["within-rest"]) if d["within-rest"] else np.nan,
                "D_task_raw": np.mean(d["within-task"]) if d["within-task"] else np.nan,
                "D_cross_raw": np.mean(d["cross"]) if d["cross"] else np.nan,
                "D_rest_norm": np.mean(d_norm["within-rest"]) if d_norm["within-rest"] else np.nan,
                "D_task_norm": np.mean(d_norm["within-task"]) if d_norm["within-task"] else np.nan,
                "D_cross_norm": np.mean(d_norm["cross"]) if d_norm["cross"] else np.nan,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
PAT_COLORS = {
    "Pat_02": "#66c2a5", "Pat_03": "#fc8d62", "Pat_05": "#8da0cb",
    "Pat_06": "#a6d854", "Pat_07": "#ffd92f", "Pat_08": "#e78ac3",
}
PAT_MARKERS = {
    "Pat_02": "o", "Pat_03": "s", "Pat_05": "D",
    "Pat_06": "^", "Pat_07": "v", "Pat_08": "P",
}


def plot_ri_summary(df):
    """Main RI summary figure — all patients."""
    # Separate RI-capable patients from Pat_06
    df_ri = df[df["n_phases"] > 2].copy()
    ri_patients = sorted(df_ri["patient"].unique())
    df_06 = df[df["patient"] == "Pat_06"]

    fig, axes = plt.subplots(1, 3, figsize=(22, 7))
    x = np.arange(len(BANDS))

    # --- Panel A: RI per band (all patients with RI) ---
    ax = axes[0]
    offsets = np.linspace(-0.2, 0.2, len(ri_patients))
    for i, pat in enumerate(ri_patients):
        pdf = df_ri[df_ri["patient"] == pat]
        ris = [pdf[pdf["band"] == b]["RI_norm"].values[0] for b in BANDS]
        ax.scatter(x + offsets[i], ris, s=80, color=PAT_COLORS[pat],
                   marker=PAT_MARKERS[pat], label=pat,
                   edgecolors="black", linewidths=0.5, zorder=5)

    # Band means + std
    for j, band in enumerate(BANDS):
        vals = df_ri[df_ri["band"] == band]["RI_norm"].values
        m, s = np.nanmean(vals), np.nanstd(vals)
        ax.errorbar(j, m, yerr=s, fmt="s", color="black", ms=8,
                    capsize=5, capthick=2, lw=2, zorder=6)

    ax.axhline(1.0, color="gray", ls="--", lw=1.5, alpha=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(BAND_TEX, fontsize=12)
    ax.set_ylabel("Reorganization Index (normalized)", fontsize=12)
    ax.set_title("(A) RI per band — all patients\n"
                 "RI = D(rest↔task) / D(within-condition)", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.2, axis="y")

    # --- Panel B: Normalized within vs cross ---
    ax = axes[1]
    width = 0.35
    wi_means, cr_means = [], []
    wi_stds, cr_stds = [], []
    for band in BANDS:
        bdf = df_ri[df_ri["band"] == band]
        # Within = mean of rest and task distances per patient
        wi = []
        for _, row in bdf.iterrows():
            vals = [v for v in [row["D_rest_norm"], row["D_task_norm"]] if not np.isnan(v)]
            if vals:
                wi.append(np.mean(vals))
        cr = bdf["D_cross_norm"].dropna().values
        wi_means.append(np.mean(wi) if wi else 0)
        wi_stds.append(np.std(wi) if wi else 0)
        cr_means.append(np.mean(cr) if len(cr) else 0)
        cr_stds.append(np.std(cr) if len(cr) else 0)

    ax.bar(x - width / 2, wi_means, width, yerr=wi_stds,
           label="Within-condition", color="#2166ac", alpha=0.8, capsize=4)
    ax.bar(x + width / 2, cr_means, width, yerr=cr_stds,
           label="Cross-condition", color="#b2182b", alpha=0.8, capsize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(BAND_TEX, fontsize=12)
    ax.set_ylabel("Normalized MS-ARI distance", fontsize=12)
    ax.set_title("(B) Within vs Cross distances\n(per-patient normalized)", fontsize=12, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(alpha=0.2, axis="y")

    # --- Panel C: Pat_06 rest_pre vs rest_post ---
    ax = axes[2]
    vals_06 = []
    for band in BANDS:
        key_raw = df_06[df_06["band"] == band]["D_rest_raw"].values
        vals_06.append(key_raw[0] if len(key_raw) else np.nan)
    bars = ax.bar(x, vals_06, color="#a6d854", edgecolor="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(BAND_TEX, fontsize=12)
    ax.set_ylabel("MS-ARI distance (rest_pre ↔ rest_post)", fontsize=12)
    ax.set_title("(C) Pat_06 — rest_pre vs rest_post\n(only rest phases available)",
                 fontsize=12, fontweight="bold")
    ax.grid(alpha=0.2, axis="y")
    for i, v in enumerate(vals_06):
        if not np.isnan(v):
            ax.text(i, v + 0.005, f"{v:.3f}", ha="center", va="bottom", fontsize=9)

    fig.suptitle(
        "Cross-Phase Reorganization — Multiscale ARI (all patients)",
        fontsize=15, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    path = OUTPUT_DIR / "reorganization_msari_summary.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


def plot_per_patient_heatmaps(raw, normed):
    """Per-patient phase distance heatmaps (normalized)."""
    for pat in PATIENTS:
        phases = PATIENT_PHASES[pat]
        n = len(phases)
        pairs = list(combinations(phases, 2))

        fig, axes = plt.subplots(2, 3, figsize=(20, 12))
        axes_flat = axes.ravel()

        for j, band in enumerate(BANDS):
            ax = axes_flat[j]
            matrix = np.zeros((n, n))
            for i1, p1 in enumerate(phases):
                for i2, p2 in enumerate(phases):
                    if i1 < i2:
                        pair = (p1, p2)
                        val = normed.get((pat, band, pair), np.nan)
                        matrix[i1, i2] = val
                        matrix[i2, i1] = val

            mask = np.eye(n, dtype=bool)
            display = np.where(mask, np.nan, matrix)
            im = ax.imshow(display, cmap="YlOrRd", vmin=0, vmax=1, aspect="equal")
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.06, shrink=0.85)

            for i1 in range(n):
                for i2 in range(n):
                    if i1 != i2:
                        val = matrix[i1, i2]
                        color = "white" if val > 0.55 else "black"
                        ax.text(i2, i1, f"{val:.2f}", ha="center", va="center",
                                fontsize=11, fontweight="bold", color=color)

            ax.set_xticks(range(n))
            ax.set_yticks(range(n))
            ax.set_xticklabels(phases, rotation=45, ha="right", fontsize=9)
            ax.set_yticklabels(phases, fontsize=9)
            ax.set_title(BRAIN_BAND_TEX_DICT.get(band, band),
                         fontsize=13, fontweight="bold")

        fig.suptitle(
            f"{pat} — Normalized MS-ARI phase distances\n"
            f"({n} phases: {', '.join(phases)})",
            fontsize=14, fontweight="bold", y=1.01,
        )
        plt.tight_layout()
        path = OUTPUT_DIR / f"heatmaps_msari_{pat}.pdf"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        fig.savefig(path.with_suffix(".png"), dpi=100, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {path}")


def plot_ri_heatmap(df):
    """Patient × Band RI heatmap."""
    df_ri = df[df["n_phases"] > 2].copy()
    ri_patients = sorted(df_ri["patient"].unique())

    fig, ax = plt.subplots(figsize=(10, 5))
    ri_matrix = np.zeros((len(ri_patients), len(BANDS)))
    for i, pat in enumerate(ri_patients):
        for j, band in enumerate(BANDS):
            row = df_ri[(df_ri["patient"] == pat) & (df_ri["band"] == band)]
            ri_matrix[i, j] = row["RI_norm"].values[0] if len(row) else np.nan

    im = ax.imshow(ri_matrix, cmap="RdBu_r", vmin=0.5, vmax=2.5, aspect="auto")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="RI")

    for i in range(len(ri_patients)):
        for j in range(len(BANDS)):
            val = ri_matrix[i, j]
            color = "white" if abs(val - 1.0) > 0.5 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=10, fontweight="bold", color=color)

    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels(BAND_TEX, fontsize=12)
    ax.set_yticks(range(len(ri_patients)))
    ax.set_yticklabels(ri_patients, fontsize=12)
    ax.set_title("Reorganization Index — MS-ARI (per-patient normalized)\n"
                 "Red = reorganization (RI > 1)  |  Blue = persistence (RI < 1)",
                 fontsize=13, fontweight="bold")

    plt.tight_layout()
    path = OUTPUT_DIR / "ri_heatmap_msari.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("Loading all LRG results...")
    lrg_data = load_all_lrg()
    print(f"  Loaded {len(lrg_data)} entries")

    print("\nComputing pairwise MS-ARI distances...")
    raw = compute_all_distances(lrg_data)

    print("\nNormalizing per patient...")
    normed = normalize_per_patient(raw)

    print("\nBuilding results table...")
    df = build_dataframe(raw, normed)
    csv_path = OUTPUT_DIR / "reorganization_msari_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")

    print("\nGenerating figures...")
    plot_ri_summary(df)
    plot_per_patient_heatmaps(raw, normed)
    plot_ri_heatmap(df)

    # Summary table
    print("\n" + "=" * 80)
    print("MS-ARI REORGANIZATION INDEX — SUMMARY")
    print("=" * 80)

    df_ri = df[df["n_phases"] > 2]
    print(f"\nPatients with RI: {sorted(df_ri['patient'].unique())}")
    print(f"Pat_06: rest-only (no RI, see D_rest column)\n")

    pivot = df_ri.pivot(index="patient", columns="band", values="RI_norm")
    pivot = pivot[BANDS]
    print("RI per patient/band (normalized):")
    print(pivot.round(3).to_string())

    print(f"\n{'Band':<12} {'Mean':>7} {'Median':>7} {'Std':>7} {'CV':>7} "
          f"{'Min':>7} {'Max':>7} {'All>1':>6}")
    print("-" * 68)
    for band in BANDS:
        vals = pivot[band].dropna().values
        m = np.nanmean(vals)
        med = np.nanmedian(vals)
        s = np.nanstd(vals)
        cv = s / m if m > 0 else 0
        print(f"{band:<12} {m:>7.3f} {med:>7.3f} {s:>7.3f} {cv:>7.3f} "
              f"{np.nanmin(vals):>7.3f} {np.nanmax(vals):>7.3f} "
              f"{str(np.all(vals > 1)):>6}")

    # Pat_06 separate
    print(f"\nPat_06 — D(rest_pre, rest_post) per band (MS-ARI distance):")
    for band in BANDS:
        row = df[(df["patient"] == "Pat_06") & (df["band"] == band)]
        d = row["D_rest_raw"].values[0] if len(row) else np.nan
        print(f"  {band:<12} {d:.4f}")


if __name__ == "__main__":
    main()
