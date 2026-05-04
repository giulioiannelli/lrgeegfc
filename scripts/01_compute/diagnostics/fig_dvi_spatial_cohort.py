#!/usr/bin/env python3
"""Cohort heatmaps of Δ_VI(spatial scale) — one per spatial summary.

Mirrors the integer-k headline figure layout but with the spatial axis.
For each (band × spatial-bin) cell, colour = cohort patient-mean Δ_VI;
black border = ≥8/10 patients sign-correct in that bin; dashed grey =
6–7/10. Probe-bias zone (< 15 mm) shown with vertical hatching.

Reads:  data/audit/dvi_spatial_scale/dvi_spatial_scale_n10_imcoh_abs.csv
Writes: data/audit/dvi_spatial_scale/figures/fig_cohort_band_<summary>.pdf
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics import SPATIAL_SUMMARY_NAMES


COHORT_N10 = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
              "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
THRESH = 8
ADVISORY = 6
N_BINS = 24
X_LO_MM = 1.0
X_HI_MM = 100.0
SINGLETON_FRACTION_NOISE = 0.80   # criterion (A)


def _cohort_floor_ceiling(df, summary_col, max_pairwise_per_patient):
    """Return (single fine_floor_mm, single implant_ceiling_mm) cohort-wide.

    fine_floor: median over all (patient, band) of the noise-zone boundary —
                largest diameter at which `n_singletons/n_clusters ≥ 0.8`.
                A single cohort-wide scalar, NOT per band.
    implant_ceiling: cohort min over patients of `max pairwise distance
                     between electrodes` — the anatomy-only upper bound.
                     Above this diameter at least one patient cannot
                     contribute regardless of partition. Band-independent.
    """
    floors = []
    for (pat, band), ps in df.groupby(["patient", "band"]):
        ps = ps.sort_values(summary_col)
        sf = ps["n_singletons_post"] / ps["n_clusters_post"]
        mask = sf >= SINGLETON_FRACTION_NOISE
        if mask.any():
            floors.append(float(ps.loc[mask, summary_col].max()))
    fine = float(np.median(floors)) if floors else float(X_LO_MM)
    ceil = float(min(max_pairwise_per_patient.values()))
    return fine, ceil


def _max_pairwise_distance_per_patient():
    """Per-patient anatomical max distance (mm) from electrode coordinates.

    Reads `implant_pat_NN.csv` and `channel_labels.csv`, normalises labels,
    drops per-patient channel rows, converts μm→mm, computes max pairwise
    distance over valid electrodes.
    """
    import re as _re
    from scipy.spatial.distance import pdist
    from lrg_eegfc.config.const import PATIENT_CHANNEL_DROP
    from lrg_eegfc.config.paths import SEEG_DATAPATH

    def _norm(s):
        s = str(s).strip().strip('"').strip().split(",")[0].strip()
        return _re.sub(r"\s+", "", s).lower()

    out = {}
    for pat in COHORT_N10:
        pat_dir = SEEG_DATAPATH / pat
        patnum = int(pat.split("_")[-1])
        implant = pd.read_csv(pat_dir / f"implant_pat_{patnum:02d}.csv")
        labels = pd.read_csv(pat_dir / "channel_labels.csv")
        drops = PATIENT_CHANNEL_DROP.get(pat, {}).get("__labels__", [])
        if drops:
            labels = labels.drop(index=list(drops)).reset_index(drop=True)
        implant["_key"] = implant["label"].apply(_norm)
        labels["_key"] = labels[labels.columns[0]].apply(_norm)
        merged = labels.merge(implant[["_key", "x", "y", "z"]],
                              on="_key", how="left")
        for col in ("x", "y", "z"):
            merged[col] = (merged[col].astype(str)
                           .str.replace(",", ".", regex=False).astype(float))
        merged[["x", "y", "z"]] = merged[["x", "y", "z"]] / 1000.0
        coords = merged[["x", "y", "z"]].to_numpy(dtype=float)
        coords = coords[np.isfinite(coords).all(axis=1)]
        out[pat] = float(pdist(coords).max()) if len(coords) > 1 else float("nan")
    return out

OUT_DIR = ROOT / "data" / "audit" / "dvi_spatial_scale" / "figures"


def _cmap_div():
    return LinearSegmentedColormap.from_list(
        "trace_div",
        [(0.05, 0.20, 0.55), (1.0, 1.0, 1.0), (0.70, 0.05, 0.05)],
        N=256,
    )


def _draw_panel(ax, df: pd.DataFrame, summary: str, bins: np.ndarray,
                fine_floor: float, implant_ceil: float):
    summary_col = f"{summary}_post"
    B = len(BRAIN_BANDS_NAMES)
    n_b = len(bins) - 1
    mean = np.full((B, n_b), np.nan)
    npos = np.zeros((B, n_b), dtype=int)
    npat = np.zeros((B, n_b), dtype=int)

    for ib, band in enumerate(BRAIN_BANDS_NAMES):
        sub = df[df["band"] == band].copy()
        sub["bin_id"] = pd.cut(sub[summary_col], bins, include_lowest=True,
                                labels=False)
        for jb in range(n_b):
            cell = sub[sub["bin_id"] == jb]
            if cell.empty:
                continue
            pp = cell.groupby("patient")["d_VI"].mean().dropna()
            if pp.empty:
                continue
            mean[ib, jb] = float(pp.mean())
            npos[ib, jb] = int((pp > 0).sum())
            npat[ib, jb] = int(pp.size)

    cmap = _cmap_div()
    vmax = np.nanmax(np.abs(mean)) if np.isfinite(mean).any() else 1.0
    norm = plt.Normalize(vmin=-vmax, vmax=+vmax)

    for ib in range(B):
        for jb in range(n_b):
            v = mean[ib, jb]
            if not np.isfinite(v):
                continue
            x_lo = bins[jb]
            x_hi = bins[jb + 1]
            ax.add_patch(Rectangle((x_lo, ib - 0.5), x_hi - x_lo, 1,
                                    facecolor=cmap(norm(v)),
                                    edgecolor="none", zorder=1))
            n = npos[ib, jb]
            np_pat = npat[ib, jb]
            if np_pat >= 7:
                if n >= THRESH:
                    ax.add_patch(Rectangle((x_lo, ib - 0.5), x_hi - x_lo, 1,
                                            facecolor="none", edgecolor="black",
                                            lw=1.4, zorder=3))
                elif n >= ADVISORY:
                    ax.add_patch(Rectangle((x_lo, ib - 0.5), x_hi - x_lo, 1,
                                            facecolor="none", edgecolor="0.4",
                                            lw=0.8, ls=":", zorder=3))
            else:
                ax.add_patch(Rectangle((x_lo, ib - 0.5), x_hi - x_lo, 1,
                                        facecolor="none", edgecolor="0.6",
                                        lw=0.4, ls="--", alpha=0.5, zorder=3))
            color = "white" if abs(v) > 0.4 * vmax else "black"
            ax.text(np.sqrt(x_lo * x_hi), ib, f"{n}/{np_pat}",
                    ha="center", va="center", fontsize=6.5,
                    color=color)

    # Single cohort-wide noise floor and implant ceiling — band-independent.
    if fine_floor > bins[0]:
        ax.add_patch(Rectangle((bins[0], -0.5),
                                fine_floor - bins[0], B,
                                facecolor="none", edgecolor="#1f77b4",
                                hatch="\\\\\\\\", alpha=0.45,
                                linewidth=0.0, zorder=4))
        ax.axvline(fine_floor, color="#1f77b4", lw=1.4, ls="--",
                    alpha=0.9, zorder=5)
    if implant_ceil < bins[-1]:
        ax.add_patch(Rectangle((implant_ceil, -0.5),
                                bins[-1] - implant_ceil, B,
                                facecolor="none", edgecolor="#d62728",
                                hatch="////", alpha=0.45,
                                linewidth=0.0, zorder=4))
        ax.axvline(implant_ceil, color="#d62728", lw=1.4, ls="--",
                    alpha=0.9, zorder=5)

    ax.set_xscale("log")
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(B - 0.5, -0.5)
    ax.set_yticks(range(B))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                        fontsize=10)
    label = {
        "mean_diameter": "mean cluster diameter (mm) — log",
        "mean_radius": "mean cluster radius (mm) — log",
        "nn_centroid": "mean NN centroid distance (mm) — log",
        "size_weighted_radius": "size-weighted radius (mm) — log",
    }[summary]
    ax.set_xlabel(label, fontsize=10)
    return cmap, norm, vmax


def main():
    df = pd.read_csv(ROOT / "data" / "audit" / "dvi_spatial_scale"
                     / "dvi_spatial_scale_n10_imcoh_abs.csv")

    # Build a single 2x2 figure with all 4 summaries side-by-side, plus the
    # integer-k headline as a 5th panel for direct comparison.
    fig, axes = plt.subplots(2, 2, figsize=(15.0, 7.0), dpi=160,
                              gridspec_kw={"hspace": 0.42, "wspace": 0.18})

    print("Computing per-patient anatomical max pairwise distance (mm)...")
    max_pw = _max_pairwise_distance_per_patient()
    for pat, d in max_pw.items():
        print(f"  {pat}: max pairwise = {d:.1f} mm")
    cohort_implant_ceiling = float(min(max_pw.values()))
    print(f"  → cohort-min anatomical max pairwise = {cohort_implant_ceiling:.1f} mm")

    last = None
    for ax, summary in zip(axes.flat, SPATIAL_SUMMARY_NAMES):
        col = f"{summary}_post"
        bins = np.linspace(X_LO_MM, X_HI_MM, N_BINS + 1)
        fine, _ = _cohort_floor_ceiling(df, col, max_pw)
        print(f"  {summary}: cohort noise floor = {fine:.2f} mm, "
              f"cohort implant ceiling = {cohort_implant_ceiling:.2f} mm")
        last = _draw_panel(ax, df, summary, bins, fine, cohort_implant_ceiling)
        ax.set_title(f"$\\Delta_{{\\mathrm{{VI}}}}$ vs {summary}",
                     fontsize=11, loc="left")

    cm, nm, _ = last
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=nm, cmap=cm),
                        ax=list(axes.flat), shrink=0.65, pad=0.02,
                        fraction=0.025)
    cbar.set_label("patient-mean $\\Delta_{\\mathrm{VI}}$ in bin", fontsize=10)
    out = OUT_DIR / "fig_cohort_band_spatial_2x2.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")

    # Single-axis bigger figure for the headline summary (mean_diameter).
    fig2, ax = plt.subplots(figsize=(11.0, 4.5), dpi=160)
    col = "mean_diameter_post"
    bins = np.linspace(X_LO_MM, X_HI_MM, N_BINS + 1)
    fine, _ = _cohort_floor_ceiling(df, col, max_pw)
    cm, nm, vmax = _draw_panel(ax, df, "mean_diameter", bins,
                                fine, cohort_implant_ceiling)
    ax.set_title("$\\Delta_{\\mathrm{VI}}$ vs mean cluster diameter — c_post-anchored, n=10",
                 fontsize=11, loc="left")
    cbar = fig2.colorbar(plt.cm.ScalarMappable(norm=nm, cmap=cm),
                         ax=ax, shrink=0.85, pad=0.02, fraction=0.04)
    cbar.set_label("patient-mean $\\Delta_{\\mathrm{VI}}$", fontsize=10)
    out2 = OUT_DIR / "fig_cohort_band_mean_diameter_headline.pdf"
    fig2.savefig(out2, bbox_inches="tight")
    plt.close(fig2)
    print(f"saved {out2}")


if __name__ == "__main__":
    main()
