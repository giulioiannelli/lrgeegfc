#!/usr/bin/env python3
"""Audit Step 15 — Anatomy-aware MSPC cohort analysis.

User concern: per-patient variability in trace counts may be explained
by *electrode placement* (each implant samples a different anatomy).
This step joins MSPC leaf-pattern labels with Desikan-Killiany regions
parsed from each patient's `implant_pat_NN.csv`, then aggregates across
the cohort to ask: do certain anatomical regions consistently host
trace structures in specific bands?

White-matter and unmapped contacts are dropped (no neural meaning).

Outputs:
  data/audit/per_patient_hierarchy_mspc/anatomy/
    leaf_x_region.csv               — leaf × region table (joined)
    cohort_lobe_band_pattern.csv    — counts per (lobe, band, pattern)
    cohort_region_band_pattern.csv  — counts per (DK region, band, pattern)
    fig_lobe_band_trace.pdf         — heatmap: lobe × band, trace fraction
    fig_lobe_band_persist.pdf       — same for persist
    fig_lobe_band_stacked.pdf       — stacked-bar per lobe across bands
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_08_per_patient_hierarchy import (  # noqa: E402
    BAND_ORDER, PATIENTS_4PHASE, PATTERN_COLOR, BRAIN_BAND_TEX_DICT,
)
from lrg_eegfc.utils.io.regions import load_channel_regions  # noqa: E402

OUT_DIR = ROOT / "data" / "audit" / "per_patient_hierarchy_mspc" / "anatomy"
LOBES_KEEP = ["frontal", "parietal", "temporal", "occipital",
              "cingulate", "insula", "subcortical"]
PATTERNS = ["trace", "persist", "reset", "rearrange"]


def join_leaves_with_regions() -> pd.DataFrame:
    """Read MSPC leaf_assignment.csv and join each row with its DK region."""
    leaves = pd.read_csv(ROOT / "data" / "audit" /
                          "per_patient_hierarchy_mspc" / "leaf_assignment.csv")
    region_dfs = []
    for pat in leaves.patient.unique():
        rdf = load_channel_regions(pat).reset_index().rename(
            columns={"index": "leaf_id"})
        rdf["patient"] = pat
        region_dfs.append(rdf[["patient", "leaf_id", "label", "region",
                                  "lobe", "hemisphere"]])
    regions = pd.concat(region_dfs, ignore_index=True)
    return leaves.merge(regions, on=["patient", "leaf_id"], how="left")


def cohort_lobe_band(joined: pd.DataFrame) -> pd.DataFrame:
    sub = joined[joined.lobe.isin(LOBES_KEEP)].copy()
    grouped = (sub.groupby(["lobe", "band", "dominant"])
                  .size().unstack(fill_value=0).reset_index())
    for p in PATTERNS:
        if p not in grouped.columns:
            grouped[p] = 0
    grouped["total"] = grouped[PATTERNS].sum(axis=1)
    for p in PATTERNS:
        grouped[f"frac_{p}"] = grouped[p] / grouped["total"].replace(0, np.nan)
    return grouped


def cohort_region_band(joined: pd.DataFrame, min_n: int = 30) -> pd.DataFrame:
    sub = joined[joined.lobe.isin(LOBES_KEEP)].copy()
    grouped = (sub.groupby(["region", "band", "dominant"])
                  .size().unstack(fill_value=0).reset_index())
    for p in PATTERNS:
        if p not in grouped.columns:
            grouped[p] = 0
    grouped["total"] = grouped[PATTERNS].sum(axis=1)
    for p in PATTERNS:
        grouped[f"frac_{p}"] = grouped[p] / grouped["total"].replace(0, np.nan)
    region_total = grouped.groupby("region")["total"].sum()
    keep = region_total[region_total >= min_n].index
    return grouped[grouped.region.isin(keep)].reset_index(drop=True)


def _patients_per_lobe(joined: pd.DataFrame) -> pd.Series:
    sub = joined[joined.lobe.isin(LOBES_KEEP)]
    return sub.groupby("lobe")["patient"].nunique()


def _heatmap(ax, M: np.ndarray, lobes: list[str], bands: list[str],
              title: str, cmap: str = "Reds", vmin: float = 0.0,
              vmax: float | None = None) -> None:
    if vmax is None:
        vmax = float(np.nanmax(M))
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0
    im = ax.imshow(M, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                    interpolation="nearest")
    im.set_rasterized(True)
    ax.set_xticks(range(len(bands)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in bands],
                        rotation=20, fontsize=9)
    ax.set_yticks(range(len(lobes)))
    ax.set_yticklabels(lobes, fontsize=9)
    ax.set_title(title, fontsize=11)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if not np.isfinite(v):
                continue
            txt = f"{v:.2f}"
            ax.text(j, i, txt, ha="center", va="center",
                     fontsize=7, color="black" if v < 0.5 * vmax else "white")
    plt.colorbar(im, ax=ax, fraction=0.04, pad=0.03)


def _matrix(grouped: pd.DataFrame, lobes: list[str], bands: list[str],
              col: str) -> np.ndarray:
    M = np.full((len(lobes), len(bands)), np.nan, dtype=float)
    for i, lobe in enumerate(lobes):
        for j, band in enumerate(bands):
            row = grouped[(grouped.lobe == lobe) & (grouped.band == band)]
            if not row.empty:
                M[i, j] = float(row[col].values[0])
    return M


def figure_lobe_band(grouped: pd.DataFrame, pat_per_lobe: pd.Series) -> None:
    lobes = [l for l in LOBES_KEEP if l in grouped.lobe.values]
    lobes_label = [f"{l} (n={pat_per_lobe.get(l, 0)})" for l in lobes]
    bands = BAND_ORDER

    for pattern in PATTERNS:
        M = _matrix(grouped, lobes, bands, f"frac_{pattern}")
        fig, ax = plt.subplots(figsize=(7.2, 4.2))
        cmap = {
            "trace": "Reds", "persist": "Greys",
            "reset": "Greens", "rearrange": "Blues",
        }[pattern]
        _heatmap(ax, M, lobes_label, bands,
                  title=f"{pattern.upper()} fraction by lobe × band "
                        f"(cohort, white-matter excluded)",
                  cmap=cmap, vmin=0.0, vmax=1.0)
        ax.set_xlabel("band"); ax.set_ylabel("anatomical lobe (cohort patients)")
        out = OUT_DIR / f"fig_lobe_band_{pattern}.pdf"
        fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def figure_lobe_stacked(grouped: pd.DataFrame, pat_per_lobe: pd.Series) -> None:
    lobes = [l for l in LOBES_KEEP if l in grouped.lobe.values]
    bands = BAND_ORDER
    n_lobes = len(lobes)
    fig, axes = plt.subplots(1, n_lobes, figsize=(3.2 * n_lobes, 3.6),
                                sharey=True)
    if n_lobes == 1:
        axes = [axes]
    for ax, lobe in zip(axes, lobes):
        bottom = np.zeros(len(bands))
        x = np.arange(len(bands))
        for pattern in PATTERNS:
            vals = np.array([
                float(grouped[(grouped.lobe == lobe) &
                                 (grouped.band == b)][f"frac_{pattern}"].values[0])
                if not grouped[(grouped.lobe == lobe) &
                                 (grouped.band == b)].empty else 0.0
                for b in bands
            ])
            ax.bar(x, vals, bottom=bottom,
                    color=PATTERN_COLOR[pattern], label=pattern.upper(),
                    edgecolor="white", linewidth=0.4)
            bottom += vals
        n_pat = pat_per_lobe.get(lobe, 0)
        ax.set_title(f"{lobe} (n={n_pat})", fontsize=10)
        ax.set_xticks(x)
        ax.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in bands],
                            rotation=30, fontsize=8)
        ax.set_ylim(0, 1.0)
        ax.set_xlabel("band", fontsize=9)
    axes[0].set_ylabel("pattern fraction", fontsize=10)
    axes[-1].legend(loc="upper right", fontsize=7, frameon=False,
                      bbox_to_anchor=(1.05, 1.05))
    fig.suptitle("MSPC pattern composition by lobe × band (cohort, "
                  "white-matter excluded)", fontsize=11)
    fig.tight_layout()
    out = OUT_DIR / "fig_lobe_band_stacked.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def figure_top_regions(grouped_region: pd.DataFrame) -> None:
    """For each band, show top-10 DK regions by trace fraction."""
    bands = BAND_ORDER
    n_bands = len(bands)
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=False)
    axes = axes.flatten()
    for ax, band in zip(axes, bands):
        sub = grouped_region[grouped_region.band == band].copy()
        if sub.empty:
            ax.set_title(f"{BRAIN_BAND_TEX_DICT.get(band, band)}: no data",
                          fontsize=10)
            ax.axis("off")
            continue
        sub = sub.sort_values("frac_trace", ascending=False).head(10)
        y = np.arange(len(sub))
        ax.barh(y, sub.frac_trace.values, color=PATTERN_COLOR["trace"])
        ax.set_yticks(y)
        ax.set_yticklabels([r.replace("ctx-", "") for r in sub.region],
                            fontsize=7)
        ax.invert_yaxis()
        ax.set_xlim(0, 1.0)
        ax.set_xlabel("trace fraction", fontsize=8)
        ax.set_title(f"{BRAIN_BAND_TEX_DICT.get(band, band)} — top regions",
                      fontsize=10)
        for yi, (frac, n) in enumerate(zip(sub.frac_trace.values,
                                                  sub.total.values)):
            ax.text(frac + 0.02, yi, f"n={int(n)}", fontsize=6, va="center")
    fig.suptitle("Top-10 DK regions by TRACE fraction per band "
                  "(cohort, regions with ≥30 leaves total)", fontsize=11)
    fig.tight_layout()
    out = OUT_DIR / "fig_top_regions_trace.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def figure_per_patient_lobe_trace(joined: pd.DataFrame) -> None:
    """For each (lobe, band), show every patient's trace fraction as a dot.

    Reveals whether the cohort-mean trace is driven by 1-2 patients or by
    a consistent majority. Black square = cohort mean; coloured dots =
    individual patients (only those with ≥5 leaves in that lobe).
    """
    sub = joined[joined.lobe.isin(LOBES_KEEP)].copy()
    g = sub.groupby(["patient", "lobe", "band"])
    totals = g.size().reset_index(name="total")
    tr = (sub[sub.dominant == "trace"]
            .groupby(["patient", "lobe", "band"]).size()
            .reset_index(name="n_trace"))
    pp = totals.merge(tr, on=["patient", "lobe", "band"], how="left").fillna(0)
    pp["frac_trace"] = pp["n_trace"] / pp["total"]
    pp = pp[pp.total >= 5]

    lobes = [l for l in LOBES_KEEP if l in pp.lobe.values]
    bands = BAND_ORDER
    fig, axes = plt.subplots(len(lobes), 1, figsize=(8.5, 1.2 * len(lobes)),
                                sharex=True)
    if len(lobes) == 1:
        axes = [axes]
    rng = np.random.default_rng(0)
    for ax, lobe in zip(axes, lobes):
        for j, band in enumerate(bands):
            s = pp[(pp.lobe == lobe) & (pp.band == band)]
            if s.empty:
                continue
            xs = j + (rng.random(len(s)) - 0.5) * 0.25
            ax.scatter(xs, s.frac_trace, color="#d62728", s=22,
                        alpha=0.7, edgecolors="white", linewidths=0.5,
                        zorder=3)
            ax.scatter(j, s.frac_trace.mean(), color="black", marker="s",
                        s=44, zorder=4)
        ax.set_xlim(-0.5, len(bands) - 0.5)
        ax.set_ylim(-0.02, 1.0)
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_ylabel(f"{lobe}\n(n_pat={pp[pp.lobe==lobe].patient.nunique()})",
                        fontsize=8)
        ax.grid(axis="y", alpha=0.25, linestyle=":")
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
    axes[-1].set_xticks(range(len(bands)))
    axes[-1].set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in bands],
                                fontsize=9)
    axes[-1].set_xlabel("band")
    fig.suptitle("Per-patient TRACE fraction (red dots) per lobe × band — "
                  "black squares = cohort mean", fontsize=10)
    fig.tight_layout()
    out = OUT_DIR / "fig_per_patient_lobe_trace.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Joining MSPC leaves with DK regions ...", flush=True)
    joined = join_leaves_with_regions()
    joined.to_csv(OUT_DIR / "leaf_x_region.csv", index=False)
    print(f"  {len(joined)} leaf rows; "
          f"{joined.region.eq('unknown').sum()} unmapped, "
          f"{joined.lobe.eq('white_matter').sum()} white-matter")

    grouped_lobe = cohort_lobe_band(joined)
    grouped_region = cohort_region_band(joined, min_n=30)
    grouped_lobe.to_csv(OUT_DIR / "cohort_lobe_band_pattern.csv", index=False)
    grouped_region.to_csv(OUT_DIR / "cohort_region_band_pattern.csv",
                            index=False)

    pat_per_lobe = _patients_per_lobe(joined)
    print("\nPatients per lobe:")
    print(pat_per_lobe.to_string())

    print("\nWriting figures ...")
    figure_lobe_band(grouped_lobe, pat_per_lobe)
    figure_lobe_stacked(grouped_lobe, pat_per_lobe)
    figure_top_regions(grouped_region)
    figure_per_patient_lobe_trace(joined)
    print(f"\nDone. See {OUT_DIR.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
