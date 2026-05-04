#!/usr/bin/env python3
"""Δ_VI(spatial scale) — replace the integer-k axis with a physical scale.

Builds 4 spatial-scale summaries of the c_post(k) partition per (patient,
band, k), then plots Δ_VI(k) against each of the 4 summaries cohort-wide.
The axis becomes physical (mm) and is biologically interpretable.

Reference partition for the spatial summary: **c_post** (per user
2026-04-27 — "the trace claim's destination phase").

Reads:
  data/raw/stereoeeg_patients/Pat_NN/implant_pat_NN.csv  (label, x, y, z)
  data/raw/stereoeeg_patients/Pat_NN/channel_labels.csv  (label)
  data/cache/imcoh_lrg/Pat_NN/<band>_<phase>_lrg_imcoh-abs.npz

Writes:
  data/audit/dvi_spatial_scale/dvi_spatial_scale_n10_imcoh_abs.csv
  data/audit/dvi_spatial_scale/coords_alignment_summary.md
  data/audit/dvi_spatial_scale/figures/fig_dvi_vs_<summary>.pdf  (×4)
  data/audit/dvi_spatial_scale/figures/fig_diameter_vs_k.pdf
  data/audit/dvi_spatial_scale/figures/SUBSTANTIATION.md

Coordinates are converted from μm → mm by dividing by 1000.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST,
    PATIENT_CHANNEL_DROP,
)
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, SEEG_DATAPATH
from lrg_eegfc.utils.metrics import (
    compute_vi, cluster_spatial_scale, SPATIAL_SUMMARY_NAMES,
)
from lrg_eegfc.workflow.lrg import load_lrg_result


COHORT_N10 = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
              "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
# Per-patient cap at N-1 — fine resolutions until every leaf is its own cluster.
# Cohort N values: 113 (Pat_10) – 122 (Pat_03).
K_MAX_GLOBAL = 122
COORD_UNIT_DIVISOR = 1000.0  # μm → mm

OUT_DIR = ROOT / "data" / "audit" / "dvi_spatial_scale"
FIG_DIR = OUT_DIR / "figures"


def _normalize_label(s: str) -> str:
    """Normalize sEEG label for cross-file merge.

    Strips whitespace, drops grid-id suffixes after a comma, lowercases.
    Examples:
        '"A 1,G2"'  → 'a1'
        'A1'        → 'a1'
        'P  3'      → 'p3'
    """
    s = str(s).strip().strip('"').strip()
    s = s.split(",")[0].strip()
    s = re.sub(r"\s+", "", s)
    return s.lower()


def load_coords_mm(patient: str) -> tuple[np.ndarray, np.ndarray, int]:
    """Return (coords_mm, valid_mask, n_dropped).

    `coords_mm` has shape (N, 3) — entries for unmatched contacts are NaN.
    `valid_mask` is True where coords are present. `n_dropped` is the
    count of NaN rows.
    """
    pat_dir = SEEG_DATAPATH / patient
    patnum = int(patient.split("_")[-1])
    implant_csv = pat_dir / f"implant_pat_{patnum:02d}.csv"
    channel_csv = pat_dir / "channel_labels.csv"
    implant = pd.read_csv(implant_csv)
    labels = pd.read_csv(channel_csv)

    # apply per-patient channel drops (e.g., Pat_10 [53,54,55])
    drops = PATIENT_CHANNEL_DROP.get(patient, {}).get("__labels__", [])
    if drops:
        labels = labels.drop(index=list(drops)).reset_index(drop=True)

    implant["_key"] = implant["label"].apply(_normalize_label)
    labels["_key"] = labels[labels.columns[0]].apply(_normalize_label)
    merged = labels.merge(implant[["_key", "x", "y", "z"]],
                          on="_key", how="left")

    for col in ("x", "y", "z"):
        merged[col] = (merged[col].astype(str).str.replace(",", ".",
                                                            regex=False)
                       .astype(float))
    merged[["x", "y", "z"]] = merged[["x", "y", "z"]] / COORD_UNIT_DIVISOR
    coords = merged[["x", "y", "z"]].to_numpy(dtype=float)
    valid = np.isfinite(coords).all(axis=1)
    return coords, valid, int((~valid).sum())


def _Z(pat: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, "imcoh_abs", IMCOH_LRG_CACHE)
    except Exception:
        return None
    return None if r is None else np.asarray(r.linkage_matrix)


def collect() -> tuple[pd.DataFrame, dict]:
    rows = []
    coords_summary: dict[str, dict] = {}
    for pat in PATIENTS_LIST:
        coords, valid, n_drop = load_coords_mm(pat)
        coords_summary[pat] = {
            "n_total": int(coords.shape[0]),
            "n_valid": int(valid.sum()),
            "n_dropped": n_drop,
            "extent_mm": (np.nanmax(coords[valid], axis=0)
                          - np.nanmin(coords[valid], axis=0)).tolist(),
        }
        if valid.sum() < 0.5 * coords.shape[0]:
            print(f"  {pat}: only {valid.sum()}/{coords.shape[0]} coords valid; skipping spatial analysis")
            continue

        for band in BRAIN_BANDS_NAMES:
            Z_pre = _Z(pat, "rest_pre", band)
            Z_test = _Z(pat, "task_test", band)
            Z_post = _Z(pat, "rest_post", band)
            if any(z is None for z in (Z_pre, Z_test, Z_post)):
                continue
            n = Z_pre.shape[0] + 1
            if Z_test.shape[0] + 1 != n or Z_post.shape[0] + 1 != n:
                continue
            if n != coords.shape[0]:
                print(f"  {pat} {band}: tree N={n} vs coords N={coords.shape[0]}; skipping")
                continue

            # Per-patient k range: [2, N-1] — go all the way to the
            # singleton-saturation regime so cluster diameters can reach 0.
            k_max_pat = min(K_MAX_GLOBAL, n - 1)
            for k in range(2, k_max_pat + 1):
                c_pre = fcluster(Z_pre, k, criterion="maxclust")
                c_test = fcluster(Z_test, k, criterion="maxclust")
                c_post = fcluster(Z_post, k, criterion="maxclust")
                d_VI = compute_vi(c_pre, c_post) - compute_vi(c_test, c_post)
                # Spatial summary computed on c_post, restricted to valid coords.
                summary_dict = cluster_spatial_scale(
                    c_post[valid], coords[valid], summary="mean_diameter",
                )
                row = {
                    "patient": pat, "band": band, "k": k,
                    "d_VI": d_VI,
                    "n_clusters_post": summary_dict["n_clusters"],
                    "n_singletons_post": summary_dict["n_singletons"],
                }
                row.update({f"{name}_post": summary_dict["all"][name]
                            for name in SPATIAL_SUMMARY_NAMES})
                rows.append(row)
    return pd.DataFrame(rows), coords_summary


# ─────────────────────────── plotting ─────────────────────────────────────


def _per_band_plot(df: pd.DataFrame, summary: str, out_pdf):
    """6-panel figure: one per band; x = spatial summary; y = Δ_VI; one curve per patient."""
    cmap = plt.get_cmap("tab10")
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.0), dpi=160,
                              sharex=False, sharey=True,
                              gridspec_kw={"hspace": 0.30, "wspace": 0.10})
    summary_col = f"{summary}_post"
    label = {
        "mean_diameter": "mean cluster diameter (mm)",
        "mean_radius": "mean cluster radius (mm)",
        "nn_centroid": "mean NN centroid distance (mm)",
        "size_weighted_radius": "size-weighted radius (mm)",
    }[summary]
    for ib, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes.flat[ib]
        sub = df[df["band"] == band]
        # Per-patient curve, sorted by spatial summary
        for ip, pat in enumerate(COHORT_N10):
            ps = sub[sub["patient"] == pat].sort_values(summary_col)
            if ps.empty:
                continue
            ax.plot(ps[summary_col], ps["d_VI"], color=cmap(ip % 10),
                    lw=1.0, alpha=0.9, marker="o", ms=3,
                    label=pat if ib == 0 else None)
        ax.axhline(0, color="0.2", lw=0.6)
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11, loc="left")
        ax.grid(alpha=0.25)
    for ax in axes[1]:
        ax.set_xlabel(label, fontsize=10)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\Delta_{\mathrm{VI}}(k)$", fontsize=10)
    axes.flat[0].legend(fontsize=7, ncol=2, loc="upper right",
                        frameon=False)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_pdf}")


def _diameter_vs_k_plot(df: pd.DataFrame, out_pdf):
    """6-panel figure: per-band, mean cluster diameter as a function of k. One curve per patient."""
    cmap = plt.get_cmap("tab10")
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.0), dpi=160,
                              sharex=True, sharey=True,
                              gridspec_kw={"hspace": 0.30, "wspace": 0.10})
    for ib, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes.flat[ib]
        sub = df[df["band"] == band]
        for ip, pat in enumerate(COHORT_N10):
            ps = sub[sub["patient"] == pat].sort_values("k")
            if ps.empty:
                continue
            ax.plot(ps["k"], ps["mean_diameter_post"],
                    color=cmap(ip % 10), lw=1.0, alpha=0.9,
                    label=pat if ib == 0 else None)
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11, loc="left")
        ax.grid(alpha=0.25)
    for ax in axes[1]:
        ax.set_xlabel("k (integer cut)", fontsize=10)
    for ax in axes[:, 0]:
        ax.set_ylabel("mean cluster diameter (mm)", fontsize=10)
    axes.flat[0].legend(fontsize=7, ncol=2, loc="upper right",
                        frameon=False)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_pdf}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("Building Δ_VI(spatial-scale) cohort dataset...")
    df, coords_summary = collect()
    csv = OUT_DIR / "dvi_spatial_scale_n10_imcoh_abs.csv"
    df.to_csv(csv, index=False)
    print(f"{len(df)} rows → {csv.relative_to(ROOT)}")

    # Coordinate alignment summary
    lines = ["# Coordinate alignment per patient", ""]
    lines.append("| patient | n_total | n_valid | n_dropped | extent_mm (Δx, Δy, Δz) |")
    lines.append("|---------|--------:|--------:|----------:|---------|")
    for pat in PATIENTS_LIST:
        s = coords_summary.get(pat, {})
        ex = s.get("extent_mm", [float("nan")] * 3)
        lines.append(
            f"| {pat} | {s.get('n_total','—')} | {s.get('n_valid','—')} | "
            f"{s.get('n_dropped','—')} | "
            f"({ex[0]:.1f}, {ex[1]:.1f}, {ex[2]:.1f}) |"
        )
    (OUT_DIR / "coords_alignment_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )

    # Plots: one figure per spatial summary
    for summary in SPATIAL_SUMMARY_NAMES:
        _per_band_plot(df, summary,
                       FIG_DIR / f"fig_dvi_vs_{summary}.pdf")
    _diameter_vs_k_plot(df, FIG_DIR / "fig_diameter_vs_k.pdf")

    md = FIG_DIR / "SUBSTANTIATION.md"
    md.write_text("""---
name: dvi-spatial-scale-figures
type: figure-index
era: COHORT_N10
created: 2026-04-27
---

# Δ_VI(spatial scale) — figures

| file | what it shows |
|---|---|
| `fig_dvi_vs_mean_diameter.pdf` | Per-band Δ_VI(k) plotted against mean cluster diameter (mm) of c_post(k). One curve per patient. |
| `fig_dvi_vs_mean_radius.pdf` | Same with mean cluster radius (mm). |
| `fig_dvi_vs_nn_centroid.pdf` | Same with mean nearest-neighbour centroid distance (mm). |
| `fig_dvi_vs_size_weighted_radius.pdf` | Same with size-weighted radius (mm). |
| `fig_diameter_vs_k.pdf` | mean cluster diameter as a function of integer k, per patient per band. The k → mm mapping. |

Reading: a vertical alignment of all 10 curves at a particular spatial
scale = cohort-wide consensus at that scale. Wide divergence = no shared
spatial scale for the trace.
""", encoding="utf-8")
    print(f"wrote {md}")


if __name__ == "__main__":
    main()
