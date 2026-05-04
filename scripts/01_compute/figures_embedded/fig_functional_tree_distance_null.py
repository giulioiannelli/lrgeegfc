#!/usr/bin/env python3
"""Cohort obs-vs-noise figure for the functional tree distance.

Reads the .npz cells produced by
    scripts/01_compute/diagnostics/diag_functional_tree_distance_null.py
under
    data/cache/functional_tree_distance_null/<patient>/<band>_null.npz
plus the observed Δ values from
    data/cache/functional_tree_distance/<patient>/<band>_<phaseA>-<phaseB>.npz

Produces a multi-page PDF:

(F1) **Contrast vs noise scatter — both metrics.** Two side-by-side panels
     (Δ_S left, Δ_P right).  X = within-phase median (noise floor);
     Y = obs contrast `C = Δ(test,post) − Δ(pre,post)`.  Negative Y =
     task-trace; amber wedge = noise envelope `|C| ≤ within-phase median`.
     One dot per (patient, band), colour-coded by patient.

(F2…F7) **One page per band.** Two stacked panels (Δ_S top, Δ_P bottom).
     X = patient (cohort N).  Per patient column: amber IQR rectangle +
     dashed median line (within-phase noise floor); grey × marks for the
     4 within-phase Δ values; blue dot for obs Δ(pre,post); green dot
     for obs Δ(test,post); connecting line between them.  Reading: green
     below blue AND outside amber → task-trace persists above noise.

Output:
    data/outputs/figures/2026-04-28_functional_tree_distance_null.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (  # noqa: E402
    BRAIN_BAND_TEX_DICT,
    BRAIN_BANDS_NAMES,
)
from lrg_eegfc.config.paths import CACHE_ROOT, FIGURES_ROOT  # noqa: E402

DEFAULT_PATIENTS: Tuple[str, ...] = (
    "Pat_02",
    "Pat_03",
    "Pat_05",
    "Pat_06",
    "Pat_07",
    "Pat_08",
    "Pat_10",
    "Pat_13",
    "Pat_14",
    "Pat_15",
)
DEFAULT_FC_METHOD = "imcoh_abs"

OBS_ROOT = CACHE_ROOT / "functional_tree_distance"
NULL_ROOT = CACHE_ROOT / "functional_tree_distance_null"
OUT_PATH = FIGURES_ROOT / "2026-04-28_functional_tree_distance_null.pdf"

# 10-patient categorical palette (Tableau 10 + extras for clarity)
PATIENT_COLOR = {
    "Pat_02": "#1f77b4",
    "Pat_03": "#d62728",
    "Pat_05": "#9467bd",
    "Pat_06": "#2ca02c",
    "Pat_07": "#ff7f0e",
    "Pat_08": "#17becf",
    "Pat_10": "#bcbd22",
    "Pat_13": "#8c564b",
    "Pat_14": "#e377c2",
    "Pat_15": "#7f7f7f",
}


def load_null_cell(patient: str, band: str) -> Optional[dict]:
    path = NULL_ROOT / patient / f"{band}_null.npz"
    if not path.exists():
        return None
    data = np.load(path, allow_pickle=True)
    return {
        "patient": str(data["patient"]),
        "band": str(data["band"]),
        "n_nodes": int(data["n_nodes"]),
        "within_phase_labels": [str(x) for x in data["within_phase_labels"]],
        "within_phase_Delta_S": np.asarray(data["within_phase_Delta_S"]),
        "within_phase_Delta_P": np.asarray(data["within_phase_Delta_P"]),
    }


def load_obs_cell(
    patient: str, band: str, phase_A: str, phase_B: str
) -> Optional[dict]:
    path = OBS_ROOT / patient / f"{band}_{phase_A}-{phase_B}.npz"
    if not path.exists():
        return None
    data = np.load(path, allow_pickle=True)
    return {
        "Delta_S": float(data["Delta_S"]),
        "Delta_P": float(data["Delta_P"]),
    }


def render_contrast_vs_noise(
    pdf: PdfPages, patients: List[str], bands: List[str]
) -> None:
    rows = []
    for patient in patients:
        for band in bands:
            null_cell = load_null_cell(patient, band)
            if null_cell is None:
                continue
            obs_pp = load_obs_cell(patient, band, "rest_pre", "rest_post")
            obs_tp = load_obs_cell(patient, band, "task_test", "rest_post")
            if obs_pp is None or obs_tp is None:
                continue
            rows.append(
                {
                    "patient": patient,
                    "band": band,
                    "within_S_med": float(np.median(null_cell["within_phase_Delta_S"])),
                    "within_P_med": float(np.median(null_cell["within_phase_Delta_P"])),
                    "C_S": obs_tp["Delta_S"] - obs_pp["Delta_S"],
                    "C_P": obs_tp["Delta_P"] - obs_pp["Delta_P"],
                }
            )
    if not rows:
        return
    df = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), dpi=160)
    metrics = [
        ("S", r"rank view  $C_S$ vs noise floor"),
        ("P", r"height view  $C_P$ vs noise floor"),
    ]
    for ax, (key, title) in zip(axes, metrics):
        within_col = f"within_{key}_med"
        contrast_col = f"C_{key}"
        x_max = max(df[within_col].max(), 0.05) * 1.1
        xs = np.array([0, x_max])
        ax.fill_between(
            xs, np.zeros_like(xs), xs, color="#ffd6a0", alpha=0.85, zorder=0,
            label=r"within-phase noise envelope ($|C| \leq $ noise)",
        )
        ax.fill_between(
            xs, -xs, np.zeros_like(xs), color="#ffd6a0", alpha=0.85, zorder=0,
        )
        ax.plot(xs, xs, color="#b07a00", lw=1.2, ls="--", zorder=1)
        ax.plot(xs, -xs, color="#b07a00", lw=1.2, ls="--", zorder=1)
        ax.axhline(0, color="black", lw=0.9, alpha=0.7, zorder=1)
        for patient in patients:
            sub = df[df["patient"] == patient]
            if sub.empty:
                continue
            ax.scatter(
                sub[within_col],
                sub[contrast_col],
                s=70,
                facecolor=PATIENT_COLOR.get(patient, "#444444"),
                edgecolor="black",
                linewidth=0.5,
                alpha=0.92,
                label=patient,
                zorder=3,
            )
            for _, r in sub.iterrows():
                ax.annotate(
                    BRAIN_BAND_TEX_DICT.get(r["band"], r["band"]),
                    (r[within_col], r[contrast_col]),
                    fontsize=6,
                    xytext=(4, 2),
                    textcoords="offset points",
                    color=PATIENT_COLOR.get(patient, "#444444"),
                )
        ax.set_xlabel("within-phase noise floor (median Δ across 4 phases)", fontsize=9)
        ax.set_ylabel(
            r"$C = \Delta(\mathrm{test, post}) - \Delta(\mathrm{pre, post})$  "
            r"(negative = task-trace persists)",
            fontsize=9,
        )
        ax.set_title(title, fontsize=10)
        ax.set_xlim(0, x_max)
        ax.tick_params(labelsize=8)
        ax.legend(loc="upper right", fontsize=6.5, frameon=False, ncol=2)
        ax.grid(True, alpha=0.25, lw=0.4)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _band_panel(
    ax,
    patients: List[str],
    band: str,
    metric_key: str,
    metric_label: str,
) -> None:
    """Draw one panel: 10 patients on x-axis, Δ values on y-axis."""
    rng = np.random.default_rng(hash(band) & 0xFFFF)

    for j, patient in enumerate(patients):
        null_cell = load_null_cell(patient, band)
        if null_cell is None:
            continue

        within_key = (
            "within_phase_Delta_S"
            if metric_key == "Delta_S"
            else "within_phase_Delta_P"
        )
        within_vals = null_cell[within_key]

        within_q1 = float(np.quantile(within_vals, 0.25))
        within_q3 = float(np.quantile(within_vals, 0.75))
        within_med = float(np.median(within_vals))

        # Per-patient amber IQR rectangle, width 0.7
        ax.add_patch(
            Rectangle(
                (j - 0.35, within_q1),
                0.70,
                within_q3 - within_q1,
                facecolor="#ffd680",
                alpha=0.55,
                edgecolor="none",
                zorder=0,
            )
        )
        ax.plot(
            [j - 0.35, j + 0.35],
            [within_med, within_med],
            color="#b07a00",
            lw=1.0,
            ls="--",
            alpha=0.85,
            zorder=1,
        )

        # Within-phase Δ values as small grey crosses, jittered
        x_within = j + rng.uniform(-0.20, 0.20, size=len(within_vals))
        ax.scatter(
            x_within,
            within_vals,
            marker="x",
            s=28,
            c="#666666",
            linewidths=1.0,
            alpha=0.85,
            zorder=2,
        )

        obs_pp = load_obs_cell(patient, band, "rest_pre", "rest_post")
        obs_tp = load_obs_cell(patient, band, "task_test", "rest_post")
        if obs_pp is None and obs_tp is None:
            continue
        if obs_pp is not None and obs_tp is not None:
            ax.plot(
                [j - 0.07, j + 0.07],
                [obs_pp[metric_key], obs_tp[metric_key]],
                color="#444444",
                lw=0.9,
                alpha=0.7,
                zorder=3,
            )
        if obs_pp is not None:
            ax.scatter(
                [j - 0.07],
                [obs_pp[metric_key]],
                marker="o",
                s=70,
                facecolor="#1f77b4",
                edgecolor="black",
                linewidth=0.6,
                zorder=4,
            )
        if obs_tp is not None:
            ax.scatter(
                [j + 0.07],
                [obs_tp[metric_key]],
                marker="o",
                s=70,
                facecolor="#2ca02c",
                edgecolor="black",
                linewidth=0.6,
                zorder=4,
            )

    ax.set_xticks(range(len(patients)))
    xtick_labels = [p.replace("Pat_", "") for p in patients]
    ax.set_xticklabels(xtick_labels, fontsize=8)
    ax.set_xlim(-0.6, len(patients) - 0.4)
    ax.set_ylim(0.0, 1.0)
    ax.tick_params(axis="y", labelsize=8)
    ax.set_ylabel(metric_label, fontsize=10)
    ax.grid(True, axis="y", alpha=0.2, lw=0.4)
    for tick_label, patient in zip(ax.get_xticklabels(), patients):
        if patient == "Pat_03":
            tick_label.set_color("darkred")


def render_band_page(
    pdf: PdfPages, patients: List[str], band: str
) -> None:
    fig, axes = plt.subplots(
        2, 1, figsize=(11.0, 6.5), dpi=150, sharex=True
    )

    tex_band = BRAIN_BAND_TEX_DICT.get(band, band)
    _band_panel(
        axes[0],
        patients,
        band,
        "Delta_S",
        rf"rank view  $\Delta_S$  ({tex_band} band)",
    )
    _band_panel(
        axes[1],
        patients,
        band,
        "Delta_P",
        rf"height view  $\Delta_P$  ({tex_band} band)",
    )
    axes[1].set_xlabel("patient (Pat_NN)", fontsize=9)

    legend_handles = [
        plt.Line2D(
            [0], [0], marker="x", color="#666666", lw=0, markersize=8,
            label=r"within-phase $\Delta$ (noise; 4 phases)",
        ),
        plt.Rectangle(
            (0, 0), 1, 1, facecolor="#ffd680", alpha=0.7,
            label="within-phase IQR (per patient)",
        ),
        plt.Line2D(
            [0], [0], color="#b07a00", lw=1.2, ls="--",
            label="within-phase median (per patient)",
        ),
        plt.Line2D(
            [0], [0], marker="o", color="#1f77b4", lw=0, markersize=8,
            markeredgecolor="black", markeredgewidth=0.5,
            label=r"obs $\Delta(\mathrm{pre},\mathrm{post})$",
        ),
        plt.Line2D(
            [0], [0], marker="o", color="#2ca02c", lw=0, markersize=8,
            markeredgecolor="black", markeredgewidth=0.5,
            label=r"obs $\Delta(\mathrm{test},\mathrm{post})$",
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=5,
        fontsize=8,
        bbox_to_anchor=(0.5, -0.02),
        frameon=False,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.99])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--patients", nargs="+", default=list(DEFAULT_PATIENTS))
    parser.add_argument("--bands", nargs="+", default=list(BRAIN_BANDS_NAMES))
    parser.add_argument("--out", type=Path, default=OUT_PATH)
    args = parser.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(args.out) as pdf:
        render_contrast_vs_noise(pdf, args.patients, args.bands)
        for band in args.bands:
            render_band_page(pdf, args.patients, band)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
