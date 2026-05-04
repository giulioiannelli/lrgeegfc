#!/usr/bin/env python3
"""Primitive figure for the functional tree distance — focused on the
task-trace contrast.

Reads .npz cells produced by
    scripts/01_compute/diagnostics/diag_functional_tree_distance.py
under
    data/cache/functional_tree_distance/<patient>/<band>_<phaseA>-<phaseB>.npz

Writes a 2-page PDF:

(F1) **HEADLINE — task-trace contrast heatmap.** Two side-by-side panels,
     rows = patient (4), cols = band (6). Cell value = the contrast
         C_X(p, b) := Δ_X(test, post; p, b) − Δ_X(pre, post; p, b)
     for X ∈ {S, P}.  Negative ⇒ post closer to test than to pre ⇒ task-trace
     persists.  Diverging green/red colormap centered at 0.  Cell numbers
     printed.

(F2) **Three-curve δ_S(τ) panels** — 4 patients × 6 bands.  Three curves
     only: (pre,post) grey, (test,post) green, (pre,test) orange.
     `fill_between` highlights the τ-window where δ(test,post) < δ(pre,post)
     (the task-trace zone).

Diagnostic τ-interval segments are produced as a separate sanity PDF.

Output:
    data/outputs/figures/2026-04-28_functional_tree_distance_primitive.pdf
    data/outputs/figures/2026-04-28_functional_tree_distance_sanity.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (  # noqa: E402
    BRAIN_BAND_TEX_DICT,
    BRAIN_BANDS_NAMES,
    PHASE_LABELS,
)
from lrg_eegfc.config.paths import CACHE_ROOT, FIGURES_ROOT  # noqa: E402

PRIMITIVE_PATIENTS: Tuple[str, ...] = ("Pat_02", "Pat_03", "Pat_06", "Pat_13")
DEFAULT_FC_METHOD = "imcoh_abs"

IN_ROOT = CACHE_ROOT / "functional_tree_distance"
OUT_HEADLINE = FIGURES_ROOT / "2026-04-28_functional_tree_distance_primitive.pdf"
OUT_SANITY = FIGURES_ROOT / "2026-04-28_functional_tree_distance_sanity.pdf"

# Three curves we plot per panel:
PRIMARY_CURVES = (
    ("rest_pre", "rest_post"),    # baseline-drift reference (grey)
    ("task_test", "rest_post"),   # task-trace claim         (green)
    ("rest_pre", "task_test"),    # task-perturbation         (orange)
)
CURVE_COLOR = {
    ("rest_pre", "rest_post"): "#808080",
    ("task_test", "rest_post"): "#2ca02c",
    ("rest_pre", "task_test"): "#ff7f0e",
}
CURVE_LABEL = {
    ("rest_pre", "rest_post"): r"$(\mathrm{pre},\,\mathrm{post})$  (baseline drift)",
    ("task_test", "rest_post"): r"$(\mathrm{test},\,\mathrm{post})$  (task-trace)",
    ("rest_pre", "task_test"): r"$(\mathrm{pre},\,\mathrm{test})$  (task perturbation)",
}

PHASE_COLORS = {
    "rest_pre": "#1f77b4",
    "task_learn": "#2ca02c",
    "task_test": "#ff7f0e",
    "rest_post": "#d62728",
}


def load_cell(patient: str, band: str, phase_A: str, phase_B: str) -> dict | None:
    path = IN_ROOT / patient / f"{band}_{phase_A}-{phase_B}.npz"
    if not path.exists():
        return None
    data = np.load(path, allow_pickle=True)
    return {
        "patient": str(data["patient"]),
        "band": str(data["band"]),
        "phase_A": str(data["phase_A"]),
        "phase_B": str(data["phase_B"]),
        "tau_grid": data["tau_grid"],
        "delta_S": data["delta_S"],
        "delta_P": data["delta_P"],
        "delta_hat_S": float(data["delta_hat_S"]),
        "delta_hat_P": float(data["delta_hat_P"]),
        "Delta_S": float(data["Delta_S"]),
        "Delta_P": float(data["Delta_P"]),
        "tau_min": float(data["tau_min"]),
        "tau_max": float(data["tau_max"]),
        "tau_star_A": float(data["tau_star_A"]),
        "tau_star_B": float(data["tau_star_B"]),
        "lambda_max_A": float(data["lambda_max_A"]),
        "lambda_max_B": float(data["lambda_max_B"]),
        "n_nodes": int(data["n_nodes"]),
        "interval_compatible": bool(int(data["interval_compatible"])),
    }


def _contrast_grid(
    patients: List[str], bands: List[str], metric_key: str
) -> np.ndarray:
    """Per-(patient, band) `Δ(test,post) − Δ(pre,post)` for `metric_key` ∈
    {"Delta_S", "Delta_P"}.  Returns shape `(n_patients, n_bands)` with NaN
    for missing cells.
    """
    grid = np.full((len(patients), len(bands)), np.nan)
    for i, pat in enumerate(patients):
        for j, band in enumerate(bands):
            tp = load_cell(pat, band, "task_test", "rest_post")
            pp = load_cell(pat, band, "rest_pre", "rest_post")
            if tp is None or pp is None:
                continue
            grid[i, j] = tp[metric_key] - pp[metric_key]
    return grid


def render_contrast_heatmap(pdf: PdfPages, patients: List[str], bands: List[str]) -> None:
    grids: Dict[str, np.ndarray] = {
        "Delta_S": _contrast_grid(patients, bands, "Delta_S"),
        "Delta_P": _contrast_grid(patients, bands, "Delta_P"),
    }

    vmax = float(np.nanmax(np.abs(np.concatenate([g.ravel() for g in grids.values()]))))
    vmax = max(0.1, vmax)
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)
    cmap = plt.get_cmap("RdYlGn_r")

    fig, axes = plt.subplots(
        1, 2, figsize=(11.0, 3.4), dpi=160, constrained_layout=True
    )
    metric_titles = {
        "Delta_S": r"rank view  $C_S = \Delta_S(\mathrm{test,post}) - \Delta_S(\mathrm{pre,post})$",
        "Delta_P": r"height view  $C_P = \Delta_P(\mathrm{test,post}) - \Delta_P(\mathrm{pre,post})$",
    }

    for ax, key in zip(axes, ["Delta_S", "Delta_P"]):
        grid = grids[key]
        im = ax.imshow(grid, cmap=cmap, norm=norm, aspect="auto")
        for i in range(len(patients)):
            for j in range(len(bands)):
                v = grid[i, j]
                if np.isnan(v):
                    continue
                txt_color = "black" if abs(v) < 0.25 else "white"
                ax.text(
                    j,
                    i,
                    f"{v:+.2f}",
                    ha="center",
                    va="center",
                    fontsize=10,
                    color=txt_color,
                )
        ax.set_xticks(range(len(bands)))
        ax.set_xticklabels(
            [BRAIN_BAND_TEX_DICT.get(b, b) for b in bands], fontsize=11
        )
        ax.set_yticks(range(len(patients)))
        ytick_labels = []
        for p in patients:
            label = p
            ytick_labels.append(label)
        ax.set_yticklabels(ytick_labels, fontsize=10)
        for tick_label, patient in zip(ax.get_yticklabels(), patients):
            if patient == "Pat_03":
                tick_label.set_color("darkred")
        ax.set_title(metric_titles[key], fontsize=10, pad=8)
        ax.set_xlabel("band", fontsize=10)

    cbar = fig.colorbar(im, ax=axes, shrink=0.85, aspect=18, pad=0.02)
    cbar.set_label(
        "task-trace contrast\n← trace persists       no signal       anti-trace →",
        fontsize=9,
    )
    cbar.ax.tick_params(labelsize=9)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def render_three_curve_grid(
    pdf: PdfPages, patients: List[str], bands: List[str], metric: str = "delta_S"
) -> None:
    n_rows = len(patients)
    n_cols = len(bands)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(2.4 * n_cols, 1.7 * n_rows),
        dpi=150,
        sharex=False,
        sharey=True,
    )
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes[np.newaxis, :]
    elif n_cols == 1:
        axes = axes[:, np.newaxis]

    for i, patient in enumerate(patients):
        for j, band in enumerate(bands):
            ax = axes[i, j]

            cells = {}
            for pair in PRIMARY_CURVES:
                cell = load_cell(patient, band, pair[0], pair[1])
                if cell is not None and cell["interval_compatible"]:
                    cells[pair] = cell

            if not cells:
                ax.text(0.5, 0.5, "no data", ha="center", va="center",
                        fontsize=8, transform=ax.transAxes)
                ax.set_axis_off()
                continue

            ref = next(iter(cells.values()))
            log_min = np.log10(ref["tau_min"])
            log_max = np.log10(ref["tau_max"])
            ax.axvspan(log_min, log_max, color="0.96", zorder=0)

            tp = cells.get(("task_test", "rest_post"))
            pp = cells.get(("rest_pre", "rest_post"))
            if tp is not None and pp is not None:
                # Both share the same τ-grid by construction; safe to align.
                logtau = np.log10(tp["tau_grid"])
                d_tp = tp[metric]
                d_pp = pp[metric]
                trace_zone = d_tp < d_pp
                ax.fill_between(
                    logtau,
                    d_tp,
                    d_pp,
                    where=trace_zone,
                    interpolate=True,
                    facecolor="#2ca02c",
                    alpha=0.18,
                    linewidth=0,
                )
                ax.fill_between(
                    logtau,
                    d_tp,
                    d_pp,
                    where=~trace_zone,
                    interpolate=True,
                    facecolor="#d62728",
                    alpha=0.18,
                    linewidth=0,
                )

            for pair in PRIMARY_CURVES:
                cell = cells.get(pair)
                if cell is None:
                    continue
                ax.plot(
                    np.log10(cell["tau_grid"]),
                    cell[metric],
                    color=CURVE_COLOR[pair],
                    linewidth=1.6 if pair == ("task_test", "rest_post") else 1.1,
                    alpha=0.95 if pair == ("task_test", "rest_post") else 0.85,
                    zorder=3 if pair == ("task_test", "rest_post") else 2,
                )

            ax.set_ylim(0.0, 1.0)
            ax.tick_params(labelsize=7)
            if i == n_rows - 1:
                ax.set_xlabel(r"$\log_{10}\tau$", fontsize=8)
            if j == 0:
                ylabel = (
                    r"$\delta_S(\tau)$"
                    if metric == "delta_S"
                    else r"$\delta_P(\tau)$"
                )
                ax.set_ylabel(f"{patient}\n{ylabel}", fontsize=8)
            tex_band = BRAIN_BAND_TEX_DICT.get(band, band)
            if i == 0:
                ax.set_title(f"${tex_band[1:-1]}$", fontsize=10)
            if patient == "Pat_03" and j == 0:
                ax.yaxis.label.set_color("darkred")

    legend_handles = [
        plt.Line2D([0], [0], color=CURVE_COLOR[p], lw=1.8, label=CURVE_LABEL[p])
        for p in PRIMARY_CURVES
    ]
    legend_handles.append(
        plt.Line2D([0], [0], marker="s", color="#2ca02c", lw=0, markersize=10,
                   alpha=0.5, label="task-trace zone (test-post below pre-post)")
    )
    legend_handles.append(
        plt.Line2D([0], [0], marker="s", color="#d62728", lw=0, markersize=10,
                   alpha=0.5, label="anti-trace zone")
    )
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=3,
        fontsize=8.5,
        bbox_to_anchor=(0.5, -0.03),
        frameon=False,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 0.99])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def render_tau_interval_diagnostic(
    pdf: PdfPages, patients: List[str], bands: List[str]
) -> None:
    """Sanity diagnostic — kept in a separate PDF, not in the headline."""
    rows = []
    for patient in patients:
        for band in bands:
            phase_to_lmax: dict[str, float] = {}
            phase_to_taustar: dict[str, float] = {}
            for phase_A in PHASE_LABELS:
                for phase_B in PHASE_LABELS:
                    if phase_A >= phase_B:
                        continue
                    cell = load_cell(patient, band, phase_A, phase_B)
                    if cell is None:
                        continue
                    phase_to_lmax[phase_A] = cell["lambda_max_A"]
                    phase_to_lmax[phase_B] = cell["lambda_max_B"]
                    phase_to_taustar[phase_A] = cell["tau_star_A"]
                    phase_to_taustar[phase_B] = cell["tau_star_B"]
            for phase in PHASE_LABELS:
                if phase not in phase_to_lmax:
                    continue
                rows.append(
                    {
                        "patient": patient,
                        "band": band,
                        "phase": phase,
                        "tau_min": 1.0 / phase_to_lmax[phase],
                        "tau_max": phase_to_taustar[phase],
                    }
                )

    fig, axes = plt.subplots(
        1,
        len(patients),
        figsize=(2.8 * len(patients), 2.6),
        dpi=150,
        sharey=True,
    )
    if len(patients) == 1:
        axes = [axes]

    for ax, patient in zip(axes, patients):
        sub = [r for r in rows if r["patient"] == patient]
        for j, band in enumerate(bands):
            for k, phase in enumerate(PHASE_LABELS):
                row = next(
                    (r for r in sub if r["band"] == band and r["phase"] == phase),
                    None,
                )
                if row is None:
                    continue
                x = j + (k - 1.5) * 0.18
                ax.plot(
                    [x, x],
                    [np.log10(row["tau_min"]), np.log10(row["tau_max"])],
                    color=PHASE_COLORS[phase],
                    lw=2.0,
                    alpha=0.85,
                    solid_capstyle="round",
                )
        ax.set_xticks(np.arange(len(bands)))
        ax.set_xticklabels(
            [BRAIN_BAND_TEX_DICT.get(b, b) for b in bands], fontsize=8
        )
        title_color = "darkred" if patient == "Pat_03" else "black"
        ax.set_title(patient, fontsize=9, color=title_color)
        ax.set_xlabel("band", fontsize=8)
        if ax is axes[0]:
            ax.set_ylabel(r"$\log_{10}\tau$", fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(True, axis="y", alpha=0.3, lw=0.4)

    legend_handles = [
        plt.Line2D([0], [0], color=c, lw=2.2, label=p)
        for p, c in PHASE_COLORS.items()
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=4,
        fontsize=8,
        bbox_to_anchor=(0.5, -0.03),
        frameon=False,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 0.99])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--patients", nargs="+", default=list(PRIMITIVE_PATIENTS))
    parser.add_argument("--bands", nargs="+", default=list(BRAIN_BANDS_NAMES))
    parser.add_argument("--out-headline", type=Path, default=OUT_HEADLINE)
    parser.add_argument("--out-sanity", type=Path, default=OUT_SANITY)
    args = parser.parse_args()

    args.out_headline.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(args.out_headline) as pdf:
        render_contrast_heatmap(pdf, args.patients, args.bands)
        render_three_curve_grid(pdf, args.patients, args.bands, metric="delta_S")
        render_three_curve_grid(pdf, args.patients, args.bands, metric="delta_P")

    with PdfPages(args.out_sanity) as pdf:
        render_tau_interval_diagnostic(pdf, args.patients, args.bands)

    print(f"Wrote {args.out_headline}")
    print(f"Wrote {args.out_sanity}")


if __name__ == "__main__":
    main()
