#!/usr/bin/env python3
"""Audit 61b — figures for the epi Grassmann embedding (Direction E).

Reads ``data/audit/epi_grassmann/E_per_patient.csv`` (produced by
``audit_61_epi_grassmann_compute.py``) and emits three publication PDFs
per variant:

- ``fig_01_cohort_heatmap_<variant>.pdf`` — 6×3 grid of small heatmaps
  (rows = bands, cols = phases). Each heatmap rows = 9 patients + cohort
  median; cols = k ∈ K_GRID; cell = z-score against strength-stratified
  null. Diverging RdBu_r, range [-3, 3].
- ``fig_02_spaghetti_<variant>.pdf`` — 6×3 grid of line panels
  (rows = bands, cols = phases). Per panel: 9 faint per-patient lines +
  thick cohort median + shaded IQR. x = k (log), y = z. Pat_03 marked
  distinctly (1024 Hz outlier).
- ``fig_03_raw_<variant>.pdf`` — same structure as fig_02 with
  y = raw observed scalar (eps_align or delta_resect) and per-patient
  null-mean reference dashed.

PDF only, full vector (no rasterisation), no fig.suptitle, no
watermark / provenance footer (file names carry metadata).

Per the scope: NO regularity / regime annotations on the figures. We
read the cohort distributions directly.
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

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT

PHASES = ("rest_pre", "task_test", "rest_post")
PHASE_TEX = {"rest_pre": "RPre", "task_test": "TT", "rest_post": "RPost"}
K_GRID = (2, 3, 5, 8, 13, 21)
PAT03 = "Pat_03"
OUT_DIR = ROOT / "data" / "audit" / "epi_grassmann"
FIG_DIR = OUT_DIR / "figures"

# Variants: (label, observed col, z col, ylabel raw, ylabel z, dir-hint)
VARIANTS = (
    dict(
        name="align",
        obs_col="eps_align_obs",
        null_mu_col="eps_align_null_mean",
        z_col="z_align",
        raw_label=r"$\varepsilon^E_{\mathrm{align}}$",
        z_label=r"$z(\varepsilon^E_{\mathrm{align}})$",
    ),
    dict(
        name="resect",
        obs_col="delta_resect_obs",
        null_mu_col="delta_resect_null_mean",
        z_col="z_resect",
        raw_label=r"$\delta^E_{\mathrm{resect}}$",
        z_label=r"$z(\delta^E_{\mathrm{resect}})$",
    ),
)


def _patients(df: pd.DataFrame) -> list[str]:
    return sorted(df["patient"].unique())


def fig_01_heatmap(df: pd.DataFrame, variant: dict, out_path: Path) -> None:
    """Cohort z-score heatmap, 6 (bands) × 3 (phases) grid."""
    z_col = variant["z_col"]
    pats = _patients(df)
    n_rows = len(pats) + 1                  # patients + cohort median
    fig, axes = plt.subplots(
        len(BRAIN_BANDS_NAMES), len(PHASES),
        figsize=(11, 13), sharex=True, sharey=False,
    )
    for bi, band in enumerate(BRAIN_BANDS_NAMES):
        for pi, phase in enumerate(PHASES):
            ax = axes[bi, pi]
            mat = np.full((n_rows, len(K_GRID)), np.nan)
            for ri, pat in enumerate(pats):
                for ki, k in enumerate(K_GRID):
                    cell = df[
                        (df.patient == pat)
                        & (df.band == band)
                        & (df.phase == phase)
                        & (df.k == k)
                    ]
                    if len(cell) == 1:
                        mat[ri, ki] = cell[z_col].iloc[0]
            # cohort-median row at the bottom
            for ki, k in enumerate(K_GRID):
                col = mat[:-1, ki]
                col = col[~np.isnan(col)]
                mat[-1, ki] = float(np.median(col)) if col.size else np.nan
            im = ax.imshow(
                mat, cmap="RdBu_r", vmin=-3, vmax=3,
                aspect="auto", interpolation="nearest",
            )
            if bi == 0:
                ax.set_title(PHASE_TEX[phase], fontsize=10)
            if pi == 0:
                ax.set_ylabel(BRAIN_BAND_TEX_DICT[band], fontsize=10)
                yt = list(range(n_rows))
                ax.set_yticks(yt)
                ax.set_yticklabels(
                    pats + ["median"], fontsize=6,
                )
            else:
                ax.set_yticks([])
            if bi == len(BRAIN_BANDS_NAMES) - 1:
                ax.set_xticks(range(len(K_GRID)))
                ax.set_xticklabels([str(k) for k in K_GRID], fontsize=8)
                ax.set_xlabel(r"$k$", fontsize=9)
            else:
                ax.set_xticks([])
            # divider line between cohort-median row and patients
            ax.axhline(n_rows - 1.5, color="k", lw=0.6)
    fig.subplots_adjust(left=0.10, right=0.92, top=0.96, bottom=0.06,
                        wspace=0.05, hspace=0.10)
    cax = fig.add_axes([0.94, 0.20, 0.012, 0.60])
    fig.colorbar(im, cax=cax, label=variant["z_label"])
    fig.savefig(out_path)
    plt.close(fig)


def _bounded_scalar(df_row_group: pd.DataFrame, variant_name: str) -> np.ndarray:
    """Return a per-row bounded scalar that's robust to null-variance collapse.

    For ``align``: ``f_align_obs − f_uniform`` where
    ``f_uniform = |E_p| / N`` (per-patient chance level). Signed; 0 = random;
    > 0 = modes lean toward epi; < 0 = away.

    For ``resect``: ``(delta_resect_obs − delta_resect_null_mean) / sqrt(k)``.
    Signed; 0 = epi-resection ≈ random-resection; > 0 = epi-resection
    rotates the eigenspace MORE than random.

    Both stay finite even when the null std collapses (low-k delocalised
    modes); both are directly interpretable on a unit-ish scale.
    """
    if variant_name == "align":
        f_uni = df_row_group["n_epi"] / df_row_group["n_nodes"]
        return (df_row_group["f_align_obs"] - f_uni).values
    elif variant_name == "resect":
        denom = np.sqrt(df_row_group["k"].astype(float))
        return ((df_row_group["delta_resect_obs"]
                 - df_row_group["delta_resect_null_mean"]) / denom).values
    else:
        raise ValueError(variant_name)


def fig_02_spaghetti(df: pd.DataFrame, variant: dict, out_path: Path) -> None:
    """Per-patient bounded-scalar spaghetti, 6×3 grid.

    Uses ``f_align − f_uniform`` (align) or ``(δ − null_mean) / sqrt(k)``
    (resect) — both bounded, signed, and finite even when the null
    variance collapses (low-k delocalised modes make z-scores blow up;
    these scalars don't).
    """
    pats = _patients(df)
    pat_labels = {
        "align": r"$f^E_{\mathrm{align}} - |E_p|/N$",
        "resect": r"$(\delta^E_{\mathrm{resect}} - \mu_{\mathrm{null}}) / \sqrt{k}$",
    }[variant["name"]]
    fig, axes = plt.subplots(
        len(BRAIN_BANDS_NAMES), len(PHASES),
        figsize=(11, 13), sharex=True, sharey="row",
    )
    for bi, band in enumerate(BRAIN_BANDS_NAMES):
        for pi, phase in enumerate(PHASES):
            ax = axes[bi, pi]
            sub = df[(df.band == band) & (df.phase == phase)]
            cohort_at_k = {k: [] for k in K_GRID}
            for pat in pats:
                pat_rows = sub[sub.patient == pat].sort_values("k")
                if pat_rows.empty:
                    continue
                vals_full = _bounded_scalar(pat_rows, variant["name"])
                ks = pat_rows["k"].values
                vals = []
                for k in K_GRID:
                    matches = np.where(ks == k)[0]
                    v = float(vals_full[matches[0]]) if matches.size else np.nan
                    vals.append(v)
                    if not np.isnan(v):
                        cohort_at_k[k].append(v)
                style = dict(color="tab:orange", lw=0.9, ls=":") \
                    if pat == PAT03 else dict(color="0.55", lw=0.6, alpha=0.8)
                ax.plot(K_GRID, vals, **style)
            med = [np.median(cohort_at_k[k]) if cohort_at_k[k] else np.nan
                   for k in K_GRID]
            q1 = [np.quantile(cohort_at_k[k], 0.25) if cohort_at_k[k] else np.nan
                  for k in K_GRID]
            q3 = [np.quantile(cohort_at_k[k], 0.75) if cohort_at_k[k] else np.nan
                  for k in K_GRID]
            ax.fill_between(K_GRID, q1, q3, color="tab:blue", alpha=0.18, lw=0)
            ax.plot(K_GRID, med, color="tab:blue", lw=1.6)
            ax.axhline(0, color="0.3", lw=0.5, ls="--")
            ax.set_xscale("log")
            ax.set_xticks(K_GRID)
            ax.set_xticklabels([str(k) for k in K_GRID], fontsize=7)
            if bi == 0:
                ax.set_title(PHASE_TEX[phase], fontsize=10)
            if pi == 0:
                ax.set_ylabel(BRAIN_BAND_TEX_DICT[band] + "\n" + pat_labels,
                              fontsize=8)
            if bi == len(BRAIN_BANDS_NAMES) - 1:
                ax.set_xlabel(r"$k$", fontsize=9)
    fig.subplots_adjust(left=0.12, right=0.97, top=0.96, bottom=0.06,
                        wspace=0.10, hspace=0.18)
    fig.savefig(out_path)
    plt.close(fig)


def fig_04_per_patient_heatmap(df: pd.DataFrame, variant: dict, out_path: Path) -> None:
    """Per-patient × k heatmap of the bounded scalar, 6×3 (band×phase) grid.

    Rows of each small heatmap = patients (sorted alphabetically). Columns
    = k. Cell = the bounded scalar (`f_align − f_uniform` for align;
    `(δ − null_mean)/sqrt(k)` for resect). Diverging colormap centred at 0.
    Lets the reader see *which patients* drive any apparent cohort
    structure — distinguishes "single-patient extreme" from "broad cohort".
    """
    pats = _patients(df)
    n_pat = len(pats)
    fig, axes = plt.subplots(
        len(BRAIN_BANDS_NAMES), len(PHASES),
        figsize=(11, 13), sharex=True, sharey=False,
    )
    # auto-scale per variant: align in [-0.3, +0.8], resect in [-0.5, +0.5]
    if variant["name"] == "align":
        vmin, vmax = -0.4, 0.4
    else:
        vmin, vmax = -0.5, 0.5
    for bi, band in enumerate(BRAIN_BANDS_NAMES):
        for pi, phase in enumerate(PHASES):
            ax = axes[bi, pi]
            mat = np.full((n_pat, len(K_GRID)), np.nan)
            for ri, pat in enumerate(pats):
                sub = df[
                    (df.band == band)
                    & (df.phase == phase)
                    & (df.patient == pat)
                ].sort_values("k")
                if sub.empty:
                    continue
                vals = _bounded_scalar(sub, variant["name"])
                ks = sub["k"].values
                for ki, k in enumerate(K_GRID):
                    matches = np.where(ks == k)[0]
                    if matches.size:
                        mat[ri, ki] = float(vals[matches[0]])
            im = ax.imshow(
                mat, cmap="RdBu_r", vmin=vmin, vmax=vmax,
                aspect="auto", interpolation="nearest",
            )
            if bi == 0:
                ax.set_title(PHASE_TEX[phase], fontsize=10)
            if pi == 0:
                ax.set_ylabel(BRAIN_BAND_TEX_DICT[band], fontsize=10)
                ax.set_yticks(range(n_pat))
                ax.set_yticklabels(pats, fontsize=6)
            else:
                ax.set_yticks([])
            if bi == len(BRAIN_BANDS_NAMES) - 1:
                ax.set_xticks(range(len(K_GRID)))
                ax.set_xticklabels([str(k) for k in K_GRID], fontsize=8)
                ax.set_xlabel(r"$k$", fontsize=9)
            else:
                ax.set_xticks([])
    fig.subplots_adjust(left=0.10, right=0.92, top=0.96, bottom=0.06,
                        wspace=0.05, hspace=0.10)
    cax = fig.add_axes([0.94, 0.20, 0.012, 0.60])
    label = (r"$f^E_{\mathrm{align}} - |E_p|/N$"
             if variant["name"] == "align"
             else r"$(\delta^E_{\mathrm{resect}} - \mu_{\mathrm{null}})/\sqrt{k}$")
    fig.colorbar(im, cax=cax, label=label)
    fig.savefig(out_path)
    plt.close(fig)


def fig_03_raw(df: pd.DataFrame, variant: dict, out_path: Path) -> None:
    """Raw observed scalar with per-patient null-mean reference, 6×3 grid."""
    obs_col = variant["obs_col"]
    null_mu_col = variant["null_mu_col"]
    pats = _patients(df)
    fig, axes = plt.subplots(
        len(BRAIN_BANDS_NAMES), len(PHASES),
        figsize=(11, 13), sharex=True, sharey="row",
    )
    for bi, band in enumerate(BRAIN_BANDS_NAMES):
        for pi, phase in enumerate(PHASES):
            ax = axes[bi, pi]
            sub = df[(df.band == band) & (df.phase == phase)]
            cohort_obs = {k: [] for k in K_GRID}
            cohort_mu = {k: [] for k in K_GRID}
            for pat in pats:
                obs = []
                mus = []
                for k in K_GRID:
                    cell = sub[(sub.patient == pat) & (sub.k == k)]
                    v_obs = cell[obs_col].iloc[0] if len(cell) == 1 else np.nan
                    v_mu = cell[null_mu_col].iloc[0] if len(cell) == 1 else np.nan
                    obs.append(v_obs)
                    mus.append(v_mu)
                    if not np.isnan(v_obs):
                        cohort_obs[k].append(v_obs)
                    if not np.isnan(v_mu):
                        cohort_mu[k].append(v_mu)
                style_obs = dict(color="tab:orange", lw=0.9, ls=":") \
                    if pat == PAT03 else dict(color="0.55", lw=0.6, alpha=0.8)
                ax.plot(K_GRID, obs, **style_obs)
            med_obs = [np.median(cohort_obs[k]) if cohort_obs[k] else np.nan
                       for k in K_GRID]
            q1 = [np.quantile(cohort_obs[k], 0.25) if cohort_obs[k] else np.nan
                  for k in K_GRID]
            q3 = [np.quantile(cohort_obs[k], 0.75) if cohort_obs[k] else np.nan
                  for k in K_GRID]
            med_mu = [np.median(cohort_mu[k]) if cohort_mu[k] else np.nan
                      for k in K_GRID]
            ax.fill_between(K_GRID, q1, q3, color="tab:blue", alpha=0.18, lw=0)
            ax.plot(K_GRID, med_obs, color="tab:blue", lw=1.6, label="observed")
            ax.plot(K_GRID, med_mu, color="tab:red", lw=1.0, ls="--",
                    label="null mean")
            ax.set_xscale("log")
            ax.set_xticks(K_GRID)
            ax.set_xticklabels([str(k) for k in K_GRID], fontsize=7)
            if bi == 0:
                ax.set_title(PHASE_TEX[phase], fontsize=10)
            if pi == 0:
                ax.set_ylabel(BRAIN_BAND_TEX_DICT[band] + "\n" + variant["raw_label"],
                              fontsize=8)
            if bi == len(BRAIN_BANDS_NAMES) - 1:
                ax.set_xlabel(r"$k$", fontsize=9)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center",
               bbox_to_anchor=(0.5, 0.005), ncol=2, frameon=False, fontsize=9)
    fig.subplots_adjust(left=0.10, right=0.97, top=0.96, bottom=0.07,
                        wspace=0.10, hspace=0.18)
    fig.savefig(out_path)
    plt.close(fig)


def main() -> None:
    csv_path = OUT_DIR / "E_per_patient.csv"
    if not csv_path.exists():
        raise SystemExit(f"missing {csv_path}; run audit_61 first")
    df = pd.read_csv(csv_path)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for variant in VARIANTS:
        v = variant["name"]
        fig_01_heatmap(df, variant, FIG_DIR / f"fig_01_cohort_heatmap_{v}.pdf")
        fig_02_spaghetti(df, variant, FIG_DIR / f"fig_02_spaghetti_{v}.pdf")
        fig_03_raw(df, variant, FIG_DIR / f"fig_03_raw_{v}.pdf")
        fig_04_per_patient_heatmap(
            df, variant, FIG_DIR / f"fig_04_per_patient_heatmap_{v}.pdf"
        )
        print(f"[audit_61b] wrote {v} figures to {FIG_DIR}")


if __name__ == "__main__":
    main()
