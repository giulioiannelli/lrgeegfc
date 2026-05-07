#!/usr/bin/env python3
"""Audit 46 — Grassmann principal-angle decomposition for §5.4.

Un-collapses the chordal distance d_G = sqrt(k - sum sigma_i^2) into its
constituent principal angles theta_i = arccos(sigma_i). For each
(patient, band, k) cell the script saves:

  - theta_pre_tt[i]  = principal angles of (V^RPre_k)^T V^TT_k
  - theta_tt_post[i] = principal angles of (V^TT_k)^T V^RPost_k
  - delta_theta[i]   = theta_pre_tt[i] - theta_tt_post[i]
                       (positive = trace direction at mode i: post-task
                        subspace is closer to task than pre-task is)

A sanity check verifies sqrt(sum sin^2(theta_i)) reproduces the chordal
distance reported in `04_grassmann_triangle/cohort_summary_k.csv` /
`section5_v2_round2/tables/grassmann_k_sweep.csv`.

Cohort aggregation: median delta_theta and n_trace per (band, k, mode_i).

Outputs
-------
``data/audit/section5_v2_round3_redo/tables/
    grassmann_principal_angles_per_patient.csv``
``data/audit/section5_v2_round3_redo/tables/
    grassmann_principal_angles_cohort.csv``
``data/audit/section5_v2_round3_redo/figures/
    grassmann_principal_angles.pdf``
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PATIENTS_LIST,
)
from lrg_eegfc.workflow.lrg import load_lrg_result

K_GRID = list(range(2, 81))  # continuous k = 2..80, matches VI(k) clean range
PHASES = ("rest_pre", "task_test", "rest_post")
OUT = ROOT / "data" / "audit" / "section5_v2_round3_redo"
TBL = OUT / "tables"
FIG = OUT / "figures"

TRACE_BANDS = {"alpha", "beta", "low_gamma"}
NULL_BANDS = {"delta", "theta", "high_gamma"}
BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]


def topk_basis(eigvecs: np.ndarray, k: int) -> np.ndarray:
    return np.ascontiguousarray(eigvecs[:, 1 : k + 1])


def principal_angles(V_a: np.ndarray, V_b: np.ndarray) -> np.ndarray:
    sigma = np.linalg.svd(V_a.T @ V_b, compute_uv=False)
    sigma = np.clip(sigma, 0.0, 1.0)
    return np.arccos(sigma)


def chordal_from_angles(theta: np.ndarray) -> float:
    return float(np.sqrt(np.sum(np.sin(theta) ** 2)))


def compute_per_patient() -> pd.DataFrame:
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            try:
                EV = {}
                for ph in PHASES:
                    res = load_lrg_result(pat, ph, band, fc_method="imcoh_abs")
                    if res is None or res.eigenvectors is None:
                        raise FileNotFoundError(
                            f"missing eigenvectors {pat} {ph} {band}"
                        )
                    EV[ph] = res.eigenvectors
            except Exception as e:
                print(f"[audit_46] WARN load {pat} {band}: {e}")
                continue

            for k in K_GRID:
                try:
                    Vpre = topk_basis(EV["rest_pre"], k)
                    Vtt = topk_basis(EV["task_test"], k)
                    Vpost = topk_basis(EV["rest_post"], k)
                    th_pre_tt = principal_angles(Vpre, Vtt)
                    th_tt_post = principal_angles(Vtt, Vpost)
                    d_pre_tt = chordal_from_angles(th_pre_tt)
                    d_tt_post = chordal_from_angles(th_tt_post)
                    for i, (a, b) in enumerate(zip(th_pre_tt, th_tt_post), start=1):
                        rows.append({
                            "patient": pat,
                            "band": band,
                            "k": k,
                            "mode_index": i,
                            "theta_pre_tt": float(a),
                            "theta_tt_post": float(b),
                            "delta_theta": float(a - b),
                            "d_pre_tt": d_pre_tt,
                            "d_tt_post": d_tt_post,
                            "T_E1": d_tt_post - d_pre_tt,
                        })
                except Exception as e:
                    print(f"[audit_46] WARN compute {pat} {band} k={k}: {e}")
    return pd.DataFrame(rows)


def aggregate_cohort(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for band in BAND_ORDER:
        for k in K_GRID:
            for i in range(1, k + 1):
                sub = df[
                    (df["band"] == band)
                    & (df["k"] == k)
                    & (df["mode_index"] == i)
                ]
                if not len(sub):
                    continue
                rows.append({
                    "band": band,
                    "k": k,
                    "mode_index": i,
                    "delta_theta_median_radians": float(sub["delta_theta"].median()),
                    "delta_theta_iqr_radians": float(
                        np.subtract(*np.percentile(sub["delta_theta"], [75, 25]))
                    ),
                    "n_patients": int(len(sub)),
                })
            sub_k = df[(df["band"] == band) & (df["k"] == k) & (df["mode_index"] == 1)]
            if not len(sub_k):
                continue
            T_vals = sub_k.groupby("patient")["T_E1"].first().values
    return pd.DataFrame(rows)


def cohort_n_trace(df: pd.DataFrame) -> dict:
    out = {}
    for band in BAND_ORDER:
        for k in K_GRID:
            sub = (
                df[(df["band"] == band) & (df["k"] == k)]
                .drop_duplicates(subset=["patient"])
            )
            if not len(sub):
                continue
            out[(band, k)] = (
                int((sub["T_E1"] < 0).sum()),
                int(len(sub)),
            )
    return out


def sanity_check_chordal(df: pd.DataFrame) -> None:
    """Verify chordal scalar from extracted angles matches existing tables."""
    canonical = pd.read_csv(
        ROOT
        / "data"
        / "reports"
        / "section_5_lrg_trace"
        / "04_grassmann_triangle"
        / "tables"
        / "Td_per_patient_per_band_k.csv"
    )
    sweep = pd.read_csv(
        ROOT
        / "data"
        / "audit"
        / "section5_v2_round2"
        / "tables"
        / "grassmann_k_sweep.csv"
    )
    print("[audit_46] sanity check: chordal scalar reproduction")
    sample = df.drop_duplicates(subset=["patient", "band", "k"])[
        ["patient", "band", "k", "d_pre_tt", "d_tt_post", "T_E1"]
    ]
    for k_test in [3, 13]:
        ours = sample[sample["k"] == k_test].head(3)
        for _, r in ours.iterrows():
            cmp = canonical[
                (canonical["patient"] == r["patient"])
                & (canonical["band"] == r["band"])
                & (canonical["k"] == k_test)
            ]
            if not len(cmp):
                continue
            ref_T = float(cmp["T_E1"].iloc[0])
            err = abs(ref_T - r["T_E1"])
            print(
                f"  {r['patient']:7s} {r['band']:11s} k={k_test:2d}  "
                f"T_ref={ref_T:+.6f}  T_ours={r['T_E1']:+.6f}  |Δ|={err:.2e}"
            )
    print()


def render_figure(cohort: pd.DataFrame, n_trace_lookup: dict) -> None:
    K_min, K_max = K_GRID[0], K_GRID[-1]

    abs_vals = cohort["delta_theta_median_radians"].abs().values
    vmax = float(np.percentile(abs_vals, 95)) if len(abs_vals) else 0.5
    if vmax < 0.05:
        vmax = 0.05

    fig = plt.figure(figsize=(13.2, 10.0))
    outer = GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.34)
    panel_order = [
        ("delta", 0, 0), ("theta", 0, 1), ("alpha", 0, 2),
        ("beta", 1, 0), ("low_gamma", 1, 1), ("high_gamma", 1, 2),
    ]

    for band, r, c in panel_order:
        inner = outer[r, c].subgridspec(2, 1, height_ratios=[3.5, 1.0], hspace=0.18)
        ax_h = fig.add_subplot(inner[0, 0])
        ax_s = fig.add_subplot(inner[1, 0], sharex=ax_h)

        H = np.full((K_max, len(K_GRID)), np.nan)
        for j, k in enumerate(K_GRID):
            sub = cohort[(cohort["band"] == band) & (cohort["k"] == k)]
            for _, row in sub.iterrows():
                i = int(row["mode_index"]) - 1
                if i < K_max:
                    H[i, j] = row["delta_theta_median_radians"]

        if band in TRACE_BANDS:
            ax_h.set_facecolor("#fff7e6")
            ax_s.set_facecolor("#fff7e6")

        im = ax_h.imshow(
            H,
            aspect="auto",
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
            origin="upper",
            interpolation="nearest",
            extent=(K_min - 0.5, K_max + 0.5, K_max + 0.5, 0.5),
        )

        ax_h.plot(
            [K_min - 0.5, K_max + 0.5],
            [K_min - 0.5, K_max + 0.5],
            color="0.55",
            lw=0.5,
            ls="--",
            zorder=4,
        )

        x_ticks = [2, 10, 20, 30, 40, 50, 60, 70, 80]
        ax_h.set_xticks(x_ticks)
        ax_h.set_xticklabels([str(k) for k in x_ticks], fontsize=8)
        y_ticks = [1, 10, 20, 30, 40, 50, 60, 70, 80]
        ax_h.set_yticks(y_ticks)
        ax_h.set_yticklabels([str(t) for t in y_ticks], fontsize=8)
        ax_h.set_ylabel("mode index $i$", fontsize=8)
        title_band = BRAIN_BAND_TEX_DICT[band]
        suffix = "  (trace)" if band in TRACE_BANDS else ""
        ax_h.set_title(rf"{title_band}{suffix}", fontsize=11, pad=4)
        for spine in ax_h.spines.values():
            spine.set_linewidth(0.6)
        plt.setp(ax_h.get_xticklabels(), visible=False)

        n_trace_vals = [n_trace_lookup.get((band, k), (0, 10))[0] for k in K_GRID]
        line_color = "#1a7c3e" if band in TRACE_BANDS else "#333333"
        ax_s.plot(
            K_GRID,
            n_trace_vals,
            color=line_color,
            lw=1.2,
        )
        ax_s.axhline(7, color="0.5", lw=0.7, ls="--", zorder=1)
        ax_s.set_xlim(K_min - 0.5, K_max + 0.5)
        ax_s.set_ylim(-0.4, 10.4)
        ax_s.set_yticks([0, 5, 10])
        ax_s.set_yticklabels(["0", "5", "10"], fontsize=8)
        ax_s.set_xticks(x_ticks)
        ax_s.set_xticklabels([str(k) for k in x_ticks], fontsize=8)
        ax_s.set_xlabel(r"spectral cutoff $k$", fontsize=8)
        ax_s.set_ylabel(r"$n_{\rm trace}/10$", fontsize=8)
        for spine in ax_s.spines.values():
            spine.set_linewidth(0.6)

    cax = fig.add_axes([0.92, 0.20, 0.014, 0.60])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(
        r"cohort-median $\Delta\theta_i$ (radians) — red = trace direction",
        fontsize=9,
    )
    cb.ax.tick_params(labelsize=8)

    fig.subplots_adjust(left=0.07, right=0.90, top=0.94, bottom=0.06)
    out_path = FIG / "grassmann_principal_angles.pdf"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[audit_46] figure saved: {out_path}")


def main() -> None:
    TBL.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    print("[audit_46] computing per-patient principal angles...")
    df = compute_per_patient()
    df.to_csv(TBL / "grassmann_principal_angles_per_patient.csv", index=False)
    print(
        f"[audit_46] per-patient table: {len(df)} rows -> "
        f"{TBL / 'grassmann_principal_angles_per_patient.csv'}"
    )

    sanity_check_chordal(df)

    cohort = aggregate_cohort(df)
    n_trace_lookup = cohort_n_trace(df)
    cohort["n_trace"] = cohort.apply(
        lambda r: n_trace_lookup.get((r["band"], r["k"]), (0, 10))[0], axis=1
    )
    cohort.to_csv(TBL / "grassmann_principal_angles_cohort.csv", index=False)
    print(
        f"[audit_46] cohort table: {len(cohort)} rows -> "
        f"{TBL / 'grassmann_principal_angles_cohort.csv'}"
    )

    render_figure(cohort, n_trace_lookup)


if __name__ == "__main__":
    main()
