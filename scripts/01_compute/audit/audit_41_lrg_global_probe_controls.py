#!/usr/bin/env python3
"""Audit 41 -- within-baseline null controls for LRG global probes.

For each of D-rank (measure 02), KC (measure 03), Grassmann (measure 04),
mirrors the substrate audit_28 / 30 / 31 controls at the LRG layer using
the half-baseline cache. Real T_d numbers are READ from the measure-02/03/04
tables; only the null is computed here.

Within-baseline null: replace the task_test phase with a second rest_pre
half. This makes a triangle of pure rest data: T_d_null =
d(RPreB, RPostA) - d(RPreA, RPreB). Under the null hypothesis "no task
involvement", the two distances should be similar and T_d_null ~ 0. If
the real T_d < T_d_null cohort-wide (Wilcoxon paired one-sided), the
real triangle's signed gap exceeds the within-baseline floor.

Pat_03 sampling-rate sensitivity (task #27): re-run cohort Wilcoxon for
KC at lambda=0 and lambda=1 with Pat_03 dropped, to verify that the
β trace finding does not collapse without the 1024 Hz outlier.

Outputs
-------
data/audit/lrg_global_probe_controls/{drank,kc,grassmann}_null.csv
data/audit/lrg_global_probe_controls/cohort_controls_summary.csv
data/audit/lrg_global_probe_controls/pat03_dropout.csv
data/outputs/figures/section_5_lrg_trace/lrg_controls/{drank,kc,grassmann}_real_vs_null.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.utils.metrics.tree_distance import kc_distance

HALVES = CACHE_ROOT / "imcoh_lrg_halves"
OUT = ROOT / "data" / "audit" / "lrg_global_probe_controls"
FIG = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "lrg_controls"

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {"delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
            "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$"}
KC_LAMBDAS = [0.0, 0.5, 1.0]
GRASSMANN_K = 13

REPORTS = ROOT / "data" / "reports" / "section_5_lrg_trace"
DRANK_REAL = REPORTS / "02_d_rank_triangle" / "tables" / "Td_per_patient_per_band.csv"
KC_REAL = REPORTS / "03_kc_lambda_triangle" / "tables" / "Td_per_patient_per_band_lambda.csv"
GRASSMANN_REAL = REPORTS / "04_grassmann_triangle" / "tables" / "Td_per_patient_per_band_k.csv"


# ---------------------------------------------------------------------------
# Distance primitives (mirror measure 02 / 03 / 04)
# ---------------------------------------------------------------------------

def d_S(u: np.ndarray, v: np.ndarray) -> float:
    """Spearman rank distance: 1 - Spearman(u, v)."""
    if u.ndim > 1:
        u = u[np.triu_indices_from(u, k=1)] if u.shape[0] == u.shape[1] else u.ravel()
    if v.ndim > 1:
        v = v[np.triu_indices_from(v, k=1)] if v.shape[0] == v.shape[1] else v.ravel()
    r, _ = spearmanr(u, v)
    return 1.0 - float(r) if not np.isnan(r) else 1.0


def d_P(u: np.ndarray, v: np.ndarray) -> float:
    if u.ndim > 1:
        u = u[np.triu_indices_from(u, k=1)] if u.shape[0] == u.shape[1] else u.ravel()
    if v.ndim > 1:
        v = v[np.triu_indices_from(v, k=1)] if v.shape[0] == v.shape[1] else v.ravel()
    r, _ = pearsonr(u, v)
    return 1.0 - float(r) if not np.isnan(r) else 1.0


def d_F(u: np.ndarray, v: np.ndarray) -> float:
    if u.ndim > 1:
        u = u[np.triu_indices_from(u, k=1)] if u.shape[0] == u.shape[1] else u.ravel()
    if v.ndim > 1:
        v = v[np.triu_indices_from(v, k=1)] if v.shape[0] == v.shape[1] else v.ravel()
    diff = np.linalg.norm(u - v)
    nrm = np.sqrt(np.linalg.norm(u) * np.linalg.norm(v))
    return float(diff / nrm) if nrm > 1e-12 else 1.0


def chordal(V_a: np.ndarray, V_b: np.ndarray) -> float:
    k = V_a.shape[1]
    M = V_a.T @ V_b
    sigma = np.linalg.svd(M, compute_uv=False)
    inner = min(float(np.sum(sigma ** 2)), float(k))
    return float(np.sqrt(max(0.0, k - inner)))


def topk_basis(V: np.ndarray, k: int) -> np.ndarray:
    return np.ascontiguousarray(V[:, 1: k + 1])


# ---------------------------------------------------------------------------
# Load helpers
# ---------------------------------------------------------------------------

def load_phase_lrg(pat, band, phase, half=False):
    if half:
        return load_lrg_result(pat, phase, band, "imcoh_abs", cache_root=HALVES)
    return load_lrg_result(pat, phase, band, "imcoh_abs")


# ---------------------------------------------------------------------------
# Null T_d on halves  (RPreA, RPreB, RPostA)  — RPreB substitutes for TT
# ---------------------------------------------------------------------------

def null_drank_for(pat, band):
    rpa = load_phase_lrg(pat, band, "rest_pre_A", half=True)
    rpb = load_phase_lrg(pat, band, "rest_pre_B", half=True)
    rpostA = load_phase_lrg(pat, band, "rest_post_A", half=True)
    if any(x is None for x in (rpa, rpb, rpostA)):
        return None
    Da, Db, Dpa = rpa.ultrametric_matrix, rpb.ultrametric_matrix, rpostA.ultrametric_matrix
    # T_d_null = d(RPreB, RPostA) - d(RPreA, RPreB)
    return dict(
        T_S=d_S(Db, Dpa) - d_S(Da, Db),
        T_P=d_P(Db, Dpa) - d_P(Da, Db),
        T_F=d_F(Db, Dpa) - d_F(Da, Db),
    )


def null_kc_for(pat, band):
    rpa = load_phase_lrg(pat, band, "rest_pre_A", half=True)
    rpb = load_phase_lrg(pat, band, "rest_pre_B", half=True)
    rpostA = load_phase_lrg(pat, band, "rest_post_A", half=True)
    if any(x is None for x in (rpa, rpb, rpostA)):
        return None
    Za, Zb, Zpa = rpa.linkage_matrix, rpb.linkage_matrix, rpostA.linkage_matrix
    out = {}
    for lam in KC_LAMBDAS:
        out[f"T_lam{lam:.1f}"] = (
            kc_distance(Zb, Zpa, lam=lam) - kc_distance(Za, Zb, lam=lam)
        )
    return out


def null_grassmann_for(pat, band):
    rpa = load_phase_lrg(pat, band, "rest_pre_A", half=True)
    rpb = load_phase_lrg(pat, band, "rest_pre_B", half=True)
    rpostA = load_phase_lrg(pat, band, "rest_post_A", half=True)
    if any(x is None or x.eigenvectors is None for x in (rpa, rpb, rpostA)):
        return None
    k = GRASSMANN_K
    if k + 1 > rpa.eigenvectors.shape[1]:
        return None
    Va = topk_basis(rpa.eigenvectors, k)
    Vb = topk_basis(rpb.eigenvectors, k)
    Vpa = topk_basis(rpostA.eigenvectors, k)
    return dict(
        T_E1=chordal(Vb, Vpa) - chordal(Va, Vb),
    )


# ---------------------------------------------------------------------------
# Real T_d numbers (read from existing CSVs)
# ---------------------------------------------------------------------------

def load_real_drank():
    df = pd.read_csv(DRANK_REAL)
    df = df[df["band"].isin(BANDS)].copy()
    return df.rename(columns={"T_S": "real_T_S", "T_P": "real_T_P", "T_F": "real_T_F"})


def load_real_kc():
    df = pd.read_csv(KC_REAL)
    df = df[df["band"].isin(BANDS)].copy()
    return df.rename(columns={"T_KC": "real_T_KC"})


def load_real_grassmann():
    df = pd.read_csv(GRASSMANN_REAL)
    df = df[(df["band"].isin(BANDS)) & (df["k"] == GRASSMANN_K)].copy()
    return df.rename(columns={"T_E1": "real_T_E1"})


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    # Load real T_d
    real_dr = load_real_drank()
    real_kc = load_real_kc()
    real_g = load_real_grassmann()

    # Compute null T_d
    rows_dr, rows_kc, rows_g = [], [], []
    for pat in PATIENTS:
        print(f"=== {pat} ===")
        for band in BANDS:
            r_dr = null_drank_for(pat, band)
            if r_dr is not None:
                rows_dr.append(dict(patient=pat, band=band, **{f"null_{k}": v for k, v in r_dr.items()}))
            r_kc = null_kc_for(pat, band)
            if r_kc is not None:
                rows_kc.append(dict(patient=pat, band=band, **{f"null_{k}": v for k, v in r_kc.items()}))
            r_g = null_grassmann_for(pat, band)
            if r_g is not None:
                rows_g.append(dict(patient=pat, band=band, **{f"null_{k}": v for k, v in r_g.items()}))
        print(f"  done")

    null_dr = pd.DataFrame(rows_dr)
    null_kc = pd.DataFrame(rows_kc)
    null_g = pd.DataFrame(rows_g)
    null_dr.to_csv(OUT / "drank_null.csv", index=False)
    null_kc.to_csv(OUT / "kc_null.csv", index=False)
    null_g.to_csv(OUT / "grassmann_null.csv", index=False)

    # Merge real + null
    dr = real_dr.merge(null_dr, on=["patient", "band"])
    kc = real_kc.merge(null_kc, on=["patient", "band"])
    g = real_g.merge(null_g, on=["patient", "band"])

    # Cohort gate per (probe, band, distance variant)
    cohort = []

    # D-rank: 3 distances
    for band in BANDS:
        for dist, real_col, null_col in [
            ("d_S", "real_T_S", "null_T_S"),
            ("d_P", "real_T_P", "null_T_P"),
            ("d_F", "real_T_F", "null_T_F"),
        ]:
            sub = dr[dr["band"] == band]
            real = sub[real_col].values
            null = sub[null_col].values
            # T_d > 0 = trace; "real > null" means real T_d more in the trace direction than null T_d.
            n_above_null = int(np.sum(real > null))
            try:
                w_p = float(wilcoxon(real, null, alternative="greater", zero_method="wilcox").pvalue)
            except ValueError:
                w_p = np.nan
            cohort.append(dict(probe="drank", band=band, variant=dist,
                               n_pat=len(sub),
                               n_real_above_null=n_above_null,
                               median_real=float(np.median(real)),
                               median_null=float(np.median(null)),
                               wilcoxon_p_real_gt_null=w_p))

    # KC: 3 lambdas
    for band in BANDS:
        for lam in KC_LAMBDAS:
            sub = kc[(kc["band"] == band) & (kc["lam"] == lam)]
            if sub.empty:
                continue
            real = sub["real_T_KC"].values
            null = sub[f"null_T_lam{lam:.1f}"].values
            n_above = int(np.sum(real > null))
            try:
                w_p = float(wilcoxon(real, null, alternative="greater", zero_method="wilcox").pvalue)
            except ValueError:
                w_p = np.nan
            cohort.append(dict(probe="kc", band=band, variant=f"lambda={lam:.1f}",
                               n_pat=len(sub),
                               n_real_above_null=n_above,
                               median_real=float(np.median(real)),
                               median_null=float(np.median(null)),
                               wilcoxon_p_real_gt_null=w_p))

    # Grassmann k=13
    for band in BANDS:
        sub = g[g["band"] == band]
        if sub.empty:
            continue
        real = sub["real_T_E1"].values
        null = sub["null_T_E1"].values
        n_above = int(np.sum(real > null))
        try:
            w_p = float(wilcoxon(real, null, alternative="greater", zero_method="wilcox").pvalue)
        except ValueError:
            w_p = np.nan
        cohort.append(dict(probe="grassmann", band=band, variant=f"k={GRASSMANN_K}",
                           n_pat=len(sub),
                           n_real_above_null=n_above,
                           median_real=float(np.median(real)),
                           median_null=float(np.median(null)),
                           wilcoxon_p_real_gt_null=w_p))

    coh_df = pd.DataFrame(cohort)
    coh_df.to_csv(OUT / "cohort_controls_summary.csv", index=False)
    print("\n=== cohort controls summary ===")
    print(coh_df.to_string(index=False))

    # Pat_03 dropout for KC β λ=0 and λ=1
    print("\n=== Pat_03 dropout (KC β at λ=0, λ=1) ===")
    drops = []
    kc_full = pd.read_csv(KC_REAL)
    for lam in [0.0, 1.0]:
        for band in ["beta"]:
            sub = kc_full[(kc_full["band"] == band) & (kc_full["lam"] == lam)]
            T_full = sub["T_KC"].values
            T_drop = sub[sub["patient"] != "Pat_03"]["T_KC"].values
            try:
                p_full = float(wilcoxon(T_full, alternative="greater", zero_method="wilcox").pvalue)
            except ValueError:
                p_full = np.nan
            try:
                p_drop = float(wilcoxon(T_drop, alternative="greater", zero_method="wilcox").pvalue)
            except ValueError:
                p_drop = np.nan
            drops.append(dict(probe="kc", band=band, variant=f"lambda={lam}",
                              n_full=len(T_full), n_drop=len(T_drop),
                              n_above_zero_full=int(np.sum(T_full > 0)),
                              n_above_zero_drop=int(np.sum(T_drop > 0)),
                              median_T_full=float(np.median(T_full)),
                              median_T_drop=float(np.median(T_drop)),
                              wilcoxon_p_full=p_full,
                              wilcoxon_p_drop=p_drop))
    drop_df = pd.DataFrame(drops)
    drop_df.to_csv(OUT / "pat03_dropout.csv", index=False)
    print(drop_df.to_string(index=False))

    # Per-band figure: real_T vs null_T scatter
    fig, axes = plt.subplots(3, len(BANDS), figsize=(2.2 * len(BANDS) + 2, 10),
                             sharex=False, sharey=False, squeeze=False)
    # rows: D-rank d_S, KC λ=0, Grassmann k=13
    rows_def = [
        ("D-rank d_S", dr, "real_T_S", "null_T_S"),
        ("KC $\\lambda=0$", kc[kc["lam"] == 0.0], "real_T_KC", "null_T_lam0.0"),
        (f"Grassmann k={GRASSMANN_K}", g, "real_T_E1", "null_T_E1"),
    ]
    for i, (label, frame, real_col, null_col) in enumerate(rows_def):
        for j, band in enumerate(BANDS):
            ax = axes[i, j]
            sub = frame[frame["band"] == band]
            if sub.empty:
                continue
            ax.scatter(sub[null_col], sub[real_col], color="#1a4f73", s=30, zorder=3)
            for _, r in sub.iterrows():
                ax.annotate(r["patient"][-2:], (r[null_col], r[real_col]),
                            fontsize=7, alpha=0.7, xytext=(3, 3),
                            textcoords="offset points")
            xy = np.concatenate([sub[real_col].values, sub[null_col].values])
            xy = xy[~np.isnan(xy)]
            if len(xy):
                lim = max(abs(xy.min()), abs(xy.max())) * 1.2
                ax.plot([-lim, lim], [-lim, lim], color="#888", lw=0.7,
                        ls="--", zorder=1, label="real = null")
                ax.fill_between([-lim, lim], [-lim, lim], -lim, color="#fff7e6",
                                alpha=0.4, zorder=0, label="real < null (trace)")
                ax.set_xlim(-lim, lim)
                ax.set_ylim(-lim, lim)
            ax.axhline(0, color="0.6", lw=0.4, ls=":")
            ax.axvline(0, color="0.6", lw=0.4, ls=":")
            if i == 0:
                ax.set_title(BAND_TEX[band])
            if j == 0:
                ax.set_ylabel(f"{label}\nreal $T_d$", fontsize=9)
            if i == 2:
                ax.set_xlabel(f"null $T_d$ (within-baseline)", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "real_vs_null_grid.pdf")
    plt.close(fig)

    print(f"\nWrote outputs to {OUT}\nFigures to {FIG}")


if __name__ == "__main__":
    main()
