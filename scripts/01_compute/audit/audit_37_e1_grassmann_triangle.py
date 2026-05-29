#!/usr/bin/env python3
"""Audit 37 — E1 Grassmann triangle (eigenspace overlap as a chained triangle).

Recasts the audit_27 pairwise frac_pos test into the same triangle form used
in measures 02 (D-rank) and 03 (KC). For each (p, b, k) compute the Grassmann
chordal distance between top-k non-trivial L̂ eigenspaces of phase pairs:

    d_chord(V^a_k, V^b_k) = sqrt(k - sum_i sigma_i^2)

where sigma_i are the singular values of (V^a_k)^T V^b_k. Eigenvectors are
columns 1..k+1 of the ascending-eigenvalue matrix (skipping the trivial zero
eigenvector). Triangle:

    T_E1(p, b, k) = d_chord(V^TT_k, V^RPost_k) - d_chord(V^RPre_k, V^TT_k)

Negative = trace direction. The pairwise frac_pos cohort verdict at audit_27
was negative; the triangle form has not been computed before.

Outputs
-------
``data/audit/e1_grassmann_triangle/Td_per_patient_per_band_k.csv``
``data/audit/e1_grassmann_triangle/cohort_summary_k.csv``
``data/outputs/figures/section_5_lrg_trace/grassmann_triangle/boxplots_per_k.pdf``
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.workflow.lrg import load_lrg_result

K_GRID = [3, 5, 8, 13]
OUT_DIR = ROOT / "data" / "audit" / "e1_grassmann_triangle"
FIG_DIR = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "grassmann_triangle"


def topk_basis(eigvecs: np.ndarray, k: int) -> np.ndarray:
    return np.ascontiguousarray(eigvecs[:, 1 : k + 1])


def chordal(V_a: np.ndarray, V_b: np.ndarray) -> float:
    """Chordal Grassmann distance: sqrt(k - sum sigma_i^2)."""
    M = V_a.T @ V_b
    sigma = np.linalg.svd(M, compute_uv=False)
    sigma = np.clip(sigma, 0.0, 1.0)
    return float(np.sqrt(max(V_a.shape[1] - float((sigma ** 2).sum()), 0.0)))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            try:
                EV = {}
                for phi in ("rest_pre", "task_test", "rest_post"):
                    res = load_lrg_result(pat, phi, band, fc_method="imcoh_abs")
                    if res is None or res.eigenvectors is None:
                        raise FileNotFoundError(f"missing eigenvectors: {pat} {phi} {band}")
                    EV[phi] = res.eigenvectors
                for k in K_GRID:
                    V_pre = topk_basis(EV["rest_pre"], k)
                    V_tt = topk_basis(EV["task_test"], k)
                    V_post = topk_basis(EV["rest_post"], k)
                    d_pre_tt = chordal(V_pre, V_tt)
                    d_tt_post = chordal(V_tt, V_post)
                    rows.append({
                        "patient": pat, "band": band, "k": k,
                        "d_pre_tt": d_pre_tt, "d_tt_post": d_tt_post,
                        # T_d > 0 = trace (rsPost closer to task than rsPre).
                        "T_E1": d_pre_tt - d_tt_post,
                    })
            except Exception as e:
                print(f"[audit_37] WARN {pat} {band}: {e}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "Td_per_patient_per_band_k.csv", index=False)

    summary_rows = []
    for b in BRAIN_BANDS_NAMES:
        for k in K_GRID:
            sub = df[(df["band"] == b) & (df["k"] == k)]
            v = sub["T_E1"].dropna().values
            n_trace = int((v > 0).sum())
            try:
                _, p = wilcoxon(v, alternative="greater")
            except Exception:
                p = float("nan")
            summary_rows.append({
                "band": b, "k": k,
                "n_trace": f"{n_trace}/{len(v)}",
                "n_trace_int": n_trace,
                "n_eligible": len(v),
                "T_median": float(np.median(v)) if len(v) else float("nan"),
                "T_iqr": float(np.subtract(*np.percentile(v, [75, 25]))) if len(v) else float("nan"),
                "wilcoxon_one_sided_p": float(p),
            })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_DIR / "cohort_summary_k.csv", index=False)
    print(summary.to_string(index=False))

    fig, axes = plt.subplots(1, len(K_GRID), figsize=(13, 3.4), sharey=False)
    for ax, k in zip(axes, K_GRID):
        data = [df[(df["band"] == b) & (df["k"] == k)]["T_E1"].dropna().values
                for b in BRAIN_BANDS_NAMES]
        ax.boxplot(data, labels=[BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                   showmeans=True, meanline=True, patch_artist=True,
                   boxprops=dict(facecolor="#e8e8e8"))
        for i, b in enumerate(BRAIN_BANDS_NAMES):
            v = df[(df["band"] == b) & (df["k"] == k)]["T_E1"].dropna().values
            jitter = (np.random.RandomState(7).rand(len(v)) - 0.5) * 0.12
            ax.scatter(np.full(len(v), i + 1) + jitter, v, s=10, color="#1f77b4", alpha=0.7, zorder=3)
        ax.axhline(0, color="0.4", lw=0.7, ls="--")
        ax.set_title(rf"$k = {k}$")
    axes[0].set_ylabel(r"$T_{E1}(p, b, k)$ -- positive = trace")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "boxplots_per_k.pdf")
    plt.close(fig)
    print(f"[audit_37] outputs at {OUT_DIR}; figure at {FIG_DIR}")


if __name__ == "__main__":
    main()
