"""Reproduce the cohort split-baseline rho_split pipeline on **raw D(τ)**.

The cached `ultrametric_matrix` field that the existing manuscript pipeline uses
is `cophenet(UPGMA(D))` — the UPGMA-cophenetic ultrametric, NOT the raw LRG
communication distance `D(τ) = 1/ρ̂_ij(τ)`. This script recomputes the same
within-baseline ρ_split statistic directly on the raw `D(τ)` at `τ = 1/λ_max`
and compares cohort-wide to the cached cophenet pipeline.

Pipeline definition (same as ctm_triangle / continuous_trace_matrix.py):
    D̃_taskT  = D(τ; task_test)
    D̃_rsPre_A = D(τ; rest_pre half A)
    D̃_rsPre_B = D(τ; rest_pre half B)
    D̃_rsPost = D(τ; rest_post)
    Δ_task = triu(D̃_taskT)  - triu(D̃_rsPre_A)
    Δ_rest = triu(D̃_rsPost) - triu(D̃_rsPre_B)
    ρ_split = Spearman(Δ_task, Δ_rest)

Comparison reference: `data/reports/imcoh_continuous_trace/per_cell_summary_split.csv`
(per-patient cophenet ρ_split) and `data/audit/ctm_triangle/cohort_summary.csv`.

Output: `data/preprint/rho_split_raw_D/beta_per_patient.csv` + a one-line cohort
summary printed to stdout.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.config.paths import IMCOH_LRG_CACHE
from lrg_eegfc.notebook import move_to_rootf

move_to_rootf(pathname="lrgeegfc")

# ----------------------------------------------------------------------------
# Cohort + paths
# ----------------------------------------------------------------------------

PATIENTS = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
BAND = "beta"
FC = "imcoh_abs"
HALVES_CACHE = Path("data/cache/imcoh_lrg_halves")
OUT_DIR = Path("data/preprint/rho_split_raw_D")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------------------------
# Core: raw D(τ) from cached eigvals / eigvecs
# ----------------------------------------------------------------------------

def D_raw_from_eig(eigvals: np.ndarray, eigvecs: np.ndarray) -> np.ndarray:
    """Build raw LRG communication distance D(τ) = 1/ρ̂_ij(τ) at τ = 1/λ_max.

    ρ̂(τ) = exp(-τL) / Tr exp(-τL) reconstructed from the Laplacian spectrum
    (eigvals, eigvecs) cached by compute_lrg_analysis.

    Returns a dense NxN array with zero diagonal and max-symmetrization
    (matching audit_63's `lrg_ultrametric_condensed` convention; symmetric to
    numerical precision in practice).
    """
    lam_max = float(np.max(eigvals))
    tau = 1.0 / lam_max
    diag_exp = np.exp(-tau * eigvals)
    Z = float(np.sum(diag_exp))
    rho = (eigvecs * diag_exp) @ eigvecs.T / Z  # NxN, symmetric, positive
    with np.errstate(divide="ignore", invalid="ignore"):
        D = 1.0 / rho
    D = np.maximum(D, D.T)
    np.fill_diagonal(D, 0.0)
    return D


def load_eig(npz_path: Path) -> tuple[np.ndarray, np.ndarray]:
    z = np.load(npz_path, allow_pickle=True)
    return np.asarray(z["eigenvalues"], dtype=float), np.asarray(z["eigenvectors"], dtype=float)


def triu_vec(M: np.ndarray) -> np.ndarray:
    n = M.shape[0]
    iu = np.triu_indices(n, k=1)
    return M[iu]


# ----------------------------------------------------------------------------
# Per-patient rho_split on raw D
# ----------------------------------------------------------------------------

def rho_split_raw(patient: str, band: str) -> dict:
    pre_A = HALVES_CACHE / patient / f"{band}_rest_pre_A_lrg_imcoh-abs.npz"
    pre_B = HALVES_CACHE / patient / f"{band}_rest_pre_B_lrg_imcoh-abs.npz"
    task  = IMCOH_LRG_CACHE / patient / f"{band}_task_test_lrg_imcoh-abs.npz"
    post  = IMCOH_LRG_CACHE / patient / f"{band}_rest_post_lrg_imcoh-abs.npz"

    for p in (pre_A, pre_B, task, post):
        if not p.exists():
            raise FileNotFoundError(p)

    eA = load_eig(pre_A);  D_preA  = D_raw_from_eig(*eA)
    eB = load_eig(pre_B);  D_preB  = D_raw_from_eig(*eB)
    et = load_eig(task);   D_task  = D_raw_from_eig(*et)
    ep = load_eig(post);   D_post  = D_raw_from_eig(*ep)

    n = D_preA.shape[0]
    if not all(D.shape == (n, n) for D in (D_preB, D_task, D_post)):
        raise RuntimeError(f"{patient} {band}: shape mismatch across phases")

    dt = triu_vec(D_task) - triu_vec(D_preA)
    dr = triu_vec(D_post) - triu_vec(D_preB)
    rho, _ = spearmanr(dt, dr)

    return {
        "patient": patient,
        "band": band,
        "N_p": n,
        "m_p": dt.size,
        "rho_split_raw_D": rho,
        "D_preA_min": float(D_preA[D_preA > 0].min()),
        "D_preA_max": float(D_preA.max()),
    }


# ----------------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------------

rows = [rho_split_raw(p, BAND) for p in PATIENTS]
df_raw = pd.DataFrame(rows)

cophe = pd.read_csv("data/reports/imcoh_continuous_trace/per_cell_summary_split.csv")
cophe_beta = cophe[cophe["band"] == BAND].rename(columns={"rho": "rho_split_cophenet"})[
    ["patient", "rho_split_cophenet"]
]
df = df_raw.merge(cophe_beta, on="patient", how="left")
df["delta_raw_minus_cophenet"] = df["rho_split_raw_D"] - df["rho_split_cophenet"]

print("\n=== beta rho_split: raw D(τ) vs cached cophenet ===\n")
print(
    df[
        [
            "patient", "N_p", "rho_split_raw_D", "rho_split_cophenet",
            "delta_raw_minus_cophenet", "D_preA_min", "D_preA_max",
        ]
    ].to_string(index=False, float_format=lambda v: f"{v:+.4f}")
)

# Cohort summary
rho_raw  = df["rho_split_raw_D"].values
rho_coph = df["rho_split_cophenet"].values
med_raw  = float(np.median(rho_raw))
med_coph = float(np.median(rho_coph))
n_pos_raw  = int((rho_raw  > 0).sum())
n_pos_coph = int((rho_coph > 0).sum())
W_raw,  p_raw  = wilcoxon(rho_raw,  alternative="greater", zero_method="wilcox")
W_coph, p_coph = wilcoxon(rho_coph, alternative="greater", zero_method="wilcox")
spearman_raw_vs_coph, _ = spearmanr(rho_raw, rho_coph)

print("\nCohort summary (n=10):")
print(f"  raw D(τ)  : median = {med_raw:+.4f}, n>0 = {n_pos_raw}/10, Wilcoxon p (one-sided) = {p_raw:.4f}")
print(f"  cophenet  : median = {med_coph:+.4f}, n>0 = {n_pos_coph}/10, Wilcoxon p (one-sided) = {p_coph:.4f}")
print(f"  per-patient Spearman(rho_raw, rho_cophenet) = {spearman_raw_vs_coph:+.4f}")

out_csv = OUT_DIR / f"{BAND}_per_patient.csv"
df.to_csv(out_csv, index=False)
print(f"\nWrote {out_csv}")
