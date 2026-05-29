"""All-bands ρ_split on raw `D(τ) = 1/ρ̂(τ_max)`, with three statistics:

    - **Spearman** (rank-based, current preprint primary)
    - **Pearson**  (linear, magnitude-sensitive)
    - **Cosine**   (vector direction in R^m, magnitude-sensitive)

Tests whether the β trace generalises to other bands when computed on the
raw LRG communication distance instead of the UPGMA-cophenetic wrap, and
whether the magnitude information that Spearman discards changes the picture.

Cohort: 10 patients (Pat_02..15). For each (patient, band):
    D_φ = 1/ρ̂(τ_max) built from cached eigvals/eigvecs in
    `data/cache/imcoh_lrg_halves/` (rest_pre_A, rest_pre_B)
    and `data/cache/imcoh_lrg/` (task_test, rest_post).

Reference: preprint_02_beta_rho_split_raw_D.py (β only, Spearman only).

Output: `data/preprint/rho_split_raw_D/all_bands_per_patient.csv`
        `data/preprint/rho_split_raw_D/all_bands_cohort.csv`
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, wilcoxon

from lrg_eegfc.config.paths import IMCOH_LRG_CACHE
from lrg_eegfc.notebook import move_to_rootf
from lrg_eegfc.visuals.styles import use_lrg_style
use_lrg_style()


move_to_rootf(pathname="lrgeegfc")


PATIENTS = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
HALVES_CACHE = Path("data/cache/imcoh_lrg_halves")
OUT_DIR = Path("data/preprint/rho_split_raw_D")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def D_raw_from_eig(eigvals: np.ndarray, eigvecs: np.ndarray) -> np.ndarray:
    lam_max = float(np.max(eigvals))
    tau = 1.0 / lam_max
    diag_exp = np.exp(-tau * eigvals)
    Z = float(np.sum(diag_exp))
    rho = (eigvecs * diag_exp) @ eigvecs.T / Z
    with np.errstate(divide="ignore", invalid="ignore"):
        D = 1.0 / rho
    D = np.maximum(D, D.T)
    np.fill_diagonal(D, 0.0)
    return D


def load_eig(p: Path) -> tuple[np.ndarray, np.ndarray]:
    z = np.load(p, allow_pickle=True)
    return np.asarray(z["eigenvalues"], dtype=float), np.asarray(z["eigenvectors"], dtype=float)


def triu_vec(M: np.ndarray) -> np.ndarray:
    return M[np.triu_indices(M.shape[0], k=1)]


def cosine(u: np.ndarray, v: np.ndarray) -> float:
    nu = float(np.linalg.norm(u))
    nv = float(np.linalg.norm(v))
    if nu == 0.0 or nv == 0.0:
        return float("nan")
    return float(np.dot(u, v) / (nu * nv))


def per_cell(patient: str, band: str) -> dict | None:
    paths = {
        "preA": HALVES_CACHE / patient / f"{band}_rest_pre_A_lrg_imcoh-abs.npz",
        "preB": HALVES_CACHE / patient / f"{band}_rest_pre_B_lrg_imcoh-abs.npz",
        "task": IMCOH_LRG_CACHE / patient / f"{band}_task_test_lrg_imcoh-abs.npz",
        "post": IMCOH_LRG_CACHE / patient / f"{band}_rest_post_lrg_imcoh-abs.npz",
    }
    for p in paths.values():
        if not p.exists():
            return None

    D = {k: D_raw_from_eig(*load_eig(p)) for k, p in paths.items()}
    n = D["preA"].shape[0]
    if not all(M.shape == (n, n) for M in D.values()):
        return None

    dt = triu_vec(D["task"]) - triu_vec(D["preA"])
    dr = triu_vec(D["post"]) - triu_vec(D["preB"])

    rho_s, _ = spearmanr(dt, dr)
    rho_p, _ = pearsonr(dt, dr)
    cos_v = cosine(dt, dr)

    return {
        "patient": patient,
        "band": band,
        "N_p": n,
        "m_p": dt.size,
        "rho_spearman": rho_s,
        "rho_pearson":  rho_p,
        "cosine":       cos_v,
        "D_preA_max":   float(D["preA"].max()),
        "dt_norm":      float(np.linalg.norm(dt)),
        "dr_norm":      float(np.linalg.norm(dr)),
    }


rows = []
for b in BANDS:
    for p in PATIENTS:
        r = per_cell(p, b)
        if r is None:
            print(f"SKIP {p} {b} (missing cache)")
        else:
            rows.append(r)

df = pd.DataFrame(rows)

# Per-patient cophenet reference
cophe = pd.read_csv("data/reports/imcoh_continuous_trace/per_cell_summary_split.csv")
cophe = cophe.rename(columns={"rho": "rho_spearman_cophenet"})[
    ["patient", "band", "rho_spearman_cophenet"]
]
df = df.merge(cophe, on=["patient", "band"], how="left")
df["delta_spearman_raw_minus_cophe"] = df["rho_spearman"] - df["rho_spearman_cophenet"]

out_per = OUT_DIR / "all_bands_per_patient.csv"
df.to_csv(out_per, index=False)

# Cohort summary
co_rows = []
for b in BANDS:
    sub = df[df["band"] == b]
    rho_s  = sub["rho_spearman"].values
    rho_p  = sub["rho_pearson"].values
    cos_v  = sub["cosine"].values
    rho_co = sub["rho_spearman_cophenet"].values

    def wilco(x):
        x = x[np.isfinite(x)]
        if x.size < 2:
            return float("nan")
        try:
            return float(wilcoxon(x, alternative="greater", zero_method="wilcox").pvalue)
        except Exception:
            return float("nan")

    co_rows.append({
        "band": b,
        "n_patients": len(sub),
        # raw D Spearman
        "median_spearman_raw":  float(np.median(rho_s)),
        "n_pos_spearman_raw":   int((rho_s > 0).sum()),
        "wilcoxon_p_spearman_raw": wilco(rho_s),
        # raw D Pearson
        "median_pearson_raw":   float(np.median(rho_p)),
        "n_pos_pearson_raw":    int((rho_p > 0).sum()),
        "wilcoxon_p_pearson_raw": wilco(rho_p),
        # raw D Cosine
        "median_cosine_raw":    float(np.median(cos_v)),
        "n_pos_cosine_raw":     int((cos_v > 0).sum()),
        "wilcoxon_p_cosine_raw": wilco(cos_v),
        # cophenet reference
        "median_spearman_cophenet":  float(np.median(rho_co)),
        "n_pos_spearman_cophenet":   int((rho_co > 0).sum()),
        "wilcoxon_p_spearman_cophenet": wilco(rho_co),
        # raw vs cophenet agreement
        "per_pat_spearman_raw_vs_cophe": float(spearmanr(rho_s, rho_co).statistic),
    })

cohort = pd.DataFrame(co_rows)
out_co = OUT_DIR / "all_bands_cohort.csv"
cohort.to_csv(out_co, index=False)

# Pretty print
print("\n=== All-bands ρ_split: raw D vs cophenet (per-band cohort summary) ===\n")
cols_to_show = [
    "band", "median_spearman_raw", "n_pos_spearman_raw", "wilcoxon_p_spearman_raw",
    "median_pearson_raw",   "n_pos_pearson_raw",   "wilcoxon_p_pearson_raw",
    "median_cosine_raw",    "n_pos_cosine_raw",    "wilcoxon_p_cosine_raw",
    "median_spearman_cophenet", "n_pos_spearman_cophenet", "wilcoxon_p_spearman_cophenet",
    "per_pat_spearman_raw_vs_cophe",
]
def fmt(v):
    if isinstance(v, float):
        return f"{v:+.4f}"
    return str(v)
print(cohort[cols_to_show].to_string(index=False, formatters={c: fmt for c in cols_to_show if c not in ("band", "n_pos_spearman_raw", "n_pos_pearson_raw", "n_pos_cosine_raw", "n_pos_spearman_cophenet")}))

print(f"\nWrote {out_per}")
print(f"Wrote {out_co}")
