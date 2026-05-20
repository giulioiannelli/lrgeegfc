"""Cohort matched-strength surrogate null on raw `D(τ) = 1/ρ̂(τ_max)`,
ALL 6 bands. Tests whether the matched-strength filter narrows the
"trace everywhere" within-baseline result on raw D down to a single
band (β) or remains permissive (multiple bands survive).

Same pipeline as `preprint_03_beta_matched_strength_raw_D.py`, looped
over 5 remaining bands (β already done; this script re-runs it for
sanity check + writes a unified all-bands CSV).

Output: `data/preprint/rho_split_raw_D/<band>_matched_strength.csv` (×6)
        `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv`
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, CACHE_ROOT
from lrg_eegfc.notebook import move_to_rootf
from lrg_eegfc.utils.surrogate.matched_strength import load_or_compute_surrogate_eigs
from lrg_eegfc.workflow.fc import load_fc_matrix

move_to_rootf(pathname="lrgeegfc")


PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["alpha", "beta", "low_gamma", "delta", "theta", "high_gamma"]
N_SURROGATES = 200
SWAP_FACTOR = 20
SEED = 20260511

HALVES_LRG_CACHE = Path("data/cache/imcoh_lrg_halves")
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"
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
    finite = np.isfinite(D)
    if not finite.all():
        cap = np.nanmax(D[finite]) if finite.any() else 1e6
        D = np.where(finite, D, cap)
    return D


def triu_vec(M): return M[np.triu_indices(M.shape[0], k=1)]


def rho_split(D_preA, D_preB, D_task, D_post):
    dt = triu_vec(D_task) - triu_vec(D_preA)
    dr = triu_vec(D_post) - triu_vec(D_preB)
    rho, _ = spearmanr(dt, dr)
    return float(rho)


def load_eig(p: Path):
    z = np.load(p, allow_pickle=True)
    return np.asarray(z["eigenvalues"], dtype=float), np.asarray(z["eigenvectors"], dtype=float)


def observed_eigs(pat: str, band: str):
    return {
        "rest_pre_A": load_eig(HALVES_LRG_CACHE / pat / f"{band}_rest_pre_A_lrg_imcoh-abs.npz"),
        "rest_pre_B": load_eig(HALVES_LRG_CACHE / pat / f"{band}_rest_pre_B_lrg_imcoh-abs.npz"),
        "task_test":  load_eig(IMCOH_LRG_CACHE / pat / f"{band}_task_test_lrg_imcoh-abs.npz"),
        "rest_post":  load_eig(IMCOH_LRG_CACHE / pat / f"{band}_rest_post_lrg_imcoh-abs.npz"),
    }


def load_W_for_surrogate(pat: str, phase: str, band: str) -> np.ndarray:
    if phase == "rest_pre_A":
        return np.asarray(np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_A_imcoh_abs.npy"), dtype=float)
    if phase == "rest_pre_B":
        return np.asarray(np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_B_imcoh_abs.npy"), dtype=float)
    if phase == "task_test":
        return np.asarray(load_fc_matrix(pat, "task_test", band, fc_method="imcoh_abs"), dtype=float)
    if phase == "rest_post":
        return np.asarray(load_fc_matrix(pat, "rest_post", band, fc_method="imcoh_abs"), dtype=float)
    raise ValueError(phase)


def per_cell(pat: str, band: str) -> dict | None:
    t0 = time.time()
    try:
        obs_eigs = observed_eigs(pat, band)
    except FileNotFoundError as e:
        print(f"[{pat} {band}] SKIP obs: {e}")
        return None

    D_obs = {p: D_raw_from_eig(*obs_eigs[p]) for p in obs_eigs}
    N = D_obs["rest_pre_A"].shape[0]
    if not all(D.shape == (N, N) for D in D_obs.values()):
        return None
    rho_obs = rho_split(D_obs["rest_pre_A"], D_obs["rest_pre_B"],
                        D_obs["task_test"], D_obs["rest_post"])

    rng = np.random.default_rng(SEED)
    surr_eigs = {}
    for phase in ("rest_pre_A", "rest_pre_B", "task_test", "rest_post"):
        W = load_W_for_surrogate(pat, phase, band)
        np.fill_diagonal(W, 0.0)
        W = np.clip(W, 0.0, 1.0)
        W = 0.5 * (W + W.T)
        evals, evecs = load_or_compute_surrogate_eigs(
            pat, band, phase, W,
            n_surr=N_SURROGATES, swap_factor=SWAP_FACTOR, seed=SEED, rng=rng,
            fc_method="imcoh_abs", verbose=False)
        surr_eigs[phase] = (evals, evecs)

    R = N_SURROGATES
    rho_surr = np.full(R, np.nan, dtype=float)
    for r in range(R):
        try:
            D_r = {}
            for phase, (evals_R, evecs_R) in surr_eigs.items():
                ev_r = evals_R[r]; vc_r = evecs_R[r]
                if not (np.all(np.isfinite(ev_r)) and np.all(np.isfinite(vc_r))):
                    raise ValueError("nan")
                D_r[phase] = D_raw_from_eig(ev_r, vc_r)
        except Exception:
            continue
        rho_surr[r] = rho_split(D_r["rest_pre_A"], D_r["rest_pre_B"],
                                 D_r["task_test"], D_r["rest_post"])

    surr = rho_surr[np.isfinite(rho_surr)]
    if surr.size == 0:
        return None
    p_one = float(np.mean(surr >= rho_obs))
    z = (rho_obs - float(np.mean(surr))) / (float(np.std(surr, ddof=1)) or float("nan"))
    out = {
        "patient": pat, "band": band, "N_nodes": int(N),
        "n_surrogates": int(surr.size),
        "obs_rho": rho_obs,
        "surr_mean_rho": float(np.mean(surr)),
        "surr_std_rho":  float(np.std(surr, ddof=1)),
        "surr_p5":  float(np.quantile(surr, 0.05)),
        "surr_p50": float(np.quantile(surr, 0.50)),
        "surr_p95": float(np.quantile(surr, 0.95)),
        "obs_z": float(z), "obs_p_one_sided": p_one,
        "dur_s": round(time.time() - t0, 1),
    }
    print(f"[{pat} {band}] obs={rho_obs:+.4f}  surr_p50={out['surr_p50']:+.4f}  "
          f"p1s={p_one:.4f}  z={out['obs_z']:+.2f}  ({out['dur_s']}s)")
    return out


cohort_rows = []
for band in BANDS:
    band_rows = []
    print(f"\n=== {band} ===\n")
    for pat in PATIENTS:
        r = per_cell(pat, band)
        if r is not None:
            band_rows.append(r)

    if not band_rows:
        continue
    df_band = pd.DataFrame(band_rows)
    df_band.to_csv(OUT_DIR / f"{band}_matched_strength.csv", index=False)

    obs = df_band["obs_rho"].values
    surr_med = df_band["surr_p50"].values
    diff = obs - surr_med
    n_above = int((df_band["obs_p_one_sided"] < 0.05).sum())
    try:
        W, p_paired = wilcoxon(diff, alternative="greater")
    except Exception:
        p_paired = float("nan")
    cohort_med_obs = float(np.median(obs))
    cohort_med_surr = float(np.median(surr_med))
    ratio = cohort_med_obs / cohort_med_surr if cohort_med_surr != 0 else float("inf")
    cohort_rows.append({
        "band": band, "n": len(df_band),
        "obs_median": cohort_med_obs,
        "surr_median": cohort_med_surr,
        "ratio_obs_over_surr": ratio,
        "wilcoxon_paired_p_one_sided": p_paired,
        "n_obs_p_lt_005": n_above,
    })

cohort = pd.DataFrame(cohort_rows)
cohort.to_csv(OUT_DIR / "all_bands_matched_strength_cohort.csv", index=False)

print("\n=== ALL-BANDS COHORT MATCHED-STRENGTH on raw D(τ_max) ===\n")
print(cohort.to_string(index=False, float_format=lambda v: f"{v:+.4f}"))
print(f"\nWrote {OUT_DIR}/all_bands_matched_strength_cohort.csv")
