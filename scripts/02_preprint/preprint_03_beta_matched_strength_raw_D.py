"""β cohort matched-strength surrogate null on raw `D(τ) = 1/ρ̂(τ_max)`.

Replicates the audit_63 split-baseline pipeline but with the distance
function swapped from `lrg_ultrametric_condensed` (UPGMA-cophenet) to
the raw LRG communication distance `D = 1/ρ̂(τ_max)`. The point is to
test whether the matched-strength cohort separation (23.7× cohort ratio
on cophenet) survives when the cophenet step is removed.

Pipeline:
    For each (patient, β) cell:
      Observed:
        - For each phase φ ∈ {rsPre_A, rsPre_B, taskT, rsPost}, build
          raw D from the cached eigvals/eigvecs in `imcoh_lrg_halves`
          / `imcoh_lrg`.
        - Compute Δ_task, Δ_rest, ρ_split = Spearman(...).
      Surrogate (R=200 replicates, seed=20260511, swap_factor=20):
        - Use cached matched-strength surrogate eigvecs at
          `data/cache/matched_strength_surrogate_lrg/` for {rsPre_A,
          taskT, rsPost}. Build rsPre_B surrogate eigs on-the-fly via
          `lrg_eegfc.utils.surrogate.matched_strength.load_or_compute_surrogate_eigs`
          (will cache).
        - For each replicate r: build D_r per phase from eigvals/eigvecs,
          compute ρ_split_r.
        - Compare observed ρ_split to surrogate distribution.

Cohort verdict: Wilcoxon paired (obs_rho − surr_p50) one-sided.

Output: `data/preprint/rho_split_raw_D/beta_matched_strength.csv`
        `data/preprint/rho_split_raw_D/beta_matched_strength_cohort.txt`
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
from lrg_eegfc.visuals.styles import use_lrg_style
use_lrg_style()


move_to_rootf(pathname="lrgeegfc")


PATIENTS = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
BAND = "beta"
N_SURROGATES = 200
SWAP_FACTOR = 20
SEED = 20260511

HALVES_LRG_CACHE = Path("data/cache/imcoh_lrg_halves")
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"
OUT_DIR = Path("data/preprint/rho_split_raw_D")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------------------------
# Raw D(τ) from eigvals / eigvecs
# ----------------------------------------------------------------------------

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


def triu_vec(M: np.ndarray) -> np.ndarray:
    return M[np.triu_indices(M.shape[0], k=1)]


def rho_split(D_preA, D_preB, D_task, D_post) -> float:
    dt = triu_vec(D_task) - triu_vec(D_preA)
    dr = triu_vec(D_post) - triu_vec(D_preB)
    rho, _ = spearmanr(dt, dr)
    return float(rho)


# ----------------------------------------------------------------------------
# Observed (real) eigvals/eigvecs loaders
# ----------------------------------------------------------------------------

def load_eig(p: Path) -> tuple[np.ndarray, np.ndarray]:
    z = np.load(p, allow_pickle=True)
    return np.asarray(z["eigenvalues"], dtype=float), np.asarray(z["eigenvectors"], dtype=float)


def observed_eigs(pat: str, band: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    return {
        "rest_pre_A": load_eig(HALVES_LRG_CACHE / pat / f"{band}_rest_pre_A_lrg_imcoh-abs.npz"),
        "rest_pre_B": load_eig(HALVES_LRG_CACHE / pat / f"{band}_rest_pre_B_lrg_imcoh-abs.npz"),
        "task_test":  load_eig(IMCOH_LRG_CACHE / pat / f"{band}_task_test_lrg_imcoh-abs.npz"),
        "rest_post":  load_eig(IMCOH_LRG_CACHE / pat / f"{band}_rest_post_lrg_imcoh-abs.npz"),
    }


def load_W_for_surrogate(pat: str, phase: str, band: str) -> np.ndarray:
    """Load adjacency matrix used to generate surrogates (only consulted on cache miss)."""
    if phase == "rest_pre_A":
        return np.asarray(np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_A_imcoh_abs.npy"), dtype=float)
    if phase == "rest_pre_B":
        return np.asarray(np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_B_imcoh_abs.npy"), dtype=float)
    if phase == "task_test":
        return np.asarray(load_fc_matrix(pat, "task_test", band, fc_method="imcoh_abs"), dtype=float)
    if phase == "rest_post":
        return np.asarray(load_fc_matrix(pat, "rest_post", band, fc_method="imcoh_abs"), dtype=float)
    raise ValueError(phase)


# ----------------------------------------------------------------------------
# Per-cell pipeline
# ----------------------------------------------------------------------------

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
        print(f"[{pat} {band}] SKIP: phase shape mismatch")
        return None

    rho_obs = rho_split(D_obs["rest_pre_A"], D_obs["rest_pre_B"],
                        D_obs["task_test"], D_obs["rest_post"])

    rng = np.random.default_rng(SEED)
    surr_eigs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for phase in ("rest_pre_A", "rest_pre_B", "task_test", "rest_post"):
        W = load_W_for_surrogate(pat, phase, band)
        # Match audit_63's safety: zero diag, clip to [0,1], symmetric
        np.fill_diagonal(W, 0.0)
        W = np.clip(W, 0.0, 1.0)
        W = 0.5 * (W + W.T)
        evals, evecs = load_or_compute_surrogate_eigs(
            pat, band, phase, W,
            n_surr=N_SURROGATES, swap_factor=SWAP_FACTOR, seed=SEED, rng=rng,
            fc_method="imcoh_abs", verbose=False,
        )
        surr_eigs[phase] = (evals, evecs)

    R = N_SURROGATES
    rho_surr = np.full(R, np.nan, dtype=float)
    for r in range(R):
        # Skip replicates where any phase has NaN (failed strength preservation)
        try:
            D_r = {}
            for phase, (evals_R, evecs_R) in surr_eigs.items():
                ev_r = evals_R[r]
                vc_r = evecs_R[r]
                if not (np.all(np.isfinite(ev_r)) and np.all(np.isfinite(vc_r))):
                    raise ValueError("nan in eigs")
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
        "patient": pat,
        "band": band,
        "N_nodes": int(N),
        "n_surrogates": int(surr.size),
        "obs_rho": rho_obs,
        "surr_mean_rho": float(np.mean(surr)),
        "surr_std_rho":  float(np.std(surr, ddof=1)),
        "surr_p5":  float(np.quantile(surr, 0.05)),
        "surr_p25": float(np.quantile(surr, 0.25)),
        "surr_p50": float(np.quantile(surr, 0.50)),
        "surr_p75": float(np.quantile(surr, 0.75)),
        "surr_p95": float(np.quantile(surr, 0.95)),
        "obs_z":   float(z),
        "obs_p_one_sided": p_one,
        "dur_s":   round(time.time() - t0, 1),
    }
    print(f"[{pat} {band}] obs={rho_obs:+.4f}  surr_p50={out['surr_p50']:+.4f}  "
          f"p1s={p_one:.4f}  z={out['obs_z']:+.2f}  ({out['dur_s']}s)")
    return out


# ----------------------------------------------------------------------------
# Run cohort
# ----------------------------------------------------------------------------

rows = []
for pat in PATIENTS:
    r = per_cell(pat, BAND)
    if r is not None:
        rows.append(r)

df = pd.DataFrame(rows)
out_csv = OUT_DIR / f"{BAND}_matched_strength.csv"
df.to_csv(out_csv, index=False)

obs   = df["obs_rho"].values
surr_med = df["surr_p50"].values
diff = obs - surr_med
n_above = int((df["obs_p_one_sided"] < 0.05).sum())

W_paired, p_paired = wilcoxon(diff, alternative="greater")
cohort_med_obs   = float(np.median(obs))
cohort_med_surr  = float(np.median(surr_med))
ratio = cohort_med_obs / cohort_med_surr if cohort_med_surr != 0 else float("inf")

summary = (
    f"\n=== β matched-strength on raw D(τ) — cohort verdict (n={len(df)}) ===\n"
    f"  obs median ρ_split           = {cohort_med_obs:+.4f}\n"
    f"  surrogate median (per-pat p50)= {cohort_med_surr:+.4f}\n"
    f"  cohort ratio obs/surr-median = {ratio:+.4g}×\n"
    f"  Wilcoxon paired (obs > surr_p50) p = {p_paired:.4f}\n"
    f"  n patients with obs_p < 0.05  = {n_above}/{len(df)}\n"
    "\n  comparison to cophenet matched-strength (reference, from audit_63):\n"
    "    obs median +0.222, surr median +0.0094, ratio 23.7×, n>0p=0.005 7/10, Wilcoxon p≈0.005\n"
)
print(summary)
(OUT_DIR / f"{BAND}_matched_strength_cohort.txt").write_text(summary)
print(f"Wrote {out_csv}")
