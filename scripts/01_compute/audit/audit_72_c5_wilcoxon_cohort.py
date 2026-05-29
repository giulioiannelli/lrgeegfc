#!/usr/bin/env python3
"""Audit 72 — C5 epi-zone exclusion Wilcoxon gate (replaces ≳80% retention).

Implements CONTROLS.md §C5 (locked 2026-05-19 Decision 10, post writing-agent
feedback). Two probes, two gates:

1. **Grassmann (all 6 bands)**: re-run audit_70 cluster-extent permutation
   on the **epi-X surrogate eigvec cache** at
   `data/cache/matched_strength_surrogate_epi_excluded_lrg/Pat_NN/`.
   Compute the resilient all-clusters cluster mass
   `T_G^*^epi-X(b) = Σ_{k : p_k^epi-X < α_k} (−log_10 p_k^epi-X)` and its
   empirical `cluster_p_mass^epi-X` from R=200 phantom permutations.
   C5 passes iff `cluster_p_mass^epi-X < 0.05`. No retention threshold.

2. **Cophenet (α only)**: one-sample one-sided Wilcoxon on per-patient
   `obs_rho^epi-X` (from `data/audit/alpha_epi_exclusion/per_patient.csv`)
   under `H_1: rho_split^epi-X > 0` (trace direction at cophenet).
   C5 passes iff `wilcoxon_one_sided_p < 0.05`.

LOO diagnostics per `feedback_no_single_patient_p_driven.md`: every gate
includes a `wilcoxon_loo_max_p_epiX` + `loo_argmax_patient` column.
Descriptive only, never a gate.

Inputs
------
- Cophenet α: `data/audit/alpha_epi_exclusion/per_patient.csv`.
- Grassmann all bands: surrogate eigvec cache
  `data/cache/matched_strength_surrogate_epi_excluded_lrg/Pat_NN/{band}_{phase}_epiX_R200_swap20_seed20260514_imcoh_abs.npz`.
- Observed T_G^epi-X per patient per band per k:
  `data/audit/grassmann_epi_exclusion/per_patient_per_band_per_k.csv`.

Outputs
-------
- `data/audit/grassmann_epi_exclusion/c5_wilcoxon_cohort.csv`
- `data/audit/alpha_epi_exclusion/c5_wilcoxon_cohort.csv`
"""
from __future__ import annotations

import gc
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

# Configuration
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
PATIENTS = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
K_GRID = list(range(2, 113))
R = 200
SWAP_FACTOR = 20
SEED_EPIX = 20260514
ALPHA_K = 0.05

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SURR_EPIX_ROOT = PROJECT_ROOT / "data/cache/matched_strength_surrogate_epi_excluded_lrg"
OBS_EPIX_CSV = PROJECT_ROOT / "data/audit/grassmann_epi_exclusion/per_patient_per_band_per_k.csv"
ALPHA_PER_PATIENT_CSV = PROJECT_ROOT / "data/audit/alpha_epi_exclusion/per_patient.csv"
BETA_PER_PATIENT_CSV = PROJECT_ROOT / "data/audit/beta_epi_exclusion/per_patient.csv"
OUT_GRASSMANN_CSV = PROJECT_ROOT / "data/audit/grassmann_epi_exclusion/c5_wilcoxon_cohort.csv"
OUT_ALPHA_CSV = PROJECT_ROOT / "data/audit/alpha_epi_exclusion/c5_wilcoxon_cohort.csv"
OUT_BETA_CSV = PROJECT_ROOT / "data/audit/beta_epi_exclusion/c5_wilcoxon_cohort.csv"


def chordal_distances_for_k_grid(eigvecs_a, eigvecs_b, k_grid):
    A = eigvecs_a[:, 1:]
    B = eigvecs_b[:, 1:]
    M = A.T @ B
    M2 = M * M
    csum = np.cumsum(np.cumsum(M2, axis=0), axis=1)
    out = np.zeros(len(k_grid))
    n_max = csum.shape[0]
    for i, k in enumerate(k_grid):
        if k > n_max:
            out[i] = np.nan; continue
        frob2 = csum[k - 1, k - 1]
        out[i] = np.sqrt(max(0.0, k - frob2))
    return out


def wilcoxon_per_k_greater(a, b):
    """Per-k paired Wilcoxon one-sided 'greater' (a > b) under project T_d > 0 = trace."""
    K = a.shape[1]
    p_vals = np.full(K, 1.0)
    for k in range(K):
        d = a[:, k] - b[:, k]
        d = d[np.isfinite(d)]
        if d.size < 3 or np.all(d == 0):
            continue
        try:
            _, p = wilcoxon(d, alternative="greater", zero_method="wilcox", correction=False)
            p_vals[k] = p
        except Exception:
            pass
    return p_vals


def cluster_mass(p_values, alpha):
    sig = (p_values < alpha) & np.isfinite(p_values)
    if not sig.any():
        return 0.0
    return float(np.sum(-np.log10(np.clip(p_values[sig], 1e-300, 1.0))))


def longest_run_below(p_values, alpha):
    sig = (p_values < alpha) & np.isfinite(p_values)
    max_run = cur = 0
    for s in sig:
        if s:
            cur += 1
            if cur > max_run: max_run = cur
        else:
            cur = 0
    return int(max_run)


def compute_surr_T_G_band_epiX(band, patients, k_grid, R):
    """Compute surrogate T_G[p, r, k] for one band under epi-X."""
    P, K = len(patients), len(k_grid)
    out = np.full((P, R, K), np.nan)
    for pi, patient in enumerate(patients):
        cache_dir = SURR_EPIX_ROOT / patient
        files = {
            "pre_A": cache_dir / f"{band}_rest_pre_A_epiX_R{R}_swap{SWAP_FACTOR}_seed{SEED_EPIX}_imcoh_abs.npz",
            "task":  cache_dir / f"{band}_task_test_epiX_R{R}_swap{SWAP_FACTOR}_seed{SEED_EPIX}_imcoh_abs.npz",
            "post":  cache_dir / f"{band}_rest_post_epiX_R{R}_swap{SWAP_FACTOR}_seed{SEED_EPIX}_imcoh_abs.npz",
        }
        for phase, path in files.items():
            if not path.exists():
                raise FileNotFoundError(f"epi-X surrogate eigvecs missing: {path}")
        eigvecs_pre = np.load(files["pre_A"])["eigvecs"]
        eigvecs_task = np.load(files["task"])["eigvecs"]
        eigvecs_post = np.load(files["post"])["eigvecs"]
        # k_max for this patient = N-1 (smallest dimension across surrogates).
        # eigvecs shape: (R, N, N). N differs per patient under epi-X.
        for r in range(R):
            d_taskpost = chordal_distances_for_k_grid(eigvecs_task[r], eigvecs_post[r], k_grid)
            d_pretask = chordal_distances_for_k_grid(eigvecs_pre[r], eigvecs_task[r], k_grid)
            out[pi, r, :] = d_pretask - d_taskpost  # T_G > 0 = trace
        del eigvecs_pre, eigvecs_task, eigvecs_post
        gc.collect()
    return out


def run_grassmann_c5():
    """C5 for Grassmann: cluster-mass test on epi-X eigvec cache."""
    print(f"\n=== Grassmann C5 (epi-X cluster-mass test) ===")
    obs_df = pd.read_csv(OBS_EPIX_CSV)
    rows = []
    for band in BANDS:
        print(f"\n  band: {band}")
        t0 = time.time()
        P, K = len(PATIENTS), len(K_GRID)

        # Observed T_G^epi-X per patient per k
        obs_T_G = np.full((P, K), np.nan)
        for pi, patient in enumerate(PATIENTS):
            sub = obs_df[(obs_df["band"] == band) & (obs_df["patient"] == patient)]
            for ki, k in enumerate(K_GRID):
                hit = sub[sub["k"] == k]
                if len(hit) == 1:
                    obs_T_G[pi, ki] = float(hit["obs_T_G"].values[0])
        # Surrogate T_G^epi-X per patient per r per k
        surr_T_G = compute_surr_T_G_band_epiX(band, PATIENTS, K_GRID, R)
        surr_mean = surr_T_G.mean(axis=1)

        # Observed cohort Wilcoxon and cluster mass
        obs_p = wilcoxon_per_k_greater(obs_T_G, surr_mean)
        obs_mass = cluster_mass(obs_p, ALPHA_K)
        obs_LR = longest_run_below(obs_p, ALPHA_K)

        # Null distribution via phantom-r tests
        null_mass = np.zeros(R)
        null_LR = np.zeros(R, dtype=int)
        for r in range(R):
            phantom = surr_T_G[:, r, :]
            mask_r = np.ones(R, dtype=bool); mask_r[r] = False
            ref = surr_T_G[:, mask_r, :].mean(axis=1)
            null_p = wilcoxon_per_k_greater(phantom, ref)
            null_mass[r] = cluster_mass(null_p, ALPHA_K)
            null_LR[r] = longest_run_below(null_p, ALPHA_K)

        cluster_p_mass = (1 + np.sum(null_mass >= obs_mass)) / (R + 1)
        cluster_p_LR = (1 + np.sum(null_LR >= obs_LR)) / (R + 1)

        # LOO diagnostic on cluster_p_mass
        print(f"    LOO ...", end="", flush=True)
        loo_p = np.zeros(P)
        for pi_drop in range(P):
            keep = np.ones(P, dtype=bool); keep[pi_drop] = False
            obs_T_G_loo = obs_T_G[keep]
            surr_T_G_loo = surr_T_G[keep]
            surr_mean_loo = surr_T_G_loo.mean(axis=1)
            obs_p_loo = wilcoxon_per_k_greater(obs_T_G_loo, surr_mean_loo)
            obs_mass_loo = cluster_mass(obs_p_loo, ALPHA_K)
            null_mass_loo = np.zeros(R)
            for r in range(R):
                phantom_loo = surr_T_G_loo[:, r, :]
                mask_r = np.ones(R, dtype=bool); mask_r[r] = False
                ref_loo = surr_T_G_loo[:, mask_r, :].mean(axis=1)
                null_p_loo = wilcoxon_per_k_greater(phantom_loo, ref_loo)
                null_mass_loo[r] = cluster_mass(null_p_loo, ALPHA_K)
            loo_p[pi_drop] = (1 + np.sum(null_mass_loo >= obs_mass_loo)) / (R + 1)
        loo_argmax = int(np.argmax(loo_p))
        loo_max_p = float(loo_p[loo_argmax])

        c5_pass = bool(cluster_p_mass < 0.05)
        elapsed = time.time() - t0
        print(f" obs_mass={obs_mass:.2f}, p_mass={cluster_p_mass:.4f}, "
              f"LOO max p={loo_max_p:.4f} ({PATIENTS[loo_argmax]}), "
              f"c5_pass={c5_pass}, elapsed={elapsed:.1f}s")

        rows.append({
            "band": band, "n_patients": P, "R": R,
            "obs_longest_run_epiX": obs_LR,
            "obs_cluster_mass_neglog10p_epiX": obs_mass,
            "cluster_p_longest_run_epiX": cluster_p_LR,
            "cluster_p_cluster_mass_epiX": cluster_p_mass,
            "cluster_p_mass_loo_max_epiX": loo_max_p,
            "cluster_p_mass_loo_argmax_patient": PATIENTS[loo_argmax],
            "c5_pass": c5_pass,
        })

    out_df = pd.DataFrame(rows)
    OUT_GRASSMANN_CSV.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT_GRASSMANN_CSV, index=False)
    print(f"\n  Wrote: {OUT_GRASSMANN_CSV}")
    print(out_df.to_string(index=False))


def _run_cophenet_c5(band: str, per_patient_csv: Path, out_csv: Path):
    """Generic cophenet C5: one-sample Wilcoxon on per-patient obs_rho^epi-X.

    Shared body for α and β cophenet C5 runs (the gate is identical;
    only the input/output paths differ). Returns the cohort row dict.
    """
    print(f"\n=== Cophenet {band} C5 (one-sample Wilcoxon on obs_rho^epi-X) ===")
    df = pd.read_csv(per_patient_csv)
    df = df.sort_values("patient").reset_index(drop=True)
    patients = df["patient"].tolist()
    obs_rho = df["obs_rho"].to_numpy()

    _, p_main = wilcoxon(obs_rho, alternative="greater",
                         zero_method="wilcox", correction=False)

    n = len(patients)
    p_loo = np.zeros(n)
    for i in range(n):
        keep = np.ones(n, dtype=bool); keep[i] = False
        try:
            _, p_loo[i] = wilcoxon(obs_rho[keep], alternative="greater",
                                    zero_method="wilcox", correction=False)
        except Exception:
            p_loo[i] = 1.0
    loo_argmax = int(np.argmax(p_loo))
    loo_max_p = float(p_loo[loo_argmax])

    med = float(np.median(obs_rho))
    c5_pass = bool(float(p_main) < 0.05)

    print(f"  obs_rho median={med:.4f}, "
          f"wilcoxon p (one-sided greater)={p_main:.4f}, "
          f"LOO max p={loo_max_p:.4f} ({patients[loo_argmax]}), "
          f"c5_pass={c5_pass}")

    row = pd.DataFrame([{
        "band": band, "probe": "cophenet",
        "n_patients": n,
        "obs_rho_median_epiX": med,
        "wilcoxon_one_sided_p_epiX": float(p_main),
        "wilcoxon_loo_max_p_epiX": loo_max_p,
        "wilcoxon_loo_argmax_patient": patients[loo_argmax],
        "c5_pass": c5_pass,
    }])
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    row.to_csv(out_csv, index=False)
    print(f"  Wrote: {out_csv}")


def run_alpha_cophenet_c5():
    """C5 for cophenet α: one-sample Wilcoxon on per-patient obs_rho^epi-X."""
    _run_cophenet_c5("alpha", ALPHA_PER_PATIENT_CSV, OUT_ALPHA_CSV)


def run_beta_cophenet_c5():
    """C5 for cophenet β: one-sample Wilcoxon on per-patient obs_rho^epi-X.

    Requires `audit_68_beta_epi_exclusion.py` to have produced
    `data/audit/beta_epi_exclusion/per_patient.csv` first.
    """
    if not BETA_PER_PATIENT_CSV.exists():
        print(f"\n=== Cophenet β C5: SKIP — {BETA_PER_PATIENT_CSV} not found ===")
        print("  Run audit_68_beta_epi_exclusion.py first.")
        return
    _run_cophenet_c5("beta", BETA_PER_PATIENT_CSV, OUT_BETA_CSV)


if __name__ == "__main__":
    run_alpha_cophenet_c5()
    run_beta_cophenet_c5()
    run_grassmann_c5()
