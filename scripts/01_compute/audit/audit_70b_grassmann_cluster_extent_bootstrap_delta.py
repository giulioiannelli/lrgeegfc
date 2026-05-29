#!/usr/bin/env python3
"""Audit 70b — δ-only cluster-extent under bootstrap-of-Wilcoxon-Z per-k input.

One-off check requested 2026-05-28: does δ's "strong trace" verdict from
audit_70 (cluster_p_mass = 0.005) survive when the per-k cohort test is
replaced by a patient-bootstrap-of-Wilcoxon-Z (outlier-robust by
construction) instead of the paired Wilcoxon? Motivation: audit_70 LOO
flagged δ's cluster_p_mass as Pat_08-leveraged (drops to 0.055 when
Pat_08 is removed). The bootstrap is meant to bake that single-patient-
robustness into the test itself rather than diagnosing it post-hoc.

Pipeline reuses audit_70 verbatim through surrogate T_G computation; the
swap is only in the per-k cohort statistic (`bootstrap_wilcoxon_z_p`)
and propagates through `cluster_mass` and `longest_run_below`
identically.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

# Reuse the heavy lifting from audit_70.
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from audit_70_grassmann_cluster_extent import (  # noqa: E402
    PATIENTS, K_GRID, R, ALPHA_K,
    OBS_CSV, OUT_DIR as AUDIT70_OUT,
    compute_surr_T_G_band, longest_run_below, cluster_mass,
)

BAND = sys.argv[1] if len(sys.argv) > 1 else "delta"
B_BOOT = 2000
SEED = 20260528
PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = PROJECT_ROOT / f"data/audit/grassmann_cluster_extent_bootstrap_{BAND}"


def bootstrap_wilcoxon_z_p(d_per_pat_per_k, B, seed):
    """Per-k bootstrap-of-Wilcoxon-Z one-sided 'greater' p-value.

    Parameters
    ----------
    d_per_pat_per_k : (P, K) array of per-patient differences (obs − ref).
    B : int          number of patient-bootstrap resamples.
    seed : int       RNG seed.

    Returns
    -------
    p_boot : (K,) array. p_boot[k] = fraction of bootstrap resamples with
        Wilcoxon signed-rank Z ≤ 0 (i.e., the resampled cohort fails to
        show the trace direction).
    """
    P, K = d_per_pat_per_k.shape
    rng = np.random.default_rng(seed)
    # Bootstrap indices once, shared across all k → patient-level resampling
    # is consistent within each bootstrap, preserving any across-k coherence
    # at the patient level.
    idx = rng.integers(0, P, size=(B, P))            # (B, P)
    p_boot = np.full(K, 1.0)
    for ki in range(K):
        d_k = d_per_pat_per_k[:, ki]
        if not np.all(np.isfinite(d_k)):
            d_k = np.nan_to_num(d_k, nan=0.0)
        D_bs = d_k[idx]                              # (B, P)
        abs_D = np.abs(D_bs)
        # Average-ties rank along the patient axis.
        ranks = rankdata(abs_D, axis=1, method="average")  # (B, P)
        sign = np.sign(D_bs)
        # Zero-handling: scipy 'wilcox' drops zeros; with bootstrapped
        # continuous data, exact zeros are essentially impossible — but
        # exclude defensively.
        w_plus = np.sum(ranks * (sign > 0), axis=1)  # (B,)
        n_active = np.sum(sign != 0, axis=1).astype(float)
        # Null mean / std of W+ under H0 (no ties correction — n=P=10 is
        # small enough that the modest tie correction is negligible).
        mean_w = n_active * (n_active + 1.0) / 4.0
        std_w = np.sqrt(n_active * (n_active + 1.0) * (2.0 * n_active + 1.0) / 24.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            z = (w_plus - mean_w) / std_w
        # Fraction of bootstraps where Z ≤ 0 (trace direction fails).
        p_boot[ki] = float(np.mean(z <= 0.0))
    return p_boot


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print(f"=== δ-only bootstrap-input cluster-extent ===")
    print(f"  output: {OUT_DIR}")
    print(f"  B_boot={B_BOOT}, seed={SEED}, R={R}, alpha_per_k={ALPHA_K}")

    # --- Observed T_G[p, k]
    obs_df = pd.read_csv(OBS_CSV)
    P, K = len(PATIENTS), len(K_GRID)
    obs_T_G = np.full((P, K), np.nan)
    for pi, patient in enumerate(PATIENTS):
        sub = obs_df[(obs_df["band"] == BAND) & (obs_df["patient"] == patient)]
        for ki, k in enumerate(K_GRID):
            hit = sub[sub["k"] == k]
            if len(hit) == 1:
                obs_T_G[pi, ki] = float(hit["obs_T_G"].values[0])
    print(f"  obs rows: {np.sum(np.all(np.isfinite(obs_T_G), axis=1))}/{P} patients")

    # --- Surrogate T_G[p, r, k]
    print(f"  computing surrogate T_G (200 surrogates × 10 patients × 111 k) ...")
    surr_T_G = compute_surr_T_G_band(BAND, PATIENTS, K_GRID, R)
    surr_mean = surr_T_G.mean(axis=1)  # (P, K)
    print(f"  surrogate done in {time.time() - t0:.1f}s")

    # --- Observed cohort bootstrap-of-Wilcoxon-Z per k
    d_obs = obs_T_G - surr_mean
    obs_p_boot = bootstrap_wilcoxon_z_p(d_obs, B_BOOT, SEED)
    obs_LR = longest_run_below(obs_p_boot, ALPHA_K)
    obs_mass = cluster_mass(obs_p_boot, ALPHA_K)
    print(f"  obs (bootstrap input): LR={obs_LR}, mass={obs_mass:.3f}")

    # --- Null distribution
    print(f"  null distribution: {R} phantom tests, bootstrap input per phantom ...")
    null_LR = np.zeros(R, dtype=int)
    null_mass = np.zeros(R)
    for r in range(R):
        phantom = surr_T_G[:, r, :]
        mask = np.ones(R, dtype=bool); mask[r] = False
        ref = surr_T_G[:, mask, :].mean(axis=1)
        d_phantom = phantom - ref
        # Use phantom-specific seed so each null trial gets an independent
        # bootstrap; keep them deterministic for reproducibility.
        null_p_boot = bootstrap_wilcoxon_z_p(d_phantom, B_BOOT, SEED + 1 + r)
        null_LR[r] = longest_run_below(null_p_boot, ALPHA_K)
        null_mass[r] = cluster_mass(null_p_boot, ALPHA_K)
        if (r + 1) % 25 == 0:
            print(f"    phantom {r+1}/{R} (elapsed {time.time() - t0:.1f}s)")

    cluster_p_LR = (1 + int(np.sum(null_LR >= obs_LR))) / (R + 1)
    cluster_p_mass = (1 + int(np.sum(null_mass >= obs_mass))) / (R + 1)
    null_p95_LR = float(np.percentile(null_LR, 95))
    null_p95_mass = float(np.percentile(null_mass, 95))

    # --- LOO under bootstrap input (descriptive)
    print(f"  computing LOO cluster_p_mass under bootstrap input ...")
    loo_p_mass = np.zeros(P)
    loo_mass = np.zeros(P)
    for pi_drop in range(P):
        keep = np.ones(P, dtype=bool); keep[pi_drop] = False
        d_obs_loo = (obs_T_G[keep] - surr_T_G[keep].mean(axis=1))
        p_boot_loo = bootstrap_wilcoxon_z_p(d_obs_loo, B_BOOT, SEED + 10000 + pi_drop)
        obs_mass_loo = cluster_mass(p_boot_loo, ALPHA_K)
        null_mass_loo = np.zeros(R)
        for r in range(R):
            phantom_loo = surr_T_G[keep, r, :]
            mask_r = np.ones(R, dtype=bool); mask_r[r] = False
            ref_loo = surr_T_G[keep][:, mask_r, :].mean(axis=1)
            d_phantom_loo = phantom_loo - ref_loo
            null_p_loo = bootstrap_wilcoxon_z_p(
                d_phantom_loo, B_BOOT, SEED + 20000 + pi_drop * R + r
            )
            null_mass_loo[r] = cluster_mass(null_p_loo, ALPHA_K)
        loo_p_mass[pi_drop] = (1 + int(np.sum(null_mass_loo >= obs_mass_loo))) / (R + 1)
        loo_mass[pi_drop] = obs_mass_loo
        print(f"    drop {PATIENTS[pi_drop]}: obs_mass_loo={obs_mass_loo:.3f}, "
              f"cluster_p_mass_loo={loo_p_mass[pi_drop]:.4f}")
    loo_argmax = int(np.argmax(loo_p_mass))
    loo_max_p = float(loo_p_mass[loo_argmax])

    if cluster_p_mass < 0.01:
        verdict = "strong"
    elif cluster_p_mass < 0.05:
        verdict = "weak"
    else:
        verdict = "no_trace"

    # --- Compare to locked audit_70 numbers
    locked = pd.read_csv(AUDIT70_OUT / "cohort_summary.csv")
    locked_row = locked[locked.band == BAND].iloc[0]
    print()
    print(f"=== δ comparison (paired Wilcoxon vs bootstrap-of-Wilcoxon-Z) ===")
    print(f"  obs_LR:                   wilcoxon={int(locked_row.obs_longest_run):3d}    bootstrap={obs_LR:3d}")
    print(f"  obs_mass:                 wilcoxon={float(locked_row.obs_cluster_mass_neglog10p):6.3f} bootstrap={obs_mass:6.3f}")
    print(f"  null_p95_LR:              wilcoxon={float(locked_row.null_p95_LR):6.3f} bootstrap={null_p95_LR:6.3f}")
    print(f"  null_p95_mass:            wilcoxon={float(locked_row.null_p95_mass):6.3f} bootstrap={null_p95_mass:6.3f}")
    print(f"  cluster_p_longest_run:    wilcoxon={float(locked_row.cluster_p_longest_run):.4f}  bootstrap={cluster_p_LR:.4f}")
    print(f"  cluster_p_cluster_mass:   wilcoxon={float(locked_row.cluster_p_cluster_mass):.4f}  bootstrap={cluster_p_mass:.4f}")
    print(f"  LOO max p_mass (worst):   wilcoxon={float(locked_row.cluster_p_mass_loo_max):.4f}  bootstrap={loo_max_p:.4f} ({PATIENTS[loo_argmax]})")
    print(f"  verdict:                  wilcoxon={locked_row.verdict_cluster_extent:>10s}  bootstrap={verdict:>10s}")

    # --- Persist
    summary = pd.DataFrame([{
        "band": BAND, "n_patients": P, "R": R, "B_boot": B_BOOT,
        "seed": SEED, "alpha_per_k": ALPHA_K,
        "obs_longest_run": obs_LR,
        "obs_cluster_mass_neglog10p": obs_mass,
        "null_mean_LR": float(null_LR.mean()),
        "null_p95_LR": null_p95_LR,
        "null_max_LR": int(null_LR.max()),
        "null_mean_mass": float(null_mass.mean()),
        "null_p95_mass": null_p95_mass,
        "cluster_p_longest_run": cluster_p_LR,
        "cluster_p_cluster_mass": cluster_p_mass,
        "cluster_p_mass_loo_max": loo_max_p,
        "cluster_p_mass_loo_argmax_patient": PATIENTS[loo_argmax],
        "verdict_cluster_extent_bootstrap": verdict,
    }])
    summary.to_csv(OUT_DIR / "cohort_summary_band.csv", index=False)
    pd.DataFrame({
        "band": BAND, "k": K_GRID, "obs_p_boot": obs_p_boot,
    }).to_csv(OUT_DIR / "per_k_obs_p_boot_band.csv", index=False)
    pd.DataFrame({
        "band": BAND,
        "surrogate_idx": np.arange(R),
        "longest_run": null_LR,
        "cluster_mass": null_mass,
    }).to_csv(OUT_DIR / "null_distribution_band.csv", index=False)
    pd.DataFrame({
        "band": BAND,
        "dropped_patient": PATIENTS,
        "obs_cluster_mass_loo": loo_mass,
        "cluster_p_mass_loo": loo_p_mass,
    }).to_csv(OUT_DIR / "loo_cluster_p_mass_band.csv", index=False)
    print()
    print(f"Wrote: {OUT_DIR / 'cohort_summary_band.csv'}")
    print(f"Total elapsed: {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
