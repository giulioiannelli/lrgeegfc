#!/usr/bin/env python3
"""Audit 70 — Grassmann cluster-extent permutation null.

Replaces the 8-cell hardcoded contiguous-run threshold (VERDICT_LEDGER.md
2026-05-18 Decision 6) with a data-derived gate using the R=200
matched-strength surrogates from audit_66.

For each band ∈ {δ, θ, α, β, γ_l, γ_h}:
  1. Read observed T_G(k) per patient per k from audit_66's CSV.
  2. Recompute surrogate T_G(k) per patient per surrogate per k from cached
     surrogate eigvecs at data/cache/matched_strength_surrogate_lrg/.
  3. Observed longest contiguous run = max contiguous k with paired
     Wilcoxon one-sided-less (obs < surr_mean across patients) p < 0.05.
  4. Null distribution: for each r ∈ 1..R, treat surrogate r as the
     "phantom observation", reference = mean of remaining R−1 surrogates,
     paired Wilcoxon per k → phantom longest run. R=200 null values.
  5. Cluster p = (1 + #(null_LR ≥ obs_LR)) / (R + 1).

Outputs
-------
data/audit/grassmann_cluster_extent/
    cohort_summary.csv     — per-band cluster_p + verdict
    null_distribution.csv  — per-band per-r phantom longest run
    per_k_obs_p.csv        — per-band per-k observed Wilcoxon p
    README.md

Method notes
------------
* Chordal Grassmann distance computed via the trace identity
  d_chord² = k − ‖Uᵀ V‖_F², no SVD needed (matches the closed form
  Σ sin²θ_i = k − Σ σ_i² where σ_i are singular values of UᵀV).
* The non-trivial mode subspace skips column 0 (λ_1 = 0) and takes
  columns [1:k+1], matching audit_66's `topk_basis()` convention.
* Per-k Wilcoxon uses one-sided alternative='less' since the trace
  direction is T_G < 0 (rsPost closer to task than rsPre is to task).
"""
from __future__ import annotations

import gc
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

# Configuration — frozen
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
PATIENTS = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
K_GRID = list(range(2, 113))  # 111 values, matches audit_66
R = 200
SWAP_FACTOR = 20
SEED = 20260511
ALPHA_K = 0.05  # per-k gate

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SURR_CACHE_ROOT = PROJECT_ROOT / "data/cache/matched_strength_surrogate_lrg"
OBS_CSV = PROJECT_ROOT / "data/audit/grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv"
OUT_DIR = PROJECT_ROOT / "data/audit/grassmann_cluster_extent"


def chordal_distances_for_k_grid(eigvecs_a, eigvecs_b, k_grid):
    """Chordal Grassmann distance d_chord(k) for each k in k_grid.

    d_chord² = k − ‖Aᵀ B‖_F²  where A, B are N×k orthonormal columns
    (eigvecs columns 1..k+1, skipping the trivial mode at column 0).

    Implementation: compute the (N−1)×(N−1) cross-correlation once and
    derive d_chord for every k via 2-D cumulative sum of the element-wise
    squared matrix.
    """
    A = eigvecs_a[:, 1:]  # (N, N-1)
    B = eigvecs_b[:, 1:]  # (N, N-1)
    M = A.T @ B           # (N-1, N-1)
    M2 = M * M
    csum = np.cumsum(np.cumsum(M2, axis=0), axis=1)  # cumsum[i, j] = sum over [0:i+1, 0:j+1]

    out = np.zeros(len(k_grid))
    n_max = csum.shape[0]
    for i, k in enumerate(k_grid):
        if k > n_max:
            out[i] = np.nan
            continue
        frob2 = csum[k - 1, k - 1]
        out[i] = np.sqrt(max(0.0, k - frob2))
    return out


def compute_surr_T_G_band(band, patients, k_grid, R):
    """Compute surrogate T_G[p, r, k] for one band, shape (P, R, K)."""
    P = len(patients)
    K = len(k_grid)
    out = np.full((P, R, K), np.nan)

    for pi, patient in enumerate(patients):
        cache_dir = SURR_CACHE_ROOT / patient
        files = {
            "pre_A": cache_dir / f"{band}_rest_pre_A_R{R}_swap{SWAP_FACTOR}_seed{SEED}_imcoh_abs.npz",
            "task":  cache_dir / f"{band}_task_test_R{R}_swap{SWAP_FACTOR}_seed{SEED}_imcoh_abs.npz",
            "post":  cache_dir / f"{band}_rest_post_R{R}_swap{SWAP_FACTOR}_seed{SEED}_imcoh_abs.npz",
        }
        for phase, path in files.items():
            if not path.exists():
                raise FileNotFoundError(f"Surrogate eigvecs missing: {path}")

        eigvecs_pre = np.load(files["pre_A"])["eigvecs"]   # (R, N, N)
        eigvecs_task = np.load(files["task"])["eigvecs"]
        eigvecs_post = np.load(files["post"])["eigvecs"]

        for r in range(R):
            d_taskpost = chordal_distances_for_k_grid(eigvecs_task[r], eigvecs_post[r], k_grid)
            d_pretask = chordal_distances_for_k_grid(eigvecs_pre[r], eigvecs_task[r], k_grid)
            out[pi, r, :] = d_taskpost - d_pretask

        del eigvecs_pre, eigvecs_task, eigvecs_post
        gc.collect()

    return out


def longest_run_below(p_values, alpha):
    """Longest contiguous run with p < alpha. Ignores NaN (treats as not-sig)."""
    sig = (p_values < alpha) & np.isfinite(p_values)
    max_run = cur = 0
    for s in sig:
        if s:
            cur += 1
            max_run = max(max_run, cur)
        else:
            cur = 0
    return int(max_run)


def wilcoxon_per_k_less(a, b):
    """Per-k paired Wilcoxon one-sided 'less' (a < b). a, b: (P, K) → p-vals (K,)."""
    K = a.shape[1]
    p_vals = np.full(K, 1.0)
    for k in range(K):
        d = a[:, k] - b[:, k]
        d = d[np.isfinite(d)]
        if d.size < 3 or np.all(d == 0):
            continue
        try:
            _, p = wilcoxon(d, alternative="less", zero_method="wilcox", correction=False)
            p_vals[k] = p
        except Exception:
            pass
    return p_vals


def cluster_mass(p_values, alpha):
    """Sum of -log10(p) over ALL k-cells with p < alpha (resilient cluster extent).

    Sums significance across every contiguous-significant cluster, not just the
    longest one. Resilient to a single missing k that fragments a long run:
    breaking a length-20 cluster into 10+10 keeps the total mass unchanged,
    whereas the longest-run-only version would halve it. The per-cell
    α threshold (alpha) is the cluster-forming threshold; the statistic is
    the sum of -log10(p_k) over the cells that pass it.
    """
    sig = (p_values < alpha) & np.isfinite(p_values)
    if not sig.any():
        return 0.0
    return float(np.sum(-np.log10(np.clip(p_values[sig], 1e-300, 1.0))))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Cluster-extent permutation — output at: {OUT_DIR}")

    obs_df = pd.read_csv(OBS_CSV)
    summary_rows, null_rows, per_k_rows, loo_rows = [], [], [], []
    overall_t0 = time.time()

    for band in BANDS:
        print(f"\n=== {band} ===")
        t0 = time.time()
        P, K = len(PATIENTS), len(K_GRID)

        # --- Observed T_G per patient per k from audit_66 CSV
        obs_T_G = np.full((P, K), np.nan)
        for pi, patient in enumerate(PATIENTS):
            sub = obs_df[(obs_df["band"] == band) & (obs_df["patient"] == patient)]
            for ki, k in enumerate(K_GRID):
                hit = sub[sub["k"] == k]
                if len(hit) == 1:
                    obs_T_G[pi, ki] = float(hit["obs_T_G"].values[0])
        n_obs_complete = np.sum(np.all(np.isfinite(obs_T_G), axis=1))
        print(f"  obs rows complete: {n_obs_complete}/{P} patients × {K} k-cells")

        # --- Surrogate T_G per patient per r per k
        print(f"  computing surrogate T_G ...")
        surr_T_G = compute_surr_T_G_band(band, PATIENTS, K_GRID, R)
        surr_mean = surr_T_G.mean(axis=1)  # (P, K)
        print(f"  surrogate done in {time.time() - t0:.1f}s")

        # --- Observed cohort Wilcoxon per k
        obs_p = wilcoxon_per_k_less(obs_T_G, surr_mean)
        obs_LR = longest_run_below(obs_p, ALPHA_K)
        obs_mass = cluster_mass(obs_p, ALPHA_K)
        for ki, k in enumerate(K_GRID):
            per_k_rows.append({"band": band, "k": k, "obs_p_one_sided_less": obs_p[ki]})

        # --- Null distribution via phantom-r tests
        print(f"  building null distribution ({R} phantom tests) ...")
        null_LR = np.zeros(R, dtype=int)
        null_mass = np.zeros(R)
        for r in range(R):
            phantom = surr_T_G[:, r, :]
            mask = np.ones(R, dtype=bool)
            mask[r] = False
            ref = surr_T_G[:, mask, :].mean(axis=1)
            null_p = wilcoxon_per_k_less(phantom, ref)
            null_LR[r] = longest_run_below(null_p, ALPHA_K)
            null_mass[r] = cluster_mass(null_p, ALPHA_K)
            null_rows.append({
                "band": band, "surrogate_idx": r,
                "longest_run": null_LR[r], "cluster_mass": null_mass[r],
            })

        # --- Cluster p-values
        cluster_p_LR = (1 + np.sum(null_LR >= obs_LR)) / (R + 1)
        cluster_p_mass = (1 + np.sum(null_mass >= obs_mass)) / (R + 1)

        # --- Leave-one-out diagnostic (per `feedback_no_single_patient_p_driven.md`):
        # for each patient i, drop them from BOTH observation and null,
        # recompute cluster_p_mass on n-1 cohort. Reports the worst-case
        # (max) LOO p-value and the patient whose removal produced it.
        # Descriptive only — never a gate. Verdict is unchanged by LOO.
        print(f"  computing LOO cluster_p_mass diagnostic ...")
        loo_p_mass = np.zeros(P)
        loo_mass = np.zeros(P)
        loo_LR = np.zeros(P, dtype=int)
        for pi_drop in range(P):
            keep = np.ones(P, dtype=bool); keep[pi_drop] = False
            obs_T_G_loo = obs_T_G[keep]
            surr_T_G_loo = surr_T_G[keep]
            surr_mean_loo = surr_T_G_loo.mean(axis=1)
            obs_p_loo = wilcoxon_per_k_less(obs_T_G_loo, surr_mean_loo)
            obs_mass_loo = cluster_mass(obs_p_loo, ALPHA_K)
            obs_LR_loo = longest_run_below(obs_p_loo, ALPHA_K)
            null_mass_loo = np.zeros(R)
            for r in range(R):
                phantom_loo = surr_T_G_loo[:, r, :]
                mask_r = np.ones(R, dtype=bool); mask_r[r] = False
                ref_loo = surr_T_G_loo[:, mask_r, :].mean(axis=1)
                null_p_loo = wilcoxon_per_k_less(phantom_loo, ref_loo)
                null_mass_loo[r] = cluster_mass(null_p_loo, ALPHA_K)
            loo_p_mass[pi_drop] = (1 + np.sum(null_mass_loo >= obs_mass_loo)) / (R + 1)
            loo_mass[pi_drop] = obs_mass_loo
            loo_LR[pi_drop] = obs_LR_loo
        loo_argmax = int(np.argmax(loo_p_mass))
        loo_max_p = float(loo_p_mass[loo_argmax])
        loo_argmax_patient = PATIENTS[loo_argmax]

        # --- Verdict mapping based on cluster_p_mass (mass-only gate,
        # locked 2026-05-19 per writing-agent feedback). cluster_p_LR
        # stays in the CSV as a descriptive co-statistic but does not
        # gate the verdict.
        if cluster_p_mass < 0.01:
            verdict = "strong"
        elif cluster_p_mass < 0.05:
            verdict = "weak"
        else:
            verdict = "no_trace"

        elapsed = time.time() - t0
        print(
            f"  obs LR: {obs_LR}, null mean LR: {null_LR.mean():.2f}, "
            f"null 95th LR: {np.percentile(null_LR, 95):.1f}, "
            f"cluster p (LR): {cluster_p_LR:.4f}, "
            f"cluster p (mass): {cluster_p_mass:.4f}, "
            f"LOO max p_mass: {loo_max_p:.4f} (drop {loo_argmax_patient}), "
            f"verdict: {verdict}, elapsed: {elapsed:.1f}s"
        )

        summary_rows.append({
            "band": band,
            "n_patients": P,
            "R": R,
            "k_grid_lo": K_GRID[0], "k_grid_hi": K_GRID[-1],
            "alpha_per_k": ALPHA_K,
            "obs_longest_run": obs_LR,
            "obs_cluster_mass_neglog10p": obs_mass,
            "null_mean_LR": float(null_LR.mean()),
            "null_median_LR": float(np.median(null_LR)),
            "null_p95_LR": float(np.percentile(null_LR, 95)),
            "null_max_LR": int(null_LR.max()),
            "null_mean_mass": float(null_mass.mean()),
            "null_p95_mass": float(np.percentile(null_mass, 95)),
            "cluster_p_longest_run": cluster_p_LR,
            "cluster_p_cluster_mass": cluster_p_mass,
            "cluster_p_mass_loo_max": loo_max_p,
            "cluster_p_mass_loo_argmax_patient": loo_argmax_patient,
            "verdict_cluster_extent": verdict,
        })

        # --- Per-patient LOO rows (one per patient dropped)
        for pi_drop in range(P):
            loo_rows.append({
                "band": band,
                "dropped_patient": PATIENTS[pi_drop],
                "obs_longest_run_loo": int(loo_LR[pi_drop]),
                "obs_cluster_mass_loo": float(loo_mass[pi_drop]),
                "cluster_p_mass_loo": float(loo_p_mass[pi_drop]),
            })

    pd.DataFrame(summary_rows).to_csv(OUT_DIR / "cohort_summary.csv", index=False)
    pd.DataFrame(null_rows).to_csv(OUT_DIR / "null_distribution.csv", index=False)
    pd.DataFrame(per_k_rows).to_csv(OUT_DIR / "per_k_obs_p.csv", index=False)
    pd.DataFrame(loo_rows).to_csv(OUT_DIR / "loo_cluster_p_mass.csv", index=False)

    readme = f"""---
name: grassmann-cluster-extent
era: IMCOH_ABS_COHORT_N10
status: complete
date: 2026-05-19
companion: ../grassmann_matched_strength_surrogate/
supersedes_in: .agents/preprint/locked/VERDICT_LEDGER.md (decision 6 — replaces 8-cell hardcoded threshold)
---

# Grassmann cluster-extent permutation null

Replaces the 8-cell contiguous-run threshold in `VERDICT_LEDGER.md` with a
data-derived gate. For each band, builds the empirical null distribution
of "longest contiguous-significant run" by re-running the matched-strength
cohort Wilcoxon `R={R}` times treating each surrogate as the phantom
observation.

## Verdict mapping (from `cluster_p_cluster_mass`, locked 2026-05-19 pm)
- `< 0.01` → strong
- `0.01 ≤ p < 0.05` → weak
- `≥ 0.05` → no trace

The mass-only gate replaces the earlier disjunctive
`min(cluster_p_LR, cluster_p_mass) < α` gate (per writing-agent
feedback 2026-05-19). The resilient all-clusters mass already
encodes both contiguity and depth; the parallel LR statistic stays
in the CSV as a descriptive co-statistic but does not gate the
verdict.

## Files
- `cohort_summary.csv` — per band: obs longest run, null distribution stats,
  cluster_p (longest-run-based and cluster-mass-based), LOO max p_mass,
  LOO argmax patient, verdict.
- `null_distribution.csv` — per band × surrogate phantom-test result.
- `per_k_obs_p.csv` — per band × k observed paired Wilcoxon p-value.
- `loo_cluster_p_mass.csv` — per band × dropped patient leave-one-out
  cluster_p_mass diagnostic (descriptive only, never a gate; per
  `feedback_no_single_patient_p_driven.md`).

## Cohort
{', '.join(PATIENTS)} (n={len(PATIENTS)}).

## Parameters
- `k_grid`: 2..112 ({len(K_GRID)} values, matches audit_66 K_GRID).
- `R`: {R} matched-strength surrogates (SWAP_FACTOR={SWAP_FACTOR}, seed={SEED}).
- `alpha_per_k`: {ALPHA_K} for the per-k paired Wilcoxon (one-sided less,
  trace direction T_G < 0).
"""
    (OUT_DIR / "README.md").write_text(readme)

    print(f"\nTotal elapsed: {(time.time() - overall_t0) / 60:.1f} min")
    print(f"Wrote: {OUT_DIR / 'cohort_summary.csv'}")
    print(f"Wrote: {OUT_DIR / 'null_distribution.csv'}")
    print(f"Wrote: {OUT_DIR / 'per_k_obs_p.csv'}")
    print(f"Wrote: {OUT_DIR / 'README.md'}")


if __name__ == "__main__":
    main()
