#!/usr/bin/env python
"""Matched-strength surrogate weight-distribution drift diagnostic.

Question: does the 4-cycle ±δ strength-preserving rewriting
(`lrg_eegfc.utils.surrogate.strength_preserving_shuffle`) leave the
edge-weight *distribution* approximately intact, or does it materially
drift the weight histogram from the observed FC?

The null is matched-strength (per-node strength exactly preserved). The
algorithm modifies individual edge weights by ±δ per swap, so the weight
distribution is only approximately preserved. This diagnostic measures
how much drift the SWAP_FACTOR=20 default produces on real β `imcoh_abs`
FCs across the n=10 cohort × 3 phases (rest_pre, task_test, rest_post).

For each (patient, phase) cell we generate R surrogates fresh, extract
upper-triangular edge weights from observed and from each surrogate,
then compare distributions via:

- Kolmogorov-Smirnov statistic `D` (max gap of empirical CDFs).
- Mean / std / 90th-percentile shifts (relative to observed).
- Pooled (cohort × phases) histogram overlay.

Outputs:
- ``data/audit/matched_strength_weight_drift/per_cell.csv`` — one row per
  (patient, phase) cell with KS + tail stats.
- ``data/audit/matched_strength_weight_drift/pooled_histogram.pdf`` —
  observed vs pooled-surrogate weight histogram.
- ``data/audit/matched_strength_weight_drift/per_cell_ks_distribution.pdf``
  — histogram of per-cell KS values.
- ``data/audit/matched_strength_weight_drift/README.md`` — auto-generated
  summary with decision-relevant numbers.

Decision rule:
- max KS < 0.05  → no material drift; matched-strength null is effectively
  weight-distribution-preserving. Stop with action (A).
- max KS in [0.05, 0.15] → mild drift; flag in methods but no action.
- max KS > 0.15 → material drift; justifies adding pure weight-permutation
  null as a complementary sanity check (action B).

Runtime: ~50 ms per surrogate × R × 10 patients × 3 phases ≈ 1–2 min for
R=50.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.surrogate import strength_preserving_shuffle, verify_strengths
from lrg_eegfc.workflow.fc import load_fc_matrix


OUT_ROOT = Path("data/audit/matched_strength_weight_drift")
PHASES = ("rest_pre", "task_test", "rest_post")
DEFAULT_BAND = "beta"
DEFAULT_R = 50
DEFAULT_SWAP_FACTOR = 20
DEFAULT_SEED = 20260528


def _triu_weights(W: np.ndarray) -> np.ndarray:
    """Upper-triangular off-diagonal weights of a symmetric matrix."""
    iu = np.triu_indices_from(W, k=1)
    return W[iu]


def _per_cell_stats(
    pat: str,
    phase: str,
    band: str,
    R: int,
    swap_factor: int,
    seed: int,
    fc_method: str,
) -> tuple[dict, np.ndarray, np.ndarray]:
    """Return (stats dict, observed triu weights, pooled-surrogate triu weights)."""
    W_obs = load_fc_matrix(pat, phase, band, fc_method)
    N = W_obs.shape[0]
    n_swaps = swap_factor * (N * (N - 1)) // 2

    obs_w = _triu_weights(W_obs)
    rng = np.random.default_rng(seed + hash((pat, phase, band)) % 2**31)
    surr_pool = []
    ks_per_surr = []
    n_ok = 0
    for r in range(R):
        W_s = strength_preserving_shuffle(W_obs, n_swaps, rng)
        if not verify_strengths(W_obs, W_s, tol=1e-4):
            continue
        surr_w = _triu_weights(W_s)
        surr_pool.append(surr_w)
        ks_per_surr.append(ks_2samp(obs_w, surr_w).statistic)
        n_ok += 1

    if n_ok == 0:
        return ({
            "patient": pat, "phase": phase, "band": band,
            "N": N, "R_ok": 0, "ks_pooled": np.nan,
            "ks_per_surr_mean": np.nan, "ks_per_surr_p95": np.nan,
            "mean_obs": np.nan, "mean_surr": np.nan, "mean_shift_rel": np.nan,
            "std_obs": np.nan, "std_surr": np.nan, "std_shift_rel": np.nan,
            "p90_obs": np.nan, "p90_surr": np.nan, "p90_shift_rel": np.nan,
        }, obs_w, np.array([]))

    surr_w_pool = np.concatenate(surr_pool)
    ks_per_surr = np.array(ks_per_surr)

    stats = {
        "patient": pat,
        "phase": phase,
        "band": band,
        "N": N,
        "R_ok": n_ok,
        "ks_pooled": ks_2samp(obs_w, surr_w_pool).statistic,
        "ks_per_surr_mean": float(ks_per_surr.mean()),
        "ks_per_surr_p95": float(np.percentile(ks_per_surr, 95)),
        "mean_obs": float(obs_w.mean()),
        "mean_surr": float(surr_w_pool.mean()),
        "mean_shift_rel": float((surr_w_pool.mean() - obs_w.mean()) / obs_w.mean()),
        "std_obs": float(obs_w.std()),
        "std_surr": float(surr_w_pool.std()),
        "std_shift_rel": float((surr_w_pool.std() - obs_w.std()) / obs_w.std()),
        "p90_obs": float(np.percentile(obs_w, 90)),
        "p90_surr": float(np.percentile(surr_w_pool, 90)),
        "p90_shift_rel": float(
            (np.percentile(surr_w_pool, 90) - np.percentile(obs_w, 90))
            / np.percentile(obs_w, 90)
        ),
    }
    return stats, obs_w, surr_w_pool


def _plot_pooled_histogram(
    obs_pool: np.ndarray, surr_pool: np.ndarray, out_path: Path
) -> None:
    import matplotlib.pyplot as plt

    from lrg_eegfc.visuals.styles import use_lrg_style
    use_lrg_style()

    fig, ax = plt.subplots(1, 1, figsize=(4.5, 3.2))
    bins = np.linspace(
        min(obs_pool.min(), surr_pool.min()),
        max(obs_pool.max(), surr_pool.max()),
        80,
    )
    ax.hist(obs_pool, bins=bins, density=True, alpha=0.55,
            color="#1f77b4", label="observed")
    ax.hist(surr_pool, bins=bins, density=True, alpha=0.55,
            color="#d62728", label="matched-strength surrogate")
    ax.set_xlabel(r"edge weight $|\mathrm{ImCoh}|$")
    ax.set_ylabel("density")
    ax.legend(frameon=False, loc="upper right")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _plot_per_cell_ks(per_cell: pd.DataFrame, out_path: Path) -> None:
    import matplotlib.pyplot as plt

    from lrg_eegfc.visuals.styles import use_lrg_style
    use_lrg_style()

    fig, ax = plt.subplots(1, 1, figsize=(4.5, 3.2))
    ax.hist(per_cell["ks_pooled"].dropna(), bins=15,
            color="#2ca02c", alpha=0.75)
    ax.axvline(0.05, color="grey", linestyle="--", linewidth=1,
               label="0.05 (no material drift)")
    ax.axvline(0.15, color="black", linestyle="--", linewidth=1,
               label="0.15 (material drift)")
    ax.set_xlabel("KS distance (observed vs pooled surrogate)")
    ax.set_ylabel("count of (patient, phase) cells")
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def _write_summary(per_cell: pd.DataFrame, band: str, R: int, swap_factor: int,
                   out_path: Path) -> None:
    ks = per_cell["ks_pooled"].dropna()
    mean_shift = per_cell["mean_shift_rel"].dropna()
    std_shift = per_cell["std_shift_rel"].dropna()
    p90_shift = per_cell["p90_shift_rel"].dropna()

    md = f"""---
name: matched-strength-weight-distribution-drift
era: IMCOH_ABS_COHORT_N10
status: current
kind: diagnostic-report
date: 2026-05-28
band: {band}
n_cells_total: {len(per_cell)}
n_cells_ok: {int((per_cell['R_ok'] > 0).sum())}
R_per_cell: {R}
swap_factor: {swap_factor}
fc_method: imcoh_abs
---

# Matched-strength weight-distribution drift on β

Question: does the 4-cycle ±δ strength-preserving rewiring leave the edge
weight *distribution* effectively intact across the n=10 cohort × 3
phases?

## Decision-relevant numbers

| Statistic | Median | p95 | max |
|---|---|---|---|
| KS distance (observed vs pooled surrogate) | {ks.median():.3f} | {np.percentile(ks, 95):.3f} | {ks.max():.3f} |
| relative mean shift | {mean_shift.median():+.3f} | {np.percentile(mean_shift, 95):+.3f} | {mean_shift.abs().max():.3f} |
| relative std shift | {std_shift.median():+.3f} | {np.percentile(std_shift, 95):+.3f} | {std_shift.abs().max():.3f} |
| relative 90th-percentile shift | {p90_shift.median():+.3f} | {np.percentile(p90_shift, 95):+.3f} | {p90_shift.abs().max():.3f} |

## Verdict

Decision thresholds:
- max KS < 0.05 → no material drift; matched-strength null is effectively
  weight-distribution-preserving.
- max KS in [0.05, 0.15] → mild drift; flag in methods but no action.
- max KS > 0.15 → material drift; justifies adding a pure
  weight-permutation null as a complementary sanity check.

Observed max KS = **{ks.max():.3f}** → {"NO MATERIAL DRIFT" if ks.max() < 0.05 else "MILD DRIFT" if ks.max() < 0.15 else "MATERIAL DRIFT"}.

## Per-cell table

See `per_cell.csv` (one row per (patient, phase) cell).

## Pooled histogram

See `pooled_histogram.pdf` — overlay of the observed and pooled-surrogate
weight distributions, aggregated across all cohort × phase cells.

## Per-cell KS distribution

See `per_cell_ks_distribution.pdf` — histogram of the {len(ks)} per-cell
KS values, with the 0.05 and 0.15 decision thresholds.

## Method

For each cell:
1. Load observed `|ImCoh|` FC adjacency from cache.
2. Generate R={R} matched-strength surrogates via
   `strength_preserving_shuffle` with `SWAP_FACTOR={swap_factor}`.
3. Extract upper-triangular off-diagonal weights from observed.
4. Pool upper-triangular off-diagonal weights from all R surrogates.
5. Compute KS statistic between observed and pooled-surrogate.

The algorithm preserves per-node strength exactly (within `1e-4`
tolerance, enforced by `verify_strengths`), so any drift in the weight
distribution is bounded by what the strength constraint permits.
"""
    out_path.write_text(md)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--band", default=DEFAULT_BAND)
    ap.add_argument("--R", type=int, default=DEFAULT_R)
    ap.add_argument("--swap-factor", type=int, default=DEFAULT_SWAP_FACTOR)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--fc-method", default="imcoh_abs")
    ap.add_argument("--out-root", type=Path, default=OUT_ROOT)
    args = ap.parse_args()

    rows = []
    obs_pool = []
    surr_pool = []
    for pat in PATIENTS_4PHASE:
        for phase in PHASES:
            print(f"[diag] {pat} {phase} ... ", end="", flush=True)
            try:
                stats, obs_w, surr_w = _per_cell_stats(
                    pat, phase, args.band, args.R, args.swap_factor,
                    args.seed, args.fc_method,
                )
            except Exception as e:
                print(f"FAIL ({e})")
                continue
            rows.append(stats)
            obs_pool.append(obs_w)
            if surr_w.size > 0:
                surr_pool.append(surr_w)
            print(f"KS = {stats['ks_pooled']:.3f}  "
                  f"Δmean = {stats['mean_shift_rel']:+.3f}  "
                  f"Δp90 = {stats['p90_shift_rel']:+.3f}  "
                  f"R_ok = {stats['R_ok']}/{args.R}")

    per_cell = pd.DataFrame(rows)
    args.out_root.mkdir(parents=True, exist_ok=True)
    per_cell.to_csv(args.out_root / "per_cell.csv", index=False)

    obs_arr = np.concatenate(obs_pool) if obs_pool else np.array([])
    surr_arr = np.concatenate(surr_pool) if surr_pool else np.array([])

    if obs_arr.size and surr_arr.size:
        _plot_pooled_histogram(obs_arr, surr_arr,
                               args.out_root / "pooled_histogram.pdf")
        _plot_per_cell_ks(per_cell,
                          args.out_root / "per_cell_ks_distribution.pdf")
        _write_summary(per_cell, args.band, args.R, args.swap_factor,
                       args.out_root / "README.md")

        pooled_ks = ks_2samp(obs_arr, surr_arr).statistic
        print(f"\n[diag] POOLED across cohort × phases: KS = {pooled_ks:.4f}")
        print(f"[diag] per-cell KS — median {per_cell['ks_pooled'].median():.3f}, "
              f"p95 {np.percentile(per_cell['ks_pooled'].dropna(), 95):.3f}, "
              f"max {per_cell['ks_pooled'].max():.3f}")
        print(f"[diag] wrote summary → {args.out_root / 'README.md'}")


if __name__ == "__main__":
    main()
