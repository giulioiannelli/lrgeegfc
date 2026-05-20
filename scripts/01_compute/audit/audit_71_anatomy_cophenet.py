#!/usr/bin/env python3
"""Audit 71 — cophenet trace anatomy under A1 + A3 (matched-strength) controls.

5-point critical preamble (per `feedback_critical_null_preamble.md`)
-------------------------------------------------------------------
1. Claim: top-decile per-pair Δ(ρ_split^coph) pairs at a given band
   concentrate disproportionately in specific DK regions vs the
   cohort-wide marginal distribution.

2. Null (A3, mandatory): R=200 strength-preserving 4-cycle ±δ surrogate
   adjacency matrices (already cached at
   ``data/cache/matched_strength_surrogate_lrg/Pat_NN/<band>_<phase>_*.npz``);
   for each surrogate r and phase, reconstruct ρ̂(τ_max), Trho, UPGMA,
   cophenet → D_coph_r; recompute per-pair Δ_task_r and Δ_rest_r;
   identify top-decile trace pairs per surrogate; aggregate region
   counts → empirical null distribution per region. Z-scores and
   p-empirical values per region.

3. Strongest plausible alternative: the trace pairs are anatomically
   distributed in proportion to baseline contact density per region
   (i.e., regions with more contacts have more trace pairs by chance
   without any genuine biological enrichment).

4. Does the null cover (3): yes for the strength-rearrangement
   component — matched-strength preserves node strength but reshuffles
   topology. The null thus tests whether the OBSERVED region-counts
   exceed what arises from rewiring under conserved strengths. What
   A3 does NOT directly control: cohort-level systematic clinical-
   implant bias (e.g., hippocampus densely sampled across patients)
   — addressed separately by A2 sampling-corrected bootstrap (NOT
   implemented here; this audit is A1 + A3 only; A2/A4 deferred).

5. Falsification: a region with observed enrichment ≥ surrogate p95
   AND z > 2 is flagged. If no region passes, anatomy is `not localized`
   for this (band, probe) cell.

Inputs
------
- ``data/reports/imcoh_continuous_trace/per_pair_split/Pat_NN_<band>.npz``
  (observed per-pair Δ_task/Δ_rest/iu_i/iu_j arrays from audit_63 era)
- ``data/cache/matched_strength_surrogate_lrg/Pat_NN/<band>_<phase>_R200_*.npz``
  (surrogate eigvals + eigvecs for 4 phases: rest_pre_A, rest_pre_B,
  task_test, rest_post)
- DK region labels via ``lrg_eegfc.utils.io.regions.load_channel_regions``

Outputs (under ``data/audit/anatomy_<band>_cophenet/``)
------------------------------------------------------
- ``cohort_summary.csv`` — per-region observed enrichment + A3 stats
- ``per_patient_trace_flags.csv`` — per-patient trace pair counts
- ``README.md`` — audit provenance

Usage
-----
    python audit_71_anatomy_cophenet.py --band beta
    python audit_71_anatomy_cophenet.py --band alpha --n-surrogates 50  # smoke test
    python audit_71_anatomy_cophenet.py --band alpha --epi-x            # C5 epi-X subset

Reads only cached artifacts; no FC recomputation.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.spatial.distance import squareform
from scipy.stats import hypergeom

from lrg_eegfc.utils.io.patient import load_epileptic_nodes
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
PHASES_SPLIT = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
TOP_DECILE = 0.10
SEED_TAG = "R200_swap20_seed20260511_imcoh_abs"

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def load_obs_pair_shifts(pat: str, band: str) -> dict:
    """Load observed per-pair Δ_task, Δ_rest, pair-endpoint indices."""
    p = ROOT / "data/reports/imcoh_continuous_trace/per_pair_split" / f"{pat}_{band}.npz"
    d = np.load(p)
    return {
        "dD_task": np.asarray(d["dD_task"]),
        "dD_rest": np.asarray(d["dD_rest"]),
        "iu_i": np.asarray(d["iu_i"]).astype(int),
        "iu_j": np.asarray(d["iu_j"]).astype(int),
    }


def trace_flag(dD_task: np.ndarray, dD_rest: np.ndarray,
               top_decile: float = TOP_DECILE) -> np.ndarray:
    """Top-decile pairs with consistent sign (dD_task · dD_rest > 0)."""
    consistent = np.sign(dD_task) == np.sign(dD_rest)
    score = np.abs(dD_task) * np.abs(dD_rest)
    score_masked = np.where(consistent, score, -np.inf)
    valid = np.isfinite(score_masked)
    if not valid.any():
        return np.zeros_like(dD_task, dtype=bool)
    threshold = np.quantile(score_masked[valid], 1.0 - top_decile)
    return score_masked >= threshold


def cophenet_from_eigs(eigvals: np.ndarray, eigvecs: np.ndarray) -> np.ndarray:
    """Reconstruct D_coph from cached eigvals/eigvecs of one surrogate.

    Steps: τ_max = 1/λ_max ⇒ ρ̂(τ_max) = V diag(exp(-τ_max λ)) V^T / Σ_k exp(-τ_max λ_k)
    ⇒ Trho = 1/ρ̂ (diag=0) ⇒ UPGMA(squareform(Trho)) ⇒ cophenet → condensed D_coph.

    Returns condensed (1D) cophenet distance array, matching the
    ``ultrametric_matrix`` field shape in cached LRG results.
    """
    # Drop the trivial zero eigenvalue (connected component); use the
    # largest eigenvalue for tau_max.
    lam = np.asarray(eigvals, dtype=np.float64)
    V = np.asarray(eigvecs, dtype=np.float64)
    lam_max = float(lam[-1])
    tau = 1.0 / lam_max
    w = np.exp(-tau * lam)
    Z = float(w.sum())
    # rho = V diag(w) V^T / Z
    rho = (V * w) @ V.T / Z
    # Numerical floor: rho can have tiny negative artefacts off-diagonal
    # from finite eigh precision. Clip to a small positive value before
    # inverting.
    rho = np.maximum(rho, 1e-12)
    Trho = 1.0 / rho
    Trho = np.maximum(Trho, Trho.T)
    np.fill_diagonal(Trho, 0.0)
    dists = squareform(Trho, checks=False)
    Zlnk = linkage(dists, method="average")
    coph = cophenet(Zlnk)
    return np.asarray(coph, dtype=np.float64)


def surrogate_pair_shifts(pat: str, band: str, n_surrogates: int = 200) -> dict:
    """Compute per-pair Δ_task_r, Δ_rest_r per surrogate r.

    Loads cached surrogate eigvals/eigvecs for 4 phases (rest_pre_A,
    rest_pre_B, task_test, rest_post), reconstructs D_coph per
    surrogate per phase, computes per-pair shifts.

    Returns
    -------
    dict with keys:
        dD_task_surr: array (R, n_pairs)
        dD_rest_surr: array (R, n_pairs)
    """
    cache_root = ROOT / "data/cache/matched_strength_surrogate_lrg" / pat
    files = {phase: cache_root / f"{band}_{phase}_{SEED_TAG}.npz" for phase in PHASES_SPLIT}
    for phase, fp in files.items():
        if not fp.exists():
            raise FileNotFoundError(f"Missing surrogate cache: {fp}")

    # Load all 4 phases' surrogate eigs
    surr_eigs = {phase: np.load(fp) for phase, fp in files.items()}
    R_avail = surr_eigs["rest_pre_A"]["eigvals"].shape[0]
    R = min(n_surrogates, R_avail)

    # Build D_coph per phase per surrogate; aggregate into per-pair shift arrays.
    # Each cophenet is condensed shape (N*(N-1)/2,).
    n_pairs_check = None
    coph = {phase: [] for phase in PHASES_SPLIT}
    for r in range(R):
        for phase in PHASES_SPLIT:
            ev = surr_eigs[phase]["eigvals"][r]
            vv = surr_eigs[phase]["eigvecs"][r]
            c = cophenet_from_eigs(ev, vv)
            if n_pairs_check is None:
                n_pairs_check = c.size
            coph[phase].append(c)
    coph_arr = {phase: np.stack(coph[phase], axis=0) for phase in PHASES_SPLIT}

    dD_task_surr = coph_arr["task_test"] - coph_arr["rest_pre_A"]
    dD_rest_surr = coph_arr["rest_post"] - coph_arr["rest_pre_B"]

    return {
        "dD_task_surr": dD_task_surr,
        "dD_rest_surr": dD_rest_surr,
        "R": R,
    }


def per_pair_region_pair(regions: pd.Series, iu_i: np.ndarray,
                         iu_j: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return arrays of DK region labels for both endpoints of each pair."""
    r = np.asarray(regions.values)
    return r[iu_i], r[iu_j]


def count_region_endpoints(region_i: np.ndarray, region_j: np.ndarray,
                           mask: np.ndarray) -> dict[str, int]:
    """Count contributions of each DK region to flagged pairs.

    A pair contributes its endpoints' DK labels to the region count
    (each pair → 2 contributions). Same-region pairs contribute twice
    to that region.
    """
    counts: dict[str, int] = {}
    for r in region_i[mask]:
        counts[r] = counts.get(r, 0) + 1
    for r in region_j[mask]:
        counts[r] = counts.get(r, 0) + 1
    return counts


def hypergeom_per_region(trace_counts: dict[str, int], all_counts: dict[str, int],
                         trace_total: int, all_total: int) -> pd.DataFrame:
    """A1 hypergeometric one-sided greater per DK region.

    Universe = all_total (counted as endpoint-instances, so each pair
    contributes 2). Trace = trace_total (top-decile endpoints). For
    each region r: K_r = all_counts[r], k_r = trace_counts[r].
    H0: k_r ~ hypergeom(N=all_total, K=K_r, n=trace_total).
    """
    rows = []
    for region, K in all_counts.items():
        k = trace_counts.get(region, 0)
        # one-sided right-tail: P(X >= k) = 1 - cdf(k-1)
        p = float(hypergeom.sf(k - 1, all_total, K, trace_total))
        expected = float(trace_total * K / all_total) if all_total > 0 else 0.0
        enrichment = float(k / expected) if expected > 0 else float("nan")
        rows.append({
            "region": region,
            "K_region": int(K),
            "k_trace": int(k),
            "expected": expected,
            "enrichment": enrichment,
            "p_hyper_one_sided": p,
        })
    df = pd.DataFrame(rows).sort_values("p_hyper_one_sided")
    # BH-FDR within band/probe
    p = df["p_hyper_one_sided"].values
    m = len(p)
    sorted_idx = np.argsort(p)
    sorted_p = p[sorted_idx]
    q_sorted = sorted_p * m / (np.arange(m) + 1)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q = np.empty(m)
    q[sorted_idx] = q_sorted
    df["q_bh"] = q
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# main per-band entry
# ---------------------------------------------------------------------------

def run_band(band: str, n_surrogates: int, epi_x: bool, out_dir: Path,
             verbose: bool = True) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # 1. Collect observed per-patient trace flags and per-pair region tags
    obs_records = []
    obs_region_counts: dict[str, int] = {}
    obs_trace_total = 0
    all_region_counts: dict[str, int] = {}
    all_total = 0
    obs_traces_by_pat: dict[str, np.ndarray] = {}
    pair_regions_by_pat: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    for pat in COHORT:
        obs = load_obs_pair_shifts(pat, band)
        regions_df = load_channel_regions(pat)
        if epi_x:
            epi = load_epileptic_nodes(pat)
            keep_mask = ~np.isin(np.arange(len(regions_df)), np.asarray(list(epi)))
            pair_keep = keep_mask[obs["iu_i"]] & keep_mask[obs["iu_j"]]
        else:
            pair_keep = np.ones_like(obs["iu_i"], dtype=bool)

        dD_task = obs["dD_task"][pair_keep]
        dD_rest = obs["dD_rest"][pair_keep]
        iu_i = obs["iu_i"][pair_keep]
        iu_j = obs["iu_j"][pair_keep]
        flag = trace_flag(dD_task, dD_rest)
        region_i, region_j = per_pair_region_pair(regions_df["region"], iu_i, iu_j)

        pat_counts = count_region_endpoints(region_i, region_j, flag)
        all_counts_pat = count_region_endpoints(region_i, region_j,
                                                np.ones_like(flag, dtype=bool))
        for r, c in pat_counts.items():
            obs_region_counts[r] = obs_region_counts.get(r, 0) + c
        for r, c in all_counts_pat.items():
            all_region_counts[r] = all_region_counts.get(r, 0) + c

        obs_trace_total += 2 * int(flag.sum())
        all_total += 2 * len(flag)
        obs_records.append({
            "patient": pat,
            "n_pairs": int(len(flag)),
            "n_trace_pairs": int(flag.sum()),
            "epi_x": bool(epi_x),
        })
        obs_traces_by_pat[pat] = flag
        pair_regions_by_pat[pat] = (region_i, region_j)
        if verbose:
            print(f"  {pat} {band}: {int(flag.sum())}/{len(flag)} trace pairs (epi_x={epi_x})")

    obs_records_df = pd.DataFrame(obs_records)
    obs_records_df.to_csv(out_dir / "per_patient_trace_flags.csv", index=False)

    # 2. A1 hypergeometric on observed trace
    a1 = hypergeom_per_region(obs_region_counts, all_region_counts,
                              obs_trace_total, all_total)
    if verbose:
        print(f"  A1 hypergeometric: {len(a1)} regions, "
              f"{int((a1.q_bh < 0.05).sum())} regions with q_BH < 0.05")

    # 3. A3 matched-strength surrogate enrichment per region
    # For each surrogate r ∈ {1..R}:
    #   For each patient: compute trace flag on (dD_task_r, dD_rest_r)
    #   Aggregate region counts cohort-wide
    # → empirical null distribution per region (R values)
    surr_counts_by_region: dict[str, list[int]] = {r: [] for r in all_region_counts}
    R_used = n_surrogates
    for ri in range(R_used):
        cohort_surr_counts: dict[str, int] = {}
        for pat in COHORT:
            if ri == 0:
                # Lazy load on first iteration
                pass
            # Compute per-pair surrogate shifts on demand
            pass

    # Loop differently: for each patient, compute ALL surrogates' Δ, then
    # tabulate per-region counts per surrogate.
    if verbose:
        print(f"  A3: computing surrogate enrichments for R={R_used} per patient")
    null_distributions = {r: np.zeros(R_used, dtype=int) for r in all_region_counts}
    for pat in COHORT:
        t_pat = time.time()
        surr = surrogate_pair_shifts(pat, band, n_surrogates=R_used)
        dD_task_surr = surr["dD_task_surr"]
        dD_rest_surr = surr["dD_rest_surr"]
        R = surr["R"]
        # Need to mask same way as obs for epi_x (the iu_i/iu_j ordering
        # matches obs since both come from squareform of an N×N matrix).
        obs = load_obs_pair_shifts(pat, band)
        regions_df = load_channel_regions(pat)
        if epi_x:
            epi = load_epileptic_nodes(pat)
            keep_mask = ~np.isin(np.arange(len(regions_df)), np.asarray(list(epi)))
            pair_keep = keep_mask[obs["iu_i"]] & keep_mask[obs["iu_j"]]
        else:
            pair_keep = np.ones_like(obs["iu_i"], dtype=bool)
        iu_i = obs["iu_i"][pair_keep]
        iu_j = obs["iu_j"][pair_keep]
        region_i, region_j = per_pair_region_pair(regions_df["region"], iu_i, iu_j)
        for ri in range(R):
            flag_r = trace_flag(dD_task_surr[ri][pair_keep], dD_rest_surr[ri][pair_keep])
            counts_r = count_region_endpoints(region_i, region_j, flag_r)
            for r, c in counts_r.items():
                null_distributions[r][ri] += c
        if verbose:
            print(f"    {pat} surrogate enrichments done in {time.time() - t_pat:.1f}s")

    # 4. Build per-region A3 summary
    a3_rows = []
    for region, K in all_region_counts.items():
        obs_k = obs_region_counts.get(region, 0)
        null = null_distributions[region]
        z = (obs_k - null.mean()) / max(null.std(ddof=0), 1e-9)
        p_emp = (1 + int((null >= obs_k).sum())) / (R_used + 1)
        a3_rows.append({
            "region": region,
            "obs_k": int(obs_k),
            "null_mean": float(null.mean()),
            "null_std": float(null.std(ddof=0)),
            "null_p95": float(np.quantile(null, 0.95)),
            "obs_z": float(z),
            "p_empirical": float(p_emp),
        })
    a3 = pd.DataFrame(a3_rows)

    # Combined cohort_summary: A1 + A3 per region
    cohort_summary = a1.merge(a3, on="region", how="outer")
    cohort_summary["passes_A1"] = cohort_summary["q_bh"] < 0.05
    cohort_summary["passes_A3"] = (cohort_summary["p_empirical"] < 0.05) & (cohort_summary["obs_z"] > 2.0)
    cohort_summary["passes_both"] = cohort_summary["passes_A1"] & cohort_summary["passes_A3"]
    cohort_summary = cohort_summary.sort_values("p_empirical")
    cohort_summary.to_csv(out_dir / "cohort_summary.csv", index=False)

    # 5. Save null distributions
    np.savez_compressed(
        out_dir / "A3_null_distributions.npz",
        **{r: null_distributions[r] for r in null_distributions},
        regions=np.array(list(null_distributions.keys())),
    )

    # 6. README
    n_a1_pass = int(cohort_summary.passes_A1.fillna(False).sum())
    n_a3_pass = int(cohort_summary.passes_A3.fillna(False).sum())
    n_both = int(cohort_summary.passes_both.fillna(False).sum())
    readme = f"""# Anatomy audit (audit_71) — {band} cophenet probe{' (C5 epi-X)' if epi_x else ''}

Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}
Compute time: {time.time() - t0:.1f} s
Cohort: {len(COHORT)} patients
Surrogates: R={R_used} (matched-strength, seed 20260511)
Top decile cutoff: {TOP_DECILE}

## Verdict counts
- Regions with A1 q_BH < 0.05: {n_a1_pass}
- Regions with A3 (p_emp < 0.05 AND z > 2): {n_a3_pass}
- Regions passing both (verdict floor): {n_both}

## Files
- cohort_summary.csv: per-region A1+A3 statistics + pass flags
- per_patient_trace_flags.csv: per-patient trace pair counts
- A3_null_distributions.npz: per-region surrogate null distributions (R={R_used} values per region)

## A1 + A3 — what they test
A1 = hypergeometric one-sided greater per DK region, BH-FDR within band.
     Tests whether the OBSERVED region endpoint count exceeds chance
     given the cohort-wide marginal contact distribution.
A3 = matched-strength surrogate null (R={R_used}). Tests whether the
     observed region endpoint count exceeds what arises from
     strength-preserving 4-cycle ±δ rewiring of the FC adjacency.

## What's NOT in this audit
- A2 (sampling-corrected bootstrap): deferred. Would correct for
  per-patient implant coverage bias not addressed by A3.
- A4 (implant-geometry regression): deferred. Would test whether
  geometric covariates (epi-distance, implant centroid distance,
  hemisphere) independently predict trace-flag at q_BH < 0.05.
"""
    (out_dir / "README.md").write_text(readme)

    if verbose:
        print(f"\n  Wrote {out_dir}/cohort_summary.csv")
        print(f"  A1 pass: {n_a1_pass}/{len(cohort_summary)}, "
              f"A3 pass: {n_a3_pass}/{len(cohort_summary)}, "
              f"both: {n_both}/{len(cohort_summary)}")
        print(f"  Total time: {time.time() - t0:.1f} s")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--band", required=True,
                        choices=["delta", "theta", "alpha", "beta",
                                 "low_gamma", "high_gamma"])
    parser.add_argument("--n-surrogates", type=int, default=200,
                        help="Number of surrogates (default 200; use 20 for smoke test)")
    parser.add_argument("--epi-x", action="store_true",
                        help="Run on C5 epi-X subset (drops epi-zone contacts)")
    args = parser.parse_args()

    epi_tag = "_epiX" if args.epi_x else ""
    out_dir = ROOT / "data" / "audit" / f"anatomy_{args.band}_cophenet{epi_tag}"
    run_band(args.band, args.n_surrogates, args.epi_x, out_dir)


if __name__ == "__main__":
    main()
