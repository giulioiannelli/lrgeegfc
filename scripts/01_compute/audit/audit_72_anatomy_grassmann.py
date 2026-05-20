#!/usr/bin/env python3
"""Audit 72 — Grassmann trace anatomy under A1 + A3 controls.

5-point critical preamble
-------------------------
1. Claim: top-decile per-node participation `p_i(k)` in the leading-`k`
   Laplacian eigenmode subspace `U_k = span{φ_2..φ_{k+1}}` (averaged
   across the band's contiguous-significant `k`-window AND across the
   four phases rest_pre_A / rest_pre_B / task_test / rest_post and the
   cohort) concentrates disproportionately in specific DK regions.

   Operational definition: per-node participation at a single `k` is
   `p_i(k) = Σ_{j∈U_k} φ_j(i)²` (probability the i-th node carries
   weight in the leading subspace). The TRACE-FLAGGED participation
   is **per-patient ranked**: a node is flagged if its mean
   participation over (k ∈ window) × (phases) ranks in the top decile
   of all nodes for that patient.

2. Null (A3, mandatory): R=200 matched-strength surrogate eigvecs
   already cached. For each surrogate r, compute the same per-node
   participation, identify top-decile nodes, aggregate region counts
   → empirical null per region.

3. Strongest plausible alternative: nodes with high per-node degree
   trivially concentrate in the leading subspace (Laplacian
   eigenmodes weighted toward high-strength nodes). Region-level
   enrichment then reflects implant coverage of high-strength nodes,
   not genuine subspace rotation between phases.

4. Does the null cover (3): YES — matched-strength surrogacy
   preserves node strengths to 10⁻⁴, so the strength-dependence of
   participation is preserved on average. A region that passes A3
   has participation above what strength-conserved rewiring produces.

5. Falsification: per-region p_emp ≥ 0.05 OR z ≤ 2 → not anatomy-
   localized at this band's Grassmann probe.

Inputs
------
- `data/cache/imcoh_lrg/Pat_NN/<band>_<phase>_lrg_imcoh-abs.npz`
  (observed eigvecs per phase)
- `data/cache/matched_strength_surrogate_lrg/Pat_NN/<band>_<phase>_R200_*.npz`
  (surrogate eigvecs per phase × R=200)
- DK region labels via load_channel_regions

Outputs (under `data/audit/anatomy_<band>_grassmann[_epiX]/`)
------------------------------------------------------------
- cohort_summary.csv — per-region A1+A3 + pass flags
- per_patient_trace_flags.csv — per-patient node flag counts
- README.md — provenance

Usage
-----
    python audit_72_anatomy_grassmann.py --band beta
    python audit_72_anatomy_grassmann.py --band delta --epi-x

Per-band k-window comes from `VERDICT_LEDGER.md`:
    β: k=27..55, α: no Grassmann trace (skip), γ_l: k=12..23,
    γ_h: k=19..27 (no-trace but interpretable), θ: no Grassmann window,
    δ: k=57..63 (full-cohort) OR k=33..39 (epi-X).
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from lrg_eegfc.utils.io.patient import load_epileptic_nodes
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
# Observed LRG cache phases (no split-baseline halves on the obs side)
OBS_PHASES = ("rest_pre", "task_test", "rest_post")
# Surrogate cache phases (split-baseline halves)
SURR_PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
TOP_DECILE = 0.10
SEED_TAG = "R200_swap20_seed20260511_imcoh_abs"

# k-windows per band (default; --k-lo/--k-hi can override)
BAND_K_WINDOWS = {
    "beta":       (27, 55),    # 29-cell, strong
    "low_gamma":  (12, 23),    # 13-cell, weak
    "delta":      (57, 63),    # 7-cell, weak (full-cohort)
    "high_gamma": (19, 27),    # 9-cell, no trace at gate but interpretable
    "alpha":      (11, 14),    # 4-cell, no trace; here for completeness
    "theta":      (10, 14),    # null; cohort-min-p region as a sham
}
# For epi-X overrides at δ: k=33..39
BAND_K_WINDOWS_EPIX = {
    "delta":      (33, 39),
}


def participation_kset(eigvecs: np.ndarray, k_iter) -> np.ndarray:
    """Per-node participation in U_k averaged over an arbitrary `k` set.

    Eigvecs are stored ascending; we skip the trivial zero mode (col 0).
    For each k, U_k = columns 1..k+1 (i.e., φ_2..φ_{k+1}).
    Per-node participation at k = Σ_{j∈U_k} φ_j(i)².
    Average over the supplied iterable of `k` values (can be a contiguous
    range or the disjoint p<0.05 cluster-extent set).
    """
    N = eigvecs.shape[0]
    weights = np.zeros(N, dtype=np.float64)
    n_k = 0
    for k in k_iter:
        if k + 1 > N:
            continue
        U = eigvecs[:, 1:k + 1]
        weights += (U * U).sum(axis=1)
        n_k += 1
    return weights / max(n_k, 1)


def participation_phase_avg(eig_per_phase: dict[str, np.ndarray],
                            k_iter,
                            phases: tuple = OBS_PHASES) -> np.ndarray:
    """Average per-node participation across phases over a `k` set."""
    k_list = list(k_iter)
    parts = [participation_kset(eig_per_phase[ph], k_list) for ph in phases]
    return np.mean(np.stack(parts, axis=0), axis=0)


def load_obs_eigvecs(pat: str, band: str) -> dict[str, np.ndarray]:
    """Load observed eigvecs for OBS_PHASES from standard imcoh_lrg cache."""
    out = {}
    cache_root = ROOT / "data" / "cache" / "imcoh_lrg"
    for ph in OBS_PHASES:
        path = cache_root / pat / f"{band}_{ph}_lrg_imcoh-abs.npz"
        if not path.exists():
            raise FileNotFoundError(f"Missing obs eigvecs: {path}")
        d = np.load(path)
        out[ph] = np.asarray(d["eigenvectors"], dtype=np.float64)
    return out


def load_surr_eigvecs(pat: str, band: str, n_surrogates: int = 200) -> dict[str, np.ndarray]:
    """Load surrogate eigvecs for SURR_PHASES; trim to common R."""
    cache_root = ROOT / "data/cache/matched_strength_surrogate_lrg" / pat
    out = {}
    R_min = n_surrogates
    for ph in SURR_PHASES:
        path = cache_root / f"{band}_{ph}_{SEED_TAG}.npz"
        if not path.exists():
            raise FileNotFoundError(f"Missing surrogate eigvecs: {path}")
        d = np.load(path)
        ev = np.asarray(d["eigvecs"], dtype=np.float64)  # (R, N, N)
        R_min = min(R_min, ev.shape[0])
        out[ph] = ev
    for ph in SURR_PHASES:
        out[ph] = out[ph][:R_min]
    return out


def trace_flag_by_participation(participation: np.ndarray,
                                top_decile: float = TOP_DECILE) -> np.ndarray:
    """Top-decile nodes by participation (within patient)."""
    thr = np.quantile(participation, 1.0 - top_decile)
    return participation >= thr


def count_region_nodes(regions: np.ndarray, mask: np.ndarray) -> dict[str, int]:
    counts: dict[str, int] = {}
    for r in regions[mask]:
        counts[r] = counts.get(r, 0) + 1
    return counts


def hypergeom_per_region(trace_counts, all_counts, trace_total, all_total) -> pd.DataFrame:
    rows = []
    for region, K in all_counts.items():
        k = trace_counts.get(region, 0)
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


def run_band(band: str, k_iter, n_surrogates: int, epi_x: bool,
             out_dir: Path, verbose: bool = True, scope_tag: str = "") -> None:
    k_list = sorted(set(int(k) for k in k_iter))
    if not k_list:
        raise ValueError("k_iter is empty")
    k_lo, k_hi = k_list[0], k_list[-1]
    n_k = len(k_list)
    contiguous = (n_k == (k_hi - k_lo + 1))
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # 1. Observed enrichment
    obs_records = []
    obs_region_counts: dict[str, int] = {}
    all_region_counts: dict[str, int] = {}
    obs_trace_total = 0
    all_total = 0
    pat_state: dict[str, dict] = {}

    for pat in COHORT:
        regions_df = load_channel_regions(pat)
        regions = regions_df["region"].values
        eig_obs = load_obs_eigvecs(pat, band)
        N = eig_obs[OBS_PHASES[0]].shape[0]

        if epi_x:
            epi = load_epileptic_nodes(pat)
            keep_mask = ~np.isin(np.arange(N), np.asarray(list(epi)))
        else:
            keep_mask = np.ones(N, dtype=bool)

        # Per-patient observed participation, masked to non-epi nodes if epi_x.
        # Computing in full-graph spectral space then masking the node set —
        # this is the operational definition used elsewhere in audit_67.
        part_obs = participation_phase_avg(eig_obs, k_list, phases=OBS_PHASES)
        part_obs_kept = part_obs[keep_mask]
        regions_kept = regions[keep_mask]

        flag = trace_flag_by_participation(part_obs_kept)
        trace_regions_obs = regions_kept[flag]
        all_regions_kept = regions_kept

        pat_counts = count_region_nodes(regions_kept, flag)
        all_counts_pat = count_region_nodes(regions_kept, np.ones_like(flag))
        for r, c in pat_counts.items():
            obs_region_counts[r] = obs_region_counts.get(r, 0) + c
        for r, c in all_counts_pat.items():
            all_region_counts[r] = all_region_counts.get(r, 0) + c

        obs_trace_total += int(flag.sum())
        all_total += len(flag)
        obs_records.append({
            "patient": pat,
            "n_nodes_kept": int(keep_mask.sum()),
            "n_trace_nodes": int(flag.sum()),
            "epi_x": bool(epi_x),
        })
        pat_state[pat] = {
            "regions_kept": regions_kept,
            "keep_mask": keep_mask,
        }
        if verbose:
            print(f"  {pat} {band}: {int(flag.sum())}/{int(keep_mask.sum())} trace nodes")

    pd.DataFrame(obs_records).to_csv(out_dir / "per_patient_trace_flags.csv", index=False)

    a1 = hypergeom_per_region(obs_region_counts, all_region_counts,
                              obs_trace_total, all_total)
    n_a1 = int((a1.q_bh < 0.05).sum())
    if verbose:
        print(f"  A1: {n_a1}/{len(a1)} regions with q_BH < 0.05")

    # 2. A3 — surrogate enrichment
    null_distributions = {r: np.zeros(n_surrogates, dtype=int) for r in all_region_counts}
    R_used = n_surrogates
    for pat in COHORT:
        t_pat = time.time()
        eig_surr = load_surr_eigvecs(pat, band, n_surrogates=n_surrogates)
        R = min(eig_surr[SURR_PHASES[0]].shape[0], R_used)
        regions_kept = pat_state[pat]["regions_kept"]
        keep_mask = pat_state[pat]["keep_mask"]
        for ri in range(R):
            eig_r = {ph: eig_surr[ph][ri] for ph in SURR_PHASES}
            part_r = participation_phase_avg(eig_r, k_list, phases=SURR_PHASES)
            part_r_kept = part_r[keep_mask]
            flag_r = trace_flag_by_participation(part_r_kept)
            counts_r = count_region_nodes(regions_kept, flag_r)
            for r, c in counts_r.items():
                if r in null_distributions:
                    null_distributions[r][ri] += c
        if verbose:
            print(f"    {pat} surrogate enrichments in {time.time() - t_pat:.1f}s")

    # 3. A3 summary
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

    cohort_summary = a1.merge(a3, on="region", how="outer")
    cohort_summary["passes_A1"] = cohort_summary["q_bh"] < 0.05
    cohort_summary["passes_A3"] = (cohort_summary["p_empirical"] < 0.05) & (cohort_summary["obs_z"] > 2.0)
    cohort_summary["passes_both"] = cohort_summary["passes_A1"] & cohort_summary["passes_A3"]
    cohort_summary = cohort_summary.sort_values("p_empirical")
    cohort_summary.to_csv(out_dir / "cohort_summary.csv", index=False)

    np.savez_compressed(
        out_dir / "A3_null_distributions.npz",
        **{r: null_distributions[r] for r in null_distributions},
        regions=np.array(list(null_distributions.keys())),
    )

    n_a1_pass = int(cohort_summary.passes_A1.fillna(False).sum())
    n_a3_pass = int(cohort_summary.passes_A3.fillna(False).sum())
    n_both = int(cohort_summary.passes_both.fillna(False).sum())
    scope_descr = (f"contiguous k=[{k_lo},{k_hi}] ({n_k} cells)" if contiguous
                   else f"cluster-extent set ({n_k} cells, span k=[{k_lo},{k_hi}])")
    readme = f"""# Anatomy audit (audit_72) — {band} Grassmann probe{' (C5 epi-X)' if epi_x else ''}

Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}
k-scope: {scope_descr}  ({scope_tag or 'default'})
k-cells: {k_list}
Compute time: {time.time() - t0:.1f} s
Cohort: {len(COHORT)} patients
Surrogates: R={R_used}
Top decile cutoff (per patient): {TOP_DECILE}

## Verdict counts
- A1 q_BH < 0.05: {n_a1_pass}/{len(cohort_summary)}
- A3 (p_emp < 0.05 AND z > 2): {n_a3_pass}/{len(cohort_summary)}
- both: {n_both}/{len(cohort_summary)}
"""
    (out_dir / "README.md").write_text(readme)
    if verbose:
        print(f"  Wrote {out_dir}/cohort_summary.csv")
        print(f"  A1={n_a1_pass}/{len(cohort_summary)}, A3={n_a3_pass}, both={n_both}")
        print(f"  Total: {time.time() - t0:.1f} s")


def _load_cluster_extent_k_set(band: str) -> list[int]:
    """All p<0.05 k cells from `grassmann_cluster_extent/per_k_obs_p.csv`.

    This is the disjunctive `cluster-mass` set: every k where the cohort
    one-sided Wilcoxon p < 0.05, regardless of contiguity. Used when the
    user asks for the cluster-extent anatomy scope rather than the
    longest-contiguous-run window.
    """
    csv = ROOT / "data" / "audit" / "grassmann_cluster_extent" / "per_k_obs_p.csv"
    df = pd.read_csv(csv)
    sub = df[(df["band"] == band) & (df["obs_p_one_sided_less"] < 0.05)]
    return sorted(int(k) for k in sub["k"].values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--band", required=True,
                        choices=["delta", "theta", "alpha", "beta",
                                 "low_gamma", "high_gamma"])
    parser.add_argument("--n-surrogates", type=int, default=200)
    parser.add_argument("--epi-x", action="store_true")
    parser.add_argument("--k-lo", type=int, default=None,
                        help="Override k-window low (default from band)")
    parser.add_argument("--k-hi", type=int, default=None,
                        help="Override k-window high (default from band)")
    parser.add_argument("--k-list", type=str, default=None,
                        help="Comma-separated explicit list of k values "
                             "(overrides --k-lo/--k-hi and disables the "
                             "contiguous-window default)")
    parser.add_argument("--cluster-extent", action="store_true",
                        help="Use all p<0.05 cells from grassmann cluster-"
                             "extent per_k_obs_p.csv (cluster-mass set) "
                             "instead of the longest-contiguous-run window")
    parser.add_argument("--out-suffix", type=str, default="",
                        help="Suffix appended to the audit output directory")
    args = parser.parse_args()

    scope_tag = ""
    if args.cluster_extent:
        k_list = _load_cluster_extent_k_set(args.band)
        if not k_list:
            raise SystemExit(f"No p<0.05 k cells found for band {args.band}")
        scope_tag = f"cluster-extent ({len(k_list)} cells)"
    elif args.k_list:
        k_list = [int(x) for x in args.k_list.split(",") if x.strip()]
        scope_tag = f"explicit k-list ({len(k_list)} cells)"
    else:
        if args.epi_x and args.band in BAND_K_WINDOWS_EPIX:
            default_lo, default_hi = BAND_K_WINDOWS_EPIX[args.band]
        else:
            default_lo, default_hi = BAND_K_WINDOWS[args.band]
        k_lo = args.k_lo if args.k_lo is not None else default_lo
        k_hi = args.k_hi if args.k_hi is not None else default_hi
        k_list = list(range(k_lo, k_hi + 1))
        scope_tag = f"contiguous [{k_lo},{k_hi}]"

    epi_tag = "_epiX" if args.epi_x else ""
    suffix = args.out_suffix or ("_clusterext" if args.cluster_extent else "")
    out_dir = ROOT / "data" / "audit" / f"anatomy_{args.band}_grassmann{epi_tag}{suffix}"
    run_band(args.band, k_list, args.n_surrogates, args.epi_x, out_dir,
             scope_tag=scope_tag)


if __name__ == "__main__":
    main()
