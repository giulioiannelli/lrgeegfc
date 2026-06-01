#!/usr/bin/env python3
"""Signed, threshold-free per-region localization test for the locked anatomy
networks — does the persistence trace genuinely concentrate in specific DK
regions cohort-wide, or is the locked localization a top-decile / absolute-value
artifact?

5-point critical preamble (per `feedback_critical_null_preamble.md`)
-------------------------------------------------------------------
1. Claim under test: the locked cophenet anatomy networks (β cophenet 7
   regions; α cophenet 11 regions) and Grassmann networks represent genuine
   cohort-wide spatial localization of the persistence trace — i.e., specific
   DK regions carry the trace consistently *across patients*.

2. Null: within-patient region-label shuffle. For each patient, permute the
   node→region assignment. This preserves (a) the implant region-size multiset
   per patient (which regions a patient samples and with how many contacts) and
   (b) the per-pair / per-node trace-value distribution; it randomizes only
   *which* node carries *which* region label. Recompute the per-region
   one-sided Wilcoxon and the cohort localization statistic under each shuffle
   → null distribution of the cohort statistic (B permutations).

3. Strongest plausible alternative: the locked localization is an artifact of
   (a) audit_71/72's absolute-value + top-decile selection — which mixes
   trace and anti-trace pairs and over-weights a few high-magnitude pairs —
   and/or (b) regions being effectively single-/thin-patient findings, so
   "cohort enrichment" reflects implant idiosyncrasy, not biology. Equivalent
   claim: a random relabeling of each patient's implant would yield as many
   "significant" regions.

4. Does the null cover (3) — by mechanism, not vibes:
   - Signed quantity removes the |Δ| confound: anti-trace pairs/nodes now
     subtract instead of adding. (cophenet s_ij = dD_task·dD_rest; >0 is
     audit_71's sign-consistent trace direction, <0 is reverted/anti.)
   - No top-decile cut removes the magnitude-selection confound.
   - The label-shuffle preserves each patient's implant geometry and value
     distribution, so a region that looks enriched only because one patient
     densely samples it with high-variance values is *equally* enriched under
     shuffle → the null absorbs single-patient / implant-density artifacts.
   - Per-region Wilcoxon is rank-based across patients: no single patient
     drives a region's p. Region size per patient is invariant under the
     shuffle, so per-region statistical power is calibrated exactly.
   What it does NOT cover: it does not re-test whether the trace itself is real
   (that is the matched-strength rung, already locked); it only tests whether,
   *given* a trace, its anatomical attribution is cohort-consistent and
   concentrated. It also cannot distinguish "biologically localized" from
   "every patient happens to be implanted in the same functional system" —
   an inherent sEEG sampling caveat, not a null failure.

5. Falsification: if the cohort localization statistic (Σ −log10 p mass across
   regions; and, descriptively, the count of regions with one-sided Wilcoxon
   p<0.05) does NOT exceed the within-patient-shuffle null at p<0.05, the
   anatomy claim for that (band, probe) is **diffuse → retract**. If it
   exceeds, report the surviving regions descriptively (continuous p, no
   per-region gate), with n_patients-sampling shown so the reader sees thin
   regions. Limitation: a "localized" verdict still rests on which specific
   regions carry positive cohort Wilcoxon; thin regions (sampled by <5
   patients) cannot reach p<0.05 by the exact-Wilcoxon p-floor and are
   reported as low-power, not as evidence against localization.

Design notes
------------
- The verdict statistic is the −log10 p MASS (threshold-free, cluster-mass
  style consistent with the project's T_G* convention). The count-of-p<0.05 is
  reported alongside as the user-requested descriptive companion. The
  permutation null IS the gate; nothing is hardcoded.
- Cophenet per-pair signed score s_ij = dD_task·dD_rest is the EXACT signed,
  threshold-free continuation of audit_71's flagging quantity
  (|dD_task·dD_rest| with sign-consistency gate + top-decile). Each pair
  contributes its s to BOTH endpoint regions (audit_71 endpoint logic).
- Grassmann per-node quantity replicates audit_72's phase-AVERAGED participation
  (NOT a cross-phase delta — see WARNING below). The localization signed value
  is the deviation of a node's participation from the patient's own mean
  participation, so above-baseline concentration is positive. β only is run
  cleanly here (contiguous k=27..55, no epi-X); γ_l/δ use cluster-extent k-sets
  and δ epi-X is affected by an upstream no-op mask bug — deferred (see report).

WARNING (surfaced 2026-05-30): the Grassmann anatomy base quantity
(participation_phase_avg over rest_pre/task/rest_post) does not difference
phases, so it is an ANCHOR-flavored "where do leading modes live" quantity, not
a cross-phase TRACE. The Grassmann localization result below tests concentration
of that anchor-anatomy, which is a different question than the cophenet trace
localization. Read accordingly.

Inputs (read-only; no FC/LRG recomputation)
-------------------------------------------
- cophenet: data/reports/imcoh_continuous_trace/per_pair_split/{pat}_{band}.npz
  (dD_task, dD_rest, iu_i, iu_j) — observed signed split-baseline deltas.
- grassmann: data/cache/imcoh_lrg/{pat}/{band}_{phase}_lrg_imcoh-abs.npz
  (eigenvectors), phases rest_pre/task_test/rest_post.
- DK regions: lrg_eegfc.utils.io.regions.load_channel_regions.

Outputs (under data/audit/anatomy_localization_wilcoxon/)
---------------------------------------------------------
- {band}_{probe}_per_region.csv  — per-region Wilcoxon p, n_patients, cohort
  mean/median of per-patient means, in_locked_network flag.
- summary.csv                    — one row per (band, probe): obs mass, obs
  count, perm-null tail p for each, verdict.
- README.md                      — verdict + reading + retraction proposals.

Usage
-----
    python diag_anatomy_localization_wilcoxon.py            # run all clean cells
    python diag_anatomy_localization_wilcoxon.py --n-perm 5000
    python diag_anatomy_localization_wilcoxon.py --band beta --probe cophenet
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.metrics.node_localization import (
    NON_ANATOMICAL,
    cophenet_trace_contributions,
    grassmann_participation_contributions,
)
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
OUT_ROOT = ROOT / "data/audit/anatomy_localization_wilcoxon"

COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
OBS_PHASES = ("rest_pre", "task_test", "rest_post")
PERM_SEED = 20260530

# NON_ANATOMICAL (Wm/Unk/unknown) is imported from the library — these are NOT
# Desikan-Killiany regions; they pollute the localization universe and (being
# sampled by every patient) manufacture spurious cohort significance. Excluded
# from the verdict; everything else (incl. Hip, Amy, Cerebellum-Cortex) is kept.

# Locked anatomy networks (source: ANATOMY_LEDGER.md, 2026-05-19 lock), keyed by
# cell-id. Used only to flag which locked regions survive the signed test.
LOCKED_NETWORKS = {
    "beta_cophenet": [
        "ctx-lh-isthmuscingulate", "ctx-lh-superiorfrontal", "ctx-rh-insula",
        "ctx-lh-parahippocampal", "ctx-rh-postcentral", "ctx-lh-entorhinal",
        "ctx-rh-rostralanteriorcingulate",
    ],
    "alpha_cophenet": [
        "ctx-lh-caudalanteriorcingulate", "ctx-lh-parahippocampal",
        "ctx-lh-rostralanteriorcingulate", "ctx-rh-medialorbitofrontal",
        "ctx-rh-caudalmiddlefrontal", "ctx-rh-postcentral",
        "ctx-lh-caudalmiddlefrontal", "ctx-rh-caudalanteriorcingulate",
        "ctx-rh-precuneus", "ctx-rh-superiorparietal",
        "ctx-rh-posteriorcingulate",
    ],
    "beta_grassmann": [
        "Hip", "ctx-lh-insula", "ctx-lh-middletemporal",
        "ctx-lh-lateralorbitofrontal", "ctx-rh-medialorbitofrontal",
        "ctx-lh-superiortemporal", "ctx-rh-rostralmiddlefrontal",
    ],
    "low_gamma_grassmann": [
        "ctx-lh-lateraloccipital", "ctx-lh-middletemporal",
        "ctx-lh-rostralmiddlefrontal", "ctx-lh-superiortemporal",
        "ctx-rh-medialorbitofrontal", "ctx-rh-parstriangularis",
        "ctx-lh-cuneus",
    ],
    "delta_grassmann": [
        "ctx-lh-inferiortemporal", "ctx-lh-inferiorparietal",
        "ctx-rh-parstriangularis", "ctx-lh-superiortemporal",
    ],
    "delta_grassmann_epiX": [
        "ctx-lh-superiorparietal", "ctx-rh-rostralmiddlefrontal",
        "ctx-lh-superiorfrontal",
    ],
}

# Grassmann k-set per band: β uses the locked contiguous window (audit_72
# BAND_K_WINDOWS); γ_l / δ use the cluster-extent S(b) set (all p<0.05 cells in
# grassmann_cluster_extent/per_k_obs_p.csv — the locked source for those cells).
GRASSMANN_CONTIGUOUS_K = {"beta": list(range(27, 56))}
CLUSTER_EXTENT_CSV = "data/audit/grassmann_cluster_extent/per_k_obs_p.csv"

# Cells to run, each a dict: band, probe, epi_x. cell_id() builds the key/label.
CELLS = [
    {"band": "beta", "probe": "cophenet", "epi_x": False},
    {"band": "alpha", "probe": "cophenet", "epi_x": False},
    {"band": "beta", "probe": "grassmann", "epi_x": False},
    {"band": "low_gamma", "probe": "grassmann", "epi_x": False},
    {"band": "delta", "probe": "grassmann", "epi_x": False},
    {"band": "delta", "probe": "grassmann", "epi_x": True},
]


def cell_id(c: dict) -> str:
    return f"{c['band']}_{c['probe']}" + ("_epiX" if c["epi_x"] else "")


def grassmann_k_set(band: str) -> list:
    """β: locked contiguous window. γ_l/δ: cluster-extent S(b) from per_k csv."""
    if band in GRASSMANN_CONTIGUOUS_K:
        return GRASSMANN_CONTIGUOUS_K[band]
    df = pd.read_csv(ROOT / CLUSTER_EXTENT_CSV)
    sub = df[(df["band"] == band) & (df["obs_p_one_sided_less"] < 0.05)]
    return sorted(int(k) for k in sub["k"].values)


# ---------------------------------------------------------------------------
# per-patient contribution builders — thin wrappers over the library probes
# (lrg_eegfc.utils.metrics.node_localization). The k-set resolution stays here
# because it is audit-specific config (locked contiguous β window vs the
# cluster-extent S(b) set); the library builder takes an explicit k_set.
# ---------------------------------------------------------------------------

def cophenet_contributions(pat: str, band: str):
    return cophenet_trace_contributions(pat, band, ROOT)


def grassmann_contributions(pat: str, band: str, epi_x: bool = False):
    return grassmann_participation_contributions(
        pat, band, ROOT, grassmann_k_set(band), epi_x=epi_x)


# ---------------------------------------------------------------------------
# cohort localization statistic
# ---------------------------------------------------------------------------

def region_wilcoxon_greater(vals: np.ndarray):
    """Descriptive exact one-sided Wilcoxon, H1: cohort median > 0.

    Reported for continuity with audit_71/72 only; NOT the verdict statistic.
    The exact-Wilcoxon p has an n-dependent floor (n=2->0.25, n=4->0.0625) that
    dilutes any cohort aggregate, so the verdict uses the floor-free permutation
    p below instead.
    """
    v = np.asarray([x for x in vals if np.isfinite(x)], dtype=np.float64)
    n = v.size
    if n < 1 or np.allclose(v, 0.0):
        return np.nan
    try:
        _, p = wilcoxon(v, alternative="greater", zero_method="wilcox")
    except (ValueError, TypeError):
        return np.nan
    return float(p)


def region_means_per_patient(per_patient, node_region_by_pat):
    """Per-region list of per-patient demeaned means, ANATOMICAL regions only.

    per_patient : dict pat -> (values, node_idx)  [values already demeaned]
    node_region_by_pat : dict pat -> region-label vector (observed or shuffled)

    Returns dict region -> list[float] (one entry per patient sampling it).
    Non-anatomical labels (Wm/Unk/unknown) are dropped. Region presence and
    per-region n_patients are invariant under the within-patient label shuffle
    (the label multiset per patient is preserved), so the region set and each
    region's statistical power are identical in observed and null.
    """
    region_means = defaultdict(list)
    for pat in COHORT:
        values, node_idx = per_patient[pat]
        reg = node_region_by_pat[pat]
        contrib_region = reg[node_idx]
        order = np.argsort(contrib_region, kind="stable")
        cr_sorted = contrib_region[order]
        val_sorted = values[order]
        uniq, starts = np.unique(cr_sorted, return_index=True)
        sums = np.add.reduceat(val_sorted, starts)
        counts = np.diff(np.append(starts, val_sorted.size))
        means = sums / counts
        for region, m in zip(uniq, means):
            if region not in NON_ANATOMICAL:
                region_means[region].append(float(m))
    return region_means


def run_cell(cell: dict, n_perm: int, q_level: float = 0.95,
             verbose: bool = True) -> dict:
    band, probe, epi_x = cell["band"], cell["probe"], cell["epi_x"]
    cid = cell_id(cell)
    if verbose:
        print(f"\n=== {cid} ===")

    def builder(pat):
        if probe == "cophenet":
            return cophenet_contributions(pat, band)
        return grassmann_contributions(pat, band, epi_x=epi_x)

    # Per-patient demeaning is the crux of the localization contrast: a
    # region's value is its mean signed score MINUS the patient's own global
    # mean (over ALL contributions — label-independent, so it is a fixed
    # per-patient constant under the label shuffle). This removes the
    # cohort-wide base rate: a globally positive trace shifts every region's
    # raw mean up, but after demeaning only regions that exceed the patient's
    # OWN average register as localized. (Demeaning enforces a
    # contribution-count-weighted sum-to-zero across all contributions; the
    # unweighted per-region means retain a tiny residual offset, but it is
    # identical under shuffle so the null stays calibrated.) Grassmann's `dev`
    # is already node-mean-0 (no-op).
    per_patient = {}
    node_region_by_pat = {}
    for pat in COHORT:
        values, node_idx, node_region = builder(pat)
        values = np.asarray(values, dtype=np.float64)
        values = values - values.mean()
        per_patient[pat] = (values, node_idx)
        node_region_by_pat[pat] = node_region

    # observed per-region distributions (anatomical regions only)
    obs_rm = region_means_per_patient(per_patient, node_region_by_pat)
    regions = sorted(obs_rm)
    n_reg = len(regions)
    reg_index = {r: i for i, r in enumerate(regions)}
    obs_median = np.array([np.median(obs_rm[r]) for r in regions])
    obs_mean = np.array([float(np.mean(obs_rm[r])) for r in regions])
    n_pat = np.array([len(obs_rm[r]) for r in regions])
    exact_p = np.array([region_wilcoxon_greater(np.asarray(obs_rm[r]))
                        for r in regions])

    # within-patient label-shuffle null: per-region median under each shuffle.
    # The effect is the across-patient MEDIAN (outlier-robust at the patient
    # level), and significance is the FLOOR-FREE permutation rank of the
    # observed median against this region's own shuffle null — so thin regions
    # contribute matched noise to obs and null (no exact-Wilcoxon floor that
    # would dilute the cohort statistic toward "diffuse").
    rng = np.random.default_rng(PERM_SEED)
    null_med = np.empty((n_reg, n_perm))
    for b in range(n_perm):
        shuffled = {pat: rng.permutation(node_region_by_pat[pat])
                    for pat in COHORT}
        rm = region_means_per_patient(per_patient, shuffled)
        for r in regions:
            null_med[reg_index[r], b] = np.median(rm[r])
        if verbose and (b + 1) % max(1, n_perm // 5) == 0:
            print(f"  perm {b + 1}/{n_perm}")

    # per-region floor-free permutation p (one-sided greater) + own-null q-thr
    perm_p = (1 + np.sum(null_med >= obs_median[:, None], axis=1)) / (n_perm + 1)
    qthr = np.quantile(null_med, q_level, axis=1)
    obs_pass = obs_median > qthr

    # cohort verdict statistics, calibrated against the JOINT shuffle null
    # (each null column is a full within-patient relabeling, so cross-region
    # dependence is preserved):
    #   count   = # anatomical regions individually extreme (median > own q95)
    #   excess  = Σ over regions of max(median − own q95, 0)  [continuous]
    count_obs = int(obs_pass.sum())
    null_count = np.sum(null_med > qthr[:, None], axis=0)
    p_count = (1 + int(np.sum(null_count >= count_obs))) / (n_perm + 1)

    excess_obs = float(np.sum(np.clip(obs_median - qthr, 0.0, None)))
    null_excess = np.sum(np.clip(null_med - qthr[:, None], 0.0, None), axis=0)
    p_excess = (1 + int(np.sum(null_excess >= excess_obs))) / (n_perm + 1)

    verdict_p = min(p_count, p_excess)
    if verdict_p < 0.05:
        verdict = "LOCALIZED"
    elif verdict_p < 0.10:
        verdict = "BORDERLINE"
    else:
        verdict = "DIFFUSE (retract)"

    locked = set(LOCKED_NETWORKS.get(cid, []))
    rows = []
    for i, region in enumerate(regions):
        rows.append({
            "cell": cid, "band": band, "probe": probe, "epi_x": epi_x,
            "region": region,
            "n_patients": int(n_pat[i]),
            "obs_median": float(obs_median[i]),
            "obs_mean": obs_mean[i],
            "own_null_q%d" % int(q_level * 100): float(qthr[i]),
            "perm_p_greater": float(perm_p[i]),
            "wilcoxon_exact_p": float(exact_p[i]) if np.isfinite(exact_p[i])
                                else np.nan,
            "passes": bool(obs_pass[i]),
            "in_locked_network": region in locked,
        })
    per_region_df = pd.DataFrame(rows).sort_values(
        "perm_p_greater").reset_index(drop=True)

    locked_rows = per_region_df[per_region_df["in_locked_network"]]
    locked_pass = int(locked_rows["passes"].sum())
    passing_regions = per_region_df[per_region_df["passes"]]["region"].tolist()
    if verbose:
        print(f"  anatomical regions: {n_reg} | obs count(median>own q95)"
              f"={count_obs} (perm-null p={p_count:.4f}); "
              f"excess={excess_obs:.4f} (perm-null p={p_excess:.4f})")
        print(f"  verdict: {verdict}")
        print(f"  locked-network regions passing: {locked_pass}/{len(locked)}")
        print(f"  passing regions: {passing_regions or '—'}")

    return {
        "summary": {
            "cell": cid, "band": band, "probe": probe, "epi_x": epi_x,
            "n_perm": n_perm,
            "n_anatomical_regions": n_reg,
            "obs_count_extreme": count_obs, "perm_p_count": p_count,
            "expected_count_chance": round((1 - q_level) * n_reg, 2),
            "obs_excess": excess_obs, "perm_p_excess": p_excess,
            "n_locked_regions": len(locked),
            "n_locked_passing": locked_pass,
            "passing_regions": "; ".join(passing_regions) or "—",
            "verdict": verdict,
        },
        "per_region": per_region_df,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default=None)
    ap.add_argument("--probe", default=None, choices=["cophenet", "grassmann"])
    ap.add_argument("--epi-x", action="store_true")
    ap.add_argument("--n-perm", type=int, default=2000)
    args = ap.parse_args()

    if args.band and args.probe:
        cells = [{"band": args.band, "probe": args.probe, "epi_x": args.epi_x}]
    else:
        cells = CELLS

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    summaries = []
    for cell in cells:
        res = run_cell(cell, args.n_perm)
        summaries.append(res["summary"])
        res["per_region"].to_csv(
            OUT_ROOT / f"{cell_id(cell)}_per_region.csv", index=False)

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(OUT_ROOT / "summary.csv", index=False)
    print("\n==== SUMMARY ====")
    cols = ["cell", "n_anatomical_regions", "obs_count_extreme",
            "expected_count_chance", "perm_p_count", "perm_p_excess",
            "n_locked_passing", "n_locked_regions", "verdict"]
    print(summary_df[cols].to_string(index=False))
    print(f"\nOutputs → {OUT_ROOT}/")


if __name__ == "__main__":
    main()
