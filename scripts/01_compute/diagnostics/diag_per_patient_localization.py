#!/usr/bin/env python3
"""Per-patient localization of the persistence trace — does each patient's trace
concentrate in a subset of its OWN implant's DK regions, and (because implants
differ) does the localized region differ across patients?

This is the strictly-weaker companion to the (retracted) cohort-localization
audit. Cohort localization — "all patients localize to the SAME regions" — was
retracted on 2026-05-30 (no DK region reaches a defensible cohort claim; max
coverage 5/10; `anatomy_localization_wilcoxon/README.md`). The open question the
user posed: maybe the trace localizes WITHIN each patient, just to a different
region per patient, consistent with the heterogeneous sEEG implants.

5-point critical preamble (per `feedback_critical_null_preamble.md`)
-------------------------------------------------------------------
1. Claim: for each patient the persistence trace concentrates in a subset of
   that patient's implanted DK regions (within-subject localization), and the
   localized region differs across patients (idiosyncratic localization).

2. Nulls (TWO; the verdict needs both):
   (a) Within-patient region-label shuffle. For one patient, permute the
       node->region label vector (preserves the per-patient node-count multiset
       of regions = the implant geometry, and the per-contribution value
       distribution; randomizes which node carries which label). Statistic:
       T_p = max over anatomical regions of the region mean signed score.
       Null: T_p under B shuffles. p_p = floor-free permutation rank. This asks
       "does this patient have a region whose trace exceeds a random size-matched
       region of its own implant?". The MAX auto-corrects for testing many
       regions and for small-region variance (the null max is over the same
       region-size multiset, so noisy small regions inflate the null too).
   (b) Idiosyncrasy null. Under (a)'s shuffles, record each patient's argmax
       (null) region -> per-patient null distribution over its own regions.
       Draw one region per patient from these implant-constrained distributions
       many times; null distribution of the cross-patient overlap statistic
       (max recurrence of any single region). Observed overlap of the real top
       regions vs this null tells shared-vs-idiosyncratic.

3. Strongest alternative: "no localization at all." Under it, (i) each patient's
   max region mean is within its own shuffle null (no within-subject
   concentration) AND (ii) top regions still differ across patients TRIVIALLY
   (random argmax over different implants). The non-overlap of top regions is
   the NULL's own prediction, so "they all localize differently" is NOT by
   itself evidence of localization — it is exactly what no-localization predicts.

4. Does the null cover (3) — by mechanism:
   - Null (a) is the hinge. It distinguishes real within-subject concentration
     from scatter: a patient counts as localizing only if its observed hotspot
     beats its own size-matched implant shuffle. Cohort tally of pass count vs
     Binomial(10, 0.05) tests whether MORE patients localize than chance.
   - Null (b) distinguishes shared cohort localization (overlap >> chance —
     would contradict the retraction) from idiosyncratic (overlap ~ chance).
   - Idiosyncrasy is a POSITIVE result ONLY when paired with (a) passing; alone
     it is the no-localization prediction (see point 3).
   - What it CANNOT do: separate "biological per-patient localization" from
     "each patient's most-sampled / highest-degree region trivially carries the
     most trace mass" beyond what the shuffle controls (the shuffle preserves
     region node-counts and randomizes degree x region coupling, so a pure
     size/degree effect is absorbed; a region x value-variance interaction is
     only partly controlled). It cannot establish cross-patient biological
     meaning of the idiosyncratic regions. sEEG sampling caveat: we only see
     implanted regions.

5. Falsification:
   - n_pass ~ chance (Binomial not significant) => per-patient localization NOT
     established; the heterogeneity of top regions is trivial scatter, not a
     result (honest retraction-style verdict).
   - n_pass >> chance AND overlap ~ chance => per-patient IDIOSYNCRATIC
     localization IS a result (the user's hypothesis).
   - n_pass >> chance AND overlap >> chance => SHARED cohort localization
     (would require re-examining the cohort retraction).
   Limitation: even a positive per-patient result is sEEG-sampling-bound and
   per-patient region counts are small; the test is descriptive of where, not
   why.

Probes
------
- cophenet trace score s_ij = dD_task*dD_rest (signed cross-phase TRACE; primary
  for all 6 bands).
- Grassmann phase-averaged participation deviation (ANCHOR-flavoured, not a
  trace — reported as a secondary "where do leading modes live per patient"
  companion, flagged accordingly).

Both probes come from `lrg_eegfc.utils.metrics.node_localization`.

Outputs (under data/audit/per_patient_localization/)
----------------------------------------------------
- {cell}_per_patient.csv  — one row per patient: top region, T_p, p_p, n_regions,
  n_nodes_anatomical, 2nd/3rd regions, demeaned top score.
- summary.csv             — one row per cell: n_pass, Binomial p, obs/null max
  recurrence, idiosyncrasy verdict, distinct-top-region count.
- README.md               — verdict + reading.

Usage
-----
    python diag_per_patient_localization.py
    python diag_per_patient_localization.py --n-perm 5000
    python diag_per_patient_localization.py --band beta --probe cophenet
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.node_localization import (
    NON_ANATOMICAL,
    cophenet_trace_contributions,
    grassmann_participation_contributions,
    grassmann_trace_contributions,
)
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
OUT_ROOT = ROOT / "data/audit/per_patient_localization"

COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
PERM_SEED = 20260601
ALPHA = 0.05                     # per-patient permutation significance level
IDIO_DRAWS = 20000               # Monte-Carlo draws for the idiosyncrasy null

# Non-anatomical labels per grouping. DK region: Wm/Unk/unknown (library const).
# Lobe: the coarse `lobe` column's non-cortical bins. The lobe grouping is the
# higher-powered sensitivity (fewer, larger groups -> the MAX-over-groups null is
# less conservative); it pre-empts the reviewer who would re-run at lobe level.
NON_ANATOMICAL_LOBE = frozenset({"white_matter", "other", "unknown"})

# Grassmann k-set per band (audit-specific config; β contiguous window, γ_l/δ
# cluster-extent S(b)). Mirrors diag_anatomy_localization_wilcoxon.
GRASSMANN_CONTIGUOUS_K = {"beta": list(range(27, 56))}
CLUSTER_EXTENT_CSV = "data/audit/grassmann_cluster_extent/per_k_obs_p.csv"

CELLS = [
    {"band": "delta", "probe": "cophenet", "epi_x": False},
    {"band": "theta", "probe": "cophenet", "epi_x": False},
    {"band": "alpha", "probe": "cophenet", "epi_x": False},
    {"band": "beta", "probe": "cophenet", "epi_x": False},
    {"band": "low_gamma", "probe": "cophenet", "epi_x": False},
    {"band": "high_gamma", "probe": "cophenet", "epi_x": False},
    {"band": "beta", "probe": "grassmann", "epi_x": False},
    {"band": "low_gamma", "probe": "grassmann", "epi_x": False},
    {"band": "delta", "probe": "grassmann", "epi_x": False},
    {"band": "delta", "probe": "grassmann", "epi_x": True},
]


def cell_id(c: dict) -> str:
    return f"{c['band']}_{c['probe']}" + ("_epiX" if c["epi_x"] else "")


def grassmann_k_set(band: str) -> list:
    if band in GRASSMANN_CONTIGUOUS_K:
        return GRASSMANN_CONTIGUOUS_K[band]
    df = pd.read_csv(ROOT / CLUSTER_EXTENT_CSV)
    sub = df[(df["band"] == band) & (df["obs_p_one_sided_less"] < 0.05)]
    return sorted(int(k) for k in sub["k"].values)


def region_to_lobe(pat: str) -> dict:
    """Per-patient DK-region -> lobe lookup (from the `lobe` column)."""
    df = load_channel_regions(pat)
    return dict(zip(df["region"].astype(str), df["lobe"].astype(str)))


def build(pat: str, cell: dict, grouping: str = "region"):
    """Return (values, node_idx, node_group). `grouping` swaps the node label
    vector between DK region (default) and the coarser lobe partition."""
    if cell["probe"] == "cophenet":
        values, node_idx, node_region = cophenet_trace_contributions(
            pat, cell["band"], ROOT)
    elif cell["probe"] == "grassmann_trace":
        values, node_idx, node_region = grassmann_trace_contributions(
            pat, cell["band"], ROOT, grassmann_k_set(cell["band"]),
            epi_x=cell["epi_x"])
    else:
        values, node_idx, node_region = grassmann_participation_contributions(
            pat, cell["band"], ROOT, grassmann_k_set(cell["band"]),
            epi_x=cell["epi_x"])
    if grouping == "lobe":
        r2l = region_to_lobe(pat)
        node_group = np.array([r2l.get(str(r), "unknown") for r in node_region])
        return values, node_idx, node_group
    return values, node_idx, node_region


def region_means_from_codes(values, contrib_codes, n_codes):
    """Per-code (region) mean of `values`, aligned to 0..n_codes-1.

    Returns (means, counts); means is NaN where count==0.
    """
    sums = np.bincount(contrib_codes, weights=values, minlength=n_codes)
    counts = np.bincount(contrib_codes, minlength=n_codes).astype(np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        means = np.where(counts > 0, sums / counts, np.nan)
    return means, counts


def per_patient_test(values, node_idx, node_region, n_perm, rng,
                     non_anatomical=NON_ANATOMICAL):
    """Within-patient max-region-mean localization test for one patient.

    Returns dict with observed top region/score, floor-free permutation p_p,
    per-region observed means (anatomical), and the null argmax-region array
    (region labels, length n_perm) for the idiosyncrasy null. Demeaning the
    values is reporting-only (the max statistic and its label-shuffle null shift
    by the same per-patient constant, so p_p is invariant to it).
    """
    values = np.asarray(values, dtype=np.float64)
    values = values - values.mean()                 # report-comparable; p-invariant
    node_region = np.asarray(node_region).astype(str)
    uniq, codes = np.unique(node_region, return_inverse=True)  # codes length N
    R = uniq.size
    anat = np.array([u not in non_anatomical for u in uniq])
    anat_idx = np.where(anat)[0]
    if anat_idx.size == 0:
        return None

    # observed
    obs_codes = codes[node_idx]
    means, counts = region_means_from_codes(values, obs_codes, R)
    valid = anat & np.isfinite(means)
    valid_idx = np.where(valid)[0]
    obs_means = means[valid_idx]
    obs_top_local = valid_idx[int(np.argmax(obs_means))]
    obs_T = float(means[obs_top_local])
    obs_top_region = str(uniq[obs_top_local])

    # within-patient label-shuffle null of the max statistic + argmax region
    null_T = np.empty(n_perm)
    null_argmax_region = np.empty(n_perm, dtype=object)
    for b in range(n_perm):
        cshuf = rng.permutation(codes)
        m, _ = region_means_from_codes(values, cshuf[node_idx], R)
        ma = m[anat_idx]
        finite = np.isfinite(ma)
        loc = int(np.argmax(np.where(finite, ma, -np.inf)))
        null_T[b] = ma[loc]
        null_argmax_region[b] = str(uniq[anat_idx[loc]])
    p_p = (1 + int(np.sum(null_T >= obs_T))) / (n_perm + 1)

    # observed region ranking (anatomical), for 2nd/3rd reporting
    order = valid_idx[np.argsort(-means[valid_idx])]
    ranked = [(str(uniq[i]), float(means[i]), int(counts[i])) for i in order]

    return {
        "top_region": obs_top_region,
        "T": obs_T,
        "p_p": p_p,
        "n_regions": int(valid.sum()),
        "ranked": ranked,
        "null_argmax_region": null_argmax_region,
    }


def run_cell(cell: dict, n_perm: int, grouping: str = "region",
             verbose: bool = True) -> dict:
    cid = cell_id(cell)
    non_anatomical = (NON_ANATOMICAL_LOBE if grouping == "lobe"
                      else NON_ANATOMICAL)
    if verbose:
        print(f"\n=== {cid} [{grouping}] ===")
    rng = np.random.default_rng(PERM_SEED)

    rows = []
    null_argmax = {}        # pat -> array of null argmax regions
    obs_top = {}            # pat -> observed top region
    p_by_pat = {}
    for pat in COHORT:
        try:
            values, node_idx, node_region = build(pat, cell, grouping)
        except FileNotFoundError as e:
            if verbose:
                print(f"  {pat}: MISSING ({e}) — skipped")
            continue
        res = per_patient_test(values, node_idx, node_region, n_perm, rng,
                               non_anatomical=non_anatomical)
        if res is None:
            continue
        ranked = res["ranked"]
        rows.append({
            "cell": cid, "band": cell["band"], "probe": cell["probe"],
            "epi_x": cell["epi_x"], "patient": pat,
            "top_region": res["top_region"],
            "top_score_demeaned": res["T"],
            "p_perm": res["p_p"],
            "localizes": bool(res["p_p"] < ALPHA),
            "n_regions_sampled": res["n_regions"],
            "region_2": ranked[1][0] if len(ranked) > 1 else "",
            "region_3": ranked[2][0] if len(ranked) > 2 else "",
            "top_n_nodes": ranked[0][2],
        })
        null_argmax[pat] = res["null_argmax_region"]
        obs_top[pat] = res["top_region"]
        p_by_pat[pat] = res["p_p"]

    per_patient_df = pd.DataFrame(rows).sort_values("p_perm").reset_index(drop=True)
    pats = list(obs_top)
    n_used = len(pats)
    n_pass = int(per_patient_df["localizes"].sum())

    # cohort tally vs Binomial(n_used, ALPHA): P(>= n_pass localize by chance)
    binom_p = float(binom.sf(n_pass - 1, n_used, ALPHA)) if n_used else np.nan

    # ---- idiosyncrasy null (max recurrence of any single top region) --------
    def max_recurrence(region_list):
        return max(Counter(region_list).values()) if region_list else 0

    # all patients
    obs_regions_all = [obs_top[p] for p in pats]
    obs_rec_all = max_recurrence(obs_regions_all)
    obs_distinct_all = len(set(obs_regions_all))
    # passing patients only (scientifically meaningful overlap)
    pass_pats = [p for p in pats if p_by_pat[p] < ALPHA]
    obs_regions_pass = [obs_top[p] for p in pass_pats]
    obs_rec_pass = max_recurrence(obs_regions_pass)
    obs_distinct_pass = len(set(obs_regions_pass))

    # null: draw one argmax region per patient from its own null distribution
    rng2 = np.random.default_rng(PERM_SEED + 1)

    def idio_null_p(pat_subset, obs_rec):
        if len(pat_subset) < 2:
            return np.nan, np.nan
        draws = np.empty(IDIO_DRAWS)
        for d in range(IDIO_DRAWS):
            picks = [str(rng2.choice(null_argmax[p])) for p in pat_subset]
            draws[d] = max_recurrence(picks)
        p = (1 + int(np.sum(draws >= obs_rec))) / (IDIO_DRAWS + 1)
        return p, float(draws.mean())

    p_idio_all, null_rec_all = idio_null_p(pats, obs_rec_all)
    p_idio_pass, null_rec_pass = idio_null_p(pass_pats, obs_rec_pass)

    # ---- verdict ------------------------------------------------------------
    if not n_used:
        verdict = "NO DATA"
    elif binom_p >= 0.05:
        verdict = ("NO per-patient localization (n_pass ~ chance; heterogeneous "
                   "top regions are trivial scatter, not a result)")
    else:
        # more patients localize than chance -> look at overlap of passers
        if np.isfinite(p_idio_pass) and p_idio_pass < 0.05:
            verdict = ("SHARED per-patient localization (passers overlap > "
                       "implant-constrained chance — re-examine cohort retraction)")
        else:
            verdict = ("IDIOSYNCRATIC per-patient localization (more patients "
                       "localize than chance; top regions ~ implant-constrained "
                       "chance overlap = different region per patient)")

    if verbose:
        print(f"  patients used: {n_used} | n_pass(p<{ALPHA}): {n_pass} "
              f"(Binomial p={binom_p:.4f}, chance~{ALPHA*n_used:.1f})")
        print(f"  top regions (all): distinct={obs_distinct_all}, "
              f"max_recurrence={obs_rec_all} (null~{null_rec_all}, p={p_idio_all})")
        print(f"  top regions (passers n={len(pass_pats)}): distinct="
              f"{obs_distinct_pass}, max_recurrence={obs_rec_pass} "
              f"(null~{null_rec_pass}, p={p_idio_pass})")
        print(f"  VERDICT: {verdict}")
        for _, r in per_patient_df.iterrows():
            star = "*" if r["localizes"] else " "
            print(f"    {star} {r['patient']}: {r['top_region']:<32s} "
                  f"p={r['p_perm']:.4f}  (n_reg={r['n_regions_sampled']})")

    summary = {
        "cell": cid, "band": cell["band"], "probe": cell["probe"],
        "epi_x": cell["epi_x"], "grouping": grouping, "n_perm": n_perm,
        "n_patients_used": n_used, "n_pass": n_pass,
        "expected_pass_chance": round(ALPHA * n_used, 2),
        "binom_p": binom_p,
        "distinct_top_all": obs_distinct_all,
        "max_recurrence_all": obs_rec_all,
        "null_max_recurrence_all": null_rec_all,
        "p_idiosyncrasy_all": p_idio_all,
        "n_passers": len(pass_pats),
        "distinct_top_passers": obs_distinct_pass,
        "max_recurrence_passers": obs_rec_pass,
        "null_max_recurrence_passers": null_rec_pass,
        "p_idiosyncrasy_passers": p_idio_pass,
        "verdict": verdict,
    }
    return {"summary": summary, "per_patient": per_patient_df}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default=None)
    ap.add_argument("--probe", default=None,
                    choices=["cophenet", "grassmann", "grassmann_trace"])
    ap.add_argument("--epi-x", action="store_true")
    ap.add_argument("--grouping", default="region", choices=["region", "lobe"],
                    help="DK region (default) or coarser lobe (power steelman)")
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=None,
                    help="override PERM_SEED (for Monte-Carlo stability checks)")
    args = ap.parse_args()

    if args.seed is not None:
        global PERM_SEED
        PERM_SEED = args.seed

    if args.band and args.probe:
        cells = [{"band": args.band, "probe": args.probe, "epi_x": args.epi_x}]
    else:
        cells = CELLS

    sfx = "" if args.grouping == "region" else f"_{args.grouping}"
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    summaries = []
    for cell in cells:
        res = run_cell(cell, args.n_perm, grouping=args.grouping)
        summaries.append(res["summary"])
        res["per_patient"].to_csv(
            OUT_ROOT / f"{cell_id(cell)}{sfx}_per_patient.csv", index=False)

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(OUT_ROOT / f"summary{sfx}.csv", index=False)
    print("\n==== SUMMARY ====")
    cols = ["cell", "grouping", "n_patients_used", "n_pass",
            "expected_pass_chance", "binom_p", "distinct_top_all",
            "max_recurrence_all", "p_idiosyncrasy_all", "verdict"]
    print(summary_df[cols].to_string(index=False))
    print(f"\nOutputs → {OUT_ROOT}/")


if __name__ == "__main__":
    main()
