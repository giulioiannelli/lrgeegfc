#!/usr/bin/env python3
"""audit_81 — Sampling-conditioned localization atlas (descriptive, multi-granularity).

Scope report: `.agents/guides/task-persistence-investigation/2026-06-05_localization-atlas.md`.

Re-examination of the over-aggressive 2026-05-30 anatomy retraction. Drops the
cohort permutation GATE (a structural false-negative at n=10 / <=5-per-region) for
a DESCRIPTIVE, sampling-conditioned reading: per anatomical unit `g`, the positive
rate `r_g = J+_g / K_g` = (patients in whom the trace concentrates in `g` above
their own baseline) / (patients IMPLANTED in `g`). The denominator is patients
implanted there, not the full cohort: Hip-beta reads 4/5, not 4/10. Four
granularities (DK region -> system -> lobe -> hemisphere); both matched-strength
TRACE probes (cophenet edge + cross-phase Grassmann node; NOT the phase-avg
anchor); epi-include/exclude split; depth-shaft adjacency control reported
alongside with a deep_target flag.

5-point critical preamble (per feedback_critical_null_preamble.md)
------------------------------------------------------------------
1. Claim: given a verified cohort trace (beta/alpha, matched-strength-locked), the
   trace concentrates in specific anatomy among patients sampled there, possibly
   at a coarser SYSTEM granularity even when no single DK region recurs.
2. Null (DEMOTED to a column, not the gate): within-patient node->unit label
   shuffle (preserves implant unit-size multiset + value distribution).
3. Strongest alternative: (a) the positive rate is a coin-flip artifact (K=5, 4/5
   has binom p~0.19 alone); (b) apparent concentration is electrode-shaft spatial
   autocorrelation, not anatomy.
4. Coverage: (a) the rate is NEVER read alone -- paired with effect-size
   distribution (median_m, median_z), the demoted perm_p, and a stated prior;
   (b) the shaft control (eta2 region-vs-shaft + NMI) is shown side by side, with
   a deep_target flag: for deep depth targets (Hip/Amy/Thal) the electrode IS the
   structure, so shaft-level concentration there is the localization, not a
   confound; the shaft control discounts only cortical-surface adjacency.
5. Falsification/limits: a unit's localization is descriptive + CONDITIONAL (holds
   only among implanted patients); it does NOT re-test trace existence (demeaning
   removes the base rate). High rate + near-zero effect, or a rate that collapses
   under epi-exclusion with epileptic support, is not a clean localization.

NB. At REGION granularity in INCLUDE mode this reproduces the locked
`anatomy_localization_wilcoxon` numbers exactly (same demeaning, same shuffle null,
node-level aggregation is algebraically identical to the contribution-level one) --
a built-in validation of the atlas against the gated audit.

Inputs (read-only; no FC/LRG recompute)
---------------------------------------
- cophenet: data/reports/imcoh_continuous_trace/per_pair_split/{pat}_{band}.npz
- grassmann: data/cache/imcoh_lrg/{pat}/{band}_{phase}_lrg_imcoh-abs.npz
- regions/systems/shaft: lrg_eegfc.utils.io.regions + node_localization lib.

Outputs (data/audit/localization_atlas/)
-----------------------------------------
- atlas_{probe}_{epi}.csv         per (band, granularity, unit) rate + effects
- per_patient_{probe}_{epi}.csv   per (band, granularity, patient) top-unit + system
- shaft_{probe}_{epi}.csv         eta2 region vs shaft + NMI (deep-target annotated)
- summary.csv                     headline rows (beta Hip/MTL highlighted)
- README.md (written separately after inspection)

Usage
-----
    python audit_81_localization_atlas.py                      # full sweep
    python audit_81_localization_atlas.py --probe cophenet --band beta
    python audit_81_localization_atlas.py --n-perm 2000
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.io.regions import (
    ANATOMICAL_SYSTEMS, DEEP_TARGET_SYSTEMS, load_channel_regions,
)
from lrg_eegfc.utils.metrics.node_localization import (
    NON_ANATOMICAL, cophenet_trace_contributions, eta2,
    epi_keep_mask, grassmann_trace_contributions, nmi, perm_p_eta2, shaft_of,
)
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
OUT_ROOT = ROOT / "data/audit/localization_atlas"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
GRANULARITIES = ["region", "system", "lobe", "hemi"]
PERM_SEED = 20260605

GRASSMANN_CONTIGUOUS_K = {"beta": list(range(27, 56))}
CLUSTER_EXTENT_CSV = "data/audit/grassmann_cluster_extent/per_k_obs_p.csv"

# Non-anatomical / non-interpretable units dropped at SCORING per granularity
# (they still enter the shuffle multiset, matching the locked audit; only the
# scored columns are filtered). region: white-matter/unknown; system: the
# catch-all; lobe: white-matter/unknown; hemi: undetermined hemisphere.
DROPSET = {
    "region": set(NON_ANATOMICAL),
    "system": {"non_anatomical"},
    "lobe": {"white_matter", "unknown"},
    "hemi": {"?"},
}
_DEEP_REGION_TOKENS = frozenset(
    {"Hip", "Amy"} | set(ANATOMICAL_SYSTEMS["subcortical_other"]))


def grassmann_k_set(band: str) -> list:
    if band in GRASSMANN_CONTIGUOUS_K:
        return GRASSMANN_CONTIGUOUS_K[band]
    df = pd.read_csv(ROOT / CLUSTER_EXTENT_CSV)
    sub = df[(df["band"] == band) & (df["obs_p_one_sided_less"] < 0.05)]
    return sorted(int(k) for k in sub["k"].values)


def contributions(probe: str, pat: str, band: str):
    """(values, node_idx, _) for a TRACE probe, NO internal epi masking."""
    if probe == "cophenet":
        return cophenet_trace_contributions(pat, band, ROOT)
    return grassmann_trace_contributions(pat, band, ROOT, grassmann_k_set(band))


def is_deep(granularity: str, unit: str) -> bool:
    if granularity == "region":
        return unit in _DEEP_REGION_TOKENS
    if granularity == "system":
        return unit in DEEP_TARGET_SYSTEMS
    return False


def unit_vectors(regions_df) -> dict:
    """Per-node unit label vector for each granularity (length N)."""
    return {
        "region": regions_df["region"].to_numpy().astype(str),
        "system": regions_df["system"].to_numpy().astype(str),
        "lobe": regions_df["lobe"].to_numpy().astype(str),
        "hemi": regions_df["hemisphere"].to_numpy().astype(str),
    }


def node_aggregates(values, node_idx, n_nodes):
    """Per-node sum/count of (demeaned) contribution values.

    Aggregating per endpoint node first makes the within-patient label shuffle a
    permutation over N~117 nodes instead of over ~10^4 contributions -- exactly
    equivalent (a unit's mean = sum of its nodes' sums / sum of its nodes' counts)
    but ~100x faster.
    """
    node_sum = np.bincount(node_idx, weights=values, minlength=n_nodes)
    node_cnt = np.bincount(node_idx, minlength=n_nodes).astype(np.float64)
    return node_sum, node_cnt


def per_patient_band(probe, pat, band, epi_excl):
    """Build the node-level aggregation for one (probe, patient, band, epi).

    Returns dict with: node_sum, node_cnt (over SAMPLED nodes), and the per-node
    unit-code vectors per granularity (over the same sampled nodes), plus the
    global-unit name lists. Demeaning baseline = the analysed (epi-kept) set.
    """
    values, node_idx, _ = contributions(probe, pat, band)
    values = np.asarray(values, dtype=np.float64)
    node_idx = np.asarray(node_idx).astype(int)
    regions_df = load_channel_regions(pat)
    n_nodes = len(regions_df)

    if epi_excl:
        keep_node = epi_keep_mask(regions_df, pat)           # label-based, fixed
        cmask = keep_node[node_idx]
        values = values[cmask]
        node_idx = node_idx[cmask]

    values = values - values.mean()                          # demean over analysed set
    node_sum, node_cnt = node_aggregates(values, node_idx, n_nodes)

    sampled = node_cnt > 0                                    # nodes with >=1 contribution
    uv = unit_vectors(regions_df)
    return {
        "node_sum": node_sum[sampled],
        "node_cnt": node_cnt[sampled],
        "units": {g: uv[g][sampled] for g in GRANULARITIES},
        "label_raw": regions_df["label_raw"].to_numpy()[sampled],
    }


def unit_means(node_sum, node_cnt, code, n_units):
    """Per-unit mean = (sum of member-node sums) / (sum of member-node counts)."""
    s = np.bincount(code, weights=node_sum, minlength=n_units)
    c = np.bincount(code, weights=node_cnt, minlength=n_units)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = np.where(c > 0, s / c, np.nan)
    return m


def run_cell(probe, band, epi_excl, granularity, n_perm, rng):
    """One (probe, band, epi, granularity) -> list of per-unit atlas rows."""
    pp = {p: per_patient_band(probe, p, band, epi_excl) for p in COHORT}

    # global unit index across cohort (scored units only)
    drop = DROPSET[granularity]
    global_units = sorted(
        {u for p in COHORT for u in set(pp[p]["units"][granularity])
         if u not in drop})
    gidx = {u: i for i, u in enumerate(global_units)}
    G = len(global_units)
    if G == 0:
        return []

    # observed per-patient per-unit means + per-node global codes (incl. drop
    # units, which map to code -1 and are ignored by the masked bincount below)
    obs = np.full((len(COHORT), G), np.nan)
    pat_codes = {}
    for pi, p in enumerate(COHORT):
        units_p = pp[p]["units"][granularity]
        code = np.array([gidx.get(u, -1) for u in units_p])  # -1 = drop unit
        pat_codes[p] = code
        keep = code >= 0
        if keep.sum() == 0:
            continue
        obs[pi] = unit_means(pp[p]["node_sum"][keep], pp[p]["node_cnt"][keep],
                             code[keep], G)

    # within-patient label-shuffle null (permute each patient's node->unit codes)
    null = np.full((n_perm, len(COHORT), G), np.nan)
    for b in range(n_perm):
        for pi, p in enumerate(COHORT):
            code = pat_codes[p]
            cs = rng.permutation(code)
            keep = cs >= 0
            if keep.sum() == 0:
                continue
            null[b, pi] = unit_means(pp[p]["node_sum"][keep],
                                     pp[p]["node_cnt"][keep], cs[keep], G)

    # cohort median null per unit (for the demoted perm_p). Non-sampler patients
    # are all-NaN slices (sampling is shuffle-invariant) -> nan-aggregates warn
    # harmlessly; suppress and guard downstream (sd>0).
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        M_null = np.nanmedian(null, axis=1)                  # (n_perm, G)
        null_mean = np.nanmean(null, axis=0)                 # (COHORT, G) for z
        null_std = np.nanstd(null, axis=0)

    rows = []
    for u, g in gidx.items():
        col = obs[:, g]
        samplers = np.where(np.isfinite(col))[0]
        K = samplers.size
        if K == 0:
            continue
        vals = col[samplers]
        M = float(np.median(vals))
        Jpos = int(np.sum(vals > 0))
        # per-patient z = (obs - null_mean)/null_std over samplers
        zs = []
        for pi in samplers:
            sd = null_std[pi, g]
            if sd > 0:
                zs.append((obs[pi, g] - null_mean[pi, g]) / sd)
        Z = float(np.median(zs)) if zs else np.nan
        Mn = M_null[:, g]
        Mn = Mn[np.isfinite(Mn)]
        perm_p = (1 + int(np.sum(Mn >= M))) / (Mn.size + 1) if Mn.size else np.nan
        rows.append({
            "probe": probe, "band": band,
            "epi": "exclude" if epi_excl else "include",
            "granularity": granularity, "unit": u,
            "K_implanted": K, "n_pos": Jpos, "n_neg": K - Jpos,
            "rate": Jpos / K, "median_m": M, "median_z": Z,
            "perm_p": perm_p, "deep_target": is_deep(granularity, u),
            "samplers": ";".join(COHORT[i].split("_")[1] for i in samplers),
        })
    return rows


def per_patient_rows(probe, band, epi_excl, granularity):
    """Per-patient top-unit (argmax mean) + its system, for region/system gran."""
    pp = {p: per_patient_band(probe, p, band, epi_excl) for p in COHORT}
    drop = DROPSET[granularity]
    rows = []
    for p in COHORT:
        units_p = pp[p]["units"][granularity]
        keep = np.array([u not in drop for u in units_p])
        if keep.sum() == 0:
            continue
        uu = units_p[keep]
        uniq, code = np.unique(uu, return_inverse=True)
        m = unit_means(pp[p]["node_sum"][keep], pp[p]["node_cnt"][keep],
                       code, uniq.size)
        top = int(np.nanargmax(m))
        top_unit = uniq[top]
        # which SYSTEM does the top unit belong to (for region granularity)
        if granularity == "region":
            sysv = pp[p]["units"]["system"][keep]
            sys_top = sysv[uu == top_unit]
            top_system = sys_top[0] if len(sys_top) else "?"
        else:
            top_system = top_unit
        rows.append({
            "probe": probe, "band": band,
            "epi": "exclude" if epi_excl else "include",
            "granularity": granularity, "patient": p.split("_")[1],
            "top_unit": top_unit, "top_m": float(m[top]),
            "top_system": top_system,
            "n_units_sampled": int(uniq.size),
        })
    return rows


def shaft_rows(probe, band, epi_excl, n_perm, rng):
    """eta2 region-vs-shaft + NMI per (probe, band, epi) -- the shaft control."""
    rows = []
    for p in COHORT:
        values, node_idx, node_region = contributions(probe, p, band)
        values = np.asarray(values, dtype=np.float64)
        node_idx = np.asarray(node_idx).astype(int)
        regions_df = load_channel_regions(p)
        if epi_excl:
            keep_node = epi_keep_mask(regions_df, p)
            cm = keep_node[node_idx]
            values, node_idx = values[cm], node_idx[cm]
        values = values - values.mean()
        node_region = np.asarray(node_region).astype(str)
        node_shaft = np.array([shaft_of(l) for l in regions_df["label_raw"]])
        anat = np.array([r not in NON_ANATOMICAL for r in node_region])
        kc = anat[node_idx]
        v, ni = values[kc], node_idx[kc]
        obs_r, p_r, n_r = perm_p_eta2(v, node_region, ni, n_perm, rng)
        obs_s, p_s, n_s = perm_p_eta2(v, node_shaft, ni, n_perm, rng)
        coll = nmi(node_region[anat], node_shaft[anat])
        rows.append({
            "probe": probe, "band": band,
            "epi": "exclude" if epi_excl else "include", "patient": p.split("_")[1],
            "eta2_region": obs_r, "p_region": p_r, "n_regions": n_r,
            "eta2_shaft": obs_s, "p_shaft": p_s, "n_shafts": n_s,
            "region_shaft_nmi": coll,
        })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", choices=["cophenet", "grassmann"], default=None)
    ap.add_argument("--band", default=None)
    ap.add_argument("--epi-mode", choices=["include", "exclude"], default=None)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--no-shaft", action="store_true")
    args = ap.parse_args()

    probes = [args.probe] if args.probe else ["cophenet", "grassmann"]
    bands = [args.band] if args.band else BANDS
    epimodes = ([args.epi_mode == "exclude"] if args.epi_mode
                else [False, True])

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    summary = []
    for probe in probes:
        for epi_excl in epimodes:
            etag = "exclude" if epi_excl else "include"
            rng = np.random.default_rng(PERM_SEED)
            atlas, prof, shaft = [], [], []
            for band in bands:
                for gran in GRANULARITIES:
                    atlas += run_cell(probe, band, epi_excl, gran,
                                      args.n_perm, rng)
                for gran in ("region", "system"):
                    prof += per_patient_rows(probe, band, epi_excl, gran)
                if not args.no_shaft:
                    shaft += shaft_rows(probe, band, epi_excl, args.n_perm, rng)
            adf = pd.DataFrame(atlas)
            adf.to_csv(OUT_ROOT / f"atlas_{probe}_{etag}.csv", index=False)
            pd.DataFrame(prof).to_csv(
                OUT_ROOT / f"per_patient_{probe}_{etag}.csv", index=False)
            if shaft:
                pd.DataFrame(shaft).to_csv(
                    OUT_ROOT / f"shaft_{probe}_{etag}.csv", index=False)

            # headline console: MTL + Hip across bands
            print(f"\n=== {probe} / epi-{etag} : MTL(system) + Hip(region) ===")
            for band in bands:
                for gran, unit in (("system", "MTL"), ("region", "Hip")):
                    r = adf[(adf.band == band) & (adf.granularity == gran)
                            & (adf.unit == unit)]
                    if len(r):
                        rr = r.iloc[0]
                        print(f"  {band:10s} {unit:4s}: {rr.n_pos}/{rr.K_implanted}"
                              f"  med_m={rr.median_m:+.4f} med_z={rr.median_z:+.2f}"
                              f"  perm_p={rr.perm_p:.3f}")
                        summary.append(rr.to_dict())

    if summary:
        pd.DataFrame(summary).to_csv(OUT_ROOT / "summary.csv", index=False)
    print(f"\nOutputs -> {OUT_ROOT}/")


if __name__ == "__main__":
    main()
