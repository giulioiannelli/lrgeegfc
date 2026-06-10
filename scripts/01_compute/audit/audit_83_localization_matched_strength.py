#!/usr/bin/env python3
"""audit_83 — matched-strength localization control for the audit_81 atlas.

THE mandatory control (feedback_matched_strength_mandatory.md): does the
sampling-conditioned localization of audit_81 (e.g. beta -> MTL/Hip/OFC) survive
a STRENGTH-PRESERVING surrogate null, or is it node-strength geometry? The
label-shuffle null in audit_81 tests "is the trace non-uniform across anatomy
given the per-patient value distribution" but does NOT rule out that high-strength
(central) nodes accumulate more trace and happen to sit in MTL. This control does.

5-point critical preamble
-------------------------
1. Claim: the beta trace's concentration in MTL/Hip/OFC is ANATOMICAL, not a
   by-product of those contacts being high-strength nodes in the FC graph.
2. Null: matched-strength 4-cycle +/-delta Laplacian surrogate (audit_63 ensemble,
   R=200, cached). Per surrogate realization r, reconstruct the cophenetic trace
   s_ij^(r) = dD_task^(r) * dD_rest^(r) on strength-preserving surrogate graphs
   (dD_task = D_coph(task) - D_coph(pre_A); dD_rest = D_coph(post) - D_coph(pre_B);
   exactly the audit_63 split-baseline construction), then run the IDENTICAL
   per-patient-demeaned per-unit aggregation as audit_81. The surrogate preserves
   each node's strength exactly, so if MTL concentration is strength-driven it is
   reproduced under the surrogate; if it is anatomy it is not.
3. Strongest alternative the null covers: MTL/Hip contacts are systematically
   higher (or lower) strength; high-strength nodes have larger |s|; demeaning
   (global, per patient) does not remove a strength->trace gradient. The surrogate
   holds strength FIXED and randomizes topology -> isolates exactly this.
4. What it cannot do: it does not re-test trace existence (locked at 5.2/5.4); it
   conditions on the observed strengths, so it cannot separate "anatomy" from "a
   strength pattern that itself is anatomically structured" if the two are
   perfectly collinear -- but then node strength itself would localize to MTL,
   which is reported as a companion descriptive (strength_localizes column).
5. Falsification: if the real per-unit cohort median M_g does NOT exceed the
   surrogate distribution {M_g^(r)} at p<0.05 for MTL/Hip/OFC at beta, the
   localization is a strength artifact and must be retracted (KC-style). If it
   does, the localization is strength-independent.

Inputs (read-only; surrogate eigs already cached at R=200)
----------------------------------------------------------
- surrogate eigs: data/cache/matched_strength_surrogate_lrg/{pat}/
  {band}_{phase}_R200_swap20_seed20260511_imcoh_abs.npz  (phases: rest_pre_A,
  rest_pre_B, task_test, rest_post), via surrogate_cache_path.
- observed trace: data/reports/imcoh_continuous_trace/per_pair_split/{pat}_{band}.npz
- regions/systems: lrg_eegfc.utils.io.regions.

Outputs (data/audit/localization_atlas/)
-----------------------------------------
- matched_strength_{epi}.csv  per (band, granularity, unit): M_obs, surrogate
  mean/median, matched-strength p (one-sided greater), n_surr_valid, strength_localizes.

Usage
-----
    python audit_83_localization_matched_strength.py            # beta, alpha, low_gamma
    python audit_83_localization_matched_strength.py --band beta --epi-mode include
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.node_localization import (
    NON_ANATOMICAL, epi_keep_mask, shaft_of,
)
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.surrogate import (
    cophenetic_condensed_from_eigs, surrogate_cache_path,
)
from lrg_eegfc.workflow.fc import load_fc_matrix

ROOT = setup_script_env()
OUT_ROOT = ROOT / "data/audit/localization_atlas"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"   # audit_63 split-half FCs

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
TARGET_BANDS = ["beta", "alpha", "low_gamma"]   # the trace bands (R200 cached)
GRANULARITIES = ["region", "system", "supersystem", "paracore"]
PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
R, SWAP, SEED = 200, 20, 20260511
#: paralimbic-ring-minus-OFC: isolates whether the pooled "limbic" survivor is
#: carried by OFC alone or whether MTL+cingulate+insula add independent signal.
PARACORE = {"MTL", "cingulate", "insula"}
DROP = {"region": set(NON_ANATOMICAL), "system": {"non_anatomical"},
        "supersystem": {"non_anatomical"},
        "paracore": {"non_anatomical", "other_np"}}


def load_surrogate_cophenet(pat, band):
    """Per-phase surrogate cophenetic condensed distances: dict phase -> (R, P)."""
    out = {}
    for ph in PHASES:
        path = surrogate_cache_path(pat, band, ph, R, SWAP, SEED, "imcoh_abs")
        if not path.exists():
            return None
        with np.load(path) as d:
            evals, evecs = d["eigvals"], d["eigvecs"]
        cond = []
        for r in range(evals.shape[0]):
            ev = evals[r]
            if not np.all(np.isfinite(ev)):
                cond.append(None)
                continue
            cond.append(cophenetic_condensed_from_eigs(ev, evecs[r]))
        out[ph] = cond
    return out


def concordance(dD_task, dD_rest):
    """Per-pair Spearman-rho_split contribution (RANK-based, scale-invariant).

    c_ij = (rank(dD_task) - (n+1)/2) * (rank(dD_rest) - (n+1)/2); sum(c)/norm =
    rho_split. Rank-based so it is immune to the surrogate cophenetic-magnitude
    explosion (1/rho on near-disconnected rewired graphs) that corrupts the raw
    product, and it is consistent with the LOCKED rho_split = Spearman trace
    measure (feedback_no_partition_metrics_use_rho_coph + audit_73: magnitude
    weighting reintroduces strength structure and fails matched-strength).
    c_ij > 0 = sign-consistent trace direction.
    """
    n = dD_task.size
    mid = (n + 1) / 2.0
    return (rankdata(dD_task) - mid) * (rankdata(dD_rest) - mid)


def _canon_cophenet(W):
    """Canonical LRG cophenetic (condensed) built with the IDENTICAL function
    the surrogate uses (``cophenetic_condensed_from_eigs`` on ``eigh(L)``), so
    observed and surrogate trace are constructed the same way."""
    W = np.asarray(W, dtype=np.float64)
    deg = W.sum(axis=1)
    ev, V = np.linalg.eigh(np.diag(deg) - W)
    return cophenetic_condensed_from_eigs(ev, V)


def obs_trace(pat, band):
    """Observed split-baseline cophenetic trace, recomputed CANONICALLY.

    Built directly from the phase FCs with the same construction as the
    matched-strength surrogate (τ=1/λ_max, T=1/ρ, average-linkage cophenet) —
    NOT read from the LRG-pipeline ``per_pair_split`` cache, whose ultrametric
    differs from the canonical/surrogate one by a per-phase ρ-trace-normalisation
    scale. That scale leaves single-phase ranks intact but perturbs the
    cross-phase DIFFERENCE ranks (Spearman ≈ 0.92–0.99 on dD, 0.85–0.97 on the
    concordance), so reading the pipeline trace would compare a pipeline-scaled
    observed against a canonical surrogate. Recomputing here makes the null
    exact and matches the locked §5.3 audit_63 methodology. (Fix 2026-06-08.)
    """
    Wt = load_fc_matrix(pat, "task_test", band, "imcoh_abs")
    Wpost = load_fc_matrix(pat, "rest_post", band, "imcoh_abs")
    WA = np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_A_imcoh_abs.npy")
    WB = np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_B_imcoh_abs.npy")
    Dt, Dpost, DA, DB = (_canon_cophenet(W) for W in (Wt, Wpost, WA, WB))
    s = concordance(Dt - DA, Dpost - DB)
    N = Wt.shape[0]
    iu = np.triu_indices(N, k=1)
    return s, iu[0].astype(int), iu[1].astype(int)


def unit_means_from_s(s, iu_i, iu_j, unit_vec, keep_node, drop, shaft_vec=None):
    """Per-unit demeaned endpoint mean of a per-pair trace vector s.

    Default: mean over all endpoint-incidences in the unit. With ``shaft_vec``
    (the SHAFT-COLLAPSED control): collapse each electrode shaft to a single
    value first (mean of its incidences in the unit), then average over the
    DISTINCT shafts in the unit. This removes within-shaft pseudoreplication /
    spatial autocorrelation — a unit sampled by 5 contacts on one shaft counts
    as ONE observation, not five — isolating multi-shaft anatomical concentration
    from "one shaft threading the structure". Applied identically to the observed
    trace and to every matched-strength surrogate realization (shaft labels are
    fixed per node, so the collapse is well-defined under the rewired surrogate).
    """
    vals = np.concatenate([s, s])
    nidx = np.concatenate([iu_i, iu_j])
    if keep_node is not None:
        m = keep_node[nidx]
        vals, nidx = vals[m], nidx[m]
    vals = vals - vals.mean()
    u = unit_vec[nidx]
    if shaft_vec is None:
        order = np.argsort(u, kind="stable")
        us, vs = u[order], vals[order]
        uniq, starts = np.unique(us, return_index=True)
        sums = np.add.reduceat(vs, starts)
        counts = np.diff(np.append(starts, vs.size))
        return {un: sm / ct for un, sm, ct in zip(uniq, sums, counts)
                if un not in drop}
    # shaft-collapsed: per (unit, shaft) mean -> mean over distinct shafts/unit
    sh = shaft_vec[nidx]
    key = np.char.add(np.char.add(u.astype(str), "\x1f"), sh.astype(str))
    ukey, inv = np.unique(key, return_inverse=True)
    sums = np.bincount(inv, weights=vals, minlength=ukey.size)
    cnts = np.bincount(inv, minlength=ukey.size).astype(float)
    shaft_means = sums / cnts
    units_of_key = np.array([k.split("\x1f", 1)[0] for k in ukey])
    out = {}
    for un in np.unique(units_of_key):
        if un in drop:
            continue
        out[un] = float(shaft_means[units_of_key == un].mean())
    return out


def run_band(band, epi_excl, verbose=True, shaft_collapse=False):
    rows = []
    # gather per-patient: obs unit-means + surrogate unit-means (per realization)
    obs_um = {g: {} for g in GRANULARITIES}        # gran -> {unit: {pat: m}}
    surr_um = {g: {} for g in GRANULARITIES}        # gran -> {unit: {pat: array(R)}}
    strength_um = {g: {} for g in GRANULARITIES}    # gran -> {unit: {pat: strength-mean}}
    for pat in COHORT:
        surco = load_surrogate_cophenet(pat, band)
        if surco is None:
            if verbose:
                print(f"  [skip] {pat} {band}: surrogate cache missing")
            continue
        try:
            s_obs, iu_i, iu_j = obs_trace(pat, band)
        except FileNotFoundError:
            if verbose:
                print(f"  [skip] {pat} {band}: half-FC missing")
            continue
        rdf = load_channel_regions(pat)
        rdf["paracore"] = rdf["system"].map(
            lambda s: "paralimbic_core" if s in PARACORE
            else ("non_anatomical" if s == "non_anatomical" else "other_np"))
        keep = epi_keep_mask(rdf, pat) if epi_excl else None
        shaft_vec = (np.array([shaft_of(l) for l in rdf["label_raw"]])
                     if shaft_collapse else None)
        # node strength (mean over the 3 main phases) for the companion column
        smat = []
        for ph in ("rest_pre", "task_test", "rest_post"):
            W = load_fc_matrix(pat, ph, band, "imcoh_abs")
            smat.append(W.sum(axis=1))
        s_node = np.mean(np.vstack(smat), axis=0)

        for g in GRANULARITIES:
            uv = rdf[g].to_numpy().astype(str)
            om = unit_means_from_s(s_obs, iu_i, iu_j, uv, keep, DROP[g],
                                   shaft_vec=shaft_vec)
            for u, m in om.items():
                obs_um[g].setdefault(u, {})[pat] = m
            # surrogate per realization
            per_unit_r = {u: [] for u in om}
            for r in range(R):
                cr = {ph: surco[ph][r] for ph in PHASES}
                if any(cr[ph] is None for ph in PHASES):
                    continue
                s_r = concordance(cr["task_test"] - cr["rest_pre_A"],
                                  cr["rest_post"] - cr["rest_pre_B"])
                sm = unit_means_from_s(s_r, iu_i, iu_j, uv, keep, DROP[g],
                                       shaft_vec=shaft_vec)
                for u in per_unit_r:
                    if u in sm:
                        per_unit_r[u].append(sm[u])
            for u in om:
                surr_um[g].setdefault(u, {})[pat] = np.array(per_unit_r[u])
            # strength per unit (demeaned, anatomical-node mean), for companion
            unode = uv
            keepn = np.ones(len(unode), bool) if keep is None else keep
            sd = s_node - s_node[keepn].mean()
            for u in om:
                sel = (unode == u) & keepn
                if sel.any():
                    strength_um[g].setdefault(u, {})[pat] = float(sd[sel].mean())

    # cohort: real median vs surrogate-median distribution
    for g in GRANULARITIES:
        for u, pm in obs_um[g].items():
            samplers = sorted(pm)
            K = len(samplers)
            if K == 0:
                continue
            M_obs = float(np.median([pm[p] for p in samplers]))
            # surrogate cohort median per realization (samplers fixed)
            Rmat = np.vstack([surr_um[g][u][p] for p in samplers
                              if surr_um[g][u].get(p) is not None
                              and len(surr_um[g][u][p]) == R]) \
                if all(len(surr_um[g][u].get(p, [])) == R for p in samplers) \
                else None
            if Rmat is None or Rmat.shape[0] != K:
                # fall back: align on the min common length
                arrs = [surr_um[g][u][p] for p in samplers]
                Lmin = min(len(a) for a in arrs) if arrs else 0
                if Lmin == 0:
                    continue
                Rmat = np.vstack([a[:Lmin] for a in arrs])
            M_surr = np.median(Rmat, axis=0)                 # (R,)
            nval = M_surr.size
            p_ms = (1 + int(np.sum(M_surr >= M_obs))) / (nval + 1)
            strv = strength_um[g].get(u, {})
            str_med = (float(np.median([strv[p] for p in samplers if p in strv]))
                       if any(p in strv for p in samplers) else np.nan)
            rows.append({
                "band": band, "epi": "exclude" if epi_excl else "include",
                "granularity": g, "unit": u, "K_implanted": K,
                "M_obs": M_obs, "surr_mean": float(np.mean(M_surr)),
                "surr_median": float(np.median(M_surr)),
                "matched_strength_p": p_ms, "n_surr": nval,
                "strength_dev_median": str_med,
            })
            if verbose and u in ("MTL", "Hip", "OFC", "limbic", "paralimbic_core"):
                print(f"  {band:10s} {g:6s} {u:5s}: M_obs={M_obs:+.4f} "
                      f"surr_med={np.median(M_surr):+.4f} "
                      f"MS_p={p_ms:.3f}  (K={K}, str_dev={str_med:+.3f})")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default=None)
    ap.add_argument("--epi-mode", choices=["include", "exclude"], default=None)
    ap.add_argument("--shaft-collapse", action="store_true",
                    help="collapse each electrode shaft to one value before "
                         "aggregating (controls within-shaft autocorrelation)")
    ap.add_argument("--R", type=int, default=200,
                    help="surrogate ensemble size to READ (cache must exist; "
                         "R!=200 writes a separate _R{R} CSV so it does not "
                         "clobber the canonical R200 results)")
    args = ap.parse_args()
    global R
    R = args.R
    bands = [args.band] if args.band else TARGET_BANDS
    epis = ([args.epi_mode == "exclude"] if args.epi_mode else [False, True])
    suffix = ("_shaftcollapsed" if args.shaft_collapse else "") \
        + (f"_R{args.R}" if args.R != 200 else "")
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    for epi_excl in epis:
        etag = "exclude" if epi_excl else "include"
        tag = f"epi-{etag}" + (" SHAFT-COLLAPSED" if args.shaft_collapse else "")
        print(f"\n==== matched-strength localization control : {tag} ====")
        allrows = []
        for band in bands:
            allrows += run_band(band, epi_excl,
                                shaft_collapse=args.shaft_collapse)
        pd.DataFrame(allrows).to_csv(
            OUT_ROOT / f"matched_strength{suffix}_{etag}.csv", index=False)
    print(f"\nOutputs -> {OUT_ROOT}/matched_strength{suffix}_*.csv")


if __name__ == "__main__":
    main()
