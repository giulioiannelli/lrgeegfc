#!/usr/bin/env python3
"""audit_110 — matched-strength VERIFICATION of the cross-phase taxonomy localization.

The decisive test the dissociation needs to be publishable. audit_107 localized
the reset / anchor channels with the GEOMETRY baseline only (count-matched node
relabel) — which controls contact count but NOT node strength, violating the
locked mandatory-matched-strength rule. This script reruns the localization of
ALL channels (trace / reset / anchor, + persist / reset_dom) against the
STRENGTH-PRESERVING surrogate (the same null behind the locked beta->OFC trace),
at R=1000 for beta, with BH-FDR within (band, channel), leave-one-patient-out,
and the shaft-collapse control.

5-point critical preamble
-------------------------
1. Claim: at beta the REVERSIBLE (reset, phi2^2) cross-phase cophenetic flow
   concentrates in LATERAL TEMPORAL cortex, anatomically distinct from the
   PERSISTENT (trace) flow in OFC; the rigid ANCHOR (-var) in insula. These are
   ANATOMICAL, not by-products of those contacts being high-strength nodes.
2. Null: matched-strength 4-cycle +/-delta Laplacian surrogate (audit_63/83/92
   ensemble, R=1000 beta cached), INDEPENDENT per phase, preserving each node's
   strength exactly. Recompute the channel signal on surrogate cophenetic
   distances, run the IDENTICAL per-patient-demeaned per-system aggregation
   (audit_83.unit_means_from_s). If a system's concentration is strength-driven,
   the surrogate reproduces it.
3. Strongest alternative it must control: lateral-temporal / insula contacts are
   systematically higher (or lower) strength; high-strength nodes accumulate more
   excursion / rigidity; the per-system mean inherits a strength->signal gradient.
   The surrogate holds strength FIXED and randomizes topology -> isolates exactly
   this. (This is precisely what the audit_107 geometry baseline did NOT do.)
4. What it cannot do: conditions on observed strengths, so it cannot separate
   "anatomy" from "a strength pattern that is itself anatomically structured" if
   the two are collinear -> then node strength localizes there too (reported as
   strength_dev companion via audit_83). It does not test that reset is "more than
   noise" (reset = the phi2 contrast by construction); it tests LOCALIZATION only.
   EXPECTATION: anchor = -var is near-definitionally tied to strength (rigid =
   stable backbone) and is the most likely to FAIL / be flagged strength-driven.
5. Falsification: if the cohort-median reset unit-mean for lateral_temporal does
   NOT exceed the matched-strength surrogate at p<0.05 (BH within beta), or fails
   LOO-patient, or fails shaft-collapse, the reset->lateral_temporal localization
   is a strength / sampling artifact and is DROPPED (KC-style). trace->OFC MUST
   reproduce (regression guard).

Outputs (data/audit/cross_phase_taxonomy/):
    matched_strength_localization{_shaftcollapsed}_R{R}.csv
      per (band, channel, system): K, M_obs, surr_median, matched_strength_p,
      q_bh, strength_dev_median, and loo_p_max / loo_robust for BH-significant
      systems.

Usage:
    python audit_110_taxonomy_matched_strength.py --bands beta --R 200      # fast check
    python audit_110_taxonomy_matched_strength.py --bands beta --R 1000     # publishable
    python audit_110_taxonomy_matched_strength.py --bands beta --R 1000 --shaft-collapse
"""
from __future__ import annotations

import argparse
import importlib.util

import numpy as np
import pandas as pd

from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask, shaft_of
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.utils.surrogate import (
    cophenetic_condensed_from_eigs, surrogate_cache_path,
)
from lrg_eegfc.workflow.fc import load_fc_matrix

ROOT = setup_script_env()


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "scripts" / "01_compute" / "audit" / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


a83 = _load("audit_83", "audit_83_localization_matched_strength.py")
a107 = _load("audit_107", "audit_107_taxonomy_localization.py")

OUT = ROOT / "data" / "audit" / "cross_phase_taxonomy"
COHORT = a83.COHORT
PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
CHANNELS = ["trace", "persist", "reset", "reset_dom", "anchor"]
SWAP, SEED = 20, 20260511
DROP_SYS = {"non_anatomical"}


# --------------------------------------------------------------------------- #
def obs_signals(pat, band):
    """Observed per-pair channel signals, canonical construction (== audit_83)."""
    Wt = load_fc_matrix(pat, "task_test", band, "imcoh_abs")
    Wpost = load_fc_matrix(pat, "rest_post", band, "imcoh_abs")
    WA = np.load(a83.HALVES_FC_CACHE / pat / f"{band}_rest_pre_A_imcoh_abs.npy")
    WB = np.load(a83.HALVES_FC_CACHE / pat / f"{band}_rest_pre_B_imcoh_abs.npy")
    Dt, Dpost, DA, DB = (a83._canon_cophenet(W) for W in (Wt, Wpost, WA, WB))
    sig = a107.channel_signals(DA, DB, Dt, Dpost)
    N = Wt.shape[0]
    iu = np.triu_indices(N, k=1)
    return sig, iu[0].astype(int), iu[1].astype(int)


def surr_cophenet(pat, band, R):
    """Per-phase surrogate cophenetic condensed distances: dict phase -> [R]."""
    out = {}
    for ph in PHASES:
        path = surrogate_cache_path(pat, band, ph, R, SWAP, SEED, "imcoh_abs")
        if not path.exists():
            return None
        with np.load(path) as d:
            evals, evecs = d["eigvals"], d["eigvecs"]
        out[ph] = [cophenetic_condensed_from_eigs(evals[r], evecs[r])
                   if np.all(np.isfinite(evals[r])) else None
                   for r in range(evals.shape[0])]
    return out


def node_strength(pat, band):
    smat = [load_fc_matrix(pat, ph, band, "imcoh_abs").sum(axis=1)
            for ph in ("rest_pre", "task_test", "rest_post")]
    return np.mean(np.vstack(smat), axis=0)


def _cohort_p(obs_pm, surr_pm, samplers, R):
    """M_obs (cohort median) vs surrogate cohort-median distribution -> p."""
    M_obs = float(np.median([obs_pm[p] for p in samplers]))
    arrs = [surr_pm[p] for p in samplers]
    Lmin = min(len(a) for a in arrs)
    if Lmin == 0:
        return M_obs, np.nan, np.nan
    Rmat = np.vstack([a[:Lmin] for a in arrs])
    M_surr = np.median(Rmat, axis=0)
    p = (1 + int(np.sum(M_surr >= M_obs))) / (M_surr.size + 1)
    return M_obs, p, float(np.median(M_surr))


# --------------------------------------------------------------------------- #
def run_band(band, R, shaft_collapse, verbose=True):
    obs_um = {ch: {} for ch in CHANNELS}       # ch -> {unit: {pat: m}}
    surr_um = {ch: {} for ch in CHANNELS}      # ch -> {unit: {pat: array}}
    strength_um = {}                            # unit -> {pat: strength dev}
    for pat in COHORT:
        surco = surr_cophenet(pat, band, R)
        if surco is None:
            if verbose:
                print(f"  [skip] {pat} {band}: surrogate R={R} cache missing")
            continue
        try:
            obs_sig, iu_i, iu_j = obs_signals(pat, band)
        except FileNotFoundError:
            if verbose:
                print(f"  [skip] {pat} {band}: FC / half-FC missing")
            continue
        rdf = load_channel_regions(pat)
        uv = rdf["system"].to_numpy().astype(str)
        shaft_vec = (np.array([shaft_of(l) for l in rdf["label_raw"]])
                     if shaft_collapse else None)
        # observed unit means per channel
        om = {ch: a83.unit_means_from_s(obs_sig[ch], iu_i, iu_j, uv, None,
                                        DROP_SYS, shaft_vec=shaft_vec)
              for ch in CHANNELS}
        for ch in CHANNELS:
            for u, m in om[ch].items():
                obs_um[ch].setdefault(u, {})[pat] = m
        # surrogate: per realization compute all channels once (shared cophenet)
        per = {ch: {u: [] for u in om[ch]} for ch in CHANNELS}
        for r in range(R):
            cr = {ph: surco[ph][r] for ph in PHASES}
            if any(cr[ph] is None for ph in PHASES):
                continue
            sig_r = a107.channel_signals(cr["rest_pre_A"], cr["rest_pre_B"],
                                         cr["task_test"], cr["rest_post"])
            for ch in CHANNELS:
                sm = a83.unit_means_from_s(sig_r[ch], iu_i, iu_j, uv, None,
                                           DROP_SYS, shaft_vec=shaft_vec)
                for u in per[ch]:
                    if u in sm:
                        per[ch][u].append(sm[u])
        for ch in CHANNELS:
            for u in om[ch]:
                surr_um[ch].setdefault(u, {})[pat] = np.array(per[ch][u])
        # node-strength deviation per system (companion)
        sd = node_strength(pat, band)
        sd = sd - sd.mean()
        for u in np.unique(uv):
            if u in DROP_SYS:
                continue
            sel = uv == u
            if sel.any():
                strength_um.setdefault(u, {})[pat] = float(sd[sel].mean())

    rows = []
    for ch in CHANNELS:
        recs = []
        for u, pm in obs_um[ch].items():
            samplers = sorted(pm)
            if not samplers:
                continue
            M_obs, p, surr_med = _cohort_p(pm, surr_um[ch][u], samplers, R)
            if not np.isfinite(p):
                continue
            strv = strength_um.get(u, {})
            str_med = (float(np.median([strv[p_] for p_ in samplers if p_ in strv]))
                       if any(p_ in strv for p_ in samplers) else np.nan)
            recs.append(dict(band=band, signal=ch, unit=u, K_implanted=len(samplers),
                             M_obs=M_obs, surr_median=surr_med,
                             matched_strength_p=p, strength_dev_median=str_med))
        if not recs:
            continue
        q = bh_fdr(np.array([r["matched_strength_p"] for r in recs]))
        for r, qq in zip(recs, q):
            r["q_bh"] = float(qq)
        # LOO for BH-significant systems
        for r in recs:
            if r["q_bh"] < 0.05:
                u = r["unit"]
                samplers = sorted(obs_um[ch][u])
                loo = []
                for drop in samplers:
                    s2 = [p_ for p_ in samplers if p_ != drop]
                    _, p2, _ = _cohort_p(obs_um[ch][u], surr_um[ch][u], s2, R)
                    loo.append(p2)
                r["loo_p_max"] = float(np.nanmax(loo))
                r["loo_robust"] = bool(np.all(np.array(loo) < 0.05))
            else:
                r["loo_p_max"] = np.nan
                r["loo_robust"] = False
        rows += recs

    if verbose:
        for ch in CHANNELS:
            sub = sorted([r for r in rows if r["signal"] == ch],
                         key=lambda r: r["matched_strength_p"])[:3]
            top = ", ".join(
                f"{r['unit']}(p={r['matched_strength_p']:.3f},q={r['q_bh']:.3f}"
                + (f",LOO≤{r['loo_p_max']:.3f}" if r["q_bh"] < 0.05 else "")
                + f",str{r['strength_dev_median']:+.2f})" for r in sub)
            print(f"  {band:9s} {ch:9s} top: {top}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=["beta"])
    ap.add_argument("--R", type=int, default=200)
    ap.add_argument("--shaft-collapse", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    suffix = ("_shaftcollapsed" if args.shaft_collapse else "") + f"_R{args.R}"
    print(f"\n==== matched-strength taxonomy localization  R={args.R}"
          f"{' SHAFT-COLLAPSED' if args.shaft_collapse else ''} ====")
    allrows = []
    for band in args.bands:
        allrows += run_band(band, args.R, args.shaft_collapse)
    df = pd.DataFrame(allrows)
    out = OUT / f"matched_strength_localization{suffix}.csv"
    df.to_csv(out, index=False)
    print(f"\nOutputs -> {out}")


if __name__ == "__main__":
    main()
