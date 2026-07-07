#!/usr/bin/env python3
"""audit_161 — length-matched inference localization null under rho_sym (R2.6 resolver).

audit_160 (rho_sym) showed the beta inference-specific component's ONLY cortical
concentration is the cingulate, clearing the FULL-LENGTH matched-strength null
(BH q=0.030 include / 0.040 exclude). But f = D_taskTest - D_taskLearn inherits the
test-vs-learn recording-length asymmetry, and under rho_split the LENGTH-MATCHED
control (audit_113b) DEMOTED the cingulate to a directional hint. This script repeats
audit_113b under rho_sym: does the cingulate concentration still beat its
strength-matched null AFTER task_test is truncated to task_learn's length?

Faithful re-run of audit_113b with exactly ONE change: the per-pair inference_pe
partial concordance becomes the SYMMETRIC average of both split-half arm assignments
(rho_sym; f is arm-invariant, e & p swap). Truncation, the length-matched surrogate
ensemble (cached lm_head, R=200), and the demeaned aggregation are audit_113b verbatim.
Does NOT edit audit_113b; new script, new output.

5-point preamble
1. Claim: after removing the test/learn length asymmetry, the beta inference-specific
   component still over-accumulates in the cingulate above its strength-matched null
   (BH q<0.05 across systems), under rho_sym.
2. Null: matched-strength surrogate of the TRUNCATED task_test + cached full-length
   surrogates of the 4 unchanged phases (audit_113b null, only the estimator symmetrized).
3. Strongest alternative: the full-length cingulate clearance was carried by task_test's
   extra data sharpening D_TT; at matched length it no longer clears (the rho_split verdict).
4. Does it address it: YES -- observed AND surrogate both built on the truncated task_test,
   so the data-amount advantage is removed from both sides. Cannot rescue power lost to
   truncation (fewer Welch segments -> noisier both sides; a modest q-rise is a power
   statement, not necessarily a confound).
5. Falsify: if the cingulate matched-LENGTH p fails BH (q>=0.05) under rho_sym, R2.6 stays
   a directional lead (duration-assisted). If it clears, the cingulate is an established
   inference location and R2.6/R2.7 upgrade hint->result.

Correctness anchor A1: arm1 (rho_split) matched-length inference_pe M_obs per system
reproduces audit_113b lenmatched_null_R200_{epi}.csv (target=inference_pe).

Output: data/audit/inference_localization_rhosym/lenmatched_null_rhosym_R{R}_{include,exclude}.csv
"""
from __future__ import annotations

import argparse
import gc
import sys

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask
from lrg_eegfc.utils.surrogate import cophenetic_condensed_from_eigs
from lrg_eegfc.utils.surrogate.matched_strength import (
    load_or_compute_eigs_at_path, surrogate_cache_path,
)

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_110_inference_mark_localization import (  # type: ignore
    COHORT, PHASES5, SEED, SWAP, _targets_from_D, load_surrogate_cophenet5,
)
from audit_63_split_baseline_surrogate import ensure_half_fcs, load_phase_fc  # type: ignore
from audit_83_localization_matched_strength import (  # type: ignore
    DROP, _canon_cophenet, unit_means_from_s,
)
from audit_113b_inference_lenmatched_null import _tt_fc_lenmatched, PLACEMENT  # type: ignore

OUT = ROOT / "data/audit/inference_localization_rhosym"
A113B = {  # arm1 anchor (rho_split, R=200)
    "include": ROOT / "data/audit/inference_localization/lenmatched_null_R200_include.csv",
    "exclude": ROOT / "data/audit/inference_localization/lenmatched_null_R200_exclude.csv",
}


def _swap_halves(D: dict) -> dict:
    Ds = dict(D)
    Ds["rest_pre_A"], Ds["rest_pre_B"] = D["rest_pre_B"], D["rest_pre_A"]
    return Ds


def _sym_ipe(D: dict):
    """rho_sym per-pair inference_pe concordance + arm1 (== audit_113b)."""
    s1 = _targets_from_D(D)["inference_pe"]
    s2 = _targets_from_D(_swap_halves(D))["inference_pe"]
    return 0.5 * (s1 + s2), s1


def run_band(band, epis, R, verbose=True):
    obs_um = {e: {} for e in epis}          # sym
    obs_um_a1 = {e: {} for e in epis}       # arm1 (anchor)
    surr_um = {e: {} for e in epis}
    for pat in COHORT:
        surco, miss = load_surrogate_cophenet5(pat, band, R)
        if surco is None:
            if verbose:
                print(f"  [skip] {pat}: full-length surrogate missing ({miss}, R={R})")
            continue
        try:
            W_ttm = _tt_fc_lenmatched(pat, band)
        except FileNotFoundError:
            if verbose:
                print(f"  [skip] {pat}: timeseries missing")
            continue
        D_obs = {ph: _canon_cophenet(load_phase_fc(pat, ph, band)) for ph in PHASES5}
        if W_ttm.shape[0] != int((1 + np.sqrt(1 + 8 * D_obs["task_test"].size)) / 2):
            if verbose:
                print(f"  [skip] {pat}: matched-length N mismatch")
            continue
        D_obs["task_test"] = _canon_cophenet(W_ttm)
        s_sym, s_arm1 = _sym_ipe(D_obs)
        N = W_ttm.shape[0]
        iu = np.triu_indices(N, k=1)
        iu_i, iu_j = iu[0].astype(int), iu[1].astype(int)

        # override task_test surrogate with the cached LENGTH-MATCHED ensemble
        path = surrogate_cache_path(pat, band, f"task_test_lm_{PLACEMENT}", R,
                                    SWAP, SEED, "imcoh_abs")
        rng = np.random.default_rng(SEED)
        evals, evecs = load_or_compute_eigs_at_path(path, W_ttm, R, SWAP, rng,
                                                    verbose=verbose)
        surco["task_test"] = [
            None if not np.all(np.isfinite(evals[r]))
            else cophenetic_condensed_from_eigs(evals[r], evecs[r]) for r in range(R)]
        sr = []
        for r in range(R):
            Dr = {ph: surco[ph][r] for ph in PHASES5}
            sr.append(None if any(Dr[ph] is None for ph in PHASES5)
                      else _sym_ipe(Dr)[0])

        rdf = load_channel_regions(pat)
        uv = rdf["system"].to_numpy().astype(str)
        for epi_excl in epis:
            keep = epi_keep_mask(rdf, pat) if epi_excl else None
            om = unit_means_from_s(s_sym, iu_i, iu_j, uv, keep, DROP["system"])
            om_a1 = unit_means_from_s(s_arm1, iu_i, iu_j, uv, keep, DROP["system"])
            for u, m in om.items():
                obs_um[epi_excl].setdefault(u, {})[pat] = m
            for u, m in om_a1.items():
                obs_um_a1[epi_excl].setdefault(u, {})[pat] = m
            per = {u: [] for u in om}
            for r in range(R):
                if sr[r] is None:
                    continue
                sm = unit_means_from_s(sr[r], iu_i, iu_j, uv, keep, DROP["system"])
                for u in per:
                    if u in sm:
                        per[u].append(sm[u])
            for u in om:
                surr_um[epi_excl].setdefault(u, {})[pat] = np.array(per[u])
        del sr, evals, evecs, surco
        gc.collect()
        if verbose:
            print(f"  [{pat}] matched-length rho_sym done (N={N}, R={R})")

    rows = []
    for epi_excl in epis:
        etag = "exclude" if epi_excl else "include"
        sys_rows = []
        for u, pm in obs_um[epi_excl].items():
            samplers = sorted(pm)
            if not samplers:
                continue
            M_obs = float(np.median([pm[p] for p in samplers]))
            arrs = [surr_um[epi_excl][u].get(p) for p in samplers]
            arrs = [a for a in arrs if a is not None and len(a) > 0]
            if not arrs:
                continue
            Lmin = min(len(a) for a in arrs)
            M_surr = np.median(np.vstack([a[:Lmin] for a in arrs]), axis=0)
            p_ms = (1 + int(np.sum(M_surr >= M_obs))) / (M_surr.size + 1)
            a1 = obs_um_a1[epi_excl].get(u, {})
            M_a1 = (float(np.median([a1[p] for p in samplers if p in a1]))
                    if any(p in a1 for p in samplers) else np.nan)
            row = {"band": band, "epi": etag, "target": "inference_pe", "unit": u,
                   "K_implanted": len(samplers), "M_obs": M_obs, "M_obs_arm1": M_a1,
                   "surr_median": float(np.median(M_surr)),
                   "matched_strength_p": p_ms, "n_surr": M_surr.size,
                   "bh_q_system": np.nan}
            rows.append(row)
            sys_rows.append(row)
        if sys_rows:
            qs = bh_fdr([r["matched_strength_p"] for r in sys_rows])
            for r, q in zip(sys_rows, qs):
                r["bh_q_system"] = q
        if verbose:
            print(f"\n  == {band} inference_pe epi-{etag} (MATCHED LENGTH, rho_sym, R={R}) ==")
            for r in sorted(sys_rows, key=lambda r: r["matched_strength_p"]):
                star = "*" if r["bh_q_system"] < 0.05 else " "
                tag = " <<< CINGULATE" if r["unit"] == "cingulate" else ""
                print(f"    {r['unit']:16s} M_sym={r['M_obs']:+.3e} (arm1={r['M_obs_arm1']:+.3e}) "
                      f"p={r['matched_strength_p']:.4f} q={r['bh_q_system']:.3f}{star}{tag}")
    return rows


def _anchor(df):
    print("\n[anchor A1] arm1 (rho_split) matched-length inference_pe vs audit_113b R=200:")
    worst = 0.0
    for etag in ("include", "exclude"):
        try:
            ref = pd.read_csv(A113B[etag])
        except FileNotFoundError:
            print(f"           [{etag}] audit_113b CSV absent — skip"); continue
        ref = ref[(ref.target == "inference_pe")].set_index("unit")["M_obs"]
        cur = df[df.epi == etag]
        for _, r in cur.iterrows():
            if r.unit in ref.index:
                worst = max(worst, abs(float(ref[r.unit]) - r.M_obs_arm1))
    print(f"           max |Δ| = {worst:.2e}  {'PASS' if worst < 1e-6 else 'CHECK'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default="beta")
    ap.add_argument("--R", type=int, default=200)
    args = ap.parse_args()
    epis = [False, True]
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"[audit_161] pre-flight: ensure rsPre half FCs cached ({args.band})")
    for pat in COHORT:
        ensure_half_fcs(pat, [args.band])
    print(f"\n#### matched-length inference null (rho_sym) : {args.band} R={args.R} "
          f"placement={PLACEMENT} ####")
    rows = run_band(args.band, epis, args.R)
    df = pd.DataFrame(rows)
    for epi_excl in epis:
        etag = "exclude" if epi_excl else "include"
        df[df.epi == etag].to_csv(
            OUT / f"lenmatched_null_rhosym_R{args.R}_{etag}.csv", index=False)
    _anchor(df)
    print(f"\nOutputs -> {OUT}/lenmatched_null_rhosym_R{args.R}_*.csv")


if __name__ == "__main__":
    main()
