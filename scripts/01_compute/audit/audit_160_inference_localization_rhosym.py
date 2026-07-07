#!/usr/bin/env python3
"""audit_160 — inference-specific (f | e) localization under rho_sym (R2.6).

Faithful re-run of audit_110's INFERENCE-MARK localizer (the demeaned per-system
concentration of the beta inference-specific persistent trace) with exactly ONE
change: the per-pair partial concordance becomes the SYMMETRIC average of both
split-half arm assignments (rho_sym = 1/2[arm_AB + arm_BA]).

The inference target is per-pair  concordance_partial(f, p | e), with
  e = D_taskLearn - D_preX (encoding),  f = D_taskTest - D_taskLearn (arm-invariant),
  p = D_restPost  - D_preY (persistence).
Arm1 (audit_110): X=A, Y=B.  Arm2: swap the two rest_pre halves (X=B, Y=A); f is
unchanged. Implemented by calling audit_110._targets_from_D on the phase dict and on
the same dict with rest_pre_A<->rest_pre_B swapped, then averaging the inference_pe
per-pair vectors. Aggregation to systems is audit_151's demeaned unit_means_from_s
(the same aggregation that reproduced audit_110's encoding M_obs bit-for-bit in
audit_158). Reuses the seed-20260511 R=1000 surrogate cophenetics (5 phases); NO
regeneration. Does NOT edit audit_110; new script, new output.

CONTEXT / ceiling: inference->cingulate is a DIRECTIONAL HINT in the rho_split era
(audit_110 BH-borderline at full length; DOWNGRADED by the matched-length control,
audit_113b). This migration asks only whether the rho_sym localizer keeps the same
DIRECTION (cingulate as the leading inference concentration); it does not upgrade a
hint to a result. Report the sign/lead + MS_p/BH honestly, significant or not.

5-point preamble
1. Claim: the beta inference-specific concentration LEANS toward cingulate under rho_sym
   (same direction as rho_split), not orbitofrontal.
2. Null: the cingulate lean was an artifact of the arbitrary A/B arm.
3. Strongest alternative: node strength. Controlled by the same matched-strength
   surrogate (strength-preserving); estimator change is orthogonal to strength.
4. Cannot: R=1000 cached surrogate; BH over a-priori systems within beta x inference_pe.
   Does NOT establish a hotspot (the matched-length duration control already demoted it);
   this is direction + magnitude under the corrected estimator only.
5. Falsify: if cingulate is no longer the leading (top upper-tail) inference system under
   rho_sym, the directional hint does not survive the estimator fix.

Correctness anchor A1: arm1 (rho_split) inference_pe M_obs per beta/system reproduces
audit_110 inference_mark_R1000_{epi}.csv (target=inference_pe, granularity=system).

Output: data/audit/inference_localization_rhosym/inference_pe_localization_rhosym_{include,exclude}.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

import sys
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.io.regions import load_channel_regions
from audit_151_localization_rhosym import (  # type: ignore
    unit_means_from_s, PARACORE, DROP, GRANULARITIES,
)
from audit_110_inference_mark_localization import (  # type: ignore
    _targets_from_D, load_surrogate_cophenet5, _canon_cophenet, load_phase_fc,
    COHORT, PHASES5,
)

BAND = "beta"
R, SEED = 1000, 20260511
OUT_ROOT = ROOT / "data/audit/inference_localization_rhosym"
A110 = {  # arm1 anchor (rho_split, R=1000)
    "include": ROOT / "data/audit/inference_localization/inference_mark_R1000_include.csv",
    "exclude": ROOT / "data/audit/inference_localization/inference_mark_R1000_exclude.csv",
}


def _swap_halves(D: dict) -> dict:
    Ds = dict(D)
    Ds["rest_pre_A"], Ds["rest_pre_B"] = D["rest_pre_B"], D["rest_pre_A"]
    return Ds


def _sym_inference_pe(D: dict):
    """rho_sym per-pair inference_pe concordance + arm1 (== audit_110)."""
    s1 = _targets_from_D(D)["inference_pe"]
    s2 = _targets_from_D(_swap_halves(D))["inference_pe"]
    return 0.5 * (s1 + s2), s1


def run(epi_excl, verbose=True):
    rows = []
    obs_um = {g: {} for g in GRANULARITIES}
    obs_um_arm1 = {g: {} for g in GRANULARITIES}
    surr_um = {g: {} for g in GRANULARITIES}
    strength_um = {g: {} for g in GRANULARITIES}
    for pat in COHORT:
        surco, miss = load_surrogate_cophenet5(pat, BAND, R)
        if surco is None:
            if verbose:
                print(f"  [skip] {pat}: surrogate missing ({miss}, R={R})")
            continue
        try:
            D = {ph: _canon_cophenet(load_phase_fc(pat, ph, BAND)) for ph in PHASES5}
        except FileNotFoundError:
            continue
        s_sym, s_arm1 = _sym_inference_pe(D)
        N = int((1 + np.sqrt(1 + 8 * D["task_test"].size)) / 2)
        iu = np.triu_indices(N, k=1)
        iu_i, iu_j = iu[0].astype(int), iu[1].astype(int)
        rdf = load_channel_regions(pat)
        rdf["paracore"] = rdf["system"].map(
            lambda s: "paralimbic_core" if s in PARACORE
            else ("non_anatomical" if s == "non_anatomical" else "other_np"))
        keep = epi_keep_mask(rdf, pat) if epi_excl else None
        s_node = np.mean(np.vstack(
            [load_phase_fc(pat, ph, BAND).sum(axis=1)
             for ph in ("rest_pre", "task_test", "rest_post")]), axis=0)
        # surrogate per-pair sym inference_pe per realization
        s_surr = []
        for r in range(R):
            Dr = {ph: surco[ph][r] for ph in PHASES5}
            s_surr.append(None if any(Dr[ph] is None for ph in PHASES5)
                          else _sym_inference_pe(Dr)[0])
        for g in GRANULARITIES:
            uv = rdf[g].to_numpy().astype(str)
            om = unit_means_from_s(s_sym, iu_i, iu_j, uv, keep, DROP[g])
            om_a1 = unit_means_from_s(s_arm1, iu_i, iu_j, uv, keep, DROP[g])
            for u, mval in om.items():
                obs_um[g].setdefault(u, {})[pat] = mval
            for u, mval in om_a1.items():
                obs_um_arm1[g].setdefault(u, {})[pat] = mval
            per_unit_r = {u: [] for u in om}
            for sr in s_surr:
                if sr is None:
                    continue
                sm = unit_means_from_s(sr, iu_i, iu_j, uv, keep, DROP[g])
                for u in per_unit_r:
                    if u in sm:
                        per_unit_r[u].append(sm[u])
            for u in om:
                surr_um[g].setdefault(u, {})[pat] = np.array(per_unit_r[u])
            keepn = np.ones(len(uv), bool) if keep is None else keep
            sd = s_node - s_node[keepn].mean()
            for u in om:
                sel = (uv == u) & keepn
                if sel.any():
                    strength_um[g].setdefault(u, {})[pat] = float(sd[sel].mean())
    for g in GRANULARITIES:
        for u, pm in obs_um[g].items():
            samplers = sorted(pm)
            K = len(samplers)
            if K == 0:
                continue
            M_obs = float(np.median([pm[p] for p in samplers]))
            arrs = [surr_um[g][u][p] for p in samplers]
            Lmin = min((len(a) for a in arrs), default=0)
            if Lmin == 0:
                continue
            Rmat = np.vstack([a[:Lmin] for a in arrs])
            M_surr = np.median(Rmat, axis=0)
            nval = M_surr.size
            p_up = (1 + int(np.sum(M_surr >= M_obs))) / (nval + 1)
            p_lo = (1 + int(np.sum(M_surr <= M_obs))) / (nval + 1)
            a1 = obs_um_arm1[g].get(u, {})
            M_arm1 = (float(np.median([a1[p] for p in samplers if p in a1]))
                      if any(p in a1 for p in samplers) else np.nan)
            strv = strength_um[g].get(u, {})
            str_med = (float(np.median([strv[p] for p in samplers if p in strv]))
                       if any(p in strv for p in samplers) else np.nan)
            rows.append({
                "band": BAND, "epi": "exclude" if epi_excl else "include",
                "target": "inference_pe", "granularity": g, "unit": u,
                "K_implanted": K, "M_obs": M_obs, "M_obs_arm1": M_arm1,
                "surr_median": float(np.median(M_surr)),
                "matched_strength_p": p_up, "matched_strength_p_lower": p_lo,
                "n_surr": nval, "strength_dev_median": str_med})
    return rows


def _anchor(df):
    print("\n[anchor A1] arm1 (rho_split) inference_pe M_obs vs audit_110 R=1000:")
    worst = 0.0
    for etag in ("include", "exclude"):
        ref = pd.read_csv(A110[etag])
        ref = ref[(ref.band == "beta") & (ref.target == "inference_pe")
                  & (ref.granularity == "system")].set_index("unit")["M_obs"]
        cur = df[(df.epi == etag) & (df.granularity == "system")]
        for _, r in cur.iterrows():
            if r.unit in ref.index:
                worst = max(worst, abs(float(ref[r.unit]) - r.M_obs_arm1))
    print(f"           max |Δ| = {worst:.2e}  {'PASS' if worst < 1e-9 else 'CHECK'}")


def main():
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    allrows = []
    for epi_excl in (False, True):
        etag = "exclude" if epi_excl else "include"
        print(f"\n==== inference_pe localization (rho_sym) : epi-{etag} (R={R}) ====")
        rows = run(epi_excl)
        pd.DataFrame(rows).to_csv(
            OUT_ROOT / f"inference_pe_localization_rhosym_{etag}.csv", index=False)
        allrows += rows
        sysdf = pd.DataFrame([r for r in rows if r["granularity"] == "system"])
        if not sysdf.empty:
            sysdf["q_up"] = bh_fdr(sysdf["matched_strength_p"].values)
            print(f"  -- beta inference_pe systems, epi-{etag}, BH over {len(sysdf)} --")
            for _, r in sysdf.sort_values("matched_strength_p").iterrows():
                tag = "LEAD" if r.q_up < 0.05 else ""
                print(f"     {r.unit:16s} M={r.M_obs:+.4f} (arm1={r.M_obs_arm1:+.4f}) "
                      f"p_up={r.matched_strength_p:.3f} q_up={r.q_up:.3f} "
                      f"strdev={r.strength_dev_median:+.2f} {tag}")
    _anchor(pd.DataFrame(allrows))
    print(f"\nOutputs -> {OUT_ROOT}/inference_pe_localization_rhosym_*.csv")


if __name__ == "__main__":
    main()
