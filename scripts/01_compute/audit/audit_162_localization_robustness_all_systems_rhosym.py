#!/usr/bin/env python3
"""audit_162 — R=1000 + LOO + shaft-collapse robustness for ALL beta systems (rho_sym).

audit_155 put ONLY OFC (upper tail) and sensorimotor/PFC (lower tail) through the full
R=1000 + leave-one-out + shaft-collapse battery. The section-1 localization figure also
shows occipital, cingulate and lateral_temporal glowing enriched at the R=200 tail — but
those were never carried through the full battery, so they were labelled "not
robustness-tested", which is ambiguous. This audit closes that gap: it runs the IDENTICAL
battery (reusing audit_155.collect -> the cached R=1000 seed-20260511 surrogates, no
regeneration) and records the upper- AND lower-tail BH q for EVERY a-priori system, at the
FULL cohort and worst single-patient LOO drop, across all 4 conditions
(epi incl/excl x contact/shaft). Each system is then classified per direction:
    robust_enriched  = worst-LOO upper-tail q < 0.05 in ALL 4 conditions
    robust_depleted  = worst-LOO lower-tail q < 0.05 in ALL 4 conditions
so a system that clears the R=200 tail but fails here is "tested, not robust", not
"untested". Does NOT edit audit_155 / audit_151.

5-point preamble
1. Claim: which beta systems (not just OFC) survive R=1000 + LOO + shaft-collapse?
2. Null: R=1000 matched-strength surrogate, symmetric split-half concordance (== audit_155).
3. Alternative controlled: node strength / single-patient leverage (LOO) / shaft
   autocorrelation (shaft-collapse) / arbitrary A-B half (rho_sym).
4. Cannot: n=10 BH family = a-priori systems within beta; LOO is single-drop (not 2-out);
   does not re-test trace existence; occipital rests on K=3 sparse nodes (small-support
   caveat stands even if it were to clear).
5. Falsify: a system's tail q fails in ANY LOO drop or ANY condition => not robust in that
   direction.

Output: data/audit/localization_atlas_rhosym/beta_all_systems_robustness_R1000.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

AUDIT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(AUDIT_DIR))
import audit_151_localization_rhosym as a151  # noqa: E402
import audit_155_beta_ofc_robustness_rhosym as a155  # noqa: E402

from lrg_eegfc.utils.metrics.hypothesis import bh_fdr  # noqa: E402

CONDITIONS = a155.CONDITIONS
OUT = a151.OUT_ROOT / "beta_all_systems_robustness_R1000.csv"


def loo_full(obs_c, surr_c, units):
    """Per LOO drop -> (q_up, q_lo) dicts over ALL units (BH within each drop)."""
    samplers = list(obs_c.keys())
    per_drop = {}
    for p_out in [None] + samplers:
        keep = [p for p in samplers if p != p_out]
        p_up, p_lo, m_med = {}, {}, {}
        for u in units:
            ks = [p for p in keep if u in obs_c[p] and len(surr_c[p].get(u, [])) > 0]
            if not ks:
                continue
            M_obs = float(np.median([obs_c[p][u] for p in ks]))
            Lmin = min(len(surr_c[p][u]) for p in ks)
            Rmat = np.vstack([surr_c[p][u][:Lmin] for p in ks])
            M_surr = np.median(Rmat, axis=0)
            p_up[u] = (1 + int(np.sum(M_surr >= M_obs))) / (M_surr.size + 1)
            p_lo[u] = (1 + int(np.sum(M_surr <= M_obs))) / (M_surr.size + 1)
            m_med[u] = M_obs
        us = list(p_up)
        q_up = dict(zip(us, bh_fdr(np.array([p_up[u] for u in us]))))
        q_lo = dict(zip(us, bh_fdr(np.array([p_lo[u] for u in us]))))
        per_drop["FULL" if p_out is None else p_out] = (q_up, q_lo, m_med)
    return per_drop


def main():
    band = "beta"
    print(f"[audit_162] all-system robustness under rho_sym: R={a155.R}, LOO, "
          f"shaft-collapse, {len(CONDITIONS)} conditions\n", flush=True)
    obs, surr, units = a155.collect(band)

    rows = []
    for c in CONDITIONS:
        tag = a155._cond_tag(*c)
        per_drop = loo_full(obs[c], surr[c], units)
        full_qup, full_qlo, full_m = per_drop["FULL"]
        loo_keys = [k for k in per_drop if k != "FULL"]
        for u in units:
            worst_qup = max(per_drop[k][0].get(u, np.nan) for k in loo_keys)
            worst_qlo = max(per_drop[k][1].get(u, np.nan) for k in loo_keys)
            rows.append(dict(
                condition=tag, unit=u, M_full=full_m.get(u, np.nan),
                full_q_up=full_qup.get(u, np.nan), worst_loo_q_up=worst_qup,
                full_q_lo=full_qlo.get(u, np.nan), worst_loo_q_lo=worst_qlo))
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)

    # classification summary across the 4 conditions (worst case)
    print(f"\n{'system':16s}{'dir(R200)':>10s}{'full q_up':>10s}"
          f"{'worstLOO up':>12s}{'worstLOO lo':>12s}  robustness")
    a200 = pd.read_csv(a151.OUT_ROOT / "matched_strength_rhosym_include.csv")
    a200 = a200[(a200.band == band) & (a200.granularity == "system")]
    dir200 = {}
    for _, r in a200.iterrows():
        dir200[r.unit] = "enriched" if r.M_obs > r.surr_median else "depleted"
    for u in sorted(set(df.unit)):
        sub = df[df.unit == u]
        # robust in a direction = worst-LOO tail q < 0.05 in EVERY condition
        rob_up = bool((sub.worst_loo_q_up < 0.05).all())
        rob_lo = bool((sub.worst_loo_q_lo < 0.05).all())
        fq_up = float(sub[sub.condition == a155._cond_tag(False, False)].full_q_up.iloc[0])
        w_up = float(sub.worst_loo_q_up.max())
        w_lo = float(sub.worst_loo_q_lo.max())
        d = dir200.get(u, "?")
        rob = ("ROBUST enriched" if (d == "enriched" and rob_up)
               else "ROBUST depleted" if (d == "depleted" and rob_lo)
               else "tested, NOT robust")
        print(f"{u:16s}{d:>10s}{fq_up:>10.4f}{w_up:>12.4f}{w_lo:>12.4f}  {rob}")
    print(f"\n[audit_162] -> {OUT}")
    print("[read] robust_enriched = worst-LOO upper-tail q<0.05 in ALL 4 conditions; "
          "robust_depleted likewise on the lower tail. Anything enriched at R=200 but not "
          "robust here is 'tested, not robust' (was mislabelled 'untested').")


if __name__ == "__main__":
    main()
