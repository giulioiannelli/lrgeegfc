#!/usr/bin/env python3
"""audit_155 — beta->OFC localization ROBUSTNESS under rho_sym: R=1000 + LOO + shaft-collapse.

Like-for-like migration of the CANONICAL beta->OFC localization lock (audit_92 R=1000
surrogates + audit_83c leave-one-patient-out + audit_83 shaft-collapse) onto the rho_sym
estimator, so the ANATOMY_LEDGER headline (R=1000 BH q=0.009-0.013, LOO- & shaft-robust,
4/4 conditions) migrates EXACTLY, not just via the R=200 audit_151 confirmation.

Single change vs the canonical pipeline: the per-pair concordance is the SYMMETRIC average
of both split-half arm assignments (a151._sym_concordance), applied identically to observed
and to every cached R=1000 matched-strength surrogate. Reuses the R=1000 seed-20260511
surrogate eig cache (surrogate_cache_path) — no regeneration. Does NOT edit audit_83/92/151.

5-point preamble
1. Claim: beta trace concentrates in OFC (upper-tail) + is depleted in sensorimotor/PFC
   (lower-tail) at R=1000, LOO-robust and shaft-collapse-robust, under rho_sym.
2. Null: R=1000 matched-strength surrogate, per-node aggregated identically, sym concordance.
3. Alternative it controls: node strength (matched-strength) OR one patient carries it (LOO)
   OR electrode-shaft autocorrelation (shaft-collapse) OR the arbitrary A/B half (rho_sym).
4. Cannot: n=10 BH family = a-priori systems within beta; LOO tests single-patient leverage
   (not 2-out); does not re-test trace existence.
5. Falsify: OFC upper-tail q fails in ANY LOO drop or under shaft-collapse => estimator-,
   patient-, or shaft-dependent.

Output: data/audit/localization_atlas_rhosym/beta_ofc_robustness_R1000.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

AUDIT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(AUDIT_DIR))
import audit_151_localization_rhosym as a151  # noqa: E402

from lrg_eegfc.utils.io.regions import load_channel_regions  # noqa: E402
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask, shaft_of  # noqa: E402
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr  # noqa: E402
from lrg_eegfc.utils.surrogate import (  # noqa: E402
    cophenetic_condensed_from_eigs, surrogate_cache_path,
)

COHORT, PHASES = a151.COHORT, a151.PHASES
SWAP, SEED = a151.SWAP, a151.SEED
R = 1000
GRAN = "system"
DROP_SYS = a151.DROP["system"]
CARRIER = "OFC"
DEPLETED = ("sensorimotor", "PFC")
# (epi_excl, shaft_collapse) — the 4 conditions the ANATOMY_LEDGER cites
CONDITIONS = [(False, False), (True, False), (False, True), (True, True)]
OUT = a151.OUT_ROOT / "beta_ofc_robustness_R1000.csv"


def _cond_tag(epi_excl, shaft):
    return f"epi-{'excl' if epi_excl else 'incl'}/{'shaft' if shaft else 'contact'}"


def load_surr_R1000(pat, band):
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
            cond.append(None if not np.all(np.isfinite(ev))
                        else cophenetic_condensed_from_eigs(ev, evecs[r]))
        out[ph] = cond
    return out


def collect(band):
    """obs[cond][pat] = {unit: M}; surr[cond][pat] = {unit: array}. One R-loop/patient."""
    obs = {c: {} for c in CONDITIONS}
    surr = {c: {} for c in CONDITIONS}
    units = set()
    for pat in COHORT:
        surco = load_surr_R1000(pat, band)
        if surco is None:
            print(f"  [skip] {pat}: R=1000 cache missing", flush=True)
            continue
        try:
            s_obs, iu_i, iu_j = a151.obs_trace_sym(pat, band)
        except FileNotFoundError:
            print(f"  [skip] {pat}: half-FC missing", flush=True)
            continue
        rdf = load_channel_regions(pat)
        uv = rdf[GRAN].to_numpy().astype(str)
        keep_epi = epi_keep_mask(rdf, pat)
        shaft_full = np.array([shaft_of(l) for l in rdf["label_raw"]])
        cfg = {}
        for (epi_excl, shaft) in CONDITIONS:
            keep = keep_epi if epi_excl else None
            sv = shaft_full if shaft else None
            om = a151.unit_means_from_s(s_obs, iu_i, iu_j, uv, keep, DROP_SYS, shaft_vec=sv)
            obs[(epi_excl, shaft)][pat] = om
            units.update(om)
            cfg[(epi_excl, shaft)] = (keep, sv, {u: [] for u in om})
        for r in range(R):
            if any(surco[ph][r] is None for ph in PHASES):
                continue
            s_r = a151._sym_concordance(surco["rest_pre_A"][r], surco["rest_pre_B"][r],
                                        surco["task_test"][r], surco["rest_post"][r])
            for (keep, sv, acc) in cfg.values():
                sm = a151.unit_means_from_s(s_r, iu_i, iu_j, uv, keep, DROP_SYS, shaft_vec=sv)
                for u in acc:
                    if u in sm:
                        acc[u].append(sm[u])
        for c, (keep, sv, acc) in cfg.items():
            surr[c][pat] = {u: np.asarray(v) for u, v in acc.items()}
        print(f"  collected {pat} {band}", flush=True)
    return obs, surr, sorted(units)


def loo(obs_c, surr_c, units):
    """LOO over patients: cohort-median per-system upper/lower-tail p + BH per drop."""
    samplers = list(obs_c.keys())
    rows = []
    for p_out in [None] + samplers:
        keep = [p for p in samplers if p != p_out]
        p_up, p_lo = {}, {}
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
        us = list(p_up)
        q_up = dict(zip(us, bh_fdr(np.array([p_up[u] for u in us]))))
        q_lo = dict(zip(us, bh_fdr(np.array([p_lo[u] for u in us]))))
        row = {"drop": "FULL" if p_out is None else p_out, "n": len(keep),
               "OFC_M": float(np.median([obs_c[p][CARRIER] for p in keep if CARRIER in obs_c[p]])),
               "OFC_p_up": p_up.get(CARRIER, np.nan), "OFC_q_up": q_up.get(CARRIER, np.nan),
               "OFC_is_top": (CARRIER == min(p_up, key=p_up.get)) if p_up else False}
        for d in DEPLETED:
            row[f"{d}_q_lo"] = q_lo.get(d, np.nan)
        rows.append(row)
    return rows


def main():
    band = "beta"
    print(f"[audit_155] beta->OFC robustness under rho_sym: R={R}, LOO, shaft-collapse, "
          f"{len(CONDITIONS)} conditions\n")
    obs, surr, units = collect(band)
    all_rows = []
    print(f"\n{'condition':22s}{'OFC full q_up':>14}{'worst-LOO q_up':>16}"
          f"{'OFC always top':>15}{'sens q_lo':>11}{'PFC q_lo':>10}  verdict")
    for c in CONDITIONS:
        rows = loo(obs[c], surr[c], units)
        for r in rows:
            r["condition"] = _cond_tag(*c)
            all_rows.append(r)
        full = next(r for r in rows if r["drop"] == "FULL")
        loo_rows = [r for r in rows if r["drop"] != "FULL"]
        worst_q = max(r["OFC_q_up"] for r in loo_rows)
        always_top = all(r["OFC_is_top"] for r in rows)
        sens_q = full.get("sensorimotor_q_lo", np.nan)
        pfc_q = full.get("PFC_q_lo", np.nan)
        verdict = "LOO-ROBUST carrier" if worst_q < 0.05 else "NOT LOO-robust"
        print(f"{_cond_tag(*c):22s}{full['OFC_q_up']:>14.4f}{worst_q:>16.4f}"
              f"{str(always_top):>15}{sens_q:>11.4f}{pfc_q:>10.4f}  {verdict}")
    df = pd.DataFrame(all_rows)
    df.to_csv(OUT, index=False)
    print(f"\n[audit_155] outputs -> {OUT}")
    print("[read] OFC_q_up = upper-tail BH q (carrier); *_q_lo = lower-tail BH q (depleted). "
          "OFC LOO-robust across all 4 conditions => estimator-, patient-, shaft-invariant.")


if __name__ == "__main__":
    main()
