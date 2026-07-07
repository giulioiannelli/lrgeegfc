#!/usr/bin/env python3
"""audit_83c — leave-one-patient-out robustness of the consistent-tier carriers.

Closes the one open item in the per-band consistency-taxonomy lock
(`.agents/reports/2026-06-25_per-band-consistency-taxonomy-lock.md`): gamma_l ->
PFC was locked on matched-strength x both-epi x shaft-collapse but had NO
leave-one-patient-out check (audit_83 carries no LOO column), so it stood one
rung below the fully-LOO-confirmed beta -> OFC. This runs LOO for both, REUSING
the audit_83 machinery (imported, not re-implemented): for each left-out patient,
recompute the cohort-median per-system trace over the remaining 9 + the
matched-strength p + BH across systems; the carrier is LOO-robust iff it stays
significant in every drop.

5-point preamble
1. Claim: gamma_l's PFC localization (and beta's OFC) is cohort-consistent, not
   carried by any single patient.
2. Null: the same R=200 matched-strength surrogate, re-aggregated per LOO subset.
3. Alternative it controls: one patient carries the whole cohort-median
   localization; dropping them collapses it.
4. Cannot: drops n=10 -> n=9; tests single-patient leverage only, not 2-out; does
   not re-test trace existence.
5. Falsify: if the carrier fails matched-strength in ANY leave-one-out subset it
   is single-patient-driven and must be downgraded from "consistent".
"""
import sys
from pathlib import Path

import numpy as np

AUDIT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(AUDIT_DIR))
import audit_83_localization_matched_strength as a83  # noqa: E402
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr  # noqa: E402

COHORT, PHASES, DROP, R = a83.COHORT, a83.PHASES, a83.DROP, a83.R
GRAN = "system"


def collect_all(band, epi_excl=False):
    """Per-patient obs + surrogate unit-means for ALL systems (one band)."""
    obs, surr, units = {}, {}, set()
    for pat in COHORT:
        surco = a83.load_surrogate_cophenet(pat, band)
        if surco is None:
            print(f"  [skip] {pat}: surrogate cache missing"); continue
        try:
            s_obs, iu_i, iu_j = a83.obs_trace(pat, band)
        except FileNotFoundError:
            print(f"  [skip] {pat}: half-FC missing"); continue
        rdf = a83.load_channel_regions(pat)
        uv = rdf[GRAN].to_numpy().astype(str)
        keep = a83.epi_keep_mask(rdf, pat) if epi_excl else None
        om = a83.unit_means_from_s(s_obs, iu_i, iu_j, uv, keep, DROP[GRAN])
        obs[pat] = om
        units.update(om.keys())
        per_r = {u: [] for u in om}
        for r in range(R):
            cr = {ph: surco[ph][r] for ph in PHASES}
            if any(cr[ph] is None for ph in PHASES):
                continue
            s_r = a83.concordance(cr["task_test"] - cr["rest_pre_A"],
                                  cr["rest_post"] - cr["rest_pre_B"])
            sm = a83.unit_means_from_s(s_r, iu_i, iu_j, uv, keep, DROP[GRAN])
            for u in per_r:
                if u in sm:
                    per_r[u].append(sm[u])
        surr[pat] = {u: np.asarray(v) for u, v in per_r.items()}
        print(f"  collected {pat} {band}")
    return obs, surr, sorted(units)


def loo_carrier(obs, surr, units, carrier):
    samplers = list(obs.keys())
    rows = []
    for p_out in [None] + samplers:
        keep = [p for p in samplers if p != p_out]
        pvals = {}
        for u in units:
            ks = [p for p in keep if u in obs[p] and len(surr[p].get(u, [])) > 0]
            if not ks:
                continue
            M_obs = float(np.median([obs[p][u] for p in ks]))
            Lmin = min(len(surr[p][u]) for p in ks)
            Rmat = np.vstack([surr[p][u][:Lmin] for p in ks])
            M_surr = np.median(Rmat, axis=0)
            pvals[u] = (1 + int(np.sum(M_surr >= M_obs))) / (M_surr.size + 1)
        us = list(pvals)
        qs = dict(zip(us, bh_fdr(np.array([pvals[u] for u in us]))))
        rows.append({
            "drop": "FULL" if p_out is None else p_out,
            "n": len(keep),
            "carrier_p": pvals.get(carrier, np.nan),
            "carrier_q": qs.get(carrier, np.nan),
            "carrier_is_top": (carrier == min(pvals, key=pvals.get)) if pvals else False,
        })
    return rows


def main():
    import pandas as pd
    out = []
    for band, carrier in [("low_gamma", "PFC"), ("beta", "OFC")]:
        print(f"\n=== {band} -> {carrier} : leave-one-patient-out matched-strength (epi-include, contact) ===")
        obs, surr, units = collect_all(band)
        rows = loo_carrier(obs, surr, units, carrier)
        loo_rows = [r for r in rows if r["drop"] != "FULL"]
        worst_p = max(r["carrier_p"] for r in loo_rows)
        worst_q = max(r["carrier_q"] for r in loo_rows)
        for r in rows:
            tag = r["drop"] if r["drop"] == "FULL" else f"-{r['drop']}"
            flag = "" if r["carrier_q"] < 0.05 else "  <-- FAILS q<.05"
            top = "" if r["carrier_is_top"] else " (not top carrier)"
            print(f"  {tag:9s} (n={r['n']}): p={r['carrier_p']:.4f} q={r['carrier_q']:.4f}{flag}{top}")
            r["band"], r["carrier"] = band, carrier
            out.append(r)
        verdict = "LOO-ROBUST" if worst_q < 0.05 else "NOT robust"
        print(f"  -> worst LOO: p={worst_p:.4f} q={worst_q:.4f}  =>  {band}->{carrier} {verdict}")
    OUT = a83.ROOT / "data/audit/localization_atlas/carrier_loo.csv"
    pd.DataFrame(out).to_csv(OUT, index=False)
    print(f"\nOutputs -> {OUT}")


if __name__ == "__main__":
    main()
