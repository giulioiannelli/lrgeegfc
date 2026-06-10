#!/usr/bin/env python3
"""audit_84 — strength-residualization control for the localization atlas.

Companion to audit_83 (matched-strength surrogate). The matched-strength control
showed beta OFC clears (MS_p=0.005) while Hip/MTL do not (MS_p~0.10). The natural
INTERPRETATION "Hip fails because it is high-strength hub geometry" must be tested
DIRECTLY, not assumed -- a Nature reviewer (and our own adversarial review,
2026-06-05) flagged exactly this. This script regresses the per-node trace
concordance on per-node FC strength WITHIN each patient and re-aggregates on the
residuals: if Hip's concentration vanishes after removing the linear
strength->trace component, it was strength geometry; if it persists, Hip's
sub-threshold matched-strength result is an effect-size/null-width phenomenon, NOT
a strength confound.

5-point critical preamble
-------------------------
1. Claim: Hip/MTL's failure to clear matched-strength is due to node-strength
   (hubness), i.e. the concentration is strength-explained.
2. Null/control: within-patient OLS residualization of per-node concordance on
   per-node strength (mean over rest_pre/task/rest_post). Re-aggregate the
   residuals per unit (per-patient demeaned, cohort median over implanted patients).
3. Strongest alternative: strength and trace are collinear within patient, so the
   raw per-unit concentration is a strength projection.
4. Coverage: if Spearman(strength, concordance) is near zero within patients AND
   residualized medians ~ raw medians, the strength-explanation is refuted by
   mechanism (there is nothing to regress out). It cannot prove anatomy, only
   remove strength as the explanation.
5. Falsification: if Hip's residualized median collapses toward 0 (or below the
   surrogate band) while OFC's persists, the "Hip = strength geometry" reading is
   SUPPORTED. If Hip persists, it is REFUTED -> Hip is a strength-independent but
   sub-threshold localization, not an artifact.

Inputs (read-only): observed trace per_pair_split + load_fc_matrix strengths.
Outputs: data/audit/localization_atlas/strength_residual_{epi}.csv
    per (band, granularity, unit): raw_median, resid_median, strength_trace_rho.

Usage
-----
    python audit_84_localization_strength_residual.py
    python audit_84_localization_strength_residual.py --band beta --epi-mode include
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.node_localization import NON_ANATOMICAL, epi_keep_mask
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.fc import load_fc_matrix

ROOT = setup_script_env()
OUT_ROOT = ROOT / "data/audit/localization_atlas"
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
TARGET_BANDS = ["beta", "alpha", "low_gamma"]
GRANULARITIES = ["region", "system"]
DROP = {"region": set(NON_ANATOMICAL), "system": {"non_anatomical"}}


def concordance(dD_task, dD_rest):
    n = dD_task.size
    mid = (n + 1) / 2.0
    return (rankdata(dD_task) - mid) * (rankdata(dD_rest) - mid)


def per_node_concordance(pat, band, n_nodes, keep_node):
    """Per-node mean concordance contribution (endpoint-aggregated)."""
    d = np.load(ROOT / "data/reports/imcoh_continuous_trace/per_pair_split"
                / f"{pat}_{band}.npz")
    c = concordance(np.asarray(d["dD_task"], float), np.asarray(d["dD_rest"], float))
    iu_i, iu_j = np.asarray(d["iu_i"]).astype(int), np.asarray(d["iu_j"]).astype(int)
    vals = np.concatenate([c, c])
    nidx = np.concatenate([iu_i, iu_j])
    if keep_node is not None:
        m = keep_node[nidx]
        vals, nidx = vals[m], nidx[m]
    s = np.bincount(nidx, weights=vals, minlength=n_nodes)
    cnt = np.bincount(nidx, minlength=n_nodes).astype(float)
    with np.errstate(invalid="ignore"):
        node_c = np.where(cnt > 0, s / cnt, np.nan)
    return node_c, cnt > 0


def node_strength(pat, band):
    smat = []
    for ph in ("rest_pre", "task_test", "rest_post"):
        W = load_fc_matrix(pat, ph, band, "imcoh_abs")
        smat.append(W.sum(axis=1))
    return np.mean(np.vstack(smat), axis=0)


def run(band, epi_excl, verbose=True):
    # per patient: node concordance, node strength, residual concordance, unit map
    rows_raw = {g: {} for g in GRANULARITIES}     # unit -> {pat: raw median}
    rows_res = {g: {} for g in GRANULARITIES}     # unit -> {pat: resid median}
    rhos = []
    for pat in COHORT:
        rdf = load_channel_regions(pat)
        n = len(rdf)
        keep = epi_keep_mask(rdf, pat) if epi_excl else None
        node_c, sampled = per_node_concordance(pat, band, n, keep)
        s_node = node_strength(pat, band)
        valid = sampled & np.isfinite(node_c) & np.isfinite(s_node)
        if keep is not None:
            valid &= keep
        cc, ss = node_c[valid], s_node[valid]
        if cc.size < 3:
            continue
        # within-patient OLS residual of concordance on strength
        b, a = np.polyfit(ss, cc, 1)
        resid_full = np.full(n, np.nan)
        resid_full[valid] = cc - (a + b * ss)
        cmean = cc.mean()
        c_full = np.full(n, np.nan)
        c_full[valid] = cc
        rho, _ = spearmanr(ss, cc)
        rhos.append(rho)
        for g in GRANULARITIES:
            uv = rdf[g if g == "system" else "region"].to_numpy().astype(str)
            # demean per patient (over valid nodes) so it is a localization contrast
            craw = c_full - np.nanmean(c_full)
            cres = resid_full - np.nanmean(resid_full)
            for u in np.unique(uv[valid]):
                if u in DROP[g]:
                    continue
                sel = valid & (uv == u)
                if sel.any():
                    rows_raw[g].setdefault(u, {})[pat] = float(np.nanmean(craw[sel]))
                    rows_res[g].setdefault(u, {})[pat] = float(np.nanmean(cres[sel]))
    out = []
    for g in GRANULARITIES:
        for u, pm in rows_raw[g].items():
            samplers = sorted(pm)
            K = len(samplers)
            raw_med = float(np.median([pm[p] for p in samplers]))
            res_med = float(np.median([rows_res[g][u][p] for p in samplers]))
            out.append({
                "band": band, "epi": "exclude" if epi_excl else "include",
                "granularity": g, "unit": u, "K_implanted": K,
                "raw_median": raw_med, "resid_median": res_med,
                "resid_over_raw": (res_med / raw_med) if raw_med else np.nan,
            })
    if verbose:
        med_rho = float(np.median(rhos)) if rhos else np.nan
        print(f"  {band}: within-patient median Spearman(strength, concordance) "
              f"= {med_rho:+.3f}")
        for g, u in (("system", "OFC"), ("region", "Hip"), ("system", "MTL")):
            r = [x for x in out if x["granularity"] == g and x["unit"] == u]
            if r:
                rr = r[0]
                print(f"    {g:6s} {u:5s}: raw={rr['raw_median']:+.0f} "
                      f"resid={rr['resid_median']:+.0f} "
                      f"(resid/raw={rr['resid_over_raw']:+.2f})")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default=None)
    ap.add_argument("--epi-mode", choices=["include", "exclude"], default=None)
    args = ap.parse_args()
    bands = [args.band] if args.band else TARGET_BANDS
    epis = ([args.epi_mode == "exclude"] if args.epi_mode else [False, True])
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    for epi_excl in epis:
        etag = "exclude" if epi_excl else "include"
        print(f"\n==== strength-residualization control : epi-{etag} ====")
        allrows = []
        for band in bands:
            allrows += run(band, epi_excl)
        pd.DataFrame(allrows).to_csv(
            OUT_ROOT / f"strength_residual_{etag}.csv", index=False)
    print(f"\nOutputs -> {OUT_ROOT}/strength_residual_*.csv")


if __name__ == "__main__":
    main()
