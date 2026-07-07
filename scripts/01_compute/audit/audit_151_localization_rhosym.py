#!/usr/bin/env python3
"""audit_151 — beta -> OFC localization under the rho_sym estimator.

The flagship's SECOND pillar (after the audit_150 gate): does the beta cophenetic
trace still CONCENTRATE in OFC when the per-node decomposition is built from the
arbitrary-half-free rho_sym estimator instead of bare rho_split?

Faithful re-run of audit_83's matched-strength localization with exactly ONE change:
the per-pair rank concordance becomes the SYMMETRIC average of both split-half arm
assignments (mirrors rho_sym = 1/2[rho_AB + rho_BA]):

    s_sym = 1/2 [ concordance(D_task-D_preA, D_post-D_preB)
                + concordance(D_task-D_preB, D_post-D_preA) ]

applied IDENTICALLY to the observed trace and to every cached matched-strength
surrogate realization. Reuses audit_83's cached surrogate cophenetic eigs (no
regeneration). Does NOT edit audit_83; new script, new output dir.

5-point preamble
1. Claim: beta trace concentrates in OFC (upper-tail matched-strength) and is
   DEPLETED in sensorimotor (lower-tail) under rho_sym, as under rho_split.
2. Null: the OFC concentration was an artifact of the arbitrary A/B arm choice.
3. Strongest alternative: it is node strength, not the estimator. Controlled by the
   same matched-strength surrogate (strength-preserving) audit_83 uses; the
   estimator change is orthogonal to strength.
4. Cannot: reuses R=200 cached surrogate (same as canonical audit_83 lock); n=10
   BH family = a-priori systems within beta. Shaft-collapse/LOO are separate runs.
5. Falsify: if OFC upper-tail MS_p (or its BH q over systems) fails, or sensorimotor
   lower-tail fails, the localization is estimator-dependent.

Output: data/audit/localization_atlas_rhosym/matched_strength_rhosym_{include,exclude}.csv
"""
from __future__ import annotations
import numpy as np, pandas as pd

from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.node_localization import (
    NON_ANATOMICAL, canonical_cophenet, epi_keep_mask, rank_concordance, shaft_of,
)
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.surrogate import cophenetic_condensed_from_eigs, surrogate_cache_path
from lrg_eegfc.workflow.fc import load_fc_matrix

ROOT = setup_script_env()
OUT_ROOT = ROOT / "data/audit/localization_atlas_rhosym"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
GRANULARITIES = ["region", "system", "supersystem", "paracore"]
PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
R, SWAP, SEED = 200, 20, 20260511
PARACORE = {"MTL", "cingulate", "insula"}
DROP = {"region": set(NON_ANATOMICAL), "system": {"non_anatomical"},
        "supersystem": {"non_anatomical"},
        "paracore": {"non_anatomical", "other_np"}}


def _sym_concordance(cA, cB, ctt, cpost):
    """rho_sym per-pair concordance: mean of both arm assignments."""
    s_ab = rank_concordance(ctt - cA, cpost - cB)
    s_ba = rank_concordance(ctt - cB, cpost - cA)
    return 0.5 * (s_ab + s_ba)


def obs_trace_sym(pat, band):
    Wt = load_fc_matrix(pat, "task_test", band, "imcoh_abs")
    Wpost = load_fc_matrix(pat, "rest_post", band, "imcoh_abs")
    WA = np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_A_imcoh_abs.npy")
    WB = np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_B_imcoh_abs.npy")
    Dt, Dpost, DA, DB = (canonical_cophenet(W) for W in (Wt, Wpost, WA, WB))
    s = _sym_concordance(DA, DB, Dt, Dpost)
    N = Wt.shape[0]; iu = np.triu_indices(N, k=1)
    return s, iu[0].astype(int), iu[1].astype(int)


def load_surrogate_cophenet(pat, band):
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


def unit_means_from_s(s, iu_i, iu_j, unit_vec, keep_node, drop, shaft_vec=None):
    vals = np.concatenate([s, s]); nidx = np.concatenate([iu_i, iu_j])
    if keep_node is not None:
        m = keep_node[nidx]; vals, nidx = vals[m], nidx[m]
    vals = vals - vals.mean(); u = unit_vec[nidx]
    if shaft_vec is None:
        order = np.argsort(u, kind="stable"); us, vs = u[order], vals[order]
        uniq, starts = np.unique(us, return_index=True)
        sums = np.add.reduceat(vs, starts)
        counts = np.diff(np.append(starts, vs.size))
        return {un: sm / ct for un, sm, ct in zip(uniq, sums, counts) if un not in drop}
    sh = shaft_vec[nidx]
    key = np.char.add(np.char.add(u.astype(str), "\x1f"), sh.astype(str))
    ukey, inv = np.unique(key, return_inverse=True)
    sums = np.bincount(inv, weights=vals, minlength=ukey.size)
    cnts = np.bincount(inv, minlength=ukey.size).astype(float)
    shaft_means = sums / cnts
    units_of_key = np.array([k.split("\x1f", 1)[0] for k in ukey])
    out = {}
    for un in np.unique(units_of_key):
        if un not in drop:
            out[un] = float(shaft_means[units_of_key == un].mean())
    return out


def run_band(band, epi_excl, shaft_collapse=False, verbose=True):
    rows = []
    obs_um = {g: {} for g in GRANULARITIES}
    surr_um = {g: {} for g in GRANULARITIES}
    strength_um = {g: {} for g in GRANULARITIES}
    for pat in COHORT:
        surco = load_surrogate_cophenet(pat, band)
        if surco is None:
            if verbose: print(f"  [skip] {pat} {band}: surrogate cache missing")
            continue
        try:
            s_obs, iu_i, iu_j = obs_trace_sym(pat, band)
        except FileNotFoundError:
            if verbose: print(f"  [skip] {pat} {band}: half-FC missing")
            continue
        rdf = load_channel_regions(pat)
        rdf["paracore"] = rdf["system"].map(
            lambda s: "paralimbic_core" if s in PARACORE
            else ("non_anatomical" if s == "non_anatomical" else "other_np"))
        keep = epi_keep_mask(rdf, pat) if epi_excl else None
        shaft_vec = (np.array([shaft_of(l) for l in rdf["label_raw"]])
                     if shaft_collapse else None)
        smat = []
        for ph in ("rest_pre", "task_test", "rest_post"):
            smat.append(load_fc_matrix(pat, ph, band, "imcoh_abs").sum(axis=1))
        s_node = np.mean(np.vstack(smat), axis=0)
        for g in GRANULARITIES:
            uv = rdf[g].to_numpy().astype(str)
            om = unit_means_from_s(s_obs, iu_i, iu_j, uv, keep, DROP[g], shaft_vec=shaft_vec)
            for u, m in om.items():
                obs_um[g].setdefault(u, {})[pat] = m
            per_unit_r = {u: [] for u in om}
            for r in range(R):
                cr = {ph: surco[ph][r] for ph in PHASES}
                if any(cr[ph] is None for ph in PHASES):
                    continue
                s_r = _sym_concordance(cr["rest_pre_A"], cr["rest_pre_B"],
                                       cr["task_test"], cr["rest_post"])
                sm = unit_means_from_s(s_r, iu_i, iu_j, uv, keep, DROP[g], shaft_vec=shaft_vec)
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
            samplers = sorted(pm); K = len(samplers)
            if K == 0:
                continue
            M_obs = float(np.median([pm[p] for p in samplers]))
            arrs = [surr_um[g][u][p] for p in samplers]
            Lmin = min((len(a) for a in arrs), default=0)
            if Lmin == 0:
                continue
            Rmat = np.vstack([a[:Lmin] for a in arrs])
            M_surr = np.median(Rmat, axis=0); nval = M_surr.size
            p_ms = (1 + int(np.sum(M_surr >= M_obs))) / (nval + 1)
            p_ms_lower = (1 + int(np.sum(M_surr <= M_obs))) / (nval + 1)
            strv = strength_um[g].get(u, {})
            str_med = (float(np.median([strv[p] for p in samplers if p in strv]))
                       if any(p in strv for p in samplers) else np.nan)
            rows.append({
                "band": band, "epi": "exclude" if epi_excl else "include",
                "granularity": g, "unit": u, "K_implanted": K, "M_obs": M_obs,
                "surr_median": float(np.median(M_surr)),
                "matched_strength_p": p_ms, "matched_strength_p_lower": p_ms_lower,
                "n_surr": nval, "strength_dev_median": str_med})
            if verbose and u in ("OFC", "MTL", "cingulate", "sensorimotor", "PFC", "limbic"):
                print(f"  {band:9s} {g:10s} {u:12s}: M_obs={M_obs:+.4f} "
                      f"surr={np.median(M_surr):+.4f} MS_p={p_ms:.3f} low={p_ms_lower:.3f} (K={K})")
    return rows


def main():
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    bands = ["beta"]
    for epi_excl in (False, True):
        etag = "exclude" if epi_excl else "include"
        print(f"\n==== rho_sym localization : epi-{etag} ====")
        allrows = []
        for band in bands:
            allrows += run_band(band, epi_excl)
        df = pd.DataFrame(allrows)
        df.to_csv(OUT_ROOT / f"matched_strength_rhosym_{etag}.csv", index=False)
        # BH-FDR over a-priori SYSTEMS within beta (upper-tail = carrier)
        sysdf = df[(df.band == "beta") & (df.granularity == "system")].copy()
        if not sysdf.empty:
            sysdf["q_upper"] = bh_fdr(sysdf["matched_strength_p"].values)
            sysdf["q_lower"] = bh_fdr(sysdf["matched_strength_p_lower"].values)
            print(f"\n  -- beta systems, epi-{etag}, BH-FDR over {len(sysdf)} systems --")
            for _, r in sysdf.sort_values("matched_strength_p").iterrows():
                tag = "CARRIER" if r.q_upper < 0.05 else ("DEPLETED" if r.q_lower < 0.05 else "")
                print(f"     {r.unit:14s} M={r.M_obs:+.4f} p_up={r.matched_strength_p:.3f} "
                      f"q_up={r.q_upper:.3f} | p_lo={r.matched_strength_p_lower:.3f} "
                      f"q_lo={r.q_lower:.3f}  {tag}")
    print(f"\nOutputs -> {OUT_ROOT}/matched_strength_rhosym_*.csv")


if __name__ == "__main__":
    main()
