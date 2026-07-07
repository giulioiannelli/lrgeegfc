#!/usr/bin/env python3
"""audit_158 — ENCODING -> OFC localization under the rho_sym estimator (R2.3).

R2.3 needs the ENCODING component of the consolidation arc localized under rho_sym:
does the learning-phase (task_learn) trace CONCENTRATE in OFC, the same hotspot the
test-phase trace uses (R1.2)? The rho_split answer (audit_110, R=1000) is OFC
MS_p=0.001, BH q=0.010 — but it is bare rho_split. This migrates it onto the
arbitrary-half-free rho_sym estimator so no rho_split number lands in R2.3.

The encoding target is the arc's `encoding` concordance: e = D_taskLearn - D_restPre
against p = D_restPost - D_restPre. This is audit_151's localizer with ONE swap: the
"task-side" cophenetic is task_learn instead of task_test. audit_151._sym_concordance
is generic in that argument, so:

    s_sym(encoding) = 1/2[ concordance(D_taskLearn-D_preA, D_post-D_preB)
                         + concordance(D_taskLearn-D_preB, D_post-D_preA) ]

applied identically to observed and to every cached matched-strength surrogate.
Reuses the canonical seed-20260511 R=1000 surrogate cophenetic (rest_pre_A/B +
rest_post from audit_155's ensemble, task_learn from audit_111) — NO regeneration.
Does NOT edit audit_110/151. band=beta, both epi modes, contact-level.

5-point preamble
1. Claim: the ENCODING (task_learn) beta trace concentrates in OFC (upper-tail
   matched-strength, BH-clearing over the 9 a-priori systems) under rho_sym, as under
   rho_split — so OFC is anchored by BOTH task phases, not the estimator's A/B half.
2. Null: the encoding->OFC concentration was an artifact of the arbitrary A/B arm.
3. Strongest alternative: (a) node strength -> controlled by the same strength-
   preserving surrogate; (b) rho_sym trivially reproduces rho_split -> controlled by
   recording arm1 and cross-checking it reproduces audit_110's encoding M_obs bit-for-
   bit. e (encoding) uses the SHORTER learn phase and carries no test/learn length
   asymmetry, so the duration confound never touches it.
4. Cannot: reuses R=1000 cached surrogate; n=10 BH family = a-priori systems within
   beta; OFC sampled in 5/10 patients (a coverage limit, not a null result).
5. Falsify: if OFC upper-tail MS_p or its BH q over systems fails under rho_sym, the
   encoding->OFC anchor is estimator-dependent and R2.3 drops the "learning sets the
   anchor" leg.

Correctness anchor (printed):
  A1  arm1 (rho_split) encoding M_obs per system == audit_110 R=1000 encoding CSV.

Output: data/audit/inference_localization_rhosym/encoding_localization_rhosym_{include,exclude}.csv
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.metrics.node_localization import (
    canonical_cophenet, epi_keep_mask, rank_concordance,
)
from lrg_eegfc.utils.surrogate import cophenetic_condensed_from_eigs, surrogate_cache_path
from lrg_eegfc.workflow.fc import load_fc_matrix

# Reuse audit_151's rho_sym localizer primitives verbatim (import, no fork).
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_151_localization_rhosym import (  # type: ignore
    COHORT, DROP, GRANULARITIES, HALVES_FC_CACHE, PARACORE,
    _sym_concordance, unit_means_from_s,
)

ENC_PHASES = ("rest_pre_A", "rest_pre_B", "task_learn", "rest_post")
R, SWAP, SEED = 1000, 20, 20260511
BAND = "beta"
OUT_ROOT = ROOT / "data/audit/inference_localization_rhosym"
# arm1 (rho_split encoding) cross-check reference:
A110_R1000 = ROOT / "data/audit/inference_localization/inference_mark_R1000_include.csv"


def obs_encoding(pat: str, band: str):
    """Observed sym + arm1 (rho_split) encoding concordance, canonical cophenetic."""
    Wl = load_fc_matrix(pat, "task_learn", band, "imcoh_abs")
    Wpost = load_fc_matrix(pat, "rest_post", band, "imcoh_abs")
    WA = np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_A_imcoh_abs.npy")
    WB = np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_B_imcoh_abs.npy")
    Dl, Dpost, DA, DB = (canonical_cophenet(W) for W in (Wl, Wpost, WA, WB))
    s_sym = _sym_concordance(DA, DB, Dl, Dpost)
    s_arm1 = rank_concordance(Dl - DA, Dpost - DB)   # == audit_110 encoding (rho_split)
    N = Wl.shape[0]
    iu = np.triu_indices(N, k=1)
    return s_sym, s_arm1, iu[0].astype(int), iu[1].astype(int)


def load_surr_encoding(pat: str, band: str):
    out = {}
    for ph in ENC_PHASES:
        path = surrogate_cache_path(pat, band, ph, R, SWAP, SEED, "imcoh_abs")
        if not path.exists():
            return None, ph
        with np.load(path) as d:
            evals, evecs = d["eigvals"], d["eigvecs"]
        cond = [None if not np.all(np.isfinite(evals[r]))
                else cophenetic_condensed_from_eigs(evals[r], evecs[r])
                for r in range(evals.shape[0])]
        out[ph] = cond
    return out, None


def run_band(band: str, epi_excl: bool, verbose: bool = True):
    rows = []
    obs_um = {g: {} for g in GRANULARITIES}
    arm1_um = {g: {} for g in GRANULARITIES}
    surr_um = {g: {} for g in GRANULARITIES}
    strength_um = {g: {} for g in GRANULARITIES}
    for pat in COHORT:
        surco, miss = load_surr_encoding(pat, band)
        if surco is None:
            if verbose:
                print(f"  [skip] {pat} {band}: surrogate cache missing ({miss}, R={R})")
            continue
        try:
            s_sym, s_arm1, iu_i, iu_j = obs_encoding(pat, band)
        except FileNotFoundError:
            if verbose:
                print(f"  [skip] {pat} {band}: half-FC / phase missing")
            continue
        rdf = load_channel_regions(pat)
        rdf["paracore"] = rdf["system"].map(
            lambda s: "paralimbic_core" if s in PARACORE
            else ("non_anatomical" if s == "non_anatomical" else "other_np"))
        keep = epi_keep_mask(rdf, pat) if epi_excl else None
        s_node = np.mean(np.vstack(
            [load_fc_matrix(pat, ph, band, "imcoh_abs").sum(axis=1)
             for ph in ("rest_pre", "task_learn", "rest_post")]), axis=0)
        for g in GRANULARITIES:
            uv = rdf[g].to_numpy().astype(str)
            om = unit_means_from_s(s_sym, iu_i, iu_j, uv, keep, DROP[g])
            om1 = unit_means_from_s(s_arm1, iu_i, iu_j, uv, keep, DROP[g])
            for u, m in om.items():
                obs_um[g].setdefault(u, {})[pat] = m
            for u, m in om1.items():
                arm1_um[g].setdefault(u, {})[pat] = m
            per = {u: [] for u in om}
            for r in range(R):
                cr = {ph: surco[ph][r] for ph in ENC_PHASES}
                if any(cr[ph] is None for ph in ENC_PHASES):
                    continue
                s_r = _sym_concordance(cr["rest_pre_A"], cr["rest_pre_B"],
                                       cr["task_learn"], cr["rest_post"])
                sm = unit_means_from_s(s_r, iu_i, iu_j, uv, keep, DROP[g])
                for u in per:
                    if u in sm:
                        per[u].append(sm[u])
            for u in om:
                surr_um[g].setdefault(u, {})[pat] = np.array(per[u])
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
            p_ms = (1 + int(np.sum(M_surr >= M_obs))) / (nval + 1)
            p_ms_lower = (1 + int(np.sum(M_surr <= M_obs))) / (nval + 1)
            strv = strength_um[g].get(u, {})
            str_med = (float(np.median([strv[p] for p in samplers if p in strv]))
                       if any(p in strv for p in samplers) else np.nan)
            arm1_med = (float(np.median([arm1_um[g][u][p] for p in samplers
                                         if p in arm1_um[g].get(u, {})]))
                        if u in arm1_um[g] else np.nan)
            rows.append({
                "band": band, "epi": "exclude" if epi_excl else "include",
                "target": "encoding", "granularity": g, "unit": u, "K_implanted": K,
                "M_obs": M_obs, "M_obs_arm1": arm1_med, "surr_median": float(np.median(M_surr)),
                "matched_strength_p": p_ms, "matched_strength_p_lower": p_ms_lower,
                "n_surr": nval, "strength_dev_median": str_med})
            if verbose and u in ("OFC", "MTL", "cingulate", "sensorimotor", "PFC"):
                print(f"  {band:9s} {g:10s} {u:12s}: M_obs={M_obs:+.4e} "
                      f"surr={np.median(M_surr):+.4e} MS_p={p_ms:.3f} (K={K})")
    return rows


def _arm1_crosscheck(all_rows: list[dict]) -> None:
    """arm1 encoding M_obs per beta/system == audit_110 R=1000 encoding (rho_split)."""
    if not A110_R1000.exists():
        print("[A1] skipped — audit_110 R=1000 CSV not found.")
        return
    ref = pd.read_csv(A110_R1000)
    ref = ref[(ref.band == "beta") & (ref.granularity == "system")
              & (ref.target == "encoding")].set_index("unit")["M_obs"]
    diffs = []
    for r in all_rows:
        if (r["epi"] == "include" and r["granularity"] == "system"
                and r["unit"] in ref.index and np.isfinite(r["M_obs_arm1"])):
            diffs.append(abs(r["M_obs_arm1"] - float(ref.loc[r["unit"]])))
    if diffs:
        md = max(diffs)
        st = "PASS" if md < 1e-6 else ("CLOSE" if md < 1e-2 else "FAIL")
        print(f"[A1] arm1 encoding M_obs == audit_110 R=1000 (rho_split): "
              f"n={len(diffs)} systems max|Δ|={md:.2e} → {st}")


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    include_rows = None
    for epi_excl in (False, True):
        etag = "exclude" if epi_excl else "include"
        print(f"\n==== rho_sym ENCODING->OFC localization : epi-{etag}  (R={R}, beta) ====")
        rows = run_band(BAND, epi_excl)
        df = pd.DataFrame(rows)
        df.to_csv(OUT_ROOT / f"encoding_localization_rhosym_{etag}.csv", index=False)
        sysdf = df[(df.band == "beta") & (df.granularity == "system")].copy()
        if not sysdf.empty:
            sysdf["q_upper"] = bh_fdr(sysdf["matched_strength_p"].values)
            sysdf["q_lower"] = bh_fdr(sysdf["matched_strength_p_lower"].values)
            print(f"\n  -- beta ENCODING systems, epi-{etag}, BH-FDR over "
                  f"{len(sysdf)} systems (upper = carrier) --")
            for _, r in sysdf.sort_values("matched_strength_p").iterrows():
                tag = ("CARRIER" if r.q_upper < 0.05
                       else ("DEPLETED" if r.q_lower < 0.05 else ""))
                print(f"     {r.unit:16s} M={r.M_obs:+.3e} p_up={r.matched_strength_p:.3f} "
                      f"q_up={r.q_upper:.3f} | q_lo={r.q_lower:.3f}  {tag}")
        if not epi_excl:
            include_rows = rows
    print("\n=== CORRECTNESS ANCHOR ===")
    if include_rows is not None:
        _arm1_crosscheck(include_rows)
    print(f"\nOutputs -> {OUT_ROOT}/encoding_localization_rhosym_*.csv")


if __name__ == "__main__":
    main()
