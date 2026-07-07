#!/usr/bin/env python3
"""audit_159 — within-system ABSOLUTE trace under the rho_sym estimator (R2.5).

Faithful re-run of audit_112 (the focal low-gamma cingulate encoding trace that
whole-brain averaging hides) with exactly ONE change: the per-pair rank concordance
that feeds every within-system rho becomes the SYMMETRIC average of both split-half
arm assignments (mirrors rho_sym = 1/2[rho_AB + rho_BA]).

audit_112's within-system statistic is a Spearman (or first-order partial Spearman)
over the pairs INCIDENT to a system S. Under rho_sym we average the two arms:

    encoding_sym(S)     = 1/2 [ rho(e_A[keep_S], p_B[keep_S]) + rho(e_B[keep_S], p_A[keep_S]) ]
    inference_pe_sym(S) = 1/2 [ pr(f[keep], p_B[keep] | e_A[keep]) + pr(f[keep], p_A[keep] | e_B[keep]) ]
    standard_sym(S)     = 1/2 [ rho(g_A[keep], p_B[keep]) + rho(g_B[keep], p_A[keep]) ]

where arm1 uses (rest_pre_A on the task/encoding side, rest_pre_B on persistence)
and arm2 SWAPS the two rest_pre halves. f = D_taskTest - D_taskLearn is arm-invariant.
Implemented by calling audit_110._targets_from_D on the phase dict and on the same
dict with rest_pre_A<->rest_pre_B swapped -> exact arm2 for every target.

Applied IDENTICALLY to observed and to every cached matched-strength surrogate
realization (canonical seed-20260511 R=200 ensemble; NO regeneration). Does NOT edit
audit_112; new script, new output.

5-point preamble
1. Claim: in low-gamma there is a genuine ABSOLUTE encoding trace on cingulate edges
   that the weak whole-brain low-gamma trace hides, and it survives rho_sym.
2. Null: the cingulate low-gamma encoding trace was an artifact of the arbitrary A/B arm.
3. Strongest alternative: it is node strength (a hub), not a real trace. Controlled by
   the same strength-preserving matched-strength surrogate audit_112 uses; the
   estimator change is orthogonal to strength.
4. Cannot: reuses R=200 cached surrogate; endpoint-incidence keep_S; BH across systems
   within band x target. This is an ABSOLUTE within-region test (no demeaning), the
   same cohort gate as the arc.
5. Falsify: if cingulate low-gamma encoding upper-tail Wilcoxon (or its BH q over
   systems) fails under rho_sym, the focal trace is estimator-dependent.

Correctness anchor A1: arm1 (rho_split) cohort-median rho per (band, target, system)
reproduces audit_112 within_system_trace_{epi}.csv (median_obs_rho).

Output: data/audit/inference_localization_rhosym/within_system_trace_rhosym_{include,exclude}.csv
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_110_inference_mark_localization import (  # type: ignore
    _targets_from_D, _partial_rho_formula, load_surrogate_cophenet5,
    _canon_cophenet, load_phase_fc, COHORT, PHASES5,
)

TARGET_BANDS = ["alpha", "beta", "low_gamma"]
TARGETS = ("standard", "encoding", "inference_pe")
MIN_PAIRS = 30
OUT_ROOT = ROOT / "data/audit/inference_localization_rhosym"
SYSTEMS = ["OFC", "cingulate", "MTL", "insula", "lateral_temporal", "PFC",
           "parietal", "sensorimotor", "occipital"]
A112 = {  # arm1 anchor reference (rho_split)
    "include": ROOT / "data/audit/inference_localization/within_system_trace_include.csv",
    "exclude": ROOT / "data/audit/inference_localization/within_system_trace_exclude.csv",
}


def _swap_halves(D: dict) -> dict:
    Ds = dict(D)
    Ds["rest_pre_A"], Ds["rest_pre_B"] = D["rest_pre_B"], D["rest_pre_A"]
    return Ds


def _trace_arm(a, b, kind, c=None):
    """audit_112's within-subset statistic (already-sliced vectors)."""
    if a.size < MIN_PAIRS:
        return np.nan
    if kind == "inference_pe":
        return _partial_rho_formula(a, b, c)
    return float(spearmanr(a, b).correlation)


def _within_sym(D: dict, keep: np.ndarray):
    """rho_sym + arm1 within-system rho per target, restricted to keep_S pairs."""
    d1 = _targets_from_D(D)["_diffs"]
    d2 = _targets_from_D(_swap_halves(D))["_diffs"]
    spec = {"standard": ("g", "p", None), "encoding": ("e", "p", None),
            "inference_pe": ("f", "p", "e")}
    sym, arm1 = {}, {}
    for t, (av, bv, cv) in spec.items():
        r1 = _trace_arm(d1[av][keep], d1[bv][keep], t,
                        d1[cv][keep] if cv else None)
        r2 = _trace_arm(d2[av][keep], d2[bv][keep], t,
                        d2[cv][keep] if cv else None)
        arm1[t] = r1
        sym[t] = np.nan if not (np.isfinite(r1) and np.isfinite(r2)) else 0.5 * (r1 + r2)
    return sym, arm1


def run_band(band, epi_excl, R, verbose=True):
    rows = []
    obs = {t: {} for t in TARGETS}          # target -> system -> {pat: rho_sym}
    obs_arm1 = {t: {} for t in TARGETS}     # arm1 (anchor)
    surr = {t: {} for t in TARGETS}         # target -> system -> {pat: array(R)}
    for pat in COHORT:
        surco, miss = load_surrogate_cophenet5(pat, band, R)
        if surco is None:
            if verbose:
                print(f"  [skip] {pat} {band}: surrogate missing ({miss}, R={R})")
            continue
        try:
            D = {ph: _canon_cophenet(load_phase_fc(pat, ph, band)) for ph in PHASES5}
        except FileNotFoundError:
            continue
        N = int((1 + np.sqrt(1 + 8 * D["task_test"].size)) / 2)
        iu = np.triu_indices(N, k=1)
        rdf = load_channel_regions(pat)
        uv = rdf["system"].to_numpy().astype(str)
        ndrop = epi_keep_mask(rdf, pat) if epi_excl else None
        si, sj = uv[iu[0]], uv[iu[1]]
        # surrogate per-realization phase dicts (None if any phase strength-violating)
        Dr_list = []
        for r in range(R):
            Dr = {ph: surco[ph][r] for ph in PHASES5}
            Dr_list.append(None if any(Dr[ph] is None for ph in PHASES5) else Dr)
        for S in SYSTEMS:
            keep = (si == S) | (sj == S)
            if ndrop is not None:
                keep = keep & ndrop[iu[0]] & ndrop[iu[1]]
            if int(keep.sum()) < MIN_PAIRS:
                continue
            osym, oarm1 = _within_sym(D, keep)
            srt = {t: [] for t in TARGETS}
            for Dr in Dr_list:
                if Dr is None:
                    continue
                ssym, _ = _within_sym(Dr, keep)
                for t in TARGETS:
                    if np.isfinite(ssym[t]):
                        srt[t].append(ssym[t])
            for t in TARGETS:
                if np.isfinite(osym[t]):
                    obs[t].setdefault(S, {})[pat] = osym[t]
                    obs_arm1[t].setdefault(S, {})[pat] = oarm1[t]
                    surr[t].setdefault(S, {})[pat] = np.asarray(srt[t])

    for t in TARGETS:
        sys_rows = []
        for S in SYSTEMS:
            pm = obs[t].get(S, {})
            samplers = sorted(pm)
            if len(samplers) < 3:
                continue
            o = np.array([pm[p] for p in samplers])
            smed = np.array([np.median(surr[t][S][p]) if len(surr[t][S][p]) else np.nan
                             for p in samplers])
            m = np.isfinite(o) & np.isfinite(smed)
            if m.sum() < 3:
                continue
            try:
                wp = wilcoxon(o[m], smed[m], alternative="greater").pvalue
            except ValueError:
                wp = np.nan
            n_sig = 0
            for p in samplers:
                a = surr[t][S][p]
                if len(a):
                    pp = (1 + int(np.sum(a >= pm[p]))) / (len(a) + 1)
                    n_sig += int(pp < 0.05)
            arm1_med = float(np.median([obs_arm1[t][S][p] for p in samplers]))
            row = {"band": band, "epi": "exclude" if epi_excl else "include",
                   "target": t, "system": S, "K_implanted": int(m.sum()),
                   "median_obs_rho": float(np.median(o[m])),
                   "median_obs_rho_arm1": arm1_med,
                   "median_surr_rho": float(np.median(smed[m])),
                   "wilcoxon_p": float(wp), "n_pat_p05": n_sig,
                   "n_pos": int((o[m] > 0).sum()), "bh_q": np.nan}
            rows.append(row)
            sys_rows.append(row)
        if sys_rows:
            qs = bh_fdr([r["wilcoxon_p"] for r in sys_rows])
            for r, q in zip(sys_rows, qs):
                r["bh_q"] = q
        if verbose:
            for r in sorted(sys_rows, key=lambda r: r["wilcoxon_p"]):
                star = "*" if r["bh_q"] < 0.05 else (" " if r["wilcoxon_p"] >= 0.05 else "·")
                print(f"  {band:9s} {t:13s} {r['system']:16s} "
                      f"obsρ_sym={r['median_obs_rho']:+.3f} (arm1={r['median_obs_rho_arm1']:+.3f}) "
                      f"surrρ={r['median_surr_rho']:+.3f} W_p={r['wilcoxon_p']:.3f} "
                      f"q={r['bh_q']:.3f}{star} ({r['n_pos']}/{r['K_implanted']} pos)")
    return rows


def _anchor(df_all):
    """A1: arm1 cohort-median rho per (band,target,system) == audit_112 median_obs_rho."""
    print("\n[anchor A1] arm1 (rho_split) vs audit_112 within_system_trace:")
    worst = 0.0
    for etag in ("include", "exclude"):
        ref = pd.read_csv(A112[etag])
        cur = df_all[df_all.epi == etag]
        for _, r in cur.iterrows():
            m = ref[(ref.band == r.band) & (ref.target == r.target)
                    & (ref.system == r.system)]
            if not m.empty:
                d = abs(float(m.iloc[0]["median_obs_rho"]) - r.median_obs_rho_arm1)
                worst = max(worst, d)
    print(f"           max |Δ| arm1 vs audit_112 = {worst:.2e}"
          f"  {'PASS' if worst < 1e-9 else 'CHECK'}")
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default=None)
    ap.add_argument("--R", type=int, default=200)
    args = ap.parse_args()
    bands = [args.band] if args.band else TARGET_BANDS
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for epi_excl in (False, True):
        etag = "exclude" if epi_excl else "include"
        print(f"\n==== within-system ABSOLUTE trace (rho_sym) : epi-{etag} (R={args.R}) ====")
        rows = []
        for band in bands:
            rows += run_band(band, epi_excl, args.R)
        pd.DataFrame(rows).to_csv(
            OUT_ROOT / f"within_system_trace_rhosym_{etag}.csv", index=False)
        all_rows += rows
    _anchor(pd.DataFrame(all_rows))
    print(f"\nOutputs -> {OUT_ROOT}/within_system_trace_rhosym_*.csv")


if __name__ == "__main__":
    main()
