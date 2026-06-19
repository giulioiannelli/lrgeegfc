#!/usr/bin/env python3
"""audit_112 — within-system ABSOLUTE trace: is there a focal trace where the
whole-brain trace is weak?

Tests the hypothesis (user, 2026-06-18): maybe α / low-γ traces are "not clear
overall but specific areas present strong traces" — i.e. a focal trace diluted to
non-significance by whole-brain averaging. The audit_110 localization CANNOT
answer this: it is per-patient DEMEANED, so a net-zero band still shows some
systems "above average" (a concentration), which is not the same as an absolute
trace. This test removes the demeaning: it restricts the cohort split-baseline
trace to the pairs INCIDENT to a system and asks whether THAT clears
matched-strength — an absolute within-region trace, no concentration confound.

Statistic (per patient, per band, per system S, per target):
  keep_S = pairs with at least one endpoint in S (endpoint incidence, matching
           the localization's aggregation).
  standard      : ρ = Spearman(g[keep], p[keep])
  encoding      : ρ = Spearman(e[keep], p[keep])
  inference_pe  : ρ = partial Spearman(f[keep], p[keep] | e[keep])
  (g=D_TT−D_preA, e=D_TL−D_preA, f=D_TT−D_TL, p=D_RP−D_preB; canonical cophenetic)
Null: same ρ on the matched-strength surrogate cophenetics restricted to keep_S
(the canonical ensemble; R=200 all bands, R=1000 β task_learn available).
Verdict: one-sided paired Wilcoxon (obs > surrogate median) across implanted
patients + per-patient n(p<0.05) + BH across systems within band×target. This is
the SAME cohort gate as the arc / audit_83 trace tests — NOT a demeaned
concentration. A system clears here ⇒ a genuine absolute trace lives on its edges,
even if the whole-brain trace does not.

Output: data/audit/inference_localization/within_system_trace_{include,exclude}.csv
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
    obs_targets, load_surrogate_cophenet5, _targets_from_D,
    _partial_rho_formula, COHORT, PHASES5,
)

TARGET_BANDS = ["alpha", "beta", "low_gamma"]
TARGETS = ("standard", "encoding", "inference_pe")
MIN_PAIRS = 30
OUT_ROOT = ROOT / "data/audit/inference_localization"
# the a-priori systems we report (drop tiny / non-anatomical)
SYSTEMS = ["OFC", "cingulate", "MTL", "insula", "lateral_temporal", "PFC",
           "parietal", "sensorimotor", "occipital"]


def _trace(a, b, kind, c=None):
    """ρ of (a,b) restricted to the already-sliced vectors; partial if kind=pe."""
    if a.size < MIN_PAIRS:
        return np.nan
    if kind == "inference_pe":
        return _partial_rho_formula(a, b, c)
    return float(spearmanr(a, b).correlation)


def _targets_within(diffs, keep):
    g, e, f, p = (diffs["g"][keep], diffs["e"][keep],
                  diffs["f"][keep], diffs["p"][keep])
    return {
        "standard":     _trace(g, p, "standard"),
        "encoding":     _trace(e, p, "encoding"),
        "inference_pe": _trace(f, p, "inference_pe", e),
    }


def run_band(band, epi_excl, R, verbose=True):
    rows = []
    # obs[target][system][pat], surr[target][system][pat] = array(R)
    obs = {t: {} for t in TARGETS}
    surr = {t: {} for t in TARGETS}
    for pat in COHORT:
        surco, miss = load_surrogate_cophenet5(pat, band, R)
        if surco is None:
            if verbose:
                print(f"  [skip] {pat} {band}: surrogate missing ({miss}, R={R})")
            continue
        try:
            t_obs, iu_i, iu_j = obs_targets(pat, band)
        except FileNotFoundError:
            continue
        diffs = t_obs["_diffs"]
        rdf = load_channel_regions(pat)
        uv = rdf["system"].to_numpy().astype(str)
        ndrop = epi_keep_mask(rdf, pat) if epi_excl else None
        si, sj = uv[iu_i], uv[iu_j]
        # precompute surrogate diffs per realization (None if strength-violating)
        srd = []
        for r in range(R):
            Dr = {ph: surco[ph][r] for ph in PHASES5}
            srd.append(None if any(Dr[ph] is None for ph in PHASES5)
                       else _targets_from_D(Dr)["_diffs"])
        for S in SYSTEMS:
            keep = (si == S) | (sj == S)
            if ndrop is not None:  # epi-exclude: drop pairs touching an epi node
                keep = keep & ndrop[iu_i] & ndrop[iu_j]
            if int(keep.sum()) < MIN_PAIRS:
                continue
            ot = _targets_within(diffs, keep)
            sr_t = {t: [] for t in TARGETS}
            for r in range(R):
                if srd[r] is None:
                    continue
                st = _targets_within(srd[r], keep)
                for t in TARGETS:
                    if np.isfinite(st[t]):
                        sr_t[t].append(st[t])
            for t in TARGETS:
                if np.isfinite(ot[t]):
                    obs[t].setdefault(S, {})[pat] = ot[t]
                    surr[t].setdefault(S, {})[pat] = np.asarray(sr_t[t])

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
            # per-patient matched-strength p
            n_sig = 0
            for p in samplers:
                a = surr[t][S][p]
                if len(a):
                    pp = (1 + int(np.sum(a >= pm[p]))) / (len(a) + 1)
                    n_sig += int(pp < 0.05)
            row = {"band": band, "epi": "exclude" if epi_excl else "include",
                   "target": t, "system": S, "K_implanted": int(m.sum()),
                   "median_obs_rho": float(np.median(o[m])),
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
                      f"obsρ={r['median_obs_rho']:+.3f} surrρ={r['median_surr_rho']:+.3f} "
                      f"W_p={r['wilcoxon_p']:.3f} q={r['bh_q']:.3f}{star} "
                      f"({r['n_pos']}/{r['K_implanted']} pos, {r['n_pat_p05']} sig)")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default=None)
    ap.add_argument("--epi-mode", choices=["include", "exclude"], default=None)
    ap.add_argument("--R", type=int, default=200)
    args = ap.parse_args()
    bands = [args.band] if args.band else TARGET_BANDS
    epis = ([args.epi_mode == "exclude"] if args.epi_mode else [False, True])
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    for epi_excl in epis:
        etag = "exclude" if epi_excl else "include"
        print(f"\n==== within-system ABSOLUTE trace : epi-{etag} (R={args.R}) ====")
        print("  (Wilcoxon obs>surr_med across implanted patients; * = BH q<0.05 "
              "within band×target; absolute trace, NOT demeaned concentration)")
        allrows = []
        for band in bands:
            allrows += run_band(band, epi_excl, args.R)
        out = OUT_ROOT / f"within_system_trace_{etag}.csv"
        pd.DataFrame(allrows).to_csv(out, index=False)
        print(f"  -> {out} ({len(allrows)} rows)")


if __name__ == "__main__":
    main()
