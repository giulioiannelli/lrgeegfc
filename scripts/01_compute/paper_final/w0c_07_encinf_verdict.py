#!/usr/bin/env python3
r"""W0-C: does learning a structure and applying it persist at different scales?

Reads the five-phase cells from ``w0c_06_encinf_scale_grid`` and answers, in
order:

1. **Is the conditional functional reportable at all?** ``T_infspec_pe`` is a
   partial correlation. Two calibration checks run before any p-value is
   quoted: the held-out-realization test (is the gate calibrated against the
   null actually in use?) and the role-permutation placebo (does the estimator
   return systematically positive values when the five observed graphs are
   assigned to the five roles arbitrarily?). If either fails, the gate withholds
   the p-value and this script says so instead of reporting it.
2. **Do the four functionals persist?** The locked margin gate, per scale, whole
   grid, plus the axis-cluster version so the multiplicity family is per (band,
   functional) rather than per (band, functional, scale).
3. **The scale dissociation.** ``T_learn(s)`` and ``T_infspec_pe(s)`` are two
   curves over the same axis for the same patient and band. The claim is not
   that one is larger -- it is that they live at different scales. Tested on the
   per-patient profile CENTROID in log s, paired within patient, so amplitude
   differences cannot produce it. Reported with its effect size and with what
   the design could have detected.

Outputs (data/paper_final/w0c_gate_tau/):
  encinf_gate.csv, encinf_placebo.csv, encinf_profiles.csv, encinf_cluster.csv
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.utils.metrics.cohort_gate import (
    DESCRIPTIVE_ALPHA, axis_cluster_gate, effective_tests, gate_grid,
    patient_margin)
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr, boot_ci_mean
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
BASE = Path(os.environ.get("W0C_BASE", ROOT / "data" / "paper_final" / "w0c_gate_tau"))
GRID = Path(os.environ.get("W0C_ENCINF_OUT", BASE / "encinf"))


def load():
    cfg = json.loads((GRID / "config.json").read_text())
    pats, bands, funcs = cfg["cohort"], cfg["bands"], cfg["funcs"]
    s = np.array(cfg["s_grid"])
    nF, nB, K, nS = len(funcs), len(bands), len(pats), s.size
    obs = np.full((nF, nB, K, nS), np.nan)
    surr = np.full((nF, nB, K, nS, cfg["R"]), np.nan)
    perm = np.full((nF, nB, K, nS, cfg["n_perm"]), np.nan)
    for ib, b in enumerate(bands):
        for ip, p in enumerate(pats):
            f = GRID / "cells" / f"{p}__{b}.npz"
            if not f.exists():
                continue
            d = np.load(f, allow_pickle=False)
            for iF in range(nF):
                obs[iF, ib, ip] = d["obs"][:, iF]
                surr[iF, ib, ip] = d["surr"][:, :, iF].T
                pm = d["perm"][:, :, iF].T                       # (nS, n_perm)
                perm[iF, ib, ip, :, :pm.shape[1]] = pm
    return dict(cfg=cfg, pats=pats, bands=bands, funcs=funcs, s=s,
                obs=obs, surr=surr, perm=perm)


def _paired(x, y, alt):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3 or np.allclose(x[m], y[m]):
        return np.nan, int(m.sum())
    return float(wilcoxon(x[m], y[m], alternative=alt).pvalue), int(m.sum())


# --------------------------------------------------------------------------- #
def placebo(D):
    """Role-permutation placebo: what each functional returns from arbitrary roles."""
    print("\n" + "=" * 78, flush=True)
    print("1. IS THE CONDITIONAL FUNCTIONAL REPORTABLE?  role-permutation placebo",
          flush=True)
    print("=" * 78, flush=True)
    print("\n  The five OBSERVED graphs, with the five role labels permuted. A", flush=True)
    print("  functional that is systematically positive here is manufacturing sign", flush=True)
    print("  from the estimator, not from consolidation.", flush=True)
    rows = []
    for iF, fn in enumerate(D["funcs"]):
        for ib, b in enumerate(D["bands"]):
            for js, sv in enumerate(D["s"]):
                v = D["perm"][iF, ib, :, js, :]                  # (K, n_perm)
                med = np.nanmedian(v)
                # one value per patient = that patient's placebo median, then a
                # cohort signed-rank against zero: is the placebo centred?
                pp = np.nanmedian(v, axis=1)
                pp = pp[np.isfinite(pp)]
                p_pos = (float(wilcoxon(pp, alternative="greater").pvalue)
                         if pp.size >= 3 and not np.allclose(pp, 0) else np.nan)
                rows.append(dict(functional=fn, band=b, s=float(sv),
                                 placebo_med=float(med), placebo_p_positive=p_pos,
                                 obs_med=float(np.nanmedian(D["obs"][iF, ib, :, js]))))
    pl = pd.DataFrame(rows)
    pl.to_csv(BASE / "encinf_placebo.csv", index=False)
    print("\n  functional      placebo median (cohort, over bands and scales)   "
          "scales where the", flush=True)
    print("                  [min, median, max]                               "
          "placebo is itself", flush=True)
    print("                                                                    "
          "significantly > 0", flush=True)
    verdict = {}
    for fn in D["funcs"]:
        x = pl[pl.functional == fn]
        n_bad = int((x.placebo_p_positive < DESCRIPTIVE_ALPHA).sum())
        verdict[fn] = n_bad
        print(f"  {fn:15s} [{x.placebo_med.min():+.3f}, {x.placebo_med.median():+.3f}, "
              f"{x.placebo_med.max():+.3f}]                      "
              f"{n_bad:3d}/{len(x)}", flush=True)
    print("\n  (a well-behaved functional sits at ~0 under permuted roles; a biased", flush=True)
    print("   one is offset. These p-values are uncorrected and descriptive -- the", flush=True)
    print("   question is whether the offset exists, not whether one cell clears.)",
          flush=True)
    return pl, verdict


def gates(D):
    print("\n" + "=" * 78, flush=True)
    print("2. DO THE FOUR FUNCTIONALS PERSIST?  locked margin gate", flush=True)
    print("=" * 78, flush=True)
    frames, clus_rows = [], []
    for iF, fn in enumerate(D["funcs"]):
        cells = [({"functional": fn, "band": b, "s": float(sv)},
                  D["obs"][iF, ib, :, js], D["surr"][iF, ib, :, js, :])
                 for ib, b in enumerate(D["bands"])
                 for js, sv in enumerate(D["s"])]
        g = gate_grid(cells, labels=D["pats"], rng=np.random.default_rng(0),
                      expect_n=len(PATIENTS_4PHASE),
                      calibrate=True, calibration_draws=200,
                      refuse_uncalibrated=True)
        frames.append(g)
        n_unc = int((~g["calibrated"]).sum())
        print(f"\n  {fn}   (whole-grid BH over {len(g)} cells; "
              f"{n_unc} withheld as uncalibrated)", flush=True)
        print("    band        cleared   q_min     s@q_min   margin_med   "
              "FPR@.05 (median)", flush=True)
        for ib, b in enumerate(D["bands"]):
            x = g[g.band == b]
            if x.q.notna().sum() == 0:
                print(f"    {b:11s}  --  (no reportable cell)", flush=True)
                continue
            r = x.loc[x.q.idxmin()]
            print(f"    {b:11s} {int((x.q < DESCRIPTIVE_ALPHA).sum()):2d}/{len(x):2d}    "
                  f"{'*' if r.q < DESCRIPTIVE_ALPHA else ' '}{float(r.q):6.4f}  "
                  f"{float(r.s):8.2f}  {float(r.margin_med):+.4f}      "
                  f"{x.fpr_05.median():.3f}", flush=True)
            # axis-cluster companion: one test for the whole scale axis
            M = patient_margin_profile(D, iF, ib)
            ac = axis_cluster_gate(M, n_perm=5000, rng=np.random.default_rng(7))
            et = effective_tests(M)
            clus_rows.append(dict(functional=fn, band=b, cluster_p=ac["p"],
                                  mass=ac["mass"],
                                  s_lo=float(D["s"][ac["cluster"][0]]) if ac["cluster"] else np.nan,
                                  s_hi=float(D["s"][ac["cluster"][1]]) if ac["cluster"] else np.nan,
                                  n_eff_pr=et["n_eff_pr"], n_eff_cn=et["n_eff_cn"],
                                  mean_offdiag=et["mean_offdiag"]))
    g_all = pd.concat(frames, ignore_index=True)
    g_all.to_csv(BASE / "encinf_gate.csv", index=False)

    cl = pd.DataFrame(clus_rows)
    if not cl.empty:
        cl["cluster_q"] = np.nan
        for fn in D["funcs"]:
            m = cl.functional == fn
            if m.any():
                cl.loc[m, "cluster_q"] = bh_fdr(cl.loc[m, "cluster_p"].to_numpy())
    cl.to_csv(BASE / "encinf_cluster.csv", index=False)
    print("\n  axis-cluster companion: ONE sign-flip cluster-mass test per (band,", flush=True)
    print("  functional), so the family is 6 bands rather than 6 x 28 cells.", flush=True)
    print("    functional      band        cluster_p  cluster_q   scale span      "
          "n_eff of the 28-point axis", flush=True)
    for _, r in cl.iterrows():
        span = (f"{r.s_lo:6.2f}-{r.s_hi:6.2f}" if np.isfinite(r.s_lo) else "     --     ")
        print(f"    {r.functional:15s} {r.band:11s} {r.cluster_p:8.4f}  "
              f"{r.cluster_q:8.4f}   {span}   "
              f"PR {r.n_eff_pr:5.2f} / CN {r.n_eff_cn:5.2f} "
              f"(mean r {r.mean_offdiag:+.2f})", flush=True)
    return g_all, cl


def patient_margin_profile(D, iF, ib):
    """``(K, nS)`` per-patient margins for one (functional, band)."""
    return patient_margin(
        D["obs"][iF, ib].ravel(),
        D["surr"][iF, ib].reshape(-1, D["surr"].shape[-1])
    ).reshape(D["obs"].shape[2], D["obs"].shape[3])


# --------------------------------------------------------------------------- #
def dissociation(D):
    print("\n" + "=" * 78, flush=True)
    print("3. SCALE DISSOCIATION -- do encoding and inference persist at different s?",
          flush=True)
    print("=" * 78, flush=True)
    ls = np.log(D["s"])
    iL = D["funcs"].index("T_learn")
    iI = D["funcs"].index("T_infspec_pe")
    iI2 = D["funcs"].index("T_infspec")
    rows = []
    for name, iF in (("T_learn", iL), ("T_infspec", iI2), ("T_infspec_pe", iI),
                     ("T_test", D["funcs"].index("T_test"))):
        for ib, b in enumerate(D["bands"]):
            M = patient_margin_profile(D, iF, ib)
            for ip, p in enumerate(D["pats"]):
                m = M[ip]
                w = np.clip(m, 0.0, None)
                cen = float(np.sum(w * ls) / np.sum(w)) if np.nansum(w) > 0 else np.nan
                has = np.isfinite(m).any()
                rows.append(dict(functional=name, band=b, patient=p,
                                 centroid_log_s=cen,
                                 s_centroid=float(np.exp(cen)) if np.isfinite(cen) else np.nan,
                                 s_at_max=float(D["s"][int(np.nanargmax(m))]) if has else np.nan,
                                 margin_max=float(np.nanmax(m)) if has else np.nan,
                                 margin_mean=float(np.nanmean(m)) if has else np.nan))
    prof = pd.DataFrame(rows)
    prof.to_csv(BASE / "encinf_profiles.csv", index=False)

    print("\n  per-patient profile centroid in s (cohort median [IQR]):", flush=True)
    print("    band        T_learn              T_infspec_pe          paired "
          "p(learn<inf)  p(learn>inf)   n", flush=True)
    res = []
    for b in D["bands"]:
        a = prof[(prof.functional == "T_learn") & (prof.band == b)].set_index("patient")
        c = prof[(prof.functional == "T_infspec_pe") & (prof.band == b)].set_index("patient")
        common = [p for p in D["pats"] if p in a.index and p in c.index]
        x = a.loc[common, "centroid_log_s"].to_numpy(float)
        y = c.loc[common, "centroid_log_s"].to_numpy(float)
        p_lt, n = _paired(x, y, "less")
        p_gt, _ = _paired(x, y, "greater")
        res.append(dict(band=b, p_less=p_lt, p_greater=p_gt, n=n,
                        learn_s=float(np.nanmedian(np.exp(x))),
                        inf_s=float(np.nanmedian(np.exp(y)))))
        print(f"    {b:11s} {np.nanmedian(np.exp(x)):6.2f} "
              f"[{np.nanpercentile(np.exp(x),25):5.2f},{np.nanpercentile(np.exp(x),75):6.2f}]"
              f"    {np.nanmedian(np.exp(y)):6.2f} "
              f"[{np.nanpercentile(np.exp(y),25):5.2f},{np.nanpercentile(np.exp(y),75):6.2f}]"
              f"      {p_lt:8.4f}      {p_gt:8.4f}   {n}", flush=True)
    rdf = pd.DataFrame(res)
    both = np.concatenate([rdf.p_less.to_numpy(), rdf.p_greater.to_numpy()])
    q = bh_fdr(both[np.isfinite(both)])
    print(f"\n    BH over the {int(np.isfinite(both).sum())} directional tests: "
          f"min q = {min(q):.4f}", flush=True)

    # effect size and detectability of the dissociation
    print("\n  effect size of the centroid difference (log s units, paired):", flush=True)
    for b in D["bands"]:
        a = prof[(prof.functional == "T_learn") & (prof.band == b)].set_index("patient")
        c = prof[(prof.functional == "T_infspec_pe") & (prof.band == b)].set_index("patient")
        common = [p for p in D["pats"] if p in a.index and p in c.index]
        d = (c.loc[common, "centroid_log_s"].to_numpy(float)
             - a.loc[common, "centroid_log_s"].to_numpy(float))
        d = d[np.isfinite(d)]
        if d.size < 3:
            continue
        mean, lo, hi = boot_ci_mean(d, rng=np.random.default_rng(0))
        # smallest paired shift a signed-rank at this n and this spread would catch
        rng = np.random.default_rng(5)
        mde = np.nan
        for eff in np.linspace(0.05, 2.0, 40):
            hits = sum(
                float(wilcoxon(rng.normal(eff * np.std(d), np.std(d), d.size),
                               alternative="greater").pvalue) < 0.05
                for _ in range(200)) / 200
            if hits >= 0.8:
                mde = eff * np.std(d)
                break
        mde_txt = f"{mde:.3f}" if np.isfinite(mde) else "> 2 SD (not reached)"
        print(f"    {b:11s} mean shift {mean:+.3f} [{lo:+.3f}, {hi:+.3f}] "
              f"(x{np.exp(mean):.2f} in s)   smallest shift detectable at 80% "
              f"power: {mde_txt}", flush=True)
    return prof


def main():
    D = load()
    print(f"[w0c-encinf] {len(D['funcs'])} functionals x {len(D['bands'])} bands x "
          f"{D['s'].size} scales x {len(D['pats'])} patients | "
          f"{D['cfg']['backbone']}@{D['cfg']['frac']} | R={D['cfg']['R']} "
          f"perm={D['cfg']['n_perm']}", flush=True)
    placebo(D)
    gates(D)
    dissociation(D)
    print(f"\n[w0c-encinf] -> {BASE}", flush=True)


if __name__ == "__main__":
    main()
