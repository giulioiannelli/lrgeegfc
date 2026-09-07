#!/usr/bin/env python3
r"""Lane E / E3 -- gate the node-loading dissociation.

Same gate as the trace verdict: knob-integrated margin against matched-strength,
one sign-flip cluster-mass test per (band, statistic) over the whole 16-point
scale axis, BH within a statistic across the bands run, leave-one-patient-out of
the verdict, and the per-scale grid reported in full.

Reads ``loading/cells/*.npz`` from ``lane_e_05_loading_dissociation.py``.
Outputs: data/paper_final/lane_e_encinf/verdict/loading_*.csv
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.utils.metrics.cohort_gate import axis_cluster_gate, cohort_margin_gate
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr, boot_ci_mean
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
BASE = Path(os.environ.get("LANE_E_OUT",
                           ROOT / "data" / "paper_final" / "lane_e_encinf"))
CELLS, OUT = BASE / "loading" / "cells", BASE / "verdict"
N_PERM = int(os.environ.get("LANE_E_NPERM", "10000"))
SEED = 20260903


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", default="delta,theta,alpha,beta,low_gamma")
    a = ap.parse_args()
    bands = a.bands.split(",")
    pats = list(PATIENTS_4PHASE)
    OUT.mkdir(parents=True, exist_ok=True)

    s = stats = None
    store = {}
    for band in bands:
        M, OK, SK, keep = [], [], [], []
        for p in pats:
            f = CELLS / f"{p}__{band}.npz"
            if not f.exists():
                continue
            d = np.load(f)
            s, stats = d["s"], [str(x) for x in d["stats"]]
            obs, surr = d["obs"], d["surr"]                # (nF,nS,nK)/(nF,R,nS,nK)
            M.append(np.nanmedian(obs - np.nanmedian(surr, axis=1), axis=0))
            OK.append(np.nanmedian(obs, axis=0))
            SK.append(np.nanmedian(surr, axis=0))
            keep.append(p)
        if keep:
            store[band] = dict(patients=keep, margin=np.array(M),
                               obs=np.array(OK), surr=np.array(SK))
    if s is None:
        raise SystemExit(f"no loading cells under {CELLS}")
    nS = s.size

    rows, grid, loo = [], [], []
    for band, st in store.items():
        for k, nm in enumerate(stats):
            Mk = st["margin"][:, :, k]
            r = axis_cluster_gate(Mk, n_perm=N_PERM,
                                  rng=np.random.default_rng(SEED))
            mean, lo, hi = boot_ci_mean(np.nanmean(Mk, axis=1), B=10000,
                                        rng=np.random.default_rng(SEED))
            rows.append(dict(band=band, stat=nm, n_pat=Mk.shape[0],
                             p_cluster=r["p"], margin_med=float(np.nanmedian(Mk)),
                             margin_mean=mean, ci_lo=lo, ci_hi=hi,
                             obs_med=float(np.nanmedian(st["obs"][:, :, k])),
                             surr_med=float(np.nanmedian(st["surr"][:, :, :, k])),
                             s_lo=s[r["cluster"][0]] if r["cluster"] else np.nan,
                             s_hi=s[r["cluster"][1]] if r["cluster"] else np.nan))
            for i, p in enumerate(st["patients"]):
                keep = np.ones(Mk.shape[0], bool)
                keep[i] = False
                loo.append(dict(band=band, stat=nm, dropped=p,
                                p_cluster=axis_cluster_gate(
                                    Mk[keep], n_perm=N_PERM,
                                    rng=np.random.default_rng(SEED))["p"]))
            for j in range(nS):
                g = cohort_margin_gate(st["obs"][:, j, k], st["surr"][:, :, j, k],
                                       labels=st["patients"], full=False,
                                       expect_n=len(pats))
                grid.append(dict(band=band, stat=nm, s=float(s[j]),
                                 p=g["p"], margin_med=g["margin_med"],
                                 obs_med=g["obs_med"], surr_med=g["surr_med"]))
    cl = pd.DataFrame(rows)
    for nm in stats:
        m = cl.stat == nm
        cl.loc[m, "q_within_stat"] = bh_fdr(cl.loc[m, "p_cluster"].values)
    cl.to_csv(OUT / "loading_cluster_gate.csv", index=False)
    pd.DataFrame(loo).to_csv(OUT / "loading_cluster_loo.csv", index=False)
    gd = pd.DataFrame(grid)
    for nm in stats:
        m = gd.stat == nm
        gd.loc[m, "q_within_stat"] = bh_fdr(gd.loc[m, "p"].values)
    gd.to_csv(OUT / "loading_per_scale.csv", index=False)

    pd.set_option("display.width", 220)
    print("\n=== E3 node-loading dissociation, margin vs matched-strength ===")
    print(cl[["band", "stat", "p_cluster", "q_within_stat", "margin_med",
              "margin_mean", "ci_lo", "ci_hi", "obs_med", "surr_med"]]
          .round(4).to_string(index=False))
    lo = pd.DataFrame(loo)
    print("\n=== LOO of the verdict (drops pushing p >= 0.05) ===")
    agg = lo.groupby(["band", "stat"]).p_cluster.agg(
        loo_min="min", loo_max="max",
        n_ge_05=lambda x: int((x >= 0.05).sum())).reset_index()
    print(agg.round(4).to_string(index=False))
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
