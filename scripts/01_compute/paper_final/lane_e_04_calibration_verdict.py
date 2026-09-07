#!/usr/bin/env python3
r"""Lane E / E2 -- is the gate calibrated for the encoding/inference functionals?

Reads the no-signal sham arcs of ``lane_e_02_sham_gate_calibration.py`` and asks
the two questions that decide whether any encoding/inference p-value in this
project is reportable.

**Q1. Does the MARGIN sit at zero on no-signal input?** W0-B measured the sham
RAW value and found all four functionals significantly positive at 16/16 scales.
The gate does not test the raw value; it tests ``obs - (that arc's own
matched-strength surrogate median)``. Both are computed here on the same arcs,
side by side, so the reader can see whether the construction offset survives the
margin subtraction. If it does, the gate is broken for these functionals and
every cell gated with it -- including the conditional one -- must be withdrawn.

**Q2. What is the measured false-positive rate?** Each shuffled sham realization
is an independent no-signal draw, so assembling one realization per patient
gives a no-signal cohort on which the entire gate can be run. Doing that many
times measures the gate's FPR directly, for the exact statistic, the exact null
and the exact cohort size in use.

The learn/test ROLE-SWAP test is calibrated on the same arcs. The sham's
``task_learn`` and ``task_test`` pseudo-phases differ in duration and position
but carry no cognitive difference, so a swap test that fires on them is
measuring block asymmetry rather than cognition -- which is exactly the
confound the real-arc swap result has to survive.

Outputs: data/paper_final/lane_e_encinf/verdict/sham_*.csv
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.utils.metrics.cohort_gate import axis_cluster_gate
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
BASE = Path(os.environ.get("LANE_E_OUT",
                           ROOT / "data" / "paper_final" / "lane_e_encinf"))
SHAM, OUT = BASE / "sham", BASE / "verdict"
FUNCS = ("T_test", "T_learn", "T_infspec", "T_infspec_pe")
N_ASSEMBLY = int(os.environ.get("LANE_E_NASSEMBLY", "100"))
N_PERM = int(os.environ.get("LANE_E_NPERM_SHAM", "1000"))
SEED = 20260831


def _p_greater(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    x = x[x != 0]
    if x.size < 3:
        return np.nan
    return float(wilcoxon(x, alternative="greater").pvalue)


def load_sham(pats, bands, sources, tag):
    """``{(band, source): dict(patients, obs, swap, margin)}`` for one construction.

    ``obs`` / ``swap`` / ``margin`` are ``(K, n_real, nS, nK)``, knob-integrated
    over the plateau fractions exactly as the real arc is.
    """
    out, s = {}, None
    for band in bands:
        for src in sources:
            O, W, M, keep = [], [], [], []
            for p in pats:
                f = SHAM / "cells" / f"{p}__{band}__{src}.npz"
                if not f.exists():
                    continue
                d = np.load(f)
                if f"obs_{tag}" not in d:
                    continue
                s = d["s"]
                o = np.nanmedian(d[f"obs_{tag}"], axis=1)          # (n_real,nS,nK)
                w = np.nanmedian(d[f"swap_{tag}"], axis=1)
                su = d[f"surr_{tag}"]                              # (n_real,nF,R,nS,nK)
                mg = np.nanmedian(d[f"obs_{tag}"] - np.nanmedian(su, axis=2), axis=1)
                O.append(o); W.append(w); M.append(mg); keep.append(p)
            if keep:
                n = min(x.shape[0] for x in O)
                out[(band, src)] = dict(
                    patients=keep,
                    obs=np.array([x[:n] for x in O]),
                    swap=np.array([x[:n] for x in W]),
                    margin=np.array([x[:n] for x in M]))
    return s, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", default="delta,theta,alpha,beta")
    ap.add_argument("--sources", default="rest_pre,rest_post")
    a = ap.parse_args()
    bands, sources = a.bands.split(","), a.sources.split(",")
    pats = list(PATIENTS_4PHASE)
    OUT.mkdir(parents=True, exist_ok=True)

    rows, fpr_rows, swap_rows = [], [], []
    for tag in ("ordered", "shuffled", "eqdur"):
        s, store = load_sham(pats, bands, sources, tag)
        if s is None:
            print(f"[lane-e/E2] no {tag} sham cells")
            continue
        nS = s.size
        for (band, src), st in store.items():
            obs, mar, swp = st["obs"], st["margin"], st["swap"]     # (K,n_real,nS,nK)
            K, n_real = obs.shape[0], obs.shape[1]
            rng = np.random.default_rng(SEED)

            # --- Q1: raw-vs-zero and margin-vs-zero, per scale, realization 0 --
            for k, fn in enumerate(FUNCS):
                for j in range(nS):
                    rows.append(dict(
                        construction=f"sham_{tag}", source=src, band=band, func=fn,
                        s=float(s[j]), n_pat=K,
                        raw_med=float(np.nanmedian(obs[:, :, j, k])),
                        raw_p_gt0=_p_greater(np.nanmedian(obs[:, :, j, k], axis=1)),
                        margin_med=float(np.nanmedian(mar[:, :, j, k])),
                        margin_p_gt0=_p_greater(np.nanmedian(mar[:, :, j, k], axis=1)),
                        surr_med=float(np.nanmedian(obs[:, :, j, k] - mar[:, :, j, k])),
                    ))

            # --- Q2: FPR of the full gate over no-signal cohort assemblies ---
            for k, fn in enumerate(FUNCS):
                pc, pw, pm = [], [], []
                n_draw = N_ASSEMBLY if n_real > 1 else 1
                for _ in range(n_draw):
                    pick = (rng.integers(0, n_real, size=K) if n_real > 1
                            else np.zeros(K, int))
                    M = mar[np.arange(K), pick, :, k]
                    D = (obs[np.arange(K), pick, :, k]
                         - swp[np.arange(K), pick, :, k])
                    pm.append(np.nanmin([_p_greater(M[:, j]) for j in range(nS)]))
                    pc.append(axis_cluster_gate(
                        M, n_perm=N_PERM, rng=np.random.default_rng(SEED))["p"])
                    pw.append(axis_cluster_gate(
                        D, n_perm=N_PERM, rng=np.random.default_rng(SEED))["p"])
                fpr_rows.append(dict(
                    construction=f"sham_{tag}", source=src, band=band, func=fn,
                    n_assembly=n_draw, n_real=n_real,
                    cluster_p_med=float(np.nanmedian(pc)),
                    fpr_cluster_05=float(np.nanmean(np.array(pc) < 0.05)),
                    fpr_cluster_10=float(np.nanmean(np.array(pc) < 0.10)),
                    swap_p_med=float(np.nanmedian(pw)),
                    fpr_swap_05=float(np.nanmean(np.array(pw) < 0.05)),
                    fpr_swap_10=float(np.nanmean(np.array(pw) < 0.10)),
                    best_scale_p_med=float(np.nanmedian(pm)),
                ))

            # --- the swap statistic's own no-signal level ---
            for k, fn in enumerate(FUNCS):
                d = np.nanmedian(obs[:, :, :, k] - swp[:, :, :, k], axis=1)
                swap_rows.append(dict(
                    construction=f"sham_{tag}", source=src, band=band, func=fn,
                    d_med=float(np.nanmedian(d)),
                    p_cluster_pos=axis_cluster_gate(
                        d, n_perm=N_PERM, rng=np.random.default_rng(SEED))["p"],
                    p_cluster_neg=axis_cluster_gate(
                        -d, n_perm=N_PERM, rng=np.random.default_rng(SEED))["p"]))

    df = pd.DataFrame(rows)
    if not df.empty:
        for key, g in df.groupby(["construction", "source", "func"]):
            df.loc[g.index, "margin_q"] = bh_fdr(g["margin_p_gt0"].values)
            df.loc[g.index, "raw_q"] = bh_fdr(g["raw_p_gt0"].values)
    df.to_csv(OUT / "sham_per_scale.csv", index=False)
    fp = pd.DataFrame(fpr_rows)
    fp.to_csv(OUT / "sham_fpr.csv", index=False)
    pd.DataFrame(swap_rows).to_csv(OUT / "sham_swap.csv", index=False)

    pd.set_option("display.width", 220)
    if not df.empty:
        print("\n=== Q1  raw vs margin on NO-SIGNAL sham arcs (medians over scales) ===")
        agg = (df.groupby(["construction", "source", "band", "func"])
                 .agg(raw_med=("raw_med", "median"),
                      raw_sig=("raw_p_gt0", lambda x: int(np.sum(np.asarray(x) < 0.05))),
                      margin_med=("margin_med", "median"),
                      margin_sig=("margin_p_gt0", lambda x: int(np.sum(np.asarray(x) < 0.05))),
                      margin_q_sig=("margin_q", lambda x: int(np.sum(np.asarray(x) < 0.05))))
                 .reset_index())
        print(agg.to_string(index=False))
        agg.to_csv(OUT / "sham_summary.csv", index=False)
    if not fp.empty:
        print("\n=== Q2  measured FPR of the whole gate on no-signal cohorts ===")
        print(fp.to_string(index=False))
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
