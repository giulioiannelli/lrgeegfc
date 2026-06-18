#!/usr/bin/env python3
"""Decisive control for the between-region CONCENTRATION statistic (eta^2).

Context. The per-patient MAX-region-mean localization test (diag_per_patient_
localization.py) is null for the cophenet trace (0-1/10 patients per band).
An adversarial verification (workflow 2026-06-01, lens A) showed that a DIFFERENT
statistic — the between-region concentration eta^2 = between-region-variance /
total-variance of the signed trace, on the SAME within-patient region-label
shuffle null — is significant for 6-8/10 patients in alpha/beta/low_gamma/
high_gamma. eta^2 does NOT identify a hotspot region (the argmax stays
idiosyncratic; MAX is null); it only says the signed trace is non-uniform across
DK regions. THE QUESTION THIS SCRIPT ANSWERS: is that non-uniformity ANATOMICAL,
or is it the trivial spatial autocorrelation of adjacent contacts on the same
sEEG electrode shaft?

5-point critical preamble
-------------------------
1. Claim under test (lens A's): the signed trace concentrates by DK region
   within patients (eta^2 >> label-shuffle null) = anatomical regional structure.
2. Null / control: re-run the IDENTICAL eta^2 + label-shuffle test grouping by
   ELECTRODE SHAFT (probe) instead of by DK region. A shaft is a purely spatial
   unit (where the surgeon placed the electrode); contacts on one shaft are
   spatially adjacent and measure near-identical signal, so their per-pair trace
   scores s_ij are autocorrelated regardless of anatomy.
3. Strongest alternative: the regional eta^2 is fully explained by same-shaft
   spatial autocorrelation — adjacent contacts -> correlated node trace -> high
   between-group variance for ANY partition that follows the shafts, and DK
   regions largely follow shafts (a shaft passes through 1-3 adjacent regions).
4. Does the control cover it: if between-SHAFT eta^2 is significant for the same
   bands and comparable in strength to between-REGION eta^2, the concentration
   lives at the spatial-sampling level, and attributing it to DK anatomy is
   unjustified (regions inherit the shaft structure). Region<->shaft collinearity
   (normalized mutual information) quantifies how non-separable the two
   partitions are. What it cannot do: if regions and shafts were perfectly
   collinear the two tests are identical by construction — then the honest
   statement is simply "regional == spatial, inseparable", which equally defuses
   the anatomical reading.
5. Falsification of lens A's anatomical reading: if shaft-eta^2 ~ region-eta^2
   (both fire on the same bands), the eta^2 concentration is spatial/sampling,
   not anatomical localization. Region-eta^2 would only support anatomy if it
   clearly EXCEEDED shaft-eta^2 (regions add structure beyond shafts).

Both groupings use the SAME anatomical-node set (Wm/Unk endpoints dropped) so the
comparison is on identical data; only the grouping label vector differs.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.node_localization import (
    NON_ANATOMICAL, cophenet_trace_contributions, eta2,
    grassmann_trace_contributions, nmi, perm_p_eta2, shaft_of,
)
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
OUT_ROOT = ROOT / "data/audit/per_patient_localization"
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
PERM_SEED = 20260601
ALPHA = 0.05
GRASSMANN_CONTIGUOUS_K = {"beta": list(range(27, 56))}
CLUSTER_EXTENT_CSV = "data/audit/grassmann_cluster_extent/per_k_obs_p.csv"


def grassmann_k_set(band: str) -> list:
    if band in GRASSMANN_CONTIGUOUS_K:
        return GRASSMANN_CONTIGUOUS_K[band]
    df = pd.read_csv(ROOT / CLUSTER_EXTENT_CSV)
    sub = df[(df["band"] == band) & (df["obs_p_one_sided_less"] < 0.05)]
    return sorted(int(k) for k in sub["k"].values)


def contributions(probe: str, pat: str, band: str):
    if probe == "cophenet":
        return cophenet_trace_contributions(pat, band, ROOT)
    return grassmann_trace_contributions(pat, band, ROOT, grassmann_k_set(band))


def run_band(band: str, n_perm: int, probe: str = "cophenet",
             verbose: bool = True):
    rng = np.random.default_rng(PERM_SEED)
    rows = []
    for pat in COHORT:
        values, node_idx, node_region = contributions(probe, pat, band)
        values = np.asarray(values, dtype=np.float64)
        values = values - values.mean()
        node_region = np.asarray(node_region).astype(str)
        regions_df = load_channel_regions(pat)
        node_shaft = np.array([shaft_of(l) for l in regions_df["label_raw"]])

        # restrict to anatomical-endpoint contributions (same node set both ways)
        anat_node = np.array([r not in NON_ANATOMICAL for r in node_region])
        keep_contrib = anat_node[node_idx]
        v = values[keep_contrib]
        ni = node_idx[keep_contrib]

        obs_r, p_r, n_r = perm_p_eta2(v, node_region, ni, n_perm, rng)
        obs_s, p_s, n_s = perm_p_eta2(v, node_shaft, ni, n_perm, rng)
        coll = nmi(node_region[anat_node], node_shaft[anat_node])
        rows.append({
            "band": band, "patient": pat,
            "eta2_region": obs_r, "p_region": p_r, "n_regions": n_r,
            "eta2_shaft": obs_s, "p_shaft": p_s, "n_shafts": n_s,
            "region_shaft_nmi": coll,
            "pass_region": p_r < ALPHA, "pass_shaft": p_s < ALPHA,
        })
    df = pd.DataFrame(rows)
    summ = {
        "band": band,
        "n_pass_region": int(df["pass_region"].sum()),
        "n_pass_shaft": int(df["pass_shaft"].sum()),
        "median_eta2_region": float(df["eta2_region"].median()),
        "median_eta2_shaft": float(df["eta2_shaft"].median()),
        "median_region_shaft_nmi": float(df["region_shaft_nmi"].median()),
    }
    if verbose:
        print(f"  {band:11s} region-eta2 n_pass={summ['n_pass_region']}/10 "
              f"(med {summ['median_eta2_region']:.3f}) | "
              f"shaft-eta2 n_pass={summ['n_pass_shaft']}/10 "
              f"(med {summ['median_eta2_shaft']:.3f}) | "
              f"region<->shaft NMI={summ['median_region_shaft_nmi']:.2f}")
    return df, summ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--band", default=None)
    ap.add_argument("--probe", default="cophenet",
                    choices=["cophenet", "grassmann_trace"])
    args = ap.parse_args()
    bands = [args.band] if args.band else BANDS
    sfx = "" if args.probe == "cophenet" else f"_{args.probe}"
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    print(f"eta^2 concentration: DK region vs electrode shaft ({args.probe})\n")
    all_rows, summ = [], []
    for b in bands:
        df, s = run_band(b, args.n_perm, probe=args.probe)
        all_rows.append(df)
        summ.append(s)
    pd.concat(all_rows, ignore_index=True).to_csv(
        OUT_ROOT / f"concentration_region_vs_shaft{sfx}_per_patient.csv", index=False)
    sdf = pd.DataFrame(summ)
    sdf.to_csv(OUT_ROOT / f"concentration_region_vs_shaft{sfx}_summary.csv", index=False)
    print("\n==== SUMMARY (region vs shaft eta^2) ====")
    print(sdf.to_string(index=False))
    print(f"\nOutputs → {OUT_ROOT}/concentration_region_vs_shaft_*.csv")


if __name__ == "__main__":
    main()
