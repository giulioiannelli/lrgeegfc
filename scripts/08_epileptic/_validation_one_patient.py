#!/usr/bin/env python3
"""Run V1+V2+V5 for a single patient; write results to CSV.

Invoked by driver script once per patient so Python releases all memory
between patients (avoids OOM).
"""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import mannwhitneyu, binomtest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io import load_epileptic_nodes
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result

BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)
OUT = FIGURES_ROOT / "epileptic_imcoh_validation"
OUT.mkdir(parents=True, exist_ok=True)


def load_ch(pat):
    p = SEEG_DATAPATH / pat / "channel_labels.csv"
    if not p.exists():
        p = SEEG_DATAPATH / pat / "channel_labels.txt"
    with open(p) as f:
        first = f.readline().strip()
    skip = 1 if first.lower() == "label" else 0
    df = pd.read_csv(p, header=None, skiprows=skip)
    return [str(l).strip('"').split(",")[0].strip().replace(" ", "") for l in df.iloc[:, 0]]


def probe(label):
    m = re.match(r"([A-Za-z]+'?)", label)
    return m.group(1) if m else label


def pair_values(M, idx_a, idx_b, probes_arr, cross_probe_only, symmetric):
    """Get upper-triangle pair values. symmetric=True for within-group (no i<=j dup)."""
    vals = []
    if symmetric:
        for ii in range(len(idx_a)):
            for jj in range(ii + 1, len(idx_a)):
                i, j = idx_a[ii], idx_a[jj]
                if cross_probe_only and probes_arr[i] == probes_arr[j]:
                    continue
                vals.append(M[i, j])
    else:
        for i in idx_a:
            for j in idx_b:
                if i >= j:
                    continue
                if cross_probe_only and probes_arr[i] == probes_arr[j]:
                    continue
                vals.append(M[i, j])
    return np.asarray(vals, dtype=np.float64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("patient")
    args = ap.parse_args()
    pat = args.patient

    epi_all = load_epileptic_nodes(pat)
    ch = load_ch(pat)
    N = len(ch)
    epi_set = set(epi_all)
    epi_mask = np.array([l in epi_set for l in ch])
    epi_idx = np.where(epi_mask)[0]
    non_idx = np.where(~epi_mask)[0]
    probes_arr = np.array([probe(l) for l in ch])

    print(f"{pat}: N={N}, epi={len(epi_idx)}, non={len(non_idx)}", flush=True)

    v1_rows, v2_rows, v5_rows = [], [], []

    for phase in PHASES:
        for band in BANDS:
            # V1: FC edge strength
            # imcoh_abs = |ImCoh| (Ewald 2012) — loader applies np.abs on the
            # signed Nolte-2004 cache, matching the paper narrative.
            A = load_fc_matrix(pat, phase, band, "imcoh_abs")
            if A is not None and A.shape[0] == N:
                A_abs = A.copy()
                np.fill_diagonal(A_abs, 0)

                ee_cp = pair_values(A_abs, epi_idx, epi_idx, probes_arr, True, True)
                nn_cp = pair_values(A_abs, non_idx, non_idx, probes_arr, True, True)
                ee_all = pair_values(A_abs, epi_idx, epi_idx, probes_arr, False, True)
                nn_all = pair_values(A_abs, non_idx, non_idx, probes_arr, False, True)

                if len(ee_cp) >= 2 and len(nn_cp) >= 2:
                    ratio_cp = float(ee_cp.mean()) / (float(nn_cp.mean()) + 1e-30)
                    ratio_all = float(ee_all.mean()) / (float(nn_all.mean()) + 1e-30) if len(nn_all) > 0 else np.nan
                    _, p_cp = mannwhitneyu(ee_cp, nn_cp, alternative="greater", method="asymptotic")
                    v1_rows.append(dict(
                        patient=pat, phase=phase, band=band,
                        ee_mean=float(ee_cp.mean()), nn_mean=float(nn_cp.mean()),
                        ratio=ratio_cp, ratio_all=ratio_all, p_value=p_cp,
                        n_ee=len(ee_cp), n_nn=len(nn_cp),
                    ))
                del A, A_abs, ee_cp, nn_cp, ee_all, nn_all

            # V2 + V5: ultrametric
            lrg = load_lrg_result(pat, phase, band, "imcoh_abs")
            if lrg is not None:
                um = squareform(lrg.ultrametric_matrix)

                # V2: cross-probe NN enrichment
                nn_is_epi = 0
                total_exp = 0.0
                valid = 0
                for i in epi_idx:
                    cp_mask = probes_arr != probes_arr[i]
                    cp_mask[i] = False
                    cp_indices = np.where(cp_mask)[0]
                    if len(cp_indices) == 0:
                        continue
                    valid += 1
                    dists = um[i, cp_indices]
                    nn_local = cp_indices[np.argmin(dists)]
                    if epi_mask[nn_local]:
                        nn_is_epi += 1
                    total_exp += int(epi_mask[cp_indices].sum()) / len(cp_indices)

                if valid > 0:
                    exp_frac = total_exp / valid
                    obs_frac = nn_is_epi / valid
                    p_val = binomtest(nn_is_epi, valid, exp_frac, alternative="greater").pvalue
                    v2_rows.append(dict(
                        patient=pat, phase=phase, band=band,
                        nn_is_epi=nn_is_epi, n_epi=valid,
                        obs_frac=obs_frac, expected_frac=exp_frac,
                        enrichment=obs_frac / (exp_frac + 1e-30),
                        p_value=p_val,
                    ))

                # V5: cross-probe ultrametric ratio
                ee_um = pair_values(um, epi_idx, epi_idx, probes_arr, True, True)
                nn_um = pair_values(um, non_idx, non_idx, probes_arr, True, True)
                if len(ee_um) > 0 and len(nn_um) > 0:
                    um_ratio = float(ee_um.mean()) / (float(nn_um.mean()) + 1e-30)
                    v5_rows.append(dict(
                        patient=pat, phase=phase, band=band, um_ratio=um_ratio,
                    ))
                del um, lrg, ee_um, nn_um

    # Append to global CSVs (locked per-patient; driver merges)
    pd.DataFrame(v1_rows).to_csv(OUT / f"_v1_{pat}.csv", index=False)
    pd.DataFrame(v2_rows).to_csv(OUT / f"_v2_{pat}.csv", index=False)
    pd.DataFrame(v5_rows).to_csv(OUT / f"_v5_{pat}.csv", index=False)
    print(f"{pat}: wrote {len(v1_rows)} V1, {len(v2_rows)} V2, {len(v5_rows)} V5 rows",
          flush=True)


if __name__ == "__main__":
    main()
