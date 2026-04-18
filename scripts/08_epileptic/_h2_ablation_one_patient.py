#!/usr/bin/env python3
"""Run ablation H2a computation for a single patient (memory-isolated)."""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io import load_epileptic_nodes
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import compute_lrg_analysis, load_lrg_result

BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)
H_GRID = np.linspace(0.01, 0.99, 60)

OUT = FIGURES_ROOT / "epileptic_h2_stratified"
OUT.mkdir(parents=True, exist_ok=True)


def load_ch(pat):
    p = SEEG_DATAPATH / pat / "channel_labels.csv"
    if not p.exists():
        p = SEEG_DATAPATH / pat / "channel_labels.txt"
    with open(p) as f:
        first = f.readline().strip()
    skip = 1 if first.lower() == "label" else 0
    df = pd.read_csv(p, header=None, skiprows=skip)
    return [str(l).strip('"').split(",")[0].strip().replace(" ", "")
            for l in df.iloc[:, 0]]


def compute_vi(lab1, lab2):
    n = len(lab1)
    if n == 0:
        return 0.0
    c1 = np.unique(lab1)
    c2 = np.unique(lab2)
    h1 = -sum(p * np.log(p) for p in (np.sum(lab1 == c) / n for c in c1) if p > 0)
    h2 = -sum(p * np.log(p) for p in (np.sum(lab2 == c) / n for c in c2) if p > 0)
    mi = 0.0
    for a in c1:
        for b in c2:
            pxy = np.sum((lab1 == a) & (lab2 == b)) / n
            if pxy > 0:
                px = np.sum(lab1 == a) / n
                py = np.sum(lab2 == b) / n
                mi += pxy * np.log(pxy / (px * py))
    return h1 + h2 - 2 * mi


def h2a_auc_from_linkages(Z_pre, Z_post, Z_tt):
    """Return integrated positive H2a = d_VI(rsPre,rsPost) - d_VI(taskTest,rsPost)."""
    vals = []
    for h in H_GRID:
        l_pre = fcluster(Z_pre, t=h, criterion="distance")
        l_post = fcluster(Z_post, t=h, criterion="distance")
        l_tt = fcluster(Z_tt, t=h, criterion="distance")
        vals.append(compute_vi(l_pre, l_post) - compute_vi(l_tt, l_post))
    vals = np.asarray(vals)
    return float(np.trapz(np.clip(vals, 0, None), H_GRID)), vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("patient")
    args = ap.parse_args()
    pat = args.patient

    ch = load_ch(pat)
    N = len(ch)
    epi_set = set(load_epileptic_nodes(pat))
    epi_mask = np.array([l in epi_set for l in ch])
    epi_idx = np.where(epi_mask)[0]

    ee_mask_full = epi_mask[:, None] & epi_mask[None, :]
    np.fill_diagonal(ee_mask_full, False)

    rows = []
    for band in BANDS:
        # Original linkages from cache
        orig = {}
        for phase in PHASES:
            lrg = load_lrg_result(pat, phase, band, "imcoh_abs")
            if lrg is None or lrg.n_nodes != N:
                continue
            orig[phase] = lrg.linkage_matrix
        if not all(p in orig for p in ("rsPre", "rsPost", "taskTest")):
            continue

        # Ablated: zero epi-epi entries, recompute LRG per phase (use_cache=False)
        abl = {}
        for phase in PHASES:
            A = load_fc_matrix(pat, phase, band, "imcoh_abs")
            if A is None or A.shape[0] != N:
                continue
            A_abl = A.copy()
            A_abl[ee_mask_full] = 0.0  # zero epi-epi edges (symmetric)
            # Compute LRG (no cache write, no read)
            res = compute_lrg_analysis(
                A_abl, pat, phase, band, "imcoh_abs",
                use_cache=False, overwrite_cache=False, verbose=False,
            )
            abl[phase] = res.linkage_matrix
        if not all(p in abl for p in ("rsPre", "rsPost", "taskTest")):
            continue

        auc_o, _ = h2a_auc_from_linkages(orig["rsPre"], orig["rsPost"], orig["taskTest"])
        auc_a, _ = h2a_auc_from_linkages(abl["rsPre"], abl["rsPost"], abl["taskTest"])

        rows.append(dict(
            patient=pat, band=band,
            h2a_auc_original=auc_o, h2a_auc_ablated=auc_a,
            delta_abs=auc_a - auc_o,
            delta_pct=(auc_a - auc_o) / abs(auc_o) * 100 if abs(auc_o) > 1e-10 else np.nan,
            n_epi=int(epi_mask.sum()),
        ))
        print(f"  {pat} {band}: H2a orig={auc_o:.4f} abl={auc_a:.4f} "
              f"Δ={(auc_a - auc_o):+.4f} ({rows[-1]['delta_pct']:+.1f}%)",
              flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / f"_ablation_{pat}.csv", index=False)
    print(f"{pat}: wrote {len(rows)} rows")


if __name__ == "__main__":
    main()
