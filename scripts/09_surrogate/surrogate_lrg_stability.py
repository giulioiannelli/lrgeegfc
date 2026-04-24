#!/usr/bin/env python3
"""LRG stability verification: does surrogate soft-masking change
multiscale VI(k) profiles enough to affect the manuscript?

Parts:
  A — VI(k) noise floor: W vs W_soft (all patients, rest_pre)
  B — VI(k) signal: cross-phase comparison (from cached LRG)
  C — Noise-to-signal ratio + decision
  D — Random perturbation control (Pat_02 only)
  E — Extended Pat_02 check (all 4 phases)

Outputs:
  data/surrogate_validation/lrg_stability_report.txt
  data/surrogate_validation/lrg_stability_vi_profiles.csv
  data/surrogate_validation/fig_vi_noise_vs_signal.pdf  (if discrepancies)
"""

from __future__ import annotations

import gc
import time
from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PHASE_LABELS,
    nperseg_for_fs,
)
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch, band_average_msc
from lrg_eegfc.utils.fc.msc.surrogates import surrogate_msc_null
from lrg_eegfc.utils.io.patient_robust import load_patient_dataset_robust
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrgsglib.core import (
    compute_laplacian_properties,
    compute_normalized_linkage,
    compute_optimal_threshold,
    get_giant_component,
)

# ── constants ────────────────────────────────────────────────────────────────
from lrg_eegfc.config.paths import SURROGATE_VALIDATION_CACHE
OUTPUT_DIR = SURROGATE_VALIDATION_CACHE
N_SURROGATES = 100
K_VALUES = [2, 3, 4, 5, 6, 8, 10, 12, 15, 20]
BAND_NAMES = BRAIN_BANDS_NAMES
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
# Pat_06 only has rest_pre + rest_post
PATIENTS_4PHASE = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
FC_METHOD = "msc"

# ── VI functions (from deep_multiscale_tree_analysis.py) ─────────────────────

def compute_vi(labels1, labels2):
    """Variation of Information between two partitions."""
    n = len(labels1)
    if n == 0:
        return 0.0
    classes1 = np.unique(labels1)
    classes2 = np.unique(labels2)
    h1 = sum(-np.sum(labels1 == c) / n * np.log(np.sum(labels1 == c) / n)
             for c in classes1 if np.sum(labels1 == c) > 0)
    h2 = sum(-np.sum(labels2 == c) / n * np.log(np.sum(labels2 == c) / n)
             for c in classes2 if np.sum(labels2 == c) > 0)
    mi = 0.0
    for c1 in classes1:
        for c2 in classes2:
            pxy = np.sum((labels1 == c1) & (labels2 == c2)) / n
            if pxy > 0:
                px = np.sum(labels1 == c1) / n
                py = np.sum(labels2 == c2) / n
                mi += pxy * np.log(pxy / (px * py))
    return h1 + h2 - 2 * mi


def normalized_vi(labels1, labels2):
    """VI normalized by joint entropy (0 = identical, 1 = maximally different)."""
    n = len(labels1)
    if n == 0:
        return 0.0
    vi = compute_vi(labels1, labels2)
    joint = {}
    for l1, l2 in zip(labels1, labels2):
        k = (l1, l2)
        joint[k] = joint.get(k, 0) + 1
    h_joint = sum(-cnt / n * np.log(cnt / n) for cnt in joint.values() if cnt > 0)
    if h_joint == 0:
        return 0.0
    return vi / h_joint


def vi_profile(Z1, Z2, n1, n2, k_values=K_VALUES):
    """Compute VI at multiple scales k (number of clusters)."""
    if n1 != n2:
        return np.full(len(k_values), np.nan)
    vis = []
    for k in k_values:
        if k > n1:
            vis.append(np.nan)
            continue
        labels1 = fcluster(Z1, k, criterion="maxclust")
        labels2 = fcluster(Z2, k, criterion="maxclust")
        vis.append(compute_vi(labels1, labels2))
    return np.array(vis)


def nvi_profile(Z1, Z2, n1, n2, k_values=K_VALUES):
    """Compute normalized VI at multiple scales k."""
    if n1 != n2:
        return np.full(len(k_values), np.nan)
    nvis = []
    for k in k_values:
        if k > n1:
            nvis.append(np.nan)
            continue
        labels1 = fcluster(Z1, k, criterion="maxclust")
        labels2 = fcluster(Z2, k, criterion="maxclust")
        nvis.append(normalized_vi(labels1, labels2))
    return np.array(nvis)


# ── helpers ──────────────────────────────────────────────────────────────────

def _upper_tri(W):
    idx = np.triu_indices_from(W, k=1)
    return W[idx]


def _load_ts(patient, phase):
    recordings = load_patient_dataset_robust(patient, phases=[phase])
    rec = recordings[phase]
    fs = float(rec.parameters.get("fs", 2048.0))
    return rec.timeseries, fs


def _compute_empirical(X, fs, nperseg):
    freqs, Coh = compute_msc_welch(X, fs, nperseg=nperseg, noverlap=nperseg // 2)
    W_bands = band_average_msc(Coh, freqs, BRAIN_BANDS)
    del Coh
    gc.collect()
    return W_bands


def _compute_null(X, fs, nperseg, n_workers=None):
    return surrogate_msc_null(
        X, fs, BRAIN_BANDS, N_SURROGATES,
        nperseg=nperseg, noverlap=nperseg // 2,
        n_workers=n_workers,
    )


def _compute_w_soft(W, W_null_band):
    """Compute soft-masked matrix from raw W and surrogate null for one band."""
    n_surr = W_null_band.shape[0]
    p_values = (1 + np.sum(W_null_band >= W[None, :, :], axis=0)) / (n_surr + 1)
    W_soft = W * (1.0 - p_values)
    np.fill_diagonal(W_soft, 0.0)
    return W_soft


def _run_lrg(W):
    """Run LRG pipeline on adjacency matrix, return (linkage, threshold, n_nodes)."""
    W = W.copy()
    np.fill_diagonal(W, 0.0)
    G = nx.from_numpy_array(W)
    giant = get_giant_component(G)
    n_nodes = giant.number_of_nodes()
    _, _, _, Trho, _ = compute_laplacian_properties(giant)
    dists = squareform(Trho)
    linkage, labels, _ = compute_normalized_linkage(dists, giant)
    threshold, *_ = compute_optimal_threshold(linkage)
    return linkage, threshold, n_nodes


# ═════════════════════════════════════════════════════════════════════════════
#  PART A: VI(k) noise floor — W vs W_soft (all patients, rest_pre)
# ═════════════════════════════════════════════════════════════════════════════

def part_a(report, n_workers=None):
    """Compute VI(k) between raw-W and W_soft dendrograms for all patients, rest_pre."""
    report.append("=" * 70)
    report.append("PART A: VI(k) noise floor — W vs W_soft (all patients, rest_pre)")
    report.append("=" * 70)
    report.append("")

    noise_rows = []

    for patient in ALL_PATIENTS:
        print(f"\n  Part A: {patient} rest_pre")
        X, fs = _load_ts(patient, "rest_pre")
        nperseg = nperseg_for_fs(fs)
        print(f"    fs={fs}, nperseg={nperseg}, shape={X.shape}")

        W_bands = _compute_empirical(X, fs, nperseg)
        print(f"    Computing {N_SURROGATES} surrogates...")
        W_null = _compute_null(X, fs, nperseg, n_workers=n_workers)
        del X
        gc.collect()

        for band in BAND_NAMES:
            W = W_bands[band].copy()
            np.fill_diagonal(W, 0.0)
            W_soft = _compute_w_soft(W, W_null[band])

            # Run LRG on both
            Z_raw, thresh_raw, n_raw = _run_lrg(W)
            Z_soft, thresh_soft, n_soft = _run_lrg(W_soft)

            # Compute VI profiles
            vis = vi_profile(Z_raw, Z_soft, n_raw, n_soft)
            nvis = nvi_profile(Z_raw, Z_soft, n_raw, n_soft)

            max_vi = np.nanmax(vis)
            mean_vi = np.nanmean(vis)
            max_nvi = np.nanmax(nvis)

            row = {
                "patient": patient, "phase": "rest_pre", "band": band,
                "type": "noise", "pair": "W_vs_Wsoft",
                "max_vi": max_vi, "mean_vi": mean_vi, "max_nvi": max_nvi,
                "n_raw": n_raw, "n_soft": n_soft,
                "thresh_raw": thresh_raw, "thresh_soft": thresh_soft,
            }
            for ki, k in enumerate(K_VALUES):
                row[f"vi_k{k}"] = vis[ki]
                row[f"nvi_k{k}"] = nvis[ki]
            noise_rows.append(row)

            line = f"    {band:>12s}: max_VI={max_vi:.4f}  mean_VI={mean_vi:.4f}  max_NVI={max_nvi:.4f}"
            print(line)
            report.append(line)

        del W_bands, W_null
        gc.collect()
        report.append("")

    return pd.DataFrame(noise_rows)


# ═════════════════════════════════════════════════════════════════════════════
#  PART B: VI(k) signal — cross-phase comparison (from cached LRG)
# ═════════════════════════════════════════════════════════════════════════════

def part_b(report):
    """Load cached LRG results, compute cross-phase VI(k) = the actual signal."""
    report.append("=" * 70)
    report.append("PART B: VI(k) signal — cross-phase comparison (cached LRG)")
    report.append("=" * 70)
    report.append("")

    signal_rows = []

    for patient in PATIENTS_4PHASE:
        phases = list(PHASE_LABELS)
        print(f"\n  Part B: {patient}")

        for band in BAND_NAMES:
            # Load all cached LRG results for this patient/band
            results = {}
            for phase in phases:
                r = load_lrg_result(patient, phase, band, FC_METHOD)
                if r is not None:
                    results[phase] = r

            if len(results) < 2:
                continue

            # Compute VI between all phase pairs
            for p1, p2 in combinations(results.keys(), 2):
                r1, r2 = results[p1], results[p2]
                vis = vi_profile(r1.linkage_matrix, r2.linkage_matrix,
                                 r1.n_nodes, r2.n_nodes)
                nvis = nvi_profile(r1.linkage_matrix, r2.linkage_matrix,
                                   r1.n_nodes, r2.n_nodes)

                row = {
                    "patient": patient, "phase": f"{p1}_vs_{p2}", "band": band,
                    "type": "signal", "pair": f"{p1}_vs_{p2}",
                    "max_vi": np.nanmax(vis), "mean_vi": np.nanmean(vis),
                    "max_nvi": np.nanmax(nvis),
                    "n_raw": r1.n_nodes, "n_soft": r2.n_nodes,
                    "thresh_raw": r1.optimal_threshold,
                    "thresh_soft": r2.optimal_threshold,
                }
                for ki, k in enumerate(K_VALUES):
                    row[f"vi_k{k}"] = vis[ki]
                    row[f"nvi_k{k}"] = nvis[ki]
                signal_rows.append(row)

        # Summary per band
        for band in BAND_NAMES:
            band_rows = [r for r in signal_rows
                         if r["patient"] == patient and r["band"] == band]
            if band_rows:
                mean_signal = np.mean([r["mean_vi"] for r in band_rows])
                report.append(f"    {patient} {band:>12s}: mean cross-phase VI = {mean_signal:.4f}")

    # Also do Pat_06 (only rest_pre vs rest_post)
    for band in BAND_NAMES:
        r1 = load_lrg_result("Pat_06", "rest_pre", band, FC_METHOD)
        r2 = load_lrg_result("Pat_06", "rest_post", band, FC_METHOD)
        if r1 is not None and r2 is not None:
            vis = vi_profile(r1.linkage_matrix, r2.linkage_matrix,
                             r1.n_nodes, r2.n_nodes)
            nvis = nvi_profile(r1.linkage_matrix, r2.linkage_matrix,
                               r1.n_nodes, r2.n_nodes)
            row = {
                "patient": "Pat_06", "phase": "rsPre_vs_rsPost", "band": band,
                "type": "signal", "pair": "rsPre_vs_rsPost",
                "max_vi": np.nanmax(vis), "mean_vi": np.nanmean(vis),
                "max_nvi": np.nanmax(nvis),
                "n_raw": r1.n_nodes, "n_soft": r2.n_nodes,
                "thresh_raw": r1.optimal_threshold,
                "thresh_soft": r2.optimal_threshold,
            }
            for ki, k in enumerate(K_VALUES):
                row[f"vi_k{k}"] = vis[ki]
                row[f"nvi_k{k}"] = nvis[ki]
            signal_rows.append(row)

    report.append("")
    return pd.DataFrame(signal_rows)


# ═════════════════════════════════════════════════════════════════════════════
#  PART C: Noise-to-signal ratio + decision
# ═════════════════════════════════════════════════════════════════════════════

def part_c(df_noise, df_signal, report):
    """Compute noise/signal ratio per (patient, band)."""
    report.append("=" * 70)
    report.append("PART C: Noise-to-signal ratio")
    report.append("=" * 70)
    report.append("")

    report.append(f"  {'Patient':>8s}  {'Band':>12s}  {'noise_maxVI':>11s}  {'signal_meanVI':>13s}  {'ratio':>7s}  {'verdict':>8s}")
    report.append("  " + "-" * 70)

    n_pass, n_flag = 0, 0
    flagged = []

    for patient in ALL_PATIENTS:
        for band in BAND_NAMES:
            noise_row = df_noise[
                (df_noise["patient"] == patient) & (df_noise["band"] == band)
            ]
            signal_rows = df_signal[
                (df_signal["patient"] == patient) & (df_signal["band"] == band)
            ]

            if noise_row.empty or signal_rows.empty:
                report.append(f"  {patient:>8s}  {band:>12s}  {'N/A':>11s}  {'N/A':>13s}  {'N/A':>7s}  {'SKIP':>8s}")
                continue

            noise_max_vi = noise_row.iloc[0]["max_vi"]
            signal_mean_vi = signal_rows["mean_vi"].mean()

            if signal_mean_vi > 0:
                ratio = noise_max_vi / signal_mean_vi
            else:
                ratio = 0.0 if noise_max_vi == 0 else np.inf

            verdict = "PASS" if ratio < 0.10 else "FLAG"
            if verdict == "PASS":
                n_pass += 1
            else:
                n_flag += 1
                flagged.append((patient, band, noise_max_vi, signal_mean_vi, ratio))

            report.append(
                f"  {patient:>8s}  {band:>12s}  {noise_max_vi:>11.4f}  "
                f"{signal_mean_vi:>13.4f}  {ratio:>7.3f}  {verdict:>8s}"
            )

    report.append("")
    report.append(f"  PASS: {n_pass}  |  FLAG: {n_flag}  |  Total: {n_pass + n_flag}")
    report.append("")

    if n_flag == 0:
        report.append("  CONCLUSION: All cells PASS (noise/signal < 0.10).")
        report.append("  Surrogate filtering has negligible impact on multiscale VI(k).")
        report.append("  The manuscript results are stable.")
    else:
        report.append(f"  WARNING: {n_flag} cells flagged (noise/signal >= 0.10):")
        for patient, band, nv, sv, r in flagged:
            report.append(f"    {patient} {band}: noise={nv:.4f}, signal={sv:.4f}, ratio={r:.3f}")

    report.append("")
    return flagged


# ═════════════════════════════════════════════════════════════════════════════
#  PART D: Random perturbation control (Pat_02 only, rest_pre)
# ═════════════════════════════════════════════════════════════════════════════

def part_d(report, n_workers=None):
    """Compare surrogate-specific VI(k) to random perturbation VI(k)."""
    report.append("=" * 70)
    report.append("PART D: Random perturbation control (Pat_02, rest_pre)")
    report.append("=" * 70)
    report.append("")

    patient = "Pat_02"
    X, fs = _load_ts(patient, "rest_pre")
    nperseg = nperseg_for_fs(fs)

    W_bands = _compute_empirical(X, fs, nperseg)
    print(f"  Part D: Computing surrogates for {patient}...")
    W_null = _compute_null(X, fs, nperseg, n_workers=n_workers)
    del X
    gc.collect()

    n_random = 10
    rng = np.random.default_rng(42)

    for band in BAND_NAMES:
        W = W_bands[band].copy()
        np.fill_diagonal(W, 0.0)
        W_soft = _compute_w_soft(W, W_null[band])

        # Frobenius norm of the surrogate perturbation
        diff = W - W_soft
        frob_diff = np.linalg.norm(_upper_tri(diff))

        # Run LRG on raw W
        Z_raw, _, n_raw = _run_lrg(W)

        # VI(k) for surrogate perturbation
        Z_soft, _, n_soft = _run_lrg(W_soft)
        vi_surr = vi_profile(Z_raw, Z_soft, n_raw, n_soft)

        # VI(k) for random perturbations of same magnitude
        vi_random = []
        for _ in range(n_random):
            # Random symmetric perturbation
            noise = rng.normal(0, 1, W.shape)
            noise = (noise + noise.T) / 2
            np.fill_diagonal(noise, 0.0)
            # Scale to match Frobenius norm of surrogate perturbation
            noise_frob = np.linalg.norm(_upper_tri(noise))
            if noise_frob > 0:
                noise = noise * (frob_diff / noise_frob)
            W_perturbed = np.clip(W + noise, 0, 1)
            np.fill_diagonal(W_perturbed, 0.0)

            Z_pert, _, n_pert = _run_lrg(W_perturbed)
            vi_rand = vi_profile(Z_raw, Z_pert, n_raw, n_pert)
            vi_random.append(vi_rand)

        vi_random = np.array(vi_random)  # (n_random, len(K_VALUES))
        vi_random_mean = np.nanmean(vi_random, axis=0)
        vi_random_max = np.nanmax(vi_random, axis=0)

        report.append(f"  {band:>12s}: ‖W-W_soft‖_F = {frob_diff:.4f}")
        report.append(f"    Surrogate VI(k):     max={np.nanmax(vi_surr):.4f}  mean={np.nanmean(vi_surr):.4f}")
        report.append(f"    Random VI(k) mean:   max={np.nanmax(vi_random_mean):.4f}  mean={np.nanmean(vi_random_mean):.4f}")
        report.append(f"    Random VI(k) worst:  max={np.nanmax(vi_random_max):.4f}  mean={np.nanmean(vi_random_max):.4f}")

        if np.nanmax(vi_surr) <= np.nanmax(vi_random_max):
            report.append(f"    → Surrogate VI within random perturbation range (generic noise)")
        else:
            report.append(f"    → Surrogate VI exceeds random perturbations (surrogate-specific)")

    del W_bands, W_null
    gc.collect()
    report.append("")


# ═════════════════════════════════════════════════════════════════════════════
#  PART E: Extended Pat_02 check (all 4 phases)
# ═════════════════════════════════════════════════════════════════════════════

def part_e(report, n_workers=None):
    """Repeat Part A for Pat_02 across all 4 phases, re-examining the 4 discrepant cells."""
    report.append("=" * 70)
    report.append("PART E: Pat_02 all phases — re-examining 4 discrepant cells")
    report.append("=" * 70)
    report.append("")

    patient = "Pat_02"
    extra_phases = ["task_learn", "task_test", "rest_post"]  # rest_pre already done in Part A
    extra_rows = []

    for phase in extra_phases:
        print(f"\n  Part E: {patient} {phase}")
        X, fs = _load_ts(patient, phase)
        nperseg = nperseg_for_fs(fs)

        W_bands = _compute_empirical(X, fs, nperseg)
        print(f"    Computing surrogates...")
        W_null = _compute_null(X, fs, nperseg, n_workers=n_workers)
        del X
        gc.collect()

        for band in BAND_NAMES:
            W = W_bands[band].copy()
            np.fill_diagonal(W, 0.0)
            W_soft = _compute_w_soft(W, W_null[band])

            Z_raw, thresh_raw, n_raw = _run_lrg(W)
            Z_soft, thresh_soft, n_soft = _run_lrg(W_soft)

            vis = vi_profile(Z_raw, Z_soft, n_raw, n_soft)
            nvis = nvi_profile(Z_raw, Z_soft, n_raw, n_soft)

            max_vi = np.nanmax(vis)
            mean_vi = np.nanmean(vis)
            max_nvi = np.nanmax(nvis)

            row = {
                "patient": patient, "phase": phase, "band": band,
                "type": "noise", "pair": "W_vs_Wsoft",
                "max_vi": max_vi, "mean_vi": mean_vi, "max_nvi": max_nvi,
                "n_raw": n_raw, "n_soft": n_soft,
                "thresh_raw": thresh_raw, "thresh_soft": thresh_soft,
            }
            for ki, k in enumerate(K_VALUES):
                row[f"vi_k{k}"] = vis[ki]
                row[f"nvi_k{k}"] = nvis[ki]
            extra_rows.append(row)

            line = f"    {phase:>12s} {band:>12s}: max_VI={max_vi:.4f}  mean_VI={mean_vi:.4f}"
            print(line)
            report.append(line)

        del W_bands, W_null
        gc.collect()

    # Highlight the 4 previously flagged cells
    report.append("")
    report.append("  Previously flagged cells (cophenetic < 0.99):")
    flagged_cells = [
        ("rest_pre", "low_gamma", 0.8706),
        ("task_learn", "theta", 0.6699),
        ("task_test", "alpha", 0.8987),
        ("rest_post", "beta", 0.9781),
    ]
    for phase, band, old_coph in flagged_cells:
        match = [r for r in extra_rows if r["phase"] == phase and r["band"] == band]
        if match:
            r = match[0]
            report.append(
                f"    {phase} {band}: old_coph={old_coph:.4f} → VI max={r['max_vi']:.4f}, NVI max={r['max_nvi']:.4f}"
            )
        else:
            # rest_pre was in Part A, not Part E
            report.append(f"    {phase} {band}: old_coph={old_coph:.4f} → (see Part A for rest_pre)")

    report.append("")
    return pd.DataFrame(extra_rows)


# ═════════════════════════════════════════════════════════════════════════════
#  PART F: FC structural stability cross-verification (all patients, all phases)
# ═════════════════════════════════════════════════════════════════════════════

def part_f(report, n_workers=None):
    """Verify that surrogate soft-masking preserves the structural backbone of FC.

    For every (patient, phase, band): compute W and W_soft, then check:
      1. Spearman rank correlation of edge weights (global)
      2. Degree sequence correlation (hub identity preserved?)
      3. Top-10% edge overlap (backbone preserved?)
      4. Giant component size preserved?
    """
    report.append("=" * 70)
    report.append("PART F: FC structural stability (all patients, all available phases)")
    report.append("=" * 70)
    report.append("")

    fc_rows = []

    for patient in ALL_PATIENTS:
        # Determine available phases
        if patient == "Pat_06":
            phases = ["rest_pre", "rest_post"]
        else:
            phases = list(PHASE_LABELS)

        for phase in phases:
            print(f"  Part F: {patient} {phase}")
            try:
                X, fs = _load_ts(patient, phase)
            except Exception as e:
                print(f"    SKIP: {e}")
                continue
            nperseg = nperseg_for_fs(fs)

            W_bands = _compute_empirical(X, fs, nperseg)
            print(f"    Computing surrogates...")
            W_null = _compute_null(X, fs, nperseg, n_workers=n_workers)
            del X
            gc.collect()

            for band in BAND_NAMES:
                W = W_bands[band].copy()
                np.fill_diagonal(W, 0.0)
                W_soft = _compute_w_soft(W, W_null[band])

                ut_W = _upper_tri(W)
                ut_Ws = _upper_tri(W_soft)

                # 1. Spearman rank correlation
                rho_spearman, _ = spearmanr(ut_W, ut_Ws)

                # 2. Degree sequence correlation (hub preservation)
                deg_W = W.sum(axis=1)
                deg_Ws = W_soft.sum(axis=1)
                rho_degree, _ = spearmanr(deg_W, deg_Ws)

                # 3. Top-10% edge overlap (backbone)
                n_edges = len(ut_W)
                top_n = max(1, n_edges // 10)
                top_W = set(np.argsort(ut_W)[-top_n:])
                top_Ws = set(np.argsort(ut_Ws)[-top_n:])
                overlap = len(top_W & top_Ws) / top_n

                # 4. Giant component size
                G_raw = nx.from_numpy_array(W)
                G_soft = nx.from_numpy_array(W_soft)
                gc_raw = max(nx.connected_components(G_raw), key=len)
                gc_soft = max(nx.connected_components(G_soft), key=len)

                fc_rows.append({
                    "patient": patient, "phase": phase, "band": band,
                    "spearman_edges": rho_spearman,
                    "spearman_degree": rho_degree,
                    "top10_overlap": overlap,
                    "gc_raw": len(gc_raw), "gc_soft": len(gc_soft),
                })

            del W_bands, W_null
            gc.collect()

    df_fc = pd.DataFrame(fc_rows)

    # Save CSV
    fc_csv = OUTPUT_DIR / "fc_structural_stability.csv"
    df_fc.to_csv(fc_csv, index=False)
    report.append(f"  Saved: {fc_csv}")
    report.append("")

    # Summary statistics
    report.append("  Summary across all (patient, phase, band) cells:")
    report.append(f"    Edge Spearman:   min={df_fc['spearman_edges'].min():.4f}  "
                  f"median={df_fc['spearman_edges'].median():.4f}  "
                  f"mean={df_fc['spearman_edges'].mean():.4f}")
    report.append(f"    Degree Spearman: min={df_fc['spearman_degree'].min():.4f}  "
                  f"median={df_fc['spearman_degree'].median():.4f}  "
                  f"mean={df_fc['spearman_degree'].mean():.4f}")
    report.append(f"    Top-10% overlap: min={df_fc['top10_overlap'].min():.4f}  "
                  f"median={df_fc['top10_overlap'].median():.4f}  "
                  f"mean={df_fc['top10_overlap'].mean():.4f}")

    # Check giant component preservation
    gc_diff = df_fc["gc_raw"] != df_fc["gc_soft"]
    if gc_diff.any():
        report.append(f"    Giant component size DIFFERS in {gc_diff.sum()} cells:")
        for _, row in df_fc[gc_diff].iterrows():
            report.append(f"      {row['patient']} {row['phase']} {row['band']}: "
                          f"raw={row['gc_raw']}, soft={row['gc_soft']}")
    else:
        report.append("    Giant component size: IDENTICAL in all cells")

    report.append("")

    # Flag any structural instability
    weak_cells = df_fc[
        (df_fc["spearman_edges"] < 0.99) |
        (df_fc["spearman_degree"] < 0.99) |
        (df_fc["top10_overlap"] < 0.90)
    ]
    if weak_cells.empty:
        report.append("  CONCLUSION: FC structure is fully preserved by surrogate masking.")
        report.append("    - All edge Spearman correlations >= 0.99")
        report.append("    - All degree sequence correlations >= 0.99 (hub identity stable)")
        report.append("    - All top-10% edge overlaps >= 0.90 (backbone preserved)")
    else:
        report.append(f"  WARNING: {len(weak_cells)} cells show structural sensitivity:")
        for _, row in weak_cells.iterrows():
            report.append(
                f"    {row['patient']} {row['phase']} {row['band']}: "
                f"ρ_edge={row['spearman_edges']:.4f}  ρ_deg={row['spearman_degree']:.4f}  "
                f"top10={row['top10_overlap']:.3f}"
            )

    report.append("")
    return df_fc


# ═════════════════════════════════════════════════════════════════════════════
#  SUMMARY FIGURE
# ═════════════════════════════════════════════════════════════════════════════

def make_summary_figure(df_noise, df_signal):
    """Noise vs signal comparison figure."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.ravel()

    for idx, band in enumerate(BAND_NAMES):
        ax = axes[idx]

        # Noise: max VI per patient
        noise_band = df_noise[df_noise["band"] == band]
        if not noise_band.empty:
            for _, row in noise_band.iterrows():
                vis = [row.get(f"vi_k{k}", np.nan) for k in K_VALUES]
                ax.plot(K_VALUES, vis, "r-", alpha=0.3, linewidth=1)
            # Plot mean noise
            noise_vis_mean = noise_band[[f"vi_k{k}" for k in K_VALUES]].mean()
            ax.plot(K_VALUES, noise_vis_mean.values, "r-", linewidth=2, label="Noise (W vs W_soft)")

        # Signal: mean cross-phase VI per patient
        signal_band = df_signal[df_signal["band"] == band]
        if not signal_band.empty:
            for _, row in signal_band.iterrows():
                vis = [row.get(f"vi_k{k}", np.nan) for k in K_VALUES]
                ax.plot(K_VALUES, vis, "b-", alpha=0.1, linewidth=0.5)
            signal_vis_mean = signal_band[[f"vi_k{k}" for k in K_VALUES]].mean()
            ax.plot(K_VALUES, signal_vis_mean.values, "b-", linewidth=2, label="Signal (cross-phase)")

        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=13)
        ax.set_xlabel("k (clusters)")
        ax.set_ylabel("VI")
        ax.legend(fontsize=8)
        ax.set_xlim(K_VALUES[0], K_VALUES[-1])

    fig.suptitle("LRG Stability: Surrogate noise (red) vs Cross-phase signal (blue)", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = OUTPUT_DIR / "fig_vi_noise_vs_signal.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  → saved {out}")


# ═════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    import argparse
    parser = argparse.ArgumentParser(description="LRG stability verification")
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    report = [
        "LRG STABILITY VERIFICATION REPORT",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"N_surrogates: {N_SURROGATES}",
        f"K_values: {K_VALUES}",
        "",
    ]

    t0 = time.time()

    # Part A: noise floor
    print("══════ PART A: VI(k) noise floor ══════")
    df_noise_a = part_a(report, n_workers=args.workers)

    # Part B: signal (from cache)
    print("\n══════ PART B: VI(k) signal ══════")
    df_signal = part_b(report)

    # Part C: noise/signal ratio
    print("\n══════ PART C: Noise/signal ratio ══════")
    flagged = part_c(df_noise_a, df_signal, report)

    # Part D: random perturbation control
    print("\n══════ PART D: Random perturbation control ══════")
    part_d(report, n_workers=args.workers)

    # Part E: extended Pat_02
    print("\n══════ PART E: Pat_02 all phases ══════")
    df_noise_e = part_e(report, n_workers=args.workers)

    # Part F: FC structural stability (all patients, all phases)
    print("\n══════ PART F: FC structural stability ══════")
    df_fc = part_f(report, n_workers=args.workers)

    # Combine all noise rows
    df_noise_all = pd.concat([df_noise_a, df_noise_e], ignore_index=True)

    # Save CSV
    df_all = pd.concat([df_noise_all, df_signal], ignore_index=True)
    csv_path = OUTPUT_DIR / "lrg_stability_vi_profiles.csv"
    df_all.to_csv(csv_path, index=False)
    print(f"\n  → saved {csv_path}")

    # Summary figure
    make_summary_figure(df_noise_all, df_signal)

    elapsed = time.time() - t0
    report.append(f"Total runtime: {elapsed / 60:.1f} min")
    print(f"\nTotal runtime: {elapsed / 60:.1f} min")

    report_path = OUTPUT_DIR / "lrg_stability_report.txt"
    report_path.write_text("\n".join(report))
    print(f"Report saved: {report_path}")


if __name__ == "__main__":
    main()
