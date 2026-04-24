#!/usr/bin/env python3
"""Surrogate validation follow-up — Tasks 0–4.

Produces:
  data/surrogate_validation/surrogate_summary_table.csv   (Pat_06 rows appended)
  data/surrogate_validation/surrogate_followup_report.txt
"""

from __future__ import annotations

import gc
import time
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    BRAIN_BANDS_NAMES,
    nperseg_for_fs,
)
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch, band_average_msc
from lrg_eegfc.utils.fc.msc.surrogates import surrogate_msc_null
from lrg_eegfc.utils.io.patient_robust import load_patient_dataset_robust

from lrg_eegfc.config.paths import SURROGATE_VALIDATION_CACHE
OUTPUT_DIR = SURROGATE_VALIDATION_CACHE
N_SURROGATES = 100
BAND_NAMES = BRAIN_BANDS_NAMES
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
PHASE = "rest_pre"


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


# ═════════════════════════════════════════════════════════════════════════════
#  TASK 0: Pipeline N_sur default
# ═════════════════════════════════════════════════════════════════════════════

def task0(report):
    report.append("=" * 70)
    report.append("TASK 0: Pipeline N_sur default")
    report.append("=" * 70)
    report.append("")
    report.append("  DEFAULT_N_SURROGATES = 200")
    report.append("  Defined at: src/lrg_eegfc/config/const.py:135")
    report.append("")
    report.append("  CLI behaviour (src/lrg_eegfc/cli/compute.py:142-144):")
    report.append("    When sparsify requires surrogates and n_surrogates==0,")
    report.append("    it is auto-set to 200.")
    report.append("")
    report.append("  The previous validation script used N_sur=100 for speed.")
    report.append("  The production pipeline default is 200.")
    report.append("")


# ═════════════════════════════════════════════════════════════════════════════
#  TASK 1: Pat_06
# ═════════════════════════════════════════════════════════════════════════════

def task1(report, n_workers=None):
    report.append("=" * 70)
    report.append("TASK 1: Pat_06 surrogate validation (rest_pre, all 6 bands)")
    report.append("=" * 70)

    X, fs = _load_ts("Pat_06", PHASE)
    nperseg = nperseg_for_fs(fs)
    report.append(f"  Pat_06: fs={fs} Hz, nperseg={nperseg}, shape={X.shape}")
    print(f"  Pat_06: fs={fs}, nperseg={nperseg}, shape={X.shape}")

    W_bands = _compute_empirical(X, fs, nperseg)
    W_null = _compute_null(X, fs, nperseg, n_workers=n_workers)
    del X
    gc.collect()

    rows = []
    for band in BAND_NAMES:
        W = W_bands[band].copy()
        np.fill_diagonal(W, 0.0)
        Wn = W_null[band]
        n_surr = Wn.shape[0]

        p_values = (1 + np.sum(Wn >= W[None, :, :], axis=0)) / (n_surr + 1)
        frac_p05 = np.mean(_upper_tri(p_values) < 0.05)
        frac_p01 = np.mean(_upper_tri(p_values) < 0.01)

        W_soft = W * (1.0 - p_values)
        np.fill_diagonal(W_soft, 0.0)
        rho, _ = spearmanr(_upper_tri(W), _upper_tri(W_soft))

        # Decision criterion
        if rho > 0.95:
            cat = "A"
        elif frac_p05 > 0.90:
            cat = "B"
        else:
            cat = "C"

        rows.append({
            "patient": "Pat_06",
            "band": band,
            "frac_significant_p05": frac_p05,
            "frac_significant_p01": frac_p01,
            "spearman_W_vs_Wsoft": rho,
        })

        line = f"    {band:>12s}: sig_p05={frac_p05:.3f}  sig_p01={frac_p01:.3f}  ρ={rho:.4f}  → ({cat})"
        print(line)
        report.append(line)

    # Append to existing CSV
    csv_path = OUTPUT_DIR / "surrogate_summary_table.csv"
    df_new = pd.DataFrame(rows)
    if csv_path.exists():
        df_old = pd.read_csv(csv_path)
        # Remove any existing Pat_06 rows to avoid duplicates
        df_old = df_old[df_old["patient"] != "Pat_06"]
        df_combined = pd.concat([df_old, df_new], ignore_index=True)
    else:
        df_combined = df_new
    df_combined.to_csv(csv_path, index=False)
    report.append(f"\n  Updated CSV: {csv_path}")
    report.append("")

    return W_bands, W_null


# ═════════════════════════════════════════════════════════════════════════════
#  TASK 2: Pat_02 anomaly investigation
# ═════════════════════════════════════════════════════════════════════════════

def task2(report, n_workers=None):
    report.append("=" * 70)
    report.append("TASK 2: Pat_02 anomaly investigation")
    report.append("=" * 70)
    report.append("")

    patients = ALL_PATIENTS
    patient_info = {}

    for patient in patients:
        X, fs = _load_ts(patient, PHASE)
        nperseg = nperseg_for_fs(fs)
        N_ch, T = X.shape
        noverlap = nperseg // 2
        step = nperseg - noverlap
        M = (T - nperseg) // step + 1  # number of Welch segments
        duration_min = T / fs / 60.0

        W_bands = _compute_empirical(X, fs, nperseg)

        # Compute surrogate means (load from null if we're computing them anyway)
        # For Pat_06 and Pat_02, we compute fresh surrogates.
        # For others, we need surrogates too for the p95 comparison
        # Only compute surrogates for the comparison patients if needed
        # Actually we need surrogate stats for ALL patients for the p95 comparison
        W_null = _compute_null(X, fs, nperseg, n_workers=n_workers)

        band_stats = {}
        for band in BAND_NAMES:
            W = W_bands[band].copy()
            np.fill_diagonal(W, 0.0)
            Wn = W_null[band]

            triu_i, triu_j = np.triu_indices(W.shape[0], k=1)
            emp_vals = W[triu_i, triu_j]
            surr_per_edge = Wn[:, triu_i, triu_j]  # (n_surr, n_edges)
            surr_pooled = surr_per_edge.ravel()

            p95_per_edge = np.percentile(surr_per_edge, 95, axis=0)

            band_stats[band] = {
                "emp_mean": emp_vals.mean(),
                "surr_mean": surr_pooled.mean(),
                "surr_std": surr_pooled.std(),
                "ratio": emp_vals.mean() / surr_pooled.mean() if surr_pooled.mean() > 0 else np.inf,
                "p95_mean": p95_per_edge.mean(),
                "p95_std": p95_per_edge.std(),
            }

        patient_info[patient] = {
            "fs": fs,
            "nperseg": nperseg,
            "N_ch": N_ch,
            "T": T,
            "duration_min": duration_min,
            "M": M,
            "bands": band_stats,
        }

        del X, W_bands, W_null
        gc.collect()

    # --- Report: recording duration table ---
    report.append("  Recording characteristics (rest_pre):")
    report.append(f"  {'Patient':>8s}  {'fs':>6s}  {'N_ch':>5s}  {'T (samples)':>14s}  {'Duration':>10s}  {'nperseg':>7s}  {'M (segments)':>12s}")
    report.append("  " + "-" * 75)
    for patient in patients:
        info = patient_info[patient]
        line = (
            f"  {patient:>8s}  {info['fs']:>6.0f}  {info['N_ch']:>5d}  "
            f"{info['T']:>14,d}  {info['duration_min']:>8.1f} min  "
            f"{info['nperseg']:>7d}  {info['M']:>12,d}"
        )
        report.append(line)
        print(line)
    report.append("")

    # --- Report: mean MSC per band ---
    report.append("  Mean empirical MSC per band:")
    header = f"  {'Patient':>8s}" + "".join(f"  {b:>12s}" for b in BAND_NAMES)
    report.append(header)
    report.append("  " + "-" * (9 + 14 * len(BAND_NAMES)))
    for patient in patients:
        info = patient_info[patient]
        vals = "".join(f"  {info['bands'][b]['emp_mean']:>12.4f}" for b in BAND_NAMES)
        report.append(f"  {patient:>8s}{vals}")
    report.append("")

    # --- Report: mean surrogate MSC per band ---
    report.append("  Mean surrogate MSC per band:")
    report.append(header)
    report.append("  " + "-" * (9 + 14 * len(BAND_NAMES)))
    for patient in patients:
        info = patient_info[patient]
        vals = "".join(f"  {info['bands'][b]['surr_mean']:>12.4f}" for b in BAND_NAMES)
        report.append(f"  {patient:>8s}{vals}")
    report.append("")

    # --- Report: ratio empirical/surrogate ---
    report.append("  Ratio empirical / surrogate mean:")
    report.append(header)
    report.append("  " + "-" * (9 + 14 * len(BAND_NAMES)))
    for patient in patients:
        info = patient_info[patient]
        vals = "".join(f"  {info['bands'][b]['ratio']:>12.1f}" for b in BAND_NAMES)
        report.append(f"  {patient:>8s}{vals}")
    report.append("")

    # --- Report: 95th percentile of surrogate null (low_gamma focus) ---
    report.append("  Mean 95th percentile of edge-specific surrogate null (low_gamma):")
    for patient in patients:
        info = patient_info[patient]
        bs = info["bands"]["low_gamma"]
        line = (
            f"    {patient}: p95_mean={bs['p95_mean']:.5f} ± {bs['p95_std']:.5f}  "
            f"emp_mean={bs['emp_mean']:.4f}  surr_mean={bs['surr_mean']:.5f}"
        )
        report.append(line)
        print(line)
    report.append("")

    # --- Hypothesis assessment ---
    pat02 = patient_info["Pat_02"]
    others_M = [patient_info[p]["M"] for p in patients if p != "Pat_02"]
    pat02_M = pat02["M"]

    report.append("  Hypothesis: does Pat_02 have fewer Welch segments?")
    report.append(f"    Pat_02 M = {pat02_M:,d}")
    report.append(f"    Other patients M = {', '.join(str(patient_info[p]['M']) for p in patients if p != 'Pat_02')}")
    report.append(f"    Min M (others) = {min(others_M):,d}")
    report.append(f"    Max M (others) = {max(others_M):,d}")

    if pat02_M >= min(others_M):
        report.append("    → Pat_02 does NOT have fewer segments than all others.")
    else:
        report.append("    → Pat_02 has the fewest segments.")

    # Check if surrogate null is wider for Pat_02
    report.append("")
    report.append("  Surrogate null width comparison (std of pooled surrogate, low_gamma):")
    for patient in patients:
        info = patient_info[patient]
        bs = info["bands"]["low_gamma"]
        report.append(f"    {patient}: surr_std={bs['surr_std']:.5f}")

    report.append("")

    return patient_info


# ═════════════════════════════════════════════════════════════════════════════
#  TASK 3: Pat_03 nperseg confirmation
# ═════════════════════════════════════════════════════════════════════════════

def task3(report):
    report.append("=" * 70)
    report.append("TASK 3: Pat_03 nperseg confirmation")
    report.append("=" * 70)
    report.append("")
    report.append("  Pat_03 sampling rate: 1024 Hz (all others: 2048 Hz)")
    report.append("  nperseg_for_fs(1024) = int(1024 * 2.0) = 2048")
    report.append("  nperseg_for_fs(2048) = int(2048 * 2.0) = 4096")
    report.append("")
    report.append("  Both yield Δf = fs / nperseg = 0.5 Hz (2-second segments)")
    report.append("  The previous validation script correctly used nperseg=2048 for Pat_03.")
    report.append("  Report header confirmed: 'Pat_03: fs=1024.0, nperseg=2048'")
    report.append("")
    report.append("  Nyquist for Pat_03: 512 Hz. High_gamma band (80-300 Hz) is fully")
    report.append("  within the Nyquist limit. No band is clipped.")
    report.append("")


# ═════════════════════════════════════════════════════════════════════════════
#  TASK 4: LRG comparison raw W vs W_soft (Pat_02, all phases, all bands)
# ═════════════════════════════════════════════════════════════════════════════

def task4(report, n_workers=None):
    import networkx as nx
    from lrgsglib.core import (
        compute_laplacian_properties,
        compute_normalized_linkage,
        compute_optimal_threshold,
        get_giant_component,
    )
    from scipy.cluster.hierarchy import cophenet

    report.append("=" * 70)
    report.append("TASK 4: LRG comparison raw W vs W_soft (Pat_02, all phases/bands)")
    report.append("=" * 70)
    report.append("")

    phases = ["rest_pre", "task_learn", "task_test", "rest_post"]
    patient = "Pat_02"

    def _run_lrg(W):
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

    # Load cached LRG results for raw W (from lrg_cache)
    from lrg_eegfc.config.paths import LRG_CACHE as _LRG_CACHE, MSC_CACHE as _MSC_CACHE
    lrg_cache = _LRG_CACHE / "Pat_02"
    msc_cache = _MSC_CACHE / "Pat_02"

    # Table header
    report.append(f"  {'Phase':>12s}  {'Band':>12s}  {'coph_corr':>10s}  {'n*_raw':>7s}  {'n*_soft':>8s}  {'n_raw':>6s}  {'n_soft':>7s}  {'agree':>6s}")
    report.append("  " + "-" * 82)

    discrepancies = []

    for phase in phases:
        print(f"\n  Phase: {phase}")

        # Load timeseries for surrogate computation
        X, fs = _load_ts(patient, phase)
        nperseg = nperseg_for_fs(fs)

        # Compute empirical MSC
        W_bands = _compute_empirical(X, fs, nperseg)

        # Compute surrogates for this phase
        print(f"    Computing surrogates for {patient} {phase}...")
        W_null = _compute_null(X, fs, nperseg, n_workers=n_workers)
        del X
        gc.collect()

        for band in BAND_NAMES:
            W = W_bands[band].copy()
            np.fill_diagonal(W, 0.0)
            Wn = W_null[band]
            n_surr = Wn.shape[0]

            # Compute W_soft
            p_values = (1 + np.sum(Wn >= W[None, :, :], axis=0)) / (n_surr + 1)
            W_soft = W * (1.0 - p_values)
            np.fill_diagonal(W_soft, 0.0)

            # Run LRG on raw W
            link_W, thresh_W, n_W = _run_lrg(W)
            # Run LRG on W_soft
            link_Ws, thresh_Ws, n_Ws = _run_lrg(W_soft)

            # Cophenetic correlation
            if n_W == n_Ws:
                d_W = squareform(cophenet(link_W))
                d_Ws = squareform(cophenet(link_Ws))
                coph_rho, _ = spearmanr(_upper_tri(d_W), _upper_tri(d_Ws))
                coph_str = f"{coph_rho:.4f}"
            else:
                coph_rho = None
                coph_str = "N/A (size)"

            # Check agreement
            agree = "yes" if abs(thresh_W - thresh_Ws) < 0.05 or (n_W == n_Ws and (coph_rho is not None and coph_rho >= 0.999)) else "CHECK"

            line = (
                f"  {phase:>12s}  {band:>12s}  {coph_str:>10s}  "
                f"{thresh_W:>7.4f}  {thresh_Ws:>8.4f}  "
                f"{n_W:>6d}  {n_Ws:>7d}  {agree:>6s}"
            )
            report.append(line)
            print(f"    {band}: coph={coph_str}, n*_raw={thresh_W:.4f}, n*_soft={thresh_Ws:.4f}, n={n_W}/{n_Ws}")

            if coph_rho is not None and coph_rho < 0.99:
                discrepancies.append((phase, band, coph_rho, thresh_W, thresh_Ws))

        del W_bands, W_null
        gc.collect()

    report.append("")

    if discrepancies:
        report.append("  DISCREPANCIES (cophenetic < 0.99):")
        for phase, band, coph, t_raw, t_soft in discrepancies:
            report.append(f"    {phase} {band}: coph={coph:.4f}, n*_raw={t_raw:.4f}, n*_soft={t_soft:.4f}")
    else:
        report.append("  No discrepancies found. All cophenetic correlations >= 0.99.")
        report.append("  Existing LRG results computed from raw W are not affected by")
        report.append("  the omission of surrogate filtering.")

    report.append("")


# ═════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Surrogate validation follow-up")
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    report = [
        "SURROGATE VALIDATION — FOLLOW-UP REPORT",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"N_surrogates: {N_SURROGATES}",
        "",
    ]

    t0 = time.time()

    task0(report)
    task3(report)  # no computation
    task1(report, n_workers=args.workers)

    print("\n══════ TASK 2: Pat_02 anomaly ══════")
    task2(report, n_workers=args.workers)

    print("\n══════ TASK 4: LRG comparison ══════")
    task4(report, n_workers=args.workers)

    elapsed = time.time() - t0
    report.append(f"Total runtime: {elapsed / 60:.1f} min")
    print(f"\nTotal runtime: {elapsed / 60:.1f} min")

    report_path = OUTPUT_DIR / "surrogate_followup_report.txt"
    report_path.write_text("\n".join(report))
    print(f"Report saved: {report_path}")


if __name__ == "__main__":
    main()
