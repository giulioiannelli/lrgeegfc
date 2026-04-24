#!/usr/bin/env python3
"""Surrogate validation follow-up: Tasks 1, 2, and 4.

Task 1: Pat_06 surrogate validation (rest_pre, 6 bands)
Task 2: Pat_02 anomaly investigation (recording durations, Welch segments)
Task 4: LRG output comparison raw W vs W_soft (Pat_02, 4 phases, 6 bands)

Outputs:
  data/surrogate_validation/surrogate_summary_table.csv  (Pat_06 rows appended)
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
from lrg_eegfc.workflow.lrg import load_lrg_result, compute_lrg_analysis
from lrg_eegfc.workflow.msc import load_msc_matrix

from lrg_eegfc.config.paths import SURROGATE_VALIDATION_CACHE
OUTPUT_DIR = SURROGATE_VALIDATION_CACHE
N_SURROGATES = 100
PHASE = "rest_pre"
BAND_NAMES = BRAIN_BANDS_NAMES
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
PHASES = ["rest_pre", "rest_post", "task_learn", "task_test"]


def _upper_tri(W):
    idx = np.triu_indices_from(W, k=1)
    return W[idx]


def _load_timeseries(patient, phase):
    recordings = load_patient_dataset_robust(patient, phases=[phase])
    rec = recordings[phase]
    fs = float(rec.parameters.get("fs", 2048.0))
    return rec.timeseries, fs


def _compute_empirical_msc(X, fs, nperseg):
    freqs, Coh = compute_msc_welch(X, fs, nperseg=nperseg, noverlap=nperseg // 2)
    W_bands = band_average_msc(Coh, freqs, BRAIN_BANDS)
    del Coh
    gc.collect()
    return W_bands


def _compute_null(X, fs, nperseg, n_surrogates=N_SURROGATES, n_workers=None):
    return surrogate_msc_null(
        X, fs, BRAIN_BANDS, n_surrogates,
        nperseg=nperseg, noverlap=nperseg // 2,
        n_workers=n_workers,
    )


# ═════════════════════════════════════════════════════════════════════════════
#  TASK 1: Pat_06 surrogate validation
# ═════════════════════════════════════════════════════════════════════════════

def task1_pat06(report_lines, n_workers=None):
    """Surrogate validation for Pat_06, rest_pre, all 6 bands."""
    print("\n══════ TASK 1: Pat_06 surrogate validation ══════")
    report_lines.append("=" * 70)
    report_lines.append("TASK 1: Pat_06 surrogate validation (rest_pre)")
    report_lines.append("=" * 70)

    X, fs = _load_timeseries("Pat_06", PHASE)
    nperseg = nperseg_for_fs(fs)
    n_channels, n_samples = X.shape
    duration_min = n_samples / fs / 60

    report_lines.append(f"  Pat_06: fs={fs} Hz, nperseg={nperseg}, shape={X.shape}")
    report_lines.append(f"  Duration: {duration_min:.2f} min, N_surrogates={N_SURROGATES}")
    report_lines.append(f"  Welch: hann, noverlap={nperseg // 2}, df={fs / nperseg:.2f} Hz")
    report_lines.append("")
    print(f"  fs={fs}, nperseg={nperseg}, shape={X.shape}, duration={duration_min:.1f} min")

    print("  Computing empirical MSC...")
    W_bands = _compute_empirical_msc(X, fs, nperseg)
    print(f"  Computing {N_SURROGATES} surrogates...")
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

        line = (
            f"  {band:>12s}: sig_p05={frac_p05:.3f}  sig_p01={frac_p01:.3f}  "
            f"rho={rho:.4f}  -> ({cat})"
        )
        print(f"  {line}")
        report_lines.append(line)

    report_lines.append("")
    del W_bands, W_null
    gc.collect()

    return pd.DataFrame(rows)


# ═════════════════════════════════════════════════════════════════════════════
#  TASK 2: Pat_02 anomaly investigation
# ═════════════════════════════════════════════════════════════════════════════

def task2_pat02_anomaly(report_lines, n_workers=None):
    """Compare recording durations and MSC statistics across patients."""
    print("\n══════ TASK 2: Pat_02 anomaly investigation ══════")
    report_lines.append("=" * 70)
    report_lines.append("TASK 2: Pat_02 anomaly investigation")
    report_lines.append("=" * 70)

    patient_info = {}

    for patient in ALL_PATIENTS:
        print(f"  Loading {patient}...")
        X, fs = _load_timeseries(patient, PHASE)
        nperseg = nperseg_for_fs(fs)
        n_channels, n_samples = X.shape
        duration_sec = n_samples / fs
        duration_min = duration_sec / 60

        # Number of Welch segments: M = floor((L - noverlap) / (nperseg - noverlap))
        noverlap = nperseg // 2
        step = nperseg - noverlap
        M = (n_samples - noverlap) // step

        # Compute empirical MSC
        print(f"    Computing empirical MSC (fs={fs}, nperseg={nperseg})...")
        W_bands = _compute_empirical_msc(X, fs, nperseg)

        mean_msc = {}
        for band in BAND_NAMES:
            W = W_bands[band].copy()
            np.fill_diagonal(W, 0.0)
            mean_msc[band] = float(_upper_tri(W).mean())

        patient_info[patient] = {
            "fs": fs,
            "nperseg": nperseg,
            "n_channels": n_channels,
            "n_samples": n_samples,
            "duration_sec": duration_sec,
            "duration_min": duration_min,
            "M_welch": M,
            "mean_msc": mean_msc,
        }

        del X, W_bands
        gc.collect()

    # Report: recording parameters
    report_lines.append("")
    report_lines.append("  Recording parameters (rest_pre):")
    header = f"  {'Patient':<10s} {'fs(Hz)':>8s} {'nperseg':>8s} {'N_ch':>6s} {'N_samples':>12s} {'Duration(min)':>14s} {'M_welch':>8s}"
    report_lines.append(header)
    report_lines.append("  " + "-" * (len(header) - 2))
    print(header)

    for patient in ALL_PATIENTS:
        info = patient_info[patient]
        line = (
            f"  {patient:<10s} {info['fs']:>8.0f} {info['nperseg']:>8d} "
            f"{info['n_channels']:>6d} {info['n_samples']:>12d} "
            f"{info['duration_min']:>14.2f} {info['M_welch']:>8d}"
        )
        print(line)
        report_lines.append(line)

    # Report: mean MSC per band
    report_lines.append("")
    report_lines.append("  Mean raw MSC per band (upper triangle, diagonal excluded):")
    band_header = f"  {'Patient':<10s}" + "".join(f" {b:>12s}" for b in BAND_NAMES)
    report_lines.append(band_header)
    report_lines.append("  " + "-" * (len(band_header) - 2))
    print("\n" + band_header)

    for patient in ALL_PATIENTS:
        vals = patient_info[patient]["mean_msc"]
        line = f"  {patient:<10s}" + "".join(f" {vals[b]:>12.4f}" for b in BAND_NAMES)
        print(line)
        report_lines.append(line)

    # Surrogate analysis for Pat_02 vs others: compute surrogate means
    report_lines.append("")
    report_lines.append("  Surrogate null analysis (mean surrogate MSC per band):")

    surr_means = {}
    for patient in ALL_PATIENTS:
        print(f"  Computing surrogates for {patient}...")
        X, fs = _load_timeseries(patient, PHASE)
        nperseg = nperseg_for_fs(fs)
        W_null = _compute_null(X, fs, nperseg, n_surrogates=N_SURROGATES, n_workers=n_workers)

        surr_means[patient] = {}
        for band in BAND_NAMES:
            Wn = W_null[band]
            triu_i, triu_j = np.triu_indices(Wn.shape[1], k=1)
            surr_vals = Wn[:, triu_i, triu_j]
            surr_means[patient][band] = {
                "mean": float(surr_vals.mean()),
                "std": float(surr_vals.std()),
                "p95": float(np.percentile(surr_vals, 95)),
                "p99": float(np.percentile(surr_vals, 99)),
            }

        del X, W_null
        gc.collect()

    surr_header = f"  {'Patient':<10s}" + "".join(f" {b:>12s}" for b in BAND_NAMES)
    report_lines.append(surr_header)
    report_lines.append("  " + "-" * (len(surr_header) - 2))
    print("\n  Surrogate means:")

    for patient in ALL_PATIENTS:
        vals = surr_means[patient]
        line = f"  {patient:<10s}" + "".join(
            f" {vals[b]['mean']:>12.5f}" for b in BAND_NAMES
        )
        print(line)
        report_lines.append(line)

    # Ratio: empirical mean / surrogate mean
    report_lines.append("")
    report_lines.append("  Ratio empirical mean / surrogate mean:")
    report_lines.append(surr_header)
    report_lines.append("  " + "-" * (len(surr_header) - 2))
    print("\n  Empirical/surrogate ratio:")

    for patient in ALL_PATIENTS:
        emp = patient_info[patient]["mean_msc"]
        surr = surr_means[patient]
        line = f"  {patient:<10s}" + "".join(
            f" {emp[b] / surr[b]['mean']:>12.1f}" for b in BAND_NAMES
        )
        print(line)
        report_lines.append(line)

    # 95th percentile of surrogate null for low_gamma (Pat_02's worst band)
    report_lines.append("")
    report_lines.append("  95th percentile of surrogate null (low_gamma — Pat_02's lowest sig band):")
    for patient in ALL_PATIENTS:
        s = surr_means[patient]["low_gamma"]
        emp_lg = patient_info[patient]["mean_msc"]["low_gamma"]
        line = (
            f"  {patient:<10s}: surr_p95={s['p95']:.5f}  surr_p99={s['p99']:.5f}  "
            f"emp_mean={emp_lg:.4f}  M_welch={patient_info[patient]['M_welch']}"
        )
        print(line)
        report_lines.append(line)

    # Hypothesis conclusion
    report_lines.append("")
    pat02_M = patient_info["Pat_02"]["M_welch"]
    other_Ms = [patient_info[p]["M_welch"] for p in ALL_PATIENTS if p != "Pat_02"]
    avg_other_M = np.mean(other_Ms)
    report_lines.append(
        f"  Hypothesis check: Pat_02 M_welch={pat02_M}, "
        f"other patients avg M_welch={avg_other_M:.0f}"
    )
    if pat02_M < avg_other_M * 0.8:
        report_lines.append(
            "  CONFIRMED: Pat_02 has substantially fewer Welch segments, producing "
            "less stable MSC estimates and wider surrogate null distributions."
        )
    else:
        report_lines.append(
            "  NOT CONFIRMED by recording duration alone. The difference in M_welch "
            "is modest. Other factors (signal quality, electrode placement) may contribute."
        )

    report_lines.append("")
    return patient_info, surr_means


# ═════════════════════════════════════════════════════════════════════════════
#  TASK 4: LRG comparison raw W vs W_soft
# ═════════════════════════════════════════════════════════════════════════════

def task4_lrg_comparison(report_lines):
    """Compare LRG outputs from raw W vs W_soft for Pat_02, all phases/bands."""
    import networkx as nx
    from lrgsglib.core import (
        compute_laplacian_properties,
        compute_normalized_linkage,
        compute_optimal_threshold,
        get_giant_component,
    )
    from scipy.cluster.hierarchy import cophenet

    print("\n══════ TASK 4: LRG comparison raw W vs W_soft ══════")
    report_lines.append("=" * 70)
    report_lines.append("TASK 4: LRG output comparison — raw W vs W_soft (Pat_02)")
    report_lines.append("=" * 70)

    patient = "Pat_02"
    nperseg = 4096

    def _run_lrg(W):
        """Run LRG on adjacency matrix, return linkage, threshold, n_nodes."""
        np.fill_diagonal(W, 0.0)
        G = nx.from_numpy_array(W)
        giant = get_giant_component(G)
        n_nodes = giant.number_of_nodes()
        _, _, _, Trho, _ = compute_laplacian_properties(giant)
        dists = squareform(Trho)
        linkage, labels, _ = compute_normalized_linkage(dists, giant)
        threshold, *_ = compute_optimal_threshold(linkage)
        return linkage, threshold, n_nodes

    results = []
    n_mismatch = 0
    n_low_coph = 0

    for phase in PHASES:
        report_lines.append(f"\n  Phase: {phase}")
        print(f"\n  Phase: {phase}")

        for band in BAND_NAMES:
            # Load raw W from MSC cache (sparsify=none)
            W_raw = load_msc_matrix(
                patient, phase, band,
                sparsify="none", nperseg=nperseg,
            )
            if W_raw is None:
                msg = f"    {band:>12s}: raw W not cached — SKIPPED"
                print(msg)
                report_lines.append(msg)
                continue

            # Load W_soft from MSC cache (sparsify=soft, nsurr=200)
            W_soft = load_msc_matrix(
                patient, phase, band,
                sparsify="soft", n_surrogates=200, nperseg=nperseg,
            )
            if W_soft is None:
                msg = f"    {band:>12s}: W_soft not cached — SKIPPED"
                print(msg)
                report_lines.append(msg)
                continue

            # Run LRG on both
            try:
                link_raw, thresh_raw, n_raw = _run_lrg(W_raw.copy())
                link_soft, thresh_soft, n_soft = _run_lrg(W_soft.copy())
            except Exception as e:
                msg = f"    {band:>12s}: LRG failed — {e}"
                print(msg)
                report_lines.append(msg)
                continue

            # Cophenetic correlation (only if same n_nodes)
            if n_raw == n_soft:
                d_raw = squareform(cophenet(link_raw))
                d_soft = squareform(cophenet(link_soft))
                coph_rho, _ = spearmanr(
                    d_raw[np.triu_indices_from(d_raw, k=1)],
                    d_soft[np.triu_indices_from(d_soft, k=1)],
                )
            else:
                coph_rho = np.nan

            # Load existing LRG result to get the cached n* (optimal threshold)
            lrg_cached = load_lrg_result(patient, phase, band, "msc")
            n_star_cached = lrg_cached.optimal_threshold if lrg_cached else thresh_raw

            # Check for meaningful disagreement
            flag = ""
            if coph_rho < 0.99 and not np.isnan(coph_rho):
                flag = " *** LOW COPHENETIC ***"
                n_low_coph += 1

            results.append({
                "phase": phase,
                "band": band,
                "coph_corr": coph_rho,
                "n_star_raw": thresh_raw,
                "n_star_soft": thresh_soft,
                "n_nodes_raw": n_raw,
                "n_nodes_soft": n_soft,
            })

            line = (
                f"    {band:>12s}: coph={coph_rho:.4f}  "
                f"n*_raw={thresh_raw:.4f}  n*_soft={thresh_soft:.4f}  "
                f"nodes_raw={n_raw}  nodes_soft={n_soft}{flag}"
            )
            print(line)
            report_lines.append(line)

    # Summary table
    report_lines.append("")
    report_lines.append("  Summary table (cophenetic_corr, n*_raw, n*_soft):")
    report_lines.append("")

    # Build pivot-style table
    header = f"  {'Phase':<12s}" + "".join(f" {b:>22s}" for b in BAND_NAMES)
    report_lines.append(header)
    report_lines.append("  " + "-" * (len(header) - 2))
    print("\n" + header)

    df_results = pd.DataFrame(results)
    for phase in PHASES:
        phase_data = df_results[df_results["phase"] == phase]
        cells = []
        for band in BAND_NAMES:
            row = phase_data[phase_data["band"] == band]
            if len(row) == 0:
                cells.append(f" {'—':>22s}")
            else:
                r = row.iloc[0]
                cell = f"{r['coph_corr']:.3f}/{r['n_star_raw']:.3f}/{r['n_star_soft']:.3f}"
                cells.append(f" {cell:>22s}")
        line = f"  {phase:<12s}" + "".join(cells)
        print(line)
        report_lines.append(line)

    # Overall assessment
    report_lines.append("")
    all_coph = [r["coph_corr"] for r in results if not np.isnan(r["coph_corr"])]
    if all_coph:
        min_coph = min(all_coph)
        mean_coph = np.mean(all_coph)
        report_lines.append(
            f"  Cophenetic correlations: min={min_coph:.4f}, mean={mean_coph:.4f}, "
            f"n_cells={len(all_coph)}"
        )
        if min_coph >= 0.99:
            report_lines.append(
                "  CONCLUSION: All cophenetic correlations >= 0.99. LRG outputs are "
                "effectively identical between raw W and W_soft. Existing results are "
                "unaffected by the choice to omit surrogate filtering."
            )
        elif n_low_coph <= 2:
            report_lines.append(
                f"  CONCLUSION: {n_low_coph} cell(s) with cophenetic < 0.99 (flagged above). "
                "Impact is minor. Existing results are largely unaffected."
            )
        else:
            report_lines.append(
                f"  WARNING: {n_low_coph} cells with cophenetic < 0.99. "
                "Some existing results may need recomputation."
            )

    # Check n* agreement
    n_star_agree = sum(
        1 for r in results
        if abs(r["n_star_raw"] - r["n_star_soft"]) < 0.05
    )
    n_star_close = sum(
        1 for r in results
        if abs(r["n_star_raw"] - r["n_star_soft"]) < 0.15
    )
    report_lines.append(
        f"  n* agreement: {n_star_agree}/{len(results)} exact (diff<0.05), "
        f"{n_star_close}/{len(results)} close (diff<0.15)"
    )
    report_lines.append("")

    return df_results


# ═════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Surrogate validation follow-up")
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    report_lines = [
        "SURROGATE VALIDATION FOLLOW-UP REPORT",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"N_surrogates: {N_SURROGATES} (validation), 200 (production pipeline)",
        f"Phase: {PHASE} (Tasks 1-2), all phases (Task 4)",
        "",
    ]

    t0 = time.time()

    # ── Task 0 ──────────────────────────────────────────────────────────────
    report_lines.append("=" * 70)
    report_lines.append("TASK 0: Production N_sur value")
    report_lines.append("=" * 70)
    report_lines.append("  Production pipeline: DEFAULT_N_SURROGATES = 200")
    report_lines.append("    Source: src/lrg_eegfc/config/const.py:135")
    report_lines.append("    CLI auto-set: src/lrg_eegfc/cli/compute.py:143")
    report_lines.append("  Validation script: N_SURROGATES = 100")
    report_lines.append("    Source: scripts/surrogate_validation.py:45")
    report_lines.append("  Manuscript text should say 200, not 100.")
    report_lines.append("")

    # ── Task 1: Pat_06 ─────────────────────────────────────────────────────
    df_pat06 = task1_pat06(report_lines, n_workers=args.workers)

    # Append to existing CSV
    csv_path = OUTPUT_DIR / "surrogate_summary_table.csv"
    df_existing = pd.read_csv(csv_path)
    # Remove any existing Pat_06 rows to avoid duplicates
    df_existing = df_existing[df_existing["patient"] != "Pat_06"]
    df_updated = pd.concat([df_existing, df_pat06], ignore_index=True)
    df_updated.to_csv(csv_path, index=False)
    report_lines.append(f"  Updated CSV: {csv_path} ({len(df_updated)} rows)")
    print(f"\n  Updated CSV with Pat_06: {csv_path}")

    # ── Task 2: Pat_02 anomaly ─────────────────────────────────────────────
    patient_info, surr_means = task2_pat02_anomaly(report_lines, n_workers=args.workers)

    # ── Task 3 ──────────────────────────────────────────────────────────────
    report_lines.append("=" * 70)
    report_lines.append("TASK 3: Pat_03 nperseg confirmation")
    report_lines.append("=" * 70)
    report_lines.append("  Pat_03 sampling rate: 1024 Hz")
    report_lines.append("  nperseg_for_fs(1024) = int(1024 * 2.0) = 2048")
    report_lines.append("  Confirmed: surrogate_validation.py:330 calls nperseg_for_fs(fs)")
    report_lines.append("  This gives 2-second Welch segments (df=0.50 Hz), matching all")
    report_lines.append("  other patients. The cached data under nperseg-2048 is correct.")
    report_lines.append("  high_gamma band (80-300 Hz) is within Pat_03 Nyquist (512 Hz).")
    report_lines.append("")

    # ── Task 4: LRG comparison ─────────────────────────────────────────────
    df_lrg = task4_lrg_comparison(report_lines)

    # ── Save report ────────────────────────────────────────────────────────
    elapsed = time.time() - t0
    report_lines.append(f"Total runtime: {elapsed / 60:.1f} min")
    print(f"\nTotal runtime: {elapsed / 60:.1f} min")

    report_path = OUTPUT_DIR / "surrogate_followup_report.txt"
    report_path.write_text("\n".join(report_lines))
    print(f"Report saved: {report_path}")


if __name__ == "__main__":
    main()
