#!/usr/bin/env python3
"""Compute time-window FC matrices (corr/MSC) and cache results."""

from __future__ import annotations

import argparse
import csv
import json
import os
import time
from pathlib import Path
from typing import Dict, Optional

os.environ.setdefault("JOBLIB_MULTIPROCESSING", "0")

import networkx as nx
import numpy as np
from scipy.spatial.distance import squareform

from lrgsglib.core import (
    compute_laplacian_properties,
    compute_normalized_linkage,
    compute_optimal_threshold,
    entropy,
    extract_ultrametric_matrix,
    get_giant_component,
)
from lrgsglib.utils.basic.signals import bandpass_sos

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS, list_patients
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.fc.corr import build_corr_network
from lrg_eegfc.utils.fc.msc import coherence_fc_pipeline
from lrg_eegfc.utils.io import load_patient_dataset_robust
from lrg_eegfc.workflow.time_windows import (
    DEFAULT_WINDOW_CYCLES,
    DEFAULT_WINDOW_MIN_SEC,
    DEFAULT_WINDOW_OVERLAP,
    compute_window_params,
    generate_window_indices,
    get_lrg_window_cache_dir,
    get_window_cache_dir,
    suggest_window_sec,
)


def _compute_lrg(
    adjacency_matrix: np.ndarray,
    *,
    entropy_steps: int,
    entropy_t1: float,
    entropy_t2: float,
    verbose: bool,
) -> Dict[str, np.ndarray]:
    graph = nx.from_numpy_array(adjacency_matrix)
    giant = get_giant_component(graph)
    if verbose:
        print(f"    LRG giant component: {giant.number_of_nodes()}/{graph.number_of_nodes()} nodes")

    sm1, spec, *_, tau = entropy(giant, steps=entropy_steps, t1=entropy_t1, t2=entropy_t2)
    _, _, _, Trho, _ = compute_laplacian_properties(giant)
    dists = squareform(Trho)
    linkage_matrix, labels, _ = compute_normalized_linkage(dists, giant)
    threshold, *_ = compute_optimal_threshold(linkage_matrix)
    ultrametric_square = extract_ultrametric_matrix(linkage_matrix, giant.number_of_nodes())
    ultrametric_matrix = squareform(ultrametric_square)

    return {
        "ultrametric_matrix": ultrametric_matrix,
        "linkage_matrix": linkage_matrix,
        "entropy_tau": tau,
        "entropy_1_minus_S": sm1,
        "entropy_C": spec,
        "optimal_threshold": float(threshold),
        "n_nodes": int(giant.number_of_nodes()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute time-window FC matrices (corr/MSC) and cache results."
    )

    parser.add_argument("--patients", nargs="+", default=None, help="Patient IDs to process")
    parser.add_argument("--phases", nargs="+", default=None, help="Phase names (default: all)")
    parser.add_argument("--bands", nargs="+", default=None, help="Band names (default: all)")
    parser.add_argument(
        "--dataset-root",
        type=Path,
        default=SEEG_DATAPATH,
        help="Root directory containing patient data",
    )

    parser.add_argument(
        "--fc-method",
        choices=["corr", "msc"],
        default="corr",
        help="FC method to compute per window",
    )

    parser.add_argument("--window-sec", type=float, default=None, help="Window length in seconds")
    parser.add_argument(
        "--overlap",
        type=float,
        default=DEFAULT_WINDOW_OVERLAP,
        help="Fractional overlap between windows (default: 0.25)",
    )
    parser.add_argument(
        "--max-windows",
        type=int,
        default=None,
        help="Limit number of windows per run (dev-only convenience)",
    )
    parser.add_argument(
        "--min-window-sec",
        type=float,
        default=DEFAULT_WINDOW_MIN_SEC,
        help="Minimum window length when auto-selecting (default: 10)",
    )
    parser.add_argument(
        "--cycles",
        type=float,
        default=DEFAULT_WINDOW_CYCLES,
        help="Minimum cycles per lowest-band frequency when auto-selecting window length",
    )
    parser.add_argument(
        "--filter-time",
        type=int,
        default=None,
        help="Limit to first N samples (dev-only convenience)",
    )

    parser.add_argument(
        "--filter-order",
        type=int,
        default=4,
        help="Bandpass filter order (corr only)",
    )
    parser.add_argument(
        "--filter-type",
        choices=["abs", "pos", "neg", "none"],
        default="abs",
        help="Correlation filter type (corr only)",
    )
    parser.add_argument(
        "--keep-diagonal",
        action="store_true",
        help="Keep diagonal in correlation matrices (corr only)",
    )

    parser.add_argument(
        "--sparsify",
        choices=["none", "soft", "fdr", "disparity", "hybrid", "ecm"],
        default="none",
        help="MSC sparsification method (msc only)",
    )
    parser.add_argument(
        "--n-surrogates",
        type=int,
        default=0,
        help="Number of surrogates (msc only, soft sparsify)",
    )
    parser.add_argument(
        "--nperseg",
        type=int,
        default=1024,
        help="Welch segment length (msc only)",
    )
    parser.add_argument(
        "--noverlap",
        type=int,
        default=None,
        help="Welch overlap length (msc only)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=64,
        help="Welch batch size (msc only)",
    )

    parser.add_argument(
        "--lrg",
        action="store_true",
        help="Compute LRG outputs per window",
    )
    parser.add_argument("--entropy-steps", type=int, default=400, help="LRG entropy steps")
    parser.add_argument("--entropy-t1", type=float, default=-3.0, help="LRG entropy t1")
    parser.add_argument("--entropy-t2", type=float, default=5.0, help="LRG entropy t2")

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Recompute even if cached",
    )
    parser.add_argument("--verbose", action="store_true", help="Verbose logging")

    args = parser.parse_args()

    patients = args.patients
    if patients is None:
        patients = list_patients(args.dataset_root)
        if not patients:
            raise SystemExit("No patients found under data/stereoeeg_patients")

    phases = args.phases or list(PHASE_LABELS)
    bands = args.bands or list(BRAIN_BANDS.keys())

    if args.window_sec is None and len(bands) > 1:
        print("Note: auto window length will be computed per band. Use --window-sec to keep windows aligned.")

    print("=" * 70)
    print("Time-window FC computation")
    print("=" * 70)
    print(f"Patients: {', '.join(patients)}")
    print(f"Phases: {', '.join(phases)}")
    print(f"Bands: {', '.join(bands)}")
    print(f"FC method: {args.fc_method}")
    print(f"Overlap: {args.overlap}")
    if args.window_sec is not None:
        print(f"Window length: {args.window_sec:.3f} sec")
    else:
        print(f"Window length: auto (min {args.min_window_sec}s, {args.cycles} cycles)")
    if args.filter_time:
        print(f"Filter time: {args.filter_time}")
    if args.max_windows:
        print(f"Max windows: {args.max_windows}")
    print("=" * 70)

    total_computed = 0
    total_skipped = 0
    total_failed = 0

    zero_diagonal = not args.keep_diagonal

    for patient in patients:
        dataset = load_patient_dataset_robust(patient, args.dataset_root, phases=phases)
        for phase in phases:
            if phase not in dataset:
                print(f"Skipping {patient} {phase}: not available")
                continue

            recording = dataset[phase]
            fs = float(recording.parameters.get("fs", 2048.0))
            data = recording.timeseries
            if args.filter_time is not None and args.filter_time > 0:
                data = data[:, : args.filter_time]

            n_channels, n_samples = data.shape
            if n_samples < 2:
                print(f"Skipping {patient} {phase}: not enough samples")
                continue

            for band in bands:
                low, high = BRAIN_BANDS[band]
                window_sec = args.window_sec
                if window_sec is None:
                    window_sec = suggest_window_sec(
                        low, min_window_sec=args.min_window_sec, cycles=args.cycles
                    )

                try:
                    window_len, step = compute_window_params(
                        n_samples, fs, window_sec, args.overlap
                    )
                except ValueError as exc:
                    print(f"Skipping {patient} {phase} {band}: {exc}")
                    continue

                ratio = window_len / max(1, n_channels)
                if ratio < 10:
                    print(
                        f"Warning: {patient} {phase} {band} window_len/channels={ratio:.2f} "
                        "may be too small for stable FC estimates."
                    )

                starts, total_windows = generate_window_indices(
                    n_samples, window_len, step, max_windows=args.max_windows
                )

                run_dir = get_window_cache_dir(
                    patient=patient,
                    phase=phase,
                    band=band,
                    fc_method=args.fc_method,
                    window_sec=window_sec,
                    overlap=args.overlap,
                    filter_time=args.filter_time,
                    filter_type=args.filter_type,
                    zero_diagonal=zero_diagonal,
                    filter_order=args.filter_order,
                    sparsify=args.sparsify,
                    n_surrogates=args.n_surrogates,
                    nperseg=args.nperseg,
                    noverlap=args.noverlap,
                )

                lrg_dir = None
                if args.lrg:
                    lrg_dir = get_lrg_window_cache_dir(
                        patient=patient,
                        phase=phase,
                        band=band,
                        fc_method=args.fc_method,
                        window_sec=window_sec,
                        overlap=args.overlap,
                        filter_time=args.filter_time,
                        filter_type=args.filter_type,
                        zero_diagonal=zero_diagonal,
                        filter_order=args.filter_order,
                        sparsify=args.sparsify,
                        n_surrogates=args.n_surrogates,
                        nperseg=args.nperseg,
                        noverlap=args.noverlap,
                    )

                window_meta_path = run_dir / "windows.csv"
                run_meta_path = run_dir / "run_meta.json"
                window_rows = []

                if args.verbose:
                    print(
                        f"\n{patient} {phase} {band}: fs={fs} n={n_channels} samples={n_samples} "
                        f"window={window_len} step={step} windows={len(starts)}/{total_windows}"
                    )

                for idx, start in enumerate(starts):
                    stop = start + window_len
                    window_id = f"win-{idx:04d}"
                    matrix_path = run_dir / f"{window_id}.npy"
                    lrg_path = lrg_dir / f"{window_id}_lrg.npz" if lrg_dir else None

                    if matrix_path.exists() and not args.overwrite:
                        total_skipped += 1
                        if args.verbose:
                            print(f"  [{idx + 1}/{len(starts)}] cache hit: {matrix_path.name}")
                        window_rows.append((idx, start, stop, start / fs, stop / fs))
                        continue

                    t0 = time.perf_counter()
                    window = data[:, start:stop]

                    try:
                        if args.fc_method == "corr":
                            filtered = bandpass_sos(window, low, high, fs, args.filter_order)
                            matrix = build_corr_network(
                                filtered,
                                filter_type=args.filter_type,
                                zero_diagonal=zero_diagonal,
                            )
                        else:
                            if window.shape[1] < args.nperseg:
                                raise ValueError(
                                    f"Window length {window.shape[1]} < nperseg {args.nperseg}"
                                )
                            matrices = coherence_fc_pipeline(
                                window,
                                fs=fs,
                                bands={band: (low, high)},
                                sparsify=args.sparsify,
                                n_surrogates=args.n_surrogates if args.sparsify == "soft" else 0,
                                nperseg=args.nperseg,
                                noverlap=args.noverlap,
                                batch_size=args.batch_size,
                                zero_diagonal=True,
                                verbose=False,
                            )
                            matrix = matrices[band]

                        np.save(matrix_path, matrix)
                        total_computed += 1

                        if args.lrg and lrg_path is not None:
                            if lrg_path.exists() and not args.overwrite:
                                if args.verbose:
                                    print(f"    LRG cache hit: {lrg_path.name}")
                            else:
                                lrg_payload = _compute_lrg(
                                    matrix,
                                    entropy_steps=args.entropy_steps,
                                    entropy_t1=args.entropy_t1,
                                    entropy_t2=args.entropy_t2,
                                    verbose=args.verbose,
                                )
                                np.savez_compressed(
                                    lrg_path,
                                    **lrg_payload,
                                    patient=patient,
                                    phase=phase,
                                    band=band,
                                    fc_method=args.fc_method,
                                    window_index=idx,
                                    window_start=int(start),
                                    window_stop=int(stop),
                                )

                        elapsed = time.perf_counter() - t0
                        if args.verbose:
                            print(
                                f"  [{idx + 1}/{len(starts)}] saved {matrix_path.name} "
                                f"({elapsed:.2f}s)"
                            )

                        window_rows.append((idx, start, stop, start / fs, stop / fs))
                    except Exception as exc:
                        total_failed += 1
                        print(f"  [{idx + 1}/{len(starts)}] failed: {exc}")

                run_meta = {
                    "patient": patient,
                    "phase": phase,
                    "band": band,
                    "fc_method": args.fc_method,
                    "fs": fs,
                    "n_channels": n_channels,
                    "n_samples": n_samples,
                    "window_sec": window_sec,
                    "window_len": window_len,
                    "overlap": args.overlap,
                    "step": step,
                    "total_windows": total_windows,
                    "saved_windows": len(window_rows),
                    "filter_time": args.filter_time,
                    "params": {
                        "filter_type": args.filter_type,
                        "zero_diagonal": zero_diagonal,
                        "filter_order": args.filter_order,
                        "sparsify": args.sparsify,
                        "n_surrogates": args.n_surrogates,
                        "nperseg": args.nperseg,
                        "noverlap": args.noverlap,
                        "batch_size": args.batch_size,
                    },
                    "lrg": {
                        "enabled": args.lrg,
                        "entropy_steps": args.entropy_steps,
                        "entropy_t1": args.entropy_t1,
                        "entropy_t2": args.entropy_t2,
                    },
                }

                with open(run_meta_path, "w", encoding="utf-8") as handle:
                    json.dump(run_meta, handle, indent=2)

                with open(window_meta_path, "w", newline="", encoding="utf-8") as handle:
                    writer = csv.writer(handle)
                    writer.writerow(["window_index", "start", "stop", "start_sec", "stop_sec"])
                    writer.writerows(window_rows)

    print("\n" + "=" * 70)
    print("Time-window computation complete")
    print("=" * 70)
    print(f"Computed: {total_computed}")
    print(f"Skipped:  {total_skipped}")
    if total_failed:
        print(f"Failed:   {total_failed}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
