#!/usr/bin/env python3
"""Compute MSC-based functional connectivity matrices.

This script computes MSC matrices for all patients, phases, and bands.
Supports dense MSC and multiple sparsification methods.

Usage:
    # Dense MSC (no validation)
    python scripts/py/compute_msc_matrices.py

    # Validated MSC with surrogates (soft)
    python scripts/py/compute_msc_matrices.py --sparsify soft --n-surrogates 200

    # FDR-corrected thresholding
    python scripts/py/compute_msc_matrices.py --sparsify fdr --n-surrogates 200

    # Disparity filter (no surrogates needed)
    python scripts/py/compute_msc_matrices.py --sparsify disparity

    # Hybrid (surrogate excess + disparity)
    python scripts/py/compute_msc_matrices.py --sparsify hybrid --n-surrogates 200

    # Specific patients
    python scripts/py/compute_msc_matrices.py --patients Pat_02 Pat_03

    # Custom parameters
    python scripts/py/compute_msc_matrices.py --nperseg 512 --verbose
"""

import argparse
import gc
import os
from pathlib import Path
from lrg_eegfc import compute_msc_for_patient
from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS, VALID_SPARSIFY_METHODS
from lrg_eegfc.config.paths import MSC_CACHE, MSC_DEV_CACHE


def main():
    parser = argparse.ArgumentParser(
        description="Compute MSC-based functional connectivity matrices"
    )

    # Patient selection
    parser.add_argument(
        "--patients",
        nargs="+",
        default=["Pat_02", "Pat_03", "Pat_05", "Pat_08"],
        help="Patient IDs to process (default: Pat_02 Pat_03 Pat_05 Pat_08)"
    )

    # MSC parameters
    parser.add_argument(
        "--sparsify",
        choices=list(VALID_SPARSIFY_METHODS),
        default="none",
        help="Sparsification method (default: none for dense MSC)"
    )
    parser.add_argument(
        "--n-surrogates",
        type=int,
        default=0,
        help="Number of surrogates for validation (default: 0)"
    )
    parser.add_argument(
        "--fdr-q",
        type=float,
        default=0.05,
        help="FDR q-value threshold (only if sparsify=fdr, default: 0.05)"
    )
    parser.add_argument(
        "--disparity-alpha",
        type=float,
        default=0.05,
        help="Disparity filter significance level (only if sparsify=disparity or hybrid, default: 0.05)"
    )
    parser.add_argument(
        "--ecm-alpha",
        type=float,
        default=0.05,
        help="ECM z-test significance level (only if sparsify=ecm, default: 0.05)"
    )
    parser.add_argument(
        "--ecm-n-ensemble",
        type=int,
        default=100,
        help="ECM ensemble size for null model sampling (only if sparsify=ecm, default: 100)"
    )
    parser.add_argument(
        "--ecm-weight-scale",
        type=int,
        default=1000,
        help="ECM weight scaling factor for float→int conversion (only if sparsify=ecm, default: 1000)"
    )
    parser.add_argument(
        "--nperseg",
        type=int,
        default=4096,
        help="Window length for Welch's method (default: 4096)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=64,
        help="Number of Welch segments to process per batch (default: 64)"
    )
    parser.add_argument(
        "--sample-rate",
        type=float,
        default=2048.0,
        help="Sampling rate in Hz (default: 2048.0)"
    )
    parser.add_argument(
        "--filter-time",
        type=int,
        default=None,
        help="Limit to first N samples (dev-only convenience)"
    )
    parser.add_argument(
        "--n-workers",
        type=int,
        default=None,
        help="Number of parallel workers for surrogate computation (default: all CPUs)"
    )

    # Subset selection
    parser.add_argument(
        "--bands",
        nargs="+",
        default=None,
        help="Band names to process (default: all bands)"
    )
    parser.add_argument(
        "--phases",
        nargs="+",
        default=None,
        help="Phase names to process (default: all phases)"
    )

    # Cache control
    parser.add_argument(
        "--cache-root",
        type=Path,
        default=MSC_CACHE,
        help="Cache directory (default: data/cache/msc)"
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Recompute even if cached"
    )

    # Output control
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed progress"
    )

    args = parser.parse_args()

    if args.filter_time is not None and args.filter_time > 0 and args.cache_root == MSC_CACHE:
        args.cache_root = MSC_DEV_CACHE
        print(f"Using dev cache root for filter_time: {args.cache_root}")

    # Validate arguments: auto-set surrogates for methods that need them
    _NEEDS_SURROGATES = {"soft", "fdr", "hybrid"}
    if args.sparsify in _NEEDS_SURROGATES and args.n_surrogates == 0:
        print(f"WARNING: sparsify='{args.sparsify}' but n_surrogates=0. Setting n_surrogates=200.")
        args.n_surrogates = 200

    # Print configuration
    print("=" * 70)
    print("MSC Functional Connectivity Matrix Computation")
    print("=" * 70)
    print(f"Patients: {', '.join(args.patients)}")
    selected_bands = args.bands or list(BRAIN_BANDS.keys())
    selected_phases = args.phases or list(PHASE_LABELS)
    print(f"Bands: {', '.join(selected_bands)}")
    print(f"Phases: {', '.join(selected_phases)}")
    print(f"Sparsify: {args.sparsify}")
    if args.sparsify in _NEEDS_SURROGATES:
        print(f"N surrogates: {args.n_surrogates}")
        n_workers = args.n_workers if args.n_workers else os.cpu_count()
        print(f"N workers: {n_workers}")
    if args.sparsify == "fdr":
        print(f"FDR q: {args.fdr_q}")
    if args.sparsify in ("disparity", "hybrid"):
        print(f"Disparity alpha: {args.disparity_alpha}")
    if args.sparsify == "ecm":
        print(f"ECM alpha: {args.ecm_alpha}")
        print(f"ECM ensemble: {args.ecm_n_ensemble}")
        print(f"ECM weight scale: {args.ecm_weight_scale}")
    print(f"nperseg: {args.nperseg}")
    print(f"batch_size: {args.batch_size}")
    if args.filter_time:
        print(f"Filter time: {args.filter_time}")
    print(f"Cache root: {args.cache_root}")
    expected = len(args.patients) * len(selected_bands) * len(selected_phases)
    print(f"Expected matrices: {expected}")
    print("=" * 70)

    total_computed = 0
    total_failed = 0

    show_totals = len(args.patients) > 1

    for patient in args.patients:
        patient_cache = args.cache_root / patient
        print(f"\nProcessing {patient}... ({patient_cache}/)")

        results = compute_msc_for_patient(
            patient,
            bands=selected_bands,
            phases=selected_phases,
            verbose=args.verbose,
            sparsify=args.sparsify,
            n_surrogates=args.n_surrogates,
            nperseg=args.nperseg,
            batch_size=args.batch_size,
            n_workers=args.n_workers,
            sample_rate=args.sample_rate,
            filter_time=args.filter_time,
            cache_root=args.cache_root,
            overwrite_cache=args.overwrite,
            return_results=False,
            fdr_q=args.fdr_q,
            disparity_alpha=args.disparity_alpha,
            ecm_alpha=args.ecm_alpha,
            ecm_n_ensemble=args.ecm_n_ensemble,
            ecm_weight_scale=args.ecm_weight_scale,
        )

        # Count successes (True or MSCResult object means success, None means failure)
        n_computed = sum(
            1 for band in results
            for phase in results[band]
            if results[band][phase] is not None
        )
        n_failed = len(selected_bands) * len(selected_phases) - n_computed

        total_computed += n_computed
        total_failed += n_failed

        print(f"  Computed: {n_computed}/{len(selected_bands) * len(selected_phases)}")
        if n_failed > 0:
            print(f"  Failed: {n_failed}")

        # Free memory between patients
        del results
        gc.collect()

    print("\n" + "=" * 70)
    print("MSC matrix computation complete!")
    print("=" * 70)
    if show_totals:
        print(f"Total computed: {total_computed}/{expected}")
        if total_failed > 0:
            print(f"Total failed: {total_failed}")
    elif total_failed > 0:
        print(f"Failed: {total_failed}")

    # Print cache information
    if args.sparsify == "none":
        filename_pattern = f"{{band}}_{{phase}}_msc_sparsify-none_nperseg-{args.nperseg}.npy"
    elif args.sparsify == "soft":
        filename_pattern = (
            f"{{band}}_{{phase}}_msc_sparsify-soft_"
            f"nsurr-{args.n_surrogates}_nperseg-{args.nperseg}.npy"
        )
    elif args.sparsify == "fdr":
        filename_pattern = (
            f"{{band}}_{{phase}}_msc_sparsify-fdr_"
            f"nsurr-{args.n_surrogates}_q-{args.fdr_q}_nperseg-{args.nperseg}.npy"
        )
    elif args.sparsify == "disparity":
        filename_pattern = (
            f"{{band}}_{{phase}}_msc_sparsify-disparity_"
            f"alpha-{args.disparity_alpha}_nperseg-{args.nperseg}.npy"
        )
    elif args.sparsify == "hybrid":
        filename_pattern = (
            f"{{band}}_{{phase}}_msc_sparsify-hybrid_"
            f"nsurr-{args.n_surrogates}_alpha-{args.disparity_alpha}_nperseg-{args.nperseg}.npy"
        )
    elif args.sparsify == "ecm":
        filename_pattern = (
            f"{{band}}_{{phase}}_msc_sparsify-ecm_"
            f"alpha-{args.ecm_alpha}_nens-{args.ecm_n_ensemble}_wscale-{args.ecm_weight_scale}_nperseg-{args.nperseg}.npy"
        )
    if args.filter_time is not None and args.filter_time > 0:
        filename_pattern = filename_pattern.replace(".npy", f"_ftime-{args.filter_time}.npy")

    print("\nCache")
    print(f"  root:        {args.cache_root}/")
    print(f"  per-patient: {args.cache_root}/{{patient}}/")
    print(f"  filename:    {filename_pattern}")


if __name__ == "__main__":
    main()
