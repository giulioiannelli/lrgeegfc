"""``lrg-eegfc compute`` subcommands."""

from __future__ import annotations

import gc
from pathlib import Path

import click

from ._common import (
    CliReporter,
    band_phase_options,
    cache_options,
    fc_method_option,
    filter_time_option,
    patient_options,
    resolve_bands,
    resolve_patients,
    resolve_phases,
    verbose_option,
)


@click.group()
def compute() -> None:
    """Compute functional connectivity matrices and analyses."""


# ---------------------------------------------------------------------------
# compute corr
# ---------------------------------------------------------------------------

@compute.command()
@patient_options()
@band_phase_options()
@cache_options("data/corr_cache")
@filter_time_option()
@click.option("--filter-type", type=click.Choice(["abs", "pos", "neg", "none"]),
              default="abs", show_default=True, help="Correlation filter type.")
@click.option("--zero-diagonal/--keep-diagonal", default=True, show_default=True)
@click.option("--filter-order", type=int, default=4, show_default=True)
@click.option("--sample-rate", type=float, default=2048.0, show_default=True)
@verbose_option()
@click.pass_context
def corr(ctx, patients, bands, band, phases, phase, cache_root, overwrite,
         filter_time, filter_type, zero_diagonal, filter_order, sample_rate,
         verbose):
    """Compute correlation-based FC matrices."""
    from lrg_eegfc.workflow.corr import compute_corr_for_patient

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    pat_list = resolve_patients(patients)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)

    rpt.header(
        "Correlation FC Matrix Computation",
        Patients=", ".join(pat_list),
        Bands=", ".join(band_list),
        Phases=", ".join(phase_list),
        filter_type=filter_type,
        zero_diagonal=zero_diagonal,
    )

    for patient in pat_list:
        rpt.info(f"\nProcessing {patient}...")
        results = compute_corr_for_patient(
            patient,
            bands=band_list,
            phases=phase_list,
            filter_type=filter_type,
            zero_diagonal=zero_diagonal,
            filter_order=filter_order,
            sample_rate=sample_rate,
            filter_time=filter_time,
            cache_root=cache_root,
            overwrite_cache=overwrite,
            verbose=rpt.verbose,
        )

        n_ok = sum(
            1 for b in results for p in results[b] if results[b][p] is not None
        )
        n_total = len(band_list) * len(phase_list)
        rpt.info(f"  {n_ok}/{n_total} computed")
        for _ in range(n_ok):
            rpt.success()
        for _ in range(n_total - n_ok):
            rpt.fail(f"{patient}: some combos failed")

        del results
        gc.collect()

    rpt.summary()


# ---------------------------------------------------------------------------
# compute msc
# ---------------------------------------------------------------------------

@compute.command()
@patient_options()
@band_phase_options()
@cache_options("data/msc_cache")
@filter_time_option()
@click.option("--sparsify",
              type=click.Choice(["none", "soft", "fdr", "disparity", "hybrid", "ecm"]),
              default="none", show_default=True, help="Sparsification method.")
@click.option("--n-surrogates", type=int, default=0, show_default=True,
              help="Number of surrogates (auto-set to 200 if needed).")
@click.option("--nperseg", type=int, default=4096, show_default=True,
              help="Welch window length.")
@click.option("--batch-size", type=int, default=64, show_default=True)
@click.option("--n-workers", type=int, default=None,
              help="Parallel workers for surrogates (default: all CPUs).")
@click.option("--sample-rate", type=float, default=2048.0, show_default=True)
@click.option("--fdr-q", type=float, default=0.05, show_default=True)
@click.option("--disparity-alpha", type=float, default=0.05, show_default=True)
@click.option("--ecm-alpha", type=float, default=0.05, show_default=True)
@click.option("--ecm-n-ensemble", type=int, default=100, show_default=True)
@click.option("--ecm-weight-scale", type=int, default=1000, show_default=True)
@verbose_option()
@click.pass_context
def msc(ctx, patients, bands, band, phases, phase, cache_root, overwrite,
        filter_time, sparsify, n_surrogates, nperseg, batch_size, n_workers,
        sample_rate, fdr_q, disparity_alpha, ecm_alpha, ecm_n_ensemble,
        ecm_weight_scale, verbose):
    """Compute MSC-based FC matrices."""
    from lrg_eegfc.workflow.msc import compute_msc_for_patient, DEFAULT_MSC_CACHE_ROOT

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    pat_list = resolve_patients(patients)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)

    # Auto-set surrogates for methods that need them
    needs_surr = {"soft", "fdr", "hybrid"}
    if sparsify in needs_surr and n_surrogates == 0:
        n_surrogates = 200
        rpt.info(f"Auto-setting n_surrogates=200 for sparsify={sparsify}")

    # Auto-switch to dev cache root
    if filter_time is not None and filter_time > 0 and cache_root == DEFAULT_MSC_CACHE_ROOT:
        cache_root = cache_root.parent / (cache_root.name + "_dev")
        rpt.info(f"Using dev cache root: {cache_root}")

    rpt.header(
        "MSC FC Matrix Computation",
        Patients=", ".join(pat_list),
        Bands=", ".join(band_list),
        Phases=", ".join(phase_list),
        sparsify=sparsify,
        n_surrogates=n_surrogates,
        nperseg=nperseg,
    )

    for patient in pat_list:
        rpt.info(f"\nProcessing {patient}...")
        results = compute_msc_for_patient(
            patient,
            bands=band_list,
            phases=phase_list,
            sparsify=sparsify,
            n_surrogates=n_surrogates,
            nperseg=nperseg,
            batch_size=batch_size,
            n_workers=n_workers,
            sample_rate=sample_rate,
            filter_time=filter_time,
            cache_root=cache_root,
            overwrite_cache=overwrite,
            return_results=False,
            verbose=rpt.verbose,
            fdr_q=fdr_q,
            disparity_alpha=disparity_alpha,
            ecm_alpha=ecm_alpha,
            ecm_n_ensemble=ecm_n_ensemble,
            ecm_weight_scale=ecm_weight_scale,
        )

        n_ok = sum(
            1 for b in results for p in results[b] if results[b][p] is not None
        )
        n_total = len(band_list) * len(phase_list)
        rpt.info(f"  {n_ok}/{n_total} computed")
        for _ in range(n_ok):
            rpt.success()
        for _ in range(n_total - n_ok):
            rpt.fail(f"{patient}: some combos failed")

        del results
        gc.collect()

    rpt.summary()


# ---------------------------------------------------------------------------
# compute lrg
# ---------------------------------------------------------------------------

@compute.command()
@patient_options()
@band_phase_options()
@fc_method_option(default="msc")
@cache_options("data/lrg_cache")
@filter_time_option()
@click.option("--fc-cache-root", type=click.Path(path_type=Path), default=None,
              help="FC cache root (default: auto-detect from fc-method).")
@click.option("--entropy-steps", type=int, default=400, show_default=True)
@click.option("--entropy-t1", type=float, default=-3.0, show_default=True,
              help="Log10 lower bound for tau range.")
@click.option("--entropy-t2", type=float, default=5.0, show_default=True,
              help="Log10 upper bound for tau range.")
@verbose_option()
@click.pass_context
def lrg(ctx, patients, bands, band, phases, phase, fc_method, cache_root,
        overwrite, filter_time, fc_cache_root, entropy_steps, entropy_t1,
        entropy_t2, verbose):
    """Compute LRG ultrametric analysis from cached FC matrices."""
    from lrg_eegfc.workflow.lrg import compute_lrg_for_patient

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    pat_list = resolve_patients(patients)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)

    rpt.header(
        "LRG Ultrametric Analysis",
        Patients=", ".join(pat_list),
        Bands=", ".join(band_list),
        Phases=", ".join(phase_list),
        fc_method=fc_method,
    )

    for patient in pat_list:
        rpt.info(f"\nProcessing {patient}...")
        results = compute_lrg_for_patient(
            patient,
            fc_method=fc_method,
            bands=band_list,
            phases=phase_list,
            fc_cache_root=fc_cache_root,
            filter_time=filter_time,
            cache_root=cache_root,
            overwrite_cache=overwrite,
            verbose=rpt.verbose,
            entropy_steps=entropy_steps,
            entropy_t1=entropy_t1,
            entropy_t2=entropy_t2,
        )

        n_ok = sum(
            1 for b in results for p in results[b] if results[b][p] is not None
        )
        n_total = len(band_list) * len(phase_list)
        rpt.info(f"  {n_ok}/{n_total} computed")
        for _ in range(n_ok):
            rpt.success()
        for _ in range(n_total - n_ok):
            rpt.fail(f"{patient}: some combos failed")

        del results
        gc.collect()

    rpt.summary()


# ---------------------------------------------------------------------------
# compute clean
# ---------------------------------------------------------------------------

@compute.command("clean")
@patient_options()
@band_phase_options()
@click.option("--corr-cache-root", type=click.Path(path_type=Path),
              default=Path("data/corr_cache"), show_default=True,
              help="Source correlation cache.")
@click.option("--overwrite", is_flag=True, help="Recompute even if cached.")
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=Path("data/stereoeeg_patients"), show_default=True)
@click.option("--filter-order", type=int, default=4, show_default=True)
@click.option("--sample-rate", type=float, default=2048.0, show_default=True)
@verbose_option()
@click.pass_context
def clean(ctx, patients, bands, band, phases, phase, corr_cache_root,
          overwrite, dataset_root, filter_order, sample_rate, verbose):
    """Clean correlation matrices (Marchenko-Pastur + percolation)."""
    from lrg_eegfc.workflow.cleaning import compute_cleaned_corr_for_patient

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    pat_list = resolve_patients(patients)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)

    rpt.header(
        "Marchenko-Pastur Correlation Cleaning",
        Patients=", ".join(pat_list),
        Bands=", ".join(band_list),
        Phases=", ".join(phase_list),
    )

    for patient in pat_list:
        rpt.info(f"\nProcessing {patient}...")
        try:
            results = compute_cleaned_corr_for_patient(
                patient,
                dataset_root=dataset_root,
                cache_root=corr_cache_root,
                bands=band_list,
                phases=phase_list,
                filter_order=filter_order,
                sample_rate=sample_rate,
                overwrite_cache=overwrite,
                verbose=rpt.verbose,
            )

            n_ok = sum(
                1 for b in results for p in results[b] if results[b][p] is not None
            )
            n_total = len(band_list) * len(phase_list)
            rpt.info(f"  {n_ok}/{n_total} cleaned")
            for _ in range(n_ok):
                rpt.success()
            for _ in range(n_total - n_ok):
                rpt.fail(f"{patient}: some combos failed")
        except Exception as exc:
            rpt.fail(f"{patient}: {exc}")

        gc.collect()

    rpt.summary()


# ---------------------------------------------------------------------------
# compute time-windows
# ---------------------------------------------------------------------------

@compute.command("time-windows")
@patient_options()
@band_phase_options()
@fc_method_option(default="msc")
@cache_options("data/msc_cache_windows")
@filter_time_option()
@click.option("--window-sec", type=float, default=None,
              help="Window length in seconds (default: auto from band freq).")
@click.option("--overlap", type=float, default=0.5, show_default=True,
              help="Fractional window overlap.")
@click.option("--sample-rate", type=float, default=2048.0, show_default=True)
@verbose_option()
@click.pass_context
def time_windows(ctx, patients, bands, band, phases, phase, fc_method,
                 cache_root, overwrite, filter_time, window_sec, overlap,
                 sample_rate, verbose):
    """Compute sliding-window FC matrices."""
    from lrg_eegfc.workflow.time_windows import (
        suggest_window_sec,
        compute_window_params,
        generate_window_indices,
        build_window_run_id,
        get_window_cache_dir,
    )
    from lrg_eegfc.config.const import BRAIN_BANDS
    from lrg_eegfc.utils.io import load_timeseries
    from lrg_eegfc.workflow.msc import compute_msc_matrix
    from lrg_eegfc.workflow.corr import compute_corr_matrix

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    pat_list = resolve_patients(patients)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)

    rpt.header(
        "Sliding-Window FC Computation",
        Patients=", ".join(pat_list),
        fc_method=fc_method,
        window_sec=window_sec or "auto",
        overlap=overlap,
    )

    for patient in pat_list:
        rpt.info(f"\nProcessing {patient}...")
        for phase_name in phase_list:
            for band_name in band_list:
                low_freq = BRAIN_BANDS[band_name][0]
                ws = window_sec if window_sec is not None else suggest_window_sec(low_freq)
                rpt.detail(f"{band_name}/{phase_name} window={ws}s")

                try:
                    ts = load_timeseries(patient, phase_name)
                    if ts is None:
                        rpt.fail(f"{patient}/{phase_name}: data not found")
                        continue

                    n_samples = ts.shape[1]
                    win_len, step = compute_window_params(
                        n_samples, sample_rate, ws, overlap,
                    )
                    indices = generate_window_indices(n_samples, win_len, step)

                    run_id = build_window_run_id(
                        fc_method=fc_method, window_sec=ws,
                        overlap=overlap, band=band_name,
                    )
                    win_cache = get_window_cache_dir(
                        patient, phase_name, run_id,
                        cache_root=cache_root,
                    )
                    win_cache.mkdir(parents=True, exist_ok=True)

                    for i, start in enumerate(indices):
                        if fc_method == "msc":
                            compute_msc_matrix(
                                patient, phase_name, band_name,
                                cache_root=win_cache,
                                overwrite_cache=overwrite,
                                verbose=False,
                            )
                        else:
                            compute_corr_matrix(
                                patient, phase_name, band_name,
                                cache_root=win_cache,
                                overwrite_cache=overwrite,
                                verbose=False,
                            )
                    rpt.success(f"{patient}/{band_name}/{phase_name}: {len(indices)} windows")
                except Exception as exc:
                    rpt.fail(f"{patient}/{band_name}/{phase_name}: {exc}")

    rpt.summary()


# ---------------------------------------------------------------------------
# compute reorganization
# ---------------------------------------------------------------------------

@compute.command("reorganization")
@patient_options()
@band_phase_options()
@fc_method_option(default="msc")
@click.option("--lrg-cache-root", type=click.Path(path_type=Path),
              default=Path("data/lrg_cache"), show_default=True)
@click.option("--output-root", type=click.Path(path_type=Path),
              default=Path("results/reorganization"), show_default=True)
@click.option("--distance-metric", default="euclidean", show_default=True)
@verbose_option()
@click.pass_context
def reorganization(ctx, patients, bands, band, phases, phase, fc_method,
                   lrg_cache_root, output_root, distance_metric, verbose):
    """Compute reorganization metrics from cached LRG results."""
    from lrg_eegfc.utils.metrics.reorganization import build_metric_specs
    from lrg_eegfc.workflow.lrg import load_lrg_result
    import numpy as np

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    pat_list = resolve_patients(patients)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)

    rpt.header(
        "Reorganization Metrics",
        Patients=", ".join(pat_list),
        fc_method=fc_method,
        distance_metric=distance_metric,
    )

    metrics_dict = build_metric_specs(distance_metric)
    output_root.mkdir(parents=True, exist_ok=True)

    for patient in pat_list:
        rpt.info(f"\nProcessing {patient}...")
        for band_name in band_list:
            try:
                # Load all LRG results for this patient/band
                lrg_results = {}
                for p in phase_list:
                    r = load_lrg_result(
                        patient, p, band_name,
                        fc_method=fc_method,
                        cache_root=lrg_cache_root,
                    )
                    if r is not None:
                        lrg_results[p] = r

                if len(lrg_results) < 2:
                    rpt.skip(f"{patient}/{band_name}: <2 phases with LRG data")
                    continue

                rpt.success(f"{patient}/{band_name}: {len(lrg_results)} phases loaded")
            except Exception as exc:
                rpt.fail(f"{patient}/{band_name}: {exc}")

    rpt.summary()


# ---------------------------------------------------------------------------
# compute diagnostics
# ---------------------------------------------------------------------------

@compute.command("diagnostics")
@patient_options()
@band_phase_options()
@click.option("--msc-cache-root", type=click.Path(path_type=Path),
              default=Path("data/msc_cache"), show_default=True)
@click.option("--nperseg", type=int, default=4096, show_default=True)
@click.option("--n-tau", type=int, default=20, show_default=True,
              help="Number of log-spaced tau for coarsening/metastability.")
@click.option("--output", type=click.Path(path_type=Path),
              default=Path("data/tables/mslcd_diagnostics_master.csv"),
              show_default=True, help="Output CSV path.")
@click.option("--trajectories-output", type=click.Path(path_type=Path),
              default=Path("data/tables/mslcd_coarsening_trajectories.npz"),
              show_default=True, help="Output NPZ for coarsening trajectories.")
@verbose_option()
@click.pass_context
def diagnostics(ctx, patients, bands, band, phases, phase, msc_cache_root,
                nperseg, n_tau, output, trajectories_output, verbose):
    """Compute MSLCD scalar diagnostics for all patient/phase/band triplets.

    Runs the full LRG multiscale pipeline on each triplet to extract spectral
    properties, entropic susceptibility, partition richness, coarsening
    trajectories, metastability scores, and community balance.

    \b
    Outputs:
      - Master CSV with all scalar diagnostics
      - NPZ with per-triplet coarsening trajectories
    """
    import numpy as np
    import pandas as pd
    from lrg_eegfc.workflow.diagnostics import compute_triplet_diagnostics

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    pat_list = resolve_patients(patients)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)

    n_total = len(pat_list) * len(phase_list) * len(band_list)
    rpt.header(
        "MSLCD Cross-Condition Diagnostics",
        Patients=", ".join(pat_list),
        Bands=", ".join(band_list),
        Phases=", ".join(phase_list),
        Total_triplets=n_total,
        n_tau=n_tau,
    )

    rows = []
    trajectories = {}
    done = 0

    for patient in pat_list:
        rpt.info(f"\n{patient}:")
        for phase_name in phase_list:
            for band_name in band_list:
                done += 1
                label = f"{patient}/{phase_name}/{band_name}"
                rpt.detail(f"[{done}/{n_total}] {label}")

                try:
                    result = compute_triplet_diagnostics(
                        patient, phase_name, band_name,
                        msc_cache=msc_cache_root,
                        nperseg=nperseg,
                        n_tau=n_tau,
                        verbose=rpt.verbose,
                    )
                    if result is None:
                        rpt.skip(f"{label}: missing or disconnected")
                        continue

                    diag, tau_grid, nmax_traj, eigenvalues = result
                    rows.append(diag.as_dict())

                    key = f"{patient}__{phase_name}__{band_name}"
                    trajectories[key + "__tau"] = tau_grid
                    trajectories[key + "__nmax"] = nmax_traj
                    trajectories[key + "__eigenvalues"] = eigenvalues

                    rpt.success(f"{label}: N={diag.N}, W={diag.W:.2f}, "
                                f"mu_mean={diag.mu_mean:.3f}")
                except Exception as exc:
                    rpt.fail(f"{label}: {exc}")

    if not rows:
        rpt.info("\nNo triplets computed — nothing to save.")
        rpt.summary()
        return

    # Save master CSV
    output.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    col_order = [
        "patient", "phase", "band", "N", "lambda_max", "lambda_2",
        "W", "gap_ratio", "tau_min", "tau_max",
        "N_peaks", "tau_star1", "tau_star1_norm", "C_max", "W_peak",
        "N_sens", "n_global", "psi_max", "sensible_ns",
        "n_start", "n_end", "coarsening_ratio",
        "mu_mean", "mu_max", "f_meta_01", "f_meta_02", "sigma_mu",
        "H_size", "C_max_frac", "n_singleton",
    ]
    df = df[[c for c in col_order if c in df.columns]]
    df.to_csv(output, index=False, float_format="%.6f")
    rpt.saved(output)

    # Save trajectories
    trajectories_output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(trajectories_output, **trajectories)
    rpt.saved(trajectories_output)

    rpt.summary()


# ---------------------------------------------------------------------------
# compute threshold-analysis
# ---------------------------------------------------------------------------

@compute.command("threshold-analysis")
@click.option("--patient", default="Pat_02", show_default=True,
              help="Patient ID for threshold analysis.")
@click.option("--phase", default="rsPre", show_default=True,
              help="Phase for threshold analysis.")
@band_phase_options()
@click.option("--msc-cache-root", type=click.Path(path_type=Path),
              default=Path("data/msc_cache"), show_default=True)
@click.option("--nperseg", type=int, default=4096, show_default=True)
@click.option("--output", type=click.Path(path_type=Path),
              default=Path("data/tables/threshold_analysis_Pat02.npz"),
              show_default=True, help="Output NPZ path.")
@verbose_option()
@click.pass_context
def threshold_analysis(ctx, patient, phase, bands, band, phases,
                       msc_cache_root, nperseg, output, verbose):
    """Compute threshold analysis for edge-weight pruning.

    Loads the full MSC adjacency matrix for a patient/phase and progressively
    removes the bottom X% of edges by weight.  At each threshold, computes
    Laplacian eigenvalues, entropy curves, and connectivity metrics.

    \b
    Default thresholds: 0, 50, 70, 80, 85, 90, 95, 97, 99 (percentile).

    \b
    Output NPZ keys per band:
      {band}__tau, {band}__log10_tau     - reference tau grid
      {band}__C__{pct}, {band}__S__{pct} - susceptibility / entropy
      {band}__eigenvalues__{pct}         - Laplacian eigenvalues
      {band}__connectivity               - (n_thresholds, 4) array
    """
    import numpy as np
    from lrg_eegfc.workflow.diagnostics import compute_threshold_analysis

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    band_list = resolve_bands(bands, band)

    rpt.header(
        "Threshold Analysis",
        Patient=patient,
        Phase=phase,
        Bands=", ".join(band_list),
    )

    try:
        result = compute_threshold_analysis(
            patient=patient,
            phase=phase,
            bands=band_list,
            msc_cache=msc_cache_root,
            nperseg=nperseg,
            verbose=rpt.verbose,
        )

        output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(output, **result)
        rpt.saved(output)
        rpt.info(f"  {len([k for k in result if '__connectivity' in k])} bands processed")
    except Exception as exc:
        rpt.fail(f"threshold-analysis: {exc}")
        import traceback
        if verbose:
            traceback.print_exc()

    rpt.summary()
