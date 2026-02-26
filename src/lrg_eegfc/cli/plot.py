"""``lrg-eegfc plot`` subcommands.

Each command loads from cache and delegates to a visuals module function.
No computation is performed here.
"""

from __future__ import annotations

from pathlib import Path

import click

from ._common import (
    CliReporter,
    band_phase_options,
    fc_method_option,
    output_options,
    resolve_bands,
    resolve_phases,
    show_option,
    single_patient_option,
    patient_options,
    resolve_patients,
    verbose_option,
)


@click.group()
def plot() -> None:
    """Generate visualizations from cached results."""


# ---------------------------------------------------------------------------
# plot corr
# ---------------------------------------------------------------------------

@plot.command()
@single_patient_option()
@band_phase_options()
@click.option("--plot-type",
              type=click.Choice(["summary", "network", "mp", "percolation", "all"]),
              default="summary", show_default=True)
@click.option("--cleaned", is_flag=True, help="Use MP-cleaned matrix.")
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/corr_cache"), show_default=True)
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=Path("data/stereoeeg_patients"), show_default=True)
@output_options("data/figures/correlation")
@show_option()
@verbose_option()
@click.pass_context
def corr(ctx, patient, band, bands, phase, phases, plot_type, cleaned,
         cache_root, dataset_root, output_dir, dpi, fmt, show_fig, verbose):
    """Visualize correlation FC matrices and networks."""
    from lrg_eegfc.visuals.correlation import (
        plot_correlation_summary,
        plot_correlation_and_network,
        plot_marchenko_pastur_comparison,
        plot_percolation_curves,
    )

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)
    out_sub = output_dir / patient
    out_sub.mkdir(parents=True, exist_ok=True)

    rpt.header(f"Correlation Visualization: {patient}", plot_type=plot_type)

    for b in band_list:
        for p in phase_list:
            try:
                if plot_type in ("summary", "all"):
                    out = out_sub / f"{b}_{p}_corr_summary.{fmt}"
                    if not out.exists():
                        plot_correlation_summary(
                            patient, p, b,
                            cache_root=cache_root,
                            dataset_root=dataset_root,
                            output_path=out,
                        )
                        rpt.saved(out, show=show_fig)
                    else:
                        rpt.skip(f"{b}/{p} summary (exists)", path=out, show=show_fig)

                if plot_type in ("network", "all"):
                    out = out_sub / f"{b}_{p}_corr_network.{fmt}"
                    if not out.exists():
                        plot_correlation_and_network(
                            patient, p, b,
                            cache_root=cache_root,
                            dataset_root=dataset_root,
                            cleaned=cleaned,
                            output_path=out,
                        )
                        rpt.saved(out, show=show_fig)
                    else:
                        rpt.skip(f"{b}/{p} network (exists)", path=out, show=show_fig)

                if plot_type in ("mp", "all"):
                    out = out_sub / f"{b}_{p}_corr_mp.{fmt}"
                    if not out.exists():
                        plot_marchenko_pastur_comparison(
                            patient, p, b,
                            cache_root=cache_root,
                            dataset_root=dataset_root,
                            output_path=out,
                        )
                        rpt.saved(out, show=show_fig)
                    else:
                        rpt.skip(f"{b}/{p} MP (exists)", path=out, show=show_fig)

                if plot_type in ("percolation", "all"):
                    out = out_sub / f"{b}_{p}_corr_percolation.{fmt}"
                    if not out.exists():
                        plot_percolation_curves(
                            patient, p, b,
                            cache_root=cache_root,
                            cleaned=cleaned,
                            dataset_root=dataset_root,
                            output_path=out,
                        )
                        rpt.saved(out, show=show_fig)
                    else:
                        rpt.skip(f"{b}/{p} percolation (exists)", path=out, show=show_fig)

            except Exception as exc:
                rpt.fail(f"{b}/{p}: {exc}")

    rpt.summary()


# ---------------------------------------------------------------------------
# plot msc
# ---------------------------------------------------------------------------

@plot.command()
@single_patient_option()
@band_phase_options()
@click.option("--plot-type",
              type=click.Choice(["network", "comparison", "summary", "all"]),
              default="summary", show_default=True)
@click.option("--sparsify", default="none", show_default=True)
@click.option("--n-surrogates", type=int, default=0, show_default=True)
@click.option("--nperseg", type=int, default=4096, show_default=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/msc_cache"), show_default=True)
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=Path("data/stereoeeg_patients"), show_default=True)
@output_options("data/figures/msc")
@show_option()
@verbose_option()
@click.pass_context
def msc(ctx, patient, band, bands, phase, phases, plot_type, sparsify,
        n_surrogates, nperseg, cache_root, dataset_root, output_dir, dpi,
        fmt, show_fig, verbose):
    """Visualize MSC FC matrices and networks."""
    from lrg_eegfc.visuals.msc import (
        plot_msc_and_network,
        plot_msc_comparison_dense_vs_validated,
        plot_msc_summary,
    )

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)
    out_sub = output_dir / patient
    out_sub.mkdir(parents=True, exist_ok=True)

    rpt.header(f"MSC Visualization: {patient}", plot_type=plot_type, sparsify=sparsify)

    for b in band_list:
        for p in phase_list:
            try:
                if plot_type in ("network", "all"):
                    out = out_sub / f"{b}_{p}_msc_network.{fmt}"
                    if not out.exists():
                        plot_msc_and_network(
                            patient, p, b,
                            cache_root=cache_root,
                            sparsify=sparsify,
                            n_surrogates=n_surrogates,
                            nperseg=nperseg,
                            dataset_root=dataset_root,
                            output_path=out,
                        )
                        rpt.saved(out, show=show_fig)
                    else:
                        rpt.skip(f"{b}/{p} network (exists)", path=out, show=show_fig)

                if plot_type in ("comparison", "all"):
                    out = out_sub / f"{b}_{p}_msc_comparison.{fmt}"
                    if not out.exists():
                        plot_msc_comparison_dense_vs_validated(
                            patient, p, b,
                            n_surrogates=n_surrogates or 200,
                            nperseg=nperseg,
                            cache_root=cache_root,
                            output_path=out,
                        )
                        rpt.saved(out, show=show_fig)
                    else:
                        rpt.skip(f"{b}/{p} comparison (exists)", path=out, show=show_fig)

                if plot_type in ("summary", "all"):
                    out = out_sub / f"{b}_{p}_msc_summary.{fmt}"
                    if not out.exists():
                        plot_msc_summary(
                            patient, p, b,
                            cache_root=cache_root,
                            nperseg=nperseg,
                            dataset_root=dataset_root,
                            output_path=out,
                        )
                        rpt.saved(out, show=show_fig)
                    else:
                        rpt.skip(f"{b}/{p} summary (exists)", path=out, show=show_fig)

            except Exception as exc:
                rpt.fail(f"{b}/{p}: {exc}")

    rpt.summary()


# ---------------------------------------------------------------------------
# plot lrg
# ---------------------------------------------------------------------------

@plot.command()
@single_patient_option()
@band_phase_options()
@fc_method_option(required=True)
@click.option("--plot-type",
              type=click.Choice(["entropy", "dendrogram", "ultrametric", "full", "all"]),
              default="full", show_default=True)
@click.option("--orientation",
              type=click.Choice(["top", "right", "bottom", "left"]),
              default="top", show_default=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/lrg_cache"), show_default=True)
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=Path("data/stereoeeg_patients"), show_default=True)
@output_options("data/figures/lrg")
@show_option()
@verbose_option()
@click.pass_context
def lrg(ctx, patient, band, bands, phase, phases, fc_method, plot_type,
        orientation, cache_root, dataset_root, output_dir, dpi, fmt,
        show_fig, verbose):
    """Visualize LRG hierarchical analysis."""
    from lrg_eegfc.visuals.lrg import (
        plot_lrg_dendrogram,
        plot_lrg_entropy_curves,
        plot_lrg_full_panel,
        plot_ultrametric_heatmap,
    )

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)
    out_sub = output_dir / patient
    out_sub.mkdir(parents=True, exist_ok=True)

    rpt.header(
        f"LRG Visualization: {patient}",
        fc_method=fc_method,
        plot_type=plot_type,
    )

    for b in band_list:
        for p in phase_list:
            try:
                if plot_type in ("entropy", "all"):
                    out = out_sub / f"{b}_{p}_lrg_{fc_method}_entropy.{fmt}"
                    if not out.exists():
                        plot_lrg_entropy_curves(
                            patient, p, b, fc_method,
                            cache_root=cache_root,
                            output_path=out,
                        )
                        rpt.saved(out, show=show_fig)
                    else:
                        rpt.skip(f"{b}/{p} entropy (exists)", path=out, show=show_fig)

                if plot_type in ("dendrogram", "all"):
                    out = out_sub / f"{b}_{p}_lrg_{fc_method}_dendrogram.{fmt}"
                    if not out.exists():
                        plot_lrg_dendrogram(
                            patient, p, b, fc_method,
                            cache_root=cache_root,
                            orientation=orientation,
                            dataset_root=dataset_root,
                            output_path=out,
                        )
                        rpt.saved(out, show=show_fig)
                    else:
                        rpt.skip(f"{b}/{p} dendrogram (exists)", path=out, show=show_fig)

                if plot_type in ("ultrametric", "all"):
                    out = out_sub / f"{b}_{p}_lrg_{fc_method}_ultrametric.{fmt}"
                    if not out.exists():
                        plot_ultrametric_heatmap(
                            patient, p, b, fc_method,
                            cache_root=cache_root,
                            dataset_root=dataset_root,
                            output_path=out,
                        )
                        rpt.saved(out, show=show_fig)
                    else:
                        rpt.skip(f"{b}/{p} ultrametric (exists)", path=out, show=show_fig)

                if plot_type in ("full", "all"):
                    out = out_sub / f"{b}_{p}_lrg_{fc_method}_full.{fmt}"
                    if not out.exists():
                        plot_lrg_full_panel(
                            patient, p, b, fc_method,
                            cache_root=cache_root,
                            dataset_root=dataset_root,
                            output_path=out,
                            verbose=rpt.verbose,
                        )
                        rpt.saved(out, show=show_fig)
                    else:
                        rpt.skip(f"{b}/{p} full panel (exists)", path=out, show=show_fig)

            except Exception as exc:
                rpt.fail(f"{b}/{p}: {exc}")

    rpt.summary()


# ---------------------------------------------------------------------------
# plot reorganization
# ---------------------------------------------------------------------------

@plot.command("reorganization")
@single_patient_option()
@band_phase_options()
@fc_method_option(required=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/lrg_cache"), show_default=True)
@output_options("data/figures/reorganization")
@show_option()
@verbose_option()
@click.pass_context
def reorganization(ctx, patient, band, bands, phase, phases, fc_method,
                   cache_root, output_dir, dpi, fmt, show_fig, verbose):
    """Visualize phase reorganization of brain networks."""
    from lrg_eegfc.visuals.reorganization import plot_phase_reorganization

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)
    out_sub = output_dir / patient
    out_sub.mkdir(parents=True, exist_ok=True)

    rpt.header(f"Phase Reorganization: {patient}", fc_method=fc_method)

    for b in band_list:
        try:
            out = out_sub / f"{b}_reorganization_{fc_method}.{fmt}"
            if not out.exists():
                plot_phase_reorganization(
                    patient, b, fc_method,
                    phases=phase_list,
                    cache_root=cache_root,
                    output_path=out,
                    verbose=rpt.verbose,
                )
                rpt.saved(out, show=show_fig)
            else:
                rpt.skip(f"{b} (exists)", path=out, show=show_fig)
        except Exception as exc:
            rpt.fail(f"{b}: {exc}")

    rpt.summary()


# ---------------------------------------------------------------------------
# plot metastable
# ---------------------------------------------------------------------------

@plot.command("metastable")
@single_patient_option()
@click.option("--band", required=True)
@click.option("--phase", required=True)
@fc_method_option(required=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/lrg_cache"), show_default=True)
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=Path("data/stereoeeg_patients"), show_default=True)
@output_options("data/figures/metastable")
@show_option()
@verbose_option()
@click.pass_context
def metastable(ctx, patient, band, phase, fc_method, cache_root,
               dataset_root, output_dir, dpi, fmt, show_fig, verbose):
    """Generate metastable nodes Sankey diagrams."""
    from lrg_eegfc.visuals.metastable import (
        compute_clustering_across_tau,
        create_sankey_diagram,
    )
    from lrg_eegfc.workflow.lrg import load_lrg_result

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    out_sub = output_dir / patient
    out_sub.mkdir(parents=True, exist_ok=True)

    rpt.header(f"Metastable Sankey: {patient}", band=band, phase=phase)

    try:
        lrg_result = load_lrg_result(
            patient, phase, band,
            fc_method=fc_method,
            cache_root=cache_root,
        )
        if lrg_result is None:
            rpt.fail(f"No LRG data for {patient}/{phase}/{band}")
            return

        partitions, tau_values = compute_clustering_across_tau(
            lrg_result.linkage_matrix,
        )

        fig = create_sankey_diagram(
            partitions, tau_values,
            title=f"{patient} {band} {phase} ({fc_method})",
        )

        out = out_sub / f"{band}_{phase}_sankey_{fc_method}.html"
        fig.write_html(str(out))
        rpt.saved(out, show=show_fig)
    except Exception as exc:
        rpt.fail(str(exc))

    rpt.summary()


# ---------------------------------------------------------------------------
# plot cleaning
# ---------------------------------------------------------------------------

@plot.command("cleaning")
@single_patient_option()
@band_phase_options()
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/corr_cache"), show_default=True)
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=Path("data/stereoeeg_patients"), show_default=True)
@output_options("data/figures/cleaning")
@show_option()
@verbose_option()
@click.pass_context
def cleaning(ctx, patient, band, bands, phase, phases, cache_root,
             dataset_root, output_dir, dpi, fmt, show_fig, verbose):
    """Plot Marchenko-Pastur cleaning diagnostics."""
    from lrg_eegfc.visuals.correlation import plot_marchenko_pastur_comparison

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)
    out_sub = output_dir / patient
    out_sub.mkdir(parents=True, exist_ok=True)

    rpt.header(f"MP Cleaning Diagnostics: {patient}")

    for b in band_list:
        for p in phase_list:
            try:
                out = out_sub / f"{b}_{p}_mp_cleaning.{fmt}"
                if not out.exists():
                    plot_marchenko_pastur_comparison(
                        patient, p, b,
                        cache_root=cache_root,
                        dataset_root=dataset_root,
                        output_path=out,
                    )
                    rpt.saved(out, show=show_fig)
                else:
                    rpt.skip(f"{b}/{p} (exists)", path=out, show=show_fig)
            except Exception as exc:
                rpt.fail(f"{b}/{p}: {exc}")

    rpt.summary()


# ---------------------------------------------------------------------------
# plot comparison
# ---------------------------------------------------------------------------

@plot.command("comparison")
@single_patient_option()
@band_phase_options()
@click.option("--corr-cache", type=click.Path(path_type=Path),
              default=Path("data/corr_cache"), show_default=True)
@click.option("--msc-cache", type=click.Path(path_type=Path),
              default=Path("data/msc_cache"), show_default=True)
@click.option("--nperseg", type=int, default=4096, show_default=True)
@output_options("data/figures/comparison")
@show_option()
@verbose_option()
@click.pass_context
def comparison(ctx, patient, band, bands, phase, phases,
               corr_cache, msc_cache, nperseg, output_dir, dpi, fmt,
               show_fig, verbose):
    """Compare correlation vs MSC FC methods visually."""
    from lrg_eegfc.visuals.compare import plot_fc_comparison

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)
    out_sub = output_dir / patient
    out_sub.mkdir(parents=True, exist_ok=True)

    rpt.header(f"FC Method Comparison: {patient}")

    for b in band_list:
        for p in phase_list:
            try:
                out = out_sub / f"{b}_{p}_fc_comparison.{fmt}"
                if not out.exists():
                    plot_fc_comparison(
                        patient, p, b,
                        cache_root_corr=corr_cache,
                        cache_root_msc=msc_cache,
                        nperseg=nperseg,
                        output_path=out,
                    )
                    rpt.saved(out, show=show_fig)
                else:
                    rpt.skip(f"{b}/{p} (exists)", path=out, show=show_fig)
            except Exception as exc:
                rpt.fail(f"{b}/{p}: {exc}")

    rpt.summary()


# ---------------------------------------------------------------------------
# plot msc-grid
# ---------------------------------------------------------------------------

@plot.command("msc-grid")
@single_patient_option()
@click.option("--sparsify", default="none", show_default=True)
@click.option("--nperseg", type=int, default=4096, show_default=True)
@click.option("--n-surrogates", type=int, default=0, show_default=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/msc_cache"), show_default=True)
@output_options("data/figures/msc")
@show_option()
@verbose_option()
@click.pass_context
def msc_grid(ctx, patient, sparsify, nperseg, n_surrogates, cache_root,
             output_dir, dpi, fmt, show_fig, verbose):
    """Plot band x phase grid of MSC matrices."""
    rpt = CliReporter.from_context(ctx, verbose=verbose)
    rpt.header(f"MSC Grid: {patient}", sparsify=sparsify)

    try:
        from lrg_eegfc.workflow.msc import load_msc_matrix
        from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
        import matplotlib.pyplot as plt

        band_list = list(BRAIN_BANDS.keys())
        phase_list = list(PHASE_LABELS)
        n_bands = len(band_list)
        n_phases = len(phase_list)

        fig, axes = plt.subplots(n_bands, n_phases, figsize=(4 * n_phases, 4 * n_bands))
        for i, b in enumerate(band_list):
            for j, p in enumerate(phase_list):
                ax = axes[i, j]
                mat = load_msc_matrix(
                    patient, p, b,
                    sparsify=sparsify,
                    n_surrogates=n_surrogates,
                    nperseg=nperseg,
                    cache_root=cache_root,
                )
                if mat is not None:
                    ax.imshow(mat, cmap="viridis", vmin=0, vmax=1, aspect="equal")
                else:
                    ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                            transform=ax.transAxes)
                ax.set_title(f"{b} / {p}", fontsize=8)
                ax.tick_params(labelbottom=False, labelleft=False)

        fig.suptitle(f"{patient} MSC Grid (sparsify={sparsify})", fontsize=14)
        fig.tight_layout()
        out = output_dir / patient / f"msc_grid_{sparsify}.{fmt}"
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        rpt.saved(out, show=show_fig)
    except Exception as exc:
        rpt.fail(str(exc))

    rpt.summary()


# ---------------------------------------------------------------------------
# plot msc-all-patients
# ---------------------------------------------------------------------------

@plot.command("msc-all-patients")
@patient_options()
@band_phase_options()
@click.option("--sparsify", default="none", show_default=True)
@click.option("--nperseg", type=int, default=4096, show_default=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/msc_cache"), show_default=True)
@output_options("data/figures/msc")
@show_option()
@verbose_option()
@click.pass_context
def msc_all_patients(ctx, patients, band, bands, phase, phases, sparsify,
                     nperseg, cache_root, output_dir, dpi, fmt, show_fig,
                     verbose):
    """Plot MSC matrices across all patients."""
    rpt = CliReporter.from_context(ctx, verbose=verbose)
    pat_list = resolve_patients(patients)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)

    rpt.header("MSC All Patients", Patients=", ".join(pat_list))

    try:
        from lrg_eegfc.workflow.msc import load_msc_matrix
        import matplotlib.pyplot as plt

        for b in band_list:
            for p in phase_list:
                fig, axes = plt.subplots(1, len(pat_list),
                                         figsize=(4 * len(pat_list), 4))
                if len(pat_list) == 1:
                    axes = [axes]
                for k, pat in enumerate(pat_list):
                    mat = load_msc_matrix(
                        pat, p, b, sparsify=sparsify, nperseg=nperseg,
                        cache_root=cache_root,
                    )
                    if mat is not None:
                        axes[k].imshow(mat, cmap="viridis", vmin=0, vmax=1)
                    else:
                        axes[k].text(0.5, 0.5, "N/A", ha="center", va="center",
                                     transform=axes[k].transAxes)
                    axes[k].set_title(pat, fontsize=9)
                    axes[k].tick_params(labelbottom=False, labelleft=False)

                fig.suptitle(f"{b} / {p} (sparsify={sparsify})")
                fig.tight_layout()
                out = output_dir / f"all_patients_{b}_{p}_{sparsify}.{fmt}"
                out.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(out, dpi=dpi, bbox_inches="tight")
                plt.close(fig)
                rpt.saved(out, show=show_fig)
    except Exception as exc:
        rpt.fail(str(exc))

    rpt.summary()


# ---------------------------------------------------------------------------
# plot lrg-phase-grid
# ---------------------------------------------------------------------------

@plot.command("lrg-phase-grid")
@single_patient_option()
@band_phase_options()
@fc_method_option(required=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/lrg_cache"), show_default=True)
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=Path("data/stereoeeg_patients"), show_default=True)
@output_options("data/figures/lrg")
@show_option()
@verbose_option()
@click.pass_context
def lrg_phase_grid(ctx, patient, band, bands, phase, phases, fc_method,
                   cache_root, dataset_root, output_dir, dpi, fmt, show_fig,
                   verbose):
    """Plot LRG full panels across phases."""
    from lrg_eegfc.visuals.lrg import plot_lrg_full_panel

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)
    out_sub = output_dir / patient
    out_sub.mkdir(parents=True, exist_ok=True)

    rpt.header(f"LRG Phase Grid: {patient}", fc_method=fc_method)

    for b in band_list:
        for p in phase_list:
            try:
                out = out_sub / f"{b}_{p}_lrg_{fc_method}_full.{fmt}"
                if not out.exists():
                    plot_lrg_full_panel(
                        patient, p, b, fc_method,
                        cache_root=cache_root,
                        dataset_root=dataset_root,
                        output_path=out,
                        verbose=rpt.verbose,
                    )
                    rpt.saved(out, show=show_fig)
                else:
                    rpt.skip(f"{b}/{p} (exists)", path=out, show=show_fig)
            except Exception as exc:
                rpt.fail(f"{b}/{p}: {exc}")

    rpt.summary()


# ---------------------------------------------------------------------------
# plot lrg-video
# ---------------------------------------------------------------------------

@plot.command("lrg-video")
@single_patient_option()
@click.option("--band", required=True)
@click.option("--phase", required=True)
@fc_method_option(required=True)
@click.option("--n-frames", type=int, default=50, show_default=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/lrg_cache"), show_default=True)
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=Path("data/stereoeeg_patients"), show_default=True)
@output_options("data/figures/lrg_video")
@verbose_option()
@click.pass_context
def lrg_video(ctx, patient, band, phase, fc_method, n_frames, cache_root,
              dataset_root, output_dir, dpi, fmt, verbose):
    """Generate ultrametric threshold animation frames."""
    rpt = CliReporter.from_context(ctx, verbose=verbose)
    rpt.header(f"LRG Video Frames: {patient}", band=band, phase=phase)
    # TODO: integrate src/visualize_lrg_video.py logic
    rpt.info("Not yet implemented in CLI. Use: python src/visualize_lrg_video.py")


# ---------------------------------------------------------------------------
# plot msc-validation
# ---------------------------------------------------------------------------

@plot.command("msc-validation")
@single_patient_option()
@band_phase_options()
@click.option("--nperseg", type=int, default=4096, show_default=True)
@click.option("--n-surrogates", type=int, default=200, show_default=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/msc_cache"), show_default=True)
@output_options("data/figures/msc_validation")
@show_option()
@verbose_option()
@click.pass_context
def msc_validation(ctx, patient, band, bands, phase, phases, nperseg,
                   n_surrogates, cache_root, output_dir, dpi, fmt, show_fig,
                   verbose):
    """Dense vs validated MSC comparison figures."""
    from lrg_eegfc.visuals.msc import plot_msc_comparison_dense_vs_validated

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)
    out_sub = output_dir / patient
    out_sub.mkdir(parents=True, exist_ok=True)

    rpt.header(f"MSC Validation: {patient}")

    for b in band_list:
        for p in phase_list:
            try:
                out = out_sub / f"{b}_{p}_msc_validation.{fmt}"
                if not out.exists():
                    plot_msc_comparison_dense_vs_validated(
                        patient, p, b,
                        n_surrogates=n_surrogates,
                        nperseg=nperseg,
                        cache_root=cache_root,
                        output_path=out,
                    )
                    rpt.saved(out, show=show_fig)
                else:
                    rpt.skip(f"{b}/{p} (exists)", path=out, show=show_fig)
            except Exception as exc:
                rpt.fail(f"{b}/{p}: {exc}")

    rpt.summary()


# ---------------------------------------------------------------------------
# plot time-windows
# ---------------------------------------------------------------------------

@plot.command("time-windows")
@single_patient_option()
@band_phase_options()
@fc_method_option(default="msc")
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/msc_cache_windows"), show_default=True)
@output_options("data/figures/time_windows")
@verbose_option()
@click.pass_context
def time_windows(ctx, patient, band, bands, phase, phases, fc_method,
                 cache_root, output_dir, dpi, fmt, verbose):
    """Visualize cached time-window FC matrices."""
    rpt = CliReporter.from_context(ctx, verbose=verbose)
    rpt.header(f"Time Window Visualization: {patient}")
    # TODO: integrate src/visualize_time_windows.py logic
    rpt.info("Not yet implemented in CLI. Use: python src/visualize_time_windows.py")


# ---------------------------------------------------------------------------
# plot reorg-metrics
# ---------------------------------------------------------------------------

@plot.command("reorg-metrics")
@single_patient_option()
@band_phase_options()
@fc_method_option(default="msc")
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/lrg_cache"), show_default=True)
@output_options("data/figures/reorganization")
@show_option()
@verbose_option()
@click.pass_context
def reorg_metrics(ctx, patient, band, bands, phase, phases, fc_method,
                  cache_root, output_dir, dpi, fmt, show_fig, verbose):
    """Plot reorganization metric matrices."""
    from lrg_eegfc.visuals.reorganization import plot_reorganization_distance_matrix

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)
    out_sub = output_dir / patient
    out_sub.mkdir(parents=True, exist_ok=True)

    rpt.header(f"Reorganization Metrics: {patient}", fc_method=fc_method)

    for b in band_list:
        try:
            out = out_sub / f"{b}_reorg_metrics_{fc_method}.{fmt}"
            if not out.exists():
                plot_reorganization_distance_matrix(
                    patient, b, fc_method,
                    phases=phase_list,
                    cache_root=cache_root,
                    output_path=out,
                    verbose=rpt.verbose,
                )
                rpt.saved(out, show=show_fig)
            else:
                rpt.skip(f"{b} (exists)", path=out, show=show_fig)
        except Exception as exc:
            rpt.fail(f"{b}: {exc}")

    rpt.summary()


# ---------------------------------------------------------------------------
# plot reorg-summary
# ---------------------------------------------------------------------------

@plot.command("reorg-summary")
@single_patient_option()
@fc_method_option(default="msc")
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/lrg_cache"), show_default=True)
@output_options("data/figures/reorganization")
@verbose_option()
@click.pass_context
def reorg_summary(ctx, patient, fc_method, cache_root, output_dir, dpi, fmt,
                  verbose):
    """Band-level reorganization metric summary."""
    rpt = CliReporter.from_context(ctx, verbose=verbose)
    rpt.header(f"Reorganization Summary: {patient}", fc_method=fc_method)
    # TODO: integrate src/visualize_reorganization_summary.py logic
    rpt.info("Not yet implemented in CLI. Use: python src/visualize_reorganization_summary.py")


# ---------------------------------------------------------------------------
# plot metric-correlation
# ---------------------------------------------------------------------------

@plot.command("metric-correlation")
@single_patient_option()
@fc_method_option(default="msc")
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/lrg_cache"), show_default=True)
@output_options("data/figures/reorganization")
@verbose_option()
@click.pass_context
def metric_correlation(ctx, patient, fc_method, cache_root, output_dir, dpi,
                       fmt, verbose):
    """Cross-metric agreement analysis."""
    rpt = CliReporter.from_context(ctx, verbose=verbose)
    rpt.header(f"Metric Correlation: {patient}", fc_method=fc_method)
    # TODO: integrate src/visualize_metric_correlation.py logic
    rpt.info("Not yet implemented in CLI. Use: python src/visualize_metric_correlation.py")


# ---------------------------------------------------------------------------
# plot cross
# ---------------------------------------------------------------------------

@plot.command("cross")
@click.option("--figure-type",
              type=click.Choice([
                  "coarsening-beta", "coarsening-patient",
                  "metastability-heatmaps", "metastability-boxplots",
                  "partition-richness", "community-balance",
                  "threshold-susceptibility", "threshold-connectivity",
                  "latex-table", "all",
              ]),
              default="all", show_default=True,
              help="Which cross-condition figure to generate.")
@click.option("--csv", "csv_path", type=click.Path(path_type=Path),
              default=Path("data/tables/mslcd_diagnostics_master.csv"),
              show_default=True, help="Master diagnostics CSV.")
@click.option("--trajectories", type=click.Path(path_type=Path),
              default=Path("data/tables/mslcd_coarsening_trajectories.npz"),
              show_default=True, help="Coarsening trajectories NPZ.")
@click.option("--threshold-data", type=click.Path(path_type=Path),
              default=Path("data/tables/threshold_analysis_Pat02.npz"),
              show_default=True, help="Threshold analysis NPZ.")
@output_options("data/figures/report_mslcd_section/cross_condition")
@click.option("--threshold-output-dir", type=click.Path(path_type=Path),
              default=Path("data/figures/report_mslcd_section/threshold_analysis"),
              show_default=True, help="Output dir for threshold figures.")
@click.option("--latex-output", type=click.Path(path_type=Path),
              default=Path("data/tables/mslcd_diagnostics_Pat02.tex"),
              show_default=True)
@show_option()
@verbose_option()
@click.pass_context
def cross(ctx, figure_type, csv_path, trajectories, threshold_data,
          output_dir, threshold_output_dir, dpi, fmt, latex_output, show_fig,
          verbose):
    """Generate cross-condition MSLCD diagnostic figures (Section 4.8).

    \b
    Reads from the master CSV produced by ``compute diagnostics`` and generates
    publication-quality figures comparing diagnostics across patients, phases,
    and frequency bands.

    \b
    Figure types:
      coarsening-beta          B1 - Coarsening trajectories (beta, all patients)
      coarsening-patient       B2 - Coarsening trajectories (Pat_02, all bands)
      metastability-heatmaps   B3 - Metastability mu_mean / mu_max heatmaps
      metastability-boxplots   B4 - Metastability boxplots with patient overlay
      partition-richness       B5 - N_sens heatmap
      community-balance        B6 - H_size heatmap
      threshold-susceptibility C1 - C(tau) at different thresholds
      threshold-connectivity   C2 - Connectivity metrics vs threshold
      latex-table              Part D - LaTeX diagnostic table
      all                      Generate everything
    """
    from lrg_eegfc.visuals.cross_condition import (
        generate_cross_condition_figure,
        FIGURE_TYPES,
    )

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    rpt.header(
        "Cross-Condition Figures",
        figure_type=figure_type,
        csv=str(csv_path),
        output_dir=str(output_dir),
    )

    # Check that required data files exist
    needs_csv = {
        "coarsening-beta", "coarsening-patient",
        "metastability-heatmaps", "metastability-boxplots",
        "partition-richness", "community-balance", "latex-table",
    }
    needs_threshold = {
        "threshold-susceptibility", "threshold-spectrum",
        "threshold-connectivity",
    }
    targets = (
        set(FIGURE_TYPES.keys()) - {"all"}
        if figure_type == "all" else {figure_type}
    )

    if (needs_csv & targets) and not csv_path.exists():
        rpt.fail(f"Master CSV not found: {csv_path}")
        rpt.info("Run first: lrg-eegfc compute diagnostics")
        return

    if (needs_threshold & targets) and not threshold_data.exists():
        rpt.fail(f"Threshold data not found: {threshold_data}")
        rpt.info("Run first: lrg-eegfc compute threshold-analysis")
        return

    try:
        paths = generate_cross_condition_figure(
            figure_type=figure_type,
            csv_path=csv_path,
            trajectories_path=trajectories,
            threshold_path=threshold_data,
            output_dir=output_dir,
            threshold_output_dir=threshold_output_dir,
            latex_output=latex_output,
            dpi=dpi,
        )
        for p in paths:
            rpt.saved(p, show=show_fig)
    except Exception as exc:
        rpt.fail(f"{figure_type}: {exc}")
        import traceback
        if verbose:
            traceback.print_exc()

    rpt.summary()
