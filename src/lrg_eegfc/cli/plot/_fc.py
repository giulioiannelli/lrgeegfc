"""``lrg-eegfc plot`` subcommands -- FC adjacency-matrix family.

Split out of cli/plot.py on 2026-05-29 (Phase 4-B split 5/7). Each
``@plot.command(...)`` here registers itself on the ``plot`` click group
defined in this package's ``__init__.py``.
"""
from __future__ import annotations

from pathlib import Path

import click

from lrg_eegfc.config.paths import (
    CORR_CACHE,
    MSC_CACHE,
    LRG_CACHE,
    SEEG_DATAPATH,
    FIGURES_ROOT,
    MSC_WINDOWS_CACHE,
    TABLES_ROOT,
)

from .._common import (
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

# The `plot` click group is defined in __init__.py; importing it here lets
# each `@plot.command(...)` decorator below register on the shared group.
from . import plot


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
              default=CORR_CACHE, show_default=True)
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=SEEG_DATAPATH, show_default=True)
@output_options(str(FIGURES_ROOT / "correlation"))
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
              default=MSC_CACHE, show_default=True)
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=SEEG_DATAPATH, show_default=True)
@output_options(str(FIGURES_ROOT / "msc"))
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
# plot msc-grid
# ---------------------------------------------------------------------------

@plot.command("msc-grid")
@single_patient_option()
@click.option("--sparsify", default="none", show_default=True)
@click.option("--nperseg", type=int, default=4096, show_default=True)
@click.option("--n-surrogates", type=int, default=0, show_default=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=MSC_CACHE, show_default=True)
@output_options(str(FIGURES_ROOT / "msc"))
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
              default=MSC_CACHE, show_default=True)
@output_options(str(FIGURES_ROOT / "msc"))
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
# plot msc-validation
# ---------------------------------------------------------------------------

@plot.command("msc-validation")
@single_patient_option()
@band_phase_options()
@click.option("--nperseg", type=int, default=4096, show_default=True)
@click.option("--n-surrogates", type=int, default=200, show_default=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=MSC_CACHE, show_default=True)
@output_options(str(FIGURES_ROOT / "msc_validation"))
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
# plot cleaning
# ---------------------------------------------------------------------------

@plot.command("cleaning")
@single_patient_option()
@band_phase_options()
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=CORR_CACHE, show_default=True)
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=SEEG_DATAPATH, show_default=True)
@output_options(str(FIGURES_ROOT / "cleaning"))
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
              default=CORR_CACHE, show_default=True)
@click.option("--msc-cache", type=click.Path(path_type=Path),
              default=MSC_CACHE, show_default=True)
@click.option("--nperseg", type=int, default=4096, show_default=True)
@output_options(str(FIGURES_ROOT / "comparison"))
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
