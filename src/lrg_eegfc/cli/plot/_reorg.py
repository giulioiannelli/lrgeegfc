"""``lrg-eegfc plot`` subcommands -- reorganization metrics family.

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
# plot reorganization
# ---------------------------------------------------------------------------

@plot.command("reorganization")
@single_patient_option()
@band_phase_options()
@fc_method_option(required=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=LRG_CACHE, show_default=True)
@output_options(str(FIGURES_ROOT / "reorganization"))
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
# plot reorg-metrics
# ---------------------------------------------------------------------------

@plot.command("reorg-metrics")
@single_patient_option()
@band_phase_options()
@fc_method_option(default="imcoh_abs")
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=LRG_CACHE, show_default=True)
@output_options(str(FIGURES_ROOT / "reorganization"))
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
@fc_method_option(default="imcoh_abs")
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=LRG_CACHE, show_default=True)
@output_options(str(FIGURES_ROOT / "reorganization"))
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
@fc_method_option(default="imcoh_abs")
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=LRG_CACHE, show_default=True)
@output_options(str(FIGURES_ROOT / "reorganization"))
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
# plot metastable
# ---------------------------------------------------------------------------

@plot.command("metastable")
@single_patient_option()
@click.option("--band", required=True)
@click.option("--phase", required=True)
@fc_method_option(required=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=LRG_CACHE, show_default=True)
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=SEEG_DATAPATH, show_default=True)
@output_options(str(FIGURES_ROOT / "metastable"))
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
