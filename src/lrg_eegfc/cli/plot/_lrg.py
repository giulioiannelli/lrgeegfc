"""``lrg-eegfc plot`` subcommands -- LRG dendrogram/panel family.

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
              default=LRG_CACHE, show_default=True)
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=SEEG_DATAPATH, show_default=True)
@output_options(str(FIGURES_ROOT / "lrg"))
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
        plot_ultrametric_heatmap,
    )
    from lrg_eegfc.visuals.lrg_panels import plot_lrg_full_panel

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
# plot lrg-phase-grid
# ---------------------------------------------------------------------------

@plot.command("lrg-phase-grid")
@single_patient_option()
@band_phase_options()
@fc_method_option(required=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=LRG_CACHE, show_default=True)
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=SEEG_DATAPATH, show_default=True)
@output_options(str(FIGURES_ROOT / "lrg"))
@show_option()
@verbose_option()
@click.pass_context
def lrg_phase_grid(ctx, patient, band, bands, phase, phases, fc_method,
                   cache_root, dataset_root, output_dir, dpi, fmt, show_fig,
                   verbose):
    """Plot LRG full panels across phases."""
    from lrg_eegfc.visuals.lrg_panels import plot_lrg_full_panel

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
              default=LRG_CACHE, show_default=True)
@click.option("--dataset-root", type=click.Path(path_type=Path),
              default=SEEG_DATAPATH, show_default=True)
@output_options(str(FIGURES_ROOT / "lrg_video"))
@verbose_option()
@click.pass_context
def lrg_video(ctx, patient, band, phase, fc_method, n_frames, cache_root,
              dataset_root, output_dir, dpi, fmt, verbose):
    """Generate ultrametric threshold animation frames."""
    rpt = CliReporter.from_context(ctx, verbose=verbose)
    rpt.header(f"LRG Video Frames: {patient}", band=band, phase=phase)
    # TODO: integrate src/visualize_lrg_video.py logic
    rpt.info("Not yet implemented in CLI. Use: python src/visualize_lrg_video.py")
