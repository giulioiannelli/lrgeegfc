"""``lrg-eegfc plot`` subcommands -- miscellaneous (cross-condition, time-windows) family.

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
              default=TABLES_ROOT / "mslcd_diagnostics_master.csv",
              show_default=True, help="Master diagnostics CSV.")
@click.option("--trajectories", type=click.Path(path_type=Path),
              default=TABLES_ROOT / "mslcd_coarsening_trajectories.npz",
              show_default=True, help="Coarsening trajectories NPZ.")
@click.option("--threshold-data", type=click.Path(path_type=Path),
              default=TABLES_ROOT / "threshold_analysis_Pat02.npz",
              show_default=True, help="Threshold analysis NPZ.")
@output_options(str(FIGURES_ROOT / "report_mslcd_section" / "cross_condition"))
@click.option("--threshold-output-dir", type=click.Path(path_type=Path),
              default=FIGURES_ROOT / "report_mslcd_section" / "threshold_analysis",
              show_default=True, help="Output dir for threshold figures.")
@click.option("--latex-output", type=click.Path(path_type=Path),
              default=TABLES_ROOT / "mslcd_diagnostics_Pat02.tex",
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


# ---------------------------------------------------------------------------
# plot time-windows
# ---------------------------------------------------------------------------

@plot.command("time-windows")
@single_patient_option()
@band_phase_options()
@fc_method_option(default="imcoh_abs")
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=MSC_WINDOWS_CACHE, show_default=True)
@output_options(str(FIGURES_ROOT / "time_windows"))
@verbose_option()
@click.pass_context
def time_windows(ctx, patient, band, bands, phase, phases, fc_method,
                 cache_root, output_dir, dpi, fmt, verbose):
    """Visualize cached time-window FC matrices."""
    rpt = CliReporter.from_context(ctx, verbose=verbose)
    rpt.header(f"Time Window Visualization: {patient}")
    # TODO: integrate src/visualize_time_windows.py logic
    rpt.info("Not yet implemented in CLI. Use: python src/visualize_time_windows.py")
