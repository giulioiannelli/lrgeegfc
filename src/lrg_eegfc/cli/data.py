"""``lrg-eegfc data`` subcommands."""

from __future__ import annotations

from pathlib import Path

import click

from ._common import (
    CliReporter,
    patient_options,
    resolve_patients,
    verbose_option,
)


@click.group()
def data() -> None:
    """Inspect and manage patient data."""


# ---------------------------------------------------------------------------
# data inspect
# ---------------------------------------------------------------------------

@data.command()
@patient_options()
@click.option("--root-path", type=click.Path(path_type=Path),
              default=Path("data/stereoeeg_patients"), show_default=True)
@verbose_option()
@click.pass_context
def inspect(ctx, patients, root_path, verbose):
    """Inspect patient data and generate report."""
    from lrg_eegfc.utils.io.inspect import inspect_patient, inspect_all_patients

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    pat_list = resolve_patients(patients, dataset_root=root_path)

    rpt.header("Patient Data Inspection", Patients=", ".join(pat_list))

    if len(pat_list) == 1:
        result = inspect_patient(pat_list[0], root_path)
        if result:
            click.echo(f"\n{result.get('patient', pat_list[0])}")
            phases_info = result.get("phases", {})
            for phase_name, info in phases_info.items():
                shape = info.get("data_shape", "?")
                fs = info.get("fs", "?")
                fmt = info.get("file_format", "?")
                click.echo(f"  {phase_name}: shape={shape}, fs={fs} Hz, format={fmt}")
        else:
            click.echo(f"  No data found for {pat_list[0]}")
    else:
        results = inspect_all_patients(root_path, pat_list)
        for pat, pat_result in results.items():
            click.echo(f"\n{pat}:")
            if pat_result is None:
                click.echo("  No data found")
                continue
            phases_info = pat_result.get("phases", pat_result)
            if isinstance(phases_info, dict):
                for phase_name, info in phases_info.items():
                    if isinstance(info, dict):
                        shape = info.get("data_shape", "?")
                        fs = info.get("fs", "?")
                        click.echo(f"  {phase_name}: shape={shape}, fs={fs}")


# ---------------------------------------------------------------------------
# data stats
# ---------------------------------------------------------------------------

@data.command()
@patient_options()
@click.option("--root-path", type=click.Path(path_type=Path),
              default=Path("data/stereoeeg_patients"), show_default=True)
@verbose_option()
@click.pass_context
def stats(ctx, patients, root_path, verbose):
    """Aggregate statistics across subjects."""
    from lrg_eegfc.utils.io.inspect import inspect_all_patients

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    pat_list = resolve_patients(patients, dataset_root=root_path)

    rpt.header("Aggregate Subject Statistics")

    results = inspect_all_patients(root_path, pat_list)
    n_ok = sum(1 for v in results.values() if v is not None)

    click.echo(f"Patients scanned: {len(pat_list)}")
    click.echo(f"With data: {n_ok}")
    click.echo(f"Missing: {len(pat_list) - n_ok}")


# ---------------------------------------------------------------------------
# data compare
# ---------------------------------------------------------------------------

@data.command()
@patient_options()
@click.option("--band", required=True, help="Frequency band.")
@click.option("--phase", required=True, help="Recording phase.")
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/lrg_cache"), show_default=True)
@verbose_option()
@click.pass_context
def compare(ctx, patients, band, phase, cache_root, verbose):
    """Compare MSC vs correlation FC methods (tabular output)."""
    from lrg_eegfc.utils.metrics.compare import compare_fc_methods

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    pat_list = resolve_patients(patients)

    rpt.header("FC Method Comparison", band=band, phase=phase)

    for pat in pat_list:
        try:
            result = compare_fc_methods(pat, phase, band, cache_root=cache_root)
            if result is not None:
                click.echo(f"\n{pat}:")
                click.echo(f"  matrix_distance:    {result.matrix_distance:.4f}")
                click.echo(f"  scaled_distance:    {result.scaled_distance:.4f}")
                click.echo(f"  rank_correlation:   {result.rank_correlation:.4f}")
                click.echo(f"  tree_similarity:    {result.tree_similarity:.4f}")
                rpt.success(pat)
            else:
                rpt.skip(f"{pat}: no comparison data")
        except Exception as exc:
            rpt.fail(f"{pat}: {exc}")

    rpt.summary()
