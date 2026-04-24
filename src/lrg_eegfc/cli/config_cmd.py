"""``lrg-eegfc config`` subcommands."""

from __future__ import annotations

import click

from lrg_eegfc.config.paths import (
    CORR_CACHE,
    MSC_CACHE,
    LRG_CACHE,
    CLEANED_CORR_CACHE,
    FIGURES_ROOT,
    SEEG_DATAPATH,
)


@click.group("config")
def config() -> None:
    """Display project configuration and verify paths."""


@config.command()
def show() -> None:
    """Display frequency bands, phases, and default parameters."""
    from lrg_eegfc.config.const import (
        BRAIN_BANDS,
        DEFAULT_DISPARITY_ALPHA,
        DEFAULT_ECM_N_ENSEMBLE,
        DEFAULT_ECM_WEIGHT_SCALE,
        DEFAULT_FDR_Q,
        DEFAULT_N_SURROGATES,
        DEFAULT_NPERSEG,
        DEFAULT_SAMPLE_RATE,
        PHASE_LABELS,
        VALID_SPARSIFY_METHODS,
    )

    click.echo("Frequency bands")
    click.echo("-" * 40)
    for band, (lo, hi) in BRAIN_BANDS.items():
        click.echo(f"  {band:12s}  {lo:6.2f} – {hi:6.2f} Hz")

    click.echo(f"\nPhases: {', '.join(PHASE_LABELS)}")

    click.echo("\nDefaults")
    click.echo("-" * 40)
    click.echo(f"  sample_rate      {DEFAULT_SAMPLE_RATE} Hz")
    click.echo(f"  nperseg          {DEFAULT_NPERSEG}")
    click.echo(f"  n_surrogates     {DEFAULT_N_SURROGATES}")
    click.echo(f"  fdr_q            {DEFAULT_FDR_Q}")
    click.echo(f"  disparity_alpha  {DEFAULT_DISPARITY_ALPHA}")
    click.echo(f"  ecm_n_ensemble   {DEFAULT_ECM_N_ENSEMBLE}")
    click.echo(f"  ecm_weight_scale {DEFAULT_ECM_WEIGHT_SCALE}")

    click.echo(f"\nSparsify methods: {', '.join(VALID_SPARSIFY_METHODS)}")


@config.command()
def paths() -> None:
    """Show data and cache paths and verify they exist."""
    from lrg_eegfc.config.const import list_patients

    _check = lambda p: "ok" if p.exists() else "MISSING"

    dirs = {
        "Raw data":       SEEG_DATAPATH,
        "Corr cache":     CORR_CACHE,
        "MSC cache":      MSC_CACHE,
        "LRG cache":      LRG_CACHE,
        "Cleaned cache":  CLEANED_CORR_CACHE,
        "Figures":        FIGURES_ROOT,
    }

    click.echo("Data paths")
    click.echo("-" * 50)
    for label, path in dirs.items():
        click.echo(f"  {label:16s}  {path}  [{_check(path)}]")

    patients = list_patients(SEEG_DATAPATH)
    if patients:
        click.echo(f"\nPatients found ({len(patients)}): {', '.join(patients)}")
    else:
        click.echo(f"\nNo patients found under {SEEG_DATAPATH}")
