"""``lrg-eegfc data`` subcommands."""

from __future__ import annotations

import re
from pathlib import Path

import click

from lrg_eegfc.config.const import DEFAULT_SAMPLE_RATE, FS_OVERRIDES
from lrg_eegfc.config.paths import SEEG_DATAPATH, LRG_CACHE

from ._common import (
    CliReporter,
    patient_options,
    resolve_patients,
    verbose_option,
)


_FS_INFO_RE = re.compile(r"Sampling\s*Rate\s*[-:]\s*([0-9]+)", re.IGNORECASE)


def _resolve_fs(patient: str, phases_info: dict, root_path: Path) -> tuple[float, str]:
    """Pick a single fs to display per patient.

    Probes (in order): .mat ``Parameters.fs`` extracted by ``inspect_patient``,
    vendor ``.info.txt`` files under ``resting/`` and ``task/``, and finally
    the canonical config (``FS_OVERRIDES``/``DEFAULT_SAMPLE_RATE``).
    """
    for info in phases_info.values():
        if isinstance(info, dict) and info.get("fs"):
            return float(info["fs"]), "from .mat"
    for info_txt in sorted((root_path / patient).rglob("*.info.txt")):
        try:
            m = _FS_INFO_RE.search(info_txt.read_text(errors="replace"))
        except OSError:
            continue
        if m:
            return float(m.group(1)), f"from {info_txt.relative_to(root_path / patient)}"
    return float(FS_OVERRIDES.get(patient, DEFAULT_SAMPLE_RATE)), "from config (no .mat/.info.txt)"


@click.group()
def data() -> None:
    """Inspect and manage patient data."""


# ---------------------------------------------------------------------------
# data inspect
# ---------------------------------------------------------------------------

@data.command()
@patient_options()
@click.option("--root-path", type=click.Path(path_type=Path),
              default=SEEG_DATAPATH, show_default=True)
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
            phases_info = result.get("phases", {})
            fs, fs_src = _resolve_fs(pat_list[0], phases_info, root_path)
            click.echo(f"\n{result.get('patient', pat_list[0])} (fs={fs:g} Hz [{fs_src}])")
            for phase_name, info in phases_info.items():
                shape = info.get("data_shape", "?")
                fmt = info.get("file_format", "?")
                click.echo(f"  {phase_name}: shape={shape}, format={fmt}")
        else:
            click.echo(f"  No data found for {pat_list[0]}")
    else:
        results = inspect_all_patients(root_path, pat_list)
        for pat, pat_result in results.get("patients", {}).items():
            if pat_result is None:
                click.echo(f"\n{pat}:")
                click.echo("  No data found")
                continue
            phases_info = pat_result.get("phases", {}) if isinstance(pat_result, dict) else {}
            fs, fs_src = _resolve_fs(pat, phases_info, root_path)
            click.echo(f"\n{pat} (fs={fs:g} Hz [{fs_src}])")
            for phase_name, info in phases_info.items():
                if isinstance(info, dict):
                    shape = info.get("data_shape", "?")
                    fmt = info.get("file_format", "?")
                    click.echo(f"  {phase_name}: shape={shape}, format={fmt}")


# ---------------------------------------------------------------------------
# data stats
# ---------------------------------------------------------------------------

@data.command()
@patient_options()
@click.option("--root-path", type=click.Path(path_type=Path),
              default=SEEG_DATAPATH, show_default=True)
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
# data normalize
# ---------------------------------------------------------------------------

@data.command()
@patient_options()
@click.option("--root-path", type=click.Path(path_type=Path),
              default=SEEG_DATAPATH, show_default=True)
@click.option("--apply", "apply_", is_flag=True, default=False,
              help="Execute the plan. Default is dry-run (print only).")
@click.option("--report", type=click.Path(path_type=Path), default=None,
              help="Write the combined plan to this markdown file.")
@verbose_option()
@click.pass_context
def normalize(ctx, patients, root_path, apply_, report, verbose):
    """Normalize patient directory layout to the canonical scheme.

    Produces an actionable plan per patient (dry-run by default). Use
    ``--apply`` to execute. The per-patient ``provenance.md`` file is
    written as the last step of each applied plan.
    """
    from lrg_eegfc.utils.io.patient_bootstrap import (
        apply_plan,
        plan_migration,
        render_plan,
    )

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    pat_list = resolve_patients(patients, dataset_root=root_path)
    rpt.header(
        "Patient Layout Normalization",
        Patients=", ".join(pat_list),
        Mode=("APPLY" if apply_ else "dry-run"),
    )

    sections: list[str] = []
    total_actions = 0
    total_warnings = 0
    for pat in pat_list:
        patient_dir = Path(root_path) / pat
        if not patient_dir.exists():
            click.echo(f"  {pat}: directory missing, skipping")
            continue
        plan = plan_migration(patient_dir)
        total_actions += len(plan.actions)
        total_warnings += len(plan.warnings)
        text = render_plan(plan)
        sections.append(text)
        click.echo("\n" + text)

        if apply_:
            results = apply_plan(plan)
            for line in results:
                click.echo(f"  {line}")
            rpt.success(pat)
        else:
            rpt.skip(f"{pat}: dry-run ({len(plan.actions)} actions, "
                     f"{len(plan.warnings)} warnings)")

    click.echo(f"\nTotal: {total_actions} actions, {total_warnings} warnings "
               f"across {len(pat_list)} patient(s).")
    if report is not None:
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text("# Patient layout normalization plan\n\n"
                          + "\n\n".join(sections) + "\n")
        click.echo(f"Wrote plan to {report}")
    rpt.summary()


# ---------------------------------------------------------------------------
# data compare
# ---------------------------------------------------------------------------

@data.command()
@patient_options()
@click.option("--band", required=True, help="Frequency band.")
@click.option("--phase", required=True, help="Recording phase.")
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=LRG_CACHE, show_default=True)
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
