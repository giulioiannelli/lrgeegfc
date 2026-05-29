"""``lrg-eegfc bundle`` subcommands."""

from __future__ import annotations

from pathlib import Path

import click

from lrg_eegfc.config.paths import FIGURES_ROOT

from ._common import CliReporter


@click.group()
def bundle() -> None:
    """Build publication bundles."""


@bundle.command()
@click.option("--patient", default=None, help="Single patient or omit for all.")
@click.option("--output-root", type=click.Path(path_type=Path),
              default=Path("outputs/overleaf"), show_default=True)
@click.option("--figures-root", type=click.Path(path_type=Path),
              default=FIGURES_ROOT, show_default=True)
@click.option("--png", "include_png", is_flag=True,
              help="Collect PNG figures instead of the default PDF.")
@click.option("--dry-run", is_flag=True, help="Show what would be copied.")
@click.pass_context
def overleaf(ctx, patient, output_root, figures_root, include_png, dry_run):
    """Collect figures into an Overleaf-ready bundle with manifest.

    Default collects ``*.pdf`` (canonical, vector, manuscript-grade).
    Pass ``--png`` to collect ``*.png`` instead (opt-in; PNG is not the
    canonical format per project rules).
    """
    import shutil

    rpt = CliReporter.from_context(ctx)
    suffix = "*.png" if include_png else "*.pdf"
    rpt.header(
        "Overleaf Bundle",
        output_root=output_root,
        format=suffix.lstrip("*."),
        dry_run=dry_run,
    )

    if not figures_root.exists():
        rpt.fail(f"Figures root not found: {figures_root}")
        return

    output_root.mkdir(parents=True, exist_ok=True)

    collected = 0
    for src in sorted(figures_root.rglob(suffix)):
        if patient and patient not in str(src):
            continue
        rel = src.relative_to(figures_root)
        dest = output_root / rel
        label = f"{'would copy' if dry_run else 'copy'}: {rel}"
        rpt.detail(label)
        if not dry_run:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)
        collected += 1

    rpt.info(f"\n{'Would collect' if dry_run else 'Collected'}: {collected} figures")
