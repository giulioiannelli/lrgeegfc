"""``lrg-eegfc bundle`` subcommands."""

from __future__ import annotations

from pathlib import Path

import click

from ._common import CliReporter


@click.group()
def bundle() -> None:
    """Build publication bundles."""


@bundle.command()
@click.option("--patient", default=None, help="Single patient or omit for all.")
@click.option("--output-root", type=click.Path(path_type=Path),
              default=Path("outputs/overleaf"), show_default=True)
@click.option("--figures-root", type=click.Path(path_type=Path),
              default=Path("data/figures"), show_default=True)
@click.option("--dry-run", is_flag=True, help="Show what would be copied.")
@click.pass_context
def overleaf(ctx, patient, output_root, figures_root, dry_run):
    """Collect figures into an Overleaf-ready bundle with manifest."""
    import shutil

    rpt = CliReporter.from_context(ctx)
    rpt.header("Overleaf Bundle", output_root=output_root, dry_run=dry_run)

    if not figures_root.exists():
        rpt.fail(f"Figures root not found: {figures_root}")
        return

    output_root.mkdir(parents=True, exist_ok=True)

    collected = 0
    for png in sorted(figures_root.rglob("*.png")):
        if patient and patient not in str(png):
            continue
        rel = png.relative_to(figures_root)
        dest = output_root / rel
        label = f"{'would copy' if dry_run else 'copy'}: {rel}"
        rpt.detail(label)
        if not dry_run:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(png, dest)
        collected += 1

    rpt.info(f"\n{'Would collect' if dry_run else 'Collected'}: {collected} figures")
