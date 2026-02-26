"""Shared CLI option factories, parameter resolvers, and output helpers.

Every reusable piece of CLI plumbing lives here so that the individual
subcommand modules (compute.py, plot.py, etc.) stay focused on delegation.
"""

from __future__ import annotations

import sys
from functools import wraps
from pathlib import Path
from typing import Callable, Optional, Sequence

import click


# ---------------------------------------------------------------------------
# Option decorator factories
# ---------------------------------------------------------------------------

def patient_options(
    required: bool = False,
    default: Optional[Sequence[str]] = None,
) -> Callable:
    """Add ``--patients`` option (multiple values, auto-discovered if empty)."""
    def decorator(f: Callable) -> Callable:
        return click.option(
            "--patients", "-p",
            multiple=True,
            default=default,
            required=required,
            help="Patient IDs to process (e.g. Pat_02 Pat_03). Default: all.",
        )(f)
    return decorator


def band_phase_options() -> Callable:
    """Add ``--bands``, ``--phases``, ``--band``, ``--phase`` options."""
    def decorator(f: Callable) -> Callable:
        f = click.option(
            "--phases", multiple=True, default=None,
            help="Phase names (default: all).",
        )(f)
        f = click.option(
            "--phase", default=None,
            help="Single phase name.",
        )(f)
        f = click.option(
            "--bands", multiple=True, default=None,
            help="Band names (default: all).",
        )(f)
        f = click.option(
            "--band", default=None,
            help="Single band name.",
        )(f)
        return f
    return decorator


def single_patient_option(required: bool = True) -> Callable:
    """Add a single ``--patient`` option (for per-patient visualisation)."""
    def decorator(f: Callable) -> Callable:
        return click.option(
            "--patient", required=required,
            help="Patient identifier (e.g. Pat_02).",
        )(f)
    return decorator


def fc_method_option(required: bool = False, default: str = "msc") -> Callable:
    """Add ``--fc-method`` option."""
    def decorator(f: Callable) -> Callable:
        return click.option(
            "--fc-method",
            type=click.Choice(["corr", "msc"]),
            required=required,
            default=default,
            show_default=True,
            help="Functional connectivity method.",
        )(f)
    return decorator


def cache_options(default_root: str = "data/corr_cache") -> Callable:
    """Add ``--cache-root`` and ``--overwrite`` options."""
    def decorator(f: Callable) -> Callable:
        f = click.option(
            "--overwrite", is_flag=True,
            help="Recompute even if cached.",
        )(f)
        f = click.option(
            "--cache-root",
            type=click.Path(path_type=Path),
            default=Path(default_root),
            show_default=True,
            help="Root directory for cache files.",
        )(f)
        return f
    return decorator


def output_options(default_dir: str = "data/figures") -> Callable:
    """Add ``--output-dir`` and ``--dpi`` options."""
    def decorator(f: Callable) -> Callable:
        f = click.option(
            "--dpi", type=int, default=150, show_default=True,
            help="Figure resolution.",
        )(f)
        f = click.option(
            "--output-dir",
            type=click.Path(path_type=Path),
            default=Path(default_dir),
            show_default=True,
            help="Directory for output figures.",
        )(f)
        return f
    return decorator


def filter_time_option() -> Callable:
    """Add ``--filter-time`` option (dev mode: limit samples)."""
    def decorator(f: Callable) -> Callable:
        return click.option(
            "--filter-time", type=int, default=None,
            help="Limit to first N samples (dev mode).",
        )(f)
    return decorator


def verbose_option() -> Callable:
    """Add local ``--verbose / -v`` flag (also inherits from root context)."""
    def decorator(f: Callable) -> Callable:
        return click.option(
            "-v", "--verbose", is_flag=True, default=False,
            help="Enable verbose output.",
        )(f)
    return decorator


def show_option() -> Callable:
    """Add ``--show / -s`` flag to open figures after saving."""
    def decorator(f: Callable) -> Callable:
        return click.option(
            "-s", "--show", "show_fig", is_flag=True, default=False,
            help="Open figure(s) after saving.",
        )(f)
    return decorator


# ---------------------------------------------------------------------------
# Parameter resolution helpers
# ---------------------------------------------------------------------------

def resolve_patients(
    patients: tuple[str, ...],
    dataset_root: Optional[Path] = None,
) -> list[str]:
    """Return an explicit patient list, auto-discovering if *patients* is empty."""
    if patients:
        return list(patients)
    from lrg_eegfc.config.const import list_patients, sEEG_DATAPATH
    root = dataset_root or sEEG_DATAPATH
    result = list_patients(root)
    if not result:
        raise click.UsageError(f"No patients found under {root}")
    return result


def resolve_bands(
    bands: tuple[str, ...] = (),
    band: Optional[str] = None,
) -> list[str]:
    """Resolve band list: ``--bands`` > ``--band`` > all."""
    if bands:
        return list(bands)
    if band:
        return [band]
    from lrg_eegfc.config.const import BRAIN_BANDS
    return list(BRAIN_BANDS.keys())


def resolve_phases(
    phases: tuple[str, ...] = (),
    phase: Optional[str] = None,
) -> list[str]:
    """Resolve phase list: ``--phases`` > ``--phase`` > all."""
    if phases:
        return list(phases)
    if phase:
        return [phase]
    from lrg_eegfc.config.const import PHASE_LABELS
    return list(PHASE_LABELS)


# ---------------------------------------------------------------------------
# Progress / output helpers
# ---------------------------------------------------------------------------

class CliReporter:
    """Structured progress reporting for CLI commands."""

    def __init__(self, *, verbose: bool = False, quiet: bool = False):
        self.verbose = verbose
        self.quiet = quiet
        self._success = 0
        self._failed = 0
        self._skipped = 0
        self._saved_paths: list[Path] = []

    # -- Output methods ---

    def header(self, title: str, **params: object) -> None:
        if self.quiet:
            return
        click.echo("=" * 70)
        click.echo(title)
        click.echo("=" * 70)
        for key, value in params.items():
            click.echo(f"{key}: {value}")
        click.echo("=" * 70)

    def info(self, msg: str) -> None:
        if not self.quiet:
            click.echo(msg)

    def detail(self, msg: str) -> None:
        if self.verbose:
            click.echo(f"  {msg}")

    def success(self, msg: str = "") -> None:
        self._success += 1
        if self.verbose and msg:
            click.echo(f"  [ok] {msg}")

    def saved(self, path: "Path", *, show: bool = False) -> None:
        """Record a saved file — always prints the path, optionally opens it."""
        self._success += 1
        self._saved_paths.append(path)
        if not self.quiet:
            click.echo(f"  saved: {path}")
        if show:
            self._open_file(path)

    def fail(self, msg: str) -> None:
        self._failed += 1
        click.echo(f"  [FAIL] {msg}", err=True)

    def skip(self, msg: str = "", *, path: "Path | None" = None,
             show: bool = False) -> None:
        self._skipped += 1
        if self.verbose and msg:
            click.echo(f"  [skip] {msg}")
        if show and path is not None and path.exists():
            self._open_file(path)

    def summary(self) -> None:
        if self.quiet:
            return
        parts = []
        if self._success:
            parts.append(f"{self._success} saved")
        if self._skipped:
            parts.append(f"{self._skipped} skipped")
        if self._failed:
            parts.append(f"{self._failed} failed")
        click.echo(f"\nDone: {', '.join(parts) or 'nothing to do'}.")

    # -- File opener ---

    @staticmethod
    def _open_file(path: Path) -> None:
        """Open a file with the platform default viewer."""
        import subprocess
        import platform

        system = platform.system()
        try:
            if system == "Darwin":
                subprocess.Popen(["open", str(path)])
            elif system == "Linux":
                subprocess.Popen(["xdg-open", str(path)])
            elif system == "Windows":
                subprocess.Popen(["start", "", str(path)], shell=True)
        except OSError:
            click.echo(f"  (could not open: {path})", err=True)

    # -- Context helpers ---

    @classmethod
    def from_context(cls, ctx: click.Context, **overrides: bool) -> "CliReporter":
        """Build a reporter from the root Click context.

        Parameters in *overrides* (e.g. ``verbose=True``) take priority over
        the root context so that per-command ``-v`` flags work.
        """
        obj = ctx.ensure_object(dict)
        verbose = overrides.get("verbose", False) or obj.get("verbose", False)
        quiet = overrides.get("quiet", False) or obj.get("quiet", False)
        return cls(verbose=verbose, quiet=quiet)
