"""Shared helpers for standalone analysis scripts.

This module eliminates the most common boilerplate duplicated across 120+
scripts in ``scripts/``.  Import it at the top of any script to get:

* :func:`setup_script_env` — configure ``matplotlib`` for headless rendering,
  add ``src/`` to ``sys.path``, and return the project root.
* :func:`iter_patient_band_phase` — iterate over (patient, band, phase) with
  standard filtering (e.g. exclude Pat_06 for cross-phase work).
* :func:`save_figure` — save a figure + companion ``.md`` description, create
  output directories, and close the figure.
* :func:`write_report` — write collected text lines to a Markdown file.
* :func:`load_all_lrg` / :func:`load_all_fc` — bulk-load cached results.

Example usage in a script::

    from lrg_eegfc.utils.scripting import setup_script_env, iter_patient_band_phase, save_figure
    ROOT = setup_script_env()

    for patient, band, phase in iter_patient_band_phase(require_all_phases=True):
        ...
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Iterator, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Script environment setup
# ---------------------------------------------------------------------------

def setup_script_env() -> Path:
    """Configure ``matplotlib`` for headless rendering, find the project root,
    and add ``src/`` to ``sys.path``.

    The project root is detected by walking upward from the calling script's
    directory until a ``pyproject.toml`` is found.  This makes the function
    independent of directory depth.

    Returns
    -------
    Path
        Absolute path to the project root directory.
    """
    import matplotlib
    matplotlib.use("Agg")

    # Walk upward to find project root (contains pyproject.toml)
    import inspect
    caller_file = inspect.stack()[1].filename
    candidate = Path(caller_file).resolve().parent
    for _ in range(10):
        if (candidate / "pyproject.toml").exists():
            src = candidate / "src"
            if str(src) not in sys.path:
                sys.path.insert(0, str(src))
            return candidate
        candidate = candidate.parent

    # Fallback: assume two levels up from script
    root = Path(caller_file).resolve().parents[2]
    src = root / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))
    return root


# ---------------------------------------------------------------------------
# Patient / band / phase iteration
# ---------------------------------------------------------------------------

def iter_patient_band_phase(
    patients: Optional[Sequence[str]] = None,
    bands: Optional[Sequence[str]] = None,
    phases: Optional[Sequence[str]] = None,
    skip_patients: Optional[Sequence[str]] = None,
    require_all_phases: bool = False,
) -> Iterator[Tuple[str, str, str]]:
    """Yield ``(patient, band, phase)`` tuples with standard filtering.

    Parameters
    ----------
    patients : sequence of str, optional
        Patient IDs to iterate over.  Defaults to
        :data:`~lrg_eegfc.config.const.PATIENTS_4PHASE`.
    bands : sequence of str, optional
        Band names.  Defaults to
        :data:`~lrg_eegfc.config.const.BRAIN_BANDS_NAMES`.
    phases : sequence of str, optional
        Phase labels.  Defaults to
        :data:`~lrg_eegfc.config.const.PHASE_LABELS`.
    skip_patients : sequence of str, optional
        Patient IDs to skip.
    require_all_phases : bool
        If ``True``, automatically exclude patients known to have fewer than
        4 phases (currently Pat_06).

    Yields
    ------
    tuple of (str, str, str)
        ``(patient, band, phase)`` combinations.
    """
    from lrg_eegfc.config.const import (
        BRAIN_BANDS_NAMES,
        PATIENTS_4PHASE,
        PHASE_LABELS,
    )

    patients = list(patients or PATIENTS_4PHASE)
    bands = list(bands or BRAIN_BANDS_NAMES)
    phases = list(phases or PHASE_LABELS)

    if require_all_phases:
        patients = [p for p in patients if p != "Pat_06"]

    if skip_patients:
        skip = set(skip_patients)
        patients = [p for p in patients if p not in skip]

    for patient in patients:
        for band in bands:
            for phase in phases:
                yield patient, band, phase


# ---------------------------------------------------------------------------
# Figure saving
# ---------------------------------------------------------------------------

def save_figure(
    fig,
    path: os.PathLike,
    *,
    what: str = "",
    proves: str = "",
    how_to_read: str = "",
    dpi: int = 300,
    fmt: str = "pdf",
    close: bool = True,
) -> Path:
    """Save a figure and optionally write a companion ``.md`` description.

    Creates parent directories automatically and closes the figure by default.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    path : path-like
        Output file path.  If it has no suffix, *fmt* is appended.
    what, proves, how_to_read : str
        When any of these is non-empty, a companion Markdown file is written
        next to the figure with structured metadata for future agents.
    dpi : int
        Resolution for raster formats.
    fmt : str
        Default format when *path* has no suffix.
    close : bool
        Close the figure after saving (default ``True``).

    Returns
    -------
    Path
        The resolved path where the figure was saved.
    """
    import matplotlib.pyplot as plt

    path = Path(path)
    if not path.suffix:
        path = path.with_suffix(f".{fmt}")
    path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(path, bbox_inches="tight", dpi=dpi)
    if close:
        plt.close(fig)
    print(f"  Saved: {path}")

    # Write companion .md when metadata is provided
    if what or proves or how_to_read:
        md_path = path.with_suffix(".md")
        sections = [f"# {path.stem}\n"]
        if what:
            sections.append(f"## What the figure shows\n{what}\n")
        if proves:
            sections.append(f"## What it proves\n{proves}\n")
        if how_to_read:
            sections.append(f"## How to read it\n{how_to_read}\n")
        md_path.write_text("\n".join(sections))

    return path


# ---------------------------------------------------------------------------
# Report writing
# ---------------------------------------------------------------------------

def write_report(
    lines: Sequence[str],
    path: os.PathLike,
    *,
    title: str = "",
) -> Path:
    """Write collected text lines to a Markdown file.

    Parameters
    ----------
    lines : sequence of str
        Lines to write (without trailing newlines — they are added).
    path : path-like
        Output file path.
    title : str
        Optional title added as ``# title`` header.

    Returns
    -------
    Path
        The resolved path where the report was saved.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    parts = []
    if title:
        parts.append(f"# {title}\n")
    parts.extend(lines)

    path.write_text("\n".join(parts) + "\n")
    print(f"  Report saved: {path}")
    return path


# ---------------------------------------------------------------------------
# Bulk data loading
# ---------------------------------------------------------------------------

def load_all_lrg(
    patients: Optional[Sequence[str]] = None,
    bands: Optional[Sequence[str]] = None,
    phases: Optional[Sequence[str]] = None,
    fc_method: str = "imcoh_abs",
    cache_root: Optional[Path] = None,
    verbose: bool = True,
) -> dict:
    """Load all cached LRG results.

    Returns a dict keyed by ``(patient, phase, band)`` → ``LRGResult``.

    Parameters
    ----------
    patients, bands, phases : sequence of str, optional
        Override the default iteration scope (see :func:`iter_patient_band_phase`).
    fc_method : str
        FC method used for LRG computation. Default ``"imcoh_abs"`` (current
        era). ``"msc"``, ``"corr"``, ``"imcoh"``, ``"imcoh_sq"`` accessible
        for diagnostics.
    cache_root : Path, optional
        Override the LRG cache directory.
    verbose : bool
        Print progress summary.
    """
    from lrg_eegfc.config.paths import lrg_cache_for
    from lrg_eegfc.workflow.lrg import load_lrg_result

    cache_root = cache_root or lrg_cache_for(fc_method)
    data: dict = {}
    n_loaded = n_missing = 0

    for patient, band, phase in iter_patient_band_phase(patients=patients, bands=bands, phases=phases):
        res = load_lrg_result(patient, phase, band, fc_method, cache_root=cache_root)
        if res is not None:
            data[(patient, phase, band)] = res
            n_loaded += 1
        else:
            n_missing += 1
            if verbose:
                print(f"  MISS: {patient} {phase} {band}")

    if verbose:
        print(f"  LRG loaded: {n_loaded}, missing: {n_missing}")
    return data


def load_all_fc(
    patients: Optional[Sequence[str]] = None,
    bands: Optional[Sequence[str]] = None,
    phases: Optional[Sequence[str]] = None,
    fc_method: str = "imcoh_abs",
    cache_root: Optional[Path] = None,
    sparsify: str = "none",
    n_surrogates: int = 0,
    verbose: bool = True,
) -> dict:
    """Load all cached FC matrices.

    Returns a dict keyed by ``(patient, phase, band)`` → ``np.ndarray``.

    Parameters
    ----------
    patients, bands, phases : sequence of str, optional
        Override the default iteration scope.
    fc_method : str
        Default ``"imcoh_abs"`` (current era). ``"msc"``, ``"corr"``,
        ``"imcoh"``, ``"imcoh_sq"`` accessible for diagnostics.
    cache_root : Path, optional
        Override the cache directory.
    sparsify : str
        Sparsification method (default ``"none"`` for raw matrices).
    n_surrogates : int
        Number of surrogates (0 for no surrogate correction).
    verbose : bool
        Print progress summary.
    """
    from lrg_eegfc.workflow.fc import load_fc_matrix as _load_fc

    def _load(pat, phase, band):
        return _load_fc(
            pat, phase, band, fc_method,
            cache_root=cache_root,
            sparsify=sparsify,
            n_surrogates=n_surrogates,
        )

    data: dict = {}
    n_loaded = n_missing = 0

    for patient, band, phase in iter_patient_band_phase(patients=patients, bands=bands, phases=phases):
        W = _load(patient, phase, band)
        if W is not None:
            data[(patient, phase, band)] = W
            n_loaded += 1
        else:
            n_missing += 1
            if verbose:
                print(f"  MISS FC: {patient} {phase} {band}")

    if verbose:
        print(f"  FC loaded: {n_loaded}, missing: {n_missing}")
    return data
