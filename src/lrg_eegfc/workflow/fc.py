"""Unified FC matrix dispatcher.

Provides a single entry point — :func:`load_fc_matrix` — that routes to the
correct loader and cache directory based on ``fc_method``.  This eliminates
per-method ``if``/``elif`` conditionals scattered throughout visualization and
analysis code.

Supported ``fc_method`` shortcuts
---------------------------------
``"corr"``
    Pearson correlation (``workflow.corr``).
``"msc"``
    Magnitude-squared coherence in ``[0, 1]`` (``workflow.msc``).
``"imcoh"``
    Signed Nolte-2004 imaginary coherency in ``[-1, 1]`` (loads from
    :data:`IMCOH_CACHE`).  Skew-symmetric; feeding this directly to LRG is
    invalid — use ``"imcoh_abs"`` or ``"imcoh_sq"`` instead.
``"imcoh_abs"``
    ``|ImCoh|`` in ``[0, 1]``, derived at load time via ``np.abs`` of the
    signed cache (no separate on-disk cache).  Ewald 2012 / Bastos-Schoffelen
    2016 convention.  This is the default for LRG.
``"imcoh_sq"``
    ``|ImCoh|^2`` in ``[0, 1]``, derived at load time via squaring the
    signed cache.

The mapping to the legacy triple ``(family, coh_method, coh_transform)`` is
exposed as :data:`FC_METHOD_SHORTCUTS` so future code can migrate to the
two-level API without touching callers.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

from lrg_eegfc.config.const import FC_METHODS, DEFAULT_NPERSEG, nperseg_for_fs
from lrg_eegfc.config.paths import (
    CORR_CACHE,
    IMCOH_CACHE,
    MSC_CACHE,
    fc_cache_for,
)

__all__ = ["load_fc_matrix", "FC_METHOD_SHORTCUTS"]

# ---------------------------------------------------------------------------
# Shortcut → (family, coh_method, coh_transform) routing table
# ---------------------------------------------------------------------------
#: Explicit mapping from the single-string ``fc_method`` shortcut to the
#: underlying (family, coherence method, magnitude transform) triple.
#: ``None`` entries mean the axis does not apply (e.g. correlation has no
#: coherence method or transform).
FC_METHOD_SHORTCUTS: dict[str, tuple[str, Optional[str], Optional[str]]] = {
    "corr":      ("correlation", None,    None),
    "msc":       ("coherence",   "msc",   "identity"),
    "imcoh":     ("coherence",   "imcoh", "identity"),
    "imcoh_abs": ("coherence",   "imcoh", "abs"),
    "imcoh_sq":  ("coherence",   "imcoh", "sq"),
}

# ---------------------------------------------------------------------------
# Patient sampling-rate map — imported from the single source of truth in
# config/const.py. Do NOT define a local FS_MAP here.
# ---------------------------------------------------------------------------
from ..config.const import FS_OVERRIDES, DEFAULT_SAMPLE_RATE


def _default_nperseg(patient: str) -> int:
    """Return the correct nperseg for *patient* (handles Pat_03 at 1024 Hz)."""
    fs = FS_OVERRIDES.get(patient, DEFAULT_SAMPLE_RATE)
    return nperseg_for_fs(fs)


def _apply_coh_transform(raw: np.ndarray, coh_transform: str) -> np.ndarray:
    """Apply the requested magnitude transform to a signed coherence matrix."""
    if coh_transform == "identity":
        return raw
    if coh_transform == "abs":
        return np.abs(raw)
    if coh_transform == "sq":
        return raw ** 2
    raise ValueError(
        f"Unknown coh_transform {coh_transform!r}; expected one of "
        "'identity', 'abs', 'sq'."
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_fc_matrix(
    patient: str,
    phase: str,
    band: str,
    fc_method: str,
    *,
    cache_root: Optional[Path] = None,
    filter_time: Optional[int] = None,
    **kwargs,
) -> Optional[np.ndarray]:
    """Load a cached FC matrix for any supported method.

    Parameters
    ----------
    patient : str
        Patient identifier (e.g. ``"Pat_02"``).
    phase : str
        Recording phase (e.g. ``"rest_pre"``).
    band : str
        Frequency band name (must exist in ``BRAIN_BANDS``).
    fc_method : str
        One of the keys of :data:`FC_METHOD_SHORTCUTS`.
    cache_root : Path, optional
        Override the default cache directory for the chosen method.
    filter_time : int, optional
        Dev-cache limit (passed to downstream loaders).
    **kwargs
        Additional keyword arguments forwarded to the method-specific loader
        (e.g. ``filter_type`` for corr, ``sparsify`` / ``n_surrogates`` /
        ``nperseg`` for msc).

    Returns
    -------
    np.ndarray or None
        The (N, N) FC matrix, or ``None`` if not cached.  For ``imcoh_abs``
        / ``imcoh_sq`` the transform is applied to the raw signed cache on
        the fly.

    Raises
    ------
    ValueError
        If *fc_method* is not recognised.
    """
    if fc_method not in FC_METHOD_SHORTCUTS:
        raise ValueError(
            f"Unknown fc_method {fc_method!r}, expected one of "
            f"{list(FC_METHOD_SHORTCUTS)}"
        )
    family, coh_method, coh_transform = FC_METHOD_SHORTCUTS[fc_method]

    if family == "correlation":
        return _load_corr(patient, phase, band, cache_root, filter_time, **kwargs)
    if family == "coherence" and coh_method == "msc":
        return _load_msc(patient, phase, band, cache_root, filter_time, **kwargs)
    if family == "coherence" and coh_method == "imcoh":
        return _load_imcoh(
            patient,
            phase,
            band,
            cache_root,
            filter_time,
            coh_transform=coh_transform,
            **kwargs,
        )
    raise ValueError(  # pragma: no cover — defensive
        f"Unroutable fc_method {fc_method!r} -> {(family, coh_method, coh_transform)}"
    )


# ---------------------------------------------------------------------------
# Private loaders
# ---------------------------------------------------------------------------


def _load_corr(
    patient: str,
    phase: str,
    band: str,
    cache_root: Optional[Path],
    filter_time: Optional[int],
    **kwargs,
) -> Optional[np.ndarray]:
    from lrg_eegfc.workflow.corr import load_corr_matrix

    root = cache_root or CORR_CACHE
    kw = {"filter_type": "abs", "zero_diagonal": True}
    kw.update(kwargs)
    return load_corr_matrix(
        patient, phase, band, cache_root=root, filter_time=filter_time, **kw
    )


def _load_msc(
    patient: str,
    phase: str,
    band: str,
    cache_root: Optional[Path],
    filter_time: Optional[int],
    **kwargs,
) -> Optional[np.ndarray]:
    from lrg_eegfc.workflow.msc import load_msc_matrix

    root = cache_root or MSC_CACHE
    kw = {"sparsify": "none", "n_surrogates": 0}
    kw.update(kwargs)
    return load_msc_matrix(
        patient, phase, band, cache_root=root, filter_time=filter_time, **kw
    )


def _load_imcoh(
    patient: str,
    phase: str,
    band: str,
    cache_root: Optional[Path],
    filter_time: Optional[int],
    *,
    coh_transform: str = "identity",
    nperseg: Optional[int] = None,
    **_kwargs,
) -> Optional[np.ndarray]:
    """Load a freq-resolved signed ImCoh cache and band-average on the fly.

    The on-disk cache stores the per-band **freq-resolved signed** ImCoh
    with shape ``(N, N, F_band)`` so that the three literature conventions
    can be derived correctly via the appropriate order of operations:

    - ``identity``: ``mean(signed, axis=F_band)`` — sign-preserving average
      (Nolte 2004). Can cancel positive/negative phase-lag contributions.
    - ``abs``:      ``mean(|signed|, axis=F_band)`` — connectivity-strength
      magnitude (Ewald 2012 / Bastos-Schoffelen 2016).
    - ``sq``:       ``mean(signed**2, axis=F_band)`` — squared imaginary
      coherence (Ewald 2012). Matches the pre-reset archive to 1e-10.

    Order matters: ``abs(mean(signed)) ≠ mean(|signed|)`` and
    ``mean(signed)**2 ≠ mean(signed**2)`` by Jensen's inequality.

    File pattern::

        {cache_root}/{patient}/{band}_{phase}_imcoh_freqresolved_nperseg-{N}.npy
    """
    root = cache_root or IMCOH_CACHE
    if nperseg is None:
        nperseg = _default_nperseg(patient)

    cache_dir = root / patient
    fname = f"{band}_{phase}_imcoh_freqresolved_nperseg-{nperseg}.npy"
    fpath = cache_dir / fname

    if not fpath.exists():
        return None

    freq_resolved = np.load(fpath).astype(np.float64, copy=False)  # (N, N, F)
    # Apply transform per-frequency-bin first, then average over band.
    if coh_transform == "identity":
        return freq_resolved.mean(axis=-1)
    if coh_transform == "abs":
        return np.abs(freq_resolved).mean(axis=-1)
    if coh_transform == "sq":
        return (freq_resolved ** 2).mean(axis=-1)
    raise ValueError(
        f"Unknown coh_transform {coh_transform!r}; expected one of "
        "'identity', 'abs', 'sq'."
    )
