"""Core constants used across :mod:`lrg_eegfc`.

The values defined in this module are intentionally lightweight and import
without touching the filesystem.  They are primarily used by the correlation
and data-loading utilities so keeping them together avoids circular
dependencies.
"""
#
from __future__ import annotations
#
from pathlib import Path
from typing import Dict, Tuple, List
from functools import lru_cache
#
__all__ = [
    'sEEG_DATAPATH',
    'PATIENTS_LIST',
    'list_patients',
    'PHASE_LABELS',
    'PARAMETER_KEYS',
    'BRAIN_BANDS_NAMES',
    'BRAIN_BANDS_FREQ',
    'BRAIN_BANDS',
    'BRAIN_BANDS_TEX_NAMES',
    'BRAIN_BAND_TEX_DICT',
    'BRAIN_BAND_LABELS',
    'DEFAULT_N_SURROGATES',
    'DEFAULT_NPERSEG',
    'DEFAULT_SAMPLE_RATE',
    'DEFAULT_FDR_Q',
    'DEFAULT_DISPARITY_ALPHA',
    'DEFAULT_ECM_N_ENSEMBLE',
    'DEFAULT_ECM_WEIGHT_SCALE',
    'VALID_SPARSIFY_METHODS',
]
#
sEEG_DATAPATH = Path('data') / 'stereoeeg_patients'
#
@lru_cache(maxsize=None)
def list_patients(dataset_root: Path = sEEG_DATAPATH) -> List[str]:
    """Return available patient IDs without scanning on import."""
    if not dataset_root.exists():
        return []
    return sorted(
        p.name
        for p in dataset_root.iterdir()
        if p.is_dir() and p.name.startswith("Pat_")
    )


class _LazyPatientsList(list):
    def __init__(self):
        super().__init__()
        self._loaded = False

    def _load(self) -> None:
        if not self._loaded:
            self[:] = list_patients()
            self._loaded = True

    def __iter__(self):
        self._load()
        return super().__iter__()

    def __len__(self) -> int:
        self._load()
        return super().__len__()

    def __getitem__(self, index):
        self._load()
        return super().__getitem__(index)

    def __repr__(self) -> str:
        self._load()
        return super().__repr__()


#: Patients available in the LRG EEG-FC dataset (lazy-loaded).
PATIENTS_LIST: List[str] = _LazyPatientsList()
#: Recording phases expected in the publicly shared SEEG datasets.
PHASE_LABELS: Tuple[str, ...] = ('rsPre', 'taskLearn', 'taskTest', 'rsPost')
#: Keys typically embedded inside the ``Parameters`` struct of the ``.mat``
#: files distributed with the LRG datasets.
PARAMETER_KEYS: Tuple[str, ...] = (
    'fs',
    'fcutHigh',
    'fcutLow',
    'filter_order',
    'NotchFilter',
    'DataDimensions',
)
#: Brain bands keys
BRAIN_BANDS_NAMES: List[str] = [
    'delta',
    'theta',
    'alpha',
    'beta',
    'low_gamma',
    'high_gamma'
]
#: Brain bands frequency tuples
BRAIN_BANDS_FREQ: List[Tuple[float, float]] = [
    (0.53, 4.0),
    (4.0, 8.0),
    (8.0, 13.0),
    (13.0, 30.0),
    (30.0, 80.0),
    (80.0, 300.0),
]
#: Canonical EEG frequency bands expressed as ``(low, high)`` Hz pairs.
BRAIN_BANDS: Dict[str, Tuple[float, float]] = {
    band: freq for band, freq in zip(BRAIN_BANDS_NAMES, BRAIN_BANDS_FREQ)
}
#: LaTeX-friendly labels for each EEG band. Useful when generating plots.
BRAIN_BANDS_TEX_NAMES: List[str] = [
    r"$\delta$",
    r"$\theta$",
    r"$\alpha$",
    r"$\beta$",
    r"$\gamma_{\mathrm{l}}$",
    r"$\gamma_{\mathrm{h}}$",
]
BRAIN_BAND_TEX_DICT: Dict[str, str] = {
    band: tex_label
    for band, tex_label in zip(BRAIN_BANDS_NAMES, BRAIN_BANDS_TEX_NAMES)
}
#: Backwards-compatible alias for LaTeX-friendly labels.
BRAIN_BAND_LABELS: Dict[str, str] = BRAIN_BAND_TEX_DICT
#: Default number of surrogates for coherence-based FC null model estimation
DEFAULT_N_SURROGATES: int = 200
#: Default Welch segment length for MSC estimation (df = fs / nperseg)
DEFAULT_NPERSEG: int = 4096
#: Default sampling rate when `fs` is missing from metadata
DEFAULT_SAMPLE_RATE: float = 2048.0
#: Default FDR q-value for FDR-corrected sparsification
DEFAULT_FDR_Q: float = 0.05
#: Default significance level for disparity filter
DEFAULT_DISPARITY_ALPHA: float = 0.05
#: Default number of ensemble samples for ECM sparsification
DEFAULT_ECM_N_ENSEMBLE: int = 100
#: Default weight scaling factor for ECM (float→int conversion)
DEFAULT_ECM_WEIGHT_SCALE: int = 1000
#: Valid sparsification methods for MSC matrices
VALID_SPARSIFY_METHODS: Tuple[str, ...] = (
    "none", "soft", "fdr", "disparity", "hybrid", "ecm",
)
