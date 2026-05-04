"""Core constants used across :mod:`lrg_eegfc`.

The values defined in this module are intentionally lightweight and import
without touching the filesystem.  They are primarily used by the correlation
and data-loading utilities so keeping them together avoids circular
dependencies.
"""
#
from __future__ import annotations
#
from itertools import combinations
from pathlib import Path
from typing import Dict, Tuple, List
from functools import lru_cache
#
__all__ = [
    'sEEG_DATAPATH',
    'PATIENTS_LIST',
    'list_patients',
    'PHASE_LABELS',
    'PHASE_SUBDIR',
    'PATIENT_CHANNEL_DROP',
    'PATIENTS_4PHASE',
    'PHASE_SHORT_LABELS',
    'ALL_PHASE_PAIRS',
    'REST_PHASES',
    'TASK_PHASES',
    'classify_pair',
    'PARAMETER_KEYS',
    'BRAIN_BANDS_NAMES',
    'BRAIN_BANDS_FREQ',
    'BRAIN_BANDS',
    'BRAIN_BANDS_TEX_NAMES',
    'BRAIN_BAND_TEX_DICT',
    'BRAIN_BAND_LABELS',
    'DEFAULT_N_SURROGATES',
    'DEFAULT_NPERSEG',
    'DEFAULT_WELCH_SEGMENT_SEC',
    'nperseg_for_fs',
    'DEFAULT_SAMPLE_RATE',
    'DEFAULT_FDR_Q',
    'DEFAULT_DISPARITY_ALPHA',
    'DEFAULT_ECM_N_ENSEMBLE',
    'DEFAULT_ECM_WEIGHT_SCALE',
    'DEFAULT_ECM_ALPHA_MIN',
    'DEFAULT_ECM_ALPHA_MAX',
    'DEFAULT_ECM_ALPHA_TOL',
    'VALID_SPARSIFY_METHODS',
    'FC_METHODS',
    'COHERENCE_METRICS',
]
#
from lrg_eegfc.config.paths import SEEG_DATAPATH as _SEEG_DATAPATH

#: Path to raw stereoEEG patient data.  Alias for :data:`paths.SEEG_DATAPATH`.
sEEG_DATAPATH = _SEEG_DATAPATH
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
#: Canonical snake_case names matching the on-disk filenames under
#: ``resting/{rest_pre,rest_post}.mat`` and ``task/{task_learn,task_test}.mat``.
PHASE_LABELS: Tuple[str, ...] = ('rest_pre', 'task_learn', 'task_test', 'rest_post')
#: Map phase → subdirectory under the per-patient root. Mirrors the canonical
#: layout enforced by ``lrg-eegfc data normalize``.
PHASE_SUBDIR: Dict[str, str] = {
    "rest_pre":  "resting",
    "rest_post": "resting",
    "task_learn": "task",
    "task_test":  "task",
}

#: Per-patient channel drops applied at load time to reconcile heterogeneous
#: channel counts across phases.  Keys: patient id.  Values: dict mapping
#: phase name → 0-based row indices to delete from the ``Data`` matrix of
#: ``channel_labels.csv`` (row indices are relative to the *full* 116-row
#: channel_labels.csv; after dropping them the remaining rows are the
#: canonical channel set for that patient across every phase).
#:
#: A special key ``"__labels__"`` under each patient lists the label-row
#: indices to drop from ``channel_labels.csv`` and the implant CSV so that
#: per-phase Data (post-drop) lines up with the metadata row-by-row.
#:
#: Pat_10 history: vendor originally shipped resting at 113 ch and task at
#: 116 ch (2026-04-22 import). The 3 task-only extras at rows [53, 54, 55]
#: (labels ``['f  3,G2', 'c  3,G2', 'o  1,G2']``) were inferred via
#: monotonic-constrained channel-fingerprint matching on rest↔task power
#: and previously dropped at load. **Vendor re-supplied uniform 113-channel
#: data on 2026-04-25**, so the mask is no longer needed: the new files
#: ship without the 3 phantom rows. The mapping is intentionally left empty
#: so future per-patient drops can re-use this hook.
PATIENT_CHANNEL_DROP: Dict[str, Dict[str, List[int]]] = {}
#: Patients with all 4 recording phases. As of 2026-04-22 the full n=10
#: roster has canonical resting + task recordings; Pat_06 was re-completed.
PATIENTS_4PHASE: List[str] = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
#: Short labels for phases (useful for compact figure annotations).
PHASE_SHORT_LABELS: Dict[str, str] = {
    "rest_pre": "Pre", "task_learn": "TL", "task_test": "TT", "rest_post": "Post",
}
#: All ordered pairs of phases (6 combinations).
ALL_PHASE_PAIRS: List[Tuple[str, str]] = list(combinations(PHASE_LABELS, 2))
#: Resting-state phases.
REST_PHASES: frozenset = frozenset({"rest_pre", "rest_post"})
#: Task phases.
TASK_PHASES: frozenset = frozenset({"task_learn", "task_test"})


def classify_pair(p1: str, p2: str) -> str:
    """Classify a phase pair as ``'within'`` or ``'cross'``.

    A pair is *within* if both phases are resting-state or both are task;
    otherwise it is *cross*.
    """
    s = {p1, p2}
    if s <= REST_PHASES or s <= TASK_PHASES:
        return "within"
    return "cross"


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
#: Default Welch segment duration in seconds for MSC estimation.
#: The actual nperseg (in samples) is computed as ``int(fs * DEFAULT_WELCH_SEGMENT_SEC)``
#: so that all patients get the same temporal resolution regardless of sampling rate.
DEFAULT_WELCH_SEGMENT_SEC: float = 2.0
#: Default Welch segment length for MSC estimation (df = fs / nperseg).
#: Legacy constant kept for backwards compatibility. New code should use
#: :func:`nperseg_for_fs` to compute the correct value for a given sampling rate.
DEFAULT_NPERSEG: int = 4096
#: Default sampling rate when `fs` is missing from metadata
DEFAULT_SAMPLE_RATE: float = 2048.0

#: Per-patient sampling-rate overrides — only deviations from
#: :data:`DEFAULT_SAMPLE_RATE` need to be listed. This is the single
#: source of truth for the project; every script/workflow needing a
#: patient's sampling rate should look it up here rather than redefine a
#: local ``FS_MAP``. Consumers typically do
#: ``fs = FS_OVERRIDES.get(patient, DEFAULT_SAMPLE_RATE)``.
#:
#: History: Pat_03 was the only patient recorded at 1024 Hz (all others
#: at 2048 Hz); see :file:`.agents/guides/03_implementation/DATA_LAYOUT.md`
#: §6 for the full per-patient audit.
FS_OVERRIDES: Dict[str, float] = {"Pat_03": 1024.0}


def nperseg_for_fs(fs: float) -> int:
    """Return the nperseg value that gives a 2-second Welch segment at the given fs.

    This ensures comparable frequency resolution (Δf = 0.5 Hz) across patients
    recorded at different sampling rates (e.g. 1024 Hz vs 2048 Hz).
    """
    return int(fs * DEFAULT_WELCH_SEGMENT_SEC)
#: Default FDR q-value for FDR-corrected sparsification
DEFAULT_FDR_Q: float = 0.05
#: Default significance level for disparity filter
DEFAULT_DISPARITY_ALPHA: float = 0.05
#: Default number of ensemble samples for ECM sparsification
DEFAULT_ECM_N_ENSEMBLE: int = 100
#: Default weight scaling factor for ECM (float->int conversion)
DEFAULT_ECM_WEIGHT_SCALE: int = 1000
#: Default lower bound for adaptive ECM alpha search
DEFAULT_ECM_ALPHA_MIN: float = 0.01
#: Default upper bound for adaptive ECM alpha search
DEFAULT_ECM_ALPHA_MAX: float = 0.50
#: Default precision (tolerance) for adaptive ECM binary search
DEFAULT_ECM_ALPHA_TOL: float = 0.01
#: Valid sparsification methods for MSC matrices
VALID_SPARSIFY_METHODS: Tuple[str, ...] = (
    "none", "soft", "fdr", "disparity", "hybrid", "ecm", "ecm_adaptive",
)
#: Supported functional-connectivity methods across the pipeline.
#:
#: ``"imcoh"`` is the signed Nolte-2004 imaginary coherency stored on disk.
#: ``"imcoh_abs"`` and ``"imcoh_sq"`` are non-negative transforms derived by
#: the loader (``np.abs`` and ``**2`` respectively) and are what feeds LRG.
FC_METHODS: Tuple[str, ...] = ("corr", "msc", "imcoh", "imcoh_abs", "imcoh_sq")
#: Coherence-based metrics (all computed from the Welch cross-spectral density).
COHERENCE_METRICS: Tuple[str, ...] = ("msc", "imcoh", "imcoh_abs", "imcoh_sq", "wpli")

#: Display label per fc_method (LaTeX-formatted; safe to use in matplotlib
#: titles, axis labels, suptitles). Single source of truth — never hardcode
#: "ImCoh" in plot scripts; always import this dict.
FC_METHOD_DISPLAY_LABELS: Dict[str, str] = {
    "corr":      r"Pearson",
    "msc":       r"MSC",
    "imcoh":     r"$\mathrm{ImCoh}$",       # signed [-1, 1]
    "imcoh_abs": r"$|\mathrm{ImCoh}|$",     # |ImCoh|, [0, 1]
    "imcoh_sq":  r"$|\mathrm{ImCoh}|^2$",   # squared, [0, 1]
}

#: Default fc_method for Section 2 narrative figures (Ewald 2012 /
#: Bastos & Schoffelen 2016 magnitude convention). Single source of truth
#: — every Section 2 script imports this rather than hardcoding the string.
SECTION2_FC_METHOD: str = "imcoh_abs"

# ---------------------------------------------------------------------------
# Visualization defaults (FC-method-aware)
# ---------------------------------------------------------------------------

#: Spring-layout repulsion constant per band.
#: Higher = sparser band → nodes need more separation.
#: Keyed by ``(fc_method, band)`` with fallback to ``("_default", band)``.
SPRING_K: Dict[str, Dict[str, float]] = {
    "msc": {
        "delta": 5, "theta": 6, "alpha": 8,
        "beta": 18, "low_gamma": 25, "high_gamma": 30,
    },
    "imcoh": {
        "delta": 3, "theta": 4, "alpha": 5,
        "beta": 8, "low_gamma": 12, "high_gamma": 18,
    },
    # imcoh_abs / imcoh_sq are monotone transforms of signed ImCoh, so
    # the same layout knobs apply; kept as independent entries so future
    # tuning per transform does not require touching the dispatcher.
    "imcoh_abs": {
        "delta": 3, "theta": 4, "alpha": 5,
        "beta": 8, "low_gamma": 12, "high_gamma": 18,
    },
    "imcoh_sq": {
        "delta": 3, "theta": 4, "alpha": 5,
        "beta": 8, "low_gamma": 12, "high_gamma": 18,
    },
    "corr": {
        "delta": 5, "theta": 6, "alpha": 8,
        "beta": 18, "low_gamma": 25, "high_gamma": 30,
    },
}

#: Edge-weight scaling mode per FC method.
#: ``"rank"`` = rank-based (uniform quantile), good for heavy-tailed MSC.
#: ``"linear"`` = linear min-max, good for ImCoh (naturally less skewed).
EDGE_WEIGHT_SCALING: Dict[str, str] = {
    "msc": "rank",
    "imcoh": "linear",
    "imcoh_abs": "linear",
    "imcoh_sq": "linear",
    "corr": "rank",
}

#: Edge visibility (alpha) range ``(min_alpha, max_alpha)`` per FC method.
EDGE_ALPHA_RANGE: Dict[str, Tuple[float, float]] = {
    "msc": (0.05, 0.9),
    "imcoh": (0.03, 0.95),
    "imcoh_abs": (0.03, 0.95),
    "imcoh_sq": (0.03, 0.95),
    "corr": (0.05, 0.9),
}

#: Edge width range ``(min_width, max_width)`` per FC method.
EDGE_WIDTH_RANGE: Dict[str, Tuple[float, float]] = {
    "msc": (0.3, 4.0),
    "imcoh": (0.15, 4.0),
    "imcoh_abs": (0.15, 4.0),
    "imcoh_sq": (0.15, 4.0),
    "corr": (0.3, 4.0),
}


def spring_k_for(fc_method: str, band: str, n_nodes: int) -> float:
    """Return the spring-layout *k* parameter scaled by network size.

    ``k = K_base / sqrt(N)`` where ``K_base`` comes from :data:`SPRING_K`.
    """
    import math
    method_k = SPRING_K.get(fc_method, SPRING_K["msc"])
    k_base = method_k.get(band, 10)
    return k_base / math.sqrt(n_nodes)


# ---------------------------------------------------------------------------
# Edge-drawing defaults per FC method (for production network figures).
# The ImCoh values were picked from the parameter sweep in
# `scripts/10_notes_imcoh/fig_helper_three_layouts.py` and produce a
# smooth rank-based gradient where only the strongest few edges stand out.
# ---------------------------------------------------------------------------

#: Rank-based gamma exponent: width = wmin + ((rank+1)/N)^γ · (wmax - wmin).
EDGE_RANK_GAMMA: Dict[str, float] = {
    "msc":       2.0,
    "imcoh":     30.0,
    "imcoh_abs": 30.0,
    "imcoh_sq":  30.0,
    "corr":      2.0,
}

#: (wmin, wmax) for edge widths.
EDGE_RANK_WIDTH: Dict[str, Tuple[float, float]] = {
    "msc":       (0.15, 4.0),
    "imcoh":     (0.0, 4.0),
    "imcoh_abs": (0.0, 4.0),
    "imcoh_sq":  (0.0, 4.0),
    "corr":      (0.15, 4.0),
}

#: (amin, amax) for edge alphas.
EDGE_RANK_ALPHA: Dict[str, Tuple[float, float]] = {
    "msc":       (0.06, 0.9),
    "imcoh":     (0.05, 0.85),
    "imcoh_abs": (0.05, 0.85),
    "imcoh_sq":  (0.05, 0.85),
    "corr":      (0.06, 0.9),
}

#: Preferred layout algorithm per FC method.
#: ``"spring"`` works for MSC (heavy-tailed, clear modules).
#: ``"kk_ultrametric"`` uses the LRG ultrametric as the pairwise distance
#: for Kamada-Kawai — this is the only recipe that reveals modular
#: structure in ImCoh's near-uniform weight distribution.
LAYOUT_METHOD: Dict[str, str] = {
    # All FC methods use classical MDS (PCoA) on the LRG ultrametric
    # distance matrix D_ij. Direct eigendecomposition of the
    # double-centered -D²/2 matrix → first 2 eigenvectors give the 2D
    # embedding. Deterministic, no iterative optimisation, distances
    # preserved in least-squares sense. Same physics for every method,
    # so layouts are directly comparable across MSC / |ImCoh| / corr.
    # Per-method choice — INTENTIONAL ASYMMETRY documenting MSC's
    # sensitivity to volume-conduction artefacts:
    #
    #  - MSC → "spring" (Fruchterman-Reingold on raw |FC| weights).
    #    This makes MSC's same-probe-inflation visible: nodes cluster
    #    by electrode shaft because the same-probe edges dominate the
    #    spring forces. The "chaotic" / "shaft-driven" appearance is
    #    not a bug — it's the point. Section 2 narrative uses this as
    #    direct evidence of why MSC is the wrong estimator for our data.
    #
    #  - ImCoh / corr → "mds_lrg_continuous" (classical MDS on the
    #    continuous LRG heat-kernel transition distance Trho at
    #    τ = 1/λ_max, winsorised at 99th percentile). Reveals real
    #    functional structure once volume-conduction is removed.
    "msc":       "spring",
    "imcoh":     "mds_lrg_continuous",
    "imcoh_abs": "mds_lrg_continuous",
    "imcoh_sq":  "mds_lrg_continuous",
    "corr":      "mds_lrg_continuous",
}

#: Default seed for reproducibility of all layouts.
LAYOUT_SEED_DEFAULT: int = 22

#: Default iterations for spring layout.
SPRING_ITER_DEFAULT: int = 10000

#: Default node marker size for network figures.
NODE_SIZE_DEFAULT: float = 40.0
