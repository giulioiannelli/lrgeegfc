"""Centralized data path defaults for the lrg_eegfc toolkit.

All default paths are relative to the project root (the directory containing
the ``data/`` folder).  Every module that needs a data path should import from
here rather than constructing its own ``Path("data/...")`` literal.

Override the base ``data/`` directory by setting the environment variable
``LRGEEGFC_DATA_ROOT`` (useful for HPC or alternative dataset locations).
"""

from __future__ import annotations

import os
from pathlib import Path

__all__ = [
    # roots
    "DATA_ROOT",
    "RAW_ROOT",
    "CACHE_ROOT",
    "REPORTS_ROOT",
    "OUTPUTS_ROOT",
    # raw data
    "SEEG_DATAPATH",
    # cache
    "CORR_CACHE",
    "CORR_DEV_CACHE",
    "MSC_CACHE",
    "MSC_DEV_CACHE",
    "LRG_CACHE",
    "LRG_DEV_CACHE",
    "CLEANED_CORR_CACHE",
    "LRG_CREMA_CACHE",
    "FC_FIG_CACHE",
    "CORR_WINDOWS_CACHE",
    "CORR_WINDOWS_DEV_CACHE",
    "MSC_WINDOWS_CACHE",
    "MSC_WINDOWS_DEV_CACHE",
    "LRG_WINDOWS_CACHE",
    "LRG_WINDOWS_DEV_CACHE",
    "METRIC_CONCORDANCE_CACHE",
    "SURROGATE_VALIDATION_CACHE",
    # imcoh (imaginary coherence — volume-conduction immune)
    "IMCOH_CACHE",
    "IMCOH_LRG_CACHE",
    "IMCOH_HALVES_CACHE",
    # experimental
    "BIPOLAR_CACHE",
    "RESCALED_CACHE",
    # outputs
    "FIGURES_ROOT",
    "TABLES_ROOT",
    # routing helpers
    "fc_cache_for",
    "lrg_cache_for",
    "lrg_filename",
]

# ---------------------------------------------------------------------------
# Base data root — overridable via environment
# ---------------------------------------------------------------------------
DATA_ROOT = Path(os.environ.get("LRGEEGFC_DATA_ROOT", "data"))

# ---------------------------------------------------------------------------
# Section 1: Raw data
# ---------------------------------------------------------------------------
RAW_ROOT = DATA_ROOT / "raw"
SEEG_DATAPATH = RAW_ROOT / "stereoeeg_patients"

# ---------------------------------------------------------------------------
# Section 2: Cache (computed results)
# ---------------------------------------------------------------------------
CACHE_ROOT = DATA_ROOT / "cache"

CORR_CACHE = CACHE_ROOT / "corr"
CORR_DEV_CACHE = CACHE_ROOT / "corr_dev"
MSC_CACHE = CACHE_ROOT / "msc"
MSC_DEV_CACHE = CACHE_ROOT / "msc_dev"
LRG_CACHE = CACHE_ROOT / "lrg"
LRG_DEV_CACHE = CACHE_ROOT / "lrg_dev"
CLEANED_CORR_CACHE = CACHE_ROOT / "cleaned_corr"
LRG_CREMA_CACHE = CACHE_ROOT / "lrg_crema"
FC_FIG_CACHE = CACHE_ROOT / "fc_fig"

CORR_WINDOWS_CACHE = CACHE_ROOT / "corr_windows"
CORR_WINDOWS_DEV_CACHE = CACHE_ROOT / "corr_windows_dev"
MSC_WINDOWS_CACHE = CACHE_ROOT / "msc_windows"
MSC_WINDOWS_DEV_CACHE = CACHE_ROOT / "msc_windows_dev"
LRG_WINDOWS_CACHE = CACHE_ROOT / "lrg_windows"
LRG_WINDOWS_DEV_CACHE = CACHE_ROOT / "lrg_windows_dev"

METRIC_CONCORDANCE_CACHE = CACHE_ROOT / "metric_concordance"
SURROGATE_VALIDATION_CACHE = CACHE_ROOT / "surrogate_validation"

# ImCoh (Imaginary Coherence — volume-conduction immune FC)
IMCOH_CACHE = CACHE_ROOT / "imcoh"
IMCOH_LRG_CACHE = CACHE_ROOT / "imcoh_lrg"
#: Split halves of a continuous phase, one file per (band, phase, half,
#: transform). These supply the within-condition baseline arms of any
#: split-half-referenced cross-phase statistic, and the estimator's own
#: test-retest. Written by ``imcoh_split_half_adjacencies``; file pattern
#: ``{patient}/{band}_{phase}_{A|B}_imcoh_{abs|sq}.npy``. The legacy
#: ``imcoh_halves_fc`` directory holds the ``abs`` arm only and is superseded.
IMCOH_HALVES_CACHE = CACHE_ROOT / "imcoh_halves"

# Experimental (failed debiasing approaches, kept for reference)
BIPOLAR_CACHE = CACHE_ROOT / "bipolar"
RESCALED_CACHE = CACHE_ROOT / "rescaled"

# ---------------------------------------------------------------------------
# Section 3: Reports (investigation outputs by research topic)
# ---------------------------------------------------------------------------
REPORTS_ROOT = DATA_ROOT / "reports"

# ---------------------------------------------------------------------------
# Section 4: General outputs
# ---------------------------------------------------------------------------
OUTPUTS_ROOT = DATA_ROOT / "outputs"
FIGURES_ROOT = OUTPUTS_ROOT / "figures"
TABLES_ROOT = OUTPUTS_ROOT / "tables"


# ---------------------------------------------------------------------------
# Section 5: FC-method → cache-directory routing
# ---------------------------------------------------------------------------

def fc_cache_for(fc_method: str) -> Path:
    """Return the FC matrix cache root for a given method.

    Note
    ----
    All three ImCoh variants (``"imcoh"``, ``"imcoh_abs"``, ``"imcoh_sq"``)
    share the same on-disk cache (the signed Nolte-2004 ImCoh).  The
    magnitude and squared-magnitude transforms are applied at load time
    by :func:`lrg_eegfc.workflow.fc.load_fc_matrix`.
    """
    _MAP = {
        "corr": CORR_CACHE,
        "msc": MSC_CACHE,
        "imcoh": IMCOH_CACHE,
        "imcoh_abs": IMCOH_CACHE,
        "imcoh_sq": IMCOH_CACHE,
    }
    if fc_method not in _MAP:
        raise ValueError(
            f"Unknown fc_method {fc_method!r}, expected one of {list(_MAP)}"
        )
    return _MAP[fc_method]


def lrg_cache_for(fc_method: str) -> Path:
    """Return the LRG result cache root for a given method.

    The signed ImCoh (``"imcoh"``) cannot feed LRG because the graph
    Laplacian requires non-negative edge weights.  Use the non-negative
    transforms ``"imcoh_abs"`` or ``"imcoh_sq"`` instead.
    """
    if fc_method == "imcoh":
        raise ValueError(
            "Signed ImCoh is in [-1, 1] and cannot feed LRG (Laplacian "
            "requires non-negative edge weights). Use fc_method='imcoh_abs' "
            "(|ImCoh|) or fc_method='imcoh_sq' (|ImCoh|^2)."
        )
    _MAP = {
        "corr": LRG_CACHE,
        "msc": LRG_CACHE,
        "imcoh_abs": IMCOH_LRG_CACHE,
        "imcoh_sq": IMCOH_LRG_CACHE,
    }
    if fc_method not in _MAP:
        raise ValueError(
            f"Unknown fc_method {fc_method!r}, expected one of {list(_MAP)}"
        )
    return _MAP[fc_method]


def lrg_filename(band: str, phase: str, fc_method: str) -> str:
    """Return the LRG ``.npz`` filename for a given band/phase/method.

    For ImCoh variants the transform is encoded in the suffix
    (``lrg_imcoh-abs.npz`` / ``lrg_imcoh-sq.npz``) so that both variants
    can coexist in :data:`IMCOH_LRG_CACHE`.  For ``corr`` and ``msc`` the
    filename retains the original plain form.
    """
    if fc_method in ("imcoh_abs", "imcoh_sq"):
        transform = fc_method.split("_", 1)[1]  # "abs" or "sq"
        return f"{band}_{phase}_lrg_imcoh-{transform}.npz"
    if fc_method == "imcoh":
        raise ValueError(
            "Signed ImCoh cannot be stored as an LRG result; use "
            "'imcoh_abs' or 'imcoh_sq'."
        )
    return f"{band}_{phase}_lrg_{fc_method}.npz"
