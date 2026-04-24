"""WP0 — Metric Exploration: shared constants, data loaders, and utilities."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    DEFAULT_NPERSEG,
    PHASE_LABELS,
)
from lrg_eegfc.config.paths import LRG_CACHE, MSC_CACHE, REPORTS_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.utils.metrics.reorganization import (
    build_metric_specs,
    compute_cluster_labels,
    cluster_swap_coefficient,
)
from sklearn.metrics import adjusted_rand_score

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
BANDS = list(BRAIN_BANDS_NAMES)
PHASES = list(PHASE_LABELS)  # rest_pre, task_learn, task_test, rest_post
FC_METHOD = "msc"
ALL_PAIRS = list(combinations(PHASES, 2))  # 6 pairs
PAIR_SHORT = {
    (a, b): f"{a[:2]}{a[-4:]}-{b[:2]}{b[-4:]}" for a, b in ALL_PAIRS
}

REST = {"rest_pre", "rest_post"}
TASK = {"task_learn", "task_test"}

OUT_ROOT = REPORTS_ROOT / "metric_exploration"


def classify_pair(p1: str, p2: str) -> str:
    """Classify a phase pair as 'within' or 'cross'."""
    s = {p1, p2}
    if s <= REST or s <= TASK:
        return "within"
    return "cross"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_all_lrg(
    patients: list[str] | None = None,
    verbose: bool = True,
) -> dict:
    """Load all cached LRG results.

    Returns dict keyed by (patient, phase, band) -> LRGResult.
    """
    patients = patients or PATIENTS
    data = {}
    n_loaded = n_missing = 0
    for pat in patients:
        for band in BANDS:
            for phase in PHASES:
                res = load_lrg_result(
                    pat, phase, band, FC_METHOD, cache_root=LRG_CACHE
                )
                if res is not None:
                    data[(pat, phase, band)] = res
                    n_loaded += 1
                else:
                    n_missing += 1
                    if verbose:
                        print(f"  MISS: {pat} {phase} {band}")
    if verbose:
        print(f"  LRG loaded: {n_loaded}, missing: {n_missing}")
    return data


def load_all_fc(
    patients: list[str] | None = None,
    verbose: bool = True,
) -> dict:
    """Load all cached raw MSC FC matrices.

    Returns dict keyed by (patient, phase, band) -> np.ndarray.
    """
    patients = patients or PATIENTS
    data = {}
    n_loaded = n_missing = 0
    for pat in patients:
        for band in BANDS:
            for phase in PHASES:
                W = load_msc_matrix(
                    pat, phase, band,
                    cache_root=MSC_CACHE,
                    sparsify="none",
                    n_surrogates=0,
                    nperseg=DEFAULT_NPERSEG,
                )
                if W is not None:
                    data[(pat, phase, band)] = W
                    n_loaded += 1
                else:
                    n_missing += 1
                    if verbose:
                        print(f"  MISS FC: {pat} {phase} {band}")
    if verbose:
        print(f"  FC loaded: {n_loaded}, missing: {n_missing}")
    return data


# ---------------------------------------------------------------------------
# Raw FC baseline metrics  (NEW — Task 2B)
# ---------------------------------------------------------------------------
def _upper_tri(W: np.ndarray) -> np.ndarray:
    """Extract upper-triangle elements (k=1)."""
    idx = np.triu_indices_from(W, k=1)
    return W[idx]


def fc_frobenius(W1: np.ndarray, W2: np.ndarray) -> float:
    """Frobenius distance of upper triangles: ||W1 - W2||_F."""
    return float(np.linalg.norm(_upper_tri(W1) - _upper_tri(W2)))


def fc_scaled_frobenius(W1: np.ndarray, W2: np.ndarray) -> float:
    """Scaled Frobenius: ||W1 - W2||_F / (||W1||_F + ||W2||_F)."""
    v1, v2 = _upper_tri(W1), _upper_tri(W2)
    diff = np.linalg.norm(v1 - v2)
    denom = np.linalg.norm(v1) + np.linalg.norm(v2)
    return float(diff / max(denom, 1e-15))


def fc_rank_distance(W1: np.ndarray, W2: np.ndarray) -> float:
    """Rank distance: 1 - Spearman(upper_tri(W1), upper_tri(W2))."""
    v1, v2 = _upper_tri(W1), _upper_tri(W2)
    corr, _ = spearmanr(v1, v2)
    return float(1.0 - corr)


def fc_mean_abs_diff(W1: np.ndarray, W2: np.ndarray) -> float:
    """Mean absolute difference of upper triangles."""
    return float(np.mean(np.abs(_upper_tri(W1) - _upper_tri(W2))))


# ---------------------------------------------------------------------------
# Metric registries
# ---------------------------------------------------------------------------
# Each entry: (label, is_distance, range_str, reference_value_note)
_METRIC_META = {
    # --- LRG ultrametric value-based ---
    "matrix_distance":     ("Ultrametric matrix distance",           True,  "[0, ∞)",   "0 for identical"),
    "scaled_distance":     ("Scaled ultrametric distance (log)",     True,  "[0, 1]",   "0 for identical"),
    "rank_distance":       ("1 − Spearman rank corr.",               True,  "[0, 2]",   "0 for identical"),
    "quantile_rmse":       ("Quantile RMSE (log)",                   True,  "[0, ∞)",   "0 for identical"),
    "permutation_robust":  ("Permutation-robust distance",           True,  "[0, ∞)",   "0 for identical"),
    # --- LRG tree-structural ---
    "tree_robinson_foulds": ("Robinson-Foulds distance",             True,  "[0, 1]",   "0 for identical topology"),
    "tree_cophenetic_corr": ("Cophenetic correlation",               False, "[-1, 1]",  "1 for identical"),
    "tree_baker_gamma":     ("Baker gamma",                          False, "[-1, 1]",  "1 for identical"),
    "tree_fowlkes_mallows": ("Fowlkes-Mallows index",               False, "[0, 1]",   "1 for identical"),
    # --- LRG partition-based ---
    "ari":                  ("Adjusted Rand Index",                  False, "[-1, 1]",  "1 for identical, ~0 for random"),
    "cluster_swap":         ("Cluster swap coefficient",             False, "[0, 1]",   "1 for identical"),
    # --- Raw FC baselines ---
    "fc_frobenius":         ("FC Frobenius distance",                True,  "[0, ∞)",   "0 for identical; unbounded"),
    "fc_scaled_frobenius":  ("FC scaled Frobenius",                  True,  "[0, 1]",   "0 for identical"),
    "fc_rank_distance":     ("FC rank distance (1 − Spearman)",      True,  "[0, 2]",   "0 for identical"),
    "fc_mean_abs_diff":     ("FC mean absolute difference",          True,  "[0, ∞)",   "0 for identical; unbounded"),
}


def _ensure_square(m: np.ndarray) -> np.ndarray:
    if m.ndim == 1:
        return squareform(m)
    return m


def get_lrg_metric_specs() -> dict:
    """Return dict of LRG metric specs: name -> {'label', 'fn', 'is_distance', 'range', 'ref'}.

    The 'fn' has signature (lrg_a, lrg_b) -> float where lrg_a/b are LRGResult objects.
    """
    raw_specs = build_metric_specs()

    def _wrap(metric_fn):
        """Wrap build_metric_specs fn(U1, U2, Z1, Z2) to accept (lrg_a, lrg_b)."""
        def wrapped(a, b):
            U1 = _ensure_square(a.ultrametric_matrix)
            U2 = _ensure_square(b.ultrametric_matrix)
            return metric_fn(U1, U2, a.linkage_matrix, b.linkage_matrix)
        return wrapped

    specs = {}
    for name, spec in raw_specs.items():
        meta = _METRIC_META[name]
        specs[name] = {
            "label": meta[0],
            "fn": _wrap(spec["fn"]),
            "is_distance": meta[1],
            "range": meta[2],
            "ref": meta[3],
            "level": "lrg",
        }

    # Partition-based metrics (ARI + cluster_swap at optimal threshold)
    def _ari(a, b):
        if a.linkage_matrix is None or b.linkage_matrix is None:
            return np.nan
        if a.n_nodes != b.n_nodes:
            return np.nan
        la = compute_cluster_labels(a.linkage_matrix, a.optimal_threshold)
        lb = compute_cluster_labels(b.linkage_matrix, b.optimal_threshold)
        return float(adjusted_rand_score(la, lb))

    def _swap(a, b):
        if a.linkage_matrix is None or b.linkage_matrix is None:
            return np.nan
        if a.n_nodes != b.n_nodes:
            return np.nan
        la = compute_cluster_labels(a.linkage_matrix, a.optimal_threshold)
        lb = compute_cluster_labels(b.linkage_matrix, b.optimal_threshold)
        return cluster_swap_coefficient(la, lb)

    meta_ari = _METRIC_META["ari"]
    specs["ari"] = {
        "label": meta_ari[0], "fn": _ari,
        "is_distance": meta_ari[1], "range": meta_ari[2], "ref": meta_ari[3],
        "level": "lrg",
    }
    meta_swap = _METRIC_META["cluster_swap"]
    specs["cluster_swap"] = {
        "label": meta_swap[0], "fn": _swap,
        "is_distance": meta_swap[1], "range": meta_swap[2], "ref": meta_swap[3],
        "level": "lrg",
    }
    return specs


def get_fc_metric_specs() -> dict:
    """Return dict of raw FC metric specs: name -> {'label', 'fn', 'is_distance', 'range', 'ref'}.

    The 'fn' has signature (W1, W2) -> float where W1/W2 are N×N arrays.
    """
    return {
        "fc_frobenius": {
            "label": _METRIC_META["fc_frobenius"][0],
            "fn": fc_frobenius,
            "is_distance": True, "range": "[0, ∞)", "ref": "0 for identical; unbounded",
            "level": "raw_fc",
        },
        "fc_scaled_frobenius": {
            "label": _METRIC_META["fc_scaled_frobenius"][0],
            "fn": fc_scaled_frobenius,
            "is_distance": True, "range": "[0, 1]", "ref": "0 for identical",
            "level": "raw_fc",
        },
        "fc_rank_distance": {
            "label": _METRIC_META["fc_rank_distance"][0],
            "fn": fc_rank_distance,
            "is_distance": True, "range": "[0, 2]", "ref": "0 for identical",
            "level": "raw_fc",
        },
        "fc_mean_abs_diff": {
            "label": _METRIC_META["fc_mean_abs_diff"][0],
            "fn": fc_mean_abs_diff,
            "is_distance": True, "range": "[0, ∞)", "ref": "0 for identical; unbounded",
            "level": "raw_fc",
        },
    }


# ---------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------
def save_fig(
    fig,
    path: Path,
    *,
    what: str,
    proves: str,
    how_to_read: str,
    dpi: int = 300,
):
    """Save figure as PDF and write companion .md file, then close fig."""
    import matplotlib.pyplot as plt

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"  Saved: {path}")

    md_path = path.with_suffix(".md")
    md_path.write_text(
        f"# {path.stem}\n\n"
        f"## What the figure shows\n{what}\n\n"
        f"## What it proves\n{proves}\n\n"
        f"## How to read it\n{how_to_read}\n"
    )


def to_distance(value: float, is_distance: bool) -> float:
    """Convert metric value to distance convention (higher = more different).

    For similarity metrics (is_distance=False), negate or invert as appropriate.
    """
    if is_distance:
        return value
    # Similarity metrics: negate so higher = more different
    return -value
