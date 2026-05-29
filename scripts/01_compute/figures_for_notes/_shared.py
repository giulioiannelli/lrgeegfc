"""Shared utilities for Section 2 (notes_imcoh) figures.

Centralizes FC-method-agnostic data loading, probe-block overlay drawing,
channel ordering, and visual constants so individual figure scripts stay
compact.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PHASE_LABELS,
    PATIENTS_4PHASE,
    nperseg_for_fs,
    SPRING_K,
    EDGE_WEIGHT_SCALING,
    EDGE_ALPHA_RANGE,
    EDGE_WIDTH_RANGE,
    EDGE_RANK_GAMMA,
    EDGE_RANK_WIDTH,
    EDGE_RANK_ALPHA,
    LAYOUT_METHOD,
    LAYOUT_SEED_DEFAULT,
    SPRING_ITER_DEFAULT,
    NODE_SIZE_DEFAULT,
    spring_k_for,
    FC_METHOD_DISPLAY_LABELS,
    SECTION2_FC_METHOD,
)

# Section 2 always contrasts MSC against the chosen ImCoh transform.
# Single source of truth for the (METHODS, METHOD_LABELS) pair used by
# every figure script — never redefine in individual scripts.
SECTION2_METHODS: list[str] = ["msc", SECTION2_FC_METHOD]
SECTION2_METHOD_LABELS: dict[str, str] = {
    m: FC_METHOD_DISPLAY_LABELS[m] for m in SECTION2_METHODS
}
from lrg_eegfc.config.paths import (
    FIGURES_ROOT,
    SEEG_DATAPATH,
)
from lrg_eegfc.utils.probe import (
    extract_probe_labels,
    build_probe_mask,
    compute_probe_weight_ratio,
    probe_weight_distributions,
    compute_enrichment_vs_scale,
)
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result

# Layout + drawing helpers were promoted to the library on 2026-05-29
# (library-naming meta-rule locked 2026-05-28: names reflect general
# graph/network concepts, never manuscript-local scope). Re-exported
# here so existing callers in figures_for_notes/ keep working without
# per-script import edits.
from lrg_eegfc.visuals.network_layouts import (
    compute_network_layout,
    compute_percolation_threshold,
)
from lrg_eegfc.visuals.network_drawing import (
    EDGE_GAMMA,
    scale_edge_weights,
    draw_network_edges,
    render_sbm_panel,
    render_lrg_panel,
)
from lrg_eegfc.visuals.network_drawing import _probe_color_map  # noqa: F401


# ---------------------------------------------------------------------------
# Patient / phase constants
# ---------------------------------------------------------------------------
ALL_PATIENTS = sorted(set(PATIENTS_4PHASE) | {"Pat_06"})
BANDS = list(BRAIN_BANDS_NAMES)
BAND_TEX = [BRAIN_BAND_TEX_DICT[b] for b in BANDS]

#: Sampling-rate overrides (only deviations from 2048 Hz).
from lrg_eegfc.config.const import FS_OVERRIDES as FS_MAP  # canonical

#: Default output root for Section 2 figures.
SECTION2_ROOT = FIGURES_ROOT / "section2"

# ---------------------------------------------------------------------------
# Visual constants
# ---------------------------------------------------------------------------
#: Per-method colours (consistent across all figures).
CLR_MSC = "#2166AC"
CLR_IMCOH = "#B2182B"
METHOD_COLORS = {"msc": CLR_MSC, "imcoh": CLR_IMCOH}

#: Same-probe / cross-probe colours.
CLR_SAME = "#E64A19"
CLR_CROSS = "#1565C0"

#: Per-band colours.
BAND_COLORS = {
    "delta": "#4C72B0", "theta": "#55A868", "alpha": "#C44E52",
    "beta": "#8172B3", "low_gamma": "#CCB974", "high_gamma": "#64B5CD",
}

#: Per-phase colours.
PHASE_COLORS = {
    "rest_pre": "#1f77b4", "task_learn": "#ff7f0e",
    "task_test": "#2ca02c", "rest_post": "#d62728",
}
PHASE_MARKERS = {"rest_pre": "o", "task_learn": "s", "task_test": "D", "rest_post": "^"}

#: Per-patient colours (tab10).
import matplotlib.pyplot as plt
_TAB10 = plt.get_cmap("tab10")
PATIENT_COLORS = {pat: _TAB10(i) for i, pat in enumerate(ALL_PATIENTS)}

#: Colormap for adjacency matrices.
CMAP_MATRIX = "magma"

# Spring layout config is now in config.const (SPRING_K, spring_k_for).

#: Community scales for enrichment analysis.
N_COMMUNITIES = [3, 5, 10, 15, 20, 30]


# ---------------------------------------------------------------------------
# Publication-quality matplotlib defaults
# ---------------------------------------------------------------------------
def apply_pub_style():
    """Set matplotlib rcParams for clean publication figures."""
    plt.rcParams.update({
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.grid": False,
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "font.family": "sans-serif",
    })


# ---------------------------------------------------------------------------
# Memory safety guards
# ---------------------------------------------------------------------------

#: Soft cap (bytes) on process virtual memory. Exceeding this raises MemoryError
#: instead of filling the system RAM. Default 12 GB; override per-script if needed.
DEFAULT_MEM_LIMIT_GB = 12.0


def set_memory_limit(limit_gb: float = DEFAULT_MEM_LIMIT_GB):
    """Cap the current process's address space.

    Call once at the top of any figure-producing script. Exceeding the cap
    will raise :class:`MemoryError` and halt the script cleanly, instead of
    swapping the whole machine.
    """
    import resource
    bytes_cap = int(limit_gb * 1024 ** 3)
    try:
        _, hard = resource.getrlimit(resource.RLIMIT_AS)
        new_hard = hard if hard > 0 else bytes_cap
        resource.setrlimit(resource.RLIMIT_AS, (bytes_cap, new_hard))
        print(f"  [mem] address space capped at {limit_gb:.1f} GB")
    except (ValueError, OSError) as exc:
        print(f"  [mem] could not set limit ({exc}); continuing unprotected")


def memory_usage_mb() -> float:
    """Return current process RSS in MB (linux-specific but portable enough)."""
    import resource
    # ru_maxrss is KB on Linux, bytes on macOS
    ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    import sys
    return ru / 1024 if sys.platform.startswith("linux") else ru / (1024 * 1024)


def checkpoint_cleanup(label: str = ""):
    """Call between iterations of a figure loop to keep RAM bounded.

    - Closes every matplotlib figure still held in pyplot's registry
    - Forces a garbage collection pass
    - Prints a one-line RSS update so runaway loops are obvious
    """
    import gc
    plt.close("all")
    gc.collect()
    if label:
        print(f"  [mem] after {label}: RSS = {memory_usage_mb():.0f} MB")


# ---------------------------------------------------------------------------
# Channel label loading
# ---------------------------------------------------------------------------

def load_channel_labels(patient: str) -> list[str]:
    """Load cleaned monopolar channel labels for *patient*.

    Handles varying CSV formats: with/without header, quoted labels,
    reference suffixes (e.g. ``,G2``).
    """
    for ext in ("csv", "txt"):
        fpath = SEEG_DATAPATH / patient / f"channel_labels.{ext}"
        if not fpath.exists():
            continue
        labels = []
        with open(fpath) as f:
            for i, line in enumerate(f):
                line = line.strip().strip('"')
                if not line:
                    continue
                if i == 0 and line.lower() == "label":
                    continue
                label = line.split(",")[0].strip().strip('"').replace(" ", "")
                if label:
                    labels.append(label)
        if labels:
            return labels
    raise FileNotFoundError(f"No channel_labels for {patient}")


# ---------------------------------------------------------------------------
# FC data loading (method-agnostic)
# ---------------------------------------------------------------------------

def load_fc_for_patient(
    patient: str,
    fc_method: str = "imcoh_abs",
    phases: Sequence[str] = PHASE_LABELS,
    bands: Sequence[str] = BANDS,
) -> dict[str, dict[str, NDArray | None]]:
    """Load FC matrices for all band x phase combinations.

    Returns nested dict ``{band: {phase: ndarray | None}}``.
    """
    data: dict[str, dict[str, NDArray | None]] = {}
    for band in bands:
        data[band] = {}
        for phase in phases:
            data[band][phase] = load_fc_matrix(patient, phase, band, fc_method)
    return data


def get_nperseg(patient: str) -> int:
    """Return correct nperseg for *patient*."""
    fs = FS_MAP.get(patient, 2048.0)
    return nperseg_for_fs(fs)


# ---------------------------------------------------------------------------
# Channel ordering by probe
# ---------------------------------------------------------------------------

def probe_sort_indices(channel_labels: Sequence[str]) -> NDArray:
    """Return indices that sort channels by probe, then by contact number."""
    probes = extract_probe_labels(channel_labels)
    unique = sorted(set(probes))
    return np.argsort([unique.index(p) * 1000 + i for i, p in enumerate(probes)])


def probe_boundaries(sorted_probes: Sequence[str]) -> list[int]:
    """Return boundary positions (first index of each new probe)."""
    return [i for i in range(1, len(sorted_probes)) if sorted_probes[i] != sorted_probes[i - 1]]


# ---------------------------------------------------------------------------
# Same-probe block outlines on heatmaps
# ---------------------------------------------------------------------------

from lrg_eegfc.visuals.fc_templates import (
    draw_probe_outlines as _lib_draw_probe_outlines,
)


def draw_probe_outlines(
    ax,
    channel_labels: Sequence[str],
    sort_idx: NDArray | None = None,
    color: str = "cyan",
    lw: float = 1.0,
    alpha: float = 0.8,
):
    """Backwards-compatible thin wrapper.

    Promoted to :func:`lrg_eegfc.visuals.draw_probe_outlines` on
    2026-05-10 (second-caller rule). New scripts should import the
    library version directly:

        from lrg_eegfc.visuals import draw_probe_outlines

    The wrapper preserves the original positional ``sort_idx`` for
    legacy callers in this folder.
    """
    _lib_draw_probe_outlines(
        ax, channel_labels,
        sort_idx=sort_idx,
        color=color, lw=lw, alpha=alpha,
    )


# ---------------------------------------------------------------------------
# Network metrics (extracted from gen_fc_figures_fast.py)
# ---------------------------------------------------------------------------

def compute_network_metrics(mat: NDArray) -> dict[str, float]:
    """Compute standard network metrics from an adjacency matrix."""
    import networkx as nx

    mat = mat.copy()
    np.fill_diagonal(mat, 0)
    triu_vals = mat[np.triu_indices_from(mat, k=1)]

    strengths = mat.sum(axis=1)
    G = nx.from_numpy_array(mat)
    cc = nx.average_clustering(G, weight="weight")

    si = strengths[np.triu_indices_from(mat, k=1)[0]]
    sj = strengths[np.triu_indices_from(mat, k=1)[1]]
    strength_assort = np.corrcoef(si, sj)[0, 1] if len(si) > 1 else 0.0

    sorted_w = np.sort(triu_vals)[::-1]
    n_top5 = max(1, int(0.05 * len(sorted_w)))
    weight_conc = sorted_w[:n_top5].sum() / (sorted_w.sum() + 1e-12)

    return {
        "mean_weight": triu_vals.mean(),
        "median_weight": np.median(triu_vals),
        "density_gt005": np.mean(triu_vals > 0.05),
        "density_gt01": np.mean(triu_vals > 0.1),
        "mean_strength": strengths.mean(),
        "std_strength": strengths.std(),
        "clustering_coeff": cc,
        "strength_assort": strength_assort,
        "weight_concentration_5pct": weight_conc,
    }


def compute_participation_coefficient(
    mat: NDArray,
    probe_labels: Sequence[str],
) -> NDArray:
    """Compute participation coefficient relative to probe membership.

    P_i = 1 - sum_s (k_{is} / k_i)^2

    where k_{is} = strength of node i to probe s, k_i = total strength.
    High P = edges spread across probes. Low P = edges concentrated on own probe.
    """
    mat = np.abs(mat.copy())
    np.fill_diagonal(mat, 0)
    N = mat.shape[0]
    probes = list(probe_labels[:N])
    unique_probes = sorted(set(probes))

    strengths = mat.sum(axis=1)
    P = np.ones(N)

    for i in range(N):
        if strengths[i] < 1e-12:
            P[i] = 0.0
            continue
        for s in unique_probes:
            mask = np.array([p == s for p in probes])
            k_is = mat[i, mask].sum()
            P[i] -= (k_is / strengths[i]) ** 2

    return P


# ---------------------------------------------------------------------------
# Laplacian eigenvalues
# ---------------------------------------------------------------------------

def compute_laplacian_eigenvalues(mat: NDArray) -> NDArray:
    """Compute sorted Laplacian eigenvalues from an adjacency matrix."""
    A = np.abs(mat.copy())
    np.fill_diagonal(A, 0)
    D = np.diag(A.sum(axis=1))
    L = D - A
    return np.sort(np.linalg.eigvalsh(L))


# ---------------------------------------------------------------------------
# Figure saving helper
# ---------------------------------------------------------------------------

def save_fig(fig, path: Path, dpi: int = 300, fmt: str | None = None):
    """Save figure, create dirs, close.

    If *fmt* is None, infers from path suffix (defaults to pdf).
    For exploration / helper figures, pass fmt="png" and dpi=100 to keep
    file size small (< 5 MB) so figures can always be embedded in replies.
    """
    path = Path(path)
    if fmt is not None:
        path = path.with_suffix(f".{fmt}")
    elif not path.suffix:
        path = path.with_suffix(".pdf")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"  Saved: {path}")
    return path


def save_helper_fig(fig, path: Path):
    """Save an exploration / sweep figure as PNG at low DPI.
    Guaranteed to stay small enough (< ~10 MB) for embedding in context.
    """
    return save_fig(fig, path, dpi=90, fmt="png")


# ---------------------------------------------------------------------------
# FC-method-aware edge scaling + drawing
# ---------------------------------------------------------------------------
# Promoted to lrg_eegfc.visuals.network_drawing on 2026-05-29:
#   EDGE_GAMMA, scale_edge_weights, draw_network_edges, _probe_color_map.
# Re-exported at the top of this module for backwards compatibility.


# ---------------------------------------------------------------------------
# Representative example sets (variety)
# ---------------------------------------------------------------------------

#: Default representative patients (those with valid coordinates + diversity).
REPR_PATIENTS = ["Pat_02", "Pat_05", "Pat_08"]

#: Default representative bands (low, mid, high frequency).
REPR_BANDS = ["alpha", "beta"]

#: Default representative phases (always include at least one task phase).
REPR_PHASES = ["rest_pre", "task_learn", "rest_post"]

#: Patients with valid MNI coordinates (for brain connectome figures).
PATIENTS_WITH_COORDS = ["Pat_02", "Pat_03", "Pat_05"]


# ---------------------------------------------------------------------------
# Canonical network layout + SBM/LRG panel rendering
# ---------------------------------------------------------------------------
# Promoted to lrg_eegfc.visuals.network_layouts on 2026-05-29:
#   compute_network_layout, compute_percolation_threshold
#   (plus internals _seeded_initial_pos, _compute_sbm_data).
# Promoted to lrg_eegfc.visuals.network_drawing on 2026-05-29:
#   render_sbm_panel, render_lrg_panel.
# Re-exported at the top of this module for backwards compatibility.
