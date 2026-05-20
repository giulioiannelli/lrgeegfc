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
# FC-method-aware edge scaling
# ---------------------------------------------------------------------------

#: Gamma exponent per FC method for edge width/alpha scaling.
#: Higher gamma → only top-weight edges are visually prominent.
EDGE_GAMMA = {
    "msc": 2.0,
    "imcoh": 6.0,
    "imcoh_abs": 6.0,
    "imcoh_sq": 6.0,
    "corr": 2.0,
}


def scale_edge_weights(
    weights: NDArray,
    fc_method: str = "imcoh_abs",
) -> tuple[NDArray, NDArray]:
    """Scale edge weights to visual width and alpha using method-aware config.

    Uses a gamma transform t = (w / w_max) ** gamma so that weak edges
    get small widths and high-weight edges stand out. Works well for
    near-uniform distributions (ImCoh) at gamma ≥ 6 and for heavy-tailed
    ones (MSC) at gamma ≈ 2.
    """
    w_min, w_max = EDGE_WIDTH_RANGE.get(fc_method, (0.15, 4.0))
    a_min, a_max = EDGE_ALPHA_RANGE.get(fc_method, (0.03, 0.9))
    gamma = EDGE_GAMMA.get(fc_method, 2.0)

    w_norm = weights / (weights.max() + 1e-30)
    t = w_norm ** gamma

    widths = w_min + t * (w_max - w_min)
    alphas = a_min + t * (a_max - a_min)
    return widths, alphas


def _probe_color_map(probes: Sequence[str]) -> dict:
    """Return {probe_id: rgba} using the same tab20 palette as the nodes."""
    uniq = sorted(set(probes))
    cmap = plt.get_cmap("tab20", len(uniq))
    return {p: cmap(i) for i, p in enumerate(uniq)}


def draw_network_edges(
    ax,
    pos_arr: NDArray,
    mat: NDArray,
    fc_method: str = "imcoh_abs",
    probe_labels: Sequence[str] | None = None,
    highlight_same_probe: bool = True,
    base_color=(0.0, 0.0, 0.0),
    gamma: float | None = None,
    wmin: float | None = None,
    wmax: float | None = None,
    amin: float | None = None,
    amax: float | None = None,
    scaling: str = "rank",
    value_range: tuple[float, float] | None = None,
):
    """Edge drawing per method (gamma, width and alpha from
    :mod:`config.const`, overridable per call).

    *scaling*:
      - ``"rank"`` (default): ``t = ((rank+1)/N)**gamma`` — scale-invariant
        across patients but visual density depends on edge count.
      - ``"value"``: per-patient min-max normalise weights to [0,1], then
        ``t = w_norm**gamma`` — same width/alpha mapping applies to the
        same *relative* weight level across patients, so the γ/wmin/wmax/
        αmin/αmax knobs produce comparable density across subjects.

    If *highlight_same_probe* is ``True`` and *probe_labels* is given,
    same-probe edges are drawn in the shaft's tab20 colour at full alpha
    (so they stand out); cross-probe edges keep *base_color* with the
    value/rank-based alpha fade.
    """
    from matplotlib.collections import LineCollection

    _g = EDGE_RANK_GAMMA.get(fc_method, 2.0)
    _wmn, _wmx = EDGE_RANK_WIDTH.get(fc_method, (0.15, 4.0))
    _amn, _amx = EDGE_RANK_ALPHA.get(fc_method, (0.06, 0.9))
    gamma = _g if gamma is None else gamma
    wmin = _wmn if wmin is None else wmin
    wmax = _wmx if wmax is None else wmax
    amin = _amn if amin is None else amin
    amax = _amx if amax is None else amax

    N = mat.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = mat[r, c]
    active = np.where(w > 0)[0]
    if active.size == 0:
        return
    w_act = w[active]
    ii, jj = r[active], c[active]

    # Same-probe mask BEFORE scaling so we can exclude same-probe from the
    # normalisation range (otherwise same-probe edges dominate the top of
    # the distribution and cross-probe edges all squash to t ≈ 0).
    if highlight_same_probe and probe_labels is not None:
        same_all = np.array(
            [probe_labels[ii[k]] == probe_labels[jj[k]]
             for k in range(len(active))],
            dtype=bool,
        )
        probe_cmap = _probe_color_map(probe_labels)
    else:
        same_all = np.zeros(len(active), dtype=bool)
        probe_cmap = {}

    cross_mask = ~same_all
    t = np.zeros_like(w_act, dtype=float)

    def _value_t(w_sub: np.ndarray) -> np.ndarray:
        if w_sub.size == 0:
            return np.zeros(0)
        lo, hi = (float(np.percentile(w_sub, 1)),
                  float(np.percentile(w_sub, 99)))
        if hi <= lo:
            return np.zeros_like(w_sub)
        wn = np.clip((w_sub - lo) / (hi - lo), 0.0, 1.0)
        return wn ** gamma

    def _rank_t(w_sub: np.ndarray) -> np.ndarray:
        if w_sub.size == 0:
            return np.zeros(0)
        ranks = np.argsort(np.argsort(w_sub)).astype(float)
        return ((ranks + 1) / len(ranks)) ** gamma

    if scaling == "value":
        if value_range is not None:
            # Caller-supplied global range applies to all edges uniformly
            lo, hi = float(value_range[0]), float(value_range[1])
            if hi > lo:
                wn = np.clip((w_act - lo) / (hi - lo), 0.0, 1.0)
                t = wn ** gamma
        else:
            # Per-population scaling: cross-probe and same-probe each get
            # their own 1st/99th-percentile window, so within each
            # population we still see relative weight differences (instead
            # of every same-probe edge clipping to t=1 = wmax).
            t[cross_mask] = _value_t(w_act[cross_mask])
            t[~cross_mask] = _value_t(w_act[~cross_mask])
    else:
        t[cross_mask] = _rank_t(w_act[cross_mask])
        t[~cross_mask] = _rank_t(w_act[~cross_mask])

    widths = wmin + t * (wmax - wmin)
    alphas = amin + t * (amax - amin)

    base_rgb = np.array(base_color, dtype=float)

    # Two passes: cross-probe first (lower zorder, painted first), then
    # same-probe on top (higher zorder, full alpha, shaft-coloured) so
    # the coloured shaft edges always sit above the black cross-probe
    # background regardless of magnitude.
    #
    # PDF size: each edge with its own (color, width, alpha) emits a
    # separate graphics state in the PDF. With ~7000 edges per panel
    # this makes vector PDFs ~2-3 MB. We bucket the styles into a small
    # set of (linewidth, alpha) bins per (rgb-color) and drop edges that
    # would render at alpha < 0.01 (invisible). Each bucket becomes one
    # LineCollection → 1-2 orders of magnitude smaller PDF.
    N_LW_BINS = 12
    N_AL_BINS = 12
    ALPHA_INVISIBLE = 0.01

    def _emit_bucketed(segs_arr, widths_arr, alphas_arr, rgb_arr, zorder):
        """Quantise (lw, alpha) per (rgb tuple) and emit one LineCollection
        per bucket. Drops edges that would render at alpha < ALPHA_INVISIBLE.
        """
        if len(segs_arr) == 0:
            return
        # Drop near-invisible edges
        keep = alphas_arr >= ALPHA_INVISIBLE
        if not keep.any():
            return
        segs_arr = segs_arr[keep]
        widths_arr = widths_arr[keep]
        alphas_arr = alphas_arr[keep]
        rgb_arr = rgb_arr[keep]
        # Quantise width and alpha
        w_lo, w_hi = float(widths_arr.min()), float(widths_arr.max())
        a_lo, a_hi = float(alphas_arr.min()), float(alphas_arr.max())
        if w_hi > w_lo:
            w_idx = np.clip(((widths_arr - w_lo) / (w_hi - w_lo)
                              * (N_LW_BINS - 1)).round().astype(int),
                            0, N_LW_BINS - 1)
            w_bins = w_lo + (w_hi - w_lo) * (np.arange(N_LW_BINS) / (N_LW_BINS - 1))
        else:
            w_idx = np.zeros(len(widths_arr), dtype=int)
            w_bins = np.array([w_lo])
        if a_hi > a_lo:
            a_idx = np.clip(((alphas_arr - a_lo) / (a_hi - a_lo)
                              * (N_AL_BINS - 1)).round().astype(int),
                            0, N_AL_BINS - 1)
            a_bins = a_lo + (a_hi - a_lo) * (np.arange(N_AL_BINS) / (N_AL_BINS - 1))
        else:
            a_idx = np.zeros(len(alphas_arr), dtype=int)
            a_bins = np.array([a_lo])
        # Quantise RGB to 8-bit so per-shaft variants collapse into one bucket
        rgb_q = np.round(rgb_arr * 255).astype(np.uint8)
        # Build bucket key
        keys = (rgb_q[:, 0].astype(np.int64) << 32
                | rgb_q[:, 1].astype(np.int64) << 24
                | rgb_q[:, 2].astype(np.int64) << 16
                | (w_idx.astype(np.int64) & 0xFF) << 8
                | (a_idx.astype(np.int64) & 0xFF))
        unique_keys, inverse = np.unique(keys, return_inverse=True)
        for b, key in enumerate(unique_keys):
            mask = inverse == b
            r = (key >> 32) & 0xFF
            g = (key >> 24) & 0xFF
            bl = (key >> 16) & 0xFF
            wi = (key >> 8) & 0xFF
            ai = key & 0xFF
            color = (float(r) / 255.0, float(g) / 255.0, float(bl) / 255.0,
                     float(a_bins[ai]))
            ax.add_collection(LineCollection(
                segs_arr[mask],
                colors=[color],
                linewidths=float(w_bins[wi]),
                zorder=zorder, rasterized=False,
            ))

    for layer_mask, layer_zorder, force_alpha in (
        (cross_mask, 1, None),     # cross-probe — alphas from t
        (~cross_mask, 2, 1.0),     # same-probe — full alpha, on top
    ):
        if not layer_mask.any():
            continue
        sub_idx = np.where(layer_mask)[0]
        sub_t = t[sub_idx]
        order = sub_idx[np.argsort(sub_t)]  # strong edges painted last
        segs = np.stack([pos_arr[ii[order]], pos_arr[jj[order]]], axis=1)
        widths_o = widths[order]
        alphas_o = (np.full(len(order), force_alpha)
                    if force_alpha is not None else alphas[order])
        rgb = np.tile(base_rgb, (len(order), 1))
        if layer_zorder == 2:  # same-probe shaft colours
            for k_pos, k in enumerate(order):
                col = probe_cmap.get(probe_labels[ii[k]])
                if col is not None:
                    rgb[k_pos] = col[:3]
        _emit_bucketed(segs, widths_o, alphas_o, rgb, layer_zorder)
        # Skip the original single-LineCollection emission below.
        continue
        colors = np.concatenate([rgb, alphas_o[:, None]], axis=1)
        ax.add_collection(LineCollection(
            segs, colors=colors, linewidths=widths_o,
            zorder=layer_zorder, rasterized=False,
        ))


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
# Canonical network layout (method-aware)
# ---------------------------------------------------------------------------

def _seeded_initial_pos(N: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    theta = 2 * np.pi * np.arange(N) / N + rng.uniform(-0.2, 0.2, N)
    pos0 = np.column_stack([np.cos(theta), np.sin(theta)])
    return {i: pos0[i] for i in range(N)}


def compute_network_layout(
    A: NDArray,
    fc_method: str,
    patient: str | None = None,
    phase: str | None = None,
    band: str | None = None,
    seed: int = LAYOUT_SEED_DEFAULT,
    spring_iter: int = SPRING_ITER_DEFAULT,
    spring_k: float | None = 0.07,
    layout_method: str | None = None,
) -> NDArray:
    """Return an (N, 2) layout matrix for an FC adjacency.

    Method chosen by :data:`LAYOUT_METHOD` per fc_method:
      - ``"mds_ultrametric"`` (default): classical MDS / PCoA on the LRG
        ultrametric distance matrix D_ij. Eigendecomposition of the
        double-centered ``-D²/2`` matrix; the first two eigenvectors
        scaled by ``√λ`` are the 2D coordinates. Deterministic, no
        iterative optimisation, distances preserved in least-squares.
      - ``"spring_lrg_distance"``: spring with edge weight ``~1/D_ij``.
      - ``"kk_ultrametric"``: Kamada-Kawai using D_ij directly.
      - ``"spring"``: spring on raw FC weights.
      - Falls back to spring on raw FC if the LRG cache is missing.
    """
    import networkx as nx
    from scipy.spatial.distance import squareform
    from lrg_eegfc.workflow.lrg import load_lrg_result

    N = A.shape[0]
    # Per-call override (e.g. fig_F forcing all methods to mds_lrg_continuous)
    # takes precedence over the per-fc_method config dispatch.
    method = (layout_method if layout_method is not None
              else LAYOUT_METHOD.get(fc_method, "mds_ultrametric"))

    needs_lrg = method in (
        "kk_ultrametric", "spring_lrg_distance",
        "mds_ultrametric", "mds_log_ultrametric",
        "mds_lrg_continuous",
    )
    if needs_lrg and all(x is not None for x in (patient, phase, band)):
        lrg = load_lrg_result(patient, phase, band, fc_method)
        if lrg is not None:
            um = lrg.ultrametric_matrix
            U = squareform(um) if um.ndim == 1 else um
            M = min(N, U.shape[0])

            if method == "mds_lrg_continuous":
                # Continuous LRG transition distance ~1/rho_ij from the
                # heat kernel at tau = susceptibility-peak τ* (where
                # C(τ) = -dS/d log τ is maximum). This is the diffusion
                # time at which the heat kernel discriminates *most*
                # between communities, so the MDS embedding sees the
                # richest multi-scale structure. The previous choice
                # `optimal_threshold` is calibrated for community
                # detection (where to cut the dendrogram), not for
                # spatial embedding — for MSC it sits past equilibrium
                # and produces a near-uniform Trho ("star plot").
                from lrgsglib.utils.lrg.infocomm import lapl_dists
                from scipy.spatial.distance import squareform as _sf
                # Reload the FC matrix that fed this LRG (|FC| ≥ 0)
                Adense = np.abs(A.copy())
                np.fill_diagonal(Adense, 0.0)
                # Restrict to the giant component used in LRG (M nodes)
                if Adense.shape[0] != M:
                    Adense = Adense[:M, :M]
                L = np.diag(Adense.sum(axis=1)) - Adense

                # τ = 1 / λ_max — inverse fastest-mode timescale of the
                # graph Laplacian. A real diffusion time (units of the
                # Laplacian's inverse). Maximises Trho coefficient-of-
                # variation across both MSC and |ImCoh| (verified
                # empirically). Dimensionless across methods (each gets
                # its own λ_max), so layouts are comparable.
                #
                # DO NOT use lrg.optimal_threshold here — that field is
                # a dendrogram cut HEIGHT in the cophenetic ultrametric,
                # NOT a diffusion time. They live in different spaces.
                eigvals_L = np.linalg.eigvalsh(L)
                lam_max = float(np.max(eigvals_L))
                if lam_max <= 0:
                    raise ValueError(
                        "Cannot pick τ for layout: Laplacian λ_max ≤ 0; "
                        "graph likely empty or disconnected."
                    )
                tau_star = 1.0 / lam_max

                Dcond = lapl_dists(L, tau=tau_star)
                # Winsorize at 99th percentile to prevent a single weakly-
                # connected node (very large Trho to everyone else) from
                # dominating the MDS eigendecomposition. The capped node
                # still ends up at the periphery, but the bulk geometry
                # is preserved instead of collapsing to a core.
                cap = float(np.percentile(Dcond, 99))
                Dcond = np.minimum(Dcond, cap)
                Dwork = _sf(Dcond)  # condensed -> square
                D2 = Dwork.astype(np.float64) ** 2
                n = M
                J = np.eye(n) - np.full((n, n), 1.0 / n)
                B = -0.5 * J @ D2 @ J
                B = 0.5 * (B + B.T)
                eigvals, eigvecs = np.linalg.eigh(B)
                idx = np.argsort(eigvals)[::-1][:2]
                lam = np.maximum(eigvals[idx], 0.0)
                coords = eigvecs[:, idx] * np.sqrt(lam)
                pos = {i: coords[i] for i in range(M)}

            elif method in ("mds_ultrametric", "mds_log_ultrametric"):
                # Classical MDS / PCoA. B = -1/2 * J D^2 J with
                # J = I - 1/n * 1·1^T (double-centering matrix).
                # log variant compresses the multi-scale span of D
                # (LRG distances span orders of magnitude — a few
                # extreme entries otherwise dominate the eigenvalues
                # and squash the bulk into a triangle).
                if method == "mds_log_ultrametric":
                    Dwork = np.log1p(U.astype(np.float64))
                else:
                    Dwork = U.astype(np.float64)
                D2 = Dwork ** 2
                n = M
                J = np.eye(n) - np.full((n, n), 1.0 / n)
                B = -0.5 * J @ D2 @ J
                # Symmetrise + eigendecompose
                B = 0.5 * (B + B.T)
                eigvals, eigvecs = np.linalg.eigh(B)
                idx = np.argsort(eigvals)[::-1][:2]  # top-2 eigenvalues
                lam = np.maximum(eigvals[idx], 0.0)  # clip negatives
                coords = eigvecs[:, idx] * np.sqrt(lam)
                pos = {i: coords[i] for i in range(M)}

            elif method == "kk_ultrametric":
                G = nx.complete_graph(M)
                dist = {i: {j: float(U[i, j]) for j in range(M) if j != i}
                        for i in range(M)}
                pos = nx.kamada_kawai_layout(
                    G, dist=dist, pos=_seeded_initial_pos(M, seed),
                )
            else:  # spring_lrg_distance
                with np.errstate(divide="ignore"):
                    W = np.where(U > 0, 1.0 / U, 0.0)
                np.fill_diagonal(W, 0.0)
                G = nx.from_numpy_array(W)
                k_val = spring_k if spring_k is not None else spring_k_for(
                    fc_method, band or "beta", M,
                )
                pos = nx.spring_layout(
                    G, k=k_val, iterations=spring_iter, seed=seed,
                    weight="weight",
                    pos=_seeded_initial_pos(M, seed),
                )
            arr = np.zeros((N, 2))
            for i in range(M):
                arr[i] = pos[i]
            return arr

    # Nested SBM layout (graph-tool): community-aware radial embedding.
    # Returns positions only; use `draw_sbm_hierarchical` for the full
    # SBM + bezier-edge pipeline.
    if method == "sbm_nested":
        data = _compute_sbm_data(A, seed=seed)
        return data["pos"]

    # Final fallback: spring on raw FC weights
    G = nx.from_numpy_array(A)
    k_val = spring_k if spring_k is not None else spring_k_for(
        fc_method, band or "beta", N,
    )
    pos = nx.spring_layout(
        G, k=k_val, iterations=spring_iter, seed=seed, weight="weight",
    )
    return np.array([pos[i] for i in range(N)])


# ---------------------------------------------------------------------------
# Nested SBM hierarchical layout (graph-tool): SBM + bezier edges
# ---------------------------------------------------------------------------

def _compute_sbm_data(A: NDArray, seed: int = LAYOUT_SEED_DEFAULT,
                     beta: float = 0.8) -> dict:
    """Nested-SBM inference + hierarchy tree + bezier control points.

    Returns a dict with keys:
        ``pos``   : (N, 2) vertex positions (leaves of the hierarchy tree).
        ``g``     : graph-tool Graph.
        ``weight``: edge property map with |FC| weights.
        ``state`` : fitted NestedBlockState.
        ``t``     : hierarchy tree (graph).
        ``tpos``  : tree-vertex positions (radial).
        ``cts``   : edge control points in the LOCAL edge frame
                    (x in [0, 1] along source→target, y is bend offset).
    """
    import graph_tool.all as gt
    gt.seed_rng(seed)
    N = A.shape[0]

    g = gt.Graph(directed=False)
    g.add_vertex(N)
    weight = g.new_edge_property("double")
    ri, ci = np.triu_indices(N, k=1)
    vals = A[ri, ci]
    mask = vals > 0
    elist = np.column_stack([ri[mask], ci[mask], vals[mask]])
    g.add_edge_list(elist, eprops=[weight])

    state = gt.minimize_nested_blockmodel_dl(
        g, state_args=dict(recs=[weight], rec_types=["real-exponential"])
    )

    t, _tb, _tpos_base = gt.get_hierarchy_tree(state)
    root = next(v for v in t.vertices() if v.in_degree() == 0)
    tpos = gt.radial_tree_layout(t, root)
    cts = gt.get_hierarchy_control_points(g, t, tpos, beta=beta)

    pos = np.array(
        [[tpos[t.vertex(int(v))][0], tpos[t.vertex(int(v))][1]]
         for v in g.vertices()]
    )
    return {"pos": pos, "g": g, "weight": weight, "state": state,
            "t": t, "tpos": tpos, "cts": cts}


def compute_percolation_threshold(A: NDArray) -> float:
    """Minimum edge weight needed to keep the (undirected) graph connected.

    Defined as the minimum weight in the *maximum* spanning tree of |A|.
    Keeping edges with weight ``>= threshold`` yields the minimally
    connected sub-graph (the classical bond-percolation threshold from
    above).
    """
    from scipy.sparse.csgraph import minimum_spanning_tree
    from scipy.sparse import csr_matrix

    A_pos = np.where(A > 0, A, 0.0)
    np.fill_diagonal(A_pos, 0.0)
    # Max spanning tree == MST on negated weights.
    mst_neg = minimum_spanning_tree(csr_matrix(-A_pos)).toarray()
    mst = np.where(mst_neg < 0, -mst_neg, 0.0)
    mst_w = mst[mst > 0]
    return float(mst_w.min()) if len(mst_w) else 0.0


def render_sbm_panel(
    output_path: Path,
    A: NDArray,
    probe_labels: list[str] | None = None,
    shaft_colors: list | None = None,
    *,
    seed: int = LAYOUT_SEED_DEFAULT,
    threshold: float | None = None,
    threshold_scale: float = 1.0,
    highlight_same_probe: bool = True,
    output_size: tuple = (600, 600),
    vertex_size: float = 10.0,
    edge_width_range: tuple = (0.3, 3.5),
    title: str | None = None,
) -> None:
    """Render ONE nested-SBM panel directly via graph-tool's native
    ``state.draw()`` to a PDF file (the graph-tool homepage look).

    Edges below the percolation threshold (min-weight edge of the maximum
    spanning tree, scaled by ``threshold_scale``) are hidden.  SBM
    inference runs on the FULL graph; thresholding is for display only.

    Compose multiple panels into a grid with pdfjam after rendering:
        pdfjam --nup 3x2 panel_*.pdf -o fig_E2_sbm.pdf
    """
    import graph_tool.all as gt
    import gc

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    gt.seed_rng(seed)
    N = A.shape[0]

    g = gt.Graph(directed=False)
    g.add_vertex(N)
    weight = g.new_edge_property("double")
    ri, ci = np.triu_indices(N, k=1)
    vals = A[ri, ci]
    mask = vals > 0
    elist = np.column_stack([ri[mask], ci[mask], vals[mask]])
    g.add_edge_list(elist, eprops=[weight])

    state = gt.minimize_nested_blockmodel_dl(
        g, state_args=dict(recs=[weight], rec_types=["real-exponential"])
    )

    theta = (threshold if threshold is not None
             else compute_percolation_threshold(A)) * threshold_scale
    efilter = g.new_edge_property("bool")
    for e in g.edges():
        efilter[e] = float(weight[e]) >= theta
    g.set_edge_filter(efilter)

    cross_color = [0.1, 0.1, 0.1, 0.5]
    ecolor = g.new_edge_property("vector<double>")
    # draw_order: cross-probe by weight, same-probe offset by a large
    # constant so they always paint on top of the cross-probe layer.
    draw_order = g.new_edge_property("double")
    w_max = float(max((float(weight[e]) for e in g.edges()), default=1.0))
    same_offset = w_max + 1.0
    for e in g.edges():
        u = int(e.source()); v = int(e.target())
        same = (highlight_same_probe and probe_labels is not None
                and probe_labels[u] == probe_labels[v])
        if same and shaft_colors is not None:
            c = shaft_colors[u]
            ecolor[e] = ([c[0], c[1], c[2], 0.95]
                         if len(c) == 3 else list(c))
            draw_order[e] = same_offset + float(weight[e])
        elif same:
            ecolor[e] = [0.85, 0.15, 0.15, 0.95]
            draw_order[e] = same_offset + float(weight[e])
        else:
            ecolor[e] = cross_color
            draw_order[e] = float(weight[e])

    vcolor = g.new_vertex_property("vector<double>")
    for v in g.vertices():
        if shaft_colors is not None:
            c = shaft_colors[int(v)]
            vcolor[v] = [c[0], c[1], c[2], 1.0] if len(c) == 3 else list(c)
        else:
            vcolor[v] = [0.5, 0.5, 0.5, 1.0]

    pen = gt.prop_to_size(
        weight, mi=edge_width_range[0], ma=edge_width_range[1],
        power=1, log=True,
    )

    state.draw(
        output=str(output_path),
        output_size=output_size,
        vertex_fill_color=vcolor,
        vertex_color=[0, 0, 0, 1],
        vertex_size=vertex_size,
        edge_color=ecolor,
        edge_pen_width=pen,
        eorder=draw_order,
        edge_gradient=[],
    )

    # Explicit cleanup to keep peak RSS low across many calls
    del state, g, weight, efilter, ecolor, vcolor, pen
    gc.collect()


def render_lrg_panel(
    output_path: Path,
    A: NDArray,
    lrg_result,
    probe_labels: list[str] | None = None,
    shaft_colors: list | None = None,
    *,
    seed: int = LAYOUT_SEED_DEFAULT,
    threshold: float | None = None,
    threshold_scale: float = 1.0,
    highlight_same_probe: bool = True,
    output_size: tuple = (600, 600),
    vertex_size: float = 10.0,
    edge_width_range: tuple = (0.03, 2.5),
    edge_alpha_range: tuple = (0.10, 0.95),
    edge_width_power: float = 4.0,
    beta: float = 0.8,
    k_clusters: int = 30,
) -> None:
    """Render ONE panel via graph-tool using the LRG dendrogram as the
    hierarchy tree (no SBM fit).

    Nodes are placed on the outer ring in dendrogram order via
    ``radial_tree_layout`` on the linkage-derived tree — same-cluster
    nodes land angularly adjacent by construction. Edges curve through
    the LRG hierarchy via ``get_hierarchy_control_points``: the graph-tool
    "homepage look" driven by LRG clustering instead of SBM.

    ``lrg_result`` must be an ``LRGResult`` (from
    :func:`lrg_eegfc.workflow.lrg.load_lrg_result`) — its
    ``linkage_matrix`` attribute is scipy-format ``Z`` built from the
    canonical ultrametric distance ``Trho = 1/rho``.

    Implementation note: graph-tool's ``get_hierarchy_control_points``
    only supports trees of depth <= 2, so we can't pass the full scipy
    dendrogram directly.  Instead we cut the LRG dendrogram at
    ``k_clusters`` communities and build a depth-2 tree: root -> one
    internal node per cluster -> the leaves of that cluster, inserted
    in scipy's dendrogram-leaf order.  Radial layout then places nodes
    on the outer ring sorted by LRG cluster, and bezier edges bundle
    through the cluster-center internal nodes.
    """
    import graph_tool.all as gt
    import gc
    from scipy.cluster.hierarchy import fcluster, leaves_list

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    gt.seed_rng(seed)

    Z = np.asarray(lrg_result.linkage_matrix)
    # Linkage was fitted on the LRG giant component (M leaves).
    M = int(Z.shape[0]) + 1
    # Restrict A + probe metadata to the giant component if necessary.
    if M != A.shape[0]:
        A = A[:M, :M]
        if probe_labels is not None:
            probe_labels = list(probe_labels)[:M]
        if shaft_colors is not None:
            shaft_colors = list(shaft_colors)[:M]
    N = M

    # Build graph from |FC|.
    g = gt.Graph(directed=False)
    g.add_vertex(N)
    weight = g.new_edge_property("double")
    ri, ci = np.triu_indices(N, k=1)
    vals = A[ri, ci]
    mask = vals > 0
    elist = np.column_stack([ri[mask], ci[mask], vals[mask]])
    g.add_edge_list(elist, eprops=[weight])

    # Cut dendrogram at k communities and order leaves per scipy dendro order.
    k = int(min(max(2, k_clusters), N - 1))
    labels = fcluster(Z, t=k, criterion="maxclust")
    dendro_order = leaves_list(Z)
    # Collect clusters in first-appearance order (preserves dendrogram layout);
    # leaves within each cluster preserve scipy's dendrogram order too.
    seen = {}
    cluster_leaves: dict = {}
    for leaf in dendro_order:
        lab = int(labels[leaf])
        if lab not in seen:
            seen[lab] = len(seen)
            cluster_leaves[lab] = []
        cluster_leaves[lab].append(int(leaf))
    ordered_labels = sorted(cluster_leaves.keys(), key=lambda x: seen[x])

    # Depth-2 tree: root (vid 2N-1+n_cl) -> cluster internals (vid N..N+n_cl-1)
    # -> leaves (vid 0..N-1). Insertion order determines radial angular order.
    n_cl = len(ordered_labels)
    t = gt.Graph(directed=True)
    t.add_vertex(N + n_cl + 1)
    root_vid = N + n_cl
    for j, lab in enumerate(ordered_labels):
        cluster_vid = N + j
        t.add_edge(t.vertex(root_vid), t.vertex(cluster_vid))
        for leaf in cluster_leaves[lab]:
            t.add_edge(t.vertex(cluster_vid), t.vertex(leaf))
    root = t.vertex(root_vid)

    tpos = gt.radial_tree_layout(t, root)
    cts = gt.get_hierarchy_control_points(g, t, tpos, beta=beta)

    # Leaf positions (first N tree vertices correspond to g's vertices).
    pos = g.new_vertex_property("vector<double>")
    for v in g.vertices():
        tv = t.vertex(int(v))
        pos[v] = [float(tpos[tv][0]), float(tpos[tv][1])]

    # Percolation threshold applies ONLY to cross-probe edges.  Same-probe
    # edges are kept in full (they're the signal we want to see) and their
    # thickness is already weight-scaled via prop_to_size below.  The
    # threshold is computed on the cross-probe subgraph so it reflects
    # cross-probe connectivity, not total graph connectivity.
    if highlight_same_probe and probe_labels is not None:
        A_cross = A.copy()
        for i in range(N):
            for j in range(i + 1, N):
                if probe_labels[i] == probe_labels[j]:
                    A_cross[i, j] = 0.0
                    A_cross[j, i] = 0.0
    else:
        A_cross = A
    theta = (threshold if threshold is not None
             else compute_percolation_threshold(A_cross)) * threshold_scale
    efilter = g.new_edge_property("bool")
    for e in g.edges():
        u = int(e.source()); v = int(e.target())
        same = (highlight_same_probe and probe_labels is not None
                and probe_labels[u] == probe_labels[v])
        if same:
            efilter[e] = True
        else:
            efilter[e] = float(weight[e]) >= theta
    g.set_edge_filter(efilter)

    # Per-edge weight^power in [0, 1] drives both width and alpha so most
    # edges are simultaneously thin AND transparent — only the heaviest
    # few percent stand out.  Normalise to the GLOBAL 95th-percentile max
    # so the absolute weight ordering is preserved (heavy same-probe
    # edges dominate over medium cross-probe ones) while a lone outlier
    # doesn't crush the entire distribution into t ≈ 0.
    edge_list = list(g.edges())
    same_flags = np.array([
        (highlight_same_probe and probe_labels is not None
         and probe_labels[int(e.source())] == probe_labels[int(e.target())])
        for e in edge_list
    ], dtype=bool)
    w_vals = np.array([float(weight[e]) for e in edge_list])
    if w_vals.size:
        w_hi = float(np.percentile(w_vals, 95))
        t_vals = (np.clip(w_vals / w_hi, 0.0, 1.0) ** edge_width_power
                  if w_hi > 0 else np.zeros_like(w_vals))
    else:
        t_vals = np.zeros(0)
    a_mi, a_ma = float(edge_alpha_range[0]), float(edge_alpha_range[1])
    alphas = a_mi + t_vals * (a_ma - a_mi)
    w_mi, w_ma = float(edge_width_range[0]), float(edge_width_range[1])
    pen_vals = w_mi + t_vals * (w_ma - w_mi)

    # Edge coloring: same-probe = shaft color, cross-probe = dark gray.
    # draw_order ensures same-probe edges always paint on top of cross-probe.
    ecolor = g.new_edge_property("vector<double>")
    draw_order = g.new_edge_property("double")
    w_max = float(w_vals.max()) if w_vals.size else 1.0
    same_offset = w_max + 1.0
    for idx, e in enumerate(edge_list):
        u = int(e.source()); v = int(e.target())
        a = float(alphas[idx])
        if same_flags[idx] and shaft_colors is not None:
            c = shaft_colors[u]
            ecolor[e] = [c[0], c[1], c[2], a] if len(c) == 3 else list(c[:3]) + [a]
            draw_order[e] = same_offset + float(weight[e])
        elif same_flags[idx]:
            ecolor[e] = [0.85, 0.15, 0.15, a]
            draw_order[e] = same_offset + float(weight[e])
        else:
            ecolor[e] = [0.1, 0.1, 0.1, a]
            draw_order[e] = float(weight[e])

    vcolor = g.new_vertex_property("vector<double>")
    for v in g.vertices():
        if shaft_colors is not None:
            c = shaft_colors[int(v)]
            vcolor[v] = [c[0], c[1], c[2], 1.0] if len(c) == 3 else list(c)
        else:
            vcolor[v] = [0.5, 0.5, 0.5, 1.0]

    # Width scaling uses the same per-class compressed t as alpha (via
    # pen_vals computed above) so the two visual channels stay in sync.
    pen = g.new_edge_property("double")
    for idx, e in enumerate(edge_list):
        pen[e] = float(pen_vals[idx])

    gt.graph_draw(
        g, pos=pos,
        output=str(output_path),
        output_size=output_size,
        vertex_fill_color=vcolor,
        vertex_color=[0, 0, 0, 1],
        vertex_size=vertex_size,
        edge_color=ecolor,
        edge_pen_width=pen,
        eorder=draw_order,
        edge_gradient=[],
        edge_control_points=cts,
    )

    del g, weight, efilter, ecolor, vcolor, pen, t, tpos, cts, pos
    gc.collect()
