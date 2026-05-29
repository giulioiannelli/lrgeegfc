"""Shared header + dendrogram helpers for the audit_round3_section5_*
per-figure scripts.

Split out of audit_round3_section5_redo.py on 2026-05-29 (Phase 4-B split
2/7). Holds the import block, output-directory constants, the asterisks
helper, and the dendrogram-helper functions used across the KC redos.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram

# Centralized font sizing — change here, propagates everywhere.
# (Per-call fontsize=... overrides still take precedence if needed.)
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
})

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PATIENTS_LIST,
)
from lrg_eegfc.utils.metrics.tree_distance import kc_vectors
from lrg_eegfc.workflow.lrg import load_lrg_result

S5 = ROOT / "data" / "reports" / "section_5_lrg_trace"
R2 = ROOT / "data" / "audit" / "section5_v2_round2"
R2_TBL = R2 / "tables"

OUT = ROOT / "data" / "audit" / "section5_v2_round3_redo"
FIG = OUT / "figures"
TBL = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TBL.mkdir(parents=True, exist_ok=True)

BAND_ORDER = list(BRAIN_BANDS_NAMES)
PHASES = ("rest_pre", "task_test", "rest_post")

# --- redo-1 dendrogram leaf-label placement -------------------------------
# Vertical position of each leaf label as a fraction of its parent merge
# height (0 < value < 1). Closer to 1.0 → closer to the merge bracket.
DEND_LABEL_Y_FRAC = 0.9
DEND_LABEL_FONTSIZE = 4.0
# Number of label-heights (visual, display-pixel space) by which the
# alternating "low-shelf" label sits below the high shelf. 1.0 = exactly
# one label-height; 1.5 = one and a half label-heights, etc.
DEND_LABEL_DROP_HEIGHTS = 2.0


def asterisks(p: float) -> str:
    if not np.isfinite(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""

# ===================================================================
# Redo 7 — KC dendrogram eye-proof (Pat_06, β trace vs θ null)
# ===================================================================
def _reorder_linkage_to_match(Z: np.ndarray, target_leaves: list[int]) -> np.ndarray:
    """Greedy permutation of merge children so the natural dendrogram leaf
    order is as close to ``target_leaves`` as possible. Topology preserved
    (only left/right swaps at each merge); KC, MC, RF distances unchanged.
    """
    target_pos = {leaf: i for i, leaf in enumerate(target_leaves)}
    n = Z.shape[0] + 1
    leaves_under: dict[int, list[int]] = {i: [i] for i in range(n)}
    Z_perm = Z.copy()
    for k in range(n - 1):
        a = int(Z[k, 0])
        b = int(Z[k, 1])
        cid = n + k
        leaves_a = leaves_under[a]
        leaves_b = leaves_under[b]
        avg_a = np.mean([target_pos.get(l, n) for l in leaves_a])
        avg_b = np.mean([target_pos.get(l, n) for l in leaves_b])
        if avg_a > avg_b:
            Z_perm[k, 0] = b
            Z_perm[k, 1] = a
            leaves_under[cid] = leaves_b + leaves_a
        else:
            leaves_under[cid] = leaves_a + leaves_b
    return Z_perm

# ===================================================================
# Redo 8 — Per-patient KC-clade dendrogram visualization (Pat_06)
# ===================================================================
def _enumerate_clades(Z: np.ndarray) -> list[dict]:
    """For every internal node of a scipy linkage Z, return
    [{'cid': int, 'leaves': frozenset[int], 'height': float}, ...]
    sorted by node id (root last). Length = N-1.
    """
    n = Z.shape[0] + 1
    leaves_under: dict[int, frozenset[int]] = {i: frozenset({i}) for i in range(n)}
    out = []
    for k in range(n - 1):
        a, b = int(Z[k, 0]), int(Z[k, 1])
        cid = n + k
        leaves_under[cid] = leaves_under[a] | leaves_under[b]
        out.append(dict(cid=cid, leaves=leaves_under[cid], height=float(Z[k, 2])))
    return out


def _max_jaccard_match(
    clades_A: list[dict], clades_B: list[dict]
) -> list[dict]:
    """For each clade in A return its best-Jaccard partner in B.
    Returns list aligned with clades_A: each row has
    {'cid_A','leaves_A','height_A','cid_B','leaves_B','height_B','jaccard'}.
    """
    rows = []
    sets_B = [c["leaves"] for c in clades_B]
    for cA in clades_A:
        sA = cA["leaves"]
        best_j = 0.0
        best_idx = -1
        for idx, sB in enumerate(sets_B):
            inter = len(sA & sB)
            if inter == 0:
                continue
            uni = len(sA | sB)
            j = inter / uni
            if j > best_j:
                best_j = j
                best_idx = idx
        if best_idx < 0:
            continue
        cB = clades_B[best_idx]
        rows.append(dict(
            cid_A=cA["cid"], leaves_A=cA["leaves"], height_A=cA["height"],
            cid_B=cB["cid"], leaves_B=cB["leaves"], height_B=cB["height"],
            jaccard=best_j,
        ))
    return rows


def _clade_leaf_xrange(Z: np.ndarray, leaf_to_x: dict[int, float],
                        cid: int) -> tuple[float, float]:
    """Min/max leaf-x position covered by the subtree rooted at cid."""
    n = Z.shape[0] + 1
    if cid < n:
        return (leaf_to_x[cid], leaf_to_x[cid])
    # Walk subtree iteratively
    stack = [cid]
    leaves = []
    while stack:
        node = stack.pop()
        if node < n:
            leaves.append(leaf_to_x[node])
        else:
            kk = node - n
            stack.append(int(Z[kk, 0]))
            stack.append(int(Z[kk, 1]))
    return (min(leaves), max(leaves))


def _draw_clade(ax, Z: np.ndarray, leaf_to_x: dict[int, float],
                cid: int, color: str, lw: float, alpha: float,
                tmin: float):
    """Draw a single horizontal bar at the clade's merge height, spanning
    the leaf-x range covered by its subtree, with two short vertical whiskers
    dropping to ``tmin`` so the eye can locate the clade's leaf base.
    Far less visual mass than recursively recoloring the whole subtree.
    """
    n = Z.shape[0] + 1
    if cid < n:
        return  # leaves are not clades
    kk = cid - n
    h = float(Z[kk, 2])
    x_lo, x_hi = _clade_leaf_xrange(Z, leaf_to_x, cid)
    # Horizontal bar at merge height
    ax.plot([x_lo, x_hi], [h, h], color=color, lw=lw, alpha=alpha,
            solid_capstyle="round", zorder=5)
    # Short whiskers (drop ~25% of the log-axis range visually)
    y_whisk = max(tmin, h * 0.5)
    ax.plot([x_lo, x_lo], [y_whisk, h], color=color, lw=lw * 0.7,
            alpha=alpha * 0.85, zorder=5)
    ax.plot([x_hi, x_hi], [y_whisk, h], color=color, lw=lw * 0.7,
            alpha=alpha * 0.85, zorder=5)

# ===================================================================
# Redo 9 — Cohort per-band dendrogram view (one PDF per band)
# ===================================================================
def _best_jaccard_partner(
    anchor: dict, candidates: list[dict]
) -> tuple[dict | None, float]:
    """For one anchor clade, return (best_match_clade, best_jaccard) over a
    list of candidate clades. Returns (None, 0.0) if no overlap.
    """
    sA = anchor["leaves"]
    best_j = 0.0
    best = None
    for cB in candidates:
        sB = cB["leaves"]
        inter = len(sA & sB)
        if inter == 0:
            continue
        j = inter / len(sA | sB)
        if j > best_j:
            best_j = j
            best = cB
    return best, best_j
