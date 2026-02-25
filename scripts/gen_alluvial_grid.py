#!/usr/bin/env python3
"""Generate alluvial community diagrams (with pin labels) for all phase/band combos.

Matplotlib PDF alluvial with side-by-side packed pin names inside community
blocks. Uses fixed (τ, n) pairs for clean merging visualization.

Run: python scripts/gen_alluvial_grid.py
"""
import sys
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from scipy.optimize import linear_sum_assignment

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PHASE_LABELS
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.visuals.lrg import _load_channel_labels
from lrg_eegfc.visuals.metastable import compute_clustering_across_tau

# ── Config ────────────────────────────────────────────────────────────
PATIENT = "Pat_02"
MSC_CACHE = ROOT / "data" / "msc_cache"
DATASET_ROOT = ROOT / "data" / "stereoeeg_patients"
OUTPUT_DIR = ROOT / "data" / "figures" / "report_mslcd_section" / PATIENT / "alluvial"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

COMBOS = [(phase, band) for phase in PHASE_LABELS for band in BRAIN_BANDS_NAMES]

TAU_NCLUST_PAIRS = [
    (0.1, 15), (0.3, 12), (0.5, 8), (0.8, 6),
    (1.0, 5), (2.0, 3), (5.0, 2),
]
ALLUVIAL_TAUS = [t for t, _ in TAU_NCLUST_PAIRS]
ALLUVIAL_NCLUST = [n for _, n in TAU_NCLUST_PAIRS]

# ── Visual parameters ────────────────────────────────────────────────
PIN_FONTSIZE = 11
ROW_HEIGHT_MIN = 1.8
CHARS_PER_LINE = 22
BOX_PAD_Y = 1.0
BAR_WIDTH = 0.75

# Kelly/Boynton maximally-distinct colors
_DISTINCT_COLORS = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#d63b00", "#f032e6", "#469990", "#9A6324", "#800000",
    "#808000", "#000075", "#e03080", "#1a7820", "#7030c0",
    "#b05010", "#206060", "#c04080", "#305090", "#704020",
]


def get_distinct_colors(n):
    if n <= len(_DISTINCT_COLORS):
        return [matplotlib.colors.to_rgba(c) for c in _DISTINCT_COLORS[:n]]
    colors = []
    for i in range(n):
        hue = (i * 0.618033988749895) % 1.0
        colors.append(matplotlib.colors.hsv_to_rgb([hue, 0.75, 0.85]))
    return colors


def align_cluster_labels(labels_list):
    """Align cluster labels across scales by maximum overlap (Hungarian)."""
    aligned = [labels_list[0].copy()]
    for s in range(1, len(labels_list)):
        prev, curr = aligned[-1], labels_list[s]
        u_prev, u_curr = np.unique(prev), np.unique(curr)
        overlap = np.zeros((len(u_prev), len(u_curr)))
        for i, lp in enumerate(u_prev):
            for j, lc in enumerate(u_curr):
                overlap[i, j] = np.sum((prev == lp) & (curr == lc))
        dim = max(len(u_prev), len(u_curr))
        cost = np.zeros((dim, dim))
        cost[:len(u_prev), :len(u_curr)] = -overlap
        row_ind, col_ind = linear_sum_assignment(cost)
        label_map = {}
        for r, c in zip(row_ind, col_ind):
            if r < len(u_prev) and c < len(u_curr):
                label_map[u_curr[c]] = u_prev[r]
        used = set(label_map.values())
        next_lbl = max(max(used) + 1, max(u_prev) + 1) if used else 1
        for lc in u_curr:
            if lc not in label_map:
                label_map[lc] = next_lbl
                next_lbl += 1
        aligned.append(np.array([label_map[l] for l in curr]))
    return aligned


def _pack_pins_into_rows(pin_names, chars_per_line=CHARS_PER_LINE):
    """Pack pin names side-by-side, wrapping when the line is full."""
    rows = []
    current = ""
    for pin in pin_names:
        candidate = f"{current}  {pin}" if current else pin
        if len(candidate) > chars_per_line and current:
            rows.append(current)
            current = pin
        else:
            current = candidate
    if current:
        rows.append(current)
    return rows


def make_alluvial(phase, band, partitions, channel_labels, N, out_path):
    """Build and save matplotlib alluvial with pin labels."""
    n_scales = len(ALLUVIAL_TAUS)
    gap = N * 0.015

    communities_raw = [partitions[t] for t in ALLUVIAL_TAUS]
    communities_aligned = align_cluster_labels(communities_raw)

    all_labels = sorted(set(l for arr in communities_aligned for l in arr))
    comm_colors = {l: get_distinct_colors(len(all_labels))[i]
                   for i, l in enumerate(all_labels)}

    # Sort nodes by trajectory
    trajectories = [tuple(communities_aligned[s][i] for s in range(n_scales))
                    for i in range(N)]
    sort_order = sorted(range(N), key=lambda i: trajectories[i])
    labels_sorted = [communities_aligned[s][sort_order] for s in range(n_scales)]

    # First pass: compute block heights
    block_info = []
    for s in range(n_scales):
        labels = labels_sorted[s]
        seen, unique = set(), []
        for l in labels:
            if l not in seen:
                unique.append(l); seen.add(l)
        info = []
        for comm in unique:
            member_mask = labels == comm
            member_sorted_idx = np.where(member_mask)[0]
            pins = [channel_labels[sort_order[si]] for si in member_sorted_idx]
            rows = _pack_pins_into_rows(pins)
            height = max(len(pins), len(rows) * ROW_HEIGHT_MIN + 2 * BOX_PAD_Y)
            info.append((comm, pins, rows, height))
        block_info.append(info)

    # Second pass: positions
    block_pos = []
    for s in range(n_scales):
        positions = {}
        y = 0
        for comm, pins, rows, height in block_info[s]:
            positions[comm] = (y, y + height, rows)
            y += height + gap
        block_pos.append(positions)

    fig, ax = plt.subplots(figsize=(max(4.5 * n_scales, 16), 14))

    # Draw blocks + text
    for s in range(n_scales):
        for comm, (y0, y1, rows) in block_pos[s].items():
            height = y1 - y0
            rect = plt.Rectangle(
                (s - BAR_WIDTH / 2, y0), BAR_WIDTH, height,
                facecolor=comm_colors[comm], edgecolor="white",
                lw=1.5, zorder=3,
            )
            ax.add_patch(rect)

            n_rows = len(rows)
            usable_h = height - 2 * BOX_PAD_Y
            row_h = usable_h / max(n_rows, 1)
            y_start = y0 + BOX_PAD_Y + row_h / 2
            for r, row_text in enumerate(rows):
                ax.text(s, y_start + r * row_h, row_text,
                        ha="center", va="center",
                        fontsize=PIN_FONTSIZE, fontweight="bold",
                        color="white", zorder=4, family="monospace")

    # Bezier flows
    for s in range(n_scales - 1):
        left, right = labels_sorted[s], labels_sorted[s + 1]
        trans = Counter(zip(left, right))
        off_left = {}
        for comm, (y0, y1, _) in block_pos[s].items():
            n_mem = int(np.sum(labels_sorted[s] == comm))
            off_left[comm] = (y0, (y1 - y0) / n_mem)
        off_right = {}
        for comm, (y0, y1, _) in block_pos[s + 1].items():
            n_mem = int(np.sum(labels_sorted[s + 1] == comm))
            off_right[comm] = (y0, (y1 - y0) / n_mem)

        for (c_l, c_r), count in sorted(trans.items(),
                                         key=lambda x: (-x[0][0], -x[1])):
            cy_l, sc_l = off_left[c_l]
            y_src_bot = cy_l
            y_src_top = cy_l + count * sc_l
            off_left[c_l] = (y_src_top, sc_l)

            cy_r, sc_r = off_right[c_r]
            y_tgt_bot = cy_r
            y_tgt_top = cy_r + count * sc_r
            off_right[c_r] = (y_tgt_top, sc_r)

            x_left = s + BAR_WIDTH / 2
            x_right = (s + 1) - BAR_WIDTH / 2
            t = np.linspace(0, 1, 80)
            x_curve = x_left + t * (x_right - x_left)
            h = 3 * t**2 - 2 * t**3
            ax.fill_between(
                x_curve,
                y_src_bot + h * (y_tgt_bot - y_src_bot),
                y_src_top + h * (y_tgt_top - y_src_top),
                color=comm_colors[c_l], alpha=0.3,
                edgecolor="none", zorder=1,
            )

    ax.set_xticks(range(n_scales))
    ax.set_xticklabels(
        [rf"$\tau={t:.1f}$" + "\n" + rf"$n={n}$"
         for t, n in TAU_NCLUST_PAIRS],
        fontsize=14,
    )
    y_top = max(bp[1] for bps in block_pos for _, bp in bps.items()
                if len(bp) >= 2)
    ax.set_xlim(-0.8, n_scales - 0.2)
    ax.set_ylim(-gap, y_top + gap * 2)
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


# ── Load channel labels (shared) ─────────────────────────────────────
channel_labels = _load_channel_labels(PATIENT, DATASET_ROOT)
N = len(channel_labels)
print(f"{PATIENT}: {N} channels, {len(COMBOS)} combos to generate\n")

# ── Generate ──────────────────────────────────────────────────────────
for phase, band in COMBOS:
    print(f"  {phase}/{band}...", end=" ")

    A = load_msc_matrix(PATIENT, phase, band,
                        cache_root=MSC_CACHE, sparsify="none",
                        n_surrogates=0, nperseg=4096)
    if A is None:
        print("SKIP (no MSC cache)")
        continue
    np.fill_diagonal(A, 0)
    assert A.shape[0] == N

    G = nx.from_numpy_array(A)
    Gcc_nodes = max(nx.connected_components(G), key=len)
    Gcc = G.subgraph(Gcc_nodes).copy()

    partitions, _ = compute_clustering_across_tau(
        None, np.array(ALLUVIAL_TAUS), Gcc,
        n_clusters_list=ALLUVIAL_NCLUST,
    )

    out = OUTPUT_DIR / f"fig_alluvial_{phase}_{band}.pdf"
    make_alluvial(phase, band, partitions, channel_labels, N, out)
    print(f"saved: {out.name}")

print(f"\nDone! Figures in {OUTPUT_DIR}")
