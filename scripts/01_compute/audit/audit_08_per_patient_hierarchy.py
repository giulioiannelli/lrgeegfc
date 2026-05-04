#!/usr/bin/env python3
"""Audit Step 8 — per-(patient, band) cross-phase dendrogram figures with
CBR-colored branches AND leaves.

For each (patient, band), render a single figure with the four phase
dendrograms (rest_pre, task_learn, task_test, rest_post) side by side.
Categories are computed on the rest_post dendrogram via the gallery's
strict CBR classifier (verbatim) and projected per-leaf with priority
trace > persist > reset > rearrange (default rearrange).

Branches are colored by the rule:
  * if every leaf below an internal node belongs to the same category,
    the link gets that category's color;
  * else the link is background gray (mixed).
Leaf stems are over-drawn in the leaf's own category color so isolated
classified leaves remain visible even when their immediate parent merges
them with leaves of a different category.

Output structure:
  data/audit/per_patient_hierarchy/{Pat_XX}/{band}_hierarchy-crossphase.pdf
  data/audit/per_patient_hierarchy/leaf_assignment.csv
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, to_tree

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE, PHASE_LABELS,
)
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE
from lrg_eegfc.utils.metrics import jaccard_leafsets


OUT_DIR = ROOT / "data" / "audit" / "per_patient_hierarchy"
BAND_ORDER = list(BRAIN_BANDS.keys())   # 6 bands
PHASE_ORDER = list(PHASE_LABELS)        # 4 phases
PHASE_TEX = {
    "rest_pre":   r"rest$_\mathrm{pre}$",
    "task_learn": r"task$_\mathrm{learn}$",
    "task_test":  r"task$_\mathrm{test}$",
    "rest_post":  r"rest$_\mathrm{post}$",
}

# Same gallery thresholds — keep classifier identical to gallery_cbr_modules.py
SIZE_MIN = 5
SIZE_MAX = 60

PATTERN_COLOR = {
    "persist":   "#000000",   # black
    "trace":     "#d62728",   # red
    "reset":     "#7f7f7f",   # gray
    "rearrange": "#1f77b4",   # blue
}
PATTERN_PRIORITY = ["trace", "persist", "reset", "rearrange"]
BACKGROUND_GRAY = "#dddddd"
STRIP_BG = "#f5f5f5"


# ── tree utilities (verbatim from gallery_cbr_modules.py) ──────────────────

def nodes_with_leaves(Z: np.ndarray) -> list[dict]:
    n_leaves = Z.shape[0] + 1
    n_int = Z.shape[0]
    root = to_tree(Z, rd=False)
    nodes: list[dict | None] = [None] * n_int

    def walk(nd) -> frozenset[int]:
        if nd.is_leaf():
            return frozenset([nd.id])
        left = walk(nd.get_left())
        right = walk(nd.get_right())
        leaves = left | right
        row = nd.id - n_leaves
        nodes[row] = {"row": row, "id": nd.id, "h": float(nd.dist),
                       "size": len(leaves), "leaves": leaves}
        return leaves

    walk(root)
    return [nd for nd in nodes if nd is not None]


def best_jaccard(leafset: frozenset[int],
                  ref_nodes: list[dict]) -> float:
    best = 0.0
    for rn in ref_nodes:
        j = jaccard_leafsets(leafset, rn["leaves"])
        if j > best:
            best = j
    return best


def classify(J_rpre: float, J_tl: float, J_tt: float) -> str | None:
    """Verbatim from gallery_cbr_modules.py — strict mutually-exclusive bins."""
    if J_rpre >= 0.80 and J_tl >= 0.80 and J_tt >= 0.80:
        return "persist"
    if J_rpre <= 0.50 and J_tl >= 0.75 and J_tt >= 0.75:
        return "trace"
    if J_rpre >= 0.75 and min(J_tl, J_tt) <= 0.50:
        return "reset"
    if J_rpre <= 0.40 and J_tl <= 0.40 and J_tt <= 0.40:
        return "rearrange"
    return None


# ── per-leaf assignment ────────────────────────────────────────────────────

def leaf_categories(Zs: dict[str, np.ndarray]) -> np.ndarray:
    nodes_rp = nodes_with_leaves(Zs["rest_post"])
    nodes_rpre = nodes_with_leaves(Zs["rest_pre"])
    nodes_tl = nodes_with_leaves(Zs["task_learn"])
    nodes_tt = nodes_with_leaves(Zs["task_test"]) if Zs["task_test"] is not None else None

    n = Zs["rest_post"].shape[0] + 1
    leaf_patterns: list[set[str]] = [set() for _ in range(n)]

    for η in nodes_rp:
        if η["size"] < SIZE_MIN or η["size"] > SIZE_MAX:
            continue
        J_rpre = best_jaccard(η["leaves"], nodes_rpre)
        J_tl = best_jaccard(η["leaves"], nodes_tl)
        J_tt = best_jaccard(η["leaves"], nodes_tt) if nodes_tt is not None else J_tl
        cls = classify(J_rpre, J_tl, J_tt)
        if cls is None:
            continue
        for lid in η["leaves"]:
            leaf_patterns[lid].add(cls)

    out = np.empty(n, dtype=object)
    for lid in range(n):
        chosen = "rearrange"
        for p in PATTERN_PRIORITY:
            if p in leaf_patterns[lid]:
                chosen = p
                break
        out[lid] = chosen
    return out


# ── dendrogram coloring ────────────────────────────────────────────────────

def leaves_under_each_row(Z: np.ndarray) -> list[frozenset[int]]:
    """For each internal-node row r, return frozenset of leaf ids below it."""
    n_leaves = Z.shape[0] + 1
    rows_leaves: list[frozenset[int]] = [frozenset()] * Z.shape[0]
    for r in range(Z.shape[0]):
        out = set()
        for col in (0, 1):
            child = int(Z[r, col])
            if child < n_leaves:
                out.add(child)
            else:
                out |= rows_leaves[child - n_leaves]
        rows_leaves[r] = frozenset(out)
    return rows_leaves


def link_color_map(Z: np.ndarray, leaf_cats: np.ndarray) -> dict[int, str]:
    """Map scipy node id -> color. Internal nodes whose leaves all share a
    category get that color; mixed nodes get BACKGROUND_GRAY."""
    n_leaves = Z.shape[0] + 1
    rows_leaves = leaves_under_each_row(Z)
    out: dict[int, str] = {}
    for r in range(Z.shape[0]):
        cats = {leaf_cats[l] for l in rows_leaves[r]}
        if len(cats) == 1:
            out[n_leaves + r] = PATTERN_COLOR[cats.pop()]
        else:
            out[n_leaves + r] = BACKGROUND_GRAY
    return out


def overdraw_leaf_stems(ax, Z: np.ndarray, dend: dict,
                          leaf_cats: np.ndarray) -> None:
    """Draw a colored vertical line for each leaf stem in its own category
    color, so isolated classified leaves are visible even when their parent
    link is gray (mixed)."""
    n_leaves = Z.shape[0] + 1
    parent_h = np.full(n_leaves, np.nan)
    for r in range(Z.shape[0]):
        h = float(Z[r, 2])
        for col in (0, 1):
            child = int(Z[r, col])
            if child < n_leaves:
                parent_h[child] = h
    leaf_order = dend["leaves"]
    for i, lid in enumerate(leaf_order):
        cat = leaf_cats[lid]
        h_top = parent_h[lid]
        if not np.isfinite(h_top):
            continue
        ax.plot([5 + 10 * i, 5 + 10 * i], [0, h_top],
                 color=PATTERN_COLOR[cat], linewidth=0.9, zorder=4,
                 solid_capstyle="butt")


# ── load helpers ──────────────────────────────────────────────────────────

def _load_Z(pat: str, phase: str, band: str) -> np.ndarray | None:
    p = IMCOH_LRG_CACHE / pat / f"{band}_{phase}_lrg_imcoh-abs.npz"
    if not p.exists():
        return None
    d = np.load(p, allow_pickle=True)
    return np.asarray(d["linkage_matrix"])


# ── plotting ──────────────────────────────────────────────────────────────

def plot_dendrogram_with_strip(ax_dend, ax_strip, Z: np.ndarray,
                                 leaf_cats: np.ndarray, title: str = "") -> None:
    lc = link_color_map(Z, leaf_cats)
    dend = dendrogram(
        Z, ax=ax_dend, no_labels=True, color_threshold=0,
        above_threshold_color=BACKGROUND_GRAY,
        link_color_func=lambda k: lc.get(k, BACKGROUND_GRAY),
    )
    overdraw_leaf_stems(ax_dend, Z, dend, leaf_cats)
    if title:
        ax_dend.set_title(title, fontsize=11)
    ax_dend.set_xticks([])
    ax_dend.tick_params(axis="y", labelsize=8)

    leaf_order = dend["leaves"]
    n = len(leaf_order)
    ax_strip.set_facecolor(STRIP_BG)
    ax_strip.set_xlim(0, 10 * n)
    ax_strip.set_ylim(0, 1)
    for i, lid in enumerate(leaf_order):
        cat = leaf_cats[lid]
        ax_strip.axvspan(10 * i, 10 * (i + 1),
                          color=PATTERN_COLOR[cat], alpha=0.95, lw=0)
    ax_strip.set_xticks([])
    ax_strip.set_yticks([])
    for spine in ax_strip.spines.values():
        spine.set_visible(False)


def figure_for_patient_band(pat: str, band: str) -> tuple[Path | None, list[dict]]:
    Zs = {ph: _load_Z(pat, ph, band) for ph in PHASE_ORDER}
    if Zs["rest_post"] is None or Zs["rest_pre"] is None or Zs["task_learn"] is None:
        return None, []

    cats = leaf_categories(Zs)

    rows_csv = [{"patient": pat, "band": band, "leaf_id": int(lid),
                  "category": cat} for lid, cat in enumerate(cats)]

    n_phases = len(PHASE_ORDER)
    fig = plt.figure(figsize=(4.5 * n_phases, 5.4))
    outer = fig.add_gridspec(1, n_phases, wspace=0.18, top=0.84, bottom=0.05)

    for ip, phase in enumerate(PHASE_ORDER):
        cell = outer[0, ip].subgridspec(2, 1, height_ratios=[8, 1], hspace=0.04)
        ax_dend = fig.add_subplot(cell[0])
        ax_strip = fig.add_subplot(cell[1])
        Z = Zs[phase]
        if Z is None:
            ax_dend.text(0.5, 0.5, "missing", ha="center", va="center",
                          transform=ax_dend.transAxes, color="#999")
            ax_dend.set_xticks([]); ax_dend.set_yticks([])
            ax_strip.axis("off")
            ax_dend.set_title(PHASE_TEX[phase], fontsize=11)
            if ip == 0:
                ax_dend.set_ylabel(BRAIN_BAND_TEX_DICT.get(band, band),
                                     fontsize=12)
            continue
        plot_dendrogram_with_strip(ax_dend, ax_strip, Z, cats,
                                     title=PHASE_TEX[phase])
        if ip == 0:
            ax_dend.set_ylabel(BRAIN_BAND_TEX_DICT.get(band, band),
                                 fontsize=12)

    handles = [Patch(color=PATTERN_COLOR[p], label=p.upper())
                for p in PATTERN_PRIORITY]
    fig.legend(handles=handles, loc="upper center", ncol=4,
                bbox_to_anchor=(0.5, 0.98), fontsize=10, frameon=False,
                title=f"{pat} — {BRAIN_BAND_TEX_DICT.get(band, band)} — "
                      f"leaf categorisation (rsPost-anchored, size ∈ [5, 60])",
                title_fontsize=11)

    out_dir = OUT_DIR / pat
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{band}_hierarchy-crossphase.pdf"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path, rows_csv


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict] = []
    n_written = 0
    for pat in PATIENTS_4PHASE:
        print(f"{pat} ...", flush=True)
        for band in BAND_ORDER:
            out, rows = figure_for_patient_band(pat, band)
            if out is None:
                print(f"  {band}: missing phase, skipped")
                continue
            all_rows.extend(rows)
            n_written += 1
            print(f"  wrote {out.relative_to(ROOT)}")

    csv_path = OUT_DIR / "leaf_assignment.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["patient", "band", "leaf_id", "category"]
        )
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nWrote {csv_path.relative_to(ROOT)}  ({len(all_rows)} rows, "
          f"{n_written} figures)")


if __name__ == "__main__":
    main()
