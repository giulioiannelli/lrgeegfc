#!/usr/bin/env python3
"""Audit Step 14 — Multi-Scale Partition Coherence (MSPC).

Per scope: .agents/guides/task-persistence-investigation/2026-04-26_multiscale-partition-coherence.md

Per-leaf-per-scale Jaccard of cluster-mate sets across phases. Two
visualisations per (patient, band):

  (1) {band}_hierarchy-crossphase.pdf — audit_08-compatible 4 phase
       dendrograms, per-leaf colour = K-averaged dominant pattern.
  (2) {band}_multiscale-strip.pdf     — N × K_max image of dominant
       pattern at every (leaf, scale), keyed to rsPost leaf order.

Outputs:
  data/audit/per_patient_hierarchy_mspc/{Pat_XX}/{band}_hierarchy-crossphase.pdf
  data/audit/per_patient_hierarchy_mspc/{Pat_XX}/{band}_multiscale-strip.pdf
  data/audit/per_patient_hierarchy_mspc/leaf_assignment.csv
  data/audit/per_patient_hierarchy_mspc/multiscale_assignment.csv
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba, ListedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, fcluster

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_08_per_patient_hierarchy import (  # noqa: E402
    BAND_ORDER, PHASE_ORDER, PHASE_TEX, BRAIN_BAND_TEX_DICT,
    PATIENTS_4PHASE, PATTERN_COLOR, BACKGROUND_GRAY, STRIP_BG,
    _load_Z,
)


OUT_DIR = ROOT / "data" / "audit" / "per_patient_hierarchy_mspc"
PATTERNS = ["trace", "persist", "reset", "rearrange"]
PATTERN_INT = {p: i for i, p in enumerate(PATTERNS)}
CORNERS = {
    "trace":     (0.0, 1.0),
    "persist":   (1.0, 1.0),
    "reset":     (1.0, 0.0),
    "rearrange": (0.0, 0.0),
}
_DIAG = float(np.sqrt(2.0))


def _affinity(jpp: np.ndarray, jtp: np.ndarray,
                corner: tuple[float, float]) -> np.ndarray:
    cx, cy = corner
    d = np.hypot(jpp - cx, jtp - cy)
    return np.maximum(0.0, 1.0 - d / _DIAG)


def _row_jaccard_same(labels_a: np.ndarray, labels_b: np.ndarray) -> np.ndarray:
    """For each leaf ℓ, return Jaccard of {j: a[j]==a[ℓ]}\\{ℓ} and
    {j: b[j]==b[ℓ]}\\{ℓ}. Returns array shape (N,)."""
    same_a = labels_a[:, None] == labels_a[None, :]
    same_b = labels_b[:, None] == labels_b[None, :]
    inter = (same_a & same_b).sum(axis=1) - 1   # subtract self-match
    union = (same_a | same_b).sum(axis=1) - 1
    inter = np.maximum(inter, 0)
    union = np.maximum(union, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        J = np.where(union > 0, inter / union, 0.0)
    return J


def mspc_for_cell(Zs: dict[str, np.ndarray]
                    ) -> tuple[np.ndarray, np.ndarray, list[int]] | None:
    """Returns (A_table, dom_at_k, k_values).

    A_table: (N, 4) per-leaf K-averaged affinities in PATTERNS order.
    dom_at_k: (N, K) integer pattern-id at each (leaf, k).
    k_values: list of k values used.
    """
    if Zs["rest_post"] is None or Zs["rest_pre"] is None or Zs["task_learn"] is None:
        return None
    has_tt = Zs["task_test"] is not None
    Z_rpost = Zs["rest_post"]
    N = Z_rpost.shape[0] + 1
    K_max = max(N // 2, 3)
    k_values = list(range(2, K_max + 1))
    K = len(k_values)

    A_table = np.zeros((N, 4), dtype=float)
    dom_at_k = np.zeros((N, K), dtype=np.int8)

    for ki, k in enumerate(k_values):
        lab_pre = fcluster(Zs["rest_pre"], k, criterion="maxclust")
        lab_tl = fcluster(Zs["task_learn"], k, criterion="maxclust")
        lab_rp = fcluster(Zs["rest_post"], k, criterion="maxclust")
        if has_tt:
            lab_tt = fcluster(Zs["task_test"], k, criterion="maxclust")

        J_pp = _row_jaccard_same(lab_pre, lab_rp)
        J_tl = _row_jaccard_same(lab_tl, lab_rp)
        if has_tt:
            J_tt = _row_jaccard_same(lab_tt, lab_rp)
            J_tp = (J_tl + J_tt) / 2.0
        else:
            J_tp = J_tl

        a_trace = _affinity(J_pp, J_tp, CORNERS["trace"])
        a_persist = _affinity(J_pp, J_tp, CORNERS["persist"])
        a_reset = _affinity(J_pp, J_tp, CORNERS["reset"])
        a_rearrange = _affinity(J_pp, J_tp, CORNERS["rearrange"])
        stack = np.stack([a_trace, a_persist, a_reset, a_rearrange], axis=1)  # (N, 4)
        dom_at_k[:, ki] = np.argmax(stack, axis=1).astype(np.int8)
        A_table += stack

    A_table /= K
    return A_table, dom_at_k, k_values


def _dominant_and_confidence(A_table: np.ndarray
                                ) -> tuple[list[str], np.ndarray]:
    """argmax of per-leaf affinity averages."""
    n = A_table.shape[0]
    dominant = []
    confidence = np.zeros(n)
    for ℓ in range(n):
        order = np.argsort(-A_table[ℓ])
        dominant.append(PATTERNS[order[0]])
        confidence[ℓ] = A_table[ℓ, order[0]] - A_table[ℓ, order[1]]
    return dominant, confidence


# ── plotting (figure 1: 4-phase dendrograms) ──────────────────────────────

def _link_colors_from_dominant(Z: np.ndarray, dominant: list[str]) -> dict[int, str]:
    n_leaves = Z.shape[0] + 1
    rows_leaves: list[set[int]] = [set() for _ in range(Z.shape[0])]
    for r in range(Z.shape[0]):
        out = set()
        for col in (0, 1):
            child = int(Z[r, col])
            if child < n_leaves:
                out.add(child)
            else:
                out |= rows_leaves[child - n_leaves]
        rows_leaves[r] = out
    out_map: dict[int, str] = {}
    for r in range(Z.shape[0]):
        labels = {dominant[l] for l in rows_leaves[r]}
        if len(labels) == 1:
            out_map[n_leaves + r] = PATTERN_COLOR[labels.pop()]
        else:
            out_map[n_leaves + r] = BACKGROUND_GRAY
    return out_map


def _plot_phase_panel(ax_dend, ax_strip, Z: np.ndarray,
                        dominant: list[str], title: str = "") -> list[int]:
    if Z is None:
        ax_dend.text(0.5, 0.5, "missing", ha="center", va="center",
                      transform=ax_dend.transAxes, color="#999")
        ax_dend.set_xticks([]); ax_dend.set_yticks([])
        ax_strip.axis("off")
        if title:
            ax_dend.set_title(title, fontsize=11)
        return []
    lc = _link_colors_from_dominant(Z, dominant)
    dend = dendrogram(
        Z, ax=ax_dend, no_labels=True, color_threshold=0,
        above_threshold_color=BACKGROUND_GRAY,
        link_color_func=lambda k: lc.get(k, BACKGROUND_GRAY),
    )
    if title:
        ax_dend.set_title(title, fontsize=11)
    ax_dend.set_xticks([])
    ax_dend.tick_params(axis="y", labelsize=8)

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
        h_top = parent_h[lid]
        if not np.isfinite(h_top):
            continue
        ax_dend.plot([5 + 10 * i, 5 + 10 * i], [0, h_top],
                       color=PATTERN_COLOR[dominant[lid]],
                       linewidth=0.9, zorder=4, solid_capstyle="butt")
    ax_strip.set_facecolor(STRIP_BG)
    n = len(leaf_order)
    ax_strip.set_xlim(0, 10 * n)
    ax_strip.set_ylim(0, 1)
    for i, lid in enumerate(leaf_order):
        ax_strip.axvspan(10 * i, 10 * (i + 1),
                          color=PATTERN_COLOR[dominant[lid]], lw=0)
    ax_strip.set_xticks([]); ax_strip.set_yticks([])
    for spine in ax_strip.spines.values():
        spine.set_visible(False)
    return leaf_order


def figure_dendrograms(pat: str, band: str,
                          Zs: dict[str, np.ndarray],
                          dominant: list[str]) -> tuple[Path, list[int]]:
    n_phases = len(PHASE_ORDER)
    fig = plt.figure(figsize=(4.5 * n_phases, 5.4))
    outer = fig.add_gridspec(1, n_phases, wspace=0.18, top=0.84, bottom=0.05)
    rpost_order: list[int] = []
    for ip, phase in enumerate(PHASE_ORDER):
        cell = outer[0, ip].subgridspec(2, 1, height_ratios=[8, 1], hspace=0.04)
        ax_dend = fig.add_subplot(cell[0])
        ax_strip = fig.add_subplot(cell[1])
        order = _plot_phase_panel(ax_dend, ax_strip, Zs[phase], dominant,
                                     title=PHASE_TEX[phase])
        if phase == "rest_post":
            rpost_order = order
        if ip == 0:
            ax_dend.set_ylabel(BRAIN_BAND_TEX_DICT.get(band, band), fontsize=12)

    handles = [Patch(color=PATTERN_COLOR[p], label=p.upper()) for p in PATTERNS]
    fig.legend(handles=handles, loc="upper center", ncol=4,
                bbox_to_anchor=(0.5, 0.98), fontsize=10, frameon=False,
                title=f"{pat} — {BRAIN_BAND_TEX_DICT.get(band, band)} — "
                      f"MSPC (multi-scale partition coherence; K-averaged dominant)",
                title_fontsize=11)
    out_dir = OUT_DIR / pat
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{band}_hierarchy-crossphase.pdf"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path, rpost_order


# ── plotting (figure 2: multiscale strip) ─────────────────────────────────

def figure_multiscale_strip(pat: str, band: str,
                                dom_at_k: np.ndarray,
                                k_values: list[int],
                                rpost_order: list[int],
                                dominant: list[str]) -> Path:
    """Image: leaves (y, ordered by rsPost dendrogram) × scale (x) ×
    dominant-at-(leaf,k) colour."""
    if not rpost_order:
        rpost_order = list(range(dom_at_k.shape[0]))
    # Reorder rows to match rsPost dendrogram display order
    matrix = dom_at_k[np.array(rpost_order), :]
    # Discrete colormap aligned to PATTERN_INT
    cmap = ListedColormap([PATTERN_COLOR[p] for p in PATTERNS])

    fig, axes = plt.subplots(
        1, 2, figsize=(11, 6),
        gridspec_kw={"width_ratios": [16, 1], "wspace": 0.05},
    )
    ax_main, ax_codebar = axes

    ax_main.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=3,
                    interpolation="nearest",
                    extent=[k_values[0] - 0.5, k_values[-1] + 0.5,
                            len(rpost_order) - 0.5, -0.5])
    ax_main.set_xlabel("scale k", fontsize=10)
    ax_main.set_ylabel("leaf (rsPost dendrogram order)", fontsize=10)
    ax_main.tick_params(labelsize=8)

    # Right codebar = K-averaged dominant per leaf (in same y order)
    cb = np.array([PATTERN_INT[dominant[ℓ]] for ℓ in rpost_order],
                    dtype=np.int8).reshape(-1, 1)
    ax_codebar.imshow(cb, aspect="auto", cmap=cmap, vmin=0, vmax=3,
                       interpolation="nearest")
    ax_codebar.set_xticks([]); ax_codebar.set_yticks([])
    ax_codebar.set_title("K-avg", fontsize=8)

    handles = [Patch(color=PATTERN_COLOR[p], label=p.upper()) for p in PATTERNS]
    fig.legend(handles=handles, loc="upper center", ncol=4,
                bbox_to_anchor=(0.5, 1.0), fontsize=10, frameon=False,
                title=f"{pat} — {BRAIN_BAND_TEX_DICT.get(band, band)} — "
                      f"MSPC multiscale strip "
                      f"(leaves × scales; rows ordered by rsPost dendrogram)",
                title_fontsize=11)

    out_dir = OUT_DIR / pat
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{band}_multiscale-strip.pdf"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ── orchestration ────────────────────────────────────────────────────────

def figure_for(pat: str, band: str
                ) -> tuple[Path | None, list[dict], list[dict]]:
    Zs = {ph: _load_Z(pat, ph, band) for ph in PHASE_ORDER}
    out = mspc_for_cell(Zs)
    if out is None:
        return None, [], []
    A_table, dom_at_k, k_values = out
    n = A_table.shape[0]
    dominant, confidence = _dominant_and_confidence(A_table)

    leaf_rows = [{
        "patient": pat, "band": band, "leaf_id": int(ℓ),
        "dominant": dominant[ℓ], "confidence": float(confidence[ℓ]),
        "A_trace": float(A_table[ℓ, 0]), "A_persist": float(A_table[ℓ, 1]),
        "A_reset": float(A_table[ℓ, 2]), "A_rearrange": float(A_table[ℓ, 3]),
    } for ℓ in range(n)]

    multiscale_rows = [{
        "patient": pat, "band": band, "leaf_id": int(ℓ),
        "k": int(k_values[ki]),
        "dominant_at_k": PATTERNS[int(dom_at_k[ℓ, ki])],
    } for ℓ in range(n) for ki in range(len(k_values))]

    out_path, rpost_order = figure_dendrograms(pat, band, Zs, dominant)
    figure_multiscale_strip(pat, band, dom_at_k, k_values, rpost_order, dominant)
    return out_path, leaf_rows, multiscale_rows


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_leaves: list[dict] = []
    all_ms: list[dict] = []
    skipped: list[str] = []
    n_written = 0
    for pat in PATIENTS_4PHASE:
        print(f"{pat} ...", flush=True)
        for band in BAND_ORDER:
            out, leaves, ms = figure_for(pat, band)
            if out is None:
                skipped.append(f"{pat}/{band}")
                continue
            all_leaves.extend(leaves)
            all_ms.extend(ms)
            n_written += 1

    # Per-leaf CSV
    leaf_csv = OUT_DIR / "leaf_assignment.csv"
    with leaf_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "patient", "band", "leaf_id", "dominant", "confidence",
            "A_trace", "A_persist", "A_reset", "A_rearrange",
        ])
        w.writeheader(); w.writerows(all_leaves)

    # Multiscale CSV (large — leaf × k rows)
    ms_csv = OUT_DIR / "multiscale_assignment.csv"
    with ms_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["patient", "band", "leaf_id",
                                              "k", "dominant_at_k"])
        w.writeheader(); w.writerows(all_ms)

    if skipped:
        (OUT_DIR / "missing.txt").write_text("\n".join(skipped) + "\n")
    print(f"\n{n_written} (patient, band) cells, "
          f"{len(skipped)} skipped, "
          f"{len(all_leaves)} leaf rows, {len(all_ms)} multiscale rows")


if __name__ == "__main__":
    main()
