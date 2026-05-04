#!/usr/bin/env python3
"""Audit Step 12 — Cohesion-CBR (variant E).

Per scope: .agents/guides/task-persistence-investigation/2026-04-25_cohesion-cbr.md

Tightness T(S, T_φ) = |S| / |leaves(LCA(S, T_φ))| replaces both Jaccard
and containment. Four soft affinities replace the binary classifier.
Per-leaf colour = dominant affinity, alpha = top1 - top2 (confidence).

Outputs:
  data/audit/per_patient_hierarchy_cohesion/{Pat_XX}/{band}_hierarchy-crossphase.pdf
  data/audit/per_patient_hierarchy_cohesion/leaf_assignment.csv
  data/audit/per_patient_hierarchy_cohesion/anchor_classification.csv
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from matplotlib.patches import Patch
from matplotlib.colors import to_rgba
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_08_per_patient_hierarchy import (  # noqa: E402
    BAND_ORDER, PHASE_ORDER, PHASE_TEX, BRAIN_BAND_TEX_DICT,
    PATIENTS_4PHASE, PATTERN_COLOR, PATTERN_PRIORITY, BACKGROUND_GRAY,
    SIZE_MIN, SIZE_MAX, STRIP_BG,
    nodes_with_leaves, _load_Z,
)


OUT_DIR = ROOT / "data" / "audit" / "per_patient_hierarchy_cohesion"
ALPHA_CLASSIFIED = 1.00     # solid colour for every classified leaf
ALPHA_UNCLASSIFIED = 0.30   # faint gray for leaves with no eligible anchor
COH_THETA = 0.8             # robust-cohesion coverage gate (≥80 % of S in one v)
NEUTRAL_GRAY = "#bbbbbb"
PATTERNS = ["trace", "persist", "reset", "rearrange"]


# ── tightness (kept as a diagnostic sidecar; not used for classification) ──

def tightness(S: frozenset[int], ref_nodes: list[dict]) -> tuple[float, int]:
    if not S:
        return 0.0, 0
    s = len(S)
    best_size = None
    for v in ref_nodes:
        if v["size"] >= s and S.issubset(v["leaves"]):
            if best_size is None or v["size"] < best_size:
                best_size = v["size"]
    if best_size is None:
        best_size = max((v["size"] for v in ref_nodes), default=s)
    return s / best_size, best_size


# ── Jaccard ───────────────────────────────────────────────────────────────

def best_jaccard_local(S: frozenset[int], ref_nodes: list[dict]) -> float:
    best = 0.0
    for v in ref_nodes:
        inter = len(S & v["leaves"])
        if inter == 0:
            continue
        un = len(S | v["leaves"])
        j = inter / un
        if j > best:
            best = j
    return best


# ── robust cohesion (gated containment) ──────────────────────────────────

def robust_cohesion(S: frozenset[int], ref_nodes: list[dict],
                     theta: float = COH_THETA) -> float:
    """Coverage-gated purity: score `|S ∩ leaves(v)| / |leaves(v)|` over only
    those reference subtrees v that already cover ≥ θ·|S| of the anchor.

    A leafset that is split across multiple distant subtrees of T_φ will
    fail the coverage gate at every reasonable v and fall back to v = root
    (low score). A leafset that is cleanly contained in some v of similar
    size scores ≈ 1. Tolerates θ·|S| boundary drift — fixes tightness's
    one-straggler-crashes-the-LCA problem.

    θ = 0.7 default: anchor must have at least 70 % of its leaves in a
    single reference subtree to be considered cohesive there.
    """
    if not S:
        return 0.0
    s = len(S)
    threshold = theta * s
    best = 0.0
    for v in ref_nodes:
        inter = len(S & v["leaves"])
        if inter < threshold:
            continue
        score = inter / v["size"]
        if score > best:
            best = score
    return best


# ── corner-distance soft affinities ───────────────────────────────────────

_SQRT2 = float(np.sqrt(2.0))


def corner_affinities(j_pre: float, j_task: float) -> dict[str, float]:
    """Distance-from-corner affinities on the unit square (j_pre, j_task).

    Corners:
      TRACE     = (0, 1)   low rsPre Jaccard, high task Jaccard
      PERSIST   = (1, 1)   high in both
      RESET     = (1, 0)   high rsPre Jaccard, low task
      REARRANGE = (0, 0)   low in both

    a_p = 1 − euclidean((j_pre, j_task), corner_p) / sqrt(2)        ∈ [0, 1].
    """
    def _aff(cx: float, cy: float) -> float:
        d = np.hypot(j_pre - cx, j_task - cy) / _SQRT2
        return float(max(0.0, 1.0 - d))
    return {
        "trace":     _aff(0.0, 1.0),
        "persist":   _aff(1.0, 1.0),
        "reset":     _aff(1.0, 0.0),
        "rearrange": _aff(0.0, 0.0),
    }


# ── per-anchor pass ───────────────────────────────────────────────────────

def per_anchor_classification(Zs: dict[str, np.ndarray]
                                ) -> tuple[list[dict], dict[str, list[dict]]]:
    """Returns (anchor_rows, nodes_dict).

    Uses Jaccard scores `J_φ = best_jaccard(S, T_φ)` per phase; collapses the
    two task phases into `J_task = mean(J_tl, J_tt)`. Soft affinities are
    corner-distance based on the (J_rpre, J_task) point in the unit square.
    Tightness is also computed and reported as a diagnostic sidecar.
    """
    nodes = {ph: nodes_with_leaves(Z) if Z is not None else None
             for ph, Z in Zs.items()}
    if nodes["rest_post"] is None or nodes["rest_pre"] is None or nodes["task_learn"] is None:
        return [], nodes
    has_tt = nodes["task_test"] is not None

    Z_rpost = Zs["rest_post"]
    dmax_rpost = float(Z_rpost[-1, 2])
    rows: list[dict] = []
    for η in nodes["rest_post"]:
        if η["size"] < SIZE_MIN or η["size"] > SIZE_MAX:
            continue
        S = η["leaves"]

        # Robust cohesion (PRIMARY classifier input).
        COH_rpre = robust_cohesion(S, nodes["rest_pre"])
        COH_tl = robust_cohesion(S, nodes["task_learn"])
        if has_tt:
            COH_tt = robust_cohesion(S, nodes["task_test"])
        else:
            COH_tt = float("nan")
        # Task cohesion: REQUIRE both task phases to cohere (use min).
        # If either task phase has the leaves scattered, this is not a trace.
        COH_task = min(COH_tl, COH_tt) if has_tt else COH_tl

        # Jaccard sidecar (diagnostic).
        J_rpre = best_jaccard_local(S, nodes["rest_pre"])
        J_tl = best_jaccard_local(S, nodes["task_learn"])
        J_tt = best_jaccard_local(S, nodes["task_test"]) if has_tt else float("nan")

        # Tightness sidecar (diagnostic).
        T_rpre, _ = tightness(S, nodes["rest_pre"])
        T_tl, _ = tightness(S, nodes["task_learn"])
        T_tt, _ = tightness(S, nodes["task_test"]) if has_tt else (float("nan"), 0)

        aff = corner_affinities(COH_rpre, COH_task)
        rows.append({
            "anchor_row": int(η["row"]),
            "size": int(η["size"]),
            "h_rel": float(η["h"]) / dmax_rpost if dmax_rpost > 0 else 0.0,
            "COH_rpre": COH_rpre, "COH_tl": COH_tl, "COH_tt": COH_tt,
            "COH_task": COH_task,
            "J_rpre": J_rpre, "J_tl": J_tl, "J_tt": J_tt,
            "T_rpre": T_rpre, "T_tl": T_tl, "T_tt": T_tt,
            "a_trace": aff["trace"], "a_persist": aff["persist"],
            "a_reset": aff["reset"], "a_rearrange": aff["rearrange"],
            "leaves": S,
        })
    return rows, nodes


# ── per-leaf projection ───────────────────────────────────────────────────

def project_to_leaves(n: int, anchor_rows: list[dict]
                        ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-leaf projection — *smallest containing anchor* rule.

    For each leaf ℓ, find the smallest rsPost anchor (size ∈ [s_min, s_max])
    whose leafset contains ℓ. The leaf inherits that anchor's affinity
    vector, so dominant = argmax pattern of that anchor and confidence =
    top1 − top2 of that anchor.

    Rationale: averaging across ancestors of all scales washes out the
    trace signal (small trace anchors are diluted by larger containing
    anchors that classify as persist because they include rsPre-shared
    material). The smallest containing anchor is the most specific
    evidence about which subtree ℓ first forms a coherent group with —
    that's what the eye reads on the dendrogram and what the leaf's
    colour should reflect.

    Returns (dominant_label, alpha, A_table, classified_mask). A_table
    has shape (n, 4) in PATTERNS order.
    """
    # Sort anchors ascending by size so smallest-containing wins.
    sorted_anchors = sorted(anchor_rows, key=lambda r: r["size"])
    A_table = np.zeros((n, 4), dtype=float)
    classified = np.zeros(n, dtype=bool)
    for row in sorted_anchors:
        a = np.array([row[f"a_{p}"] for p in PATTERNS], dtype=float)
        for lid in row["leaves"]:
            if classified[lid]:
                continue          # smaller anchor already wrote this leaf
            A_table[lid] = a
            classified[lid] = True

    dominant = np.empty(n, dtype=object)
    alpha = np.zeros(n, dtype=float)
    for ℓ in range(n):
        if not classified[ℓ]:
            dominant[ℓ] = "unclassified"
            alpha[ℓ] = ALPHA_UNCLASSIFIED
            continue
        order = np.argsort(-A_table[ℓ])
        dominant[ℓ] = PATTERNS[order[0]]
        # Solid colour for any classified leaf — no confidence-based fading
        # (geometry of the corner-distance score systematically gives PERSIST
        #  larger gaps than the other categories, which made trace fade out).
        alpha[ℓ] = ALPHA_CLASSIFIED
    return dominant, alpha, A_table, classified


# ── plotting ──────────────────────────────────────────────────────────────

def _leaf_rgba(label: str, alpha: float) -> tuple[float, float, float, float]:
    if label == "unclassified":
        rgba = list(to_rgba(NEUTRAL_GRAY))
    else:
        rgba = list(to_rgba(PATTERN_COLOR[label]))
    rgba[3] = alpha
    return tuple(rgba)


def link_colors_for_phase(Z: np.ndarray, dominant: np.ndarray,
                            alpha: np.ndarray) -> dict[int, str]:
    """Internal-link colour as a hex string (scipy requires this).

    Per scope: dominant colour iff every leaf below shares the same dominant
    category AND mean alpha ≥ 0.5; else background gray. Per-link alpha is
    not encoded (scipy limitation); the saturation signal is carried by leaf
    stems and codebar instead.
    """
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
        leaves = rows_leaves[r]
        labels = {dominant[l] for l in leaves}
        labels.discard("unclassified")
        if len(labels) == 1:
            cat = labels.pop()
            mean_alpha = float(np.mean([alpha[l] for l in leaves]))
            if mean_alpha >= 0.5:
                out_map[n_leaves + r] = PATTERN_COLOR[cat]
                continue
        out_map[n_leaves + r] = BACKGROUND_GRAY
    return out_map


def plot_phase_panel(ax_dend, ax_strip, Z: np.ndarray,
                       dominant: np.ndarray, alpha: np.ndarray,
                       title: str = "") -> None:
    if Z is None:
        ax_dend.text(0.5, 0.5, "missing", ha="center", va="center",
                      transform=ax_dend.transAxes, color="#999")
        ax_dend.set_xticks([]); ax_dend.set_yticks([])
        ax_strip.axis("off")
        if title:
            ax_dend.set_title(title, fontsize=11)
        return

    lc = link_colors_for_phase(Z, dominant, alpha)
    dend = dendrogram(
        Z, ax=ax_dend, no_labels=True, color_threshold=0,
        above_threshold_color=BACKGROUND_GRAY,
        link_color_func=lambda k: lc.get(k, BACKGROUND_GRAY),
    )
    if title:
        ax_dend.set_title(title, fontsize=11)
    ax_dend.set_xticks([])
    ax_dend.tick_params(axis="y", labelsize=8)

    # Overdraw leaf stems with the leaf's (hue, alpha)
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
        rgba = _leaf_rgba(dominant[lid], alpha[lid])
        ax_dend.plot([5 + 10 * i, 5 + 10 * i], [0, h_top],
                       color=rgba, linewidth=0.9, zorder=4,
                       solid_capstyle="butt")

    # Codebar
    n = len(leaf_order)
    ax_strip.set_facecolor(STRIP_BG)
    ax_strip.set_xlim(0, 10 * n)
    ax_strip.set_ylim(0, 1)
    for i, lid in enumerate(leaf_order):
        rgba = _leaf_rgba(dominant[lid], alpha[lid])
        ax_strip.axvspan(10 * i, 10 * (i + 1), color=rgba, lw=0)
    ax_strip.set_xticks([])
    ax_strip.set_yticks([])
    for spine in ax_strip.spines.values():
        spine.set_visible(False)


def figure_for(pat: str, band: str
                ) -> tuple[Path | None, list[dict], list[dict]]:
    Zs = {ph: _load_Z(pat, ph, band) for ph in PHASE_ORDER}
    if Zs["rest_post"] is None or Zs["rest_pre"] is None or Zs["task_learn"] is None:
        return None, [], []
    n = Zs["rest_post"].shape[0] + 1
    anchor_rows, _ = per_anchor_classification(Zs)
    dominant, alpha, A_table, classified = project_to_leaves(n, anchor_rows)

    leaf_rows = []
    for lid in range(n):
        leaf_rows.append({
            "patient": pat, "band": band, "leaf_id": int(lid),
            "dominant": dominant[lid],
            "confidence": float((A_table[lid].max() - np.partition(A_table[lid], -2)[-2])
                                  if classified[lid] else 0.0),
            "A_trace": float(A_table[lid, 0]),
            "A_persist": float(A_table[lid, 1]),
            "A_reset": float(A_table[lid, 2]),
            "A_rearrange": float(A_table[lid, 3]),
        })
    for r in anchor_rows:
        r2 = {k: v for k, v in r.items() if k != "leaves"}
        r2["patient"] = pat; r2["band"] = band
        # T_tt may be NaN — preserve
    anchor_rows_clean = [{k: v for k, v in r.items() if k != "leaves"} | {"patient": pat, "band": band}
                            for r in anchor_rows]

    n_phases = len(PHASE_ORDER)
    fig = plt.figure(figsize=(4.5 * n_phases, 5.4))
    outer = fig.add_gridspec(1, n_phases, wspace=0.18, top=0.84, bottom=0.05)
    for ip, phase in enumerate(PHASE_ORDER):
        cell = outer[0, ip].subgridspec(2, 1, height_ratios=[8, 1], hspace=0.04)
        ax_dend = fig.add_subplot(cell[0])
        ax_strip = fig.add_subplot(cell[1])
        plot_phase_panel(ax_dend, ax_strip, Zs[phase], dominant, alpha,
                          title=PHASE_TEX[phase])
        if ip == 0:
            ax_dend.set_ylabel(BRAIN_BAND_TEX_DICT.get(band, band), fontsize=12)

    handles = [Patch(color=PATTERN_COLOR[p], label=p.upper()) for p in PATTERNS]
    handles.append(Patch(color=NEUTRAL_GRAY, label="UNCLASSIFIED"))
    fig.legend(handles=handles, loc="upper center", ncol=5,
                bbox_to_anchor=(0.5, 0.98), fontsize=10, frameon=False,
                title=f"{pat} — {BRAIN_BAND_TEX_DICT.get(band, band)} — "
                      f"Cohesion-CBR (rsPost-anchored, size ∈ [{SIZE_MIN}, {SIZE_MAX}]; "
                      f"α = top1 − top2)",
                title_fontsize=11)

    out_dir = OUT_DIR / pat
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{band}_hierarchy-crossphase.pdf"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path, leaf_rows, anchor_rows_clean


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_leaves: list[dict] = []
    all_anchors: list[dict] = []
    skipped: list[str] = []
    n_written = 0

    for pat in PATIENTS_4PHASE:
        print(f"{pat} ...", flush=True)
        for band in BAND_ORDER:
            out, leaves, anchors = figure_for(pat, band)
            if out is None:
                skipped.append(f"{pat}/{band}")
                continue
            all_leaves.extend(leaves)
            all_anchors.extend(anchors)
            n_written += 1

    leaf_csv = OUT_DIR / "leaf_assignment.csv"
    with leaf_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "patient", "band", "leaf_id", "dominant", "confidence",
            "A_trace", "A_persist", "A_reset", "A_rearrange",
        ])
        w.writeheader(); w.writerows(all_leaves)

    anchor_csv = OUT_DIR / "anchor_classification.csv"
    with anchor_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "patient", "band", "anchor_row", "size", "h_rel",
            "COH_rpre", "COH_tl", "COH_tt", "COH_task",
            "J_rpre", "J_tl", "J_tt",
            "T_rpre", "T_tl", "T_tt",
            "a_trace", "a_persist", "a_reset", "a_rearrange",
        ])
        w.writeheader(); w.writerows(all_anchors)

    if skipped:
        (OUT_DIR / "missing.txt").write_text("\n".join(skipped) + "\n")

    print(f"\n{n_written} figures written, {len(skipped)} skipped")
    print(f"Wrote {leaf_csv.relative_to(ROOT)} ({len(all_leaves)} leaves)")
    print(f"Wrote {anchor_csv.relative_to(ROOT)} ({len(all_anchors)} anchors)")


if __name__ == "__main__":
    main()
