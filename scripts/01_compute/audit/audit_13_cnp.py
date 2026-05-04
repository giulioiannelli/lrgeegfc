#!/usr/bin/env python3
"""Audit Step 13 — Cophenetic-Neighbourhood Persistence (CNP).

Per scope: .agents/guides/task-persistence-investigation/2026-04-25_cophenetic-neighbourhood.md

Per-leaf hierarchy-level classifier. For each leaf ℓ and phase φ, build
the cophenetic-distance vector D^φ_ℓ. Compute Spearman correlations
between phase-vector pairs (excluding self), assign 4 corner-distance
soft affinities on (ρ_pre_post, ρ_task_post), color leaves by dominant.

Outputs:
  data/audit/per_patient_hierarchy_cnp/{Pat_XX}/{band}_hierarchy-crossphase.pdf
  data/audit/per_patient_hierarchy_cnp/leaf_assignment.csv
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import cophenet, dendrogram
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_08_per_patient_hierarchy import (  # noqa: E402
    BAND_ORDER, PHASE_ORDER, PHASE_TEX, BRAIN_BAND_TEX_DICT,
    PATIENTS_4PHASE, PATTERN_COLOR, BACKGROUND_GRAY, STRIP_BG,
    _load_Z,
)


OUT_DIR = ROOT / "data" / "audit" / "per_patient_hierarchy_cnp"
PATTERNS = ["trace", "persist", "reset", "rearrange"]
CORNERS = {
    "trace":     (-1.0,  1.0),
    "persist":   ( 1.0,  1.0),
    "reset":     ( 1.0, -1.0),
    "rearrange": (-1.0, -1.0),
}
_DIAG = float(2 * np.sqrt(2))    # corner-to-opposite-corner distance on [-1,1]²


def cophenet_matrix(Z: np.ndarray) -> np.ndarray:
    """Full N × N cophenetic distance matrix."""
    return squareform(cophenet(Z))


def corner_affinity(rho_pp: float, rho_tp: float, corner: tuple[float, float]
                     ) -> float:
    cx, cy = corner
    d = np.hypot(rho_pp - cx, rho_tp - cy)
    return float(max(0.0, 1.0 - d / _DIAG))


def per_leaf_classification(Zs: dict[str, np.ndarray]) -> list[dict]:
    """Returns one dict per leaf with the CNP fields."""
    if Zs["rest_post"] is None or Zs["rest_pre"] is None or Zs["task_learn"] is None:
        return []

    D_rpre = cophenet_matrix(Zs["rest_pre"])
    D_tl = cophenet_matrix(Zs["task_learn"])
    D_rpost = cophenet_matrix(Zs["rest_post"])
    has_tt = Zs["task_test"] is not None
    if has_tt:
        D_tt = cophenet_matrix(Zs["task_test"])
        D_task = (D_tl + D_tt) / 2.0
    else:
        D_task = D_tl

    n = D_rpost.shape[0]
    out: list[dict] = []
    for ℓ in range(n):
        mask = np.ones(n, dtype=bool); mask[ℓ] = False
        v_rpre = D_rpre[ℓ, mask]
        v_task = D_task[ℓ, mask]
        v_rpost = D_rpost[ℓ, mask]

        rho_pp, _ = spearmanr(v_rpre, v_rpost)
        rho_tp, _ = spearmanr(v_task, v_rpost)
        rho_pt, _ = spearmanr(v_rpre, v_task)
        # Spearman returns NaN if a vector has zero variance (constant); guard.
        rho_pp = 0.0 if np.isnan(rho_pp) else float(rho_pp)
        rho_tp = 0.0 if np.isnan(rho_tp) else float(rho_tp)
        rho_pt = 0.0 if np.isnan(rho_pt) else float(rho_pt)

        affs = {p: corner_affinity(rho_pp, rho_tp, CORNERS[p]) for p in PATTERNS}
        ordered = sorted(affs.items(), key=lambda kv: -kv[1])
        dominant = ordered[0][0]
        confidence = float(ordered[0][1] - ordered[1][1])
        out.append({
            "leaf_id": int(ℓ),
            "rho_pre_post": rho_pp,
            "rho_task_post": rho_tp,
            "rho_pre_task": rho_pt,
            "a_trace": affs["trace"], "a_persist": affs["persist"],
            "a_reset": affs["reset"], "a_rearrange": affs["rearrange"],
            "dominant": dominant, "confidence": confidence,
        })
    return out


# ── plotting ──────────────────────────────────────────────────────────────

def link_colors_from_dominant(Z: np.ndarray, dominant: list[str]) -> dict[int, str]:
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


def plot_panel(ax_dend, ax_strip, Z: np.ndarray,
                 dominant: list[str], title: str = "") -> None:
    if Z is None:
        ax_dend.text(0.5, 0.5, "missing", ha="center", va="center",
                      transform=ax_dend.transAxes, color="#999")
        ax_dend.set_xticks([]); ax_dend.set_yticks([])
        ax_strip.axis("off")
        if title:
            ax_dend.set_title(title, fontsize=11)
        return
    lc = link_colors_from_dominant(Z, dominant)
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


def figure_for(pat: str, band: str
                ) -> tuple[Path | None, list[dict]]:
    Zs = {ph: _load_Z(pat, ph, band) for ph in PHASE_ORDER}
    leaves = per_leaf_classification(Zs)
    if not leaves:
        return None, []
    n = len(leaves)
    dominant = [r["dominant"] for r in leaves]

    rows_csv = [{"patient": pat, "band": band, **r} for r in leaves]

    n_phases = len(PHASE_ORDER)
    fig = plt.figure(figsize=(4.5 * n_phases, 5.4))
    outer = fig.add_gridspec(1, n_phases, wspace=0.18, top=0.84, bottom=0.05)
    for ip, phase in enumerate(PHASE_ORDER):
        cell = outer[0, ip].subgridspec(2, 1, height_ratios=[8, 1], hspace=0.04)
        ax_dend = fig.add_subplot(cell[0])
        ax_strip = fig.add_subplot(cell[1])
        plot_panel(ax_dend, ax_strip, Zs[phase], dominant, title=PHASE_TEX[phase])
        if ip == 0:
            ax_dend.set_ylabel(BRAIN_BAND_TEX_DICT.get(band, band), fontsize=12)

    handles = [Patch(color=PATTERN_COLOR[p], label=p.upper()) for p in PATTERNS]
    fig.legend(handles=handles, loc="upper center", ncol=4,
                bbox_to_anchor=(0.5, 0.98), fontsize=10, frameon=False,
                title=f"{pat} — {BRAIN_BAND_TEX_DICT.get(band, band)} — "
                      f"CNP (per-leaf cophenetic-neighbourhood Spearman; "
                      f"corner-distance affinities)",
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
    skipped: list[str] = []
    n_written = 0
    for pat in PATIENTS_4PHASE:
        print(f"{pat} ...", flush=True)
        for band in BAND_ORDER:
            out, rows = figure_for(pat, band)
            if out is None:
                skipped.append(f"{pat}/{band}")
                continue
            all_rows.extend(rows)
            n_written += 1

    csv_path = OUT_DIR / "leaf_assignment.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "patient", "band", "leaf_id",
            "dominant", "confidence",
            "rho_pre_post", "rho_task_post", "rho_pre_task",
            "a_trace", "a_persist", "a_reset", "a_rearrange",
        ])
        w.writeheader(); w.writerows(all_rows)

    if skipped:
        (OUT_DIR / "missing.txt").write_text("\n".join(skipped) + "\n")
    print(f"\n{n_written} figures written, {len(skipped)} skipped")
    print(f"Wrote {csv_path.relative_to(ROOT)} ({len(all_rows)} leaves)")


if __name__ == "__main__":
    main()
