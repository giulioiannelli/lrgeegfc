#!/usr/bin/env python3
"""Audit Step 9 — Task-anchored CBR (variant A).

Per scope: .agents/guides/task-persistence-investigation/2026-04-25_task-anchored-cbr.md

Same Jaccard machinery and thresholds as the legacy CBR, but anchors on
`task_test` internal nodes instead of `rest_post` ones. Surfaces task
subtrees that fragment in rsPost (failure mode #1 of the legacy procedure).

Outputs:
  data/audit/per_patient_hierarchy_taskanchored/{Pat_XX}/{band}_hierarchy-crossphase.pdf
  data/audit/per_patient_hierarchy_taskanchored/leaf_assignment.csv
  data/audit/per_patient_hierarchy_taskanchored/anchor_classification.csv
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

# Reuse audit_08's plotting / dendrogram helpers verbatim (same figure layout).
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_08_per_patient_hierarchy import (  # noqa: E402
    BAND_ORDER, PHASE_ORDER, PHASE_TEX, BRAIN_BAND_TEX_DICT,
    PATIENTS_4PHASE, PATTERN_COLOR, PATTERN_PRIORITY, BACKGROUND_GRAY,
    SIZE_MIN, SIZE_MAX,
    nodes_with_leaves, best_jaccard, _load_Z,
    plot_dendrogram_with_strip,
)


OUT_DIR = ROOT / "data" / "audit" / "per_patient_hierarchy_taskanchored"
TAU_MATCH = 0.75
TAU_DISSIM = 0.50


def classify_TACBR(J_rpre: float, J_rpost: float) -> str | None:
    """Per scope §2.2 — mutually exclusive priority chain."""
    if J_rpre <= TAU_DISSIM and J_rpost >= TAU_MATCH:
        return "trace"
    if J_rpre >= TAU_MATCH and J_rpost >= TAU_MATCH:
        return "persist"
    if J_rpre >= TAU_MATCH and J_rpost <= TAU_DISSIM:
        return "reset"
    if J_rpre <= TAU_DISSIM and J_rpost <= TAU_DISSIM:
        return "rearrange"
    return None


def classify_anchors_and_project(Zs: dict[str, np.ndarray]
                                   ) -> tuple[np.ndarray, list[dict]]:
    """Returns (leaf_cats, per_anchor_rows).

    Anchor universe = internal nodes of T_taskT with size ∈ [SIZE_MIN, SIZE_MAX].
    """
    if Zs.get("task_test") is None:
        return np.empty(0, dtype=object), []

    nodes_tt = nodes_with_leaves(Zs["task_test"])
    nodes_rpre = nodes_with_leaves(Zs["rest_pre"])
    nodes_rpost = nodes_with_leaves(Zs["rest_post"])
    nodes_tl = nodes_with_leaves(Zs["task_learn"])

    n = Zs["task_test"].shape[0] + 1
    leaf_patterns: list[set[str]] = [set() for _ in range(n)]
    rows: list[dict] = []
    dmax_tt = float(Zs["task_test"][-1, 2])

    for η in nodes_tt:
        if η["size"] < SIZE_MIN or η["size"] > SIZE_MAX:
            continue
        J_rpre = best_jaccard(η["leaves"], nodes_rpre)
        J_rpost = best_jaccard(η["leaves"], nodes_rpost)
        J_tl = best_jaccard(η["leaves"], nodes_tl)
        label = classify_TACBR(J_rpre, J_rpost)
        rows.append({
            "anchor_row": int(η["row"]),
            "size": int(η["size"]),
            "h_rel": float(η["h"]) / dmax_tt if dmax_tt > 0 else 0.0,
            "J_rpre": J_rpre, "J_rpost": J_rpost, "J_tl": J_tl,
            "label": label if label is not None else "",
        })
        if label is None:
            continue
        for lid in η["leaves"]:
            leaf_patterns[lid].add(label)

    out = np.empty(n, dtype=object)
    for lid in range(n):
        chosen = "rearrange"
        for p in PATTERN_PRIORITY:
            if p in leaf_patterns[lid]:
                chosen = p
                break
        out[lid] = chosen
    return out, rows


def figure_for_patient_band(pat: str, band: str
                              ) -> tuple[Path | None, list[dict], list[dict]]:
    Zs = {ph: _load_Z(pat, ph, band) for ph in PHASE_ORDER}
    if Zs["rest_post"] is None or Zs["rest_pre"] is None or Zs["task_learn"] is None or Zs["task_test"] is None:
        return None, [], []

    cats, anchor_rows = classify_anchors_and_project(Zs)

    leaf_rows = [{"patient": pat, "band": band, "leaf_id": int(lid),
                   "category": cat} for lid, cat in enumerate(cats)]
    for r in anchor_rows:
        r["patient"] = pat
        r["band"] = band

    n_phases = len(PHASE_ORDER)
    fig = plt.figure(figsize=(4.5 * n_phases, 5.4))
    outer = fig.add_gridspec(1, n_phases, wspace=0.18, top=0.84, bottom=0.05)

    for ip, phase in enumerate(PHASE_ORDER):
        cell = outer[0, ip].subgridspec(2, 1, height_ratios=[8, 1], hspace=0.04)
        ax_dend = fig.add_subplot(cell[0])
        ax_strip = fig.add_subplot(cell[1])
        Z = Zs[phase]
        plot_dendrogram_with_strip(ax_dend, ax_strip, Z, cats,
                                     title=PHASE_TEX[phase])
        if ip == 0:
            ax_dend.set_ylabel(BRAIN_BAND_TEX_DICT.get(band, band), fontsize=12)

    handles = [Patch(color=PATTERN_COLOR[p], label=p.upper())
                for p in PATTERN_PRIORITY]
    fig.legend(handles=handles, loc="upper center", ncol=4,
                bbox_to_anchor=(0.5, 0.98), fontsize=10, frameon=False,
                title=f"{pat} — {BRAIN_BAND_TEX_DICT.get(band, band)} — "
                      f"TA-CBR (task_test-anchored, size ∈ [{SIZE_MIN}, {SIZE_MAX}])",
                title_fontsize=11)

    out_dir = OUT_DIR / pat
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{band}_hierarchy-crossphase.pdf"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path, leaf_rows, anchor_rows


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_leaves: list[dict] = []
    all_anchors: list[dict] = []
    skipped: list[str] = []
    n_written = 0

    for pat in PATIENTS_4PHASE:
        print(f"{pat} ...", flush=True)
        for band in BAND_ORDER:
            out, leaves, anchors = figure_for_patient_band(pat, band)
            if out is None:
                skipped.append(f"{pat}/{band}")
                print(f"  {band}: missing phase (likely Pat_14 task_test), skipped")
                continue
            all_leaves.extend(leaves)
            all_anchors.extend(anchors)
            n_written += 1
            print(f"  wrote {out.relative_to(ROOT)}")

    leaf_csv = OUT_DIR / "leaf_assignment.csv"
    with leaf_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["patient", "band", "leaf_id", "category"])
        w.writeheader(); w.writerows(all_leaves)

    anchor_csv = OUT_DIR / "anchor_classification.csv"
    with anchor_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "patient", "band", "anchor_row", "size", "h_rel",
            "J_rpre", "J_rpost", "J_tl", "label",
        ])
        w.writeheader(); w.writerows(all_anchors)

    if skipped:
        (OUT_DIR / "missing.txt").write_text("\n".join(skipped) + "\n")

    print(f"\nWrote {leaf_csv.relative_to(ROOT)}  ({len(all_leaves)} rows)")
    print(f"Wrote {anchor_csv.relative_to(ROOT)}  ({len(all_anchors)} anchors)")
    print(f"{n_written} figures, {len(skipped)} (patient, band) skipped")


if __name__ == "__main__":
    main()
