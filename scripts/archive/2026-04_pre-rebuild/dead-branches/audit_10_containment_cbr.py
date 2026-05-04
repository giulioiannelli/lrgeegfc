#!/usr/bin/env python3
"""Audit Step 10 — Containment-CBR (variant B).

Per scope: .agents/guides/task-persistence-investigation/2026-04-25_containment-cbr.md

Replaces symmetric Jaccard with directional containment
    C(S, T) = max_{v: |leaves(v)| ≤ κ} |S ∩ leaves(v)| / |S|.
Restores sensitivity to small anchors that are clean subsets of larger
reference subtrees. Both anchor variants (rsPost-anchored, task-anchored)
emitted side by side.

Outputs:
  data/audit/per_patient_hierarchy_containment/{rpost_anchor,task_anchor}/{Pat_XX}/{band}_hierarchy-crossphase.pdf
  data/audit/per_patient_hierarchy_containment/{rpost_anchor,task_anchor}/leaf_assignment.csv
  data/audit/per_patient_hierarchy_containment/{rpost_anchor,task_anchor}/anchor_classification.csv
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

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_08_per_patient_hierarchy import (  # noqa: E402
    BAND_ORDER, PHASE_ORDER, PHASE_TEX, BRAIN_BAND_TEX_DICT,
    PATIENTS_4PHASE, PATTERN_COLOR, PATTERN_PRIORITY,
    SIZE_MIN, SIZE_MAX,
    nodes_with_leaves, _load_Z,
    plot_dendrogram_with_strip,
)


OUT_DIR = ROOT / "data" / "audit" / "per_patient_hierarchy_containment"
TAU_MATCH = 0.75
TAU_DISSIM = 0.50


def containment(S: frozenset[int], ref_nodes: list[dict],
                 kappa: int) -> tuple[float, bool]:
    """Returns (best containment, cap_was_binding) for S against ref tree
    capped at size ≤ kappa. cap_binding = True iff the unconstrained max
    would have been higher than the capped max (the cap clipped the search)."""
    if not S:
        return 0.0, False
    capped_best = 0.0
    uncapped_best = 0.0
    sz_S = len(S)
    for rn in ref_nodes:
        score = len(S & rn["leaves"]) / sz_S
        if score > uncapped_best:
            uncapped_best = score
        if rn["size"] <= kappa and score > capped_best:
            capped_best = score
    return capped_best, uncapped_best > capped_best + 1e-12


def classify_C_rsPost(C_rpre: float, C_tl: float, C_tt: float) -> str | None:
    if C_rpre <= TAU_DISSIM and min(C_tl, C_tt) >= TAU_MATCH:
        return "trace"
    if min(C_rpre, C_tl, C_tt) >= TAU_MATCH:
        return "persist"
    if C_rpre >= TAU_MATCH and min(C_tl, C_tt) <= TAU_DISSIM:
        return "reset"
    if max(C_rpre, C_tl, C_tt) <= TAU_DISSIM:
        return "rearrange"
    return None


def classify_C_taskT(C_rpre: float, C_rpost: float) -> str | None:
    if C_rpre <= TAU_DISSIM and C_rpost >= TAU_MATCH:
        return "trace"
    if C_rpre >= TAU_MATCH and C_rpost >= TAU_MATCH:
        return "persist"
    if C_rpre >= TAU_MATCH and C_rpost <= TAU_DISSIM:
        return "reset"
    if C_rpre <= TAU_DISSIM and C_rpost <= TAU_DISSIM:
        return "rearrange"
    return None


def build_categories(Zs: dict[str, np.ndarray], variant: str
                       ) -> tuple[np.ndarray, list[dict]]:
    """variant ∈ {'rpost', 'task'}."""
    if variant == "rpost":
        anchor_Z = Zs["rest_post"]
    else:  # 'task'
        anchor_Z = Zs.get("task_test")
        if anchor_Z is None:
            return np.empty(0, dtype=object), []

    nodes_anchor = nodes_with_leaves(anchor_Z)
    nodes_rpre = nodes_with_leaves(Zs["rest_pre"])
    nodes_rpost = nodes_with_leaves(Zs["rest_post"])
    nodes_tl = nodes_with_leaves(Zs["task_learn"])
    nodes_tt = nodes_with_leaves(Zs["task_test"]) if Zs["task_test"] is not None else None

    n = anchor_Z.shape[0] + 1
    kappa = n // 2
    leaf_patterns: list[set[str]] = [set() for _ in range(n)]
    rows: list[dict] = []
    dmax_anchor = float(anchor_Z[-1, 2])

    for η in nodes_anchor:
        if η["size"] < SIZE_MIN or η["size"] > SIZE_MAX:
            continue
        S = η["leaves"]
        C_rpre, cb_rpre = containment(S, nodes_rpre, kappa)
        C_rpost, cb_rpost = containment(S, nodes_rpost, kappa)
        C_tl, cb_tl = containment(S, nodes_tl, kappa)
        C_tt = float("nan"); cb_tt = False
        if nodes_tt is not None:
            C_tt, cb_tt = containment(S, nodes_tt, kappa)

        if variant == "rpost":
            label = (classify_C_rsPost(C_rpre, C_tl, C_tt)
                       if not np.isnan(C_tt) else None)
        else:
            label = classify_C_taskT(C_rpre, C_rpost)

        rows.append({
            "anchor_row": int(η["row"]),
            "size": int(η["size"]),
            "h_rel": float(η["h"]) / dmax_anchor if dmax_anchor > 0 else 0.0,
            "C_rpre": C_rpre, "C_rpost": C_rpost,
            "C_tl": C_tl, "C_tt": C_tt,
            "cap_bound_rpre": int(cb_rpre),
            "cap_bound_rpost": int(cb_rpost),
            "cap_bound_tl": int(cb_tl),
            "cap_bound_tt": int(cb_tt),
            "label": label if label is not None else "",
            "variant": variant,
        })
        if label is None:
            continue
        for lid in S:
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


def figure(pat: str, band: str, variant: str
            ) -> tuple[Path | None, list[dict], list[dict]]:
    Zs = {ph: _load_Z(pat, ph, band) for ph in PHASE_ORDER}
    if Zs["rest_post"] is None or Zs["rest_pre"] is None or Zs["task_learn"] is None:
        return None, [], []
    if variant == "task" and Zs["task_test"] is None:
        return None, [], []

    cats, anchor_rows = build_categories(Zs, variant)
    if cats.size == 0:
        return None, [], []

    leaf_rows = [{"patient": pat, "band": band, "leaf_id": int(lid),
                   "category": cat, "variant": variant}
                  for lid, cat in enumerate(cats)]
    for r in anchor_rows:
        r["patient"] = pat; r["band"] = band

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
        else:
            plot_dendrogram_with_strip(ax_dend, ax_strip, Z, cats,
                                         title=PHASE_TEX[phase])
        if ip == 0:
            ax_dend.set_ylabel(BRAIN_BAND_TEX_DICT.get(band, band), fontsize=12)

    handles = [Patch(color=PATTERN_COLOR[p], label=p.upper())
                for p in PATTERN_PRIORITY]
    var_label = "rsPost-anchored" if variant == "rpost" else "task_test-anchored"
    n = Zs["rest_post"].shape[0] + 1
    fig.legend(handles=handles, loc="upper center", ncol=4,
                bbox_to_anchor=(0.5, 0.98), fontsize=10, frameon=False,
                title=f"{pat} — {BRAIN_BAND_TEX_DICT.get(band, band)} — "
                      f"Containment-CBR ({var_label}, κ = ⌊N/2⌋ = {n // 2})",
                title_fontsize=11)

    sub = "rpost_anchor" if variant == "rpost" else "task_anchor"
    out_dir = OUT_DIR / sub / pat
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{band}_hierarchy-crossphase.pdf"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path, leaf_rows, anchor_rows


def write_csvs(variant_label: str, leaves: list[dict], anchors: list[dict]) -> None:
    sub = OUT_DIR / variant_label
    sub.mkdir(parents=True, exist_ok=True)
    leaf_csv = sub / "leaf_assignment.csv"
    with leaf_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["patient", "band", "leaf_id",
                                            "category", "variant"])
        w.writeheader(); w.writerows(leaves)
    anchor_csv = sub / "anchor_classification.csv"
    with anchor_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "patient", "band", "anchor_row", "size", "h_rel",
            "C_rpre", "C_rpost", "C_tl", "C_tt",
            "cap_bound_rpre", "cap_bound_rpost", "cap_bound_tl", "cap_bound_tt",
            "label", "variant",
        ])
        w.writeheader(); w.writerows(anchors)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    by_variant_leaves: dict[str, list[dict]] = {"rpost": [], "task": []}
    by_variant_anchors: dict[str, list[dict]] = {"rpost": [], "task": []}
    skipped: list[str] = []
    n_written = 0

    for pat in PATIENTS_4PHASE:
        print(f"{pat} ...", flush=True)
        for band in BAND_ORDER:
            for variant in ("rpost", "task"):
                out, leaves, anchors = figure(pat, band, variant)
                if out is None:
                    skipped.append(f"{variant}/{pat}/{band}")
                    continue
                by_variant_leaves[variant].extend(leaves)
                by_variant_anchors[variant].extend(anchors)
                n_written += 1
        print(f"  done {pat}")

    write_csvs("rpost_anchor", by_variant_leaves["rpost"], by_variant_anchors["rpost"])
    write_csvs("task_anchor", by_variant_leaves["task"], by_variant_anchors["task"])

    if skipped:
        (OUT_DIR / "missing.txt").write_text("\n".join(skipped) + "\n")

    print(f"\n{n_written} figures written across both anchor variants")
    print(f"{len(skipped)} (variant, patient, band) skipped")


if __name__ == "__main__":
    main()
