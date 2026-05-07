#!/usr/bin/env python3
"""Audit 51b — Section 5.6 Figure 2: per-patient beta 4x3 dendrogram composite.

Rows = taxonomy class (trace / reset / rearrange / anchor).
Cols = phase (rest_pre / task_test / rest_post).
Each cell is a Pat_NN beta dendrogram for that phase, with the row
class's modules highlighted as colored leaf bars + colored internal
links. Within-class shades distinguish individual modules.

Module sets are recomputed at lambda=0 (size>=3 stratum for rearrange)
using the canonical detection functions in audits 47/48/49/50.
Highlighting uses the strict intersection of per-class leaf sets:
  trace     = leaves_tt & leaves_post
  reset     = leaves_pre & leaves_post
  rearrange = leaves_post − strict_leaves(trace ∪ reset ∪ anchor)
  anchor    = leaves_pre & leaves_tt & leaves_post

Outputs:
  data/audit/section5_v2_round3_redo/figures/taxonomy_pat<NN>_beta_composite.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import numpy as np
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_47_kc_trace_network_view import find_trace_pairs
from audit_48_kc_reset_module_view import find_reset_pairs
from audit_49_kc_rearrangement_module_view import find_rearrangement_modules
from audit_50_kc_anchor_module_view import find_anchor_triples

PATIENTS = (
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
)
BAND = "beta"
PHASES = ("rest_pre", "task_test", "rest_post")
PHASE_LABELS = {p: p for p in PHASES}
CLASSES = ["trace", "reset", "rearrange", "anchor", "diffuse"]
CLASS_COLOR = {
    "trace":     "#2c7fb8",
    "reset":     "#41ab5d",
    "rearrange": "#fc9272",
    "anchor":    "#807dba",
    "diffuse":   "#8c8c8c",
}
CLASS_CMAPS = {
    "trace":     "Blues",
    "reset":     "Greens",
    "rearrange": "Oranges",
    "anchor":    "Purples",
    "diffuse":   "Greys",
}
GRAY_LEAF = "#c8c8c8"
GRAY_LINK = "#bababa"

OUT_DIR = ROOT / "data/audit/section5_v2_round3_redo/figures"


def shades_for_class(cls: str, n: int) -> list:
    if n == 0:
        return []
    cmap = plt.get_cmap(CLASS_CMAPS[cls])
    if n == 1:
        return [to_hex(cmap(0.75))]
    return [to_hex(cmap(0.95 - 0.45 * i / (n - 1))) for i in range(n)]


def load_phase(patient: str, phase: str):
    A = load_fc_matrix(patient, phase, BAND, fc_method="imcoh_abs")
    res = load_lrg_result(patient, phase, BAND, fc_method="imcoh_abs")
    if A is None or res is None or res.linkage_matrix is None:
        raise FileNotFoundError(f"{patient} {phase} {BAND}")
    return np.asarray(A), np.asarray(res.linkage_matrix)


def detect_modules(Z_pre, Z_tt, Z_post, n):
    out = {}
    out["trace"]     = find_trace_pairs(Z_tt, Z_post, Z_pre, n)[0.0]
    out["reset"]     = find_reset_pairs(Z_pre, Z_post, Z_tt, n)[0.0]
    out["rearrange"] = [m for m in find_rearrangement_modules(Z_pre, Z_tt, Z_post, n)
                        if m["size"] >= 3]
    out["anchor"]    = find_anchor_triples(Z_pre, Z_tt, Z_post, n)[0.0]
    return out


def strict_module_leaves(cls: str, module: dict) -> set:
    if cls == "trace":
        return module["leaves_tt"] & module["leaves_post"]
    if cls == "reset":
        return module["leaves_pre"] & module["leaves_post"]
    if cls == "anchor":
        return module["leaves_pre"] & module["leaves_tt"] & module["leaves_post"]
    if cls == "diffuse":
        return module["leaves"]
    return module.get("leaves_post_residual", module["leaves_post"])


def class_leaves_for_phase(cls: str, module: dict, phase: str) -> set:
    return strict_module_leaves(cls, module)


def assign_leaf_for_cell(cls: str, modules: list, phase: str, n: int) -> list:
    out = [None] * n
    for m in sorted(modules, key=lambda x: x.get("size", len(class_leaves_for_phase(cls, x, phase)))):
        leaves = class_leaves_for_phase(cls, m, phase)
        for lf in leaves:
            if out[lf] is None:
                out[lf] = m
    return out


def make_link_color_func(Z, modules, cls, phase, n):
    leaf_m = assign_leaf_for_cell(cls, modules, phase, n)
    sorted_mods = sorted(modules, key=lambda x: x.get("size", 0))

    subtree_leaves = {i: {i} for i in range(n)}
    link_color = {}
    for i, row in enumerate(Z):
        a, b = int(row[0]), int(row[1])
        leaves = subtree_leaves[a] | subtree_leaves[b]
        new_id = n + i
        subtree_leaves[new_id] = leaves
        link_color[new_id] = GRAY_LINK
        for m in sorted_mods:
            mod_leaves = class_leaves_for_phase(cls, m, phase)
            if mod_leaves and leaves <= mod_leaves:
                link_color[new_id] = m["_shade"]
                break

    leaf_color = {}
    for lf in range(n):
        leaf_color[lf] = leaf_m[lf]["_shade"] if leaf_m[lf] is not None else GRAY_LEAF

    def fn(node_id):
        if node_id < n:
            return leaf_color[node_id]
        return link_color.get(node_id, GRAY_LINK)

    return fn, leaf_m


def draw_cell(ax, Z, modules, cls, phase, n):
    fn, leaf_m = make_link_color_func(Z, modules, cls, phase, n)
    R = dendrogram(
        Z, ax=ax, no_labels=True,
        color_threshold=0,
        above_threshold_color=GRAY_LINK,
        link_color_func=fn,
    )
    merge_heights = sorted(set(Z[:, 2]))
    tmin = 0.0
    tmax = merge_heights[-1] * 1.05
    ax.set_ylim(tmin, tmax)
    ax.set_xticks([])
    ax.tick_params(axis="y", labelsize=6)
    for sn in ("top", "right", "bottom"):
        ax.spines[sn].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)

    leaf_parent_h = {}
    for i, row in enumerate(Z):
        a, b, h = int(row[0]), int(row[1]), float(row[2])
        if a < n:
            leaf_parent_h.setdefault(a, h)
        if b < n:
            leaf_parent_h.setdefault(b, h)

    leaf_order = R["leaves"]
    for i, lf_idx in enumerate(leaf_order):
        m = leaf_m[lf_idx]
        if m is None:
            continue
        x = 5.0 + 10.0 * i
        h_parent = leaf_parent_h.get(lf_idx, tmax * 0.05)
        ax.vlines(x, tmin, h_parent, color=m["_shade"], lw=1.4, zorder=4)


def render_for_patient(patient: str) -> dict:
    A_by_phase = {}
    Z_by_phase = {}
    for ph in PHASES:
        A, Z = load_phase(patient, ph)
        A_by_phase[ph] = A
        Z_by_phase[ph] = Z
    n = A_by_phase["rest_pre"].shape[0]
    for ph in PHASES:
        if A_by_phase[ph].shape[0] != n:
            n = min(A_by_phase[ph].shape[0], n)

    modules_by_class = detect_modules(
        Z_by_phase["rest_pre"], Z_by_phase["task_test"],
        Z_by_phase["rest_post"], n,
    )
    raw_counts = {cls: len(modules_by_class[cls])
                  for cls in ("trace", "reset", "rearrange", "anchor")}

    non_rearrange_leaves: set = set()
    for cls in ("trace", "reset", "anchor"):
        for m in modules_by_class[cls]:
            non_rearrange_leaves |= strict_module_leaves(cls, m)
    residual_modules = []
    for m in modules_by_class["rearrange"]:
        residual = m["leaves_post"] - non_rearrange_leaves
        if len(residual) >= 3:
            m = {**m, "leaves_post_residual": residual, "size": len(residual)}
            residual_modules.append(m)
    modules_by_class["rearrange"] = residual_modules

    classified_leaves: set = set()
    for cls in ("trace", "reset", "rearrange", "anchor"):
        for m in modules_by_class[cls]:
            classified_leaves |= strict_module_leaves(cls, m)
    diffuse_leaves = set(range(n)) - classified_leaves
    modules_by_class["diffuse"] = (
        [{"leaves": diffuse_leaves, "size": len(diffuse_leaves)}]
        if diffuse_leaves else []
    )

    rendered_counts = {cls: len(modules_by_class[cls])
                       for cls in ("trace", "reset", "rearrange", "anchor")}
    rendered_counts["diffuse_leaves"] = len(diffuse_leaves)

    for cls in CLASSES:
        mods = modules_by_class[cls]
        mods_by_size = sorted(mods, key=lambda m: -len(strict_module_leaves(cls, m)))
        shades = shades_for_class(cls, len(mods_by_size))
        for m, sh in zip(mods_by_size, shades):
            m["_shade"] = sh

    fig = plt.figure(figsize=(7.0, 6.5))
    gs = fig.add_gridspec(
        5, 3,
        left=0.085, right=0.99, top=0.94, bottom=0.085,
        hspace=0.18, wspace=0.04,
    )
    row_axes: dict = {}
    for r, cls in enumerate(CLASSES):
        modules = modules_by_class[cls]
        for c, phase in enumerate(PHASES):
            sharey = row_axes.get(r)
            ax = fig.add_subplot(gs[r, c], sharey=sharey)
            if sharey is None:
                row_axes[r] = ax
            draw_cell(ax, Z_by_phase[phase], modules, cls, phase, n)
            if c > 0:
                ax.tick_params(axis="y", labelleft=False)
            if r == 0:
                ax.set_title(PHASE_LABELS[phase], fontsize=10, pad=4)
            if c == 0:
                ax.text(
                    -0.20, 0.5, cls, transform=ax.transAxes,
                    rotation=90, ha="center", va="center",
                    fontsize=11, color=CLASS_COLOR[cls],
                    fontweight="bold",
                )

    line1 = (
        f"{patient}  $\\beta$:   trace {raw_counts['trace']}  "
        f"·  reset {raw_counts['reset']}  "
        f"·  rearrange {rendered_counts['rearrange']}/{raw_counts['rearrange']} (residual)  "
        f"·  anchor {raw_counts['anchor']}  "
        f"·  diffuse {rendered_counts['diffuse_leaves']}/{n} leaves"
    )
    line2 = "highlight = strict intersection of per-class leaf sets"
    fig.text(0.5, 0.040, line1, ha="center", va="center", fontsize=8)
    fig.text(0.5, 0.018, line2, ha="center", va="center", fontsize=7,
             color="#666")

    pat_lower = patient.replace("Pat_", "pat").lower()
    out_path = OUT_DIR / f"taxonomy_{pat_lower}_beta_composite.pdf"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    return {
        "patient": patient,
        "raw_counts": raw_counts,
        "rendered_counts": rendered_counts,
        "out_path": out_path,
    }


def main() -> int:
    rows = []
    for patient in PATIENTS:
        info = render_for_patient(patient)
        rc, rd = info["raw_counts"], info["rendered_counts"]
        print(f"{patient}: trace={rc['trace']:2d} reset={rc['reset']:2d} "
              f"rearrange={rd['rearrange']:2d}/{rc['rearrange']:2d} (residual) "
              f"anchor={rc['anchor']:2d} diffuse={rd['diffuse_leaves']:3d} leaves "
              f"-> {info['out_path'].name}")
        rows.append({
            "patient": patient, "band": BAND,
            "trace": rc["trace"], "reset": rc["reset"],
            "rearrange_raw": rc["rearrange"],
            "rearrange_residual": rd["rearrange"],
            "anchor": rc["anchor"],
            "diffuse_leaves": rd["diffuse_leaves"],
        })

    import pandas as pd
    out_csv = ROOT / "data/audit/section5_v2_round3_redo/tables/taxonomy_per_patient_beta_counts.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
