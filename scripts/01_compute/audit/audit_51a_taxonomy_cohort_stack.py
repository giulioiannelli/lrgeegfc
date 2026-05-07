#!/usr/bin/env python3
"""Audit 51a — Section 5.6 Figure 1: cohort-total taxonomy bar chart.

Single-panel grouped bar chart. 6 bands x 4 classes (trace / reset /
rearrange / anchor) per band. Module counts read from the canonical
CSVs of audits 47/48/49/50 for trace/reset/anchor; rearrange uses
RESIDUAL counts (leaves_post minus strict-intersection leaves of
trace/reset/anchor modules at the same cell). Residuals are computed
inline by re-running detection at lambda=0 for the full 10 x 6 = 60
cohort cells (~5 min).

Inputs (canonical):
  data/audit/kc_trace_network_view/trace_subtrees_summary.csv
  data/audit/kc_reset_module_view/reset_subtrees_summary.csv
  data/audit/kc_rearrangement_module_view/rearrangement_subtrees_summary.csv
  data/audit/kc_anchor_module_view/anchor_subtrees_summary.csv

Outputs:
  data/audit/section5_v2_round3_redo/figures/taxonomy_cohort_stack.pdf
  data/audit/section5_v2_round3_redo/tables/taxonomy_cohort_stack_data.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.workflow.lrg import load_lrg_result

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_47_kc_trace_network_view import find_trace_pairs
from audit_48_kc_reset_module_view import find_reset_pairs
from audit_49_kc_rearrangement_module_view import find_rearrangement_modules
from audit_50_kc_anchor_module_view import find_anchor_triples

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
CLASSES = ["trace", "reset", "rearrange", "anchor", "diffuse"]
CLASS_COLORS = {
    "trace":     "#2c7fb8",
    "reset":     "#41ab5d",
    "rearrange": "#fc9272",
    "anchor":    "#807dba",
    "diffuse":   "#8c8c8c",
}
CLASS_LABEL = {
    "trace":     "trace (modules)",
    "reset":     "reset (modules)",
    "rearrange": "rearrange residual (modules)",
    "anchor":    "anchor (modules)",
    "diffuse":   "diffuse (leaves)",
}
COHORT_N9 = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_10",
             "Pat_13", "Pat_14", "Pat_15"]
COHORT_N10 = COHORT_N9 + ["Pat_08"]
COHORT_FOR = {
    "trace":     COHORT_N9,
    "reset":     COHORT_N10,
    "rearrange": COHORT_N10,
    "anchor":    COHORT_N10,
    "diffuse":   COHORT_N10,
}
COHORT_N_FOR = {"trace": 9, "reset": 10, "rearrange": 10, "anchor": 10, "diffuse": 10}
CSV_BY_CLASS = {
    "trace":     ROOT / "data/audit/kc_trace_network_view/trace_subtrees_summary.csv",
    "reset":     ROOT / "data/audit/kc_reset_module_view/reset_subtrees_summary.csv",
    "rearrange": ROOT / "data/audit/kc_rearrangement_module_view/rearrangement_subtrees_summary.csv",
    "anchor":    ROOT / "data/audit/kc_anchor_module_view/anchor_subtrees_summary.csv",
}
COUNT_COL = {
    "trace":     "n_modules_lambda0",
    "reset":     "n_modules_lambda0",
    "rearrange": "n_modules",
    "anchor":    "n_modules_lambda0",
}

EXPECTED_TOTAL_TRA_RES_ANC = {
    "trace":  {"delta": 15, "theta": 15, "alpha": 23, "beta": 23, "low_gamma": 21, "high_gamma": 10},
    "reset":  {"delta": 10, "theta": 15, "alpha": 13, "beta": 18, "low_gamma": 26, "high_gamma":  1},
    "anchor": {"delta": 11, "theta": 12, "alpha": 16, "beta": 43, "low_gamma": 45, "high_gamma": 17},
}

OUT_FIG = ROOT / "data/audit/section5_v2_round3_redo/figures/taxonomy_cohort_stack.pdf"
OUT_CSV = ROOT / "data/audit/section5_v2_round3_redo/tables/taxonomy_cohort_stack_data.csv"
RESIDUAL_CACHE = ROOT / "data/audit/section5_v2_round3_redo/tables/rearrange_residual_per_cell.csv"


def strict_leaves_for_class(cls: str, m: dict) -> set:
    if cls == "trace":
        return m["leaves_tt"] & m["leaves_post"]
    if cls == "reset":
        return m["leaves_pre"] & m["leaves_post"]
    if cls == "anchor":
        return m["leaves_pre"] & m["leaves_tt"] & m["leaves_post"]
    return m["leaves_post"]


def load_Z(patient: str, phase: str, band: str) -> np.ndarray | None:
    res = load_lrg_result(patient, phase, band, fc_method="imcoh_abs")
    if res is None or res.linkage_matrix is None:
        return None
    return np.asarray(res.linkage_matrix)


def compute_rearrange_residual_for_cell(patient: str, band: str) -> tuple[int, int, int, int]:
    Zs = {ph: load_Z(patient, ph, band) for ph in ("rest_pre", "task_test", "rest_post")}
    if any(Z is None for Z in Zs.values()):
        return 0, 0, 0, 0
    n = int(Zs["rest_pre"].shape[0]) + 1
    n = min(n, int(Zs["task_test"].shape[0]) + 1, int(Zs["rest_post"].shape[0]) + 1)
    Z_pre, Z_tt, Z_post = Zs["rest_pre"], Zs["task_test"], Zs["rest_post"]

    trace_mods  = find_trace_pairs(Z_tt, Z_post, Z_pre, n)[0.0]
    reset_mods  = find_reset_pairs(Z_pre, Z_post, Z_tt, n)[0.0]
    anchor_mods = find_anchor_triples(Z_pre, Z_tt, Z_post, n)[0.0]
    rearr_mods  = [m for m in find_rearrangement_modules(Z_pre, Z_tt, Z_post, n)
                   if m["size"] >= 3]

    excluded: set = set()
    for m in trace_mods:  excluded |= strict_leaves_for_class("trace", m)
    for m in reset_mods:  excluded |= strict_leaves_for_class("reset", m)
    for m in anchor_mods: excluded |= strict_leaves_for_class("anchor", m)

    classified: set = set(excluded)
    residual_count = 0
    for m in rearr_mods:
        residual_leaves = m["leaves_post"] - excluded
        if len(residual_leaves) >= 3:
            residual_count += 1
            classified |= residual_leaves
    diffuse_leaves = n - len(classified)
    return len(rearr_mods), residual_count, n, diffuse_leaves


def compute_residuals_all_cells(force: bool = False) -> pd.DataFrame:
    if not force and RESIDUAL_CACHE.exists():
        df = pd.read_csv(RESIDUAL_CACHE)
        print(f"loaded residuals cache: {RESIDUAL_CACHE.relative_to(ROOT)} "
              f"({len(df)} rows)")
        return df

    rows = []
    for patient in COHORT_N10:
        for band in BANDS:
            try:
                raw, residual, n, diffuse = compute_rearrange_residual_for_cell(patient, band)
            except Exception as e:
                print(f"[skip] {patient} {band}: {e}")
                continue
            print(f"{patient:8s} {band:11s} n={n:3d}  rearrange raw={raw:2d}  "
                  f"residual={residual:2d}  diffuse_leaves={diffuse:3d}")
            rows.append({"patient": patient, "band": band,
                         "raw": raw, "residual": residual,
                         "n_nodes": n, "diffuse_leaves": diffuse})
    df = pd.DataFrame(rows)
    RESIDUAL_CACHE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESIDUAL_CACHE, index=False)
    print(f"wrote residuals cache: {RESIDUAL_CACHE.relative_to(ROOT)}")
    return df


def collect(residuals: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cls in CLASSES:
        if cls == "rearrange":
            for band in BANDS:
                sub = residuals[residuals["band"] == band]
                total = int(sub["residual"].sum())
                cells = int((sub["residual"] > 0).sum())
                rows.append({
                    "band": band, "class": cls,
                    "n_modules": total, "n_cells_with_module": cells,
                    "cohort_n": COHORT_N_FOR[cls],
                })
            continue
        if cls == "diffuse":
            for band in BANDS:
                sub = residuals[residuals["band"] == band]
                total = int(sub["diffuse_leaves"].sum())
                cells = int((sub["diffuse_leaves"] > 0).sum())
                rows.append({
                    "band": band, "class": cls,
                    "n_modules": total, "n_cells_with_module": cells,
                    "cohort_n": COHORT_N_FOR[cls],
                })
            continue
        df = pd.read_csv(CSV_BY_CLASS[cls])
        df = df[df["patient"].isin(COHORT_FOR[cls])].copy()
        col = COUNT_COL[cls]
        for band in BANDS:
            sub = df[df["band"] == band]
            total = int(sub[col].fillna(0).sum())
            cells = int((sub[col].fillna(0) > 0).sum())
            rows.append({
                "band": band, "class": cls,
                "n_modules": total, "n_cells_with_module": cells,
                "cohort_n": COHORT_N_FOR[cls],
            })
    return pd.DataFrame(rows)


def sanity_check(df: pd.DataFrame) -> bool:
    rows = []
    pass_all = True
    for cls in ("trace", "reset", "anchor"):
        for band in BANDS:
            r = df[(df["class"] == cls) & (df["band"] == band)].iloc[0]
            tot = int(r["n_modules"])
            exp = EXPECTED_TOTAL_TRA_RES_ANC[cls][band]
            ok = tot == exp
            if not ok:
                pass_all = False
            rows.append({"class": cls, "band": band, "total": tot, "exp": exp,
                         "status": "PASS" if ok else "FAIL"})
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"\nFigure 1 sanity (trace/reset/anchor only): {'PASS' if pass_all else 'FAIL'}")
    rearr_rows = df[df["class"] == "rearrange"].sort_values("band")
    print("\nrearrange residual per band:")
    for _, r in rearr_rows.iterrows():
        print(f"  {r['band']:11s}  total={int(r['n_modules']):3d}  "
              f"cells={int(r['n_cells_with_module']):2d}/{int(r['cohort_n'])}")
    return pass_all


def render(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.375))
    n_classes = len(CLASSES)
    bar_w = 0.16
    group_centers = np.arange(len(BANDS))
    offsets = (np.arange(n_classes) - (n_classes - 1) / 2.0) * bar_w

    module_max = int(df[df["class"] != "diffuse"]["n_modules"].max())
    ymax = max(module_max * 1.22, 110)

    cohort_totals = {}
    for ci, cls in enumerate(CLASSES):
        heights = []
        cells_list = []
        cohorts = []
        for band in BANDS:
            r = df[(df["class"] == cls) & (df["band"] == band)].iloc[0]
            heights.append(int(r["n_modules"]))
            cells_list.append(int(r["n_cells_with_module"]))
            cohorts.append(int(r["cohort_n"]))
        cohort_totals[cls] = sum(heights)
        xs = group_centers + offsets[ci]
        hatch = "///" if cls == "diffuse" else None
        if cls == "diffuse":
            clipped_heights = [min(h, ymax) for h in heights]
            ax.bar(
                xs, clipped_heights, width=bar_w,
                color=CLASS_COLORS[cls], edgecolor="white", linewidth=0.4,
                label=f"{CLASS_LABEL[cls]} ({cohort_totals[cls]})",
                hatch=hatch, zorder=3,
            )
            for x, h_actual in zip(xs, heights):
                ax.text(
                    x, ymax * 0.985, f"↑{h_actual}",
                    ha="center", va="top", fontsize=7, color="#333",
                    rotation=90, zorder=5,
                    bbox=dict(facecolor="white", edgecolor="none",
                              alpha=0.85, pad=0.6),
                )
        else:
            ax.bar(
                xs, heights, width=bar_w,
                color=CLASS_COLORS[cls], edgecolor="white", linewidth=0.4,
                label=f"{CLASS_LABEL[cls]} ({cohort_totals[cls]})",
                zorder=3,
            )
            for x, h, c, n in zip(xs, heights, cells_list, cohorts):
                ax.text(
                    x, h + 1.5, f"{c}/{n}",
                    ha="center", va="bottom", fontsize=7, color="#333",
                    rotation=90, zorder=4,
                )

    ax.set_xticks(group_centers)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=11)
    ax.set_ylabel("module count (4 classes)  ·  leaf count (diffuse, ↑clipped)",
                  fontsize=9)
    ax.tick_params(axis="y", labelsize=9)
    ax.yaxis.set_major_locator(plt.MultipleLocator(25))
    ax.grid(axis="y", color="#dcdcdc", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for sn in ("top", "right"):
        ax.spines[sn].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.set_xlim(-0.5, len(BANDS) - 0.5)
    ax.set_ylim(0, ymax)

    ax.legend(loc="upper right", frameon=False, fontsize=8, handlelength=1.4, ncol=1)

    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT_FIG.relative_to(ROOT)}")


def main() -> int:
    residuals = compute_residuals_all_cells(force=False)
    df = collect(residuals)
    if not sanity_check(df):
        print("sanity check failed; not rendering")
        return 1
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV.relative_to(ROOT)}")
    render(df)
    print(f"figure dimensions: 7.0 x 4.375 in (1.6:1 aspect)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
