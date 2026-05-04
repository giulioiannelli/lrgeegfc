#!/usr/bin/env python3
"""Audit Step 16 — Restricted Δ_VI(k) on a leaf subset.

Reconciliation primitive for the MSPC ↔ canonical Δ_VI mismatch:

  Whole-brain Δ_VI says γ_l silent, MSPC anatomy says γ_l strongest.
  α/β go the other way (Δ_VI sees, MSPC argmax doesn't).

This script restricts the partition VI computation to a leaf subset M
and reports Δ_VI(k) on M only:

    Δ_VI(M, k) = VI(P_pre|M, P_post|M) − VI(P_test|M, P_post|M)

where the partitions are induced by the WHOLE-tree fcluster at scale k
(not refit). High Δ_VI(M, k) means P_test|M is closer to P_post|M than
P_pre|M is — the trace pattern. The subset M is fed by ``--mask``:

  --mask mspc_trace : leaves with MSPC dominant=trace (per pat × band)
                       → tests whether MSPC's per-leaf trace claim
                         creates a partition shift on those leaves.
  --mask lobe       : leaves grouped by Desikan-Killiany lobe
                       → tests whether canonical Δ_VI α/β trace
                         localises to specific anatomy.

Outputs:
  data/audit/per_patient_hierarchy_mspc/restricted_vi/{mask}/
    per_cell.csv                    — pat × band [× lobe] × k rows
    cohort_summary.csv              — ≥8/10 cohort runs per band [× lobe]
    fig_delta_vi_curves_{band}.pdf  — per-band Δ_VI(k) cohort + ribbon
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_08_per_patient_hierarchy import (  # noqa: E402
    BAND_ORDER, PHASE_ORDER, PATIENTS_4PHASE, BRAIN_BAND_TEX_DICT,
    _load_Z,
)
from lrg_eegfc.utils.metrics.tree import partition_vi_on_subset  # noqa: E402
from lrg_eegfc.utils.io.regions import load_channel_regions  # noqa: E402

OUT_BASE = ROOT / "data" / "audit" / "per_patient_hierarchy_mspc" / "restricted_vi"
LOBES_KEEP = ["frontal", "parietal", "temporal", "occipital",
              "cingulate", "insula", "subcortical"]
N_MIN_SUBSET = 5     # subsets smaller than this give degenerate VI
COHORT_GATE = 0.8    # ≥8/10 patients


def _build_mspc_trace_masks() -> dict[tuple[str, str], np.ndarray]:
    """Return mask[(pat, band)] → boolean array over leaves (True for trace)."""
    leaves = pd.read_csv(ROOT / "data" / "audit" /
                            "per_patient_hierarchy_mspc" / "leaf_assignment.csv")
    masks: dict[tuple[str, str], np.ndarray] = {}
    for (pat, band), grp in leaves.groupby(["patient", "band"]):
        n = grp.leaf_id.max() + 1
        m = np.zeros(n, dtype=bool)
        m[grp.loc[grp.dominant == "trace", "leaf_id"].values] = True
        masks[(pat, band)] = m
    return masks


def _build_lobe_masks() -> dict[tuple[str, str], np.ndarray]:
    """Return mask[(pat, lobe)] → boolean array over leaves (True for that lobe)."""
    masks: dict[tuple[str, str], np.ndarray] = {}
    for pat in PATIENTS_4PHASE:
        rdf = load_channel_regions(pat)
        for lobe in LOBES_KEEP:
            m = (rdf["lobe"].values == lobe)
            if m.sum() >= N_MIN_SUBSET:
                masks[(pat, lobe)] = m
    return masks


def _delta_vi_curve(Zs: dict[str, np.ndarray], mask: np.ndarray,
                       k_values: list[int]) -> list[dict]:
    """For each k, compute Δ_VI on the subset AND on the whole brain.
    Uses task_test as canonical `task` panel."""
    if Zs.get("rest_pre") is None or Zs.get("rest_post") is None:
        return []
    if Zs.get("task_test") is None:
        return []
    rows: list[dict] = []
    full_mask = np.ones_like(mask, dtype=bool)
    for k in k_values:
        lab_pre = fcluster(Zs["rest_pre"], k, criterion="maxclust")
        lab_test = fcluster(Zs["task_test"], k, criterion="maxclust")
        lab_post = fcluster(Zs["rest_post"], k, criterion="maxclust")

        d_pp_sub = partition_vi_on_subset(lab_pre, lab_post, mask)
        d_tp_sub = partition_vi_on_subset(lab_test, lab_post, mask)
        d_pp_full = partition_vi_on_subset(lab_pre, lab_post, full_mask)
        d_tp_full = partition_vi_on_subset(lab_test, lab_post, full_mask)

        rows.append({
            "k": k,
            "n_subset": d_pp_sub["n_subset"],
            "vi_pre_post": d_pp_sub["vi"],
            "vi_test_post": d_tp_sub["vi"],
            "delta_vi": d_pp_sub["vi"] - d_tp_sub["vi"],
            "vi_pre_post_full": d_pp_full["vi"],
            "vi_test_post_full": d_tp_full["vi"],
            "delta_vi_full": d_pp_full["vi"] - d_tp_full["vi"],
            "max_cluster_fraction_pre":
                d_pp_sub["max_cluster_fraction_1"],
            "max_cluster_fraction_test":
                d_tp_sub["max_cluster_fraction_1"],
            "max_cluster_fraction_post":
                d_pp_sub["max_cluster_fraction_2"],
        })
    return rows


def run_mspc_trace(out_dir: Path) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    masks = _build_mspc_trace_masks()
    rows: list[dict] = []
    for pat in PATIENTS_4PHASE:
        print(f"{pat} ...", flush=True)
        Zs_all = {(pat, b): {ph: _load_Z(pat, ph, b) for ph in PHASE_ORDER}
                  for b in BAND_ORDER}
        for band in BAND_ORDER:
            mask = masks.get((pat, band))
            if mask is None or mask.sum() < N_MIN_SUBSET:
                continue
            Zs = Zs_all[(pat, band)]
            if Zs["rest_post"] is None:
                continue
            N = Zs["rest_post"].shape[0] + 1
            k_values = list(range(2, max(N // 2, 3) + 1))
            curve = _delta_vi_curve(Zs, mask, k_values)
            for r in curve:
                rows.append({"patient": pat, "band": band, **r})
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "per_cell.csv", index=False)
    return df


def run_lobe(out_dir: Path) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    masks = _build_lobe_masks()
    rows: list[dict] = []
    for pat in PATIENTS_4PHASE:
        print(f"{pat} ...", flush=True)
        for band in BAND_ORDER:
            Zs = {ph: _load_Z(pat, ph, band) for ph in PHASE_ORDER}
            if Zs["rest_post"] is None:
                continue
            N = Zs["rest_post"].shape[0] + 1
            k_values = list(range(2, max(N // 2, 3) + 1))
            for lobe in LOBES_KEEP:
                mask = masks.get((pat, lobe))
                if mask is None or mask.sum() < N_MIN_SUBSET:
                    continue
                curve = _delta_vi_curve(Zs, mask, k_values)
                for r in curve:
                    rows.append({"patient": pat, "band": band,
                                  "lobe": lobe, **r})
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "per_cell.csv", index=False)
    return df


def cohort_summary(df: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    """For each (band [× lobe]) and k: count patients with delta_vi > 0."""
    grouped = df.groupby(by + ["k"]).agg(
        n_pat=("patient", "nunique"),
        n_pos=("delta_vi", lambda s: int((s > 0).sum())),
        median_delta_vi=("delta_vi", "median"),
        q25_delta_vi=("delta_vi", lambda s: float(s.quantile(0.25))),
        q75_delta_vi=("delta_vi", lambda s: float(s.quantile(0.75))),
        median_delta_vi_full=("delta_vi_full", "median"),
        q25_delta_vi_full=("delta_vi_full", lambda s: float(s.quantile(0.25))),
        q75_delta_vi_full=("delta_vi_full", lambda s: float(s.quantile(0.75))),
    ).reset_index()
    grouped["frac_pos"] = grouped["n_pos"] / grouped["n_pat"].replace(0, np.nan)
    grouped["cohort_pos"] = (grouped["frac_pos"] >= COHORT_GATE).astype(int)
    return grouped


def figure_per_band(summary: pd.DataFrame, out_dir: Path,
                       group_col: str | None = None) -> None:
    """Δ_VI(k) median + IQR per band; trace-restricted vs whole-brain.
    Filled circles = ≥80% cohort positivity (frac_pos ≥ 0.8 over patients
    with subsets ≥ 5 leaves)."""
    bands = BAND_ORDER
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.5), sharex=False)
    axes = axes.flatten()
    for ax, band in zip(axes, bands):
        sub = summary[summary["band"] == band]
        if sub.empty:
            ax.set_title(f"{BRAIN_BAND_TEX_DICT.get(band, band)}: no data",
                          fontsize=10)
            ax.axis("off")
            continue
        if group_col is None or group_col not in sub.columns:
            grp_iter = [(None, sub)]
            colors = ["#d62728"]
        else:
            unique_groups = sorted(sub[group_col].unique())
            grp_iter = [(g, sub[sub[group_col] == g]) for g in unique_groups]
            cmap = plt.get_cmap("tab10")
            colors = [cmap(i % 10) for i in range(len(unique_groups))]
        for (g, gdf), color in zip(grp_iter, colors):
            gdf = gdf.sort_values("k")
            # restricted (subset)
            ax.fill_between(gdf["k"], gdf["q25_delta_vi"], gdf["q75_delta_vi"],
                              color=color, alpha=0.18, linewidth=0)
            ax.plot(gdf["k"], gdf["median_delta_vi"], color=color,
                     label=(f"{g} (subset)" if g is not None else "subset"),
                     linewidth=1.6)
            cohort = gdf[gdf["cohort_pos"] == 1]
            if not cohort.empty:
                ax.scatter(cohort["k"], cohort["median_delta_vi"],
                            color=color, edgecolors="black", linewidths=0.6,
                            s=22, zorder=5)
            # whole-brain baseline (only once per ax for top-level mask)
            if group_col is None:
                ax.fill_between(gdf["k"], gdf["q25_delta_vi_full"],
                                  gdf["q75_delta_vi_full"],
                                  color="#888", alpha=0.18, linewidth=0)
                ax.plot(gdf["k"], gdf["median_delta_vi_full"],
                         color="#444", linestyle="--", linewidth=1.2,
                         label="whole-brain")
        ax.axhline(0.0, color="#888", linewidth=0.6, linestyle=":")
        ax.set_title(BRAIN_BAND_TEX_DICT.get(band, band), fontsize=11)
        ax.set_xlabel("k", fontsize=9); ax.set_ylabel("Δ_VI (nats)", fontsize=9)
        ax.legend(loc="best", fontsize=7, frameon=False)
    fig.suptitle("Restricted Δ_VI(k) on MSPC trace subset (red) vs whole brain "
                  "(grey dashed) — median + IQR", fontsize=11)
    fig.tight_layout()
    out = out_dir / "fig_delta_vi_curves.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--mask", choices=["mspc_trace", "lobe"],
                     default="mspc_trace")
    args = ap.parse_args()

    out_dir = OUT_BASE / args.mask
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[audit_16] mask={args.mask}; writing to {out_dir.relative_to(ROOT)}/")
    if args.mask == "mspc_trace":
        df = run_mspc_trace(out_dir)
        summary = cohort_summary(df, by=["band"])
        summary.to_csv(out_dir / "cohort_summary.csv", index=False)
        figure_per_band(summary, out_dir, group_col=None)
    else:
        df = run_lobe(out_dir)
        summary = cohort_summary(df, by=["band", "lobe"])
        summary.to_csv(out_dir / "cohort_summary.csv", index=False)
        figure_per_band(summary, out_dir, group_col="lobe")

    print(f"\n[audit_16] {len(df)} cell-rows; "
          f"{int(summary['cohort_pos'].sum())} cohort-positive cells "
          f"(≥{int(COHORT_GATE*10)}/10)")
    print(f"Wrote {out_dir.relative_to(ROOT)}/per_cell.csv, "
          f"cohort_summary.csv, fig_delta_vi_curves.pdf")


if __name__ == "__main__":
    main()
