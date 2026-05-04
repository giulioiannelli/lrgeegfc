#!/usr/bin/env python3
"""Task-anchored module (TAM) persistence sweep.

A TAM is a leafset that is *jointly* a coherent subtree in BOTH
``task_learn`` and ``task_test``. Concretely:

    Candidate v ∈ Internal(T_task_learn) with k_min ≤ size(v) ≤ N − 1.
    leaves(v) is a TAM iff
        J*(leaves(v), T_task_test) ≥ J_anchor.

The candidate set is therefore defined ONLY by task data — no
selection bias from rest_post structure.

For each TAM with leafset S we compute

    J_pre  = J*(S, T_rest_pre)        (was the module there beforehand?)
    J_post = J*(S, T_rest_post)       (does it persist after task?)
    Δ(S)   = J_post − J_pre           (persistence asymmetry)

A TAM is a "trace TAM" iff J_pre < J_low AND J_post ≥ J_high
(default 0.5 / 0.7 — discrete rule).

Per (patient, band) outputs:
    n_TAMs        — count of TAMs in that cell
    mean_delta    — average J_post − J_pre across TAMs (continuous)
    n_trace       — count of TAMs satisfying the discrete rule
    frac_trace    — n_trace / max(n_TAMs, 1)

Long-format output (one row per TAM): ``trace_tam_per_module.csv``.
Cell summary: ``trace_tam_summary.csv``.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, PATIENTS_4PHASE
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.utils.metrics.tree import (
    tree_internal_nodes,
    jaccard_leafsets,
)
from lrg_eegfc.workflow.lrg import load_lrg_result


def load_Z(pat: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, fc_method="imcoh_abs",
                            cache_root=IMCOH_LRG_CACHE)
    except Exception:
        return None
    if r is None:
        return None
    return np.asarray(r.linkage_matrix)


def best_match(S: frozenset[int], V_other: list[dict]) -> float:
    if not V_other:
        return 0.0
    best_j = 0.0
    for u in V_other:
        j = jaccard_leafsets(S, u["leaves"])
        if j > best_j:
            best_j = j
            if best_j >= 1.0:
                return 1.0
    return best_j


def find_tams(Z_TL: np.ndarray, Z_TT: np.ndarray,
              k_min: int, J_anchor: float) -> list[dict]:
    V_TL = [v for v in tree_internal_nodes(Z_TL)
            if k_min <= v["size"] <= Z_TL.shape[0]]
    V_TT = tree_internal_nodes(Z_TT)
    tams = []
    for v in V_TL:
        j_anchor = best_match(v["leaves"], V_TT)
        if j_anchor >= J_anchor:
            tams.append({
                "leaves": v["leaves"],
                "size": int(v["size"]),
                "h_rel_TL": float(v["h_rel"]),
                "J_anchor": float(j_anchor),
            })
    return tams


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", default=None)
    ap.add_argument("--bands", default=None)
    ap.add_argument("--k-min", type=int, default=5)
    ap.add_argument("--j-anchor", type=float, default=0.7,
                    help="J_TL↔TT match threshold for TAM (default 0.7).")
    ap.add_argument("--j-low", type=float, default=0.5,
                    help="rest_pre absence threshold (default 0.5).")
    ap.add_argument("--j-high", type=float, default=0.7,
                    help="rest_post presence threshold (default 0.7).")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    patients = (args.patients.split(",") if args.patients
                else list(PATIENTS_4PHASE))
    bands = args.bands.split(",") if args.bands else list(BRAIN_BANDS_NAMES)
    k_min = int(args.k_min)
    J_anchor = float(args.j_anchor)
    J_low = float(args.j_low)
    J_high = float(args.j_high)

    out_dir = REPORTS_ROOT / "imcoh_mrl"
    out_dir.mkdir(parents=True, exist_ok=True)

    per_module: list[dict] = []
    summary: list[dict] = []

    for pat in patients:
        for band in bands:
            Z_TL = load_Z(pat, "task_learn", band)
            Z_TT = load_Z(pat, "task_test",  band)
            Z_pre  = load_Z(pat, "rest_pre",  band)
            Z_post = load_Z(pat, "rest_post", band)
            if any(z is None for z in (Z_TL, Z_TT, Z_pre, Z_post)):
                continue
            n_post = Z_post.shape[0] + 1
            if not (Z_TL.shape[0] + 1 == Z_TT.shape[0] + 1 == n_post
                    == Z_pre.shape[0] + 1):
                continue

            tams = find_tams(Z_TL, Z_TT, k_min, J_anchor)
            V_pre  = tree_internal_nodes(Z_pre)
            V_post = tree_internal_nodes(Z_post)

            mods = []
            for tam in tams:
                S = tam["leaves"]
                J_pre  = best_match(S, V_pre)
                J_post = best_match(S, V_post)
                row = {
                    "patient": pat, "band": band,
                    "size": tam["size"], "h_rel_TL": tam["h_rel_TL"],
                    "J_anchor": tam["J_anchor"],
                    "J_pre": J_pre, "J_post": J_post,
                    "delta": J_post - J_pre,
                    "is_trace": (J_pre < J_low) and (J_post >= J_high),
                    "leaves": json.dumps(sorted(int(x) for x in S)),
                }
                per_module.append(row)
                mods.append(row)

            n_tams = len(mods)
            mean_d = float(np.mean([m["delta"] for m in mods])) if mods else np.nan
            mean_pre = float(np.mean([m["J_pre"] for m in mods])) if mods else np.nan
            mean_post = float(np.mean([m["J_post"] for m in mods])) if mods else np.nan
            n_trace = int(sum(1 for m in mods if m["is_trace"]))
            frac_trace = (n_trace / n_tams) if n_tams else np.nan
            summary.append({
                "patient": pat, "band": band, "n_TAMs": n_tams,
                "mean_J_pre": mean_pre, "mean_J_post": mean_post,
                "mean_delta": mean_d, "n_trace": n_trace,
                "frac_trace": frac_trace,
            })
            if args.verbose:
                print(f"[tam] {pat:7s} {band:10s}  n_TAMs={n_tams:3d}  "
                      f"mean_pre={mean_pre:.2f}  mean_post={mean_post:.2f}  "
                      f"mean_Δ={mean_d:+.3f}  trace={n_trace}/{n_tams}")

    pm_df = pd.DataFrame(per_module)
    sm_df = pd.DataFrame(summary)
    pm_path = out_dir / "trace_tam_per_module.csv"
    sm_path = out_dir / "trace_tam_summary.csv"
    pm_df.to_csv(pm_path, index=False)
    sm_df.to_csv(sm_path, index=False)
    print()
    print(f"[tam] wrote {pm_path}  ({len(pm_df)} TAMs)")
    print(f"[tam] wrote {sm_path}  ({len(sm_df)} cells)")
    print()
    print(f"[tam] thresholds: J_anchor={J_anchor}, J_low={J_low}, "
          f"J_high={J_high}, k_min={k_min}")
    print()
    print("[tam] per-band cohort summary:")
    print(f"  {'band':<10s}  {'cells':>5s}  {'TAMs':>5s}  "
          f"{'med_Δ':>7s}  {'mean_Δ':>7s}  "
          f"{'n_pat_Δ>0':>10s}  {'mean_frac_trace':>15s}")
    for band in bands:
        sub = sm_df[sm_df["band"] == band]
        if sub.empty:
            print(f"  {band:<10s}  0")
            continue
        n_cells = len(sub)
        n_tams_total = int(sub["n_TAMs"].sum())
        med = float(sub["mean_delta"].median())
        mn  = float(sub["mean_delta"].mean())
        n_pos = int((sub["mean_delta"] > 0).sum())
        m_frac = float(sub["frac_trace"].mean())
        print(f"  {band:<10s}  {n_cells:>5d}  {n_tams_total:>5d}  "
              f"{med:>+.3f}  {mn:>+.3f}  {n_pos:>4d}/{n_cells:<3d}     "
              f"{m_frac:>+.3f}")


if __name__ == "__main__":
    main()
