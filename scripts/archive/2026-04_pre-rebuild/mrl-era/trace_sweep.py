#!/usr/bin/env python3
"""Step 2 — sweep all (patient, band) cells, one trace candidate each.

For every (patient, band) under ``imcoh_abs``:
  1. Enumerate internal nodes ``v ∈ T_rest_post`` with ``size(v) ≥ k_min``.
  2. For each ``v``: compute
        J_test(v) = J*(leaves(v), T_task_test)   (ImCoh task_test tree)
        J_pre(v)  = J*(leaves(v), T_rest_pre)
        score(v)  = J_test(v) − J_pre(v)
  3. Pick the candidate ``v*`` as the LARGEST subtree with
     ``score ≥ score_floor``; fall back to global argmax if none.

Output: ``data/reports/imcoh_mrl/trace_sweep.csv`` — one row per cell.

This is the cohort-level data product that Step 3 (per-band dotplot)
and Step 4 (multi-patient dendrogram reprojection) read from.
"""
from __future__ import annotations

import argparse
import json
from typing import Iterable

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


PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]


def load_Z(patient: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(patient, phase, band, fc_method="imcoh_abs",
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


def find_candidate(Z_post: np.ndarray, Z_test: np.ndarray, Z_pre: np.ndarray,
                   k_min: int, score_floor: float) -> dict | None:
    V_post = [v for v in tree_internal_nodes(Z_post)
              if k_min <= v["size"] <= Z_post.shape[0]]
    if not V_post:
        return None
    V_test = tree_internal_nodes(Z_test)
    V_pre  = tree_internal_nodes(Z_pre)

    rows = []
    for v in V_post:
        j_test = best_match(v["leaves"], V_test)
        j_pre  = best_match(v["leaves"], V_pre)
        rows.append({
            "v": v, "size": v["size"], "h_rel": v["h_rel"],
            "j_test": j_test, "j_pre": j_pre,
            "score": j_test - j_pre,
        })

    above = [r for r in rows if r["score"] >= score_floor]
    chosen = (max(above, key=lambda r: (r["size"], r["score"])) if above
              else max(rows, key=lambda r: r["score"]))
    chosen["passed_floor"] = bool(above)
    return chosen


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", default=None,
                    help="Comma-sep IDs (default PATIENTS_4PHASE).")
    ap.add_argument("--bands", default=None,
                    help="Comma-sep band names (default all 6).")
    ap.add_argument("--k-min", type=int, default=5)
    ap.add_argument("--score-floor", type=float, default=0.5)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    patients = (args.patients.split(",") if args.patients
                else list(PATIENTS_4PHASE))
    bands = args.bands.split(",") if args.bands else list(BRAIN_BANDS_NAMES)
    k_min = int(args.k_min)
    floor = float(args.score_floor)

    out_dir = REPORTS_ROOT / "imcoh_mrl"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    skipped: list[tuple[str, str]] = []
    for pat in patients:
        for band in bands:
            Zs = {phase: load_Z(pat, phase, band) for phase in PHASES}
            if any(z is None for z in Zs.values()):
                skipped.append((pat, band))
                if args.verbose:
                    print(f"[sweep] SKIP {pat:7s} {band:10s} — incomplete LRG cache")
                continue
            n_post = Zs["rest_post"].shape[0] + 1
            if not all(z.shape[0] + 1 == n_post for z in Zs.values()):
                skipped.append((pat, band))
                continue

            cand = find_candidate(Zs["rest_post"], Zs["task_test"],
                                  Zs["rest_pre"], k_min=k_min,
                                  score_floor=floor)
            if cand is None:
                skipped.append((pat, band))
                continue
            v = cand["v"]
            rows.append({
                "patient": pat,
                "band": band,
                "n_leaves": n_post,
                "size": int(cand["size"]),
                "h_rel": float(cand["h_rel"]),
                "j_test": float(cand["j_test"]),
                "j_pre": float(cand["j_pre"]),
                "score": float(cand["score"]),
                "passed_floor": cand["passed_floor"],
                "leaves": json.dumps(sorted(int(x) for x in v["leaves"])),
            })
            if args.verbose:
                print(f"[sweep]   {pat:7s} {band:10s}  size={cand['size']:3d}  "
                      f"h_rel={cand['h_rel']:.3f}  J_test={cand['j_test']:.2f}  "
                      f"J_pre={cand['j_pre']:.2f}  score={cand['score']:+.3f}  "
                      f"floor={'yes' if cand['passed_floor'] else 'no'}")

    df = pd.DataFrame(rows)
    out_path = out_dir / "trace_sweep.csv"
    df.to_csv(out_path, index=False)
    print()
    print(f"[sweep] wrote {out_path}  ({len(df)} cells)")
    if skipped:
        print(f"[sweep] skipped cells: {skipped}")
    print()
    print("[sweep] per-band summary (n=cells, median score, frac score≥floor):")
    print(f"[sweep]   {'band':<10s}  {'n':>3s}  {'med_score':>9s}  "
          f"{'%≥floor':>8s}  {'med_size':>8s}")
    for band in bands:
        sub = df[df["band"] == band]
        if sub.empty:
            print(f"[sweep]   {band:<10s}  0  --")
            continue
        med_score = sub["score"].median()
        frac_above = (sub["passed_floor"]).mean()
        med_size = sub["size"].median()
        print(f"[sweep]   {band:<10s}  {len(sub):>3d}  {med_score:>+9.3f}  "
              f"{frac_above*100:>7.1f}%  {med_size:>8.1f}")


if __name__ == "__main__":
    main()
