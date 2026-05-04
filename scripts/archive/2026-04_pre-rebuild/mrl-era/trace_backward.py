#!/usr/bin/env python3
"""Backward control for the trace effect.

Forward (already done):
    candidate v* ∈ T_rest_post
    score_forward = J*(v*, T_task_test) − J*(v*, T_rest_pre)

Backward (this script):
    candidate v_b* ∈ T_rest_pre
    score_backward = J*(v_b*, T_task_learn) − J*(v_b*, T_rest_post)

Symmetry of the two scores in expectation:
  - If the effect is **task-induced asymmetric**: rest_post inherits
    task structure that rest_pre was untouched by, so forward >> 0
    while backward ≈ 0 (pre subtrees are anatomical baselines, equally
    distant from learn and from post).
  - If the effect is **temporal/session-order**: temporally-adjacent
    phases are simply more similar; both scores should be of similar
    magnitude.

Output: ``data/reports/imcoh_mrl/trace_backward_sweep.csv`` and the
500-perm size-matched null in ``trace_backward_null.csv``.
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


def find_candidate(Z_anchor: np.ndarray,
                   Z_close: np.ndarray, Z_far: np.ndarray,
                   k_min: int, score_floor: float) -> dict | None:
    """Pick subtree v ∈ T_anchor maximising
        score(v) = J*(v, T_close) − J*(v, T_far)
    using the same "largest with score≥floor; else argmax" logic as
    the forward sweep."""
    V_anchor = [v for v in tree_internal_nodes(Z_anchor)
                if k_min <= v["size"] <= Z_anchor.shape[0]]
    if not V_anchor:
        return None
    V_close = tree_internal_nodes(Z_close)
    V_far   = tree_internal_nodes(Z_far)

    rows = []
    for v in V_anchor:
        j_close = best_match(v["leaves"], V_close)
        j_far   = best_match(v["leaves"], V_far)
        rows.append({
            "v": v, "size": v["size"], "h_rel": v["h_rel"],
            "j_close": j_close, "j_far": j_far,
            "score": j_close - j_far,
        })

    above = [r for r in rows if r["score"] >= score_floor]
    chosen = (max(above, key=lambda r: (r["size"], r["score"])) if above
              else max(rows, key=lambda r: r["score"]))
    chosen["passed_floor"] = bool(above)
    return chosen


def perm_null(S_size: int, n_leaves: int, V_close: list[dict],
              V_far: list[dict], n_perm: int, rng: np.random.Generator
              ) -> np.ndarray:
    """500 random size-matched leafsets → null score distribution."""
    out = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        idx = rng.choice(n_leaves, size=S_size, replace=False)
        S = frozenset(int(x) for x in idx)
        j_c = best_match(S, V_close)
        j_f = best_match(S, V_far)
        out[i] = j_c - j_f
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", default=None)
    ap.add_argument("--bands", default=None)
    ap.add_argument("--k-min", type=int, default=5)
    ap.add_argument("--score-floor", type=float, default=0.5)
    ap.add_argument("--n-perm", type=int, default=500)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    patients = (args.patients.split(",") if args.patients
                else list(PATIENTS_4PHASE))
    bands = args.bands.split(",") if args.bands else list(BRAIN_BANDS_NAMES)
    k_min = int(args.k_min)
    floor = float(args.score_floor)
    n_perm = int(args.n_perm)
    rng = np.random.default_rng(args.seed)

    out_dir = REPORTS_ROOT / "imcoh_mrl"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    null_rows: list[dict] = []
    for pat in patients:
        for band in bands:
            Zs = {phase: load_Z(pat, phase, band) for phase in
                  ["rest_pre", "task_learn", "task_test", "rest_post"]}
            if any(z is None for z in Zs.values()):
                continue
            n_post = Zs["rest_post"].shape[0] + 1
            if not all(z.shape[0] + 1 == n_post for z in Zs.values()):
                continue

            # Backward: anchor = rest_pre; close = task_learn; far = rest_post.
            cand = find_candidate(Zs["rest_pre"], Zs["task_learn"],
                                  Zs["rest_post"], k_min, floor)
            if cand is None:
                continue
            v = cand["v"]
            row = {
                "patient": pat, "band": band, "n_leaves": n_post,
                "size": int(cand["size"]), "h_rel": float(cand["h_rel"]),
                "j_close": float(cand["j_close"]),  # = J*(pre_subtree, learn)
                "j_far":   float(cand["j_far"]),    # = J*(pre_subtree, post)
                "score": float(cand["score"]),
                "passed_floor": cand["passed_floor"],
                "leaves": json.dumps(sorted(int(x) for x in v["leaves"])),
            }
            rows.append(row)

            # Null on this cell.
            V_close = tree_internal_nodes(Zs["task_learn"])
            V_far   = tree_internal_nodes(Zs["rest_post"])
            null_scores = perm_null(int(cand["size"]), n_post, V_close,
                                    V_far, n_perm, rng)
            null_rows.append({
                "patient": pat, "band": band,
                "size": int(cand["size"]), "n_leaves": n_post,
                "observed_score": float(cand["score"]),
                "null_mean": float(null_scores.mean()),
                "null_std":  float(null_scores.std()),
                "null_p95":  float(np.quantile(null_scores, 0.95)),
                "p_empirical": float((null_scores >= cand["score"]).sum() + 1)
                                / (n_perm + 1),
                "z_observed": float((cand["score"] - null_scores.mean())
                                    / (null_scores.std() + 1e-12)),
            })
            if args.verbose:
                print(f"[back] {pat:7s} {band:10s}  size={int(cand['size']):3d}  "
                      f"score={cand['score']:+.3f}  "
                      f"null_mean={null_scores.mean():+.3f}  "
                      f"p={null_rows[-1]['p_empirical']:.3f}  "
                      f"z={null_rows[-1]['z_observed']:+.2f}")

    df = pd.DataFrame(rows)
    df_null = pd.DataFrame(null_rows)
    sweep_path = out_dir / "trace_backward_sweep.csv"
    null_path = out_dir / "trace_backward_null.csv"
    df.to_csv(sweep_path, index=False)
    df_null.to_csv(null_path, index=False)
    print()
    print(f"[back] wrote {sweep_path}  ({len(df)} cells)")
    print(f"[back] wrote {null_path}   ({len(df_null)} cells, {n_perm} perms each)")
    print()
    print("[back] per-band summary:")
    print(f"[back]   {'band':<10s}  {'n':>3s}  {'med_obs':>7s}  "
          f"{'med_null':>8s}  {'med_z':>6s}  {'p<.05':>6s}  {'p<.01':>6s}")
    for band in bands:
        sub = df_null[df_null["band"] == band]
        if sub.empty:
            print(f"[back]   {band:<10s}  0")
            continue
        med_obs = sub["observed_score"].median()
        med_null = sub["null_mean"].median()
        med_z = sub["z_observed"].median()
        n_05 = int((sub["p_empirical"] < 0.05).sum())
        n_01 = int((sub["p_empirical"] < 0.01).sum())
        print(f"[back]   {band:<10s}  {len(sub):>3d}  {med_obs:>+.3f}  "
              f"{med_null:>+.3f}  {med_z:>+.2f}  "
              f"{n_05}/{len(sub):<3d}  {n_01}/{len(sub):<3d}")


if __name__ == "__main__":
    main()
