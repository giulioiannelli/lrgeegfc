#!/usr/bin/env python3
"""(b) Size-matched permutation null on trace scores.

For each row of ``trace_sweep.csv`` (one (patient, band) cell), draw
``--n-perm`` (default 500) random subsets of the *same size* as the
observed candidate from the patient's leaf set. For each random
subset ``S_rand``:

    score_rand = J*(S_rand, T_task_test) − J*(S_rand, T_rest_pre)

The empirical null is the distribution of ``score_rand`` over the
permutations. Per-cell p-value = fraction of null draws ≥ observed
score. The observed candidate's leafset was filtered by score-floor
and size-argmax; the random size-matched subsets ignore that filter
and ask "what's the score of an arbitrary group of leaves of the
same size?". A high observed score relative to this null = the
specific contacts cluster together in task and post in a way that
random groupings of the same size do not.

Output: ``data/reports/imcoh_mrl/trace_null.csv``
"""
from __future__ import annotations

import argparse

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


PHASES = ["rest_pre", "task_test"]


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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=500)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    n_perm = int(args.n_perm)

    sweep_path = REPORTS_ROOT / "imcoh_mrl" / "trace_sweep.csv"
    sweep = pd.read_csv(sweep_path)

    out_rows: list[dict] = []
    for _, row in sweep.iterrows():
        pat = str(row["patient"])
        band = str(row["band"])
        size = int(row["size"])
        n_leaves = int(row["n_leaves"])
        observed = float(row["score"])

        Z_test = load_Z(pat, "task_test", band)
        Z_pre  = load_Z(pat, "rest_pre",  band)
        if Z_test is None or Z_pre is None:
            continue
        V_test = tree_internal_nodes(Z_test)
        V_pre  = tree_internal_nodes(Z_pre)

        leaves = np.arange(n_leaves)
        scores_null = np.empty(n_perm, dtype=float)
        for i in range(n_perm):
            idx = rng.choice(n_leaves, size=size, replace=False)
            S = frozenset(int(x) for x in idx)
            j_test = best_match(S, V_test)
            j_pre  = best_match(S, V_pre)
            scores_null[i] = j_test - j_pre

        null_mean = float(scores_null.mean())
        null_std  = float(scores_null.std())
        null_p95  = float(np.quantile(scores_null, 0.95))
        p_emp     = float((scores_null >= observed).sum() + 1) / (n_perm + 1)
        z_obs     = (observed - null_mean) / (null_std + 1e-12)

        out_rows.append({
            "patient": pat, "band": band,
            "size": size, "n_leaves": n_leaves,
            "observed_score": observed,
            "null_mean": null_mean,
            "null_std": null_std,
            "null_p95": null_p95,
            "p_empirical": p_emp,
            "z_observed": z_obs,
        })
        if args.verbose:
            print(f"[null] {pat:7s} {band:10s}  size={size:3d}  "
                  f"obs={observed:+.3f}  null_mean={null_mean:+.3f}  "
                  f"null_p95={null_p95:+.3f}  p={p_emp:.3f}  z={z_obs:+.2f}")

    df = pd.DataFrame(out_rows)
    out_path = REPORTS_ROOT / "imcoh_mrl" / "trace_null.csv"
    df.to_csv(out_path, index=False)
    print()
    print(f"[null] wrote {out_path}  ({len(df)} cells, {n_perm} perms each)")
    print()
    print("[null] per-band summary:")
    print(f"[null]   {'band':<10s}  {'n':>3s}  {'med_obs':>7s}  "
          f"{'med_null':>8s}  {'med_z':>6s}  {'p<.05':>6s}  {'p<.01':>6s}")
    for band in BRAIN_BANDS_NAMES:
        sub = df[df["band"] == band]
        if sub.empty:
            print(f"[null]   {band:<10s}  0")
            continue
        med_obs = sub["observed_score"].median()
        med_null = sub["null_mean"].median()
        med_z = sub["z_observed"].median()
        n_05 = int((sub["p_empirical"] < 0.05).sum())
        n_01 = int((sub["p_empirical"] < 0.01).sum())
        print(f"[null]   {band:<10s}  {len(sub):>3d}  {med_obs:>+.3f}  "
              f"{med_null:>+.3f}  {med_z:>+.2f}  "
              f"{n_05}/{len(sub):<3d}  {n_01}/{len(sub):<3d}")


if __name__ == "__main__":
    main()
