#!/usr/bin/env python3
"""MRL — Module-Retention Landscape.

For every internal node ``v`` of the pooled task-phase trees
``V_task = V(T_task_learn) ∪ V(T_task_test)`` we record:

  - retained(v; J_min) := (J*(leaves(v), T_rest_pre)  < J_min)
                         ∧ (J*(leaves(v), T_rest_post) ≥ J_min)
  - ΔJ(v) := J*(leaves(v), T_rest_post) − J*(leaves(v), T_rest_pre)

per (patient, band).  No synthetic query grid, no tolerance — every
node sits at its native fractional merge height ``h_rel(v) = h(v) /
dmax(T)``.

Cohort views aggregate to log-spaced ``h_rel`` bins **at plot time only**
(this script also writes the binned cohort CSVs as a convenience for
``fig_mrl_landscape.py``).

See the scope report for full mathematical detail:
``.agents/guides/task-persistence-investigation/2026-04-25_module-retention-landscape.md``

Outputs (under ``data/reports/imcoh_mrl/``):
  - ``mrl_per_node.csv``       — long format, every (p, b, v) row.
  - ``mrl_cohort_step.csv``    — cohort `M̄_step` and `Π` per
                                 (band, bin, J_min).
  - ``mrl_cohort_smooth.csv``  — cohort `M̄_smooth` per (band, bin)
                                 (no J_min dependence).
"""
from __future__ import annotations

import argparse
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
    h_log_grid,
)
from lrg_eegfc.workflow.lrg import load_lrg_result


# ─────────────────────────── helpers ───────────────────────────


def load_Z(pat: str, phase: str, band: str) -> np.ndarray | None:
    """Load LRG linkage for (pat, phase, band) under imcoh_abs.

    Returns ``None`` on any load failure (e.g. missing cache for a
    vendor-corrupt phase). Callers must skip the (p, b) cell if any
    of the four phases is unavailable.
    """
    try:
        r = load_lrg_result(pat, phase, band, fc_method="imcoh_abs",
                            cache_root=IMCOH_LRG_CACHE)
    except Exception:
        return None
    if r is None:
        return None
    return np.asarray(r.linkage_matrix)


def best_match(S: frozenset[int], others: list[dict]) -> float:
    """Return ``J*(S, T) = max_u J(S, leaves(u))`` over ``others``.

    ``others`` is the list returned by ``tree_internal_nodes`` for the
    target tree, optionally pre-filtered by ``size``. Returns ``0.0``
    if ``others`` is empty (no candidate to match).
    """
    if not others:
        return 0.0
    best = 0.0
    for u in others:
        j = jaccard_leafsets(S, u["leaves"])
        if j > best:
            best = j
            if best >= 1.0:
                return 1.0
    return best


def pooled_task_nodes(Z_TL: np.ndarray, Z_TT: np.ndarray, k_min: int,
                      n_leaves: int) -> list[dict]:
    """Return ``V_task = V(T_TL) ∪ V(T_TT)`` as a multiset.

    Each entry carries the standard ``tree_internal_nodes`` keys plus
    ``"source" ∈ {"TL", "TT"}``. Filter by
    ``k_min ≤ size ≤ n_leaves − 1`` (excludes leaf pairs and root).
    """
    out: list[dict] = []
    for Z, src in ((Z_TL, "TL"), (Z_TT, "TT")):
        for v in tree_internal_nodes(Z):
            if k_min <= v["size"] <= n_leaves - 1:
                v_out = {**v, "source": src}
                out.append(v_out)
    return out


def filter_by_size(nodes: list[dict], k_min: int) -> list[dict]:
    """Return only the internal nodes with ``size >= k_min``."""
    return [u for u in nodes if u["size"] >= k_min]


# ─────────────────────────── per-patient pass ───────────────────────────


def per_patient_rows(pat: str, band: str, k_min: int) -> list[dict]:
    """Compute one row per task internal node for (pat, band).

    Returns ``[]`` if any of the four phase trees fails to load (e.g.
    Pat_14 task_test before vendor cache rebuild).
    """
    Z_TL   = load_Z(pat, "task_learn", band)
    Z_TT   = load_Z(pat, "task_test",  band)
    Z_pre  = load_Z(pat, "rest_pre",   band)
    Z_post = load_Z(pat, "rest_post",  band)
    if any(z is None for z in (Z_TL, Z_TT, Z_pre, Z_post)):
        return []

    n_TL  = Z_TL.shape[0] + 1
    n_TT  = Z_TT.shape[0] + 1
    n_pre = Z_pre.shape[0] + 1
    n_post = Z_post.shape[0] + 1
    if not (n_TL == n_TT == n_pre == n_post):
        # Different leaf counts across phases — skip cell. Pat_10 task↔rest
        # is reconciled by PATIENT_CHANNEL_DROP at FC load time, so this
        # should never trigger in normal operation; guard anyway.
        return []
    n_leaves = n_TL

    V_pre  = filter_by_size(tree_internal_nodes(Z_pre),  k_min)
    V_post = filter_by_size(tree_internal_nodes(Z_post), k_min)
    V_task = pooled_task_nodes(Z_TL, Z_TT, k_min, n_leaves)

    rows: list[dict] = []
    for v in V_task:
        S = v["leaves"]
        J_pre  = best_match(S, V_pre)
        J_post = best_match(S, V_post)
        rows.append({
            "patient": pat,
            "band": band,
            "source": v["source"],
            "h_rel": v["h_rel"],
            "size": v["size"],
            "J_pre": J_pre,
            "J_post": J_post,
            "delta_J": J_post - J_pre,
        })
    return rows


# ─────────────────────────── cohort aggregation ───────────────────────────


def make_bin_edges(K: int, h_min: float) -> np.ndarray:
    """Log-spaced bin edges in ``(0, 1]`` of length K+1.

    Reuses ``h_log_grid`` (returns ``n`` log-spaced points from h_min to 1);
    we request ``K+1`` points → K bins.
    """
    return h_log_grid([1.0], n=K + 1, h_min_floor=h_min)


def assign_bin(h_rel: float, edges: np.ndarray) -> int:
    """Return the bin index ``k ∈ {0, …, K-1}`` of ``h_rel`` in ``edges``.

    Edges are half-open ``[edges[k], edges[k+1])`` for ``k < K-1`` and
    closed on the right at the last edge: ``[edges[K-1], edges[K]]``.
    Returns ``-1`` if ``h_rel`` lies below ``edges[0]``.
    """
    K = len(edges) - 1
    if h_rel < edges[0]:
        return -1
    # Right-closed last bin so h_rel = 1 lands in bin K-1.
    if h_rel >= edges[K]:
        return K - 1
    # np.searchsorted: index i s.t. edges[i-1] <= h_rel < edges[i].
    return int(np.searchsorted(edges, h_rel, side="right") - 1)


def aggregate_cohort(per_node: pd.DataFrame, j_min_sweep: list[float],
                     K: int, h_min: float) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    """Aggregate per-node rows into cohort fields binned along h_rel.

    Returns:
      - ``cohort_step``  — long: (band, bin, J_min, M_step, Pi, n_rows, n_pat)
      - ``cohort_smooth`` — long: (band, bin, M_smooth, n_rows, n_pat)
      - ``edges``        — bin edges array of length K+1
    """
    edges = make_bin_edges(K, h_min)
    df = per_node.copy()
    df["bin"] = df["h_rel"].map(lambda x: assign_bin(x, edges))
    df = df[df["bin"] >= 0].copy()  # drop rows below h_min

    bands = BRAIN_BANDS_NAMES
    step_rows: list[dict] = []
    smooth_rows: list[dict] = []
    for band in bands:
        for k in range(K):
            cell = df[(df["band"] == band) & (df["bin"] == k)]
            n_rows = len(cell)
            n_pat = cell["patient"].nunique()
            if n_rows == 0:
                M_smooth = np.nan
            else:
                M_smooth = float(np.mean(np.maximum(cell["delta_J"].to_numpy(), 0.0)))
            smooth_rows.append({
                "band": band, "bin": k,
                "h_lo": float(edges[k]), "h_hi": float(edges[k + 1]),
                "h_mid_log": float(np.exp(0.5 * (np.log(edges[k]) + np.log(edges[k + 1])))),
                "M_smooth": M_smooth,
                "n_rows": n_rows, "n_pat": n_pat,
            })
            for J_min in j_min_sweep:
                if n_rows == 0:
                    M_step = np.nan
                    Pi = np.nan
                    pat_with = 0
                else:
                    retained = ((cell["J_pre"] < J_min) & (cell["J_post"] >= J_min))
                    M_step = float(retained.mean())
                    pats_in = cell["patient"].unique()
                    pats_with_retained = {p for p in pats_in
                                          if retained[cell["patient"] == p].any()}
                    pat_with = len(pats_with_retained)
                    Pi = pat_with / max(n_pat, 1) if n_pat > 0 else np.nan
                step_rows.append({
                    "band": band, "bin": k,
                    "h_lo": float(edges[k]), "h_hi": float(edges[k + 1]),
                    "h_mid_log": float(np.exp(0.5 * (np.log(edges[k]) + np.log(edges[k + 1])))),
                    "J_min": J_min,
                    "M_step": M_step,
                    "Pi": Pi,
                    "pat_with_retained": pat_with,
                    "n_rows": n_rows, "n_pat": n_pat,
                })
    return pd.DataFrame(step_rows), pd.DataFrame(smooth_rows), edges


# ─────────────────────────── main ───────────────────────────


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="MRL — Module-Retention Landscape.")
    ap.add_argument("--patients", type=str, default=None,
                    help="Comma-sep patient IDs (default: PATIENTS_4PHASE).")
    ap.add_argument("--bands", type=str, default=None,
                    help="Comma-sep band names (default: BRAIN_BANDS_NAMES).")
    ap.add_argument("--j-min-sweep", type=str, default="0.85,0.9,0.95",
                    help="Comma-sep J_min values (default 0.85,0.9,0.95).")
    ap.add_argument("--k-min", type=int, default=3,
                    help="Subtree-size floor (default 3).")
    ap.add_argument("--n-bins", type=int, default=12,
                    help="Number of log-spaced h_rel bins for cohort viz "
                         "(default 12).")
    ap.add_argument("--h-min", type=float, default=0.05,
                    help="Lower edge of h_rel bin grid (default 0.05).")
    ap.add_argument("-v", "--verbose", action="store_true")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    patients = (args.patients.split(",") if args.patients
                else list(PATIENTS_4PHASE))
    bands = args.bands.split(",") if args.bands else list(BRAIN_BANDS_NAMES)
    j_min_sweep = [float(x) for x in args.j_min_sweep.split(",")]
    k_min = int(args.k_min)
    K = int(args.n_bins)
    h_min = float(args.h_min)
    v = bool(args.verbose)

    out_dir = REPORTS_ROOT / "imcoh_mrl"
    out_dir.mkdir(parents=True, exist_ok=True)

    if v:
        print(f"[mrl] patients = {patients}")
        print(f"[mrl] bands = {bands}")
        print(f"[mrl] J_min sweep = {j_min_sweep}")
        print(f"[mrl] k_min = {k_min}, K = {K}, h_min = {h_min}")
        print(f"[mrl] output dir = {out_dir}")

    rows: list[dict] = []
    skipped: list[tuple[str, str]] = []
    for pat in patients:
        for band in bands:
            cell_rows = per_patient_rows(pat, band, k_min)
            if not cell_rows:
                skipped.append((pat, band))
                if v:
                    print(f"[mrl] SKIP {pat:7s} {band:10s} — incomplete LRG cache")
                continue
            rows.extend(cell_rows)
            if v:
                print(f"[mrl]   {pat:7s} {band:10s}  V_task = {len(cell_rows):4d}")

    if not rows:
        raise SystemExit("[mrl] no rows produced — check LRG caches.")

    per_node = pd.DataFrame(rows)
    per_node_path = out_dir / "mrl_per_node.csv"
    per_node.to_csv(per_node_path, index=False)
    if v:
        print(f"[mrl] wrote {per_node_path}  ({len(per_node)} rows)")

    cohort_step, cohort_smooth, edges = aggregate_cohort(
        per_node, j_min_sweep, K, h_min,
    )
    cohort_step_path = out_dir / "mrl_cohort_step.csv"
    cohort_smooth_path = out_dir / "mrl_cohort_smooth.csv"
    cohort_step.to_csv(cohort_step_path, index=False)
    cohort_smooth.to_csv(cohort_smooth_path, index=False)
    np.savetxt(out_dir / "mrl_bin_edges.txt", edges, fmt="%.6f")
    if v:
        print(f"[mrl] wrote {cohort_step_path}    ({len(cohort_step)} rows)")
        print(f"[mrl] wrote {cohort_smooth_path}  ({len(cohort_smooth)} rows)")

    # Quick stdout summary at primary J_min = 0.9.
    primary = 0.9 if 0.9 in j_min_sweep else j_min_sweep[len(j_min_sweep) // 2]
    s = cohort_step[cohort_step["J_min"] == primary].copy()
    print()
    print(f"[mrl] === summary at J_min = {primary} ===")
    print(f"[mrl] patients with full cell coverage: "
          f"{len([(p, b) for p in patients for b in bands]) - len(skipped)} / "
          f"{len(patients) * len(bands)}")
    if skipped:
        print(f"[mrl] skipped (p, b) pairs: {skipped}")
    print()
    print(f"[mrl] M_step  per band (max over bins, mean over bins):")
    for band in bands:
        bs = s[s["band"] == band]
        if bs.empty or bs["M_step"].dropna().empty:
            print(f"[mrl]   {band:10s}  no data")
            continue
        m_max = float(np.nanmax(bs["M_step"].to_numpy()))
        m_mean = float(np.nanmean(bs["M_step"].to_numpy()))
        pi_max = float(np.nanmax(bs["Pi"].to_numpy()))
        cells_cohort = int((bs["Pi"] >= 7/9).fillna(False).sum())
        print(f"[mrl]   {band:10s}  max={m_max:.3f}  mean={m_mean:.3f}  "
              f"Π_max={pi_max:.3f}  cohort-wide cells (Π≥7/9)={cells_cohort}")

    print()
    print(f"[mrl] M_smooth  per band (max over bins, mean over bins):")
    for band in bands:
        bs = cohort_smooth[cohort_smooth["band"] == band]
        if bs.empty or bs["M_smooth"].dropna().empty:
            print(f"[mrl]   {band:10s}  no data")
            continue
        m_max = float(np.nanmax(bs["M_smooth"].to_numpy()))
        m_mean = float(np.nanmean(bs["M_smooth"].to_numpy()))
        print(f"[mrl]   {band:10s}  max={m_max:.3f}  mean={m_mean:.3f}")


if __name__ == "__main__":
    main()
