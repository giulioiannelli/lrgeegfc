#!/usr/bin/env python3
"""Trace-modules — module identification at significant `(band, k)` cells.

For each `(band, (k_lo, k_hi))` cell flagged sig in the headline figure
(default: δ k=23–31, α k=2–4 from the n=9 cluster-perm pass), enumerate
candidate task subtrees `V_task = V(T_TL) ∪ V(T_TT)` per patient and emit
one CSV row per T-regime subtree (`¬present_pre ∧ present_post` at
`J_min = 0.9`, matches MRL).

Output:
  data/audit/trace_modules/trace_subtrees_n10_imcoh_abs.csv
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, to_tree

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, SEEG_DATAPATH
from lrg_eegfc.utils.metrics.tree import (
    dmax_from_Z, jaccard_leafsets, tree_internal_nodes,
)
from lrg_eegfc.workflow.lrg import load_lrg_result


COHORT_N10 = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
              "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
DEFAULT_K_BINS = {
    "delta": [(23, 31)],
    "alpha": [(2, 4)],
}
DEFAULT_J_MIN = 0.9
DEFAULT_K_MIN = 3


def _load_Z(pat: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, "imcoh_abs", IMCOH_LRG_CACHE)
    except Exception:
        return None
    return np.asarray(r.linkage_matrix) if r is not None else None


def _load_channel_labels(pat: str) -> list[str]:
    p = Path(SEEG_DATAPATH) / pat / "channel_labels.csv"
    if not p.exists():
        return []
    df = pd.read_csv(p)
    return df.iloc[:, 0].astype(str).tolist()


def _internal_nodes_at_k(Z: np.ndarray, k: int, k_min: int) -> list[dict]:
    """Internal nodes whose leafset corresponds to a cluster at scale `k`.

    We use scipy's `fcluster(maxclust)` to get the partition at `k` and then
    pull the leaf-sets of those clusters as candidate subtrees. Sizes < k_min
    are filtered.

    For trees where internal nodes don't naturally correspond to fcluster
    partitions, we still want subtrees at the right scale: this returns each
    cluster's leafset with its `h_rel` set to the merge height inside that
    cluster's subtree (taken as the deepest internal node fully contained
    in the cluster's leaf list).
    """
    labels = fcluster(Z, k, criterion="maxclust")
    leafsets: dict[int, frozenset[int]] = {}
    for leaf_idx, lab in enumerate(labels):
        leafsets.setdefault(int(lab), set()).add(int(leaf_idx))
    leafsets = {k_: frozenset(v) for k_, v in leafsets.items()}

    # Look up each cluster's internal-node merge height by finding the node
    # whose leaves match.
    nodes_all = tree_internal_nodes(Z)
    leafset_to_h: dict[frozenset[int], tuple[float, float]] = {}
    for nd in nodes_all:
        leafset_to_h[nd["leaves"]] = (nd["h"], nd["h_rel"])

    out = []
    for cl_id, leafset in leafsets.items():
        if len(leafset) < k_min:
            continue
        h_pair = leafset_to_h.get(leafset)
        if h_pair is None:
            # Singleton clusters or odd partitions; assign h=0
            h_h_rel = (0.0, 0.0)
        else:
            h_h_rel = h_pair
        out.append({
            "leaves": leafset,
            "size": len(leafset),
            "h": h_h_rel[0],
            "h_rel": h_h_rel[1],
            "cluster_id_at_k": cl_id,
        })
    return out


def _Vpost_candidates(Z: np.ndarray, k_min: int) -> list[frozenset[int]]:
    return [nd["leaves"] for nd in tree_internal_nodes(Z) if nd["size"] >= k_min]


def _jstar(S: frozenset[int], V: list[frozenset[int]]) -> float:
    if not V:
        return 0.0
    return max(jaccard_leafsets(S, leaves) for leaves in V)


def _classify(present_pre: bool, present_post: bool) -> str:
    """task is implicit (candidate drawn from task tree)."""
    if present_pre and present_post:
        return "P"
    if not present_pre and present_post:
        return "T"
    if not present_pre and not present_post:
        return "R"
    return "RA"  # present_pre and not present_post


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bands", default=",".join(DEFAULT_K_BINS.keys()),
                    help="comma-separated bands (default: delta,alpha)")
    ap.add_argument("--k-bins", default=None,
                    help="overrides DEFAULT_K_BINS; format 'lo-hi,lo-hi' "
                         "(applied identically to every band in --bands)")
    ap.add_argument("--j-min", type=float, default=DEFAULT_J_MIN)
    ap.add_argument("--k-min", type=int, default=DEFAULT_K_MIN)
    ap.add_argument("--j-min-sweep", default=None,
                    help="comma-separated extra J_min values; emits one CSV "
                         "per value (e.g. '0.7,0.85,0.95'). The primary --j-min "
                         "still drives the default output filename.")
    ap.add_argument("--output", default="data/audit/trace_modules/trace_subtrees_n10_imcoh_abs.csv")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    bands = args.bands.split(",")
    if args.k_bins:
        kbin_pairs = [tuple(int(x) for x in s.split("-"))
                      for s in args.k_bins.split(",")]
        k_bins_per_band = {b: kbin_pairs for b in bands}
    else:
        k_bins_per_band = {b: DEFAULT_K_BINS.get(b, [(2, 49)]) for b in bands}

    # Build the J_min sweep list; primary first, then extras.
    j_min_values = [args.j_min]
    if args.j_min_sweep:
        for s in args.j_min_sweep.split(","):
            v = float(s)
            if v not in j_min_values:
                j_min_values.append(v)

    # rows[(j_min)] -> list[dict]
    rows_per_jmin: dict[float, list[dict]] = {jm: [] for jm in j_min_values}
    # Always also dump every (J_pre, J_post) candidate row so downstream can
    # inspect the full distribution (matches the MRL "long format" convention).
    cand_rows: list[dict] = []

    for band in bands:
        for (k_lo, k_hi) in k_bins_per_band[band]:
            k_rep = int(round(math.exp((math.log(k_lo) + math.log(k_hi)) / 2.0)))
            if args.verbose:
                print(f"== {band} ({k_lo}–{k_hi}, rep k={k_rep}) ==")
            for pat in COHORT_N10:
                Z_pre  = _load_Z(pat, "rest_pre",  band)
                Z_TL   = _load_Z(pat, "task_learn", band)
                Z_TT   = _load_Z(pat, "task_test",  band)
                Z_post = _load_Z(pat, "rest_post", band)
                if any(Z is None for Z in (Z_pre, Z_TL, Z_TT, Z_post)):
                    if args.verbose:
                        print(f"  SKIP {pat} (missing Z)")
                    continue
                ch = _load_channel_labels(pat)
                V_pre  = _Vpost_candidates(Z_pre,  args.k_min)
                V_post = _Vpost_candidates(Z_post, args.k_min)
                V_task = []
                for src, Z_src in (("TL", Z_TL), ("TT", Z_TT)):
                    for nd in _internal_nodes_at_k(Z_src, k_rep, args.k_min):
                        nd["source"] = src
                        V_task.append(nd)
                emitted_per_jm = {jm: 0 for jm in j_min_values}
                for nd in V_task:
                    S = nd["leaves"]
                    J_pre  = _jstar(S, V_pre)
                    J_post = _jstar(S, V_post)
                    leaf_idx_sorted = sorted(S)
                    leaves_named = (
                        ";".join(ch[i] for i in leaf_idx_sorted if i < len(ch))
                        if ch else ""
                    )
                    cand_rows.append({
                        "patient": pat, "band": band,
                        "k_rep": k_rep, "k_lo": k_lo, "k_hi": k_hi,
                        "source": nd["source"],
                        "h_rel": float(nd["h_rel"]),
                        "size": int(nd["size"]),
                        "J_pre": float(J_pre),
                        "J_post": float(J_post),
                        "leaf_indices": ";".join(str(i) for i in leaf_idx_sorted),
                        "leaves_named": leaves_named,
                    })
                    for jm in j_min_values:
                        regime = _classify(J_pre >= jm, J_post >= jm)
                        if regime != "T":
                            continue
                        rows_per_jmin[jm].append({
                            "patient": pat, "band": band,
                            "k_rep": k_rep, "k_lo": k_lo, "k_hi": k_hi,
                            "j_min": jm,
                            "source": nd["source"],
                            "h_rel": float(nd["h_rel"]),
                            "size": int(nd["size"]),
                            "J_pre": float(J_pre),
                            "J_post": float(J_post),
                            "regime": regime,
                            "leaf_indices": ";".join(str(i) for i in leaf_idx_sorted),
                            "leaves_named": leaves_named,
                            "pointer_to_MRL_node_id": "",
                            "pointer_to_CBR_module_id": "",
                        })
                        emitted_per_jm[jm] += 1
                if args.verbose:
                    summary = " ".join(f"j={jm}:{c}" for jm, c in emitted_per_jm.items())
                    print(f"  {pat:8s} {band:10s} k={k_rep}: {len(V_task)} candidates → {summary} T-modules")

    out_path = Path(args.output)
    if not out_path.is_absolute():
        out_path = ROOT / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_primary = pd.DataFrame(rows_per_jmin[args.j_min])
    df_primary.to_csv(out_path, index=False)
    print(f"\nprimary J_min={args.j_min}: {len(df_primary)} T-rows → {out_path}")

    # Extra J_min files
    for jm in j_min_values[1:]:
        suffix = f"_jmin{int(jm*100):02d}"
        extra_path = out_path.with_name(out_path.stem + suffix + out_path.suffix)
        df_jm = pd.DataFrame(rows_per_jmin[jm])
        df_jm.to_csv(extra_path, index=False)
        print(f"sweep J_min={jm}: {len(df_jm)} T-rows → {extra_path}")

    # Always emit the candidate-distribution dump for downstream/diagnostic.
    cand_path = out_path.with_name(out_path.stem + "_candidates" + out_path.suffix)
    pd.DataFrame(cand_rows).to_csv(cand_path, index=False)
    print(f"all candidates: {len(cand_rows)} rows → {cand_path}")

    if not df_primary.empty:
        print("\nT-modules per (band, k_rep, patient) at primary J_min:")
        print(df_primary.groupby(["band", "k_rep", "patient"]).size().unstack(fill_value=0))


if __name__ == "__main__":
    main()
