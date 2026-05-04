#!/usr/bin/env python3
"""Audit Step 14b — Trace-coherence diagnostic for MSPC.

User concern: red ("trace") leaves on the dendrogram look scattered and
might not actually co-cluster across `task` ↔ `rest_post`. MSPC is per-leaf
by construction (each leaf labelled trace iff its OWN cluster-mates are
similar in task and post), so it does NOT enforce that all trace leaves
belong to ONE persistent cluster. This script checks groupwise coherence:

  For each (patient, band) cell, take the K-averaged-dominant=trace leaves T.
  At each scale k, compute the fcluster ids restricted to T in task and post.
  Report the maximum-cluster purity (largest task-cluster ∩ T) / |T| and
  same for post. A trace claim is "groupwise" iff *both* purities are high
  at the same scale.

Outputs:
  data/audit/per_patient_hierarchy_mspc/trace_coherence_check.csv
  data/audit/per_patient_hierarchy_mspc/trace_coherence_check.md  (top-line)
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_08_per_patient_hierarchy import (  # noqa: E402
    BAND_ORDER, PHASE_ORDER, PATIENTS_4PHASE, _load_Z,
)
from audit_14_mspc import mspc_for_cell, _dominant_and_confidence  # noqa: E402


OUT_DIR = ROOT / "data" / "audit" / "per_patient_hierarchy_mspc"


def _restricted_purity(labels: np.ndarray, T: np.ndarray) -> tuple[float, int]:
    """For label vector and trace-leaf index array T, return:
        - max_cluster_purity = max_c |labels[T]==c| / |T|
        - n_clusters_containing_T = number of distinct labels in labels[T]
    """
    if T.size == 0:
        return 0.0, 0
    sub = labels[T]
    vals, counts = np.unique(sub, return_counts=True)
    return float(counts.max() / T.size), int(vals.size)


def coherence_for_cell(pat: str, band: str) -> dict | None:
    """Test whether trace-AT-k leaves form one coherent task↔post group.

    At each scale k where ≥3 leaves are MSPC-classified `trace`, ask: do
    those leaves cluster together in `task` AND in `rest_post`? Purity is
    measured only at *meaningful* scales (k ≥ 5 to avoid the trivial
    root-split where everything lives in one big cluster). The reported
    cell is the scale that maximises min(task_purity, post_purity); if no
    such k exists the row gets `n_trace=0`.
    """
    Zs = {ph: _load_Z(pat, ph, band) for ph in PHASE_ORDER}
    out = mspc_for_cell(Zs)
    if out is None:
        return None
    A_table, dom_at_k, k_values = out
    N = A_table.shape[0]
    has_tt = Zs["task_test"] is not None
    K_MIN = 5  # below this the partition is too coarse to mean anything
    # Find scale with the MOST trace leaves (the strongest evidence scale).
    # Tie-break by k (prefer the lowest k, i.e. coarser partition).
    best = None  # (n_trace, k, tp, pp, tnc, pnc)

    for ki, k in enumerate(k_values):
        if k < K_MIN:
            continue
        T = np.where(dom_at_k[:, ki] == 0)[0]
        if T.size < 3:
            continue
        lab_tl = fcluster(Zs["task_learn"], k, criterion="maxclust")
        lab_rp = fcluster(Zs["rest_post"], k, criterion="maxclust")
        if has_tt:
            lab_tt = fcluster(Zs["task_test"], k, criterion="maxclust")
            tp_tt, tnc_tt = _restricted_purity(lab_tt, T)
            tp_tl, tnc_tl = _restricted_purity(lab_tl, T)
            tp, tnc = (tp_tt, tnc_tt) if tp_tt >= tp_tl else (tp_tl, tnc_tl)
        else:
            tp, tnc = _restricted_purity(lab_tl, T)
        pp, pnc = _restricted_purity(lab_rp, T)
        if best is None or T.size > best[0]:
            best = (int(T.size), k, tp, pp, tnc, pnc)

    if best is None:
        return {
            "patient": pat, "band": band, "n_trace": 0,
            "best_k": -1, "task_purity": float("nan"),
            "post_purity": float("nan"),
            "task_n_clusters": 0, "post_n_clusters": 0,
            "min_purity": float("nan"),
        }
    n_trace_at_k, k, tp, pp, tnc, pnc = best
    return {
        "patient": pat, "band": band, "n_trace": int(n_trace_at_k),
        "best_k": int(k), "task_purity": float(tp), "post_purity": float(pp),
        "task_n_clusters": int(tnc), "post_n_clusters": int(pnc),
        "min_purity": float(min(tp, pp)),
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for pat in PATIENTS_4PHASE:
        print(f"{pat} ...", flush=True)
        for band in BAND_ORDER:
            r = coherence_for_cell(pat, band)
            if r is not None:
                rows.append(r)

    csv_path = OUT_DIR / "trace_coherence_check.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "patient", "band", "n_trace", "best_k",
            "task_purity", "post_purity", "min_purity",
            "task_n_clusters", "post_n_clusters",
        ])
        w.writeheader(); w.writerows(rows)

    # Top-line markdown
    n_cells = sum(1 for r in rows if r["n_trace"] > 0)
    n_groupwise = sum(1 for r in rows
                      if r["n_trace"] > 0 and r["min_purity"] >= 0.7)
    n_scattered = sum(1 for r in rows
                      if r["n_trace"] > 0 and r["min_purity"] < 0.4)
    md = [
        "# MSPC trace-coherence diagnostic\n",
        "Tests whether per-leaf trace labels cluster TOGETHER in `task` and "
        "`rest_post`, vs being scattered across many small persistent groups.\n",
        f"- Cells with ≥1 trace leaf: **{n_cells}**",
        f"- Cells where trace leaves form one coherent group "
        f"(min_purity ≥ 0.70): **{n_groupwise}** ({n_groupwise/n_cells*100:.1f}%)",
        f"- Cells where trace is scattered (min_purity < 0.40): "
        f"**{n_scattered}** ({n_scattered/n_cells*100:.1f}%)\n",
        "min_purity = min over phases of (largest cluster ∩ trace-leaves) / |trace-leaves|, "
        "maximised over scale k.\n",
        "See `trace_coherence_check.csv` for the full per-cell table.",
    ]
    (OUT_DIR / "trace_coherence_check.md").write_text("\n".join(md) + "\n")
    print(f"\n{n_cells} cells with trace leaves; "
          f"{n_groupwise} groupwise (≥0.70), {n_scattered} scattered (<0.40)")
    print(f"Wrote {csv_path.relative_to(ROOT)} and trace_coherence_check.md")


if __name__ == "__main__":
    main()
