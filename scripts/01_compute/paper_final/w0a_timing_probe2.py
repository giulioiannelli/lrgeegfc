#!/usr/bin/env python3
"""Timing probe 2: per-backbone cost + connectivity, so the A2 grid can be sized."""
from __future__ import annotations
import time, sys
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import load_phase
from lrg_eegfc.utils.fc.backbone import (
    select_backbone, backbone_density, mst_union_top_fraction,
    tmfg_backbone, pmfg_backbone, disparity_backbone, percolation_backbone,
    maximum_spanning_tree,
)

W = load_phase("Pat_02", "task_test", "beta")
N = W.shape[0]
print(f"N={N}", flush=True)


def ncomp(A):
    return connected_components(csr_matrix(A > 0), directed=False)[0]


def thresh_backbone(W, frac):
    A = np.asarray(W, float).copy()
    A = np.maximum(A, 0.0); np.fill_diagonal(A, 0.0); A = 0.5 * (A + A.T)
    r, c = np.triu_indices(N, 1)
    w = A[r, c]
    k = max(1, int(np.ceil(frac * w.size)))
    thr = np.partition(w, -k)[-k]
    return np.where(A >= thr, A, 0.0)


cases = [("mst", lambda: maximum_spanning_tree(W)),
         ("tmfg", lambda: tmfg_backbone(W)),
         ("pmfg", lambda: pmfg_backbone(W)),
         ("perc", lambda: percolation_backbone(W)[0])]
for f in (0.02, 0.05, 0.10, 0.20, 0.40, 1.00):
    cases.append((f"mst_union@{f}", lambda f=f: mst_union_top_fraction(W, f)))
for f in (0.05, 0.10, 0.20, 0.40):
    cases.append((f"thresh@{f}", lambda f=f: thresh_backbone(W, f)))
for a in (0.01, 0.05, 0.10, 0.20, 0.50):
    cases.append((f"disp@{a}", lambda a=a: disparity_backbone(W, alpha=a, ensure_connected=True)))
    cases.append((f"dispRAW@{a}", lambda a=a: disparity_backbone(W, alpha=a, ensure_connected=False)))

for name, fn in cases:
    t0 = time.time(); B = fn(); dt = (time.time() - t0) * 1000
    print(f"  {name:16s} dens={backbone_density(B):.4f} ncomp={ncomp(B):3d} "
          f"wfrac={B.sum()/W.sum():.3f} {dt:8.1f} ms", flush=True)
