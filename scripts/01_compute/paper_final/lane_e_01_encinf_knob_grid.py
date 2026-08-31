#!/usr/bin/env python3
r"""Lane E / E1 -- the five-phase encoding-vs-inference arc, knob-integrated.

Re-derives the encoding / inference decomposition on the LOCKED substrate of
``.agents/preprint/locked/PIPELINE_CONTRACT.md``: ``imcoh_abs``, mst-union
backbone **integrated over the plateau fractions** ``f in {0.07, 0.10, 0.14,
0.20}`` rather than reported at one, the locked 16-point scale grid, and
matched-strength on the dense FC before sparsification.

Per (patient, band) this writes, for every plateau fraction:

  obs   (nF, nS, 4)      the four functionals on the real arc
  swap  (nF, nS, 4)      the same four with ``task_learn`` and ``task_test``
                         EXCHANGED -- the construction-preserving role null
  surr  (nF, R,  nS, 4)  matched-strength; the SAME shuffled dense matrix feeds
                         all four fractions, so fraction-to-fraction differences
                         are attributable to the backbone and not to surrogate
                         noise

The swap costs nothing: the cophenetic distances are already computed and the
library's ``cross_phase_functionals`` takes the role map as data, so exchanging
two roles is a different ``roles`` dict, not different code.

5-point critical preamble: see ``lane_e_00_preregistration.md``, frozen before
any number in this lane. In one line -- the claim is that the encoding and the
inference-specific reorganisation each persist into ``rest_post`` (C1) and that
the two are distinguishable (C2); the null for C1 is matched-strength, which
matches the arc's construction and therefore its construction baseline; the null
for C2 is the role swap, which is IDENTICAL in construction and differs only in
which task block is called the premises block; neither null reaches the
coherency estimator, session nonstationarity or lag-specificity, and the swap
null cannot separate cognitive role from block order or block duration.

Outputs: data/paper_final/lane_e_encinf/grid/cells/<patient>__<band>.npz
"""
from __future__ import annotations

import argparse
import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, PATIENTS_4PHASE
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import (
    CROSS_PHASE_ROLES,
    cophenetic_at_scale,
    cross_phase_functionals,
    laplacian_eig,
)
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.substrate import (
    CANONICAL,
    CANONICAL_PHASES,
    canonical_graph,
    canonical_scale_grid,
)

ROOT = setup_script_env()

R = int(os.environ.get("LANE_E_R", "200"))
WORKERS = int(os.environ.get("LANE_E_WORKERS", "8"))
OUT = Path(os.environ.get("LANE_E_OUT", ROOT / "data" / "paper_final" / "lane_e_encinf")) / "grid"
RESUME = os.environ.get("LANE_E_RESUME", "1") not in ("0", "false", "False")

FUNCS = ("T_test", "T_learn", "T_infspec", "T_infspec_pe")
SGRID = canonical_scale_grid()
FRACS = tuple(float(f) for f in CANONICAL.plateau_fracs)
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260831

#: The role map with ``encode`` and ``probe`` exchanged. Same five graphs, same
#: split-half baseline structure, same durations -- only the semantic assignment
#: of the two task blocks moves.
SWAP_ROLES = dict(CROSS_PHASE_ROLES,
                  encode=CROSS_PHASE_ROLES["probe"],
                  probe=CROSS_PHASE_ROLES["encode"])


# --------------------------------------------------------------------------- #
def _funcs_over_scales(Ws: dict, roles_list) -> list[np.ndarray]:
    """``[(nS, 4)]`` per role map, from one set of five graphs.

    The eigendecomposition and the cophenetic distances are computed once and
    shared across the role maps, so the role swap is free.
    """
    eig = {ph: laplacian_eig(Ws[ph]) for ph in CANONICAL_PHASES}
    out = [np.full((SGRID.size, len(FUNCS)), np.nan) for _ in roles_list]
    for j, s in enumerate(SGRID):
        try:
            D = {ph: cophenetic_at_scale(*eig[ph], s) for ph in CANONICAL_PHASES}
            if any((not np.all(np.isfinite(v))) or np.std(v) == 0 for v in D.values()):
                continue
        except Exception:                                          # noqa: BLE001
            continue
        for i, roles in enumerate(roles_list):
            try:
                v = cross_phase_functionals(D, roles)
            except Exception:                                      # noqa: BLE001
                continue
            out[i][j] = [v[k] for k in FUNCS]
    return out


def per_cell(job):
    idx, pat, band = job
    t0 = time.time()
    cell = OUT / "cells" / f"{pat}__{band}.npz"
    if cell.exists() and RESUME:
        return dict(patient=pat, band=band, reused=True, elapsed_s=0.0)
    try:
        dense = {ph: canonical_graph(pat, ph, band, dense=True)
                 for ph in CANONICAL_PHASES}
    except Exception as exc:                                       # noqa: BLE001
        return dict(patient=pat, band=band, error=f"{type(exc).__name__}: {exc}")
    N = dense["A"].shape[0]

    nF, nS, nK = len(FRACS), SGRID.size, len(FUNCS)
    obs = np.full((nF, nS, nK), np.nan)
    swap = np.full((nF, nS, nK), np.nan)
    for i, fr in enumerate(FRACS):
        Ws = {ph: select_backbone(dense[ph], CANONICAL.backbone, frac=fr)
              for ph in CANONICAL_PHASES}
        obs[i], swap[i] = _funcs_over_scales(Ws, [CROSS_PHASE_ROLES, SWAP_ROLES])

    rng = np.random.default_rng([BASE_SEED, idx])
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    surr = np.full((nF, R, nS, nK), np.nan)
    for r in range(R):
        sh = {ph: matched_strength_shuffle(dense[ph], n_swaps, rng, W_MAX)
              for ph in CANONICAL_PHASES}
        for i, fr in enumerate(FRACS):
            Ws = {ph: select_backbone(sh[ph], CANONICAL.backbone, frac=fr)
                  for ph in CANONICAL_PHASES}
            surr[i, r] = _funcs_over_scales(Ws, [CROSS_PHASE_ROLES])[0]

    cell.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cell, s=SGRID, fracs=np.array(FRACS), funcs=np.array(FUNCS),
                        obs=obs, swap=swap, surr=surr, N=np.array([N]))
    return dict(patient=pat, band=band, N=N, reused=False, elapsed_s=time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", default="")
    ap.add_argument("--bands", default="")
    ap.add_argument("--R", type=int, default=0)
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    global R, WORKERS
    R = a.R or R
    WORKERS = a.workers or WORKERS
    pats = a.patients.split(",") if a.patients else list(PATIENTS_4PHASE)
    bands = a.bands.split(",") if a.bands else list(BRAIN_BANDS_NAMES)
    (OUT / "cells").mkdir(parents=True, exist_ok=True)
    matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))   # jit warm

    jobs = [(i, p, b) for i, (b, p) in
            enumerate((b, p) for b in bands for p in pats)]
    print(f"[lane-e/E1] {len(jobs)} cells | 5 phases | {CANONICAL.label()} "
          f"| fracs={FRACS} | R={R} | {SGRID.size} scales | {WORKERS} workers "
          f"-> {OUT}", flush=True)
    t0 = time.time()
    rows = []
    with Pool(WORKERS) as pool:
        for i, res in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            rows.append(res)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {res.get('patient','?')}/{res.get('band','?'):11s}"
                  f" {el:6.0f}s ETA {el/i*(len(jobs)-i):6.0f}s"
                  + (f"  ERROR {res['error'][:70]}" if res.get("error") else ""),
                  flush=True)
    pd.DataFrame(rows).to_csv(OUT / "descriptive.csv", index=False)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="lane_e_encinf_knob_grid", substrate=CANONICAL.label(),
        transform=CANONICAL.transform, backbone=CANONICAL.backbone,
        fracs=list(FRACS), plateau=list(CANONICAL.plateau or ()),
        R=R, swap_factor=SWAP_FACTOR, s_grid=[float(x) for x in SGRID],
        funcs=list(FUNCS), phases=list(CANONICAL_PHASES), cohort=pats, bands=bands,
        base_seed=BASE_SEED,
        nulls="matched-strength on dense FC then sparsify (shared draw across fracs); "
              "learn/test role swap on the observed graphs",
    ), indent=2))
    print(f"\n[lane-e/E1] done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
