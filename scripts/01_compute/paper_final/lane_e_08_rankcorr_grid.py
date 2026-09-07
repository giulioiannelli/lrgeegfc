#!/usr/bin/env python3
r"""Lane E / E5-E8 -- the real five-phase arc, stored as a rank-correlation tensor.

Same substrate, same seeds, same scale grid and same null as
``lane_e_01_encinf_knob_grid.py``. One thing changes, and it changes what the
lane can be asked afterwards: instead of four collapsed numbers per (fraction,
scale, realization) this writes the ``(8, 8)`` Spearman matrix of the eight
reorganisation vectors ``e, e2, f, g, g2, p, p2, d``.

That matrix is the sufficient statistic for the entire cross-phase family. The
four canonical functionals, the encode/probe role swap, the drift-controlled
variants and any contrast between them are exact functions of it -- verified to
3e-17 by ``lane_e_verify_rankcorr.py``. The practical consequence is that the
next time the estimator is questioned, the answer costs a CSV read rather than
three hours of eigendecomposition. Given that this lane has now revised its
estimator twice in response to what the nulls showed, that is the difference
between being able to keep questioning it and not.

The new eighth vector is the reason for the re-run. ``d = D_B - D_A`` is the
TASK-FREE drift: the two baselines are contiguous halves of one pre-task
recording, so ``d`` is the reorganisation the brain performs over a comparable
timespan with no task in it. The ordered sham showed that a five-phase arc laid
out in time inherits a positive cross-phase signal from adjacency alone, and the
symmetric A/B average cannot remove it (both arms pick up drift with the same
sign). ``d`` is the covariate that reaches that alternative, and it was always
implied by the arc -- it just was not being stored.

5-point critical preamble: ``lane_e_07_preregistration_addendum.md``, frozen
before any number from these estimators existed.

Outputs: data/paper_final/lane_e_encinf/rankcorr_grid/cells/<patient>__<band>.npz
  Robs  (nF, nS, 8, 8)      real arc
  Rsurr (nF, R, nS, 8, 8)   matched-strength on dense FC, one draw shared across
                            fractions so fraction-to-fraction differences are
                            attributable to the backbone, not to surrogate noise
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
    CROSS_PHASE_VECTOR_NAMES,
    cross_phase_rank_corr_over_scales,
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
WORKERS = int(os.environ.get("LANE_E_WORKERS", "6"))
OUT = Path(os.environ.get("LANE_E_OUT",
                          ROOT / "data" / "paper_final" / "lane_e_encinf")) / "rankcorr_grid"
RESUME = os.environ.get("LANE_E_RESUME", "1") not in ("0", "false", "False")

SGRID = canonical_scale_grid()
FRACS = tuple(float(f) for f in CANONICAL.plateau_fracs)
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260831           # identical to E1: same surrogate stream per cell


def _tensor(Ws: dict) -> np.ndarray:
    eig = {ph: laplacian_eig(Ws[ph]) for ph in CANONICAL_PHASES}
    return cross_phase_rank_corr_over_scales(eig, SGRID)


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
    nF, nS, nV = len(FRACS), SGRID.size, len(CROSS_PHASE_VECTOR_NAMES)

    Robs = np.full((nF, nS, nV, nV), np.nan, np.float32)
    for i, fr in enumerate(FRACS):
        Robs[i] = _tensor({ph: select_backbone(dense[ph], CANONICAL.backbone, frac=fr)
                           for ph in CANONICAL_PHASES})

    rng = np.random.default_rng([BASE_SEED, idx])
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    Rsurr = np.full((nF, R, nS, nV, nV), np.nan, np.float32)
    for r in range(R):
        sh = {ph: matched_strength_shuffle(dense[ph], n_swaps, rng, W_MAX)
              for ph in CANONICAL_PHASES}
        for i, fr in enumerate(FRACS):
            Rsurr[i, r] = _tensor(
                {ph: select_backbone(sh[ph], CANONICAL.backbone, frac=fr)
                 for ph in CANONICAL_PHASES})

    cell.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cell, s=SGRID, fracs=np.array(FRACS),
                        vectors=np.array(CROSS_PHASE_VECTOR_NAMES),
                        Robs=Robs, Rsurr=Rsurr, N=np.array([N]))
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
    print(f"[lane-e/E5-8] {len(jobs)} cells | rank-corr tensor "
          f"{len(CROSS_PHASE_VECTOR_NAMES)}x{len(CROSS_PHASE_VECTOR_NAMES)} | "
          f"{CANONICAL.label()} fracs={FRACS} | R={R} | {SGRID.size} scales | "
          f"{WORKERS} workers -> {OUT}", flush=True)
    t0, rows = time.time(), []
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
        deliverable="lane_e_rankcorr_grid", substrate=CANONICAL.label(),
        transform=CANONICAL.transform, backbone=CANONICAL.backbone,
        fracs=list(FRACS), plateau=list(CANONICAL.plateau or ()), R=R,
        swap_factor=SWAP_FACTOR, s_grid=[float(x) for x in SGRID],
        vectors=list(CROSS_PHASE_VECTOR_NAMES), phases=list(CANONICAL_PHASES),
        cohort=pats, bands=bands, base_seed=BASE_SEED,
        nulls="matched-strength on dense FC then sparsify (shared draw across fracs); "
              "role swap and drift control are exact post-hoc derivations",
    ), indent=2))
    print(f"\n[lane-e/E5-8] done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
