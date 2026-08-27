#!/usr/bin/env python3
r"""W0-C: the five-phase arc across diffusion scales -- learning a structure vs applying it.

The cross-phase object has FIVE graphs, not four. ``task_learn`` presents the
ordered premises; ``task_test`` demands judgements about non-adjacent pairs that
only an inferred order can answer. Separating the two is the difference between
a cognition result and a generic task-effect, so every functional here is
computed at every diffusion scale rather than at one.

Per (patient, band, scale), on cophenetic distances of the five phases:

  e  = D_learn - D_A          e2 = D_learn - D_B        encoding
  f  = D_test  - D_learn                                inference-specific
                                                        (arm-invariant: never
                                                        references rest_pre)
  g  = D_test  - D_A          g2 = D_test  - D_B        whole task
  p  = D_post  - D_B          p2 = D_post  - D_A        persistence

  T_test       = 1/2[ rho(g, p)      + rho(g2, p2)      ]
  T_learn      = 1/2[ rho(e, p)      + rho(e2, p2)      ]
  T_infspec    = 1/2[ rho(f, p)      + rho(f,  p2)      ]
  T_infspec_pe = 1/2[ pr(f, p | e)   + pr(f,  p2 | e2)  ]

Two nulls, both computed here, because a CONDITIONAL statistic needs more than a
graph-level shuffle:

  matched-strength  the incumbent. Shuffle each dense FC strength-preservingly,
                    push all five through the identical pipeline.
  role permutation  the five OBSERVED graphs, with the five role labels
                    permuted. Under this null the estimator sees real cophenetic
                    geometries in an arbitrary arrangement, so any systematic
                    positivity it returns is the estimator's own bias rather
                    than consolidation. This exists because a partial
                    correlation has a documented failure of exactly that kind in
                    this project: a sham arc built inside pre-task rest returned
                    a significantly positive conditional trace (p = 0.007) that
                    exceeded the real value. It costs almost nothing: the
                    cophenetic distances are already computed, and permuting
                    roles only recombines them.

5-point critical preamble
1. Claim. (i) The whole-task trace, the encoding echo and the inference-specific
   component each persist into rest_post; (ii) encoding and inference persist at
   DIFFERENT diffusion scales -- a scale dissociation between learning a
   structure and applying it, which no single-scale method can state.
2. Null. Matched-strength on all five graphs for the trace claims; role
   permutation for the estimator-bias question. The cohort decision is the
   locked margin gate.
3. Strongest alternative each null must control for. For the trace: node
   strength alone. For the conditional functional: the estimator manufactures
   positive values from conditioning, independent of any consolidation -- the
   role-permutation null attacks precisely this, since it holds the graphs fixed
   and destroys only the assignment of roles. For the scale-dissociation claim:
   the two curves differ only in amplitude, not in where they live, so the
   dissociation is tested on scale POSITION (per-patient profile centroid),
   paired within patient, not on which curve is larger.
4. What the nulls cannot do. Neither reaches the connectivity estimator, the
   band split or session nonstationarity (lane W0-B). Role permutation treats
   the five phases as exchangeable, which they are not -- it bounds the
   estimator's structural bias, it does not test the science. n = 10.
5. Falsification. If T_infspec_pe is significantly positive under role
   permutation, the conditional functional is not reportable at all and the gate
   must refuse it. If the per-patient scale centroids of T_learn and
   T_infspec_pe do not separate, the scale dissociation is a negative and is
   reported as one.

Environment: same injection points as w0c_01 (W0C_FC_METHOD / W0C_BACKBONE /
W0C_FRAC / W0C_R / W0C_WORKERS), plus W0C_ENCINF_OUT and W0C_NPERM.

Outputs (data/paper_final/w0c_gate_tau/encinf/):
  cells/<patient>__<band>.npz   obs, matched-strength ensemble, role-permutation
                                ensemble, all (functional x scale)
  descriptive.csv, config.json
"""
from __future__ import annotations

import argparse
import json
import os
import time
from itertools import permutations
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import cophenetic_at_scale, laplacian_eig
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import phase_graphs

ROOT = setup_script_env()

FC_METHOD = os.environ.get("W0C_FC_METHOD", "imcoh_abs")
BACKBONE = os.environ.get("W0C_BACKBONE", "mst020")
FRAC = float(os.environ.get("W0C_FRAC", "0.20"))
R = int(os.environ.get("W0C_R", "200"))
N_PERM = int(os.environ.get("W0C_NPERM", "200"))
WORKERS = int(os.environ.get("W0C_WORKERS", "4"))
OUT = Path(os.environ.get("W0C_ENCINF_OUT",
                          ROOT / "data" / "paper_final" / "w0c_gate_tau" / "encinf"))
RESUME = os.environ.get("W0C_RESUME", "1") not in ("0", "false", "False")

PHASES5 = ("A", "B", "task_learn", "task_test", "rest_post")
FUNCS = ("T_test", "T_learn", "T_infspec", "T_infspec_pe")
SGRID = np.logspace(np.log10(0.05), np.log10(180.0), 28)
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260826


# --------------------------------------------------------------------------- #
def _ranks_centred(x):
    from scipy.stats import rankdata
    r = rankdata(x)
    return r - r.mean()


def _corr_matrix(vecs):
    """Spearman correlation matrix of a list of vectors, computed once."""
    Rk = np.array([_ranks_centred(v) for v in vecs])
    n = np.linalg.norm(Rk, axis=1, keepdims=True)
    Rk = np.divide(Rk, n, out=np.zeros_like(Rk), where=n > 0)
    return Rk @ Rk.T


def _partial_from(C, i, j, k):
    """Partial Spearman of ``i`` and ``j`` given ``k`` from a correlation matrix."""
    rij, rik, rjk = C[i, j], C[i, k], C[j, k]
    den = np.sqrt(max(0.0, (1 - rik ** 2) * (1 - rjk ** 2)))
    return (rij - rik * rjk) / den if den > 0 else np.nan


def functionals(D: dict) -> np.ndarray:
    """The four symmetric functionals from one dict of five cophenetic vectors.

    All six difference vectors go into a single Spearman correlation matrix, so
    each functional -- including the two partials -- is read off it rather than
    recomputing ranks. Returns an array in :data:`FUNCS` order.
    """
    e = D["task_learn"] - D["A"]
    e2 = D["task_learn"] - D["B"]
    f = D["task_test"] - D["task_learn"]
    g = D["task_test"] - D["A"]
    g2 = D["task_test"] - D["B"]
    p = D["rest_post"] - D["B"]
    p2 = D["rest_post"] - D["A"]
    C = _corr_matrix([e, e2, f, g, g2, p, p2])       # indices 0..6
    E, E2, F, G, G2, P, P2 = range(7)
    return np.array([
        0.5 * (C[G, P] + C[G2, P2]),
        0.5 * (C[E, P] + C[E2, P2]),
        0.5 * (C[F, P] + C[F, P2]),
        0.5 * (_partial_from(C, F, P, E) + _partial_from(C, F, P2, E2)),
    ])


def _backbone(W):
    return W if BACKBONE == "dense" else select_backbone(W, BACKBONE, frac=FRAC)


def coph_by_scale(Ws: dict) -> list[dict]:
    """Cophenetic distances of every phase at every scale: ``[{phase: vec}, ...]``."""
    eig = {ph: laplacian_eig(_backbone(Ws[ph])) for ph in PHASES5}
    out = []
    for s in SGRID:
        try:
            D = {ph: cophenetic_at_scale(*eig[ph], s) for ph in PHASES5}
            if any(not np.all(np.isfinite(v)) or np.std(v) == 0 for v in D.values()):
                D = None
        except Exception:                                          # noqa: BLE001
            D = None
        out.append(D)
    return out


def realization(Ws: dict) -> np.ndarray:
    """``(nS, nFunc)`` functionals for one set of five graphs."""
    vals = np.full((SGRID.size, len(FUNCS)), np.nan)
    for j, D in enumerate(coph_by_scale(Ws)):
        if D is not None:
            vals[j] = functionals(D)
    return vals


def role_permutation_null(Ds: list, n_perm: int, rng) -> np.ndarray:
    """``(n_perm, nS, nFunc)`` from permuting which observed graph plays which role.

    The cophenetic geometries are the real ones; only the mapping of graphs to
    the roles {A, B, task_learn, task_test, rest_post} is scrambled. The identity
    permutation is excluded. This isolates whatever the estimator returns from an
    arbitrary arrangement of real data -- the failure mode a partial correlation
    is prone to.
    """
    perms = [q for q in permutations(range(5)) if q != (0, 1, 2, 3, 4)]
    pick = rng.choice(len(perms), size=min(n_perm, len(perms)), replace=False)
    out = np.full((len(pick), SGRID.size, len(FUNCS)), np.nan)
    for i, ip in enumerate(pick):
        q = perms[ip]
        for j, D in enumerate(Ds):
            if D is None:
                continue
            Dp = {PHASES5[a]: D[PHASES5[q[a]]] for a in range(5)}
            out[i, j] = functionals(Dp)
    return out


# --------------------------------------------------------------------------- #
def per_cell(job):
    idx, pat, band = job
    t0 = time.time()
    try:
        Ws = phase_graphs(pat, band, fc_method=FC_METHOD, phases=PHASES5)
    except Exception as exc:                                       # noqa: BLE001
        return dict(patient=pat, band=band, error=str(exc))
    N = Ws["A"].shape[0]
    cell_path = OUT / "cells" / f"{pat}__{band}.npz"
    if cell_path.exists() and RESUME:
        return dict(patient=pat, band=band, N=N, reused=True,
                    elapsed_s=time.time() - t0)

    Ds = coph_by_scale(Ws)
    obs = np.full((SGRID.size, len(FUNCS)), np.nan)
    for j, D in enumerate(Ds):
        if D is not None:
            obs[j] = functionals(D)

    rng = np.random.default_rng(BASE_SEED + idx)
    perm = role_permutation_null(Ds, N_PERM, rng)

    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    surr = np.full((R, SGRID.size, len(FUNCS)), np.nan)
    for r in range(R):
        Wsh = {ph: matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX) for ph in PHASES5}
        try:
            surr[r] = realization(Wsh)
        except Exception:                                          # noqa: BLE001
            continue

    cell_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cell_path, s=SGRID, funcs=np.array(FUNCS), obs=obs,
                        surr=surr, perm=perm, N=np.array([N]))
    return dict(patient=pat, band=band, N=N, reused=False, elapsed_s=time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", type=str, default="")
    ap.add_argument("--bands", type=str, default="")
    ap.add_argument("--R", type=int, default=0)
    a = ap.parse_args()
    global R
    if a.R:
        R = a.R
    pats = a.patients.split(",") if a.patients else list(PATIENTS_4PHASE)
    bands = a.bands.split(",") if a.bands else list(BRAIN_BANDS_NAMES)
    OUT.mkdir(parents=True, exist_ok=True)
    matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))

    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in bands for p in pats)]
    print(f"[w0c-encinf] {len(jobs)} cells | 5 phases | fc={FC_METHOD} "
          f"backbone={BACKBONE}@{FRAC} | R={R} perm={N_PERM} | {SGRID.size} scales "
          f"| {WORKERS} workers -> {OUT}", flush=True)
    t0 = time.time()
    rows = []
    with Pool(WORKERS) as pool:
        for i, res in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if res:
                rows.append(res)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {res.get('patient','?')}/{res.get('band','?'):11s} "
                  f"{el:6.0f}s ETA {el/i*(len(jobs)-i):6.0f}s"
                  + (f"  ERROR {res['error'][:60]}" if res.get("error") else ""),
                  flush=True)
    pd.DataFrame(rows).to_csv(OUT / "descriptive.csv", index=False)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="w0c_encinf_scale_grid", fc_method=FC_METHOD, backbone=BACKBONE,
        frac=FRAC, R=R, n_perm=N_PERM, swap_factor=SWAP_FACTOR,
        s_grid=[float(x) for x in SGRID], funcs=list(FUNCS), phases=list(PHASES5),
        cohort=pats, bands=bands, base_seed=BASE_SEED,
        nulls="matched-strength (all 5 graphs) + role permutation (observed graphs)",
    ), indent=2))
    print(f"\n[w0c-encinf] done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
