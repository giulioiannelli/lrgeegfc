#!/usr/bin/env python3
r"""Lane E / E3 -- do encoding and inference persistence load on DIFFERENT contacts?

The band axis and the scale axis have both been tested and neither separates
encoding from inference on this substrate. The remaining axis named in the lane
brief is spatial: *which contact pairs carry each persistence*. This script
tests it, and it is the one separability question that does not go through the
conditional statistic at all.

CONSTRUCTION. A Spearman correlation is a sum of per-pair terms, so the trace
statistics decompose exactly over contact pairs. With centred ranks (``~``),

  c_enc = e~ * p~     the per-pair contribution to the ENCODING persistence
  c_inf = f~ * p~     the per-pair contribution to the INFERENCE-specific one

and ``sum_pairs c = rho`` in each case. Aggregating each to contacts with the
library's endpoint-incidence mean gives two node profiles, ``L_enc`` and
``L_inf``: how much each contact contributes to each persistence.

THE STATISTIC, and why it needs a surrogate rather than a zero.
``L_enc`` and ``L_inf`` share the factor ``p~``, so they correlate across nodes
even when nothing is going on -- a raw correlation between them is therefore
uninterpretable against zero, exactly as the conditional statistic is. What is
interpretable is whether they are LESS alike than a null that inherits the same
shared factor. So the reported quantity is

  S_dissoc = -Spearman_nodes(L_enc, L_inf)

margin-tested against matched-strength, which builds its five phases the same
way and therefore carries the same ``p~``-sharing. A positive margin means the
two persistences distribute over contacts more differently than the shared
construction alone produces.

A second, more directly anatomical pair of readouts uses magnitudes rather than
signed contributions, and leans on the fact that ``f = D_test - D_learn`` never
references ``rest_pre`` and so is structurally immune to split-half baseline
artifacts:

  rho_PE = Spearman_nodes(|p~| loading, |e~| loading)
  rho_PF = Spearman_nodes(|p~| loading, |f~| loading)
  delta  = rho_PF - rho_PE

margin-tested the same way: positive means what persists is spatially closer to
what the PROBE moved than to what the PREMISES moved.

5-point critical preamble (the lane pre-registration is
``lane_e_00_preregistration.md``; this is the E3 instance of it).
1. CLAIM: encoding persistence and inference-specific persistence load on
   different contacts.
2. NULL: matched-strength on the dense FC of all five phases, sparsified
   afterwards, one draw shared across plateau fractions -- the same null the
   trace gate uses, so it inherits the arc's construction including the shared
   ``p~`` factor that makes a zero reference meaningless here.
3. STRONGEST ALTERNATIVE: that any apparent spatial difference is the arithmetic
   of two difference vectors sharing a term, plus per-contact differences in
   estimator noise (contact strength, probe membership), rather than anything
   about cognition.
4. DOES THE NULL CONTROL FOR IT: for the shared-term arithmetic, yes, by
   construction -- the surrogate arc has the identical algebra. For per-contact
   noise, partly: matched-strength preserves each node's strength, which is the
   dominant driver of per-contact estimator noise, but not the probe geometry.
   It CANNOT reach the coherency estimator, session nonstationarity, or the fact
   that ``task_test`` is the longer and later block in all 10 patients -- a
   contact whose signal drifts across the session will look "probe-like" for
   reasons that have nothing to do with inference. NO spatial-proximity or
   node-label permutation null is used; this project forbids them on |ImCoh|.
5. FALSIFICATION: if the margin of ``S_dissoc`` and of ``delta`` are both
   indistinguishable from the matched-strength null across the scale axis, the
   spatial axis does not separate encoding from inference and E3 is negative on
   its last available axis. That is reported as the result.

Outputs: data/paper_final/lane_e_encinf/loading/cells/<patient>__<band>.npz
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
from scipy.stats import rankdata, spearmanr

from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import cophenetic_at_scale, laplacian_eig
from lrg_eegfc.utils.metrics.node_localization import node_incidence_mean
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.substrate import (
    CANONICAL,
    CANONICAL_PHASES,
    canonical_graph,
    canonical_scale_grid,
)

ROOT = setup_script_env()

R = int(os.environ.get("LANE_E_LOAD_R", "100"))
WORKERS = int(os.environ.get("LANE_E_WORKERS", "4"))
OUT = Path(os.environ.get("LANE_E_OUT",
                          ROOT / "data" / "paper_final" / "lane_e_encinf")) / "loading"
RESUME = os.environ.get("LANE_E_RESUME", "1") not in ("0", "false", "False")

STATS = ("S_dissoc", "rho_PE", "rho_PF", "delta_PF_PE")
SGRID = canonical_scale_grid()
FRACS = tuple(float(f) for f in CANONICAL.plateau_fracs)
SWAP_FACTOR, W_MAX = 20, 1.0
BASE_SEED = 20260903


def _cr(x):
    """Centred ranks -- Pearson on these is Spearman."""
    r = rankdata(x)
    return r - r.mean()


def loading_stats(D: dict, ni, nj, n_nodes) -> np.ndarray:
    """The four node-loading statistics from one set of five cophenetic vectors."""
    e = _cr(D["task_learn"] - D["A"])
    f = _cr(D["task_test"] - D["task_learn"])
    p = _cr(D["rest_post"] - D["B"])
    L_enc = node_incidence_mean(e * p, ni, nj, n_nodes)
    L_inf = node_incidence_mean(f * p, ni, nj, n_nodes)
    P = node_incidence_mean(np.abs(p), ni, nj, n_nodes)
    E = node_incidence_mean(np.abs(e), ni, nj, n_nodes)
    F = node_incidence_mean(np.abs(f), ni, nj, n_nodes)
    ok = np.isfinite(L_enc) & np.isfinite(L_inf)
    s_d = -spearmanr(L_enc[ok], L_inf[ok]).statistic if ok.sum() > 4 else np.nan
    r_pe = spearmanr(P[ok], E[ok]).statistic if ok.sum() > 4 else np.nan
    r_pf = spearmanr(P[ok], F[ok]).statistic if ok.sum() > 4 else np.nan
    return np.array([s_d, r_pe, r_pf, r_pf - r_pe], float)


def stats_over_scales(Ws: dict, ni, nj, n_nodes) -> np.ndarray:
    eig = {ph: laplacian_eig(Ws[ph]) for ph in CANONICAL_PHASES}
    out = np.full((SGRID.size, len(STATS)), np.nan)
    for j, s in enumerate(SGRID):
        try:
            D = {ph: cophenetic_at_scale(*eig[ph], s) for ph in CANONICAL_PHASES}
            if any((not np.all(np.isfinite(v))) or np.std(v) == 0 for v in D.values()):
                continue
            out[j] = loading_stats(D, ni, nj, n_nodes)
        except Exception:                                          # noqa: BLE001
            continue
    return out


def per_cell(job):
    idx, pat, band = job
    t0 = time.time()
    cell = OUT / "cells" / f"{pat}__{band}.npz"
    if cell.exists() and RESUME:
        return dict(patient=pat, band=band, reused=True)
    try:
        dense = {ph: canonical_graph(pat, ph, band, dense=True)
                 for ph in CANONICAL_PHASES}
    except Exception as exc:                                       # noqa: BLE001
        return dict(patient=pat, band=band, error=f"{type(exc).__name__}: {exc}")
    N = dense["A"].shape[0]
    ni, nj = np.triu_indices(N, 1)

    nF, nS, nK = len(FRACS), SGRID.size, len(STATS)
    obs = np.full((nF, nS, nK), np.nan)
    for i, fr in enumerate(FRACS):
        Ws = {ph: select_backbone(dense[ph], CANONICAL.backbone, frac=fr)
              for ph in CANONICAL_PHASES}
        obs[i] = stats_over_scales(Ws, ni, nj, N)

    rng = np.random.default_rng([BASE_SEED, idx])
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    surr = np.full((nF, R, nS, nK), np.nan)
    for r in range(R):
        sh = {ph: matched_strength_shuffle(dense[ph], n_swaps, rng, W_MAX)
              for ph in CANONICAL_PHASES}
        for i, fr in enumerate(FRACS):
            Ws = {ph: select_backbone(sh[ph], CANONICAL.backbone, frac=fr)
                  for ph in CANONICAL_PHASES}
            surr[i, r] = stats_over_scales(Ws, ni, nj, N)

    cell.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cell, s=SGRID, fracs=np.array(FRACS),
                        stats=np.array(STATS), obs=obs, surr=surr, N=np.array([N]))
    return dict(patient=pat, band=band, N=N, reused=False, elapsed_s=time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", default="")
    ap.add_argument("--bands", default="theta,alpha,beta")
    ap.add_argument("--R", type=int, default=0)
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    global R, WORKERS
    R = a.R or R
    WORKERS = a.workers or WORKERS
    pats = a.patients.split(",") if a.patients else list(PATIENTS_4PHASE)
    bands = a.bands.split(",")
    (OUT / "cells").mkdir(parents=True, exist_ok=True)
    matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))

    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in bands for p in pats)]
    print(f"[lane-e/E3-loading] {len(jobs)} cells | {CANONICAL.label()} "
          f"fracs={FRACS} | R={R} | {SGRID.size} scales | {WORKERS} workers -> {OUT}",
          flush=True)
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
        deliverable="lane_e_loading_dissociation", substrate=CANONICAL.label(),
        fracs=list(FRACS), R=R, stats=list(STATS), bands=bands, cohort=pats,
        s_grid=[float(x) for x in SGRID], base_seed=BASE_SEED,
        null="matched-strength on dense FC then sparsify, shared draw across fracs",
    ), indent=2))
    print(f"\n[lane-e/E3-loading] done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
