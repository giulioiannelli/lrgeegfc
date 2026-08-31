#!/usr/bin/env python3
r"""Lane E / E2a -- does the GATE calibrate on a no-signal five-phase arc?

W0-B established that all four cross-phase functionals return significantly
positive values on a sham arc carved out of a single resting recording -- 16 of
16 scales for alpha and beta even with temporal order (hence drift) destroyed,
and for beta ``T_infspec_pe`` the sham value (+0.119) EXCEEDED the real one
(+0.090). That was measured **against zero**. The project's gate has never been
"T > 0": it is "observed > that patient's own surrogate median". Whether the
gate calibrates is therefore an open and separate question, and it is the one
that decides whether any encoding/inference p-value is reportable.

This script answers it by construction: it builds the same no-signal arcs and
pushes them through the ENTIRE gate, matched-strength surrogates included. If
the margin on a no-signal arc sits at zero and the measured false-positive rate
is nominal, the gate is calibrated for these functionals and W0-C's cells stand;
if the margin is positive, the gate is broken for them and no p-value from it --
including the conditional one -- may be quoted.

The learn/test ROLE-SWAP test (lane E's construction-preserving null for the
separability claim) is calibrated on exactly the same arcs, because the sham's
``task_learn`` and ``task_test`` pseudo-phases differ in duration and position
but not in cognitive content -- which is precisely the confound the swap test
must be shown not to fire on.

CONSTRUCTIONS (five pseudo-phases from ONE resting recording; no task, hence no
consolidation, can exist inside any of them):

  sham_ordered    five contiguous windows in true temporal order, sized in
                  proportion to the real phase durations. Retains drift.
  sham_shuffled   the same five window sizes assembled from randomly reassigned
                  30 s blocks. Destroys temporal order, hence drift.

Run on ``rest_pre`` and, independently, on ``rest_post``, so a quirk of one
recording cannot masquerade as a property of the gate.

One deviation from W0-B's calibration, and it is a fix: the real pipeline
estimates the two ``rest_pre`` halves at ``nperseg // 2`` and the other three
phases at the full ``nperseg``. W0-B's sham used the full ``nperseg`` for all
five and flagged the inconsistency. Here the sham mirrors the real pipeline
exactly -- the A/B pseudo-phases are estimated at ``nperseg // 2`` -- so the
calibration is of the gate as actually used.

5-point critical preamble: ``lane_e_00_preregistration.md`` (frozen). In one
line -- the claim under test is that the MARGIN, not the raw value, of these
functionals is zero on no-signal input; the null is the sham arc itself; the
strongest alternative is that the construction offset survives the margin
subtraction because the matched-strength surrogate does not reproduce it; the
ordered/shuffled contrast separates drift from construction; and if the sham
margin is significantly positive the lane's own headline is that the gate does
not calibrate.

Outputs: data/paper_final/lane_e_encinf/sham/cells/<pat>__<band>__<source>.npz
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

from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    DEFAULT_SAMPLE_RATE,
    FS_OVERRIDES,
    PATIENTS_4PHASE,
    nperseg_for_fs,
)
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import (
    CROSS_PHASE_ROLES,
    cophenetic_at_scale,
    cross_phase_functionals,
    laplacian_eig,
)
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.utils.surrogate.timeseries_nulls import (
    coherency_from_csd,
    csd_from_segment_subset,
    imcoh_abs_from_coherency,
    segment_fft,
)
from lrg_eegfc.workflow.substrate import (
    CANONICAL,
    CANONICAL_PHASES,
    canonical_scale_grid,
)

ROOT = setup_script_env()

R = int(os.environ.get("LANE_E_SHAM_R", "100"))
N_SHUF = int(os.environ.get("LANE_E_SHAM_NSHUF", "3"))
WORKERS = int(os.environ.get("LANE_E_WORKERS", "8"))
OUT = Path(os.environ.get("LANE_E_OUT",
                          ROOT / "data" / "paper_final" / "lane_e_encinf")) / "sham"
RESUME = os.environ.get("LANE_E_RESUME", "1") not in ("0", "false", "False")

FUNCS = ("T_test", "T_learn", "T_infspec", "T_infspec_pe")
SGRID = canonical_scale_grid()
FRACS = tuple(float(f) for f in CANONICAL.plateau_fracs)
SWAP_FACTOR, W_MAX = 20, 1.0
BLOCK_SECONDS = 30.0
BASE_SEED = 20260831
SWAP_ROLES = dict(CROSS_PHASE_ROLES,
                  encode=CROSS_PHASE_ROLES["probe"],
                  probe=CROSS_PHASE_ROLES["encode"])
#: pseudo-phases estimated at ``nperseg // 2``, mirroring the real split-half arms
HALF_NPERSEG_PHASES = ("A", "B")

_DP = ROOT / "data" / "paper_final" / "w0b_nulls" / "durations.json"
DURS = json.loads(_DP.read_text()) if _DP.exists() else None


# --------------------------------------------------------------------------- #
def _funcs_over_scales(Ws: dict, roles_list) -> list[np.ndarray]:
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
                out[i][j] = [cross_phase_functionals(D, roles)[k] for k in FUNCS]
            except Exception:                                      # noqa: BLE001
                continue
    return out


def _arc_scores(dense: dict, rng, n_surr: int):
    """``(obs, swap, surr)`` for one five-phase set of DENSE matrices.

    Identical code path to the real arc: sparsify at every plateau fraction, and
    draw the matched-strength surrogate on the dense matrix before sparsifying,
    with one draw shared across the fractions.
    """
    nF, nS, nK = len(FRACS), SGRID.size, len(FUNCS)
    obs = np.full((nF, nS, nK), np.nan)
    swap = np.full((nF, nS, nK), np.nan)
    for i, fr in enumerate(FRACS):
        Ws = {ph: select_backbone(dense[ph], CANONICAL.backbone, frac=fr)
              for ph in CANONICAL_PHASES}
        obs[i], swap[i] = _funcs_over_scales(Ws, [CROSS_PHASE_ROLES, SWAP_ROLES])
    N = dense["A"].shape[0]
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    surr = np.full((nF, n_surr, nS, nK), np.nan)
    for r in range(n_surr):
        sh = {ph: matched_strength_shuffle(dense[ph], n_swaps, rng, W_MAX)
              for ph in CANONICAL_PHASES}
        for i, fr in enumerate(FRACS):
            Ws = {ph: select_backbone(sh[ph], CANONICAL.backbone, frac=fr)
                  for ph in CANONICAL_PHASES}
            surr[i, r] = _funcs_over_scales(Ws, [CROSS_PHASE_ROLES])[0]
    return obs, swap, surr


def _block_of_segment(n_seg: int, nperseg: int, blk: int) -> np.ndarray:
    """Block index of each Welch segment; ``-1`` for segments straddling a block."""
    seg_start = np.arange(n_seg) * (nperseg // 2)
    b0 = seg_start // blk
    b1 = (seg_start + nperseg - 1) // blk
    return np.where(b0 == b1, b0, -1).astype(np.int64)


def sham_dense(pat: str, band: str, source: str, mode: str, rng):
    """Five DENSE sham FC matrices carved out of ONE resting recording.

    ``mode`` is ``"identity"`` (contiguous windows in true order -- drift
    retained) or ``"free"`` (blocks randomly reassigned -- drift destroyed).
    Window sizes are proportional to the real phase durations, so the sham
    reproduces the real arc's duration profile and spectral degrees of freedom.
    """
    fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
    nps = nperseg_for_fs(fs)
    nps_h = max(256, nps // 2)
    bnd = BRAIN_BANDS[band]

    X = np.asarray(load_timeseries(pat, source, SEEG_DATAPATH), float)
    if X.shape[0] > X.shape[1]:
        X = X.T
    grids = {}
    for tag, nperseg in (("full", nps), ("half", nps_h)):
        _, F, sc = segment_fft(X, fs, nperseg, band=bnd)
        grids[tag] = (F, sc, _block_of_segment(F.shape[1], nperseg, int(round(BLOCK_SECONDS * fs))))
    del X

    d = DURS[pat]
    prop = np.array([d["rest_pre"] / 2, d["rest_pre"] / 2, d["task_learn"],
                     d["task_test"], d["rest_post"]], float)
    prop /= prop.sum()
    nb = int(min(g[2].max() for g in grids.values())) + 1
    sizes = np.maximum(2, np.floor(prop * nb).astype(int))
    while sizes.sum() > nb:
        sizes[int(np.argmax(sizes))] -= 1
    while sizes.sum() < nb:
        sizes[int(np.argmin(sizes))] += 1

    order = np.arange(nb) if mode == "identity" else rng.permutation(nb)
    dense, k = {}, 0
    for ph, n in zip(CANONICAL_PHASES, sizes):
        blocks = order[k:k + n]
        k += n
        F, sc, bos = grids["half" if ph in HALF_NPERSEG_PHASES else "full"]
        idx = np.sort(np.concatenate([np.flatnonzero(bos == b) for b in blocks]))
        if idx.size < 8:
            return None, sizes.tolist(), nb
        W = imcoh_abs_from_coherency(
            coherency_from_csd(csd_from_segment_subset(F, idx, sc)))
        W = np.asarray(W, float)
        np.fill_diagonal(W, 0.0)
        dense[ph] = 0.5 * (np.clip(W, 0.0, 1.0) + np.clip(W, 0.0, 1.0).T)
    return dense, sizes.tolist(), nb


def per_cell(job):
    idx, pat, band, source = job
    t0 = time.time()
    cell = OUT / "cells" / f"{pat}__{band}__{source}.npz"
    if cell.exists() and RESUME:
        return dict(patient=pat, band=band, source=source, reused=True)
    rng = np.random.default_rng([BASE_SEED, idx])
    store = {}
    try:
        for tag, mode, n_rep in (("ordered", "identity", 1), ("shuffled", "free", N_SHUF)):
            obs_l, swp_l, sur_l = [], [], []
            for _ in range(n_rep):
                dense, sizes, nb = sham_dense(pat, band, source, mode, rng)
                if dense is None:
                    continue
                o, w, s = _arc_scores(dense, rng, R)
                obs_l.append(o); swp_l.append(w); sur_l.append(s)
            if not obs_l:
                return dict(patient=pat, band=band, source=source,
                            error="no usable sham realization")
            store[f"obs_{tag}"] = np.array(obs_l)
            store[f"swap_{tag}"] = np.array(swp_l)
            store[f"surr_{tag}"] = np.array(sur_l)
            store[f"sizes_{tag}"] = np.array(sizes)
            store[f"nblocks_{tag}"] = np.array([nb])
    except Exception as exc:                                       # noqa: BLE001
        return dict(patient=pat, band=band, source=source,
                    error=f"{type(exc).__name__}: {exc}")
    cell.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cell, s=SGRID, fracs=np.array(FRACS),
                        funcs=np.array(FUNCS), **store)
    return dict(patient=pat, band=band, source=source, reused=False,
                elapsed_s=time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", default="")
    ap.add_argument("--bands", default="delta,theta,alpha,beta")
    ap.add_argument("--sources", default="rest_pre,rest_post")
    ap.add_argument("--R", type=int, default=0)
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    global R, WORKERS
    R = a.R or R
    WORKERS = a.workers or WORKERS
    if DURS is None:
        raise SystemExit(f"missing {_DP}; run null_ladder/09_cache_durations.py")
    pats = a.patients.split(",") if a.patients else list(PATIENTS_4PHASE)
    bands, srcs = a.bands.split(","), a.sources.split(",")
    (OUT / "cells").mkdir(parents=True, exist_ok=True)
    matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))

    jobs = [(i, p, b, s) for i, (s, b, p) in
            enumerate((s, b, p) for s in srcs for b in bands for p in pats)]
    print(f"[lane-e/E2a] {len(jobs)} cells | sham arcs, ordered + {N_SHUF} shuffled "
          f"| {CANONICAL.label()} fracs={FRACS} | R={R} | {WORKERS} workers -> {OUT}",
          flush=True)
    t0, rows = time.time(), []
    with Pool(WORKERS) as pool:
        for i, res in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            rows.append(res)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {res.get('patient','?')}/{res.get('band','?'):11s}"
                  f"/{res.get('source','?'):9s} {el:6.0f}s ETA {el/i*(len(jobs)-i):6.0f}s"
                  + (f"  ERROR {res['error'][:70]}" if res.get("error") else ""),
                  flush=True)
    pd.DataFrame(rows).to_csv(OUT / "descriptive.csv", index=False)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="lane_e_sham_gate_calibration", substrate=CANONICAL.label(),
        fracs=list(FRACS), R=R, n_shuffled=N_SHUF, block_seconds=BLOCK_SECONDS,
        s_grid=[float(x) for x in SGRID], funcs=list(FUNCS),
        phases=list(CANONICAL_PHASES), cohort=pats, bands=bands, sources=srcs,
        half_nperseg_phases=list(HALF_NPERSEG_PHASES), base_seed=BASE_SEED,
    ), indent=2))
    print(f"\n[lane-e/E2a] done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
