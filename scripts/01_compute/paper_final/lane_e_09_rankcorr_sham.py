#!/usr/bin/env python3
r"""Lane E / E5-E8 -- the no-signal sham arcs, stored as rank-correlation tensors.

The reference distribution for every claim this lane makes. Five pseudo-phases
are carved out of ONE resting recording -- so no task, and therefore no
consolidation, can exist inside any of them -- and pushed through the entire
pipeline the real arc goes through: the same backbone at the same four plateau
fractions, the same scale grid, and the same matched-strength null on top.

Three arms, and the contrast between them is a mechanism ladder rather than a
robustness check:

  ordered    five contiguous windows in TRUE temporal order, sized in proportion
             to the real phase durations. Retains session drift and phase
             adjacency. This is the arm that fires.
  shuffled   the same five window sizes reassembled from randomly reassigned
             30 s blocks. Destroys temporal order, hence drift.
  eqdur      shuffled, plus the two task pseudo-phases given equal duration.
             ``task_test`` is longer than ``task_learn`` in all 10 patients
             (ratio 0.40-0.74), so the two task blocks are not exchangeable on
             duration alone; this arm isolates that mechanism because duration
             is the only thing it changes.

Reading the ladder: an effect present in ``ordered`` and absent in ``shuffled``
is temporal order. An effect present in ``shuffled`` and absent in ``eqdur`` is
the duration asymmetry. An effect present in all three is the construction
itself. An effect in none of them is a candidate for being real.

Difference from ``lane_e_02_sham_gate_calibration.py``, which produced the
4-band cells this supersedes: it stores the ``(8, 8)`` rank-correlation tensor
rather than four collapsed functionals, so the drift-controlled estimators and
any future one are recoverable from the sham exactly as they are from the real
arc -- which is what makes the two comparable at all. It also covers all six
bands: without low gamma and high gamma no claim that a dissociation is
beta-specific or delta-specific is admissible. To buy those two bands in the
same wall-clock the shuffled and equal-duration arms drop to 2 and 1
realizations; the ORDERED arm, the load-bearing one, is unchanged, and the
previous cells remain on disk as an independent check of the arms that shrank.

The recording is loaded as float32. The rFFT coefficients are stored as
``complex64`` by ``segment_fft`` regardless, so the raw timeseries is the memory
peak and halving it is what keeps concurrent workers inside the box;
``--verify-dtype`` measures the resulting difference in the dense FC rather than
assuming it is negligible.

5-point critical preamble: ``lane_e_07_preregistration_addendum.md``.

Outputs: data/paper_final/lane_e_encinf/rankcorr_sham/cells/<pat>__<band>__<src>.npz
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
    BRAIN_BANDS_NAMES,
    DEFAULT_SAMPLE_RATE,
    FS_OVERRIDES,
    PATIENTS_4PHASE,
    nperseg_for_fs,
)
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import (
    CROSS_PHASE_VECTOR_NAMES,
    cross_phase_rank_corr_over_scales,
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
WORKERS = int(os.environ.get("LANE_E_WORKERS", "3"))
STAGGER_S = float(os.environ.get("LANE_E_STAGGER_S", "60"))
TS_DTYPE = np.dtype(os.environ.get("LANE_E_TS_DTYPE", "float32"))
OUT = Path(os.environ.get("LANE_E_OUT",
                          ROOT / "data" / "paper_final" / "lane_e_encinf")) / "rankcorr_sham"
RESUME = os.environ.get("LANE_E_RESUME", "1") not in ("0", "false", "False")

SGRID = canonical_scale_grid()
FRACS = tuple(float(f) for f in CANONICAL.plateau_fracs)
SWAP_FACTOR, W_MAX = 20, 1.0
BLOCK_SECONDS = 30.0
BASE_SEED = 20260903
HALF_NPERSEG_PHASES = ("A", "B")     # mirrors the real pipeline's split-half arms

#: (tag, block mode, n realizations, equalise task durations)
ARMS = (("ordered", "identity", 1, False),
        ("shuffled", "free", 2, False),
        ("eqdur", "free", 1, True))

_DP = ROOT / "data" / "paper_final" / "w0b_nulls" / "durations.json"
DURS = json.loads(_DP.read_text()) if _DP.exists() else None


# --------------------------------------------------------------------------- #
def _tensor(Ws: dict) -> np.ndarray:
    eig = {ph: laplacian_eig(Ws[ph]) for ph in CANONICAL_PHASES}
    return cross_phase_rank_corr_over_scales(eig, SGRID)


def _arc_tensors(dense: dict, rng, n_surr: int):
    """``(Robs, Rsurr)`` for one five-phase set of DENSE matrices.

    Identical code path to the real arc: sparsify at every plateau fraction, and
    draw the matched-strength surrogate on the DENSE matrix before sparsifying,
    with one draw shared across the fractions.
    """
    nF, nS, nV = len(FRACS), SGRID.size, len(CROSS_PHASE_VECTOR_NAMES)
    Robs = np.full((nF, nS, nV, nV), np.nan, np.float32)
    for i, fr in enumerate(FRACS):
        Robs[i] = _tensor({ph: select_backbone(dense[ph], CANONICAL.backbone, frac=fr)
                           for ph in CANONICAL_PHASES})
    N = dense["A"].shape[0]
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    Rsurr = np.full((nF, n_surr, nS, nV, nV), np.nan, np.float32)
    for r in range(n_surr):
        sh = {ph: matched_strength_shuffle(dense[ph], n_swaps, rng, W_MAX)
              for ph in CANONICAL_PHASES}
        for i, fr in enumerate(FRACS):
            Rsurr[i, r] = _tensor(
                {ph: select_backbone(sh[ph], CANONICAL.backbone, frac=fr)
                 for ph in CANONICAL_PHASES})
    return Robs, Rsurr


def _block_of_segment(n_seg: int, nperseg: int, blk: int) -> np.ndarray:
    """Block index of each Welch segment; ``-1`` for segments straddling a block."""
    seg_start = np.arange(n_seg) * (nperseg // 2)
    b0 = seg_start // blk
    b1 = (seg_start + nperseg - 1) // blk
    return np.where(b0 == b1, b0, -1).astype(np.int64)


def band_grids(pat: str, source: str, bands, dtype=None) -> dict:
    """``{band: {"full"/"half": (F, scale, block_of_segment)}}`` from ONE load.

    The recording is the expensive object and the segment rFFT over it is the
    dominant cost, so both are paid once: one transform per ``nperseg`` covers
    the union of the requested bands and each band is a bin slice of it -- six
    bands cost two passes, not twelve. Both ``nperseg`` grids are built because
    the real pipeline estimates the split-half arms at ``nperseg // 2``.
    """
    fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
    nps = nperseg_for_fs(fs)
    nps_h = max(256, nps // 2)
    blk = int(round(BLOCK_SECONDS * fs))
    lo = min(BRAIN_BANDS[b][0] for b in bands)
    hi = max(BRAIN_BANDS[b][1] for b in bands)
    X = np.asarray(load_timeseries(pat, source, SEEG_DATAPATH), dtype or TS_DTYPE)
    if X.shape[0] > X.shape[1]:
        X = X.T
    out = {b: {} for b in bands}
    for tag, nperseg in (("full", nps), ("half", nps_h)):
        freqs, F, sc = segment_fft(X, fs, nperseg, band=(lo, hi))
        bos = _block_of_segment(F.shape[1], nperseg, blk)
        for band in bands:
            f0, f1 = BRAIN_BANDS[band]
            m = (freqs >= f0) & (freqs <= f1)
            if not m.any():
                raise ValueError(f"no rFFT bins in {band} at nperseg={nperseg}")
            out[band][tag] = (np.ascontiguousarray(F[:, :, m]), sc, bos)
        del F
    del X
    return out


def sham_dense(pat: str, grids: dict, mode: str, rng,
               equal_task_durations: bool = False):
    """Five DENSE sham FC matrices carved out of ONE resting recording."""
    d = DURS[pat]
    prop = np.array([d["rest_pre"] / 2, d["rest_pre"] / 2, d["task_learn"],
                     d["task_test"], d["rest_post"]], float)
    if equal_task_durations:
        prop[2] = prop[3] = 0.5 * (prop[2] + prop[3])
    prop /= prop.sum()
    nb = int(min(int(g[2].max()) for g in grids.values())) + 1
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
        W = np.asarray(imcoh_abs_from_coherency(
            coherency_from_csd(csd_from_segment_subset(F, idx, sc))), float)
        np.fill_diagonal(W, 0.0)
        W = np.clip(W, 0.0, 1.0)
        dense[ph] = 0.5 * (W + W.T)
    return dense, sizes.tolist(), nb


def per_cell(job):
    """One (patient, source): read the recording once, score every band's arcs.

    Worker starts are staggered because ``systemd-oomd`` kills the whole compute
    slice on memory *pressure*, not on a hard limit, and simultaneous first-wave
    loads are exactly the spike that triggers it.
    """
    idx, pat, source, bands = job
    t0 = time.time()
    todo = [b for b in bands
            if not ((OUT / "cells" / f"{pat}__{b}__{source}.npz").exists() and RESUME)]
    if not todo:
        return dict(patient=pat, source=source, bands=len(bands), reused=True)
    if idx < WORKERS and STAGGER_S:
        time.sleep(STAGGER_S * idx)
    rng = np.random.default_rng([BASE_SEED, idx])
    try:
        grids = band_grids(pat, source, todo)
    except Exception as exc:                                       # noqa: BLE001
        return dict(patient=pat, source=source,
                    error=f"load/fft {type(exc).__name__}: {exc}")
    done = []
    for band in todo:
        store = {}
        try:
            for tag, mode, n_rep, eqd in ARMS:
                obs_l, sur_l = [], []
                for _ in range(n_rep):
                    dense, sizes, nb = sham_dense(pat, grids[band], mode, rng, eqd)
                    if dense is None:
                        continue
                    o, s = _arc_tensors(dense, rng, R)
                    obs_l.append(o)
                    sur_l.append(s)
                if not obs_l:
                    raise RuntimeError(f"no usable {tag} realization")
                store[f"Robs_{tag}"] = np.array(obs_l)
                store[f"Rsurr_{tag}"] = np.array(sur_l)
                store[f"sizes_{tag}"] = np.array(sizes)
                store[f"nblocks_{tag}"] = np.array([nb])
        except Exception as exc:                                   # noqa: BLE001
            return dict(patient=pat, source=source, band=band,
                        error=f"{type(exc).__name__}: {exc}")
        cell = OUT / "cells" / f"{pat}__{band}__{source}.npz"
        cell.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cell, s=SGRID, fracs=np.array(FRACS),
                            vectors=np.array(CROSS_PHASE_VECTOR_NAMES), **store)
        done.append(band)
        del grids[band]
    return dict(patient=pat, source=source, bands=len(done), reused=False,
                elapsed_s=time.time() - t0)


# --------------------------------------------------------------------------- #
def verify_dtype(pat: str, source: str, band: str) -> int:
    """Measure what loading the recording as float32 costs, rather than assume.

    Builds the ordered sham's five DENSE matrices from a float64 load and from a
    float32 load with the same block layout, and reports the largest absolute
    disagreement. ``segment_fft`` stores its coefficients as ``complex64``
    whatever the input, so the expectation is a difference at that precision.
    """
    outs = {}
    for dt in (np.float64, np.float32):
        g = band_grids(pat, source, [band], dtype=dt)[band]
        dense, _, _ = sham_dense(pat, g, "identity", np.random.default_rng(0))
        outs[dt] = dense
        del g
    worst = max(float(np.abs(outs[np.float64][ph] - outs[np.float32][ph]).max())
                for ph in CANONICAL_PHASES)
    rel = max(float(np.abs(outs[np.float64][ph]).max()) for ph in CANONICAL_PHASES)
    print(f"[dtype] {pat}/{source}/{band}: max |float64 - float32| = {worst:.3e} "
          f"on entries up to {rel:.3f}  (complex64 eps ~ 1.2e-7)", flush=True)
    return 0 if worst < 1e-5 else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", default="")
    ap.add_argument("--bands", default="")
    ap.add_argument("--sources", default="rest_pre,rest_post")
    ap.add_argument("--R", type=int, default=0)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--verify-dtype", action="store_true")
    a = ap.parse_args()
    global R, WORKERS
    R = a.R or R
    WORKERS = a.workers or WORKERS
    if DURS is None:
        raise SystemExit(f"missing {_DP}; run null_ladder/09_cache_durations.py")
    pats = a.patients.split(",") if a.patients else list(PATIENTS_4PHASE)
    bands = a.bands.split(",") if a.bands else list(BRAIN_BANDS_NAMES)
    srcs = a.sources.split(",")
    if a.verify_dtype:
        raise SystemExit(verify_dtype(pats[0], srcs[0], bands[0]))

    (OUT / "cells").mkdir(parents=True, exist_ok=True)
    matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))

    jobs = [(i, p, s, bands) for i, (s, p) in
            enumerate((s, p) for s in srcs for p in pats)]
    arms = " + ".join(f"{n} {t}" for t, _, n, _ in ARMS)
    print(f"[lane-e/sham] {len(jobs)} (patient, source) jobs x {len(bands)} bands "
          f"| arms: {arms} | {CANONICAL.label()} fracs={FRACS} | R={R} "
          f"| ts dtype {TS_DTYPE} | {WORKERS} workers -> {OUT}", flush=True)
    t0, rows = time.time(), []
    with Pool(WORKERS) as pool:
        for i, res in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            rows.append(res)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {res.get('patient','?')}"
                  f"/{res.get('source','?'):9s} {el:6.0f}s ETA {el/i*(len(jobs)-i):6.0f}s"
                  + (f"  ERROR {res['error'][:70]}" if res.get("error") else ""),
                  flush=True)
    pd.DataFrame(rows).to_csv(OUT / "descriptive.csv", index=False)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="lane_e_rankcorr_sham", substrate=CANONICAL.label(),
        fracs=list(FRACS), R=R, block_seconds=BLOCK_SECONDS,
        arms={t: dict(mode=m, n_realizations=n, equal_task_durations=e)
              for t, m, n, e in ARMS},
        s_grid=[float(x) for x in SGRID], vectors=list(CROSS_PHASE_VECTOR_NAMES),
        phases=list(CANONICAL_PHASES), cohort=pats, bands=bands, sources=srcs,
        half_nperseg_phases=list(HALF_NPERSEG_PHASES), ts_dtype=str(TS_DTYPE),
        base_seed=BASE_SEED,
    ), indent=2))
    print(f"\n[lane-e/sham] done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
