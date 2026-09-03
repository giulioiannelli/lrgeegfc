#!/usr/bin/env python3
"""Do different diffusion scales build the trace out of the SAME contact pairs?

The scale axis can be flat in two completely different ways, and the incumbent
statistic cannot tell them apart because it only ever reports a total.

  degenerate  -- the same pairs carry the trace at every scale, so every scale is
                 scored on the same information and the flat curve is a tautology
                 about the readout;
  pervasive   -- different pairs carry it at different scales, and the totals
                 happen to coincide, which is a real multiscale statement (same
                 strength, different carriers) and NOT a failure of the readout.

Spearman is an inner product of normalised centred ranks, so the trace at scale
``s`` owns an exact per-pair contribution vector ``c(s)`` with ``sum_e c_e(s) =
T(s)``. The overlap between ``c(s_i)`` and ``c(s_j)`` is therefore the direct
measurement: cosine similarity (carrier overlap in the estimator's own inner
product), Spearman (rank-robust companion), and the fraction of pairs whose
contribution keeps its sign.

Reference, not a test. This describes the estimator's internal structure rather
than asserting an effect, so it is not gated -- but it is still referenced to the
matched-strength null (a smaller ensemble, ``W0S_CARRIER_R`` realizations), so
"the carriers turn over" can be compared with how much turnover a
strength-matched graph produces on its own. Nothing here is compared to zero.

Observed pass and null pass both knob-integrate over the locked plateau
fractions; the per-scale contribution vectors are computed at each fraction and
the similarity matrices are the median across fractions.
"""
from __future__ import annotations

import os
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet
from scipy.stats import rankdata

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, PATIENTS_4PHASE
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, linkage_at_scale
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle
from lrg_eegfc.utils.metrics.tree import pair_merge_level
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.substrate import CANONICAL, canonical_graph

ROOT = setup_script_env()

from w0s_01_scale_locality_grid import (                      # noqa: E402
    BACKBONE, FRACS, NQ, PHASES, SGRID, SWAP_FACTOR, W_MAX, _contrib, _sp,
    _strata,
)

R_NULL = int(os.environ.get("W0S_CARRIER_R", "20"))
WORKERS = int(os.environ.get("W0S_WORKERS", "8"))
OUT = Path(os.environ.get(
    "W0S_OUT_ANALYSIS", ROOT / "data" / "paper_final" / "lane_s_scale"))
BASE_SEED = 20260831


def _carrier_matrices(Ws: dict) -> tuple:
    """``(cosine, spearman, sign_agreement, split_half_reliability)``, knob-median.

    The last one is the yardstick the cross-scale numbers need. Asking whether a
    Spearman over ~7000 pairs "changes with scale" is meaningless without knowing
    how much the same quantity changes for reasons that have nothing to do with
    scale. ``rho_S(C_A(s), C_B(s))`` is exactly that: the same cophenetic vector,
    same scale, computed on the two contiguous halves of the same pre-task rest
    recording, so the only difference between them is measurement. If moving four
    decades along the scale axis perturbs the vector less than splitting one rest
    recording in half does, the axis is not delivering distinguishable objects to
    the readout, whatever the hierarchy is doing.
    """
    nS = SGRID.size
    acc = np.full((len(FRACS), 3, nS, nS), np.nan)
    #: three rows: the cophenetic vector over ALL pairs (what T correlates), the
    #: same restricted to one tree-level stratum (what Tloc/Qloc correlate,
    #: median over quintiles), and the merge-height profile (what Thei
    #: correlates). One reliability per candidate family, same scale, two halves
    #: of the same rest recording.
    rel = np.full((len(FRACS), 3, nS), np.nan)
    for fi, f in enumerate(FRACS):
        eig = {ph: laplacian_eig(select_backbone(Ws[ph], BACKBONE, frac=float(f)))
               for ph in PHASES}
        C_by_s = []
        for j, s in enumerate(SGRID):
            try:
                Z = {ph: linkage_at_scale(*eig[ph], s) for ph in PHASES}
                C = {ph: cophenet(Z[ph]) for ph in PHASES}
            except Exception:                                    # noqa: BLE001
                C_by_s.append(None)
                continue
            if any((not np.all(np.isfinite(v))) or np.std(v) == 0
                   for v in C.values()):
                C_by_s.append(None)
                continue
            rel[fi, 0, j] = _sp(C["A"], C["B"])
            lev = pair_merge_level(Z["A"])
            _, qlab = _strata(lev, Ws["A"].shape[0])
            qq = [_sp(C["A"][qlab == g], C["B"][qlab == g])
                  for g in range(1, NQ + 1) if int((qlab == g).sum()) >= 30]
            rel[fi, 1, j] = float(np.median(qq)) if qq else np.nan
            rel[fi, 2, j] = _sp(np.sort(Z["A"][:, 2]), np.sort(Z["B"][:, 2]))
            c = 0.5 * (_contrib(C["task_test"] - C["A"], C["rest_post"] - C["B"])
                       + _contrib(C["task_test"] - C["B"], C["rest_post"] - C["A"]))
            C_by_s.append(c)
        ok = [i for i, c in enumerate(C_by_s) if c is not None]
        if len(ok) < 2:
            continue
        X = np.vstack([C_by_s[i] for i in ok])
        Xn = X / np.linalg.norm(X, axis=1, keepdims=True)
        Rk = np.vstack([rankdata(x) for x in X])
        Rk = Rk - Rk.mean(axis=1, keepdims=True)
        Rn = Rk / np.linalg.norm(Rk, axis=1, keepdims=True)
        Sg = np.sign(X)
        idx = np.ix_(ok, ok)
        acc[fi, 0][idx] = Xn @ Xn.T
        acc[fi, 1][idx] = Rn @ Rn.T
        acc[fi, 2][idx] = (Sg @ Sg.T) / X.shape[1] * 0.5 + 0.5
    m = np.nanmedian(acc, axis=0)
    return m[0], m[1], m[2], np.nanmedian(rel, axis=0)  # rel: (3, nS)


def per_cell(job):
    idx, pat, band = job
    t0 = time.time()
    done = OUT / "carrier" / f"{pat}__{band}.npz"
    if done.exists():
        z = np.load(done)
        iu = np.triu_indices(z["s"].size, 1)
        rel = z["split_half_reliability"]
        return dict(patient=pat, band=band, reused=True, elapsed_s=0.0,
                    cos_mean_offdiag=float(np.nanmean(z["cos"][iu])),
                    cos_first_last=float(z["cos"][0, -1]),
                    spearman_mean_offdiag=float(np.nanmean(z["spearman"][iu])),
                    sign_mean_offdiag=float(np.nanmean(z["sign"][iu])),
                    cos_null_mean_offdiag=float(np.nanmean(
                        np.nanmedian(z["cos_null"], axis=0)[iu])),
                    coph_split_half_med=float(np.nanmedian(rel[0])),
                    coph_split_half_min=float(np.nanmin(rel[0])),
                    stratum_split_half_med=float(np.nanmedian(rel[1])),
                    heights_split_half_med=float(np.nanmedian(rel[2])))
    try:
        Ws = {ph: canonical_graph(pat, ph, band, dense=True) for ph in PHASES}
    except Exception as exc:                                     # noqa: BLE001
        return dict(patient=pat, band=band, error=str(exc))
    N = Ws["A"].shape[0]
    cos_o, spr_o, sgn_o, rel_o = _carrier_matrices(Ws)

    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    cos_n = np.full((R_NULL, SGRID.size, SGRID.size), np.nan, dtype=np.float32)
    for r in range(R_NULL):
        Wsh = {ph: matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX)
               for ph in PHASES}
        try:
            cos_n[r] = _carrier_matrices(Wsh)[0]                  # noqa: PLW2901
        except Exception:                                        # noqa: BLE001
            continue
    d = OUT / "carrier"
    d.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(d / f"{pat}__{band}.npz", s=SGRID, cos=cos_o,
                        spearman=spr_o, sign=sgn_o, cos_null=cos_n,
                        split_half_reliability=rel_o)
    iu = np.triu_indices(SGRID.size, 1)
    return dict(patient=pat, band=band, N=N, elapsed_s=time.time() - t0,
                cos_mean_offdiag=float(np.nanmean(cos_o[iu])),
                cos_first_last=float(cos_o[0, -1]),
                spearman_mean_offdiag=float(np.nanmean(spr_o[iu])),
                sign_mean_offdiag=float(np.nanmean(sgn_o[iu])),
                cos_null_mean_offdiag=float(np.nanmean(
                    np.nanmedian(cos_n, axis=0)[iu])),
                coph_split_half_med=float(np.nanmedian(rel_o[0])),
                coph_split_half_min=float(np.nanmin(rel_o[0])),
                stratum_split_half_med=float(np.nanmedian(rel_o[1])),
                heights_split_half_med=float(np.nanmedian(rel_o[2])))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    jobs = [(i, p, b) for i, (b, p) in
            enumerate((b, p) for b in BRAIN_BANDS_NAMES for p in PATIENTS_4PHASE)]
    print(f"[w0s-carrier] {len(jobs)} cells | {CANONICAL.label()} "
          f"| R_null={R_NULL} | {WORKERS} workers", flush=True)
    t0 = time.time()
    rows = []
    with Pool(WORKERS) as pool:
        for i, res in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            rows.append(res)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {res.get('patient','?')}/"
                  f"{str(res.get('band','?')):11s} {el:6.0f}s "
                  f"ETA {el/i*(len(jobs)-i):6.0f}s", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "carrier_overlap_cells.csv", index=False)
    summ = (df.groupby("band")[["cos_mean_offdiag", "cos_first_last",
                                "spearman_mean_offdiag", "sign_mean_offdiag",
                                "cos_null_mean_offdiag", "coph_split_half_med",
                                "coph_split_half_min", "stratum_split_half_med",
                                "heights_split_half_med"]]
            .median().reset_index())
    summ.to_csv(OUT / "carrier_overlap_summary.csv", index=False)
    print("\n-- carrier overlap across scales (cohort median) --")
    print(summ.to_string(index=False))
    print(f"\n[w0s-carrier] done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
