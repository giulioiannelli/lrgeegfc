#!/usr/bin/env python3
r"""Lane W0-S master grid: is the diffusion-scale axis flat, or is the readout blind?

Pre-registration (written before any number, with the full 5-point preamble,
notation, definitions, decision rules and caveats):
``.agents/guides/task-persistence-investigation/2026-08-31_scale-local-trace-readouts.md``.
Nothing in that document may be edited in response to a result.

One pass produces every number the lane reports, so the incumbent statistic, the
dilution diagnostics and all candidate scale-local readouts run on exactly the
same graphs, the same scale grid, the same knob fractions and the same
matched-strength draws.

Per (patient, band), on the four cross-phase graphs {A, B, task_test, rest_post},
at every fraction of the locked plateau f in {0.07, 0.10, 0.14, 0.20} and every
scale s of a 32-point log grid spanning 0.02 -> 180:

  T          incumbent rho_sym / T_probe                          (the control)
  Ccon[o]    additive contribution of tree-OCTAVE o to T; sums to T exactly (R1)
  Tloc[o]    within-octave trace, re-ranked inside the stratum            (R2)
  Qcon[q]    the same contribution decomposition on EQUAL-COUNT strata    (R1')
  Qloc[q]    the same local trace on equal-count strata                   (R2')
  Thei       trace on the log merge-height profile, not pairwise at all   (R3)
  *_diag     each stratified readout evaluated at the stratum the diffusion is
             currently resolving  (derived from N_eff(s), not tuned)

A pair's tree octave is floor(log2 k), where k is the number of clusters present
just before the merge that joins it (k = N before the first merge, k = 2 at the
root): octave 1 is the coarsest split of the implant, the top octave is
hierarchical near-neighbours. Boundaries are fixed powers of two. Because a UPGMA
tree puts most pairs near the root, the top octaves are pair-starved -- so the
identical decomposition is also run on EQUAL-COUNT quintiles of merge level,
which are never starved. Both stratifications are pre-registered; neither was
chosen after seeing a result. The stratification is read off the BASELINE arm's
tree (A for arm 1, B for arm 2) so both members of a cross-phase pair are scored
on the same index set.

Each readout gets the identical matched-strength ensemble, pipeline-identical:
shuffle the DENSE FC at fixed node strength, then sparsify at each fraction, then
Laplacian -> heat kernel -> UPGMA -> stratify -> rank-correlate. One shuffled
dense matrix serves all four fractions of a realization, so fraction-to-fraction
differences are attributable to the knob and not to surrogate noise.

Observed-only diagnostics (no null; they describe the readout, they do not test
anything) answering S1 -- how much of the statistic is the same information at
every scale:

  xs_coph    32x32 cross-scale Spearman of the baseline cophenetic vector
  xs_reorg   the same on the symmetric task reorganisation -- the actual input
  n_distinct effective number of distinct cophenetic values a Spearman over
             ~7000 pairs actually sees (exp of the entropy of the pair-count
             distribution over merge heights)
  pairfrac   share of pairs in each octave
  n_eff, m   the locked scale units, recomputed fresh, never sourced

5-point critical preamble (full version in the scope report)
1. Claim. The near-constancy of the cross-phase trace along the scale axis is a
   property of the READOUT -- a Spearman over all ~7000 pairs dominated by a
   global tree ordering that barely changes with s -- not of the phenomenon.
2. Null. Matched-strength shuffle of the dense FC through the identical
   pipeline; the gate is a one-sided signed-rank on the per-patient margin
   obs - median(own surrogates). NOTHING is tested against zero: W0-B showed all
   four cross-phase functionals are significantly positive at 16/16 scales on a
   no-task sham arc with temporal order destroyed.
3. Strongest alternative. "A scale-local readout is just a noisier readout."
   For the p-value this is controlled by construction (the surrogate passes
   through the same stratification with the same pair counts). For the
   scale-structure claim it is NOT: independent noise decorrelates the columns
   of the margin matrix, so n_eff RISES with noise and an n_eff criterion would
   be satisfied by simply making the statistic worse. Hence n_eff must be
   null-referenced against held-out surrogate realizations of the SAME readout
   -- which is why the full (R, n_frac, nS, n_readout) ensemble is stored here
   rather than collapsed to a median.
4. What the null cannot do. Matched-strength conditions on the finished N x N
   matrix: every verdict is "given this FC estimate" (W0-B's timeseries ladder
   is the live dependency). And a flat result means "no scale structure
   detectable at n = 10", never "scale-invariant" -- equivalence is not
   establishable at this n, and cross-patient inconsistency is indistinguishable
   from within-patient flatness under a cohort test.
5. Falsification. The dilution hypothesis is WRONG if the stratum contribution
   shares are flat (no stratum above 2x the median stratum's share) at most
   scales AND the per-stratum local traces are as cross-scale-correlated as the
   global one. The lane's positive claim FAILS if no readout satisfies all of
   (a) reproduces the incumbent when aggregated, (b) carries more independent
   cross-scale information than its own noise floor, (c) is calibrated.

Environment (the injection points)
  W0S_R           default 100     matched-strength realizations per (cell, frac)
  W0S_WORKERS     default 8
  W0S_OUT         default data/paper_final/lane_s_scale/grid
  W0S_RESUME      default 1       reuse an existing per-cell ensemble file

Outputs (under W0S_OUT):
  cells/<patient>__<band>.npz   obs + full (R, n_frac, nS, n_readout) ensemble
                                + observed-only diagnostics
  descriptive.csv               per-cell summary and timings
  config.json
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
from scipy.cluster.hierarchy import cophenet
from scipy.stats import rankdata

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, PATIENTS_4PHASE
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import (
    communication_neighbourhood_size,
    effective_cluster_count,
    laplacian_eig,
    linkage_at_scale,
)
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle
from lrg_eegfc.utils.metrics.tree import (
    merge_height_pair_counts,
    pair_merge_level,
)
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.substrate import CANONICAL, canonical_graph

ROOT = setup_script_env()

R = int(os.environ.get("W0S_R", "100"))
WORKERS = int(os.environ.get("W0S_WORKERS", "8"))
OUT = Path(os.environ.get(
    "W0S_OUT", ROOT / "data" / "paper_final" / "lane_s_scale" / "grid"))
RESUME = os.environ.get("W0S_RESUME", "1") not in ("0", "false", "False")

#: Four-phase cross-phase arc. The five-phase encoding/inference split is W0-C's;
#: the scale-axis question is asked on the standard trace so it stays comparable
#: to every incumbent number.
PHASES = ("A", "B", "task_test", "rest_post")

#: 32 log-spaced scales. Extends BELOW the locked contract grid (which starts at
#: s = 1) into the near-local regime, because that is the half of the axis the
#: contract grid never samples. As tau -> 0 the kernel is I + tau*W, so D ~ 1/W
#: and the tree is UPGMA on a monotone transform of the raw FC -- the fine end is
#: the raw-FC limit by construction and is reported as a diagnostic of the axis.
SGRID = np.logspace(np.log10(0.02), np.log10(180.0), 32)

FRACS = tuple(CANONICAL.plateau_fracs)
BACKBONE = CANONICAL.backbone
FC_METHOD = CANONICAL.fc_method

#: Octave slots. floor(log2 N) is 6 at the cohort's implant sizes; 8 slots leave
#: headroom and unused slots stay NaN rather than being renormalised away.
NOCT = 8
#: Equal-count strata of merge level -- the starvation-free companion to octaves.
NQ = 5
#: A within-stratum rank correlation on fewer pairs than this is not estimated.
#: The additive contribution is always kept: it is a sum, exact at any count.
MIN_PAIRS_LOCAL = 30

SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260831

READOUTS = (["T"]
            + [f"Ccon_o{o}" for o in range(1, NOCT + 1)]
            + [f"Tloc_o{o}" for o in range(1, NOCT + 1)]
            + [f"Qcon_q{q}" for q in range(1, NQ + 1)]
            + [f"Qloc_q{q}" for q in range(1, NQ + 1)]
            + ["Thei", "Ccon_diag", "Tloc_diag", "Qcon_diag", "Qloc_diag"])
NREAD = len(READOUTS)
IX = {k: i for i, k in enumerate(READOUTS)}


# --------------------------------------------------------------------------- #
# rank primitives (Spearman as an inner product of normalised centred ranks)
# --------------------------------------------------------------------------- #
def _cr(x: np.ndarray) -> np.ndarray:
    r = rankdata(x)
    return r - r.mean()


def _sp(a: np.ndarray, b: np.ndarray) -> float:
    ca, cb = _cr(a), _cr(b)
    na, nb = np.sqrt(ca @ ca), np.sqrt(cb @ cb)
    return float(ca @ cb / (na * nb)) if na > 0 and nb > 0 else np.nan


def _contrib(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Per-element additive contributions to ``Spearman(x, y)``; they sum to it."""
    a, b = _cr(x), _cr(y)
    na, nb = np.sqrt(a @ a), np.sqrt(b @ b)
    if na <= 0 or nb <= 0:
        return np.full(x.shape, np.nan)
    return (a * b) / (na * nb)


def _shannon_neff(counts: np.ndarray) -> float:
    """``exp`` of the Shannon entropy of a count vector -- effective support size."""
    c = np.asarray(counts, float)
    c = c[c > 0]
    if c.size == 0:
        return np.nan
    p = c / c.sum()
    return float(np.exp(-np.sum(p * np.log(p))))


# --------------------------------------------------------------------------- #
# stratifications of the pair set, both read off one baseline tree
# --------------------------------------------------------------------------- #
def _strata(lev: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    """``(octave, quintile)`` labels for every pair, from its merge level.

    ``lev`` is the 1-based merge index of each pair (:func:`pair_merge_level`).
    ``k = n - lev + 1`` is the number of clusters standing just before that
    merge, so the octave ``floor(log2 k)`` runs from 1 (root-adjacent pairs) to
    ``floor(log2 n)`` (hierarchical near-neighbours). The quintile is an
    equal-count split of the same ordering, so every quintile holds ``M/5``
    pairs whatever the tree's shape.
    """
    k = np.maximum(n - lev + 1, 2)
    oct_ = np.clip(np.floor(np.log2(k)).astype(np.int64), 1, NOCT)
    order = np.argsort(rankdata(lev, method="ordinal"))
    q = np.empty(lev.size, dtype=np.int64)
    q[order] = (np.arange(lev.size) * NQ // lev.size) + 1
    return oct_, q


def _active_stratum(lev: np.ndarray, n: int, n_eff: float) -> tuple[int, int]:
    """The octave and quintile the diffusion is currently resolving.

    A diffusion resolving ``n_eff`` components is separating exactly the pairs
    whose merge stands at ``k ~ n_eff`` clusters, i.e. merge level
    ``l* = n - n_eff + 1``. The octave is ``floor(log2 n_eff)``; the quintile is
    the one holding the pairs at ``l*``. Both are consequences of the locked
    scale unit ``N_eff(s)``, not fitted windows.
    """
    ne = float(np.clip(n_eff, 2.0, n))
    o = int(np.clip(np.floor(np.log2(ne)), 1, NOCT))
    lstar = n - ne + 1.0
    p = float(np.mean(lev <= lstar))
    q = int(np.clip(np.ceil(p * NQ), 1, NQ))
    return o, q


def _arm(dx: np.ndarray, dy: np.ndarray, lab: np.ndarray, nlab: int) -> tuple:
    """Contributions and within-stratum traces for one split-half arm.

    ``dx`` is the task-referenced change, ``dy`` the cross-baseline persistence,
    ``lab`` the stratum label of every pair read off this arm's baseline tree.
    Returns ``(T, con[nlab], loc[nlab], counts[nlab])`` with ``sum(con) == T``.
    """
    c = _contrib(dx, dy)
    T = float(np.nansum(c))
    con = np.full(nlab, np.nan)
    loc = np.full(nlab, np.nan)
    cnt = np.zeros(nlab, dtype=np.int64)
    for i in range(nlab):
        m = lab == (i + 1)
        k = int(m.sum())
        cnt[i] = k
        if k == 0:
            continue
        con[i] = float(np.nansum(c[m]))
        if k >= MIN_PAIRS_LOCAL:
            loc[i] = _sp(dx[m], dy[m])
    return T, con, loc, cnt


# --------------------------------------------------------------------------- #
# one realization: all fractions x all scales
# --------------------------------------------------------------------------- #
def realization(Ws: dict, want_diag: bool = False):
    """All readouts for one set of four dense phase matrices.

    Returns ``(vals, diag)`` with ``vals`` of shape ``(n_frac, nS, NREAD)``.
    ``diag`` is ``None`` unless ``want_diag`` (the observed pass only).
    """
    nF, nS = len(FRACS), SGRID.size
    n = Ws["A"].shape[0]
    vals = np.full((nF, nS, NREAD), np.nan)
    diag = None
    if want_diag:
        diag = dict(
            n_eff=np.full((nF, nS), np.nan), m_comm=np.full((nF, nS), np.nan),
            n_distinct=np.full((nF, nS), np.nan),
            pairfrac=np.full((nF, nS, NOCT), np.nan),
            ostar=np.full((nF, nS), np.nan), qstar=np.full((nF, nS), np.nan),
            xs_coph=np.full((nF, nS, nS), np.nan),
            xs_reorg=np.full((nF, nS, nS), np.nan),
        )

    for fi, f in enumerate(FRACS):
        G = {ph: select_backbone(Ws[ph], BACKBONE, frac=float(f)) for ph in PHASES}
        eig = {ph: laplacian_eig(G[ph]) for ph in PHASES}
        nbar = np.array([0.5 * (effective_cluster_count(eig["A"][0], s)
                                + effective_cluster_count(eig["B"][0], s))
                         for s in SGRID])

        rk_coph, rk_reorg = [], []
        for j, s in enumerate(SGRID):
            try:
                Z = {ph: linkage_at_scale(*eig[ph], s) for ph in PHASES}
            except Exception:                                       # noqa: BLE001
                continue
            C = {ph: cophenet(Z[ph]) for ph in PHASES}
            if any((not np.all(np.isfinite(v))) or np.std(v) == 0
                   for v in C.values()):
                continue
            lev = {ph: pair_merge_level(Z[ph]) for ph in ("A", "B")}
            lab = {ph: _strata(lev[ph], n) for ph in ("A", "B")}
            act = {ph: _active_stratum(lev[ph], n, nbar[j]) for ph in ("A", "B")}

            dxy = {"A": (C["task_test"] - C["A"], C["rest_post"] - C["B"]),
                   "B": (C["task_test"] - C["B"], C["rest_post"] - C["A"])}
            res = {}
            for si, (nlab, tag) in enumerate(((NOCT, "o"), (NQ, "q"))):
                acc_T, acc_con, acc_loc, acc_cnt, acc_dc, acc_dl = [], [], [], [], [], []
                for ph in ("A", "B"):
                    T_, con_, loc_, cnt_ = _arm(*dxy[ph], lab[ph][si], nlab)
                    acc_T.append(T_)
                    acc_con.append(con_)
                    acc_loc.append(loc_)
                    acc_cnt.append(cnt_)
                    idx = act[ph][si] - 1
                    acc_dc.append(con_[idx])
                    acc_dl.append(loc_[idx])
                res[tag] = (float(np.mean(acc_T)),
                            0.5 * (acc_con[0] + acc_con[1]),
                            0.5 * (acc_loc[0] + acc_loc[1]),
                            acc_cnt[0],
                            0.5 * (acc_dc[0] + acc_dc[1]),
                            0.5 * (acc_dl[0] + acc_dl[1]))

            vals[fi, j, IX["T"]] = res["o"][0]
            for o in range(NOCT):
                vals[fi, j, IX[f"Ccon_o{o+1}"]] = res["o"][1][o]
                vals[fi, j, IX[f"Tloc_o{o+1}"]] = res["o"][2][o]
            for q in range(NQ):
                vals[fi, j, IX[f"Qcon_q{q+1}"]] = res["q"][1][q]
                vals[fi, j, IX[f"Qloc_q{q+1}"]] = res["q"][2][q]
            vals[fi, j, IX["Ccon_diag"]] = res["o"][4]
            vals[fi, j, IX["Tloc_diag"]] = res["o"][5]
            vals[fi, j, IX["Qcon_diag"]] = res["q"][4]
            vals[fi, j, IX["Qloc_diag"]] = res["q"][5]

            H = {ph: np.log(np.maximum(np.sort(Z[ph][:, 2]), 1e-300))
                 for ph in PHASES}
            vals[fi, j, IX["Thei"]] = 0.5 * (
                _sp(H["task_test"] - H["A"], H["rest_post"] - H["B"])
                + _sp(H["task_test"] - H["B"], H["rest_post"] - H["A"]))

            if want_diag:
                diag["n_eff"][fi, j] = nbar[j]
                diag["m_comm"][fi, j] = communication_neighbourhood_size(*eig["A"], s)
                diag["n_distinct"][fi, j] = _shannon_neff(
                    merge_height_pair_counts(Z["A"]))
                tot = res["o"][3].sum()
                if tot:
                    diag["pairfrac"][fi, j] = res["o"][3] / tot
                diag["ostar"][fi, j] = act["A"][0]
                diag["qstar"][fi, j] = act["A"][1]
                rk_coph.append(_cr(C["A"]))
                rk_reorg.append(_cr(0.5 * ((C["task_test"] - C["A"])
                                           + (C["task_test"] - C["B"]))))

        if want_diag and len(rk_coph) == nS:
            for key, mats in (("xs_coph", rk_coph), ("xs_reorg", rk_reorg)):
                Xr = np.vstack(mats)
                Xn = Xr / np.linalg.norm(Xr, axis=1, keepdims=True)
                diag[key][fi] = Xn @ Xn.T
    return vals, diag


# --------------------------------------------------------------------------- #
# one cell
# --------------------------------------------------------------------------- #
def per_cell(job):
    idx, pat, band = job
    t0 = time.time()
    cell_path = OUT / "cells" / f"{pat}__{band}.npz"
    if cell_path.exists() and RESUME:
        return dict(patient=pat, band=band, reused=True,
                    elapsed_s=time.time() - t0)
    try:
        Ws = {ph: canonical_graph(pat, ph, band, dense=True) for ph in PHASES}
    except Exception as exc:                                        # noqa: BLE001
        return dict(patient=pat, band=band, error=str(exc))
    N = Ws["A"].shape[0]

    obs, diag = realization(Ws, want_diag=True)

    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    surr = np.full((R, len(FRACS), SGRID.size, NREAD), np.nan, dtype=np.float32)
    for r in range(R):
        Wsh = {ph: matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX)
               for ph in PHASES}
        try:
            surr[r] = realization(Wsh)[0]
        except Exception:                                           # noqa: BLE001
            continue

    cell_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cell_path, s=SGRID, fracs=np.array(FRACS), readouts=np.array(READOUTS),
        obs=obs, surr=surr, N=np.array([N]), **diag)
    return dict(patient=pat, band=band, N=N, reused=False,
                elapsed_s=time.time() - t0,
                obs_T_med=float(np.nanmedian(obs[:, :, IX["T"]])))


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

    jobs = [(i, p, b) for i, (b, p) in
            enumerate((b, p) for b in bands for p in pats)]
    print(f"[w0s-grid] {len(jobs)} cells | {FC_METHOD} x {BACKBONE} "
          f"f={FRACS} | R={R} | {SGRID.size} scales "
          f"s={SGRID[0]:.3g}..{SGRID[-1]:.0f} | {NREAD} readouts "
          f"| {WORKERS} workers -> {OUT}", flush=True)
    t0 = time.time()
    rows = []
    with Pool(WORKERS) as pool:
        for i, res in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            rows.append(res)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {res.get('patient','?')}/"
                  f"{str(res.get('band','?')):11s} {el:6.0f}s "
                  f"ETA {el/i*(len(jobs)-i):6.0f}s", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "descriptive.csv", index=False)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="w0s_scale_locality_grid",
        scope=".agents/guides/task-persistence-investigation/"
              "2026-08-31_scale-local-trace-readouts.md",
        substrate=CANONICAL.label(), fc_method=FC_METHOD, backbone=BACKBONE,
        fracs=[float(f) for f in FRACS], plateau=list(CANONICAL.plateau),
        R=R, swap_factor=SWAP_FACTOR, phases=list(PHASES),
        s_grid=[float(x) for x in SGRID], readouts=list(READOUTS),
        n_octave_slots=NOCT, n_quantile_strata=NQ,
        min_pairs_local=MIN_PAIRS_LOCAL, cohort=pats, bands=bands,
        base_seed=BASE_SEED,
        null="matched-strength on dense FC -> sparsify(each f) -> LRG -> "
             "cophenetic -> stratify -> rank-correlate",
    ), indent=2))
    print(f"\n[w0s-grid] done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
