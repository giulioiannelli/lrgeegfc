#!/usr/bin/env python3
r"""W0-C master grid: one pass that feeds the cohort gate and all three tau claims.

Every number this lane reports comes from this one compute, so the gate, the
detection test, the characterization test and the selection test all run on
exactly the same graphs, the same scale grid and the same surrogate draws. The
substrate (FC transform x backbone) and the null are injected via environment
variables so the grid can be re-run on the locked pipeline contract without
editing anything.

Per (patient, band) cell, on the four cross-phase graphs {A, B, task_test,
rest_post}, at every scale s in a 28-point log grid spanning 0.05 -> 180:

  coph          rho_sym of the UPGMA cophenetic geometry             (the trace)
  coph_res_lin  rho_sym of cophenetic RANK-LINEARLY residualized on raw edges
  coph_res_np   rho_sym of cophenetic residualized on raw by BINNED conditional
                mean -- removes ANY function of raw, not only monotone ones
  raw_res_lin   rho_sym of raw edges rank-linearly residualized on cophenetic

and once per cell, scale-free:

  raw           rho_sym of the dense upper-triangular FC (the s -> 0 limit)

Each of the five gets the identical matched-strength surrogate ensemble
(R realizations, pipeline-identical: shuffle the DENSE FC, then sparsify, then
LRG, then residualize), so the gate sees (obs, ensemble) for every cell.

Descriptive, observed-only, no null (they are mechanism probes, not tests):

  gini_raw, top{1,5}pct_raw   concentration of the per-pair contributions to the
                              raw trace -- the edge-locality hypothesis
  abl_*                       raw trace after deleting the top-k fraction of
                              contributing pairs (k = 0.1% .. 10%)
  absorb_r2(s)                squared Spearman between the raw per-pair task
                              reorganization and the cophenetic one at s: how
                              much of what raw sees the hierarchy can express
  n_eff(s), m_comm(s)         the resolution unit of the scale axis

5-point critical preamble
1. Claim. (a) DETECTION: no scale of the hierarchy carries cross-phase trace
   that the dense raw FC does not already carry. (b) CHARACTERIZATION: the shape
   of the trace-vs-scale profile differs by band. (c) SELECTION: the hierarchy
   accepts alpha/beta and rejects low_gamma, which raw accepts.
2. Null. For every gated quantity: the matched-strength shuffle of the dense FC
   pushed through the identical sparsify -> Laplacian -> heat kernel -> UPGMA ->
   residualize pipeline. The cohort gate is a signed-rank test on the per-patient
   margin obs - median(surrogates), never on the raw observed value.
3. Strongest alternative each null must control for. For the trace: node strength
   alone reproduces the geometry (matched-strength addresses this by construction
   -- and only this; it cannot reach the connectivity estimator, the band split,
   or session nonstationarity, which is lane W0-B's job). For the DETECTION
   negative: the residualization is too weak to remove raw's contribution, so
   coph|raw "clears" spuriously -- controlled by running the null through the
   same residualizer. For the DETECTION negative the risk runs the other way:
   the residualization is too STRONG or too NARROW, so a real incremental
   contribution is destroyed -- controlled by running two residualizers of very
   different strength (rank-linear removes monotone dependence only; binned
   conditional mean removes every measurable function of raw). If the negative
   survives the weaker one it is not a residualization artifact.
4. What the nulls cannot do. Matched-strength conditions on the observed FC
   matrix: every verdict here is "given this connectivity estimate". The
   cophenetic geometry is a deterministic function of the graph, so coph|raw is
   the hierarchy's own non-monotone content and not an independent measurement;
   a null result there means "no NEW information", never "no structure". n = 10
   caps the signed-rank floor at p = 1/2^10; equivalence claims at this n can
   bound only large effects, which is stated rather than worked around.
5. Falsification. DETECTION is closed negative if coph|raw fails to clear at
   every scale under BOTH residualizers while raw|coph clears. CHARACTERIZATION
   fails if per-patient scale profiles do not separate by band above a
   within-band/across-band permutation null. SELECTION fails if low_gamma clears
   the hierarchy at any scale, or if the edge-locality probe shows low_gamma's
   raw trace is no more concentrated than beta's.

Environment (the injection points)
  W0C_FC_METHOD   default imcoh_abs
  W0C_BACKBONE    default mst020        (any select_backbone name, or "dense")
  W0C_FRAC        default 0.20
  W0C_R           default 200           surrogate realizations
  W0C_WORKERS     default 12
  W0C_OUT         default data/paper_final/w0c_gate_tau/grid
  W0C_RESUME      default 1             reuse existing per-cell surrogate files

Outputs (under W0C_OUT):
  cells/<patient>__<band>.npz   obs + full (R, nS) surrogate ensembles
  descriptive.csv               per-cell mechanism probes and resolution units
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
from scipy.stats import rankdata

from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import (
    cophenetic_at_scale, effective_cluster_count,
    communication_neighbourhood_size, laplacian_eig)
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import CROSS_PHASE_PHASES, phase_graphs

ROOT = setup_script_env()

FC_METHOD = os.environ.get("W0C_FC_METHOD", "imcoh_abs")
BACKBONE = os.environ.get("W0C_BACKBONE", "mst020")
FRAC = float(os.environ.get("W0C_FRAC", "0.20"))
R = int(os.environ.get("W0C_R", "200"))
WORKERS = int(os.environ.get("W0C_WORKERS", "12"))
OUT = Path(os.environ.get("W0C_OUT", ROOT / "data" / "paper_final" / "w0c_gate_tau" / "grid"))
#: Reuse an existing cell file instead of redrawing its surrogate ensemble. The
#: draws are seeded per cell, so a resumed run is identical to an unbroken one.
RESUME = os.environ.get("W0C_RESUME", "1") not in ("0", "false", "False")

PHASES = CROSS_PHASE_PHASES
SGRID = np.logspace(np.log10(0.05), np.log10(180.0), 28)
MEASURES = ("coph", "coph_res_lin", "coph_res_np", "raw_res_lin")
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260825
N_BINS = 64                      # quantile bins of the non-parametric residualizer
ABLATE = (0.001, 0.005, 0.01, 0.02, 0.05, 0.10)

_IU: dict[int, tuple] = {}


# --------------------------------------------------------------------------- #
# fast rank primitives (Spearman without scipy's p-value machinery)
# --------------------------------------------------------------------------- #
def _cr(x: np.ndarray) -> np.ndarray:
    """Mean-centred ranks."""
    r = rankdata(x)
    return r - r.mean()


def _sp_centred(a: np.ndarray, b: np.ndarray) -> float:
    """Pearson of two already-mean-centred vectors (= Spearman if they are ranks)."""
    na, nb = np.sqrt(a @ a), np.sqrt(b @ b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else np.nan


def _sp(a: np.ndarray, b: np.ndarray) -> float:
    return _sp_centred(_cr(a), _cr(b))


def trace_sym(vA, vB, vt, vp) -> float:
    """Symmetric split-half cross-phase trace on four per-pair vectors."""
    return 0.5 * (_sp(vt - vA, vp - vB) + _sp(vt - vB, vp - vA))


def _triu(M: np.ndarray) -> np.ndarray:
    N = M.shape[0]
    if N not in _IU:
        _IU[N] = np.triu_indices(N, 1)
    return M[_IU[N]]


# --------------------------------------------------------------------------- #
# residualizers
# --------------------------------------------------------------------------- #
def resid_rank_linear(x: np.ndarray, c: np.ndarray) -> np.ndarray:
    """rank(x) with its rank-linear (monotone) dependence on ``c`` projected out."""
    rx, rc = _cr(x), _cr(c)
    d = float(rc @ rc)
    return rx - (float(rx @ rc) / d) * rc if d > 0 else rx


def bin_index(c: np.ndarray, n_bins: int = N_BINS) -> np.ndarray:
    """Equal-count quantile bin labels of ``c`` (ties broken by rank order)."""
    order = np.argsort(np.argsort(c))          # 0..M-1 rank positions
    return (order * n_bins // c.size).astype(np.int64)


def resid_binned(x: np.ndarray, bins: np.ndarray, n_bins: int = N_BINS) -> np.ndarray:
    """rank(x) minus its conditional mean within each quantile bin of the covariate.

    Removes *any* function of the covariate resolvable at ``n_bins`` quantiles,
    not merely the monotone part that :func:`resid_rank_linear` removes. The
    strictly stronger of the two residualizers, so the pair brackets the
    question "is the negative an artifact of how we residualized?".
    """
    rx = _cr(x)
    cnt = np.bincount(bins, minlength=n_bins).astype(float)
    tot = np.bincount(bins, weights=rx, minlength=n_bins)
    mean = np.divide(tot, cnt, out=np.zeros_like(tot), where=cnt > 0)
    return rx - mean[bins]


# --------------------------------------------------------------------------- #
# the per-realization computation
# --------------------------------------------------------------------------- #
def _backbone(W: np.ndarray) -> np.ndarray:
    return W if BACKBONE == "dense" else select_backbone(W, BACKBONE, frac=FRAC)


def realization(Ws: dict) -> tuple[np.ndarray, float]:
    """All scale-resolved measures + the raw reference for one set of four graphs.

    Returns ``(vals: (nS, n_measures), raw_rho)``.
    """
    A = {ph: _triu(Ws[ph]) for ph in PHASES}
    eig = {ph: laplacian_eig(_backbone(Ws[ph])) for ph in PHASES}
    bins = {ph: bin_index(A[ph]) for ph in PHASES}

    raw_rho = trace_sym(A["A"], A["B"], A["task_test"], A["rest_post"])

    vals = np.full((SGRID.size, len(MEASURES)), np.nan)
    for j, s in enumerate(SGRID):
        try:
            C = {ph: cophenetic_at_scale(*eig[ph], s) for ph in PHASES}
        except Exception:
            continue
        if any(not np.all(np.isfinite(v)) or np.std(v) == 0 for v in C.values()):
            continue
        cl = {ph: resid_rank_linear(C[ph], A[ph]) for ph in PHASES}
        cn = {ph: resid_binned(C[ph], bins[ph]) for ph in PHASES}
        al = {ph: resid_rank_linear(A[ph], C[ph]) for ph in PHASES}
        for k, d in enumerate((C, cl, cn, al)):
            vals[j, k] = trace_sym(d["A"], d["B"], d["task_test"], d["rest_post"])
    return vals, raw_rho


# --------------------------------------------------------------------------- #
# mechanism probes (observed only, no null)
# --------------------------------------------------------------------------- #
def pair_contributions(vA, vB, vt, vp) -> np.ndarray:
    """Per-pair additive contributions to ``trace_sym``; they sum to the trace.

    Spearman is an inner product of normalised centred ranks, so each pair owns
    a signed share of it. Averaged over the two split-half arms, exactly as the
    estimator is.
    """
    out = np.zeros_like(np.asarray(vA, float))
    for x, y in ((vt - vA, vp - vB), (vt - vB, vp - vA)):
        a, b = _cr(x), _cr(y)
        na, nb = np.sqrt(a @ a), np.sqrt(b @ b)
        if na > 0 and nb > 0:
            out += 0.5 * (a * b) / (na * nb)
    return out


def gini(x: np.ndarray) -> float:
    """Gini concentration of the non-negative part of ``x`` (0 = flat, 1 = one pair)."""
    v = np.sort(np.clip(np.asarray(x, float), 0.0, None))
    n = v.size
    tot = v.sum()
    if n == 0 or tot <= 0:
        return np.nan
    idx = np.arange(1, n + 1)
    return float((2.0 * (idx * v).sum()) / (n * tot) - (n + 1.0) / n)


def edge_locality(A: dict) -> dict:
    """Concentration of the raw trace, and how fast it dies as top pairs are deleted.

    The edge-locality hypothesis says a band whose trace rides on a few strong
    pairs is one that hierarchical coarse-graining will absorb. ``gini`` and the
    top-share are the static form of that; the ablation curve is the dynamic
    form -- delete the top ``k`` fraction of *contributing* pairs and recompute.
    A distributed trace decays gently; a concentrated one collapses.
    """
    c = pair_contributions(A["A"], A["B"], A["task_test"], A["rest_post"])
    M = c.size
    order = np.argsort(-c)
    pos = np.clip(c, 0.0, None)
    tot = pos.sum()
    out = dict(gini_raw=gini(c),
               top1pct_raw=float(pos[order[:max(1, M // 100)]].sum() / tot) if tot > 0 else np.nan,
               top5pct_raw=float(pos[order[:max(1, M // 20)]].sum() / tot) if tot > 0 else np.nan)
    full = trace_sym(A["A"], A["B"], A["task_test"], A["rest_post"])
    out["abl_full"] = full
    for k in ABLATE:
        drop = order[:max(1, int(round(k * M)))]
        keep = np.ones(M, bool)
        keep[drop] = False
        out[f"abl_{k:g}"] = trace_sym(A["A"][keep], A["B"][keep],
                                      A["task_test"][keep], A["rest_post"][keep])
    return out


def reorg(d: dict) -> np.ndarray:
    """Symmetric per-pair task reorganization ``1/2[(t-A) + (t-B)]``."""
    return 0.5 * ((d["task_test"] - d["A"]) + (d["task_test"] - d["B"]))


# --------------------------------------------------------------------------- #
# one cell
# --------------------------------------------------------------------------- #
def per_cell(job):
    idx, pat, band = job
    t0 = time.time()
    try:
        Ws = phase_graphs(pat, band, fc_method=FC_METHOD)
    except Exception as exc:                                   # noqa: BLE001
        return dict(patient=pat, band=band, error=str(exc))
    N = Ws["A"].shape[0]

    # Resume: the surrogate ensemble is ~99% of the cost, so an existing cell
    # file is reused as-is. The descriptive probes below are recomputed either
    # way (they take about a second and are needed for the summary table).
    cell_path = OUT / "cells" / f"{pat}__{band}.npz"
    have = cell_path.exists() and RESUME

    obs_vals, obs_raw = realization(Ws)

    # --- descriptive probes on the observed graphs only -------------------- #
    A = {ph: _triu(Ws[ph]) for ph in PHASES}
    eig = {ph: laplacian_eig(_backbone(Ws[ph])) for ph in PHASES}
    desc = dict(patient=pat, band=band, N=N, obs_raw=obs_raw)
    desc.update(edge_locality(A))
    u_raw = reorg(A)
    absorb = np.full(SGRID.size, np.nan)
    n_eff = np.full(SGRID.size, np.nan)
    m_comm = np.full(SGRID.size, np.nan)
    ev_p, V_p = eig["rest_post"]
    for j, s in enumerate(SGRID):
        try:
            C = {ph: cophenetic_at_scale(*eig[ph], s) for ph in PHASES}
            absorb[j] = _sp(u_raw, reorg(C)) ** 2
        except Exception:
            pass
        n_eff[j] = effective_cluster_count(ev_p, s)
        m_comm[j] = communication_neighbourhood_size(ev_p, V_p, s)

    # --- matched-strength ensemble ----------------------------------------- #
    if have:
        desc["elapsed_s"] = time.time() - t0
        desc["reused"] = True
        return desc
    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    surr = np.full((R, SGRID.size, len(MEASURES)), np.nan)
    surr_raw = np.full(R, np.nan)
    for r in range(R):
        Wsh = {ph: matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX) for ph in PHASES}
        try:
            surr[r], surr_raw[r] = realization(Wsh)
        except Exception:                                       # noqa: BLE001
            continue

    cell_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cell_path,
        s=SGRID, measures=np.array(MEASURES), obs=obs_vals, surr=surr,
        obs_raw=np.array([obs_raw]), surr_raw=surr_raw,
        absorb_r2=absorb, n_eff=n_eff, m_comm=m_comm, N=np.array([N]))
    desc["elapsed_s"] = time.time() - t0
    desc["reused"] = False
    return desc


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
    matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))   # warm numba

    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in bands for p in pats)]
    print(f"[w0c-grid] {len(jobs)} cells | fc={FC_METHOD} backbone={BACKBONE}@{FRAC} "
          f"| R={R} | {SGRID.size} scales s={SGRID[0]:.3g}..{SGRID[-1]:.0f} "
          f"| {WORKERS} workers -> {OUT}", flush=True)
    t0 = time.time()
    rows = []
    with Pool(WORKERS) as pool:
        for i, res in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if res:
                rows.append(res)
            el = time.time() - t0
            tag = f"{res.get('patient','?')}/{res.get('band','?')}" if res else "?"
            print(f"[{i}/{len(jobs)}] {tag:22s} {el:6.0f}s ETA {el/i*(len(jobs)-i):6.0f}s",
                  flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "descriptive.csv", index=False)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="w0c_gate_and_tau_grid", fc_method=FC_METHOD, backbone=BACKBONE,
        frac=FRAC, R=R, swap_factor=SWAP_FACTOR, n_bins=N_BINS,
        s_grid=[float(x) for x in SGRID], measures=list(MEASURES),
        ablate=list(ABLATE), cohort=pats, bands=bands, base_seed=BASE_SEED,
        null="matched-strength on dense FC -> sparsify -> LRG -> cophenetic -> residualize",
    ), indent=2))
    print(f"\n[w0c-grid] done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
