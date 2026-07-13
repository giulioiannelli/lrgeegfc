#!/usr/bin/env python3
"""Controls ladder on mst@0.20 -- what can simple / multiscale-blind descriptors see?

Runs the SAME cross-phase rho_sym functionals (T_test whole-task, T_learn encoding,
T_infspec inference-specific, T_infspec|e inference|encoding) through a ladder of
representations, all under the IDENTICAL matched-strength null, so we can locate WHERE
each signal is readable:

  raw_fc      dense upper-tri edges                      simple pairwise
  strength    node strength (degree)                     simplest node scalar
  clustering  weighted clustering (Onnela)               basic network metric
  geodesic    shortest-path distance on mst@0.20         topological, single-scale
  resistance  effective resistance = L^+ (mst@0.20)      SPECTRAL, multiscale-blind
  coph_taumin cophenetic UPGMA at s=1 (tau_min)          LRG propagator, single-scale
  coph_meso   cophenetic UPGMA at s=5.6 (mesoscale)      LRG propagator, multiscale

Two claims:
  (A) band-selective TRACE (T_test alpha/beta, others null) is UNIQUE to the LRG
      cophenetic -- raw/strength/clustering/geodesic/resistance do not resolve it;
  (B) ENCODING (T_learn) is readable from a basic network metric while INFERENCE
      (T_infspec) is readable only from the multiscale LRG cophenetic -- the
      encoding/inference dissociation reframed as a representation dissociation.

Null: matched-strength shuffle of the DENSE FC (audit_150 bit-identical), then the
identical descriptor pipeline, per cell, R=200, shared draws across descriptors.
Reads per functional; cohort Wilcoxon(obs - surr_p50) one-sided greater, n=10.

Outputs (data/sparsified_arc/controls_ladder/):
  per_cell.csv       per (descriptor, patient, band): obs + surrogate p, 4 functionals
  cohort_gate.csv    per (descriptor, band, functional): Wilcoxon gate
  config.json
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import load_phase, COHORT, BANDS
from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction, geodesic_distance
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, cophenetic_at_scale
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle
from lrg_eegfc.utils.metrics.graph_descriptors import (
    raw_edges, node_strength, weighted_clustering_onnela,
)

OUT = ROOT / "data" / "sparsified_arc" / "controls_ladder"
PHASES5 = ("A", "B", "task_learn", "task_test", "rest_post")
FUNCS = ("T_test", "T_learn", "T_infspec", "T_infspec_pe")
DESCS = ("raw_fc", "strength", "clustering", "geodesic", "resistance",
         "coph_taumin", "coph_meso")
FRAC = 0.20
S_MESO = 5.6
R = 200
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260712
_TRIU = {}


def _triu(M):
    N = M.shape[0]
    key = N
    if key not in _TRIU:
        _TRIU[key] = np.triu_indices(N, k=1)
    return M[_TRIU[key]]


def _resistance_triu(ev, V):
    """Effective-resistance (commute-time) distances = Laplacian pseudoinverse."""
    tol = 1e-10 * ev[-1]
    inv = np.where(ev > tol, 1.0 / np.where(ev > tol, ev, 1.0), 0.0)
    Lp = (V * inv) @ V.T
    d = np.diag(Lp)
    Rm = d[:, None] + d[None, :] - 2.0 * Lp
    return _triu(Rm)


def descriptor_vectors(W):
    """All ladder descriptors for one phase FC matrix W, read on the SAME mst@0.20
    backbone B -- apples-to-apples, so the read-out is the ONLY variable.

    Fixed 2026-07-13: the earlier version read raw_fc/strength/clustering on the
    dense W and geodesic/resistance/cophenetic on B, confounding graph with
    read-out. Now every descriptor sees B. Note the matched-strength null preserves
    *dense* node strength, so strength-on-B is no longer null-invariant (the old
    ``p=1 by construction`` was an artifact of reading strength on the dense graph);
    it is now an honest low-order backbone descriptor.
    """
    B = mst_union_top_fraction(W, FRAC)
    ev, V = laplacian_eig(B)
    out = {"raw_fc": raw_edges(B), "strength": node_strength(B),
           "clustering": weighted_clustering_onnela(B)}
    out["geodesic"] = _triu(geodesic_distance(B))
    out["resistance"] = _resistance_triu(ev, V)
    out["coph_taumin"] = cophenetic_at_scale(ev, V, 1.0)
    out["coph_meso"] = cophenetic_at_scale(ev, V, S_MESO)
    return out


def _rho(a, b):
    r, _ = spearmanr(a, b)
    return float(r)


def _partial(a, b, c):
    rab, rac, rbc = _rho(a, b), _rho(a, c), _rho(b, c)
    den = np.sqrt(max(0.0, (1 - rac ** 2) * (1 - rbc ** 2)))
    return (rab - rac * rbc) / den if den > 0 else np.nan


def _arm(D_TL, D_TT, D_RP, D_bt, D_br):
    e = D_TL - D_bt; g = D_TT - D_bt; f = D_TT - D_TL; p = D_RP - D_br
    return {"T_test": _rho(g, p), "T_learn": _rho(e, p),
            "T_infspec": _rho(f, p), "T_infspec_pe": _partial(f, p, e)}


def _sym(D):
    a1 = _arm(D["task_learn"], D["task_test"], D["rest_post"], D["A"], D["B"])
    a2 = _arm(D["task_learn"], D["task_test"], D["rest_post"], D["B"], D["A"])
    return {k: 0.5 * (a1[k] + a2[k]) for k in FUNCS}


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES5}
    except Exception:
        return None
    N = Ws["A"].shape[0]
    obs_vec = {ph: descriptor_vectors(Ws[ph]) for ph in PHASES5}
    obs_ff = {d: _sym({ph: obs_vec[ph][d] for ph in PHASES5}) for d in DESCS}

    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    surr = {d: {f: np.full(R, np.nan) for f in FUNCS} for d in DESCS}
    for r in range(R):
        sv = {ph: descriptor_vectors(matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX))
              for ph in PHASES5}
        for d in DESCS:
            ff = _sym({ph: sv[ph][d] for ph in PHASES5})
            for f in FUNCS:
                surr[d][f][r] = ff[f]

    rows = []
    for d in DESCS:
        row = {"patient": pat, "band": band, "descriptor": d, "N": int(N)}
        for f in FUNCS:
            s = surr[d][f][np.isfinite(surr[d][f])]; o = obs_ff[d][f]
            row[f"{f}_obs"] = float(o) if np.isfinite(o) else np.nan
            row[f"{f}_surr_p50"] = float(np.median(s)) if s.size else np.nan
            row[f"{f}_p"] = float(np.mean(s >= o)) if s.size and np.isfinite(o) else np.nan
        rows.append(row)
    return rows


def cohort_gate(df):
    out = []
    for d in DESCS:
        for band in BANDS:
            x = df[(df.descriptor == d) & (df.band == band)]
            if x.empty:
                continue
            for f in FUNCS:
                z = x.dropna(subset=[f"{f}_obs", f"{f}_surr_p50"])
                if len(z) < 5:
                    continue
                diff = z[f"{f}_obs"].values - z[f"{f}_surr_p50"].values
                try:
                    p = float(wilcoxon(diff, alternative="greater")[1]) if np.any(diff != 0) else 1.0
                except Exception:
                    p = np.nan
                out.append(dict(descriptor=d, band=band, functional=f, gate_p=p,
                                med_obs=float(z[f"{f}_obs"].median()),
                                med_surr=float(z[f"{f}_surr_p50"].median()),
                                n_above=int((z[f"{f}_p"] < 0.05).sum()), n_pat=len(z)))
    return pd.DataFrame(out)


def main():
    global R, OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--bands", type=str, default="")
    ap.add_argument("--R", type=int, default=0)
    ap.add_argument("--tag", type=str, default="")
    a = ap.parse_args()
    if a.R:
        R = a.R
    if a.tag:
        OUT = OUT.parent / f"controls_ladder_{a.tag}"
    OUT.mkdir(parents=True, exist_ok=True)
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    pats = COHORT[:a.limit] if a.limit else COHORT
    bands = a.bands.split(",") if a.bands else BANDS
    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in bands for p in pats)]
    ncpu = int(os.environ.get("SA_WORKERS", 12))
    print(f"[ladder] {len(jobs)} cells, {len(DESCS)} descriptors, R={R}, {ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, rl in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if rl:
                rows.extend(rl)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s", flush=True)
    df = pd.DataFrame(rows)
    if df.empty:
        print("[ladder] no cells", flush=True); return
    df.to_csv(OUT / "per_cell.csv", index=False)
    gate = cohort_gate(df); gate.to_csv(OUT / "cohort_gate.csv", index=False)
    print(f"\n[ladder] {len(df)} rows in {time.time()-t0:.0f}s -> {OUT}\n", flush=True)
    if gate.empty or "functional" not in gate.columns:
        print("[ladder] cohort gate needs >=5 patients (timing/subset run)", flush=True); return
    for f in FUNCS:
        piv = gate[gate.functional == f].pivot(index="descriptor", columns="band",
                                               values="gate_p").reindex(DESCS)[BANDS]
        print(f"=== {f}: cohort gate_p (descriptor x band) ===", flush=True)
        print(piv.round(3).to_string(), flush=True)
        print("", flush=True)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="controls_ladder_mst020", descriptors=list(DESCS), functionals=list(FUNCS),
        frac=FRAC, s_meso=S_MESO, R=R, cohort=pats, bands=bands,
        null="matched-strength on dense FC -> identical descriptor pipeline, per cell"), indent=2))
    print(f"[ladder] done -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
