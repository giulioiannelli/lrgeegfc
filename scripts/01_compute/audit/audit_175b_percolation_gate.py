#!/usr/bin/env python3
"""audit_175b — cohort matched-strength trace gate on the PARAMETER-FREE percolation backbone.

The principled sparsification (user, 2026-07-11): keep every edge >= the connectivity
bottleneck (weakest max-spanning-tree edge) -- the largest global weight threshold at
which P_inf stays 1 (single connected component). Density is SET BY THE DATA (surviving
edge fraction), not chosen; the graph is connected and cycle-rich (not a tree). Proximity
is KEPT (imcoh is volume-conduction-immune -> near-neighbour coupling is genuine signal).

Same harness as audit_175 (reuses diff_coph/geo_coph/gate machinery) but the ONE backbone
is percolation_backbone, applied identically to the observed FC and to each matched-strength
surrogate (its own bottleneck). Reports per band: diffusion vs geodesic gate at tau_min,
data-driven density, C(tau) peak count, and the observed trace at its best scale s.

Outputs:
    data/audit/sparse_backbone_propagator/percolation_gate.csv
"""
from __future__ import annotations
import os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import (
    load_phase, shuffle, rho_sym_split, COHORT, BANDS, N_SURROGATES, SWAP_FACTOR, BASE_SEED,
)
from audit_175_sparse_backbone_propagator import diff_coph, geo_coph, _gate_pack, PHASES, RHO_FLOOR
from lrg_eegfc.utils.fc.backbone import percolation_backbone, backbone_density
from lrg_eegfc.workflow.diagnostics import compute_susceptibility_diagnostics
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import cophenet, linkage

OUT = ROOT / "data" / "audit" / "sparse_backbone_propagator"
S_GRID = np.logspace(0.0, np.log10(60.0), 12)     # scale sweep for best-scale trace
R = N_SURROGATES


def _npeaks(W):
    deg = W.sum(1); ev = np.maximum(np.linalg.eigvalsh(np.diag(deg) - W), 0.0)
    try:
        return int(compute_susceptibility_diagnostics(ev)["N_peaks"])
    except Exception:
        return -1


def _coph_s(ev, V, s):
    tau = s / ev[-1]
    rho = (V * np.exp(-tau * ev)) @ V.T; rho /= np.trace(rho)
    rho = np.where(rho > RHO_FLOOR, rho, RHO_FLOOR)
    T = 1.0 / rho; np.fill_diagonal(T, 0.0); T = np.maximum(T, T.T)
    return cophenet(linkage(squareform(T, checks=False), method="average"))


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception:
        return None
    N = Ws["A"].shape[0]
    if any(v.shape[0] != N for v in Ws.values()):
        return None
    Wb = {ph: percolation_backbone(Ws[ph])[0] for ph in PHASES}
    Dd = {ph: diff_coph(Wb[ph]) for ph in PHASES}
    Dg = {ph: geo_coph(Wb[ph]) for ph in PHASES}
    obs_d, _ = rho_sym_split(Dd["A"], Dd["B"], Dd["task_test"], Dd["rest_post"])
    obs_g, _ = rho_sym_split(Dg["A"], Dg["B"], Dg["task_test"], Dg["rest_post"])

    # best-scale observed trace (result at different tau)
    eig = {}
    for ph in PHASES:
        deg = Wb[ph].sum(1); ev, V = np.linalg.eigh(np.diag(deg) - Wb[ph]); eig[ph] = (np.maximum(ev, 0), V)
    rho_s = []
    for s in S_GRID:
        D = {ph: _coph_s(*eig[ph], s) for ph in PHASES}
        rho_s.append(rho_sym_split(D["A"], D["B"], D["task_test"], D["rest_post"])[0])
    rho_s = np.array(rho_s); i_best = int(np.nanargmax(rho_s))

    # matched-strength null: percolation-sparsify the SAME surrogate draw
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    rng = np.random.default_rng(BASE_SEED + idx)
    surr_d = np.empty(R); surr_g = np.empty(R)
    for r in range(R):
        Wsurr = {ph: shuffle(Ws[ph], n_swaps, rng) for ph in PHASES}
        Wsb = {ph: percolation_backbone(Wsurr[ph])[0] for ph in PHASES}
        Sd = {ph: diff_coph(Wsb[ph]) for ph in PHASES}
        Sg = {ph: geo_coph(Wsb[ph]) for ph in PHASES}
        surr_d[r], _ = rho_sym_split(Sd["A"], Sd["B"], Sd["task_test"], Sd["rest_post"])
        surr_g[r], _ = rho_sym_split(Sg["A"], Sg["B"], Sg["task_test"], Sg["rest_post"])

    rec = {"patient": pat, "band": band, "N_nodes": int(N),
           "density": float(backbone_density(Wb["rest_post"])),
           "npeaks": _npeaks(Wb["A"]),
           "treesim_diff_geo": float(spearmanr(Dd["rest_post"], Dg["rest_post"]).statistic),
           "best_s": float(S_GRID[i_best]), "rho_at_best_s": float(rho_s[i_best]),
           "rho_at_tau_min": float(rho_s[0])}
    rec.update(_gate_pack(obs_d, surr_d, "diff"))
    rec.update(_gate_pack(obs_g, surr_g, "geo"))
    return rec


def cohort_gate(df):
    out = []
    for band in BANDS:
        x = df[df.band == band]
        if x.empty:
            continue
        rec = {"band": band, "n_patients": len(x),
               "median_density": float(x.density.median()),
               "npeaks_median": float(x.npeaks.median()),
               "median_best_s": float(x.best_s.median()),
               "treesim_diff_geo": float(x.treesim_diff_geo.median())}
        for arm in ("diff", "geo"):
            try:
                _, p = wilcoxon(x[f"obs_{arm}"].values - x[f"surr_p50_{arm}"].values, alternative="greater")
            except Exception:
                p = float("nan")
            rec[f"gate_p_{arm}"] = float(p)
            rec[f"obs_median_{arm}"] = float(x[f"obs_{arm}"].median())
            rec[f"n_above_{arm}"] = int((x[f"p_{arm}"] < 0.05).sum())
            rec[f"verdict_{arm}"] = "CLEAR" if p < 0.05 else "fail"
        out.append(rec)
    return pd.DataFrame(out)


def main():
    global R
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--R", type=int, default=N_SURROGATES)
    ap.add_argument("--limit", type=int, default=0); a = ap.parse_args(); R = a.R
    OUT.mkdir(parents=True, exist_ok=True)
    _ = shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in BANDS for p in COHORT)]
    if a.limit:
        jobs = jobs[:a.limit]
    ncpu = min(14, (os.cpu_count() or 4) - 2)
    print(f"[audit_175b] {len(jobs)} cells, R={R}, PERCOLATION backbone, {ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, r in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if r:
                rows.append(r)
                el = time.time() - t0
                print(f"[{i}/{len(jobs)}] {r['patient']}/{r['band']} d={r['density']:.2f} "
                      f"diff={r['obs_diff']:+.3f}(p{r['p_diff']:.3f}) geo={r['obs_geo']:+.3f}(p{r['p_geo']:.3f}) "
                      f"best_s={r['best_s']:.1f}(ρ{r['rho_at_best_s']:+.2f}) pk={r['npeaks']} "
                      f"[{el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s]", flush=True)
    df = pd.DataFrame(rows).sort_values(["band", "patient"]); df.to_csv(OUT / "percolation_per_cell.csv", index=False)
    gate = cohort_gate(df); gate.to_csv(OUT / "percolation_gate.csv", index=False)
    print(f"\n[audit_175b] {len(rows)} cells in {time.time()-t0:.0f}s")
    print("\n=== PERCOLATION-BACKBONE GATE (matched-strength, data-driven density) ===")
    print(gate[["band", "median_density", "npeaks_median", "obs_median_diff", "gate_p_diff",
                "verdict_diff", "gate_p_geo", "verdict_geo", "median_best_s"]].to_string(index=False))
    print(f"[audit_175b] outputs -> {OUT}")


if __name__ == "__main__":
    main()
