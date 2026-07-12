#!/usr/bin/env python3
"""audit_175c — scale-MAX matched-strength gate on the percolation backbone.

audit_175b showed the trace on the percolation backbone does NOT live at tau_min
(alpha/beta fail there) but at a mesoscale (best_s ~ 6-14). This is the user's
thesis: "result at a different tau." The honest test of that is a MAX-over-scales
statistic with the null treated identically -- the matched-strength surrogate is
ALSO free to pick its best scale, controlling for the winner's-curse of scale
selection (a cluster/max-statistic permutation test).

Per cell, diffusion arm on the percolation backbone:
  obs_max  = max_s rho_sym_obs(s)                      over a pre-committed s-grid
  surr_max = max_s rho_sym_surr(s)   for each surrogate (its OWN best scale)
  p = mean(surr_max >= obs_max)
Cohort gate: Wilcoxon(obs_max - surr_max_p50, greater). Reports the argmax scale.

Outputs: data/audit/sparse_backbone_propagator/percolation_scalemax_gate.csv
"""
from __future__ import annotations
import os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import wilcoxon
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import cophenet, linkage

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import (
    load_phase, shuffle, rho_sym_split, COHORT, BANDS, N_SURROGATES, SWAP_FACTOR, BASE_SEED,
)
from audit_175_sparse_backbone_propagator import PHASES, RHO_FLOOR
from lrg_eegfc.utils.fc.backbone import percolation_backbone, backbone_density

OUT = ROOT / "data" / "audit" / "sparse_backbone_propagator"
S_GRID = np.logspace(0.0, np.log10(50.0), 9)     # pre-committed scale grid, s = tau*lambda_max
R = N_SURROGATES


def rho_sym_over_scales(Wb):
    """rho_sym(s) for s in S_GRID on the 4 phase backbones (eigh once per phase)."""
    eig = {}
    for ph in PHASES:
        deg = Wb[ph].sum(1); ev, V = np.linalg.eigh(np.diag(deg) - Wb[ph]); eig[ph] = (ev, V)
    out = np.empty(len(S_GRID))
    for i, s in enumerate(S_GRID):
        D = {}
        for ph in PHASES:
            ev, V = eig[ph]; tau = s / ev[-1]
            rho = (V * np.exp(-tau * ev)) @ V.T; rho /= np.trace(rho)
            rho = np.where(rho > RHO_FLOOR, rho, RHO_FLOOR)
            T = 1.0 / rho; np.fill_diagonal(T, 0.0); T = np.maximum(T, T.T)
            D[ph] = cophenet(linkage(squareform(T, checks=False), method="average"))
        out[i], _ = rho_sym_split(D["A"], D["B"], D["task_test"], D["rest_post"])
    return out


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
    obs_s = rho_sym_over_scales(Wb)
    obs_max = float(np.nanmax(obs_s)); i_arg = int(np.nanargmax(obs_s))

    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    rng = np.random.default_rng(BASE_SEED + idx)
    surr_max = np.empty(R)
    for r in range(R):
        Wsurr = {ph: shuffle(Ws[ph], n_swaps, rng) for ph in PHASES}
        Wsb = {ph: percolation_backbone(Wsurr[ph])[0] for ph in PHASES}
        s_s = rho_sym_over_scales(Wsb)
        surr_max[r] = np.nanmax(s_s)
    sm = surr_max[np.isfinite(surr_max)]
    return {"patient": pat, "band": band,
            "density": float(backbone_density(Wb["rest_post"])),
            "obs_max": obs_max, "argmax_s": float(S_GRID[i_arg]),
            "obs_at_tau_min": float(obs_s[0]),
            "surr_max_p50": float(np.quantile(sm, 0.50)) if sm.size else float("nan"),
            "surr_max_p95": float(np.quantile(sm, 0.95)) if sm.size else float("nan"),
            "p_scalemax": float(np.mean(sm >= obs_max)) if sm.size else float("nan")}


def cohort_gate(df):
    out = []
    for band in BANDS:
        x = df[df.band == band]
        if x.empty:
            continue
        try:
            _, p = wilcoxon(x["obs_max"].values - x["surr_max_p50"].values, alternative="greater")
        except Exception:
            p = float("nan")
        out.append({"band": band, "n_patients": len(x),
                    "obs_max_median": float(x.obs_max.median()),
                    "surr_max_median": float(x.surr_max_p50.median()),
                    "median_argmax_s": float(x.argmax_s.median()),
                    "obs_taumin_median": float(x.obs_at_tau_min.median()),
                    "n_above": int((x.p_scalemax < 0.05).sum()),
                    "gate_p_scalemax": float(p),
                    "verdict": "CLEAR" if p < 0.05 else "fail"})
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
    print(f"[audit_175c] {len(jobs)} cells, R={R}, scale-MAX gate, percolation backbone, {ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, r in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if r:
                rows.append(r)
                el = time.time() - t0
                print(f"[{i}/{len(jobs)}] {r['patient']}/{r['band']} obs_max={r['obs_max']:+.3f}@s{r['argmax_s']:.1f} "
                      f"(τmin {r['obs_at_tau_min']:+.2f}) surr_p50={r['surr_max_p50']:+.3f} p={r['p_scalemax']:.3f} "
                      f"[{el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s]", flush=True)
    df = pd.DataFrame(rows).sort_values(["band", "patient"]); df.to_csv(OUT / "percolation_scalemax_per_cell.csv", index=False)
    gate = cohort_gate(df); gate.to_csv(OUT / "percolation_scalemax_gate.csv", index=False)
    print(f"\n[audit_175c] {len(rows)} cells in {time.time()-t0:.0f}s")
    print("\n=== SCALE-MAX matched-strength gate (percolation backbone; null also maxes over scale) ===")
    print(gate[["band", "obs_max_median", "surr_max_median", "median_argmax_s",
                "obs_taumin_median", "n_above", "gate_p_scalemax", "verdict"]].to_string(index=False))
    print(f"[audit_175c] outputs -> {OUT}")


if __name__ == "__main__":
    main()
