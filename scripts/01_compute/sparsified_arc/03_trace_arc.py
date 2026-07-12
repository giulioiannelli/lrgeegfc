#!/usr/bin/env python3
"""D2 -- tau-resolved cross-phase trace arc on the percolation backbone.

The headline. For every (patient, band) cell:
  - percolation backbone per phase (A,B,task_test,rest_post), eigendecompose once;
  - observed rho_sym(s) swept over the dimensionless scale grid s = tau*lambda_max
    (s=1 == tau_min == the old single-scale anchor);
  - matched-strength null: R surrogates, each shuffled (audit_150 rng/order),
    percolation-sparsified the SAME way, rho_sym(s) recomputed -> the null is free
    to express whatever trace strength-preservation allows at EVERY scale.

Per-scale gate p(s) = mean(surr(s) >= obs(s)). Scale-MAX gate (controls scale
selection): obs_max = max_s obs(s); each surrogate maxes over its OWN best scale;
p = mean(surr_max >= obs_max). Cohort gate per band = Wilcoxon(obs-surr_p50,
greater) at each s, plus the scale-max Wilcoxon. Band trace NATURE classified
from the cohort fire-set F = {s : gate_p(s) < 0.05}:
  none / single-scale (<1/3 grid, one band) / multiscale (>=2 bands)
  / continuous-multiscale (>=2/3 grid).

s=1 reproduces the audit_150 gate (anchor) up to the backbone (percolation vs full).

Outputs:
    data/sparsified_arc/trace_arc/{band}/{patient}.npz    (s, obs, surr_p50/p95, p_of_s)
    data/sparsified_arc/trace_arc/per_cell.csv
    data/sparsified_arc/trace_arc/cohort_gate_vs_s.csv     (band x s: gate_p, medians)
    data/sparsified_arc/trace_arc/band_nature.csv          (band: class, gates)
    data/sparsified_arc/trace_arc/config.json
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import (
    load_phase, COHORT, BANDS, N_SURROGATES, SWAP_FACTOR, BASE_SEED,
)
from lrg_eegfc.utils.fc.backbone import percolation_backbone, mst_union_top_fraction
from lrg_eegfc.utils.fc.heat_multiscale import (
    laplacian_eig, scale_grid, rho_sym_over_scales,
)
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle

PHASES = ("A", "B", "task_test", "rest_post")
N_S = 16                      # scale-grid resolution
R = N_SURROGATES

# backbone selection (percolation = primary; fixed densities = robustness cross-check)
BACKBONES = {
    "percolation": lambda W: percolation_backbone(W)[0],
    "d0.15": lambda W: mst_union_top_fraction(W, 0.15),
    "d0.20": lambda W: mst_union_top_fraction(W, 0.20),
}
BACKBONE_FN = BACKBONES["percolation"]
OUT = ROOT / "data" / "sparsified_arc" / "trace_arc"


def _eig_backbone(W):
    return laplacian_eig(BACKBONE_FN(W))


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception:
        return None
    N = Ws["A"].shape[0]
    if any(v.shape[0] != N for v in Ws.values()):
        return None
    eig = {ph: _eig_backbone(Ws[ph]) for ph in PHASES}
    s_grid = scale_grid(eig["A"][0], n=N_S)          # per-cell, from phase-A spectrum
    obs = rho_sym_over_scales(eig, s_grid)

    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    rng = np.random.default_rng(BASE_SEED + idx)     # same seed/order as audit_150
    surr = np.full((R, N_S), np.nan)
    surr_max = np.full(R, np.nan)
    for r in range(R):
        eig_s = {}
        for ph in PHASES:
            Wsurr = matched_strength_shuffle(Ws[ph], n_swaps, rng)
            eig_s[ph] = _eig_backbone(Wsurr)
        srow = rho_sym_over_scales(eig_s, s_grid)
        surr[r] = srow
        surr_max[r] = np.nanmax(srow) if np.isfinite(srow).any() else np.nan

    p_of_s = np.array([np.mean(surr[np.isfinite(surr[:, j]), j] >= obs[j])
                       if np.isfinite(obs[j]) and np.isfinite(surr[:, j]).any() else np.nan
                       for j in range(N_S)])
    surr_p50 = np.nanquantile(surr, 0.50, axis=0)
    surr_p95 = np.nanquantile(surr, 0.95, axis=0)

    obs_max = float(np.nanmax(obs)); i_arg = int(np.nanargmax(obs))
    sm = surr_max[np.isfinite(surr_max)]
    p_scalemax = float(np.mean(sm >= obs_max)) if sm.size else np.nan
    j1 = int(np.argmin(np.abs(s_grid - 1.0)))        # s=1 anchor

    cell = OUT / band
    cell.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cell / f"{pat}.npz", s=s_grid, obs=obs,
                        surr_p50=surr_p50, surr_p95=surr_p95, p_of_s=p_of_s,
                        surr_max_p50=float(np.nanquantile(sm, 0.5)) if sm.size else np.nan)
    return dict(patient=pat, band=band, N_nodes=int(N),
                obs_at_s1=float(obs[j1]), p_at_s1=float(p_of_s[j1]),
                obs_max=obs_max, argmax_s=float(s_grid[i_arg]),
                surr_max_p50=float(np.nanquantile(sm, 0.5)) if sm.size else np.nan,
                p_scalemax=p_scalemax,
                s_grid=";".join(f"{x:.3f}" for x in s_grid),
                obs_s=";".join(f"{x:.4f}" for x in obs),
                surr_p50_s=";".join(f"{x:.4f}" for x in surr_p50))


def cohort_analysis(df):
    """Per band: per-s Wilcoxon gate, scale-max gate, and nature classification."""
    vs_s, nature = [], []
    for band in BANDS:
        x = df[df.band == band]
        if x.empty:
            continue
        S = np.array([[float(v) for v in s.split(";")] for s in x.s_grid])
        OBS = np.array([[float(v) for v in s.split(";")] for s in x.obs_s])
        P50 = np.array([[float(v) for v in s.split(";")] for s in x.surr_p50_s])
        s_axis = np.nanmedian(S, axis=0)
        gate_p = np.full(N_S, np.nan)
        for j in range(N_S):
            d = OBS[:, j] - P50[:, j]
            d = d[np.isfinite(d)]
            if d.size >= 5 and np.any(d != 0):
                try:
                    _, gate_p[j] = wilcoxon(d, alternative="greater")
                except Exception:
                    pass
            vs_s.append(dict(band=band, s=float(s_axis[j]), gate_p=float(gate_p[j]),
                             obs_median=float(np.nanmedian(OBS[:, j])),
                             surr_median=float(np.nanmedian(P50[:, j]))))
        # scale-max cohort gate
        try:
            _, p_sm = wilcoxon(x.obs_max.values - x.surr_max_p50.values, alternative="greater")
        except Exception:
            p_sm = np.nan
        # s=1 cohort gate (anchor)
        try:
            j1 = int(np.argmin(np.abs(s_axis - 1.0)))
            _, p_s1 = wilcoxon(OBS[:, j1] - P50[:, j1], alternative="greater")
        except Exception:
            p_s1 = np.nan
        fire = np.isfinite(gate_p) & (gate_p < 0.05)
        frac_fire = float(fire.mean())
        # contiguous bands of firing scales
        n_bands = 0; prev = False
        for f in fire:
            if f and not prev:
                n_bands += 1
            prev = f
        if frac_fire == 0:
            cls = "none"
        elif frac_fire >= 2 / 3:
            cls = "continuous-multiscale"
        elif n_bands >= 2:
            cls = "multiscale"
        else:
            cls = "single-scale"
        fire_s = s_axis[fire]
        nature.append(dict(band=band, n_patients=len(x), classification=cls,
                           frac_scales_firing=round(frac_fire, 3),
                           n_firing_bands=n_bands,
                           fire_s_min=float(fire_s.min()) if fire_s.size else np.nan,
                           fire_s_max=float(fire_s.max()) if fire_s.size else np.nan,
                           gate_p_at_s1=float(p_s1), gate_p_scalemax=float(p_sm),
                           obs_max_median=float(x.obs_max.median()),
                           argmax_s_median=float(x.argmax_s.median())))
    return pd.DataFrame(vs_s), pd.DataFrame(nature)


def main():
    global R, BACKBONE_FN, OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=N_SURROGATES)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--backbone", choices=list(BACKBONES), default="percolation")
    a = ap.parse_args(); R = a.R
    BACKBONE_FN = BACKBONES[a.backbone]
    OUT = ROOT / "data" / "sparsified_arc" / (
        "trace_arc" if a.backbone == "percolation" else f"trace_arc_{a.backbone.replace('.', '')}")
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"[D2-trace] backbone={a.backbone} -> {OUT.name}", flush=True)
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))  # warm numba
    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in BANDS for p in COHORT)]
    if a.limit:
        jobs = jobs[:a.limit]
    ncpu = min(14, (os.cpu_count() or 4) - 2)
    print(f"[D2-trace] {len(jobs)} cells, R={R}, {N_S} scales, percolation backbone, {ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, r in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if r:
                rows.append(r)
                el = time.time() - t0
                print(f"[{i}/{len(jobs)}] {r['patient']}/{r['band']} "
                      f"s1={r['obs_at_s1']:+.3f}(p{r['p_at_s1']:.2f}) "
                      f"max={r['obs_max']:+.3f}@s{r['argmax_s']:.1f}(p{r['p_scalemax']:.3f}) "
                      f"[{el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s]", flush=True)
    df = pd.DataFrame(rows).sort_values(["band", "patient"])
    df.to_csv(OUT / "per_cell.csv", index=False)
    vs_s, nature = cohort_analysis(df)
    vs_s.to_csv(OUT / "cohort_gate_vs_s.csv", index=False)
    nature.to_csv(OUT / "band_nature.csv", index=False)
    rt = time.time() - t0
    print(f"\n[D2-trace] {len(rows)} cells in {rt:.0f}s", flush=True)
    print("\n=== BAND TRACE NATURE (percolation backbone, matched-strength, tau-resolved) ===", flush=True)
    print(nature[["band", "classification", "frac_scales_firing", "gate_p_at_s1",
                  "gate_p_scalemax", "argmax_s_median", "obs_max_median"]].to_string(index=False), flush=True)
    if not a.limit:
        (OUT / "config.json").write_text(json.dumps(dict(
            deliverable="D2_trace_arc", backbone="percolation (theta*)",
            estimator="rho_sym (symmetric split-half cophenetic)",
            null="matched-strength 4-cycle +/-delta, percolation-sparsified same way",
            R=R, n_scales=N_S, scale="s=tau*lambda_max in [1, s(Shat=0.02)] capped 200",
            cohort=COHORT, bands=BANDS, seed=BASE_SEED, swap_factor=SWAP_FACTOR,
            anchor="s=1 reproduces audit_150 gate up to percolation vs full graph",
            note="drift null is a SEPARATE secondary pass (pending); matched-strength is the mandatory gate.",
        ), indent=2))
    print(f"[D2-trace] outputs -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
