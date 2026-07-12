#!/usr/bin/env python3
"""D3 -- tau-resolved encoding vs inference arc on the percolation backbone.

Adds the encoding phase (task_learn) to the trace arc and computes, at every
diffusion scale s = tau*lambda_max, the four cophenetic cross-phase functionals
(symmetric split-half over the A/B rest_pre halves):

  D_x(s) = cophenetic UPGMA distances of phase x at scale s, x in
           {A, B, task_learn, task_test, rest_post}
  e = D_learn - D_A     (encoding;  e2 = D_learn - D_B)
  f = D_test  - D_learn  (inference-specific; arm-invariant)
  g = D_test  - D_A     (standard;   g2 = D_test - D_B)
  p = D_post  - D_B     (persistence; p2 = D_post - D_A)

  T_test(s)      = 1/2[ rho(g,p)  + rho(g2,p2) ]      standard trace (== D2)
  T_learn(s)     = 1/2[ rho(e,p2) + rho(e2,p)  ]      encoding echo
  T_infspec(s)   = 1/2[ rho(f,p)  + rho(f,p2)  ]      inference-specific
  T_infspec_pe(s)= 1/2[ pr(f,p|e) + pr(f,p2|e2) ]     inference-specific | encoding
                                                       (the beta-only headline)

Matched-strength null over ALL FIVE phases (shuffled with the audit_150 rng/order,
percolation-sparsified the same way), recomputed at every scale. Scale-MAX gate +
per-scale cohort Wilcoxon, per functional per band. This is the tau-resolved
generalisation of audit_152 / audit_103.

Outputs:
    data/sparsified_arc/enc_inf_arc/{band}/{patient}.npz   (s + 4 functionals x obs/surr)
    data/sparsified_arc/enc_inf_arc/per_cell.csv
    data/sparsified_arc/enc_inf_arc/cohort_gate_vs_s.csv
    data/sparsified_arc/enc_inf_arc/band_nature.csv
    data/sparsified_arc/enc_inf_arc/config.json
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import (
    load_phase, COHORT, BANDS, N_SURROGATES, SWAP_FACTOR, BASE_SEED,
)
from lrg_eegfc.utils.fc.backbone import percolation_backbone, mst_union_top_fraction
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, scale_grid, cophenetic_at_scale
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle

PHASES5 = ("A", "B", "task_learn", "task_test", "rest_post")
FUNCS = ("T_test", "T_learn", "T_infspec", "T_infspec_pe")
# Backbone selected by env: BACKBONE=mst020 (recovered scheme) | perc (legacy). mst@0.20
# recovered clean alpha/beta band-selectivity under matched-strength (13_matched_strength_mst020).
BACKBONE = os.environ.get("SA_BACKBONE", "mst020")
FRAC = 0.20
OUT = ROOT / "data" / "sparsified_arc" / ("enc_inf_arc_mst020" if BACKBONE == "mst020" else "enc_inf_arc")
N_S = 16
R = N_SURROGATES


def _rho(a, b):
    r, _ = spearmanr(a, b)
    return r


def _partial(a, b, c):
    rab, rac, rbc = _rho(a, b), _rho(a, c), _rho(b, c)
    den = np.sqrt(max(0.0, (1 - rac ** 2) * (1 - rbc ** 2)))
    return (rab - rac * rbc) / den if den > 0 else np.nan


def functionals_over_scales(eig5, s_grid):
    """Return dict func -> array(len s_grid) of the 4 symmetric functionals."""
    out = {k: np.full(len(s_grid), np.nan) for k in FUNCS}
    for i, s in enumerate(s_grid):
        try:
            D = {ph: cophenetic_at_scale(*eig5[ph], s) for ph in PHASES5}
        except Exception:
            continue
        e = D["task_learn"] - D["A"]; e2 = D["task_learn"] - D["B"]
        f = D["task_test"] - D["task_learn"]
        g = D["task_test"] - D["A"]; g2 = D["task_test"] - D["B"]
        p = D["rest_post"] - D["B"]; p2 = D["rest_post"] - D["A"]
        out["T_test"][i] = 0.5 * (_rho(g, p) + _rho(g2, p2))
        out["T_learn"][i] = 0.5 * (_rho(e, p2) + _rho(e2, p))
        out["T_infspec"][i] = 0.5 * (_rho(f, p) + _rho(f, p2))
        out["T_infspec_pe"][i] = 0.5 * (_partial(f, p, e) + _partial(f, p2, e2))
    return out


def _eig(W):
    B = mst_union_top_fraction(W, FRAC) if BACKBONE == "mst020" else percolation_backbone(W)[0]
    return laplacian_eig(B)


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES5}
    except Exception:
        return None
    N = Ws["A"].shape[0]
    if any(v.shape[0] != N for v in Ws.values()):
        return None
    eig5 = {ph: _eig(Ws[ph]) for ph in PHASES5}
    s_grid = scale_grid(eig5["A"][0], n=N_S)
    obs = functionals_over_scales(eig5, s_grid)

    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    rng = np.random.default_rng(BASE_SEED + idx)
    surr = {k: np.full((R, N_S), np.nan) for k in FUNCS}
    for r in range(R):
        eig_s = {ph: _eig(matched_strength_shuffle(Ws[ph], n_swaps, rng)) for ph in PHASES5}
        so = functionals_over_scales(eig_s, s_grid)
        for k in FUNCS:
            surr[k][r] = so[k]

    cell = OUT / band
    cell.mkdir(parents=True, exist_ok=True)
    save = {"s": s_grid}
    rec = {"patient": pat, "band": band, "N_nodes": int(N)}
    j1 = int(np.argmin(np.abs(s_grid - 1.0)))
    for k in FUNCS:
        o = obs[k]; S = surr[k]
        p_of_s = np.array([np.mean(S[np.isfinite(S[:, j]), j] >= o[j])
                           if np.isfinite(o[j]) and np.isfinite(S[:, j]).any() else np.nan
                           for j in range(N_S)])
        p50 = np.nanquantile(S, 0.50, axis=0)
        smax_obs = float(np.nanmax(o)); i_arg = int(np.nanargmax(o))
        srmax = np.nanmax(S, axis=1); srmax = srmax[np.isfinite(srmax)]
        save[f"{k}__obs"] = o
        save[f"{k}__surr_p50"] = p50
        save[f"{k}__surr_p95"] = np.nanquantile(S, 0.95, axis=0)
        save[f"{k}__p_of_s"] = p_of_s
        rec[f"{k}__obs_s1"] = float(o[j1]); rec[f"{k}__p_s1"] = float(p_of_s[j1])
        rec[f"{k}__obs_max"] = smax_obs; rec[f"{k}__argmax_s"] = float(s_grid[i_arg])
        rec[f"{k}__surr_max_p50"] = float(np.nanquantile(srmax, 0.5)) if srmax.size else np.nan
        rec[f"{k}__p_scalemax"] = float(np.mean(srmax >= smax_obs)) if srmax.size else np.nan
    np.savez_compressed(cell / f"{pat}.npz", **save)
    return rec


def cohort_analysis(df):
    vs_s, nature = [], []
    for band in BANDS:
        x = df[df.band == band]
        if x.empty:
            continue
        # load per-cell npz arrays
        cells = []
        for pat in x.patient:
            f = OUT / band / f"{pat}.npz"
            if f.exists():
                cells.append(np.load(f))
        if not cells:
            continue
        s_axis = cells[0]["s"]
        for k in FUNCS:
            OBS = np.array([c[f"{k}__obs"] for c in cells])
            P50 = np.array([c[f"{k}__surr_p50"] for c in cells])
            gate_p = np.full(N_S, np.nan)
            for j in range(N_S):
                d = OBS[:, j] - P50[:, j]; d = d[np.isfinite(d)]
                if d.size >= 5 and np.any(d != 0):
                    try:
                        _, gate_p[j] = wilcoxon(d, alternative="greater")
                    except Exception:
                        pass
                vs_s.append(dict(band=band, functional=k, s=float(s_axis[j]),
                                 gate_p=float(gate_p[j]),
                                 obs_median=float(np.nanmedian(OBS[:, j])),
                                 surr_median=float(np.nanmedian(P50[:, j]))))
            xo = x[f"{k}__obs_max"].values; xs = x[f"{k}__surr_max_p50"].values
            try:
                _, p_sm = wilcoxon(xo - xs, alternative="greater")
            except Exception:
                p_sm = np.nan
            j1 = int(np.argmin(np.abs(s_axis - 1.0)))
            try:
                _, p_s1 = wilcoxon(OBS[:, j1] - P50[:, j1], alternative="greater")
            except Exception:
                p_s1 = np.nan
            fire = np.isfinite(gate_p) & (gate_p < 0.05)
            frac = float(fire.mean())
            n_bands = 0; prev = False
            for fv in fire:
                if fv and not prev:
                    n_bands += 1
                prev = fv
            cls = ("none" if frac == 0 else "continuous-multiscale" if frac >= 2 / 3
                   else "multiscale" if n_bands >= 2 else "single-scale")
            fs = s_axis[fire]
            nature.append(dict(band=band, functional=k, classification=cls,
                               frac_scales_firing=round(frac, 3),
                               gate_p_at_s1=float(p_s1), gate_p_scalemax=float(p_sm),
                               argmax_s_median=float(x[f"{k}__argmax_s"].median()),
                               obs_max_median=float(x[f"{k}__obs_max"].median()),
                               n_pos=int((x[f"{k}__p_scalemax"] < 0.05).sum()),
                               fire_s_min=float(fs.min()) if fs.size else np.nan,
                               fire_s_max=float(fs.max()) if fs.size else np.nan))
    return pd.DataFrame(vs_s), pd.DataFrame(nature)


def main():
    global R
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=N_SURROGATES)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(); R = a.R
    OUT.mkdir(parents=True, exist_ok=True)
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in BANDS for p in COHORT)]
    if a.limit:
        jobs = jobs[:a.limit]
    ncpu = min(14, (os.cpu_count() or 4) - 2)
    print(f"[D3-encinf] {len(jobs)} cells, R={R}, {N_S} scales, 5 phases, {ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, r in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if r:
                rows.append(r)
                el = time.time() - t0
                print(f"[{i}/{len(jobs)}] {r['patient']}/{r['band']} "
                      f"Ttest_max={r['T_test__obs_max']:+.2f}(p{r['T_test__p_scalemax']:.2f}) "
                      f"Tinfpe_max={r['T_infspec_pe__obs_max']:+.2f}@s{r['T_infspec_pe__argmax_s']:.0f}"
                      f"(p{r['T_infspec_pe__p_scalemax']:.2f}) [{el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s]",
                      flush=True)
    df = pd.DataFrame(rows).sort_values(["band", "patient"])
    df.to_csv(OUT / "per_cell.csv", index=False)
    vs_s, nature = cohort_analysis(df)
    vs_s.to_csv(OUT / "cohort_gate_vs_s.csv", index=False)
    nature.to_csv(OUT / "band_nature.csv", index=False)
    print(f"\n[D3-encinf] {len(rows)} cells in {time.time()-t0:.0f}s", flush=True)
    print("\n=== ENC/INF FUNCTIONALS: scale-max gate per band (matched-strength) ===", flush=True)
    piv = nature.pivot(index="band", columns="functional", values="gate_p_scalemax")
    print(piv.reindex(BANDS)[list(FUNCS)].round(3).to_string(), flush=True)
    print("\n  T_infspec_pe detail (the beta-only headline):", flush=True)
    print(nature[nature.functional == "T_infspec_pe"][
        ["band", "classification", "gate_p_at_s1", "gate_p_scalemax",
         "argmax_s_median", "n_pos"]].to_string(index=False), flush=True)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="D3_enc_inf_arc", backbone="percolation (theta*)",
        phases=list(PHASES5), functionals=list(FUNCS), R=R, n_scales=N_S,
        cohort=COHORT, bands=BANDS, seed=BASE_SEED,
        headline="T_infspec_pe = partial rho(f,p|e); prior full-graph rho_sym p=0.0098 beta-only",
        note="drift null separate/secondary; matched-strength is the mandatory gate.",
    ), indent=2))
    print(f"[D3-encinf] outputs -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
