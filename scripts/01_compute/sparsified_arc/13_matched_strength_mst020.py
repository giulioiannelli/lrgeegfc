#!/usr/bin/env python3
"""Matched-strength null for the mst@0.20 recovered trace, tau-resolved, per patient.

The ONLY meaningful null (drift discarded: task is directional -> a trace is
monotonic by construction). Pipeline-identical to the observed measure:
  matched-strength shuffle the DENSE imcoh_abs FC (audit_150's exact bit-identical
  4-cycle +/-delta, strength-preserving) -> sparsify @ frac=0.20 (mst_union_top_
  fraction) -> LRG cophenetic at scale s -> rho_sym. Swept over s = tau*lambda_max.
Null H0: the recovered cophenetic trace is explained by node strength alone (same
sparsify+LRG pipeline on a strength-matched graph).

Reports EVERYTHING per scale (no scalar collapse; feedback_report_tau_dependence):
  - per (patient, band, s): obs_rho, surr p50/p95, one-sided p = mean(surr >= obs)
  - cohort gate_p(s) = Wilcoxon(obs - surr_p50, greater) over 10 patients
  - per-patient: which patients beat their own null, in which bands, over how many scales

Outputs (data/sparsified_arc/ms_mst020/):
  per_patient_scale.csv  : patient, band, s, obs_rho, surr_p50, surr_p95, p, N
  cohort_gate.csv        : band, s, gate_p, n_above, obs_med, surr_med
  config.json
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS, load_phase
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, rho_sym_over_scales
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle

# Backbone via env: mst020 (recovery baseline) | tmfg | pmfg | disparity.
# SA_FRAC sets the mst-union fraction (single-fraction unification test: 0.05/0.10/0.20).
BACKBONE = os.environ.get("SA_BACKBONE", "mst020")
DISP_ALPHA = float(os.environ.get("SA_DISP_ALPHA", "0.20"))
FRAC = float(os.environ.get("SA_FRAC", "0.20"))
if BACKBONE == "mst020":
    _OUTNAME = f"ms_mst{int(round(FRAC * 100)):03d}"          # ms_mst020 / ms_mst010 / ms_mst005
elif BACKBONE == "disparity":
    _OUTNAME = f"ms_disparity_a{DISP_ALPHA:g}"
else:
    _OUTNAME = f"ms_{BACKBONE}"
OUT = ROOT / "data" / "sparsified_arc" / _OUTNAME
PHASES = ("A", "B", "task_test", "rest_post")
SGRID = np.logspace(0.0, np.log10(180.0), 16)
R = 200
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260712


def eig_mst(W):
    return laplacian_eig(select_backbone(W, BACKBONE, frac=FRAC, disparity_alpha=DISP_ALPHA))


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception:
        return None
    N = Ws["A"].shape[0]
    eig_obs = {ph: eig_mst(Ws[ph]) for ph in PHASES}
    obs = rho_sym_over_scales(eig_obs, SGRID)                       # (nS,)
    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    surr = np.full((R, SGRID.size), np.nan)
    for r in range(R):
        eig_s = {ph: eig_mst(matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX))
                 for ph in PHASES}
        surr[r] = rho_sym_over_scales(eig_s, SGRID)
    rows = []
    for j, s in enumerate(SGRID):
        col = surr[:, j]; col = col[np.isfinite(col)]
        o = obs[j]
        p = float(np.mean(col >= o)) if col.size and np.isfinite(o) else np.nan
        rows.append(dict(patient=pat, band=band, s=float(s), N=int(N),
                         obs_rho=float(o) if np.isfinite(o) else np.nan,
                         surr_p50=float(np.nanpercentile(col, 50)) if col.size else np.nan,
                         surr_p95=float(np.nanpercentile(col, 95)) if col.size else np.nan,
                         p=p, n_surr=int(col.size)))
    return rows


def cohort_gate(df):
    out = []
    for band in BANDS:
        for s in np.sort(df.s.unique()):
            x = df[(df.band == band) & np.isclose(df.s, s)].dropna(subset=["obs_rho", "surr_p50"])
            if len(x) < 5:
                continue
            d = x.obs_rho.values - x.surr_p50.values
            try:
                p = float(wilcoxon(d, alternative="greater")[1]) if np.any(d != 0) else np.nan
            except Exception:
                p = np.nan
            out.append(dict(band=band, s=float(s), gate_p=p,
                            n_above=int((x.p < 0.05).sum()), n_pat=len(x),
                            obs_med=float(x.obs_rho.median()),
                            surr_med=float(x.surr_p50.median())))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--bands", type=str, default="")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))  # warm numba
    pats = COHORT[:a.limit] if a.limit else COHORT
    bands = a.bands.split(",") if a.bands else BANDS
    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in bands for p in pats)]
    ncpu = int(os.environ.get("SA_WORKERS", 12))
    print(f"[ms-{BACKBONE}] {len(jobs)} cells, R={R}, backbone={BACKBONE}"
          f"{f'(a={DISP_ALPHA:g})' if BACKBONE=='disparity' else ''}, "
          f"{len(SGRID)} scales, {ncpu} workers -> {OUT.name}", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, rl in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if rl:
                rows.extend(rl)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s", flush=True)
    df = pd.DataFrame(rows)
    if df.empty:
        print("[ms-mst020] no cells", flush=True); return
    df.to_csv(OUT / "per_patient_scale.csv", index=False)
    gate = cohort_gate(df); gate.to_csv(OUT / "cohort_gate.csv", index=False)
    print(f"\n[ms-mst020] {len(df)} rows in {time.time()-t0:.0f}s -> {OUT}\n", flush=True)
    if gate.empty or "band" not in gate.columns:
        print("[ms-mst020] cohort gate needs >=5 patients (timing/subset run)", flush=True)
        return
    # best (min gate_p) scale per band
    print("=== mst@0.20 matched-strength gate: best scale per band ===", flush=True)
    for band in bands:
        g = gate[gate.band == band].dropna(subset=["gate_p"])
        if g.empty:
            continue
        b = g.loc[g.gate_p.idxmin()]
        star = "CLEAR" if b.gate_p < 0.05 else "fail"
        print(f"  {band:11s} s*={b.s:6.1f}  gate_p={b.gate_p:.4f}  "
              f"n_above={int(b.n_above)}/{int(b.n_pat)}  obs_med={b.obs_med:+.3f}  "
              f"surr_med={b.surr_med:+.3f}  [{star}]", flush=True)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="matched_strength_trace", backbone=BACKBONE,
        disparity_alpha=DISP_ALPHA if BACKBONE == "disparity" else None,
        frac=FRAC, R=R, swap_factor=SWAP_FACTOR,
        s_grid=[float(x) for x in SGRID], cohort=pats, bands=bands,
        null=f"matched-strength on dense FC -> sparsify({BACKBONE}) -> LRG -> rho_sym, per scale",
        base_seed=BASE_SEED), indent=2))
    print(f"\n[ms-mst020] done -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
