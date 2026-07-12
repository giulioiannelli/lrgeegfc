#!/usr/bin/env python3
"""Minimal-modification sparsification recovery of the rho_sym band-selective trace.

QUESTION (user, 2026-07-12): the old FULLY-CONNECTED scheme (audit_150) gave a
band-selective cophenetic trace (alpha, beta CLEAR; others fail) but its LRG
propagator is DEGENERATE at tau_min (D=1/K ~ D=1/A, audit_174). Sparsification
makes the propagator non-trivial/multiscale (Villegas 2025), but percolation broke
band-selectivity. So: starting from the OLD scheme, which sparsification -- at some
scale s = tau*lambda_max -- best REPRODUCES the old alpha/beta selectivity, with the
result depending as LITTLE as possible on the choice of sparsifier?

This is the minimal-modification sweep. It reproduces audit_150 EXACTLY at
(method=dense, s=1): dense = mst_union_top_fraction(W, 1.0) = the cleaned full graph,
and cophenetic_at_scale(.,.,1.0) = tau=1/lambda_max = audit_150 `ultra`. Then it
sweeps sparsifier x density x scale and reports rho_sym per (patient, band).

NO drift null (invalid: the task is directional -> a trace is monotonic by
construction, so removing monotonic structure removes the signal). The ONLY
meaningful null is matched-strength -- run SEPARATELY on the shortlist (Stage 2);
this Stage-1 sweep is EXPLORATORY triage to locate the (method, scale) sweet spot.
Read per-scale curves, not a scalar (feedback_report_tau_dependence_no_scalar_collapse).

Sparsifiers (all library, all connected+spanning so per-pair rho_sym is preserved):
  dense (frac=1.0=old), mst_union_top_fraction @ frac in {0.7..0.05},
  percolation_backbone (parameter-free, current), tmfg_backbone (parameter-free).

Outputs (data/sparsified_arc/sparsification_recovery/):
  rho_sym_sweep.csv   : patient, band, method, frac, density, s, rho_sym
  multiscale.csv      : patient, band, method, density, n_peaks, dominant_s  (task_test)
  cohort_best.csv     : band x method -> best-scale cohort-median rho_sym + s
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS, load_phase   # exact old-scheme loader
from lrg_eegfc.utils.fc.backbone import (
    mst_union_top_fraction, percolation_backbone, tmfg_backbone, backbone_density,
)
from lrg_eegfc.utils.fc.heat_multiscale import (
    laplacian_eig, cophenetic_at_scale, rho_sym,
    entropy_specific_heat, specific_heat_peaks,
)

OUT = ROOT / "data" / "sparsified_arc" / "sparsification_recovery"
PHASES = ("A", "B", "task_test", "rest_post")
FRACS = [1.00, 0.70, 0.50, 0.40, 0.30, 0.20, 0.15, 0.10, 0.05]
SGRID = np.logspace(0.0, np.log10(180.0), 16)     # s = tau*lambda_max, fixed & comparable


def backbones(W):
    """name -> (sparsified W, nominal frac). Connected + spanning."""
    out = {}
    for f in FRACS:
        out[f"mst{f:.2f}"] = (mst_union_top_fraction(W, f), f)
    out["perc"] = (percolation_backbone(W)[0], np.nan)
    out["tmfg"] = (tmfg_backbone(W), np.nan)
    return out


def per_cell(job):
    pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception:
        return [], []
    N = Ws["A"].shape[0]
    rows, mrows = [], []
    for meth, (_, frac) in backbones(Ws["A"]).items():
        # per-phase backbone (spanning -> all N nodes kept -> pair vectors aligned)
        Bs, eig, dens = {}, {}, []
        ok = True
        for ph in PHASES:
            W = mst_union_top_fraction(Ws[ph], frac) if meth.startswith("mst") else (
                percolation_backbone(Ws[ph])[0] if meth == "perc" else tmfg_backbone(Ws[ph]))
            Bs[ph] = W
            dens.append(backbone_density(W))
            eig[ph] = laplacian_eig(W)
        density = float(np.mean(dens))
        for s in SGRID:
            try:
                D = {ph: cophenetic_at_scale(*eig[ph], float(s)) for ph in PHASES}
                if len({d.size for d in D.values()}) != 1:
                    raise ValueError("size mismatch")
                r, _ = rho_sym(D["A"], D["B"], D["task_test"], D["rest_post"])
            except Exception:
                r = np.nan
            rows.append(dict(patient=pat, band=band, method=meth, frac=frac,
                             density=density, s=float(s), rho_sym=float(r)))
        # multiscale indicator on task_test (graph property)
        try:
            res = entropy_specific_heat(eig["task_test"][0])
            pk = specific_heat_peaks(res)
            mrows.append(dict(patient=pat, band=band, method=meth, frac=frac,
                              density=density, n_peaks=pk["n_peaks"],
                              dominant_s=pk["dominant_s"], N=int(N)))
        except Exception:
            pass
    return rows, mrows


def cohort_best(df):
    """Per (band, method): the scale with the largest cohort-median rho_sym."""
    out = []
    for band in BANDS:
        for meth in df.method.unique():
            x = df[(df.band == band) & (df.method == meth)]
            if x.empty:
                continue
            g = x.groupby("s").rho_sym.median()
            g = g[np.isfinite(g.values)]
            if g.empty:
                continue
            s_best = float(g.idxmax()); r_best = float(g.max())
            r_s1 = float(x[np.isclose(x.s, 1.0)].rho_sym.median())
            out.append(dict(band=band, method=meth,
                            density=float(x.density.median()),
                            rho_s1=r_s1, rho_best=r_best, s_best=s_best,
                            n_pat=int(x.patient.nunique())))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="first K patients (timing)")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    pats = COHORT[:a.limit] if a.limit else COHORT
    jobs = [(p, b) for b in BANDS for p in pats]
    ncpu = int(os.environ.get("SA_WORKERS", 6))
    print(f"[recovery] {len(jobs)} (pt,band) cells | {len(FRACS)}+2 sparsifiers x "
          f"{len(SGRID)} scales | dense@s=1 == audit_150 | {ncpu} workers", flush=True)
    t0 = time.time(); rows, mrows = [], []
    with Pool(ncpu) as pool:
        for i, (rl, ml) in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            rows.extend(rl); mrows.extend(ml)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s", flush=True)
    df = pd.DataFrame(rows); md = pd.DataFrame(mrows)
    if df.empty:
        print("[recovery] no cells", flush=True); return
    df.to_csv(OUT / "rho_sym_sweep.csv", index=False)
    md.to_csv(OUT / "multiscale.csv", index=False)
    cb = cohort_best(df); cb.to_csv(OUT / "cohort_best.csv", index=False)

    print(f"\n[recovery] {len(df)} rows in {time.time()-t0:.0f}s -> {OUT}\n", flush=True)
    # --- audit_150 reproduction check: dense @ s=1 cohort-median rho_sym per band ---
    d1 = df[(df.method == "mst1.00") & np.isclose(df.s, 1.0)]
    print("=== DENSE @ s=1  (must reproduce audit_150 alpha/beta selectivity) ===", flush=True)
    print(d1.groupby("band").rho_sym.median().reindex(BANDS).round(3).to_string(), flush=True)
    # --- best-scale cohort-median rho_sym: band x method ---
    piv = cb.pivot(index="band", columns="method", values="rho_best").reindex(BANDS)
    cols = ["mst1.00", "mst0.50", "mst0.30", "mst0.20", "mst0.10", "perc", "tmfg"]
    cols = [c for c in cols if c in piv.columns]
    print("\n=== BEST-SCALE cohort-median rho_sym  (band x sparsifier) ===", flush=True)
    print(piv[cols].round(3).to_string(), flush=True)
    print("\n=== multiscale: mean C(tau) peak count per sparsifier (task_test) ===", flush=True)
    print(md.groupby("method").agg(n_peaks=("n_peaks", "mean"),
                                    density=("density", "mean")).round(2).to_string(), flush=True)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="sparsification_recovery", phases=list(PHASES), fracs=FRACS,
        s_grid=[float(x) for x in SGRID], cohort=pats, bands=BANDS,
        note="Stage-1 exploratory triage; matched-strength null on shortlist = Stage-2.",
        dense_equals="audit_150 ultra (frac=1.0, s=1)"), indent=2))
    print(f"\n[recovery] done -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
