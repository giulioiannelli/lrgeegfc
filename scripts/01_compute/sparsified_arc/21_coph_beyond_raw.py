#!/usr/bin/env python3
"""Non-circular incremental test: does the cophenetic trace carry persistence BEYOND raw FC?

The 'scale-invariance is only-multiscale' claim is partly self-referential (a scale axis is a
multiscale construct). This asks the information question instead: residualize each phase's
cophenetic geometry on its raw FC edges (coph ⊥ raw = the part of the hierarchy NOT linearly
in the pairwise edges — the higher-order / multi-path structure), and test whether THAT still
traces (rest_post re-approaches task). Then the reverse (raw ⊥ coph). Asymmetry is the result:
  coph|raw traces  AND  raw|coph does not  ⇒ the hierarchy SUBSUMES the pairwise trace and adds
  persistent higher-order structure the edges lack — raw FC is insufficient, non-circularly.

Per (patient, band), mesoscale s=5.6, mst@0.20 backbone:
  A_x   = dense |ImCoh| upper-tri (raw pairwise edges), phase x in {A,B,task_test,rest_post}
  C_x   = cophenetic UPGMA distances at s=5.6
  cr_x  = rank(C_x) residualized on rank(A_x)   (cophenetic ⊥ raw)
  ar_x  = rank(A_x) residualized on rank(C_x)   (raw ⊥ cophenetic)
  trace_sym(v) = 1/2[ spearman(v_t - v_A, v_p - v_B) + spearman(v_t - v_B, v_p - v_A) ]
  T_coph|raw = trace_sym(cr) ;  T_raw|coph = trace_sym(ar)
Matched-strength null (shuffle each dense FC -> recompute A,C,residuals,traces), R=200.
Cohort Wilcoxon(obs - surr_p50, greater) per band + BH.

5-point preamble
1. Claim: the cophenetic trace persists after removing the raw-FC-explainable part (coph|raw>0),
   while the raw trace does NOT persist after removing the cophenetic part (raw|coph~0).
2. Null: the residual 'trace' is a node-strength artifact of the shuffle-then-residualize pipeline.
3. Strongest alt it must beat: node strength -> IDENTICAL matched-strength surrogate through the
   SAME residualization pipeline (so any residualization bias is in the null too).
4. Cannot: residualization is rank-linear (removes monotone raw dependence, not arbitrary
   nonlinear); C_x is a deterministic function of A_x, so coph|raw = the higher-order transform's
   own persistent content (that is the point, not a confound). n=10 Wilcoxon.
5. Falsify: if coph|raw does not clear, the hierarchy adds no persistent structure beyond raw FC;
   if raw|coph ALSO clears, they carry complementary (not subsumed) trace -> weaker claim.

Output: data/sparsified_arc/coph_beyond_raw/{per_patient,cohort}.csv
"""
from __future__ import annotations
import os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import rankdata, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, load_phase            # noqa: E402
from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction    # noqa: E402
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, cophenetic_at_scale  # noqa: E402
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle  # noqa: E402
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr             # noqa: E402

PHASES = ("A", "B", "task_test", "rest_post")
BANDS = ["delta", "alpha", "beta", "low_gamma"]
FRAC = 0.20
S_MESO = float(os.environ.get("SA_SCALE", "5.6"))
R = 200
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260713
LIMIT = int(os.environ.get("SA_LIMIT", "0"))


def _cr(x):
    r = rankdata(x); return r - r.mean()


def resid(x, c):
    """rank(x) residualized on rank(c): remove the monotone dependence of x on c."""
    rx, rc = _cr(x), _cr(c)
    rc2 = float(np.dot(rc, rc))
    return rx - (np.dot(rx, rc) / rc2) * rc if rc2 > 0 else rx


def _sp(a, b):
    ra, rb = rankdata(a), rankdata(b)
    ra = ra - ra.mean(); rb = rb - rb.mean()
    da, db = np.sqrt(np.dot(ra, ra)), np.sqrt(np.dot(rb, rb))
    return float(np.dot(ra, rb) / (da * db)) if da > 0 and db > 0 else np.nan


def trace_sym(vA, vB, vt, vp):
    return 0.5 * (_sp(vt - vA, vp - vB) + _sp(vt - vB, vp - vA))


def coph_condensed(W):
    ev, V = laplacian_eig(mst_union_top_fraction(W, FRAC))
    return cophenetic_at_scale(ev, V, S_MESO)


def two_traces(Ws):
    """Return (T_coph|raw, T_raw|coph) for a dict of 4 phase FC matrices."""
    N = Ws["A"].shape[0]
    iu = np.triu_indices(N, 1)
    A = {ph: Ws[ph][iu] for ph in PHASES}          # raw edges (condensed)
    C = {ph: coph_condensed(Ws[ph]) for ph in PHASES}
    cr = {ph: resid(C[ph], A[ph]) for ph in PHASES}   # coph ⊥ raw
    ar = {ph: resid(A[ph], C[ph]) for ph in PHASES}   # raw ⊥ coph
    t_cr = trace_sym(cr["A"], cr["B"], cr["task_test"], cr["rest_post"])
    t_ar = trace_sym(ar["A"], ar["B"], ar["task_test"], ar["rest_post"])
    return t_cr, t_ar


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception:
        return None
    N = Ws["A"].shape[0]
    if any(v.shape[0] != N for v in Ws.values()):
        return None
    o_cr, o_ar = two_traces(Ws)
    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    s_cr = np.full(R, np.nan); s_ar = np.full(R, np.nan)
    for r in range(R):
        Wsh = {ph: matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX) for ph in PHASES}
        try:
            s_cr[r], s_ar[r] = two_traces(Wsh)
        except Exception:
            continue

    def rec(name, obs, surr):
        col = surr[np.isfinite(surr)]
        p = float((np.sum(col >= obs) + 1) / (col.size + 1)) if col.size else np.nan
        return dict(patient=pat, band=band, measure=name, obs=float(obs),
                    surr_p50=float(np.nanmedian(col)) if col.size else np.nan, p=p)
    return [rec("coph|raw", o_cr, s_cr), rec("raw|coph", o_ar, s_ar)]


def cohort(df):
    out = []
    for band in BANDS:
        for meas in ("coph|raw", "raw|coph"):
            x = df[(df.band == band) & (df.measure == meas)].dropna(subset=["obs", "surr_p50"])
            if len(x) < 5:
                continue
            diff = x.obs.values - x.surr_p50.values
            try:
                _, p = wilcoxon(diff, alternative="greater")
            except ValueError:
                p = 1.0
            out.append(dict(band=band, measure=meas, n=len(x), obs_med=float(x.obs.median()),
                            surr_med=float(x.surr_p50.median()), n_pos=int((diff > 0).sum()),
                            p=float(p)))
    c = pd.DataFrame(out)
    if not c.empty:
        for meas in ("coph|raw", "raw|coph"):
            m = c.measure == meas
            c.loc[m, "q"] = bh_fdr(c.loc[m, "p"].values)
    return c


def main():
    OUT = ROOT / "data" / "sparsified_arc" / f"coph_beyond_raw_s{S_MESO:04.1f}"
    OUT.mkdir(parents=True, exist_ok=True)
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    jobs = [(i, p, b) for i, (p, b) in enumerate((p, b) for b in BANDS for p in COHORT)]
    if LIMIT:
        jobs = jobs[:LIMIT]
    print(f"[coph-beyond-raw] s={S_MESO} R={R} {len(jobs)} cells", flush=True)
    t0 = time.time(); recs = []
    nproc = min(12, max(1, (os.cpu_count() or 2) - 2))
    with Pool(nproc) as pool:
        for k, rl in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if rl:
                recs.extend(rl)
            el = time.time() - t0
            print(f"  [{k}/{len(jobs)}] {el:.0f}s ETA {el/k*(len(jobs)-k):.0f}s", flush=True)
    per_pat = pd.DataFrame(recs); per_pat.to_csv(OUT / "per_patient.csv", index=False)
    c = cohort(per_pat); c.to_csv(OUT / "cohort.csv", index=False)
    print("\n=== INCREMENTAL TRACE (upper-tail matched-strength; coph⊥raw vs raw⊥coph) ===", flush=True)
    if c.empty:
        print("  [cohort empty — need >=5 patients]", flush=True)
        print(f"\n[done] {time.time()-t0:.0f}s -> {OUT}", flush=True)
        return
    for band in BANDS:
        b = c[c.band == band]
        if b.empty:
            continue
        cr = b[b.measure == "coph|raw"]; ar = b[b.measure == "raw|coph"]
        cs = f"coph|raw p={cr.p.iloc[0]:.3f} q={cr.q.iloc[0]:.3f} ({cr.n_pos.iloc[0]}/{cr.n.iloc[0]}, obs {cr.obs_med.iloc[0]:+.3f})" if not cr.empty else ""
        as_ = f"raw|coph p={ar.p.iloc[0]:.3f} q={ar.q.iloc[0]:.3f} ({ar.n_pos.iloc[0]}/{ar.n.iloc[0]}, obs {ar.obs_med.iloc[0]:+.3f})" if not ar.empty else ""
        print(f"  {band:10s} {cs:58s} | {as_}", flush=True)
    print(f"\n[done] {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
