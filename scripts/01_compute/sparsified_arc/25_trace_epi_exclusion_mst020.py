#!/usr/bin/env python3
"""Does the mst@0.20 trace hold with the SEIZURE-ONSET ZONE physically removed?

The β trace must be a COGNITIVE result, not one riding on epileptic tissue. Stronger
than the pair-class restriction (script 19, which keeps SOZ nodes in the graph and reads
only non-SOZ pairs): here the SOZ contacts are DROPPED from the FC, the mst@0.20 backbone
is rebuilt on the reduced graph, and the ρ_sym trace is recomputed under matched-strength.
Pipeline-identical to the trace gate (script 13) except each phase's FC is first masked to
non-SOZ nodes.

Size-matched DECIMATION control (locked rule feedback_decimation_control_for_subset_exclusion):
for each patient we ALSO drop a RANDOM set of non-SOZ nodes of the SAME size as its SOZ set
(10 seeds) and recompute the observed trace — so "trace survives SOZ removal" is read against
"trace survives an equivalent random node removal", not against the full graph (which has
more nodes). If SOZ-removal preserves the trace as well as random removal, the trace does not
depend on the seizure zone.

5-point preamble
1. Claim: on mst@0.20 the β trace clears matched-strength with the SOZ removed (β is cognitive,
   SOZ-independent); low_γ does NOT (its trace needs the SOZ).
2. Null: matched-strength R=200 on the SOZ-excluded reduced graph (same 4-cycle ±δ → sparsify
   @0.20 → LRG).
3. Strongest alternative it must beat: (a) node strength (matched-strength); (b) the effect of
   merely REDUCING N (→ size-matched random-decimation control, 10 seeds).
4. Cannot: SOZ from the label-based epi_keep_mask (audit_71 pattern); global-rewire surrogate;
   R=200. Patients with 0 SOZ contacts contribute unchanged.
5. Falsify: if β does not clear MS after SOZ removal, or clears no better than / worse than the
   decimation band, the β trace is SOZ-dependent (confounded with epilepsy).

Output: data/sparsified_arc/trace_epi_exclusion_mst020/{per_patient,cohort,decim}.csv
Run: [SA_LIMIT=n] python 25_trace_epi_exclusion_mst020.py
"""
from __future__ import annotations
import os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, load_phase                 # noqa: E402
from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction         # noqa: E402
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, rho_sym_over_scales  # noqa: E402
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle  # noqa: E402
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask     # noqa: E402
from lrg_eegfc.utils.io.regions import load_channel_regions            # noqa: E402

PHASES = ("A", "B", "task_test", "rest_post")
BANDS = ["delta", "alpha", "beta", "low_gamma"]
FRAC = 0.20
SGRID = np.logspace(0.0, np.log10(180.0), 16)
R = 200
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260713
N_DECIM = 10
LIMIT = int(os.environ.get("SA_LIMIT", "0"))
OUT = ROOT / "data" / "sparsified_arc" / "trace_epi_exclusion_mst020"


def eig_mst(W):
    return laplacian_eig(mst_union_top_fraction(W, FRAC))


def sub(W, keep):
    return np.ascontiguousarray(W[np.ix_(keep, keep)])


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception:
        return None
    N = Ws["A"].shape[0]
    rdf = load_channel_regions(pat)
    if len(rdf) != N:
        return None
    keep = epi_keep_mask(rdf, pat)                 # True = non-SOZ
    n_soz = int((~keep).sum())
    kept_idx = np.where(keep)[0]

    # --- SOZ-EXCLUDED trace + matched-strength null ---
    Wx = {ph: sub(Ws[ph], kept_idx) for ph in PHASES}
    Nx = kept_idx.size
    eig_obs = {ph: eig_mst(Wx[ph]) for ph in PHASES}
    obs = rho_sym_over_scales(eig_obs, SGRID)
    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * Nx * (Nx - 1) // 2
    surr = np.full((R, SGRID.size), np.nan)
    for r in range(R):
        eig_s = {ph: eig_mst(matched_strength_shuffle(Wx[ph], n_swaps, rng, W_MAX))
                 for ph in PHASES}
        surr[r] = rho_sym_over_scales(eig_s, SGRID)

    rows = []
    for j, s in enumerate(SGRID):
        col = surr[:, j][np.isfinite(surr[:, j])]; o = obs[j]
        p = float(np.mean(col >= o)) if col.size and np.isfinite(o) else np.nan
        rows.append(dict(patient=pat, band=band, s=float(s), N_full=int(N), N_kept=int(Nx),
                         n_soz=n_soz, obs_rho=float(o), surr_p50=float(np.nanmedian(col)),
                         p=p))

    # --- size-matched random DECIMATION (observed only, 10 seeds) ---
    drows = []
    if n_soz > 0:
        drng = np.random.default_rng(BASE_SEED + 777 + idx)
        decim = np.full((N_DECIM, SGRID.size), np.nan)
        for d in range(N_DECIM):
            drop = drng.choice(kept_idx, size=min(n_soz, Nx - 3), replace=False)
            dkeep = np.setdiff1d(np.arange(N), drop)
            Wd = {ph: sub(Ws[ph], dkeep) for ph in PHASES}
            try:
                decim[d] = rho_sym_over_scales({ph: eig_mst(Wd[ph]) for ph in PHASES}, SGRID)
            except Exception:
                pass
        for j, s in enumerate(SGRID):
            c = decim[:, j][np.isfinite(decim[:, j])]
            drows.append(dict(patient=pat, band=band, s=float(s),
                              decim_mean=float(np.mean(c)) if c.size else np.nan,
                              decim_p05=float(np.percentile(c, 5)) if c.size else np.nan,
                              decim_p95=float(np.percentile(c, 95)) if c.size else np.nan,
                              excl_rho=float(obs[j])))
    return rows, drows


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
            out.append(dict(band=band, s=float(s), gate_p=p, n_above=int((x.p < 0.05).sum()),
                            n_pat=len(x), obs_med=float(x.obs_rho.median())))
    return pd.DataFrame(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in BANDS for p in COHORT)]
    if LIMIT:
        jobs = jobs[:LIMIT]
    print(f"[trace epi-exclusion mst@0.20] {len(jobs)} cells, R={R}, +{N_DECIM}-seed decimation", flush=True)
    t0 = time.time(); rows = []; drows = []
    nproc = min(12, max(1, (os.cpu_count() or 2) - 2))
    with Pool(nproc) as pool:
        for k, res in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if res:
                rows.extend(res[0]); drows.extend(res[1])
            el = time.time() - t0
            print(f"  [{k}/{len(jobs)}] elapsed {el:5.0f}s ETA {el/k*(len(jobs)-k):5.0f}s", flush=True)
    pp = pd.DataFrame(rows); pp.to_csv(OUT / "per_patient.csv", index=False)
    dd = pd.DataFrame(drows); dd.to_csv(OUT / "decim.csv", index=False)
    gate = cohort_gate(pp); gate.to_csv(OUT / "cohort.csv", index=False)

    print("\n=== β TRACE with SOZ REMOVED — cohort gate (Wilcoxon obs-surr, K=10) ===", flush=True)
    for band in BANDS:
        g = gate[gate.band == band]
        if g.empty:
            continue
        nclear = int((g.gate_p < 0.05).sum()); best = g.loc[g.gate_p.idxmin()]
        print(f"  {band:9s}: {nclear}/{len(g)} scales clear MS | best gate_p={best.gate_p:.3f} @s{best.s:.1f} "
              f"(obs_med ρ={best.obs_med:.3f})", flush=True)
    # decimation comparison for beta: is SOZ-exclusion within the random-decimation band?
    print("\n=== DECIMATION control (β): SOZ-excl ρ vs random same-size-drop ρ band, per scale ===", flush=True)
    db = dd[dd.band == "beta"]
    if not db.empty:
        for s in np.sort(db.s.unique())[::3]:
            x = db[np.isclose(db.s, s)]
            excl = x.excl_rho.median(); dmean = x.decim_mean.median()
            lo = x.decim_p05.median(); hi = x.decim_p95.median()
            inside = "within" if lo <= excl <= hi else ("ABOVE" if excl > hi else "below")
            print(f"  s{s:5.1f}: SOZ-excl ρ={excl:.3f} vs decim {dmean:.3f} [{lo:.3f},{hi:.3f}] -> {inside}", flush=True)
    print(f"\n[done] {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
