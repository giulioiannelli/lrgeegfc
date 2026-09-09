#!/usr/bin/env python3
"""Lane T1 — the trace as a function of diffusion time: rho_sym(s), s = tau*lambda_max, with a matched-strength
null at EVERY s, per patient x band, and the cross-patient alignment of the tau window.
This is objective 1 as posed (scale = tau). The hierarchy-level construction (lane_m1) is a single-tau object and
is NOT a scale axis (user, 2026-09-09).

Preamble.
(1) Claim: for some band, the task_test -> rest_post trace measured on the LRG tree at scale s occupies a specific
    window of s, consistently across patients (peak positions aligned, window narrower than the grid), and is not
    present uniformly across all s nor only at the raw-FC limit s -> 0.
(2) Null: matched-strength (4-cycle +-delta) surrogates of the four dense matrices A, B, task_test, rest_post,
    backbone + eigendecomposition + tree recomputed, the full s profile recomputed per surrogate (R per cell).
    The excess profile e(s) = rho_obs(s) - mean_null(s) and z(s) are the objects tested.
(3) Strongest alternatives: (a) the profile is flat in s and only its level is nonzero (trace at all scales = not
    scale-specific); (b) the profile peaks where the tree is most reliable (split-half rel(s)), i.e. a reliability
    shape, not a trace shape; (c) peak positions scatter across patients as much as under the null (no alignment).
(4) The per-s null removes strength-sequence structure at every s and gives a null distribution of PEAK positions
    (each surrogate profile has an argmax), which is the reference for (c). It does not remove the reliability
    shape (b); rel(s) is stored and the excess profile is compared with it. It cannot reject "any state change".
(5) Falsifiers: cohort Wilcoxon on e(s) significant at (nearly) all s with a flat median (alternative a); peak
    s* IQR across patients not narrower than the null s* IQR (c); e(s) shape equal to rel(s) shape (b);
    instability between backbone fractions 0.14 and 0.20.
Readouts: 'coph' = cophenetic distances of the UPGMA tree of 1/rho(s) (the canonical LRG pipeline);
'wave' = graph-wavelet pair vector -tau d(rho)/d(tau) at s (band-pass in the Laplacian spectrum, non-cumulative).
Outputs: data/paper_final/lane_t1_tau/{cells/*.npz, profiles.csv, cohort_by_s.csv, peaks.csv}.
Run: PYTHONPATH=src <lapbrain python> -u <this file> [--R 200 --workers 8]
"""
import argparse, time, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import spearmanr, wilcoxon
from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, cophenetic_at_scale, rho_sym, rho_sym_over_scales
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import load_phase_fc
from lrg_eegfc.workflow.substrate import CANONICAL
ROOT = setup_script_env(); OUT = ROOT / "data" / "paper_final" / "lane_t1_tau"
PH = ("A", "B", "task_test", "rest_post"); FRACS = (0.14, 0.20)
S = np.logspace(-1, 2, 25); SWAP_FACTOR = 20; W_MAX = 1.0; SEED = 20260909; R = 200; WORKERS = 8

def wave_vec(ev, V, s, iu):
    tau = s / ev[-1]; w = np.exp(-tau * ev); Z = w.sum(); lam = (ev * w).sum() / Z
    return (tau * ((V * ((ev - lam) * w)) @ V.T) / Z)[iu]

def profiles(eig, iu):
    """-> coph (nS,), wave (nS,), rel_coph (nS,)"""
    coph = rho_sym_over_scales(eig, S)
    wave = np.full(S.size, np.nan); rel = np.full(S.size, np.nan)
    for i, s in enumerate(S):
        Wv = {ph: wave_vec(*eig[ph], s, iu) for ph in PH}
        wave[i] = rho_sym(Wv["A"], Wv["B"], Wv["task_test"], Wv["rest_post"])[0]
        try: rel[i] = spearmanr(cophenetic_at_scale(*eig["A"], s), cophenetic_at_scale(*eig["B"], s))[0]
        except Exception: pass
    return coph, wave, rel

def per_cell(job):
    i, pat, band, nsurr = job
    dense = {ph: load_phase_fc(pat, ph, band) for ph in PH}; N = dense["A"].shape[0]; iu = np.triu_indices(N, 1)
    rng = np.random.default_rng(SEED + 1000 * i); nsw = SWAP_FACTOR * N * (N - 1) // 2
    rows = []; out = {}
    raw = rho_sym(dense["A"][iu], dense["B"][iu], dense["task_test"][iu], dense["rest_post"][iu])[0]
    for frac in FRACS:
        eig = {ph: laplacian_eig(select_backbone(dense[ph], CANONICAL.backbone, frac=frac)) for ph in PH}
        coph, wave, rel = profiles(eig, iu)
        surr = np.full((nsurr, 2, S.size), np.nan); t0 = time.time()
        for rr in range(nsurr):
            sh = {ph: matched_strength_shuffle(dense[ph], nsw, rng, W_MAX) for ph in PH}
            eg = {ph: laplacian_eig(select_backbone(sh[ph], CANONICAL.backbone, frac=frac)) for ph in PH}
            c, w, _ = profiles(eg, iu); surr[rr, 0] = c; surr[rr, 1] = w
            if rr == 1 and frac == FRACS[0]: print(f"    {pat}/{band}: {(time.time()-t0)/2:.2f} s/surrogate", flush=True)
        out[f"surr_{frac}"] = surr; out[f"coph_{frac}"] = coph; out[f"wave_{frac}"] = wave; out[f"rel_{frac}"] = rel
        for k, (name, obs) in enumerate((("coph", coph), ("wave", wave))):
            nm = np.nanmean(surr[:, k], 0); ns = np.nanstd(surr[:, k], 0)
            for si, s in enumerate(S):
                rows.append(dict(patient=pat, band=band, frac=frac, readout=name, s=s, obs=obs[si], null_mean=nm[si], null_sd=ns[si],
                                 excess=obs[si] - nm[si], z=(obs[si] - nm[si]) / (ns[si] + 1e-12), p=float(np.nanmean(surr[:, k, si] >= obs[si])),
                                 rel=rel[si], rho_raw=raw))
    np.savez(OUT / "cells" / f"{pat}_{band}.npz", S=S, **out)
    return rows

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--patients", default=""); ap.add_argument("--bands", default="")
    ap.add_argument("--R", type=int, default=R); ap.add_argument("--workers", type=int, default=WORKERS); a = ap.parse_args()
    pats = a.patients.split(",") if a.patients else list(PATIENTS_4PHASE); bands = a.bands.split(",") if a.bands else list(BRAIN_BANDS_NAMES)
    (OUT / "cells").mkdir(parents=True, exist_ok=True); matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    jobs = [(i, p, b, a.R) for i, (b, p) in enumerate((b, p) for b in bands for p in pats)]
    print(f"[T1] {len(jobs)} cells | R={a.R} | fracs={FRACS} | {S.size} scales s in [{S[0]:.2f},{S[-1]:.0f}] | {a.workers} workers -> {OUT}", flush=True)
    t0 = time.time(); rows = []
    with Pool(a.workers) as pool:
        for i, res in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            rows += res; el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {res[0]['patient']}/{res[0]['band']} elapsed {el/60:.1f} min ETA {el/i*(len(jobs)-i)/60:.1f} min", flush=True)
    df = pd.DataFrame(rows); df.to_csv(OUT / "profiles.csv", index=False)
    # cohort per s: Wilcoxon (two-sided) on excess across patients + LOO; per-patient surrogate p count
    cs = []
    for (band, frac, rd, s), g in df.groupby(["band", "frac", "readout", "s"]):
        e = g["excess"].to_numpy(float); ok = ~np.isnan(e); e = e[ok]
        if e.size < 5: continue
        p = wilcoxon(e).pvalue; loo = max(wilcoxon(np.delete(e, j)).pvalue for j in range(e.size))
        cs.append(dict(band=band, frac=frac, readout=rd, s=s, n=e.size, obs_med=g["obs"].median(), null_med=g["null_mean"].median(), excess_med=np.median(e),
                       n_pos=int((e > 0).sum()), p_wilcoxon=p, p_loo_max=loo, n_pat_p05=int((g["p"] < 0.05).sum()), rel_med=g["rel"].median()))
    C = pd.DataFrame(cs); C.to_csv(OUT / "cohort_by_s.csv", index=False)
    # peaks: per patient argmax of excess, half-max window; null peak dispersion from surrogate profiles
    pk = []
    for band in bands:
        for frac in FRACS:
            for k, rd in enumerate(("coph", "wave")):
                obs_pk, null_pk, widths = [], [], []
                for p_ in pats:
                    z = np.load(OUT / "cells" / f"{p_}_{band}.npz"); surr = z[f"surr_{frac}"][:, k]; obs = z[f"{rd}_{frac}"]
                    ex = obs - np.nanmean(surr, 0)
                    if np.all(np.isnan(ex)): continue
                    j = int(np.nanargmax(ex)); obs_pk.append(np.log10(S[j]))
                    above = np.where(ex >= 0.5 * ex[j])[0]; widths.append(np.log10(S[above.max()]) - np.log10(S[above.min()]))
                    sp = surr - np.nanmean(surr, 0); null_pk += [np.log10(S[int(np.nanargmax(row))]) for row in sp if not np.all(np.isnan(row))]
                obs_pk = np.array(obs_pk); null_pk = np.array(null_pk)
                iqr = lambda x: np.percentile(x, 75) - np.percentile(x, 25)
                # null IQR of 10 peaks: resample 10 surrogate peaks 2000x
                rng = np.random.default_rng(0); niqr = np.array([iqr(rng.choice(null_pk, obs_pk.size, replace=False)) for _ in range(2000)])
                pk.append(dict(band=band, frac=frac, readout=rd, n=obs_pk.size, s_peak_median=10 ** np.median(obs_pk), s_peak_q25=10 ** np.percentile(obs_pk, 25),
                               s_peak_q75=10 ** np.percentile(obs_pk, 75), iqr_log10=iqr(obs_pk), null_iqr_log10=np.median(niqr), p_aligned=float(np.mean(niqr <= iqr(obs_pk))),
                               halfmax_width_log10_med=np.median(widths), grid_span_log10=3.0))
    P = pd.DataFrame(pk); P.to_csv(OUT / "peaks.csv", index=False)
    pd.set_option("display.width", 240, "display.max_rows", 600)
    for rd in ("coph", "wave"):
        print(f"\n== {rd}, frac 0.14: cohort median excess e(s) [Wilcoxon p, n_pos] by band x s ==")
        g = C[(C.readout == rd) & (C.frac == 0.14)]
        tab = g.pivot(index="band", columns="s", values="excess_med"); tab.columns = [f"{c:.2g}" for c in tab.columns]
        print(tab.to_string(float_format=lambda x: f"{x:.2f}"))
        pt = g.pivot(index="band", columns="s", values="p_wilcoxon"); pt.columns = [f"{c:.2g}" for c in pt.columns]
        print("p_wilcoxon:\n" + pt.to_string(float_format=lambda x: f"{x:.2f}"))
        npos = g.pivot(index="band", columns="s", values="n_pos"); npos.columns = [f"{c:.2g}" for c in npos.columns]
        print("n_pos:\n" + npos.to_string())
    print("\n== reliability rel(s) median, coph, frac 0.14 ==")
    g = C[(C.readout == "coph") & (C.frac == 0.14)].pivot(index="band", columns="s", values="rel_med"); g.columns = [f"{c:.2g}" for c in g.columns]
    print(g.to_string(float_format=lambda x: f"{x:.2f}"))
    print("\n== peak alignment across patients ==\n", P.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print(f"\nwall {(time.time()-t0)/60:.1f} min -> {OUT}", flush=True)

if __name__ == "__main__":
    main()
