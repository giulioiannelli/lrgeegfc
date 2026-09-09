#!/usr/bin/env python3
"""Lane M1 — hierarchy-level decomposition of the cross-phase trace, with matched-strength null.
Scope report (5-point preamble in full): .agents/guides/task-persistence-investigation/2026-09-08_hierarchy-level-decomposition.md

Compact preamble.
(1) Claim: the task_test -> rest_post trace (rho_sym) is carried by pairs that merge at FINE levels of the
    rest_pre diffusion hierarchy (k >= 5 clusters remaining) and is absent at the coarsest levels (k <= 4),
    consistently across patients; i.e. the trace is scale-specific on the hierarchy axis.
(2) Null: matched-strength (4-cycle +-delta) surrogates of all five dense |ImCoh| matrices (A, B, task_test,
    rest_post, rest_pre), the whole pipeline re-run (reference tree, backbones, level binning, readouts).
(3) Strongest alternative: the level profile is a reliability gradient (fine-level pairs are within-module,
    strongly connected, hence reproducible in every phase) or a strength-sequence artefact.
(4) The null removes strength-sequence structure and leaves the level populations; it does NOT remove the
    reliability gradient, which is why the split-half reliability rel_l is stored and the profile is also
    tested as rho/rel. It cannot reject "fine pairs are same-module pairs that any state change preserves";
    the xshaft mask and the per-level sham (windowed cache, later) address that.
(5) Falsifiers: flat profile under any knob; trend p > 0.05 on Wilcoxon or LOO; surrogate cohort p > 0.05;
    profile not stable across backbone fraction {0.07, 0.14, 0.20} or reference scale {1, 4.72}.
Outputs: data/paper_final/lane_m1_level/{cells/*.npz, level_profiles.csv, knob_summary.csv, cohort_null.csv}.
Run: PYTHONPATH=src <lapbrain python> -u scripts/01_compute/paper_final/lane_m1_level_decomposition.py [--R 100 --workers 4]
"""
import argparse, time, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.cluster.hierarchy import cophenet
from scipy.stats import spearmanr, wilcoxon
from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, linkage_at_scale, rho_sym
from lrg_eegfc.utils.io.patient import load_channel_labels
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle
from lrg_eegfc.utils.probe import build_probe_mask
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import load_phase_fc
from lrg_eegfc.workflow.substrate import CANONICAL
ROOT = setup_script_env(); OUT = ROOT / "data" / "paper_final" / "lane_m1_level"
PH4 = ("A", "B", "task_test", "rest_post"); PH5 = PH4 + ("rest_pre",)
FRACS = (0.07, 0.14, 0.20); SREFS = (1.0, 4.72); S_READ = 1.0
EDGES = [1, 2, 4, 8, 16, 32, 64, 10**6]; NLEV = 7; FINE0 = 2       # levels idx>=2 <=> k>=5
LEV = ["k2", "k3-4", "k5-8", "k9-16", "k17-32", "k33-64", "k>64"]
MASKS = ("xshaft", "all"); READ = ("rawW", "heat1"); SWAP_FACTOR = 20; W_MAX = 1.0; MIN_PAIRS = 30
R = 100; WORKERS = 4; SEED = 20260909

def levels_from_tree(Wpre, frac, s_ref):
    B = select_backbone(Wpre, CANONICAL.backbone, frac=frac); ev, V = laplacian_eig(B); N = B.shape[0]
    Z = linkage_at_scale(ev, V, s_ref); h = cophenet(Z)
    k = np.maximum(N - np.searchsorted(np.sort(Z[:, 2]), h, side="right"), 1)
    return np.digitize(k, EDGES[1:], right=True)

def readouts(dense, frac, iu):
    Rd = {"rawW": {ph: dense[ph][iu] for ph in PH4}, "heat1": {}}
    for ph in PH4:
        B = select_backbone(dense[ph], CANONICAL.backbone, frac=frac); ev, V = laplacian_eig(B)
        w = np.exp(-(S_READ / ev[-1]) * ev); K = (V * w) @ V.T; Rd["heat1"][ph] = (K / np.trace(K))[iu]
    return Rd

def profiles(Rd, lev, keep):
    """-> rho[readout] (NLEV,), rel[readout] (NLEV,), npairs (NLEV,)"""
    rho = {r: np.full(NLEV, np.nan) for r in READ}; rel = {r: np.full(NLEV, np.nan) for r in READ}; npair = np.zeros(NLEV, int)
    for li in range(NLEV):
        m = keep & (lev == li); npair[li] = m.sum()
        if npair[li] < MIN_PAIRS: continue
        for r in READ:
            X = Rd[r]; rho[r][li] = rho_sym(X["A"][m], X["B"][m], X["task_test"][m], X["rest_post"][m])[0]
            rel[r][li] = spearmanr(X["A"][m], X["B"][m])[0]
    return rho, rel, npair

def stats(p):
    ok = ~np.isnan(p); idx = np.arange(NLEV)
    t = spearmanr(idx[ok], p[ok])[0] if ok.sum() >= 3 else np.nan
    c = np.nanmean(p[FINE0:]) - np.nanmean(p[:FINE0]) if (ok[FINE0:].any() and ok[:FINE0].any()) else np.nan
    return t, c

def per_cell(job):
    i, pat, band, nsurr = job
    dense = {ph: load_phase_fc(pat, ph, band) for ph in PH5}; N = dense["rest_pre"].shape[0]; iu = np.triu_indices(N, 1)
    labels = load_channel_labels(pat); same = build_probe_mask(list(labels))[iu]
    keep = {"xshaft": ~same, "all": np.ones(iu[0].size, bool)}
    rows = []
    for frac in FRACS:
        Rd = readouts(dense, frac, iu)
        for s_ref in SREFS:
            lev = levels_from_tree(dense["rest_pre"], frac, s_ref)
            for mk in MASKS:
                rho, rel, npair = profiles(Rd, lev, keep[mk])
                for r in READ:
                    t, c = stats(rho[r]); tn, cn = stats(rho[r] / rel[r])
                    rows.append(dict(patient=pat, band=band, frac=frac, s_ref=s_ref, mask=mk, readout=r, trend=t, contrast=c,
                                     trend_norm=tn, contrast_norm=cn, **{f"rho_{LEV[l]}": rho[r][l] for l in range(NLEV)},
                                     **{f"rel_{LEV[l]}": rel[r][l] for l in range(NLEV)}, **{f"n_{LEV[l]}": npair[l] for l in range(NLEV)}))
    # matched-strength null at the canonical knob (frac=CANONICAL.frac, s_ref=1)
    rng = np.random.default_rng(SEED + 1000 * i); nsw = SWAP_FACTOR * N * (N - 1) // 2
    surr = np.full((nsurr, len(MASKS), len(READ), 4), np.nan)           # stats: trend, contrast, trend_norm, contrast_norm
    surr_prof = np.full((nsurr, len(MASKS), len(READ), NLEV), np.nan)
    surr_tree = np.full_like(surr, np.nan)                                  # tree-only null: observed readouts, surrogate reference tree
    Rd_obs = readouts(dense, CANONICAL.frac, iu); t0 = time.time()
    for rr in range(nsurr):
        sh = {ph: matched_strength_shuffle(dense[ph], nsw, rng, W_MAX) for ph in PH5}
        if rr == 0:
            assert max(np.abs(sh[ph].sum(1) - dense[ph].sum(1)).max() for ph in PH5) < 1e-6, "strengths not preserved"
            assert np.abs(sh["rest_pre"] - dense["rest_pre"]).max() > 1e-3, "surrogate identical to data"
        Rd = readouts(sh, CANONICAL.frac, iu); lev = levels_from_tree(sh["rest_pre"], CANONICAL.frac, 1.0)
        lev_t = levels_from_tree(matched_strength_shuffle(dense["rest_pre"], nsw, rng, W_MAX), CANONICAL.frac, 1.0)
        for mi, mk in enumerate(MASKS):
            rho, rel, _ = profiles(Rd, lev, keep[mk]); rho_t, rel_t, _ = profiles(Rd_obs, lev_t, keep[mk])
            for ri, r in enumerate(READ):
                surr[rr, mi, ri, :2] = stats(rho[r]); surr[rr, mi, ri, 2:] = stats(rho[r] / rel[r]); surr_prof[rr, mi, ri] = rho[r]
                surr_tree[rr, mi, ri, :2] = stats(rho_t[r]); surr_tree[rr, mi, ri, 2:] = stats(rho_t[r] / rel_t[r])
        if rr == 1: print(f"    cell {pat}/{band}: {(time.time()-t0)/2:.2f} s/surrogate", flush=True)
    np.savez(OUT / "cells" / f"{pat}_{band}.npz", surr=surr, surr_tree=surr_tree, surr_prof=surr_prof, masks=MASKS, readouts=READ, levels=LEV)
    return rows

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--patients", default=""); ap.add_argument("--bands", default="")
    ap.add_argument("--R", type=int, default=R); ap.add_argument("--workers", type=int, default=WORKERS); a = ap.parse_args()
    pats = a.patients.split(",") if a.patients else list(PATIENTS_4PHASE); bands = a.bands.split(",") if a.bands else list(BRAIN_BANDS_NAMES)
    (OUT / "cells").mkdir(parents=True, exist_ok=True); matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    jobs = [(i, p, b, a.R) for i, (b, p) in enumerate((b, p) for b in bands for p in pats)]
    print(f"[M1] {len(jobs)} cells | R={a.R} | fracs={FRACS} s_ref={SREFS} | {a.workers} workers -> {OUT}", flush=True)
    t0 = time.time(); rows = []
    with Pool(a.workers) as pool:
        for i, res in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            rows += res; el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {res[0]['patient']}/{res[0]['band']} elapsed {el/60:.1f} min ETA {el/i*(len(jobs)-i)/60:.1f} min", flush=True)
    df = pd.DataFrame(rows); df.to_csv(OUT / "level_profiles.csv", index=False)
    # knob summary: cohort Wilcoxon (two-sided) + LOO on trend / contrast, all knobs
    ks = []
    for key, g in df.groupby(["band", "frac", "s_ref", "mask", "readout"]):
        for st in ("trend", "contrast", "trend_norm", "contrast_norm"):
            x = g[st].to_numpy(float); x = x[~np.isnan(x)]
            if x.size < 5: continue
            p = wilcoxon(x).pvalue; loo = max(wilcoxon(np.delete(x, j)).pvalue for j in range(x.size))
            ks.append(dict(zip(["band", "frac", "s_ref", "mask", "readout"], key), stat=st, n=x.size, mean=x.mean(), median=np.median(x),
                           n_pos=int((x > 0).sum()), p_wilcoxon=p, p_loo_max=loo))
    K = pd.DataFrame(ks); K.to_csv(OUT / "knob_summary.csv", index=False)
    # matched-strength cohort null at the canonical knob
    cn = []
    for band in bands:
      for null_name, key_ in (("full", "surr"), ("tree", "surr_tree")):
        cells = [np.load(OUT / "cells" / f"{p}_{band}.npz")[key_] for p in pats]
        S = np.stack(cells)                                    # (P, R, M, RD, 4)
        obs = df[(df.band == band) & (df.frac == CANONICAL.frac) & (df.s_ref == 1.0)]
        for mi, mk in enumerate(MASKS):
            for ri, r in enumerate(READ):
                for si, st in enumerate(("trend", "contrast", "trend_norm", "contrast_norm")):
                    o = np.array([obs[(obs.patient == p) & (obs["mask"] == mk) & (obs.readout == r)][st].iloc[0] for p in pats], float)
                    s = S[:, :, mi, ri, si]                    # (P, R)
                    ok = ~np.isnan(o); om = np.nanmean(o); sm = np.nanmean(s[ok], axis=0)
                    pp = np.array([np.mean(s[j][~np.isnan(s[j])] >= o[j]) if ok[j] else np.nan for j in range(len(pats))])
                    cn.append(dict(band=band, null=null_name, mask=mk, readout=r, stat=st, obs_mean=om, null_mean=np.nanmean(sm), null_sd=np.nanstd(sm),
                                   p_cohort_surr=float(np.mean(sm >= om)), z_cohort=(om - np.nanmean(sm)) / (np.nanstd(sm) + 1e-12),
                                   n_pat_p05=int(np.nansum(pp < 0.05)), n_pat=int(ok.sum())))
    C = pd.DataFrame(cn); C.to_csv(OUT / "cohort_null.csv", index=False)
    pd.set_option("display.width", 220, "display.max_rows", 500)
    v = C[(C["mask"] == "xshaft") & (C.stat.isin(["trend", "contrast"]))].merge(
        K[(K.frac == CANONICAL.frac) & (K.s_ref == 1.0) & (K["mask"] == "xshaft")][["band", "readout", "stat", "n_pos", "p_wilcoxon", "p_loo_max"]],
        on=["band", "readout", "stat"])
    print("\n== VERDICT (xshaft, canonical knob): observed vs matched-strength null ==\n",
          v[["band", "null", "readout", "stat", "obs_mean", "null_mean", "null_sd", "z_cohort", "p_cohort_surr", "n_pat_p05", "n_pos", "p_wilcoxon", "p_loo_max"]]
          .sort_values(["stat", "readout", "band", "null"]).to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    kk = K[(K["mask"] == "xshaft") & (K.stat == "contrast")].pivot_table(index=["band", "readout"], columns=["frac", "s_ref"], values="p_wilcoxon")
    print("\n== knob stability: Wilcoxon p of fine-coarse contrast (xshaft) by frac x s_ref ==\n", kk.to_string(float_format=lambda x: f"{x:.3f}"))
    prof = df[(df.frac == CANONICAL.frac) & (df.s_ref == 1.0) & (df["mask"] == "xshaft")].groupby(["band", "readout"])[[f"rho_{l}" for l in LEV]].median()
    print("\n== median level profile rho_sym (xshaft, canonical) ==\n", prof.to_string(float_format=lambda x: f"{x:.2f}"))
    print(f"\nwall {(time.time()-t0)/60:.1f} min -> {OUT}", flush=True)

if __name__ == "__main__":
    main()
