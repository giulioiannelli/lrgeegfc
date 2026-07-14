#!/usr/bin/env python3
"""Localization ARC on mst@0.20 — WHERE the trace, encoding and inference live,
across ALL 16 scales, under one matched-strength null.

Supersedes 16_localization_mst020 (trace only, 3 scales). Ports the whole-graph
per-system localizers — audit_151 (trace), audit_158 (encoding->OFC), audit_160
(inference-specific | encoding -> cingulate) — onto the recovered mst@0.20 backbone,
and sweeps the FULL dimensionless scale grid s = tau*lambda_max in [1, 180] (16 pts,
identical to 13/16). One eigendecomposition per phase per realization; cophenetic
evaluated at every scale. The decisive question for R2 anatomy: does encoding->OFC
and inference->cingulate survive on the non-degenerate backbone at ANY scale, or do
they delocalize like the beta TRACE did (16_localization: nothing clears BH)?

Three per-pair concordance targets (all rank-based centered-rank products; the
inference target residualizes the same ranks -> one consistent primitive):
  trace       s = 1/2[ conc(cT-cA, cP-cB) + conc(cT-cB, cP-cA) ]          (audit_151)
  encoding    s = 1/2[ conc(cL-cA, cP-cB) + conc(cL-cB, cP-cA) ]          (audit_158)
  inference   s = 1/2[ partial(cT-cL, cP-cB | cL-cA)                        (audit_160)
                     + partial(cT-cL, cP-cA | cL-cB) ]  (inference-specific | encoding)
with cA,cB = rest_pre split halves, cL = task_learn, cT = task_test, cP = rest_post
cophenetics at scale s. Endpoint-incidence demeaned per-unit means -> per-system +
binary SOZ. Matched-strength surrogate IDENTICAL to 13/16 (4-cycle +/-delta over all
5 phases -> sparsify@0.20 -> LRG), R=200. Cohort Wilcoxon(obs-surr_p50, greater) + BH
per (band, grouping, scale).

5-point preamble
1. Claim: on mst@0.20 the encoding trace concentrates in OFC and the inference-specific
   trace in cingulate, at some scale, above a strength-preserving null.
2. Null: concentration is an artifact of node strength / the A|B arm choice.
3. Strongest alternative it must beat: node strength -> the IDENTICAL matched-strength
   surrogate (preserves every node's total coupling); arm choice -> symmetrized over A|B.
4. Cannot: R=200; a-priori atlas + binary SOZ; BH over systems within band x scale.
   OFC/cingulate sampled in 5/8 patients (coverage limit). Does not upgrade a hint to a
   hotspot; reports sign+lead+p/q honestly per scale.
5. Falsify: if OFC (encoding) / cingulate (inference) never clears the upper tail at any
   scale, R2 loses its anatomy and rests on scale-signature + higher-order + band
   dissociation (same fate as the beta TRACE -> delocalized).

Output: data/sparsified_arc/localization_arc_mst020/{system,soz}_enrichment.csv (all
targets x bands x scales) + per_patient.csv.
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
from lrg_eegfc.utils.metrics.node_localization import rank_concordance, epi_keep_mask  # noqa: E402
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr             # noqa: E402
from lrg_eegfc.utils.io.regions import load_channel_regions       # noqa: E402

PHASES = ("A", "B", "task_learn", "task_test", "rest_post")
BANDS = ["delta", "alpha", "beta", "low_gamma"]
TARGETS = ("trace", "encoding", "inference")
FRAC = 0.20
SGRID = np.logspace(0.0, np.log10(180.0), 16)    # dimensionless s = tau*lambda_max
R = 200
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260713
DROP_SYS = {"non_anatomical"}
LIMIT = int(os.environ.get("SA_LIMIT", "0"))     # >0 -> time a few cells then stop


def _rc(x):
    return rankdata(x) - (x.size + 1) / 2.0


def concordance_partial(a, b, c):
    """Per-pair contribution to first-order partial Spearman rho(a,b|c) — audit_110."""
    ra, rb, rc = _rc(a), _rc(b), _rc(c)
    rc2 = float(np.dot(rc, rc))
    if rc2 == 0.0:
        return ra * rb
    resid_a = ra - (np.dot(ra, rc) / rc2) * rc
    resid_b = rb - (np.dot(rb, rc) / rc2) * rc
    return resid_a * resid_b


def sym_trace(cA, cB, cT, cP):
    return 0.5 * (rank_concordance(cT - cA, cP - cB) + rank_concordance(cT - cB, cP - cA))


def sym_encoding(cA, cB, cL, cP):
    return 0.5 * (rank_concordance(cL - cA, cP - cB) + rank_concordance(cL - cB, cP - cA))


def sym_inference(cA, cB, cL, cT, cP):
    """inference-specific | encoding, symmetric over the A|B arm (f = cT-cL is arm-invariant)."""
    f = cT - cL
    return 0.5 * (concordance_partial(f, cP - cB, cL - cA)
                  + concordance_partial(f, cP - cA, cL - cB))


def concordances_at_scale(C):
    """C = {phase: cophenetic}; return {target: per-pair vector} at one scale."""
    cA, cB, cL, cT, cP = C["A"], C["B"], C["task_learn"], C["task_test"], C["rest_post"]
    return {"trace": sym_trace(cA, cB, cT, cP),
            "encoding": sym_encoding(cA, cB, cL, cP),
            "inference": sym_inference(cA, cB, cL, cT, cP)}


def unit_means(s, nidx_i, nidx_j, unit_vec, drop=frozenset()):
    """Endpoint-incidence, demeaned per-unit mean of a per-pair concordance vector."""
    vals = np.concatenate([s, s]); nidx = np.concatenate([nidx_i, nidx_j])
    vals = vals - vals.mean(); u = unit_vec[nidx]
    order = np.argsort(u, kind="stable"); us, vs = u[order], vals[order]
    uniq, starts = np.unique(us, return_index=True)
    sums = np.add.reduceat(vs, starts)
    counts = np.diff(np.append(starts, vs.size))
    return {un: sm / ct for un, sm, ct in zip(uniq, sums, counts) if un not in drop}


def eig_backbone(W):
    return laplacian_eig(mst_union_top_fraction(W, FRAC))


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception:
        return None
    N = Ws["A"].shape[0]
    if any(v.shape[0] != N for v in Ws.values()):
        return None
    rdf = load_channel_regions(pat)
    if len(rdf) != N:
        print(f"  [skip] {pat} {band}: regions {len(rdf)} != FC {N}", flush=True)
        return None
    sys_vec = rdf["system"].to_numpy().astype(str)
    soz_vec = np.where(epi_keep_mask(rdf, pat), "non_soz", "soz")
    iu = np.triu_indices(N, 1); ii, jj = iu[0].astype(int), iu[1].astype(int)
    NS = len(SGRID)

    def all_unit_means(eig5):
        """{target: {grouping: {scale_idx: {unit: mean}}}} for one FC set (obs or surr)."""
        res = {t: {"system": [None] * NS, "soz": [None] * NS} for t in TARGETS}
        for si, s in enumerate(SGRID):
            C = {ph: cophenetic_at_scale(*eig5[ph], s) for ph in PHASES}
            conc = concordances_at_scale(C)
            for t in TARGETS:
                res[t]["system"][si] = unit_means(conc[t], ii, jj, sys_vec, drop=DROP_SYS)
                res[t]["soz"][si] = unit_means(conc[t], ii, jj, soz_vec)
        return res

    eig_obs = {ph: eig_backbone(Ws[ph]) for ph in PHASES}
    obs = all_unit_means(eig_obs)

    # surrogate accumulation: surr[t][grp][si][unit] -> list over R
    surr = {t: {grp: [dict() for _ in range(NS)] for grp in ("system", "soz")} for t in TARGETS}
    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    for r in range(R):
        eig_sh = {ph: eig_backbone(matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX))
                  for ph in PHASES}
        try:
            sm = all_unit_means(eig_sh)
        except Exception:
            continue
        for t in TARGETS:
            for grp in ("system", "soz"):
                for si in range(NS):
                    for u, val in sm[t][grp][si].items():
                        surr[t][grp][si].setdefault(u, []).append(val)

    rows = []
    for t in TARGETS:
        for grp in ("system", "soz"):
            for si in range(NS):
                for u, o in obs[t][grp][si].items():
                    col = np.array(surr[t][grp][si].get(u, []), float)
                    col = col[np.isfinite(col)]
                    p = float((np.sum(col >= o) + 1) / (col.size + 1)) if col.size else np.nan
                    rows.append(dict(patient=pat, band=band, target=t, grouping=grp,
                                     scale_idx=si, s=float(SGRID[si]), unit=u, obs=float(o),
                                     surr_p50=float(np.nanmedian(col)) if col.size else np.nan,
                                     p=p))
    return rows


def cohort(df, grouping):
    """Per (target, band, scale, unit): Wilcoxon(obs-surr_p50, greater) + BH per (band,scale)."""
    out = []
    d = df[df.grouping == grouping]
    for t in TARGETS:
        dt = d[d.target == t]
        for band in BANDS:
            db = dt[dt.band == band]
            for si in sorted(db.scale_idx.unique()):
                ds = db[db.scale_idx == si]
                rows = []
                for u in sorted(ds.unit.unique()):
                    x = ds[ds.unit == u].dropna(subset=["obs", "surr_p50"])
                    if len(x) < 5:
                        continue
                    diff = x.obs.values - x.surr_p50.values
                    try:
                        _, p = wilcoxon(diff, alternative="greater")
                    except ValueError:
                        p = 1.0
                    rows.append(dict(target=t, band=band, scale_idx=int(si),
                                     s=float(ds.s.iloc[0]), unit=u, n=len(x),
                                     obs_med=float(x.obs.median()),
                                     n_pos=int((diff > 0).sum()), p=float(p)))
                if rows:
                    rdf = pd.DataFrame(rows)
                    rdf["q"] = bh_fdr(rdf.p.values)
                    out.append(rdf)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def _summarize(sys_c, soz_c):
    """Print the decisive read: for each target/band, the best (min-p) system per scale,
    and whether OFC/cingulate/insula/PFC ever clear BH."""
    if sys_c.empty:
        print("\n[summary skipped — cohort frame empty (need >=5 patients/unit)]", flush=True)
        return
    KEY = {"trace": None, "encoding": "OFC", "inference": "cingulate"}
    for t in TARGETS:
        print(f"\n########## TARGET = {t} ##########", flush=True)
        st = sys_c[sys_c.target == t]
        for band in BANDS:
            b = st[st.band == band]
            if b.empty:
                continue
            # best system at each scale
            best = b.loc[b.groupby("scale_idx").p.idxmin()].sort_values("scale_idx")
            leads = ", ".join(f"s{r.s:05.1f}:{r.unit}(p{r.p:.3f}{'*' if r.q<0.05 else ''})"
                              for _, r in best.iterrows())
            nclear = int((b.q < 0.05).sum())
            print(f"  [{band}] BH-clearing cells: {nclear}", flush=True)
            print(f"      per-scale lead: {leads}", flush=True)
            key = KEY[t]
            if key is not None:
                kk = b[b.unit == key]
                if not kk.empty:
                    kbest = kk.loc[kk.p.idxmin()]
                    print(f"      {key}: best p={kbest.p:.3f} q={kbest.q:.3f} @s{kbest.s:.1f} "
                          f"(obs_med {kbest.obs_med:+.0f}, {kbest.n_pos}/{kbest.n}); "
                          f"clears BH at {int((kk.q<0.05).sum())}/{len(kk)} scales", flush=True)
    print("\n=== SOZ enrichment (soz vs non_soz), best scale per target/band ===", flush=True)
    for t in TARGETS:
        st = soz_c[(soz_c.target == t) & (soz_c.unit == "soz")]
        for band in BANDS:
            b = st[st.band == band]
            if b.empty:
                continue
            kb = b.loc[b.p.idxmin()]
            print(f"  {t:10s} {band:10s} soz best p={kb.p:.3f} q={kb.q:.3f} @s{kb.s:.1f} "
                  f"({kb.n_pos}/{kb.n})  clears {int((b.q<0.05).sum())}/{len(b)} scales", flush=True)


def main():
    OUT = ROOT / "data" / "sparsified_arc" / "localization_arc_mst020"
    OUT.mkdir(parents=True, exist_ok=True)
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))  # warm numba
    jobs = [(i, pat, band) for i, (pat, band) in
            enumerate((p, b) for b in BANDS for p in COHORT)]
    if LIMIT:
        jobs = jobs[:LIMIT]
    print(f"[localization-arc mst@0.20] R={R}  {len(SGRID)} scales  {len(jobs)} cells "
          f"({len(BANDS)}b x {len(COHORT)}p, 5 phases, targets={TARGETS})", flush=True)
    t0 = time.time(); recs = []
    nproc = min(12, max(1, (os.cpu_count() or 2) - 2))
    with Pool(nproc) as pool:
        for k, rl in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if rl:
                recs.extend(rl)
            el = time.time() - t0
            print(f"  [{k}/{len(jobs)}] elapsed {el:6.0f}s  ETA {el/k*(len(jobs)-k):6.0f}s",
                  flush=True)
    per_pat = pd.DataFrame(recs)
    per_pat.to_csv(OUT / "per_patient.csv", index=False)
    sys_c = cohort(per_pat, "system"); sys_c.to_csv(OUT / "system_enrichment.csv", index=False)
    soz_c = cohort(per_pat, "soz");   soz_c.to_csv(OUT / "soz_enrichment.csv", index=False)
    _summarize(sys_c, soz_c)
    print(f"\n[done] {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
