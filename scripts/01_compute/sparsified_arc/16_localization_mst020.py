#!/usr/bin/env python3
"""Localization of the mst@0.20 cophenetic trace — WHERE each band lives.

Ports the whole-graph per-system localization (audit_151 / audit_158) onto the
recovered mst@0.20 multiscale backbone, at the joint report scale s=S_REPORT. This
is the non-circular arbiter for the band story: a *cognitive* trace should
concentrate in task-relevant cortex (predict beta -> OFC); a *disease* band's
cross-phase structure should concentrate in the seizure-onset zone (predict
delta/low_gamma -> SOZ); a *placeless* band concentrates nowhere (predict alpha).

5-point preamble
1. Claim: on mst@0.20, the per-pair rho_sym trace CONCENTRATES by band in a way that
   separates cognitive (beta->OFC) from disease (delta->SOZ) tissue.
2. Null: the concentration is an artifact of node strength / the A|B arm choice.
3. Strongest alternative it must beat: node strength. Controlled by the IDENTICAL
   matched-strength surrogate script 13 uses (shuffle each phase's dense FC ->
   sparsify@0.20 -> LRG), which preserves every node's total coupling.
4. Cannot: R=200 (matches 13's lock); systems = a-priori anatomical atlas + a binary
   SOZ split; BH over the anatomical family per band. Shaft-collapse/LOO are separate.
5. Falsify: if no system clears the upper tail (BH q) for beta, or if delta does not
   concentrate in the SOZ, the "cognitive vs disease band" reading is unsupported.

Null is IDENTICAL to 13_matched_strength_mst020 (bit-compatible 4-cycle +/-delta).
Output: data/sparsified_arc/localization_mst020/{system,soz}_enrichment.csv
"""
from __future__ import annotations
import os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import wilcoxon

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

PHASES = ("A", "B", "task_test", "rest_post")
BANDS = ["delta", "alpha", "beta", "low_gamma"]
FRAC = 0.20
SGRID = np.logspace(0.0, np.log10(180.0), 16)
S_REPORT = float(os.environ.get("SA_SCALE", SGRID[5]))   # default 5.646; override to sweep
OUT = None                                               # set in main() from S_REPORT
R = 200
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260713
DROP_SYS = {"non_anatomical"}


def coph(W):
    ev, V = laplacian_eig(mst_union_top_fraction(W, FRAC))
    return cophenetic_at_scale(ev, V, S_REPORT)          # condensed (triu, k=1)


def sym_concordance(cA, cB, ctt, cpost):
    """rho_sym per-pair concordance = mean of both A|B arm assignments."""
    return 0.5 * (rank_concordance(ctt - cA, cpost - cB)
                  + rank_concordance(ctt - cB, cpost - cA))


def unit_means(s, iu_i, iu_j, unit_vec, keep_node=None, drop=frozenset()):
    """Endpoint-incidence, demeaned per-unit mean of a per-pair concordance vector."""
    vals = np.concatenate([s, s]); nidx = np.concatenate([iu_i, iu_j])
    if keep_node is not None:
        m = keep_node[nidx]; vals, nidx = vals[m], nidx[m]
    vals = vals - vals.mean(); u = unit_vec[nidx]
    order = np.argsort(u, kind="stable"); us, vs = u[order], vals[order]
    uniq, starts = np.unique(us, return_index=True)
    sums = np.add.reduceat(vs, starts)
    counts = np.diff(np.append(starts, vs.size))
    return {un: sm / ct for un, sm, ct in zip(uniq, sums, counts) if un not in drop}


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception:
        return None
    N = Ws["A"].shape[0]
    rdf = load_channel_regions(pat)
    if len(rdf) != N:                                    # alignment guard
        print(f"  [skip] {pat} {band}: regions {len(rdf)} != FC {N}", flush=True)
        return None
    sys_vec = rdf["system"].to_numpy().astype(str)
    soz_vec = np.where(epi_keep_mask(rdf, pat), "non_soz", "soz")
    iu = np.triu_indices(N, 1); ii, jj = iu[0].astype(int), iu[1].astype(int)

    def concordance(Wdict):
        C = {ph: coph(Wdict[ph]) for ph in PHASES}
        return sym_concordance(C["A"], C["B"], C["task_test"], C["rest_post"])

    s_obs = concordance(Ws)
    obs_sys = unit_means(s_obs, ii, jj, sys_vec, drop=DROP_SYS)
    obs_soz = unit_means(s_obs, ii, jj, soz_vec)

    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    surr_sys = {u: np.full(R, np.nan) for u in obs_sys}
    surr_soz = {u: np.full(R, np.nan) for u in obs_soz}
    for r in range(R):
        Wsh = {ph: matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX) for ph in PHASES}
        try:
            s_r = concordance(Wsh)
        except Exception:
            continue
        m_sys = unit_means(s_r, ii, jj, sys_vec, drop=DROP_SYS)
        m_soz = unit_means(s_r, ii, jj, soz_vec)
        for u in surr_sys:
            if u in m_sys: surr_sys[u][r] = m_sys[u]
        for u in surr_soz:
            if u in m_soz: surr_soz[u][r] = m_soz[u]

    rows = []
    for grp, obs_d, surr_d in [("system", obs_sys, surr_sys), ("soz", obs_soz, surr_soz)]:
        for u, o in obs_d.items():
            col = surr_d[u][np.isfinite(surr_d[u])]
            p = float((np.sum(col >= o) + 1) / (col.size + 1)) if col.size else np.nan
            rows.append(dict(patient=pat, band=band, grouping=grp, unit=u, obs=float(o),
                             surr_p50=float(np.nanmedian(col)) if col.size else np.nan, p=p))
    return rows


def cohort(df, grouping):
    """Per (band, unit): Wilcoxon(obs - surr_p50 over patients, greater) + BH per band."""
    out = []
    d = df[df.grouping == grouping]
    for band in BANDS:
        rows = []
        db = d[d.band == band]
        for u in sorted(db.unit.unique()):
            x = db[db.unit == u].dropna(subset=["obs", "surr_p50"])
            if len(x) < 5:
                continue
            diff = x.obs.values - x.surr_p50.values
            try:
                _, p = wilcoxon(diff, alternative="greater")
            except ValueError:
                p = 1.0
            rows.append(dict(band=band, unit=u, n=len(x), obs_med=float(x.obs.median()),
                             n_pos=int((diff > 0).sum()), p=float(p)))
        if rows:
            rdf = pd.DataFrame(rows)
            rdf["q"] = bh_fdr(rdf.p.values)
            out.append(rdf)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def main():
    global OUT
    OUT = ROOT / "data" / "sparsified_arc" / f"localization_mst020_s{S_REPORT:04.1f}"
    OUT.mkdir(parents=True, exist_ok=True)
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))  # warm numba
    jobs = [(i, pat, band) for i, (pat, band) in
            enumerate((p, b) for b in BANDS for p in COHORT)]
    print(f"[localization mst@0.20] s={S_REPORT:.3f}  R={R}  {len(jobs)} cells "
          f"({len(BANDS)} bands x {len(COHORT)} pat)", flush=True)
    t0 = time.time(); recs = []
    nproc = min(12, max(1, (os.cpu_count() or 2) - 2))
    with Pool(nproc) as pool:
        for k, rl in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if rl:
                recs.extend(rl)
            el = time.time() - t0
            print(f"  [{k}/{len(jobs)}] elapsed {el:5.0f}s  ETA {el/k*(len(jobs)-k):5.0f}s",
                  flush=True)
    per_pat = pd.DataFrame(recs)
    per_pat.to_csv(OUT / "per_patient.csv", index=False)
    sys_c = cohort(per_pat, "system"); sys_c.to_csv(OUT / "system_enrichment.csv", index=False)
    soz_c = cohort(per_pat, "soz");   soz_c.to_csv(OUT / "soz_enrichment.csv", index=False)

    print(f"\n=== SYSTEM enrichment (upper-tail, BH per band) — top hits ===", flush=True)
    if not sys_c.empty:
        for band in BANDS:
            b = sys_c[sys_c.band == band].sort_values("p").head(3)
            for _, r in b.iterrows():
                flag = " <== " if r.q < 0.05 else ("  ~  " if r.p < 0.05 else "     ")
                print(f"  {band:10s} {r.unit:22s} obs={r.obs_med:+.3f} n+={int(r.n_pos)}/"
                      f"{int(r.n)} p={r.p:.3f} q={r.q:.3f}{flag}", flush=True)
    print(f"\n=== SOZ enrichment (soz vs non_soz) ===", flush=True)
    if not soz_c.empty:
        for _, r in soz_c[soz_c.unit == "soz"].iterrows():
            flag = " <== SOZ" if r.p < 0.05 else ""
            print(f"  {r.band:10s} soz obs={r.obs_med:+.3f} n+={int(r.n_pos)}/{int(r.n)} "
                  f"p={r.p:.3f} q={r.q:.3f}{flag}", flush=True)
    print(f"\n[done] {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
