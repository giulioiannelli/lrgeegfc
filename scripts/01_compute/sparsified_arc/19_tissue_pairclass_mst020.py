#!/usr/bin/env python3
"""Tissue / SOZ pair-class stratification of the trace on mst@0.20, scale-swept.

Ports audit_156 (WM + epi pair-class rho_sym) onto the recovered mst@0.20 backbone,
across the full dimensionless scale grid s = tau*lambda_max in [1,180]. Settles R1.9
(the beta trace is carried by GRAY-GRAY coupling, not white matter) and R1.11 (over the
seizure core beta SPARES while alpha RECRUITS), on the non-degenerate propagator and as
a function of scale.

Pair-class rho_sym restricts the symmetric split-half trace to pairs whose BOTH endpoints
fall in a tissue class (audit_156):

  rho_sym(class) = 1/2[ spearman((cT-cA)[m], (cP-cB)[m]) + spearman((cT-cB)[m], (cP-cA)[m]) ]

with m = the pair mask for the class. Classes:
  WM  : gray_gray (both region != 'Wm'), cross, wm_wm
  EPI : nonepi_nonepi (both kept by epi_keep_mask), cross_epi, epi_epi (both SOZ)

Matched-strength null IDENTICAL to 13/16/17 (4-cycle +/-delta over all 4 phases ->
sparsify@0.20 -> LRG), R=200, recomputed at every scale. Cohort Wilcoxon(obs-surr_p50,
greater) per (band, class, scale).

5-point preamble
1. Claim: on mst@0.20 the beta trace concentrates in gray_gray coupling and NOT in the
   seizure core (epi_epi spared), while the alpha trace RECRUITS epi_epi.
2. Null: the class contrast is a node-strength artifact (SOZ contacts are high-degree).
3. Strongest alternative it must beat: node strength -> IDENTICAL matched-strength
   surrogate (strength-preserving). Small-class pair-count is a caveat (fewer pairs ->
   noisier Spearman), reported as n_pairs.
4. Cannot: R=200; classes from a-priori WM (Desikan argmax 'Wm') + binary SOZ; global-
   rewire surrogate (not within-class). Depleted (lower-tail) reported separately.
5. Falsify: if beta gray_gray does not exceed the null, or beta epi_epi is not spared
   (>= null) while alpha epi_epi does not recruit, R1.9/R1.11 do not hold on mst@0.20.

Output: data/sparsified_arc/tissue_pairclass_mst020/{per_patient,cohort}.csv
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
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask   # noqa: E402
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr             # noqa: E402
from lrg_eegfc.utils.io.regions import load_channel_regions       # noqa: E402

PHASES = ("A", "B", "task_test", "rest_post")
BANDS = ["delta", "alpha", "beta", "low_gamma"]
CLASSES = ("gray_gray", "cross_wm", "wm_wm", "nonepi_nonepi", "cross_epi", "epi_epi")
FRAC = 0.20
SGRID = np.logspace(0.0, np.log10(180.0), 16)
R = 200
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260713
MIN_PAIRS = 15
LIMIT = int(os.environ.get("SA_LIMIT", "0"))


def _spear(a, b):
    """Spearman via centered-rank correlation (nan if degenerate)."""
    if a.size < 3:
        return np.nan
    ra, rb = rankdata(a), rankdata(b)
    ra = ra - ra.mean(); rb = rb - rb.mean()
    da, db = np.sqrt(np.dot(ra, ra)), np.sqrt(np.dot(rb, rb))
    return float(np.dot(ra, rb) / (da * db)) if da > 0 and db > 0 else np.nan


def sym_rho_class(cA, cB, cT, cP, mask):
    """Class-restricted symmetric split-half trace (audit_156)."""
    if mask.sum() < MIN_PAIRS:
        return np.nan
    t1, p1 = (cT - cA)[mask], (cP - cB)[mask]
    t2, p2 = (cT - cB)[mask], (cP - cA)[mask]
    return 0.5 * (_spear(t1, p1) + _spear(t2, p2))


def pair_masks(N, wm_node, epi_node):
    iu = np.triu_indices(N, 1); i, j = iu[0], iu[1]
    wi, wj = wm_node[i], wm_node[j]
    ei, ej = epi_node[i], epi_node[j]
    return {
        "gray_gray": (~wi) & (~wj),
        "cross_wm": (wi ^ wj),
        "wm_wm": wi & wj,
        "nonepi_nonepi": (~ei) & (~ej),
        "cross_epi": (ei ^ ej),
        "epi_epi": ei & ej,
    }


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
    wm_node = (rdf["region"].to_numpy().astype(str) == "Wm")
    epi_node = ~epi_keep_mask(rdf, pat)          # True = SOZ
    masks = pair_masks(N, wm_node, epi_node)
    npairs = {c: int(m.sum()) for c, m in masks.items()}
    NS = len(SGRID)

    def rho_over_scales(eig4):
        out = {c: np.full(NS, np.nan) for c in CLASSES}
        for si, s in enumerate(SGRID):
            C = {ph: cophenetic_at_scale(*eig4[ph], s) for ph in PHASES}
            cA, cB, cT, cP = C["A"], C["B"], C["task_test"], C["rest_post"]
            for c in CLASSES:
                out[c][si] = sym_rho_class(cA, cB, cT, cP, masks[c])
        return out

    eig_obs = {ph: eig_backbone(Ws[ph]) for ph in PHASES}
    obs = rho_over_scales(eig_obs)

    surr = {c: np.full((R, NS), np.nan) for c in CLASSES}
    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    for r in range(R):
        eig_sh = {ph: eig_backbone(matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX))
                  for ph in PHASES}
        try:
            so = rho_over_scales(eig_sh)
        except Exception:
            continue
        for c in CLASSES:
            surr[c][r] = so[c]

    rows = []
    for c in CLASSES:
        for si in range(NS):
            o = obs[c][si]
            col = surr[c][:, si]; col = col[np.isfinite(col)]
            p_up = float((np.sum(col >= o) + 1) / (col.size + 1)) if (col.size and np.isfinite(o)) else np.nan
            p_lo = float((np.sum(col <= o) + 1) / (col.size + 1)) if (col.size and np.isfinite(o)) else np.nan
            rows.append(dict(patient=pat, band=band, pair_class=c, scale_idx=si,
                             s=float(SGRID[si]), n_pairs=npairs[c], obs=float(o) if np.isfinite(o) else np.nan,
                             surr_p50=float(np.nanmedian(col)) if col.size else np.nan,
                             p_upper=p_up, p_lower=p_lo))
    return rows


def cohort(df):
    out = []
    for band in BANDS:
        for c in CLASSES:
            for si in sorted(df.scale_idx.unique()):
                x = df[(df.band == band) & (df.pair_class == c) & (df.scale_idx == si)]
                x = x.dropna(subset=["obs", "surr_p50"])
                if len(x) < 5:
                    continue
                diff = x.obs.values - x.surr_p50.values
                try:
                    _, p_up = wilcoxon(diff, alternative="greater")
                except ValueError:
                    p_up = 1.0
                try:
                    _, p_lo = wilcoxon(diff, alternative="less")
                except ValueError:
                    p_lo = 1.0
                out.append(dict(band=band, pair_class=c, scale_idx=int(si), s=float(x.s.iloc[0]),
                                n=len(x), obs_med=float(x.obs.median()),
                                surr_med=float(x.surr_p50.median()), n_pos=int((diff > 0).sum()),
                                p_upper=float(p_up), p_lower=float(p_lo),
                                n_pairs_med=float(x.n_pairs.median())))
    cdf = pd.DataFrame(out)
    if not cdf.empty:
        for band in BANDS:
            m = cdf.band == band
            cdf.loc[m, "q_upper"] = bh_fdr(cdf.loc[m, "p_upper"].values)
    return cdf


def _summary(cdf):
    print("\n=== TISSUE / SOZ pair-class trace on mst@0.20 (upper-tail = concentrated) ===", flush=True)
    for band in BANDS:
        print(f"\n[{band}]", flush=True)
        b = cdf[cdf.band == band]
        for c in CLASSES:
            bc = b[b.pair_class == c]
            if bc.empty:
                continue
            best = bc.loc[bc.p_upper.idxmin()]
            nclear = int((bc.q_upper < 0.05).sum())
            depl = int((bc.p_lower < 0.05).sum())
            print(f"  {c:16s} best p_up={best.p_upper:.3f} q={best.q_upper:.3f} @s{best.s:5.1f} "
                  f"(obs {best.obs_med:+.3f} vs surr {best.surr_med:+.3f}, {best.n_pos}/{best.n}); "
                  f"BH-clear {nclear}/{len(bc)} scales; depleted {depl} | ~{best.n_pairs_med:.0f} pairs",
                  flush=True)


def main():
    OUT = ROOT / "data" / "sparsified_arc" / "tissue_pairclass_mst020"
    OUT.mkdir(parents=True, exist_ok=True)
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    jobs = [(i, pat, band) for i, (pat, band) in
            enumerate((p, b) for b in BANDS for p in COHORT)]
    if LIMIT:
        jobs = jobs[:LIMIT]
    print(f"[tissue pair-class mst@0.20] R={R} {len(SGRID)} scales {len(jobs)} cells", flush=True)
    t0 = time.time(); recs = []
    nproc = min(12, max(1, (os.cpu_count() or 2) - 2))
    with Pool(nproc) as pool:
        for k, rl in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if rl:
                recs.extend(rl)
            el = time.time() - t0
            print(f"  [{k}/{len(jobs)}] elapsed {el:6.0f}s ETA {el/k*(len(jobs)-k):6.0f}s", flush=True)
    per_pat = pd.DataFrame(recs); per_pat.to_csv(OUT / "per_patient.csv", index=False)
    cdf = cohort(per_pat); cdf.to_csv(OUT / "cohort.csv", index=False)
    if not cdf.empty:
        _summary(cdf)
    print(f"\n[done] {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
