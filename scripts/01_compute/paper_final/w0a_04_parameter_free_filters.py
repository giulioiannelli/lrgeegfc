#!/usr/bin/env python3
"""W0-A / A3 — can a genuinely parameter-free filter replace the free fraction?

Two questions the A2 sweep cannot answer on its own:

(i)  PMFG is the principled parent of TMFG but costs ~10 s per graph, so it
     cannot sit inside a 200-surrogate null. The standing recommendation
     (sparsification-choice.md) is to verify PMFG ~ TMFG on the OBSERVED
     readout and then use TMFG as the null-compatible proxy. That verification
     has never been run. It is run here, per-scale, on all four cross-phase
     functionals -- not on a single scalar.

(ii) The 2026-07-15 suite reported that planarity KILLS the alpha trace
     (12/16 -> 0/16). That is re-derived fresh here rather than assumed, and if
     it holds, the MECHANISM is characterised: planarity caps the edge budget at
     3N-6, so the question is which edges it must discard. We measure, at
     matched edge budget, what TMFG keeps that a plain weight-ranked filter of
     the SAME size does not -- i.e. whether the loss is the budget (density) or
     the topological prior (planarity) itself. Those two are confounded in any
     bare TMFG-vs-mst@0.20 comparison and the earlier verdict did not separate
     them.

Also reports the percolation backbone's density DISTRIBUTION across the cohort:
percolation is parameter-free only in the sense that the data sets the
threshold, so if its density swings wildly across patients and bands it cannot
sit inside a single plateau, and that has to be measured, not assumed.

5-POINT CRITICAL PREAMBLE
1. Claim. A parameter-free filter (TMFG/PMFG at 3N-6, or percolation at the
   connectivity bottleneck) can serve as the substrate, and PMFG ~ TMFG so the
   fast chordal computation is not itself the effect.
2. Null. H0_a: PMFG and TMFG give different cross-phase readouts, so any TMFG
   result is a greedy-triangulation artifact rather than the planar-filter
   principle. H0_b: TMFG's behaviour relative to mst@f is explained entirely by
   its edge BUDGET (3N-6 ~ 5% density), with planarity contributing nothing.
3. Strongest plausible alternative. That "planarity kills alpha" is really
   "5% density kills alpha" -- the planar filter is simply sparser than the
   incumbent 20%, and the density confound has been read as a topology result.
4. Does the null control for it, by mechanism? H0_b is controlled by matching
   the edge budget exactly: mst-union is evaluated at the fraction that
   reproduces TMFG's own density per (patient, band), so the two graphs differ
   in WHICH edges they keep and not in HOW MANY. What this CANNOT reject: the
   two filters also differ in connectivity policy (both are spanning) and in
   the local vs global nature of the selection, so a residual difference is
   "planar selection vs weight-ranked selection at matched size", not planarity
   in isolation. PMFG vs TMFG is observed-only -- there is no surrogate null on
   the PMFG arm, so it establishes agreement of the READOUT, not that PMFG
   clears a null.
5. Falsification / limits. (i) is falsified if the per-scale PMFG and TMFG
   functionals disagree beyond the scale-to-scale noise. (ii) is falsified if
   TMFG and a budget-matched mst-union give the same verdict, which would mean
   the planarity prior contributes nothing and the earlier "planarity kills
   alpha" statement was a density statement. n=10 and observed-only: this
   section describes graphs and readouts, it does not gate a claim.

Outputs (data/paper_final/w0a_substrate/a3_parameter_free/):
  pmfg_vs_tmfg.csv        per patient/band/scale/functional, both filters
  budget_matched.csv      TMFG vs mst-union at TMFG's own density, per cell
  edge_composition.csv    what planarity discards (rank/weight profile of the
                          edges TMFG drops that the budget-matched filter keeps)
  percolation_density.csv per patient/band density of the percolation backbone
"""
from __future__ import annotations

import json
import os
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, PATIENTS_4PHASE
from lrg_eegfc.utils.fc.backbone import (
    backbone_density, backbone_structure, mst_union_top_fraction,
    percolation_backbone, pmfg_backbone, tmfg_backbone,
)
from lrg_eegfc.utils.fc.heat_multiscale import (
    cross_phase_functionals_over_scales, laplacian_eig,
)
from lrg_eegfc.workflow.substrate import (
    CANONICAL_PHASES, CANONICAL_SCALES, canonical_graph, canonical_phase_graphs,
)

COHORT = list(PATIENTS_4PHASE)
# Default band set REDUCED to the two signal + two null bands under the machine
# resource cap: PMFG is O(N^3) planarity testing (~10 s/graph x 5 phases), so the
# full six-band grid is ~1.1 h at 4 workers. delta and high_gamma are dropped
# from THIS observed-only section only; the A2 stability sweep keeps all six.
DEFAULT_BANDS = ("theta", "alpha", "beta", "low_gamma")
BANDS = os.environ.get("W0A_BANDS", "").split(",") if os.environ.get("W0A_BANDS") \
    else list(DEFAULT_BANDS)
TRANSFORM = os.environ.get("W0A_TRANSFORM", "abs")
FUNCTIONALS = ("T_probe", "T_encode", "T_probespec", "T_probespec_pe")
OUT = ROOT / "data" / "paper_final" / "w0a_substrate" / "a3_parameter_free"


def _eigs(graphs):
    return {ph: laplacian_eig(W) for ph, W in graphs.items()}


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = canonical_phase_graphs(pat, band, dense=True, transform=TRANSFORM)
    except Exception as exc:                                       # pragma: no cover
        print(f"  [{pat}/{band}] load failed: {exc}", flush=True)
        return None
    N = Ws["A"].shape[0]
    s = CANONICAL_SCALES

    B_tmfg = {ph: tmfg_backbone(W) for ph, W in Ws.items()}
    B_pmfg = {ph: pmfg_backbone(W) for ph, W in Ws.items()}
    # budget-matched weight-ranked filter: mst-union at TMFG's OWN density,
    # per phase, so the comparison isolates WHICH edges, not HOW MANY.
    f_match = {ph: backbone_density(B_tmfg[ph]) for ph in Ws}
    B_match = {ph: mst_union_top_fraction(Ws[ph], f_match[ph]) for ph in Ws}
    B_perc = {ph: percolation_backbone(Ws[ph])[0] for ph in Ws}

    fun = {name: cross_phase_functionals_over_scales(_eigs(B), s)
           for name, B in (("tmfg", B_tmfg), ("pmfg", B_pmfg),
                           ("match", B_match), ("perc", B_perc))}

    rows_f, rows_d, rows_e = [], [], []
    for j, sv in enumerate(s):
        for k in FUNCTIONALS:
            rows_f.append(dict(patient=pat, band=band, s=float(sv), functional=k,
                               tmfg=float(fun["tmfg"][k][j]), pmfg=float(fun["pmfg"][k][j]),
                               budget_matched=float(fun["match"][k][j]),
                               perc=float(fun["perc"][k][j])))
    for ph in Ws:
        st_t = backbone_structure(B_tmfg[ph]); st_p = backbone_structure(B_pmfg[ph])
        st_m = backbone_structure(B_match[ph]); st_c = backbone_structure(B_perc[ph])
        rows_d.append(dict(patient=pat, band=band, phase=ph, N=int(N),
                           tmfg_density=st_t["density"], pmfg_density=st_p["density"],
                           match_density=st_m["density"], perc_density=st_c["density"],
                           tmfg_clustering=st_t["clustering"], match_clustering=st_m["clustering"],
                           perc_clustering=st_c["clustering"],
                           tmfg_gap=st_t["spectral_gap"], match_gap=st_m["spectral_gap"],
                           perc_gap=st_c["spectral_gap"],
                           tmfg_wfrac=float(B_tmfg[ph].sum() / Ws[ph].sum()),
                           match_wfrac=float(B_match[ph].sum() / Ws[ph].sum()),
                           perc_wfrac=float(B_perc[ph].sum() / Ws[ph].sum())))

        # --- what planarity discards, at matched budget --------------------- #
        r, c = np.triu_indices(N, 1)
        w = Ws[ph][r, c]
        rank = np.argsort(np.argsort(-w)) / w.size          # 0 = strongest edge
        et = B_tmfg[ph][r, c] > 0
        em = B_match[ph][r, c] > 0
        only_match = em & ~et                                # kept by weight, dropped by planarity
        only_tmfg = et & ~em                                 # kept by planarity, not by weight
        rows_e.append(dict(patient=pat, band=band, phase=ph,
                           jaccard=float((et & em).sum() / max(1, (et | em).sum())),
                           n_only_match=int(only_match.sum()), n_only_tmfg=int(only_tmfg.sum()),
                           rank_only_match=float(np.mean(rank[only_match])) if only_match.any() else np.nan,
                           rank_only_tmfg=float(np.mean(rank[only_tmfg])) if only_tmfg.any() else np.nan,
                           w_only_match=float(np.mean(w[only_match])) if only_match.any() else np.nan,
                           w_only_tmfg=float(np.mean(w[only_tmfg])) if only_tmfg.any() else np.nan,
                           w_med_all=float(np.median(w))))
    return rows_f, rows_d, rows_e


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Capped at 4 (shared box; see w0a_02 for the resource note).
    ncpu = int(os.environ.get("SA_WORKERS", 4))
    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in BANDS for p in COHORT)]
    print(f"[a3] {len(jobs)} cells (PMFG is O(N^3) planarity testing: ~10 s/graph x "
          f"{len(CANONICAL_PHASES)} phases/cell), {ncpu} workers -> {OUT}", flush=True)
    t0 = time.time()
    F, D, E = [], [], []
    with Pool(ncpu) as pool:
        for i, res in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if res:
                F.extend(res[0]); D.extend(res[1]); E.extend(res[2])
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s", flush=True)
    dfF, dfD, dfE = pd.DataFrame(F), pd.DataFrame(D), pd.DataFrame(E)
    dfF.to_csv(OUT / "pmfg_vs_tmfg.csv", index=False)
    dfD.to_csv(OUT / "backbone_densities.csv", index=False)
    dfE.to_csv(OUT / "edge_composition.csv", index=False)

    print(f"\n=== (i) PMFG vs TMFG on the OBSERVED readout, per functional (all "
          f"{dfF.patient.nunique()} patients x {dfF.band.nunique()} bands x 16 scales) ===",
          flush=True)
    for k in FUNCTIONALS:
        x = dfF[dfF.functional == k].dropna(subset=["tmfg", "pmfg"])
        if x.empty:
            continue
        d = (x.tmfg - x.pmfg).abs()
        rr = float(np.corrcoef(x.tmfg, x.pmfg)[0, 1])
        print(f"  {k:15s} pearson(tmfg,pmfg)={rr:.4f}  median|diff|={d.median():.4f}  "
              f"p95|diff|={d.quantile(.95):.4f}  sd(tmfg)={x.tmfg.std():.4f}", flush=True)

    print(f"\n=== (ii) planarity vs a BUDGET-MATCHED weight-ranked filter "
          f"(mst-union at TMFG's own density) ===", flush=True)
    print(f"  TMFG density: median {dfD.tmfg_density.median():.4f} "
          f"[{dfD.tmfg_density.min():.4f}, {dfD.tmfg_density.max():.4f}]   "
          f"budget-matched density: median {dfD.match_density.median():.4f}", flush=True)
    for k in FUNCTIONALS:
        x = dfF[dfF.functional == k].dropna(subset=["tmfg", "budget_matched"])
        if x.empty:
            continue
        print(f"  {k:15s} per-band median obs (tmfg | budget-matched | perc):", flush=True)
        for b in BANDS:
            y = x[x.band == b]
            if y.empty:
                continue
            print(f"    {b:11s} {y.tmfg.median():+.4f} | {y.budget_matched.median():+.4f} | "
                  f"{y.perc.median():+.4f}   corr(tmfg,matched)="
                  f"{np.corrcoef(y.tmfg, y.budget_matched)[0,1]:.3f}", flush=True)

    print(f"\n=== what planarity discards at matched budget ===", flush=True)
    print(f"  edge-set Jaccard(TMFG, budget-matched) = {dfE.jaccard.median():.3f} "
          f"[{dfE.jaccard.min():.3f}, {dfE.jaccard.max():.3f}]", flush=True)
    print(f"  edges kept by WEIGHT but dropped by planarity: n={dfE.n_only_match.median():.0f}/phase, "
          f"mean normalised weight-rank {dfE.rank_only_match.median():.4f} "
          f"(0 = strongest edge), mean weight {dfE.w_only_match.median():.4f}", flush=True)
    print(f"  edges kept by PLANARITY but not by weight: n={dfE.n_only_tmfg.median():.0f}/phase, "
          f"mean normalised weight-rank {dfE.rank_only_tmfg.median():.4f}, "
          f"mean weight {dfE.w_only_tmfg.median():.4f}  "
          f"(median edge weight overall {dfE.w_med_all.median():.4f})", flush=True)

    print(f"\n=== percolation backbone: is its 'parameter-free' density stable "
          f"across the cohort? ===", flush=True)
    print(f"  overall median {dfD.perc_density.median():.4f}, "
          f"range [{dfD.perc_density.min():.4f}, {dfD.perc_density.max():.4f}], "
          f"IQR [{dfD.perc_density.quantile(.25):.4f}, {dfD.perc_density.quantile(.75):.4f}]",
          flush=True)
    for b in BANDS:
        y = dfD[dfD.band == b]
        print(f"    {b:11s} median {y.perc_density.median():.4f}  "
              f"range [{y.perc_density.min():.4f}, {y.perc_density.max():.4f}]", flush=True)
    dfD[["patient", "band", "phase", "perc_density", "tmfg_density"]].to_csv(
        OUT / "percolation_density.csv", index=False)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="a3_parameter_free", transform=TRANSFORM, cohort=COHORT,
        bands=BANDS, phases=list(CANONICAL_PHASES),
        s_grid=[float(x) for x in CANONICAL_SCALES],
        note="observed-only (no surrogate null on the PMFG arm; PMFG is O(N^3))",
    ), indent=2))
    print(f"\n[a3] done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
