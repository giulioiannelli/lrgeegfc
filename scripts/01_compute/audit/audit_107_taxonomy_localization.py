#!/usr/bin/env python3
"""audit_107 — cross-phase taxonomy localization (geometry baseline).

STEP 3 of the cross-phase taxonomy (scope 2026-06-12): WHERE do the channels
live? Localizes each per-pair channel signal to a-priori DK anatomical systems
and tests enrichment against the IMPLANT-GEOMETRY baseline (count-matched random
node-label permutation) — the spatial null chosen for the similarity channels
(the coupled temporal null failed its validation gate, audit_106).

Channel signals localized (per pair, from lrg_eegfc.utils.metrics.cross_phase):
  trace      = concordance(D_task-D_preA, D_post-D_preB)   [audit_83 signal;
               REFERENCE — must reproduce beta -> OFC]
  persist    = phi1**2     (persistence energy = trace carrier)
  reset      = phi2**2     (excursion energy)
  reset_dom  = phi2**2 - phi1**2   (excursion minus persistence = reset-specific;
               positive where movement reverts rather than persists)
  anchor     = -var        (rigidity: where the tree barely moves)

Geometry baseline null: per patient, permute the node->system labels R times
(count-matched random draw, KC audit_60 style), recompute the per-system demeaned
endpoint mean, build the cohort-median null distribution, one-sided enrichment p,
BH-FDR over systems WITHIN (band, signal). Reproduces the trace->OFC sanity check
on the SAME pipeline so reset/anchor numbers are trustworthy. Optional
--shaft-collapse robustness (audit_83 lesson: distributed claims need it).

NO tree-cutting; cophenetic built with the canonical _canon_cophenet shared with
audit_83.

Usage:
    python audit_107_taxonomy_localization.py --bands beta alpha low_gamma --R 500
    python audit_107_taxonomy_localization.py --bands beta --shaft-collapse
Outputs data/audit/cross_phase_taxonomy/localization{suffix}.csv
"""
from __future__ import annotations

import argparse
import importlib.util

import numpy as np
import pandas as pd

from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.cross_phase import phase_contrasts
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask, shaft_of
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
_spec = importlib.util.spec_from_file_location(
    "audit_83", ROOT / "scripts" / "01_compute" / "audit"
    / "audit_83_localization_matched_strength.py")
a83 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a83)
_spec63 = importlib.util.spec_from_file_location(
    "audit_63", ROOT / "scripts" / "01_compute" / "audit"
    / "audit_63_split_baseline_surrogate.py")
a63 = importlib.util.module_from_spec(_spec63)
_spec63.loader.exec_module(a63)

OUT = ROOT / "data" / "audit" / "cross_phase_taxonomy"
PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
GRANS = ["system", "supersystem"]
SEED = 20260612


def channel_signals(DA, DB, DT, DP) -> dict:
    """Per-pair channel signals for localization (all aligned, condensed)."""
    phi1, phi2, var = phase_contrasts(DA, DT, DP)
    return {
        "trace": a83.concordance(DT - DA, DP - DB),
        "persist": phi1 ** 2,
        "reset": phi2 ** 2,
        "reset_dom": phi2 ** 2 - phi1 ** 2,
        "anchor": -var,
    }


def geom_null_unit_means(s, iu_i, iu_j, unit_vec, keep, drop, shaft_vec,
                         R, rng) -> dict:
    """R count-matched node-label permutations -> per-unit means (dict unit->array)."""
    acc: dict[str, list] = {}
    n = unit_vec.size
    for _ in range(R):
        perm = rng.permutation(n)
        uv_perm = unit_vec[perm]
        sv_perm = shaft_vec[perm] if shaft_vec is not None else None
        um = a83.unit_means_from_s(s, iu_i, iu_j, uv_perm, keep, drop,
                                   shaft_vec=sv_perm)
        for u, m in um.items():
            acc.setdefault(u, []).append(m)
    return {u: np.array(v) for u, v in acc.items()}


def run_band(band, epi_excl, R, shaft_collapse, verbose=True):
    rng = np.random.default_rng(SEED)
    sig_names = ["trace", "persist", "reset", "reset_dom", "anchor"]
    obs = {(sg, g): {} for sg in sig_names for g in GRANS}      # -> {unit:{pat:m}}
    null = {(sg, g): {} for sg in sig_names for g in GRANS}     # -> {unit:{pat:arr}}
    for pat in a63.COHORT:
        try:
            Ws = {ph: a63.load_phase_fc(pat, ph, band) for ph in PHASES}
        except Exception as e:
            if verbose:
                print(f"  [skip] {pat} {band}: {e}")
            continue
        D = {ph: a83._canon_cophenet(Ws[ph]) for ph in PHASES}
        sigs = channel_signals(D["rest_pre_A"], D["rest_pre_B"],
                               D["task_test"], D["rest_post"])
        N = Ws[PHASES[0]].shape[0]
        iu_i, iu_j = np.triu_indices(N, k=1)
        rdf = load_channel_regions(pat)
        keep = epi_keep_mask(rdf, pat) if epi_excl else None
        shaft_vec = (np.array([shaft_of(l) for l in rdf["label_raw"]])
                     if shaft_collapse else None)
        for g in GRANS:
            uv = rdf[g].to_numpy().astype(str)
            for sg in sig_names:
                s = sigs[sg]
                om = a83.unit_means_from_s(s, iu_i, iu_j, uv, keep, a83.DROP[g],
                                           shaft_vec=shaft_vec)
                for u, m in om.items():
                    obs[(sg, g)].setdefault(u, {})[pat] = m
                nm = geom_null_unit_means(s, iu_i, iu_j, uv, keep, a83.DROP[g],
                                          shaft_vec, R, rng)
                for u in om:
                    if u in nm:
                        null[(sg, g)].setdefault(u, {})[pat] = nm[u]

    rows = []
    for (sg, g), pm_all in obs.items():
        recs = []
        for u, pm in pm_all.items():
            pats = sorted(pm)
            K = len(pats)
            if K == 0:
                continue
            M_obs = float(np.median([pm[p] for p in pats]))
            arrs = [null[(sg, g)].get(u, {}).get(p) for p in pats]
            arrs = [a for a in arrs if a is not None]
            if len(arrs) < K:
                continue
            Lmin = min(len(a) for a in arrs)
            Rmat = np.vstack([a[:Lmin] for a in arrs])
            M_null = np.median(Rmat, axis=0)
            p = (1 + int(np.sum(M_null >= M_obs))) / (M_null.size + 1)
            recs.append(dict(band=band, epi="exclude" if epi_excl else "include",
                             signal=sg, granularity=g, unit=u, K_implanted=K,
                             M_obs=M_obs, null_median=float(np.median(M_null)),
                             p_geom=p))
        if recs:
            q = bh_fdr(np.array([r["p_geom"] for r in recs]))
            for r, qq in zip(recs, q):
                r["q_geom"] = float(qq)
            rows += recs
    if verbose:
        for sg in sig_names:
            sub = [r for r in rows if r["signal"] == sg and r["granularity"] == "system"]
            sub = sorted(sub, key=lambda r: r["p_geom"])[:3]
            top = ", ".join(f"{r['unit']}(p={r['p_geom']:.3f},q={r['q_geom']:.3f})"
                            for r in sub)
            print(f"  {band:10s} {sg:10s} top systems: {top}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=["beta", "alpha", "low_gamma"])
    ap.add_argument("--epi-mode", choices=["include", "exclude"], default="include")
    ap.add_argument("--R", type=int, default=500)
    ap.add_argument("--shaft-collapse", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    epi_excl = args.epi_mode == "exclude"
    suffix = ("_shaftcollapsed" if args.shaft_collapse else "") \
        + (f"_{args.epi_mode}")

    print(f"\n==== taxonomy localization (geometry baseline) : epi-{args.epi_mode}"
          f"{' SHAFT-COLLAPSED' if args.shaft_collapse else ''}  R={args.R} ====")
    for pat in a63.COHORT:
        a63.ensure_half_fcs(pat, args.bands)
    allrows = []
    for band in args.bands:
        allrows += run_band(band, epi_excl, args.R, args.shaft_collapse)
    df = pd.DataFrame(allrows)
    out = OUT / f"localization{suffix}.csv"
    df.to_csv(out, index=False)
    print(f"\nOutputs -> {out}")


if __name__ == "__main__":
    main()
