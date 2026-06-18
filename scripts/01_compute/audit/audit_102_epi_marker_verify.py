#!/usr/bin/env python3
"""Audit 102 — verify the audit_101 winner: does the δ slow heat-kernel really mark
DISTANT epileptic nodes, or is it a 64-cell multiple-comparison fluke?

Winner (audit_101): heat kernel e^{-tL} at slow tau (t4/t5), δ band, strength-residual
leave-one-SHAFT-out distant-discovery AUC ~0.72-0.75, 7-8/10 patients. Slow diffusion
reaching distant nodes is the mechanistically PRE-SPECIFIED choice for distant discovery,
not cherry-picked — but we seal it three ways:

  (A) LABEL-SHUFFLE NULL: relabel n_epi random contacts (spanning >=2 shafts) as fake-SOZ,
      recompute the cohort-median LOSO AUC; 500 draws -> empirical p. If real >> null, the
      distant recovery is specific to the TRUE SOZ being a co-diffusing community.
  (B) NESTED LOPO marker selection: pick the best marker on 9 patients (by mean resid AUC),
      score it on the held-out 10th; the mean held-out AUC is selection-bias-free.
  (C) report it across the whole heat family + both Laplacian normalisations so it is a
      coherent ridge, not one lucky tau.

No spin: if the null is not cleared (p>=0.05) or LOPO collapses to 0.5, it is reported flat.
Outputs: data/audit/epi_marker_library/verify_{nulls,lopo}.csv + console verdict.
"""
from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.utils.io.patient import build_epi_masks

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI  # type: ignore
from audit_101_epi_marker_library import build_operators, loso_auc  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_marker_library"
N_SHUF = 500
HERO = ["heat_t4", "heat_t5", "heat_t3", "heatN_t4"]   # the δ winners to verify


def _prep(pat, band):
    try:
        W = load_phase_fc(pat, "rest_post", band)
    except Exception:
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    probes = np.asarray(pm.probes, object)
    epi_idx = np.where(epi)[0]
    healthy = np.where(~epi)[0]
    if epi_idx.size < max(MIN_EPI, 4) or healthy.size < 5:
        return None
    if len(set(probes[epi_idx])) < 2:
        return None
    ops, _intr, strength = build_operators(W)
    return dict(N=N, probes=probes, epi_idx=epi_idx, healthy=healthy,
                strength=strength, ops=ops)


def _auc_for(op, epi_idx, healthy, probes, strength):
    def marker(seeds, O=op):
        return O[:, seeds].mean(axis=1)
    _raw, res, _n = loso_auc(marker, epi_idx, healthy, probes, strength)
    return res


def main():
    band = "delta"
    t0 = time.time()
    P = {pat: _prep(pat, band) for pat in es.COHORT}
    P = {k: v for k, v in P.items() if v is not None}
    print(f"[audit_102] {band}: {len(P)} patients with >=2 SOZ shafts")

    # ---- real cohort AUC for each hero marker ----
    real = {}
    for name in HERO:
        v = [_auc_for(P[p]["ops"][name], P[p]["epi_idx"], P[p]["healthy"],
                      P[p]["probes"], P[p]["strength"]) for p in P]
        v = np.array([x for x in v if np.isfinite(x)])
        real[name] = v
        print(f"  REAL {name}: median={np.median(v):.3f} mean={np.mean(v):.3f} "
              f"{int((v>0.5).sum())}/{v.size}>0.5")

    # ---- (A) label-shuffle null on the median, per hero marker ----
    rng = np.random.default_rng(20260612)
    null_rows = []
    for name in HERO:
        real_med = float(np.median(real[name]))
        null_meds = []
        for _ in range(N_SHUF):
            per = []
            for p, d in P.items():
                N, probes = d["N"], d["probes"]
                n_epi = d["epi_idx"].size
                # fake SOZ: n_epi random contacts spanning >=2 shafts
                for _try in range(8):
                    fake = rng.choice(N, size=n_epi, replace=False)
                    if len(set(probes[fake])) >= 2:
                        break
                fepi = np.sort(fake)
                fhealthy = np.array([j for j in range(N) if j not in set(fepi.tolist())])
                a = _auc_for(d["ops"][name], fepi, fhealthy, probes, d["strength"])
                if np.isfinite(a):
                    per.append(a)
            if per:
                null_meds.append(float(np.median(per)))
        null_meds = np.array(null_meds)
        p_emp = float((null_meds >= real_med).mean())
        null_rows.append({"marker": name, "real_median_auc": real_med,
                          "null_median_mean": float(null_meds.mean()),
                          "null_median_p95": float(np.percentile(null_meds, 95)),
                          "p_empirical": p_emp, "n_shuf": null_meds.size})
        print(f"  NULL {name}: real_med={real_med:.3f}  null_med~{null_meds.mean():.3f} "
              f"(p95={np.percentile(null_meds,95):.3f})  p_emp={p_emp:.3f}")
    pd.DataFrame(null_rows).to_csv(OUT / "verify_nulls.csv", index=False)

    # ---- (B) nested LOPO marker selection over the full heat family ----
    fam = [m for m in P[next(iter(P))]["ops"] if m.startswith(("heat_t", "heatN_t"))]
    pats = list(P)
    auc_tab = {m: {p: _auc_for(P[p]["ops"][m], P[p]["epi_idx"], P[p]["healthy"],
                               P[p]["probes"], P[p]["strength"]) for p in pats}
               for m in fam}
    held = []
    for hp in pats:
        # select marker by mean resid-AUC on the OTHER patients
        best, bestval = None, -1
        for m in fam:
            vv = [auc_tab[m][q] for q in pats if q != hp and np.isfinite(auc_tab[m][q])]
            if vv and np.mean(vv) > bestval:
                bestval, best = np.mean(vv), m
        hv = auc_tab[best][hp]
        if np.isfinite(hv):
            held.append({"held_patient": hp, "selected_marker": best, "held_auc": hv})
    hd = pd.DataFrame(held)
    hd.to_csv(OUT / "verify_lopo.csv", index=False)
    print(f"\n[audit_102] NESTED-LOPO held-out AUC: median={hd.held_auc.median():.3f} "
          f"mean={hd.held_auc.mean():.3f}  {int((hd.held_auc>0.5).sum())}/{len(hd)}>0.5  "
          f"(selected: {hd.selected_marker.value_counts().to_dict()})")
    print(f"[audit_102] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
