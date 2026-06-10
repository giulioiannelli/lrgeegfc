#!/usr/bin/env python3
"""Audit 95 — the decisive 'mark FEW, discover MANY?' test: recovery as a function
of the number of seed SOZ contacts.

audit_92 validated recovery at HIDE_FRAC=0.5 (mark ~half the SOZ, recover the other
half, AUC 0.72–0.82). The clinically exciting claim is stronger and untested: mark
only a FEW seeds (k=2,3,5,...) and discover the MANY remaining SOZ. That regime is
HARDER (fewer seeds → noisier community definition), so it must be measured, not
assumed. This script sweeps the seed count k and reports both a ranking metric (AUC
of remaining-SOZ vs healthy) and a hard discovery metric (precision@top-m, m = number
of remaining SOZ — i.e. if you flag the top-m contacts, how many are truly SOZ).

Critical preamble
=================
(1) Claim: from a FEW seed SOZ, the relational marker recovers the MANY remaining
    SOZ with useful AUC AND precision@top-m well above prevalence.
(2) Null: prevalence (precision of a random flag) and chance AUC 0.5; node-strength
    baseline.
(3) Strongest alternative: it only works with heavy seeding (k≈n_epi/2); at k=2–3 it
    collapses to chance / prevalence.
(4) Reach: scores are strength-residualised affinity-to-seeds (strength-orthogonal);
    seeds are drawn at random many times per k; precision@top-m is the honest
    "discover many" metric (not just a ranking AUC).
(5) Falsification: if AUC → 0.5 and precision@top-m → prevalence as k shrinks, then
    'mark few discover many' is FALSE — the marker needs substantial seeding.

Outputs (``data/audit/epi_marker_relational/``)
    seed_curve_per_patient.csv   per (patient, band, k)
    seed_curve_cohort.csv        per (band, k)
    README_seed_curve.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.patient import build_epi_masks

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI, N_TAU  # type: ignore
from audit_92_epi_masked_recovery_relational import (  # type: ignore
    _aff_to_set, _auc, _resid_on_strength)

OUT = ROOT / "data" / "audit" / "epi_marker_relational"
SEED_COUNTS = (2, 3, 5, 8)
N_DRAW = 150


def _precision_at_m(score, target_idx, healthy_idx, m):
    """Flag the top-m non-seed contacts; fraction that are true (remaining) SOZ."""
    cand = np.r_[target_idx, healthy_idx]
    s = score[cand]
    is_t = np.r_[np.ones(target_idx.size, bool), np.zeros(healthy_idx.size, bool)]
    ok = np.isfinite(s)
    cand, is_t = cand[ok], is_t[ok]
    s = s[ok]
    if s.size == 0 or m <= 0:
        return np.nan
    top = np.argsort(-s)[:m]
    return float(is_t[top].mean())


def per_patient_band(pat, band, rng, verbose):
    try:
        W = load_phase_fc(pat, "rest_post", band)
    except Exception:
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    if epi.sum() < MIN_EPI or (~epi).sum() < 3:
        return None
    lam_o, U_o = np.linalg.eigh(np.diag(W.sum(1)) - W)
    if float(lam_o[-1]) <= 0:
        return None
    taus = np.geomspace(1.0 / lam_o[-1], 10.0 / lam_o[-1], N_TAU)
    strength = W.sum(1)
    epi_idx = np.where(epi)[0]
    healthy = np.where(~epi)[0]
    n_epi = epi_idx.size
    rows = []
    for k in SEED_COUNTS:
        if k >= n_epi:
            continue
        aucs, precs, sprecs = [], [], []
        m = n_epi - k                      # remaining SOZ to discover
        prevalence = m / (m + healthy.size)
        for _ in range(N_DRAW):
            seeds = rng.choice(epi_idx, size=k, replace=False)
            targets = np.array([e for e in epi_idx if e not in set(seeds)])
            f_raw, f_res = _aff_to_set(lam_o, U_o, taus, seeds, strength)
            aucs.append(_auc(f_res, targets, healthy))
            precs.append(_precision_at_m(f_res, targets, healthy, m))
            sprecs.append(_precision_at_m(strength, targets, healthy, m))
        rows.append({
            "patient": pat, "band": band, "k_seeds": k, "n_epi": n_epi,
            "n_remaining": m, "prevalence": prevalence,
            "auc_resid": float(np.nanmean(aucs)),
            "precision_at_m": float(np.nanmean(precs)),
            "precision_strength": float(np.nanmean(sprecs)),
            "lift": float(np.nanmean(precs) / max(prevalence, 1e-9)),
        })
    return pd.DataFrame(rows)


def cohort(df):
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(df.band)]:
        for k in SEED_COUNTS:
            d = df[(df.band == band) & (df.k_seeds == k)]
            if d.empty:
                continue
            out.append({
                "band": band, "k_seeds": k, "n_patients": len(d),
                "med_auc_resid": float(d.auc_resid.median()),
                "med_precision_at_m": float(d.precision_at_m.median()),
                "med_prevalence": float(d.prevalence.median()),
                "med_lift": float(d.lift.median()),
                "med_precision_strength": float(d.precision_strength.median()),
            })
    return pd.DataFrame(out)


def write_readme(coh, runtime):
    L = ["---", "name: epi_seed_curve",
         "scope: mark_few_discover_many_seed_count_sweep",
         "era: COHORT_N10 / IMCOH_ABS",
         f"date: {time.strftime('%Y-%m-%d')}",
         "build_script: scripts/01_compute/audit/audit_95_epi_seed_curve.py",
         "---", "",
         "# 'Mark few, discover many?' — recovery vs number of seed SOZ", "",
         "**Head.** Mark k seed SOZ, score the rest by strength-residualised affinity "
         "to the seeds, measure recovery of the remaining SOZ. `auc_resid` = ranking "
         "(remaining SOZ vs healthy); `precision_at_m` = flag the top-m contacts "
         "(m = remaining SOZ), fraction truly SOZ; `lift` = precision / prevalence. "
         "The honest test of whether few seeds suffice.", "",
         "| band | k | AUC | precision@top-m | prevalence | lift | strength prec |",
         "|---|---|---|---|---|---|---|"]
    for band in [b for b in es.ALL_BANDS if b in set(coh.band)]:
        for k in SEED_COUNTS:
            r = coh[(coh.band == band) & (coh.k_seeds == k)]
            if r.empty:
                continue
            r = r.iloc[0]
            L.append(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {k} "
                     f"| {r.med_auc_resid:.3f} | {r.med_precision_at_m:.3f} "
                     f"| {r.med_prevalence:.3f} | {r.med_lift:.1f}× "
                     f"| {r.med_precision_strength:.3f} |")
    L += ["", "## Reading",
          "- If AUC stays high and precision@top-m stays well above prevalence as k "
          "drops to 2–3, 'mark few → discover many' holds. If they collapse toward "
          "chance/prevalence, the marker needs substantial seeding.",
          f"- {N_DRAW} random seed draws per k; wall-clock {runtime:.1f}s"]
    (OUT / "README_seed_curve.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=["delta", "beta", "low_gamma", "alpha"])
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    frames = []
    for band in args.bands:
        for pat in args.patients:
            rng = np.random.default_rng(20260610 + sum(ord(c) for c in pat + band))
            res = per_patient_band(pat, band, rng, args.verbose)
            if res is not None:
                frames.append(res)
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(OUT / "seed_curve_per_patient.csv", index=False)
    coh = cohort(df)
    coh.to_csv(OUT / "seed_curve_cohort.csv", index=False)
    runtime = time.time() - t0
    write_readme(coh, runtime)
    print("[audit_95] recovery vs #seeds (median): AUC | precision@top-m (prev) | lift")
    for band in [b for b in es.ALL_BANDS if b in set(coh.band)]:
        print(f"  {band}:")
        for k in SEED_COUNTS:
            r = coh[(coh.band == band) & (coh.k_seeds == k)]
            if r.empty:
                continue
            r = r.iloc[0]
            print(f"    k={k}: AUC={r.med_auc_resid:.2f} "
                  f"prec@m={r.med_precision_at_m:.2f} (prev {r.med_prevalence:.2f}, "
                  f"{r.med_lift:.1f}×) str_prec={r.med_precision_strength:.2f}")
    print(f"[audit_95] done in {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
