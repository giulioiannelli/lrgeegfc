#!/usr/bin/env python3
"""Audit 98 — DEPLOYMENT pipeline: given a few labelled SOZ, rank the remaining
contacts and predict the other SOZ. This is the *deployment* view of the marker
validated in audit_91 (relational features) / audit_92 (masking recovery) /
audit_95 (few-seed curve) — it forks NO math (imports the validated helpers) and
adds only the human-readable prediction read-out (contact label, Desikan-Killiany
region, distance to nearest seed, score, HIT/miss).

It answers exactly: "can the propagator measure predict other SOZ nodes given a few
labelled?" — and SHOWS the prediction, not just an AUC.

Pipeline (per patient × band)
=============================
    1. W   = |ImCoh| FC, rest_post, band          (load_phase_fc)
    2. L   = D − W ;  λ,U = eigh(L)               (graph Laplacian eigos)
    3. τ   = geomspace(1/λmax, 10/λmax, N_TAU)    (multiscale diffusion times)
    4. seeds S  = k labelled SOZ contacts          (the "few labelled")
    5. f_aff(i) = mean_τ  ρ(i, S∖i)  with ρ=U e^{−τλ} Uᵀ / Z   (affinity to seeds)
    6. score(i) = f_aff(i) − linfit(strength→f_aff)(i)   (STRENGTH-ORTHOGONAL residual)
    7. rank non-seed contacts by score ↓ ; flag top-m (m = #remaining SOZ)
    8. read out predictions + precision@top-m, AUC(remaining-SOZ vs healthy), lift

Critical preamble
=================
(1) Claim: from k labelled SOZ, the strength-residual propagator-affinity score ranks
    the *remaining* SOZ above healthy contacts well enough to be a useful shortlist.
(2) Null: prevalence (random-flag precision) and chance AUC 0.5; node-strength baseline.
(3) Strongest alternative: the score is just re-reading hubness / spatial proximity to
    seeds (find-near-the-seed). Controlled by the strength residual (step 6) and by
    reporting distance-to-seed alongside each hit so proximity is visible, not hidden.
(4) Reach: this is SEED-BASED (semi-supervised) — it needs ≥ a few labelled SOZ to
    define S; it cannot run on a zero-label patient. precision@top-m is the honest hard
    metric; AUC is the ranking metric; both averaged over N_DRAW random seed draws.
(5) Falsification: if precision@top-m → prevalence and AUC → 0.5, the propagator adds
    nothing over chance and the marker is not predictive.

This is a demonstration / deployment script (no new statistical claim beyond audit_92/95);
it writes a per-(patient,band,k) prediction table for inspection.

Usage
-----
    python audit_98_epi_seed_prediction_pipeline.py --patients Pat_08 --bands delta --k 3 --show
    python audit_98_epi_seed_prediction_pipeline.py            # full cohort, default bands/k
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.utils.io.regions import load_channel_regions

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI, N_TAU  # type: ignore
from audit_92_epi_masked_recovery_relational import (  # type: ignore
    _aff_to_set, _auc, _resid_on_strength)

OUT = ROOT / "data" / "audit" / "epi_seed_prediction"
DEFAULT_BANDS = ["delta", "beta", "low_gamma"]
N_DRAW = 200


def _laplacian_eigs(W):
    lam, U = np.linalg.eigh(np.diag(W.sum(1)) - W)
    return lam, U


def _taus(lam):
    return np.geomspace(1.0 / lam[-1], 10.0 / lam[-1], N_TAU)


def predict_from_seeds(W, seeds, strength, lam=None, U=None):
    """Core deployment call: given the FC graph and a seed index array, return a
    strength-orthogonal SOZ-affinity score for every node (np.nan for seeds)."""
    if lam is None:
        lam, U = _laplacian_eigs(W)
    taus = _taus(lam)
    _f_raw, f_res = _aff_to_set(lam, U, taus, np.asarray(seeds), strength)
    score = f_res.copy()
    score[np.asarray(seeds)] = np.nan          # don't rank the seeds themselves
    return score


def _nearest_seed_mm(coords, i, seeds):
    d = np.linalg.norm(coords[seeds] - coords[i], axis=1)
    return float(np.nanmin(d)) / 1000.0        # implant coords are micrometres → mm


def showcase(pat, band, k, rng):
    """One concrete, human-readable prediction run on real seeds."""
    W = load_phase_fc(pat, "rest_post", band)
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N or pm.epi_mask.sum() < max(MIN_EPI, k + 1):
        return None
    reg = load_channel_regions(pat)
    coords = reg[["x", "y", "z"]].to_numpy(float)
    labels = reg["label_raw"].to_numpy()
    regions = reg["region"].to_numpy()
    epi = np.asarray(pm.epi_mask, bool)
    epi_idx = np.where(epi)[0]
    healthy = np.where(~epi)[0]
    strength = W.sum(1)
    lam, U = _laplacian_eigs(W)

    # Deterministic showcase: seeds = k highest-strength SOZ (the "clinically salient
    # onset contacts you already found"); targets = the remaining SOZ to discover.
    order = epi_idx[np.argsort(-strength[epi_idx])]
    seeds = order[:k]
    targets = np.array([e for e in epi_idx if e not in set(seeds)])
    m = targets.size
    score = predict_from_seeds(W, seeds, strength, lam, U)

    # rank all non-seed contacts; top-m are the prediction
    cand = np.r_[targets, healthy]
    cs = score[cand]
    ok = np.isfinite(cs)
    cand, cs = cand[ok], cs[ok]
    order_pred = cand[np.argsort(-cs)]
    topm = order_pred[:m]
    tset = set(targets.tolist())
    hits = [i for i in topm if i in tset]

    rows = []
    for rank, i in enumerate(order_pred[:max(m, 10)], 1):
        rows.append({
            "rank": rank, "label": labels[i], "region": regions[i],
            "dist_to_seed_mm": round(_nearest_seed_mm(coords, i, seeds), 1),
            "score": round(float(score[i]), 4),
            "is_SOZ": bool(i in tset), "in_topm": rank <= m,
        })
    show = pd.DataFrame(rows)
    auc = _auc(score, targets, healthy)
    prec = len(hits) / m if m else np.nan
    prevalence = m / (m + healthy.size)
    seed_lab = ", ".join(labels[seeds])
    return dict(pat=pat, band=band, k=k, n_epi=int(epi.sum()), m=m,
                seeds=seed_lab, auc=auc, precision=prec, prevalence=prevalence,
                lift=prec / max(prevalence, 1e-9), table=show,
                n_hits=len(hits), hit_labels=", ".join(labels[hits]))


def cohort_metrics(pat, band, k, rng):
    """Honest metrics: average over N_DRAW random seed draws (not strength-picked)."""
    W = load_phase_fc(pat, "rest_post", band)
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N or pm.epi_mask.sum() < max(MIN_EPI, k + 1):
        return None
    epi = np.asarray(pm.epi_mask, bool)
    epi_idx = np.where(epi)[0]
    healthy = np.where(~epi)[0]
    strength = W.sum(1)
    lam, U = _laplacian_eigs(W)
    taus = _taus(lam)
    aucs, precs, sprecs = [], [], []
    m = epi_idx.size - k
    prevalence = m / (m + healthy.size)
    for _ in range(N_DRAW):
        seeds = rng.choice(epi_idx, size=k, replace=False)
        targets = np.array([e for e in epi_idx if e not in set(seeds)])
        _fr, f_res = _aff_to_set(lam, U, taus, seeds, strength)
        # precision@top-m for the residual score and for the strength baseline
        for sc, bucket in ((f_res, precs), (strength, sprecs)):
            s = np.r_[sc[targets], sc[healthy]]
            is_t = np.r_[np.ones(m, bool), np.zeros(healthy.size, bool)]
            good = np.isfinite(s)
            s, is_t = s[good], is_t[good]
            bucket.append(float(is_t[np.argsort(-s)[:m]].mean()) if m else np.nan)
        aucs.append(_auc(f_res, targets, healthy))
    return dict(patient=pat, band=band, k=k, n_epi=int(epi.sum()), m=m,
                prevalence=prevalence, auc=float(np.nanmean(aucs)),
                precision_at_m=float(np.nanmean(precs)),
                precision_strength=float(np.nanmean(sprecs)),
                lift=float(np.nanmean(precs) / max(prevalence, 1e-9)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--bands", nargs="+", default=DEFAULT_BANDS)
    ap.add_argument("--k", type=int, default=3, help="number of labelled seed SOZ")
    ap.add_argument("--show", action="store_true",
                    help="print the concrete prediction table for each patient×band")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    rows = []
    showcases = []
    for band in args.bands:
        for pat in args.patients:
            rng = np.random.default_rng(20260611 + sum(ord(c) for c in pat + band))
            mc = cohort_metrics(pat, band, args.k, rng)
            if mc is not None:
                rows.append(mc)
            if args.show:
                sc = showcase(pat, band, args.k, rng)
                if sc is None:
                    continue
                tbl = sc["table"].copy()
                tbl.insert(0, "patient", pat)
                tbl.insert(1, "band", band)
                tbl.insert(2, "m", sc["m"])
                tbl.insert(3, "seeds", sc["seeds"])
                showcases.append(tbl)
                print(f"\n=== {pat} / {band} / k={args.k} seeds  "
                      f"(n_epi={sc['n_epi']}, discover m={sc['m']}) ===")
                print(f"  seeds (labelled SOZ given): {sc['seeds']}")
                print(f"  AUC={sc['auc']:.2f}  precision@top{sc['m']}={sc['precision']:.2f} "
                      f"(prevalence {sc['prevalence']:.2f}, {sc['lift']:.1f}× lift)")
                print(f"  HITS in top-{sc['m']}: {sc['n_hits']}/{sc['m']} -> {sc['hit_labels']}")
                with pd.option_context("display.max_rows", None, "display.width", 140):
                    print(sc["table"].to_string(index=False))
    if showcases:
        pd.concat(showcases, ignore_index=True).to_csv(
            OUT / f"seed_prediction_showcase_k{args.k}.csv", index=False)

    df = pd.DataFrame(rows)
    if not df.empty:
        df.to_csv(OUT / f"seed_prediction_per_patient_k{args.k}.csv", index=False)
        coh = (df.groupby("band")
                 .agg(n_patients=("patient", "size"),
                      med_auc=("auc", "median"),
                      med_precision_at_m=("precision_at_m", "median"),
                      med_prevalence=("prevalence", "median"),
                      med_lift=("lift", "median"),
                      med_precision_strength=("precision_strength", "median"))
                 .reindex([b for b in es.ALL_BANDS if b in set(df.band)]))
        coh.to_csv(OUT / f"seed_prediction_cohort_k{args.k}.csv")
        print(f"\n[audit_98] cohort median (k={args.k} seeds): "
              f"AUC | precision@top-m (prev) | lift | strength-prec")
        for band, r in coh.iterrows():
            print(f"  {band:10s}: AUC={r.med_auc:.2f}  "
                  f"prec={r.med_precision_at_m:.2f} (prev {r.med_prevalence:.2f}, "
                  f"{r.med_lift:.1f}×)  str_prec={r.med_precision_strength:.2f}")
    print(f"[audit_98] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
