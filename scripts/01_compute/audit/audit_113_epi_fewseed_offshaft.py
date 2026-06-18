#!/usr/bin/env python3
"""Audit 113 — honest OFF-SHAFT few-seed reconstruction curve: mark only k SOZ contacts,
rank every contact on no-seed electrodes by diffusion affinity, recover the REMAINING
DISTANT SOZ. The "mark few, discover many" question, evaluated with proximity removed.

Replaces the withdrawn all-contacts few-seed curve (audit_95, proximity-confounded — a
trivial nearest-contact baseline beats it, audit_99). Marker = delta heat kernel
e^{-tau L} at slow tau (heat_t5), strength-residualised (the audit_101/102 winner).

This is the *few-seed* regime (k=2,3,5 seeds -> recover the many remaining), the
complement of audit_101's leave-one-shaft-out *many-seed* regime (all-but-one-shaft ->
recover one held-out shaft). The honest expectation is a MONOTONE seed-count curve that
climbs toward the ~0.72 leave-one-shaft-out value as k grows.

Critical preamble
=================
(1) Claim: from only k=2,3,5 known SOZ seeds, the propagator ranks the REMAINING SOZ on
    no-seed electrodes above healthy contacts (AUC > 0.5) and concentrates them in the
    top-m (lift > 1), beating a node-strength ranking.
(2) Null: few-seed AUC = 0.5 / lift = 1 (remaining SOZ not recoverable from k seeds).
(3) Strongest alternative: hubness (top contacts are strong) -> killed by the strength
    residual AND the explicit strength baseline; proximity (near a seed) -> killed by
    scoring ONLY contacts on electrodes with no seed (off-shaft, exactly as audit_101).
(4) Reach: k seeds drawn at random per patient (NDRAW draws), averaged; the off-shaft
    pool is rebuilt per draw; identical residual + AUC concordance as audit_101; the two
    hub-patients (Pat_10/15) are reported, not hidden.
(5) Falsification: if few-seed AUC ~ 0.5 cohort-wide, the marker cannot reconstruct from
    a few seeds and is reported flat. No spin.

Outputs (data/audit/epi_marker_fewseed/)
    fewseed_per_patient.csv ; fewseed_cohort.csv ; README.md
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

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI  # type: ignore
from audit_101_epi_marker_library import build_operators, _resid, _conc  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_marker_fewseed"
ML = ROOT / "data" / "audit" / "epi_marker_library"
BAND = "delta"
MARKER = "heat_t5"
KSEEDS = [2, 3, 5]
MIN_OFF_CTRL = 3        # need a few healthy controls off-shaft to score
RESP_THR = 0.568        # heat_t5 label-shuffle null p95 (audit_102) -> responder threshold
SEED = 20260618


def _auc(case, ctrl):
    c, t = _conc(np.asarray(case, float), np.asarray(ctrl, float))
    return c / t if t else np.nan


def per_patient(pat, ndraw, rng):
    try:
        W = load_phase_fc(pat, "rest_post", BAND)
    except Exception:
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    probes = np.asarray(pm.probes, object)
    epi_idx = np.where(epi)[0]
    if epi_idx.size < max(MIN_EPI, 4) or len(set(probes[epi_idx])) < 2:
        return None
    ops, _intr, strength = build_operators(W)
    O = ops[MARKER]
    epi_set = set(epi_idx.tolist())

    rows = []
    for k in KSEEDS:
        if epi_idx.size <= k:
            continue
        a_prop, a_str, lifts = [], [], []
        for _ in range(ndraw):
            seeds = rng.choice(epi_idx, size=k, replace=False)
            seed_pr = set(probes[seeds])
            seed_set = set(seeds.tolist())
            off = np.array([j for j in range(N) if probes[j] not in seed_pr])
            tgt = np.array([j for j in off if epi[j] and j not in seed_set])
            ctl = np.array([j for j in off if not epi[j]])
            if tgt.size < 1 or ctl.size < MIN_OFF_CTRL:
                continue
            aff = _resid(O[:, seeds].mean(axis=1), strength)        # strength-orthogonal
            a_prop.append(_auc(aff[tgt], aff[ctl]))
            a_str.append(_auc(strength[tgt], strength[ctl]))         # hubness baseline
            m = tgt.size
            top = off[np.argsort(-aff[off])][:m]                     # precision@top-m
            hit = sum(1 for j in top if j in set(tgt.tolist()))
            prev = tgt.size / off.size
            lifts.append((hit / m) / prev if prev > 0 else np.nan)
        if a_prop:
            rows.append(dict(
                patient=pat, k=k, n_draws=len(a_prop),
                auc_prop=float(np.nanmean(a_prop)),
                auc_strength=float(np.nanmean(a_str)),
                lift_topm=float(np.nanmean(lifts))))
    return pd.DataFrame(rows) if rows else None


def cohort(df, resp):
    out = []
    for k in KSEEDS:
        d = df[df.k == k]
        if d.empty:
            continue
        dr = d[d.patient.isin(resp)]
        out.append(dict(
            k=k, n_patients=int(d.patient.nunique()),
            med_auc_prop=float(d.auc_prop.median()),
            mean_auc_prop=float(d.auc_prop.mean()),
            n_above_half=int((d.auc_prop > 0.5).sum()),
            med_auc_strength=float(d.auc_strength.median()),
            med_lift_topm=float(d.lift_topm.median()),
            med_lift_topm_responders=float(dr.lift_topm.median()) if len(dr) else np.nan))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ndraw", type=int, default=300)
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    t0 = time.time()

    # responder set = audit_101 validated heat_t5 LOSO AUC >= RESP_THR (consistent definition)
    try:
        ml = pd.read_csv(ML / "marker_library_per_patient.csv")
        mld = ml[(ml.band == BAND) & (ml.marker == MARKER)]
        auc_val = dict(zip(mld.patient, mld.auc_resid))
    except Exception:
        auc_val = {}
    resp = [p for p in args.patients if auc_val.get(p, 0) >= RESP_THR]

    frames = [per_patient(p, args.ndraw, rng) for p in args.patients]
    df = pd.concat([f for f in frames if f is not None], ignore_index=True)
    df.to_csv(OUT / "fewseed_per_patient.csv", index=False)
    coh = cohort(df, resp)
    coh.to_csv(OUT / "fewseed_cohort.csv", index=False)

    print(f"[audit_113] {BAND} {MARKER} off-shaft few-seed ({args.ndraw} draws, "
          f"{df.patient.nunique()} patients, {time.time()-t0:.1f}s)\n")
    print("  k   propAUC(med)  n>0.5   strAUC(med)  lift@top-m(all/resp)")
    for _, r in coh.iterrows():
        print(f"  {int(r.k)}     {r.med_auc_prop:.3f}      {int(r.n_above_half)}/{int(r.n_patients)}    "
              f"{r.med_auc_strength:.3f}       {r.med_lift_topm:.2f} / {r.med_lift_topm_responders:.2f}")
    print("\n  per-patient propagator few-seed AUC:")
    piv = df.pivot(index="patient", columns="k", values="auc_prop")
    print(piv.round(3).to_string())
    print(f"\n[audit_113] -> {OUT}")


if __name__ == "__main__":
    main()
