#!/usr/bin/env python3
"""Audit 103 — turn the verified delta distant-SOZ marker into a DEPLOYABLE tool and report
its HONEST precision: (A) leave-one-shaft-out precision@k / recall@k / lift for recovering
hidden SOZ, and (B) per-patient occult-candidate shortlists (rank the UNMARKED contacts by
distant diffusion affinity to the marked SOZ).

Replaces the withdrawn proximity-confounded deployment numbers (audit_98) with the clean
off-shaft version. Marker = delta heat kernel e^{-tau L} at slow tau (heat_t5),
strength-residualised (the audit_101/102 winner).

Critical preamble
=================
(1) Claim: among contacts on shafts FAR from the seeds, the top-k by marker score are
    enriched in true SOZ several-fold over the off-shaft base rate (lift > 1), beating a
    node-strength ranking.
(2) Null: lift@k = 1 (top-k no richer than a random off-shaft contact). Strength baseline
    ranks the same off pool by node strength.
(3) Strongest alternative: hubness (top contacts are just strong) -> killed by the strength
    residual AND shown explicitly via the strength-baseline ranking; proximity -> killed by
    the off-shaft restriction (targets and controls all far from every seed).
(4) Reach: precision@k is per leave-one-shaft-out fold (relevant = hidden SOZ on the held-out
    shaft), averaged per patient; lift normalises the small-targets-per-fold precision cap;
    the two hub-patients (Pat_10/15) are reported but flagged (marker known to fail there).
(5) Falsification: lift@k ~ 1 cohort-wide => the shortlist is not enriched; reported flat.
    NOTE on (B): candidate lists are HYPOTHESES — never-marked contacts have NO ground truth.
    (A) LOSO precision is the hit-rate for KNOWN hidden SOZ and is an UPPER bound on
    occult-discovery precision (some high-scoring unmarked contacts are genuinely healthy).

Outputs (data/audit/epi_marker_precision/)
    loso_precision_per_patient.csv ; loso_precision_cohort.csv ; candidates_per_patient.csv ; README.md
"""
from __future__ import annotations

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
from audit_85_epi_propagator_recovery import MIN_EPI  # type: ignore
from audit_101_epi_marker_library import build_operators, _resid  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_marker_precision"
ML = ROOT / "data" / "audit" / "epi_marker_library"
BAND = "delta"
MARKER = "heat_t5"          # the verified slow heat kernel
K_LIST = [1, 3, 5, 10]
MIN_OFF = 10                # need a reasonable off-shaft pool to rank
RESP_THR = 0.568           # heat_t5 label-shuffle null p95: responder if validated AUC >= this
TOP_CAND = 12              # candidates listed per patient
DIST_DISTANT = 20.0        # mm: a candidate this far from every seed is a non-trivial (distant) hit


def _auc_lookup():
    """Per-patient validated heat_t5 delta strength-residual LOSO AUC (audit_101)."""
    df = pd.read_csv(ML / "marker_library_per_patient.csv")
    d = df[(df.band == BAND) & (df.marker == MARKER)]
    return dict(zip(d.patient, d.auc_resid))


def process_patient(pat):
    try:
        W = load_phase_fc(pat, "rest_post", BAND)
    except Exception:
        return None, None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None, None
    epi = np.asarray(pm.epi_mask, bool)
    probes = np.asarray(pm.probes, object)
    channels = np.asarray(pm.channels, object)
    epi_idx = np.where(epi)[0]
    if epi_idx.size < 1:
        return None, None
    ops, _intr, strength = build_operators(W)
    O = ops[MARKER]

    # ---- coords + region for interpretability (mm; Desikan-Killiany) ----
    try:
        reg = load_channel_regions(pat)
        coords = reg[["x", "y", "z"]].to_numpy(float) / 1000.0
        region = reg["region"].to_numpy(object) if "region" in reg.columns else np.array(["?"] * N)
    except Exception:
        coords = np.full((N, 3), np.nan)
        region = np.array(["?"] * N)

    # =================== (A) leave-one-shaft-out precision@k ===================
    prec_rows = []
    if epi_idx.size >= max(MIN_EPI, 4) and len(set(probes[epi_idx])) >= 2:
        for p in sorted(set(probes[epi_idx])):
            seeds = epi_idx[probes[epi_idx] != p]
            seed_pr = set(probes[seeds])
            off = np.array([j for j in range(N) if probes[j] not in seed_pr])
            relevant = set(epi_idx[probes[epi_idx] == p].tolist()) & set(off.tolist())
            if off.size < MIN_OFF or len(relevant) < 1:
                continue
            m = O[:, seeds].mean(axis=1)
            mr = _resid(m, strength)
            order_m = off[np.argsort(-mr[off])]          # rank off pool by marker
            order_s = off[np.argsort(-strength[off])]    # rank off pool by node strength
            prev = len(relevant) / off.size
            for k in K_LIST:
                kk = min(k, off.size)
                hit_m = sum(1 for j in order_m[:kk] if j in relevant)
                hit_s = sum(1 for j in order_s[:kk] if j in relevant)
                prec_rows.append(dict(
                    patient=pat, shaft=str(p), k=k, n_off=int(off.size),
                    n_targets=len(relevant), prevalence=prev,
                    prec_marker=hit_m / kk, prec_strength=hit_s / kk,
                    recall_marker=hit_m / len(relevant),
                    lift_marker=(hit_m / kk) / prev if prev > 0 else np.nan,
                    lift_strength=(hit_s / kk) / prev if prev > 0 else np.nan))

    # =================== (B) deployment candidate shortlist ====================
    m = O[:, epi_idx].mean(axis=1)               # seeds = ALL marked SOZ
    mr = _resid(m, strength)
    seedco = coords[epi_idx]
    dist = np.array([
        float(np.nanmin(np.linalg.norm(seedco - coords[j], axis=1)))
        if np.isfinite(coords[j]).all() and np.isfinite(seedco).any() else np.nan
        for j in range(N)])
    unmarked = np.array([j for j in range(N) if not epi[j]])
    order = unmarked[np.argsort(-mr[unmarked])]
    # percentile of each unmarked score within the unmarked pool
    ranks = {j: r for r, j in enumerate(unmarked[np.argsort(mr[unmarked])])}
    cand_rows = []
    for rank, j in enumerate(order[:TOP_CAND], 1):
        cand_rows.append(dict(
            patient=pat, rank=rank, contact=str(channels[j]), region=str(region[j]),
            shaft=str(probes[j]), score_resid=float(mr[j]),
            score_pct=round(100 * ranks[j] / max(1, len(unmarked) - 1), 1),
            dist_to_nearest_seed_mm=round(float(dist[j]), 1) if np.isfinite(dist[j]) else np.nan,
            distant=bool(np.isfinite(dist[j]) and dist[j] >= DIST_DISTANT)))
    return (pd.DataFrame(prec_rows) if prec_rows else None,
            pd.DataFrame(cand_rows) if cand_rows else None)


def cohort_summary(pp, auc, label, patients):
    out = []
    d0 = pp[pp.patient.isin(patients)]
    for k in K_LIST:
        d = d0[d0.k == k]
        if d.empty:
            continue
        out.append(dict(
            group=label, k=k, n_patients=int(d.patient.nunique()),
            med_prevalence=float(d.prevalence.median()),
            med_prec_marker=float(d.prec_marker.median()),
            med_prec_strength=float(d.prec_strength.median()),
            med_recall_marker=float(d.recall_marker.median()),
            med_lift_marker=float(d.lift_marker.median()),
            med_lift_strength=float(d.lift_strength.median())))
    return pd.DataFrame(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    auc = _auc_lookup()
    t0 = time.time()
    prec_frames, cand_frames = [], []
    for pat in es.COHORT:
        pr, ca = process_patient(pat)
        if pr is not None:
            prec_frames.append(pr)
        if ca is not None:
            ca = ca.copy()
            ca.insert(1, "validation_auc", round(float(auc.get(pat, np.nan)), 3))
            ca.insert(2, "responder", bool(auc.get(pat, 0) >= RESP_THR))
            cand_frames.append(ca)

    prec = pd.concat(prec_frames, ignore_index=True)
    # per patient: mean over folds
    pp = prec.groupby(["patient", "k"], as_index=False).agg(
        prec_marker=("prec_marker", "mean"), prec_strength=("prec_strength", "mean"),
        recall_marker=("recall_marker", "mean"), lift_marker=("lift_marker", "mean"),
        lift_strength=("lift_strength", "mean"), prevalence=("prevalence", "mean"))
    pp.to_csv(OUT / "loso_precision_per_patient.csv", index=False)

    responders = [p for p in pp.patient.unique() if auc.get(p, 0) >= RESP_THR]
    coh = pd.concat([
        cohort_summary(pp, auc, "all", list(pp.patient.unique())),
        cohort_summary(pp, auc, "responders", responders),
    ], ignore_index=True)
    coh.to_csv(OUT / "loso_precision_cohort.csv", index=False)

    cand = pd.concat(cand_frames, ignore_index=True)
    cand.to_csv(OUT / "candidates_per_patient.csv", index=False)

    # ----------------------------------------------------------------- console
    print(f"[audit_103] {BAND} {MARKER}: {pp.patient.nunique()} patients with LOSO precision; "
          f"{cand.patient.nunique()} candidate lists. ({time.time()-t0:.1f}s)")
    print("\n(A) HONEST LOSO precision/recall/lift (cohort median):")
    print("    group        k  prevalence  prec(marker)  prec(strength)  recall  lift(marker)  lift(strength)")
    for _, r in coh.iterrows():
        print(f"    {r.group:10s} {int(r.k):2d}   {r.med_prevalence:6.3f}      "
              f"{r.med_prec_marker:6.3f}        {r.med_prec_strength:6.3f}      "
              f"{r.med_recall_marker:5.2f}     {r.med_lift_marker:5.2f}        {r.med_lift_strength:5.2f}")

    print("\n(A) per-patient lift@5 (marker vs strength), validated AUC:")
    p5 = pp[pp.k == 5].sort_values("lift_marker", ascending=False)
    for _, r in p5.iterrows():
        tag = "" if auc.get(r.patient, 0) >= RESP_THR else "  <-- hub-patient"
        print(f"    {r.patient:7s} auc={auc.get(r.patient,float('nan')):.2f}  "
              f"lift@5 marker={r.lift_marker:4.2f} strength={r.lift_strength:4.2f}{tag}")

    print("\n(B) candidate shortlist — top DISTANT (>=20mm) unmarked hypotheses per responder:")
    for pat in [p for p in es.COHORT if auc.get(p, 0) >= RESP_THR]:
        d = cand[(cand.patient == pat) & (cand.distant)].head(2)
        for _, r in d.iterrows():
            print(f"    {pat:7s} #{int(r['rank']):2d} {r.region:24s} shaft {r.shaft:4s} "
                  f"score_pct={r.score_pct:5.1f}  dist={r.dist_to_nearest_seed_mm:4.0f}mm")
    print(f"\n[audit_103] -> {OUT}")


if __name__ == "__main__":
    main()
