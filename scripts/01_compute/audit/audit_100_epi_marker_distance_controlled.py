#!/usr/bin/env python3
"""Audit 100 — THE decisive test: does propagator affinity-to-seeds mark SOZ AFTER you
control away the two trivial confounds — distance-to-nearest-seed (focal anatomy) and
node strength (hubness)?

Motivation (honest). Earlier marker scores (audit_92/98: AUC ~0.8, "beats strength")
never controlled for spatial proximity. audit_99 showed a trivial "flag the nearest
contacts" baseline beats the marker ~2x overall. So every prior "good" number is
confounded by the fact that the remaining SOZ sit next to the marked ones. This script
asks the only question that survives that: among contacts AT THE SAME DISTANCE from the
seeds AND THE SAME STRENGTH, do the true SOZ still have higher propagator affinity?

Method — conditional (matched-pair) AUC, leave-one-SOZ-out.
    For each patient/band, hide each SOZ in turn; seeds = all other SOZ. For the hidden
    SOZ ("case") and every healthy contact, compute RAW affinity-to-seeds, distance to
    nearest seed (mm), and node strength. Accumulate all (case, healthy) pairs across
    folds, then compute concordance (marker_case > marker_healthy):
      - A_raw      : all pairs                       (the inflated, uncontrolled number)
      - A_dist     : pairs matched on distance       (|Δdist| <= DELTA_MM)  <- discounts anatomy
      - A_diststr  : matched on distance AND strength (also |Δstrength| <= S_TOL)
    A_dist (and A_diststr) > 0.5 with a cohort sign test = the propagator carries SOZ
    information BEYOND proximity (and hubness). A_dist ~ 0.5 = it does not; the marker was
    proximity all along. Distance's own conditional AUC is reported as a sanity check
    (must be ~0.5 — confirms the matching neutralises distance).

Critical preamble
=================
(1) Claim: propagator affinity-to-seeds discriminates SOZ from healthy at MATCHED
    distance-to-seed and matched strength.
(2) Null: A_dist = 0.5 (no signal once proximity is matched).
(3) Strongest alternative: the whole effect is focal proximity (+ hubness) — exactly
    what the matching removes.
(4) Reach: matched-pair conditional AUC is the textbook way to partial out a confounder
    without a parametric residual; DELTA_MM swept (5/10/15) for robustness; raw affinity
    (not strength-residual) is used so strength is controlled explicitly by matching, not
    pre-removed.
(5) Falsification: A_dist <= ~0.5 cohort-wide (sign test n.s.) => the propagator does NOT
    mark SOZ beyond trivial anatomy. Reported as such, no spin.

Outputs (data/audit/epi_marker_distance_controlled/)
    distance_controlled_per_patient.csv ; distance_controlled_cohort.csv ; README.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.utils.io.regions import load_channel_regions

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI, N_TAU  # type: ignore
from audit_91_epi_relational_features import relational_feats  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_marker_distance_controlled"
DELTA_MM = 10.0           # distance-match tolerance (mm)
S_TOL_FRAC = 0.5          # strength-match tolerance, as a fraction of strength std


def _raw_affinity(lam, U, taus, seeds, N):
    """τ-mean RAW affinity-to-seeds for every node (no strength residual)."""
    e = np.zeros(N, bool)
    e[seeds] = True
    fr = relational_feats(lam, U, e, taus)
    return np.nanmean([fr[f"f_aff_t{ti}"] for ti in range(N_TAU)], axis=0)


def per_patient_band(pat, band):
    try:
        W = load_phase_fc(pat, "rest_post", band)
    except Exception:
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    epi_idx = np.where(epi)[0]
    healthy = np.where(~epi)[0]
    if epi_idx.size < max(MIN_EPI, 4) or healthy.size < 5:
        return None
    lam, U = np.linalg.eigh(np.diag(W.sum(1)) - W)
    if float(lam[-1]) <= 0:
        return None
    taus = np.geomspace(1.0 / lam[-1], 10.0 / lam[-1], N_TAU)
    strength = W.sum(1)
    s_tol = S_TOL_FRAC * float(np.std(strength))
    coords = load_channel_regions(pat)[["x", "y", "z"]].to_numpy(float) / 1000.0  # mm

    # accumulate concordance counts over all leave-one-out (case, healthy) pairs
    cnt = {kk: [0.0, 0.0] for kk in ("raw", "dist", "diststr", "dist_ofdist")}  # [conc,total]
    raw_dist_closer = [0.0, 0.0]  # fraction(case closer to seed than healthy) — how focal
    for i in epi_idx:
        seeds = epi_idx[epi_idx != i]
        f = _raw_affinity(lam, U, taus, seeds, N)
        d = np.full(N, np.inf)
        cs = coords[seeds]
        for j in np.r_[i, healthy]:
            if np.isfinite(coords[j]).all():
                d[j] = float(np.nanmin(np.linalg.norm(cs - coords[j], axis=1)))
        mi, di, si = f[i], d[i], strength[i]
        if not np.isfinite(mi) or not np.isfinite(di):
            continue
        for h in healthy:
            mh, dh, sh = f[h], d[h], strength[h]
            if not np.isfinite(mh) or not np.isfinite(dh):
                continue
            conc = 1.0 if mi > mh else (0.5 if mi == mh else 0.0)
            cnt["raw"][0] += conc; cnt["raw"][1] += 1
            raw_dist_closer[0] += 1.0 if di < dh else (0.5 if di == dh else 0.0)
            raw_dist_closer[1] += 1
            if abs(di - dh) <= DELTA_MM:
                cnt["dist"][0] += conc; cnt["dist"][1] += 1
                # sanity: distance's own conditional AUC (case lower-dist?) -> must be ~0.5
                cnt["dist_ofdist"][0] += 1.0 if di < dh else (0.5 if di == dh else 0.0)
                cnt["dist_ofdist"][1] += 1
                if abs(si - sh) <= s_tol:
                    cnt["diststr"][0] += conc; cnt["diststr"][1] += 1

    def auc(k):
        c, t = cnt[k]
        return (c / t) if t >= 20 else np.nan
    return {
        "patient": pat, "band": band, "n_epi": int(epi_idx.size),
        "auc_raw": auc("raw"),
        "auc_dist_controlled": auc("dist"),
        "auc_dist_strength_controlled": auc("diststr"),
        "sanity_distance_conditional": auc("dist_ofdist"),   # must be ~0.5
        "frac_case_closer": (raw_dist_closer[0] / raw_dist_closer[1]
                             if raw_dist_closer[1] else np.nan),
        "n_pairs_dist_matched": int(cnt["dist"][1]),
    }


def cohort(df):
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(df.band)]:
        d = df[df.band == band]
        row = {"band": band, "n_patients": len(d)}
        for col in ("auc_raw", "auc_dist_controlled", "auc_dist_strength_controlled",
                    "sanity_distance_conditional", "frac_case_closer"):
            v = d[col].to_numpy(float)
            v = v[np.isfinite(v)]
            row[f"med_{col}"] = float(np.median(v)) if v.size else np.nan
            if col in ("auc_dist_controlled", "auc_dist_strength_controlled") and v.size >= 5:
                row[f"n_above_half_{col}"] = int((v > 0.5).sum())
                try:
                    row[f"p_{col}"] = float(wilcoxon(v - 0.5, alternative="greater")[1])
                except ValueError:
                    row[f"p_{col}"] = np.nan
        out.append(row)
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=["delta", "beta", "low_gamma", "alpha"])
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rows = []
    for band in args.bands:
        for pat in args.patients:
            r = per_patient_band(pat, band)
            if r is not None:
                rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "distance_controlled_per_patient.csv", index=False)
    coh = cohort(df)
    coh.to_csv(OUT / "distance_controlled_cohort.csv", index=False)
    rt = time.time() - t0
    print("[audit_100] conditional AUC (median over patients): "
          "raw -> distance-matched -> distance+strength-matched | sanity(dist)=~0.5")
    for _, r in coh.iterrows():
        nd = r.get("n_above_half_auc_dist_controlled", np.nan)
        pd_ = r.get("p_auc_dist_controlled", np.nan)
        print(f"  {r.band:10s}: raw {r.med_auc_raw:.2f} -> dist {r.med_auc_dist_controlled:.2f} "
              f"-> dist+str {r.med_auc_dist_strength_controlled:.2f} "
              f"| dist>0.5 in {nd}/{int(r.n_patients)} (p={pd_:.3f}) "
              f"| sanity {r.med_sanity_distance_conditional:.2f} "
              f"| case-closer {r.med_frac_case_closer:.2f}")
    print(f"[audit_100] done in {rt:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
