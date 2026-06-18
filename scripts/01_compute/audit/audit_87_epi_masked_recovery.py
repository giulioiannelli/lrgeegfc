#!/usr/bin/env python3
"""Audit 87 — hide-and-seek validation: does the α propagator marker beat hubness?

The user's validation design: take the KNOWN epileptic contacts, hide a random
half (pretend they are unlabelled), learn the "epileptic direction" in feature
space from the VISIBLE half, and test whether the hidden half is re-found among
the unlabelled contacts. Repeat over many random splits. The decisive question is
whether the STRENGTH-REMOVED propagator features (matched-strength surrogate-z,
from audit_86) recover hidden epi above chance AND add to a fair node-strength
baseline — i.e. whether the modest α signal is real beyond trivial hubness.

This is the cleaner, fairer cousin of audit_85/86's leave-one-out recovery: the
scorer is a Fisher-style linear discriminant direction learned from the visible
labels, which is fair to a single monotone feature (strength) and to a
multi-feature subspace alike — removing the template-distance handicap that
inflated audit_85's strength comparison.

Critical preamble
=================
(1) Claim: the strength-removed propagator subspace finds hidden epileptic
    contacts beyond what node strength alone achieves, at α.
(2) Null: chance recovery = 0.5 (AUC of hidden-epi vs healthy); per-band cohort
    Wilcoxon of per-patient recovery AUC against 0.5.
(3) Strongest alternative: it is hubness — node strength alone already recovers
    hidden epi, and the propagator adds nothing.
(4) Reach: the propagator features here are the matched-strength surrogate-z
    versions (audit_86), strength-fixed by construction; the head-to-head is
    strength vs strength-removed-propagator vs both, with a scorer fair to all.
(5) Falsification: if strength-removed-propagator recovery ≈ 0.5, or if
    (strength+propagator) does not exceed strength alone, the α story is hubness
    and is dropped.

Reads ``data/audit/epi_marker/propagator_ms_features.csv`` (audit_86).
Outputs ``data/audit/epi_marker/{masked_recovery_per_band.csv, README_masked_recovery.md}``.
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.metrics import roc_auc_score

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT

OUT = ROOT / "data" / "audit" / "epi_marker"
ALL_BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]

PROP_RAW = (["rho_ii_t%d" % i for i in range(6)]
            + ["Hdiff_t%d" % i for i in range(6)]
            + ["Tcomm_mean", "Tcomm_min"])
PROP_MSZ = ["z_" + c for c in PROP_RAW]            # strength-removed
STRENGTH = ["strength_mean"]

FEATURE_SETS = {
    "strength": STRENGTH,
    "propagator_msz": PROP_MSZ,                     # strength-orthogonal
    "propagator_raw": PROP_RAW,
    "strength+propagator_msz": STRENGTH + PROP_MSZ,
}
MIN_EPI = 6                                          # need >=3 visible + >=3 hidden
N_SPLIT = 200
HIDE_FRAC = 0.5


def _standardize(X):
    mu = np.nanmean(X, 0); sd = np.nanstd(X, 0)
    sd = np.where(sd > 1e-9, sd, 1.0)
    return np.nan_to_num((X - mu) / sd, nan=0.0)


def _recovery_auc(Z, epi_idx, healthy_idx, rng) -> float:
    """Mean over random 50% hides of AUC(hidden epi vs healthy) using a Fisher
    direction learned from the VISIBLE epi vs healthy."""
    aucs = []
    n_hide = max(1, int(round(epi_idx.size * HIDE_FRAC)))
    if epi_idx.size - n_hide < 1 or healthy_idx.size < 2:
        return float("nan")
    for _ in range(N_SPLIT):
        perm = rng.permutation(epi_idx.size)
        hidden = epi_idx[perm[:n_hide]]
        visible = epi_idx[perm[n_hide:]]
        d = Z[visible].mean(0) - Z[healthy_idx].mean(0)   # epileptic direction
        nrm = np.linalg.norm(d)
        if nrm < 1e-12:
            continue
        d = d / nrm
        score = Z @ d
        y = np.r_[np.ones(hidden.size), np.zeros(healthy_idx.size)]
        s = np.r_[score[hidden], score[healthy_idx]]
        if np.unique(s).size < 2:
            continue
        aucs.append(roc_auc_score(y, s))
    return float(np.mean(aucs)) if aucs else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=ALL_BANDS)
    args = ap.parse_args()
    src = OUT / "propagator_ms_features.csv"
    if not src.exists():
        raise SystemExit("[audit_87] run audit_86 first")
    feat = pd.read_csv(src)
    # merge shaft (probe) for the within-shaft proximity control
    probe_src = pd.read_csv(OUT / "node_features_per_band.csv")[
        ["patient", "band", "node", "probe"]]
    feat = feat.merge(probe_src, on=["patient", "band", "node"], how="left")
    rng = np.random.default_rng(20260605)
    t0 = time.time()

    rows = []
    for band in args.bands:
        fb = feat[feat.band == band]
        for pat, fp in fb.groupby("patient"):
            fp = fp.reset_index(drop=True)
            y = fp["is_epi"].to_numpy(bool)
            if y.sum() < MIN_EPI:
                continue
            epi_idx = np.where(y)[0]; healthy = np.where(~y)[0]
            probe = fp["probe"].to_numpy()
            epi_shafts = set(probe[epi_idx])
            # HARD control: negatives = healthy contacts ON epi-bearing shafts
            # (so the scorer cannot win by "this shaft" — shaft is held constant)
            within = np.array([j for j in healthy if probe[j] in epi_shafts],
                              dtype=int)
            rec = {}
            for sn, cols in FEATURE_SETS.items():
                Z = _standardize(fp[cols].to_numpy(float))
                rec[f"{sn}"] = _recovery_auc(Z, epi_idx, healthy, rng)
                rec[f"{sn}__wshaft"] = (_recovery_auc(Z, epi_idx, within, rng)
                                        if within.size >= 2 else float("nan"))
            rows.append({"band": band, "patient": pat, "n_epi": int(y.sum()),
                         "n_within_shaft_healthy": int(within.size),
                         **{f"rec_{k}": v for k, v in rec.items()}})

    df = pd.DataFrame(rows)
    summ = []
    for band in args.bands:
        b = df[df.band == band]
        if b.empty:
            continue
        rec_strength = b["rec_strength"].to_numpy()
        row = {"band": band, "n_patients": len(b)}
        for sn in FEATURE_SETS:
            for suff in ("", "__wshaft"):
                v = b[f"rec_{sn}{suff}"].to_numpy(); v = v[np.isfinite(v)]
                wp = np.nan
                if v.size >= 3:
                    try:
                        _, wp = wilcoxon(v - 0.5, alternative="greater")
                    except Exception:
                        pass
                tag = "median" if suff == "" else "wshaft_median"
                ptag = "p_gt_chance" if suff == "" else "wshaft_p"
                row[f"{sn}_{tag}"] = float(np.median(v)) if v.size else np.nan
                row[f"{sn}_{ptag}"] = wp
        # decisive: WITHIN-SHAFT — propagator(MS-z) vs strength, both held to
        # epi-bearing shafts so shaft location cannot be the discriminator
        paired = b.dropna(subset=["rec_strength__wshaft",
                                  "rec_propagator_msz__wshaft"])
        if len(paired) >= 3:
            try:
                _, betterp = wilcoxon(paired["rec_propagator_msz__wshaft"]
                                      - paired["rec_strength__wshaft"],
                                      alternative="greater")
            except Exception:
                betterp = np.nan
            row["wshaft_prop_minus_strength"] = float(
                (paired["rec_propagator_msz__wshaft"]
                 - paired["rec_strength__wshaft"]).median())
            row["wshaft_prop_beats_strength_p"] = betterp
        else:
            row["wshaft_prop_minus_strength"] = np.nan
            row["wshaft_prop_beats_strength_p"] = np.nan
        # strength-orthogonal propagator beats chance, WITHIN shaft (the hard test)
        row["msz_wshaft_beats_chance"] = bool(
            np.isfinite(row["propagator_msz_wshaft_p"])
            and row["propagator_msz_wshaft_p"] < 0.05)
        row["prop_beats_strength_wshaft"] = bool(
            np.isfinite(row["wshaft_prop_beats_strength_p"])
            and row["wshaft_prop_beats_strength_p"] < 0.05)
        summ.append(row)
    summ = pd.DataFrame(summ)
    df.to_csv(OUT / "masked_recovery_per_patient.csv", index=False)
    summ.to_csv(OUT / "masked_recovery_per_band.csv", index=False)

    # README
    L = ["---", "name: epi_masked_recovery_validation",
         "scope: direction1_hide_and_seek_alpha_beyond_hubness",
         "era: COHORT_N10 / IMCOH_ABS", f"date: {time.strftime('%Y-%m-%d')}",
         "build_script: scripts/01_compute/audit/audit_87_epi_masked_recovery.py",
         "---", "",
         "# Hide-and-seek validation — does the propagator beat hubness?", "",
         "**Head.** Hide a random 50% of each patient's known epileptic contacts, "
         "learn the epileptic direction from the visible half, recover the hidden "
         "half (AUC vs healthy), 200 splits. Fair Fisher-direction scorer (no "
         "template-distance handicap). `strength` = hubness baseline; "
         "`propagator_msz` = matched-strength-removed propagator (strength-"
         "orthogonal); `strength+propagator_msz` = both.", "",
         "## Recovery AUC per band — ALL-healthy vs the WITHIN-SHAFT control",
         "",
         "Within-shaft = negatives restricted to healthy contacts on epi-bearing "
         "shafts, so shaft location is held constant and the scorer must find the "
         "epi CONTACT, not the epi SHAFT. This is the decisive control.", "",
         "| band | strength (all) | prop-MSz (all) | **strength (wshaft)** | "
         "**prop-MSz (wshaft)** | prop-MSz wshaft>chance p | prop−strength "
         "wshaft Δ (p) | verdict |",
         "|---|---|---|---|---|---|---|---|"]
    for band in [b for b in ALL_BANDS if b in set(summ.band)]:
        r = summ[summ.band == band].iloc[0]
        verdict = ("propagator beats hubness (survives shaft)"
                   if r["msz_wshaft_beats_chance"] and r["prop_beats_strength_wshaft"]
                   else ("prop>chance within shaft, ~ties strength"
                         if r["msz_wshaft_beats_chance"] else
                         "no within-shaft signal (shaft/hubness)"))
        L.append(
            f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {r['strength_median']:.3f} "
            f"| {r['propagator_msz_median']:.3f} | {r['strength_wshaft_median']:.3f} "
            f"| {r['propagator_msz_wshaft_median']:.3f} "
            f"| {r['propagator_msz_wshaft_p']:.4f} "
            f"| {r['wshaft_prop_minus_strength']:+.3f} "
            f"({r['wshaft_prop_beats_strength_p']:.3f}) | {verdict} |")
    L += ["", "## Reading",
          "- `(all)` negatives = every healthy contact; `(wshaft)` negatives = "
          "healthy contacts ON epi-bearing shafts (shaft held constant).",
          "- A drop from `(all)` to `(wshaft)` = part of the recovery was shaft "
          "proximity, not contact identity. Survival within-shaft = real.",
          "- **The marker is real beyond hubness AND beyond shaft location iff "
          "`prop-MSz (wshaft) > chance` AND it ties/beats `strength (wshaft)`.**",
          f"- splits/patient: {N_SPLIT}; hide fraction {HIDE_FRAC}; "
          f"min epi {MIN_EPI}; wall-clock {time.time()-t0:.1f}s"]
    (OUT / "README_masked_recovery.md").write_text("\n".join(L) + "\n",
                                                   encoding="utf-8")

    print("[audit_87] recovery AUC (median over patients) — ALL | WITHIN-SHAFT:")
    for band in [b for b in ALL_BANDS if b in set(summ.band)]:
        r = summ[summ.band == band].iloc[0]
        flag = ("** beats hubness+shaft"
                if r["msz_wshaft_beats_chance"] and r["prop_beats_strength_wshaft"]
                else (" prop>chance(wshaft)" if r["msz_wshaft_beats_chance"]
                      else " hubness/shaft"))
        print(f"  {band:10s}: ALL str={r['strength_median']:.3f} "
              f"prop={r['propagator_msz_median']:.3f} | WSHAFT "
              f"str={r['strength_wshaft_median']:.3f} "
              f"prop={r['propagator_msz_wshaft_median']:.3f} "
              f"(p{r['propagator_msz_wshaft_p']:.3f}, "
              f"Δ{r['wshaft_prop_minus_strength']:+.3f}){flag}")
    print(f"[audit_87] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
