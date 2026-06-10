#!/usr/bin/env python3
"""Audit 92 — VII1: within-patient masking recovery of held-out SOZ with the
relational diffusion-affinity marker.

audit_91 showed the matched-strength-z relational features separate epi from non
(AUC ~0.82, strength-orthogonal). This is the rigorous within-patient validation
the directive asks for: HIDE a fraction of each patient's labelled SOZ, pretend
they are unlabelled, learn nothing but "affinity to the VISIBLE SOZ", and try to
recover the hidden SOZ among the healthy contacts. The held-out epi never inform
their own feature (affinity is to the visible set only) — so this is an honest
recovery test, harder than audit_91's leave-self-out screen.

Three scores, head to head:
  - `raw`      affinity to visible SOZ (relational, strength NOT removed),
  - `resid`    that affinity residualised on node strength (linear) — the
               strength-orthogonal relational score (in-loop control; audit_91's
               matched-strength-z is the global certificate),
  - `strength` node strength (the hubness baseline that failed in audit_86/87).

Negatives: all healthy contacts AND a within-shaft control (healthy contacts on
epi-bearing shafts, so the scorer must find the epi CONTACT, not the epi SHAFT).
Two regimes: random 50%-hide (200 splits) and leave-one-epi-out (fine-grained).

Critical preamble
=================
(1) Claim: affinity to the visible SOZ recovers held-out SOZ above chance and
    above node strength, within patient, in the community bands (δ/β/low-γ; α).
(2) Null: chance AUC 0.5 (cohort Wilcoxon of per-patient recovery AUC); strength
    baseline; within-shaft negatives (shaft held constant).
(3) Strongest alternative: hubness (strength) or shaft proximity recovers hidden
    epi and affinity adds nothing.
(4) Reach: `resid` removes the linear strength channel; within-shaft holds shaft
    constant; the head-to-head `resid` vs `strength` is the decisive paired test.
(5) Falsification: if `resid` recovery ≈ 0.5 or does not beat `strength`, the
    marker is hubness/shaft and there is no deployable per-node recovery.

Outputs (``data/audit/epi_marker_relational/``)
    masked_recovery_per_patient.csv   per (patient, band): AUC/PR-AUC per score+regime
    masked_recovery_cohort.csv        per (band)
    README_masked_recovery.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.metrics import average_precision_score, roc_auc_score

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.patient import build_epi_masks

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI, N_TAU  # type: ignore
from audit_91_epi_relational_features import relational_feats  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_marker_relational"
OUT.mkdir(parents=True, exist_ok=True)

N_SPLIT = 200
HIDE_FRAC = 0.5
SCORES = ("raw", "resid", "strength")


def _resid_on_strength(f, s):
    """Linear-residualise feature f on node strength s (strength-orthogonal)."""
    f = np.asarray(f, float)
    s = np.asarray(s, float)
    m = np.isfinite(f) & np.isfinite(s)
    if m.sum() < 3 or np.nanstd(s[m]) < 1e-12:
        return f - np.nanmean(f[m]) if m.any() else f
    b = np.polyfit(s[m], f[m], 1)
    return f - np.polyval(b, s)


def _auc(score, pos_idx, neg_idx):
    s = np.r_[score[pos_idx], score[neg_idx]]
    y = np.r_[np.ones(pos_idx.size), np.zeros(neg_idx.size)]
    ok = np.isfinite(s)
    if np.unique(s[ok]).size < 2 or y[ok].sum() < 1 or (y[ok] == 0).sum() < 1:
        return np.nan
    return float(roc_auc_score(y[ok], s[ok]))


def _ap(score, pos_idx, neg_idx):
    s = np.r_[score[pos_idx], score[neg_idx]]
    y = np.r_[np.ones(pos_idx.size), np.zeros(neg_idx.size)]
    ok = np.isfinite(s)
    if y[ok].sum() < 1 or (y[ok] == 0).sum() < 1:
        return np.nan
    return float(average_precision_score(y[ok], s[ok]))


def _aff_to_set(lam_o, U_o, taus, set_mask, strength):
    """τ-mean affinity-to-set for all nodes (raw) + its strength residual."""
    e = np.zeros(U_o.shape[0], bool)
    e[set_mask] = True
    fr = relational_feats(lam_o, U_o, e, taus)
    f_aff = np.nanmean([fr[f"f_aff_t{ti}"] for ti in range(N_TAU)], axis=0)
    return f_aff, _resid_on_strength(f_aff, strength)


def per_patient_band(pat, band, rng, verbose):
    try:
        W = load_phase_fc(pat, "rest_post", band)
    except Exception as e:
        if verbose:
            print(f"[audit_92] SKIP {pat}/{band}: {e}")
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    probe = np.asarray(pm.probes, object)
    if epi.sum() < MIN_EPI or (~epi).sum() < 3:
        return None
    lam_o, U_o = np.linalg.eigh(np.diag(W.sum(1)) - W)
    if float(lam_o[-1]) <= 0:
        return None
    taus = np.geomspace(1.0 / lam_o[-1], 10.0 / lam_o[-1], N_TAU)
    strength = W.sum(1)

    epi_idx = np.where(epi)[0]
    healthy = np.where(~epi)[0]
    epi_shafts = set(probe[epi_idx])
    within = np.array([j for j in healthy if probe[j] in epi_shafts], dtype=int)
    n_hide = max(1, int(round(epi_idx.size * HIDE_FRAC)))
    if epi_idx.size - n_hide < 1:
        return None

    # --- random 50%-hide recovery ---
    acc = {f"{sc}_{neg}": [] for sc in SCORES for neg in ("all", "wshaft")}
    ap_acc = {f"{sc}_all": [] for sc in ("raw", "resid")}
    for _ in range(N_SPLIT):
        perm = rng.permutation(epi_idx.size)
        hidden = epi_idx[perm[:n_hide]]
        visible = epi_idx[perm[n_hide:]]
        f_raw, f_res = _aff_to_set(lam_o, U_o, taus, visible, strength)
        scval = {"raw": f_raw, "resid": f_res, "strength": strength}
        for sc in SCORES:
            acc[f"{sc}_all"].append(_auc(scval[sc], hidden, healthy))
            if within.size >= 2:
                acc[f"{sc}_wshaft"].append(_auc(scval[sc], hidden, within))
        for sc in ("raw", "resid"):
            ap_acc[f"{sc}_all"].append(_ap(scval[sc], hidden, healthy))

    # --- leave-one-epi-out (fine-grained) ---
    loeo = {sc: [] for sc in SCORES}
    for h in range(epi_idx.size):
        hidden = epi_idx[[h]]
        visible = np.delete(epi_idx, h)
        f_raw, f_res = _aff_to_set(lam_o, U_o, taus, visible, strength)
        scval = {"raw": f_raw, "resid": f_res, "strength": strength}
        for sc in SCORES:
            loeo[sc].append(_auc(scval[sc], hidden, healthy))

    row = {"patient": pat, "band": band, "n_epi": int(epi_idx.size),
           "n_within_shaft": int(within.size)}
    for k, v in acc.items():
        v = [x for x in v if np.isfinite(x)]
        row[f"auc_{k}"] = float(np.mean(v)) if v else np.nan
    for k, v in ap_acc.items():
        v = [x for x in v if np.isfinite(x)]
        row[f"ap_{k}"] = float(np.mean(v)) if v else np.nan
    for sc in SCORES:
        v = [x for x in loeo[sc] if np.isfinite(x)]
        row[f"loeo_{sc}"] = float(np.mean(v)) if v else np.nan
    row["prevalence"] = float(epi_idx.size / N)
    return row


def cohort(df):
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(df.band)]:
        d = df[df.band == band]
        row = {"band": band, "n_patients": len(d)}
        for col in [c for c in d.columns if c.startswith(("auc_", "ap_", "loeo_"))]:
            v = d[col].to_numpy()
            v = v[np.isfinite(v)]
            row[f"med_{col}"] = float(np.median(v)) if v.size else np.nan
            # Wilcoxon vs chance for AUC/LOEO columns
            if col.startswith(("auc_", "loeo_")) and v.size >= 3:
                try:
                    _, p = wilcoxon(v - 0.5, alternative="greater")
                except Exception:
                    p = np.nan
                row[f"p_{col}"] = p
        # decisive paired: resid beats strength (all-healthy)
        pr = d.dropna(subset=["auc_resid_all", "auc_strength_all"])
        if len(pr) >= 3:
            try:
                _, p = wilcoxon(pr["auc_resid_all"] - pr["auc_strength_all"],
                                alternative="greater")
            except Exception:
                p = np.nan
            row["resid_minus_strength_med"] = float(
                (pr["auc_resid_all"] - pr["auc_strength_all"]).median())
            row["resid_beats_strength_p"] = p
        out.append(row)
    return pd.DataFrame(out)


def write_readme(coh, runtime):
    L = ["---", "name: epi_masked_recovery_relational",
         "scope: VII1_within_patient_masking_recovery_relational_affinity_marker",
         "era: COHORT_N10 / IMCOH_ABS",
         f"date: {time.strftime('%Y-%m-%d')}",
         "build_script: scripts/01_compute/audit/audit_92_epi_masked_recovery_relational.py",
         "---", "",
         "# VII1 — masking recovery of held-out SOZ (relational affinity marker)",
         "",
         "**Head.** Hide 50% of each patient's labelled SOZ, recover them among "
         "healthy contacts using affinity to the VISIBLE SOZ (200 splits) + "
         "leave-one-epi-out. `resid` = affinity residualised on strength "
         "(strength-orthogonal); `strength` = hubness baseline. `wshaft` = "
         "negatives restricted to healthy contacts on epi-bearing shafts.", "",
         "## Recovery AUC (median over patients)", "",
         "| band | resid (all) | p | strength (all) | resid−strength (p) "
         "| resid (wshaft) | LOEO resid | PR-AUC resid (prev) |",
         "|---|---|---|---|---|---|---|---|"]
    for band in [b for b in es.ALL_BANDS if b in set(coh.band)]:
        r = coh[coh.band == band]
        if r.empty:
            continue
        r = r.iloc[0]
        prev = r.get("med_ap_resid_all", np.nan)
        L.append(
            f"| {BRAIN_BAND_TEX_DICT.get(band, band)} "
            f"| {r.get('med_auc_resid_all', np.nan):.3f} "
            f"| {r.get('p_auc_resid_all', np.nan):.4f} "
            f"| {r.get('med_auc_strength_all', np.nan):.3f} "
            f"| {r.get('resid_minus_strength_med', np.nan):+.3f} "
            f"({r.get('resid_beats_strength_p', np.nan):.3f}) "
            f"| {r.get('med_auc_resid_wshaft', np.nan):.3f} "
            f"| {r.get('med_loeo_resid', np.nan):.3f} "
            f"| {prev:.3f} |")
    L += ["", "## Reading",
          "- `resid (all) > 0.5` (small p) ⇒ strength-orthogonal affinity recovers "
          "hidden SOZ. `resid−strength > 0` ⇒ it beats the hubness baseline. "
          "`resid (wshaft)` surviving ⇒ not shaft proximity.",
          "- Numbers per patient in `masked_recovery_per_patient.csv`; per band in "
          "`masked_recovery_cohort.csv`.",
          f"- splits {N_SPLIT}, hide {HIDE_FRAC}, min epi {MIN_EPI}; "
          f"wall-clock {runtime:.1f}s"]
    (OUT / "README_masked_recovery.md").write_text("\n".join(L) + "\n",
                                                   encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    rows = []
    for band in args.bands:
        for pat in args.patients:
            rng = np.random.default_rng(20260608 + sum(ord(c) for c in pat + band))
            tc = time.time()
            r = per_patient_band(pat, band, rng, args.verbose)
            if r is not None:
                rows.append(r)
                print(f"[audit_92] {pat}/{band}: n_epi={r['n_epi']} "
                      f"resid={r['auc_resid_all']:.3f} str={r['auc_strength_all']:.3f} "
                      f"({time.time()-tc:.1f}s)")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "masked_recovery_per_patient.csv", index=False)
    coh = cohort(df)
    coh.to_csv(OUT / "masked_recovery_cohort.csv", index=False)
    runtime = time.time() - t0
    write_readme(coh, runtime)

    print("\n[audit_92] recovery AUC (median): resid | strength | resid-str p")
    for band in [b for b in es.ALL_BANDS if b in set(coh.band)]:
        r = coh[coh.band == band].iloc[0]
        print(f"  {band:10s}: resid={r.get('med_auc_resid_all', np.nan):.3f} "
              f"(p{r.get('p_auc_resid_all', np.nan):.3f}) "
              f"str={r.get('med_auc_strength_all', np.nan):.3f} "
              f"Δ={r.get('resid_minus_strength_med', np.nan):+.3f} "
              f"(p{r.get('resid_beats_strength_p', np.nan):.3f})")
    print(f"[audit_92] done in {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
