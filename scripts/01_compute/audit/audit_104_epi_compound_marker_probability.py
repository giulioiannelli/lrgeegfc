#!/usr/bin/env python3
"""Audit 104 — compound (propagator + node-strength) SOZ marker, and a leave-one-patient-out
CALIBRATED probability P(node = SOZ). Does combining the two complementary signals beat the
propagator alone, and can we emit an honest per-node probability?

Motivation. audit_103 showed the two are complementary: the δ heat-kernel propagator carries
the RESPONDERS, node strength carries the HUB-patients (Pat_10/15). A compound should cover
both. But the two populations have OPPOSITE useful features (affinity for responders, strength
for hubs, where affinity even ANTI-predicts) — so a single fixed linear model may not transfer
(this is exactly why the earlier per-node classifier audit_80 failed cross-patient). We test
this head-on, OFF-SHAFT, and report whatever happens.

Compounds (all per-node, scored on the SAME leave-one-shaft-out off-shaft distant-discovery
protocol as audit_101/103 — proximity discarded, never an all-contacts metric):
  1. affinity        — δ heat_t5 affinity-to-seeds, strength-residual   (the verified marker)
  2. strength        — node strength                                     (hub baseline)
  3. union_max       — max(pct(affinity), pct(strength))                 (cover both pops)
  4. mean            — 0.5 pct(affinity) + 0.5 pct(strength)
  5. adaptive        — (1-r) pct(affinity) + r pct(strength), r = seed strength percentile
                       (seed-based regime detector: hub-like seeds -> weight strength)
  6. logistic (LOPO) — P(SOZ)=sigma(b0+b1 aff_z+b2 str_z), trained on 9 patients, scored on the
                       10th; yields the calibrated probability.

Critical preamble
=================
(1) Claim: a propagator+strength compound raises off-shaft distant-discovery above the
    propagator alone (more patients > 0.5 / higher cohort lift), and a LOPO logistic emits a
    calibrated P(SOZ).
(2) Null: compound AUC = affinity AUC (no gain); logistic held-out AUC = 0.5 (no transfer);
    calibration no better than predicting prevalence (Brier >= prevalence Brier).
(3) Strongest alternative: the "gain" is just strength leaking proximity/hubness back in ->
    controlled because strength is scored on the SAME off-shaft pool (distant) and reported
    explicitly as its own baseline; affinity is strength-residual.
(4) Reach: every compound uses identical LOSO folds; logistic is LEAVE-ONE-PATIENT-OUT with
    within-fold z-scoring (so per-patient scale can't leak); calibration assessed on held-out
    predictions only (reliability deciles + Brier).
(5) Falsification: if no compound beats affinity AND the logistic does not transfer (held-out
    AUC ~0.5) or is mis-calibrated, we report that the probability is NOT reliable cross-patient
    and the honest deliverable stays the rank-based shortlist. No spin.

Outputs (data/audit/epi_marker_compound/)
    compound_auc_{per_patient,cohort}.csv ; logistic_lopo_per_patient.csv ;
    calibration_curve.csv ; candidates_with_prob.csv ; README.md
"""
from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.utils.io.regions import load_channel_regions

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI  # type: ignore
from audit_101_epi_marker_library import build_operators, _resid  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_marker_compound"
ML = ROOT / "data" / "audit" / "epi_marker_library"
BAND = "delta"
MARKER = "heat_t5"
K_LIST = [5, 10]
MIN_OFF = 10
RESP_THR = 0.568
TOP_CAND = 12
R_THR = 0.70            # seed strength percentile above which a patient is treated as hub-regime
COMPOUNDS = ["affinity", "strength", "union_max", "mean", "adaptive", "switch"]


def _pct(x):
    return (rankdata(x) - 0.5) / len(x)


def _conc(case, ctrl):
    c = 0.0
    for v in case:
        c += np.sum(v > ctrl) + 0.5 * np.sum(v == ctrl)
    return c, case.size * ctrl.size


def fit_logistic(X, y, iters=60, l2=1e-2):
    """IRLS logistic with L2; X already has an intercept column. Returns weights."""
    w = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-(X @ w)))
        Wd = p * (1 - p) + 1e-9
        grad = X.T @ (p - y) + l2 * w
        H = (X * Wd[:, None]).T @ X + l2 * np.eye(X.shape[1])
        try:
            w = w - np.linalg.solve(H, grad)
        except np.linalg.LinAlgError:
            break
    return w


def patient_pack(pat):
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
    if epi_idx.size < 1:
        return None
    ops, _i, strength = build_operators(W)
    return dict(pat=pat, N=N, epi=epi, probes=probes, epi_idx=epi_idx,
                strength=strength, O=ops[MARKER], channels=np.asarray(pm.channels, object))


def fold_rows(d):
    """Yield (seeds, off, target_mask_off, aff_off, str_off, r_seed) per held-out SOZ shaft."""
    epi_idx, probes, strength, O, N = (d["epi_idx"], d["probes"], d["strength"], d["O"], d["N"])
    str_pct_all = _pct(strength)
    for p in sorted(set(probes[epi_idx])):
        seeds = epi_idx[probes[epi_idx] != p]
        seed_pr = set(probes[seeds])
        off = np.array([j for j in range(N) if probes[j] not in seed_pr])
        if off.size < MIN_OFF or seeds.size < 1:
            continue
        relevant = set(epi_idx[probes[epi_idx] == p].tolist()) & set(off.tolist())
        if len(relevant) < 1:
            continue
        aff = _resid(O[:, seeds].mean(axis=1), strength)
        tmask = np.array([j in relevant for j in off])
        r_seed = float(np.mean(str_pct_all[seeds]))      # regime: are seeds hubs?
        yield seeds, off, tmask, aff[off], strength[off], r_seed


def compound_scores(aff_off, str_off, r_seed):
    pa, psr = _pct(aff_off), _pct(str_off)
    return {
        "affinity": aff_off,
        "strength": str_off,
        "union_max": np.maximum(pa, psr),
        "mean": 0.5 * pa + 0.5 * psr,
        "adaptive": (1 - r_seed) * pa + r_seed * psr,
        "switch": (str_off if r_seed >= R_THR else aff_off),   # hard regime gate (no dilution)
    }


def eval_compounds(packs):
    rows = []
    for d in packs.values():
        if d["epi_idx"].size < max(MIN_EPI, 4) or len(set(d["probes"][d["epi_idx"]])) < 2:
            continue
        acc = {c: [0.0, 0] for c in COMPOUNDS}             # concordance, pairs
        lift = {c: {k: [] for k in K_LIST} for c in COMPOUNDS}
        r_vals = []
        for _s, off, tmask, aff_off, str_off, r in fold_rows(d):
            sc = compound_scores(aff_off, str_off, r)
            r_vals.append(r)
            prev = tmask.sum() / tmask.size
            for c in COMPOUNDS:
                s = sc[c]
                cse, tot = _conc(s[tmask], s[~tmask])
                acc[c][0] += cse
                acc[c][1] += tot
                order = np.argsort(-s)
                for k in K_LIST:
                    kk = min(k, s.size)
                    hit = tmask[order[:kk]].sum()
                    lift[c][k].append((hit / kk) / prev if prev > 0 else np.nan)
        for c in COMPOUNDS:
            if acc[c][1] == 0:
                continue
            row = dict(patient=d["pat"], compound=c, auc=acc[c][0] / acc[c][1],
                       seed_strength_pct=float(np.mean(r_vals)))
            for k in K_LIST:
                row[f"lift_at_{k}"] = float(np.nanmean(lift[c][k]))
            rows.append(row)
    return pd.DataFrame(rows)


def logistic_lopo(packs, auc_lookup):
    """Pool per-fold off-node rows; LOPO-train logistic; return held-out AUC + calibration."""
    perpat = {}                                            # pat -> (X[n,2 z-scored], y)
    for d in packs.values():
        if d["epi_idx"].size < max(MIN_EPI, 4) or len(set(d["probes"][d["epi_idx"]])) < 2:
            continue
        Xs, ys = [], []
        for _s, off, tmask, aff_off, str_off, _r in fold_rows(d):
            az = (aff_off - aff_off.mean()) / (aff_off.std() + 1e-12)
            sz = (str_off - str_off.mean()) / (str_off.std() + 1e-12)
            Xs.append(np.column_stack([az, sz]))
            ys.append(tmask.astype(float))
        if Xs:
            perpat[d["pat"]] = (np.vstack(Xs), np.concatenate(ys))
    pats = list(perpat)
    held_rows, cal_P, cal_y = [], [], []
    for hp in pats:
        Xtr = np.vstack([perpat[q][0] for q in pats if q != hp])
        ytr = np.concatenate([perpat[q][1] for q in pats if q != hp])
        w = fit_logistic(np.column_stack([np.ones(len(Xtr)), Xtr]), ytr)
        Xte, yte = perpat[hp]
        P = 1.0 / (1.0 + np.exp(-(np.column_stack([np.ones(len(Xte)), Xte]) @ w)))
        cse, tot = _conc(P[yte == 1], P[yte == 0])
        held_rows.append(dict(patient=hp, held_auc=cse / tot if tot else np.nan,
                              validation_auc=round(float(auc_lookup.get(hp, np.nan)), 3),
                              n_target=int(yte.sum()), n_off=int(len(yte)),
                              w_aff=float(w[1]), w_str=float(w[2])))
        cal_P.append(P)
        cal_y.append(yte)
    P = np.concatenate(cal_P)
    y = np.concatenate(cal_y)
    # reliability deciles + Brier
    bins = np.quantile(P, np.linspace(0, 1, 11))
    bins[-1] += 1e-9
    cal = []
    for b in range(10):
        m = (P >= bins[b]) & (P < bins[b + 1])
        if m.sum() >= 5:
            cal.append(dict(bin=b, p_pred_mean=float(P[m].mean()),
                            soz_observed=float(y[m].mean()), n=int(m.sum())))
    brier = float(np.mean((P - y) ** 2))
    brier_base = float(np.mean((y.mean() - y) ** 2))
    return pd.DataFrame(held_rows), pd.DataFrame(cal), brier, brier_base, float(y.mean())


def deployment_probabilities(packs, auc_lookup):
    """Train logistic on the OTHER patients (LOPO), emit P(SOZ) for each unmarked contact."""
    # rebuild pooled per-patient training rows (same as logistic_lopo)
    perpat = {}
    for d in packs.values():
        if d["epi_idx"].size < max(MIN_EPI, 4) or len(set(d["probes"][d["epi_idx"]])) < 2:
            continue
        Xs, ys = [], []
        for _s, off, tmask, aff_off, str_off, _r in fold_rows(d):
            az = (aff_off - aff_off.mean()) / (aff_off.std() + 1e-12)
            sz = (str_off - str_off.mean()) / (str_off.std() + 1e-12)
            Xs.append(np.column_stack([az, sz]))
            ys.append(tmask.astype(float))
        if Xs:
            perpat[d["pat"]] = (np.vstack(Xs), np.concatenate(ys))
    out = []
    for d in packs.values():
        pat = d["pat"]
        train = [q for q in perpat if q != pat]
        if not train:
            continue
        Xtr = np.vstack([perpat[q][0] for q in train])
        ytr = np.concatenate([perpat[q][1] for q in train])
        w = fit_logistic(np.column_stack([np.ones(len(Xtr)), Xtr]), ytr)
        epi, strength, O, N = d["epi"], d["strength"], d["O"], d["N"]
        unmarked = np.array([j for j in range(N) if not epi[j]])
        aff = _resid(O[:, d["epi_idx"]].mean(axis=1), strength)
        au, su = aff[unmarked], strength[unmarked]
        az = (au - au.mean()) / (au.std() + 1e-12)
        sz = (su - su.mean()) / (su.std() + 1e-12)
        P = 1.0 / (1.0 + np.exp(-(np.column_stack([np.ones(len(au)), az, sz]) @ w)))
        try:
            reg = load_channel_regions(pat)["region"].to_numpy(object)
        except Exception:
            reg = np.array(["?"] * N)
        order = np.argsort(-P)
        for rank, ii in enumerate(order[:TOP_CAND], 1):
            j = unmarked[ii]
            out.append(dict(patient=pat, validation_auc=round(float(auc_lookup.get(pat, np.nan)), 3),
                            responder=bool(auc_lookup.get(pat, 0) >= RESP_THR), rank=rank,
                            contact=str(d["channels"][j]), region=str(reg[j]) if j < len(reg) else "?",
                            shaft=str(d["probes"][j]), P_soz=round(float(P[ii]), 3)))
    return pd.DataFrame(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ml = pd.read_csv(ML / "marker_library_per_patient.csv")
    mld = ml[(ml.band == BAND) & (ml.marker == MARKER)]
    auc_lookup = dict(zip(mld.patient, mld.auc_resid))
    t0 = time.time()

    packs = {p: patient_pack(p) for p in es.COHORT}
    packs = {k: v for k, v in packs.items() if v is not None}

    comp = eval_compounds(packs)
    comp.to_csv(OUT / "compound_auc_per_patient.csv", index=False)
    coh_rows = []
    for c in COMPOUNDS:
        d = comp[comp.compound == c]
        resp = d[d.patient.map(lambda p: auc_lookup.get(p, 0) >= RESP_THR)]
        coh_rows.append(dict(compound=c, n=len(d),
                             med_auc=float(d.auc.median()), n_above_half=int((d.auc > 0.5).sum()),
                             med_lift5=float(d.lift_at_5.median()),
                             med_lift5_responders=float(resp.lift_at_5.median()) if len(resp) else np.nan))
    coh = pd.DataFrame(coh_rows)
    coh.to_csv(OUT / "compound_auc_cohort.csv", index=False)

    held, cal, brier, brier_base, prev = logistic_lopo(packs, auc_lookup)
    held.to_csv(OUT / "logistic_lopo_per_patient.csv", index=False)
    cal.to_csv(OUT / "calibration_curve.csv", index=False)
    cand = deployment_probabilities(packs, auc_lookup)
    cand.to_csv(OUT / "candidates_with_prob.csv", index=False)

    # ----------------------------------------------------------------- console
    print(f"[audit_104] {BAND} {MARKER}: {len(packs)} patients ({time.time()-t0:.1f}s)\n")
    print("(1-5) COMPOUND off-shaft distant-discovery (cohort):")
    print("    compound     med_AUC  n>0.5   lift@5(all)  lift@5(resp)")
    for _, r in coh.iterrows():
        print(f"    {r.compound:11s}  {r.med_auc:.3f}   {int(r.n_above_half):2d}/{int(r.n)}    "
              f"{r.med_lift5:5.2f}        {r.med_lift5_responders:5.2f}")
    print("\n    per-patient regime (seed strength pct) + affinity / strength / switch AUC:")
    piv = comp.pivot(index="patient", columns="compound", values="auc")
    rseed = comp.groupby("patient").seed_strength_pct.first()
    for pat in piv.index:
        gate = "HUB->strength" if rseed[pat] >= R_THR else "community->affinity"
        print(f"      {pat:7s} r={rseed[pat]:.2f} [{gate:20s}]  "
              f"aff={piv.loc[pat,'affinity']:.3f} str={piv.loc[pat,'strength']:.3f} "
              f"switch={piv.loc[pat,'switch']:.3f}")
    sw = comp[comp.compound == "switch"]
    print(f"    SWITCH cohort: med_AUC={sw.auc.median():.3f} {int((sw.auc>0.5).sum())}/{len(sw)}>0.5 "
          f"med_lift@5={sw.lift_at_5.median():.2f}")

    print("\n(6) LOGISTIC leave-one-PATIENT-out P(SOZ):")
    print(f"    held-out AUC: median={held.held_auc.median():.3f} mean={held.held_auc.mean():.3f} "
          f"{int((held.held_auc>0.5).sum())}/{len(held)}>0.5")
    print(f"    mean weights: w_affinity={held.w_aff.mean():+.3f}  w_strength={held.w_str.mean():+.3f}")
    print(f"    calibration: Brier={brier:.4f} vs prevalence-baseline={brier_base:.4f} "
          f"(prevalence={prev:.3f})  {'BETTER' if brier < brier_base else 'NOT better'}")
    print("    reliability (predicted P -> observed SOZ fraction):")
    for _, r in cal.iterrows():
        print(f"      P~{r.p_pred_mean:.3f}  observed={r.soz_observed:.3f}  (n={int(r.n)})")
    print("\n    per-patient held-out AUC (logistic) vs validated affinity AUC:")
    for _, r in held.sort_values("held_auc", ascending=False).iterrows():
        print(f"      {r.patient:7s} logistic={r.held_auc:.3f}  affinity={r.validation_auc:.2f}")
    print(f"\n    top P(SOZ) candidates (LOPO), highest per responder:")
    for pat in [p for p in es.COHORT if auc_lookup.get(p, 0) >= RESP_THR and p in set(cand.patient)]:
        r = cand[cand.patient == pat].iloc[0]
        print(f"      {pat:7s} {r.region:26s} shaft {r.shaft:4s} P(SOZ)={r.P_soz:.3f}")
    print(f"\n[audit_104] -> {OUT}")


if __name__ == "__main__":
    main()
