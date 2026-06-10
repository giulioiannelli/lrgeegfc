#!/usr/bin/env python3
"""Audit 93 — VI2/VI3 + VII2/VII4: interpretable per-node SOZ classifier, its
nulls and calibration, and discovery of occult (unlabelled) SOZ candidates.

Builds on the validated relational marker (audit_91 features, audit_92 within-
patient recovery). Here the marker becomes an INTERPRETABLE cross-patient
classifier that assigns every contact a calibrated P(SOZ | node):

  - Inputs: the matched-strength-z relational features (`z_f_aff_mean` affinity to
    the SOZ community, `z_f_seg_mean` interface segregation) — both strength-
    orthogonal and within-patient-normalised, so they pool across patients.
  - Model: logistic regression (interpretable; coefficients reported as odds
    ratios), leave-ONE-PATIENT-out (LOPO) so a patient's own nodes never train
    its own scores → out-of-fold P(SOZ) for every contact.
  - Baselines: strength-only (hubness) and strength+relational, to show the
    relational features carry the signal and strength adds little.
  - Nulls (VI3): within-patient label permutation (R_PERM) → null LOPO AUC → p.
  - Calibration: reliability bins on pooled out-of-fold P(SOZ).

Discovery (VII2/VII4): apply the LOPO model to a patient's UNLABELLED contacts,
rank by P(SOZ), and propose top candidates as occult-SOZ HYPOTHESES (no resection/
outcome ground truth exists in this cohort). Each candidate carries: per-band
P(SOZ), cross-band agreement over the airtight bands (δ/β/low-γ), distance to the
nearest labelled SOZ and anatomical region (plausibility), and the feature values
that made it fire (WHY — affinity vs interface), so it is explainable.

Critical preamble
=================
(1) Claim: an interpretable logistic on relational propagator features predicts
    SOZ membership cross-patient (LOPO) above chance and above strength, and its
    P(SOZ) is calibrated enough to rank occult candidates.
(2) Null: within-patient label-permutation LOPO AUC (R_PERM) → p; strength-only
    baseline; chance 0.5.
(3) Strongest alternative: cross-patient transfer is hubness/geometry (audit_80
    showed raw cross-patient markers fail) — strength alone transfers.
(4) Reach: features are matched-strength-z (strength-orthogonal) and within-patient
    normalised (so cross-patient pooling is fair); LOPO prevents within-patient
    leakage; the permutation null breaks the label↔feature link while keeping the
    feature geometry.
(5) Falsification: if LOPO AUC ≈ permutation null or ≈ strength-only, the marker is
    within-patient only — discovery then falls back to within-patient ranks (stated
    honestly), not a calibrated cross-patient P(SOZ).

Reads ``data/audit/epi_marker_relational/node_relational_features.csv`` (audit_91).
Outputs (``data/audit/epi_marker_relational/``)
    classifier_lopo.csv         per (band, model): AUC, PR-AUC, OR, perm-p
    classifier_calibration.csv  per (band, bin): mean P vs observed frequency
    occult_candidates.csv       per (patient, node): per-band + cross-band P(SOZ)
    README_classifier.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.regions import load_channel_regions

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
import _epi_stratify as es  # type: ignore

SRC = ROOT / "data" / "audit" / "epi_marker_relational" / "node_relational_features.csv"
OUT = ROOT / "data" / "audit" / "epi_marker_relational"

REL_FEATS = ["z_f_aff_mean", "z_f_seg_mean"]      # interpretable relational marker
MODELS = {
    "strength": ["z_strength"],
    "relational": REL_FEATS,
    "strength+relational": ["z_strength"] + REL_FEATS,
}
AIRTIGHT = ["delta", "beta", "low_gamma"]          # C3+C5-survived bands
R_PERM = 300
N_CAL_BINS = 8
TOP_K = 8                                           # occult candidates per patient


def _zscore_within(df, col):
    out = np.full(len(df), np.nan)
    for _, idx in df.groupby("patient").groups.items():
        v = df.loc[idx, col].to_numpy(float)
        sd = np.nanstd(v)
        out[df.index.get_indexer(idx)] = (v - np.nanmean(v)) / (sd if sd > 1e-9 else 1.0)
    return out


def lopo_proba(X, y, groups):
    """Out-of-fold P(SOZ) (leave-one-patient-out) + per-fold coefficients."""
    P = np.full(len(y), np.nan)
    coefs = []
    for g in np.unique(groups):
        tr = groups != g
        te = groups == g
        yt = y[tr].astype(bool)
        if yt.sum() < 2 or (~yt).sum() < 2:
            continue
        clf = LogisticRegression(max_iter=3000, C=1.0)
        clf.fit(X[tr], y[tr])
        P[te] = clf.predict_proba(X[te])[:, 1]
        coefs.append(clf.coef_[0])
    return P, (np.array(coefs) if coefs else np.zeros((0, X.shape[1])))


def _pooled_auc(y, P):
    ok = np.isfinite(P)
    if y[ok].sum() < 1 or (y[ok] == 0).sum() < 1 or np.unique(P[ok]).size < 2:
        return np.nan
    return float(roc_auc_score(y[ok], P[ok]))


def fit_band(fb, band, rng):
    fb = fb.reset_index(drop=True)
    y = fb["is_epi"].to_numpy(int)
    groups = fb["patient"].to_numpy()
    res = {"band": band, "n_nodes": len(fb), "n_epi": int(y.sum()),
           "prevalence": float(y.mean())}
    P_store = {}
    for name, cols in MODELS.items():
        X = fb[cols].to_numpy(float)
        X = np.nan_to_num(X, nan=0.0)
        P, coefs = lopo_proba(X, y, groups)
        P_store[name] = P
        res[f"auc_{name}"] = _pooled_auc(y, P)
        ok = np.isfinite(P)
        res[f"prauc_{name}"] = (float(average_precision_score(y[ok], P[ok]))
                                if ok.sum() and y[ok].sum() else np.nan)
        if coefs.size:
            for j, c in enumerate(cols):
                res[f"OR_{name}_{c}"] = float(np.exp(np.median(coefs[:, j])))
    # permutation null for the relational model
    Xr = np.nan_to_num(fb[REL_FEATS].to_numpy(float), nan=0.0)
    obs = res["auc_relational"]
    null = []
    for _ in range(R_PERM):
        yp = y.copy()
        for g in np.unique(groups):
            m = groups == g
            yp[m] = rng.permutation(yp[m])
        Pp, _ = lopo_proba(Xr, yp, groups)
        a = _pooled_auc(yp, Pp)
        if np.isfinite(a):
            null.append(a)
    null = np.array(null)
    res["perm_p_relational"] = (float((null >= obs).mean())
                                if null.size and np.isfinite(obs) else np.nan)
    res["perm_null_med"] = float(np.median(null)) if null.size else np.nan
    return res, P_store["relational"], y


def calibration_bins(y, P, band):
    ok = np.isfinite(P)
    y, P = y[ok], P[ok]
    if P.size < N_CAL_BINS:
        return []
    edges = np.quantile(P, np.linspace(0, 1, N_CAL_BINS + 1))
    edges[-1] += 1e-9
    rows = []
    for b in range(N_CAL_BINS):
        m = (P >= edges[b]) & (P < edges[b + 1])
        if m.sum() < 1:
            continue
        rows.append({"band": band, "bin": b, "n": int(m.sum()),
                     "mean_pred": float(P[m].mean()),
                     "obs_freq": float(y[m].mean())})
    return rows


def bootstrap_candidate_ci(feat, rng, B=120):
    """Cluster-bootstrap (resample TRAINING PATIENTS with replacement) the LOPO
    P(SOZ) → per-node cross-band P(SOZ) CI + top-K rank stability. The relevant
    uncertainty for cross-patient transfer is *which patients* trained the model,
    so we resample patients, not nodes."""
    epi_map, band_data = {}, {}
    for band in AIRTIGHT:
        fb = feat[feat.band == band].reset_index(drop=True)
        band_data[band] = (fb["patient"].to_numpy(), fb["node"].to_numpy(int),
                           fb["is_epi"].to_numpy(int),
                           np.nan_to_num(fb[REL_FEATS].to_numpy(float), nan=0.0))
        for p, n, e in zip(fb.patient, fb.node, fb.is_epi):
            epi_map[(p, int(n))] = bool(e)
    pats = sorted(set(feat.patient))
    samples, topk = {}, {}
    for _ in range(B):
        cross = {}
        for band in AIRTIGHT:
            P_arr, N_arr, y_arr, X = band_data[band]
            for tgt in pats:
                tr_pats = [p for p in pats if p != tgt]
                boot = rng.choice(tr_pats, size=len(tr_pats), replace=True)
                tr = np.concatenate([np.where(P_arr == p)[0] for p in boot])
                te = np.where(P_arr == tgt)[0]
                ytr = y_arr[tr]
                if ytr.sum() < 2 or (ytr == 0).sum() < 2 or te.size == 0:
                    continue
                clf = LogisticRegression(max_iter=3000, C=1.0)
                clf.fit(X[tr], ytr)
                Pte = clf.predict_proba(X[te])[:, 1]
                for j, idx in enumerate(te):
                    cross.setdefault((tgt, int(N_arr[idx])), []).append(Pte[j])
        meanP = {k: float(np.mean(v)) for k, v in cross.items()}
        for k, mp in meanP.items():
            samples.setdefault(k, []).append(mp)
        for tgt in pats:
            items = sorted(((k, meanP[k]) for k in meanP
                            if k[0] == tgt and not epi_map.get(k, False)),
                           key=lambda x: -x[1])
            for k, _ in items[:TOP_K]:
                topk[k] = topk.get(k, 0) + 1
    out = {}
    for k, v in samples.items():
        v = np.array(v)
        out[k] = (float(np.percentile(v, 5)), float(np.percentile(v, 95)),
                  topk.get(k, 0) / B)
    return out


def discover(feat, P_by_band, rng):
    """Rank each patient's UNLABELLED contacts by cross-band P(SOZ)."""
    # attach per-band out-of-fold P(SOZ) to the node table
    base = feat[["patient", "band", "node", "is_epi", "probe",
                 "z_f_aff_mean", "z_f_seg_mean"]].copy()
    base["P_soz"] = np.nan
    for band, (P, idx) in P_by_band.items():
        base.loc[idx, "P_soz"] = P
    # wide: one row per (patient,node), columns = P per airtight band
    rows = []
    for pat, gp in base.groupby("patient"):
        reg = load_channel_regions(pat)
        coords = reg[["x", "y", "z"]].to_numpy(float)
        epi_nodes = sorted(set(gp[gp.is_epi].node))
        epi_xyz = coords[epi_nodes] if epi_nodes else None
        per_node = {}
        for nd, gnd in gp.groupby("node"):
            d = {b: float(gnd[gnd.band == b].P_soz.mean()) for b in AIRTIGHT}
            vals = [d[b] for b in AIRTIGHT if np.isfinite(d[b])]
            is_epi = bool(gnd.is_epi.iloc[0])
            nd = int(nd)
            # distance to nearest labelled SOZ (exclude self)
            dist = np.nan
            if epi_xyz is not None and np.isfinite(coords[nd]).all():
                dd = np.sqrt(((epi_xyz - coords[nd]) ** 2).sum(1))
                dd = dd[dd > 1e-6]
                dist = float(dd.min()) if dd.size else np.nan
            per_node[nd] = {
                "patient": pat, "node": nd, "is_epi": is_epi,
                "probe": gnd.probe.iloc[0],
                "region": reg.region.iloc[nd] if nd < len(reg) else "?",
                "hemisphere": reg.hemisphere.iloc[nd] if nd < len(reg) else "?",
                "P_mean_airtight": float(np.mean(vals)) if vals else np.nan,
                "n_bands_top_decile": 0,
                "dist_nearest_soz_um": dist,
                "z_f_aff_mean": float(gnd.z_f_aff_mean.mean()),
                "z_f_seg_mean": float(gnd.z_f_seg_mean.mean()),
                **{f"P_{b}": d[b] for b in AIRTIGHT}}
        # cross-band agreement: count bands where node is in patient top decile
        for b in AIRTIGHT:
            pv = np.array([per_node[n][f"P_{b}"] for n in per_node])
            nds = list(per_node)
            ok = np.isfinite(pv)
            if ok.sum() < 5:
                continue
            thr = np.quantile(pv[ok], 0.9)
            for n, val in zip(nds, pv):
                if np.isfinite(val) and val >= thr:
                    per_node[n]["n_bands_top_decile"] += 1
        rows.extend(per_node.values())
    df = pd.DataFrame(rows)
    # rank UNLABELLED candidates within patient by cross-band mean P(SOZ)
    df["occult_rank"] = np.nan
    cand = df[~df.is_epi].copy()
    for pat, g in cand.groupby("patient"):
        order = g.sort_values("P_mean_airtight", ascending=False)
        df.loc[order.index, "occult_rank"] = np.arange(1, len(order) + 1)
    df["is_top_candidate"] = (~df.is_epi) & (df.occult_rank <= TOP_K)
    return df.sort_values(["patient", "is_epi", "occult_rank"])


def write_readme(clf, runtime):
    L = ["---", "name: epi_soz_classifier_discovery",
         "scope: VI2_VI3_VII2_interpretable_SOZ_classifier_and_occult_discovery",
         "era: COHORT_N10 / IMCOH_ABS",
         f"date: {time.strftime('%Y-%m-%d')}",
         "build_script: scripts/01_compute/audit/audit_93_epi_soz_classifier_discovery.py",
         "---", "",
         "# Interpretable SOZ classifier P(SOZ|node) + occult-candidate discovery",
         "",
         "**Head.** Logistic regression (LOPO, interpretable) on the strength-"
         "orthogonal relational features `z_f_aff_mean` (affinity to SOZ community) "
         "and `z_f_seg_mean` (interface segregation) → out-of-fold calibrated "
         "P(SOZ|node). Reported vs a strength-only baseline and a within-patient "
         "label-permutation null. Discovery: rank each patient's UNLABELLED "
         "contacts by cross-band P(SOZ) (δ/β/low-γ) → occult-SOZ hypotheses "
         "(no ground truth → hypotheses, see `occult_candidates.csv`).", "",
         "## Cross-patient LOPO classifier", "",
         "| band | AUC relational | PR-AUC (prev) | perm-p | AUC strength | "
         "AUC str+rel | OR aff | OR seg |",
         "|---|---|---|---|---|---|---|---|"]
    for band in [b for b in es.ALL_BANDS if b in set(clf.band)]:
        r = clf[clf.band == band].iloc[0]
        L.append(
            f"| {BRAIN_BAND_TEX_DICT.get(band, band)} "
            f"| {r.auc_relational:.3f} "
            f"| {r.prauc_relational:.3f} ({r.prevalence:.2f}) "
            f"| {r.perm_p_relational:.3f} "
            f"| {r.auc_strength:.3f} | {r['auc_strength+relational']:.3f} "
            f"| {r.get('OR_relational_z_f_aff_mean', np.nan):.2f} "
            f"| {r.get('OR_relational_z_f_seg_mean', np.nan):.2f} |")
    L += ["", "## Reading",
          "- AUC relational > perm null (perm-p small) AND ≥ AUC strength ⇒ the "
          "relational marker transfers cross-patient and is not hubness. Odds "
          "ratios > 1 ⇒ higher affinity / segregation raise P(SOZ).",
          "- Calibration in `classifier_calibration.csv` (mean predicted vs "
          "observed SOZ frequency per bin).",
          "- Occult candidates in `occult_candidates.csv`: per-patient top-"
          f"{TOP_K} unlabelled contacts by cross-band P(SOZ), with distance to "
          "nearest labelled SOZ, region, and the feature values explaining the "
          "call. **Hypotheses only — no resection/outcome ground truth.**",
          f"- Permutations {R_PERM}; wall-clock {runtime:.1f}s"]
    (OUT / "README_classifier.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    args = ap.parse_args()
    if not SRC.exists():
        raise SystemExit("[audit_93] run audit_91 first")
    feat = pd.read_csv(SRC)
    feat["z_strength"] = _zscore_within(feat, "strength_mean")
    rng = np.random.default_rng(20260608)
    t0 = time.time()

    clf_rows, cal_rows = [], []
    P_by_band = {}
    for band in args.bands:
        fb = feat[feat.band == band].reset_index()      # 'index' = row in feat
        if fb.empty:
            continue
        res, P_rel, y = fit_band(fb, band, rng)
        clf_rows.append(res)
        cal_rows.extend(calibration_bins(y, P_rel, band))
        P_by_band[band] = (P_rel, fb["index"].to_numpy())
        print(f"[audit_93] {band}: AUC rel={res['auc_relational']:.3f} "
              f"(perm-p {res['perm_p_relational']:.3f}) "
              f"str={res['auc_strength']:.3f} "
              f"str+rel={res['auc_strength+relational']:.3f}")

    clf = pd.DataFrame(clf_rows)
    clf.to_csv(OUT / "classifier_lopo.csv", index=False)
    pd.DataFrame(cal_rows).to_csv(OUT / "classifier_calibration.csv", index=False)

    cand = discover(feat, P_by_band, rng)
    ci = bootstrap_candidate_ci(feat, rng)
    cand["P_mean_lo"] = [ci.get((r.patient, int(r.node)), (np.nan,) * 3)[0]
                         for r in cand.itertuples()]
    cand["P_mean_hi"] = [ci.get((r.patient, int(r.node)), (np.nan,) * 3)[1]
                         for r in cand.itertuples()]
    cand["topk_stability"] = [ci.get((r.patient, int(r.node)), (np.nan,) * 3)[2]
                              for r in cand.itertuples()]
    cand.to_csv(OUT / "occult_candidates.csv", index=False)
    runtime = time.time() - t0
    write_readme(clf, runtime)

    print("\n[audit_93] top occult candidates (cross-band P, low-epi patients):")
    low_epi = ["Pat_03", "Pat_07", "Pat_08", "Pat_06"]
    top = cand[(cand.is_top_candidate) & (cand.patient.isin(low_epi))]
    for _, r in top.sort_values("P_mean_airtight", ascending=False).head(10).iterrows():
        print(f"  {r.patient} node{int(r.node)} {r.probe:>4s} {str(r.region)[:22]:22s} "
              f"P={r.P_mean_airtight:.2f} [{r.P_mean_lo:.2f},{r.P_mean_hi:.2f}] "
              f"stab={r.topk_stability:.2f} dist={r.dist_nearest_soz_um/1000:.0f}mm "
              f"aff={r.z_f_aff_mean:+.1f} seg={r.z_f_seg_mean:+.1f}")
    print(f"[audit_93] done in {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
