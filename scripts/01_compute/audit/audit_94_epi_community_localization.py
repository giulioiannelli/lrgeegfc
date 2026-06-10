#!/usr/bin/env python3
"""Audit 94 — LABEL-FREE cross-patient SOZ localization at the COMMUNITY level.

The marker so far (audit_91–93) is seed-based: its features are affinity to the
KNOWN SOZ, so it cannot run on a patient with zero labels. The user's question:
given what we learn from single patients, is anything possible cross-patient
WITHOUT the target patient's labels?

The answer can only live at the COMMUNITY level, not the node level. Node-intrinsic
SOZ-ness does not transfer (audit_80; strength at chance in audit_93). What we PROVED
transfers is a collective property: the SOZ is a strength-independent diffusion
community. So the transferable invariant is "what a SOZ community looks like", and
the pipeline is:

  (1) Build a LABEL-FREE beyond-strength co-diffusion matrix A^z: per pair (i,j),
      z-score of ρ_ij(τ) against the matched-strength surrogate ensemble (the
      surrogate is the strength reference — no labels). A^z>0 = co-diffuse more than
      strength forces.
  (2) Detect communities UNSUPERVISED (spectral clustering on the positive part of
      A^z) — communities of beyond-strength co-diffusion, NOT strength/probe blocks.
  (3) Per community, LABEL-FREE features: internal coherence (mean A^z within),
      interface segregation (internal − external A^z), size, spatial compactness,
      mean node strength (to show the rule is NOT picking strength-blobs), n shafts.
  (4) Cross-patient LOPO: train a community-level logistic (features → P(this is the
      SOZ community)) on 8 patients, apply LABEL-FREE to the held-out patient's
      communities, rank, take the top community, measure SOZ overlap.

Verdict = does the top-ranked community capture the SOZ above a random-community
baseline, cross-patient, with zero labels in the target? Bounded by the responder
subset (a patient whose SOZ is NOT a diffusion community — e.g. Pat_10, epi-as-hubs —
will be missed by ANY diffusion-community method; that is a real false-negative
floor, reported honestly).

Critical preamble
=================
(1) Claim: a community-level rule learned from single patients localizes the SOZ in
    a held-out patient WITHOUT its labels, above chance.
(2) Null: random-community baseline (precision = SOZ prevalence; SOZ-community rank
    uniform in 1..K) + the high-γ internal negative.
(3) Strongest alternative: the rule just picks the highest-strength or largest
    community (trivial), or the clustering leaks labels.
(4) Reach: A^z is matched-strength-z (label-free strength control); the classifier
    carries node strength as a feature so its weight is auditable; the SOZ labels
    enter ONLY in training targets and final eval, never in the target patient's
    clustering or features.
(5) Falsification: if the top-ranked community's SOZ precision ≈ prevalence and the
    SOZ-community rank ≈ uniform, label-free cross-patient localization fails — the
    marker stays seed-based and we say so.

Outputs (``data/audit/epi_community_localization/``)
    community_features.csv      per (patient, band, community): label-free feats + SOZ overlap
    localization_per_patient.csv per (patient, band): top-community precision/recall, SOZ rank
    localization_cohort.csv      per (band)
    README.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.cluster import SpectralClustering
from sklearn.linear_model import LogisticRegression

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.surrogate import load_or_compute_eigs_at_path

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI, N_TAU  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_community_localization"
OUT.mkdir(parents=True, exist_ok=True)

K_COMM = 8            # unsupervised communities per patient (fixed a priori)
TAU_IDXS = (2, 3, 4)  # mid-range τ averaged for A^z (fixed a priori; signal is broad)
MIN_COMM = 3          # min nodes to score a community
MIN_OVERLAP = 3       # min SOZ nodes for a community to be the "SOZ community" target
CLF_FEATS = ["internal_z", "seg_z", "size_frac", "compactness", "strength_z", "shaft_frac"]


def rho_mean_over_tau(lam, U, taus):
    """Mean propagator over a τ subset (one (N,N) matrix)."""
    lam = np.clip(np.asarray(lam, float), 0.0, None)
    acc = None
    for t in taus:
        w = np.exp(-t * lam)
        Z = w.sum()
        rho = (U * w) @ U.T / Z
        acc = rho if acc is None else acc + rho
    return acc / len(taus)


def co_diffusion_z(W, evals, evecs, taus):
    """Label-free beyond-strength co-diffusion matrix A^z (per-pair z of ρ vs the
    matched-strength surrogate ensemble)."""
    lam_o, U_o = np.linalg.eigh(np.diag(W.sum(1)) - W)
    obs = rho_mean_over_tau(lam_o, U_o, taus)
    N = W.shape[0]
    s1 = np.zeros((N, N))
    s2 = np.zeros((N, N))
    n = 0
    for r in range(evals.shape[0]):
        if not np.isfinite(evals[r]).all():
            continue
        rr = rho_mean_over_tau(evals[r], evecs[r], taus)
        s1 += rr
        s2 += rr * rr
        n += 1
    if n < 5:
        return None
    mu = s1 / n
    var = np.maximum(s2 / n - mu * mu, 0.0)
    sd = np.sqrt(var)
    with np.errstate(divide="ignore", invalid="ignore"):
        Az = np.where(sd > 0, (obs - mu) / sd, 0.0)
    np.fill_diagonal(Az, 0.0)
    return Az


def detect_communities(Az, k):
    """Unsupervised spectral clustering on the positive part of A^z."""
    aff = np.maximum(Az, 0.0)
    aff = 0.5 * (aff + aff.T)
    k = min(k, Az.shape[0] - 1)
    try:
        labels = SpectralClustering(
            n_clusters=k, affinity="precomputed", assign_labels="discretize",
            random_state=0).fit_predict(aff)
    except Exception:
        return None
    return labels


def community_features(Az, labels, epi, strength, coords, probe):
    """Per-community label-free features + SOZ overlap (label used only downstream)."""
    N = Az.shape[0]
    srank = np.argsort(np.argsort(strength)) / max(1, N - 1)   # within-patient
    rows = []
    for c in np.unique(labels):
        idx = np.where(labels == c)[0]
        if idx.size < MIN_COMM:
            continue
        ext = np.where(labels != c)[0]
        iu = np.triu_indices(idx.size, 1)
        internal = float(Az[np.ix_(idx, idx)][iu].mean()) if iu[0].size else 0.0
        external = float(Az[np.ix_(idx, ext)].mean()) if ext.size else 0.0
        ok = np.isfinite(coords[idx]).all(1)
        if ok.sum() >= 2:
            cc = coords[idx][ok]
            compact = float(np.sqrt(((cc - cc.mean(0)) ** 2).sum(1).mean()))
        else:
            compact = np.nan
        rows.append({
            "community": int(c), "size": int(idx.size),
            "internal_z": internal, "seg_z": internal - external,
            "size_frac": idx.size / N,
            "compactness": compact,
            "strength_z": float(srank[idx].mean()),
            "shaft_frac": len(set(probe[idx])) / idx.size,
            "n_soz": int(epi[idx].sum()),
            "soz_precision": float(epi[idx].mean()),
        })
    return rows


def per_patient_band(pat, band, n_surr, swap, verbose):
    try:
        W = load_phase_fc(pat, "rest_post", band)
    except Exception as e:
        if verbose:
            print(f"[audit_94] SKIP {pat}/{band}: {e}")
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    probe = np.asarray(pm.probes, object)
    if epi.sum() < MIN_EPI or (~epi).sum() < 2:
        return None
    reg = load_channel_regions(pat)
    if len(reg) != N:
        return None
    coords = reg[["x", "y", "z"]].to_numpy(float)
    lam_o = np.linalg.eigvalsh(np.diag(W.sum(1)) - W)
    lmax = float(lam_o[-1])
    if lmax <= 0:
        return None
    taus = np.geomspace(1.0 / lmax, 10.0 / lmax, N_TAU)[list(TAU_IDXS)]

    evals, evecs = load_or_compute_eigs_at_path(
        es.surr_eig_path("full", pat, band, "rest_post", n_surr, swap),
        W, n_surr, swap, es.cell_rng(pat, band, "rest_post", "full"),
        verbose=verbose)
    Az = co_diffusion_z(W, evals, evecs, taus)
    if Az is None:
        return None
    labels = detect_communities(Az, K_COMM)
    if labels is None:
        return None
    feats = community_features(Az, labels, epi, W.sum(1), coords, probe)
    if not feats:
        return None
    df = pd.DataFrame(feats)
    df["patient"] = pat
    df["band"] = band
    df["n_soz_total"] = int(epi.sum())
    df["prevalence"] = float(epi.mean())
    return df


def lopo_localize(comm):
    """Cross-patient LOPO: train community-level SOZ classifier on 8 patients,
    apply LABEL-FREE to the held-out patient's communities, rank, score top."""
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(comm.band)]:
        cb = comm[comm.band == band].copy()
        # training target: the SOZ community = max n_soz in patient, if >= MIN_OVERLAP
        cb["is_soz_comm"] = 0
        for pat, g in cb.groupby("patient"):
            if g.n_soz.max() >= MIN_OVERLAP:
                cb.loc[g.index[g.n_soz.values.argmax()], "is_soz_comm"] = 1
        pats = sorted(set(cb.patient))
        for tgt in pats:
            tr = cb[cb.patient != tgt]
            te = cb[cb.patient == tgt].copy()
            if tr.is_soz_comm.sum() < 2 or te.empty or te.n_soz.max() < MIN_OVERLAP:
                continue
            Xtr = np.nan_to_num(tr[CLF_FEATS].to_numpy(float), nan=0.0)
            Xte = np.nan_to_num(te[CLF_FEATS].to_numpy(float), nan=0.0)
            clf = LogisticRegression(max_iter=3000, C=1.0,
                                     class_weight="balanced")
            clf.fit(Xtr, tr.is_soz_comm.to_numpy(int))
            te["P_soz_comm"] = clf.predict_proba(Xte)[:, 1]
            te = te.sort_values("P_soz_comm", ascending=False).reset_index(drop=True)
            top = te.iloc[0]
            # rank the genuine SOZ community (max n_soz) in the classifier order
            soz_row = te.n_soz.values.argmax()
            K = len(te)
            out.append({
                "patient": tgt, "band": band, "K": K,
                "prevalence": float(te.prevalence.iloc[0]),
                "top_precision": float(top.soz_precision),
                "top_recall": float(top.n_soz / te.n_soz_total.iloc[0]),
                "top_is_soz_comm": int(soz_row == 0),
                "soz_comm_rank": int(soz_row + 1),
                "rank_frac": float((soz_row + 1) / K),
                "top_size": int(top["size"]),
                "lift": float(top.soz_precision / max(te.prevalence.iloc[0], 1e-6)),
            })
    return pd.DataFrame(out)


def cohort(loc):
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(loc.band)]:
        d = loc[loc.band == band]
        prec = d.top_precision.to_numpy()
        prev = d.prevalence.to_numpy()
        wp = np.nan
        if len(d) >= 3:
            try:
                _, wp = wilcoxon(prec - prev, alternative="greater")
            except Exception:
                pass
        out.append({
            "band": band, "n_patients": len(d),
            "med_top_precision": float(np.median(prec)),
            "med_prevalence": float(np.median(prev)),
            "med_lift": float(np.median(d.lift)),
            "med_top_recall": float(np.median(d.top_recall)),
            "frac_soz_top1": float(d.top_is_soz_comm.mean()),
            "frac_soz_top2": float((d.soz_comm_rank <= 2).mean()),
            "med_soz_rank": float(np.median(d.soz_comm_rank)),
            "p_prec_gt_prev": wp,
        })
    return pd.DataFrame(out)


def write_readme(coh, runtime):
    L = ["---", "name: epi_community_localization",
         "scope: label_free_cross_patient_SOZ_community_localization",
         "era: COHORT_N10 / IMCOH_ABS",
         f"date: {time.strftime('%Y-%m-%d')}",
         "build_script: scripts/01_compute/audit/audit_94_epi_community_localization.py",
         "---", "",
         "# Label-free cross-patient SOZ localization at the community level", "",
         "**Head.** Detect beyond-strength co-diffusion communities UNSUPERVISED in a "
         "held-out patient, score each with a community-level SOZ classifier trained "
         "on the OTHER patients, take the top community — using ZERO labels in the "
         "target. `top_precision` = fraction of the top community that is SOZ; "
         "`lift` = precision / prevalence; `frac_soz_top1` = fraction of patients "
         "where the true SOZ community is ranked #1. Bounded by the responder subset "
         "(non-community SOZ, e.g. Pat_10, missed by construction).", "",
         "## Cohort (label-free, leave-one-patient-out)", "",
         "| band | top precision | prevalence | lift | recall | SOZ#1 | SOZ≤2 | "
         "med rank/K | p(prec>prev) |",
         "|---|---|---|---|---|---|---|---|---|"]
    for band in [b for b in es.ALL_BANDS if b in set(coh.band)]:
        r = coh[coh.band == band].iloc[0]
        L.append(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {r.med_top_precision:.3f} "
                 f"| {r.med_prevalence:.3f} | {r.med_lift:.2f}× | {r.med_top_recall:.2f} "
                 f"| {r.frac_soz_top1:.2f} | {r.frac_soz_top2:.2f} "
                 f"| {r.med_soz_rank:.0f} | {r.p_prec_gt_prev:.4f} |")
    L += ["", "## Reading",
          "- `lift > 1` and `p(prec>prev) < 0.05` ⇒ the top label-free community is "
          "SOZ-enriched above chance. `SOZ#1` high ⇒ the true SOZ community is the "
          "classifier's top pick. This is the ONLY cross-patient, label-free signal "
          "in the investigation; everything else is seed-based.",
          "- The classifier carries node strength as a feature — inspect its weight "
          "(`community_features.csv` + the model) to confirm it is NOT picking "
          "strength/size blobs.",
          "- Per patient in `localization_per_patient.csv` (Pat_10-type non-"
          "responders read directly). high-γ = internal negative.",
          f"- Communities K={K_COMM}, τ idx {TAU_IDXS}; wall-clock {runtime:.1f}s"]
    (OUT / "README.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--n-surrogates", type=int, default=es.N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=es.SWAP_FACTOR)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    frames = []
    for band in args.bands:
        for pat in args.patients:
            tc = time.time()
            df = per_patient_band(pat, band, args.n_surrogates, args.swap_factor,
                                  args.verbose)
            if df is not None:
                frames.append(df)
                soz = df.loc[df.n_soz.idxmax()]
                print(f"[audit_94] {pat}/{band}: K={df.community.nunique()} "
                      f"bestSOZcomm prec={soz.soz_precision:.2f} "
                      f"n_soz={int(soz.n_soz)}/{int(soz.n_soz_total)} ({time.time()-tc:.1f}s)")
    comm = pd.concat(frames, ignore_index=True)
    comm.to_csv(OUT / "community_features.csv", index=False)
    loc = lopo_localize(comm)
    loc.to_csv(OUT / "localization_per_patient.csv", index=False)
    coh = cohort(loc)
    coh.to_csv(OUT / "localization_cohort.csv", index=False)
    runtime = time.time() - t0
    write_readme(coh, runtime)

    print("\n[audit_94] LABEL-FREE cross-patient localization (top community):")
    for band in [b for b in es.ALL_BANDS if b in set(coh.band)]:
        r = coh[coh.band == band].iloc[0]
        print(f"  {band:10s}: precision={r.med_top_precision:.2f} "
              f"(prev {r.med_prevalence:.2f}, lift {r.med_lift:.1f}×) "
              f"SOZ#1={r.frac_soz_top1:.0%} SOZ≤2={r.frac_soz_top2:.0%} "
              f"p={r.p_prec_gt_prev:.3f}")
    print(f"[audit_94] done in {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
