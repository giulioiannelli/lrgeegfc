#!/usr/bin/env python3
"""Audit 91 — interpretable per-node RELATIONAL diffusion features for the SOZ marker.

The whole epi investigation's lesson (audit_85/86/87 node strength & recovery,
audit_88 signed Laplacian): NODE-INTRINSIC / strength features fail beyond hubness.
The audit_89/90 win is RELATIONAL — epileptic contacts form a strength-independent,
spatially-irreducible diffusion COMMUNITY. So the marker must be built on relational
propagator features, not node strength.

For each node i and the known epileptic set E, read the LRG propagator ρ(τ)=e^{-τL}/Z
relationally:

    f_aff  = mean_{j in E\\{i}} ρ_ij               affinity to the SOZ community
    f_seg  = mean_{j in E\\i} ρ_ij − mean_{j in N\\i} ρ_ij    interface segregation
    f_part = aff_E / (aff_E + aff_N)               share of diffusion going to SOZ

(N = non-epi.) An UNLABELED epileptic contact should diffuse preferentially into the
SOZ block (high f_aff), couple LESS to healthy tissue than its strength forces
(high f_seg — the V3 interface signal: audit_89 z_EN<0), and route a large fraction
of its diffusion to the SOZ block (high f_part). Every feature is z-scored against the
matched-strength surrogate ensemble (R=200, cached rest_post eigs, seed 20260511),
which fixes each node's strength exactly — so the MS-z feature is STRENGTH-ORTHOGONAL
by construction (a high-affinity hub is normalized away; only beyond-strength affinity
survives). Computed at all 6 τ via O(N·K) eigen-projections.

This script produces the FEATURES and a per-band epi-vs-non separation SCREEN
(within-patient leave-self-out: each node's own label is removed from E, so f_aff(i)
never sees i's label). The rigorous hide-out validation (mask a fraction of E, recover
the hidden) and the calibrated classifier live in audit_92.

Critical preamble
=================
(1) Claim: relational propagator features carry strength-ORTHOGONAL, per-node SOZ
    information (epi nodes are more diffusion-coupled to the SOZ community than
    healthy nodes, beyond strength).
(2) Null: matched-strength surrogate-z (per feature) + chance AUC 0.5 (cohort
    Wilcoxon of within-patient MW AUC).
(3) Strongest alternative: hubness — epi are high-strength so they diffuse to
    everything, including E. Also leakage: f_aff(i) seeing i's own label.
(4) Reach: MS-z fixes strength exactly (kills hubness); leave-self-out removes i's
    own label from E (kills self-leakage). The remaining f_aff is beyond-strength
    affinity to the OTHER labeled epi — exactly the deployable signal.
(5) Falsification: if no MS-z relational feature separates epi above chance in any
    band (cohort Wilcoxon), the relational bet fails like the node-intrinsic one,
    and there is no per-node marker.

Outputs (``data/audit/epi_marker_relational/``)
    node_relational_features.csv   per (patient, band, node): raw + MS-z f_aff/f_seg/f_part
    relational_separation.csv      per (band, feature): within-patient AUC, MS-z AUC, p
    README_relational.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.utils.surrogate import load_or_compute_eigs_at_path

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI, N_TAU  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_marker_relational"
OUT.mkdir(parents=True, exist_ok=True)

FEATS = ("f_aff", "f_seg", "f_part")
REL_COLS = [f"{f}_t{ti}" for f in FEATS for ti in range(N_TAU)]


def relational_feats(lam, U, epi_mask, taus):
    """Per-node relational propagator features at each τ, from a Laplacian
    eigendecomposition (identical call for observed graph and each surrogate).

    ρ = U diag(e^{-τλ}) U^T / Z (heat kernel, entrywise ≥ 0 since -L is Metzler).
    Affinity to a set S is computed as (U @ (w * U^T 1_S))/Z, O(N·K) per τ.
    """
    lam = np.clip(np.asarray(lam, float), 0.0, None)
    if lam.size == 0 or not np.isfinite(lam).all() or not np.isfinite(U).all():
        return None
    N = U.shape[0]
    e = epi_mask.astype(float)
    nE = float(e.sum())
    nN = float(N - nE)
    pE = U.T @ e                       # (K,)
    pAll = U.T @ np.ones(N)            # (K,)
    pN = pAll - pE
    nE_i = nE - e                      # epi count excluding self
    nN_i = nN - (1.0 - e)             # non-epi count excluding self
    out: dict[str, np.ndarray] = {}
    for ti, t in enumerate(taus):
        w = np.exp(-t * lam)
        Z = w.sum()
        rho_ii = (U * U * w).sum(1) / Z
        sumE = (U @ (w * pE)) / Z
        sumN = (U @ (w * pN)) / Z
        affE = (sumE - rho_ii * e) / np.where(nE_i > 0, nE_i, np.nan)
        affN = (sumN - rho_ii * (1.0 - e)) / np.where(nN_i > 0, nN_i, np.nan)
        denom = affE + affN
        out[f"f_aff_t{ti}"] = affE
        out[f"f_seg_t{ti}"] = affE - affN
        out[f"f_part_t{ti}"] = np.where(denom > 0, affE / denom, np.nan)
    return out


def _auc(a, b):
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return float("nan")
    return float(mannwhitneyu(a, b, alternative="two-sided").statistic
                 / (a.size * b.size))


def per_patient_band(pat, band, n_surr, swap, verbose):
    try:
        W = load_phase_fc(pat, "rest_post", band)
    except Exception as e:
        if verbose:
            print(f"[audit_91] SKIP {pat}/{band}: {e}")
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    probe = np.asarray(pm.probes, object)
    if epi.sum() < MIN_EPI or (~epi).sum() < 2:
        return None
    lam_o, U_o = np.linalg.eigh(np.diag(W.sum(1)) - W)
    lmax = float(lam_o[-1])
    if lmax <= 0:
        return None
    taus = np.geomspace(1.0 / lmax, 10.0 / lmax, N_TAU)
    obs = relational_feats(lam_o, U_o, epi, taus)
    if obs is None:
        return None
    strength = W.sum(1)

    evals, evecs = load_or_compute_eigs_at_path(
        es.surr_eig_path("full", pat, band, "rest_post", n_surr, swap),
        W, n_surr, swap, es.cell_rng(pat, band, "rest_post", "full"),
        verbose=verbose)
    R = evals.shape[0]
    surr = {c: np.full((R, N), np.nan) for c in REL_COLS}
    for r in range(R):
        fr = relational_feats(evals[r], evecs[r], epi, taus)
        if fr is None:
            continue
        for c in REL_COLS:
            surr[c][r] = fr[c]

    rows = []
    for i in range(N):
        row = {"patient": pat, "band": band, "node": i,
               "is_epi": bool(epi[i]), "probe": probe[i],
               "strength_mean": float(strength[i])}
        for c in REL_COLS:
            o = float(obs[c][i])
            s = surr[c][:, i]
            s = s[np.isfinite(s)]
            row[c] = o
            row[f"z_{c}"] = ((o - s.mean()) / s.std()
                             if s.size >= 5 and s.std() > 0 else 0.0)
        # τ-averaged MS-z summaries (interpretable headline features)
        for f in FEATS:
            zc = [row[f"z_{f}_t{ti}"] for ti in range(N_TAU)]
            row[f"z_{f}_mean"] = float(np.nanmean(zc))
        rows.append(row)
    return pd.DataFrame(rows)


def separation(feat):
    out = []
    summ_cols = [f"z_{f}_mean" for f in FEATS]
    for band in [b for b in es.ALL_BANDS if b in set(feat.band)]:
        fb = feat[feat.band == band]
        # fair strength baseline (oriented within-patient MW AUC)
        s_auc = []
        for _, fp in fb.groupby("patient"):
            a = _auc(fp[fp.is_epi]["strength_mean"].values,
                     fp[~fp.is_epi]["strength_mean"].values)
            if np.isfinite(a):
                s_auc.append(max(a, 1 - a))
        strength_fair = float(np.median(s_auc)) if s_auc else np.nan
        for col in summ_cols + REL_COLS:
            # per-τ cols get raw + MS-z; the z_<f>_mean summaries are MS-z only
            if col in REL_COLS:
                variants = [("raw", col), ("msz", f"z_{col}")]
            else:
                variants = [("msz", col)]
            for kind, c in variants:
                au = []
                for _, fp in fb.groupby("patient"):
                    a = _auc(fp[fp.is_epi][c].values, fp[~fp.is_epi][c].values)
                    if np.isfinite(a):
                        au.append(a)
                au = np.array(au)
                if au.size >= 3:
                    try:
                        _, wp = wilcoxon(au - 0.5, alternative="greater")
                    except Exception:
                        wp = np.nan
                    out.append({
                        "band": band, "feature": col, "kind": kind,
                        "n_pat": int(au.size),
                        "median_auc": float(np.median(au)),
                        "n_above_0.5": int((au > 0.5).sum()),
                        "wilcoxon_p_gt_chance": wp,
                        "strength_fair_auc": strength_fair})
    return pd.DataFrame(out)


def write_readme(sep, runtime):
    L = ["---", "name: epi_relational_features",
         "scope: VI1_interpretable_relational_propagator_features_for_SOZ_marker",
         "era: COHORT_N10 / IMCOH_ABS",
         f"date: {time.strftime('%Y-%m-%d')}",
         "build_script: scripts/01_compute/audit/audit_91_epi_relational_features.py",
         "---", "",
         "# Interpretable relational propagator features for the SOZ marker", "",
         "**Head.** Per-node, strength-orthogonal (matched-strength-z) relational "
         "diffusion features: `f_aff` affinity to the SOZ community, `f_seg` "
         "interface segregation (couples to SOZ more than to healthy, beyond "
         "strength — the V3 signal), `f_part` share of diffusion routed to SOZ. "
         "Screen = within-patient leave-self-out epi-vs-non AUC (epi high), τ-mean "
         "MS-z. Rigorous hide-out validation + calibrated classifier in audit_92.",
         "",
         "## τ-mean MS-z separation (within-patient, epi vs non)", "",
         "| band | feature | median AUC | n>0.5/N | Wilcoxon p | fair strength AUC |",
         "|---|---|---|---|---|---|"]
    summ = sep[sep.feature.str.endswith("_mean")]
    for band in [b for b in es.ALL_BANDS if b in set(summ.band)]:
        for f in FEATS:
            r = summ[(summ.band == band) & (summ.feature == f"z_{f}_mean")]
            if r.empty:
                continue
            r = r.iloc[0]
            L.append(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {f} "
                     f"| {r.median_auc:.3f} | {int(r['n_above_0.5'])}/{int(r.n_pat)} "
                     f"| {r.wilcoxon_p_gt_chance:.4f} | {r.strength_fair_auc:.3f} |")
    L += ["", "## Reading",
          "- MS-z AUC > 0.5 with small Wilcoxon p ⇒ epi nodes are more "
          "diffusion-coupled to the SOZ community than healthy nodes, BEYOND "
          "strength (hubness is fixed in the null). Compare to `fair strength AUC` "
          "(the hubness baseline that failed in audit_86/87).",
          "- Per-τ detail + raw (non-MS-z) AUCs in `relational_separation.csv`; "
          "per-node features in `node_relational_features.csv`.",
          "- Lead bands per audit_89/90: δ, β, low-γ (real community); α tested; "
          "high-γ expected null.",
          f"- Wall-clock: {runtime:.1f}s"]
    (OUT / "README_relational.md").write_text("\n".join(L) + "\n", encoding="utf-8")


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
            res = per_patient_band(pat, band, args.n_surrogates,
                                   args.swap_factor, args.verbose)
            if res is not None:
                frames.append(res)
                print(f"[audit_91] {pat}/{band}: {len(res)} nodes "
                      f"({int(res.is_epi.sum())} epi) ({time.time()-tc:.1f}s)")
    feat = pd.concat(frames, ignore_index=True)
    feat.to_csv(OUT / "node_relational_features.csv", index=False)
    sep = separation(feat)
    sep.to_csv(OUT / "relational_separation.csv", index=False)
    runtime = time.time() - t0
    write_readme(sep, runtime)

    print("\n[audit_91] τ-mean MS-z separation (epi vs non, within-patient):")
    summ = sep[sep.feature.str.endswith("_mean")]
    for band in [b for b in es.ALL_BANDS if b in set(summ.band)]:
        parts = []
        for f in FEATS:
            r = summ[(summ.band == band) & (summ.feature == f"z_{f}_mean")]
            if not r.empty:
                r = r.iloc[0]
                parts.append(f"{f}={r.median_auc:.3f}(p{r.wilcoxon_p_gt_chance:.3f})")
        sf = summ[summ.band == band].strength_fair_auc
        sfv = float(sf.iloc[0]) if not sf.empty else np.nan
        print(f"  {band:10s}: " + " ".join(parts) + f"  [strength {sfv:.3f}]")
    print(f"[audit_91] done in {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
