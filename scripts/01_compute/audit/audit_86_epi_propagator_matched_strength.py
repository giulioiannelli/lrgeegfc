#!/usr/bin/env python3
"""Audit 86 — matched-strength CERTIFICATION of the propagator-subspace epi marker.

The audit_85 within-patient recovery showed the LRG propagator ρ(τ) subspace
recovering held-out epi nodes (β full AUC 0.715) and apparently beating a
strength baseline (0.622). TWO things must be checked before that is called real,
because the project's whole history says "it's hubness" is the default:

  (1) FAIR strength baseline. audit_85 scored strength by template-DISTANCE,
      which handicaps a single monotone feature (it flags nodes NEAR the epi
      mean, penalising nodes more extreme than epi). The fair strength baseline
      is the within-patient Mann-Whitney AUC of strength (epi rank high). This
      script reports it.
  (2) MATCHED-STRENGTH surrogate-z. Recompute every propagator feature on the
      node's R=200 matched-strength surrogate ensemble (the cached rest_post
      eigs, seed 20260511) and z-score the observation against it. The surrogate
      fixes each node's strength exactly, so any epi separation left in the
      surrogate-z feature is STRENGTH-ORTHOGONAL by construction. If the
      surrogate-z propagator features still separate / recover epi, the marker is
      pattern, not hubness. If they collapse to chance, the audit_85 "beat" was
      a framing artifact (template-distance) on top of hubness.

Critical preamble
=================
(1) Claim: the propagator subspace carries strength-ORTHOGONAL epi information.
(2) Null: matched-strength 4-cycle±δ surrogate (R=200, cached) → per-feature z.
(3) Strongest alternative: it is all node strength (return prob, communication
    distance are strength-correlated).
(4) Reach: surrogate-z removes the strength-magnitude channel exactly; the fair
    Mann-Whitney strength baseline removes the template-distance handicap. Both
    are reported so neither can hide the truth.
(5) Falsification: if no surrogate-z propagator feature separates epi above
    chance (cohort Wilcoxon) AND the surrogate-z subspace recovery ≈ 0.5, the
    marker is hubness — report that and retract the "beats strength" framing.

Outputs (``data/audit/epi_marker/``)
------------------------------------
    propagator_ms_features.csv     per (patient, band, node): raw + surrogate-z propagator feats
    propagator_ms_certify.csv      per (band, feature): raw AUC, MSz AUC, fair-strength AUC, recovery
    README_ms_certify.md
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
from audit_63_split_baseline_surrogate import ensure_half_fcs, load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import (  # type: ignore
    MIN_EPI, N_TAU, PROP_COLS, _loo_recovery, _standardize)

OUT = ROOT / "data" / "audit" / "epi_marker"
OUT.mkdir(parents=True, exist_ok=True)


def diffusion_from_eigs(lam: np.ndarray, U: np.ndarray) -> dict | None:
    """Propagator ρ(τ) per-node features from a Laplacian eigendecomposition
    (works identically for the observed graph and each surrogate)."""
    lam = np.clip(np.asarray(lam, float), 0.0, None)
    if lam.size == 0 or not np.isfinite(lam).all() or not np.isfinite(U).all():
        return None
    lmax = float(lam[-1])
    if lmax <= 0:
        return None
    taus = np.geomspace(1.0 / lmax, 10.0 / lmax, N_TAU)
    f: dict[str, np.ndarray] = {}
    for ti, t in enumerate(taus):
        w = np.exp(-t * lam); Z = w.sum()
        f[f"rho_ii_t{ti}"] = (U ** 2 * w).sum(1) / Z
        K = np.clip((U * w) @ U.T, 1e-300, None)
        P = K / K.sum(1, keepdims=True)
        f[f"Hdiff_t{ti}"] = -(P * np.log(P)).sum(1)
    w0 = np.exp(-(1.0 / lmax) * lam); Z0 = w0.sum()
    rho0 = np.clip(((U * w0) @ U.T) / Z0, 1e-300, None)
    Tr = 1.0 / rho0; np.fill_diagonal(Tr, np.nan)
    f["Tcomm_mean"] = np.nanmean(Tr, 1); f["Tcomm_min"] = np.nanmin(Tr, 1)
    return f


def _auc(a: np.ndarray, b: np.ndarray) -> float:
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return float("nan")
    return float(mannwhitneyu(a, b, alternative="two-sided").statistic
                 / (a.size * b.size))


def per_patient_band(pat: str, band: str, n_surr: int, swap: int,
                     verbose: bool) -> tuple[pd.DataFrame, dict] | None:
    try:
        W = load_phase_fc(pat, "rest_post", band)
    except Exception as e:
        if verbose:
            print(f"[audit_86] SKIP {pat}/{band}: {e}")
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    if epi.sum() < MIN_EPI:
        return None
    # observed propagator features
    lam_o, U_o = np.linalg.eigh(np.diag(W.sum(1)) - W)
    obs = diffusion_from_eigs(lam_o, U_o)
    if obs is None:
        return None
    strength = W.sum(1)
    # surrogate ensemble (cached rest_post matched-strength eigs)
    evals, evecs = load_or_compute_eigs_at_path(
        es.surr_eig_path("full", pat, band, "rest_post", n_surr, swap),
        W, n_surr, swap, es.cell_rng(pat, band, "rest_post", "full"),
        verbose=verbose)
    R = evals.shape[0]
    surr = {c: np.full((R, N), np.nan) for c in PROP_COLS}
    for r in range(R):
        fr = diffusion_from_eigs(evals[r], evecs[r])
        if fr is None:
            continue
        for c in PROP_COLS:
            surr[c][r] = fr[c]
    rows = []
    for i in range(N):
        row = {"patient": pat, "band": band, "node": i, "is_epi": bool(epi[i]),
               "strength_mean": float(strength[i])}
        for c in PROP_COLS:
            o = float(obs[c][i])
            s = surr[c][:, i]; s = s[np.isfinite(s)]
            row[c] = o
            row[f"z_{c}"] = ((o - s.mean()) / s.std()
                             if s.size >= 5 and s.std() > 0 else 0.0)
        rows.append(row)
    df = pd.DataFrame(rows)
    return df, {}


def certify(feat: pd.DataFrame) -> pd.DataFrame:
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(feat.band)]:
        fb = feat[feat.band == band]
        # --- fair strength baseline (within-patient MW AUC, oriented) ---
        s_auc = []
        for _, fp in fb.groupby("patient"):
            a = _auc(fp[fp.is_epi]["strength_mean"].values,
                     fp[~fp.is_epi]["strength_mean"].values)
            if np.isfinite(a):
                s_auc.append(max(a, 1 - a))            # oriented = best direction
        strength_fair = float(np.median(s_auc)) if s_auc else np.nan
        # --- per-feature raw + MSz within-patient AUC ---
        for c in PROP_COLS:
            for kind, col in (("raw", c), ("msz", f"z_{c}")):
                au = []
                for _, fp in fb.groupby("patient"):
                    a = _auc(fp[fp.is_epi][col].values,
                             fp[~fp.is_epi][col].values)
                    if np.isfinite(a):
                        au.append(a)
                au = np.array(au)
                if au.size >= 3:
                    try:
                        _, wp = wilcoxon(au - 0.5, alternative="two-sided")
                    except Exception:
                        wp = np.nan
                    dev = np.array([max(a, 1 - a) for a in au])
                    out.append({
                        "band": band, "feature": c, "kind": kind,
                        "n_pat": int(au.size),
                        "median_auc": float(np.median(au)),
                        "median_auc_oriented": float(np.median(dev)),
                        "wilcoxon_p": wp, "strength_fair_auc": strength_fair})
        # --- subspace template recovery: raw vs MSz vs strength-monotone ---
        rec = {"diffusion_raw": [], "diffusion_msz": [], "strength_mono": []}
        for _, fp in fb.groupby("patient"):
            fp = fp.reset_index(drop=True)
            y = fp["is_epi"].to_numpy(bool)
            if y.sum() < MIN_EPI or (~y).sum() < 2:
                continue
            epi_idx = np.where(y)[0]; neg = np.where(~y)[0]
            Zraw = _standardize(fp[PROP_COLS].to_numpy(float))
            Zz = _standardize(fp[[f"z_{c}" for c in PROP_COLS]].to_numpy(float))
            rec["diffusion_raw"].append(_loo_recovery(Zraw, epi_idx, neg))
            rec["diffusion_msz"].append(_loo_recovery(Zz, epi_idx, neg))
            # strength as ORIENTED monotone (fair), same LOO metric
            s = fp["strength_mean"].to_numpy(float)
            d = np.sign(s[epi_idx].mean() - s.mean()) or 1.0
            sc = d * s
            wins = 0.0; tot = 0
            for e in epi_idx:
                wins += float((sc[e] > sc[neg]).sum()) + \
                    0.5 * float((sc[e] == sc[neg]).sum())
                tot += neg.size
            rec["strength_mono"].append(wins / tot if tot else np.nan)
        for sn, vals in rec.items():
            v = np.array([x for x in vals if np.isfinite(x)])
            wp = np.nan
            if v.size >= 3:
                try:
                    _, wp = wilcoxon(v - 0.5, alternative="greater")
                except Exception:
                    pass
            out.append({"band": band, "feature": f"__recovery_{sn}",
                        "kind": "recovery", "n_pat": int(v.size),
                        "median_auc": float(np.median(v)) if v.size else np.nan,
                        "median_auc_oriented": float(np.median(v)) if v.size else np.nan,
                        "wilcoxon_p": wp, "strength_fair_auc": strength_fair})
    return pd.DataFrame(out)


def write_readme(cert: pd.DataFrame, runtime: float) -> None:
    L = []; a = L.append
    a("---")
    a("name: epi_propagator_matched_strength_certification")
    a("scope: direction2c_certify_propagator_subspace_is_strength_orthogonal")
    a("era: COHORT_N10 / IMCOH_ABS")
    a(f"date: {time.strftime('%Y-%m-%d')}")
    a("status: certification")
    a("build_script: scripts/01_compute/audit/audit_86_epi_propagator_matched_strength.py")
    a("---")
    a("")
    a("# Matched-strength certification of the propagator-subspace epi marker")
    a("")
    a("**Head.** Is the audit_85 within-patient propagator recovery STRENGTH-"
      "ORTHOGONAL or hubness? Two controls: a FAIR strength baseline "
      "(within-patient Mann-Whitney AUC, not template-distance), and "
      "matched-strength surrogate-z of every propagator feature (R=200 cached "
      "rest_post ensemble) — the surrogate fixes each node's strength exactly, so "
      "surrogate-z separation is strength-orthogonal by construction.")
    a("")
    rec = cert[cert.kind == "recovery"]
    if not rec.empty:
        def recval(band, sn):
            r = rec[(rec.band == band) & (rec.feature == f"__recovery_{sn}")]
            return float(r.median_auc.iloc[0]) if not r.empty else np.nan
        a("## Bottom line — recovery AUC, fair baseline vs strength-orthogonal")
        a("")
        a("| band | strength (fair mono) | diffusion raw | diffusion MS-z "
          "(strength-removed) | verdict |")
        a("|---|---|---|---|---|")
        for band in [b for b in es.ALL_BANDS if b in set(rec.band)]:
            sm = recval(band, "strength_mono")
            dr = recval(band, "diffusion_raw")
            dz = recval(band, "diffusion_msz")
            # certified iff MS-z recovery clears chance AND >= fair strength
            cert_ok = (np.isfinite(dz) and dz > 0.55
                       and np.isfinite(sm) and dz >= sm - 0.02)
            v = ("CERTIFIED strength-orthogonal" if cert_ok
                 else ("partial (MS-z>chance, ≤strength)"
                       if np.isfinite(dz) and dz > 0.55 else
                       "NOT certified (collapses to chance → hubness)"))
            a(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {sm:.3f} | {dr:.3f} "
              f"| {dz:.3f} | {v} |")
        a("")
        a("**Reading.** `strength (fair mono)` is the honest baseline (the "
          "audit_85 template-distance strength was unfairly low). `diffusion "
          "MS-z` is recovery using ONLY the strength-orthogonal part of each "
          "propagator feature. If MS-z recovery stays well above chance and ≳ "
          "the fair strength baseline, the marker is real pattern. If MS-z "
          "collapses toward 0.5, the audit_85 advantage was template-distance "
          "framing on top of hubness.")
        a("")
    a("## Per-feature strength-orthogonal separation (MS-z within-patient AUC)")
    a("")
    a("| band | feature | raw AUC (oriented) | MS-z AUC (oriented) | MS-z "
      "Wilcoxon p | fair strength AUC |")
    a("|---|---|---|---|---|---|")
    perf = cert[cert.kind.isin(["raw", "msz"])]
    for band in [b for b in es.ALL_BANDS if b in set(perf.band)]:
        for c in PROP_COLS:
            raw = perf[(perf.band == band) & (perf.feature == c)
                       & (perf.kind == "raw")]
            msz = perf[(perf.band == band) & (perf.feature == c)
                       & (perf.kind == "msz")]
            if raw.empty or msz.empty:
                continue
            a(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {c} "
              f"| {raw.median_auc_oriented.iloc[0]:.3f} "
              f"| {msz.median_auc_oriented.iloc[0]:.3f} "
              f"| {msz.wilcoxon_p.iloc[0]:.4f} "
              f"| {raw.strength_fair_auc.iloc[0]:.3f} |")
    a("")
    a(f"- Wall-clock: {runtime:.1f}s")
    (OUT / "README_ms_certify.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--n-surrogates", type=int, default=es.N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=es.SWAP_FACTOR)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    print("[audit_86] pre-flight: ensure rsPre half FCs cached")
    for pat in args.patients:
        ensure_half_fcs(pat, args.bands)

    t0 = time.time()
    frames = []
    for band in args.bands:
        for pat in args.patients:
            tc = time.time()
            res = per_patient_band(pat, band, args.n_surrogates,
                                   args.swap_factor, args.verbose)
            if res is not None:
                frames.append(res[0])
                ne = int(res[0].is_epi.sum())
                print(f"[audit_86] {pat}/{band}: {len(res[0])} nodes "
                      f"({ne} epi) ({time.time()-tc:.1f}s)")
    feat = pd.concat(frames, ignore_index=True)
    feat.to_csv(OUT / "propagator_ms_features.csv", index=False)
    cert = certify(feat)
    cert.to_csv(OUT / "propagator_ms_certify.csv", index=False)
    runtime = time.time() - t0
    write_readme(cert, runtime)

    # console verdict
    rec = cert[cert.kind == "recovery"]
    print("\n[audit_86] CERTIFICATION (recovery AUC):")
    for band in [b for b in es.ALL_BANDS if b in set(rec.band)]:
        r = {row.feature.replace("__recovery_", ""): row.median_auc
             for row in rec[rec.band == band].itertuples()}
        print(f"  {band}: strength_mono={r.get('strength_mono', np.nan):.3f} "
              f"diffusion_raw={r.get('diffusion_raw', np.nan):.3f} "
              f"diffusion_MSz={r.get('diffusion_msz', np.nan):.3f}")
    print(f"[audit_86] done in {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
