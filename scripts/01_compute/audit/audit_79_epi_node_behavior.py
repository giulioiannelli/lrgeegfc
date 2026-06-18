#!/usr/bin/env python3
"""Audit 79 — per-node epileptic behaviour under cross-phase trace measures
(Direction-2 DESCRIPTIVE step; gates the classifier).

Lifts the validated per-pair cross-phase statistics to PER-NODE features, each
expressed as excess over the node's own matched-strength surrogate (so it is
strength-independent by construction), and asks per band whether epileptic
nodes separate from non-epileptic ones — and whether any trace feature beats
the strength and probe baselines (the artifact detectors).

Scope report: `.agents/guides/task-persistence-investigation/2026-06-05_epi-node-trace-marker.md`.

This pass computes the three lowest-cost, rank/sign/mass features (F1, F2, F7);
the magnitude features (F3–F5), LOO-resect (F8), and the LOPO classifier
(audit_80) are DEFERRED and gated on this descriptive result.

Critical preamble (per CLAUDE.md rule)
======================================
(1) **Claim.** Epileptic nodes occupy a distinct, band-specific region of the
    per-node trace-feature space (F1 ρ_split_node, F2 coherence_node, F7
    Grassmann mode-load), distinct beyond their matched-strength surrogate.
(2) **Null.** Each per-node feature recomputed on the node's R=200
    matched-strength surrogate ensemble (full-graph cache) → z-score / upper-p.
(3) **Strongest alternative.** The separation is a re-read of (a) node strength
    (per-leaf cophenetic ρ is global-position dominated) or (b) sEEG-shaft
    spatial autocorrelation (epi contacts cluster on probes).
(4) **Mechanical reach.** Surrogate-z fixes per-node strength exactly → removes
    (a) by construction. It does NOT remove (b); the probe baselines
    (n_same_shaft_epi, min_contact_dist_to_epi) are reported alongside and the
    trace features must out-separate them. Whole-shaft masking (audit_80,
    deferred) is the decisive (b) control.
(5) **Falsification.** If no band's trace-feature epi-vs-non-epi AUC exceeds
    0.5 beyond surrogate spread AND beats the strength/probe baselines, the
    marker is a strength/probe re-read → no classifier; report that verdict.

Outputs
-------
``data/audit/epi_marker/``
    node_features_per_band.csv   one row per (patient, band, node)
    epi_behavior_per_band.csv    per (band, feature) epi-vs-non-epi AUC + cohort test
    README.md
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.patient import build_epi_masks, parse_seeg_label
from lrg_eegfc.utils.surrogate import (
    cophenetic_condensed_from_eigs,
    load_or_compute_eigs_at_path,
)

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import (  # type: ignore
    ensure_half_fcs,
    load_phase_fc,
)
import _epi_stratify as es  # type: ignore


OUT = ROOT / "data" / "audit" / "epi_marker"
OUT.mkdir(parents=True, exist_ok=True)

K_LOAD = 10                       # slowest K nontrivial modes for F7 mode-load
TRACE_FEATURES = ["z_rho_split", "z_coherence", "z_grass_modeload"]
BASELINE_FEATURES = ["strength_mean", "n_same_shaft_epi", "min_contact_dist_to_epi"]


# ---------------------------------------------------------------------------
# Vectorized tie-aware Spearman over rows (cophenetic D is heavily tied)
# ---------------------------------------------------------------------------
def _rank_rows(M: np.ndarray) -> np.ndarray:
    try:
        return rankdata(M, axis=1)
    except TypeError:                                  # older scipy
        return np.vstack([rankdata(M[r]) for r in range(M.shape[0])])


def _spearman_rows(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Row-wise Spearman of (R, n) vs (R, n) → (R,), tie-average ranks."""
    ra, rb = _rank_rows(A), _rank_rows(B)
    ra = ra - ra.mean(axis=1, keepdims=True)
    rb = rb - rb.mean(axis=1, keepdims=True)
    num = (ra * rb).sum(axis=1)
    den = np.sqrt((ra ** 2).sum(axis=1) * (rb ** 2).sum(axis=1))
    return np.where(den > 0, num / den, np.nan)


def _coh_rows(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Row-wise sign-agreement fraction of (R, n) vs (R, n) → (R,)."""
    return (np.sign(A) == np.sign(B)).mean(axis=1)


def _z(obs: float, surr: np.ndarray) -> tuple[float, float, int]:
    s = surr[np.isfinite(surr)]
    if s.size < 5:
        return float("nan"), float("nan"), int(s.size)
    sd = s.std(ddof=1)
    z = (obs - s.mean()) / sd if sd > 0 else float("nan")
    p = (1 + int((s >= obs).sum())) / (s.size + 1)      # Phipson–Smyth upper
    return float(z), float(p), int(s.size)


# ---------------------------------------------------------------------------
# Per (patient, band) feature rows
# ---------------------------------------------------------------------------
def per_patient_band(pat: str, band: str, n_surr: int, swap_factor: int,
                     verbose: bool) -> list[dict]:
    try:
        Ws = {ph: load_phase_fc(pat, ph, band) for ph in es.PHASES_4}
    except Exception as e:
        if verbose:
            print(f"[audit_79] SKIP {pat}/{band}: {e}")
        return []
    N = Ws["rest_pre_A"].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        return []
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        if verbose:
            print(f"[audit_79] SKIP {pat}/{band}: channels {len(pm.channels)} != N {N}")
        return []
    epi = np.asarray(pm.epi_mask, dtype=bool)
    probes = np.asarray(pm.probes, dtype=object)
    contacts = np.array([parse_seeg_label(c)[1] for c in pm.channels], dtype=object)

    iu_i, iu_j = es.upper_tri_indices(N)
    # per-node condensed indices of incident pairs
    node_pairs: list[list[int]] = [[] for _ in range(N)]
    for c in range(iu_i.size):
        node_pairs[iu_i[c]].append(c)
        node_pairs[iu_j[c]].append(c)
    node_pairs = [np.asarray(p, dtype=int) for p in node_pairs]

    # observed eigs (one eigh/phase) → cophenetic + post eigvecs
    obs_eig = {ph: np.linalg.eigh(np.diag(Ws[ph].sum(1)) - Ws[ph])
               for ph in es.PHASES_4}
    D_obs = {ph: cophenetic_condensed_from_eigs(*obs_eig[ph])
             for ph in es.PHASES_4}
    dT_o = D_obs["task_test"] - D_obs["rest_pre_A"]
    dR_o = D_obs["rest_post"] - D_obs["rest_pre_B"]
    Vpost_obs = obs_eig["rest_post"][1]
    F7_obs = (Vpost_obs[:, 1:K_LOAD + 1] ** 2).sum(axis=1)        # (N,)

    # surrogate eigs (cached) → per-r cophenetic stacks + post mode-load
    surr_coph: dict[str, np.ndarray] = {}
    Vpost_surr = None
    for ph in es.PHASES_4:
        evals, evecs = load_or_compute_eigs_at_path(
            es.surr_eig_path("full", pat, band, ph, n_surr, swap_factor),
            Ws[ph], n_surr, swap_factor,
            es.cell_rng(pat, band, ph, "full"), verbose=verbose)
        R = evals.shape[0]
        S = np.full((R, iu_i.size), np.nan)
        for r in range(R):
            if np.isfinite(evals[r]).all():
                S[r] = cophenetic_condensed_from_eigs(evals[r], evecs[r])
        surr_coph[ph] = S
        if ph == "rest_post":
            F7s = np.full((R, N), np.nan)
            for r in range(R):
                if np.isfinite(evecs[r]).all():
                    F7s[r] = (evecs[r][:, 1:K_LOAD + 1] ** 2).sum(axis=1)
            Vpost_surr = F7s
    dT_s = surr_coph["task_test"] - surr_coph["rest_pre_A"]        # (R, n_pairs)
    dR_s = surr_coph["rest_post"] - surr_coph["rest_pre_B"]

    strengths = {ph: Ws[ph].sum(axis=1) for ph in es.PHASES_4}
    smat = np.vstack([strengths[ph] for ph in es.PHASES_4])         # (4, N)
    s_mean = smat.mean(axis=0)
    s_cv = np.where(s_mean > 0, smat.std(axis=0) / s_mean, np.nan)

    rows = []
    for i in range(N):
        ci = node_pairs[i]
        valid = np.isfinite(dT_s[:, ci]).all(axis=1) & \
            np.isfinite(dR_s[:, ci]).all(axis=1)
        # F1 rho_split_node
        f1_obs = _spearman_rows(dT_o[ci][None, :], dR_o[ci][None, :])[0]
        f1_surr = _spearman_rows(dT_s[valid][:, ci], dR_s[valid][:, ci]) \
            if valid.any() else np.array([])
        z1, p1, n1 = _z(f1_obs, f1_surr)
        # F2 coherence_node
        f2_obs = _coh_rows(dT_o[ci][None, :], dR_o[ci][None, :])[0]
        f2_surr = _coh_rows(dT_s[valid][:, ci], dR_s[valid][:, ci]) \
            if valid.any() else np.array([])
        z2, p2, _ = _z(f2_obs, f2_surr)
        # F7 grass_modeload
        f7_obs = float(F7_obs[i])
        f7_surr = Vpost_surr[np.isfinite(Vpost_surr[:, i]), i] \
            if Vpost_surr is not None else np.array([])
        z7, p7, _ = _z(f7_obs, f7_surr)

        # probe baselines (descriptive: full labels; leakage matters only in clf)
        same = (probes == probes[i]); same[i] = False
        n_shaft_epi = int((same & epi).sum())
        probe_size = int(same.sum() + 1)
        ci_num = contacts[i]
        if ci_num is not None:
            dists = [abs(int(ci_num) - int(contacts[j]))
                     for j in np.where(same & epi)[0]
                     if contacts[j] is not None]
            min_dist = float(min(dists)) if dists else np.nan
            cidx_norm = (float(ci_num) /
                         max([int(contacts[j]) for j in np.where(same)[0]
                              if contacts[j] is not None] + [int(ci_num)]))
        else:
            min_dist = np.nan; cidx_norm = np.nan

        rows.append({
            "patient": pat, "band": band, "node": i,
            "channel": pm.channels[i], "probe": str(probes[i]),
            "contact": (int(ci_num) if ci_num is not None else -1),
            "is_epi": bool(epi[i]),
            "rho_split_node": f1_obs, "z_rho_split": z1, "p_rho_split": p1,
            "coherence_node": f2_obs, "z_coherence": z2,
            "grass_modeload": f7_obs, "z_grass_modeload": z7,
            "n_surr_finite": n1,
            "strength_mean": float(s_mean[i]), "strength_cv": float(s_cv[i]),
            "probe_has_other_epi": int(n_shaft_epi > 0),
            "n_same_shaft_epi": n_shaft_epi,
            "min_contact_dist_to_epi": min_dist,
            "probe_size": probe_size, "contact_index_norm": cidx_norm,
        })
    return rows


# ---------------------------------------------------------------------------
# epi-vs-non-epi separation per band per feature
# ---------------------------------------------------------------------------
def _auc(epi_vals: np.ndarray, non_vals: np.ndarray) -> float:
    """P(feature_epi > feature_nonepi) via Mann-Whitney U."""
    epi_vals = epi_vals[np.isfinite(epi_vals)]
    non_vals = non_vals[np.isfinite(non_vals)]
    if epi_vals.size < 2 or non_vals.size < 2:
        return float("nan")
    u = mannwhitneyu(epi_vals, non_vals, alternative="two-sided").statistic
    return float(u / (epi_vals.size * non_vals.size))


def behaviour_table(feat: pd.DataFrame) -> pd.DataFrame:
    out = []
    features = TRACE_FEATURES + BASELINE_FEATURES
    for band in [b for b in es.ALL_BANDS if b in set(feat.band)]:
        fb = feat[feat.band == band]
        for f in features:
            aucs, pats = [], []
            for pat, fp in fb.groupby("patient"):
                a = _auc(fp[fp.is_epi][f].values, fp[~fp.is_epi][f].values)
                if np.isfinite(a):
                    aucs.append(a); pats.append(pat)
            aucs = np.array(aucs)
            if aucs.size >= 3:
                try:
                    _, wp = wilcoxon(aucs - 0.5, alternative="two-sided")
                except Exception:
                    wp = float("nan")
                med = float(np.median(aucs))
                n_gt = int((aucs > 0.5).sum())
            else:
                wp = float("nan"); med = float("nan"); n_gt = int((aucs > 0.5).sum())
            out.append({
                "band": band, "feature": f,
                "kind": "trace" if f in TRACE_FEATURES else "baseline",
                "n_patients": int(aucs.size),
                "median_auc": med,
                "abs_median_auc_dev": (abs(med - 0.5) if np.isfinite(med) else np.nan),
                "n_patients_auc_gt_0.5": n_gt,
                "cohort_wilcoxon_p": wp,
            })
    return pd.DataFrame(out)


def write_readme(feat: pd.DataFrame, beh: pd.DataFrame, runtime_s: float) -> None:
    L: list[str] = []
    a = L.append
    a("---")
    a("name: epi_marker_node_behaviour")
    a("scope: direction2_descriptive_per_node_epi_vs_nonepi_trace_behaviour")
    a("era: COHORT_N10 / IMCOH_ABS")
    a(f"date: {time.strftime('%Y-%m-%d')}")
    a("status: descriptive_gate")
    a("build_script: scripts/01_compute/audit/audit_79_epi_node_behavior.py")
    a("scope_report: .agents/guides/task-persistence-investigation/2026-06-05_epi-node-trace-marker.md")
    a("---")
    a("")
    a("# Per-node epi-vs-non-epi trace behaviour (descriptive gate)")
    a("")
    a("**Head.** Per-node surrogate-z trace features (F1 ρ_split_node, F2 "
      "coherence_node, F7 Grassmann mode-load) vs node-strength and probe "
      "baselines, per band, as epi-vs-non-epi separation AUC "
      "(P(feature_epi > feature_nonepi), Mann-Whitney, per patient; cohort "
      "Wilcoxon of AUC−0.5). **A trace feature is a marker candidate only "
      "where its |AUC−0.5| exceeds the strength and probe baselines.** This "
      "table GATES whether the audit_80 classifier is worth building.")
    a("")
    a("## Separation AUC — band × feature (trace vs baseline)")
    a("")
    a("| band | feature | kind | n_pat | median AUC | |AUC−.5| | n>0.5 | Wilcoxon p |")
    a("|---|---|---|---|---|---|---|---|")
    for band in [b for b in es.ALL_BANDS if b in set(beh.band)]:
        for _, r in beh[beh.band == band].iterrows():
            med = r["median_auc"]
            a(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {r['feature']} "
              f"| {r['kind']} | {r['n_patients']} "
              f"| {med:.3f} | {r['abs_median_auc_dev']:.3f} "
              f"| {r['n_patients_auc_gt_0.5']} | {r['cohort_wilcoxon_p']:.4f} |")
    a("")
    a("## Reading")
    a("- AUC = P(epi node has higher feature than non-epi node), per patient, "
      "median over patients. 0.5 = no separation.")
    a("- **Gate:** a band where a trace feature's |AUC−0.5| > the best baseline "
      "|AUC−0.5| is a marker candidate. If every trace feature is dominated by "
      "strength_mean or the probe baselines, the marker is an artifact — the "
      "classifier (audit_80) becomes a negative control only.")
    a("- Pat_15 (0 epi) excluded from AUC (no positives).")
    a("")
    a("## Provenance")
    a("- Features: F1 rho_split_node, F2 coherence_node (per-node row "
      "restriction of audit_63/76 pair stats), F7 grass_modeload "
      f"(slowest {K_LOAD} nontrivial rest_post modes), each surrogate-z'd "
      "against the full-graph R=200 matched-strength ensemble (seed 20260511).")
    a("- Baselines: strength_mean (per-node, 4-phase mean), n_same_shaft_epi, "
      "min_contact_dist_to_epi (probe spatial autocorrelation detectors).")
    a("- DEFERRED (gated): F3–F5 magnitude features, F8 LOO-resect, and the "
      "audit_80 LOPO classifier + whole-shaft masking + candidate discovery.")
    a(f"- Wall-clock: {runtime_s:.1f} s")
    a("")
    a("## Files")
    a("- `node_features_per_band.csv` — one row per (patient, band, node)")
    a("- `epi_behavior_per_band.csv` — per (band, feature) AUC + cohort test")
    (OUT / "README.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--n-surrogates", type=int, default=es.N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=es.SWAP_FACTOR)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    print("[audit_79] pre-flight: ensure rsPre half FCs cached")
    for pat in args.patients:
        ensure_half_fcs(pat, args.bands)

    t0 = time.time()
    all_rows: list[dict] = []
    for band in args.bands:
        for pat in args.patients:
            tc = time.time()
            rows = per_patient_band(pat, band, args.n_surrogates,
                                    args.swap_factor, args.verbose)
            all_rows.extend(rows)
            if rows:
                ne = sum(r["is_epi"] for r in rows)
                print(f"[audit_79] {pat}/{band}: {len(rows)} nodes ({ne} epi) "
                      f"({time.time()-tc:.1f}s)")
    runtime = time.time() - t0

    feat = pd.DataFrame(all_rows)
    feat.to_csv(OUT / "node_features_per_band.csv", index=False)
    beh = behaviour_table(feat)
    beh.to_csv(OUT / "epi_behavior_per_band.csv", index=False)
    write_readme(feat, beh, runtime)
    print(f"[audit_79] {len(feat)} node-rows, {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
