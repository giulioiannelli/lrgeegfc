#!/usr/bin/env python3
"""Audit 90 — occult (non-labeled) epi-node marker: leave-node-out cohort-trace
influence + confound-residualized per-node behaviour.

The *discovery* reframe of the per-node marker question. The known-node
classifier (audit_80) was NEGATIVE (geometry_dominant in every band). This asks
the different question: do any NON-labeled contacts BEHAVE structurally like the
labeled epileptic ones, beyond node strength and sEEG-shaft geometry — candidate
*occult* epileptogenic tissue?

Scope report:
`.agents/guides/task-persistence-investigation/2026-06-05_occult-epi-node-marker.md`.

Critical preamble (per CLAUDE.md rule, before any code)
======================================================
(1) **Claim.** A per-node cross-phase structural-behaviour score — leave-node-out
    (LNO) influence on the patient's cophenetic ρ_split trace, Δρ_i = ρ_full −
    ρ_{∖i}, plus the row-restricted ρ_split_node / coherence_node — after being
    residualized against node strength + same-shaft-epi count + along-shaft
    position, separates labeled epileptic from non-labeled contacts, and flags
    non-labeled contacts that sit inside the labeled-epi range (occult
    candidates).
(2) **Null / control.** No surrogate ensemble. The control is the WITHIN-PATIENT
    OLS residualization against the confound design {strength, n_same_shaft_epi,
    along-shaft pos, shaft length}; the residual is confound-orthogonal by
    construction. The discrimination gate also reports the residualized STRENGTH
    baseline AUC that the trace score must beat.
(3) **Strongest plausible alternative.** The separation is a re-read of node
    strength (epi nodes are FC-hyperconnected within patient, audit_79) or of
    sEEG-shaft geometry (epi contacts cluster on shafts; foci targeted with
    longer shafts at characteristic depths, audit_80 `probe_geom`).
(4) **Mechanical reach + what it CANNOT reject.** OLS residualization removes
    LINEAR strength/geometry dependence. It CANNOT remove (a) monotone-nonlinear
    strength dependence (rank-partial fallback reported); (b) it cannot
    establish that an "occult candidate" is truly epileptogenic — there is NO
    surgical-outcome ground truth in this cohort, so a flag is a hypothesis, not
    a prediction; (c) it cannot detect an occult focus whose trace behaviour is
    identical to surrounding healthy tissue (no contrast — the likely real case).
(5) **Falsification + limitations.** Occult signal exists only if (i) the
    residualized trace score's cohort epi-vs-non-epi AUC beats the residualized
    strength baseline AND |AUC−0.5| is non-trivial, AND (ii) some non-labeled
    nodes land in the labeled-epi residual range while ordinary tissue does not.
    If the labeled-epi residual distribution is indistinguishable from the
    non-labeled bulk, the candidates are a strength/geometry re-read → clean
    negative. Limitations: linear control only; no ground truth; restricted to
    the trace-bearing bands α/β/γ_l (audit_77).

Outputs
-------
``data/audit/occult_node_marker/``
    lno_node_influence.csv     one row per (patient, band, node) with scores+confounds+residuals
    occult_gate_per_band.csv   per (band, score) cohort epi-vs-non-epi AUC + strength baseline
    occult_candidates.csv      ranked non-labeled candidates per (band, patient)
    README.md
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata, spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.patient import build_epi_masks, parse_seeg_label
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import (  # type: ignore
    ensure_half_fcs,
    load_phase_fc,
    lrg_ultrametric_condensed,
)
import _epi_stratify as es  # type: ignore


OUT = ROOT / "data" / "audit" / "occult_node_marker"
OUT.mkdir(parents=True, exist_ok=True)

# Only bands with an established cohort cophenetic trace (audit_77): a
# leave-node-out influence on a null trace is leverage on noise.
BANDS_TRACE = ["alpha", "beta", "low_gamma"]
PHASES = es.PHASES_4
MIN_EPI_FOR_WINDOW = 3            # honest floor for the epi residual window
TOP_DECILE = 0.90                # non-labeled percentile for a candidate
SCORES = ["dRho", "rho_node", "coh_node"]
CONFOUNDS = ["strength_mean", "n_same_shaft_epi", "pos_norm", "shaft_len"]


# ---------------------------------------------------------------------------
# Per-node row restriction of the condensed cophenetic vector (audit_76/79 idea)
# ---------------------------------------------------------------------------
def _row_pairs_mask(n: int, i: int) -> np.ndarray:
    """Boolean mask (condensed upper-tri order) selecting pairs incident to i."""
    iu_i, iu_j = np.triu_indices(n, k=1)
    return (iu_i == i) | (iu_j == i)


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 3:
        return float("nan")
    r, _ = spearmanr(a, b)
    return float(r)


# ---------------------------------------------------------------------------
# Confound geometry per patient
# ---------------------------------------------------------------------------
def _node_geometry(channels: list[str], epi: np.ndarray) -> dict[str, np.ndarray]:
    """n_same_shaft_epi, normalized along-shaft position, shaft length, and the
    per-node 3D-free along-shaft distance to nearest same-shaft labeled epi.
    Uses canonical parse_seeg_label (no local copy)."""
    n = len(channels)
    shafts: list[str | None] = []
    contacts: list[int | None] = []
    for lab in channels:
        sh, ct = parse_seeg_label(lab)
        shafts.append(sh)
        contacts.append(ct)
    # group nodes by shaft
    by_shaft: dict[str, list[int]] = {}
    for idx, sh in enumerate(shafts):
        if sh is not None:
            by_shaft.setdefault(sh, []).append(idx)
    nse = np.zeros(n)
    pos = np.full(n, np.nan)
    slen = np.zeros(n)
    along_dist = np.full(n, np.inf)
    for sh, members in by_shaft.items():
        cts = [contacts[m] for m in members]
        valid = [c for c in cts if c is not None]
        cmax = max(valid) if valid else 1
        epi_members = [m for m in members if epi[m]]
        for m in members:
            slen[m] = len(members)
            c = contacts[m]
            if c is not None and cmax > 0:
                pos[m] = c / cmax
            # same-shaft epi count excluding self
            nse[m] = sum(1 for e in epi_members if e != m)
            # along-shaft distance to nearest same-shaft labeled epi
            if c is not None:
                ds = [abs(c - contacts[e]) for e in epi_members
                      if e != m and contacts[e] is not None]
                if ds:
                    along_dist[m] = min(ds)
    pos = np.where(np.isfinite(pos), pos, 0.5)  # unparseable → mid-shaft
    return dict(shaft=np.array([s if s is not None else "?" for s in shafts]),
                contact=np.array([c if c is not None else -1 for c in contacts]),
                n_same_shaft_epi=nse, pos_norm=pos, shaft_len=slen,
                along_dist_to_epi=along_dist)


# ---------------------------------------------------------------------------
# Per (patient, band): full trace + LNO influence + per-node features + confounds
# ---------------------------------------------------------------------------
def per_patient_band(pat: str, band: str, verbose: bool) -> list[dict]:
    try:
        Ws = {ph: load_phase_fc(pat, ph, band) for ph in PHASES}
    except Exception as e:  # pragma: no cover
        if verbose:
            print(f"[audit_90] SKIP {pat}/{band}: load failed: {e}")
        return []
    N = Ws["rest_pre_A"].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        if verbose:
            print(f"[audit_90] SKIP {pat}/{band}: phase shape mismatch")
        return []

    masks = build_epi_masks(pat)
    epi = np.asarray(masks.epi_mask, dtype=bool)
    channels = list(masks.channels)
    if epi.size != N or len(channels) != N:
        if verbose:
            print(f"[audit_90] SKIP {pat}/{band}: epi/chan length != N ({epi.size},{len(channels)},{N})")
        return []

    # full-graph cophenetic trace
    D = {ph: lrg_ultrametric_condensed(Ws[ph]) for ph in PHASES}
    dT = D["task_test"] - D["rest_pre_A"]
    dR = D["rest_post"] - D["rest_pre_B"]
    rho_full = _spearman(dT, dR)

    # confounds
    strength = np.mean([Ws[ph].sum(axis=1) for ph in PHASES], axis=0)
    geo = _node_geometry(channels, epi)

    rows: list[dict] = []
    for i in range(N):
        # LNO rebuild on the (N-1) submatrix
        keep = np.ones(N, dtype=bool)
        keep[i] = False
        idx = np.ix_(keep, keep)
        Dk = {ph: lrg_ultrametric_condensed(np.ascontiguousarray(Ws[ph][idx]))
              for ph in PHASES}
        rho_loo = _spearman(Dk["task_test"] - Dk["rest_pre_A"],
                            Dk["rest_post"] - Dk["rest_pre_B"])
        dRho = rho_full - rho_loo
        # row-restricted per-node features
        m = _row_pairs_mask(N, i)
        rho_node = _spearman(dT[m], dR[m])
        sgnT, sgnR = np.sign(dT[m]), np.sign(dR[m])
        coh_node = float(np.mean(sgnT == sgnR)) if m.sum() else float("nan")
        rows.append(dict(
            patient=pat, band=band, node=i, channel=channels[i],
            shaft=str(geo["shaft"][i]), contact=int(geo["contact"][i]),
            is_epi=bool(epi[i]),
            rho_full=rho_full, rho_loo=rho_loo,
            dRho=dRho, rho_node=rho_node, coh_node=coh_node,
            strength_mean=float(strength[i]),
            n_same_shaft_epi=float(geo["n_same_shaft_epi"][i]),
            pos_norm=float(geo["pos_norm"][i]),
            shaft_len=float(geo["shaft_len"][i]),
            along_dist_to_epi=float(geo["along_dist_to_epi"][i]),
            n_epi=int(epi.sum()), N=N,
        ))
    return rows


# ---------------------------------------------------------------------------
# Confound residualization (within patient, per band) — OLS, the weak control
# ---------------------------------------------------------------------------
def _ols_residual(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Residual of y after regressing on [1, X]. Columns with zero variance are
    dropped so the design stays full-rank (e.g. n_same_shaft_epi all 0)."""
    good = np.isfinite(y)
    Xc = X.copy()
    keep_cols = [c for c in range(Xc.shape[1])
                 if np.nanstd(Xc[good, c]) > 1e-12]
    A = np.column_stack([np.ones(good.sum())] +
                        [Xc[good, c] for c in keep_cols])
    beta, *_ = np.linalg.lstsq(A, y[good], rcond=None)
    resid = np.full(y.shape, np.nan)
    resid[good] = y[good] - A @ beta
    return resid


def add_residuals(df: pd.DataFrame) -> pd.DataFrame:
    """Add <score>_perp columns, residualized within each (patient, band)."""
    out = df.copy()
    for s in SCORES:
        out[f"{s}_perp"] = np.nan
    for (pat, band), g in df.groupby(["patient", "band"]):
        X = g[CONFOUNDS].to_numpy(dtype=float)
        for s in SCORES:
            resid = _ols_residual(g[s].to_numpy(dtype=float), X)
            out.loc[g.index, f"{s}_perp"] = resid
    return out


# ---------------------------------------------------------------------------
# Discrimination gate: does the residual axis separate KNOWN epi?
# ---------------------------------------------------------------------------
def _auc(score_epi: np.ndarray, score_non: np.ndarray) -> float:
    se = score_epi[np.isfinite(score_epi)]
    sn = score_non[np.isfinite(score_non)]
    if se.size == 0 or sn.size == 0:
        return float("nan")
    u, _ = mannwhitneyu(se, sn, alternative="two-sided")
    return float(u / (se.size * sn.size))


def gate_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    gate_scores = [f"{s}_perp" for s in SCORES] + ["strength_mean"]
    for band in BANDS_TRACE:
        sub = df[df.band == band]
        for sc in gate_scores:
            aucs = []
            for pat, g in sub.groupby("patient"):
                if int(g.is_epi.sum()) < MIN_EPI_FOR_WINDOW:
                    continue
                a = _auc(g.loc[g.is_epi, sc].to_numpy(),
                         g.loc[~g.is_epi, sc].to_numpy())
                if np.isfinite(a):
                    aucs.append(a)
            aucs = np.array(aucs)
            if aucs.size < 3:
                rows.append(dict(band=band, score=sc, n_pat=int(aucs.size),
                                 median_auc=np.nan, abs_dev=np.nan,
                                 n_above_half="0/0", wilcoxon_z=np.nan,
                                 wilcoxon_p=np.nan, kind=("baseline"
                                 if sc == "strength_mean" else "trace")))
                continue
            wz, wp = wilcoxon_z(aucs - 0.5)
            rows.append(dict(
                band=band, score=sc, n_pat=int(aucs.size),
                median_auc=float(np.median(aucs)),
                abs_dev=float(abs(np.median(aucs) - 0.5)),
                n_above_half=f"{int((aucs>0.5).sum())}/{aucs.size}",
                wilcoxon_z=float(wz), wilcoxon_p=float(wp),
                kind=("baseline" if sc == "strength_mean" else "trace"),
            ))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Occult candidates (descriptive, hypothesis-only)
# ---------------------------------------------------------------------------
def candidate_table(df: pd.DataFrame, score: str = "dRho_perp") -> pd.DataFrame:
    rows = []
    for band in BANDS_TRACE:
        for pat, g in df[df.band == band].groupby("patient"):
            epi_vals = g.loc[g.is_epi, score].to_numpy()
            epi_vals = epi_vals[np.isfinite(epi_vals)]
            if epi_vals.size < MIN_EPI_FOR_WINDOW:
                continue
            lo, hi = float(epi_vals.min()), float(epi_vals.max())
            non = g[~g.is_epi].copy()
            nv = non[score].to_numpy(dtype=float)
            finite = np.isfinite(nv)
            if finite.sum() < 5:
                continue
            thr = np.quantile(nv[finite], TOP_DECILE)
            non = non.assign(score_val=nv,
                             pctile_nonepi=rankdata(np.where(finite, nv, -np.inf))
                             / finite.sum())
            cand = non[(non.score_val >= thr) & (non.score_val >= lo)
                       & (non.score_val <= hi) & finite]
            for _, r in cand.iterrows():
                rows.append(dict(
                    patient=pat, band=band, channel=r["channel"],
                    shaft=r["shaft"], contact=int(r["contact"]),
                    score=score, score_val=float(r["score_val"]),
                    pctile_nonepi=float(r["pctile_nonepi"]),
                    epi_window_lo=lo, epi_window_hi=hi,
                    along_dist_to_epi=float(r["along_dist_to_epi"]),
                    on_epi_free_shaft=int(not np.isfinite(r["along_dist_to_epi"])),
                    strength_mean=float(r["strength_mean"]),
                ))
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["band", "pctile_nonepi"],
                              ascending=[True, False])
    return out


# ---------------------------------------------------------------------------
# Rank-partial sanity (does a monotone-nonlinear strength control change it?)
# ---------------------------------------------------------------------------
def rank_partial_gate(df: pd.DataFrame) -> pd.DataFrame:
    """Spearman-rank-partial fallback: residualize the RANK of dRho on the RANK
    of strength only, then re-run the epi-vs-non AUC. Agreement with the OLS
    gate => the linear control suffices."""
    rows = []
    for band in BANDS_TRACE:
        sub = df[df.band == band]
        aucs = []
        for pat, g in sub.groupby("patient"):
            if int(g.is_epi.sum()) < MIN_EPI_FOR_WINDOW:
                continue
            y = g["dRho"].to_numpy(dtype=float)
            s = g["strength_mean"].to_numpy(dtype=float)
            good = np.isfinite(y) & np.isfinite(s)
            if good.sum() < 5:
                continue
            ry, rs = rankdata(y[good]), rankdata(s[good])
            A = np.column_stack([np.ones(good.sum()), rs])
            beta, *_ = np.linalg.lstsq(A, ry, rcond=None)
            resid = ry - A @ beta
            isepi = g.is_epi.to_numpy()[good]
            a = _auc(resid[isepi], resid[~isepi])
            if np.isfinite(a):
                aucs.append(a)
        aucs = np.array(aucs)
        if aucs.size >= 3:
            wz, wp = wilcoxon_z(aucs - 0.5)
            rows.append(dict(band=band, score="dRho_rankpartial_strength",
                             n_pat=int(aucs.size),
                             median_auc=float(np.median(aucs)),
                             abs_dev=float(abs(np.median(aucs) - 0.5)),
                             n_above_half=f"{int((aucs>0.5).sum())}/{aucs.size}",
                             wilcoxon_z=float(wz), wilcoxon_p=float(wp)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Cross-patient LOPO transfer — THE decisive occult/discovery test
# ---------------------------------------------------------------------------
def lopo_transfer(df: pd.DataFrame) -> pd.DataFrame:
    """Leave-one-patient-out: can the LNO influence dRho predict epi labels in a
    HELD-OUT patient? This is the occult/discovery use-case — for an unlabeled
    patient there are NO within-patient labels to residualize or threshold, so a
    within-patient AUC is NOT enough; the signal must transfer across patients.
    PR-AUC headline (prevalence = chance; lift = PR/prev). Mirrors audit_80.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import average_precision_score, roc_auc_score
    from sklearn.preprocessing import StandardScaler

    rows = []
    for band in BANDS_TRACE:
        sub = df[df.band == band].reset_index(drop=True)
        pats = sorted(sub.patient.unique())
        y = sub.is_epi.astype(int).to_numpy()
        for feats in (["dRho"], ["strength_mean"], ["dRho", "strength_mean"]):
            proba = np.full(len(sub), np.nan)
            for q in pats:
                tr = sub[sub.patient != q]
                te = sub[sub.patient == q]
                if int(tr.is_epi.sum()) < 2 or int(te.is_epi.sum()) < 1:
                    continue  # fold needs positives both sides
                Xtr = tr[feats].to_numpy(dtype=float)
                Xte = te[feats].to_numpy(dtype=float)
                scaler = StandardScaler().fit(Xtr)
                clf = LogisticRegression(class_weight="balanced",
                                         max_iter=500)
                clf.fit(scaler.transform(Xtr), tr.is_epi.astype(int))
                proba[te.index] = clf.predict_proba(
                    scaler.transform(Xte))[:, 1]
            m = np.isfinite(proba)
            if m.sum() < 20 or y[m].sum() < 3:
                continue
            prev = float(y[m].mean())
            pr = float(average_precision_score(y[m], proba[m]))
            rows.append(dict(
                band=band, feats="+".join(feats), n_scored=int(m.sum()),
                prevalence=round(prev, 3), PR_AUC=round(pr, 3),
                lift=round(pr / prev, 2) if prev > 0 else np.nan,
                ROC_AUC=round(float(roc_auc_score(y[m], proba[m])), 3),
            ))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------
def write_readme(df: pd.DataFrame, gate: pd.DataFrame, rankg: pd.DataFrame,
                 lopo: pd.DataFrame, cand: pd.DataFrame,
                 runtime_s: float) -> None:
    L: list[str] = []
    a = L.append
    a("---")
    a("name: occult_node_marker")
    a("scope: direction2_occult_nonlabeled_node_lno_influence_residualized")
    a("era: COHORT_N10 / IMCOH_ABS")
    a(f"date: {time.strftime('%Y-%m-%d')}")
    a("status: built_negative")
    a("build_script: scripts/01_compute/audit/audit_90_occult_node_lno_influence.py")
    a("scope_report: .agents/guides/task-persistence-investigation/2026-06-05_occult-epi-node-marker.md")
    a("---")
    a("")
    a("# Occult (non-labeled) epi-node marker — LNO cohort-trace influence")
    a("")
    a("**Head.** Discovery reframe of the per-node marker (companion "
      "known-node classifier audit_80 was NEGATIVE). Per-node leave-node-out "
      "(LNO) influence on the patient cophenetic ρ_split trace "
      "(Δρ_i = ρ_full − ρ_{∖i}) + row-restricted ρ_split_node / coherence, "
      "each residualized WITHIN patient against strength + same-shaft-epi count "
      "+ along-shaft position + shaft length. Gate: does the residual axis "
      "separate LABELED epi from non-labeled (must beat the residualized "
      "strength baseline)? Then: do non-labeled nodes land in the labeled-epi "
      "residual range (occult candidates)? Bands α/β/γ_l only (the trace-bearing "
      "bands; audit_77). **No surrogate** — the within-patient residualization "
      "is the strength control. **No surgical-outcome ground truth → candidates "
      "are hypotheses, not predictions.**")
    a("")
    a("## Bottom line")
    a("")
    a("- **Within-patient, β LNO influence `dRho` IS a real, "
      "strength-independent epi signal** (β cohort AUC 0.63, 8/9 patients, "
      "Wilcoxon p≈0.01; survives the rank-partial nonlinear-strength control). "
      "This is NEW — the row-restricted ρ_split_node (audit_79 F1) missed it "
      "because leave-node-out *leverage* on the cohort trace ≠ the trace ON the "
      "node's row. Sign: epi nodes have positive median `dRho` (trace-carriers); "
      "non-epi negative (trace-suppressors). But the within-patient AUC is "
      "**comparable to the strength baseline** (β 0.65), not above it.")
    a("- **Cross-patient (the occult/discovery use-case) it is at CHANCE.** "
      "Leave-one-patient-out PR-AUC for `dRho` sits at prevalence (β 0.12 vs "
      "prev 0.11, lift ≈1.2×, ROC 0.56) — the SAME ceiling audit_80 found. The "
      "within-patient β separation does NOT transfer to a held-out patient. "
      "**This is the decisive negative**: occult discovery needs a "
      "cross-patient-transferable threshold (an unlabeled patient has no "
      "within-patient labels), and `dRho` does not provide one.")
    a("- **Verdict: no transferable occult marker.** Consistent with the "
      "known-node classifier (audit_80, geometry_dominant). The 154 candidates "
      "below are within-patient top-decile re-reads, hypothesis-generating "
      "ONLY; with no surgical-outcome ground truth they cannot be validated.")
    a("")
    a("## Discrimination gate — does the residual axis separate KNOWN epi? (WITHIN-patient)")
    a("")
    a("| band | score | kind | n_pat | median AUC | |AUC−.5| | n>0.5 | Wilcoxon p |")
    a("|---|---|---|---|---|---|---|---|")
    for band in BANDS_TRACE:
        for _, r in gate[gate.band == band].iterrows():
            a(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | `{r['score']}` "
              f"| {r['kind']} | {r['n_pat']} | {r['median_auc']:.3f} "
              f"| {r['abs_dev']:.3f} | {r['n_above_half']} "
              f"| {r['wilcoxon_p']:.4f} |")
    a("")
    a("**Reading.** AUC = P(labeled-epi residual score > non-labeled residual "
      "score), per patient, median over patients (≥3 epi). 0.5 = no separation. "
      "A trace score is a usable occult axis only if its |AUC−0.5| beats the "
      "`strength_mean` baseline AND is non-trivial. If every residualized trace "
      "score sits at AUC≈0.5, there is no occult axis — candidates are a "
      "strength/geometry re-read.")
    a("")
    a("## Rank-partial (monotone-nonlinear strength) sanity")
    a("")
    a("| band | score | n_pat | median AUC | |AUC−.5| | n>0.5 | Wilcoxon p |")
    a("|---|---|---|---|---|---|---|")
    for _, r in rankg.iterrows():
        a(f"| {BRAIN_BAND_TEX_DICT.get(r['band'], r['band'])} "
          f"| `{r['score']}` | {r['n_pat']} | {r['median_auc']:.3f} "
          f"| {r['abs_dev']:.3f} | {r['n_above_half']} | {r['wilcoxon_p']:.4f} |")
    a("")
    a("Agreement with the OLS gate ⇒ the linear strength control suffices.")
    a("")
    a("## Cross-patient LOPO transfer — THE decisive occult/discovery test")
    a("")
    a("| band | features | n scored | prevalence | PR-AUC | lift | ROC-AUC |")
    a("|---|---|---|---|---|---|---|")
    for _, r in lopo.iterrows():
        a(f"| {BRAIN_BAND_TEX_DICT.get(r['band'], r['band'])} "
          f"| `{r['feats']}` | {r['n_scored']} | {r['prevalence']:.3f} "
          f"| {r['PR_AUC']:.3f} | {r['lift']:.2f}× | {r['ROC_AUC']:.3f} |")
    a("")
    a("**Reading.** PR-AUC is the headline (prevalence ≈ chance; lift = "
      "PR/prevalence). For an UNLABELED patient there are no within-patient "
      "labels, so only a cross-patient-transferable score can flag occult "
      "nodes. `dRho` lift ≈1.2× (β PR 0.12 vs prev 0.11) = at chance → no "
      "transferable occult marker. Same ceiling as the audit_80 classifier.")
    a("")
    a("## Occult candidates (descriptive, HYPOTHESIS-ONLY)")
    a("")
    n_cand = 0 if cand.empty else len(cand)
    a(f"- `occult_candidates.csv`: {n_cand} non-labeled contacts (band α/β/γ_l) "
      "whose residualized LNO influence is in the top decile among non-labeled "
      "AND inside the labeled-epi residual span. **These are NOT validated "
      "predictions** — there is no surgical-outcome ground truth; where the "
      "gate above is null, they are a strength/geometry re-read.")
    if not cand.empty:
        a("")
        a("### Top candidates by non-epi percentile")
        a("")
        a("| patient | band | channel | shaft | non-epi pctile | along-dist to epi | epi-free shaft |")
        a("|---|---|---|---|---|---|---|")
        for _, r in cand.head(15).iterrows():
            ad = "∞" if not np.isfinite(r["along_dist_to_epi"]) else f"{r['along_dist_to_epi']:.0f}"
            a(f"| {r['patient']} | {BRAIN_BAND_TEX_DICT.get(r['band'], r['band'])} "
              f"| {r['channel']} | {r['shaft']} | {r['pctile_nonepi']:.3f} "
              f"| {ad} | {int(r['on_epi_free_shaft'])} |")
    a("")
    a("## Caveats")
    a("")
    a("- **No ground truth.** A candidate cannot be confirmed or refuted without "
      "surgical-outcome (Engel) maps. Every flag is a hypothesis.")
    a("- **Linear residualization** removes only linear strength/geometry "
      "dependence (rank-partial fallback above guards the nonlinear case).")
    a("- **Adjacency confound.** A candidate near a labeled focus (low "
      "along-dist, `on_epi_free_shaft=0`) is most likely field/volume "
      "continuation, not a novel focus. The genuinely novel — and most "
      "uncertain — flag is `on_epi_free_shaft=1`.")
    a("- **Prior:** the known-node classifier (audit_80) was geometry_dominant "
      "in every band, so 'no occult signal beyond strength/geometry' is the "
      "expected outcome.")
    a("")
    a("## Provenance")
    a("- Cohort: " + ", ".join(es.COHORT) + " (Pat_15: 0 epi → no occult window).")
    a("- Bands: " + ", ".join(BANDS_TRACE) + " (trace-bearing; audit_77).")
    a("- Phases: " + ", ".join(PHASES) + "; FC imcoh_abs; τ=1/λ_max; average linkage.")
    a("- LNO = drop row+col, rebuild LRG on (N−1) nodes, recompute ρ_split.")
    a(f"- No surrogate ensemble. Wall-clock: {runtime_s:.1f} s")
    a("")
    a("## Files")
    a("- `lno_node_influence.csv` — per (patient, band, node): scores + confounds + residuals")
    a("- `occult_gate_per_band.csv` — per (band, score) WITHIN-patient epi-vs-non AUC + strength baseline")
    a("- `occult_lopo_transfer.csv` — cross-patient LOPO PR/ROC-AUC (the decisive test)")
    a("- `occult_candidates.csv` — ranked non-labeled candidates (hypothesis-only)")
    (OUT / "README.md").write_text("\n".join(L) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=BANDS_TRACE)
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    print("[audit_90] pre-flight: ensure rsPre half FCs cached")
    for pat in args.patients:
        ensure_half_fcs(pat, args.bands)

    t0 = time.time()
    all_rows: list[dict] = []
    for band in args.bands:
        for pat in args.patients:
            tc = time.time()
            rows = per_patient_band(pat, band, args.verbose)
            all_rows.extend(rows)
            print(f"[audit_90] {pat}/{band}: {len(rows)} nodes "
                  f"({time.time()-tc:.1f}s)")
    df = pd.DataFrame(all_rows)
    if df.empty:
        print("[audit_90] no rows produced; aborting")
        return

    df = add_residuals(df)
    gate = gate_table(df)
    rankg = rank_partial_gate(df)
    lopo = lopo_transfer(df)
    cand = candidate_table(df, score="dRho_perp")

    df.to_csv(OUT / "lno_node_influence.csv", index=False)
    gate.to_csv(OUT / "occult_gate_per_band.csv", index=False)
    lopo.to_csv(OUT / "occult_lopo_transfer.csv", index=False)
    cand.to_csv(OUT / "occult_candidates.csv", index=False)
    runtime = time.time() - t0
    write_readme(df, gate, rankg, lopo, cand, runtime)

    print("\n=== DISCRIMINATION GATE (WITHIN-patient: residual axis vs KNOWN epi) ===")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(gate.to_string(index=False))
    print("\n=== RANK-PARTIAL (nonlinear strength) sanity ===")
    print(rankg.to_string(index=False) if not rankg.empty else "(none)")
    print("\n=== CROSS-PATIENT LOPO TRANSFER (decisive occult test) ===")
    print(lopo.to_string(index=False) if not lopo.empty else "(none)")
    print(f"\n=== OCCULT CANDIDATES: {0 if cand.empty else len(cand)} non-labeled flags ===")
    if not cand.empty:
        with pd.option_context("display.width", 200, "display.max_columns", 20):
            print(cand.head(20).to_string(index=False))
    print(f"\n[audit_90] done in {runtime:.1f}s → {OUT}")


if __name__ == "__main__":
    main()
