#!/usr/bin/env python3
"""audit_144 — per-node trace decomposition + anti-trace node characterization.

Bottom-up explanator for the per-patient heterogeneity of the cophenetic trace.
Decomposes rho_split into per-node scores T_i (endpoint-incidence mean of the
rank-concordance; mean_i T_i is proportional to rho_split, so a low-rho_split
patient is one whose node-mean T_i is dragged down by anti-trace nodes), calibrates
each node against its OWN matched-strength surrogate -> carrier / neutral / anti,
pools nodes across the cohort, and characterizes the ANTI nodes (the negative space
audit_83 never described) by a-priori properties.

Scope report (read first):
  .agents/guides/task-persistence-investigation/2026-06-25_per-node-trace-decomposition.md

5-point critical preamble: see scope report §1. In short — claim: the heterogeneity
is carried by anti-trace nodes sharing an a-priori property P; null (polarity):
matched-strength 4-cycle ±δ surrogate (audit_63 R=200 ensemble, seed 20260511),
per-node aggregated identically to the observed; null (closing test, separate
script): size-matched random-node decimation (audit_85 --mode node). This script
produces the per-node scores + the discovery-step characterization ONLY; it does NOT
run the closing test (that is gated on a candidate P emerging here) and does NOT
re-test trace existence (locked §5.2/5.4).

Reuses (library-first):
  node_localization.{rank_concordance, canonical_cophenet, node_incidence_mean,
                     epi_keep_mask, shaft_of, NON_ANATOMICAL}
  surrogate.{cophenetic_condensed_from_eigs, surrogate_cache_path}
  io.regions.load_channel_regions ; io.patient.load_epileptic_nodes ; workflow.fc.load_fc_matrix

Outputs (data/audit/per_node_trace_decomposition/)
  per_node.csv               one row per (patient, band, node): T_raw, t_rel,
                             p_plus, p_minus, polarity, pi_tilde + properties X_i
  per_patient_summary.csv    per (patient, band): rho_split + carrier/anti/neutral counts
  characterization_{band}.csv per (property, level): anti/carrier fractions +
                             within-patient-shuffle reference (the candidate-P shortlist)

Usage
  python audit_144_per_node_trace_decomposition.py
  python audit_144_per_node_trace_decomposition.py --bands beta --patients Pat_06
"""
from __future__ import annotations

import argparse
import re

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.io.patient import load_epileptic_nodes
from lrg_eegfc.utils.io.regions import _normalise_label, load_channel_regions
from lrg_eegfc.utils.metrics.node_localization import (
    canonical_cophenet, epi_keep_mask, node_incidence_mean, rank_concordance,
    shaft_of,
)
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.utils.surrogate import (
    cophenetic_condensed_from_eigs, surrogate_cache_path,
)
from lrg_eegfc.workflow.fc import load_fc_matrix

ROOT = setup_script_env()
OUT = ROOT / "data/audit/per_node_trace_decomposition"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
TARGET_BANDS = ["beta", "alpha", "low_gamma"]   # surrogate cache exists at R=200
PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
R, SWAP, SEED = 200, 20, 20260511

#: categorical / continuous node properties characterized in the discovery step
CAT_PROPS = ["system", "tissue", "is_epi", "hemisphere", "lobe", "supersystem"]
CONT_PROPS = ["depth_norm", "strength"]
N_SHUFFLE = 2000


# ---------------------------------------------------------------------------
# observed + surrogate cophenetic trace (identical construction to audit_83)
# ---------------------------------------------------------------------------
def obs_cophenet_trace(pat, band):
    """Per-pair rank-concordance c_p (length P) + upper-tri indices + N."""
    Wt = load_fc_matrix(pat, "task_test", band, "imcoh_abs")
    Wpost = load_fc_matrix(pat, "rest_post", band, "imcoh_abs")
    WA = np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_A_imcoh_abs.npy")
    WB = np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_B_imcoh_abs.npy")
    Dt, Dpost, DA, DB = (canonical_cophenet(W) for W in (Wt, Wpost, WA, WB))
    c = rank_concordance(Dt - DA, Dpost - DB)
    rho = float(spearmanr(Dt - DA, Dpost - DB).statistic)
    N = Wt.shape[0]
    iu = np.triu_indices(N, 1)
    return c, iu[0].astype(int), iu[1].astype(int), N, rho


def load_surrogate_cophenet(pat, band):
    """Per-phase surrogate cophenetic condensed distances: dict phase -> list[R]."""
    out = {}
    for ph in PHASES:
        path = surrogate_cache_path(pat, band, ph, R, SWAP, SEED, "imcoh_abs")
        if not path.exists():
            return None
        with np.load(path) as d:
            evals, evecs = d["eigvals"], d["eigvecs"]
        cond = []
        for r in range(evals.shape[0]):
            ev = evals[r]
            cond.append(cophenetic_condensed_from_eigs(ev, evecs[r])
                        if np.all(np.isfinite(ev)) else None)
        out[ph] = cond
    return out


# ---------------------------------------------------------------------------
# node properties (read-only; all a-priori / trace-independent)
# ---------------------------------------------------------------------------
_NUM_RE = re.compile(r"[A-Za-z']+\s*0*([0-9]+)")


def _contact_number(label_raw: str):
    m = _NUM_RE.match(str(label_raw))
    return int(m.group(1)) if m else np.nan


def node_properties(pat, band) -> pd.DataFrame:
    df = load_channel_regions(pat).copy()
    # de-duplicate the doubled 'label' column quirk
    df = df.loc[:, ~df.columns.duplicated()]
    N = len(df)
    df["node_idx"] = np.arange(N)
    df["tissue"] = np.where(df["system"].to_numpy() == "non_anatomical",
                            "nonanat", "gray")
    keep = epi_keep_mask(df, pat)              # True = non-epi
    df["is_epi"] = ~keep
    df["shaft"] = [shaft_of(l) for l in df["label_raw"]]
    df["contact_num"] = [_contact_number(l) for l in df["label_raw"]]
    # within-shaft normalized depth (rank by contact number); singletons -> 0.5
    depth = np.full(N, 0.5)
    for sh, idx in df.groupby("shaft").groups.items():
        idx = list(idx)
        if len(idx) > 1:
            order = df.loc[idx, "contact_num"].rank(method="first").to_numpy()
            depth[idx] = (order - 1) / (len(idx) - 1)
    df["depth_norm"] = depth
    # node strength (hubness): mean over the 3 main phases
    smat = []
    for ph in ("rest_pre", "task_test", "rest_post"):
        W = load_fc_matrix(pat, ph, band, "imcoh_abs")
        smat.append(W.sum(axis=1))
    df["strength"] = np.mean(np.vstack(smat), axis=0)
    return df


# ---------------------------------------------------------------------------
# per-cell decomposition
# ---------------------------------------------------------------------------
def per_cell(pat, band, verbose=True):
    surco = load_surrogate_cophenet(pat, band)
    if surco is None:
        if verbose:
            print(f"  [skip] {pat} {band}: surrogate cache missing")
        return None
    try:
        c, iu_i, iu_j, N, rho = obs_cophenet_trace(pat, band)
    except FileNotFoundError:
        if verbose:
            print(f"  [skip] {pat} {band}: half-FC missing")
        return None

    T = node_incidence_mean(c, iu_i, iu_j, N, demean=False)
    t = node_incidence_mean(c, iu_i, iu_j, N, demean=True)

    # surrogate per-node scores
    Ts = []
    for r in range(R):
        if any(surco[ph][r] is None for ph in PHASES):
            continue
        cr = rank_concordance(surco["task_test"][r] - surco["rest_pre_A"][r],
                              surco["rest_post"][r] - surco["rest_pre_B"][r])
        Ts.append(node_incidence_mean(cr, iu_i, iu_j, N, demean=False))
    Ts = np.vstack(Ts) if Ts else np.full((1, N), np.nan)   # (Rv, N)
    Rv = Ts.shape[0]
    Qplus = np.nanquantile(Ts, 0.95, axis=0)
    Qminus = np.nanquantile(Ts, 0.05, axis=0)
    pplus = (1 + np.sum(Ts >= T[None, :], axis=0)) / (Rv + 1)
    pminus = (1 + np.sum(Ts <= T[None, :], axis=0)) / (Rv + 1)
    polarity = np.where(T >= Qplus, "carrier",
                        np.where(T <= Qminus, "anti", "neutral"))
    pi_tilde = -np.log10(pplus) + np.log10(pminus)

    props = node_properties(pat, band)
    props["band"] = band
    props["patient"] = pat
    props["rho_split"] = rho
    props["T_raw"] = T
    props["t_rel"] = t
    props["surr_mean_T"] = np.nanmean(Ts, axis=0)
    props["p_plus"] = pplus
    props["p_minus"] = pminus
    props["polarity"] = polarity
    props["pi_tilde"] = pi_tilde
    props["n_surr_valid"] = Rv
    if verbose:
        nc = int((polarity == "carrier").sum())
        na = int((polarity == "anti").sum())
        print(f"  {pat} {band:9s}: rho={rho:+.3f}  N={N}  "
              f"carrier={nc} anti={na} neutral={N-nc-na}  (Rv={Rv})")
    return props


# ---------------------------------------------------------------------------
# discovery step: characterize anti vs carrier by a-priori property
# ---------------------------------------------------------------------------
def characterize(pool: pd.DataFrame, band: str, rng) -> pd.DataFrame:
    """Per property level: anti/carrier fractions + within-patient-shuffle p.

    Within-patient shuffle permutes the polarity labels inside each patient
    (preserving each patient's carrier/anti/neutral counts), so a level that
    collects more anti nodes than its per-patient base-rate share is flagged.
    """
    sub = pool[pool.band == band].copy()
    pats = sub.patient.to_numpy()
    pol = sub.polarity.to_numpy()
    is_anti = (pol == "anti").astype(float)
    is_carr = (pol == "carrier").astype(float)
    # precompute within-patient shuffles of the polarity index
    pat_codes, pat_idx = np.unique(pats, return_inverse=True)
    groups = [np.where(pat_idx == g)[0] for g in range(len(pat_codes))]
    shuffles = []
    for _ in range(N_SHUFFLE):
        perm = np.arange(len(sub))
        for g in groups:
            perm[g] = rng.permutation(g)
        shuffles.append(perm)
    shuffles = np.array(shuffles)                       # (N_SHUFFLE, n)

    rows = []
    for prop in CAT_PROPS:
        levels = sub[prop].astype(str).to_numpy()
        for lv in pd.unique(levels):
            sel = levels == lv
            n = int(sel.sum())
            if n < 5:
                continue
            obs_anti = is_anti[sel].mean()
            obs_carr = is_carr[sel].mean()
            sh_anti = is_anti[shuffles[:, sel]].mean(axis=1)
            sh_carr = is_carr[shuffles[:, sel]].mean(axis=1)
            p_anti = (1 + np.sum(sh_anti >= obs_anti)) / (N_SHUFFLE + 1)
            p_carr = (1 + np.sum(sh_carr >= obs_carr)) / (N_SHUFFLE + 1)
            rows.append({
                "band": band, "property": prop, "level": lv, "n_nodes": n,
                "anti_frac": obs_anti, "carrier_frac": obs_carr,
                "mean_pi_tilde": float(sub.loc[sel, "pi_tilde"].mean()),
                "shuffle_p_anti": p_anti, "shuffle_p_carrier": p_carr,
            })
    # continuous: anti vs carrier mean contrast (within-patient-shuffle on the gap)
    for prop in CONT_PROPS:
        v = sub[prop].to_numpy(dtype=float)
        ok = np.isfinite(v)
        gap_obs = (np.nanmean(v[(pol == "anti") & ok])
                   - np.nanmean(v[(pol == "carrier") & ok]))
        sh_gap = []
        for perm in shuffles:
            pp = pol[perm]
            sh_gap.append(np.nanmean(v[(pp == "anti") & ok])
                          - np.nanmean(v[(pp == "carrier") & ok]))
        sh_gap = np.array(sh_gap)
        p_hi = (1 + np.sum(sh_gap >= gap_obs)) / (N_SHUFFLE + 1)
        p_lo = (1 + np.sum(sh_gap <= gap_obs)) / (N_SHUFFLE + 1)
        rows.append({
            "band": band, "property": prop, "level": "anti_minus_carrier_mean",
            "n_nodes": int(ok.sum()), "anti_frac": float(np.nanmean(v[(pol == "anti") & ok])),
            "carrier_frac": float(np.nanmean(v[(pol == "carrier") & ok])),
            "mean_pi_tilde": gap_obs,
            "shuffle_p_anti": p_hi, "shuffle_p_carrier": p_lo,
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=TARGET_BANDS)
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--seed", type=int, default=20260625)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    all_nodes = []
    for band in args.bands:
        print(f"== band {band} ==")
        for pat in args.patients:
            r = per_cell(pat, band)
            if r is not None:
                all_nodes.append(r)
    pool = pd.concat(all_nodes, ignore_index=True)

    keep_cols = ["patient", "band", "node_idx", "label_raw", "region", "lobe",
                 "system", "supersystem", "hemisphere", "tissue", "is_epi",
                 "shaft", "contact_num", "depth_norm", "strength", "x", "y", "z",
                 "rho_split", "T_raw", "t_rel", "surr_mean_T", "p_plus",
                 "p_minus", "polarity", "pi_tilde", "n_surr_valid"]
    keep_cols = [c for c in keep_cols if c in pool.columns]
    pool[keep_cols].to_csv(OUT / "per_node.csv", index=False)

    # per-patient summary
    summ = (pool.groupby(["patient", "band"])
            .agg(N=("node_idx", "size"),
                 rho_split=("rho_split", "first"),
                 n_carrier=("polarity", lambda s: (s == "carrier").sum()),
                 n_anti=("polarity", lambda s: (s == "anti").sum()),
                 n_neutral=("polarity", lambda s: (s == "neutral").sum()))
            .reset_index())
    summ.to_csv(OUT / "per_patient_summary.csv", index=False)

    # characterization per band
    for band in args.bands:
        ch = characterize(pool, band, rng)
        ch = ch.sort_values("shuffle_p_anti")
        ch.to_csv(OUT / f"characterization_{band}.csv", index=False)
        print(f"\n=== anti-node characterization : {band} "
              f"(top levels by within-patient-shuffle p_anti) ===")
        show = ch[ch.shuffle_p_anti < 0.2].head(12)
        if len(show):
            print(show[["property", "level", "n_nodes", "anti_frac",
                        "carrier_frac", "shuffle_p_anti",
                        "shuffle_p_carrier"]].to_string(index=False))
        else:
            print("  (no level with shuffle_p_anti < 0.2)")

    print(f"\nOutputs -> {OUT}")
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
