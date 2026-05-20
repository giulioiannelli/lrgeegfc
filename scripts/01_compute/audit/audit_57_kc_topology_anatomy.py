#!/usr/bin/env python3
"""Audit 57 — KC λ=0 (topology) trace anatomical-distribution audit.

Parallel of `audit_56_anatomy_distribution.py` on a *structurally distinct*
object: the §5.2 KC λ=0 per-leaf topological-memory catalog. Asks whether
the β KC topology trace admits an anatomical address (single-region,
lobar, hemispheric) or spreads across cortex.

Predicate (NOT redefined here — reused from §5.2 / audit_53)
------------------------------------------------------------
Per pair (i, j) within a (patient, band) cell, on the KC m-vectors of
imcoh_abs LRG dendrograms at the three phases:

    T(i, j) = ( |m_tt − m_post| ≤ δ ) AND ( |m_pre − m_tt| > δ )
              AND ( m_tt ≥ min_m AND m_post ≥ min_m )
              AND ( size_tt ≤ ⌈max_mrca_frac · n⌉
                    AND size_post ≤ ⌈max_mrca_frac · n⌉ )

Per leaf:
    f_i        = | { j : T(i, j) = 1 } | / (n − 1)

The within-patient null is the *anti-trace* predicate (rspre↔rspost
agreement, task differing) with all other guards mirrored:

    T_anti(i, j) = ( |m_pre − m_post| ≤ δ ) AND ( |m_tt − m_pre| > δ )
                   AND nontrivial-m and size-cap (all three phases)

    f_i^anti   = | { j : T_anti(i, j) = 1 } | / (n − 1)

Topology-trace flag at (patient, band, leaf):
    leaf is FLAGGED iff f_i > p95(f_i^anti at that cell) AND f_i > 0
    (one-sided upper tail; ε guard against degenerate p95=0 cells).

Defaults match audit_53: δ = 1, min_m = 2, max_mrca_frac = 0.50.

Tests (mirrors audit_56)
------------------------
1. Single-region hypergeometric (m = 48).
2. Hemisphere (m_hemi = 18) and lobar (m_lobe = 42) aggregation.
3. Concentration scalars (G, N_eff, s_1, s_3) + coverage-weighted
   bootstrap null.
4. Patient base distribution.
5. Focal-leaf anatomy: report Desikan-Killiany region + lobe +
   hemisphere for the §5.2-cited focal leaves (Pat_07 β=62,
   Pat_05 β=66, Pat_10 β=37, Pat_14 high-γ=62).

Sensitivity regimes: full / no_pat03 / pro_cohort (drop Pat_07+Pat_15).

Bands (Option A — KC λ=0-aligned only): β (primary), α (control).
γ_low is excluded (KC heights band, not topology); a parallel λ=1 audit
is documented as a deferred Option B in the README.

Outputs
-------
data/kc_topology_anatomy/region_count_table.csv
data/kc_topology_anatomy/concentration_scalars.csv
data/kc_topology_anatomy/hemisphere_lobe_enrichment.csv
data/kc_topology_anatomy/patient_base_distribution.csv
data/kc_topology_anatomy/focal_leaf_anatomy.csv
data/kc_topology_anatomy/per_topology_trace_leaf.csv
data/kc_topology_anatomy/figures/topology_trace_violins.pdf
data/kc_topology_anatomy/figures/topology_lobar_heatmap.pdf
data/kc_topology_anatomy/README.md
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics.tree_distance import (
    _ancestor_chain,
    _build_parent_height_depth,
    kc_vectors,
)
from lrg_eegfc.workflow.lrg import load_lrg_result

plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
})

# --------------------------------------------------------------------------
# Cohort + parameters
# --------------------------------------------------------------------------
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
TARGET_BANDS = ["alpha", "beta"]
PHASES = ("rest_pre", "task_test", "rest_post")

# §5.2 / audit_53 predicate parameters
DELTA = 1
MIN_M = 2
MAX_MRCA_FRAC = 0.50

# Threshold quantile (one-sided upper tail of f_i^anti)
NULL_QUANTILE = 0.95
EPS = 1e-12

BONFERRONI_M = 48      # per spec — matches §5.5 manuscript denominator
B_BOOT = 1000
RNG_SEED = 20260508
SENSITIVITIES = ["full", "no_pat03", "pro_cohort"]
PRO_COHORT_DROP = {"Pat_07", "Pat_15"}

# §5.2-cited focal leaves (manuscript Fig. manuscript_kc_leaf_topology.pdf)
FOCAL_CITED = [
    ("Pat_07", "beta", 62),
    ("Pat_05", "beta", 66),
    ("Pat_10", "beta", 37),
    ("Pat_14", "high_gamma", 62),  # §5.2 reference cell (sentinel: empty predicate)
]

OUT = ROOT / "data" / "kc_topology_anatomy"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

# Lobe map (copied verbatim from audit_56 to keep the two audits comparable)
LOBE_MAP = {
    "ctx-lh-superiorfrontal": "frontal", "ctx-rh-superiorfrontal": "frontal",
    "ctx-lh-rostralmiddlefrontal": "frontal", "ctx-rh-rostralmiddlefrontal": "frontal",
    "ctx-lh-caudalmiddlefrontal": "frontal", "ctx-rh-caudalmiddlefrontal": "frontal",
    "ctx-lh-parsopercularis": "frontal", "ctx-rh-parsopercularis": "frontal",
    "ctx-lh-parstriangularis": "frontal", "ctx-rh-parstriangularis": "frontal",
    "ctx-lh-parsorbitalis": "frontal", "ctx-rh-parsorbitalis": "frontal",
    "ctx-lh-lateralorbitofrontal": "frontal", "ctx-rh-lateralorbitofrontal": "frontal",
    "ctx-lh-medialorbitofrontal": "frontal", "ctx-rh-medialorbitofrontal": "frontal",
    "ctx-lh-precentral": "frontal", "ctx-rh-precentral": "frontal",
    "ctx-lh-paracentral": "frontal", "ctx-rh-paracentral": "frontal",
    "ctx-lh-frontalpole": "frontal", "ctx-rh-frontalpole": "frontal",
    "ctx-lh-superiorparietal": "parietal", "ctx-rh-superiorparietal": "parietal",
    "ctx-lh-inferiorparietal": "parietal", "ctx-rh-inferiorparietal": "parietal",
    "ctx-lh-supramarginal": "parietal", "ctx-rh-supramarginal": "parietal",
    "ctx-lh-postcentral": "parietal", "ctx-rh-postcentral": "parietal",
    "ctx-lh-precuneus": "parietal", "ctx-rh-precuneus": "parietal",
    "ctx-lh-superiortemporal": "temporal", "ctx-rh-superiortemporal": "temporal",
    "ctx-lh-middletemporal": "temporal", "ctx-rh-middletemporal": "temporal",
    "ctx-lh-inferiortemporal": "temporal", "ctx-rh-inferiortemporal": "temporal",
    "ctx-lh-fusiform": "temporal", "ctx-rh-fusiform": "temporal",
    "ctx-lh-transversetemporal": "temporal", "ctx-rh-transversetemporal": "temporal",
    "ctx-lh-bankssts": "temporal", "ctx-rh-bankssts": "temporal",
    "ctx-lh-entorhinal": "temporal", "ctx-rh-entorhinal": "temporal",
    "ctx-lh-parahippocampal": "temporal", "ctx-rh-parahippocampal": "temporal",
    "ctx-lh-temporalpole": "temporal", "ctx-rh-temporalpole": "temporal",
    "ctx-lh-lateraloccipital": "occipital", "ctx-rh-lateraloccipital": "occipital",
    "ctx-lh-lingual": "occipital", "ctx-rh-lingual": "occipital",
    "ctx-lh-cuneus": "occipital", "ctx-rh-cuneus": "occipital",
    "ctx-lh-pericalcarine": "occipital", "ctx-rh-pericalcarine": "occipital",
    "ctx-lh-isthmuscingulate": "limbic", "ctx-rh-isthmuscingulate": "limbic",
    "ctx-lh-caudalanteriorcingulate": "limbic",
    "ctx-rh-caudalanteriorcingulate": "limbic",
    "ctx-lh-posteriorcingulate": "limbic", "ctx-rh-posteriorcingulate": "limbic",
    "ctx-lh-rostralanteriorcingulate": "limbic",
    "ctx-rh-rostralanteriorcingulate": "limbic",
    "Hip": "limbic", "Amy": "limbic",
    "ctx-lh-insula": "insular", "ctx-rh-insula": "insular",
}


def hemisphere_of(region: str) -> str:
    if region.startswith("ctx-lh-"):
        return "left"
    if region.startswith("ctx-rh-"):
        return "right"
    if region in ("Hip", "Amy"):
        return "left"
    return "midline"


def lobe_of(region: str) -> str:
    return LOBE_MAP.get(region, "other")


# --------------------------------------------------------------------------
# Anatomy + cortical pool — per-leaf, parallel to §5.5 catalog
# --------------------------------------------------------------------------
from lrg_eegfc.utils.io.patient import (
    PATIENT_CHANNEL_DROP,  # for label drops (Pat_10 rows 53,54,55)
)


def _normalize_label(raw: str) -> str:
    """'A 1,G2' -> 'A1' (drop reference suffix + collapse spaces). Mirrors
    `audit_44_anatomical_mapping.normalize_label`."""
    return str(raw).split(",")[0].replace(" ", "").upper()


def _parse_dk_dominant(dk) -> str:
    """First entry of the comma-list is the dominant region."""
    if pd.isna(dk):
        return "Unk"
    return str(dk).split(",")[0].strip()


def load_per_leaf_anatomy(patient: str) -> pd.DataFrame:
    """Return per-leaf anatomy with columns
    ``[patient, leaf_id, channel, channel_raw, x, y, z, region]``.
    ``leaf_id`` is the canonical post-drop row index, matching the LRG
    linkage leaves. Implant lookup is by **normalized** channel label
    (matches `audit_44_anatomical_mapping`).
    """
    pat_dir = ROOT / "data" / "raw" / "stereoeeg_patients" / patient
    labels = pd.read_csv(pat_dir / "channel_labels.csv")
    impl = pd.read_csv(pat_dir / f"implant_pat_{patient[-2:]}.csv")
    dk_col = next(c for c in impl.columns if c.startswith("Desikan"))

    label_drops = PATIENT_CHANNEL_DROP.get(patient, {}).get("__labels__", [])
    if label_drops:
        labels = labels.drop(index=list(label_drops)).reset_index(drop=True)

    impl["norm_label"] = impl["label"].astype(str).apply(_normalize_label)
    impl_map = impl.set_index("norm_label")

    rows = []
    for leaf_id, raw in enumerate(labels["label"].astype(str).tolist()):
        ch = _normalize_label(raw)
        rec = {"patient": patient, "leaf_id": int(leaf_id),
               "channel": ch, "channel_raw": raw,
               "x": np.nan, "y": np.nan, "z": np.nan, "region": "Unk"}
        if ch in impl_map.index:
            r = impl_map.loc[ch]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            rec["x"] = float(str(r["x"]).replace(",", "."))
            rec["y"] = float(str(r["y"]).replace(",", "."))
            rec["z"] = float(str(r["z"]).replace(",", "."))
            rec["region"] = _parse_dk_dominant(r[dk_col])
        rows.append(rec)
    return pd.DataFrame(rows)


def load_pool() -> pd.DataFrame:
    df = pd.concat([load_per_leaf_anatomy(p) for p in PATIENTS],
                   ignore_index=True)
    df = df[~df.region.isin(["Wm", "Unk"])].reset_index(drop=True)
    df["hemisphere"] = df.region.map(hemisphere_of)
    df["lobe"] = df.region.map(lobe_of)
    return df


# --------------------------------------------------------------------------
# §5.2 predicate machinery (mirrors audit_53)
# --------------------------------------------------------------------------
def pair_index(n: int, i: int, j: int) -> int:
    if i > j:
        i, j = j, i
    return n * i - (i * (i + 1)) // 2 + (j - i - 1)


def subtree_sizes(Z: np.ndarray) -> np.ndarray:
    n = Z.shape[0] + 1
    sizes = np.ones(2 * n - 1, dtype=np.int64)
    for k in range(n - 1):
        a, b = int(Z[k, 0]), int(Z[k, 1])
        sizes[n + k] = sizes[a] + sizes[b]
    return sizes


def pair_mrca_subtree_sizes(Z: np.ndarray) -> np.ndarray:
    n = Z.shape[0] + 1
    parent, _, _ = _build_parent_height_depth(Z)
    sizes = subtree_sizes(Z)
    chains = [_ancestor_chain(i, parent) for i in range(n)]
    chain_len = np.array([len(c) for c in chains], dtype=np.int64)
    npairs = n * (n - 1) // 2
    out = np.zeros(npairs, dtype=np.int64)
    idx = 0
    for i in range(n):
        ci = chains[i]
        for j in range(i + 1, n):
            cj = chains[j]
            limit = min(chain_len[i], chain_len[j])
            mrca = ci[0]
            for k in range(limit):
                if ci[k] == cj[k]:
                    mrca = ci[k]
                else:
                    break
            out[idx] = sizes[mrca]
            idx += 1
    return out


def expand_per_leaf_score(predicate_vec: np.ndarray, n: int) -> np.ndarray:
    f = np.zeros(n, dtype=float)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if predicate_vec[idx]:
                f[i] += 1.0
                f[j] += 1.0
            idx += 1
    f /= max(n - 1, 1)
    return f


def topology_trace_flags(patient: str, band: str) -> dict:
    """Return per-leaf topology-trace flags + diagnostics for a (patient, band)
    cell. Returns ``None`` keys ``flags``, ``f``, ``f_anti``, ``threshold``
    or raises FileNotFoundError if any phase is missing."""
    Z_by_phase = {}
    for ph in PHASES:
        res = load_lrg_result(patient, ph, band, fc_method="imcoh_abs")
        if res is None or res.linkage_matrix is None:
            raise FileNotFoundError(f"{patient} {band} {ph}: missing LRG")
        Z_by_phase[ph] = np.asarray(res.linkage_matrix)
    n = int(Z_by_phase["rest_pre"].shape[0]) + 1
    for ph in PHASES:
        if int(Z_by_phase[ph].shape[0]) + 1 != n:
            raise ValueError(f"{patient} {band}: phase leaf-count mismatch")

    m_pre = kc_vectors(Z_by_phase["rest_pre"])[0]
    m_tt = kc_vectors(Z_by_phase["task_test"])[0]
    m_post = kc_vectors(Z_by_phase["rest_post"])[0]
    size_pre = pair_mrca_subtree_sizes(Z_by_phase["rest_pre"])
    size_tt = pair_mrca_subtree_sizes(Z_by_phase["task_test"])
    size_post = pair_mrca_subtree_sizes(Z_by_phase["rest_post"])

    s_max = int(np.ceil(MAX_MRCA_FRAC * n))

    # Trace predicate (audit_53 verbatim)
    matched_tt_post = np.abs(m_tt - m_post) <= DELTA
    differs_pre_tt = np.abs(m_pre - m_tt) > DELTA
    nontrivial_m_trace = (m_tt >= MIN_M) & (m_post >= MIN_M)
    fine_grained_trace = (size_tt <= s_max) & (size_post <= s_max)
    T = matched_tt_post & differs_pre_tt & nontrivial_m_trace & fine_grained_trace

    # Anti-trace predicate (rspre <-> rspost agreement, task differing)
    matched_pre_post = np.abs(m_pre - m_post) <= DELTA
    differs_tt_pre = np.abs(m_tt - m_pre) > DELTA
    nontrivial_m_anti = (m_pre >= MIN_M) & (m_post >= MIN_M)
    fine_grained_anti = (size_pre <= s_max) & (size_post <= s_max)
    T_anti = matched_pre_post & differs_tt_pre & nontrivial_m_anti & fine_grained_anti

    f = expand_per_leaf_score(T, n)
    f_anti = expand_per_leaf_score(T_anti, n)

    # Threshold = NULL_QUANTILE percentile of f_anti distribution (within-cell);
    # flag iff f > thresh (strict) AND f > 0.
    thresh = float(np.quantile(f_anti, NULL_QUANTILE))
    flags = (f > thresh + EPS) & (f > 0.0)

    return dict(n=n, f=f, f_anti=f_anti, threshold=thresh, flags=flags,
                n_trace_pairs=int(T.sum()),
                n_anti_pairs=int(T_anti.sum()))


# --------------------------------------------------------------------------
# Build per-topology-trace-leaf catalog (parallel of §5.5 per_trace_leaf.csv)
# --------------------------------------------------------------------------
def build_catalog(anatomy_by_pat: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for pat in PATIENTS:
        anat = anatomy_by_pat[pat]
        for band in TARGET_BANDS:
            try:
                d = topology_trace_flags(pat, band)
            except (FileNotFoundError, ValueError) as e:
                print(f"  skip {pat} {band}: {e}")
                continue
            n = d["n"]
            if len(anat) != n:
                print(f"  WARN {pat} {band}: anatomy rows {len(anat)} vs n {n}")
                continue
            for leaf_id in range(n):
                if not d["flags"][leaf_id]:
                    continue
                r = anat.iloc[leaf_id]
                rows.append({
                    "patient": pat, "band": band,
                    "leaf_id": leaf_id,
                    "channel": r["channel"],
                    "region": r["region"],
                    "x": r["x"], "y": r["y"], "z": r["z"],
                    "f_i": float(d["f"][leaf_id]),
                    "f_i_anti": float(d["f_anti"][leaf_id]),
                    "threshold": float(d["threshold"]),
                })
    cat = pd.DataFrame(rows)
    cat = cat[~cat.region.isin(["Wm", "Unk"])].reset_index(drop=True)
    cat["hemisphere"] = cat.region.map(hemisphere_of)
    cat["lobe"] = cat.region.map(lobe_of)
    return cat


# --------------------------------------------------------------------------
# Anatomy distribution pipeline (functions ported from audit_56)
# --------------------------------------------------------------------------
def restrict(df: pd.DataFrame, regime: str) -> pd.DataFrame:
    if regime == "full":
        return df
    if regime == "no_pat03":
        return df[df.patient != "Pat_03"]
    if regime == "pro_cohort":
        return df[~df.patient.isin(PRO_COHORT_DROP)]
    raise ValueError(regime)


def per_region_table(traces: pd.DataFrame, pool: pd.DataFrame,
                     band: str, regime: str) -> pd.DataFrame:
    sub_traces = restrict(traces[traces.band == band], regime)
    sub_pool = restrict(pool, regime)
    K = len(sub_traces)
    N_total = len(sub_pool)
    base_per = sub_pool.groupby("region").size().rename("n_contacts").reset_index()
    if K == 0:
        return pd.DataFrame()
    base_rate = K / N_total
    per = (sub_traces.groupby("region")
           .agg(n_trace=("leaf_id", "count"),
                n_contributing_patients=("patient", "nunique"))
           .reset_index())
    m = per.merge(base_per, on="region", how="left")
    m["trace_rate"] = m.n_trace / m.n_contacts
    m["enrichment"] = m.trace_rate / base_rate
    m["p_hyper"] = m.apply(
        lambda r: float(hypergeom.sf(r.n_trace - 1, N_total, K, int(r.n_contacts))),
        axis=1)
    m["p_bonferroni"] = (m["p_hyper"] * BONFERRONI_M).clip(upper=1.0)
    m["band"] = band
    m["sensitivity_regime"] = regime
    m["K"] = K
    m["N_total"] = N_total
    m["baseline_rate"] = base_rate
    m["eligible"] = (m.n_contacts >= 5) & (m.n_contributing_patients >= 2)
    return m


def gini(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    if x.sum() == 0 or x.size == 0:
        return float("nan")
    x = np.sort(x)
    n = x.size
    cum = np.cumsum(x)
    return float((n + 1 - 2 * np.sum(cum) / cum[-1]) / n)


def n_eff(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    if x.sum() == 0 or x.size == 0:
        return float("nan")
    p = x / x.sum()
    p = p[p > 0]
    return float(np.exp(-np.sum(p * np.log(p))))


def top_k_share(x: np.ndarray, k: int) -> float:
    x = np.asarray(x, dtype=float)
    if x.sum() == 0 or x.size == 0:
        return float("nan")
    return float(np.sum(np.sort(x)[::-1][:k]) / x.sum())


def concentration_scalars(counts: np.ndarray) -> dict:
    return dict(G=gini(counts), N_eff=n_eff(counts),
                s_1=top_k_share(counts, 1), s_3=top_k_share(counts, 3))


def bootstrap_null(eligible_regions: list[str], pool: pd.DataFrame,
                   K: int, B: int, rng: np.random.Generator) -> dict:
    pool_arr = pool["region"].to_numpy()
    n_pool = len(pool_arr)
    rec = {key: np.empty(B) for key in ("G", "N_eff", "s_1", "s_3")}
    for b in range(B):
        idx = rng.choice(n_pool, size=K, replace=False)
        labels = pool_arr[idx]
        counts = np.array([np.sum(labels == r) for r in eligible_regions],
                          dtype=float)
        s = concentration_scalars(counts)
        for k in rec:
            rec[k][b] = s[k]
    return rec


def aggregate_bootstrap(observed: dict, null_arrays: dict) -> dict:
    out = dict(observed)
    G_b = null_arrays["G"]
    s1_b = null_arrays["s_1"]
    Neff_b = null_arrays["N_eff"]
    out["p_G_null"] = float((G_b <= observed["G"]).mean())
    out["p_s1_null"] = float((s1_b >= observed["s_1"]).mean())
    out["N_eff_null_lo"] = float(np.quantile(Neff_b, 0.025))
    out["N_eff_null_hi"] = float(np.quantile(Neff_b, 0.975))
    out["s_1_null_lo"] = float(np.quantile(s1_b, 0.025))
    out["s_1_null_hi"] = float(np.quantile(s1_b, 0.975))
    return out


def coarse_table(traces: pd.DataFrame, pool: pd.DataFrame,
                 band: str, regime: str, level: str) -> pd.DataFrame:
    sub_traces = restrict(traces[traces.band == band], regime)
    sub_pool = restrict(pool, regime)
    K = len(sub_traces)
    N_total = len(sub_pool)
    if K == 0:
        return pd.DataFrame()
    col = level
    base_per = sub_pool.groupby(col).size().rename("n_contacts").reset_index()
    per = (sub_traces.groupby(col)
           .agg(n_trace=("leaf_id", "count"),
                n_contributing_patients=("patient", "nunique"))
           .reset_index())
    m = per.merge(base_per, on=col, how="left")
    m["trace_rate"] = m.n_trace / m.n_contacts
    m["baseline_rate"] = K / N_total
    m["enrichment"] = m.trace_rate / m.baseline_rate
    m["p_hyper"] = m.apply(
        lambda r: float(hypergeom.sf(r.n_trace - 1, N_total, K, int(r.n_contacts))),
        axis=1)
    m["band"] = band
    m["group_level"] = level
    m["group"] = m[col]
    m["sensitivity_regime"] = regime
    return m[["band", "group_level", "group", "sensitivity_regime",
              "n_trace", "n_contacts", "trace_rate", "baseline_rate",
              "enrichment", "p_hyper", "n_contributing_patients"]]


# --------------------------------------------------------------------------
# Verdict
# --------------------------------------------------------------------------
def verdict(c: pd.Series, p: pd.Series, region_df: pd.DataFrame,
            band: str, regime: str) -> str:
    full = region_df[(region_df.band == band)
                     & (region_df.sensitivity_regime == regime)]
    n_elig = int(full.eligible.sum())
    G = c["G"]
    Neff = c["N_eff"]
    s_1 = c["s_1"]
    p_G = c["p_G_null"]
    p_s1 = c["p_s1_null"]
    n_bonf = int(((full.p_bonferroni < 0.05) & full.eligible).sum())
    median_pat = p["median_n_contributing"]
    frac_eq1 = p["frac_eq1"]

    if n_bonf > 0:
        return "single-region"
    distributed = (G <= 0.4 and Neff >= 0.6 * n_elig and s_1 <= 0.20
                   and p_G <= 0.05
                   and (median_pat >= 3 if not np.isnan(median_pat) else False))
    if distributed:
        return "distributed-coverage-consistent"
    fragmented = ((not np.isnan(median_pat)) and median_pat <= 1.0
                  and (not np.isnan(frac_eq1)) and frac_eq1 >= 0.7)
    if fragmented:
        return "patient-fragmented"
    return "inconclusive"


def lobar_verdict(coarse_df: pd.DataFrame, band: str, regime: str) -> str | None:
    sub = coarse_df[(coarse_df.band == band)
                    & (coarse_df.group_level == "lobe")
                    & (coarse_df.sensitivity_regime == regime)
                    & (coarse_df.p_bonferroni < 0.05)]
    if sub.empty:
        return None
    return ", ".join(sub["group"].tolist())


def hemi_verdict(coarse_df: pd.DataFrame, band: str, regime: str) -> str | None:
    sub = coarse_df[(coarse_df.band == band)
                    & (coarse_df.group_level == "hemisphere")
                    & (coarse_df.sensitivity_regime == regime)
                    & (coarse_df.p_bonferroni < 0.05)]
    if sub.empty:
        return None
    return ", ".join(sub["group"].tolist())


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------
def _short(region: str) -> str:
    return region.replace("ctx-lh-", "L:").replace("ctx-rh-", "R:")


def make_violin_figure(region_df, conc_df, pool, traces, rng,
                       focal_anatomy: pd.DataFrame):
    full = region_df[region_df.sensitivity_regime == "full"]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), sharey=True)
    for ax, band in zip(axes, TARGET_BANDS):
        sub = full[(full.band == band) & full.eligible].copy()
        if len(sub) == 0:
            ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} — no eligible regions")
            continue
        K = int(sub["K"].iloc[0])
        N_total = int(sub["N_total"].iloc[0])
        baseline = sub["baseline_rate"].iloc[0]
        sub_pool = restrict(pool, "full")
        elig = sub["region"].tolist()
        contacts = sub.set_index("region").loc[elig, "n_contacts"]
        B_fig = 300
        pool_arr = sub_pool["region"].to_numpy()
        n_pool = len(pool_arr)
        null_rates = np.zeros((B_fig, len(elig)))
        for b in range(B_fig):
            idx = rng.choice(n_pool, size=K, replace=False)
            labs = pool_arr[idx]
            for j, r in enumerate(elig):
                nt = int(np.sum(labs == r))
                null_rates[b, j] = nt / contacts[r]
        x = np.arange(len(elig))
        parts = ax.violinplot(null_rates, positions=x, widths=0.7,
                              showmeans=False, showextrema=False, showmedians=False)
        for pc in parts["bodies"]:
            pc.set_facecolor("#cfd8e3")
            pc.set_edgecolor("none")
            pc.set_alpha(0.7)
        obs_rate = sub["trace_rate"].to_numpy()
        is_bonf = (sub["p_bonferroni"] < 0.05).to_numpy()
        ax.scatter(x[~is_bonf], obs_rate[~is_bonf],
                   s=42, color="#1f3d6e", zorder=3, edgecolor="white", lw=0.6,
                   label="observed")
        ax.scatter(x[is_bonf], obs_rate[is_bonf],
                   s=80, marker="*", color="#c0392b", zorder=4,
                   edgecolor="white", lw=0.6,
                   label=f"$p_{{\\rm Bonf}}<0.05$ (m={BONFERRONI_M})")
        ax.axhline(baseline, color="#888888", lw=0.7, ls=":",
                   label=f"baseline {baseline:.3f}")
        # Mark §5.2-cited focal-leaf regions if they fall on an eligible region
        focal_b = focal_anatomy[focal_anatomy.band == band]
        for _, fl in focal_b.iterrows():
            if fl["region"] in elig:
                jj = elig.index(fl["region"])
                ax.scatter([jj], [obs_rate[jj]], s=140, marker="o",
                           facecolor="none", edgecolor="#e89c40", lw=2.0,
                           zorder=5,
                           label=("§5.2 focal-leaf region"
                                  if jj == elig.index(focal_b.iloc[0]["region"])
                                  and fl.equals(focal_b.iloc[0]) else None))
        ax.set_xticks(x)
        ax.set_xticklabels([_short(r) for r in elig], rotation=55,
                           ha="right", fontsize=8)
        ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]}  (K={K}, m={len(elig)})")
        if ax is axes[0]:
            ax.set_ylabel("topology-trace leaf rate")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "topology_trace_violins.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIG / 'topology_trace_violins.pdf'}")


def make_lobar_heatmap(coarse_df: pd.DataFrame):
    full = coarse_df[(coarse_df.sensitivity_regime == "full")
                     & (coarse_df.group_level == "lobe")
                     & coarse_df.band.isin(TARGET_BANDS)]
    if full.empty:
        return
    bands = TARGET_BANDS
    lobes = sorted(full["group"].unique(),
                   key=lambda x: ("frontal", "parietal", "temporal",
                                  "occipital", "limbic", "insular",
                                  "other").index(x)
                   if x in ("frontal", "parietal", "temporal", "occipital",
                            "limbic", "insular", "other") else 99)
    M = np.full((len(bands), len(lobes)), np.nan)
    sig = np.zeros_like(M, dtype=bool)
    n_trace = np.zeros_like(M, dtype=int)
    n_contacts = np.zeros_like(M, dtype=int)
    for i, b in enumerate(bands):
        for j, lobe in enumerate(lobes):
            row = full[(full.band == b) & (full.group == lobe)]
            if not row.empty:
                M[i, j] = float(row["enrichment"].iloc[0])
                sig[i, j] = bool(row["p_bonferroni"].iloc[0] < 0.05)
                n_trace[i, j] = int(row["n_trace"].iloc[0])
                n_contacts[i, j] = int(row["n_contacts"].iloc[0])
    fig, ax = plt.subplots(figsize=(0.9 * len(lobes) + 2.5,
                                    0.9 * len(bands) + 2.0))
    vmax = max(2.0, np.nanmax(M)) if np.isfinite(np.nanmax(M)) else 2.0
    im = ax.imshow(M, cmap="RdBu_r", vmin=0, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(lobes)))
    ax.set_xticklabels(lobes, rotation=30, ha="right")
    ax.set_yticks(range(len(bands)))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in bands])
    for i in range(len(bands)):
        for j in range(len(lobes)):
            if np.isnan(M[i, j]):
                continue
            txt = f"{M[i, j]:.2f}\n{n_trace[i, j]}/{n_contacts[i, j]}"
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=8, color="#222")
            if sig[i, j]:
                ax.text(j + 0.30, i - 0.30, "*", ha="left", va="top",
                        fontsize=14, fontweight="bold", color="#c0392b")
    fig.colorbar(im, ax=ax, label="enrichment ratio")
    fig.tight_layout()
    fig.savefig(FIG / "topology_lobar_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIG / 'topology_lobar_heatmap.pdf'}")


# --------------------------------------------------------------------------
# Focal-leaf anatomy
# --------------------------------------------------------------------------
def focal_leaf_anatomy(anatomy_by_pat: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for pat, band, leaf_id in FOCAL_CITED:
        anat = anatomy_by_pat[pat]
        if leaf_id < 0 or leaf_id >= len(anat):
            rows.append({"patient": pat, "band": band, "leaf_id": leaf_id,
                         "channel": "<oob>", "region": "<oob>",
                         "lobe": "<oob>", "hemisphere": "<oob>",
                         "x": np.nan, "y": np.nan, "z": np.nan})
            continue
        r = anat.iloc[leaf_id]
        rows.append({
            "patient": pat, "band": band, "leaf_id": leaf_id,
            "channel": r["channel"], "region": r["region"],
            "lobe": lobe_of(r["region"]),
            "hemisphere": hemisphere_of(r["region"]),
            "x": r["x"], "y": r["y"], "z": r["z"],
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# README
# --------------------------------------------------------------------------
def write_readme(catalog: pd.DataFrame, region_df: pd.DataFrame,
                 conc_df: pd.DataFrame, coarse_df: pd.DataFrame,
                 pat_df: pd.DataFrame, focal_df: pd.DataFrame,
                 n_hemi: int, n_lobe: int) -> None:
    full_region = region_df[region_df.sensitivity_regime == "full"]
    lines = [
        "---",
        "date: 2026-05-08",
        "era: COHORT_N10 / IMCOH_ABS",
        "status: current",
        "type: audit-report",
        "scope: section_5_2_kc_lambda0_topology_anatomy",
        "predicate-source: .agents/guides/task-persistence-investigation/2026-05-07_kc-leaf-topological-memory.md",
        "---",
        "",
        "# KC λ=0 (topology) trace — anatomical-distribution audit",
        "",
        "**Head.** Verdict per band ∈ {β (primary), α (control)} on whether "
        "the §5.2 KC λ=0 per-leaf topology-trace flag admits an anatomical "
        "address (single-region / lobar / hemispheric) or spreads across "
        "cortex. Predicate is **reused verbatim** from §5.2 / `audit_53` "
        "(no redefinition); the within-patient null is the anti-trace "
        "predicate (rspre↔rspost agreement, task differing). γ_low is "
        "**excluded** as it is the KC heights band, not topology — Option B "
        "(parallel λ=1 audit on a heights-aligned per-leaf predicate) is "
        "documented as a deferred extension.",
        "",
        "## Catalog provenance and counters",
        "",
        f"- Predicate: §5.2 / `audit_53` (δ=1, min_m=2, max_mrca_frac=0.50). "
        f"Topology-trace flag = `f_i > p{int(NULL_QUANTILE*100)}(f_i^anti) "
        f"AND f_i > 0`.",
        f"- Cortical pool: drop Wm/Unk → "
        f"N_total = {int(conc_df['N_total'].iloc[0]) if 'N_total' in conc_df.columns else len(catalog)} "
        f"(see region table for full-regime counts).",
        f"- Eligibility filter (per band): `n_contacts ≥ 5 AND "
        f"n_contributing_patients ≥ 2`.",
        f"- Bonferroni denominators: m_single = {BONFERRONI_M} (per spec, "
        f"matches §5.5); m_hemi = {n_hemi}; m_lobe = {n_lobe} (auto-derived "
        f"from groups present × 6 bands, matching audit_56).",
        "",
        "**K(band) verification (full regime; topology-trace leaves):**",
        "",
    ]
    for band in TARGET_BANDS:
        sub = full_region[full_region.band == band]
        K = int(sub["K"].iloc[0]) if not sub.empty else 0
        lines.append(f"- K({band}) = {K}")
    lines += [
        "",
        "**Per-patient flag count (full regime, post-Wm/Unk filter):**",
        "",
        "| patient | α | β |",
        "|---|---|---|",
    ]
    for pat in PATIENTS:
        a = int(((catalog.patient == pat) & (catalog.band == "alpha")).sum())
        b = int(((catalog.patient == pat) & (catalog.band == "beta")).sum())
        lines.append(f"| {pat} | {a} | {b} |")
    lines += [
        "",
        "**Disclosure (threshold asymmetry).** The within-patient null is the "
        "anti-trace predicate (rspre↔rspost agreement, task differing). When "
        "anti-trace pairs dominate trace pairs at a (patient, band) cell "
        "(e.g. Pat_14 at β: 360 anti vs 82 trace pairs), the 95th-percentile "
        "threshold of the f_i^anti distribution is so high that no leaf's "
        "f_i exceeds it — that patient contributes zero to K. This is a "
        "literal reading of the spec ('one-sided upper tail of f_i^anti'); "
        "we did not switch to the rspre split-half fallback. The implication "
        "is that K under-counts patients with anti-direction-dominated phase "
        "profiles; verdicts below are conditional on this asymmetry.",
        "",
        "## Per-band verdict (target bands, full regime)",
        "",
        "| band | verdict | G | N_eff | n_eligible | s_1 | s_3 | p_G_null | p_s1_null | n_Bonf | median n_pat | frac_eq1 |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for band in TARGET_BANDS:
        c = conc_df[(conc_df.band == band)
                    & (conc_df.sensitivity_regime == "full")]
        if c.empty:
            continue
        c = c.iloc[0]
        p = pat_df[(pat_df.band == band)
                   & (pat_df.sensitivity_regime == "full")].iloc[0]
        v = verdict(c, p, region_df, band, "full")
        lines.append(
            f"| {band} | **{v}** | {c['G']:.3f} | {c['N_eff']:.2f} | "
            f"{int(c['n_eligible_regions'])} | {c['s_1']:.3f} | {c['s_3']:.3f} "
            f"| {c['p_G_null']:.3f} | {c['p_s1_null']:.3f} | "
            f"{int(c['n_regions_bonferroni'])} | "
            f"{p['median_n_contributing']:.1f} | "
            f"{p['frac_eq1']:.2f} |"
        )
    lines += [
        "",
        "**Verdict labels (KC-adapted):**",
        "- **single-region**: any region survives Bonferroni at m=48.",
        "- **lobar**: any lobe survives Bonferroni at m_lobe; reported "
        "separately from the verdict column when the single-region test "
        "fails (see hemisphere/lobe table below).",
        "- **hemispheric**: any hemisphere survives Bonferroni at m_hemi; "
        "reported separately.",
        "- **distributed-coverage-consistent**: G ≤ 0.4 AND N_eff ≥ 0.6 × "
        "n_eligible AND s_1 ≤ 0.20 AND p_G_null ≤ 0.05 AND median "
        "n_contributing_patients ≥ 3.",
        "- **patient-fragmented**: median n_contributing_patients ≤ 1 AND "
        "frac_eq1 ≥ 0.7 (the §5.6 β CTM finding pattern).",
        "- **inconclusive**: anything else.",
        "",
        "## Single-region descriptive counts (full regime)",
        "",
        "| band | n_eligible | n_with_trace | n_at_or_above_baseline | n_uncorr_p<0.05 | n_Bonf_m=48 |",
        "|---|---|---|---|---|---|",
    ]
    for band in TARGET_BANDS:
        c = conc_df[(conc_df.band == band)
                    & (conc_df.sensitivity_regime == "full")]
        if c.empty:
            continue
        c = c.iloc[0]
        lines.append(
            f"| {band} | {int(c['n_eligible_regions'])} | "
            f"{int(c['n_regions_with_trace'])} | "
            f"{int(c['n_regions_at_or_above_baseline'])} | "
            f"{int(c['n_regions_uncorrected_p005'])} | "
            f"{int(c['n_regions_bonferroni'])} |"
        )
    # Coarser-grain Bonferroni survivors at full regime
    lines += [
        "",
        "## Coarser-grain Bonferroni survivors (full regime)",
        "",
        "Reported alongside the single-region verdict so the lobar / "
        "hemispheric structure is visible even when a single region survives.",
        "",
        "| band | single-region (m=48) | lobe (m=" + str(n_lobe) + ") | hemisphere (m=" + str(n_hemi) + ") |",
        "|---|---|---|---|",
    ]
    for band in TARGET_BANDS:
        full = region_df[(region_df.band == band)
                         & (region_df.sensitivity_regime == "full")
                         & region_df.eligible]
        sr = full[full.p_bonferroni < 0.05]
        sr_str = (", ".join(f"{r['region']} (×{r['enrichment']:.2f}, "
                            f"p_Bonf={r['p_bonferroni']:.3f})"
                            for _, r in sr.iterrows())
                  if not sr.empty else "—")
        lobe_s = lobar_verdict(coarse_df, band, "full") or "—"
        if lobe_s != "—":
            sub = coarse_df[(coarse_df.band == band)
                            & (coarse_df.group_level == "lobe")
                            & (coarse_df.sensitivity_regime == "full")
                            & (coarse_df.p_bonferroni < 0.05)]
            lobe_s = ", ".join(f"{r['group']} (×{r['enrichment']:.2f}, "
                               f"p_Bonf={r['p_bonferroni']:.3f})"
                               for _, r in sub.iterrows())
        hemi_s = hemi_verdict(coarse_df, band, "full") or "—"
        if hemi_s != "—":
            sub = coarse_df[(coarse_df.band == band)
                            & (coarse_df.group_level == "hemisphere")
                            & (coarse_df.sensitivity_regime == "full")
                            & (coarse_df.p_bonferroni < 0.05)]
            hemi_s = ", ".join(f"{r['group']} (×{r['enrichment']:.2f}, "
                               f"p_Bonf={r['p_bonferroni']:.3f})"
                               for _, r in sub.iterrows())
        lines.append(f"| {band} | {sr_str} | {lobe_s} | {hemi_s} |")
    lines += [
        "",
        "## Sensitivity (verdict per regime)",
        "",
        "| band | full | no_pat03 | pro_cohort |",
        "|---|---|---|---|",
    ]
    for band in TARGET_BANDS:
        cells = []
        for regime in SENSITIVITIES:
            c = conc_df[(conc_df.band == band)
                        & (conc_df.sensitivity_regime == regime)]
            p = pat_df[(pat_df.band == band)
                       & (pat_df.sensitivity_regime == regime)]
            if c.empty or p.empty:
                cells.append("(empty)")
                continue
            cells.append(verdict(c.iloc[0], p.iloc[0], region_df, band, regime))
        lines.append(f"| {band} | {cells[0]} | {cells[1]} | {cells[2]} |")
    lines += [
        "",
        "## Hemisphere + lobe coarse aggregation (full regime)",
        "",
        f"Bonferroni denominators: m_hemisphere = {n_hemi}, m_lobe = {n_lobe}.",
        "",
        "**Hemisphere (band × group):**",
        "",
        "| band | group | n_trace | n_contacts | enrichment | p_hyper | p_Bonf |",
        "|---|---|---|---|---|---|---|",
    ]
    hemi = coarse_df[(coarse_df.group_level == "hemisphere")
                     & (coarse_df.sensitivity_regime == "full")
                     & coarse_df.band.isin(TARGET_BANDS)]
    for _, r in hemi.iterrows():
        lines.append(
            f"| {r['band']} | {r['group']} | {int(r['n_trace'])} | "
            f"{int(r['n_contacts'])} | {r['enrichment']:.2f} | "
            f"{r['p_hyper']:.3f} | {r['p_bonferroni']:.3f} |"
        )
    lines += [
        "",
        "**Lobe (band × group):**",
        "",
        "| band | group | n_trace | n_contacts | enrichment | p_hyper | p_Bonf |",
        "|---|---|---|---|---|---|---|",
    ]
    lobe = coarse_df[(coarse_df.group_level == "lobe")
                     & (coarse_df.sensitivity_regime == "full")
                     & coarse_df.band.isin(TARGET_BANDS)]
    for _, r in lobe.iterrows():
        lines.append(
            f"| {r['band']} | {r['group']} | {int(r['n_trace'])} | "
            f"{int(r['n_contacts'])} | {r['enrichment']:.2f} | "
            f"{r['p_hyper']:.3f} | {r['p_bonferroni']:.3f} |"
        )
    lines += [
        "",
        "## Patient base (full regime)",
        "",
        "| band | median n_contributing | frac ≥ 3 | frac = 1 | n_regions_with_trace |",
        "|---|---|---|---|---|",
    ]
    for band in TARGET_BANDS:
        p = pat_df[(pat_df.band == band)
                   & (pat_df.sensitivity_regime == "full")]
        if p.empty:
            continue
        p = p.iloc[0]
        lines.append(
            f"| {band} | {p['median_n_contributing']:.1f} | "
            f"{p['frac_geq3']:.2f} | {p['frac_eq1']:.2f} | "
            f"{int(p['n_regions_with_trace'])} |"
        )
    lines += [
        "",
        "## §5.2 focal-leaf anatomy",
        "",
        "Independent structural-illustration check. Regions of the focal "
        "leaves cited in §5.2 (Fig. `manuscript_kc_leaf_topology.pdf`):",
        "",
        "| patient | band | leaf | channel | region | lobe | hemisphere |",
        "|---|---|---|---|---|---|---|",
    ]
    for _, r in focal_df.iterrows():
        lines.append(
            f"| {r['patient']} | {r['band']} | {int(r['leaf_id'])} | "
            f"{r['channel']} | {r['region']} | {r['lobe']} | {r['hemisphere']} |"
        )
    # Consistency line
    beta_focals = focal_df[focal_df.band == "beta"]
    if not beta_focals.empty:
        regions = beta_focals["region"].tolist()
        lobes = beta_focals["lobe"].tolist()
        hemis = beta_focals["hemisphere"].tolist()
        same_region = (len(set(regions)) == 1)
        same_lobe = (len(set(lobes)) == 1)
        same_hemi = (len(set(hemis)) == 1)
        lines += [
            "",
            "**β focal-leaf consistency:** "
            f"shared region={same_region} ({sorted(set(regions))}); "
            f"shared lobe={same_lobe} ({sorted(set(lobes))}); "
            f"shared hemisphere={same_hemi} ({sorted(set(hemis))}).",
        ]
    lines += [
        "",
        "## Option B (deferred)",
        "",
        "A parallel λ=1 (heights-aligned) audit at γ_low would require a "
        "new per-leaf predicate (MRCA height agreement instead of MRCA depth "
        "agreement). The §5.2 predicate is depth-only and the existing "
        "audit_53 catalog is depth-only; defining and validating the height "
        "predicate is left as future work.",
        "",
        "## Files",
        "",
        "- `data/kc_topology_anatomy/per_topology_trace_leaf.csv` "
        "(catalog parallel to §5.5 `per_trace_leaf.csv`)",
        "- `data/kc_topology_anatomy/region_count_table.csv`",
        "- `data/kc_topology_anatomy/concentration_scalars.csv`",
        "- `data/kc_topology_anatomy/hemisphere_lobe_enrichment.csv`",
        "- `data/kc_topology_anatomy/patient_base_distribution.csv`",
        "- `data/kc_topology_anatomy/focal_leaf_anatomy.csv`",
        "- `data/kc_topology_anatomy/figures/topology_trace_violins.pdf`",
        "- `data/kc_topology_anatomy/figures/topology_lobar_heatmap.pdf`",
        "- Build script: `scripts/01_compute/audit/audit_57_kc_topology_anatomy.py`",
        "",
    ]
    (OUT / "README.md").write_text("\n".join(lines))
    print(f"wrote {OUT / 'README.md'}")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main() -> None:
    rng = np.random.default_rng(RNG_SEED)

    print("loading per-patient anatomy ...")
    anatomy_by_pat = {p: load_per_leaf_anatomy(p) for p in PATIENTS}

    pool = load_pool()
    print(f"cortical pool: {len(pool)} contacts, "
          f"{pool.region.nunique()} regions, "
          f"{pool.patient.nunique()} patients")

    print("building topology-trace catalog (KC λ=0, target bands) ...")
    catalog = build_catalog(anatomy_by_pat)
    catalog.to_csv(OUT / "per_topology_trace_leaf.csv", index=False)
    print(f"wrote {OUT / 'per_topology_trace_leaf.csv'} "
          f"({len(catalog)} rows)")

    # Per-cell flag counts diagnostic
    print()
    print("topology-trace flag counts per (patient, band):")
    diag = (catalog.groupby(["patient", "band"]).size()
            .rename("n_flagged").reset_index())
    print(diag.pivot(index="patient", columns="band",
                     values="n_flagged").fillna(0).astype(int).to_string())

    region_rows = []
    patient_base_rows = []
    concentration_rows = []
    coarse_rows = []

    for regime in SENSITIVITIES:
        for band in TARGET_BANDS:
            tab = per_region_table(catalog, pool, band, regime)
            if tab.empty:
                continue
            region_rows.append(tab)

            elig = tab[tab.eligible]["region"].tolist()
            n_elig = len(elig)
            n_with_trace = int((tab.n_trace > 0).sum())
            n_at_or_above = int((tab.eligible & (tab.trace_rate >= tab.baseline_rate)).sum())
            n_uncorr_005 = int((tab.eligible & (tab.p_hyper < 0.05)).sum())
            n_bonf = int((tab.eligible & (tab.p_bonferroni < 0.05)).sum())

            elig_counts = (tab.set_index("region").loc[elig, "n_trace"]
                           .to_numpy(dtype=float)
                           if n_elig > 0 else np.array([]))
            obs = concentration_scalars(elig_counts) if n_elig > 0 else dict(
                G=float("nan"), N_eff=float("nan"), s_1=float("nan"),
                s_3=float("nan"))
            K = int(tab["K"].iloc[0])
            N_total = int(tab["N_total"].iloc[0])

            sub_pool = restrict(pool, regime)
            if n_elig > 0 and K > 0:
                null = bootstrap_null(elig, sub_pool, K, B_BOOT, rng)
                agg = aggregate_bootstrap(obs, null)
            else:
                agg = dict(obs)
                for k in ("p_G_null", "p_s1_null",
                          "N_eff_null_lo", "N_eff_null_hi",
                          "s_1_null_lo", "s_1_null_hi"):
                    agg[k] = float("nan")

            concentration_rows.append({
                "band": band, "sensitivity_regime": regime,
                **agg,
                "K_total": K,
                "N_total": N_total,
                "n_eligible_regions": n_elig,
                "n_regions_with_trace": n_with_trace,
                "n_regions_at_or_above_baseline": n_at_or_above,
                "n_regions_uncorrected_p005": n_uncorr_005,
                "n_regions_bonferroni": n_bonf,
            })

            with_trace = tab[tab.n_trace > 0]
            ncontrib = with_trace["n_contributing_patients"].to_numpy()
            patient_base_rows.append({
                "band": band, "sensitivity_regime": regime,
                "median_n_contributing": (float(np.median(ncontrib))
                                          if len(ncontrib) else float("nan")),
                "frac_geq3": (float((ncontrib >= 3).mean())
                              if len(ncontrib) else float("nan")),
                "frac_eq1": (float((ncontrib == 1).mean())
                             if len(ncontrib) else float("nan")),
                "n_regions_with_trace": int(len(ncontrib)),
            })

            for level in ("hemisphere", "lobe"):
                ct = coarse_table(catalog, pool, band, regime, level)
                if not ct.empty:
                    coarse_rows.append(ct)

    if not region_rows:
        print("no eligible (band, region) cells — aborting")
        return

    region_df = pd.concat(region_rows, ignore_index=True)
    region_out = region_df[[
        "band", "region", "sensitivity_regime",
        "n_trace", "n_contacts", "trace_rate", "baseline_rate",
        "enrichment", "p_hyper", "p_bonferroni",
        "n_contributing_patients", "eligible", "K", "N_total",
    ]]
    region_out.to_csv(OUT / "region_count_table.csv", index=False)
    print(f"wrote {OUT / 'region_count_table.csv'} ({len(region_out)} rows)")

    conc_df = pd.DataFrame(concentration_rows)
    conc_df.to_csv(OUT / "concentration_scalars.csv", index=False)
    print(f"wrote {OUT / 'concentration_scalars.csv'}")

    pat_df = pd.DataFrame(patient_base_rows)
    pat_df.to_csv(OUT / "patient_base_distribution.csv", index=False)
    print(f"wrote {OUT / 'patient_base_distribution.csv'}")

    coarse_df = pd.concat(coarse_rows, ignore_index=True)
    n_hemi = coarse_df[coarse_df.group_level == "hemisphere"]["group"].nunique() * 6
    n_lobe = coarse_df[coarse_df.group_level == "lobe"]["group"].nunique() * 6
    coarse_df["m_bonferroni"] = coarse_df.group_level.map(
        {"hemisphere": n_hemi, "lobe": n_lobe})
    coarse_df["p_bonferroni"] = (coarse_df.p_hyper * coarse_df.m_bonferroni
                                 ).clip(upper=1.0)
    coarse_df.to_csv(OUT / "hemisphere_lobe_enrichment.csv", index=False)
    print(f"wrote {OUT / 'hemisphere_lobe_enrichment.csv'} "
          f"(m_hemi={n_hemi}, m_lobe={n_lobe})")

    focal_df = focal_leaf_anatomy(anatomy_by_pat)
    focal_df.to_csv(OUT / "focal_leaf_anatomy.csv", index=False)
    print(f"wrote {OUT / 'focal_leaf_anatomy.csv'}")

    make_violin_figure(region_df, conc_df, pool, catalog, rng, focal_df)
    make_lobar_heatmap(coarse_df)

    write_readme(catalog, region_df, conc_df, coarse_df, pat_df, focal_df,
                 n_hemi, n_lobe)


if __name__ == "__main__":
    main()
