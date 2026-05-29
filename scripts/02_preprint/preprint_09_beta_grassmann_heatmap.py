#!/usr/bin/env python3
"""Grassmann principal-angle heatmap (mode i × cutoff k) — per band.

Default: loops over all six bands and writes one PDF per band to
``data/preprint/figures/<band>/grassmann/fig_<band>_grassmann_heatmap.pdf``.
Run with one or more band names on the command line to restrict, e.g.

    python preprint_09_beta_grassmann_heatmap.py beta
    python preprint_09_beta_grassmann_heatmap.py alpha beta low_gamma

Per-cell principal-angle decomposition of the matched-strength-controlled
Grassmann probe. The cohort-median per-mode persistence

    Δθ_i(k) = θ_i^{rsPre, taskT} − θ_i^{taskT, rsPost}   (rad)

is rendered as a diverging RdBu_r heatmap on the triangular region i ≤ k
(positive Δθ ⇒ rsPost subspace closer to taskT than rsPre ⇒ trace ⇒ red).

Per-`k` cohort-paired matched-strength Wilcoxon `p<0.05` is drawn as a
two-row sig-bar above the heatmap (full FC + epi-X) and as Wilcoxon-tick
strips below. No 7/10 reference line, no contiguous-window pre-frame —
the cohort gate is per-`k` Wilcoxon p<0.05 only.

Right panel: per-mode profiles at six representative cutoffs k spanning
2..112, plasma-coloured small→large; each profile naturally truncates at
i = k.

Inputs
------
data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv
data/audit/grassmann_epi_exclusion/cohort_summary.csv
data/cache/imcoh_lrg/Pat_NN/<band>_<phase>_lrg_imcoh-abs.npz
  (loaded via workflow.lrg.load_lrg_result)

Outputs (PDF only, no PNG sibling; one per band)
------------------------------------------------
data/preprint/figures/<band>/grassmann/fig_<band>_grassmann_heatmap.pdf
data/preprint/cache/grassmann_principal_angles_<band>_full_k.csv  (cache)
"""
from __future__ import annotations

import functools
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.spatial import (
    load_spatial_metadata,
    prepare_spatial_coordinates,
)
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result

ROOT = setup_script_env()
use_lrg_style()

# ---------------------------------------------------------------------------
# Centralised label sizing. The lrg_eegfc mplstyle baseline is tuned for
# compact 2x2 / 1x4 publication figures (font.size=9). This figure is a
# single full-page composite; bump every text element by LABEL_SCALE so
# labels remain legible at the intended print/screen size. Edit only this
# block to retune the whole figure — never hardcode `fontsize=N` below.
# ---------------------------------------------------------------------------
LABEL_SCALE = 2.0
plt.rcParams.update({
    "font.size":             9 * LABEL_SCALE,
    "axes.titlesize":        10 * LABEL_SCALE,
    "axes.labelsize":        9 * LABEL_SCALE,
    "xtick.labelsize":       7 * LABEL_SCALE,
    "ytick.labelsize":       7 * LABEL_SCALE,
    "legend.fontsize":       7 * LABEL_SCALE,
    "legend.title_fontsize": 8 * LABEL_SCALE,
    "figure.titlesize":      10 * LABEL_SCALE,
})


ALL_BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
PHASES = ("rest_pre", "task_test", "rest_post")

GRASS_DIR = ROOT / "data" / "audit" / "grassmann_matched_strength_surrogate"
GRASS_EPIX_DIR = ROOT / "data" / "audit" / "grassmann_epi_exclusion"
CLUSTER_DIR = ROOT / "data" / "audit" / "grassmann_cluster_extent"

# Half-FC + surrogate-eigvec caches drive the rest_pre_A observed
# eigendecomposition and the per-mode surrogate principal-angle compute,
# matching the construction the cluster-mass verdict is built on.
HALVES_FC_CACHE = ROOT / "data" / "cache" / "imcoh_halves_fc"
SURROGATE_LRG_CACHE = ROOT / "data" / "cache" / "matched_strength_surrogate_lrg"
SURROGATE_SEED_TAG = "R200_swap20_seed20260511_imcoh_abs"
N_SURROGATES_PA = 200

GRASS_PER_PAT_CSV = "per_patient_per_band_per_k.csv"

CACHE_DIR = ROOT / "data" / "preprint" / "cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# Full k axis: matched-strength cohort summary covers k=2..112.
K_MIN, K_MAX = 2, 112
I_MAX = K_MAX

CLR_FULL = "#1d1d1d"
CLR_EPIX = "#7f3b2c"
CLR_SIG = "#1f7a1f"
CLR_SURR = "#7d7d7d"           # T_G surrogate cohort-median line
CLR_SURR_FILL = "#dcdcdc"      # T_G surrogate IQR envelope fill
CLR_TG_OBS = "#1f3d6e"         # T_G observed (full FC)
CLR_TG_EPIX = CLR_EPIX         # T_G observed (epi-X)

# Cluster-extent verdict colour scheme (locked 2026-05-19, VERDICT_LEDGER.md).
CLR_VERDICT_STRONG = "#0d5c1c"     # deep green: cluster_p < 0.01
CLR_VERDICT_WEAK = "#a48b22"        # ochre:     0.01 ≤ cluster_p < 0.05
CLR_VERDICT_NONE = "#9a9a9a"        # grey:      cluster_p ≥ 0.05

VERDICT_COLOR = {
    "strong": CLR_VERDICT_STRONG,
    "weak": CLR_VERDICT_WEAK,
    "no_trace": CLR_VERDICT_NONE,
}

K_PROFILES = [10, 30, 50, 70, 90, 110]


# ---------------------------------------------------------------------------
# Compute helpers
# ---------------------------------------------------------------------------
def _topk_basis(eigvecs: np.ndarray, k: int) -> np.ndarray:
    return np.ascontiguousarray(eigvecs[:, 1 : k + 1])


def _principal_angles(V_a: np.ndarray, V_b: np.ndarray) -> np.ndarray:
    sigma = np.linalg.svd(V_a.T @ V_b, compute_uv=False)
    sigma = np.clip(sigma, 0.0, 1.0)
    return np.arccos(sigma)


def _load_obs_eigvecs_rest_pre_A(pat: str, band: str) -> np.ndarray:
    """Eigenvectors of the combinatorial Laplacian at rest_pre_A.

    Loads the cached half-FC matrix written by audit_66 / audit_63,
    builds L̂ = D̂ − W, and returns the eigvecs (ascending eigenvalue
    order, same convention as `load_lrg_result`).
    """
    path = HALVES_FC_CACHE / pat / f"{band}_rest_pre_A_imcoh_abs.npy"
    W = np.asarray(np.load(path), dtype=float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    W = 0.5 * (W + W.T)
    deg = W.sum(axis=1)
    L = np.diag(deg) - W
    _eigvals, eigvecs = np.linalg.eigh(L)
    return eigvecs


def compute_principal_angles_full_k(band: str) -> pd.DataFrame:
    """Cohort principal-angle Δθ_i(k) over k=2..112 for the given band.

    Convention: Δθ_i = θ_i^{rsPre_A,taskT} − θ_i^{taskT,rsPost}. Positive
    ⇒ rsPost subspace sits closer to taskT than rsPre_A does along mode
    i ⇒ trace direction.

    Phase choice: rest_pre_A (first half of rest_pre) matches the
    surrogate construction in audit_66 and the cluster-mass verdict in
    audit_70. The matched-strength surrogate ensemble at rest_pre is
    only cached for rest_pre_A / rest_pre_B, so for the residual
    heatmap to be self-consistent we anchor observed eigvecs on
    rest_pre_A as well.
    """
    cache = CACHE_DIR / f"grassmann_principal_angles_obs_rest_pre_A_{band}_full_k.csv"
    if cache.exists():
        return pd.read_csv(cache)

    per_pat = []
    for pat in COHORT:
        EV = {
            "rest_pre_A": _load_obs_eigvecs_rest_pre_A(pat, band),
            "task_test": load_lrg_result(
                pat, "task_test", band, fc_method="imcoh_abs").eigenvectors,
            "rest_post": load_lrg_result(
                pat, "rest_post", band, fc_method="imcoh_abs").eigenvectors,
        }
        N = EV["rest_pre_A"].shape[0]
        k_hi = min(K_MAX, N - 1)
        for k in range(K_MIN, k_hi + 1):
            Vpre = _topk_basis(EV["rest_pre_A"], k)
            Vtt = _topk_basis(EV["task_test"], k)
            Vpost = _topk_basis(EV["rest_post"], k)
            th_pre_tt = _principal_angles(Vpre, Vtt)
            th_tt_post = _principal_angles(Vtt, Vpost)
            for i, (a, b) in enumerate(zip(th_pre_tt, th_tt_post), start=1):
                per_pat.append(dict(
                    patient=pat, k=k, mode_index=i,
                    delta_theta=float(a - b),
                ))
    pp = pd.DataFrame(per_pat)
    rows = []
    for (k, i), g in pp.groupby(["k", "mode_index"]):
        vals = g["delta_theta"].values
        rows.append(dict(
            k=int(k),
            mode_index=int(i),
            n_patients=int(g["patient"].nunique()),
            delta_theta_median_radians=float(np.median(vals)),
            delta_theta_iqr_radians=float(
                np.subtract(*np.percentile(vals, [75, 25]))),
        ))
    cohort = pd.DataFrame(rows).sort_values(["k", "mode_index"])
    cohort.to_csv(cache, index=False)
    print(f"  [{band}] cached observed principal-angle decomposition → {cache.name}")
    return cohort


def compute_surrogate_principal_angles_cohort(band: str) -> pd.DataFrame:
    """Cohort-median per-mode surrogate Δθ_i(k) under matched-strength.

    For each patient, load the cached matched-strength surrogate
    eigvec ensembles at rest_pre_A, task_test, rest_post (R=200 each).
    For each surrogate replicate r, compute principal-angle vectors
    θ_i^{rsPre_A^r, taskT^r} and θ_i^{taskT^r, rsPost^r}, take their
    difference Δθ_i^r(k). Per-patient summary = median across r.
    Cohort summary = median across patients (matches the construction
    of the observed cohort-median in compute_principal_angles_full_k).
    """
    cache = CACHE_DIR / f"grassmann_principal_angles_surr_cohort_{band}_full_k.csv"
    if cache.exists():
        return pd.read_csv(cache)

    per_pat_med = []  # one row per (patient, k, mode_index)
    for pat in COHORT:
        E = {}
        for ph_key, ph_file in (("rest_pre_A", "rest_pre_A"),
                                  ("task_test", "task_test"),
                                  ("rest_post", "rest_post")):
            path = (SURROGATE_LRG_CACHE / pat
                    / f"{band}_{ph_file}_{SURROGATE_SEED_TAG}.npz")
            E[ph_key] = np.load(path)["eigvecs"]  # (R, N, N)

        R = min(E["rest_pre_A"].shape[0], E["task_test"].shape[0],
                E["rest_post"].shape[0], N_SURROGATES_PA)
        N = E["rest_pre_A"].shape[1]
        k_hi = min(K_MAX, N - 1)
        # Storage: (R, mode_index_max=k_hi, n_k) — NaN where mode_index > k
        delta_rkm = np.full((R, k_hi, k_hi - K_MIN + 1), np.nan)
        for r in range(R):
            Epre = E["rest_pre_A"][r]
            Ett = E["task_test"][r]
            Epost = E["rest_post"][r]
            for k_idx, k in enumerate(range(K_MIN, k_hi + 1)):
                Vpre = _topk_basis(Epre, k)
                Vtt = _topk_basis(Ett, k)
                Vpost = _topk_basis(Epost, k)
                th_pre_tt = _principal_angles(Vpre, Vtt)
                th_tt_post = _principal_angles(Vtt, Vpost)
                delta_rkm[r, :k, k_idx] = th_pre_tt - th_tt_post

        med_per_pat = np.nanmedian(delta_rkm, axis=0)  # (k_hi, n_k)
        for k_idx, k in enumerate(range(K_MIN, k_hi + 1)):
            for i in range(1, k + 1):
                per_pat_med.append(dict(
                    patient=pat, k=int(k), mode_index=int(i),
                    surr_delta_theta_median=float(med_per_pat[i - 1, k_idx]),
                ))
        print(f"  [{band}] surrogate principal angles done: {pat}")

    pp = pd.DataFrame(per_pat_med)
    rows = []
    for (k, i), g in pp.groupby(["k", "mode_index"]):
        rows.append(dict(
            k=int(k),
            mode_index=int(i),
            surr_delta_theta_cohort_median=float(np.nanmedian(
                g["surr_delta_theta_median"].values)),
        ))
    cohort = pd.DataFrame(rows).sort_values(["k", "mode_index"])
    cohort.to_csv(cache, index=False)
    print(f"  [{band}] cached surrogate principal-angle cohort median → {cache.name}")
    return cohort


def load_cluster_mass_null(band: str
                            ) -> tuple[np.ndarray, float, float, float]:
    """Return (null_array, obs_cluster_mass_neglog10p, p_mass, p_LR).

    Reads audit_70 `null_distribution.csv` for the per-band R=200 null
    cluster-mass values, and `cohort_summary.csv` for the observed
    cluster-mass `−log10 p_k` aggregate + the cluster-mass p-value
    (gate quantity) + longest-run p (companion).
    """
    nd = pd.read_csv(CLUSTER_DIR / "null_distribution.csv")
    cs = pd.read_csv(CLUSTER_DIR / "cohort_summary.csv")
    null_vals = nd[nd.band == band]["cluster_mass"].values.astype(float)
    row = cs[cs.band == band].iloc[0]
    return (null_vals,
            float(row.obs_cluster_mass_neglog10p),
            float(row.cluster_p_cluster_mass),
            float(row.cluster_p_longest_run))


def build_heatmap(pa_obs: pd.DataFrame, pa_surr: pd.DataFrame
                   ) -> tuple[np.ndarray, np.ndarray]:
    """Residual heatmap: obs cohort median − surrogate cohort median.

    H[i, k] = Δθ_i(k)_obs_cohort_median − Δθ_i(k)_surr_cohort_median.
    Trace direction = positive (rsPost closer to taskT than rsPre_A in
    excess of matched-strength expectation). α-vs-β discrimination:
    α residual ≈ 0 (rotation explained by strength), β residual > 0
    in the k ≈ 20..60 manuscript window (genuine reorganization).
    """
    merged = pa_obs.merge(pa_surr, on=["k", "mode_index"], how="left")
    ks = list(range(K_MIN, K_MAX + 1))
    H = np.full((I_MAX, len(ks)), np.nan)
    for _, row in merged.iterrows():
        k = int(row.k)
        i = int(row.mode_index) - 1
        if k < K_MIN or k > K_MAX or i >= I_MAX or i < 0:
            continue
        col = ks.index(k)
        obs = float(row.delta_theta_median_radians)
        surr = (float(row.surr_delta_theta_cohort_median)
                if pd.notna(row.surr_delta_theta_cohort_median) else 0.0)
        H[i, col] = obs - surr
    return H, np.array(ks)


def cohort_strip(ms: pd.DataFrame, band: str
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sub = ms[(ms.band == band) & (ms.k >= K_MIN) & (ms.k <= K_MAX)]
    sub = sub.sort_values("k")
    # Sign-convention adapter: post-2026-05-26 the matched-strength CSV
    # uses "above_own_surrogate" (trace = positive); the epi-X CSV still
    # uses "below_own_surrogate". Accept either to remain compatible
    # while the two caches are aligned.
    if "n_patients_below_own_surrogate" in sub.columns:
        n_col = sub.n_patients_below_own_surrogate.values.astype(float)
    else:
        n_col = sub.n_patients_above_own_surrogate.values.astype(float)
    return (sub.k.values.astype(int),
            n_col / 10.0,
            sub.paired_wilcoxon_p.values.astype(float))


def longest_run_span(ks: np.ndarray, p_vals: np.ndarray
                     ) -> tuple[int, int, int]:
    """Return (length, k_start, k_end) of the longest contiguous p<0.05 run."""
    sig = (p_vals < 0.05).astype(int)
    best_len = 0
    best_s = -1
    best_e = -1
    cur_len = 0
    cur_s = -1
    for idx, s in enumerate(sig):
        if s:
            if cur_len == 0:
                cur_s = idx
            cur_len += 1
            if cur_len > best_len:
                best_len = cur_len
                best_s = cur_s
                best_e = idx
        else:
            cur_len = 0
    if best_len == 0:
        return 0, -1, -1
    return best_len, int(ks[best_s]), int(ks[best_e])


def cohort_surrogate_envelope(per_pat: pd.DataFrame, band: str
                               ) -> pd.DataFrame:
    """Per-k cohort-median of per-patient surrogate quantiles for T_G.

    Returns columns: k, surr_p25, surr_p50, surr_p75, in the Methods
    convention (positive = trace).

    Convention auto-detect: post-2026-05-26 the matched-strength
    per-patient CSV writes `obs_p_one_sided_upper` (positive=trace);
    the legacy epi-X CSV still writes `obs_p_one_sided_lower`
    (negative=trace). When the legacy convention is detected the
    surrogate percentiles are sign-flipped and p25/p75 swapped (since
    percentile order reverses under a sign flip). Once the epi-X
    cache is regenerated under the new convention this branch becomes
    a no-op and can be deleted.
    """
    sub = per_pat[per_pat.band == band]
    legacy = "obs_p_one_sided_lower" in sub.columns
    rows = []
    for k, g in sub.groupby("k"):
        if legacy:
            rows.append(dict(
                k=int(k),
                surr_p25=float(-np.nanmedian(g.surr_T_G_p75.values)),
                surr_p50=float(-np.nanmedian(g.surr_T_G_p50.values)),
                surr_p75=float(-np.nanmedian(g.surr_T_G_p25.values)),
            ))
        else:
            rows.append(dict(
                k=int(k),
                surr_p25=float(np.nanmedian(g.surr_T_G_p25.values)),
                surr_p50=float(np.nanmedian(g.surr_T_G_p50.values)),
                surr_p75=float(np.nanmedian(g.surr_T_G_p75.values)),
            ))
    return pd.DataFrame(rows).sort_values("k").reset_index(drop=True)


def compute_bootstrap_cohort_p(band: str, B: int = 2000,
                                seed: int = 20260528,
                                source: str = "full") -> pd.DataFrame:
    """Patient-bootstrap one-sided cohort p-value per k (trace direction).

    Methodology (parallel to paired Wilcoxon, robust by construction):
        d_i(k) = obs_T_G_i(k) − surr_T_G_p50_i(k)           per patient
        S(k)   = median over patients of d_i(k)              cohort statistic
        Bootstrap: B = 2000 resamples of 10 patients with replacement,
                   recompute S*(k) for each resample.
        p_boot(k) = mean over resamples of [ S*(k) ≤ 0 ].

    A signal driven by 1–2 patients widens the bootstrap distribution
    of S*(k) (∼37 % of resamples don't contain any given patient by
    chance), so the cohort claim self-protects against single-patient
    leverage — no explicit count threshold needed, no LOO sensitivity
    required as a downstream sanity check.

    `source ∈ {"full", "epix"}` picks matched-strength full FC vs the
    epileptic-zone-exclusion variant. Sign convention auto-detected via
    column name (`obs_p_one_sided_upper` = new positive-trace; legacy
    `obs_p_one_sided_lower` flips the sign of d).
    """
    if source == "full":
        csv_path = GRASS_DIR / GRASS_PER_PAT_CSV
    elif source == "epix":
        csv_path = GRASS_EPIX_DIR / GRASS_PER_PAT_CSV
    else:
        raise ValueError(f"unknown source: {source!r}")

    cache = (CACHE_DIR
             / f"grassmann_bootstrap_p_{band}_{source}_B{B}_seed{seed}.csv")
    if cache.exists():
        return pd.read_csv(cache)

    per_pat = pd.read_csv(csv_path)
    sub = per_pat[per_pat.band == band].copy()
    legacy = "obs_p_one_sided_lower" in sub.columns
    sign = -1.0 if legacy else 1.0
    sub["d"] = sign * (sub["obs_T_G"].values - sub["surr_T_G_p50"].values)

    pivot = sub.pivot(index="k", columns="patient", values="d").sort_index()
    ks = pivot.index.values.astype(int)
    D = pivot.values    # (n_k, n_patients), NaN where patient/k missing

    rng = np.random.default_rng(seed)
    p_boot = np.full(len(ks), np.nan)
    S_obs = np.full(len(ks), np.nan)
    S_lo = np.full(len(ks), np.nan)
    S_hi = np.full(len(ks), np.nan)

    for ki in range(len(ks)):
        d_k = D[ki, :]
        valid = ~np.isnan(d_k)
        d_valid = d_k[valid]
        n = len(d_valid)
        if n < 3:
            continue
        idx = rng.integers(0, n, size=(B, n))
        S_star = np.median(d_valid[idx], axis=1)
        p_boot[ki] = float(np.mean(S_star <= 0))
        S_obs[ki] = float(np.median(d_valid))
        S_lo[ki] = float(np.percentile(S_star, 2.5))
        S_hi[ki] = float(np.percentile(S_star, 97.5))

    out = pd.DataFrame({
        "k": ks,
        "S_obs": S_obs,
        "S_boot_2_5": S_lo,
        "S_boot_97_5": S_hi,
        "p_boot": p_boot,
    })
    out.to_csv(cache, index=False)
    print(f"  [{band}/{source}] bootstrap p_boot(k) cached → {cache.name}")
    return out


def cohort_strip_bootstrap(band: str, source: str = "full", B: int = 2000
                            ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bootstrap analog of `cohort_strip` — returns (k, S_obs, p_boot)."""
    df = compute_bootstrap_cohort_p(band, B=B, source=source)
    df = df[(df.k >= K_MIN) & (df.k <= K_MAX)].sort_values("k")
    return (df.k.values.astype(int),
            df.S_obs.values.astype(float),
            df.p_boot.values.astype(float))


def cluster_verdict_row(band: str) -> dict:
    """Read the locked cluster-extent verdict for *band* from audit_70."""
    df = pd.read_csv(CLUSTER_DIR / "cohort_summary.csv")
    row = df[df.band == band]
    if row.empty:
        raise RuntimeError(f"no cluster-extent row for band={band}")
    r = row.iloc[0]
    return dict(
        verdict=str(r.verdict_cluster_extent),
        cluster_p_LR=float(r.cluster_p_longest_run),
        cluster_p_mass=float(r.cluster_p_cluster_mass),
        obs_LR=int(r.obs_longest_run),
        null_mean=float(r.null_mean_LR),
        null_p95=float(r.null_p95_LR),
        null_max=float(r.null_max_LR),
    )


# ---------------------------------------------------------------------------
# Cohort glass-brain inset helpers
# ---------------------------------------------------------------------------
def _primary_region(dk_str) -> str | None:
    """First (highest-weight) DK region in the implant string.

    Mirror of the parser in `spatial._estimate_mni_transform` /
    `preprint_11_beta_anatomy_brain._primary_region`. Returns the raw
    region string (e.g. `"Hip"`, `"ctx-lh-fusiform"`) as it appears in
    the implant CSV and in the audit anatomy `cohort_summary.csv`
    `region` column — no normalisation.
    """
    if pd.isna(dk_str):
        return None
    parts = str(dk_str).strip().split(",")
    if not parts:
        return None
    return parts[0].strip().strip('"').strip()


def _a3_trace_regions(band: str) -> set[str]:
    """A3-passing DK regions for `band` (cluster-extent preferred).

    Reads `data/audit/anatomy_<band>_grassmann_clusterext/cohort_summary.csv`
    when present, falling back to `anatomy_<band>_grassmann/cohort_summary.csv`.
    Returns the set of region strings with `passes_A3 == True`.
    """
    base = ROOT / "data" / "audit"
    cand_clusterext = base / f"anatomy_{band}_grassmann_clusterext" / "cohort_summary.csv"
    cand_plain = base / f"anatomy_{band}_grassmann" / "cohort_summary.csv"
    if cand_clusterext.exists():
        csv = cand_clusterext
    elif cand_plain.exists():
        csv = cand_plain
    else:
        print(f"  [{band}] no anatomy CSV — A3 set empty")
        return set()
    df = pd.read_csv(csv)
    if "passes_A3" not in df.columns:
        return set()
    a3 = set(df.loc[df["passes_A3"] == True, "region"].astype(str).tolist())
    return a3


# Locked manuscript palette for β (matches preprint_11_beta_anatomy_brain
# so the same region carries the same colour across both figures).
BETA_GRASS_REGIONS = [
    "Hip",
    "ctx-rh-rostralmiddlefrontal",
    "ctx-lh-lateralorbitofrontal",
    "ctx-rh-medialorbitofrontal",
    "ctx-lh-middletemporal",
    "ctx-lh-superiortemporal",
    "ctx-lh-insula",
]
PALETTE_BETA_GRASS = [
    "#cb4335",  # Hippocampus            — brick red
    "#f1c40f",  # R rostral mid frontal  — golden yellow
    "#e84393",  # L lateral OFC          — magenta
    "#16a085",  # R medial OFC           — jade
    "#5b2c6f",  # L middle temporal      — plum
    "#239b56",  # L superior temporal    — emerald
    "#85929e",  # L insula               — cool grey
]
NONSIG_COLOR = "#c8c8c8"


def _band_a3_regions_palette(band: str) -> tuple[list[str], list[str]]:
    """Ordered list of A3-passing trace regions + matching colour palette.

    For β the manuscript palette in preprint_11_beta_anatomy_brain is
    re-used verbatim (same region → same colour across figures). For
    other bands the A3 set is auto-discovered from the anatomy CSV and
    coloured from a categorical map. `Unk` is always excluded.
    """
    if band == "beta":
        return list(BETA_GRASS_REGIONS), list(PALETTE_BETA_GRASS)
    a3 = sorted(r for r in _a3_trace_regions(band) if r and r != "Unk")
    if not a3:
        return [], []
    cmap = plt.get_cmap("tab10")
    palette = [matplotlib.colors.to_hex(cmap(i % 10)) for i in range(len(a3))]
    return a3, palette


def _enrichment_by_region(band: str, regions: list[str]) -> dict[str, float]:
    """Per-region enrichment ratio (k_trace / expected) for *regions*.

    Reads the same anatomy CSV used by `_a3_trace_regions`. Regions not
    found in the CSV (or with NaN enrichment) default to 1.0.
    """
    base = ROOT / "data" / "audit"
    cand_clusterext = base / f"anatomy_{band}_grassmann_clusterext" / "cohort_summary.csv"
    cand_plain = base / f"anatomy_{band}_grassmann" / "cohort_summary.csv"
    csv = cand_clusterext if cand_clusterext.exists() else cand_plain
    if not csv.exists():
        return {r: 1.0 for r in regions}
    df = pd.read_csv(csv)
    out: dict[str, float] = {}
    for r in regions:
        sub = df[df.region == r]
        if sub.empty:
            out[r] = 1.0
            continue
        e = float(sub.iloc[0]["enrichment"])
        out[r] = e if np.isfinite(e) else 1.0
    return out


@functools.lru_cache(maxsize=None)
def _cohort_contacts(band: str) -> pd.DataFrame:
    """Cohort-aggregated contact DataFrame: patient, x, y, z (MNI mm), region.

    Identical pooling protocol to preprint_11_beta_anatomy_brain.
    Cached on `band` so the figure-pair (annotated + visual_only)
    only pays the per-patient MNI fit once.
    """
    rows = []
    for pat in COHORT:
        meta = load_spatial_metadata(pat, SEEG_DATAPATH)
        meta = meta.dropna(subset=["x", "y", "z"]).reset_index(drop=True)
        coords = prepare_spatial_coordinates(meta, scale="mm",
                                              center=False, to_mni=True)
        regions = meta["Desikan-Killany"].apply(_primary_region).values
        for i, r in enumerate(regions):
            rows.append(dict(
                patient=pat, x=coords[i, 0], y=coords[i, 1],
                z=coords[i, 2], region=r if r is not None else "unknown"))
    return pd.DataFrame(rows)


def _per_patient_distribution(contacts: pd.DataFrame, band: str,
                                regions: list[str]) -> None:
    """Print per-region distinct-patient counts so the user can verify
    the trace anatomy is cohort-distributed, not single-patient.
    """
    print(f"  [{band}] per-region patient distribution (distinct patients / 10):")
    for r in regions:
        sub = contacts[contacts.region == r]
        pats = sub["patient"].unique()
        n_p = len(pats)
        n_e = len(sub)
        marker = "" if n_p >= 3 else "  *single-patient warning*"
        print(f"    {r:<40s}  n_pat={n_p:>2d}  n_contacts={n_e:>3d}{marker}")


def _add_glassbrain_inset(ax_h, band: str, visual_only: bool) -> None:
    """Embed a unified T-shape glass-brain inset in the i > k triangle of `ax_h`.

    All three anatomical views (sagittal/coronal/axial) PLUS the marker
    size legend are rendered in **one** throw-away figure at IDENTICAL
    mm/inch scale by sizing each sub-axes proportionally to its brain
    template's MNI bounding box (sagittal/axial are 195×156 mm, coronal
    is 156×156 mm). The whole T-shape is saved as a single PNG and
    embedded as a single `inset_axes` in `ax_h`.

    Layout is fully relative: no per-view bbox tuning, no hardcoded
    placement of the size legend. The only knobs are
    `INSET_W` / `INSET_H` (axes-fraction footprint of the inset inside
    ax_h's i > k triangle), constrained by `INSET_W + INSET_H ≤ 1`
    (the i = k diagonal in axes-frac).
    """
    regions, palette = _band_a3_regions_palette(band)
    enrichment = _enrichment_by_region(band, regions)
    contacts = _cohort_contacts(band)
    _per_patient_distribution(contacts, band, regions)

    n_total = len(contacts)
    n_trace = int(contacts.region.isin(regions).sum())
    print(f"  [{band}] glass-brain: {n_trace}/{n_total} contacts in "
          f"{len(regions)} A3 trace region(s)")

    if enrichment:
        e_vals = list(enrichment.values())
        e_lo, e_hi = float(min(e_vals)), float(max(e_vals))
    else:
        e_lo, e_hi = 1.0, 2.0
    S_SIG_LO, S_SIG_HI = 60.0, 160.0
    S_BG = 2.5

    def _region_marker_size(region: str) -> float:
        e = enrichment.get(region, 1.0)
        if e_hi == e_lo:
            return 0.5 * (S_SIG_LO + S_SIG_HI)
        t = (e - e_lo) / (e_hi - e_lo)
        return S_SIG_LO + t * (S_SIG_HI - S_SIG_LO)

    from io import BytesIO
    from nilearn.plotting import plot_glass_brain

    # ---- one consolidated tmp_fig at uniform mm/inch scale ---------------
    # Brain MNI extents per display mode (long / short = 195/156 mm):
    #   "x" sagittal : 195 mm horiz × 156 mm vert
    #   "y" coronal  : 156 mm horiz × 156 mm vert
    #   "z" axial    : 156 mm horiz × 195 mm vert
    # If each axes has aspect (W, H) proportional to its brain's
    # (W_mm, H_mm), nilearn's `aspect='equal'` autofit gives EVERY
    # brain the same mm/inch scale. Stacking them in a T (sagittal +
    # axial in the bottom row, coronal centred on top, size legend in
    # the empty top-right quadrant) yields a single, consistent PNG.
    LONG_OVER_SHORT = 195.0 / 156.0
    SHORT_IN = 2.0
    LONG_IN = LONG_OVER_SHORT * SHORT_IN

    brain_w_in = LONG_IN + SHORT_IN
    brain_h_in = LONG_IN + SHORT_IN  # symmetric square outer canvas
    brain_fig = plt.figure(figsize=(brain_w_in, brain_h_in), dpi=200)

    # Bottom row: sagittal on the left (LONG × SHORT), axial on the
    # right (SHORT × LONG). Top row: coronal centred (SHORT × SHORT).
    ax_x = brain_fig.add_axes((
        0.0, 0.0,
        LONG_IN / brain_w_in,
        SHORT_IN / brain_h_in,
    ))
    ax_z = brain_fig.add_axes((
        LONG_IN / brain_w_in, 0.0,
        SHORT_IN / brain_w_in,
        LONG_IN / brain_h_in,
    ))
    y_w_frac = SHORT_IN / brain_w_in
    y_h_frac = SHORT_IN / brain_h_in
    ax_y = brain_fig.add_axes((
        0.5 * (1.0 - y_w_frac),
        LONG_IN / brain_h_in,
        y_w_frac, y_h_frac,
    ))

    def _render_brain_in_axes(ax, display_mode):
        disp = plot_glass_brain(
            None, figure=brain_fig, axes=ax,
            display_mode=display_mode,
            annotate=False, plot_abs=False, colorbar=False, black_bg=False,
        )
        sig_mask = contacts.region.isin(regions).values
        bg_coords = contacts.loc[~sig_mask, ["x", "y", "z"]].values
        if bg_coords.size > 0:
            disp.add_markers(marker_coords=bg_coords,
                              marker_color=NONSIG_COLOR,
                              marker_size=S_BG, alpha=0.50,
                              edgecolors="white", linewidths=0.15)
        for region, color in zip(regions, palette):
            sel = contacts.loc[contacts.region == region,
                                ["x", "y", "z"]].values
            if sel.size == 0:
                continue
            disp.add_markers(marker_coords=sel, marker_color=color,
                              marker_size=_region_marker_size(region),
                              alpha=0.92,
                              edgecolors="0.15", linewidths=0.5)

    _render_brain_in_axes(ax_x, "x")
    _render_brain_in_axes(ax_y, "y")
    _render_brain_in_axes(ax_z, "z")

    # (No size legend / no enrichment-ratio text — only axis labels and
    # tick labels are allowed as in-plot text; enrichment scaling is
    # documented in the figure caption.)

    buf = BytesIO()
    brain_fig.savefig(buf, format="png", dpi=200,
                       bbox_inches="tight", pad_inches=0.02,
                       transparent=True)
    plt.close(brain_fig)
    buf.seek(0)
    brain_png = plt.imread(buf)

    # ---- embed ONE consolidated PNG as ONE inset in ax_h -----------------
    # Inset bbox is allowed past the i = k diagonal because the T-shape
    # leaves the PNG's top-right corner empty (transparent background),
    # so the heatmap shows through under the diagonal there. The
    # *brain content* itself stays below the diagonal: the coronal's
    # top-right is at (PAD + 0.722*W, PAD + H) → 0.02 + 0.722*0.55 +
    # 0.02 + 0.55 = 0.987 < 1, axial's top-right at (PAD + W, PAD +
    # 0.556*H) = 0.57 + 0.326 = 0.896 < 1.
    INSET_W = 0.55
    INSET_H = 0.55
    INSET_PAD = 0.02
    ax_brain = inset_axes(
        ax_h,
        width="100%", height="100%",
        bbox_to_anchor=(INSET_PAD, INSET_PAD, INSET_W, INSET_H),
        bbox_transform=ax_h.transAxes,
        loc="lower left", borderpad=0,
    )
    ax_brain.imshow(brain_png, aspect="equal", interpolation="bilinear")
    ax_brain.set_axis_off()


# ---------------------------------------------------------------------------
# Plot driver
# ---------------------------------------------------------------------------
def build_figure(band: str, vmax_global: float | None = None,
                 visual_only: bool = False,
                 weight_method: str = "wilcoxon") -> Path:
    """Render the per-band Grassmann heatmap figure.

    visual_only=True (default in `main`): heatmap + per-mode profiles +
    T_G subplot + brain inset, no titles / verdict box / annotations.
    Output filename pattern:

        weight_method="wilcoxon" → fig_<band>_grassmann_heatmap_visual.pdf
        weight_method="bootstrap" → fig_<band>_grassmann_heatmap_visual_bootstrap.pdf

    `weight_method` controls which per-k cohort p-value drives the
    heatmap weighting (and the significance overlays):
      - "wilcoxon" reads `paired_wilcoxon_p` from the matched-strength
        cohort summary (current default; admits single-patient drivers
        at n=10).
      - "bootstrap" uses `compute_bootstrap_cohort_p` — patient
        resampling with replacement, B=2000, one-sided p on the cohort
        median of d_i = obs − surr_p50. Self-protects against single-
        patient leverage.
    """
    band_tex = BRAIN_BAND_TEX_DICT[band]
    out_dir = ROOT / "data" / "preprint" / "figures" / band / "grassmann"
    out_dir.mkdir(parents=True, exist_ok=True)

    pa = compute_principal_angles_full_k(band)
    pa_surr = compute_surrogate_principal_angles_cohort(band)
    ms_full = pd.read_csv(GRASS_DIR / "cohort_summary.csv")
    ms_epix = pd.read_csv(GRASS_EPIX_DIR / "cohort_summary.csv")

    H, ks = build_heatmap(pa, pa_surr)
    if weight_method == "wilcoxon":
        ks_full, n_full, p_full = cohort_strip(ms_full, band)
        ks_epix, n_epix, p_epix = cohort_strip(ms_epix, band)
    elif weight_method == "bootstrap":
        ks_full, _, p_full = cohort_strip_bootstrap(band, source="full")
        ks_epix, _, p_epix = cohort_strip_bootstrap(band, source="epix")
        n_full = np.zeros_like(p_full)
        n_epix = np.zeros_like(p_epix)
    else:
        raise ValueError(f"unknown weight_method: {weight_method!r}")

    # Locked cluster-extent verdict for this band (VERDICT_LEDGER.md).
    verdict = cluster_verdict_row(band)
    run_len, run_k_lo, run_k_hi = longest_run_span(ks_full, p_full)
    verdict_clr = VERDICT_COLOR[verdict["verdict"]]

    finite = H[np.isfinite(H)]
    if vmax_global is not None:
        vmax = vmax_global
    else:
        vmax = float(np.percentile(np.abs(finite), 95))
        vmax = max(vmax, 0.05)

    # `layout="constrained"` lets matplotlib auto-reserve space for every
    # tick label, axis label, suptitle, and colorbar label. No hardcoded
    # `left/right/top/bottom/wspace/hspace` — the only knobs are the
    # GridSpec width/height ratios (purely relative).
    fig = plt.figure(figsize=(12.0, 10.5), layout="constrained")
    outer = GridSpec(
        2, 3, figure=fig,
        width_ratios=[1.0, 3.6, 0.12],
        height_ratios=[3.6, 1.0],
    )
    ax_p = fig.add_subplot(outer[0, 0])
    ax_h = fig.add_subplot(outer[0, 1])
    cax = fig.add_subplot(outer[0, 2])
    ax_legend_tg = fig.add_subplot(outer[1, 0])
    ax_legend_tg.axis("off")
    ax_tg = fig.add_subplot(outer[1, 1], sharex=ax_h)

    # Pre-compute per-k significance masks for the overlay (drawn AFTER
    # the heatmap so the stripes sit above the data).
    sig_full_ks = np.array([int(kk) for kk, pp in zip(ks_full, p_full)
                             if K_MIN <= kk <= K_MAX and pp < 0.05])
    sig_epix_ks = np.array([int(kk) for kk, pp in zip(ks_epix, p_epix)
                             if K_MIN <= kk <= K_MAX and pp < 0.05])
    if len(ks_epix):
        epix_k_max = int(ks_epix.max())
    else:
        epix_k_max = K_MIN - 1

    # ---- bridged cluster-mass support (columns inside the cluster) ------
    # The cluster-mass null distribution (audit_70) bridges single empty
    # k cells inside a contiguous-significant run. Compute the bridged
    # cluster-support set FIRST so the heatmap can use it to fade out
    # columns that do not carry the cohort-level trace claim.
    sig_mask_bar = p_full < 0.05
    sig_indices_b = np.where(sig_mask_bar)[0]
    bridged_clusters: list[tuple[int, int]] = []
    if len(sig_indices_b):
        c_start = sig_indices_b[0]
        c_prev = sig_indices_b[0]
        for ix in sig_indices_b[1:]:
            if ix - c_prev - 1 > 1:        # gap > 1 non-sig cell → new cluster
                bridged_clusters.append((c_start, c_prev))
                c_start = ix
            c_prev = ix
        bridged_clusters.append((c_start, c_prev))
    support_k_set: set[int] = set()
    for r_lo, r_hi in bridged_clusters:
        for r in range(r_lo, r_hi + 1):
            support_k_set.add(int(ks_full[r]))

    # ---- heatmap: cluster-mass-weighted residual, continuous BWR cmap ----
    # The reader's "is this a trace?" question is now answered by
    # *weighting the data* by the cohort-level evidence at each k,
    # rather than thresholding the colormap. The weight at column k is
    #
    #   w(k) = clip(-log10(p_full[k]) / -log10(0.05), 0, 1)
    #
    # so columns with p < 0.05 get full weight (1.0) and those with
    # p → 1 ramp smoothly to zero. The plotted quantity is
    # H_weighted[i,k] = H[i,k] · w(k). The colormap is a plain
    # continuous blue → pure-white → red diverging map. Visually:
    #   - β trace cluster (many k with p < 0.05, large |H|) saturates red.
    #   - α (few or no k with p < 0.05) fades smoothly to white at
    #     non-significant k regardless of the unweighted residual sign.
    log_05 = -np.log10(0.05)
    neglog_p = -np.log10(np.maximum(p_full, 1e-12))
    weight_per_k = np.clip(neglog_p / log_05, 0.0, 1.0)
    if len(weight_per_k) != H.shape[1]:
        # Defensive: ks_full and the heatmap k axis should match in
        # length; if they ever drift, broadcast to the heatmap shape.
        weight_per_k = np.interp(ks, ks_full, weight_per_k)
    # Band-level cluster-extent verdict (audit_70). Strong-trace bands keep
    # full weight; no-trace bands are muted so the figure cannot claim
    # evidence the band-level verdict already rejected.
    _ce_csv = Path("data/audit/grassmann_cluster_extent/cohort_summary.csv")
    _ce = pd.read_csv(_ce_csv)
    _p_band = float(_ce.loc[_ce.band == band, "cluster_p_cluster_mass"].iloc[0])
    if   _p_band < 0.01: band_scale = 1.00  # strong → unchanged
    elif _p_band < 0.05: band_scale = 0.50  # weak
    else:                band_scale = 0.15  # no_trace → muted
    weight_per_k = band_scale * weight_per_k
    W = np.tile(weight_per_k, (I_MAX, 1))
    H_weighted = H * W

    TRACE_POS_COLOR = "#b8001f"      # deep red — trace
    TRACE_NEG_COLOR = "#001fb8"      # deep blue — anti-trace
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "trace_continuous",
        [(0.00, TRACE_NEG_COLOR),
         (0.50, "#ffffff"),           # pure white midpoint
         (1.00, TRACE_POS_COLOR)])
    cmap.set_bad("white")   # i > k triangle: white background

    H_masked = np.ma.masked_invalid(H_weighted)
    im = ax_h.imshow(
        H_masked, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax,
        origin="upper", interpolation="nearest",
        extent=(ks.min() - 0.5, ks.max() + 0.5, I_MAX + 0.5, 0.5),
    )
    ax_h.plot([ks.min() - 0.5, ks.max() + 0.5],
              [ks.min() - 0.5, ks.max() + 0.5],
              color="0.40", lw=0.7, ls="--", zorder=4)
    # The lower-left triangle (i > k) is no longer empty: a top-down
    # glass-brain inset renders the cohort sEEG electrodes coloured by
    # whether their DK region passes the per-band A3 trace anatomy test
    # (audit_72 cluster-extent set when available, else plain
    # `anatomy_<band>_grassmann`).
    _add_glassbrain_inset(ax_h, band, visual_only=visual_only)

    ax_h.set_ylim(I_MAX + 0.5, 0.5)
    x_ticks = [2, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 112]
    y_ticks = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 112]
    ax_h.set_xticks(x_ticks)
    ax_h.set_xticklabels([str(t) for t in x_ticks])
    plt.setp(ax_h.get_xticklabels(), visible=False)
    plt.setp(ax_h.get_yticklabels(), visible=False)
    ax_h.spines[["top", "right"]].set_visible(False)

    # ---- shared colorbar (GridSpec cell — height matches ax_h) ------------
    cb = fig.colorbar(im, cax=cax)
    cb.set_ticks([-vmax, -0.5 * vmax, 0.0, 0.5 * vmax, vmax])
    cb.set_ticklabels([f"{-vmax:.02f}",
                         f"{-0.5*vmax:.02f}",
                         "0",
                         f"{0.5*vmax:.02f}",
                         f"{vmax:.02f}"])
    cb.set_label(r"$w(k)\,\Delta\theta_i(k)$  (rad)")

    # ---- T_G(k) panel: observed cohort median vs matched-strength ---------
    # Load per-patient T_G to build the surrogate IQR envelope (median
    # over patients of each patient's own p25/p50/p75 across R=200
    # surrogate replicates). Plot observed cohort-median T_G(k) for both
    # full FC and epi-X overlaid; shade the surrogate IQR (full FC) so
    # the cohort signal vs strength-controlled null is visible by eye.
    per_pat_full = pd.read_csv(GRASS_DIR / GRASS_PER_PAT_CSV)
    per_pat_epix = pd.read_csv(GRASS_EPIX_DIR / GRASS_PER_PAT_CSV)
    env_full = cohort_surrogate_envelope(per_pat_full, band)
    env_epix = cohort_surrogate_envelope(per_pat_epix, band)

    cohort_full_df = ms_full[ms_full.band == band].sort_values("k")
    cohort_epix_df = ms_epix[ms_epix.band == band].sort_values("k")

    # Sign convention: auto-detect per CSV via the count-column name.
    # New convention (positive = trace, locked 2026-05-26) writes
    # `n_patients_above_own_surrogate`; legacy convention writes
    # `n_patients_below_own_surrogate` and requires a sign flip to
    # arrive at the Methods convention. Once both caches are
    # regenerated under the new convention this auto-detect becomes a
    # no-op and the `sign_*` constants can be deleted.
    sign_full = (1.0
                  if "n_patients_above_own_surrogate" in cohort_full_df.columns
                  else -1.0)
    sign_epix = (1.0
                  if "n_patients_above_own_surrogate" in cohort_epix_df.columns
                  else -1.0)
    tg_obs_full = sign_full * cohort_full_df.obs_median_T_G.values.astype(float)
    tg_obs_epix = sign_epix * cohort_epix_df.obs_median_T_G.values.astype(float)
    tg_surr_full = (sign_full
                     * cohort_full_df.surr_median_T_G_per_patient_median.values.astype(float))

    h_iqr_full = ax_tg.fill_between(
        env_full.k.values, env_full.surr_p25.values, env_full.surr_p75.values,
        color=CLR_SURR_FILL, alpha=0.95, zorder=1)
    h_iqr_epix = ax_tg.fill_between(
        env_epix.k.values, env_epix.surr_p25.values, env_epix.surr_p75.values,
        color=CLR_EPIX, alpha=0.12, zorder=1)
    h_surr_full, = ax_tg.plot(ks_full, tg_surr_full,
                               color=CLR_SURR, lw=1.8, zorder=2)
    h_obs_full, = ax_tg.plot(ks_full, tg_obs_full,
                              color=CLR_TG_OBS, lw=2.6, zorder=4)
    h_obs_epix, = ax_tg.plot(ks_epix, tg_obs_epix,
                              color=CLR_TG_EPIX, lw=2.6, zorder=4)
    ax_tg.axhline(0, color="0.4", lw=0.6, ls="--", zorder=1)
    ax_tg.set_xlim(ks.min() - 0.5, ks.max() + 0.5)
    ax_tg.set_xticks(x_ticks)
    ax_tg.set_xticklabels([str(t) for t in x_ticks])
    ax_tg.set_ylabel(r"$T_G(k)$")
    ax_tg.spines[["top", "right"]].set_visible(False)

    # T_G legend stays in its dedicated col-0 row-1 cell (one column,
    # stacked entries, no frame). With `layout="constrained"` the cell
    # is sized and separated from ax_tg's y-axis decoration
    # automatically — no manual bbox_to_anchor.
    ax_legend_tg.legend(
        handles=[
            mlines.Line2D([], [], color=CLR_TG_OBS, lw=2.6,
                          label="obs full FC"),
            mlines.Line2D([], [], color=CLR_TG_EPIX, lw=2.6,
                          label="obs epi-X"),
            mlines.Line2D([], [], color=CLR_SURR, lw=1.8,
                          label="surr median"),
            mpatches.Patch(facecolor=CLR_SURR_FILL,
                           label="surr IQR"),
        ],
        loc="center left", frameon=False, handlelength=1.4,
    )

    # The bottom-row plot is now ax_tg (T_G line). Show its x-axis
    # tick labels and the k label, replacing what used to be ax_n.
    plt.setp(ax_tg.get_xticklabels(), visible=True)
    ax_tg.set_xlabel(r"$k$")

    # ---- right-side companion: per-mode profiles at representative k -----
    # Clipped magma (0.12 -> 0.78) — keeps the warm sequential feel of
    # the original plasma palette while excluding the near-white yellow
    # end that made the largest-k profile unresolvable on a white
    # background.
    cmap_profiles = plt.get_cmap("magma")
    profile_fracs = np.linspace(0.12, 0.78, len(K_PROFILES))

    def _profile(ax, k, color):
        col = pa[pa.k == k].sort_values("mode_index")
        if col.empty:
            return
        y = col.mode_index.values.astype(int)
        dth = col.delta_theta_median_radians.values.astype(float)
        ax.plot(dth, y, color=color, lw=2.4, zorder=2,
                label=f"{k}")

    for k, frac in zip(K_PROFILES, profile_fracs):
        clr = cmap_profiles(frac)
        _profile(ax_p, k, clr)

    ax_p.axvline(0, color="0.4", lw=0.6, ls="--", zorder=1)
    ax_p.invert_yaxis()
    ax_p.invert_xaxis()
    ax_p.set_ylim(I_MAX + 0.5, 0.5)
    ax_p.set_yticks(y_ticks)
    ax_p.set_yticklabels([str(t) for t in y_ticks])
    ax_p.spines[["top", "right"]].set_visible(False)
    ax_p.set_xlabel(r"$\Delta\theta_i(k)$  (rad)")
    ax_p.set_ylabel(r"$i$")
    ax_p.legend(loc="lower left", frameon=False,
                 ncol=2, handlelength=1.4, title=r"$k$")

    suffix = "_visual" if visual_only else ""
    if weight_method == "bootstrap":
        suffix += "_bootstrap"
    out = out_dir / f"fig_{band}_grassmann_heatmap{suffix}.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"  [{band}/{weight_method}{' visual' if visual_only else ''}] "
          f"saved → {out}")
    return out


def main() -> list[Path]:
    bands = sys.argv[1:] if len(sys.argv) > 1 else ALL_BANDS
    bad = [b for b in bands if b not in ALL_BANDS]
    if bad:
        raise SystemExit(
            f"unknown band(s): {bad}; choose from {ALL_BANDS}")

    # Fixed global v-max for the residual heatmap. β is the strongest
    # band (95th-pct |residual| ≈ 0.06); fixing vmax = 0.07 saturates β
    # near the colormap extremes and pushes weaker bands (α, γ_l, γ_h,
    # δ, θ) into the muted middle range. Trace bands therefore read as
    # strongly red where the trace lives; no-trace bands read as washed
    # out — the magnitude difference between bands now carries through
    # the heatmap.
    VMAX_FIXED = 0.07
    # Render both the locked Wilcoxon-weighted figure (current
    # preprint default) AND the new patient-bootstrap variant in
    # parallel. The bootstrap PDF is suffixed `_visual_bootstrap` —
    # NEVER overwrites the Wilcoxon variant. Once the bootstrap is
    # validated against the Wilcoxon results, the integration into
    # the preprint statistics pipeline becomes a follow-up step.
    outs = []
    for b in bands:
        outs.append(build_figure(b, vmax_global=VMAX_FIXED,
                                  visual_only=True,
                                  weight_method="wilcoxon"))
        outs.append(build_figure(b, vmax_global=VMAX_FIXED,
                                  visual_only=True,
                                  weight_method="bootstrap"))
    return outs


if __name__ == "__main__":
    main()
