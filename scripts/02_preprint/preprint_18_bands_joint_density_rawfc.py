#!/usr/bin/env python3
"""Six-band cohort joint density at the per-pair substrate layer.

``--layer {raw,coph}``  ×  ``--mode {copula,drift,gate,median_excess}``

Substrate: ``raw`` = raw ``|ImCoh|`` adjacency ``A``; ``coph`` = LRG cophenetic
distance ``D``. Two floor recipes, identical panel style so every pairing is a
clean swap. Output names encode both, so each is its OWN figure file:

    data/preprint/figures/all_bands/fig_bands_joint_density_{copula,drift}_{rawfc,coph}.pdf

mode = copula  (headline selectivity figure)
------------------------------------------------
Gaussian-copula render. The heatmap is the signed excess of
``copula(cohort-median ρ_obs)`` over ``copula(surr_median_rho_p95)`` (the
matched-strength gate's p95 null). Both are parameterized by a single scalar ρ,
so this render is **tie-immune and robust to per-patient heterogeneity** — and
it reproduces the matched-strength GATE verdict exactly: ``raw`` lights every
band (non-specific), ``coph`` only β (strongest) and α.

WHY a copula model and not the raw empirical density: the cophenetic distance is
a near-ultrametric — ``dD`` is ~91% tied — so the *empirical* pooled-rank density
piles spurious mass on the diagonal in EVERY band regardless of trace, and the
cohort-mean is further inflated by a few high-ρ patients. Spearman ρ (hence the
copula) is tie-corrected and median-based, so it is immune to both. This panel
is therefore a faithful MODEL render of the per-band scalar (ρ_obs vs p95-null),
NOT the raw joint density. (raw ``|ImCoh|`` is continuous/tie-free, so its copula
and empirical renders agree — the divergence is a cophenetic-substrate effect.)

mode = drift  (separate empirical companion)
------------------------------------------------
Fully empirical. The heatmap is the signed excess of the observed per-patient-KDE
density over the WITHIN-SESSION DRIFT density — rest_pre halves vs rest_post
halves, no task: ``dD_task→D^preB−D^preA``, ``dD_rest→D^postB−D^postA``. The
drift null is the ONLY empirical floor that is tie- AND heterogeneity-matched to
the observed density (it shares the same cophenetic ties and per-patient
structure, differing only by the absence of task), so for a null band obs≈drift
→ cream and for a trace band obs>drift → green. Shows the drift-cleared set
(α/β/γ_l); the matched-strength gate (α/β) is the copula panel + the
null-triangle (preprint_21).

Construction (matches continuous_trace_matrix._per_cell_split exactly):
  obs:   dD_task = (D_test − D_pre_A)[iu],  dD_rest = (D_post − D_pre_B)[iu]
         (raw: A^tt/A^post full; coph: full task/post from IMCOH_LRG_CACHE;
          half rest_pre baselines from the halves cache)
  drift: dD_task = (D_pre_B − D_pre_A)[iu], dD_rest = (D_post_B − D_post_A)[iu]
         (all halves)

mode = median_excess  (per-patient-subtracted empirical field)
------------------------------------------------
Fully empirical, but the per-patient subtraction handles the cophenetic tie
artifact by cancellation instead of by random tie-breaking (the ``gate`` mode's
route). For each patient: a fixed-bandwidth histogram density of its midrank
pairs minus the per-pixel median (or p95, ``--floor``) of ITS OWN matched-
strength surrogate density fields; the cohort field is the per-pixel median of
those per-patient excesses. Because observed and surrogate share the same
ultrametric per patient, the diagonal tie-pile cancels where the geometry
matches — a pooled obs−surrogate cannot do this (piles never line up across
patients). Two renders: ``(u,v)`` (diagonal trace, with the same conditional-
split marginals) and ``(s,d)`` with ``s=(u+v)/2, d=u−v`` (trace = green band at
``d≈0``, inside the valid diamond ``|d|≤2·min(s,1−s)``; outside = white).
Companion data in ``data/preprint/median_excess/``: ``medexc_diagonal_mass.csv``
(central-strip positive mass vs off-diagonal) and ``medexc_tie_fractions.csv``
(obs vs surrogate tie fractions per patient — validates the cancellation; flags
gaps > 0.05). Smoothing is a single fixed σ for obs/surrogate/all patients/all
bands — NOT a free knob, so α/β selectivity cannot be manufactured.

Nested HDR contours (0.10…0.99) of the heatmap density overlay the copula/drift/
gate panels; the two marginals show the conditional split (distribution of ``u``
for top- vs bottom-half ``v`` and vice versa). n = 10. Layout 2×3, band order
δ θ α / β γ_l γ_h.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import squareform
from scipy.stats import gaussian_kde, norm, rankdata, spearmanr

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import CACHE_ROOT, IMCOH_LRG_CACHE
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.utils.surrogate import (
    adjacency_from_laplacian_eigs,
    cophenetic_condensed_from_eigs,
    surrogate_cache_path,
)
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result

ROOT = setup_script_env()
use_lrg_style()

# Per-script font bump. All non-band-label text inherits these; band labels
# stay fixed at 22 pt as the dominant element.
matplotlib.rcParams.update({
    "font.size": 12.5,
    "axes.labelsize": 12.5,
    "axes.titlesize": 12.5,
    "xtick.labelsize": 11.5,
    "ytick.labelsize": 11.5,
    "legend.fontsize": 11.5,
})


COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]

RAW_MS_DIR = ROOT / "data" / "audit" / "raw_fc_matched_strength"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"          # raw |ImCoh| halves
LRG_HALVES_CACHE = CACHE_ROOT / "imcoh_lrg_halves"        # cophenetic halves
COPH_PAIR_SPLIT_DIR = (
    ROOT / "data" / "reports" / "imcoh_continuous_trace" / "per_pair_split"
)
COPH_MS_DIR = ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"

# Cached matched-strength surrogate eigendecompositions (split-baseline: half
# rest_pre A/B, full task_test / rest_post). One ensemble drives both layers.
GATE_FLOOR_CACHE = CACHE_ROOT / "gate_floor_density_scratch"
SURR_R, SURR_SWAP, SURR_SEED = 200, 20, 20260511
SURR_PHASES = ("task_test", "rest_pre_A", "rest_post", "rest_pre_B")
TIEBREAK_SEED = 20260603

# Diverging trace cmap (green diagonal = trace, red anti-diagonal = anti, cream
# = at/below noise floor). Same as preprint_09 / preprint_17.
TRACE_CMAP = LinearSegmentedColormap.from_list(
    "trace_rwg",
    [
        "#4a0e0e", "#a02525", "#d8584c", "#f0a89e",
        "#fdf2ef",
        "#bfddc7", "#5fac7e", "#1a7c3e", "#0a3a1d",
    ],
    N=512,
)

GRID_N = 120
KDE_BW = 0.22
# median_excess mode: per-patient histogram density, single fixed Gaussian
# smoothing (NO per-band/per-panel bandwidth tuning — see honesty note in the
# mode docstring), subtracted against each patient's OWN matched-strength
# surrogate field before any cohort aggregation.
GRID_MEDEXC = 64
SIGMA_MEDEXC = 1.0                      # Gaussian smoothing, in bins
MEDEXC_CACHE = CACHE_ROOT / "medexc_perpatient_fields"
MEDEXC_DATA = ROOT / "data" / "preprint" / "median_excess"
# Colour scale. copula deviations are small (≈0.025, as preprint_09); empirical
# KDE deviations are O(1), so the drift mode sets vabs from the data (robust pct).
VABS_COPULA = 0.025
HDR_COVERAGES = (0.10, 0.25, 0.40, 0.55, 0.68, 0.80, 0.90, 0.95, 0.99)


# ---------------------------------------------------------------------------
# Condensed-vector helpers (raw |ImCoh| adjacency; audit_67 convention)
# ---------------------------------------------------------------------------
def _clean(W: np.ndarray) -> np.ndarray:
    """Diagonal-zeroed, [0,1]-clipped, symmetrized (audit_67.load_phase_fc)."""
    W = np.asarray(W, dtype=float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def _condensed(W: np.ndarray) -> np.ndarray:
    W = W.copy()
    np.fill_diagonal(W, 0.0)
    W = np.maximum(W, W.T)
    return squareform(W, checks=False)


def _raw_half_path(pat: str, band: str, phase: str, half: str) -> Path:
    return HALVES_FC_CACHE / pat / f"{band}_{phase}_{half}_imcoh_abs.npy"


def load_raw_half(pat: str, band: str, phase: str, half: str) -> np.ndarray:
    p = _raw_half_path(pat, band, phase, half)
    if not p.exists():
        raise FileNotFoundError(
            f"missing raw half-FC {p}\n  run audit_67/audit_74 to populate the "
            f"imcoh_halves_fc cache (rest_pre AND rest_post halves)"
        )
    return _clean(np.load(p))


def load_raw_full(pat: str, band: str, phase: str) -> np.ndarray:
    return _clean(load_fc_matrix(pat, phase, band, fc_method="imcoh_abs"))


# ---------------------------------------------------------------------------
# Cophenetic-distance helpers (D = squareform(ultrametric_matrix), audit_63 /
# continuous_trace_matrix convention). Full task/post from IMCOH_LRG_CACHE;
# halves from imcoh_lrg_halves.
# ---------------------------------------------------------------------------
def load_D_square(pat: str, phase: str, band: str, cache_root) -> np.ndarray:
    r = load_lrg_result(pat, phase, band, fc_method="imcoh_abs",
                        cache_root=cache_root)
    if r is None:
        raise FileNotFoundError(f"missing LRG result {pat}/{phase}/{band} "
                                f"in {cache_root}")
    um = np.asarray(r.ultrametric_matrix)
    return squareform(um) if um.ndim == 1 else um


def _coph_condensed(pat: str, phase: str, band: str, cache_root) -> np.ndarray:
    D = load_D_square(pat, phase, band, cache_root)
    iu = np.triu_indices(D.shape[0], k=1)
    return D[iu]


# ---------------------------------------------------------------------------
# OBSERVED per-pair shift vectors (one (dD_task, dD_rest) per patient)
# ---------------------------------------------------------------------------
def raw_obs_pairs(band: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    pairs = {}
    for pat in COHORT:
        A_preA = load_raw_half(pat, band, "rest_pre", "A")
        A_preB = load_raw_half(pat, band, "rest_pre", "B")
        A_tt = load_raw_full(pat, band, "task_test")
        A_post = load_raw_full(pat, band, "rest_post")
        if len({M.shape for M in (A_preA, A_preB, A_tt, A_post)}) != 1:
            raise ValueError(f"{pat}/{band}: raw phase shape mismatch")
        pairs[pat] = (_condensed(A_tt) - _condensed(A_preA),
                      _condensed(A_post) - _condensed(A_preB))
    return pairs


def coph_obs_pairs(band: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Cached continuous_trace split-baseline shift vectors (D_test/D_post full,
    rest_pre halves) — identical to continuous_trace_matrix._per_cell_split."""
    pairs = {}
    for pat in COHORT:
        p = COPH_PAIR_SPLIT_DIR / f"{pat}_{band}.npz"
        if not p.exists():
            raise FileNotFoundError(f"missing cophenetic per-pair-split {p}")
        d = np.load(p)
        pairs[pat] = (np.asarray(d["dD_task"]), np.asarray(d["dD_rest"]))
    return pairs


# ---------------------------------------------------------------------------
# DRIFT per-pair shift vectors (within-session, no task; tie-/heterogeneity-
# matched to the observed pairs). dD_task→preB−preA, dD_rest→postB−postA.
# ---------------------------------------------------------------------------
def raw_drift_pairs(band: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    pairs = {}
    for pat in COHORT:
        A_preA = load_raw_half(pat, band, "rest_pre", "A")
        A_preB = load_raw_half(pat, band, "rest_pre", "B")
        A_postA = load_raw_half(pat, band, "rest_post", "A")
        A_postB = load_raw_half(pat, band, "rest_post", "B")
        if len({M.shape for M in (A_preA, A_preB, A_postA, A_postB)}) != 1:
            raise ValueError(f"{pat}/{band}: raw drift shape mismatch")
        pairs[pat] = (_condensed(A_preB) - _condensed(A_preA),
                      _condensed(A_postB) - _condensed(A_postA))
    return pairs


def coph_drift_pairs(band: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    pairs = {}
    for pat in COHORT:
        dA = _coph_condensed(pat, "rest_pre_A", band, LRG_HALVES_CACHE)
        dB = _coph_condensed(pat, "rest_pre_B", band, LRG_HALVES_CACHE)
        pA = _coph_condensed(pat, "rest_post_A", band, LRG_HALVES_CACHE)
        pB = _coph_condensed(pat, "rest_post_B", band, LRG_HALVES_CACHE)
        pairs[pat] = (dB - dA, pB - pA)
    return pairs


# ---------------------------------------------------------------------------
# Rank pooling + empirical density + per-patient ρ
# ---------------------------------------------------------------------------
def pool_ranks(pairs, rng=None) -> tuple[np.ndarray, np.ndarray]:
    t_all, r_all = [], []
    for dD_task, dD_rest in pairs.values():
        n = len(dD_task)
        if n < 2:
            continue
        if rng is None:
            t_all.append((rankdata(dD_task, method="ordinal") - 0.5) / n)
            r_all.append((rankdata(dD_rest, method="ordinal") - 0.5) / n)
        else:
            t_all.append(rank_random_tiebreak(dD_task, rng))
            r_all.append(rank_random_tiebreak(dD_rest, rng))
    return np.concatenate(t_all), np.concatenate(r_all)


def cohort_average_empirical_density(pairs, grid_n: int, bw: float) -> np.ndarray:
    """Per-patient KDE on within-patient rank pairs, then cohort mean."""
    xy = np.linspace(0.0, 1.0, grid_n)
    xx, yy = np.meshgrid(xy, xy)
    grid_pts = np.vstack([xx.ravel(), yy.ravel()])
    per_pat = []
    for dD_task, dD_rest in pairs.values():
        n = len(dD_task)
        if n < 2:
            continue
        t = (rankdata(dD_task, method="ordinal") - 0.5) / n
        r = (rankdata(dD_rest, method="ordinal") - 0.5) / n
        k = gaussian_kde(np.vstack([t, r]), bw_method=bw)
        per_pat.append(k(grid_pts).reshape(grid_n, grid_n))
    return np.mean(per_pat, axis=0)


def per_patient_rhos(pairs) -> np.ndarray:
    rhos = []
    for dD_task, dD_rest in pairs.values():
        if len(dD_task) < 2:
            continue
        rhos.append(float(spearmanr(dD_task, dD_rest).statistic))
    return np.asarray(rhos)


# ---------------------------------------------------------------------------
# Robust empirical render (mode = gate): random tie-breaking + cohort-median +
# empirical matched-strength p95 floor. The two robustifications remove the two
# inflations of the naive empirical KDE, so it converges onto the copula while
# staying a genuine KDE of the data.
# ---------------------------------------------------------------------------
def rank_random_tiebreak(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Uniform (0,1) ranks with ties broken RANDOMLY and INDEPENDENTLY.

    `rankdata(method="ordinal")` breaks ties by index order — identically on
    both axes — so the ~91% co-tied cophenetic pairs (same value in task AND
    rest) get the same rank on both axes and pile exactly on u=v, manufacturing
    a green diagonal in every band. Random independent tie-breaking is the
    unbiased rank (the randomized midrank): co-tied pairs scatter as uniform
    background, and in expectation the dependence equals Spearman ρ — which is
    why the copula never shows the artifact.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    order = np.lexsort((rng.random(n), x))
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(n)
    return (ranks + 0.5) / n


def _kde_on_grid(t, r, grid_pts, grid_n, bw):
    return gaussian_kde(np.vstack([t, r]),
                        bw_method=bw)(grid_pts).reshape(grid_n, grid_n)


def cohort_median_empirical_density(pairs, grid_n, bw, rng):
    """Per-patient KDE on randomly-tie-broken ranks, then cohort MEDIAN per cell.

    Median (not mean) is the robust cohort centre — one ρ≈0.8 patient (δ Pat_06)
    cannot move it, matching the copula's use of ``obs_median_rho``.
    """
    xy = np.linspace(0.0, 1.0, grid_n)
    xx, yy = np.meshgrid(xy, xy)
    grid_pts = np.vstack([xx.ravel(), yy.ravel()])
    per_pat = []
    for dD_task, dD_rest in pairs.values():
        if len(dD_task) < 2:
            continue
        per_pat.append(_kde_on_grid(rank_random_tiebreak(dD_task, rng),
                                    rank_random_tiebreak(dD_rest, rng),
                                    grid_pts, grid_n, bw))
    return np.median(np.stack(per_pat, axis=0), axis=0)


def _load_surr_eigs(pat, band, phase):
    path = surrogate_cache_path(pat, band, phase, SURR_R, SURR_SWAP, SURR_SEED,
                                "imcoh_abs")
    if not path.exists():
        raise FileNotFoundError(f"missing matched-strength surrogate eigs {path}")
    with np.load(path) as d:
        return d["eigvals"], d["eigvecs"]


def _surr_condensed(layer, ev_r, ec_r):
    if layer == "coph":
        return cophenetic_condensed_from_eigs(ev_r, ec_r)
    return _condensed(adjacency_from_laplacian_eigs(ev_r, ec_r))


def _surrogate_rho_and_pairs(pat, band, layer):
    """Per-surrogate split-baseline (dD_task, dD_rest) and ρ for one patient.

    dD_task = D(task) − D(rest_pre_A),  dD_rest = D(rest_post) − D(rest_pre_B),
    each phase its OWN matched-strength realization r — same split-baseline
    construction as the observed pairs. coph → cophenetic distance, raw →
    adjacency recovered from the same eigendecomposition.
    """
    eigs = {ph: _load_surr_eigs(pat, band, ph) for ph in SURR_PHASES}
    R = eigs["task_test"][0].shape[0]
    dts, drs, rhos = [], [], []
    for r in range(R):
        if not all(np.isfinite(eigs[ph][0][r]).all() for ph in SURR_PHASES):
            continue
        cond = {ph: _surr_condensed(layer, eigs[ph][0][r], eigs[ph][1][r])
                for ph in SURR_PHASES}
        dt = cond["task_test"] - cond["rest_pre_A"]
        dr = cond["rest_post"] - cond["rest_pre_B"]
        dts.append(dt)
        drs.append(dr)
        rhos.append(float(spearmanr(dt, dr).statistic))
    return dts, drs, np.asarray(rhos)


def _floor_cache_path(layer, band, grid_n, bw):
    return (GATE_FLOOR_CACHE /
            f"{layer}_{band}_g{grid_n}_bw{bw:.3f}_seed{TIEBREAK_SEED}_p95floor.npz")


def cohort_median_p95_floor(band, layer, grid_n, bw, rng):
    """Empirical matched-strength GATE floor at the conservative 95th percentile.

    Per patient: the upper-tail surrogate realizations (ρ ≥ that patient's p95)
    give the dependence strength-matched noise reaches at its 95th percentile;
    KDE each (tie-broken), average → patient floor; cohort MEDIAN per cell. This
    is the band-structured floor (high for δ/γ where matched-strength
    manufactures rank-correlation, low for α/β) — the empirical twin of
    ``copula(surr_median_rho_p95)``. Cached: the per-surrogate cophenetic tail is
    the only expensive step.
    """
    cache = _floor_cache_path(layer, band, grid_n, bw)
    if cache.exists():
        with np.load(cache) as d:
            return d["floor"]
    xy = np.linspace(0.0, 1.0, grid_n)
    xx, yy = np.meshgrid(xy, xy)
    grid_pts = np.vstack([xx.ravel(), yy.ravel()])
    per_pat = []
    for pat in COHORT:
        dts, drs, rhos = _surrogate_rho_and_pairs(pat, band, layer)
        if rhos.size == 0 or not np.isfinite(rhos).any():
            continue
        # the single realization at this patient's 95th-percentile ρ (the
        # dependence strength-matched noise reaches at its conservative tail);
        # cohort-median over patients smooths the single-draw noise.
        p95 = np.nanpercentile(rhos, 95.0)
        r_star = int(np.nanargmin(np.abs(rhos - p95)))
        per_pat.append(_kde_on_grid(rank_random_tiebreak(dts[r_star], rng),
                                    rank_random_tiebreak(drs[r_star], rng),
                                    grid_pts, grid_n, bw))
    floor = np.median(np.stack(per_pat, axis=0), axis=0)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache, floor=floor)
    return floor


def gaussian_copula_density(rho: float, grid_n: int) -> np.ndarray:
    """Analytic bivariate Gaussian copula density on [0,1]² (== preprint_09)."""
    eps = 1.5e-3
    u = np.linspace(eps, 1.0 - eps, grid_n)
    uu, vv = np.meshgrid(u, u)
    nu = norm.ppf(uu)
    nv = norm.ppf(vv)
    rho2 = rho * rho
    norm_factor = 1.0 / np.sqrt(1.0 - rho2)
    expo = -0.5 / (1.0 - rho2) * (rho2 * (nu**2 + nv**2) - 2.0 * rho * nu * nv)
    return norm_factor * np.exp(expo)


def signed_excess_field(z_main: np.ndarray, z_floor: np.ndarray) -> np.ndarray:
    """sign(z_main−1)·max(|z_main−1| − |z_floor−1|, 0)."""
    dev_obs = z_main - 1.0
    dev_floor = z_floor - 1.0
    abs_excess = np.maximum(np.abs(dev_obs) - np.abs(dev_floor), 0.0)
    return np.sign(dev_obs) * abs_excess


# ---------------------------------------------------------------------------
# mode = median_excess  (per-patient-subtracted, cohort-median empirical field)
# ---------------------------------------------------------------------------
# Each patient's observed rank-pair density is differenced against ITS OWN
# matched-strength surrogate density (where observed and null share the same
# ultrametric / tie geometry), THEN aggregated by per-pixel cohort median. The
# per-patient subtraction is what makes the cophenetic tie-pile cancel: a tied
# block lands on the diagonal in both the observed and the surrogate field of
# the SAME patient, so it cancels in E_p where the geometry actually matches —
# unlike a pooled obs − pooled surrogate, where the per-patient piles never line
# up. Fixed histogram density + single Gaussian smoothing (σ = SIGMA_MEDEXC,
# applied identically to obs, every surrogate, every patient): bandwidth is NOT
# a free knob, so the render cannot be tuned to manufacture α/β selectivity.
# Floor: per-pixel surrogate MEDIAN (default, "excess over typical null") or p95
# (conservative, the gate-faithful floor). midranks (not random tie-break) — the
# cancellation, not the rank jitter, removes the tie artifact.
def _midrank_unit(x: np.ndarray) -> np.ndarray:
    """Within-vector midranks mapped to (0,1): rank/(M+1), ties → average."""
    x = np.asarray(x, dtype=float)
    return rankdata(x, method="average") / (x.size + 1)


def _density_field(dD_task: np.ndarray, dD_rest: np.ndarray,
                   grid_n: int, sigma: float) -> np.ndarray:
    """Smoothed 2D rank histogram density on [0,1]², indexed [iu, iv] (∫ = 1)."""
    u = _midrank_unit(dD_task)
    v = _midrank_unit(dD_rest)
    edges = np.linspace(0.0, 1.0, grid_n + 1)
    H, _, _ = np.histogram2d(u, v, bins=[edges, edges], density=True)
    H = gaussian_filter(H, sigma=sigma, mode="nearest")
    cell = (1.0 / grid_n) ** 2
    s = H.sum() * cell
    if s > 0:
        H = H / s
    return H


def _tie_frac(x: np.ndarray) -> float:
    """Fraction of entries whose value is shared with ≥1 other entry."""
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return float("nan")
    _vals, counts = np.unique(x, return_counts=True)
    return float(counts[counts > 1].sum() / x.size)


def _pair_tie_frac(dD_task: np.ndarray, dD_rest: np.ndarray) -> float:
    return 0.5 * (_tie_frac(dD_task) + _tie_frac(dD_rest))


def _medexc_cache_path(layer: str, band: str, grid_n: int, sigma: float) -> Path:
    return (MEDEXC_CACHE /
            f"{layer}_{band}_g{grid_n}_s{sigma:.2f}_perpatient.npz")


def compute_perpatient_fields(layer: str, band: str, grid_n: int, sigma: float):
    """Per-patient (obs, surrogate-median, surrogate-p95) density fields + tie
    fractions for *band*/*layer*. Cached — the per-surrogate cophenetic
    reconstruction is the only expensive step (shared with the gate mode).

    Returns dict with arrays shaped (P, grid_n, grid_n): ``obs``, ``null_median``,
    ``null_p95``; and (P,) vectors ``tie_obs``, ``tie_surr``; plus ``patients``.
    """
    cache = _medexc_cache_path(layer, band, grid_n, sigma)
    if cache.exists():
        with np.load(cache, allow_pickle=True) as d:
            return {k: d[k] for k in d.files}

    cfg = LAYERS[layer]
    obs_pairs = cfg["obs"](band)
    obs, null_median, null_p95 = [], [], []
    tie_obs, tie_surr, patients = [], [], []
    for pat in COHORT:
        dD_task, dD_rest = obs_pairs[pat]
        obs.append(_density_field(dD_task, dD_rest, grid_n, sigma))
        tie_obs.append(_pair_tie_frac(dD_task, dD_rest))

        dts, drs, _rhos = _surrogate_rho_and_pairs(pat, band, layer)
        surr_stack = np.stack(
            [_density_field(dt, dr, grid_n, sigma) for dt, dr in zip(dts, drs)],
            axis=0,
        )
        null_median.append(np.median(surr_stack, axis=0))
        null_p95.append(np.percentile(surr_stack, 95.0, axis=0))
        tie_surr.append(float(np.mean(
            [_pair_tie_frac(dt, dr) for dt, dr in zip(dts, drs)])))
        patients.append(pat)
        print(f"    [{layer}/{band}] {pat}: R={len(dts)} surrogates · "
              f"tie_obs={tie_obs[-1]:.3f} tie_surr={tie_surr[-1]:.3f}")

    out = dict(
        obs=np.stack(obs, axis=0),
        null_median=np.stack(null_median, axis=0),
        null_p95=np.stack(null_p95, axis=0),
        tie_obs=np.asarray(tie_obs),
        tie_surr=np.asarray(tie_surr),
        patients=np.asarray(patients),
    )
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache, **out)
    return out


def medexc_star(layer: str, band: str, floor: str,
                grid_n: int, sigma: float):
    """Cohort signed-excess field E_star[iu, iv] and its per-patient fields."""
    fields = compute_perpatient_fields(layer, band, grid_n, sigma)
    null = fields["null_median"] if floor == "median" else fields["null_p95"]
    E_p = fields["obs"] - null               # (P, g, g), ∫ ≈ 0 per patient
    E_star = np.median(E_p, axis=0)          # per-pixel cohort median
    return E_star, fields


def diagonal_mass(E_star: np.ndarray, strip: float = 0.1) -> tuple[float, float, float]:
    """Central-strip positive mass vs off-diagonal |mass| (E_star is [iu, iv])."""
    g = E_star.shape[0]
    centers = (np.arange(g) + 0.5) / g
    U, V = np.meshgrid(centers, centers, indexing="ij")
    duv = U - V
    cell = (1.0 / g) ** 2
    on = np.abs(duv) < strip
    diag = float(np.maximum(E_star[on], 0.0).sum() * cell)
    off = float(np.abs(E_star[~on]).sum() * cell)
    ratio = diag / off if off > 0 else float("nan")
    return diag, off, ratio


def to_sd_field(E_star: np.ndarray, n_s: int = 96, n_d: int = 192):
    """Rasterize E_star[iu, iv] onto straightened (s, d) coordinates.

    s = (u+v)/2, d = u−v. Outside the valid diamond |d| ≤ 2·min(s, 1−s) the
    field is NaN (rendered white). A trace (diagonal concentration) becomes a
    green horizontal band at d ≈ 0.
    """
    g = E_star.shape[0]
    centers = (np.arange(g) + 0.5) / g
    interp = RegularGridInterpolator(
        (centers, centers), E_star, bounds_error=False, fill_value=np.nan)
    s_ax = np.linspace(0.0, 1.0, n_s)
    d_ax = np.linspace(-1.0, 1.0, n_d)
    S, D = np.meshgrid(s_ax, d_ax)           # (n_d, n_s): rows d, cols s
    U = S + 0.5 * D
    V = S - 0.5 * D
    valid = (U >= 0.0) & (U <= 1.0) & (V >= 0.0) & (V <= 1.0)
    pts = np.stack([np.clip(U, 0.0, 1.0), np.clip(V, 0.0, 1.0)], axis=-1)
    Z = interp(pts.reshape(-1, 2)).reshape(D.shape)
    Z[~valid] = np.nan
    return Z, s_ax, d_ax


# ---------------------------------------------------------------------------
# One band panel — heatmap (precomputed signed_excess) + 2 conditional marginals
# ---------------------------------------------------------------------------
def plot_band_density(fig, outer_spec, band: str, t: np.ndarray, r: np.ndarray,
                      z_contour: np.ndarray, signed_excess: np.ndarray,
                      vabs: float, sup: str, draw_contours: bool = True):
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
    sub = outer_spec.subgridspec(
        2, 2, height_ratios=[1.0, 3.6], width_ratios=[3.6, 1.0],
        hspace=0.04, wspace=0.04,
    )
    ax_main = fig.add_subplot(sub[1, 0])
    ax_top = fig.add_subplot(sub[0, 0], sharex=ax_main)
    ax_right = fig.add_subplot(sub[1, 1], sharey=ax_main)

    # Band label in the otherwise-empty top-right corner cell.
    ax_label = fig.add_subplot(sub[0, 1])
    ax_label.set_xticks([]); ax_label.set_yticks([])
    for s in ax_label.spines.values():
        s.set_visible(False)
    ax_label.text(0.5, 0.5, band_tex, transform=ax_label.transAxes,
                  ha="center", va="center",
                  fontsize=22, fontweight="bold", color="0.18")

    ax_top.set_autoscalex_on(False)
    ax_right.set_autoscaley_on(False)
    ax_main.set_autoscalex_on(False)
    ax_main.set_autoscaley_on(False)
    ax_main.set_box_aspect(1.0)
    ax_top.set_box_aspect(1.0 / 3.6)
    ax_right.set_box_aspect(3.6)

    g = signed_excess.shape[0]
    xy = np.linspace(0.0, 1.0, g)
    xx, yy = np.meshgrid(xy, xy)

    im = ax_main.imshow(
        signed_excess, origin="lower", extent=[0.0, 1.0, 0.0, 1.0],
        cmap=TRACE_CMAP, vmin=-vabs, vmax=vabs,
        aspect="auto", interpolation="bilinear",
    )

    # Reference axes: main diagonal = trace axis (green here = trace),
    # anti-diagonal = reversal axis (the θ anti-trace reads against it). Their
    # crossing marks the independence centre.
    ax_main.plot([0, 1], [0, 1], color="white", lw=1.6, alpha=0.75, zorder=4)
    ax_main.plot([0, 1], [1, 0], color="white", lw=0.9, ls="--",
                 alpha=0.35, zorder=4)

    if draw_contours:
        # Nested HDR contours of the EMPIRICAL density — informative about the
        # data's shape. Skipped for the analytic copula, where the density is a
        # deterministic function of the single ρ (the contours add no
        # information, and being level sets of z_main rather than the plotted
        # signed-excess they would not even trace the coloured region).
        dA = 1.0 / (g * g)

        def _hdr_level(zf: np.ndarray, coverage: float) -> float:
            zs = np.sort(zf.ravel())[::-1]
            cm = np.cumsum(zs) * dA
            idx = min(int(np.searchsorted(cm, coverage)), zs.size - 1)
            return float(zs[idx])

        hdr_levels = sorted({_hdr_level(z_contour, c) for c in HDR_COVERAGES})
        ax_main.contour(xx, yy, z_contour, levels=hdr_levels,
                        colors="0.20", linewidths=1.0, alpha=0.40, zorder=4.5)
        ax_main.contour(xx, yy, z_contour, levels=hdr_levels,
                        colors="white", linewidths=0.5, alpha=0.85, zorder=5)
        ax_main.plot(0.5, 0.5, "o", markersize=5, color="white",
                     markeredgecolor="0.10", markeredgewidth=1.0, zorder=6)

    ax_main.set_xticks([0.0, 0.5, 1.0])
    ax_main.set_yticks([0.0, 0.5, 1.0])
    ax_main.set_xticklabels(["0", "½", "1"])
    ax_main.set_yticklabels(["0", "½", "1"])
    if band == BAND_ORDER[0]:
        ax_main.set_xlabel(r"$u$ = within-patient rank of "
                           rf"$\Delta_{{\mathrm{{task}}}}^{{\mathrm{{{sup}}}}}(i,j)$",
                           color="0.20")
        ax_main.set_ylabel(r"$v$ = within-patient rank of "
                           rf"$\Delta_{{\mathrm{{rest}}}}^{{\mathrm{{{sup}}}}}(i,j)$",
                           color="0.20")
    else:
        ax_main.set_xlabel(r"$u$", color="0.20")
        ax_main.set_ylabel(r"$v$", color="0.20")

    # Conditional-split marginals.
    n_bins = 28
    bin_edges = np.linspace(0.0, 1.0, n_bins + 1)
    mask_high_r = r > 0.5
    hist_top_hi, _ = np.histogram(t[mask_high_r], bins=bin_edges, density=True)
    hist_top_lo, _ = np.histogram(t[~mask_high_r], bins=bin_edges, density=True)
    yn_top = max(hist_top_hi.max(), hist_top_lo.max(), 1e-9)
    hist_top_hi_n = 0.85 * hist_top_hi / yn_top
    hist_top_lo_n = 0.85 * hist_top_lo / yn_top
    mask_high_t = t > 0.5
    hist_rt_hi, _ = np.histogram(r[mask_high_t], bins=bin_edges, density=True)
    hist_rt_lo, _ = np.histogram(r[~mask_high_t], bins=bin_edges, density=True)
    yn_rt = max(hist_rt_hi.max(), hist_rt_lo.max(), 1e-9)
    hist_rt_hi_n = 0.85 * hist_rt_hi / yn_rt
    hist_rt_lo_n = 0.85 * hist_rt_lo / yn_rt

    ax_top.stairs(hist_top_lo_n, bin_edges, baseline=0.0, fill=True,
                  facecolor="#c0392b", alpha=0.40,
                  edgecolor="#7f1d1d", linewidth=0.9, zorder=2)
    ax_top.stairs(hist_top_hi_n, bin_edges, baseline=0.0, fill=True,
                  facecolor="#1a7c3e", alpha=0.45,
                  edgecolor="#0a3a1d", linewidth=0.9, zorder=2)
    ax_top.axvline(0.5, color="0.35", lw=0.7, alpha=0.6, zorder=1)
    ax_top.set_ylim(0, 1)
    ax_top.tick_params(axis="both", which="both",
                       bottom=False, top=False, left=False, right=False,
                       labelbottom=False, labeltop=False,
                       labelleft=False, labelright=False)
    ax_top.spines[["top", "right", "left"]].set_visible(False)
    ax_top.spines["bottom"].set_color("0.70")
    if band == BAND_ORDER[0]:
        ax_top.text(0.03, 0.10, r"$P(u)$ split by $v$ half",
                    transform=ax_top.transAxes, ha="left", va="bottom",
                    color="0.20", fontstyle="italic",
                    bbox=dict(facecolor="white", alpha=0.65, edgecolor="none",
                              boxstyle="round,pad=0.25"), zorder=10)

    ax_right.stairs(hist_rt_lo_n, bin_edges, baseline=0.0, fill=True,
                    orientation="horizontal",
                    facecolor="#c0392b", alpha=0.40,
                    edgecolor="#7f1d1d", linewidth=0.9, zorder=2)
    ax_right.stairs(hist_rt_hi_n, bin_edges, baseline=0.0, fill=True,
                    orientation="horizontal",
                    facecolor="#1a7c3e", alpha=0.45,
                    edgecolor="#0a3a1d", linewidth=0.9, zorder=2)
    ax_right.axhline(0.5, color="0.35", lw=0.7, alpha=0.6, zorder=1)
    ax_right.set_xlim(0, 1)
    ax_right.tick_params(axis="both", which="both",
                         bottom=False, top=False, left=False, right=False,
                         labelbottom=False, labeltop=False,
                         labelleft=False, labelright=False)
    ax_right.spines[["top", "right", "bottom"]].set_visible(False)
    ax_right.spines["left"].set_color("0.70")
    if band == BAND_ORDER[0]:
        ax_right.text(0.22, 0.97, r"$P(v)$ split by $u$ half",
                      transform=ax_right.transAxes, ha="left", va="top",
                      rotation=270, rotation_mode="anchor",
                      color="0.20", fontstyle="italic",
                      bbox=dict(facecolor="white", alpha=0.65, edgecolor="none",
                                boxstyle="round,pad=0.15"),
                      clip_on=False, zorder=10)

    for spine in ax_main.spines.values():
        spine.set_edgecolor("0.55")
        spine.set_linewidth(0.9)

    ax_main.set_xlim(0.0, 1.0, auto=False)
    ax_main.set_ylim(0.0, 1.0, auto=False)
    return im


LAYERS = {
    "raw": dict(
        obs=raw_obs_pairs, drift=raw_drift_pairs,
        cohort_csv=RAW_MS_DIR / "cohort_summary_all_bands.csv",
        sup="raw", tag="rawfc",
    ),
    "coph": dict(
        obs=coph_obs_pairs, drift=coph_drift_pairs,
        cohort_csv=COPH_MS_DIR / "cohort_summary.csv",
        sup="coph", tag="coph",
    ),
}


def main(layer: str, mode: str) -> Path:
    cfg = LAYERS[layer]
    cohort = pd.read_csv(cfg["cohort_csv"])
    rng_obs = np.random.default_rng(TIEBREAK_SEED)
    rng_floor = np.random.default_rng(TIEBREAK_SEED + 1)

    pooled, z_contour, sx = {}, {}, {}
    for band in BAND_ORDER:
        obs_pairs = cfg["obs"](band)
        row = cohort[cohort.band == band].iloc[0]
        rho_obs = float(row["obs_median_rho"])
        rho_p95 = float(row["surr_median_rho_p95"])

        if mode == "copula":
            pooled[band] = pool_ranks(obs_pairs)
            z_main = gaussian_copula_density(rho_obs, GRID_N)
            z_floor = gaussian_copula_density(rho_p95, GRID_N)
            floor_label = "matched-strength p95"
        elif mode == "drift":
            pooled[band] = pool_ranks(obs_pairs)
            z_main = cohort_average_empirical_density(obs_pairs, GRID_N, KDE_BW)
            z_floor = cohort_average_empirical_density(
                cfg["drift"](band), GRID_N, KDE_BW)
            floor_label = "within-session drift"
        else:  # gate — robust empirical KDE + empirical matched-strength p95 floor
            pooled[band] = pool_ranks(obs_pairs, rng_obs)
            z_main = cohort_median_empirical_density(
                obs_pairs, GRID_N, KDE_BW, rng_obs)
            z_floor = cohort_median_p95_floor(
                band, layer, GRID_N, KDE_BW, rng_floor)
            floor_label = "matched-strength p95 (gate)"
        z_contour[band] = z_main
        sx[band] = signed_excess_field(z_main, z_floor)

        rhos = per_patient_rhos(obs_pairs)
        print(f"  [{layer}/{mode}] {band:11s}: rho_obs={rho_obs:+.3f} · "
              f"surr_p95={rho_p95:+.3f} · gate_p={float(row['paired_wilcoxon_p']):.3f}"
              f" · n_pos={int((rhos > 0).sum())}/{rhos.size}")

    vabs = (VABS_COPULA if mode == "copula"
            else float(np.percentile(np.abs(np.concatenate(
                [s.ravel() for s in sx.values()])), 99.0)))
    print(f"  [{layer}/{mode}] vabs={vabs:.4f}  floor={floor_label}")

    fig = plt.figure(figsize=(17.0, 11.0))
    outer = fig.add_gridspec(2, 3, wspace=0.20, hspace=0.34,
                             left=0.05, right=0.92, top=0.94, bottom=0.07)
    last_im = None
    for i, band in enumerate(BAND_ORDER):
        rr, cc = i // 3, i % 3
        last_im = plot_band_density(
            fig, outer[rr, cc], band, pooled[band][0], pooled[band][1],
            z_contour[band], sx[band], vabs, cfg["sup"],
            draw_contours=(mode != "copula"))

    cbar_ax = fig.add_axes([0.94, 0.16, 0.015, 0.62])
    cb = fig.colorbar(last_im, cax=cbar_ax, orientation="vertical", extend="both")
    cb.set_label(rf"signed enrichment over {floor_label} floor",
                 rotation=270, labelpad=18)

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"fig_bands_joint_density_{mode}_{cfg['tag']}.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  [{layer}/{mode}] Saved: {out}")
    return out


# ---------------------------------------------------------------------------
# median_excess renders + driver
# ---------------------------------------------------------------------------
def plot_band_sd(fig, ax, band: str, Z: np.ndarray, vabs: float) -> object:
    """Straightened (s, d) panel: trace = green horizontal band at d ≈ 0."""
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
    cmap = TRACE_CMAP.copy()
    cmap.set_bad("white")                          # outside the valid diamond
    im = ax.imshow(np.ma.masked_invalid(Z), origin="lower",
                   extent=[0.0, 1.0, -1.0, 1.0], cmap=cmap,
                   vmin=-vabs, vmax=vabs, aspect="auto",
                   interpolation="bilinear")
    ax.axhline(0.0, color="white", lw=1.4, alpha=0.80, zorder=4)   # trace axis
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-1.0, 1.0)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_xticklabels(["0", "½", "1"])
    ax.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    ax.set_yticklabels(["−1", "−½", "0", "½", "1"])
    for sp in ax.spines.values():
        sp.set_edgecolor("0.55")
        sp.set_linewidth(0.9)
    ax.set_box_aspect(1.0)
    if band == BAND_ORDER[0]:
        ax.set_xlabel(r"$s=(u+v)/2$", color="0.20")
        ax.set_ylabel(r"$d=u-v$  (trace at $d\simeq 0$)", color="0.20")
    else:
        ax.set_xlabel(r"$s$", color="0.20")
        ax.set_ylabel(r"$d$", color="0.20")
    ax.text(0.97, 0.05, band_tex, transform=ax.transAxes, ha="right",
            va="bottom", fontsize=22, fontweight="bold", color="0.18")
    return im


def main_medexc(layer: str, floor: str,
                grid_n: int = GRID_MEDEXC, sigma: float = SIGMA_MEDEXC):
    """Render the (u,v) and (s,d) median-excess figures for one layer; return
    the per-band scalar-mass rows and per-patient tie rows for the CSVs."""
    cfg = LAYERS[layer]
    E_stars, pooled, tie_rows = {}, {}, []
    for band in BAND_ORDER:
        E_star, fields = medexc_star(layer, band, floor, grid_n, sigma)
        E_stars[band] = E_star
        obs_pairs = cfg["obs"](band)
        pooled[band] = (
            np.concatenate([_midrank_unit(dt) for dt, _ in obs_pairs.values()]),
            np.concatenate([_midrank_unit(dr) for _, dr in obs_pairs.values()]),
        )
        for pat, to, ts in zip(fields["patients"], fields["tie_obs"],
                               fields["tie_surr"]):
            gap = float(abs(float(to) - float(ts)))
            tie_rows.append(dict(layer=layer, band=band, patient=str(pat),
                                 tie_frac_obs=float(to), tie_frac_surr=float(ts),
                                 gap=gap, flag=int(gap > 0.05)))

    # Symmetric colour limit L = p99(|E_star|) pooled over the six panels.
    vabs = float(np.percentile(
        np.concatenate([np.abs(E_stars[b]).ravel() for b in BAND_ORDER]), 99.0))
    print(f"  [{layer}/medexc/{floor}] L (colour limit) = {vabs:.5f}")

    mass_rows = []
    for band in BAND_ORDER:
        diag, off, ratio = diagonal_mass(E_stars[band])
        mass_rows.append(dict(layer=layer, band=band, floor=floor,
                              diagonal_excess=diag, offdiagonal_excess=off,
                              ratio=ratio, L=vabs))
        print(f"    {band:11s} diag+={diag:.5f} off|.|={off:.5f} ratio={ratio:.3f}")

    floor_label = ("matched-strength surrogate median" if floor == "median"
                   else "matched-strength surrogate p95")
    suffix = "" if floor == "median" else "_p95"
    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Render A: (u, v) — reuse the conditional-split-marginals panel ----
    figA = plt.figure(figsize=(17.0, 11.0))
    outerA = figA.add_gridspec(2, 3, wspace=0.20, hspace=0.34,
                               left=0.05, right=0.92, top=0.94, bottom=0.07)
    lastA = None
    for i, band in enumerate(BAND_ORDER):
        rr, cc = i // 3, i % 3
        lastA = plot_band_density(
            figA, outerA[rr, cc], band, pooled[band][0], pooled[band][1],
            None, E_stars[band].T, vabs, cfg["sup"], draw_contours=False)
    cbA = figA.add_axes([0.94, 0.16, 0.015, 0.62])
    cb = figA.colorbar(lastA, cax=cbA, orientation="vertical", extend="both")
    cb.set_label(rf"per-patient signed excess over {floor_label} (cohort median)",
                 rotation=270, labelpad=18)
    outA = out_dir / f"fig_bands_medexc_uv_{cfg['tag']}{suffix}.pdf"
    figA.savefig(outA, bbox_inches="tight", transparent=True)
    plt.close(figA)
    print(f"  [{layer}/medexc/{floor}] Saved: {outA}")

    # ---- Render B: (s, d) straightened ----
    figB = plt.figure(figsize=(16.0, 11.0))
    outerB = figB.add_gridspec(2, 3, wspace=0.24, hspace=0.30,
                               left=0.06, right=0.91, top=0.95, bottom=0.08)
    lastB = None
    for i, band in enumerate(BAND_ORDER):
        rr, cc = i // 3, i % 3
        ax = figB.add_subplot(outerB[rr, cc])
        Z, _s_ax, _d_ax = to_sd_field(E_stars[band])
        lastB = plot_band_sd(figB, ax, band, Z, vabs)
    cbB = figB.add_axes([0.925, 0.18, 0.014, 0.60])
    cb = figB.colorbar(lastB, cax=cbB, orientation="vertical", extend="both")
    cb.set_label(rf"per-patient signed excess over {floor_label} (cohort median)",
                 rotation=270, labelpad=22)
    outB = out_dir / f"fig_bands_medexc_sd_{cfg['tag']}{suffix}.pdf"
    figB.savefig(outB, bbox_inches="tight", transparent=True)
    plt.close(figB)
    print(f"  [{layer}/medexc/{floor}] Saved: {outB}")

    return mass_rows, tie_rows


def run_medexc(layers, floor: str):
    """Drive median_excess over *layers* and write the two companion CSVs."""
    MEDEXC_DATA.mkdir(parents=True, exist_ok=True)
    mass_all, tie_all = [], []
    for layer in layers:
        m, t = main_medexc(layer, floor)
        mass_all += m
        tie_all += t

    mass_df = pd.DataFrame(mass_all)
    suffix = "" if floor == "median" else "_p95"
    mass_path = MEDEXC_DATA / f"medexc_diagonal_mass{suffix}.csv"
    mass_df.to_csv(mass_path, index=False)
    print(f"  [medexc] wrote {mass_path}")

    # Tie fractions are floor-independent; write once (overwrite is idempotent).
    tie_df = pd.DataFrame(tie_all)
    coh = (tie_df.groupby(["layer", "band"], as_index=False)
           [["tie_frac_obs", "tie_frac_surr", "gap"]].mean())
    coh["patient"] = "cohort_mean"
    coh["flag"] = (coh["gap"] > 0.05).astype(int)
    tie_df = pd.concat([tie_df, coh[tie_df.columns]], ignore_index=True)
    tie_path = MEDEXC_DATA / "medexc_tie_fractions.csv"
    tie_df.to_csv(tie_path, index=False)
    print(f"  [medexc] wrote {tie_path}")
    flagged = tie_df[(tie_df.patient != "cohort_mean") & (tie_df.flag == 1)]
    if len(flagged):
        print(f"  [medexc] WARNING: {len(flagged)} (band,patient) tie-gap > 0.05 "
              f"— cancellation incomplete there:")
        for _, rrow in flagged.iterrows():
            print(f"      {rrow.layer}/{rrow.band}/{rrow.patient}: "
                  f"|{rrow.tie_frac_obs:.3f}−{rrow.tie_frac_surr:.3f}|"
                  f"={rrow.gap:.3f}")
    else:
        print("  [medexc] tie cancellation OK: all gaps ≤ 0.05")
    return mass_path, tie_path


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--layer", choices=list(LAYERS), default="raw")
    ap.add_argument("--mode",
                    choices=["copula", "drift", "gate", "median_excess"],
                    default="copula")
    ap.add_argument("--floor", choices=["median", "p95"], default="median",
                    help="median_excess surrogate floor: median (default, "
                         "'excess over typical null') or p95 (gate-faithful)")
    ap.add_argument("--both-layers", action="store_true",
                    help="render raw and coph for the chosen --mode")
    ap.add_argument("--all", action="store_true",
                    help="render all density modes (copula/drift/gate) × layers")
    args = ap.parse_args()
    if args.mode == "median_excess":
        layers = list(LAYERS) if (args.both_layers or args.all) else [args.layer]
        run_medexc(layers, args.floor)
    elif args.all:
        for m in ("copula", "drift", "gate"):
            for lyr in ("raw", "coph"):
                main(lyr, m)
    elif args.both_layers:
        for lyr in ("raw", "coph"):
            main(lyr, args.mode)
    else:
        main(args.layer, args.mode)
