#!/usr/bin/env python3
"""β per-pair persistence figure — TEST2 design.

Two landscape panels (the brain-anatomy panel was deferred):

Panel (a)  per-patient z above own matched-strength surrogate (kept
           from v1, polished — full Pat_XX labels, micro-bar for
           ``n_above_p95_local``, stems coloured by trace sign).

Panel (b)  exemplar (default Pat_05) trace + anti edges embedded in
           the LRG hierarchy-bundled chord layout (graph-tool
           ``get_hierarchy_control_points`` over an
           ``fcluster(linkage, K, 'maxclust')`` community cut of the
           exemplar's β rest_pre LRG).  Reveals the modular topology
           of the preserved-shift subgraph; epi-zone contacts are
           red in the outer label ring.

Usage
-----
    python preprint_07_beta_rho_split_figure_test2.py             # beta, Pat_05
    python preprint_07_beta_rho_split_figure_test2.py --band alpha
    python preprint_07_beta_rho_split_figure_test2.py --exemplar Pat_03

Inputs
------
data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv
data/reports/imcoh_continuous_trace/per_pair_split/<Pat>_<band>.npz
data/cache/imcoh/<Pat>/<band>_rest_pre_imcoh_freqresolved_nperseg-*.npy   (via load_fc_matrix)
data/cache/imcoh_lrg/<Pat>/<band>_rest_pre_lrg_imcoh-abs.npz              (via load_lrg_result)
data/raw/stereoeeg_patients/<Pat>/implant_pat_NN.csv                       (MNI coords)
data/raw/stereoeeg_patients/<Pat>/channel_labels.csv                       (shaft labels)

Output (PDF only, no PNG sibling)
---------------------------------
data/preprint/figures/<band>/per_pair_trace/fig_<band>_rho_split_test2.pdf
"""
from __future__ import annotations

import argparse
import math
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.image import imread
from scipy.cluster.hierarchy import cophenet, fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata, spearmanr

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.utils.io.patient import load_epileptic_nodes, load_patient_metadata
from lrg_eegfc.utils.io.regions import _normalise_label
from lrg_eegfc.visuals.network_templates import (
    _load_inputs,
    shaft_colors,
)
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result

ROOT = setup_script_env()
use_lrg_style()


COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]

LRG_CTM_DIR = ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"
CTM_TRI_DIR = ROOT / "data" / "audit" / "ctm_triangle"
PAIR_SPLIT_DIR = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "per_pair_split"

# Cached surrogate eigendecomps for the §5.2 / §5.3 matched-strength
# pipeline.  One file per (patient, band, phase) at
# {SURR_LRG_CACHE}/Pat_NN/{band}_{phase}_R200_swap20_seed20260511_imcoh_abs.npz
# containing `eigvals` (R×N) and `eigvecs` (R×N×N).
SURR_LRG_CACHE = ROOT / "data" / "cache" / "matched_strength_surrogate_lrg"
# Per-pair p-value cache (this script writes here on first run, then
# reuses).  One file per (patient, band) keyed by surrogate params.
PER_PAIR_CACHE = ROOT / "data" / "cache" / "matched_strength_per_pair"
SURROGATE_R = 200
SURROGATE_SWAP = 20
SURROGATE_SEED = 20260511

# Trace-direction palette: warm natural pairing (grass + brick) so the
# trace/anti edges read against the black null backbone in panel (b)
# as found-object colours, not synthetic CSS swatches.
CLR_TRACE = "#5fa844"   # grass green  — trace concordant
CLR_ANTI = "#c0392b"    # brick red    — anti concordant
CLR_NODE = "#2c3e50"    # node default for panel (c) when shaft palette off

# Edge classification gate.  Pair is classified trace (green) iff
# its matched-strength empirical two-sided p ≤ ALPHA AND c_obs > 0;
# anti (red) iff p ≤ ALPHA AND c_obs < 0.  No multiple-testing
# correction across pairs — the 6.9k pairs are the per-pair
# decomposition of a single patient's verdict, not 6.9k independent
# claims.  The per-patient gate is panel (a) §5.3 cohort+per-patient
# Wilcoxon, already published.
ALPHA_DEFAULT = 0.05
# Restrict the rendered edge set to cross-probe pairs (probes[i] !=
# probes[j]).  Same-shaft pairs are dominated by zero-phase-lag
# leakage (probe-bias guide); the §3.2 cross-probe restriction is
# the matching robustness check in the manuscript text.
CROSS_PROBE_ONLY_DEFAULT = True
# Phase whose ImCoh|·| sets the edge magnitude shown in panel (b).
# Default rest_post = "what the FC looks like after task"; can flip
# to rest_pre or task_test via CLI.
IMCOH_PHASE_DEFAULT = "rest_post"

# graph-tool hierarchy-bundled chord (panel c) settings — mirrors the
# `render_lrg_panel` recipe in scripts/01_compute/figures_for_notes/_shared.py
# which produced data/outputs/figures/section2/fig_E/fig_E2_sbm_*.pdf.
# K_CLUSTERS is the depth-2 cut of the LRG dendrogram (graph-tool's
# get_hierarchy_control_points caps at depth 2); leaves are placed on
# the outer ring in scipy `leaves_list(Z)` order so same-cluster nodes
# are angularly adjacent.  BETA controls bezier bundle tightness:
# 0 = chord, 1 = fully bundled through the cluster centres.
LRG_K_CLUSTERS = 7      # fewer = larger bundles (was 12 → bundles too thin)
LRG_BETA = 0.92         # higher = tighter bundling through cluster centres
LRG_RENDER_PX = 2200    # bigger raster for crisp embedding


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
def load_pair_data(pat: str, band: str) -> dict:
    npz_path = PAIR_SPLIT_DIR / f"{pat}_{band}.npz"
    d = np.load(npz_path)
    return dict(
        dD_task=np.asarray(d["dD_task"]),
        dD_rest=np.asarray(d["dD_rest"]),
        sigma=np.asarray(d["sigma"]).astype(int),
        iu_i=np.asarray(d["iu_i"]).astype(int),
        iu_j=np.asarray(d["iu_j"]).astype(int),
    )


def per_contact_displayed_trace_weight(pair: dict, edges: dict,
                                        weight_per_pair: np.ndarray,
                                        N: int) -> np.ndarray:
    """Per-contact incidence weight in the DISPLAYED green edges.

    w_trace(c) = Σ ``weight_per_pair[k]`` over the σ=+1 cross-probe
    pairs ``k`` that touch contact ``c`` AND clear the matched-strength
    FDR gate (i.e. that ended up in ``edges['trace_idx']``).
    Used only for the radial-sort ranking.
    """
    trace_idx = np.asarray(edges["trace_idx"], dtype=int)
    if trace_idx.size == 0:
        return np.zeros(N)
    iu_i = pair["iu_i"][trace_idx]
    iu_j = pair["iu_j"][trace_idx]
    w_k = np.asarray(weight_per_pair)[trace_idx]
    in_range = (iu_i < N) & (iu_j < N)

    w = np.zeros(N)
    np.add.at(w, iu_i[in_range], w_k[in_range])
    np.add.at(w, iu_j[in_range], w_k[in_range])
    return w


def select_significant_edges(pair: dict, p_pair: np.ndarray,
                              c_obs: np.ndarray,
                              probes: list[str], alpha: float,
                              cross_probe_only: bool = True) -> dict:
    """Classify each pair as trace / anti / null using the per-pair
    matched-strength two-sided p-value at threshold ``alpha``.

    Pair is **trace** (green) iff ``p_pair ≤ alpha`` AND ``c_obs > 0``
    (pair sits in the same tail of dD_task and dD_rest);  **anti**
    (red) iff ``p_pair ≤ alpha`` AND ``c_obs < 0`` (opposite tails).
    Everything else falls into the grey null.

    No multiple-testing correction across pairs.  The 6.9k pairs
    are the per-pair decomposition of one patient's signal — the
    per-patient gate is the §5.3 cohort+per-patient Wilcoxon
    summarised in panel (a) — so each pair simply needs to clear
    its own matched-strength null at α.
    """
    iu_i = pair["iu_i"]
    iu_j = pair["iu_j"]
    N_probes = len(probes)

    in_range = (iu_i < N_probes) & (iu_j < N_probes)
    if cross_probe_only:
        probes_arr = np.asarray(probes)
        xp = np.zeros(len(iu_i), dtype=bool)
        xp[in_range] = (probes_arr[iu_i[in_range]]
                        != probes_arr[iu_j[in_range]])
        keep = in_range & xp
    else:
        keep = in_range

    sig = (p_pair <= alpha) & keep
    trace_idx = np.where(sig & (c_obs > 0))[0]
    anti_idx = np.where(sig & (c_obs < 0))[0]

    return dict(
        trace_idx=trace_idx,
        anti_idx=anti_idx,
        n_tested=int(keep.sum()),
    )


# ---------------------------------------------------------------------------
# Per-pair matched-strength null
# ---------------------------------------------------------------------------
def _surrogate_d_coph(eigvals_r: np.ndarray,
                      eigvecs_r: np.ndarray) -> np.ndarray | None:
    """Reconstruct one surrogate's cophenet D from its cached Laplacian
    eigendecomp.  Returns ``None`` if the shuffle failed (NaN entries)
    or if τ_max/ρ̂ degenerate.  Mirrors ``lrg_linkage`` in
    ``scripts/01_compute/audit/audit_65_kc_matched_strength_surrogate``.
    """
    if not np.isfinite(eigvals_r).all() or not np.isfinite(eigvecs_r).all():
        return None
    lam_max = float(eigvals_r[-1])
    if lam_max <= 0:
        return None
    tau = 1.0 / lam_max
    diag_exp = np.exp(-tau * eigvals_r)
    rho = (eigvecs_r * diag_exp) @ eigvecs_r.T
    tr = float(np.trace(rho))
    if not np.isfinite(tr) or tr <= 0:
        return None
    rho = rho / tr
    with np.errstate(divide="ignore", invalid="ignore"):
        D = 1.0 / rho
    D = np.maximum(D, D.T)
    np.fill_diagonal(D, 0.0)
    finite = np.isfinite(D)
    if not finite.all():
        if not finite.any():
            return None
        cap = float(np.nanmax(D[finite]))
        D = np.where(finite, D, cap)
    try:
        Z_surr = linkage(squareform(D, checks=False), method="average")
        D_coph = squareform(cophenet(Z_surr))
    except Exception:
        return None
    return D_coph


def _per_pair_concordance(d_task: np.ndarray,
                           d_rest: np.ndarray) -> np.ndarray:
    """Per-pair contribution to Spearman ρ between ``d_task`` and
    ``d_rest``.

    .. math::

        c_k = \\bigl(\\mathrm{rank}_{\\!t}(k) - \\bar r\\bigr)
              \\bigl(\\mathrm{rank}_{\\!r}(k) - \\bar r\\bigr),
        \\qquad \\bar r = \\tfrac{n+1}{2}

    Summed over ``k``, this gives the (unnormalised) Spearman
    numerator — so the per-pair contribution is the natural rank-based
    statistic for "does this pair drive ρ_split^coph".  Sign tells
    direction: ``c > 0`` ⇔ pair sits in the same tail of both phases
    (trace); ``c < 0`` ⇔ opposite tails (anti).
    """
    n = len(d_task)
    rt = rankdata(d_task)
    rr = rankdata(d_rest)
    return (rt - (n + 1) * 0.5) * (rr - (n + 1) * 0.5)


def compute_or_load_per_pair_pvals(patient: str, band: str, *,
                                    n_surr: int = SURROGATE_R,
                                    swap_factor: int = SURROGATE_SWAP,
                                    seed: int = SURROGATE_SEED,
                                    fc_method: str = "imcoh_abs",
                                    verbose: bool = True) -> dict:
    """Empirical two-sided p-values for the per-pair Spearman-ρ
    contribution ``c_k = (rank_t - (n+1)/2)·(rank_r - (n+1)/2)`` against
    the R=200 matched-strength surrogate null.

    A two-sided gate on ``|c|`` is the right test against this null:
    matched-strength inflates per-phase ``|Δ_D_coph|`` magnitudes by
    construction (random shuffles → wider tree spread), so per-pair
    magnitude tests are biased the wrong way.  The §5.3 cohort test
    uses Spearman precisely because it is **rank-based** and immune to
    this magnitude inflation.

    P-value is the **empirical** count
    ``(1 + #{|c_surr| ≥ |c_obs|}) / (1 + R)`` — minimum is 1/(R+1).
    No multiple-testing correction across pairs; the per-patient gate
    is panel (a) §5.3 cohort+per-patient Wilcoxon.

    Cache: ``data/cache/matched_strength_per_pair/<Pat>/
    {band}_R{R}_swap{S}_seed{seed}_{fc_method}.npz``.  Pairing of
    surrogates across phases is by index — same convention as §5.3.
    """
    cache_dir = PER_PAIR_CACHE / patient
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = (
        cache_dir
        / f"{band}_R{n_surr}_swap{swap_factor}_seed{seed}_{fc_method}.npz"
    )
    if cache_path.exists():
        d = np.load(cache_path)
        out = dict(
            iu_i=np.asarray(d["iu_i"]).astype(int),
            iu_j=np.asarray(d["iu_j"]).astype(int),
            p_pair=np.asarray(d["p_pair"]).astype(float),
            c_obs=np.asarray(d["c_obs"]).astype(float),
            valid_R=int(d["valid_R"]),
        )
        if verbose:
            print(f"  loaded per-pair p-values from {cache_path.name} "
                  f"(valid R={out['valid_R']})")
        return out

    # Observed per-pair Spearman-ρ contribution from per_pair_split NPZ
    obs = load_pair_data(patient, band)
    iu_i = obs["iu_i"]; iu_j = obs["iu_j"]
    c_obs = _per_pair_concordance(obs["dD_task"], obs["dD_rest"])

    # Load surrogate eigendecomps for the four phases used by the
    # split-baseline per-pair NPZ:
    #     Δ_task = D_test  − D_pre_A      (baseline = first half)
    #     Δ_rest = D_post  − D_pre_B      (baseline = second half)
    # — matches `_per_cell_split` in
    # ``scripts/01_compute/hypothesis_tests/continuous_trace_matrix.py``.
    phases = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
    surr_paths = {
        ph: (SURR_LRG_CACHE / patient
             / f"{band}_{ph}_R{n_surr}_swap{swap_factor}_seed{seed}_{fc_method}.npz")
        for ph in phases
    }
    surr_eigs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for ph, p in surr_paths.items():
        if not p.exists():
            raise FileNotFoundError(
                f"Surrogate eigendecomp cache missing for {patient}/{band}/{ph}: {p}"
            )
        with np.load(p) as data:
            surr_eigs[ph] = (np.asarray(data["eigvals"]).astype(np.float64),
                              np.asarray(data["eigvecs"]).astype(np.float64))

    R = int(surr_eigs["rest_pre_A"][0].shape[0])
    N_full = int(surr_eigs["rest_pre_A"][0].shape[1])
    pair_count = len(iu_i)

    # Accumulator: surrogate per-pair concordance values.
    c_surr_all = np.full((R, pair_count), np.nan, dtype=np.float64)
    if verbose:
        print(f"  computing per-pair null for {patient}/{band} "
              f"(R={R}, N={N_full}, pairs={pair_count})")
    valid_R = 0
    for r in range(R):
        D_phase: dict[str, np.ndarray] = {}
        ok = True
        for ph in phases:
            evals, evecs = surr_eigs[ph]
            D_c = _surrogate_d_coph(evals[r], evecs[r])
            if D_c is None:
                ok = False
                break
            D_phase[ph] = D_c
        if not ok:
            continue
        d_task = (D_phase["task_test"] - D_phase["rest_pre_A"])[iu_i, iu_j]
        d_rest = (D_phase["rest_post"] - D_phase["rest_pre_B"])[iu_i, iu_j]
        c_surr_all[r, :] = _per_pair_concordance(d_task, d_rest)
        valid_R += 1

    if valid_R == 0:
        raise RuntimeError(
            f"All R={R} surrogates failed for {patient}/{band}"
        )
    valid_mask = np.isfinite(c_surr_all).all(axis=1)
    c_surr_valid = c_surr_all[valid_mask]
    n_above = (np.abs(c_surr_valid) >= np.abs(c_obs)[None, :]).sum(axis=0)
    p_pair = (1 + n_above) / (1 + valid_R)

    np.savez_compressed(
        cache_path,
        iu_i=iu_i, iu_j=iu_j,
        p_pair=p_pair, c_obs=c_obs,
        valid_R=valid_R, R=R, seed=seed, swap_factor=swap_factor,
    )
    if verbose:
        print(f"  cached per-pair p-values → {cache_path.name} "
              f"(valid R={valid_R}, empirical two-sided)")
    return dict(
        iu_i=iu_i, iu_j=iu_j,
        p_pair=p_pair, c_obs=c_obs,
        valid_R=valid_R,
    )


# ---------------------------------------------------------------------------
# Panel (a) — per-patient z-stems
# ---------------------------------------------------------------------------
def short(pat: str) -> str:
    """Compact label P02, P03, ..."""
    return pat.replace("Pat_", "P")


def panel_a_zstrip(ax: plt.Axes, per_pat: pd.DataFrame, cohort_row: pd.Series,
                   band: str, band_tex: str) -> None:
    """Per-patient z above own matched-strength surrogate (R=200,
    4-cycle ±δ).  Stems sorted by z descending; green if z > +2,
    red if z < −2, grey otherwise.  Each patient gets a thin
    right-edge ``frac_above_surr_local`` micro-bar that converts the
    cohort ``n_above_p95`` summary into a per-row visual.
    """
    sub = per_pat[per_pat.band == band].set_index("patient").loc[COHORT]
    z = sub["obs_z"].values
    # Fraction of surrogates the patient beat — derived from
    # ``obs_p_one_sided`` (= proportion of R=200 surrogates at-or-above
    # the observed ρ).  Robust to column-rename variants.
    if "obs_p_one_sided" in sub.columns:
        frac = 1.0 - sub["obs_p_one_sided"].values
    elif "obs_p_emp" in sub.columns:
        frac = 1.0 - sub["obs_p_emp"].values
    elif "frac_above_surr" in sub.columns:
        frac = sub["frac_above_surr"].values
    else:
        frac = np.full(len(sub), np.nan)

    order = np.argsort(z)[::-1]
    patients_sorted = [COHORT[i] for i in order]
    z_sorted = z[order]
    frac_sorted = frac[order]

    n = len(z_sorted)
    y = np.arange(n)[::-1]

    ax.axvspan(-2, 2, color="0.92", alpha=0.55, zorder=0)
    ax.axvline(0, color="0.35", lw=0.8, zorder=1)

    for yi, zi in zip(y, z_sorted):
        if zi > 2:
            clr = CLR_TRACE
        elif zi < -2:
            clr = CLR_ANTI
        else:
            clr = "0.55"
        ax.plot([0, zi], [yi, yi], color=clr, lw=1.9, alpha=0.85,
                zorder=2, solid_capstyle="round")
        ax.scatter([zi], [yi], s=110, color=clr, edgecolor="white",
                   linewidth=1.2, zorder=3)

    ax.set_yticks(y)
    ax.set_yticklabels(patients_sorted, fontsize=9.5)
    ax.set_ylim(-0.6, n - 0.4)

    z_max = max(float(z_sorted.max()) + 1.2, 3.5)
    z_min = min(float(z_sorted.min()) - 1.0, -3.0)
    ax.set_xlim(z_min, z_max)
    ax.set_xlabel(r"$z$ above own matched-strength surrogate",
                  fontsize=10)

    if np.isfinite(frac_sorted).any():
        # micro-bar at the right edge: width ∝ frac, height fixed.
        x0 = z_max - 0.35 * (z_max - z_min) * 0.04
        bar_max_w = (z_max - z_min) * 0.08
        for yi, fr in zip(y, frac_sorted):
            if not np.isfinite(fr):
                continue
            fr_clip = float(np.clip(fr, 0.0, 1.0))
            bar_clr = CLR_TRACE if fr_clip >= 0.95 else "0.55"
            ax.add_patch(plt.Rectangle(
                (x0 - bar_max_w, yi - 0.18),
                bar_max_w * fr_clip, 0.36,
                facecolor=bar_clr, edgecolor="none",
                alpha=0.75, zorder=2,
            ))
            ax.add_patch(plt.Rectangle(
                (x0 - bar_max_w, yi - 0.18),
                bar_max_w, 0.36,
                facecolor="none", edgecolor="0.65",
                lw=0.5, zorder=2.5,
            ))

    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("(a)", fontsize=12.0, loc="left", pad=4, fontweight="bold")


# ---------------------------------------------------------------------------
# Patient label helper
# ---------------------------------------------------------------------------
def _pat_display(pat: str) -> str:
    """Display label for a patient — matplotlib's plain text renderer
    handles ``_`` fine outside mathtext, so no escape is needed.
    Returns ``Pat_05`` verbatim.
    """
    return pat


# ---------------------------------------------------------------------------
# Panel (c) — LRG hierarchy-bundled chord network (graph-tool)
# ---------------------------------------------------------------------------
RGB_BACKBONE = (0.55, 0.72, 0.92)   # soft sky blue — hierarchy tree
RGB_NULL = (0.0, 0.0, 0.0)          # pure black — non-trace pairs


def _parse_contact_label(raw: str) -> str:
    """Compact ``<letter><number>`` label.

    ``"F 1,G2"`` → ``"F1"``;  ``"G' 4,X1"`` → ``"G4"``;
    ``"Pì 12,..."`` (UTF-8-mojibake of ``P' 12``) → ``"P12"``.

    The pipeline: take the head before the first comma (drops the
    "groups" annotation), then keep only ASCII alphanumeric characters
    so apostrophes, spaces and mojibake garbage disappear.  Primed and
    unprimed shafts (e.g. ``G`` and ``G'``) collide on the resulting
    label by design; the shaft colour disambiguates them visually.
    """
    head = str(raw).split(",")[0]
    return "".join(c for c in head if c.isascii() and c.isalnum())


def load_contact_labels(patient: str) -> list[str]:
    """Per-contact display labels (``F 1``, ``G' 4``, …) in the same
    order as ``load_probe_labels``.  Falls back to numeric ids if the
    metadata is unavailable.
    """
    meta = load_patient_metadata(patient, SEEG_DATAPATH)
    if meta is None or "label" not in meta.columns:
        return []
    return [_parse_contact_label(s) for s in meta["label"].tolist()]


def render_lrg_hierarchy_chord(out_path: Path, lrg, pair: dict, probes: list[str],
                                edges: dict, contact_labels: list[str],
                                imcoh_weights: np.ndarray, *,
                                k_clusters: int = LRG_K_CLUSTERS,
                                beta: float = LRG_BETA,
                                output_size: int = LRG_RENDER_PX,
                                fit_view: float = 0.92) -> dict:
    """Render the trace+anti subgraph as a hierarchy-bundled chord —
    same recipe as ``render_lrg_panel`` in
    ``scripts/01_compute/figures_for_notes/_shared.py`` — with:

    1.  every non-selected pair drawn underneath as a faint grey
        bundled edge — the "null backbone" that gives the trace +
        anti edges their context;
    2.  larger leaf nodes coloured by shaft;
    3.  contact labels (``F 1``, ``G' 4``, …) overlaid in matplotlib
        on the outer ring, rotated radially.
    """
    import graph_tool.all as gt

    Z = np.asarray(lrg.linkage_matrix)
    M = int(Z.shape[0]) + 1
    N = min(M, len(probes))

    # ------------------------------------------------------------------
    # 1.  Build the depth-2 hierarchy tree (root → clusters → leaves).
    # Cluster MEMBERSHIP comes from the LRG dendrogram (``fcluster``
    # cut of ``Z``); the radial ORDER is hierarchical-by-trace:
    #
    #   • clusters are placed around the ring in order of MEAN
    #     per-contact trace incidence (highest-trace cluster near
    #     12 o'clock);
    #   • within each cluster, leaves are placed by their own per-
    #     contact trace incidence, descending.
    #
    # This keeps each cluster's leaves angularly contiguous so the
    # cluster centroid (vector mean of leaf angles) sits BETWEEN its
    # leaves — bundled edges then exit each leaf at a consistent angle
    # toward the centroid, giving the clean chord aesthetic.
    # ------------------------------------------------------------------
    k = int(min(max(2, k_clusters), N - 1))
    labels = fcluster(Z, t=k, criterion="maxclust")[:N]

    w_trace = per_contact_displayed_trace_weight(
        pair, edges, imcoh_weights, N,
    )

    # Group leaves by LRG-cluster label; sort each cluster's leaves by
    # trace incidence; then sort the clusters themselves by their MEAN
    # trace incidence so the trace-rich cluster sits near the top.
    cluster_leaves: dict = {}
    for leaf in range(N):
        cluster_leaves.setdefault(int(labels[leaf]), []).append(leaf)
    for lab in cluster_leaves:
        cluster_leaves[lab] = sorted(
            cluster_leaves[lab], key=lambda l: -w_trace[l]
        )

    def _cluster_score(lab: int) -> float:
        leaves = cluster_leaves[lab]
        return float(np.mean([w_trace[l] for l in leaves])) if leaves else 0.0

    ordered_labels = sorted(cluster_leaves.keys(),
                            key=lambda x: -_cluster_score(x))
    radial_order: list[int] = []
    for lab in ordered_labels:
        radial_order.extend(cluster_leaves[lab])

    n_cl = len(ordered_labels)
    t = gt.Graph(directed=True)
    t.add_vertex(N + n_cl + 1)
    root_vid = N + n_cl
    for j, lab in enumerate(ordered_labels):
        cluster_vid = N + j
        t.add_edge(t.vertex(root_vid), t.vertex(cluster_vid))
        for leaf in cluster_leaves[lab]:
            t.add_edge(t.vertex(cluster_vid), t.vertex(leaf))
    root = t.vertex(root_vid)

    # Custom radial layout: leaves equally spaced on outer ring (R=1)
    # in trace-score-descending order, cluster centroids at R=0.55
    # (vector mean of their leaves' angles), root at origin.
    R_OUTER = 1.0
    R_CLUSTER = 0.55
    n_leaves_in_order = len(radial_order)
    theta_of_leaf = np.zeros(N)
    for pos_in_order, leaf_id in enumerate(radial_order):
        # Subtract pi/2 so the slot-0 leaf (highest s) sits at 12 o'clock.
        theta_of_leaf[leaf_id] = (
            2 * math.pi * pos_in_order / n_leaves_in_order - math.pi / 2
        )

    tpos = t.new_vertex_property("vector<double>")
    for leaf_id in range(N):
        th = float(theta_of_leaf[leaf_id])
        tpos[t.vertex(leaf_id)] = [
            R_OUTER * math.cos(th), R_OUTER * math.sin(th),
        ]
    for j, lab in enumerate(ordered_labels):
        sum_cos = sum(math.cos(theta_of_leaf[l]) for l in cluster_leaves[lab])
        sum_sin = sum(math.sin(theta_of_leaf[l]) for l in cluster_leaves[lab])
        avg_th = math.atan2(sum_sin, sum_cos)
        tpos[t.vertex(N + j)] = [
            R_CLUSTER * math.cos(avg_th), R_CLUSTER * math.sin(avg_th),
        ]
    tpos[t.vertex(root_vid)] = [0.0, 0.0]

    # ------------------------------------------------------------------
    # 2.  Build the rendering graph ``g`` — leaf vertices only.  The
    # hierarchy backbone (root + cluster T-junctions + spokes) is
    # rendered separately as a matplotlib overlay on the imshow'd PNG,
    # because adding root/cluster edges to ``g`` confuses graph-tool's
    # ``get_hierarchy_control_points`` (directed tree-path lookup
    # fails on the backbone single-step edges).
    # ------------------------------------------------------------------
    g = gt.Graph(directed=False)
    g.add_vertex(N)

    sigma_p = g.new_edge_property("int")
    mag_p = g.new_edge_property("double")
    cat_p = g.new_edge_property("int")  # 0=null, +1=trace, -1=anti

    pair_i = pair["iu_i"]
    pair_j = pair["iu_j"]
    sigma_arr = pair["sigma"]
    weight_all = np.asarray(imcoh_weights)
    trace_set = set(int(x) for x in edges["trace_idx"])
    anti_set = set(int(x) for x in edges["anti_idx"])
    probes_arr = np.asarray(probes)

    # Cross-probe-only rendering: the classification family is
    # cross-probe (§3.2 probe-bias guide), so the visual must match.
    # Including same-shaft pairs would inflate ``mag_max`` (intra-probe
    # ImCoh|·| is 2-8× cross-probe ImCoh|·|), crushing every cross-probe
    # null edge into the invisible tail of the width scaling.
    n_pair = len(pair_i)
    for kk in range(n_pair):
        i, j = int(pair_i[kk]), int(pair_j[kk])
        if i >= N or j >= N or i == j:
            continue
        if probes_arr[i] == probes_arr[j]:
            continue  # same-shaft pair → not in cross-probe family
        e = g.add_edge(i, j)
        if kk in trace_set:
            cat_p[e] = 1
        elif kk in anti_set:
            cat_p[e] = -1
        else:
            cat_p[e] = 0
        sigma_p[e] = int(sigma_arr[kk])
        # Edge magnitude now driven by rest_post ImCoh|·| (the actual
        # FC weight of this connection), not by |Δ_D_coph|·|Δ_D_coph|.
        # Color carries the LRG-measure trace classification.
        mag_p[e] = float(weight_all[kk])

    # ------------------------------------------------------------------
    # 3.  Layout — copy tpos into the rendering graph's pos property.
    # ------------------------------------------------------------------
    pos = g.new_vertex_property("vector<double>")
    for v in g.vertices():
        tv = t.vertex(int(v))
        pos[v] = [float(tpos[tv][0]), float(tpos[tv][1])]

    cts = gt.get_hierarchy_control_points(g, t, tpos, beta=beta)

    # ------------------------------------------------------------------
    # 4.  Edge appearance — null grey backbone first, then anti, then
    #     trace on top, with the hierarchy backbone painted just above
    #     the null so the T-junctions read through the grey wash.
    # ------------------------------------------------------------------
    rgb_trace = matplotlib.colors.to_rgb(CLR_TRACE)
    rgb_anti = matplotlib.colors.to_rgb(CLR_ANTI)
    edge_list = list(g.edges())

    # Edge appearance — width/alpha scale with rest_post ImCoh|·|
    # (raw, not rank).  Null edges are black with width and alpha
    # both modulated by the per-pair FC weight; weak FC pairs
    # disappear via the superlinear (^1.4) width exponent, surviving
    # null edges are dark, sublinear-alpha (^0.5) so the moderately
    # strong ones are clearly opaque without being uniformly black.
    all_mags = [float(mag_p[e]) for e in edge_list]
    mag_max = max(all_mags) if all_mags else 1.0

    ecolor = g.new_edge_property("vector<double>")
    pen = g.new_edge_property("double")
    draw_order = g.new_edge_property("double")

    for e in edge_list:
        cat = int(cat_p[e])
        m = float(mag_p[e])
        t_w_null = (m / max(mag_max, 1e-12)) ** 1.4
        t_w_sig = (m / max(mag_max, 1e-12)) ** 0.45
        t_a = (m / max(mag_max, 1e-12)) ** 0.5
        if cat == 0:   # null backbone — black, width carries FC pattern
            alpha = 0.30 + t_a * 0.65
            width = 0.02 + t_w_null * 3.0
            ecolor[e] = [RGB_NULL[0], RGB_NULL[1], RGB_NULL[2], alpha]
            pen[e] = width
            draw_order[e] = m
        else:          # +1 (trace) or -1 (anti)
            alpha = 0.55 + t_w_sig * 0.40
            width = 1.6 + t_w_sig * 6.4
            rgb = rgb_trace if cat > 0 else rgb_anti
            ecolor[e] = [rgb[0], rgb[1], rgb[2], alpha]
            pen[e] = width
            draw_order[e] = (2000.0 + m) if cat > 0 else (1000.0 + m)

    # ------------------------------------------------------------------
    # 5.  Vertex appearance — leaves get shaft hue + large circle +
    #     radial text label; cluster T-junctions get light-blue square;
    #     root gets a slightly larger light-blue square at the centre.
    # ------------------------------------------------------------------
    shaft_clr = shaft_colors(probes[:N])
    vcolor = g.new_vertex_property("vector<double>")
    vsize = g.new_vertex_property("double")
    vshape = g.new_vertex_property("string")
    vtext = g.new_vertex_property("string")
    vtext_pos = g.new_vertex_property("double")
    vtext_color = g.new_vertex_property("vector<double>")
    vfont_size = g.new_vertex_property("double")

    for v in g.vertices():
        vid = int(v)
        c = shaft_clr[vid]
        vcolor[v] = ([c[0], c[1], c[2], 1.0] if len(c) == 3 else list(c))
        vsize[v] = 24.0
        vshape[v] = "circle"
        # Labels are drawn by the matplotlib overlay (panel_c_hierarchy)
        # — graph-tool's vertex_text was unreadable at font 7.5 on a
        # 2200-px raster.
        vtext[v] = ""
        vtext_pos[v] = -1
        vtext_color[v] = [0.15, 0.15, 0.15, 1.0]
        vfont_size[v] = 0.0

    gt.graph_draw(
        g, pos=pos,
        output=str(out_path),
        output_size=(output_size, output_size),
        vertex_fill_color=vcolor,
        vertex_color=[0.10, 0.10, 0.10, 0.95],
        vertex_size=vsize,
        vertex_shape=vshape,
        vertex_text=vtext,
        vertex_text_position=vtext_pos,
        vertex_text_color=vtext_color,
        vertex_font_size=vfont_size,
        edge_color=ecolor,
        edge_pen_width=pen,
        eorder=draw_order,
        edge_gradient=[],
        edge_control_points=cts,
        bg_color=[1, 1, 1, 1],
        fit_view=fit_view,
    )

    # Return layout metadata so the caller can overlay per-leaf
    # contact labels in the same image-pixel coordinate space.
    return dict(
        N=N,
        theta_of_leaf=theta_of_leaf,
        R_outer=R_OUTER,
        output_size=output_size,
        fit_view=fit_view,
    )


def panel_c_hierarchy(ax: plt.Axes, exemplar: str, band: str, band_tex: str,
                      edges: dict, pair: dict, probes: list[str],
                      contact_labels: list[str],
                      imcoh_weights: np.ndarray,
                      lrg, tmp_dir: Path) -> None:
    """Render the graph-tool hierarchy-bundled chord to a temp PNG and
    embed it in ``ax`` via ``imshow``; the layout metadata returned by
    the renderer is reused to overlay the per-leaf contact labels in
    the same image-pixel coordinate space.
    """
    png_path = tmp_dir / f"panel_c_{exemplar}_{band}.png"
    meta = render_lrg_hierarchy_chord(png_path, lrg, pair, probes, edges,
                                       contact_labels, imcoh_weights)
    img = imread(png_path)
    ax.imshow(img, interpolation="lanczos")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    # ------------------------------------------------------------------
    # Coordinate transform — local (-1, 1) box → image-pixel space.
    # graph-tool's ``fit_view=fv`` scales the data bbox (here roughly
    # [-1, 1] × [-1, 1]) to occupy ``fv * H`` pixels, centered on the
    # canvas.  Image coords have origin top-left, y growing DOWN —
    # hence the y-flip below.
    # ------------------------------------------------------------------
    N = meta["N"]
    theta_of_leaf = meta["theta_of_leaf"]
    R_outer_render = meta["R_outer"]
    H = meta["output_size"]
    fv = meta["fit_view"]

    half_span = R_outer_render  # bbox extends from -R to +R in local space
    scale_px = fv * H / (2 * half_span)
    cx, cy = H / 2.0, H / 2.0

    def _to_px_xy(x, y):
        return cx + x * scale_px, cy - y * scale_px

    # ------------------------------------------------------------------
    # Per-leaf contact labels around the outer ring — rotated
    # tangentially so adjacent labels don't collide and text on the
    # left half flips for readability (standard chord-diagram idiom).
    # Epi-zone contacts (red font in the implant Excel) are coloured
    # red; matching uses the canonical ``_normalise_label`` which
    # preserves apostrophes so primed/unprimed shafts don't collide.
    # ------------------------------------------------------------------
    contact_labels_local = load_contact_labels(exemplar)
    epi_set = {_normalise_label(s) for s in load_epileptic_nodes(exemplar)}
    meta_lbl = load_patient_metadata(exemplar, SEEG_DATAPATH)
    if meta_lbl is not None and "label" in meta_lbl.columns:
        is_epi = [_normalise_label(s) in epi_set
                  for s in meta_lbl["label"].tolist()]
    else:
        is_epi = [False] * N
    offset_px = H * 0.015
    for i in range(N):
        if i >= len(contact_labels_local):
            break
        # Leaf in image-pixel coords (at r = 1.0 in local space).
        lx, ly = math.cos(theta_of_leaf[i]), math.sin(theta_of_leaf[i])
        lpx, lpy = _to_px_xy(lx, ly)
        # Outward unit vector in image coords (y-flip).
        ux, uy = (lpx - cx), (lpy - cy)
        norm = math.hypot(ux, uy)
        if norm < 1e-6:
            continue
        ux /= norm; uy /= norm
        tx, ty = lpx + ux * offset_px, lpy + uy * offset_px
        angle_deg = math.degrees(math.atan2(-(lpy - cy), lpx - cx))
        if -90 < angle_deg <= 90:
            rotation, ha = angle_deg, "left"
        else:
            rotation, ha = angle_deg + 180, "right"
        is_epi_i = i < len(is_epi) and is_epi[i]
        label_color = CLR_ANTI if is_epi_i else "0.18"
        ax.text(tx, ty, contact_labels_local[i],
                rotation=rotation, rotation_mode="anchor",
                ha=ha, va="center",
                fontsize=4.8, color=label_color,
                fontweight=("bold" if is_epi_i else "normal"),
                zorder=8, clip_on=False)
    # Title placed via figure-level text in main().


# ---------------------------------------------------------------------------
# Panel (b) — 3-phase circular dendrograms with continuous task-order
# rainbow leaf colouring (no partition cut).  Cross-phase similarity
# is the per-pair cophenet distance correlation
#   ρ^coph(A, B) = Spearman(D_coph_A[triu], D_coph_B[triu])
# matching the §5.3 ρ_split^coph metric family.  K-partition / ARI
# metrics are NOT used: at any single K the partition is dominated by
# one large cluster and the metric is noise.
# ---------------------------------------------------------------------------
PERSIST_PHASES = ("rest_pre", "task_test", "rest_post")
PERSIST_PHASE_TEX = {
    "rest_pre":  r"rest$_{\mathrm{pre}}$",
    "task_test": r"task$_{\mathrm{test}}$",
    "rest_post": r"rest$_{\mathrm{post}}$",
}


def _draw_circular_dendrogram_panel(
    ax: plt.Axes, Z: np.ndarray, *,
    leaf_colors=None, leaf_labels=None,
    branch_colors=None,
    h_max_override: float | None = None,
    node_r_override: np.ndarray | None = None,
    r_outer: float = 1.0, r_inner: float = 0.10,
    line_width: float = 0.85,
    node_size: float = 22.0,
    label_fontsize: float = 4.5, label_offset: float = 1.05,
    branch_desaturate: float = 0.0,
) -> np.ndarray:
    """Standalone circular dendrogram with subtree-coloured branches
    and leaf-coloured radial labels.

    Each branch's colour is the size-weighted RGB mean of its
    descendant leaves' colours.  Subtrees whose leaves share similar
    colours keep a vibrant branch hue; subtrees that mix far-apart
    colours wash out to gray — so persistence reads visually as
    branch vibrancy along the tree.

    Equiangular leaves in ``leaves_list(Z)`` order around the unit
    circle; internal nodes at ``r = R_outer · (1 − h / h_max)``; each
    merge is two radial spokes from children to merge radius plus one
    CCW arc joining the children at that radius.
    """
    import math
    from matplotlib.colors import to_rgb
    from scipy.cluster.hierarchy import leaves_list

    Z = np.asarray(Z)
    N = int(Z.shape[0]) + 1
    h_max = float(h_max_override) if h_max_override is not None else float(Z[-1, 2])
    n_total = 2 * N - 1

    leaf_order = [int(x) for x in leaves_list(Z)]
    n_leaves = len(leaf_order)
    pos_in_order = {lf: idx for idx, lf in enumerate(leaf_order)}

    theta_of_leaf = np.zeros(N)
    for slot, lf in enumerate(leaf_order):
        theta_of_leaf[lf] = 2 * math.pi * slot / n_leaves - math.pi / 2

    node_r = np.zeros(n_total)
    node_theta = np.zeros(n_total)
    idx_min = np.zeros(n_total, dtype=int)
    idx_max = np.zeros(n_total, dtype=int)

    for i in range(N):
        node_r[i] = r_outer
        node_theta[i] = theta_of_leaf[i]
        idx_min[i] = pos_in_order[i]
        idx_max[i] = pos_in_order[i]

    def _angle(slot_f: float) -> float:
        return 2 * math.pi * slot_f / n_leaves - math.pi / 2

    for r in range(Z.shape[0]):
        node_id = N + r
        h = float(Z[r, 2])
        c1, c2 = int(Z[r, 0]), int(Z[r, 1])
        if node_r_override is not None:
            node_r[node_id] = float(node_r_override[node_id])
        else:
            node_r[node_id] = r_outer - (h / max(h_max, 1e-12)) * (r_outer - r_inner)
        idx_min[node_id] = min(idx_min[c1], idx_min[c2])
        idx_max[node_id] = max(idx_max[c1], idx_max[c2])
        node_theta[node_id] = _angle(0.5 * (idx_min[node_id] + idx_max[node_id]))

    # Per-node subtree colour (size-weighted RGB mean of descendants).
    # This is what makes the persistence visible — subtrees of similar-
    # task-slot leaves keep a vibrant hue, mixed-slot subtrees wash to gray.
    if leaf_colors is None:
        leaf_rgb_arr = np.full((N, 3), 0.35)
    else:
        leaf_rgb_arr = np.array([to_rgb(c) for c in leaf_colors])
    node_rgb = np.zeros((n_total, 3))
    node_size_w = np.zeros(n_total)
    for i in range(N):
        node_rgb[i] = leaf_rgb_arr[i]
        node_size_w[i] = 1.0
    for r in range(Z.shape[0]):
        node_id = N + r
        c1, c2 = int(Z[r, 0]), int(Z[r, 1])
        s1, s2 = float(node_size_w[c1]), float(node_size_w[c2])
        node_rgb[node_id] = (s1 * node_rgb[c1] + s2 * node_rgb[c2]) / (s1 + s2)
        node_size_w[node_id] = s1 + s2

    def _mix_with_white(rgb, frac):
        if frac <= 0:
            return rgb
        return rgb + (np.array([1.0, 1.0, 1.0]) - rgb) * frac

    # Branches — colour source:
    #   • if ``branch_colors`` is given (one RGB(A) per merge in Z),
    #     every part of merge ``r``'s branch (both radial spokes + arc)
    #     uses ``branch_colors[r]``.  Used for the "task-cophenet at
    #     this merge" encoding so radial colour vs radial position
    #     reveals cross-phase agreement directly.
    #   • else the size-weighted subtree mean of leaf colours (legacy).
    have_branch_colors = branch_colors is not None
    for r in range(Z.shape[0]):
        node_id = N + r
        c1, c2 = int(Z[r, 0]), int(Z[r, 1])
        if idx_min[c1] <= idx_min[c2]:
            c_left, c_right = c1, c2
        else:
            c_left, c_right = c2, c1
        merge_r = node_r[node_id]
        t_left = node_theta[c_left]
        t_right = node_theta[c_right]
        if have_branch_colors:
            seg_rgb_left = seg_rgb_right = arc_rgb = tuple(
                branch_colors[r][:3]
            )
        for c, side in ((c_left, "left"), (c_right, "right")):
            cr = node_r[c]; th = node_theta[c]
            if have_branch_colors:
                rgb_seg = (seg_rgb_left if side == "left" else seg_rgb_right)
            else:
                rgb_seg = tuple(
                    _mix_with_white(node_rgb[c], branch_desaturate)
                )
            ax.plot([cr * math.cos(th), merge_r * math.cos(th)],
                    [cr * math.sin(th), merge_r * math.sin(th)],
                    color=rgb_seg, lw=line_width, zorder=3,
                    solid_capstyle="round")
        dth = (t_right - t_left) % (2 * math.pi)
        arc_th = t_left + np.linspace(0.0, dth, 48)
        if have_branch_colors:
            rgb_arc = arc_rgb
        else:
            rgb_arc = tuple(
                _mix_with_white(node_rgb[node_id], branch_desaturate)
            )
        ax.plot(merge_r * np.cos(arc_th), merge_r * np.sin(arc_th),
                color=rgb_arc, lw=line_width, zorder=3,
                solid_capstyle="round")

    # Leaf markers
    leaf_x = r_outer * np.cos(theta_of_leaf)
    leaf_y = r_outer * np.sin(theta_of_leaf)
    if leaf_colors is None:
        leaf_colors = ["0.35"] * N
    ax.scatter(leaf_x, leaf_y, c=leaf_colors,
               s=node_size, edgecolors="white", linewidths=0.55,
               zorder=5)

    # Radial labels — coloured by the leaf hue
    if leaf_labels is not None:
        for i in range(min(N, len(leaf_labels))):
            if not leaf_labels[i]:
                continue
            th = float(theta_of_leaf[i])
            lr = r_outer * label_offset
            x = lr * math.cos(th); y = lr * math.sin(th)
            ang = math.degrees(th)
            if -90 < ang <= 90:
                rot, ha = ang, "left"
            else:
                rot, ha = ang + 180, "right"
            ax.text(x, y, leaf_labels[i],
                    rotation=rot, rotation_mode="anchor",
                    ha=ha, va="center",
                    fontsize=label_fontsize,
                    color=leaf_colors[i],
                    zorder=6, clip_on=False)

    ax.set_xlim(-r_outer * 1.13, r_outer * 1.13)
    ax.set_ylim(-r_outer * 1.13, r_outer * 1.13)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    return theta_of_leaf


def panel_b_phase_persistence(ax_pre: plt.Axes, ax_task: plt.Axes,
                               ax_post: plt.Axes,
                               exemplar: str, band: str) -> dict:
    """Render 3 circular dendrograms (rest_pre, task_test, rest_post)
    for one exemplar patient, with leaves coloured continuously by
    their slot index in task_test's ``leaves_list`` order (a 1D rainbow
    on contact identity, fixed across phases).

    Persistence reads visually as the angular SMOOTHNESS of the rainbow
    around the rest_post ring (adjacent leaves carry adjacent task
    slots → task hierarchy preserved).  rest_pre tends to scramble the
    rainbow more than rest_post does.

    Annotations are the cross-phase per-pair cophenet correlation
    ``ρ^coph(A, B) = Spearman(D_coph_A[triu], D_coph_B[triu])`` — the
    project canonical metric (continuous, not partition-based).
    """
    from scipy.cluster.hierarchy import (
        leaves_list, cophenet, optimal_leaf_ordering,
    )
    from scipy.spatial.distance import squareform

    lrgs = {}
    for ph in PERSIST_PHASES:
        lrg = load_lrg_result(exemplar, ph, band, "imcoh_abs")
        if lrg is None or lrg.linkage_matrix is None:
            raise FileNotFoundError(
                f"No LRG cache for {exemplar}/{band}/{ph}/imcoh_abs"
            )
        lrgs[ph] = lrg
    Zs = {ph: np.asarray(lrgs[ph].linkage_matrix, dtype=float)
          for ph in PERSIST_PHASES}

    # (1) Optimal-leaf-order task_test against its OWN cophenet — gives
    # the colour-defining ring a smooth task-hierarchy progression
    # without altering tree topology.
    Zs["task_test"] = optimal_leaf_ordering(
        Zs["task_test"], cophenet(Zs["task_test"]),
    )
    task_leaf_order = [int(x) for x in leaves_list(Zs["task_test"])]
    N = len(task_leaf_order)
    task_slot_of = np.zeros(N, dtype=int)
    for slot, c in enumerate(task_leaf_order):
        if c < N:
            task_slot_of[c] = slot

    # (2) Optimal-leaf-order rest_pre / rest_post AGAINST the task slot
    # distance ``d(i, j) = |slot_task(i) − slot_task(j)|`` so that, at
    # every merge, scipy picks the child swap that puts similar-task-slot
    # leaves angularly adjacent in this phase's ring.  Tree topology is
    # untouched; only child orderings at each internal node move.  This
    # is the canonical "remove crossings" reordering — same idea that
    # ``seriation`` / Bar-Joseph 2003 use.
    task_dist_2d = np.abs(
        task_slot_of[:, None] - task_slot_of[None, :]
    ).astype(float)
    np.fill_diagonal(task_dist_2d, 0.0)
    y_task_slot = squareform(task_dist_2d, checks=False)
    for ph in ("rest_pre", "rest_post"):
        Zs[ph] = optimal_leaf_ordering(Zs[ph], y_task_slot)

    # (3) Continuous rainbow cmap — custom fixed-luminance HLS sweep so
    # every hue lands at the SAME dark-mid luminance (L=0.42, S=0.90).
    # In particular the yellow band (H≈0.17) becomes dark olive
    # ≈ #c1c10a instead of the near-white #ffff80 that bleached
    # leaves out against the paper background.  No matplotlib stock
    # cmap (turbo, viridis-family, Spectral) avoids the bright-yellow
    # region cleanly — fixed-luminance generation is the cleanest fix.
    import colorsys as _cs
    from matplotlib.colors import LinearSegmentedColormap as _LSC
    _anchors_no_yellow = [
        _cs.hls_to_rgb(h, 0.42, 0.90)
        for h in np.linspace(0.0, 0.83, 24)
    ]
    cmap_seq = _LSC.from_list("rainbow_fixed_L", _anchors_no_yellow, N=256)
    leaf_colors_task = [cmap_seq(task_slot_of[c] / max(N - 1, 1))
                        for c in range(N)]
    GRAY_RGBA = (0.80, 0.80, 0.80, 1.0)
    PERSIST_THRESHOLD = 0.35  # per-leaf ρ_c gate (target ~60% kept in post)

    # (4) Per-leaf cophenet correlation gate.  For each leaf ``c`` and
    # each rest phase X, ρ_c(X) = Spearman over ``j ≠ c`` of
    # ``D_coph_task[c, j]`` vs ``D_coph_X[c, j]``.  Leaves whose
    # neighbourhood in task is preserved in X (ρ_c > 0.5) keep their
    # turbo task-slot colour; the rest are greyed.  rest_post should
    # therefore carry MORE coloured leaves than rest_pre — that is
    # the persistence story made local.
    Dc = {ph: squareform(cophenet(Zs[ph])) for ph in PERSIST_PHASES}
    leaf_colors_per_phase: dict = {"task_test": leaf_colors_task}
    n_colored_per_phase: dict = {"task_test": N}
    for ph in ("rest_pre", "rest_post"):
        colors_ph = []
        n_kept = 0
        for c in range(N):
            mask = np.ones(N, dtype=bool)
            mask[c] = False
            rho_c, _ = spearmanr(Dc["task_test"][c, mask],
                                 Dc[ph][c, mask])
            if np.isfinite(rho_c) and rho_c > PERSIST_THRESHOLD:
                colors_ph.append(leaf_colors_task[c])
                n_kept += 1
            else:
                colors_ph.append(GRAY_RGBA)
        leaf_colors_per_phase[ph] = colors_ph
        n_colored_per_phase[ph] = n_kept

    contact_labels = load_contact_labels(exemplar)[:N]

    # Unified h_max across the three phases so radial scales are
    # directly comparable (a merge at the same cophenet height lands
    # at the same radius in any of the three panels).
    h_max_unified = float(max(Zs[ph][-1, 2] for ph in PERSIST_PHASES))

    ax_map = {"rest_pre": ax_pre, "task_test": ax_task,
              "rest_post": ax_post}
    for ph in PERSIST_PHASES:
        _draw_circular_dendrogram_panel(
            ax_map[ph], Zs[ph],
            leaf_colors=leaf_colors_per_phase[ph],
            leaf_labels=contact_labels,
            h_max_override=h_max_unified,
            r_outer=1.0, r_inner=0.08,
            line_width=0.85,
            node_size=14.0,
            label_fontsize=3.6,
            label_offset=1.04,
        )
        ax_map[ph].set_title(PERSIST_PHASE_TEX[ph],
                             fontsize=10.5, pad=2)

    # Per-pair cophenet distance correlation between phases
    iu = np.triu_indices(N, k=1)
    rho_pre, _ = spearmanr(Dc["task_test"][iu], Dc["rest_pre"][iu])
    rho_post, _ = spearmanr(Dc["task_test"][iu], Dc["rest_post"][iu])

    n_pre = n_colored_per_phase["rest_pre"]
    n_post = n_colored_per_phase["rest_post"]

    # Headline numbers go to stdout / caption ledger only — the figure
    # itself shows the contrast visually (more coloured leaves in post
    # than in pre).  No in-axes annotation: it stole vertical space
    # from the dendrograms.
    print(f"  panel_b ρ^coph: (task, pre) = {rho_pre:+.4f}, "
          f"(task, post) = {rho_post:+.4f}")
    print(f"  per-leaf ρ_c > {PERSIST_THRESHOLD}: "
          f"pre = {n_pre}/{N}, post = {n_post}/{N}")
    return dict(rho_pre=rho_pre, rho_post=rho_post,
                n_pre_kept=n_pre, n_post_kept=n_post, N=N)


# ---------------------------------------------------------------------------
# Panel (c, archived) — coarse joint-rank deviation heatmap
# ---------------------------------------------------------------------------
N_HEATMAP_BINS = 20           # bin count per side for the joint-rank heatmap
HEATMAP_VMAX = 0.60           # symmetric cmap clip on (P_obs/P_null) − 1


def panel_c_deviation_heatmap(ax: plt.Axes, band: str,
                              n_bins: int = N_HEATMAP_BINS,
                              vmax: float = HEATMAP_VMAX) -> None:
    r"""Cohort joint-rank-density deviation from the uniform null.

    Pool within-patient rank pairs ``(u, v) = (rank Δ_task, rank Δ_rest)``
    across the ten cohort patients, bin onto an ``n_bins × n_bins``
    grid on ``[0, 1]^2``, and colour each cell by the relative
    deviation from the uniform null:

    .. math::

        z(u_{\rm bin}, v_{\rm bin}) =
            \frac{P_{\rm obs}(u_{\rm bin}, v_{\rm bin})}{P_{\rm null}}
            - 1, \qquad P_{\rm null} = 1\ \text{on}\ [0,1]^2.

    Why this version and not the raw 60×60 cividis density:  at
    60 bins the per-cell SE swamps a +0.22 Spearman signal (cell-SE
    ~13%, diagonal-corner excess ~40% only at the very tail).  At
    10 bins each cell pools ~640 pairs (SE ~4%), so the diagonal-
    corner excess (~30–50%) and anti-corner depletion (~15–25%)
    that the Gaussian-copula model predicts for ρ ≈ 0.22 rise to
    ~10σ per cell.  The diagonal stripe is visible at a glance.
    """
    # Pool within-patient ranks across the cohort
    us, vs = [], []
    for pat in COHORT:
        pair = load_pair_data(pat, band)
        n = len(pair["dD_task"])
        if n < 2:
            continue
        u = (rankdata(pair["dD_task"], method="ordinal") - 0.5) / n
        v = (rankdata(pair["dD_rest"], method="ordinal") - 0.5) / n
        us.append(u); vs.append(v)
    u_pool = np.concatenate(us)
    v_pool = np.concatenate(vs)
    n_total = u_pool.size

    edges = np.linspace(0.0, 1.0, n_bins + 1)
    H, _, _ = np.histogram2d(u_pool, v_pool, bins=[edges, edges])
    expected = n_total / (n_bins * n_bins)
    z = H / expected - 1.0  # relative deviation from uniform null

    # Diverging cmap anchored at the trace/anti palette (grass-green +,
    # brick-red −, near-white at 0).
    from matplotlib.colors import LinearSegmentedColormap
    cmap_dev = LinearSegmentedColormap.from_list(
        "trace_dev",
        [CLR_ANTI, "#f7f3eb", CLR_TRACE],
        N=256,
    )

    im = ax.imshow(z.T, origin="lower", extent=[0.0, 1.0, 0.0, 1.0],
                   cmap=cmap_dev, vmin=-vmax, vmax=+vmax,
                   aspect="equal", interpolation="nearest", zorder=1)

    # Identity diagonal (concordance) and anti-diagonal (anti) refs
    ax.plot([0, 1], [0, 1], color="0.18", lw=1.4, ls="--",
            alpha=0.55, zorder=3)
    ax.plot([0, 1], [1, 0], color="0.18", lw=0.9, ls=":",
            alpha=0.35, zorder=3)

    # Cell-value annotations (% deviation) — threshold scales with bin
    # count so 20×20 only labels the strongest cells (otherwise the grid
    # is overrun with tiny digits).
    annot_threshold = 0.10 if n_bins <= 12 else 0.30
    for i in range(n_bins):
        for j in range(n_bins):
            zi = float(z[i, j])
            if abs(zi) < annot_threshold:
                continue
            x = (i + 0.5) / n_bins
            y = (j + 0.5) / n_bins
            text_clr = "white" if abs(zi) > vmax * 0.55 else "0.08"
            ax.text(x, y, rf"{zi * 100:+.0f}",
                    ha="center", va="center",
                    fontsize=6.5 if n_bins > 12 else 7.5,
                    color=text_clr,
                    fontweight="bold", zorder=4)

    ax.set_xlim(0.0, 1.0); ax.set_ylim(0.0, 1.0)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_xticklabels(["0", r"$\frac{1}{2}$", "1"])
    ax.set_yticklabels(["0", r"$\frac{1}{2}$", "1"])
    ax.set_xlabel(r"within-patient rank of "
                  r"$\Delta_{\mathrm{task}}(i,j)$  ($u$)",
                  fontsize=9.5, labelpad=4)
    ax.set_ylabel(r"within-patient rank of "
                  r"$\Delta_{\mathrm{rest}}(i,j)$  ($v$)",
                  fontsize=9.5, labelpad=4)
    ax.tick_params(labelsize=9.5)

    # Colorbar — tight inset on the right
    cb = ax.figure.colorbar(
        im, ax=ax, fraction=0.046, pad=0.022, shrink=0.88,
    )
    cb.set_label(r"$(P_{\rm obs} / P_{\rm null}) - 1$",
                 fontsize=8.5, labelpad=4)
    cb.ax.tick_params(labelsize=8.0)
    cb.set_ticks([-vmax, -vmax/2, 0.0, +vmax/2, +vmax])
    cb.set_ticklabels([rf"$-{vmax:.2f}$",
                       rf"$-{vmax/2:.2f}$", "$0$",
                       rf"$+{vmax/2:.2f}$",
                       rf"$+{vmax:.2f}$"])

    # Pooled Spearman is printed to stdout for the caller; no in-axes
    # annotation (kept off the figure by user request).
    pooled_rho, _ = spearmanr(u_pool, v_pool)
    print(f"  panel_c heatmap: pooled ρ = {pooled_rho:+.4f} "
          f"({n_total:,} pairs, 10 patients)")

    ax.set_title("(c)", loc="left", fontweight="bold",
                 fontsize=12.0, pad=4)


# ---------------------------------------------------------------------------
# Panel (c, archived) — cohort conditional sign-concordance curve
# ---------------------------------------------------------------------------
N_QUANTILE_BINS = 10


def _per_patient_concordance_curve(pat: str, band: str,
                                    n_bins: int = N_QUANTILE_BINS) -> np.ndarray:
    """Conditional sign-concordance of ``(Δ_task, Δ_rest)`` as a function of
    the within-patient quantile of ``|Δ_task|``.

    For each bin ``b ∈ [0, 1]`` of the per-pair quantile rank of
    ``|Δ_task|``, returns
    ``f(b) = P(sign(Δ_rest) = sign(Δ_task) | |Δ_task| in bin b)``.
    Under the matched-strength null this is 0.5 in every bin; under a
    cross-phase memory trace it rises with the quantile because the
    pairs that move the most under task tend to move in the same
    direction during rest.
    """
    pair = load_pair_data(pat, band)
    dt = pair["dD_task"]
    dr = pair["dD_rest"]
    n = len(dt)
    if n < n_bins * 10:
        return np.full(n_bins, np.nan)
    # within-patient quantile of |Δ_task|; ordinal ranks → exact splits
    q = (rankdata(np.abs(dt), method="ordinal") - 0.5) / n
    concord = (np.sign(dt) == np.sign(dr)).astype(float)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    curve = np.full(n_bins, np.nan)
    for k in range(n_bins):
        mask = (q >= edges[k]) & (q < edges[k + 1])
        if mask.sum() >= 5:
            curve[k] = float(concord[mask].mean())
    return curve


def panel_c_concordance(ax: plt.Axes, band: str,
                         n_bins: int = N_QUANTILE_BINS) -> None:
    """Cohort conditional sign-concordance curve at this band.

    Per-patient curves drawn faint grey; cohort mean overlaid in
    grass-green to match the panel (b) trace colour.  Null line at
    ``y = 0.5`` is the matched-strength expectation under a
    strength-only model (no cross-phase signal).
    """
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    curves = np.full((len(COHORT), n_bins), np.nan, dtype=float)
    for pi, pat in enumerate(COHORT):
        curves[pi] = _per_patient_concordance_curve(pat, band, n_bins)
    cohort_mean = np.nanmean(curves, axis=0)

    # Null shading (grey band) + axis line at 0.5
    ax.axhspan(0.495, 0.505, facecolor="0.85", edgecolor="none", zorder=0)
    ax.axhline(0.5, color="0.45", lw=1.0, ls="--", zorder=1)

    # Per-patient curves (thin grey)
    for pi in range(len(COHORT)):
        if np.all(np.isnan(curves[pi])):
            continue
        ax.plot(centers, curves[pi], color="0.58", lw=0.95, alpha=0.55,
                marker="o", ms=2.4, zorder=2)

    # Cohort mean — bold grass-green
    ax.plot(centers, cohort_mean, color=CLR_TRACE, lw=2.8,
            marker="o", ms=7.0, markerfacecolor="white",
            markeredgecolor=CLR_TRACE, markeredgewidth=2.0,
            zorder=5, label="cohort mean")

    # Reference: identity-tilt line for ρ=1 (every pair concordant)
    # only as a faint hint at top-right; skip if it overlaps too much.

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0.40, 0.80)
    ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "¼", "½", "¾", "1"])
    ax.set_yticks([0.4, 0.5, 0.6, 0.7, 0.8])
    ax.tick_params(labelsize=9.5)
    ax.set_xlabel(r"within-patient quantile of $|\Delta_{\mathrm{task}}(i,j)|$",
                  fontsize=10, labelpad=4)
    ax.set_ylabel(r"$P(\,\mathrm{sign}\,\Delta_{\mathrm{rest}}"
                  r" = \mathrm{sign}\,\Delta_{\mathrm{task}}\,)$",
                  fontsize=10, labelpad=4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.015, 0.502, "null = 0.5", transform=ax.get_yaxis_transform(),
            ha="left", va="bottom", fontsize=8.5, color="0.40",
            style="italic", zorder=6)

    # Annotation: bottom-decile vs top-decile concordance, cohort mean.
    bot_c = float(cohort_mean[0])
    top_c = float(cohort_mean[-1])
    delta = top_c - bot_c
    ax.annotate(
        rf"$f_{{\rm top}} - f_{{\rm bot}} = {delta:+.3f}$",
        xy=(centers[-1], top_c),
        xytext=(0.20, 0.745),
        fontsize=9.5, fontweight="bold", color=CLR_TRACE,
        arrowprops=dict(arrowstyle="->", color=CLR_TRACE, lw=1.1,
                        connectionstyle="arc3,rad=-0.18"),
        zorder=7,
    )
    ax.text(0.97, 0.06,
            rf"$f_{{\rm top}} = {top_c:.3f}$" "\n"
            rf"$f_{{\rm bot}} = {bot_c:.3f}$" "\n"
            rf"$n_{{\rm pat}} = 10$",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8.5, color="0.15",
            bbox=dict(boxstyle="round,pad=0.28",
                      facecolor="white", edgecolor="#aaa",
                      alpha=0.92), zorder=8)

    ax.set_title("(c)", loc="left", fontweight="bold",
                 fontsize=12.0, pad=4)


# ---------------------------------------------------------------------------
# Figure assembly
# ---------------------------------------------------------------------------
def main(band: str = "beta", exemplar: str = "Pat_05",
         alpha: float = ALPHA_DEFAULT,
         imcoh_phase: str = IMCOH_PHASE_DEFAULT) -> Path:
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")

    per_pat = pd.read_csv(LRG_CTM_DIR / "per_patient_per_band.csv")
    cohort = pd.read_csv(LRG_CTM_DIR / "cohort_summary.csv")
    cohort_row = cohort[cohort.band == band].iloc[0]

    _, probes = _load_inputs(exemplar, band, "rest_pre", "imcoh_abs")
    contact_labels = load_contact_labels(exemplar)

    pair = load_pair_data(exemplar, band)

    # Per-pair matched-strength p-values (cached at PER_PAIR_CACHE)
    pvals = compute_or_load_per_pair_pvals(exemplar, band)
    # FDR-gated trace/anti classification over the cross-probe family
    edges = select_significant_edges(
        pair, pvals["p_pair"], pvals["c_obs"], probes=probes,
        alpha=alpha,
        cross_probe_only=CROSS_PROBE_ONLY_DEFAULT,
    )

    # rest_post (or chosen phase) ImCoh|·| per pair → edge magnitude
    fc_imcoh = load_fc_matrix(exemplar, imcoh_phase, band,
                              fc_method="imcoh_abs")
    iu_i = pair["iu_i"]; iu_j = pair["iu_j"]
    N_fc = fc_imcoh.shape[0]
    imcoh_weights = np.zeros(len(iu_i), dtype=float)
    in_range = (iu_i < N_fc) & (iu_j < N_fc)
    imcoh_weights[in_range] = fc_imcoh[iu_i[in_range], iu_j[in_range]]

    print(f"  matched-strength gated edges (p ≤ {alpha}): "
          f"{len(edges['trace_idx'])} trace + "
          f"{len(edges['anti_idx'])} anti "
          f"out of {edges['n_tested']} cross-probe pairs tested "
          f"(valid R={pvals['valid_R']})")
    lrg = load_lrg_result(exemplar, "rest_pre", band, "imcoh_abs")
    if lrg is None or lrg.linkage_matrix is None:
        raise FileNotFoundError(
            f"No LRG result cached for {exemplar} / {band} / rest_pre / imcoh_abs"
        )

    out_dir = ROOT / "data" / "preprint" / "figures" / band / "per_pair_trace"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Single-row 2-panel layout — (a) per-patient z-stems on the left,
    # (b) the three-phase circular-dendrogram persistence triptych
    # (rest_pre, task_test, rest_post) on the right.  Figure dimensions
    # chosen so each sub-dendrogram (square aspect = "equal") has its
    # side length equal to panel (a)'s vertical axis height — both are
    # visually 4.5 inches tall.
    fig = plt.figure(figsize=(21.0, 5.5))
    gs = fig.add_gridspec(
        1, 2, width_ratios=[1.0, 2.55],
        wspace=0.10,
        left=0.030, right=0.99, top=0.93, bottom=0.10,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    gs_b = gs[0, 1].subgridspec(1, 3, wspace=0.04)
    ax_b_pre = fig.add_subplot(gs_b[0, 0])
    ax_b_task = fig.add_subplot(gs_b[0, 1])
    ax_b_post = fig.add_subplot(gs_b[0, 2])

    panel_a_zstrip(ax_a, per_pat, cohort_row, band, band_tex)
    persist_meta = panel_b_phase_persistence(
        ax_b_pre, ax_b_task, ax_b_post,
        exemplar=exemplar, band=band,
    )

    fig.text(
        ax_b_pre.get_position().x0, ax_b_pre.get_position().y1 + 0.020,
        "(b)", ha="left", va="bottom", fontsize=12.0, fontweight="bold",
    )

    out = out_dir / f"fig_{band}_rho_split_test2.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--band", default="beta",
                        choices=["delta", "theta", "alpha", "beta",
                                 "low_gamma", "high_gamma"])
    parser.add_argument("--exemplar", default="Pat_05",
                        choices=COHORT)
    parser.add_argument("--alpha", type=float, default=ALPHA_DEFAULT,
                        help="per-pair matched-strength p threshold (default 0.05)")
    parser.add_argument("--imcoh-phase", default=IMCOH_PHASE_DEFAULT,
                        choices=["rest_pre", "task_test", "rest_post"],
                        help="Phase whose ImCoh|·| weighting sets edge magnitude")
    args = parser.parse_args()
    main(args.band, args.exemplar,
         alpha=args.alpha, imcoh_phase=args.imcoh_phase)
