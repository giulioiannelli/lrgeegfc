#!/usr/bin/env python3
"""Audit 25 — raw-FC phase-distance (step 0, pre-LRG).

Implements the scope at
``.agents/guides/task-persistence-investigation/2026-04-28_raw-fc-phase-distance.md``.

Single substrate: raw ``A^φ`` from ``imcoh_abs``, restricted to the
giant-component intersection ``V*(p, b) = ⋂_φ GC(|A^φ|)``. Three
distances on the substrate ladder:

    d_P(A, B) = 1 − corr_P(triu(A, k=1), triu(B, k=1))
    d_S(A, B) = 1 − corr_S(triu(A, k=1), triu(B, k=1))
    d_F(A, B) = ‖A − B‖_F / √(‖A‖_F · ‖B‖_F)

Within-``rest_pre`` split-half null (``n_split = 50``) yields a robust
z-score per phase pair. Triangle scalar
``T_d = d(task, post) − d(pre, task)`` operationalises persistence at
the FC-edge level.

Usage
-----
    # Single patient (default mode; Pat_06 unless overridden)
    python audit_25_raw_fc_phase_distance.py
    python audit_25_raw_fc_phase_distance.py --patient Pat_02
    # Cohort sweep (all 10 PATIENTS_4PHASE)
    python audit_25_raw_fc_phase_distance.py --cohort

Outputs
-------
``data/audit/raw_fc_phase_distance/<Patient>/audit.csv``
``data/audit/raw_fc_phase_distance/<Patient>/null.npz``
``data/audit/raw_fc_phase_distance/<Patient>/diagnostic.pdf``
``data/audit/raw_fc_phase_distance/cohort_summary.csv`` (cohort mode)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import pearsonr, spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    DEFAULT_SAMPLE_RATE,
    FS_OVERRIDES,
    PATIENTS_4PHASE,
    nperseg_for_fs,
)
from lrg_eegfc.config.paths import DATA_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.visuals.correlation import imshow_colorbar_caxdivider
from lrg_eegfc.workflow.fc import load_fc_matrix
from scipy.signal import get_window

PHASES: Tuple[str, ...] = ("rest_pre", "task_test", "rest_post")
PHASE_PAIRS: Tuple[Tuple[str, str], ...] = (
    ("rest_pre", "task_test"),
    ("rest_pre", "rest_post"),
    ("task_test", "rest_post"),
)
DISTANCES: Tuple[str, ...] = ("P", "S", "F")
OUT_BASE = DATA_ROOT / "audit" / "raw_fc_phase_distance"


# ---------------------------------------------------------------------------
# Distance functions
# ---------------------------------------------------------------------------


def _triu_offdiag(A: np.ndarray) -> np.ndarray:
    iu = np.triu_indices_from(A, k=1)
    return A[iu]


def d_P(A: np.ndarray, B: np.ndarray) -> float:
    a, b = _triu_offdiag(A), _triu_offdiag(B)
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    r, _ = pearsonr(a, b)
    return float(1.0 - r)


def d_S(A: np.ndarray, B: np.ndarray) -> float:
    a, b = _triu_offdiag(A), _triu_offdiag(B)
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    r = spearmanr(a, b).statistic
    return float(1.0 - r)


def d_F(A: np.ndarray, B: np.ndarray) -> float:
    nA = np.linalg.norm(A, ord="fro")
    nB = np.linalg.norm(B, ord="fro")
    if nA == 0 or nB == 0:
        return float("nan")
    return float(np.linalg.norm(A - B, ord="fro") / np.sqrt(nA * nB))


def all_distances(A: np.ndarray, B: np.ndarray) -> Dict[str, float]:
    return {"P": d_P(A, B), "S": d_S(A, B), "F": d_F(A, B)}


# ---------------------------------------------------------------------------
# Substrate construction
# ---------------------------------------------------------------------------


def _common_giant_indices(adj_phases: Dict[str, np.ndarray]) -> np.ndarray:
    """Intersection of giant-component node sets across phases."""
    sets: List[set] = []
    for adj in adj_phases.values():
        G = nx.from_numpy_array(np.abs(adj))
        comps = list(nx.connected_components(G))
        if not comps:
            sets.append(set())
            continue
        sets.append(set(max(comps, key=len)))
    if not sets:
        return np.array([], dtype=int)
    common = set.intersection(*sets)
    return np.array(sorted(common), dtype=int)


def _segment_ffts(
    ts: np.ndarray, fs: float, nperseg: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Hann-windowed segment FFTs with **non-overlapping** segments.

    Returns ``(freqs, fft_seg)`` where ``fft_seg`` has shape
    ``(N, n_seg, F)``. Non-overlapping segments guarantee clean
    independence between any two segment subsets (split-half null), at the
    cost of ~half the segments compared to Welch's standard 50%-overlap
    convention. The cached observed FC uses the standard convention; this
    introduces a small bias (null distances slightly noisier than observed)
    that is conservative for the audit verdict.
    """
    N, L = ts.shape
    n_seg = L // nperseg
    if n_seg < 4:
        return np.array([]), np.empty((N, 0, 0))
    ts_trim = ts[:, : n_seg * nperseg].astype(np.float64, copy=False)
    segs = ts_trim.reshape(N, n_seg, nperseg)
    window = get_window("hann", nperseg).astype(np.float64)
    block = segs * window[None, None, :]
    fft_seg = np.fft.rfft(block, n=nperseg, axis=2)  # (N, n_seg, F)
    freqs = np.fft.rfftfreq(nperseg, 1.0 / fs)
    return freqs, fft_seg


def _imcoh_abs_from_subset(
    fft_seg: np.ndarray, subset: np.ndarray, freq_mask: np.ndarray
) -> np.ndarray:
    """Compute band-averaged ``|ImCoh|`` from a subset of segment FFTs."""
    sub = fft_seg[:, subset, :][:, :, freq_mask]  # (N, k, F_band)
    sub_conj = np.conj(sub)
    cross = np.einsum("nsf,msf->nmf", sub, sub_conj)  # (N, N, F_band)
    auto = np.real(np.einsum("nsf,nsf->nf", sub, sub_conj))  # (N, F_band)
    denom = np.sqrt(auto[:, None, :] * auto[None, :, :])
    with np.errstate(invalid="ignore", divide="ignore"):
        signed = np.where(denom > 0, np.imag(cross) / denom, 0.0)
    A = np.abs(signed).mean(axis=-1)
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    return A


def _split_half_null_from_ffts(
    fft_seg: np.ndarray,
    freqs: np.ndarray,
    fmin: float,
    fmax: float,
    common_idx: np.ndarray,
    n_split: int,
    rng: np.random.Generator,
) -> Dict[str, np.ndarray]:
    """Split-half null distances using precomputed segment FFTs.

    Random partition of segments into two equal halves; per split, both
    halves go through ``_imcoh_abs_from_subset`` and the three distances
    are computed on the ``V*``-restricted matrices.
    """
    if fft_seg.size == 0:
        return {d: np.full(n_split, np.nan) for d in DISTANCES}
    n_seg = fft_seg.shape[1]
    half = n_seg // 2
    if half < 4:
        return {d: np.full(n_split, np.nan) for d in DISTANCES}
    freq_mask = (freqs >= fmin) & (freqs < fmax)
    if not freq_mask.any():
        idx = int(np.clip(np.searchsorted(freqs, fmin), 0, len(freqs) - 1))
        freq_mask = np.zeros_like(freqs, dtype=bool)
        freq_mask[idx] = True

    null = {d: np.empty(n_split) for d in DISTANCES}
    for i in range(n_split):
        order = rng.permutation(n_seg)
        a_idx = order[:half]
        b_idx = order[half : 2 * half]
        A_a = _imcoh_abs_from_subset(fft_seg, a_idx, freq_mask)
        A_b = _imcoh_abs_from_subset(fft_seg, b_idx, freq_mask)
        A_a = A_a[np.ix_(common_idx, common_idx)]
        A_b = A_b[np.ix_(common_idx, common_idx)]
        for d_label, val in all_distances(A_a, A_b).items():
            null[d_label][i] = val
    return null


# ---------------------------------------------------------------------------
# Per-patient pipeline
# ---------------------------------------------------------------------------


def run_patient(
    patient: str,
    n_split: int,
    seed: int,
    out_dir: Path,
    *,
    reuse_null: bool = True,
) -> pd.DataFrame:
    """Per-patient audit pipeline.

    The expensive step is the within-`rest_pre` split-half null
    (50 imcoh-from-FFT-subset computations × 6 bands). We cache it on
    disk as ``<out_dir>/null.npz``; subsequent runs reuse it band-by-band
    and skip the ``rest_pre`` time-series load + segment-FFT precompute
    entirely if every (band, distance) entry is present. Pass
    ``reuse_null=False`` to force recompute.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    fs = FS_OVERRIDES.get(patient, DEFAULT_SAMPLE_RATE)
    nperseg = nperseg_for_fs(fs)

    print(f"[{patient}] fs={fs:.1f} Hz, nperseg={nperseg}")

    # Try to load cached null
    cached_null: Dict[str, np.ndarray] = {}
    null_path = out_dir / "null.npz"
    if reuse_null and null_path.exists():
        try:
            with np.load(null_path) as data:
                cached_null = {k: np.asarray(data[k]) for k in data.files}
            print(f"[{patient}] reusing cached null ({len(cached_null)} entries) "
                  f"from {null_path.name}")
        except Exception as exc:  # pragma: no cover
            print(f"[{patient}] failed to reload {null_path.name}: {exc}; recomputing")
            cached_null = {}

    def _have_cached(band: str) -> bool:
        return all(f"{band}_{d}" in cached_null for d in DISTANCES)

    need_fft = any(not _have_cached(b) for b in BRAIN_BANDS_NAMES)

    # Conditionally load time-series + precompute segment FFTs
    fft_seg = np.empty((0, 0, 0))
    freqs = np.array([])
    if need_fft:
        try:
            ts_pre = load_timeseries(patient, "rest_pre", SEEG_DATAPATH)
        except Exception as exc:  # pragma: no cover
            print(f"[{patient}] FAILED to load rest_pre time-series: {exc}")
            return pd.DataFrame()
        freqs, fft_seg = _segment_ffts(ts_pre, fs, nperseg)
        print(
            f"[{patient}] rest_pre time-series shape={ts_pre.shape}, "
            f"n_segments={fft_seg.shape[1]} (non-overlap)"
        )
    else:
        print(f"[{patient}] all null entries cached — skipping FFT precompute")

    rows: List[Dict] = []
    null_store: Dict[str, np.ndarray] = dict(cached_null)
    fc_store: Dict[Tuple[str, str], np.ndarray] = {}  # (band, phase) -> A_restricted

    for band in BRAIN_BANDS_NAMES:
        fmin, fmax = BRAIN_BANDS[band]
        adj_phases: Dict[str, np.ndarray] = {}
        missing = False
        for phase in PHASES:
            A = load_fc_matrix(patient, phase, band, "imcoh_abs")
            if A is None:
                print(f"[{patient}/{band}/{phase}] FC missing — skipping band")
                missing = True
                break
            adj_phases[phase] = np.asarray(A, dtype=np.float64)
        if missing:
            continue

        common_idx = _common_giant_indices(adj_phases)
        if common_idx.size < 10:
            print(
                f"[{patient}/{band}] |V*|={common_idx.size} < 10 — ineligible"
            )
            continue

        # Restrict observed FC to V*
        A_phase = {
            phase: adj_phases[phase][np.ix_(common_idx, common_idx)]
            for phase in PHASES
        }
        for phase in PHASES:
            fc_store[(band, phase)] = A_phase[phase]

        # Within-rest_pre split-half null — reuse cache if all 3 distances present
        if _have_cached(band):
            null = {d: cached_null[f"{band}_{d}"] for d in DISTANCES}
            print(f"[{patient}/{band}] null reused from cache")
        else:
            null = _split_half_null_from_ffts(
                fft_seg, freqs, fmin, fmax, common_idx, n_split, rng
            )
        for d_label in DISTANCES:
            null_store[f"{band}_{d_label}"] = null[d_label]

        # Observed phase-pair distances + Z
        for phase_pair in PHASE_PAIRS:
            phi_A, phi_B = phase_pair
            d_obs = all_distances(A_phase[phi_A], A_phase[phi_B])
            for d_label in DISTANCES:
                arr = null[d_label]
                arr_finite = arr[np.isfinite(arr)]
                if arr_finite.size < 3:
                    z = float("nan"); med = float("nan"); mad = float("nan")
                    flag = "degenerate_null"
                else:
                    med = float(np.median(arr_finite))
                    mad = float(np.median(np.abs(arr_finite - med)))
                    if mad == 0:
                        z = float("nan"); flag = "degenerate_null"
                    else:
                        z = (d_obs[d_label] - med) / mad
                        flag = ""
                rows.append({
                    "patient": patient,
                    "band": band,
                    "distance": d_label,
                    "phase_A": phi_A,
                    "phase_B": phi_B,
                    "d_obs": d_obs[d_label],
                    "null_median": med,
                    "null_mad": mad,
                    "z": z,
                    "n_common": int(common_idx.size),
                    "flag": flag,
                })

        # Triangle scalar T_d for this band
        for d_label in DISTANCES:
            d_taskpre = all_distances(A_phase["rest_pre"], A_phase["task_test"])[d_label]
            d_taskpost = all_distances(A_phase["task_test"], A_phase["rest_post"])[d_label]
            T = d_taskpost - d_taskpre
            rows.append({
                "patient": patient,
                "band": band,
                "distance": d_label,
                "phase_A": "TRIANGLE",
                "phase_B": "TRIANGLE",
                "d_obs": T,
                "null_median": float("nan"),
                "null_mad": float("nan"),
                "z": float("nan"),
                "n_common": int(common_idx.size),
                "flag": "T_d",
            })

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "audit.csv", index=False)
    if null_store:
        np.savez_compressed(out_dir / "null.npz", **null_store)

    if df.empty:
        print(f"[{patient}] no eligible bands — skipping diagnostic PDF")
        return df

    _render_diagnostic(patient, fc_store, null_store, df, out_dir)
    return df


# ---------------------------------------------------------------------------
# Diagnostic PDF
# ---------------------------------------------------------------------------


def _render_diagnostic(
    patient: str,
    fc_store: Dict[Tuple[str, str], np.ndarray],
    null_store: Dict[str, np.ndarray],
    df: pd.DataFrame,
    out_dir: Path,
) -> None:
    pdf_path = out_dir / "diagnostic.pdf"
    bands = [b for b in BRAIN_BANDS_NAMES if (b, "rest_pre") in fc_store]
    if not bands:
        return

    nb = len(bands)
    col_w = 2.0  # per-column width (inches)

    with PdfPages(pdf_path) as pdf:
        # ---- Page 1: imshow A^φ — rows=phases (3), cols=bands (≤6) ----
        # Per-row vmax (max across bands within each phase row) so a single
        # colorbar at the row's right end honestly represents every panel.
        nrows = len(PHASES)
        fig, axes = plt.subplots(
            nrows, nb,
            figsize=(col_w * nb + 1.0, col_w * nrows + 0.4),
            squeeze=False,
        )
        for r, phase in enumerate(PHASES):
            row_vmax = max(fc_store[(b, phase)].max() for b in bands)
            last_im = None
            for c, band in enumerate(bands):
                ax = axes[r, c]
                M = fc_store[(band, phase)]
                last_im = ax.imshow(M, cmap="magma", vmin=0, vmax=row_vmax,
                                    aspect="equal")
                ax.set_xticks([]); ax.set_yticks([])
                if r == 0:
                    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
                if c == 0:
                    ax.set_ylabel(phase, fontsize=9)
            if last_im is not None:
                _, _, clb = imshow_colorbar_caxdivider(
                    last_im, axes[r, nb - 1], size="4%", pad=0.06,
                )
                clb.ax.tick_params(labelsize=6)
        fig.text(0.99, 0.01, f"audit_25 P1 {patient}", ha="right", fontsize=6, color="gray")
        fig.tight_layout()
        pdf.savefig(fig, dpi=200); plt.close(fig)

        # ---- Page 2: difference matrices — rows=2 (task−pre, post−pre), cols=bands ----
        # Per-row symmetric vmax for the same reason as Page 1.
        nrows = 2
        fig, axes = plt.subplots(
            nrows, nb,
            figsize=(col_w * nb + 1.0, col_w * nrows + 0.4),
            squeeze=False,
        )
        diff_labels = ["task − pre", "post − pre"]
        diff_funcs = [
            lambda b: fc_store[(b, "task_test")] - fc_store[(b, "rest_pre")],
            lambda b: fc_store[(b, "rest_post")] - fc_store[(b, "rest_pre")],
        ]
        for r, (label, dfunc) in enumerate(zip(diff_labels, diff_funcs)):
            row_vmax = max(
                float(np.percentile(np.abs(dfunc(b)), 99)) for b in bands
            ) or 1e-9
            last_im = None
            for c, band in enumerate(bands):
                ax = axes[r, c]
                D = dfunc(band)
                last_im = ax.imshow(D, cmap="RdBu_r", vmin=-row_vmax,
                                    vmax=row_vmax, aspect="equal")
                ax.set_xticks([]); ax.set_yticks([])
                if r == 0:
                    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
                if c == 0:
                    ax.set_ylabel(label, fontsize=9)
            if last_im is not None:
                _, _, clb = imshow_colorbar_caxdivider(
                    last_im, axes[r, nb - 1], size="4%", pad=0.06,
                )
                clb.ax.tick_params(labelsize=6)
        fig.text(0.99, 0.01, f"audit_25 P2 {patient}", ha="right", fontsize=6, color="gray")
        fig.tight_layout()
        pdf.savefig(fig, dpi=200); plt.close(fig)

        # ---- Page 3: null vs observed — rows=3 distances (P/S/F), cols=bands ----
        nrows = len(DISTANCES)
        fig, axes = plt.subplots(
            nrows, nb,
            figsize=(col_w * nb + 0.6, 1.6 * nrows + 0.4),
            squeeze=False,
        )
        pair_colors = {("rest_pre", "task_test"): "tab:red",
                       ("rest_pre", "rest_post"): "tab:blue",
                       ("task_test", "rest_post"): "tab:green"}
        for c, band in enumerate(bands):
            for r, d_label in enumerate(DISTANCES):
                ax = axes[r, c]
                arr = null_store.get(f"{band}_{d_label}", np.array([]))
                arr_finite = arr[np.isfinite(arr)]
                if arr_finite.size > 0:
                    ax.hist(arr_finite, bins=20, color="lightgray", edgecolor="k", linewidth=0.3)
                sub = df[(df.band == band) & (df.distance == d_label) & (df.flag != "T_d")]
                for _, row in sub.iterrows():
                    pair = (row.phase_A, row.phase_B)
                    ax.axvline(
                        row.d_obs,
                        color=pair_colors.get(pair, "k"),
                        lw=1.2,
                        label=(
                            f"{pair[0][:4]}–{pair[1][:4]}: Z={row.z:.1f}"
                            if not np.isnan(row.z)
                            else f"{pair[0][:4]}–{pair[1][:4]}"
                        ),
                    )
                ax.tick_params(labelsize=6)
                if r == 0:
                    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
                if c == 0:
                    ax.set_ylabel(f"d_{d_label}", fontsize=9)
                if r == 0 and c == nb - 1:
                    ax.legend(fontsize=5, loc="upper right")
        fig.text(0.99, 0.01, f"audit_25 P3 {patient}", ha="right", fontsize=6, color="gray")
        fig.tight_layout()
        pdf.savefig(fig, dpi=200); plt.close(fig)

    print(f"[{patient}] wrote {pdf_path}")


# ---------------------------------------------------------------------------
# Cohort summary
# ---------------------------------------------------------------------------


def _render_cohort_diagnostic(
    per_patient: Dict[str, pd.DataFrame],
    summary: pd.DataFrame,
    out_path: Path,
) -> None:
    """Cohort-level Pages 4–5: n_+ bars and triangle scatter."""
    pooled = pd.concat(per_patient.values(), ignore_index=True)
    obs = pooled[pooled.flag != "T_d"].copy()
    tri = pooled[pooled.flag == "T_d"].copy()

    pair_label = {
        ("rest_pre", "task_test"): "pre→task",
        ("rest_pre", "rest_post"): "pre→post",
        ("task_test", "rest_post"): "task→post",
    }
    pair_color = {
        ("rest_pre", "task_test"): "tab:red",
        ("rest_pre", "rest_post"): "tab:blue",
        ("task_test", "rest_post"): "tab:green",
    }

    nb = len(BRAIN_BANDS_NAMES)
    bar_w = 0.25

    with PdfPages(out_path) as pdf:
        # ---- Page 4: cohort n_+ bars (rows = phase pairs, cols = bands) ----
        nrows = len(PHASE_PAIRS)
        fig, axes = plt.subplots(
            nrows, nb,
            figsize=(2.0 * nb + 0.6, 1.8 * nrows + 0.4),
            squeeze=False,
            sharey=True,
        )
        for r, phase_pair in enumerate(PHASE_PAIRS):
            phi_A, phi_B = phase_pair
            for c, band in enumerate(BRAIN_BANDS_NAMES):
                ax = axes[r, c]
                sub = summary[
                    (summary.band == band)
                    & (summary.phase_A == phi_A)
                    & (summary.phase_B == phi_B)
                ]
                xs = np.arange(len(DISTANCES))
                vals = [
                    int(sub[sub.distance == d].n_plus.iloc[0])
                    if not sub[sub.distance == d].empty
                    else 0
                    for d in DISTANCES
                ]
                bars = ax.bar(xs, vals,
                              color=pair_color.get(phase_pair, "k"),
                              edgecolor="k", linewidth=0.4, alpha=0.85)
                ax.axhline(8, color="k", lw=0.6, ls="--", alpha=0.7)
                ax.set_ylim(0, 9.3)
                ax.set_xticks(xs); ax.set_xticklabels(DISTANCES, fontsize=7)
                for v, b in zip(vals, bars):
                    ax.text(b.get_x() + b.get_width() / 2, v + 0.15, f"{v}",
                            ha="center", va="bottom", fontsize=6)
                if r == 0:
                    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
                if c == 0:
                    ax.set_ylabel(f"n+ — {pair_label[phase_pair]}", fontsize=8)
                ax.tick_params(labelsize=6)
        fig.text(0.99, 0.01, "audit_25 P4 cohort", ha="right", fontsize=6, color="gray")
        fig.tight_layout()
        pdf.savefig(fig, dpi=200); plt.close(fig)

        # ---- Page 5: cohort triangle scatter (rows = distances, cols = bands) ----
        nrows = len(DISTANCES)
        fig, axes = plt.subplots(
            nrows, nb,
            figsize=(2.0 * nb + 0.6, 1.9 * nrows + 0.4),
            squeeze=False,
        )
        for r, d_label in enumerate(DISTANCES):
            for c, band in enumerate(BRAIN_BANDS_NAMES):
                ax = axes[r, c]
                # x = d(pre, task), y = d(task, post), per patient
                d_pre_task = obs[
                    (obs.band == band) & (obs.distance == d_label)
                    & (obs.phase_A == "rest_pre") & (obs.phase_B == "task_test")
                ][["patient", "d_obs"]].rename(columns={"d_obs": "x"})
                d_task_post = obs[
                    (obs.band == band) & (obs.distance == d_label)
                    & (obs.phase_A == "task_test") & (obs.phase_B == "rest_post")
                ][["patient", "d_obs"]].rename(columns={"d_obs": "y"})
                merged = d_pre_task.merge(d_task_post, on="patient", how="inner")
                if merged.empty:
                    ax.set_xticks([]); ax.set_yticks([]); continue

                # All 10 patients in pool. Pat_03 highlighted as orange
                # triangle so the 1024 Hz case stays visually identifiable
                # without affecting the cohort count.
                in_pool = merged[merged.patient != "Pat_03"]
                pat03 = merged[merged.patient == "Pat_03"]

                ax.scatter(in_pool.x, in_pool.y, color="tab:blue",
                           edgecolor="k", linewidth=0.3, s=24, zorder=3)
                if not pat03.empty:
                    ax.scatter(pat03.x, pat03.y, marker="^",
                               color="tab:orange", edgecolor="k", linewidth=0.3,
                               s=30, zorder=4, label="Pat_03")
                # Identity line
                lo = float(min(merged.x.min(), merged.y.min()))
                hi = float(max(merged.x.max(), merged.y.max()))
                ax.plot([lo, hi], [lo, hi], color="k", lw=0.6, ls="--", alpha=0.6)
                ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
                ax.tick_params(labelsize=6)
                # Annotate fraction below diagonal (T_d < 0)
                n_below = int((merged.y < merged.x).sum())
                n_total = int(merged.shape[0])
                ax.text(0.04, 0.92, f"{n_below}/{n_total}",
                        transform=ax.transAxes, fontsize=7,
                        ha="left", va="top",
                        bbox=dict(facecolor="white", edgecolor="0.6",
                                  boxstyle="round,pad=0.15", alpha=0.85))
                if r == 0:
                    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
                if c == 0:
                    ax.set_ylabel(f"d_{d_label}\nd(task,post)", fontsize=8)
                if r == nrows - 1:
                    ax.set_xlabel("d(pre,task)", fontsize=7)
        fig.text(0.99, 0.01, "audit_25 P5 cohort", ha="right", fontsize=6, color="gray")
        fig.tight_layout()
        pdf.savefig(fig, dpi=200); plt.close(fig)


def _cohort_summary(per_patient: Dict[str, pd.DataFrame], out_path: Path) -> None:
    rows: List[Dict] = []
    pooled = pd.concat(per_patient.values(), ignore_index=True)
    pooled = pooled[pooled.flag != "T_d"]  # exclude triangle rows for n_+
    # Pat_03 included: Z-score is dimensionless w.r.t. fs, so the 1024 Hz
    # outlier is comparable to in-pool patients on this analysis.
    in_pool = pooled

    for band in BRAIN_BANDS_NAMES:
        for d_label in DISTANCES:
            for phase_pair in PHASE_PAIRS:
                phi_A, phi_B = phase_pair
                sub = in_pool[
                    (in_pool.band == band)
                    & (in_pool.distance == d_label)
                    & (in_pool.phase_A == phi_A)
                    & (in_pool.phase_B == phi_B)
                ]
                z = sub.z.dropna()
                n_eligible = z.size
                n_plus = int((z > 2).sum())
                rows.append({
                    "band": band,
                    "distance": d_label,
                    "phase_A": phi_A,
                    "phase_B": phi_B,
                    "n_eligible": n_eligible,
                    "n_plus": n_plus,
                    "verdict": "positive" if n_plus >= 8 else "negative",
                })

    # Triangle metrics — two complementary counts:
    #   TRIANGLE     : T_d < 0 alone (geometric persistence).
    #   TRIANGLE_SIG : T_d < 0 AND Z(pre,post) > 2 (geometric + significant
    #                  pre→post move). Stricter, can suppress real cells
    #                  where pre/post sit close to each other but post is
    #                  still on the task side of pre.
    pooled_T = pd.concat(per_patient.values(), ignore_index=True)
    for band in BRAIN_BANDS_NAMES:
        for d_label in DISTANCES:
            tri = pooled_T[
                (pooled_T.flag == "T_d")
                & (pooled_T.band == band)
                & (pooled_T.distance == d_label)
            ]
            zpost = pooled_T[
                (pooled_T.flag != "T_d")
                & (pooled_T.band == band)
                & (pooled_T.distance == d_label)
                & (pooled_T.phase_A == "rest_pre")
                & (pooled_T.phase_B == "rest_post")
            ][["patient", "z"]].rename(columns={"z": "z_post"})
            merged = tri.merge(zpost, on="patient", how="left")
            n_neg = int((merged.d_obs < 0).sum())
            n_persist = int(((merged.d_obs < 0) & (merged.z_post > 2)).sum())
            n_eligible = int(merged.shape[0])
            rows.append({
                "band": band,
                "distance": d_label,
                "phase_A": "TRIANGLE",
                "phase_B": "TRIANGLE",
                "n_eligible": n_eligible,
                "n_plus": n_neg,
                "verdict": "positive" if n_neg >= 6 else "negative",
            })
            rows.append({
                "band": band,
                "distance": d_label,
                "phase_A": "TRIANGLE_SIG",
                "phase_B": "TRIANGLE_SIG",
                "n_eligible": n_eligible,
                "n_plus": n_persist,
                "verdict": "positive" if n_persist >= 6 else "negative",
            })

    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(out_path, index=False)
    print(f"wrote {out_path}")

    cohort_pdf = out_path.parent / "cohort_diagnostic.pdf"
    _render_cohort_diagnostic(per_patient, summary_df, cohort_pdf)
    print(f"wrote {cohort_pdf}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", default="Pat_06")
    ap.add_argument("--n-split", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cohort", action="store_true",
                    help="Run all 10 PATIENTS_4PHASE patients sequentially.")
    ap.add_argument("--summarize-only", action="store_true",
                    help="Skip per-patient compute; aggregate existing "
                         "<Patient>/audit.csv into cohort_summary.csv.")
    ap.add_argument("--no-reuse-null", action="store_true",
                    help="Force recompute of the within-rest_pre null even if "
                         "<Patient>/null.npz exists. Default: reuse cache.")
    args = ap.parse_args()

    OUT_BASE.mkdir(parents=True, exist_ok=True)

    if args.summarize_only:
        per_patient: Dict[str, pd.DataFrame] = {}
        for p in PATIENTS_4PHASE:
            csv = OUT_BASE / p / "audit.csv"
            if csv.exists():
                per_patient[p] = pd.read_csv(csv)
        if per_patient:
            _cohort_summary(per_patient, OUT_BASE / "cohort_summary.csv")
        return 0

    if args.cohort:
        targets = list(PATIENTS_4PHASE)
    else:
        targets = [args.patient]

    per_patient = {}
    for p in targets:
        out_dir = OUT_BASE / p
        df = run_patient(
            p, args.n_split, args.seed, out_dir,
            reuse_null=not args.no_reuse_null,
        )
        if not df.empty:
            per_patient[p] = df

    if args.cohort and per_patient:
        _cohort_summary(per_patient, OUT_BASE / "cohort_summary.csv")

    return 0


if __name__ == "__main__":
    sys.exit(main())
