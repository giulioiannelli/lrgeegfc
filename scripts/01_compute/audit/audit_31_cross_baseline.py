#!/usr/bin/env python3
"""Audit 31 — cross-baseline distance (Control 3).

Maximum-time-gap pseudo-cross-phase distance between the **late half** of
``rest_pre`` and the **early half** of ``rest_post`` — no task between
the halves. Under monotone within-session drift this distance should be
comparable to ``d(rest_pre, rest_post)``; under task-specific
reorganization it should be substantially **smaller** than
``d(task_test, rest_post)``.

For each (patient, band, distance):
- Split ``rest_pre`` segment-FFTs into early/late halves by index.
- Same for ``rest_post``.
- Compute |ImCoh| for each half restricted to
  ``V*(p, b) = ⋂_φ GC(|A^φ|)``.
- Compute:
    * ``d_pre_late_post_early`` = d(rest_pre.late, rest_post.early)
    * ``d_pre_early_pre_late``  = d(rest_pre.early, rest_pre.late)
    * ``d_post_early_post_late`` = d(rest_post.early, rest_post.late)
- Reload cached observed full-phase distances:
    * ``d_pre_full_post_full`` = d(rest_pre, rest_post)
    * ``d_taskTest_rsPost`` = d(task_test, rest_post)

Decision rule: per band, count patients with
``d(taskTest, rsPost) < d(rsPre.late, rsPost.early)``. The trace is
task-specific in bands where this count is ≥ 7/10.

Outputs
-------
``data/audit/raw_fc_phase_distance/cohort_cross_baseline.csv``
``data/audit/raw_fc_phase_distance/cross_baseline_scatter.pdf``
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
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch
from scipy.signal import get_window
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
from lrg_eegfc.workflow.fc import load_fc_matrix

PHASES_FULL: Tuple[str, ...] = ("rest_pre", "task_test", "rest_post")
DISTANCES: Tuple[str, ...] = ("P", "S", "F")
OUT_BASE = DATA_ROOT / "audit" / "raw_fc_phase_distance"
TRACE_BANDS = ("delta", "alpha", "beta", "low_gamma")


def _triu_offdiag(A: np.ndarray) -> np.ndarray:
    return A[np.triu_indices_from(A, k=1)]


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


_DIST_FUNCS = {"P": d_P, "S": d_S, "F": d_F}


def _common_giant_indices(adj_phases: Dict[str, np.ndarray]) -> np.ndarray:
    sets: List[set] = []
    for adj in adj_phases.values():
        G = nx.from_numpy_array(np.abs(adj))
        comps = list(nx.connected_components(G))
        sets.append(set(max(comps, key=len)) if comps else set())
    if not sets:
        return np.array([], dtype=int)
    common = set.intersection(*sets)
    return np.array(sorted(common), dtype=int)


def _segment_ffts(
    ts: np.ndarray, fs: float, nperseg: int
) -> Tuple[np.ndarray, np.ndarray]:
    N, L = ts.shape
    n_seg = L // nperseg
    if n_seg < 4:
        return np.array([]), np.empty((N, 0, 0))
    ts_trim = ts[:, : n_seg * nperseg].astype(np.float64, copy=False)
    segs = ts_trim.reshape(N, n_seg, nperseg)
    window = get_window("hann", nperseg).astype(np.float64)
    block = segs * window[None, None, :]
    fft_seg = np.fft.rfft(block, n=nperseg, axis=2)
    freqs = np.fft.rfftfreq(nperseg, 1.0 / fs)
    return freqs, fft_seg


def _imcoh_abs_from_subset(
    fft_seg: np.ndarray, subset: np.ndarray, freq_mask: np.ndarray
) -> np.ndarray:
    sub = fft_seg[:, subset, :][:, :, freq_mask]
    sub_conj = np.conj(sub)
    cross = np.einsum("nsf,msf->nmf", sub, sub_conj)
    auto = np.real(np.einsum("nsf,nsf->nf", sub, sub_conj))
    denom = np.sqrt(auto[:, None, :] * auto[None, :, :])
    with np.errstate(invalid="ignore", divide="ignore"):
        signed = np.where(denom > 0, np.imag(cross) / denom, 0.0)
    A = np.abs(signed).mean(axis=-1)
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    return A


def _v_star_per_band(patient: str) -> Tuple[Dict[str, np.ndarray], Dict[str, Dict[str, np.ndarray]]]:
    """Return per-band V* indices and the cached full-phase FCs (V*-restricted)."""
    v_star: Dict[str, np.ndarray] = {}
    fc_full: Dict[str, Dict[str, np.ndarray]] = {}  # band -> phase -> A_restricted
    for band in BRAIN_BANDS_NAMES:
        adj_phases: Dict[str, np.ndarray] = {}
        skip = False
        for phase in PHASES_FULL:
            A = load_fc_matrix(patient, phase, band, "imcoh_abs")
            if A is None:
                skip = True; break
            adj_phases[phase] = np.asarray(A, dtype=np.float64)
        if skip:
            continue
        idx = _common_giant_indices(adj_phases)
        if idx.size < 10:
            continue
        v_star[band] = idx
        fc_full[band] = {
            phase: adj_phases[phase][np.ix_(idx, idx)] for phase in PHASES_FULL
        }
    return v_star, fc_full


def run_patient(patient: str) -> pd.DataFrame:
    fs = FS_OVERRIDES.get(patient, DEFAULT_SAMPLE_RATE)
    nperseg = nperseg_for_fs(fs)
    print(f"[{patient}] fs={fs:.1f} Hz, nperseg={nperseg}")

    v_star, fc_full = _v_star_per_band(patient)
    if not v_star:
        print(f"[{patient}] no eligible bands")
        return pd.DataFrame()

    # Load and FFT rest_pre, rest_post, task_test; split each into halves.
    # task_test is included so we can build a noise-symmetric counterpart
    # to d(rsPre.late, rsPost.early): d(taskTest.late, rsPost.early) uses
    # the same half-segment budget on both sides.
    half_ffts: Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]] = {}
    for phase in ("rest_pre", "rest_post", "task_test"):
        try:
            ts = load_timeseries(patient, phase, SEEG_DATAPATH)
        except Exception as exc:
            print(f"[{patient}/{phase}] load failed: {exc}")
            return pd.DataFrame()
        freqs, fft_seg = _segment_ffts(ts, fs, nperseg)
        if fft_seg.size == 0:
            print(f"[{patient}/{phase}] no segments")
            return pd.DataFrame()
        n_seg = fft_seg.shape[1]
        half = n_seg // 2
        if half < 4:
            print(f"[{patient}/{phase}] insufficient segments (n_seg={n_seg})")
            return pd.DataFrame()
        early_idx = np.arange(0, half)
        late_idx = np.arange(n_seg - half, n_seg)
        half_ffts[phase] = {
            "early": (fft_seg, freqs, early_idx),
            "late": (fft_seg, freqs, late_idx),
        }
        print(f"[{patient}/{phase}] n_seg={n_seg}, half={half}")

    rows: List[Dict] = []
    for band in BRAIN_BANDS_NAMES:
        if band not in v_star:
            continue
        common_idx = v_star[band]
        fmin, fmax = BRAIN_BANDS[band]

        # Build half-FCs for all 3 phases
        half_FC: Dict[str, Dict[str, np.ndarray]] = {}  # phase -> {early|late}: A
        for phase in ("rest_pre", "rest_post", "task_test"):
            half_FC[phase] = {}
            for tag, (fft_seg, freqs, seg_idx) in half_ffts[phase].items():
                freq_mask = (freqs >= fmin) & (freqs < fmax)
                if not freq_mask.any():
                    idx = int(np.clip(np.searchsorted(freqs, fmin), 0, len(freqs) - 1))
                    freq_mask = np.zeros_like(freqs, dtype=bool); freq_mask[idx] = True
                A = _imcoh_abs_from_subset(fft_seg, seg_idx, freq_mask)
                half_FC[phase][tag] = A[np.ix_(common_idx, common_idx)]

        for d_label, fn in _DIST_FUNCS.items():
            # Half-pair distances (all use the same half-segment noise budget)
            d_pre_late_post_early = fn(
                half_FC["rest_pre"]["late"], half_FC["rest_post"]["early"]
            )
            d_taskLate_post_early = fn(
                half_FC["task_test"]["late"], half_FC["rest_post"]["early"]
            )
            d_pre_within = fn(
                half_FC["rest_pre"]["early"], half_FC["rest_pre"]["late"]
            )
            d_post_within = fn(
                half_FC["rest_post"]["early"], half_FC["rest_post"]["late"]
            )
            d_task_within = fn(
                half_FC["task_test"]["early"], half_FC["task_test"]["late"]
            )
            # Full-phase distances (lower noise; reference only)
            d_pre_full_post_full = fn(fc_full[band]["rest_pre"], fc_full[band]["rest_post"])
            d_task_post = fn(fc_full[band]["task_test"], fc_full[band]["rest_post"])
            d_pre_task = fn(fc_full[band]["rest_pre"], fc_full[band]["task_test"])
            rows.append({
                "patient": patient,
                "band": band,
                "distance": d_label,
                "n_common": int(common_idx.size),
                # Symmetric-noise half-pair comparison (LOAD-BEARING for the
                # decision rule). Both sides are half-segments → same noise floor.
                "d_pre_late_post_early": d_pre_late_post_early,
                "d_taskLate_post_early": d_taskLate_post_early,
                # Within-baseline half-pair (estimation noise floor per phase)
                "d_pre_within": d_pre_within,
                "d_post_within": d_post_within,
                "d_task_within": d_task_within,
                # Full-phase distances (reference; do not use for the decision)
                "d_pre_full_post_full": d_pre_full_post_full,
                "d_taskTest_rsPost_full": d_task_post,
                "d_rsPre_taskTest_full": d_pre_task,
                # Symmetric decision: task half-pair < baseline-baseline half-pair
                "task_specific_symmetric": int(
                    d_taskLate_post_early < d_pre_late_post_early
                ),
                # Legacy asymmetric flag (kept for diagnostic continuity; biased
                # by noise asymmetry — full-phase d(TT, RPost) is less noisy
                # than half-half d(RPre.late, RPost.early), so this flag inflates)
                "task_specific_asymmetric_legacy": int(
                    d_task_post < d_pre_late_post_early
                ),
            })
    return pd.DataFrame(rows)


def _render_scatter(merged: pd.DataFrame, out_pdf: Path) -> None:
    """Per-band scatter on d_S (load-bearing): x = d(rsPre.late, rsPost.early),
    y = d(taskTest, rsPost). Below identity = task-specific trace direction.
    Single row, one panel per band. Style matches the bundle Fig 1.
    """
    bands = list(BRAIN_BANDS_NAMES)
    nb = len(bands)
    n_cohort = len(PATIENTS_4PHASE)
    sub_S = merged[merged.distance == "S"]

    fig, axes = plt.subplots(1, nb, figsize=(2.6 * nb + 0.0, 2.95),
                             squeeze=False)
    for c, band in enumerate(bands):
        ax = axes[0, c]
        sub = sub_S[sub_S.band == band].copy()
        if sub.empty:
            ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
            ax.set_xticks([]); ax.set_yticks([])
            continue
        pat_codes = sub.patient.values
        # SYMMETRIC half-pair comparison: same noise budget on both sides
        x = sub.d_pre_late_post_early.values.astype(float)
        y = sub.d_taskLate_post_early.values.astype(float)
        finite = np.isfinite(x) & np.isfinite(y)
        pat_codes = pat_codes[finite]; x = x[finite]; y = y[finite]
        if x.size == 0:
            ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
            continue

        span = max(x.max(), y.max()) - min(x.min(), y.min())
        lo = min(x.min(), y.min()) - 0.05 * span
        hi = max(x.max(), y.max()) + 0.05 * span

        # Trace zone (below identity) — light red shade
        poly_x = [lo, hi, hi]; poly_y = [lo, lo, hi]
        ax.fill(poly_x, poly_y, color="#d62728", alpha=0.06, zorder=0)

        # Identity
        ax.plot([lo, hi], [lo, hi], color="k", lw=0.7, ls="--",
                alpha=0.55, zorder=1)

        # Cohort 1σ covariance ellipse + mean star
        if x.size >= 2:
            mu_x, mu_y = float(x.mean()), float(y.mean())
            cov = np.cov(x, y)
            vals, vecs = np.linalg.eigh(cov)
            order = np.argsort(vals)[::-1]
            vals = vals[order]; vecs = vecs[:, order]
            angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
            width, height = 2.0 * np.sqrt(np.maximum(vals, 0.0))
            if width > 0 and height > 0:
                ell = Ellipse(
                    (mu_x, mu_y), width, height, angle=angle,
                    facecolor="#888", alpha=0.18, edgecolor="#444",
                    linewidth=0.7, zorder=2,
                )
                ax.add_patch(ell)
            ax.scatter([mu_x], [mu_y], marker="*", s=180, c="#2ca02c",
                       edgecolor="k", linewidth=0.9, zorder=5)

        # Per-patient circles; Pat_03 distinct triangle
        is_p03 = pat_codes == "Pat_03"
        ax.scatter(x[~is_p03], y[~is_p03], c="white",
                   edgecolor="#1f77b4", linewidth=1.1, s=34, zorder=3)
        if is_p03.any():
            ax.scatter(x[is_p03], y[is_p03], marker="^", c="white",
                       edgecolor="#ff7f0e", linewidth=1.2, s=46, zorder=4)

        # n_trace annotation
        n_trace = int((y < x).sum())
        n_total = int(x.size)
        ax.text(0.05, 0.95, f"trace: {n_trace}/{n_total}",
                transform=ax.transAxes, fontsize=8.5, va="top",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="#999", alpha=0.92))

        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
        ax.set_xlabel(r"$d_S$(RPre.late, RPost.early)", fontsize=8.5)
        if c == 0:
            ax.set_ylabel(r"$d_S$(TT.late, RPost.early)", fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(color="#eee", lw=0.4, zorder=0)
        ax.set_axisbelow(True)

    handles = [
        Line2D([0], [0], marker="*", color="none",
               markerfacecolor="#2ca02c", markeredgecolor="k",
               markersize=12, label="cohort mean"),
        Patch(facecolor="#888", alpha=0.18, edgecolor="#444",
              label=r"cohort 1$\sigma$ ellipse"),
        Patch(facecolor="#d62728", alpha=0.10,
              label="task-specific zone (below identity)"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="white", markeredgecolor="#1f77b4",
               markersize=8, label=f"patient (n={n_cohort - 1})"),
        Line2D([0], [0], marker="^", color="none",
               markerfacecolor="white", markeredgecolor="#ff7f0e",
               markersize=9, label="Pat_03 (1024 Hz)"),
    ]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.03), ncol=len(handles),
               frameon=False, fontsize=8.2)
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


def cohort_summary(per_patient: Dict[str, pd.DataFrame], out_base: Path) -> None:
    pooled = pd.concat(per_patient.values(), ignore_index=True)
    out_csv = out_base / "cohort_cross_baseline.csv"
    pooled.to_csv(out_csv, index=False)
    print(f"wrote {out_csv}")

    rows = []
    for band in BRAIN_BANDS_NAMES:
        for d_label in DISTANCES:
            sub = pooled[(pooled.band == band) & (pooled.distance == d_label)]
            if sub.empty:
                continue
            n_tot = int(sub.shape[0])
            n_task_below_baseline = int(sub.task_specific_symmetric.sum())
            rows.append({
                "band": band,
                "distance": d_label,
                "n_patients": n_tot,
                "n_task_below_baseline_baseline": n_task_below_baseline,
                "task_specific_pass": int(n_task_below_baseline >= 7),
                "median_d_pre_late_post_early": round(float(sub.d_pre_late_post_early.median()), 4),
                "median_d_taskLate_post_early": round(float(sub.d_taskLate_post_early.median()), 4),
                "median_d_pre_within": round(float(sub.d_pre_within.median()), 4),
                "median_d_post_within": round(float(sub.d_post_within.median()), 4),
                "median_d_task_within": round(float(sub.d_task_within.median()), 4),
                "median_d_pre_full_post_full": round(float(sub.d_pre_full_post_full.median()), 4),
                "median_d_taskTest_rsPost_full": round(float(sub.d_taskTest_rsPost_full.median()), 4),
                "median_d_rsPre_taskTest_full": round(float(sub.d_rsPre_taskTest_full.median()), 4),
                "n_task_below_baseline_baseline_legacy_asym": int(sub.task_specific_asymmetric_legacy.sum()),
            })
    band_summary = pd.DataFrame(rows)
    band_csv = out_base / "cohort_cross_baseline_band_summary.csv"
    band_summary.to_csv(band_csv, index=False)
    print(f"wrote {band_csv}")
    print(f"\n=== Cross-baseline cohort summary (n={len(per_patient)} patients) ===")
    print(band_summary.to_string(index=False))

    _render_scatter(pooled, out_base / "cross_baseline_scatter.pdf")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", default=None)
    ap.add_argument("--cohort", action="store_true")
    args = ap.parse_args()
    OUT_BASE.mkdir(parents=True, exist_ok=True)

    targets = list(PATIENTS_4PHASE) if args.cohort else [args.patient or PATIENTS_4PHASE[0]]
    per_patient: Dict[str, pd.DataFrame] = {}
    for p in targets:
        df = run_patient(p)
        if not df.empty:
            per_patient[p] = df

    if args.cohort and per_patient:
        cohort_summary(per_patient, OUT_BASE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
