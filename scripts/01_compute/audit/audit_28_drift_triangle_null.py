#!/usr/bin/env python3
"""Audit 28 — drift-triangle null for raw-FC phase distances.

Implements the scope at
``.agents/guides/task-persistence-investigation/2026-04-29_drift-triangle-null.md``.

For each patient × band × distance:

1. Split ``rest_pre`` into K consecutive non-overlapping temporal chunks;
   same for ``rest_post``.
2. Compute |ImCoh| per chunk per band, restricted to the per-patient
   ``V*(p, b) = ⋂_{φ ∈ {rest_pre, task_test, rest_post}} GC(|A^φ|)``.
3. Compute pairwise within-phase distances ``D^φ[k,l]`` and temporal
   centroid gaps ``gap^φ[k,l]``.
4. Fit ``D ≈ α + β·gap`` on pooled (rest_pre, rest_post) points.
5. Predict ``d_drift_pred`` at the inter-rest centroid gap.
6. Compute residual = observed − predicted for ``(pre, task)`` and
   ``(pre, post)`` pairs.

Usage
-----
    python audit_28_drift_triangle_null.py                       # Pat_06 default
    python audit_28_drift_triangle_null.py --patient Pat_02
    python audit_28_drift_triangle_null.py --cohort
    python audit_28_drift_triangle_null.py --K 6                 # K = chunks per phase

Outputs
-------
``data/audit/raw_fc_phase_distance/<Patient>/drift_triangle.csv``
``data/audit/raw_fc_phase_distance/<Patient>/drift_triangle.pdf``
``data/audit/raw_fc_phase_distance/cohort_drift_triangle.csv`` (cohort)
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

PHASES_REST: Tuple[str, ...] = ("rest_pre", "rest_post")
PHASES_FULL: Tuple[str, ...] = ("rest_pre", "task_test", "rest_post")
DISTANCES: Tuple[str, ...] = ("P", "S", "F")
OUT_BASE = DATA_ROOT / "audit" / "raw_fc_phase_distance"


# ---------------------------------------------------------------------------
# Distance functions (same as audit_25)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Substrate (V*) — same as audit_25 (intersection of giant components)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# FFT + ImCoh on chunked time-series
# ---------------------------------------------------------------------------

def _segment_ffts(
    ts: np.ndarray, fs: float, nperseg: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Hann-windowed non-overlapping segment FFTs. Returns (freqs, fft_seg)
    where fft_seg has shape (N, n_seg, F)."""
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


def _chunk_segment_indices(n_seg: int, K: int) -> List[np.ndarray]:
    """Split n_seg consecutive segments into K equal consecutive chunks
    by integer division; remainder is dropped from the tail."""
    per = n_seg // K
    if per < 4:
        return []
    return [np.arange(k * per, (k + 1) * per) for k in range(K)]


# ---------------------------------------------------------------------------
# Per-patient drift-triangle pipeline
# ---------------------------------------------------------------------------

def _load_full_phase_FCs(patient: str) -> Dict[str, np.ndarray]:
    out: Dict[str, np.ndarray] = {}
    for phase in PHASES_FULL:
        for band in BRAIN_BANDS_NAMES:
            A = load_fc_matrix(patient, phase, band, "imcoh_abs")
            if A is None:
                continue
            out[f"{phase}_{band}"] = np.asarray(A, dtype=np.float64)
    return out


def _common_v_star_per_band(patient: str) -> Dict[str, np.ndarray]:
    """V*(p, b) per band — intersection of GC across all 3 phases at that band."""
    v_star: Dict[str, np.ndarray] = {}
    for band in BRAIN_BANDS_NAMES:
        adj_phases = {}
        skip = False
        for phase in PHASES_FULL:
            A = load_fc_matrix(patient, phase, band, "imcoh_abs")
            if A is None:
                skip = True; break
            adj_phases[phase] = np.asarray(A, dtype=np.float64)
        if skip:
            continue
        v_star[band] = _common_giant_indices(adj_phases)
    return v_star


def run_patient(patient: str, K: int, out_dir: Path) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    fs = FS_OVERRIDES.get(patient, DEFAULT_SAMPLE_RATE)
    nperseg = nperseg_for_fs(fs)

    print(f"[{patient}] fs={fs:.1f} Hz, nperseg={nperseg}, K={K}")

    # Per-band V* (intersection of giant components across rest_pre/task_test/rest_post)
    v_star = _common_v_star_per_band(patient)
    if not v_star:
        print(f"[{patient}] no eligible bands (V* missing)")
        return pd.DataFrame()

    # Load rest_pre and rest_post time-series and compute segment FFTs
    fft_per_phase: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    duration_s: Dict[str, float] = {}
    for phase in PHASES_REST:
        try:
            ts = load_timeseries(patient, phase, SEEG_DATAPATH)
        except Exception as exc:
            print(f"[{patient}/{phase}] load failed: {exc}")
            return pd.DataFrame()
        duration_s[phase] = ts.shape[1] / fs
        freqs, fft_seg = _segment_ffts(ts, fs, nperseg)
        if fft_seg.size == 0:
            print(f"[{patient}/{phase}] no segments")
            return pd.DataFrame()
        fft_per_phase[phase] = (freqs, fft_seg)
        print(f"[{patient}/{phase}] T={duration_s[phase]:.1f}s, "
              f"n_seg={fft_seg.shape[1]}")

    # Load inter-rest gap (task_learn + task_test duration)
    inter_rest_gap_s = 0.0
    for tphase in ("task_learn", "task_test"):
        try:
            ts_t = load_timeseries(patient, tphase, SEEG_DATAPATH)
            inter_rest_gap_s += ts_t.shape[1] / fs
        except Exception as exc:
            print(f"[{patient}/{tphase}] load failed: {exc}")

    # Build chunk segment indices per phase
    chunk_idx: Dict[str, List[np.ndarray]] = {}
    chunk_T_s: Dict[str, float] = {}
    for phase in PHASES_REST:
        freqs, fft_seg = fft_per_phase[phase]
        n_seg = fft_seg.shape[1]
        chunks = _chunk_segment_indices(n_seg, K)
        if not chunks:
            print(f"[{patient}/{phase}] insufficient segments per chunk (n_seg={n_seg}, K={K})")
            return pd.DataFrame()
        chunk_idx[phase] = chunks
        # Each segment is nperseg/fs seconds; chunk duration = per_chunk_segs * nperseg/fs
        chunk_T_s[phase] = (n_seg // K) * (nperseg / fs)
        print(f"[{patient}/{phase}] chunk_T={chunk_T_s[phase]:.1f}s, "
              f"K={K}, segs/chunk={n_seg // K}")

    # Per-band: compute chunk FCs, within-phase D[k,l], drift fit, predicted, residuals
    rows: List[Dict] = []
    diagnostic_packs: Dict[str, Dict] = {}  # for PDF

    for band in BRAIN_BANDS_NAMES:
        if band not in v_star:
            continue
        common_idx = v_star[band]
        if common_idx.size < 10:
            continue
        fmin, fmax = BRAIN_BANDS[band]

        # Compute per-chunk FC restricted to V*
        A_chunk: Dict[str, List[np.ndarray]] = {}
        for phase in PHASES_REST:
            freqs, fft_seg = fft_per_phase[phase]
            freq_mask = (freqs >= fmin) & (freqs < fmax)
            if not freq_mask.any():
                idx = int(np.clip(np.searchsorted(freqs, fmin), 0, len(freqs) - 1))
                freq_mask = np.zeros_like(freqs, dtype=bool); freq_mask[idx] = True
            A_chunk[phase] = []
            for ck in chunk_idx[phase]:
                A_full = _imcoh_abs_from_subset(fft_seg, ck, freq_mask)
                A_chunk[phase].append(A_full[np.ix_(common_idx, common_idx)])

        # Pairwise within-phase D and gaps; pool rest_pre + rest_post
        gap_pre = chunk_T_s["rest_pre"]
        gap_post = chunk_T_s["rest_post"]
        for d_label in DISTANCES:
            df_func = _DIST_FUNCS[d_label]
            xs: List[float] = []
            ys: List[float] = []
            phase_tag: List[str] = []
            for phase, gap_s in (("rest_pre", gap_pre), ("rest_post", gap_post)):
                K_phase = len(A_chunk[phase])
                for k in range(K_phase):
                    for l in range(k + 1, K_phase):
                        gap = (l - k) * gap_s  # centroid gap
                        d = df_func(A_chunk[phase][k], A_chunk[phase][l])
                        xs.append(gap); ys.append(d); phase_tag.append(phase)
            xs_arr = np.asarray(xs); ys_arr = np.asarray(ys)
            mask = np.isfinite(ys_arr)
            if mask.sum() < 4:
                continue
            xv = xs_arr[mask]; yv = ys_arr[mask]
            # Linear fit D = a + b*gap
            try:
                slope, intercept = np.polyfit(xv, yv, 1)
                ypred = intercept + slope * xv
                ss_res = np.sum((yv - ypred) ** 2)
                ss_tot = np.sum((yv - np.mean(yv)) ** 2)
                R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
            except Exception:
                slope = float("nan"); intercept = float("nan"); R2 = float("nan")

            # Inter-rest centroid gap = T_pre/2 + inter_rest_gap + T_post/2
            gap_cross = duration_s["rest_pre"] / 2 + inter_rest_gap_s + duration_s["rest_post"] / 2
            # Pre→task gap = T_pre/2 + T_learn + T_test/2 ≈ T_pre/2 + inter/2 (use centroid of task_test)
            # We don't have task_test duration here easily; approximate gap_pre_task as T_pre/2 + T_learn + T_test/2.
            # Use Q5: load task durations.
            # (Computed below from time-series)
            d_drift_at_cross = intercept + slope * gap_cross

            # Observed cross-phase distances from cached full-phase FCs
            A_pre_full = load_fc_matrix(patient, "rest_pre", band, "imcoh_abs")
            A_task_full = load_fc_matrix(patient, "task_test", band, "imcoh_abs")
            A_post_full = load_fc_matrix(patient, "rest_post", band, "imcoh_abs")
            A_pre_full = np.asarray(A_pre_full, dtype=np.float64)[np.ix_(common_idx, common_idx)]
            A_task_full = np.asarray(A_task_full, dtype=np.float64)[np.ix_(common_idx, common_idx)]
            A_post_full = np.asarray(A_post_full, dtype=np.float64)[np.ix_(common_idx, common_idx)]
            d_obs_pre_post = df_func(A_pre_full, A_post_full)
            d_obs_pre_task = df_func(A_pre_full, A_task_full)
            d_obs_task_post = df_func(A_task_full, A_post_full)

            # Compute pre→task gap: end of rest_pre to centroid of task_test ≈ T_pre/2 + T_learn + T_test/2
            try:
                ts_learn = load_timeseries(patient, "task_learn", SEEG_DATAPATH)
                T_learn = ts_learn.shape[1] / fs
                ts_test = load_timeseries(patient, "task_test", SEEG_DATAPATH)
                T_test = ts_test.shape[1] / fs
            except Exception:
                T_learn = inter_rest_gap_s / 2; T_test = inter_rest_gap_s / 2

            gap_pre_task = duration_s["rest_pre"] / 2 + T_learn + T_test / 2
            gap_task_post = T_test / 2 + duration_s["rest_post"] / 2
            d_drift_pre_task = intercept + slope * gap_pre_task
            d_drift_task_post = intercept + slope * gap_task_post

            res_pre_post = d_obs_pre_post - d_drift_at_cross
            res_pre_task = d_obs_pre_task - d_drift_pre_task
            res_task_post = d_obs_task_post - d_drift_task_post

            rows.append({
                "patient": patient,
                "band": band,
                "distance": d_label,
                "K": K,
                "n_pairs_pooled": int(mask.sum()),
                "drift_intercept": float(intercept),
                "drift_slope_per_s": float(slope),
                "drift_R2": float(R2),
                "gap_cross_s": float(gap_cross),
                "gap_pre_task_s": float(gap_pre_task),
                "gap_task_post_s": float(gap_task_post),
                "d_drift_pred_pre_post": float(d_drift_at_cross),
                "d_drift_pred_pre_task": float(d_drift_pre_task),
                "d_drift_pred_task_post": float(d_drift_task_post),
                "d_obs_pre_post": float(d_obs_pre_post),
                "d_obs_pre_task": float(d_obs_pre_task),
                "d_obs_task_post": float(d_obs_task_post),
                "residual_pre_post": float(res_pre_post),
                "residual_pre_task": float(res_pre_task),
                "residual_task_post": float(res_task_post),
                "T_d_obs": float(d_obs_task_post - d_obs_pre_task),
                "T_d_drift": float(d_drift_task_post - d_drift_pre_task),
            })

            diagnostic_packs[f"{band}_{d_label}"] = {
                "xv": xv, "yv": yv, "phase_tag": np.array(phase_tag)[mask],
                "slope": slope, "intercept": intercept, "R2": R2,
                "gap_cross": gap_cross, "gap_pre_task": gap_pre_task,
                "gap_task_post": gap_task_post,
                "d_obs_pre_post": d_obs_pre_post,
                "d_obs_pre_task": d_obs_pre_task,
                "d_obs_task_post": d_obs_task_post,
            }

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "drift_triangle.csv", index=False)
    if diagnostic_packs:
        _render_drift_pdf(patient, diagnostic_packs, out_dir / "drift_triangle.pdf")
    return df


def _render_drift_pdf(
    patient: str,
    packs: Dict[str, Dict],
    out_path: Path,
) -> None:
    bands = [b for b in BRAIN_BANDS_NAMES
             if any(k.startswith(f"{b}_") for k in packs)]
    if not bands:
        return
    nb = len(bands)
    nrows = len(DISTANCES)
    with PdfPages(out_path) as pdf:
        fig, axes = plt.subplots(
            nrows, nb,
            figsize=(2.4 * nb + 0.6, 2.0 * nrows + 0.4),
            squeeze=False,
        )
        for c, band in enumerate(bands):
            for r, d_label in enumerate(DISTANCES):
                ax = axes[r, c]
                pk = packs.get(f"{band}_{d_label}")
                if pk is None:
                    ax.set_xticks([]); ax.set_yticks([]); continue
                xv = pk["xv"]; yv = pk["yv"]; tag = pk["phase_tag"]
                # Drift points
                m_pre = tag == "rest_pre"
                m_post = tag == "rest_post"
                ax.scatter(xv[m_pre], yv[m_pre], s=14, c="tab:gray",
                           label="pre", edgecolor="k", linewidth=0.3, alpha=0.7)
                ax.scatter(xv[m_post], yv[m_post], s=14, c="tab:cyan",
                           label="post", edgecolor="k", linewidth=0.3, alpha=0.7)
                # Drift fit line — extend to gap_cross
                xline = np.linspace(0, max(xv.max(), pk["gap_cross"]), 100)
                yline = pk["intercept"] + pk["slope"] * xline
                ax.plot(xline, yline, color="k", lw=0.8, alpha=0.6,
                        label=f"R²={pk['R2']:.2f}")
                # Observed cross-phase distances at the corresponding gaps
                ax.scatter(pk["gap_pre_task"], pk["d_obs_pre_task"],
                           marker="X", s=46, c="tab:red", edgecolor="k",
                           linewidth=0.4, label="pre–task obs", zorder=5)
                ax.scatter(pk["gap_cross"], pk["d_obs_pre_post"],
                           marker="X", s=46, c="tab:blue", edgecolor="k",
                           linewidth=0.4, label="pre–post obs", zorder=5)
                ax.scatter(pk["gap_task_post"], pk["d_obs_task_post"],
                           marker="X", s=46, c="tab:green", edgecolor="k",
                           linewidth=0.4, label="task–post obs", zorder=5)
                ax.tick_params(labelsize=6)
                if r == 0:
                    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
                if c == 0:
                    ax.set_ylabel(f"d_{d_label}", fontsize=9)
                if r == nrows - 1:
                    ax.set_xlabel("gap (s)", fontsize=8)
                if r == 0 and c == nb - 1:
                    ax.legend(fontsize=5, loc="upper left")
        fig.text(0.99, 0.01, f"audit_28 {patient}", ha="right", fontsize=6, color="gray")
        fig.tight_layout()
        pdf.savefig(fig, dpi=200); plt.close(fig)
    print(f"[{patient}] wrote {out_path}")


# ---------------------------------------------------------------------------
# Cohort summary
# ---------------------------------------------------------------------------

def _cohort_summary(per_patient: Dict[str, pd.DataFrame], out_path: Path) -> None:
    pooled = pd.concat(per_patient.values(), ignore_index=True)
    in_pool = pooled  # n=10; Pat_03 included

    rows = []
    for band in BRAIN_BANDS_NAMES:
        for d_label in DISTANCES:
            sub = in_pool[(in_pool.band == band) & (in_pool.distance == d_label)]
            if sub.empty:
                continue
            n_obs_lt_drift_pre_post = int((sub.residual_pre_post < 0).sum())
            n_obs_lt_drift_pre_task = int((sub.residual_pre_task < 0).sum())
            n_obs_gt_drift_pre_task = int((sub.residual_pre_task > 0).sum())
            rows.append({
                "band": band,
                "distance": d_label,
                "n_patients": int(sub.shape[0]),
                "median_R2": round(float(sub.drift_R2.median()), 3),
                "median_d_obs_pre_post": round(float(sub.d_obs_pre_post.median()), 4),
                "median_d_drift_pre_post": round(float(sub.d_drift_pred_pre_post.median()), 4),
                "median_residual_pre_post": round(float(sub.residual_pre_post.median()), 4),
                "median_d_obs_pre_task": round(float(sub.d_obs_pre_task.median()), 4),
                "median_d_drift_pre_task": round(float(sub.d_drift_pred_pre_task.median()), 4),
                "median_residual_pre_task": round(float(sub.residual_pre_task.median()), 4),
                "median_T_d_obs": round(float(sub.T_d_obs.median()), 4),
                "median_T_d_drift": round(float(sub.T_d_drift.median()), 4),
                "n_obs_below_drift_pre_post": n_obs_lt_drift_pre_post,
                "n_obs_above_drift_pre_task": n_obs_gt_drift_pre_task,
            })
    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    print(f"wrote {out_path}")
    n_cohort = len(per_patient)
    print(f"\n=== Cohort drift-triangle summary (n={n_cohort} patients) ===")
    print(df.to_string(index=False))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", default="Pat_06")
    ap.add_argument("--K", type=int, default=4, help="chunks per phase (≥3)")
    ap.add_argument("--cohort", action="store_true")
    args = ap.parse_args()
    OUT_BASE.mkdir(parents=True, exist_ok=True)

    targets = list(PATIENTS_4PHASE) if args.cohort else [args.patient]
    per_patient: Dict[str, pd.DataFrame] = {}
    for p in targets:
        out_dir = OUT_BASE / p
        df = run_patient(p, args.K, out_dir)
        if not df.empty:
            per_patient[p] = df
            if not args.cohort:
                print(f"\n=== {p} drift-triangle table ===")
                show = df[[
                    "band", "distance", "drift_R2", "drift_slope_per_s",
                    "d_obs_pre_post", "d_drift_pred_pre_post", "residual_pre_post",
                    "d_obs_pre_task", "d_drift_pred_pre_task", "residual_pre_task",
                    "T_d_obs", "T_d_drift",
                ]].copy()
                for col in show.columns:
                    if show[col].dtype.kind == "f":
                        show[col] = show[col].round(4)
                print(show.to_string(index=False))

    if args.cohort and per_patient:
        _cohort_summary(per_patient, OUT_BASE / "cohort_drift_triangle.csv")

    return 0


if __name__ == "__main__":
    sys.exit(main())
