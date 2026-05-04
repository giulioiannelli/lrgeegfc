#!/usr/bin/env python3
"""Audit 30 — within-``rest_post`` split-half null (Control 2).

Mirror of ``audit_25``'s within-``rest_pre`` null but applied to
``rest_post``. Used to confirm rsPost is itself stationary at the
50-split scale, ruling out a transient post-task baseline anomaly
inflating ``d(task_test, rest_post)`` for non-trace reasons.

Per (patient, band, distance):
- 50 random equal-half partitions of the rest_post non-overlapping
  segment FFTs.
- For each split, compute |ImCoh| on each half restricted to
  ``V*(p, b) = ⋂_φ GC(|A^φ|)`` (same V* as audit_25 — intersection
  across rest_pre, task_test, rest_post).
- Three distances ``d_P, d_S, d_F`` between halves.

Outputs
-------
``data/audit/raw_fc_phase_distance/<Patient>/null_rsPost.npz``
``data/audit/raw_fc_phase_distance/cohort_null_rsPost.csv`` (cohort)

The CSV columns mirror the spec in
``.agents/reports/2026-04-30_drift-controls-spec.md``:
``patient, band, distance, null_median, null_mad, null_q1, null_q3,
n_splits``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import networkx as nx
import numpy as np
import pandas as pd
from scipy.signal import get_window
from scipy.stats import pearsonr, spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    BRAIN_BANDS_NAMES,
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


def _split_half_null(
    fft_seg: np.ndarray,
    freqs: np.ndarray,
    fmin: float,
    fmax: float,
    common_idx: np.ndarray,
    n_split: int,
    rng: np.random.Generator,
) -> Dict[str, np.ndarray]:
    if fft_seg.size == 0:
        return {d: np.full(n_split, np.nan) for d in DISTANCES}
    n_seg = fft_seg.shape[1]
    half = n_seg // 2
    if half < 4:
        return {d: np.full(n_split, np.nan) for d in DISTANCES}
    freq_mask = (freqs >= fmin) & (freqs < fmax)
    if not freq_mask.any():
        idx = int(np.clip(np.searchsorted(freqs, fmin), 0, len(freqs) - 1))
        freq_mask = np.zeros_like(freqs, dtype=bool); freq_mask[idx] = True
    null = {d: np.empty(n_split) for d in DISTANCES}
    for i in range(n_split):
        order = rng.permutation(n_seg)
        a_idx = order[:half]
        b_idx = order[half : 2 * half]
        A_a = _imcoh_abs_from_subset(fft_seg, a_idx, freq_mask)
        A_b = _imcoh_abs_from_subset(fft_seg, b_idx, freq_mask)
        A_a = A_a[np.ix_(common_idx, common_idx)]
        A_b = A_b[np.ix_(common_idx, common_idx)]
        for d_label, fn in _DIST_FUNCS.items():
            null[d_label][i] = fn(A_a, A_b)
    return null


def _v_star_per_band(patient: str) -> Dict[str, np.ndarray]:
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
        idx = _common_giant_indices(adj_phases)
        if idx.size >= 10:
            v_star[band] = idx
    return v_star


def run_patient(patient: str, n_split: int, seed: int, out_dir: Path) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    fs = FS_OVERRIDES.get(patient, DEFAULT_SAMPLE_RATE)
    nperseg = nperseg_for_fs(fs)
    print(f"[{patient}] fs={fs:.1f} Hz, nperseg={nperseg}")

    v_star = _v_star_per_band(patient)
    if not v_star:
        print(f"[{patient}] no eligible bands")
        return pd.DataFrame()

    try:
        ts = load_timeseries(patient, "rest_post", SEEG_DATAPATH)
    except Exception as exc:
        print(f"[{patient}] FAILED to load rest_post: {exc}")
        return pd.DataFrame()
    freqs, fft_seg = _segment_ffts(ts, fs, nperseg)
    print(f"[{patient}/rest_post] shape={ts.shape}, n_seg={fft_seg.shape[1]}")

    rows: List[Dict] = []
    null_store: Dict[str, np.ndarray] = {}
    for band in BRAIN_BANDS_NAMES:
        if band not in v_star:
            continue
        fmin, fmax = BRAIN_BANDS[band]
        null = _split_half_null(
            fft_seg, freqs, fmin, fmax, v_star[band], n_split, rng
        )
        for d_label in DISTANCES:
            arr = null[d_label]
            null_store[f"{band}_{d_label}"] = arr
            arr_finite = arr[np.isfinite(arr)]
            if arr_finite.size < 3:
                med = float("nan"); mad = float("nan")
                q1 = float("nan"); q3 = float("nan")
            else:
                med = float(np.median(arr_finite))
                mad = float(np.median(np.abs(arr_finite - med)))
                q1 = float(np.quantile(arr_finite, 0.25))
                q3 = float(np.quantile(arr_finite, 0.75))
            rows.append({
                "patient": patient,
                "band": band,
                "distance": d_label,
                "null_median": med,
                "null_mad": mad,
                "null_q1": q1,
                "null_q3": q3,
                "n_splits": int(arr_finite.size),
                "n_common": int(v_star[band].size),
            })

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "null_rsPost.csv", index=False)
    np.savez_compressed(out_dir / "null_rsPost.npz", **null_store)
    print(f"[{patient}] wrote null_rsPost.{{csv,npz}}")
    return df


def _load_rsPre_null_summary(patient: str, out_base: Path) -> pd.DataFrame:
    """Reload audit_25 null.npz to get per-cell rsPre null median/mad/IQR."""
    p_dir = out_base / patient
    null_path = p_dir / "null.npz"
    if not null_path.exists():
        return pd.DataFrame()
    rows: List[Dict] = []
    with np.load(null_path) as data:
        for key in data.files:
            band, d_label = key.rsplit("_", 1)
            if d_label not in DISTANCES:
                continue
            arr = np.asarray(data[key])
            arr_finite = arr[np.isfinite(arr)]
            if arr_finite.size < 3:
                med = float("nan"); mad = float("nan")
                q1 = float("nan"); q3 = float("nan")
            else:
                med = float(np.median(arr_finite))
                mad = float(np.median(np.abs(arr_finite - med)))
                q1 = float(np.quantile(arr_finite, 0.25))
                q3 = float(np.quantile(arr_finite, 0.75))
            rows.append({
                "patient": patient, "band": band, "distance": d_label,
                "rsPre_null_median": med, "rsPre_null_mad": mad,
                "rsPre_null_q1": q1, "rsPre_null_q3": q3,
            })
    return pd.DataFrame(rows)


def cohort_summary(per_patient: Dict[str, pd.DataFrame], out_base: Path) -> None:
    pooled = pd.concat(per_patient.values(), ignore_index=True)
    pooled = pooled.rename(columns={
        "null_median": "rsPost_null_median",
        "null_mad": "rsPost_null_mad",
        "null_q1": "rsPost_null_q1",
        "null_q3": "rsPost_null_q3",
        "n_splits": "rsPost_n_splits",
    })

    pre_rows = [_load_rsPre_null_summary(p, out_base) for p in per_patient]
    pre_df = pd.concat([r for r in pre_rows if not r.empty], ignore_index=True)
    merged = pooled.merge(pre_df, on=["patient", "band", "distance"], how="left")
    merged["mad_ratio_post_over_pre"] = merged["rsPost_null_mad"] / merged["rsPre_null_mad"]
    merged["iqr_post"] = merged["rsPost_null_q3"] - merged["rsPost_null_q1"]
    merged["iqr_pre"] = merged["rsPre_null_q3"] - merged["rsPre_null_q1"]
    merged["iqr_ratio_post_over_pre"] = merged["iqr_post"] / merged["iqr_pre"]

    out_csv = out_base / "cohort_null_rsPost.csv"
    merged.to_csv(out_csv, index=False)
    print(f"wrote {out_csv}")

    # Per-(band, distance) cohort summary
    summary_rows = []
    for band in BRAIN_BANDS_NAMES:
        for d_label in DISTANCES:
            sub = merged[(merged.band == band) & (merged.distance == d_label)]
            if sub.empty:
                continue
            summary_rows.append({
                "band": band,
                "distance": d_label,
                "n_patients": int(sub.shape[0]),
                "median_mad_ratio": round(float(sub.mad_ratio_post_over_pre.median()), 3),
                "median_iqr_ratio": round(float(sub.iqr_ratio_post_over_pre.median()), 3),
                "n_post_wider_2x": int((sub.mad_ratio_post_over_pre > 2.0).sum()),
                "n_pre_wider_2x": int((sub.mad_ratio_post_over_pre < 0.5).sum()),
                "median_rsPost_mad": round(float(sub.rsPost_null_mad.median()), 5),
                "median_rsPre_mad": round(float(sub.rsPre_null_mad.median()), 5),
            })
    band_summary = pd.DataFrame(summary_rows)
    band_summary_csv = out_base / "cohort_null_rsPost_band_summary.csv"
    band_summary.to_csv(band_summary_csv, index=False)
    print(f"wrote {band_summary_csv}")
    print(f"\n=== Cohort rsPost null vs rsPre null (n={len(per_patient)} patients) ===")
    print(band_summary.to_string(index=False))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", default=None)
    ap.add_argument("--cohort", action="store_true")
    ap.add_argument("--n-split", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    OUT_BASE.mkdir(parents=True, exist_ok=True)

    targets = list(PATIENTS_4PHASE) if args.cohort else [args.patient or PATIENTS_4PHASE[0]]
    per_patient: Dict[str, pd.DataFrame] = {}
    for p in targets:
        out_dir = OUT_BASE / p
        df = run_patient(p, args.n_split, args.seed, out_dir)
        if not df.empty:
            per_patient[p] = df

    if args.cohort and per_patient:
        cohort_summary(per_patient, OUT_BASE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
