#!/usr/bin/env python3
"""Audit 26 — 4-phase FC geometry on raw `imcoh_abs` matrices.

Companion to ``audit_25_raw_fc_phase_distance.py``. Where audit_25 looks
at three phases (`rest_pre`, `task_test`, `rest_post`) against a within-
`rest_pre` null, this audit adds `task_learn` and computes the full 6
pairwise distances on cached observed FC. Cheap because no null is
recomputed — every observed FC matrix is already cached.

The point: contextualize the persistence claim with positive controls.

- `d(task_learn, task_test)` — within-task baseline. Tasks are
  cognitively similar, so their FC should be **close** on rank
  distances. This is the "task cluster width."
- `d(rest_pre, rest_post)` — within-rest "baseline." If rest does NOT
  drift across the session, this should also be small. If it IS large,
  rest_pre and rest_post ARE different — and the question is *whether
  the difference is task-shaped*.
- `d(task_test, rest_post)` vs `d(task_test, rest_pre)` — the
  persistence claim. If `rest_post` sits closer to the task cluster
  than `rest_pre` does, post is task-shaped.

Outputs
-------
``data/audit/fc_phase_geometry/<Patient>/distance_4phase.csv``
    Long-format per-patient table:
    (band, distance, phase_A, phase_B, d_obs, n_common).
``data/audit/fc_phase_geometry/<Patient>/distance_4phase.pdf``
    Per-patient diagnostic: rows = distances (P/S/F), cols = bands,
    each panel a 4×4 heatmap of inter-phase distances on `V*`.
``data/audit/fc_phase_geometry/cohort_geometry_4phase_summary.csv``
    Cohort-pooled table:
    (band, distance, phase_A, phase_B, median, q1, q3, n).
``data/audit/fc_phase_geometry/cohort_geometry.pdf``
    Cohort diagnostic: rows = distances, cols = bands, each panel a
    4×4 median-distance heatmap with patient-spread annotation, and a
    bar plot summarising the within-cluster vs cross-cluster gap.

Usage
-----
    python audit_26_fc_phase_geometry.py --patient Pat_06   # one patient
    python audit_26_fc_phase_geometry.py --cohort           # all 10
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
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PATIENTS_4PHASE,
)
from lrg_eegfc.config.paths import DATA_ROOT
from lrg_eegfc.visuals.correlation import imshow_colorbar_caxdivider
from lrg_eegfc.workflow.fc import load_fc_matrix

PHASES_4: Tuple[str, ...] = ("rest_pre", "task_learn", "task_test", "rest_post")
PHASE_SHORT = {"rest_pre": "RPre", "task_learn": "TL",
               "task_test": "TT", "rest_post": "RPost"}
DISTANCES: Tuple[str, ...] = ("P", "S", "F")
OUT_BASE = DATA_ROOT / "audit" / "raw_fc_phase_distance"


# ---------------------------------------------------------------------------
# Distance functions (mirrors audit_25)
# ---------------------------------------------------------------------------


def _triu(A: np.ndarray) -> np.ndarray:
    return A[np.triu_indices_from(A, k=1)]


def d_P(A: np.ndarray, B: np.ndarray) -> float:
    a, b = _triu(A), _triu(B)
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    r, _ = pearsonr(a, b)
    return float(1.0 - r)


def d_S(A: np.ndarray, B: np.ndarray) -> float:
    a, b = _triu(A), _triu(B)
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(1.0 - spearmanr(a, b).statistic)


def d_F(A: np.ndarray, B: np.ndarray) -> float:
    nA = np.linalg.norm(A, ord="fro")
    nB = np.linalg.norm(B, ord="fro")
    if nA == 0 or nB == 0:
        return float("nan")
    return float(np.linalg.norm(A - B, ord="fro") / np.sqrt(nA * nB))


def all_distances(A: np.ndarray, B: np.ndarray) -> Dict[str, float]:
    return {"P": d_P(A, B), "S": d_S(A, B), "F": d_F(A, B)}


# ---------------------------------------------------------------------------
# Substrate construction (V* across all 4 phases)
# ---------------------------------------------------------------------------


def _common_giant_indices(adj_phases: Dict[str, np.ndarray]) -> np.ndarray:
    sets: List[set] = []
    for adj in adj_phases.values():
        G = nx.from_numpy_array(np.abs(adj))
        comps = list(nx.connected_components(G))
        if not comps:
            sets.append(set()); continue
        sets.append(set(max(comps, key=len)))
    if not sets:
        return np.array([], dtype=int)
    common = set.intersection(*sets)
    return np.array(sorted(common), dtype=int)


# ---------------------------------------------------------------------------
# Per-patient pipeline
# ---------------------------------------------------------------------------


def run_patient(patient: str, out_dir: Path) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[{patient}] computing 4-phase distances")

    rows: List[Dict] = []
    fc_store: Dict[Tuple[str, str], np.ndarray] = {}

    for band in BRAIN_BANDS_NAMES:
        adj_phases: Dict[str, np.ndarray] = {}
        missing = False
        for phase in PHASES_4:
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
            print(f"[{patient}/{band}] |V*|={common_idx.size} < 10 — skipping band")
            continue

        A_phase = {
            phase: adj_phases[phase][np.ix_(common_idx, common_idx)]
            for phase in PHASES_4
        }
        for phase in PHASES_4:
            fc_store[(band, phase)] = A_phase[phase]

        # All ordered pairs — store as upper-triangle of the 4×4
        for i, phi_A in enumerate(PHASES_4):
            for j, phi_B in enumerate(PHASES_4):
                if j <= i:
                    continue
                d = all_distances(A_phase[phi_A], A_phase[phi_B])
                for d_label in DISTANCES:
                    rows.append({
                        "patient": patient,
                        "band": band,
                        "distance": d_label,
                        "phase_A": phi_A,
                        "phase_B": phi_B,
                        "d_obs": d[d_label],
                        "n_common": int(common_idx.size),
                    })

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "distance_4phase.csv", index=False)
    if not df.empty:
        _render_per_patient(patient, df, out_dir)
    return df


# ---------------------------------------------------------------------------
# Per-patient PDF: 4×4 heatmap per band per distance
# ---------------------------------------------------------------------------


def _build_4x4(df_band_dist: pd.DataFrame) -> np.ndarray:
    """Build a 4×4 symmetric distance matrix for one (patient, band, distance)."""
    M = np.full((4, 4), np.nan)
    np.fill_diagonal(M, 0.0)
    for _, r in df_band_dist.iterrows():
        i = PHASES_4.index(r.phase_A)
        j = PHASES_4.index(r.phase_B)
        M[i, j] = r.d_obs
        M[j, i] = r.d_obs
    return M


def _render_per_patient(patient: str, df: pd.DataFrame, out_dir: Path) -> None:
    pdf_path = out_dir / "distance_4phase.pdf"
    bands = [b for b in BRAIN_BANDS_NAMES if b in df.band.unique()]
    if not bands:
        return
    nb = len(bands)
    nrows = len(DISTANCES)
    short_labels = [PHASE_SHORT[p] for p in PHASES_4]

    with PdfPages(pdf_path) as pdf:
        fig, axes = plt.subplots(
            nrows, nb,
            figsize=(2.0 * nb + 1.2, 1.9 * nrows + 0.4),
            squeeze=False,
        )
        for r, d_label in enumerate(DISTANCES):
            row_max = 0.0
            mats = {}
            for band in bands:
                sub = df[(df.band == band) & (df.distance == d_label)]
                M = _build_4x4(sub)
                mats[band] = M
                if np.isfinite(M).any():
                    row_max = max(row_max, float(np.nanmax(M)))
            row_max = row_max or 1.0
            last_im = None
            for c, band in enumerate(bands):
                ax = axes[r, c]
                im = ax.imshow(mats[band], cmap="viridis",
                               vmin=0, vmax=row_max, aspect="equal")
                last_im = im
                ax.set_xticks(range(4)); ax.set_yticks(range(4))
                ax.set_xticklabels(short_labels, fontsize=6, rotation=45)
                ax.set_yticklabels(short_labels, fontsize=6)
                # Annotate cells with values
                for i in range(4):
                    for j in range(4):
                        v = mats[band][i, j]
                        if np.isnan(v):
                            continue
                        ax.text(j, i, f"{v:.2f}",
                                ha="center", va="center",
                                fontsize=5,
                                color="white" if v > row_max * 0.5 else "black")
                if r == 0:
                    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
                if c == 0:
                    ax.set_ylabel(f"d_{d_label}", fontsize=9)
            if last_im is not None:
                imshow_colorbar_caxdivider(
                    last_im, axes[r, nb - 1], size="4%", pad=0.06,
                )
        fig.text(0.99, 0.01, f"audit_26 P1 {patient}",
                 ha="right", fontsize=6, color="gray")
        fig.tight_layout()
        pdf.savefig(fig, dpi=200); plt.close(fig)
    print(f"[{patient}] wrote {pdf_path}")


# ---------------------------------------------------------------------------
# Cohort summary + diagnostic
# ---------------------------------------------------------------------------


def _cohort_summary(per_patient: Dict[str, pd.DataFrame], out_path: Path) -> None:
    pooled = pd.concat(per_patient.values(), ignore_index=True)
    rows: List[Dict] = []
    for band in BRAIN_BANDS_NAMES:
        for d_label in DISTANCES:
            for i, phi_A in enumerate(PHASES_4):
                for j, phi_B in enumerate(PHASES_4):
                    if j <= i:
                        continue
                    sub = pooled[
                        (pooled.band == band)
                        & (pooled.distance == d_label)
                        & (pooled.phase_A == phi_A)
                        & (pooled.phase_B == phi_B)
                    ]
                    vals = sub.d_obs.dropna().values
                    if vals.size == 0:
                        continue
                    rows.append({
                        "band": band,
                        "distance": d_label,
                        "phase_A": phi_A,
                        "phase_B": phi_B,
                        "median": float(np.median(vals)),
                        "q1": float(np.percentile(vals, 25)),
                        "q3": float(np.percentile(vals, 75)),
                        "n": int(vals.size),
                    })
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"wrote {out_path}")


def _render_cohort(per_patient: Dict[str, pd.DataFrame],
                   summary: pd.DataFrame, out_path: Path) -> None:
    pooled = pd.concat(per_patient.values(), ignore_index=True)
    bands = [b for b in BRAIN_BANDS_NAMES]
    nb = len(bands)
    nrows = len(DISTANCES)
    short_labels = [PHASE_SHORT[p] for p in PHASES_4]

    with PdfPages(out_path) as pdf:
        # ---- Page 1: cohort-median 4×4 heatmaps ----
        fig, axes = plt.subplots(
            nrows, nb,
            figsize=(2.0 * nb + 1.2, 1.9 * nrows + 0.4),
            squeeze=False,
        )
        for r, d_label in enumerate(DISTANCES):
            mats = {}
            row_max = 0.0
            for band in bands:
                M = np.full((4, 4), np.nan); np.fill_diagonal(M, 0.0)
                for _, row in summary[
                    (summary.band == band) & (summary.distance == d_label)
                ].iterrows():
                    i = PHASES_4.index(row.phase_A)
                    j = PHASES_4.index(row.phase_B)
                    M[i, j] = row["median"]
                    M[j, i] = row["median"]
                mats[band] = M
                if np.isfinite(M).any():
                    row_max = max(row_max, float(np.nanmax(M)))
            row_max = row_max or 1.0
            last_im = None
            for c, band in enumerate(bands):
                ax = axes[r, c]
                im = ax.imshow(mats[band], cmap="viridis",
                               vmin=0, vmax=row_max, aspect="equal")
                last_im = im
                ax.set_xticks(range(4)); ax.set_yticks(range(4))
                ax.set_xticklabels(short_labels, fontsize=6, rotation=45)
                ax.set_yticklabels(short_labels, fontsize=6)
                for i in range(4):
                    for j in range(4):
                        v = mats[band][i, j]
                        if np.isnan(v): continue
                        ax.text(j, i, f"{v:.2f}",
                                ha="center", va="center", fontsize=5,
                                color="white" if v > row_max * 0.5 else "black")
                if r == 0:
                    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
                if c == 0:
                    ax.set_ylabel(f"d_{d_label} (median)", fontsize=9)
            if last_im is not None:
                imshow_colorbar_caxdivider(
                    last_im, axes[r, nb - 1], size="4%", pad=0.06,
                )
        fig.text(0.99, 0.01, "audit_26 cohort medians",
                 ha="right", fontsize=6, color="gray")
        fig.tight_layout()
        pdf.savefig(fig, dpi=200); plt.close(fig)

        # ---- Page 2: within-task vs within-rest vs persistence comparison ----
        # Per (band, distance) plot a tiny bar chart with 4 bars:
        #   d(TL, TT)        — within-task
        #   d(RPre, RPost)   — within-rest (drift floor)
        #   d(RPre, TT)      — pre→task cross-cluster
        #   d(TT, RPost)     — task→post cross-cluster
        bar_pairs = [
            ("task_learn", "task_test", "TL↔TT", "tab:purple"),
            ("rest_pre", "rest_post", "RPre↔RPost", "tab:gray"),
            ("rest_pre", "task_test", "RPre↔TT", "tab:red"),
            ("task_test", "rest_post", "TT↔RPost", "tab:green"),
        ]
        fig, axes = plt.subplots(
            nrows, nb,
            figsize=(2.0 * nb + 1.0, 1.9 * nrows + 0.4),
            squeeze=False, sharey="row",
        )
        for r, d_label in enumerate(DISTANCES):
            for c, band in enumerate(bands):
                ax = axes[r, c]
                xs = np.arange(len(bar_pairs))
                vals_med = []
                vals_q1 = []
                vals_q3 = []
                colors = []
                labels = []
                for phi_A, phi_B, lab, col in bar_pairs:
                    sub = pooled[
                        (pooled.band == band)
                        & (pooled.distance == d_label)
                        & (
                            ((pooled.phase_A == phi_A) & (pooled.phase_B == phi_B))
                            | ((pooled.phase_A == phi_B) & (pooled.phase_B == phi_A))
                        )
                    ]
                    v = sub.d_obs.dropna().values
                    if v.size == 0:
                        vals_med.append(0); vals_q1.append(0); vals_q3.append(0)
                    else:
                        vals_med.append(float(np.median(v)))
                        vals_q1.append(float(np.percentile(v, 25)))
                        vals_q3.append(float(np.percentile(v, 75)))
                    colors.append(col); labels.append(lab)
                vals_med = np.asarray(vals_med)
                yerr_lower = vals_med - np.asarray(vals_q1)
                yerr_upper = np.asarray(vals_q3) - vals_med
                ax.bar(xs, vals_med, yerr=[yerr_lower, yerr_upper],
                       color=colors, edgecolor="k", linewidth=0.4,
                       capsize=2, alpha=0.85)
                ax.set_xticks(xs)
                ax.set_xticklabels(labels, fontsize=5, rotation=45, ha="right")
                ax.tick_params(axis="y", labelsize=6)
                if r == 0:
                    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
                if c == 0:
                    ax.set_ylabel(f"d_{d_label}\nmedian (Q1, Q3)",
                                  fontsize=8)
        fig.text(0.99, 0.01, "audit_26 cohort cluster bars",
                 ha="right", fontsize=6, color="gray")
        fig.tight_layout()
        pdf.savefig(fig, dpi=200); plt.close(fig)
    print(f"wrote {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", default="Pat_06")
    ap.add_argument("--cohort", action="store_true")
    args = ap.parse_args()

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    targets = list(PATIENTS_4PHASE) if args.cohort else [args.patient]

    per_patient: Dict[str, pd.DataFrame] = {}
    for p in targets:
        out_dir = OUT_BASE / p
        df = run_patient(p, out_dir)
        if not df.empty:
            per_patient[p] = df

    if args.cohort and per_patient:
        _cohort_summary(per_patient, OUT_BASE / "cohort_geometry_4phase_summary.csv")
        # Also load summary (just produced) and render cohort PDF
        summary = pd.read_csv(OUT_BASE / "cohort_geometry_4phase_summary.csv")
        _render_cohort(per_patient, summary, OUT_BASE / "cohort_geometry_4phase.pdf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
