#!/usr/bin/env python3
"""Audit Step 4 — per-patient scale-resolved H2a task-trace matrix.

Recomputes VI(k) from the rebuilt linkage matrices for all 6 phase pairs,
then computes the canonical repo H2a per (patient, band, k):

    H2a(p, b, k) = VI(rest_pre, rest_post)(p, b, k)
                 − VI(task_test, rest_post)(p, b, k)

Positive = trace detected.

Outputs:
  data/audit/task_trace_per_patient/Pat_XX.png   (6 bands x K_max heatmap)
  data/audit/task_trace_per_patient/task_trace_arrays.npz
  data/audit/vi_raw_profiles_audit.csv           (reference)
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE, PHASE_LABELS
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE
from lrg_eegfc.utils.metrics.vi import compute_vi


OUT_DIR = ROOT / "data" / "audit" / "task_trace_per_patient"
VI_CSV = ROOT / "data" / "audit" / "vi_raw_profiles_audit.csv"
NPZ = OUT_DIR / "task_trace_arrays.npz"

BAND_ORDER = list(BRAIN_BANDS.keys())  # 6 bands
PHASE_ORDER = list(PHASE_LABELS)       # 4 phases
ALL_PAIRS = [
    ("rest_pre", "task_learn"),
    ("rest_pre", "task_test"),
    ("rest_pre", "rest_post"),
    ("task_learn", "task_test"),
    ("task_learn", "rest_post"),
    ("task_test", "rest_post"),
]


def _load_linkage(patient: str, phase: str, band: str):
    p = IMCOH_LRG_CACHE / patient / f"{band}_{phase}_lrg_imcoh-abs.npz"
    if not p.exists():
        return None, None
    d = np.load(p, allow_pickle=True)
    return np.asarray(d["linkage_matrix"]), int(d["n_nodes"])


def compute_all_vi() -> list[dict]:
    rows: list[dict] = []
    for pat in PATIENTS_4PHASE:
        phase_data = {}
        for phase in PHASE_ORDER:
            for band in BAND_ORDER:
                Z, n = _load_linkage(pat, phase, band)
                if Z is not None:
                    phase_data[(phase, band)] = (Z, n)
        if not phase_data:
            continue
        # per-patient, per-band k-range from min n across phases
        for band in BAND_ORDER:
            ns = [phase_data[(ph, band)][1] for ph in PHASE_ORDER
                  if (ph, band) in phase_data]
            if not ns:
                continue
            n_min = min(ns)
            k_range = list(range(2, max(n_min - 1, 3) + 1))
            for pa, pb in ALL_PAIRS:
                if (pa, band) not in phase_data or (pb, band) not in phase_data:
                    continue
                Za, na = phase_data[(pa, band)]
                Zb, nb = phase_data[(pb, band)]
                if na != nb:
                    continue
                for k in k_range:
                    la = fcluster(Za, k, criterion="maxclust")
                    lb = fcluster(Zb, k, criterion="maxclust")
                    rows.append({
                        "patient": pat, "band": band,
                        "phase_a": pa, "phase_b": pb,
                        "k": k, "vi": compute_vi(la, lb),
                    })
        print(f"  {pat}: {len(rows)} VI rows so far")
    return rows


def build_h2a_array(vi_rows: list[dict]) -> tuple[np.ndarray, list[int]]:
    """Returns H2a[patient_idx, band_idx, k_idx], and the k_values list."""
    # Organize by (pat, band) -> dict of pair -> k -> vi
    vi_lookup: dict = {}
    for r in vi_rows:
        key = (r["patient"], r["band"])
        vi_lookup.setdefault(key, {})
        # use frozenset so order-agnostic
        pair = frozenset({r["phase_a"], r["phase_b"]})
        vi_lookup[key].setdefault(pair, {})
        vi_lookup[key][pair][int(r["k"])] = float(r["vi"])

    # Determine max k across all (patient, band)
    all_k: set[int] = set()
    for d in vi_lookup.values():
        for pair_d in d.values():
            all_k.update(pair_d.keys())
    k_values = sorted(all_k)
    k_to_idx = {k: i for i, k in enumerate(k_values)}

    n_pat = len(PATIENTS_4PHASE)
    n_band = len(BAND_ORDER)
    n_k = len(k_values)
    H2a = np.full((n_pat, n_band, n_k), np.nan)

    pair_rr = frozenset({"rest_pre", "rest_post"})
    pair_tp = frozenset({"task_test", "rest_post"})

    for ip, pat in enumerate(PATIENTS_4PHASE):
        for ib, band in enumerate(BAND_ORDER):
            d = vi_lookup.get((pat, band))
            if d is None:
                continue
            vi_rr = d.get(pair_rr, {})
            vi_tp = d.get(pair_tp, {})
            for k, v_rr in vi_rr.items():
                v_tp = vi_tp.get(k)
                if v_tp is None:
                    continue
                H2a[ip, ib, k_to_idx[k]] = v_rr - v_tp
    return H2a, k_values


def figure_for_patient(pat_idx: int, patient: str, H2a: np.ndarray,
                         k_values: list[int]) -> Path | None:
    panel = H2a[pat_idx]  # (n_band, n_k)
    if np.all(np.isnan(panel)):
        return None
    # clip to the valid k-range: trim trailing NaNs per band; use common extent
    finite_cols = np.isfinite(panel).any(axis=0)
    k_mask = np.where(finite_cols)[0]
    if k_mask.size == 0:
        return None
    k_lo, k_hi = int(k_mask.min()), int(k_mask.max())
    panel_t = panel[:, k_lo:k_hi + 1]
    k_vals = k_values[k_lo:k_hi + 1]

    vmax = np.nanmax(np.abs(panel_t)) if np.any(np.isfinite(panel_t)) else 1.0
    if vmax <= 0:
        vmax = 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(max(7, len(k_vals) * 0.06), 3.8))
    im = ax.imshow(panel_t, aspect="auto", cmap="RdBu_r", norm=norm,
                   interpolation="nearest",
                   extent=[k_vals[0] - 0.5, k_vals[-1] + 0.5,
                           len(BAND_ORDER) - 0.5, -0.5])
    ax.set_yticks(range(len(BAND_ORDER)))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BAND_ORDER])
    ax.set_xlabel("k (number of clusters)")
    ax.set_ylabel("band")
    ax.set_title(
        f"{patient} — H2a = VI(rsPre,rsPost) − VI(taskTest,rsPost); + = trace"
    )
    fig.colorbar(im, ax=ax, label="H2a (nats)")
    out = OUT_DIR / f"{patient}.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Computing VI(k) for all phase pairs ...")
    vi_rows = compute_all_vi()
    print(f"Total VI rows: {len(vi_rows)}")

    # Write VI CSV
    with VI_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["patient", "band", "phase_a", "phase_b", "k", "vi"]
        )
        writer.writeheader()
        for r in vi_rows:
            writer.writerow({
                "patient": r["patient"], "band": r["band"],
                "phase_a": r["phase_a"], "phase_b": r["phase_b"],
                "k": int(r["k"]), "vi": f"{r['vi']:.6g}",
            })
    print(f"Wrote {VI_CSV.relative_to(ROOT)}")

    # Build H2a arrays
    H2a, k_values = build_h2a_array(vi_rows)
    print(f"H2a shape: {H2a.shape}  k range: [{min(k_values)}, {max(k_values)}]")

    # Per-patient figures
    n_with = 0
    for ip, pat in enumerate(PATIENTS_4PHASE):
        p = figure_for_patient(ip, pat, H2a, k_values)
        if p is not None:
            n_with += 1
            print(f"  wrote {p.relative_to(ROOT)}")
    print(f"{n_with} per-patient heatmaps written")

    # Save arrays
    np.savez_compressed(
        NPZ,
        H2a=H2a,
        patient_ids=np.array(PATIENTS_4PHASE),
        band_names=np.array(BAND_ORDER),
        k_values=np.array(k_values, dtype=int),
    )
    print(f"Wrote {NPZ.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
