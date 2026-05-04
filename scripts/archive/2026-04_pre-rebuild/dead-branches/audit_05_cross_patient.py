#!/usr/bin/env python3
"""Audit Step 5 — cross-patient H2a summary (descriptive; no p-values).

Reads data/audit/task_trace_per_patient/task_trace_arrays.npz
Writes:
  data/audit/task_trace_cross_patient.png  (four 6 x K_max panels)
  data/audit/task_trace_cross_patient.npz  (mean, median, n_pos, sign_frac)

Panels:
  (a) mean H2a per (band, k)
  (b) median H2a per (band, k)
  (c) number of patients with H2a > 0  per (band, k)   integer 0..N
  (d) sign consistency — fraction of patients with the modal sign per (band, k)

NaN-aware: patients with missing H2a at a given (band, k) (notably Pat_14
which lacks task_test) are excluded from that cell's aggregation.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT

IN_NPZ = ROOT / "data" / "audit" / "task_trace_per_patient" / "task_trace_arrays.npz"
OUT_PNG = ROOT / "data" / "audit" / "task_trace_cross_patient.png"
OUT_NPZ = ROOT / "data" / "audit" / "task_trace_cross_patient.npz"

BAND_ORDER = list(BRAIN_BANDS.keys())


def aggregate(H2a: np.ndarray) -> dict:
    """H2a shape (n_patients, n_bands, n_k). Returns per-(band, k) aggregates."""
    n_pat, n_band, n_k = H2a.shape
    finite = np.isfinite(H2a)
    n_contrib = finite.sum(axis=0)  # (n_band, n_k)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.nanmean(H2a, axis=0)
        median = np.nanmedian(H2a, axis=0)
    n_pos = (H2a > 0).sum(axis=0)   # (n_band, n_k)
    n_neg = (H2a < 0).sum(axis=0)
    # modal sign count: max(n_pos, n_neg), then divided by n_contrib
    modal = np.maximum(n_pos, n_neg)
    with np.errstate(invalid="ignore", divide="ignore"):
        sign_frac = np.where(n_contrib > 0, modal / n_contrib, np.nan)
    return {
        "mean": mean,
        "median": median,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_contrib": n_contrib,
        "sign_frac": sign_frac,
    }


def plot(agg: dict, k_values: np.ndarray, n_patients_total: int,
          out_path: Path, title_suffix: str = "") -> None:
    fig, axes = plt.subplots(4, 1, figsize=(12, 12), sharex=True,
                               gridspec_kw={"hspace": 0.28})

    extent = [k_values[0] - 0.5, k_values[-1] + 0.5, len(BAND_ORDER) - 0.5, -0.5]

    def _diverging(mat):
        finite = mat[np.isfinite(mat)]
        if finite.size == 0:
            return TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
        vmax = max(abs(float(finite.min())), abs(float(finite.max())))
        if vmax <= 0:
            vmax = 1.0
        return TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    # panel (a): mean
    ax = axes[0]
    im = ax.imshow(agg["mean"], aspect="auto", cmap="RdBu_r",
                   norm=_diverging(agg["mean"]), interpolation="nearest",
                   extent=extent)
    ax.set_title(f"(a) mean H2a{title_suffix}")
    fig.colorbar(im, ax=ax, label="nats", shrink=0.85)

    # panel (b): median
    ax = axes[1]
    im = ax.imshow(agg["median"], aspect="auto", cmap="RdBu_r",
                   norm=_diverging(agg["median"]), interpolation="nearest",
                   extent=extent)
    ax.set_title(f"(b) median H2a{title_suffix}")
    fig.colorbar(im, ax=ax, label="nats", shrink=0.85)

    # panel (c): n_pos (integer 0..N)
    ax = axes[2]
    im = ax.imshow(agg["n_pos"], aspect="auto", cmap="viridis",
                   vmin=0, vmax=n_patients_total, interpolation="nearest",
                   extent=extent)
    ax.set_title(
        f"(c) # patients with H2a > 0  (of {n_patients_total} contributing)"
    )
    fig.colorbar(im, ax=ax, label="count", shrink=0.85)

    # panel (d): sign consistency (0.5..1.0)
    ax = axes[3]
    im = ax.imshow(agg["sign_frac"], aspect="auto", cmap="magma",
                   vmin=0.5, vmax=1.0, interpolation="nearest",
                   extent=extent)
    ax.set_title("(d) sign consistency (fraction on modal sign)")
    fig.colorbar(im, ax=ax, label="fraction", shrink=0.85)

    for ax in axes:
        ax.set_yticks(range(len(BAND_ORDER)))
        ax.set_yticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BAND_ORDER])
        ax.set_ylabel("band")
    axes[-1].set_xlabel("k (number of clusters)")

    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    d = np.load(IN_NPZ, allow_pickle=True)
    H2a = d["H2a"]  # (n_pat, n_band, n_k)
    pats = list(d["patient_ids"])
    bands = list(d["band_names"])
    k_values = np.asarray(d["k_values"], dtype=int)

    # Step 5 says "across all 9 patients". Pat_14 is all-NaN (no task_test),
    # so nan-aware aggregation naturally excludes it. Keep all 10 in the
    # array and track contributing-patient count in panel (c).
    finite_any = np.isfinite(H2a).any(axis=(1, 2))
    pats_contrib = [p for p, k in zip(pats, finite_any) if k]
    n_contrib_total = len(pats_contrib)

    print(f"Contributing patients: {n_contrib_total}  "
          f"({', '.join(pats_contrib)})")

    agg = aggregate(H2a)
    plot(agg, k_values, n_contrib_total, OUT_PNG,
          title_suffix=f"  (n={n_contrib_total})")

    np.savez_compressed(
        OUT_NPZ,
        mean=agg["mean"], median=agg["median"],
        n_pos=agg["n_pos"], n_neg=agg["n_neg"],
        n_contrib=agg["n_contrib"], sign_frac=agg["sign_frac"],
        band_names=np.array(bands), k_values=k_values,
        patient_ids_contributing=np.array(pats_contrib),
    )
    print(f"Wrote {OUT_PNG.relative_to(ROOT)}")
    print(f"Wrote {OUT_NPZ.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
