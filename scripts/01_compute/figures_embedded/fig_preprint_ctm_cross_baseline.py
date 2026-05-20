"""CTM cohort fingerprint — cross-baseline scatter (port of audit_31 idiom).

For each band, plot per-patient (rho_drift, rho_split) as a 2D scatter:
  x = rho_null_drift  (within-session drift estimate)
  y = rho_split       (controlled real CTM correlation)
  identity (y = x) line, trace zone = above identity (real beats drift)
  cohort 1σ ellipse + green mean star
  per-patient circles (Pat_03 = orange triangle)
  n_trace annotation (# patients above identity)

1×6 panels (one per band), shared layout. β panel should show the cohort
cluster sitting above the identity line; θ / γ_l straddle.

Data source: data/audit/ctm_triangle/Td_per_patient_per_band.csv
Output:      data/reports/preprint/cohort_measures/ctm_cross_baseline.pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style

use_lrg_style()


ROOT = Path(__file__).resolve().parents[3]
SRC_CSV = ROOT / "data/audit/ctm_triangle/Td_per_patient_per_band.csv"
OUT_DIR = ROOT / "data/reports/preprint/cohort_measures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_PDF = OUT_DIR / "ctm_cross_baseline.pdf"

X_COL = "rho_null_drift"
Y_COL = "rho_split"


def _ellipse_from_xy(x: np.ndarray, y: np.ndarray) -> Ellipse | None:
    if x.size < 2:
        return None
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals = vals[order]; vecs = vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    width, height = 2.0 * np.sqrt(np.maximum(vals, 0.0))
    if width <= 0 or height <= 0:
        return None
    return Ellipse(
        (x.mean(), y.mean()), width, height, angle=angle,
        facecolor="#888", alpha=0.18, edgecolor="#444", lw=0.7, zorder=2,
    )


def render(df: pd.DataFrame, out_pdf: Path) -> None:
    bands = list(BRAIN_BANDS_NAMES)
    nb = len(bands)
    fig, axes = plt.subplots(1, nb, figsize=(2.55 * nb, 3.0), squeeze=False)
    for c, band in enumerate(bands):
        ax = axes[0, c]
        sub = df[df.band == band].copy()
        if sub.empty:
            ax.set_title(BRAIN_BAND_TEX_DICT[band])
            ax.set_xticks([]); ax.set_yticks([])
            continue
        pat = sub.patient.values
        x = sub[X_COL].astype(float).values
        y = sub[Y_COL].astype(float).values
        keep = np.isfinite(x) & np.isfinite(y)
        pat, x, y = pat[keep], x[keep], y[keep]

        span = max(x.max(), y.max()) - min(x.min(), y.min())
        lo = min(x.min(), y.min()) - 0.10 * span
        hi = max(x.max(), y.max()) + 0.10 * span

        # Trace zone (above identity, y > x) — light red shade
        ax.fill([lo, hi, lo], [lo, hi, hi], color="#d62728", alpha=0.06, zorder=0)
        # Zero reference (rho = 0 horizontal + vertical)
        ax.axhline(0, color="0.7", lw=0.4, ls=":", zorder=1)
        ax.axvline(0, color="0.7", lw=0.4, ls=":", zorder=1)
        # Identity
        ax.plot([lo, hi], [lo, hi], color="k", lw=0.7, ls="--", alpha=0.55, zorder=1)

        ell = _ellipse_from_xy(x, y)
        if ell is not None:
            ax.add_patch(ell)
        ax.scatter([x.mean()], [y.mean()], marker="*", s=180, c="#2ca02c",
                   edgecolor="k", linewidth=0.9, zorder=5)

        is_p03 = pat == "Pat_03"
        ax.scatter(x[~is_p03], y[~is_p03], c="white",
                   edgecolor="#1f77b4", linewidth=1.1, s=34, zorder=3)
        if is_p03.any():
            ax.scatter(x[is_p03], y[is_p03], marker="^", c="white",
                       edgecolor="#ff7f0e", linewidth=1.2, s=46, zorder=4)

        n_trace = int((y > x).sum())
        n_total = int(x.size)
        ax.text(0.05, 0.95, f"trace: {n_trace}/{n_total}",
                transform=ax.transAxes, va="top", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor="#999", alpha=0.92))

        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_title(BRAIN_BAND_TEX_DICT[band])
        ax.set_xlabel(r"$\rho_{\mathrm{drift}}$")
        if c == 0:
            ax.set_ylabel(r"$\rho_{\mathrm{split}}$")
        ax.grid(color="#eee", lw=0.4, zorder=0)
        ax.set_axisbelow(True)

    handles = [
        Line2D([0], [0], marker="*", color="none",
               markerfacecolor="#2ca02c", markeredgecolor="k",
               markersize=12, label="cohort mean"),
        Patch(facecolor="#888", alpha=0.18, edgecolor="#444",
              label=r"cohort 1$\sigma$ ellipse"),
        Patch(facecolor="#d62728", alpha=0.10,
              label=r"trace zone ($\rho_{\mathrm{split}} > \rho_{\mathrm{drift}}$)"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="white", markeredgecolor="#1f77b4",
               markersize=7, markeredgewidth=1.1, label="patient"),
        Line2D([0], [0], marker="^", color="none",
               markerfacecolor="white", markeredgecolor="#ff7f0e",
               markersize=8, markeredgewidth=1.2, label="Pat_03 (1024 Hz)"),
    ]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.06), ncol=5, frameon=False)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


def main() -> None:
    df = pd.read_csv(SRC_CSV)
    render(df, OUT_PDF)


if __name__ == "__main__":
    main()
