"""KC cohort fingerprint — 10×6 sign-flag matrix.

Per (patient, band) cell, encode the trace state at λ=0 and λ=1:

  state = 2 (dark red)   both real_T_KC < null_T_KC  (full trace)
  state = 1 (light red)  exactly one of λ=0, λ=1 below null (partial)
  state = 0 (grey)       neither below null (no trace)

Scale-free: cell color depends only on the SIGN of `real - null`, never
on its magnitude — so per-patient outliers cannot dominate the visual.

Reader scans columns: β reads as a single solid stripe (10 dark cells);
δ / θ / γ_l have mixed columns.

Data source: data/audit/section5_v2_kc_controls/per_patient_table.csv
Output:      data/reports/preprint/cohort_measures/kc_sign_flag_matrix.pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style

use_lrg_style()


ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "data/audit/section5_v2_kc_controls/per_patient_table.csv"
OUT_DIR = ROOT / "data/reports/preprint/cohort_measures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_PDF = OUT_DIR / "kc_sign_flag_matrix.pdf"

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]


def build_state_matrix(df: pd.DataFrame) -> np.ndarray:
    n_pat = len(PATIENTS)
    n_band = len(BAND_ORDER)
    M = np.zeros((n_pat, n_band), dtype=int)
    for j, band in enumerate(BAND_ORDER):
        for i, pat in enumerate(PATIENTS):
            row0 = df[(df["band"] == band) & (df["patient"] == pat) & (df["lam"] == 0.0)]
            row1 = df[(df["band"] == band) & (df["patient"] == pat) & (df["lam"] == 1.0)]
            if not len(row0) or not len(row1):
                M[i, j] = -1
                continue
            b0 = bool(row0.iloc[0]["real_below_null"])
            b1 = bool(row1.iloc[0]["real_below_null"])
            M[i, j] = int(b0) + int(b1)
    return M


def render(M: np.ndarray, out_pdf: Path) -> None:
    cmap = ListedColormap(["#dddddd", "#f4a8a0", "#a82a1a"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

    fig, ax = plt.subplots(figsize=(5.4, 5.4))
    im = ax.imshow(M, cmap=cmap, norm=norm, aspect="equal",
                   interpolation="nearest")

    # Cell annotations: number of λ-axes below null (0/1/2)
    for (i, j), v in np.ndenumerate(M):
        if v < 0:
            ax.text(j, i, "—", ha="center", va="center",
                    fontsize=9, color="#999")
        else:
            ax.text(j, i, str(int(v)), ha="center", va="center",
                    fontsize=9, color="white" if v >= 1 else "#444",
                    fontweight="bold" if v == 2 else "normal")

    # Column counts (n_full_trace, n_partial) per band, footer row
    n_full = (M == 2).sum(axis=0)
    n_partial = (M == 1).sum(axis=0)
    for j, band in enumerate(BAND_ORDER):
        ax.text(j, len(PATIENTS) - 0.4 + 0.95, f"{n_full[j]}",
                ha="center", va="top",
                fontsize=11, fontweight="bold", color="#a82a1a")
        ax.text(j, len(PATIENTS) - 0.4 + 1.55, f"+{n_partial[j]}",
                ha="center", va="top",
                fontsize=8, color="#a82a1a")

    ax.text(-0.7, len(PATIENTS) - 0.4 + 0.95, r"$n_{\mathrm{full}}$",
            ha="right", va="top", fontsize=9, color="#a82a1a")
    ax.text(-0.7, len(PATIENTS) - 0.4 + 1.55, r"$+n_{\mathrm{partial}}$",
            ha="right", va="top", fontsize=8, color="#a82a1a")

    ax.set_xticks(range(len(BAND_ORDER)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER], fontsize=12)
    ax.set_yticks(range(len(PATIENTS)))
    ax.set_yticklabels(PATIENTS)
    ax.tick_params(axis="x", which="both", length=0, pad=4)
    ax.tick_params(axis="y", which="both", length=0)
    ax.set_xlabel("band")
    ax.set_ylabel("patient")
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax.set_xlim(-0.5, len(BAND_ORDER) - 0.5)
    ax.set_ylim(len(PATIENTS) + 1.6, -0.5)

    handles = [
        Patch(facecolor="#a82a1a", edgecolor="0.4", label="2 = full trace (both λ)"),
        Patch(facecolor="#f4a8a0", edgecolor="0.4", label="1 = partial (one λ)"),
        Patch(facecolor="#dddddd", edgecolor="0.4", label="0 = no trace"),
    ]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.02), ncol=3, frameon=False)

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


def main() -> None:
    M = build_state_matrix(pd.read_csv(SRC))
    render(M, OUT_PDF)


if __name__ == "__main__":
    main()
