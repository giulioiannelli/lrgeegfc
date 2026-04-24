#!/usr/bin/env python3
"""Evidence-chain binary heatmap — the honest-story figure.

For every patient × band, show the SIGN of each test result on a single
common scale (z-score-like, capped). Tests run left-to-right across the
evidence chain:

    [H1 task stab] [H2a-raw direct sim] [H2c residual drift]
    [H2d block persist] [H2e above drift]

Green = trace-consistent (positive contrast in the expected direction).
Red   = trace-opposite.
White = no data (e.g. Pat_14 ttest).

The visual claim:
    * H1 row universally green — task is a distinguishable state.
    * H2a-raw row universally red or near-zero — direct global
      ultrametric similarity does NOT increase toward task.
    * H2c, H2d, H2e rows universally green — residual / conditional /
      drift-floored tests do.

That is the evidence structure the paper defends.

Reads (all n=9 raw CSVs):
    h2c_ultrametric_drift_raw.csv, h2d_persistence_raw.csv,
    h2_raw_matrix_correlation_raw.csv, h2e_split_half_rho_raw.csv,
    rigorous_tests.csv   (for H1 per-patient contrast)

Writes: data/reports/imcoh_vi/figures/evidence_chain.{pdf,png}
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST,
)
from lrg_eegfc.config.paths import REPORTS_ROOT


IN_DIR = REPORTS_ROOT / "imcoh_vi"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Tests to stack, in left-to-right order on the figure. Each tuple:
#   (title, loader_fn, expected_sign_of_trace)
# Expected sign +1 means "positive value = trace-consistent".
TESTS: list[tuple[str, str]] = [
    ("H1\ntask\nstab",            "+"),
    ("H2a-raw\ndirect\nsim",      "+"),
    ("H2c\nresidual\ndrift",      "+"),
    ("H2d\nblock\npersist",       "+"),
    ("H2e\nabove\ndrift",         "+"),
]


# ─────────────────────────── per-test loaders ───────────────────────────


def _h1_values() -> pd.DataFrame:
    """Per-patient × band H1 contrast = ρ(tl,tt) − mean(cross-type), from h2_raw."""
    raw = pd.read_csv(IN_DIR / "h2_raw_matrix_correlation_raw.csv")
    PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")
    TYPE = {"rest_pre": "rest", "rest_post": "rest",
            "task_learn": "task", "task_test": "task"}
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            sub = raw[(raw["patient"] == pat) & (raw["band"] == band)]
            if sub.empty:
                rows.append({"patient": pat, "band": band, "val": np.nan})
                continue
            # ρ(task_learn, task_test): the two-phase within-task pair
            mask_task = ((sub["p1"] == "task_learn") & (sub["p2"] == "task_test")) \
                      | ((sub["p1"] == "task_test") & (sub["p2"] == "task_learn"))
            task_rho = sub[mask_task]["rho"].to_numpy()
            # Cross-type: phases with different TYPE
            cross = []
            for _, r in sub.iterrows():
                if TYPE[r["p1"]] != TYPE[r["p2"]]:
                    cross.append(r["rho"])
            cross = [c for c in cross if np.isfinite(c)]
            val = (float(task_rho[0]) - float(np.mean(cross))
                   if task_rho.size and cross else np.nan)
            rows.append({"patient": pat, "band": band, "val": val})
    return pd.DataFrame(rows)


def _h2a_raw_values() -> pd.DataFrame:
    raw = pd.read_csv(IN_DIR / "h2_raw_matrix_correlation_raw.csv")
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            sub = raw[(raw["patient"] == pat) & (raw["band"] == band)]
            m_tt   = ((sub["p1"] == "rest_post") & (sub["p2"] == "task_test")) \
                   | ((sub["p1"] == "task_test") & (sub["p2"] == "rest_post"))
            m_pre  = ((sub["p1"] == "rest_pre")  & (sub["p2"] == "rest_post")) \
                   | ((sub["p1"] == "rest_post") & (sub["p2"] == "rest_pre"))
            rho_tt  = sub[m_tt]["rho"].to_numpy()
            rho_pre = sub[m_pre]["rho"].to_numpy()
            val = (float(rho_tt[0]) - float(rho_pre[0])
                   if rho_tt.size and rho_pre.size else np.nan)
            rows.append({"patient": pat, "band": band, "val": val})
    return pd.DataFrame(rows)


def _h2c_values() -> pd.DataFrame:
    raw = pd.read_csv(IN_DIR / "h2c_ultrametric_drift_raw.csv")
    raw = raw[["patient", "band", "rho_task"]].rename(columns={"rho_task": "val"})
    return raw


def _h2d_values() -> pd.DataFrame:
    raw = pd.read_csv(IN_DIR / "h2d_persistence_raw.csv")
    return (raw.groupby(["patient", "band"])["delta_rho"].mean()
            .reset_index().rename(columns={"delta_rho": "val"}))


def _h2e_values() -> pd.DataFrame:
    raw = pd.read_csv(IN_DIR / "h2e_split_half_rho_raw.csv")
    raw["val"] = raw["rho_cross"] - raw["rho_null_drift"]
    return raw[["patient", "band", "val"]]


LOADERS = [_h1_values, _h2a_raw_values, _h2c_values, _h2d_values, _h2e_values]


# ─────────────────────────── figure ───────────────────────────


def _rescale(x: pd.DataFrame) -> pd.DataFrame:
    """Rescale by abs-max of finite values so every test lives on [-1, +1]."""
    vals = x["val"].to_numpy(dtype=float)
    finite = vals[np.isfinite(vals)]
    if finite.size == 0:
        return x
    scale = np.max(np.abs(finite))
    if scale == 0:
        return x
    x = x.copy()
    x["val"] = vals / scale
    return x


def main() -> None:
    test_frames = [_rescale(fn()) for fn in LOADERS]
    titles = [t[0] for t in TESTS]

    patients = list(PATIENTS_LIST)
    P, B, T = len(patients), len(BRAIN_BANDS_NAMES), len(test_frames)

    fig, axes = plt.subplots(
        1, T, figsize=(2.6 * T + 1.4, 0.45 * P + 1.6), dpi=160,
        gridspec_kw={"wspace": 0.12},
    )
    norm = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
    cmap = plt.get_cmap("RdYlGn")

    for ti, (ax, tdf, title) in enumerate(zip(axes, test_frames, titles)):
        mat = np.full((P, B), np.nan)
        for ip, pat in enumerate(patients):
            for ib, band in enumerate(BRAIN_BANDS_NAMES):
                hit = tdf[(tdf["patient"] == pat) & (tdf["band"] == band)]
                if len(hit) and np.isfinite(hit["val"].iloc[0]):
                    mat[ip, ib] = float(hit["val"].iloc[0])
        im = ax.imshow(mat, cmap=cmap, norm=norm, aspect="auto")
        ax.set_xticks(range(B))
        ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                            fontsize=10)
        if ti == 0:
            ax.set_yticks(range(P))
            ax.set_yticklabels(patients, fontsize=9)
        else:
            ax.set_yticks([])
        ax.set_title(title, fontsize=10, pad=10)
        # Count pos/neg
        finite = mat[np.isfinite(mat)]
        n_pos = int((finite > 0).sum())
        n_neg = int((finite < 0).sum())
        ax.set_xlabel(f"+{n_pos} / −{n_neg}", fontsize=9, labelpad=4)

    cbar = fig.colorbar(im, ax=axes, orientation="vertical",
                         pad=0.01, shrink=0.78, fraction=0.04)
    cbar.set_label("trace-consistency (per-test rescaled to [−1, +1])", fontsize=9)
    fig.suptitle(
        "Evidence chain — per-patient, per-band sign of each contrast.\n"
        "Green = trace-consistent; red = trace-opposite. "
        "Direct tests (H2a-raw) are mostly red; residual / conditional / "
        "drift-floored tests (H2c, H2d, H2e) are universally green.",
        fontsize=10, y=1.02,
    )
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"evidence_chain.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
