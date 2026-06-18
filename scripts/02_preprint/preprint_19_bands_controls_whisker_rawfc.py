#!/usr/bin/env python3
"""Raw-FC per-pair trace vs its two nulls — one-row cohort whisker, by band.

Companion to ``preprint_18_bands_joint_density_rawfc`` (Fig A, "every band
looks like a trace"). This is Fig B, the bridge that makes the substrate's
instability explicit: for each band it shows the cohort distribution of the
observed ``ρ_split^raw`` next to its two nulls —

  * ``ρ_drift^raw`` (within-session drift floor, audit_74) and
  * the matched-strength surrogate (per-patient surrogate median, audit_67) —

as dodged whiskers (median + IQR) with the per-patient points overlaid.

Reading: the matched-strength surrogate is tight near zero, but the observed
``ρ_split^raw`` spread crosses zero and overlaps the drift null, whose large
positive values fall in exactly the patients whose observed coupling is
negative (Pat_10 / Pat_13 at α/β). At the substrate the trace clears neither
check at cohort level — the apparent diagonal of Fig A cannot be separated
from within-session drift or from a strength-only rewiring.

Layout: 1 row × 6 bands (δ θ α β γ_l γ_h). Compact, not tall.

Inputs
------
data/audit/raw_fc_matched_strength/per_patient_per_band_all_bands.csv
    (obs_rho, surr_p50 — matched-strength surrogate per-patient median)
data/audit/raw_fc_matched_strength/drift_floor_per_patient.csv
    (rho_null_drift — within-session drift floor)
data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv
    (paired_wilcoxon_p — matched-strength gate, stdout only)
data/audit/raw_fc_matched_strength/drift_floor_band_stats.md  (drift p, stdout)

Output (PDF only, full vector, transparent)
-------------------------------------------
data/preprint/figures/all_bands/fig_bands_controls_whisker_rawfc.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
RAW_MS_DIR = ROOT / "data" / "audit" / "raw_fc_matched_strength"

# Three series: observed signal + the two nulls. Categorical, distinct,
# no near-white (the obs/null contrast must read on a white page).
SERIES = [
    ("obs",   r"observed $\rho_{\mathrm{split}}^{\mathrm{raw}}$", "#1f3d6e"),
    ("drift", r"drift null $\rho_{\mathrm{drift}}^{\mathrm{raw}}$", "#e08214"),
    ("ms",    r"matched-strength surrogate",                       "#888888"),
]
DX = 0.26          # horizontal dodge between the three series within a band
JIT = 0.055        # per-patient point jitter half-width


def load_long() -> pd.DataFrame:
    """Per (patient, band): obs, drift null, matched-strength surrogate median."""
    ms = pd.read_csv(RAW_MS_DIR / "per_patient_per_band_all_bands.csv")[
        ["patient", "band", "obs_rho", "surr_p50"]
    ]
    drift = pd.read_csv(RAW_MS_DIR / "drift_floor_per_patient.csv")[
        ["patient", "band", "rho_null_drift"]
    ]
    df = ms.merge(drift, on=["patient", "band"], how="inner")
    return df.rename(columns={"obs_rho": "obs", "rho_null_drift": "drift",
                              "surr_p50": "ms"})


def main() -> Path:
    df = load_long()
    cohort = pd.read_csv(RAW_MS_DIR / "cohort_summary_all_bands.csv")[
        ["band", "paired_wilcoxon_p"]
    ].set_index("band")["paired_wilcoxon_p"].to_dict()

    rng = np.random.default_rng(0)

    fig, ax = plt.subplots(figsize=(15.0, 4.4))

    for bi, band in enumerate(BAND_ORDER):
        b = df[df.band == band]
        for si, (key, _lbl, color) in enumerate(SERIES):
            vals = b[key].to_numpy()
            vals = vals[np.isfinite(vals)]
            x0 = bi + (si - 1) * DX
            med = float(np.median(vals))
            q1, q3 = np.percentile(vals, [25, 75])
            # per-patient points (jittered)
            jit = rng.uniform(-JIT, JIT, size=vals.size)
            ax.scatter(np.full(vals.size, x0) + jit, vals, s=16,
                       color=color, alpha=0.45, edgecolor="none", zorder=2)
            # IQR bar
            ax.plot([x0, x0], [q1, q3], color=color, lw=5.0, alpha=0.55,
                    solid_capstyle="round", zorder=3)
            # median tick
            ax.plot([x0 - 0.085, x0 + 0.085], [med, med], color=color, lw=2.6,
                    zorder=4)
            ax.plot(x0, med, "o", ms=7, color=color, mec="white", mew=1.1,
                    zorder=5)
        # stdout summary (no on-figure text per project convention)
        ms_p = cohort.get(band, np.nan)
        d_diff = (b["obs"] - b["drift"]).to_numpy()
        _, d_p = wilcoxon_z(d_diff[np.isfinite(d_diff)])
        print(f"  {band:11s} med obs={b['obs'].median():+.3f} "
              f"drift={b['drift'].median():+.3f} ms={b['ms'].median():+.3f}  "
              f"| gate_p(MS)={ms_p:.3f}  drift_p={d_p:.3f}")

    ax.axhline(0.0, color="0.55", lw=0.9, ls="--", zorder=1)
    ax.set_xticks(range(len(BAND_ORDER)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER], fontsize=14)
    ax.set_xlim(-0.6, len(BAND_ORDER) - 0.4)
    ax.set_ylabel(r"per-pair cross-phase rank coupling $\rho$")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", labelsize=11)

    handles = [Line2D([0], [0], marker="o", color=c, lw=0, mec="white", mew=1.0,
                      ms=8, label=lbl) for _key, lbl, c in SERIES]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.02), ncol=len(SERIES), frameon=False)
    fig.tight_layout(rect=(0, 0.05, 1, 1))

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_bands_controls_whisker_rawfc.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
