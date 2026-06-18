#!/usr/bin/env python3
"""Raw-FC per-pair trace vs its noise floor — immediate one-row read, by band.

Fig B (immediate variant). Companion to ``preprint_18`` (Fig A, joint
density: "every band looks like a trace") and the richer fallback
``preprint_19`` (three-distribution whiskers). The point of this version is
a single glanceable comparison per band:

    observed effect  vs  noise-floor zone

For each band the observed cohort coupling ``ρ_split^raw`` is drawn as one
clean 95%-bootstrap-CI bar; behind it a shaded gray zone spans the
noise floor — the union of the 95% CIs of the two nulls, the within-session
drift floor ``ρ_drift^raw`` and the matched-strength surrogate. The two null
means are drawn as horizontal lines inside the zone so their near-equality
is visible (the drift bias is no smaller than the strength bias).

Immediate read: when the observed CI bar sits inside / overlaps the gray
zone, the trace is not separated from what drift or strength-only rewiring
already produce — i.e. at the raw substrate the apparent diagonal of Fig A
cannot be ruled a task-specific reorganization. (Contrast: at the LRG
cophenetic layer the corresponding bar clears the zone at α/β.)

Layout: 1 row × 6 bands (δ θ α β γ_l γ_h). Compact.

Inputs
------
data/audit/raw_fc_matched_strength/per_patient_per_band_all_bands.csv
    (obs_rho, surr_p50 — matched-strength surrogate per-patient median)
data/audit/raw_fc_matched_strength/drift_floor_per_patient.csv
    (rho_null_drift — within-session drift floor)

Output (PDF only, full vector, transparent)
-------------------------------------------
data/preprint/figures/all_bands/fig_bands_signal_vs_noise_rawfc.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics.hypothesis import boot_ci_mean
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
RAW_MS_DIR = ROOT / "data" / "audit" / "raw_fc_matched_strength"

C_OBS = "#1f3d6e"     # observed
C_DRIFT = "#e08214"   # drift null
C_MS = "#5a5a5a"      # matched-strength null
C_ZONE = "#c4c4c4"    # noise-floor zone fill

HALF_W = 0.30         # half-width of a band's zone / bar cluster
JIT = 0.05


def load_long() -> pd.DataFrame:
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
    rng = np.random.default_rng(0)

    fig, ax = plt.subplots(figsize=(14.0, 4.6))

    for bi, band in enumerate(BAND_ORDER):
        b = df[df.band == band]
        obs = b["obs"].to_numpy()
        drift = b["drift"].to_numpy()
        ms = b["ms"].to_numpy()

        o_m, o_lo, o_hi = boot_ci_mean(obs, rng=rng)
        d_m, d_lo, d_hi = boot_ci_mean(drift, rng=rng)
        m_m, m_lo, m_hi = boot_ci_mean(ms, rng=rng)

        # Noise-floor zone = union of the two nulls' 95% CIs.
        z_lo = min(d_lo, m_lo)
        z_hi = max(d_hi, m_hi)
        ax.fill_between([bi - HALF_W, bi + HALF_W], z_lo, z_hi,
                        color=C_ZONE, alpha=0.55, lw=0, zorder=1)
        # The two null means as horizontal lines inside the zone — their
        # proximity shows the drift bias ≈ the strength bias.
        ax.plot([bi - HALF_W, bi + HALF_W], [d_m, d_m], color=C_DRIFT,
                lw=2.4, alpha=0.95, zorder=3)
        ax.plot([bi - HALF_W, bi + HALF_W], [m_m, m_m], color=C_MS,
                lw=2.0, ls=(0, (4, 2)), alpha=0.95, zorder=3)

        # Faint per-patient observed points (texture + the inversions).
        ax.scatter(np.full(obs.size, bi) + rng.uniform(-JIT, JIT, obs.size),
                   obs, s=14, color=C_OBS, alpha=0.22, edgecolor="none",
                   zorder=2)
        # Observed mean ± 95% CI — the one clean bar.
        ax.errorbar(bi, o_m, yerr=[[o_m - o_lo], [o_hi - o_m]],
                    fmt="o", color=C_OBS, ecolor=C_OBS, elinewidth=2.4,
                    capsize=6, capthick=2.4, ms=10, mfc=C_OBS, mec="white",
                    mew=1.4, zorder=5)

        clears = o_lo > z_hi
        print(f"  {band:11s} obs={o_m:+.3f} [{o_lo:+.3f},{o_hi:+.3f}]  "
              f"zone=[{z_lo:+.3f},{z_hi:+.3f}]  "
              f"drift={d_m:+.3f} ms={m_m:+.3f}  clears_zone={clears}")

    ax.axhline(0.0, color="0.6", lw=0.9, ls=":", zorder=0)
    ax.set_xticks(range(len(BAND_ORDER)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER], fontsize=15)
    ax.set_xlim(-0.6, len(BAND_ORDER) - 0.4)
    ax.set_ylabel(r"per-pair cross-phase coupling $\rho$")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", labelsize=11)

    handles = [
        Line2D([0], [0], marker="o", color=C_OBS, lw=0, mec="white", mew=1.2,
               ms=10, label=r"observed $\rho_{\mathrm{split}}^{\mathrm{raw}}$ (95% CI)"),
        Patch(facecolor=C_ZONE, alpha=0.6, label=r"noise floor (95% CI, drift $\cup$ matched-strength)"),
        Line2D([0], [0], color=C_DRIFT, lw=2.4,
               label=r"drift-null mean $\rho_{\mathrm{drift}}^{\mathrm{raw}}$"),
        Line2D([0], [0], color=C_MS, lw=2.0, ls=(0, (4, 2)),
               label="matched-strength surrogate mean"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.04),
               ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 1))

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_bands_signal_vs_noise_rawfc.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
