#!/usr/bin/env python3
"""Per-pair "null triangle" by band — signal vs its nulls, per-patient
connectors, for the raw substrate AND the LRG cophenetic layer (one script,
``--layer {raw,coph}``, identical styling so the two figures are a clean swap).

In the spirit of panel (c) of the cophenetic ``fig_<band>_rho_split`` figures.
For each band, three boxes —

    ρ_split    (observed cross-phase coupling)
    ρ_drift    (within-session drift null, no task)
    ρ_xprobe   (cross-probe restriction of ρ_split)

— with a gray shaded **matched-strength null floor** (cohort p5–p95 of the
strength-preserving surrogate) behind them, and each patient's three values
joined by a connector line.

CONNECTOR ENCODING (magnitude, not just sign). The connector is red when the
drift null exceeds that patient's signal (ρ_drift > ρ_split = an inversion),
and its **thickness + saturation scale with how far drift exceeds the signal**
(ρ_drift − ρ_split). Grey, thin = signal above its drift null. The magnitude
scale (MARGIN_REF) is shared across layers, so the raw→LRG difference is the
one the eye reads: the inverting patients largely persist in NUMBER between
raw and cophenetic, but their inversions **collapse in magnitude** at the LRG
layer (thick red spikes → thin) — which is why the cohort drift test clears at
α/β/γ_l for cophenetic and at no band for raw. (Do NOT describe this as "red
turns gray": the dissenting patients still dissent, only more mildly.)

Layers
------
raw  : raw |ImCoh| adjacency. split/xprobe from audit_74 cross_probe_per_patient,
       drift from audit_74 drift_floor_per_patient, MS floor from audit_67
       per_patient_per_band_all_bands (surr p5/p95).
coph : LRG cophenetic D. split/drift/xprobe from imcoh_continuous_trace
       controls_summary (Run A/C/B), MS floor from audit_63
       matched_strength_surrogate_split_baseline per_patient_per_band.
       (Per GAP-2, the headline obs number is audit_63's; here split/drift/
       xprobe all come from controls_summary so the per-patient triple is
       internally consistent — the ~0.01 α difference is immaterial visually.)

Layout: 1 row × 6 band panels (δ θ α β γ_l γ_h), shared fixed y-range.

Output (PDF only, full vector, transparent)
-------------------------------------------
data/preprint/figures/all_bands/fig_bands_null_triangle_<rawfc|coph>.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]

COLS = {"split": "#2f4b7c", "drift": "#b3992d", "xprobe": "#6d519c"}
COL_ORDER = ["split", "drift", "xprobe"]
C_INVERT = "#c0392b"   # connector: drift null beats the signal
C_KEEP = "#9aa0a6"     # connector: signal above its drift null

# Shared magnitude scale: an inversion margin (ρ_drift − ρ_split) of MARGIN_REF
# maps to the thickest / most saturated red. Identical across layers so the
# raw→cophenetic collapse in inversion magnitude is visible at one glance.
MARGIN_REF = 0.6
YLIM = (-0.6, 1.0)

RAW_DIR = ROOT / "data" / "audit" / "raw_fc_matched_strength"
COPH_CT = ROOT / "data" / "reports" / "imcoh_continuous_trace"
COPH_MS = ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"


# ---------------------------------------------------------------------------
# Layer loaders → (triples df[patient,band,split,drift,xprobe], zones, gate p)
# ---------------------------------------------------------------------------
def _zone_from_quantiles(df: pd.DataFrame) -> dict:
    return {b: (float(df[df.band == b]["surr_p5"].median()),
               float(df[df.band == b]["surr_p95"].median()))
            for b in BAND_ORDER}


def load_raw():
    cp = pd.read_csv(RAW_DIR / "cross_probe_per_patient.csv")[
        ["patient", "band", "rho_full", "rho_cross"]]
    dr = pd.read_csv(RAW_DIR / "drift_floor_per_patient.csv")[
        ["patient", "band", "rho_null_drift"]]
    df = cp.merge(dr, on=["patient", "band"]).rename(
        columns={"rho_full": "split", "rho_cross": "xprobe",
                 "rho_null_drift": "drift"})
    zones = _zone_from_quantiles(
        pd.read_csv(RAW_DIR / "per_patient_per_band_all_bands.csv"))
    gate = pd.read_csv(RAW_DIR / "cohort_summary_all_bands.csv").set_index(
        "band")["paired_wilcoxon_p"].to_dict()
    return df, zones, gate


def load_coph():
    c = pd.read_csv(COPH_CT / "controls_summary.csv")[
        ["patient", "band", "rho_split", "rho_null_drift", "rho_split_cross_probe"]]
    df = c.rename(columns={"rho_split": "split", "rho_null_drift": "drift",
                           "rho_split_cross_probe": "xprobe"}).dropna(
        subset=["split", "drift", "xprobe"])
    zones = _zone_from_quantiles(
        pd.read_csv(COPH_MS / "per_patient_per_band.csv"))
    gate = pd.read_csv(COPH_MS / "cohort_summary.csv").set_index(
        "band")["paired_wilcoxon_p"].to_dict()
    return df, zones, gate


LAYERS = {
    "raw":  dict(load=load_raw,  sup="raw",
                 ylab=r"per-pair cross-phase coupling $\rho$ (raw $|\mathrm{ImCoh}|$)",
                 out="fig_bands_null_triangle_rawfc.pdf"),
    "coph": dict(load=load_coph, sup="coph",
                 ylab=r"per-pair cross-phase coupling $\rho$ (LRG cophenetic $D$)",
                 out="fig_bands_null_triangle_coph.pdf"),
}


def connector_style(margin: float):
    """(color, lw, alpha) for a patient connector. margin = ρ_drift − ρ_split.

    margin > 0 → inversion: red, thickness + alpha grow with margin.
    margin ≤ 0 → signal above drift: thin neutral grey.
    """
    if margin > 0:
        t = min(margin / MARGIN_REF, 1.0)
        return C_INVERT, 1.0 + 4.0 * t, 0.45 + 0.50 * t
    return C_KEEP, 0.7, 0.30


def _style_box(bp, color):
    for patch in bp["boxes"]:
        patch.set(facecolor=color, alpha=0.22, edgecolor=color, linewidth=1.4)
    for med in bp["medians"]:
        med.set(color=color, linewidth=2.4)
    for el in ("whiskers", "caps"):
        for art in bp[el]:
            art.set(color=color, linewidth=1.2)


def main(layer: str) -> Path:
    cfg = LAYERS[layer]
    df, zones, gate = cfg["load"]()
    sup = cfg["sup"]
    box_label = {k: rf"$\rho_{{\mathrm{{{k}}}}}^{{\mathrm{{{sup}}}}}$"
                 for k in COL_ORDER}
    rng = np.random.default_rng(0)

    fig, axes = plt.subplots(1, len(BAND_ORDER), figsize=(18.0, 4.4),
                             sharey=True)

    for ax, band in zip(axes, BAND_ORDER):
        b = df[df.band == band].reset_index(drop=True)
        data = [b[k].to_numpy() for k in COL_ORDER]
        offs = rng.uniform(-0.12, 0.12, size=len(b))

        z_lo, z_hi = zones[band]
        ax.axhspan(z_lo, z_hi, color="#cfcfcf", alpha=0.6, lw=0, zorder=0)

        n_inv = 0
        for pi in range(len(b)):
            row = b.iloc[pi]
            margin = float(row["drift"] - row["split"])
            n_inv += int(margin > 0)
            color, lw, alpha = connector_style(margin)
            ax.plot([0 + offs[pi], 1 + offs[pi], 2 + offs[pi]],
                    [row["split"], row["drift"], row["xprobe"]],
                    color=color, lw=lw, alpha=alpha, zorder=2,
                    solid_capstyle="round")

        bp = ax.boxplot(data, positions=[0, 1, 2], widths=0.48,
                        patch_artist=True, showcaps=True, showfliers=False,
                        zorder=3)
        for k, patch_i in zip(COL_ORDER, range(3)):
            sub = {key: [bp[key][patch_i]] if key in ("boxes", "medians")
                   else bp[key][2 * patch_i: 2 * patch_i + 2]
                   for key in ("boxes", "medians", "whiskers", "caps")}
            _style_box(sub, COLS[k])

        for k_i, k in enumerate(COL_ORDER):
            ax.scatter(k_i + offs, b[k].to_numpy(), s=22, facecolor="white",
                       edgecolor=COLS[k], linewidth=1.1, alpha=0.95, zorder=4)

        ax.axhline(0.0, color="0.55", lw=0.9, ls=":", zorder=1)
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=18, pad=6)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels([box_label[k] for k in COL_ORDER], fontsize=11)
        ax.set_xlim(-0.6, 2.6)
        ax.set_ylim(*YLIM)
        ax.spines[["top", "right"]].set_visible(False)

        d_diff = (b["split"] - b["drift"]).to_numpy()
        inv_sum = float((b["drift"] - b["split"])[b["drift"] > b["split"]].sum())
        _, d_p = wilcoxon_z(d_diff[np.isfinite(d_diff)])
        print(f"  [{layer}] {band:11s} inversions={n_inv}/{len(b)} "
              f"inv_magnitude_sum={inv_sum:+.3f}  med split={b['split'].median():+.3f} "
              f"drift={b['drift'].median():+.3f} xprobe={b['xprobe'].median():+.3f}  "
              f"drift_p={d_p:.3f} gate_p(MS)={gate.get(band, np.nan):.3f}")

    axes[0].set_ylabel(cfg["ylab"])
    axes[0].tick_params(axis="y", labelsize=11)

    handles = [
        Patch(facecolor="#cfcfcf", alpha=0.75,
              label="matched-strength null (p5–p95)"),
        Line2D([0], [0], color=C_INVERT, lw=3.2,
               label=r"inversion $\rho_{\mathrm{drift}}>\rho_{\mathrm{split}}$ "
                     r"(thicker = drift exceeds signal by more)"),
        Line2D([0], [0], color=C_KEEP, lw=1.2,
               label=r"signal above its drift null"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.03),
               ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 1))

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / cfg["out"]
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  [{layer}] Saved: {out}")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--layer", choices=list(LAYERS), default="raw")
    ap.add_argument("--both", action="store_true", help="render raw and coph")
    args = ap.parse_args()
    if args.both:
        for lyr in ("raw", "coph"):
            main(lyr)
    else:
        main(args.layer)
