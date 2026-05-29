#!/usr/bin/env python3
"""Section 6 — Final radar chart of scalar VI contrast (H2a).

N=5 patients (era artifact). Pat_03 is rendered with the same line
style and marker as every other patient — the earlier "1024 Hz
outlier" framing was retired 2026-05-18 (see
``feedback_pat03_no_dropout``; Pat_03 is handled at the config layer
only).
Output: data/figures/section6/radar_chart_final.pdf + .md
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style
use_lrg_style()

ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.config.paths import FIGURES_ROOT

OUTDIR = FIGURES_ROOT / "section6"
OUTDIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
})

# ── Data (H2a Scalar A, all 5 patients) ──────────────────────────────
df = pd.read_csv("data/wp_scalar_vi/scalar_table_full.csv")
df = df[df["hypothesis"] == "H2a"]

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_LABELS = [r"$\delta$", r"$\theta$", r"$\alpha$", r"$\beta$",
               r"$\gamma_l$", r"$\gamma_h$"]
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
PAT_COLORS = {
    "Pat_02": "#1f77b4",
    "Pat_03": "#ff7f0e",
    "Pat_05": "#2ca02c",
    "Pat_07": "#d62728",
    "Pat_08": "#9467bd",
}
PAT_SHORT = {"Pat_02": "P2", "Pat_03": "P3", "Pat_05": "P5",
             "Pat_07": "P7", "Pat_08": "P8"}

data = np.zeros((len(PATIENTS), len(BANDS)))
for ip, pat in enumerate(PATIENTS):
    for ib, band in enumerate(BANDS):
        row = df[(df["patient"] == pat) & (df["band"] == band)]
        data[ip, ib] = row["scalar_A_mean"].values[0]

mean_vals = data.mean(axis=0)
vmax = np.ceil(np.max(np.abs(data)) * 10) / 10  # round up to nearest 0.1

# ── Radar chart ──────────────────────────────────────────────────────
n_bands = len(BANDS)
angles = np.linspace(0, 2 * np.pi, n_bands, endpoint=False).tolist()
angles_closed = angles + [angles[0]]

fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(polar=True))

ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)

# Anti-trace zone: gray fill inside zero circle
theta_fill = np.linspace(0, 2 * np.pi, 300)
ax.fill(theta_fill, np.full_like(theta_fill, 0), color="#ededed", zorder=0)

# Zero circle — prominent dashed black
ax.plot(theta_fill, np.zeros_like(theta_fill), color="black",
        linestyle="--", lw=1.8, alpha=0.7, zorder=2)

# Radial range and gridlines
ax.set_rlim(-vmax, vmax)
rticks = np.arange(-vmax, vmax + 0.05, 0.1).round(1).tolist()
ax.set_rticks(rticks)
ax.set_yticklabels([f"{v:+.1f}" if v != 0 else "0" for v in rticks],
                    fontsize=7, color="#999999")
ax.set_rlabel_position(22)

# Minimal grid
ax.grid(color="#e0e0e0", linewidth=0.4, linestyle="-")
ax.spines["polar"].set_visible(False)

# Patient polygons — Pat_03 uniformly rendered with the other patients
# (Pat_03 outlier framing retired 2026-05-18; OUTLIER kept as dead
# constant for downstream readability of the legacy script).
for ip, pat in enumerate(PATIENTS):
    vals = data[ip].tolist() + [data[ip, 0]]
    ls = "-"
    marker = "o"
    lbl = PAT_SHORT[pat]
    ax.plot(angles_closed, vals, ls, marker=marker, color=PAT_COLORS[pat],
            lw=1.2, markersize=3.5, alpha=0.85, label=lbl, zorder=3)
    ax.fill(angles_closed, vals, color=PAT_COLORS[pat], alpha=0.08, zorder=1)

# Mean polygon
mean_closed = mean_vals.tolist() + [mean_vals[0]]
ax.plot(angles_closed, mean_closed, "-s", color="black", lw=2.5,
        markersize=5, label="Mean", zorder=4)

# Band labels
ax.set_xticks(angles)
ax.set_xticklabels(BAND_LABELS, fontsize=15, fontweight="bold")

# Alpha annotation
alpha_idx = BANDS.index("alpha")
alpha_angle = angles[alpha_idx]
ax.annotate(
    "5/5\n$p$ = 0.031",
    xy=(alpha_angle, mean_vals[alpha_idx]),
    xytext=(alpha_angle, vmax * 0.82),
    fontsize=8.5, fontweight="bold", ha="center", va="center",
    color="#b8860b",
    arrowprops=dict(arrowstyle="-", color="#b8860b", lw=0.8),
    zorder=5,
)

# Legend
leg = ax.legend(loc="upper right", bbox_to_anchor=(1.28, 1.08),
                fontsize=9, frameon=True, framealpha=0.95,
                edgecolor="#cccccc", handlelength=1.5)

fig.savefig(OUTDIR / "radar_chart_final.pdf", bbox_inches="tight", dpi=200)
plt.close(fig)
print("Saved radar_chart_final.pdf")

# ── Markdown companion ───────────────────────────────────────────────
md = """\
# radar_chart_final

## What the figure shows
Radar chart of H2a scalar VI contrast (Scalar A) across 6 frequency bands
for N=5 patients (era-pinned cohort; not current n=10). Each polygon =
one patient's profile, rendered with the same line style. Thick black =
mean. Dashed circle = zero (no trace baseline). Outside = positive task
trace; gray fill inside = anti-trace.

## Key result
Alpha is the only band where all 5 patients show positive task trace (5/5,
binomial p = 0.031).

## How to read it
Vertices outside the dashed zero circle = positive trace (rest_post closer to task
than rest_pre). Inside gray zone = anti-trace. The alpha vertex consistently
points outward for all patients. Other bands show mixed profiles.

## Data
H2a Scalar A from data/wp_scalar_vi/scalar_table_full.csv.
"""
(OUTDIR / "radar_chart_final.md").write_text(md)
print("Saved radar_chart_final.md")
