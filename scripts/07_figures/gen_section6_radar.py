#!/usr/bin/env python3
"""Section 6 — Figure 1: Radar chart of scalar VI contrast (H2a).

Produces three variants:
  v1: Standard radar (polar polygon)
  v2: Half-radar with radial bars and patient dots
  v3: Circular bar chart (sectors per band, ring per patient)

Data source: data/wp_scalar_vi/scalar_table_full.csv (H2a, scalar_A_mean)
Output:      figures/section6/radar_chart_v{1,2,3}.pdf + .md
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

from lrg_eegfc.config.paths import FIGURES_ROOT

OUTDIR = FIGURES_ROOT / "section6"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ── Data ─────────────────────────────────────────────────────────────
df = pd.read_csv("data/wp_scalar_vi/scalar_table_full.csv")
df = df[df["hypothesis"] == "H2a"]

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_LABELS = [r"$\delta$", r"$\theta$", r"$\alpha$", r"$\beta$",
               r"$\gamma_l$", r"$\gamma_h$"]
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
PAT_COLORS = {
    "Pat_02": "#1f77b4",  # blue
    "Pat_03": "#ff7f0e",  # orange
    "Pat_05": "#2ca02c",  # green
    "Pat_07": "#d62728",  # red
    "Pat_08": "#9467bd",  # purple
}
PAT_SHORT = {p: p.replace("Pat_0", "P") for p in PATIENTS}

# Build data matrix: patients × bands
data = np.zeros((len(PATIENTS), len(BANDS)))
for ip, pat in enumerate(PATIENTS):
    for ib, band in enumerate(BANDS):
        row = df[(df["patient"] == pat) & (df["band"] == band)]
        data[ip, ib] = row["scalar_A_mean"].values[0]

mean_vals = data.mean(axis=0)

# Data range for axis limits
vmax = np.ceil(np.max(np.abs(data)) * 20) / 20  # round up to nearest 0.05
vmax = max(vmax, 0.30)

# ═════════════════════════════════════════════════════════════════════
# VARIANT 1: Standard radar chart
# ═════════════════════════════════════════════════════════════════════
print("Generating v1: Standard radar chart...")

n_bands = len(BANDS)
angles = np.linspace(0, 2 * np.pi, n_bands, endpoint=False).tolist()
# Close the polygon
angles_closed = angles + [angles[0]]

fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(polar=True))

# Start from top (90°), go clockwise
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)

# Gray fill inside the zero circle to mark "anti-trace" region
theta_fill = np.linspace(0, 2 * np.pi, 200)
ax.fill(theta_fill, np.full_like(theta_fill, 0), color="#f0f0f0", zorder=0)

# Dashed zero circle
ax.plot(theta_fill, np.full_like(theta_fill, 0), "k--", lw=1.2,
        alpha=0.6, zorder=2, label="_zero")

# Radial gridlines — light gray
ax.set_rlim(-vmax, vmax)
rticks = np.arange(-0.4, 0.5, 0.1)
rticks = rticks[(rticks >= -vmax - 0.01) & (rticks <= vmax + 0.01)]
ax.set_rticks(rticks)
ax.set_yticklabels([f"{v:.1f}" if v != 0 else "0" for v in rticks],
                    fontsize=7, color="#888888")
ax.tick_params(axis="y", labelsize=7, colors="#aaaaaa")
ax.set_rlabel_position(30)

# Grid styling
ax.grid(color="#dddddd", linewidth=0.5, linestyle="-")
ax.spines["polar"].set_visible(False)

# Patient polygons
for ip, pat in enumerate(PATIENTS):
    vals = data[ip].tolist() + [data[ip, 0]]
    ax.plot(angles_closed, vals, "-o", color=PAT_COLORS[pat], lw=1.3,
            markersize=4, alpha=0.85, label=PAT_SHORT[pat], zorder=3)
    ax.fill(angles_closed, vals, color=PAT_COLORS[pat], alpha=0.10, zorder=1)

# Mean polygon — thick black
mean_closed = mean_vals.tolist() + [mean_vals[0]]
ax.plot(angles_closed, mean_closed, "-s", color="black", lw=2.5,
        markersize=5, label="Mean", zorder=4)

# Band labels at axis tips
ax.set_xticks(angles)
ax.set_xticklabels(BAND_LABELS, fontsize=14, fontweight="bold")

# Statistical annotation on alpha axis
alpha_idx = BANDS.index("alpha")
alpha_angle = angles[alpha_idx]
alpha_r = vmax * 0.88
ax.annotate("5/5\n$p$=0.031", xy=(alpha_angle, mean_vals[alpha_idx]),
            xytext=(alpha_angle, alpha_r),
            fontsize=8, fontweight="bold", ha="center", va="center",
            color="#b8860b",
            arrowprops=dict(arrowstyle="-", color="#b8860b", lw=0.8),
            zorder=5)

# Legend
leg = ax.legend(loc="upper right", bbox_to_anchor=(1.28, 1.08),
                fontsize=9, frameon=True, framealpha=0.9,
                edgecolor="#cccccc", handlelength=1.5)

fig.savefig(OUTDIR / "radar_chart_v1.pdf", bbox_inches="tight", dpi=200)
plt.close(fig)
print("  Saved radar_chart_v1.pdf")


# ═════════════════════════════════════════════════════════════════════
# VARIANT 2: Half-radar with radial bars + patient dots
# ═════════════════════════════════════════════════════════════════════
print("Generating v2: Half-radar with bars...")

fig, ax = plt.subplots(figsize=(9, 5.5), subplot_kw=dict(polar=True))

# Semicircle: bands spread over π (top half)
half_angles = np.linspace(np.pi, 0, n_bands, endpoint=True)

ax.set_theta_offset(0)
ax.set_theta_direction(1)
ax.set_thetamin(0)
ax.set_thetamax(180)

# Zero circle (semicircle)
th = np.linspace(0, np.pi, 200)
ax.fill_between(th, -vmax, 0, color="#f5f5f5", zorder=0)
ax.plot(th, np.zeros_like(th), "k--", lw=1.2, alpha=0.6, zorder=2)

# Radial gridlines
ax.set_rlim(-vmax, vmax)
rticks = np.arange(-0.4, 0.5, 0.1)
rticks = rticks[(rticks >= -vmax - 0.01) & (rticks <= vmax + 0.01)]
ax.set_rticks(rticks)
ax.set_yticklabels([f"{v:.1f}" if v != 0 else "0" for v in rticks],
                    fontsize=7, color="#888888")
ax.set_rlabel_position(90)
ax.grid(color="#dddddd", linewidth=0.5)
ax.spines["polar"].set_visible(False)

# Bar width
bar_w = (np.pi / n_bands) * 0.6

# Mean bars as background
for ib in range(n_bands):
    color = "#333333"
    ax.bar(half_angles[ib], mean_vals[ib], width=bar_w, bottom=0,
           color=color, alpha=0.12, edgecolor="none", zorder=1)

# Patient dots along each bar
offsets = np.linspace(-bar_w * 0.35, bar_w * 0.35, len(PATIENTS))
for ip, pat in enumerate(PATIENTS):
    for ib in range(n_bands):
        ax.scatter(half_angles[ib] + offsets[ip], data[ip, ib],
                   c=PAT_COLORS[pat], s=40, zorder=4, edgecolors="white",
                   linewidths=0.5)
    # Invisible line for legend
    ax.plot([], [], "o", color=PAT_COLORS[pat], markersize=5,
            label=PAT_SHORT[pat])

# Mean markers
for ib in range(n_bands):
    ax.scatter(half_angles[ib], mean_vals[ib], c="black", s=60,
               marker="D", zorder=5, edgecolors="white", linewidths=0.8)
ax.plot([], [], "D", color="black", markersize=6, label="Mean")

# Band labels
ax.set_xticks(half_angles)
ax.set_xticklabels(BAND_LABELS, fontsize=14, fontweight="bold")

# Alpha annotation
alpha_ib = BANDS.index("alpha")
ax.annotate("5/5, $p$=0.031",
            xy=(half_angles[alpha_ib], mean_vals[alpha_ib]),
            xytext=(half_angles[alpha_ib], vmax * 0.85),
            fontsize=8, fontweight="bold", ha="center", color="#b8860b",
            arrowprops=dict(arrowstyle="-", color="#b8860b", lw=0.8),
            zorder=6)

leg = ax.legend(loc="upper right", bbox_to_anchor=(1.22, 1.05),
                fontsize=9, frameon=True, framealpha=0.9,
                edgecolor="#cccccc")

fig.savefig(OUTDIR / "radar_chart_v2.pdf", bbox_inches="tight", dpi=200)
plt.close(fig)
print("  Saved radar_chart_v2.pdf")


# ═════════════════════════════════════════════════════════════════════
# VARIANT 3: Circular bar chart (sectors per band, stacked rings)
# ═════════════════════════════════════════════════════════════════════
print("Generating v3: Circular bar chart...")

fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))

ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)

sector_width = 2 * np.pi / n_bands * 0.75
sector_angles = np.linspace(0, 2 * np.pi, n_bands, endpoint=False)

# Gray anti-trace fill
theta_fill = np.linspace(0, 2 * np.pi, 200)
ax.fill(theta_fill, np.full_like(theta_fill, 0), color="#f0f0f0", zorder=0)
ax.plot(theta_fill, np.zeros_like(theta_fill), "k--", lw=1.0, alpha=0.5,
        zorder=2)

ax.set_rlim(-vmax, vmax)
rticks = np.arange(-0.4, 0.5, 0.1)
rticks = rticks[(rticks >= -vmax - 0.01) & (rticks <= vmax + 0.01)]
ax.set_rticks(rticks)
ax.set_yticklabels([f"{v:.1f}" if v != 0 else "0" for v in rticks],
                    fontsize=7, color="#888888")
ax.set_rlabel_position(30)
ax.grid(color="#dddddd", linewidth=0.5)
ax.spines["polar"].set_visible(False)

# Sub-bar width per patient within each sector
sub_w = sector_width / len(PATIENTS)

for ip, pat in enumerate(PATIENTS):
    for ib in range(n_bands):
        theta = sector_angles[ib] - sector_width / 2 + sub_w * ip + sub_w / 2
        ax.bar(theta, data[ip, ib], width=sub_w * 0.85, bottom=0,
               color=PAT_COLORS[pat], alpha=0.7, edgecolor="white",
               linewidth=0.3, zorder=3)
    ax.bar(0, 0, color=PAT_COLORS[pat], label=PAT_SHORT[pat])

# Mean ring — line connecting mean values
mean_closed = mean_vals.tolist() + [mean_vals[0]]
angles_closed_v3 = sector_angles.tolist() + [sector_angles[0]]
ax.plot(angles_closed_v3, mean_closed, "k-s", lw=2, markersize=4,
        zorder=5, label="Mean")

# Band labels
ax.set_xticks(sector_angles)
ax.set_xticklabels(BAND_LABELS, fontsize=14, fontweight="bold")

# Alpha annotation
alpha_ib = BANDS.index("alpha")
ax.annotate("5/5\n$p$=0.031",
            xy=(sector_angles[alpha_ib], mean_vals[alpha_ib]),
            xytext=(sector_angles[alpha_ib], vmax * 0.88),
            fontsize=8, fontweight="bold", ha="center", color="#b8860b",
            arrowprops=dict(arrowstyle="-", color="#b8860b", lw=0.8),
            zorder=6)

leg = ax.legend(loc="upper right", bbox_to_anchor=(1.28, 1.08),
                fontsize=9, frameon=True, framealpha=0.9,
                edgecolor="#cccccc")

fig.savefig(OUTDIR / "radar_chart_v3.pdf", bbox_inches="tight", dpi=200)
plt.close(fig)
print("  Saved radar_chart_v3.pdf")


# ═════════════════════════════════════════════════════════════════════
# Markdown companions
# ═════════════════════════════════════════════════════════════════════
for v, desc in [
    ("v1", "Standard radar: one polygon per patient connecting H2a scalar VI contrast across 6 bands. "
           "Gray fill = anti-trace zone (inside zero circle). All 5 patients extend outside zero at alpha."),
    ("v2", "Half-radar: semicircular layout with patient dots along radial axes per band. "
           "Diamond = mean. Gray bars = mean direction. Same data as v1, different geometry."),
    ("v3", "Circular bar chart: each band is a sector, sub-bars per patient. "
           "Black line connects mean values. All alpha bars extend outward."),
]:
    md = (
        f"# radar_chart_{v}\n\n"
        f"## What the figure shows\n{desc}\n\n"
        f"## Data\nH2a Scalar A (mean VI contrast) from scalar_table_full.csv.\n"
        f"Positive = task trace (rest_post closer to task than rest_pre). "
        f"Alpha is the only band where all 5 patients are positive (5/5, p=0.031).\n\n"
        f"## How to read it\nOutside zero circle = positive trace. "
        f"Inside = anti-trace. Alpha consistently outside for all patients.\n"
    )
    (OUTDIR / f"radar_chart_{v}.md").write_text(md)

# Recommendation
rec = (
    "# Recommendation\n\n"
    "**v1 (standard radar)** is the clearest for the paper.\n\n"
    "Reasons:\n"
    "- The polygon shape immediately conveys each patient's profile\n"
    "- The zero circle is the natural baseline — inside/outside is intuitive\n"
    "- Alpha unanimity is visible at a glance (all 5 vertices outside)\n"
    "- The format is familiar to readers (radar/spider charts are standard)\n\n"
    "v2 loses the polygon shape that makes patient profiles readable.\n"
    "v3 is cluttered — 5 sub-bars per sector are hard to parse.\n"
)
(OUTDIR / "recommendation.md").write_text(rec)

print("\nDone!")
