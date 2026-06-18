#!/usr/bin/env python3
"""audit_108 — brain-state TRAJECTORY visualization of the trace (all bands).

The catchy figure the taxonomy was really for. Each phase (rest_pre, task_test,
rest_post) is a state; the distance between two states is d = 1 - ρ^coph (Spearman
of their LRG cophenetic vectors). The three states form a TRIANGLE that embeds
exactly in 2D. Canonical frame: RPre at the origin, Task at (1,0) (÷ L=d(RPre,Task)).
RPost lands at (x,y), and its position is the behaviour:

    Task corner (1,0) = TRACE     RPre corner (0,0) = RESET/ANCHOR     up = REORGANIZE

The canonical **trace** is T_d = d(RPre,Task) − d(Task,RPost) > 0 (locked sign
convention): RPost is closer to Task than the brain's own pre→task move — i.e.
RPost lands inside the **trace zone**, the disk of radius L around Task (which
touches RPre). The trajectory RPre→Task→RPost is the mind-state path; an
*incomplete return* (RPost stuck toward Task) is the trace.

Honest note: the triangles are "fat" (each phase differs from each other by a
similar amount), so the trace is a SUBTLE, consistent pull — not an edge-collapse.
It is real and cohort-coherent at β (8/10 T_d>0, mean T_d +0.09, the largest of
any band; corr +0.47 with the verified ρ_split). We show it as the cohort-mean
trajectory with a bootstrap CI, with per-patient RPost colored by T_d sign.

Outputs (data/audit/cross_phase_taxonomy/):
  figures/state_trajectories.pdf       2×3 band grid (taxonomy plane + mean path)
  figures/state_trajectories_3d.html   Plotly 3D: mean trajectories stacked by band
  state_trajectories.csv               per (patient, band): triangle d's, x, y, T_d

Reuses audit_63 FC→cophenet pipeline verbatim. NO tree-cutting.
Usage: python audit_108_state_trajectory.py
"""
from __future__ import annotations

import importlib.util

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse
from matplotlib.lines import Line2D
from scipy.stats import rankdata, spearmanr

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
_spec = importlib.util.spec_from_file_location(
    "audit_63", ROOT / "scripts" / "01_compute" / "audit"
    / "audit_63_split_baseline_surrogate.py")
a63 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a63)

OUT = ROOT / "data" / "audit" / "cross_phase_taxonomy"
FIG = OUT / "figures"
BANDS = list(BRAIN_BANDS)
C_TRACE, C_RESET = "#1f78b4", "#e08214"     # T_d>0 / T_d<0
NBOOT = 2000
RNG = np.random.default_rng(20260612)


def d_rho(Da, Db):
    return 1.0 - float(spearmanr(Da, Db).statistic)


def embed(d_pt, d_tr, d_pr):
    """Canonical (normalised) RPost coords: RPre=(0,0), Task=(1,0). Returns (x,y)."""
    L = max(d_pt, 1e-9)
    X = (L * L + d_pr * d_pr - d_tr * d_tr) / (2.0 * L)
    Y = np.sqrt(max(0.0, d_pr * d_pr - X * X))
    return X / L, Y / L


def compute() -> pd.DataFrame:
    rows = []
    for pat in a63.COHORT:
        a63.ensure_half_fcs(pat, BANDS)
        for band in BANDS:
            try:
                D = {ph: a63.lrg_ultrametric_condensed(a63.load_phase_fc(pat, ph, band))
                     for ph in ("rest_pre", "task_test", "rest_post",
                                "rest_pre_A", "rest_pre_B")}
            except Exception as e:
                print(f"[skip] {pat}/{band}: {e}")
                continue
            d_pt = d_rho(D["rest_pre"], D["task_test"])
            d_tr = d_rho(D["task_test"], D["rest_post"])
            d_pr = d_rho(D["rest_pre"], D["rest_post"])
            base = d_rho(D["rest_pre_A"], D["rest_pre_B"])
            x, y = embed(d_pt, d_tr, d_pr)
            rows.append(dict(patient=pat, band=band, d_pre_task=d_pt,
                             d_task_post=d_tr, d_pre_post=d_pr, baseline=base,
                             x=x, y=y, T_d=d_pt - d_tr))
    return pd.DataFrame(rows)


def cohort_mean_post(sub, boot=False):
    """Embed the cohort-mean triangle -> mean RPost (x,y). boot=resample patients."""
    if boot:
        idx = RNG.integers(0, len(sub), len(sub))
        sub = sub.iloc[idx]
    return embed(sub.d_pre_task.mean(), sub.d_task_post.mean(), sub.d_pre_post.mean())


def fig_grid(df):
    use_lrg_style()
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.6), sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, band in zip(axes, BANDS):
        sub = df[df.band == band]
        # trace zone (T_d>0): RPost within distance L of Task = disk radius 1 @ Task.
        # outline only, so blue (inside) vs orange (outside) dots read cleanly.
        ax.add_patch(Circle((1, 0), 1.0, fill=False, ec=C_TRACE, lw=1.0,
                            ls=(0, (5, 4)), alpha=0.7, zorder=0))
        ax.text(0.93, 1.06, "trace zone", color=C_TRACE, fontsize=7, alpha=0.85,
                ha="right")
        # per-patient RPost + faint trajectory (blue inside zone = trace, T_d>0)
        for _, r in sub.iterrows():
            col = C_TRACE if r.T_d > 0 else C_RESET
            ax.plot([0, 1, r.x], [0, 0, r.y], color=col, lw=0.5, alpha=0.22, zorder=1)
            ax.scatter(r.x, r.y, s=28, color=col, ec="white", lw=0.5, zorder=3)
        # cohort-mean trajectory + bootstrap CI ellipse
        mx, my = cohort_mean_post(sub)
        boots = np.array([cohort_mean_post(sub, boot=True) for _ in range(NBOOT)])
        cov = np.cov(boots.T)
        vals, vecs = np.linalg.eigh(cov)
        ang = np.degrees(np.arctan2(vecs[1, np.argmax(vals)], vecs[0, np.argmax(vals)]))
        w, h = 2 * 2 * np.sqrt(np.maximum(vals, 0))   # 2σ
        ax.add_patch(Ellipse((mx, my), w, h, angle=ang, color="#222",
                             alpha=0.18, zorder=4))
        ax.plot([0, 1, mx], [0, 0, my], color="#222", lw=2.2, zorder=5)
        ax.scatter(mx, my, s=230, marker="*", color="#222", ec="white", lw=1.2,
                   zorder=6)
        ax.scatter([0, 1], [0, 0], s=95, marker="s", color="black", zorder=6)
        ax.annotate("", xy=(1, 0), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color="0.4", lw=1.0), zorder=2)
        ax.text(0.0, -0.13, "RPre", ha="center", fontsize=8, fontweight="bold")
        ax.text(1.0, -0.13, "Task", ha="center", fontsize=8, fontweight="bold")
        ntr = int((sub.T_d > 0).sum())
        ax.set_title(f"{BRAIN_BAND_TEX_DICT.get(band, band)}    "
                     f"trace {ntr}/10   $\\overline{{T_d}}$={sub.T_d.mean():+.3f}",
                     fontsize=10.5)
        ax.set_xlim(-0.45, 1.55)
        ax.set_ylim(-0.22, 1.5)
        ax.set_aspect("equal")
    for ax in axes[3:]:
        ax.set_xlabel("RPre (reset) ← progress → Task (trace)")
    for ax in (axes[0], axes[3]):
        ax.set_ylabel("reorganize (excursion)")
    handles = [Line2D([], [], marker="o", ls="", mec="white", color=C_TRACE,
                      label="RPost, $T_d>0$ (trace)"),
               Line2D([], [], marker="o", ls="", mec="white", color=C_RESET,
                      label="RPost, $T_d<0$ (reset)"),
               Line2D([], [], color="#222", lw=2.2, marker="*", mec="white",
                      label="cohort-mean path + 2σ CI"),
               Line2D([], [], color=C_TRACE, ls=(0, (5, 4)), lw=1.0,
                      label="trace zone (RPost closer to Task than RPre→Task)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "state_trajectories.pdf", transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIG/'state_trajectories.pdf'}")


def fig_3d(df):
    try:
        import plotly.graph_objects as go
    except Exception as e:
        print(f"[3d skip] plotly unavailable: {e}")
        return
    figg = go.Figure()
    for zi, band in enumerate(BANDS):
        sub = df[df.band == band]
        figg.add_trace(go.Scatter3d(
            x=[0, 1], y=[0, 0], z=[zi, zi], mode="lines+markers+text",
            line=dict(color="black", width=5), marker=dict(size=4, color="black"),
            text=["RPre", "Task"], textposition="top center",
            showlegend=False, hoverinfo="skip"))
        for _, r in sub.iterrows():
            col = C_TRACE if r.T_d > 0 else C_RESET
            figg.add_trace(go.Scatter3d(
                x=[r.x], y=[r.y], z=[zi], mode="markers",
                marker=dict(size=4, color=col, line=dict(color="white", width=1)),
                showlegend=False, hoverinfo="text",
                hovertext=f"{r.patient} {band} T_d={r.T_d:+.3f}"))
        mx, my = cohort_mean_post(sub)
        figg.add_trace(go.Scatter3d(
            x=[0, 1, mx], y=[0, 0, my], z=[zi, zi, zi], mode="lines+markers",
            line=dict(color="#d62728" if band == "beta" else "#444",
                      width=7 if band == "beta" else 4),
            marker=dict(size=[3, 3, 9], color="#d62728" if band == "beta" else "#444"),
            name=f"{band} mean", showlegend=False,
            hovertext=f"{band} mean RPost", hoverinfo="text"))
    figg.update_layout(
        scene=dict(xaxis_title="RPre(reset) → Task(trace)", yaxis_title="reorganize",
                   zaxis=dict(title="band", tickvals=list(range(len(BANDS))),
                              ticktext=[BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS]),
                   aspectmode="manual", aspectratio=dict(x=1.3, y=1.3, z=2.4)),
        title="Mind-state trajectories RPre→Task→RPost across frequency "
              "(β mean path in red; RPost dots blue=trace)",
        margin=dict(l=0, r=0, t=42, b=0))
    out = FIG / "state_trajectories_3d.html"
    figg.write_html(str(out))
    print(f"wrote {out}")


def main():
    df = compute()
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "state_trajectories.csv", index=False)
    print("\nband        n(T_d>0)/10   mean T_d    cohort-mean RPost (x,y)")
    for band in BANDS:
        sub = df[df.band == band]
        mx, my = cohort_mean_post(sub)
        print(f"{band:11s}    {int((sub.T_d>0).sum()):2d}        {sub.T_d.mean():+.3f}"
              f"      ({mx:.2f}, {my:.2f})")
    fig_grid(df)
    fig_3d(df)


if __name__ == "__main__":
    main()
