#!/usr/bin/env python3
r"""Per-pair cophenetic trace in 3D: density on the floor, MEANINGFULNESS as height.

Rendered with Plotly (not matplotlib): matplotlib's mplot3d draws whole artists
front-to-back and cannot clip a plane against a surface, so the significance
sheet floats wrong and the floor can end up above the terrain. Plotly composites
in depth, so the significance plane correctly cuts the terrain and the figure is
rotatable.

The all-bands movers density (preprint_24) shows WHERE the per-pair co-movement
sits — but it is a pooled magnitude, blind to the patient index, so it cannot
tell a cohort-coherent trace (β: 8/10) from a two-patient diagonal (γ_l: 6/4).
The consistency t-map (preprint_25) is that missing dimension: per cell, the
cohort t of (obs movers density − the patient's OWN matched-strength surrogate)
— cross-patient AGREEMENT on the per-pair excess over the strength SURROGATE.

This figure stacks them, one 3-D scene per band, on the movers' within-patient
rank square (u = rank Δ^coph_task, v = rank Δ^coph_rest):

  * FLOOR (the 2-D KDE)        — the pooled movers-density excess (density − 1).
  * HEIGHT + COLOUR (meaning)  — the composite meaningfulness
    M(u,v) = ρ_band · t(u,v) · max(0, 2·frac_pos − 1), fusing the three factors
    behind the verdict: the gate EFFECT MAGNITUDE ρ_band (β 0.22 ≫ α 0.10),
    the per-cell cohort CLEARANCE t over the matched-strength surrogate
    (mean ÷ between-patient SE), and the cohort AGREEMENT (sign consensus).
    Why all three: a plain t divides magnitude away (α, being tighter, then
    out-ranks β) and is lifted by 2–3 strong patients in a split cohort (γ_l,
    carried by Pat_05/02, scored like β's 8/10). Re-injecting ρ restores β as
    the tallest; the agreement factor zeroes split cells so γ_l collapses.
    Diagonal M: β 0.186 > α 0.111 > γ_l 0.033 > γ_h/θ/δ ≈ 0.

Read across: β is the tallest, most diagonal-coherent ridge; α second; γ_l is
clearly demoted; θ/δ/γ_h are flat. This is a descriptive meaningfulness score;
the formal per-band verdict is the paired-Wilcoxon gate / forest (preprint_23).

Reuses preprint_24 (pooled_excess) + preprint_25 (movers/surrogate density
helpers), caching the per-patient excess fields E once at
data/cache/perpair_consistency_E/. Cophenetic only (the surrogate is the gate's).
Outputs (data/preprint/figures/all_bands/):
    fig_bands_perpair_meaningfulness_3d_coph.html   (interactive — rotate to read)
    fig_bands_perpair_meaningfulness_3d_coph.png    (static, via kaleido)
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.ndimage import gaussian_filter, zoom

from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
sys.path.insert(0, str(Path(__file__).resolve().parent))
p18 = importlib.import_module("preprint_18_bands_joint_density_rawfc")
p24 = importlib.import_module("preprint_24_bands_perpair_movers_density_coph")
p25 = importlib.import_module("preprint_25_bands_perpair_consistency_tmap_coph")

BAND_ORDER = p18.BAND_ORDER
ZOOM, SMOOTH = 10, 2.2              # 6x6 grid → smooth 60x60 terrain
FLOOR_DZ = 0.45                    # movers-density floor sits this far below z=0
ECACHE = CACHE_ROOT / "perpair_consistency_E"
# Band effect magnitude = the gate's cohort-median per-patient ρ (the trusted
# headline effect size; β 0.22 ≫ α 0.10). The coarse movers density washes this
# out and the t-statistic divides it away, so we re-inject it as a band factor.
GATE_RHO = (pd.read_csv(p18.LAYERS["coph"]["cohort_csv"])
            .set_index("band")["obs_median_rho"].to_dict())

# Plotly renders no LaTeX in subplot titles — use unicode band symbols.
BAND_UNICODE = {"delta": "δ", "theta": "θ", "alpha": "α", "beta": "β",
                "low_gamma": "γₗ", "high_gamma": "γₕ"}


def _plotly_scale(cmap, n=33):
    return [[i / (n - 1),
             "rgb(%d,%d,%d)" % tuple(int(255 * c) for c in cmap(i / (n - 1))[:3])]
            for i in range(n)]


TRACE_SCALE = _plotly_scale(p18.TRACE_CMAP)
FLOOR_SCALE = _plotly_scale(p18.TRACE_CMAP)        # same theme, muted via opacity


def _band_E(band: str) -> np.ndarray:
    """Per-patient excess fields E[p, iu, iv] = obs movers density − the patient's
    OWN matched-strength surrogate-median movers density. Cached once (the
    surrogate reconstruction is the only expensive step)."""
    cache = ECACHE / f"{band}_G{p25.GRID}_NS{p25.N_SURR}_seed{p25.SEED}.npy"
    if cache.exists():
        return np.load(cache)
    rng = np.random.default_rng(p25.SEED)
    edges = np.linspace(0.0, 1.0, p25.GRID + 1)
    pairs = p18.coph_obs_pairs(band)
    E = []
    for pat in p18.COHORT:
        dt, dr = pairs[pat]
        obs = p25._density(*p25._movers_uv(dt, dr, rng), edges)
        dts, drs, _ = p18._surrogate_rho_and_pairs(pat, band, "coph")
        surr = [p25._density(*p25._movers_uv(s_dt, s_dr, rng), edges)
                for s_dt, s_dr in list(zip(dts, drs))[:p25.N_SURR]]
        E.append(obs - np.median(np.stack(surr), axis=0))
    E = np.stack(E)
    ECACHE.mkdir(parents=True, exist_ok=True)
    np.save(cache, E)
    return E


def band_meaningfulness(band: str) -> np.ndarray:
    """Composite meaningfulness M[iu, iv] = ρ_band · t · max(0, 2·frac_pos − 1).

    Fuses the three factors behind the verdict, so β (the most significant band)
    is tallest:
      • ρ_band  — the gate effect magnitude (β 0.22 ≫ α 0.10); re-injected
                  because the coarse density washes it out and t divides it away;
      • t       — per-cell cohort clearance over the matched-strength surrogate
                  (mean ÷ between-patient SE);
      • agree   — sign-consensus, which zeroes split cells so a few strong
                  patients cannot carry a band (demotes γ_l, a 6/4 split).
    Diagonal means: β 0.186 > α 0.111 > γ_l 0.033 > γ_h/θ/δ ≈ 0. Negative-ρ
    bands (θ) clip to 0. This is a descriptive meaningfulness score, not a formal
    test — the verdict is the paired-Wilcoxon gate / forest (preprint_23).
    """
    E = _band_E(band)
    n = E.shape[0]
    t = E.mean(0) / (E.std(0, ddof=1) / np.sqrt(n) + 1e-12)
    agree = np.clip(2.0 * np.mean(E > 0, axis=0) - 1.0, 0.0, 1.0)
    rho = max(float(GATE_RHO.get(band, 0.0)), 0.0)
    return rho * t * agree                                  # [iu, iv]


def _m_surface(band: str):
    M = gaussian_filter(zoom(band_meaningfulness(band).T, ZOOM, order=1),
                        sigma=SMOOTH)                       # [iv, iu] → x=u, y=v
    g = (np.arange(M.shape[0]) + 0.5) / M.shape[0]
    return g, g, M


def _floor(band: str, rng):
    excess, _f, _d = p24.pooled_excess(band, rng)           # [v, u] 120x120
    excess = zoom(excess, 0.5, order=1)                     # 60x60 for lighter html
    g = (np.arange(excess.shape[0]) + 0.5) / excess.shape[0]
    return g, g, excess


def main() -> Path:
    rng = np.random.default_rng(p24.SEED)
    fig = make_subplots(
        rows=2, cols=3,
        specs=[[{"type": "surface"}] * 3] * 2,
        subplot_titles=[f"<b>{BAND_UNICODE[b]}</b>" for b in BAND_ORDER],
        horizontal_spacing=0.01, vertical_spacing=0.06,
    )

    fvabs = 0.0
    floors, surfs = {}, {}
    for b in BAND_ORDER:
        surfs[b] = _m_surface(b)
        floors[b] = _floor(b, rng)
        fvabs = max(fvabs, float(np.percentile(np.abs(floors[b][2]), 98.0)))
    mvabs = float(np.percentile(
        np.concatenate([np.abs(M).ravel() for *_g, M in surfs.values()]), 99.0))
    zrange = (-FLOOR_DZ - 0.25, mvabs + 0.3)

    print("Per-pair meaningfulness 3D (Plotly; height = agreement-weighted M):")
    for idx, band in enumerate(BAND_ORDER):
        r, c = idx // 3 + 1, idx % 3 + 1
        gx, gy, M = surfs[band]
        fx, fy, F = floors[band]

        # meaningfulness terrain — only positive M rises (flat plain elsewhere)
        fig.add_trace(go.Surface(
            x=gx, y=gy, z=np.clip(M, 0.0, mvabs), surfacecolor=M,
            colorscale=TRACE_SCALE, cmin=-mvabs, cmax=mvabs, showscale=(idx == 0),
            colorbar=dict(title="M", len=0.5, x=1.005, y=0.78, thickness=12),
            lighting=dict(ambient=0.85, diffuse=0.4, specular=0.05),
            contours={"z": {"show": True, "start": 0.0, "end": mvabs, "size": 0.4,
                            "color": "rgba(40,40,40,0.25)"}},
            name="M"), row=r, col=c)
        # floor: movers density excess (the 2-D KDE)
        fig.add_trace(go.Surface(
            x=fx, y=fy, z=np.full_like(F, -FLOOR_DZ), surfacecolor=F,
            colorscale=FLOOR_SCALE, cmin=-fvabs, cmax=fvabs, showscale=False,
            opacity=0.92, name="density"), row=r, col=c)

        diag = float(np.mean(np.diag(M)))
        print(f"  {band:11s} max M={M.max():+.2f}  diag-mean M={diag:+.2f}")

    scene = dict(
        xaxis=dict(title="u", range=[0, 1], tickvals=[0, 0.5, 1],
                   backgroundcolor="rgba(0,0,0,0)"),
        yaxis=dict(title="v", range=[0, 1], tickvals=[0, 0.5, 1],
                   backgroundcolor="rgba(0,0,0,0)"),
        zaxis=dict(title="meaningfulness M", range=list(zrange),
                   backgroundcolor="rgba(0,0,0,0)"),
        aspectratio=dict(x=1, y=1, z=1.05),
        camera=dict(eye=dict(x=1.55, y=-1.55, z=0.85)),
    )
    fig.update_layout(
        **{f"scene{'' if i == 0 else i + 1}": scene for i in range(6)},
        height=900, width=1500, margin=dict(l=0, r=40, t=46, b=0),
        paper_bgcolor="white",
        title_text=("Per-pair cophenetic trace · floor = movers density (the 2-D "
                    "KDE) · height M = effect ρ × surrogate clearance × patient "
                    "agreement (β tallest, then α; γ_l demoted)"),
        title_x=0.5, title_font_size=13,
    )

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_html = out_dir / "fig_bands_perpair_meaningfulness_3d_coph.html"
    fig.write_html(str(out_html), include_plotlyjs="cdn")
    print(f"  Saved (interactive): {out_html}")
    out_png = out_dir / "fig_bands_perpair_meaningfulness_3d_coph.png"
    try:
        fig.write_image(str(out_png), scale=2)
        print(f"  Saved (static): {out_png}")
    except Exception as e:  # noqa: BLE001
        print(f"  [kaleido] static export skipped: {e}")
    return out_html


if __name__ == "__main__":
    main()
