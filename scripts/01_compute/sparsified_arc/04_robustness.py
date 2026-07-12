#!/usr/bin/env python3
"""D2 robustness -- is each band's tau-trace robust, or backbone/patient-specific?

Reads the per_cell.csv from every trace-arc backbone that exists (percolation
primary; mst-union d0.15 / d0.20 cross-checks) and, per band, reports the
matched-strength scale-max gate under three stress tests the user demanded
(never trust a single p):

  1. cohort scale-max gate p            (Wilcoxon obs_max - surr_max_p50)
  2. LOO robustness: max gate p over all 10 leave-one-patient-out refits
     (a band that dies when one patient is dropped is single-patient-driven)
  3. per-patient consensus: how many patients are positive (obs_max > null) and
     how many clear their own cell gate (p_scalemax < 0.05)
  4. cross-backbone agreement: a band is ROBUST only if it clears in EVERY backbone.

Discriminator for the "continuous-multiscale" suspects (delta/gamma_high): a
genuine multiscale trace clears at a coherent scale in EVERY backbone and
survives LOO; a density-dependent null artifact clears only under percolation
(heterogeneous density) and/or collapses under LOO.

Outputs:
    data/sparsified_arc/trace_arc/robustness_table.csv
    data/sparsified_arc/trace_arc/figs/robustness_heatmap.pdf
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import BANDS
from lrg_eegfc.visuals.styles import use_lrg_style

SA = ROOT / "data" / "sparsified_arc"
BACKBONES = {"percolation": SA / "trace_arc",
             "d0.15": SA / "trace_arc_d015",
             "d0.20": SA / "trace_arc_d020"}


def _gate_p(d):
    d = d[np.isfinite(d)]
    if d.size < 5 or not np.any(d != 0):
        return np.nan
    try:
        return float(wilcoxon(d, alternative="greater")[1])
    except Exception:
        return np.nan


def _loo_max_p(d):
    d = d[np.isfinite(d)]
    if d.size < 6:
        return np.nan
    ps = [_gate_p(np.delete(d, i)) for i in range(d.size)]
    ps = [p for p in ps if np.isfinite(p)]
    return float(max(ps)) if ps else np.nan


def analyse():
    rows = []
    for name, path in BACKBONES.items():
        f = path / "per_cell.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f)
        for band in BANDS:
            x = df[df.band == band]
            if x.empty:
                continue
            d = (x.obs_max - x.surr_max_p50).values
            rows.append(dict(
                backbone=name, band=band, n=len(x),
                gate_p=_gate_p(d), loo_max_p=_loo_max_p(d),
                n_pos_sign=int((x.obs_max > x.surr_max_p50).sum()),
                n_cell_clear=int((x.p_scalemax < 0.05).sum()),
                argmax_s_med=float(x.argmax_s.median()),
                obs_max_med=float(x.obs_max.median())))
    return pd.DataFrame(rows)


def verdict(tab):
    """Robust = clears (gate<0.05) AND LOO-robust in EVERY available backbone."""
    out = []
    bbs = tab.backbone.unique()
    for band in BANDS:
        t = tab[tab.band == band]
        if t.empty:
            continue
        clears = {r.backbone: (r.gate_p < 0.05) for r in t.itertuples()}
        loo_ok = {r.backbone: (np.isfinite(r.loo_max_p) and r.loo_max_p < 0.05) for r in t.itertuples()}
        all_clear = all(clears.get(b, False) for b in bbs)
        all_loo = all(loo_ok.get(b, False) for b in bbs)
        n_clear = sum(clears.values())
        v = ("ROBUST" if all_clear and all_loo else
             "robust-gate-not-LOO" if all_clear else
             f"partial({n_clear}/{len(bbs)})" if n_clear else "none")
        out.append(dict(band=band, verdict=v,
                        clears_in="/".join(b for b in bbs if clears.get(b)),
                        **{f"gate_{b}": float(t[t.backbone == b].gate_p.iloc[0])
                           for b in bbs if (t.backbone == b).any()},
                        **{f"loo_{b}": float(t[t.backbone == b].loo_max_p.iloc[0])
                           for b in bbs if (t.backbone == b).any()}))
    return pd.DataFrame(out)


def heatmap(tab):
    use_lrg_style()
    bbs = list(tab.backbone.unique())
    M = np.full((len(BANDS), len(bbs) * 2), np.nan)
    cols = []
    for bi, b in enumerate(bbs):
        cols += [f"{b}\ngate", f"{b}\nLOO"]
    for i, band in enumerate(BANDS):
        for bi, b in enumerate(bbs):
            r = tab[(tab.band == band) & (tab.backbone == b)]
            if r.empty:
                continue
            M[i, 2 * bi] = -np.log10(max(r.gate_p.iloc[0], 1e-4))
            M[i, 2 * bi + 1] = -np.log10(max(r.loo_max_p.iloc[0], 1e-4)) if np.isfinite(r.loo_max_p.iloc[0]) else np.nan
    fig, ax = plt.subplots(figsize=(1.4 * len(cols) + 2, 4.2))
    im = ax.imshow(M, aspect="auto", cmap="viridis", vmin=0, vmax=3)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, fontsize=8)
    ax.set_yticks(range(len(BANDS))); ax.set_yticklabels(BANDS)
    for i in range(len(BANDS)):
        for j in range(len(cols)):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f"{10**(-M[i,j]):.2f}", ha="center", va="center",
                        color="w" if M[i, j] < 1.5 else "k", fontsize=7)
    cb = fig.colorbar(im, ax=ax, label=r"$-\log_{10}p$")
    cb.ax.axhline(-np.log10(0.05), color="r", lw=1.5)
    ax.set_title("τ-trace robustness: scale-max gate & LOO across backbones\n(cell text = p; p<0.05 bright)",
                 fontweight="bold", fontsize=10)
    fig.tight_layout()
    (SA / "figures" / "robustness").mkdir(parents=True, exist_ok=True)
    out = SA / "figures" / "robustness" / "robustness_heatmap.pdf"
    fig.savefig(out, transparent=True); plt.close(fig); print(f"[fig] {out}")


def main():
    tab = analyse()
    tab.to_csv(SA / "trace_arc" / "robustness_table.csv", index=False)
    vt = verdict(tab)
    vt.to_csv(SA / "trace_arc" / "robustness_verdict.csv", index=False)
    print("=== PER-BACKBONE per-band scale-max gate + LOO ===")
    print(tab[["backbone", "band", "gate_p", "loo_max_p", "n_pos_sign",
               "n_cell_clear", "argmax_s_med"]].to_string(index=False))
    print("\n=== ROBUSTNESS VERDICT (robust = clears + LOO in every backbone) ===")
    print(vt.to_string(index=False))
    heatmap(tab)
    print(f"[D2-robustness] -> {SA/'trace_arc'}")


if __name__ == "__main__":
    main()
