#!/usr/bin/env python3
"""H2-TOPO — directed topology-level trace tests.

Uses the bipartition-overlap data from `diag_topology_vs_rho.csv` to
run the DIRECTED H2-family tests at the tree-topology level:

  H2a-topo:  bip_overlap(Z^rpost, Z^ttest)  >  bip_overlap(Z^rpost, Z^rpre)
             — "is rpost's tree more topologically similar to task's
             than to rest_pre's?"
  H2b-topo:  bip_overlap(Z^rpre, Z^ttest)   >  bip_overlap(Z^ttest, Z^rpost)
             — approach > exit at topology level.
  H1-topo:   bip_overlap(Z^tlearn, Z^ttest) − mean cross-type overlap
  H3-topo:   mean bip_overlap(within-type)  −  mean bip_overlap(cross-type)

Per-band paired Wilcoxon across patients, FDR-BH m=6, 10k bootstrap CI,
rank-biserial. Plus per-patient × per-band contrast heatmap and
band-unanimity count (how many of 9 patients have positive contrast
in each band).

Reads:  data/reports/imcoh_vi/diag_topology_vs_rho.csv
Writes: data/reports/imcoh_vi/h2_topology_directed.{md,csv}
        data/reports/imcoh_vi/figures/h2_topology_directed.{pdf,png}
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.config.paths import REPORTS_ROOT

from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, rank_biserial, boot_ci_mean, bh_fdr


PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")
TYPE = {"rest_pre": "rest", "rest_post": "rest",
        "task_learn": "task", "task_test": "task"}


def _lookup(df: pd.DataFrame, pat: str, band: str, p1: str, p2: str) -> float:
    base = df[(df["patient"] == pat) & (df["band"] == band)]
    hit = base[(base["p1"] == p1) & (base["p2"] == p2)]
    if len(hit) == 0:
        hit = base[(base["p1"] == p2) & (base["p2"] == p1)]
    return float(hit["bip_overlap"].iloc[0]) if len(hit) else np.nan


def _h2a(df, pat, band):
    return (_lookup(df, pat, band, "rest_post", "task_test")
            - _lookup(df, pat, band, "rest_post", "rest_pre"))


def _h2b(df, pat, band):
    return (_lookup(df, pat, band, "rest_pre", "task_test")
            - _lookup(df, pat, band, "rest_post", "task_test"))


def _h1(df, pat, band):
    r_tlt = _lookup(df, pat, band, "task_learn", "task_test")
    crosses = []
    for p1 in PHASES:
        for p2 in PHASES:
            if p1 < p2 and TYPE[p1] != TYPE[p2]:
                v = _lookup(df, pat, band, p1, p2)
                if np.isfinite(v):
                    crosses.append(v)
    return (r_tlt - float(np.mean(crosses))) if crosses and np.isfinite(r_tlt) else np.nan


def _h3(df, pat, band):
    within, cross = [], []
    for p1 in PHASES:
        for p2 in PHASES:
            if p1 >= p2:
                continue
            v = _lookup(df, pat, band, p1, p2)
            if not np.isfinite(v):
                continue
            (within if TYPE[p1] == TYPE[p2] else cross).append(v)
    if not within or not cross:
        return np.nan
    return float(np.mean(within) - np.mean(cross))


def per_band_test(df, fn, patients):
    out = []
    for band in BRAIN_BANDS_NAMES:
        vals = np.array([fn(df, pat, band) for pat in patients], dtype=float)
        pats = [p for p, v in zip(patients, vals) if np.isfinite(v)]
        vals = vals[np.isfinite(vals)]
        if len(vals) < 3:
            out.append({"band": band, "n": len(vals)})
            continue
        z, p = wilcoxon_z(vals)
        r_rb = rank_biserial(vals)
        mean, lo, hi = boot_ci_mean(vals)
        n_pos = int((vals > 0).sum())
        out.append({"band": band, "n": len(vals),
                    "mean": mean, "ci_lo": lo, "ci_hi": hi,
                    "r_rb": r_rb, "z": z, "p": p,
                    "n_pos": n_pos, "values": vals.tolist(),
                    "patients": pats})
    ps = [r.get("p", np.nan) for r in out]
    finite = [np.isfinite(pp) for pp in ps]
    qs = bh_fdr([pp for pp, ok in zip(ps, finite) if ok])
    j = 0
    for r, ok in zip(out, finite):
        r["q"] = qs[j] if ok else np.nan
        if ok:
            j += 1
    return out


def _table(rows, tex, metric_name):
    L = [f"| band | n | mean {metric_name} | 95 % CI | pos | r_rb | p | q (BH) |",
         "|------|--:|--------:|:-------|----:|-----:|--:|-------:|"]
    for r in rows:
        tx = tex[r["band"]]
        if "mean" not in r:
            L.append(f"| {tx} | {r['n']} | — | — | — | — | — | — |")
            continue
        star = "★" if r["q"] < 0.05 else ""
        L.append(
            f"| {tx} | {r['n']} | {r['mean']:+.4f} | "
            f"[{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}] | {r['n_pos']}/{r['n']} | "
            f"{r['r_rb']:+.3f} | {r['p']:.4f} | {r['q']:.4f}{star} |"
        )
    return L


# ─────────────────────────── figure ───────────────────────────

BAND_COLORS = {
    "delta":      "#3b6e9c", "theta":      "#c44e4e",
    "alpha":      "#5ea85e", "beta":       "#b28ad1",
    "low_gamma":  "#e5a24b", "high_gamma": "#7d5a50",
}


def figure(df, h1, h2a, h2b, h3):
    patients = list(PATIENTS_LIST)
    fig = plt.figure(figsize=(13.0, 8.0), dpi=160)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.1],
                           height_ratios=[1.0, 1.15],
                           hspace=0.45, wspace=0.35)

    # (a) — Band × patient heatmap of H2a-topo contrast
    ax = fig.add_subplot(gs[0, 0])
    mat = np.full((len(patients), len(BRAIN_BANDS_NAMES)), np.nan)
    for ip, pat in enumerate(patients):
        for ib, band in enumerate(BRAIN_BANDS_NAMES):
            mat[ip, ib] = _h2a(df, pat, band)
    vmax = float(np.nanmax(np.abs(mat)))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    im = ax.imshow(mat, cmap="RdYlGn", norm=norm, aspect="auto")
    ax.set_xticks(range(len(BRAIN_BANDS_NAMES)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                       fontsize=11)
    ax.set_yticks(range(len(patients)))
    ax.set_yticklabels(patients, fontsize=9)
    for i, pat in enumerate(patients):
        for j, band in enumerate(BRAIN_BANDS_NAMES):
            v = mat[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.02f}",
                        ha="center", va="center", fontsize=8,
                        color="white" if abs(v) > vmax * 0.6 else "black")
    # Count unanimity per band at top
    pos_counts = np.nansum(mat > 0, axis=0).astype(int)
    tot_counts = np.sum(np.isfinite(mat), axis=0).astype(int)
    for j, (p, t) in enumerate(zip(pos_counts, tot_counts)):
        ax.text(j, -0.8, f"{p}/{t}", ha="center", va="bottom",
                fontsize=9, fontweight="bold")
    ax.set_title(
        "(a) H2a-topo: bip_overlap(rpost,tt) − bip_overlap(rpost,rpre)\n"
        "per patient × band. Top numbers: patients with positive contrast / n",
        fontsize=10, loc="left",
    )
    fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)

    # (b) — Band ranking by H2a-topo mean contrast + H3-topo Δ
    ax = fig.add_subplot(gs[0, 1])
    h2a_map = {r["band"]: r for r in h2a if "mean" in r}
    h3_map  = {r["band"]: r for r in h3  if "mean" in r}
    x = np.arange(len(BRAIN_BANDS_NAMES))
    w = 0.38
    means_h2a = [h2a_map[b]["mean"] for b in BRAIN_BANDS_NAMES]
    means_h3  = [h3_map[b]["mean"] for b in BRAIN_BANDS_NAMES]
    ax.bar(x - w/2, means_h2a, w,
           color=[BAND_COLORS[b] for b in BRAIN_BANDS_NAMES],
           edgecolor="black", lw=0.6, label="H2a-topo (trace)")
    ax.bar(x + w/2, means_h3, w,
           color=[BAND_COLORS[b] for b in BRAIN_BANDS_NAMES],
           edgecolor="black", lw=0.6, alpha=0.55, label="H3-topo (within−cross)")
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                       fontsize=11)
    ax.set_ylabel("mean contrast on bipartition overlap", fontsize=10)
    for xi, r in zip(x, [h2a_map[b] for b in BRAIN_BANDS_NAMES]):
        star = "★" if np.isfinite(r.get("q", np.nan)) and r["q"] < 0.05 else ""
        ax.text(xi - w/2, r["mean"] + 0.003, star,
                ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.set_title("(b) Band ranking on H2a-topo and H3-topo",
                  fontsize=11, loc="left")
    ax.grid(axis="y", alpha=0.3)

    # (c) — Per-patient H2a-topo value lines
    ax = fig.add_subplot(gs[1, :])
    for band in BRAIN_BANDS_NAMES:
        vals = [_h2a(df, pat, band) for pat in patients]
        ax.plot(range(len(patients)), vals, marker="o", lw=1.3,
                color=BAND_COLORS[band], label=BRAIN_BAND_TEX_DICT[band],
                alpha=0.85)
    ax.axhline(0, color="#888", lw=0.8, linestyle="--")
    ax.set_xticks(range(len(patients)))
    ax.set_xticklabels(patients, fontsize=9)
    ax.set_ylabel("H2a-topo contrast per patient", fontsize=10)
    ax.set_title(
        "(c) Per-patient H2a-topo contrast per band. "
        "Bands that lie consistently above zero carry a topology-level trace.",
        fontsize=10, loc="left",
    )
    ax.legend(ncol=6, loc="lower right", frameon=False, fontsize=9)
    ax.grid(alpha=0.3)

    fig.suptitle(
        "H2-TOPO: directed trace tests on tree bipartition overlap (n=9)",
        fontsize=12, y=1.00,
    )
    out = REPORTS_ROOT / "imcoh_vi" / "figures" / "h2_topology_directed.png"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    print(f"saved {out}")
    plt.close(fig)


def main() -> None:
    df = pd.read_csv(REPORTS_ROOT / "imcoh_vi" / "diag_topology_vs_rho.csv")
    patients = list(PATIENTS_LIST)

    h1  = per_band_test(df, _h1,  patients)
    h2a = per_band_test(df, _h2a, patients)
    h2b = per_band_test(df, _h2b, patients)
    h3  = per_band_test(df, _h3,  patients)

    # CSV
    records = []
    for tag, rows in [("H1_topo", h1), ("H2a_topo", h2a),
                       ("H2b_topo", h2b), ("H3_topo", h3)]:
        for r in rows:
            if "mean" not in r:
                continue
            records.append({"hypothesis": tag, **{k: r.get(k) for k in
                ("band", "n", "n_pos", "mean", "ci_lo", "ci_hi",
                 "r_rb", "z", "p", "q")}})
    csv = REPORTS_ROOT / "imcoh_vi" / "h2_topology_directed.csv"
    pd.DataFrame(records).to_csv(csv, index=False)
    print(f"wrote {csv}")

    lines: list[str] = []
    ap = lines.append
    ap("# H2-TOPO — directed topology-level trace tests (n=9)")
    ap("")
    ap("Metric per (patient, band, phase-pair): fraction of internal tree "
       "bipartitions shared between `Z^p1` and `Z^p2` (exact set-of-leaves "
       "match, out of N−1 internal nodes). Multiscale by construction — "
       "bipartitions span every merge height.")
    ap("")
    ap("## H1-topo — task phases share more bipartitions than cross-type")
    lines.extend(_table(h1, BRAIN_BAND_TEX_DICT, "Δ"))
    ap("")
    ap("## H2a-topo — rpost's tree more similar to task's than to rpre's")
    ap("`bip_overlap(rpost, tt) − bip_overlap(rpost, rpre)`")
    lines.extend(_table(h2a, BRAIN_BAND_TEX_DICT, "Δ"))
    ap("")
    ap("## H2b-topo — approach > exit")
    lines.extend(_table(h2b, BRAIN_BAND_TEX_DICT, "Δ"))
    ap("")
    ap("## H3-topo — within-type > cross-type")
    lines.extend(_table(h3, BRAIN_BAND_TEX_DICT, "Δ"))
    ap("")

    md = REPORTS_ROOT / "imcoh_vi" / "h2_topology_directed.md"
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {md}")

    # Console summary
    for tag, rows in [("H1_topo", h1), ("H2a_topo", h2a),
                       ("H2b_topo", h2b), ("H3_topo", h3)]:
        n_pass = sum(1 for r in rows
                     if np.isfinite(r.get("q", np.nan)) and r["q"] < 0.05)
        total_pos = sum(r.get("n_pos", 0) for r in rows if "n_pos" in r)
        total_n   = sum(r.get("n", 0) for r in rows if "n_pos" in r)
        print(f"  {tag}: {n_pass}/6 bands pass FDR; "
              f"{total_pos}/{total_n} patient×band cells positive")

    figure(df, h1, h2a, h2b, h3)


if __name__ == "__main__":
    main()
