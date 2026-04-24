#!/usr/bin/env python3
"""H2-PART — multiscale partition-level cross-phase comparison.

For every (patient, band, k ∈ [2, 49]) build the partitions
P^rpre_k, P^tt_k, P^rpost_k by cutting the LRG dendrogram at `k`
clusters, and compute four DIRECTED partition-level contrasts:

  Δ_VI(k)   = VI(P^rpre, P^rpost)   − VI(P^tt, P^rpost)
              ("rpost closer to tt than rpre at partition level k")

  Δ_ARI(k)  = ARI(P^tt, P^rpost)    − ARI(P^rpre, P^rpost)
              (ARI version — more sensitive to small-cluster agreement)

  Δ_H(k)    = H(P^rpost|P^rpre)     − H(P^rpost|P^tt)
              (directed conditional-entropy; tt reduces rpost's
              residual uncertainty more than rpre does)

  Δ_NMI(k)  = NMI(P^tt, P^rpost)    − NMI(P^rpre, P^rpost)
              (NMI version — normalizes MI by entropy)

All are POSITIVE when P^rpost is closer to P^tt than to P^rpre at
scale k. Trace hypothesis: positive across patients, at least in
some k-range, for some bands.

Per (band, k): fraction of patients with positive contrast.
Per band: longest contiguous k-run at ≥ 7/9 and 9/9 patient unanimity.
FDR-BH across 6 bands on per-band mean contrasts at each k.

Reads:  LRG imcoh_abs caches.
Writes: data/reports/imcoh_vi/h2_partition_multiscale_raw.csv
        data/reports/imcoh_vi/h2_partition_multiscale.md
        data/reports/imcoh_vi/figures/h2_partition_multiscale.{pdf,png}
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
)

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.utils.metrics import compute_vi, conditional_entropy
from lrg_eegfc.workflow.lrg import load_lrg_result

from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, bh_fdr


K_RANGE = list(range(2, 50))


def _load_Z(pat: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, "imcoh_abs", IMCOH_LRG_CACHE)
    except Exception:
        return None
    return np.asarray(r.linkage_matrix) if r is not None else None


def _partitions_per_k(Z_rpre, Z_tt, Z_rpost):
    """Yield (k, P_rpre, P_tt, P_rpost) with matching N."""
    if Z_rpre is None or Z_tt is None or Z_rpost is None:
        return
    n = Z_rpre.shape[0] + 1
    if Z_tt.shape[0] + 1 != n or Z_rpost.shape[0] + 1 != n:
        return
    for k in K_RANGE:
        yield (k,
               fcluster(Z_rpre,  k, criterion="maxclust"),
               fcluster(Z_tt,    k, criterion="maxclust"),
               fcluster(Z_rpost, k, criterion="maxclust"))


def collect_raw() -> pd.DataFrame:
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            Z_rpre  = _load_Z(pat, "rest_pre",  band)
            Z_tt    = _load_Z(pat, "task_test", band)
            Z_rpost = _load_Z(pat, "rest_post", band)
            for k, P_rpre, P_tt, P_rpost in _partitions_per_k(Z_rpre, Z_tt, Z_rpost):
                vi_rpre_rpost = compute_vi(P_rpre, P_rpost)
                vi_tt_rpost   = compute_vi(P_tt,   P_rpost)
                ari_tt_rpost   = adjusted_rand_score(P_tt,   P_rpost)
                ari_rpre_rpost = adjusted_rand_score(P_rpre, P_rpost)
                h_rpost_given_tt   = conditional_entropy(P_rpost, P_tt)
                h_rpost_given_rpre = conditional_entropy(P_rpost, P_rpre)
                nmi_tt_rpost   = normalized_mutual_info_score(P_tt,   P_rpost)
                nmi_rpre_rpost = normalized_mutual_info_score(P_rpre, P_rpost)
                rows.append({
                    "patient": pat, "band": band, "k": k,
                    "d_VI":  vi_rpre_rpost - vi_tt_rpost,
                    "d_ARI": ari_tt_rpost - ari_rpre_rpost,
                    "d_H":   h_rpost_given_rpre - h_rpost_given_tt,
                    "d_NMI": nmi_tt_rpost - nmi_rpre_rpost,
                })
    return pd.DataFrame(rows)


# ─────────────────────────── band-level stats ───────────────────────────


def _unanimity_map(df: pd.DataFrame, contrast: str) -> tuple[np.ndarray, np.ndarray]:
    """Return (fraction_positive, mean_contrast) per (band, k)."""
    B, K = len(BRAIN_BANDS_NAMES), len(K_RANGE)
    frac = np.full((B, K), np.nan)
    mean = np.full((B, K), np.nan)
    for ib, band in enumerate(BRAIN_BANDS_NAMES):
        bdf = df[df["band"] == band]
        for jk, k in enumerate(K_RANGE):
            sub = bdf[bdf["k"] == k][contrast].to_numpy()
            sub = sub[np.isfinite(sub)]
            if sub.size == 0:
                continue
            frac[ib, jk] = (sub > 0).sum() / sub.size
            mean[ib, jk] = sub.mean()
    return frac, mean


def _longest_run(mask: np.ndarray) -> tuple[int, int]:
    best_len = best_start = best_end = 0
    cur_start = None
    for i, v in enumerate(mask):
        if v and cur_start is None:
            cur_start = i
        elif not v and cur_start is not None:
            length = i - cur_start
            if length > best_len:
                best_len, best_start, best_end = length, cur_start, i - 1
            cur_start = None
    if cur_start is not None:
        length = len(mask) - cur_start
        if length > best_len:
            best_len, best_start, best_end = length, cur_start, len(mask) - 1
    return (best_start, best_end) if best_len > 0 else (-1, -1)


# ─────────────────────────── figure ───────────────────────────


def _cmap():
    return LinearSegmentedColormap.from_list(
        "trace",
        [(1.0, 1.0, 1.0), (1.0, 0.95, 0.7),
         (1.0, 0.6, 0.2), (0.7, 0.05, 0.05)],
        N=256,
    )


def draw_landscape(ax, frac: np.ndarray, title: str, *, show_yticks=True):
    cmap = _cmap()
    im = ax.imshow(frac, aspect="auto", cmap=cmap, vmin=0.0, vmax=1.0,
                    extent=[K_RANGE[0] - 0.5, K_RANGE[-1] + 0.5,
                            len(BRAIN_BANDS_NAMES) - 0.5, -0.5],
                    interpolation="nearest")
    # Dots at ≥ 7/9 and ≥ 9/9
    ys, xs = np.where(frac >= 7/9)
    if ys.size:
        ax.scatter(np.array(K_RANGE)[xs], ys, s=1.2, c="black", alpha=0.4)
    ys, xs = np.where(frac >= 0.999)
    if ys.size:
        ax.scatter(np.array(K_RANGE)[xs], ys, s=4, c="black")
    if show_yticks:
        ax.set_yticks(range(len(BRAIN_BANDS_NAMES)))
        ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                            fontsize=11)
    else:
        ax.set_yticks([])
    ax.set_xlabel("k", fontsize=10)
    ax.set_title(title, fontsize=11, loc="left")
    return im


def _band_summary(frac: np.ndarray) -> list[dict]:
    """Per-band longest runs at 7/9 and 9/9 unanimity."""
    out = []
    for ib, band in enumerate(BRAIN_BANDS_NAMES):
        mask_7 = frac[ib] >= 7/9
        mask_9 = frac[ib] >= 0.999
        s7, e7 = _longest_run(mask_7)
        s9, e9 = _longest_run(mask_9)
        out.append({
            "band": band,
            "mean_frac": float(np.nanmean(frac[ib])),
            "len7": int(mask_7.sum()) if mask_7.any() else 0,
            "run7_k": (K_RANGE[s7], K_RANGE[e7]) if s7 >= 0 else None,
            "run7_len": (e7 - s7 + 1) if s7 >= 0 else 0,
            "len9": int(mask_9.sum()) if mask_9.any() else 0,
            "run9_k": (K_RANGE[s9], K_RANGE[e9]) if s9 >= 0 else None,
            "run9_len": (e9 - s9 + 1) if s9 >= 0 else 0,
        })
    return out


# ─────────────────────────── main ───────────────────────────


def main() -> None:
    out_dir = REPORTS_ROOT / "imcoh_vi"
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    df = collect_raw()
    df.to_csv(out_dir / "h2_partition_multiscale_raw.csv", index=False)
    print(f"{len(df)} raw rows")

    contrasts = [
        ("d_VI",  r"$\Delta_{VI}(k)$ — VI(rpre,rpost) − VI(tt,rpost)"),
        ("d_ARI", r"$\Delta_{ARI}(k)$ — ARI(tt,rpost) − ARI(rpre,rpost)"),
        ("d_H",   r"$\Delta_H(k)$ — H(rpost|rpre) − H(rpost|tt)"),
        ("d_NMI", r"$\Delta_{NMI}(k)$ — NMI(tt,rpost) − NMI(rpre,rpost)"),
    ]

    # ── Build figure: 4 unanimity maps + band-summary strip  ──
    fig, axes = plt.subplots(
        2, 2, figsize=(14.0, 7.0), dpi=160, sharex=True,
        gridspec_kw={"hspace": 0.48, "wspace": 0.15},
    )
    summaries: dict[str, list[dict]] = {}
    for ax, (col, title) in zip(axes.ravel(), contrasts):
        frac, _ = _unanimity_map(df, col)
        im = draw_landscape(ax, frac, title)
        summaries[col] = _band_summary(frac)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(),
                         shrink=0.85, pad=0.02, fraction=0.04)
    cbar.set_label("fraction of patients with Δ > 0", fontsize=10)
    fig.suptitle(
        "Multiscale partition-level directed contrasts per (band, k).\n"
        "Dots: ≥ 7/9 (faint) and 9/9 (black) patient unanimity. "
        "Positive contrast = P^rpost closer to P^ttest than to P^rpre at scale k.",
        fontsize=11, y=1.00,
    )
    for ax in axes[-1]:
        ax.set_xlabel("k (dendrogram cut scale)", fontsize=10)
    for ext in ("pdf", "png"):
        out = fig_dir / f"h2_partition_multiscale.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"saved {out}")
    plt.close(fig)

    # ── Markdown report ─────────────────────────────────────
    lines: list[str] = []
    ap = lines.append
    ap("# H2-PART — multiscale partition-level directed contrasts")
    ap("")
    ap(f"For every (patient, band, k ∈ [{K_RANGE[0]}, {K_RANGE[-1]}]) and "
       "four partition-distance contrasts (VI, ARI, conditional H, NMI), "
       "we count per-k unanimity (fraction of 9 patients with positive "
       "contrast) and report the longest contiguous k-run at ≥ 7/9 and "
       "9/9 levels per band.")
    ap("")
    for col, title in contrasts:
        rows = summaries[col]
        ap(f"## {title}")
        ap("")
        ap("| band | mean fraction | longest ≥7/9 run | longest 9/9 run |")
        ap("|------|-------------:|:-----------------|:----------------|")
        for r in rows:
            run7 = (f"k={r['run7_k'][0]}–{r['run7_k'][1]} (L={r['run7_len']})"
                    if r["run7_k"] else "—")
            run9 = (f"k={r['run9_k'][0]}–{r['run9_k'][1]} (L={r['run9_len']})"
                    if r["run9_k"] else "—")
            ap(f"| {BRAIN_BAND_TEX_DICT[r['band']]} | "
               f"{r['mean_frac']:.3f} | {run7} | {run9} |")
        ap("")

    md = out_dir / "h2_partition_multiscale.md"
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {md}")


if __name__ == "__main__":
    main()
