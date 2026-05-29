#!/usr/bin/env python3
"""Definitive synthesis figure: best-balanced reorganization metric.

Compares mean_VI (best composite) and mean_ARI (best H3) side by side,
showing ALL hypotheses including band-dependent task trace.

Output: data/figures/metric_exploration/final_synthesis/
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.config.paths import FIGURES_ROOT

PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
PATIENTS_4PH = PATIENTS_4PHASE
BANDS = BRAIN_BANDS_NAMES
OUT_DIR = FIGURES_ROOT / "metric_exploration" / "final_synthesis"

ALL_PAIRS = [
    ("rest_pre", "rest_post"), ("task_learn", "task_test"),
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"),
    ("task_learn", "rest_post"), ("task_test", "rest_post"),
]
PS = {  # short labels
    ("rest_pre", "rest_post"): "Pre↔Post", ("task_learn", "task_test"): "TL↔TT",
    ("rest_pre", "task_learn"): "Pre↔TL", ("rest_pre", "task_test"): "Pre↔TT",
    ("task_learn", "rest_post"): "TL↔Post", ("task_test", "rest_post"): "TT↔Post",
}
PCAT = {
    ("rest_pre", "rest_post"): "W", ("task_learn", "task_test"): "W",
    ("rest_pre", "task_learn"): "C", ("rest_pre", "task_test"): "C",
    ("task_learn", "rest_post"): "C", ("task_test", "rest_post"): "C",
}

def bl(b):
    return BRAIN_BAND_TEX_DICT.get(b, b)

def gv(df, pat, band, p1, p2, col):
    r = df[(df["patient"] == pat) & (df["band"] == band) &
           (df["phase1"] == p1) & (df["phase2"] == p2)]
    return r[col].values[0] if len(r) > 0 else np.nan


def make_similarity(df, col, invert=False):
    """Convert to similarity: if invert=True, negate (for VI-type metrics)."""
    vals = {}
    for _, r in df.iterrows():
        v = r[col]
        if invert:
            v = -v
        vals[(r["patient"], r["band"], r["phase1"], r["phase2"])] = v
    return vals


def compute_contrasts(sim_vals):
    """Compute all contrasts per (patient, band)."""
    results = {}
    for pat in PATIENTS_4PH:
        for band in BANDS:
            pv = {}
            for pair in ALL_PAIRS:
                k = (pat, band, pair[0], pair[1])
                if k in sim_vals:
                    pv[pair] = sim_vals[k]
            if len(pv) != 6:
                continue

            within = [pv[p] for p in pv if PCAT[p] == "W"]
            cross = [pv[p] for p in pv if PCAT[p] == "C"]
            tt = pv[("task_learn", "task_test")]
            rr = pv[("rest_pre", "rest_post")]
            tt_post = pv[("task_test", "rest_post")]
            tl_post = pv[("task_learn", "rest_post")]
            pre_tl = pv[("rest_pre", "task_learn")]
            pre_tt = pv[("rest_pre", "task_test")]

            results[(pat, band)] = {
                "within_cross": np.mean(within) - np.mean(cross),
                "task_stability": tt - np.mean(cross),
                "rest_stability": rr - np.mean(cross),
                "task_minus_rest": tt - rr,
                "tt_post_vs_pre_post": tt_post - rr,
                "tl_post_vs_pre_post": tl_post - rr,
                "tt_post_vs_pre_tt": tt_post - pre_tt,
                "tt_post_vs_pre_tl": tt_post - pre_tl,
                "pair_values": pv,
            }
    return results


def unanimity(vals):
    """Check if all values have same sign. Return (n_agree, N, direction, unanimous, relaxed)."""
    n_pos = sum(v > 0 for v in vals)
    N = len(vals)
    n_agree = max(n_pos, N - n_pos)
    direction = "+" if n_pos > N - n_pos else "−"
    return n_agree, N, direction, n_agree == N, n_agree >= N - 1


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(FIGURES_ROOT / "metric_exploration" / "partition_multiscale" / "results.csv")
    print(f"Loaded {len(df)} rows")

    # Build similarity dicts
    ari_sim = make_similarity(df, "mean_ARI", invert=False)
    vi_sim = make_similarity(df, "mean_VI", invert=True)  # lower VI = more similar

    # Compute contrasts
    ari_c = compute_contrasts(ari_sim)
    vi_c = compute_contrasts(vi_sim)

    metrics = [
        ("mean ARI", ari_sim, ari_c),
        ("mean VI", vi_sim, vi_c),
    ]

    # ================================================================
    # FIGURE 1: Main comparison — 2 rows (ARI, VI) × 4 panels
    # ================================================================
    fig = plt.figure(figsize=(24, 16))
    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.35, wspace=0.3)

    for m_idx, (mname, sim_vals, contrasts) in enumerate(metrics):

        # --- Panel A: Raw similarity profiles ---
        ax_a = fig.add_subplot(gs[m_idx, 0])
        for pair in ALL_PAIRS:
            p1, p2 = pair
            cat = PCAT[pair]
            ls = '-' if cat == "W" else '--'
            lw = 2.5 if cat == "W" else 1.0
            ms = 7 if cat == "W" else 4
            mk = 's' if cat == "W" else 'o'
            means = []
            for band in BANDS:
                vals = [sim_vals.get((pat, band, p1, p2), np.nan) for pat in PATIENTS_4PH]
                means.append(np.nanmean(vals))
            ax_a.plot(range(len(BANDS)), means, f'{mk}{ls}', linewidth=lw,
                      markersize=ms, label=f"{PS[pair]}[{cat}]")

        ax_a.set_xticks(range(len(BANDS)))
        ax_a.set_xticklabels([bl(b) for b in BANDS], fontsize=9, rotation=20)
        ax_a.set_ylabel("Similarity")
        ax_a.set_title(f"{mname}: pair profiles", fontsize=10, fontweight='bold')
        ax_a.legend(fontsize=6, ncol=2, loc='lower right')

        # --- Panel B: Within-Cross gap heatmap ---
        ax_b = fig.add_subplot(gs[m_idx, 1])
        gap_mat = np.full((len(BANDS), len(PATIENTS_4PH)), np.nan)
        for b_idx, band in enumerate(BANDS):
            for p_idx, pat in enumerate(PATIENTS_4PH):
                c = contrasts.get((pat, band))
                if c:
                    gap_mat[b_idx, p_idx] = c["within_cross"]

        vmax = max(0.01, np.nanmax(np.abs(gap_mat)) * 1.1)
        im = ax_b.imshow(gap_mat, cmap='Greens', vmin=0, vmax=vmax, aspect='auto')
        for i in range(len(BANDS)):
            for j in range(len(PATIENTS_4PH)):
                v = gap_mat[i, j]
                if not np.isnan(v):
                    color = 'white' if v > vmax * 0.6 else 'black'
                    ax_b.text(j, i, f"{v:.3f}", ha='center', va='center',
                              fontsize=8, fontweight='bold', color=color)
        ax_b.set_xticks(range(len(PATIENTS_4PH)))
        ax_b.set_xticklabels([p[-2:] for p in PATIENTS_4PH], fontsize=9)
        ax_b.set_yticks(range(len(BANDS)))
        ax_b.set_yticklabels([bl(b) for b in BANDS], fontsize=10)
        all_pos = np.all(gap_mat > 0)
        n_pos_bands = sum(np.all(gap_mat[i, :] > 0) for i in range(len(BANDS)))
        ax_b.set_title(f"H3: Within−Cross gap\n({n_pos_bands}/6 bands ALL+{'★' if all_pos else ''})",
                        fontsize=10, fontweight='bold')

        # --- Panel C: Contrasts dot plot ---
        ax_c = fig.add_subplot(gs[m_idx, 2])
        contrast_keys = [
            ("within_cross", "Within−Cross"),
            ("task_stability", "TL↔TT−Cross"),
            ("task_minus_rest", "TL↔TT−Pre↔Post"),
        ]
        y_pos = np.arange(len(BANDS))
        colors_c = ['#2196F3', '#4CAF50', '#9C27B0']

        for c_idx, (ckey, clabel) in enumerate(contrast_keys):
            x_offset = c_idx * 0.12
            for b_idx, band in enumerate(BANDS):
                vals = [contrasts[(pat, band)][ckey] for pat in PATIENTS_4PH
                        if (pat, band) in contrasts]
                if not vals:
                    continue
                m = np.mean(vals)
                na, N, d, unan, relax = unanimity(vals)
                ax_c.scatter(m, b_idx + x_offset - 0.12, marker='D', s=50,
                             color=colors_c[c_idx], edgecolors='black', linewidths=0.3)
                if unan:
                    ax_c.annotate("★", (m, b_idx + x_offset - 0.12),
                                  textcoords="offset points", xytext=(6, 0),
                                  fontsize=9, color='gold', fontweight='bold', va='center')

        ax_c.axvline(0, color='gray', linestyle='--', linewidth=0.8)
        ax_c.set_yticks(y_pos)
        ax_c.set_yticklabels([bl(b) for b in BANDS], fontsize=10)
        ax_c.invert_yaxis()
        ax_c.set_xlabel("Contrast value")

        from matplotlib.lines import Line2D
        legend_els = [Line2D([0], [0], marker='D', color=colors_c[i], linestyle='None',
                              markersize=6, label=clabel)
                       for i, (_, clabel) in enumerate(contrast_keys)]
        ax_c.legend(handles=legend_els, fontsize=7, loc='lower left')
        ax_c.set_title(f"H1+H3: contrasts\n(★ = unanimous)", fontsize=10, fontweight='bold')

        # --- Panel D: Band-dependent task trace (H2) ---
        ax_d = fig.add_subplot(gs[m_idx, 3])
        trace_keys = [
            ("tt_post_vs_pre_post", "TT↔Post − Pre↔Post"),
            ("tt_post_vs_pre_tt", "TT↔Post − Pre↔TT"),
        ]
        colors_t = ['#E91E63', '#FF5722']

        for t_idx, (tkey, tlabel) in enumerate(trace_keys):
            x_offset = t_idx * 0.15
            for b_idx, band in enumerate(BANDS):
                vals = [contrasts[(pat, band)][tkey] for pat in PATIENTS_4PH
                        if (pat, band) in contrasts]
                if not vals:
                    continue

                # Plot individual patients
                for v in vals:
                    clr = '#1565C0' if v > 0 else '#C62828'
                    ax_d.scatter(v, b_idx + x_offset - 0.075, s=15, alpha=0.5,
                                 color=clr, edgecolors='none')

                m = np.mean(vals)
                na, N, d, unan, relax = unanimity(vals)
                ax_d.scatter(m, b_idx + x_offset - 0.075, marker='D', s=50,
                             color=colors_t[t_idx], edgecolors='black', linewidths=0.4)
                if unan:
                    symbol = "★+" if d == "+" else "★−"
                    color = '#1565C0' if d == "+" else '#C62828'
                    ax_d.annotate(symbol, (m, b_idx + x_offset - 0.075),
                                  textcoords="offset points", xytext=(8, 0),
                                  fontsize=9, color=color, fontweight='bold', va='center')
                elif relax:
                    ax_d.annotate("◆", (m, b_idx + x_offset - 0.075),
                                  textcoords="offset points", xytext=(8, 0),
                                  fontsize=8, color='orange', va='center')

        ax_d.axvline(0, color='gray', linestyle='--', linewidth=0.8)
        ax_d.set_yticks(y_pos)
        ax_d.set_yticklabels([bl(b) for b in BANDS], fontsize=10)
        ax_d.invert_yaxis()
        ax_d.set_xlabel("Contrast value")

        legend_els2 = [Line2D([0], [0], marker='D', color=colors_t[i], linestyle='None',
                               markersize=6, label=tlabel)
                        for i, (_, tlabel) in enumerate(trace_keys)]
        ax_d.legend(handles=legend_els2, fontsize=7, loc='lower left')
        ax_d.set_title(f"H2: task trace in rest_post\n(★+=persists, ★−=recovers)", fontsize=10,
                        fontweight='bold')

    fig.suptitle("Hierarchical reorganization: mean ARI (top) vs mean VI (bottom)\n"
                 "LRG dendrograms on MSC — flat clustering k=2..20 — 4-phase patients",
                 fontsize=14, fontweight='bold', y=0.99)
    fig.savefig(OUT_DIR / "final_comparison.pdf", bbox_inches='tight', dpi=150)
    plt.close(fig)
    print("Saved final_comparison.pdf")

    # ================================================================
    # FIGURE 2: Per-band phase×phase heatmaps for BOTH metrics
    # ================================================================
    phases_4 = ["rest_pre", "task_learn", "task_test", "rest_post"]
    phase_labels = ["rest_pre", "tLearn", "tTest", "rest_post"]

    for mname, sim_vals, _ in metrics:
        fig2, axes2 = plt.subplots(2, 3, figsize=(16, 10))
        fig2.suptitle(f"{mname}: phase×phase similarity per band\n"
                      "(averaged over 4 patients)", fontsize=13, fontweight='bold')
        for b_idx, band in enumerate(BANDS):
            ax = axes2[b_idx // 3, b_idx % 3]
            mat = np.eye(4)
            for i, pi in enumerate(phases_4):
                for j, pj in enumerate(phases_4):
                    if i >= j:
                        continue
                    vals = [sim_vals.get((pat, band, pi, pj), np.nan) for pat in PATIENTS_4PH]
                    vals = [v for v in vals if not np.isnan(v)]
                    if vals:
                        m = np.mean(vals)
                        mat[i, j] = m
                        mat[j, i] = m

            vmin_hm = np.min(mat[mat < 1]) - 0.02
            im = ax.imshow(mat, cmap='YlOrRd', vmin=vmin_hm, vmax=1.0, aspect='equal')
            for i in range(4):
                for j in range(4):
                    color = 'white' if mat[i, j] > 0.75 else 'black'
                    ax.text(j, i, f"{mat[i,j]:.3f}", ha='center', va='center',
                            fontsize=10, fontweight='bold', color=color)
            ax.set_xticks(range(4))
            ax.set_xticklabels(phase_labels, fontsize=9, rotation=30)
            ax.set_yticks(range(4))
            ax.set_yticklabels(phase_labels, fontsize=9)
            ax.set_title(bl(band), fontsize=12, fontweight='bold')

        fig2.tight_layout(rect=[0, 0, 0.92, 0.93])
        cbar_ax = fig2.add_axes([0.93, 0.15, 0.02, 0.7])
        fig2.colorbar(im, cax=cbar_ax, label="Similarity")
        safe_name = mname.replace(" ", "_")
        fig2.savefig(OUT_DIR / f"phase_heatmaps_{safe_name}.pdf", bbox_inches='tight', dpi=150)
        plt.close(fig2)
        print(f"Saved phase_heatmaps_{safe_name}.pdf")

    # ================================================================
    # FIGURE 3: All 6 patients for the best metric (mean_VI)
    # ================================================================
    fig3, axes3 = plt.subplots(2, 3, figsize=(18, 10))
    fig3.suptitle("mean VI (inverted): per-patient similarity profiles\n"
                  "All 6 patients — solid=within, dashed=cross", fontsize=13, fontweight='bold')

    pair_colors = {
        ("rest_pre", "rest_post"): '#6A1B9A',
        ("task_learn", "task_test"): '#2E7D32',
        ("rest_pre", "task_learn"): '#E65100',
        ("rest_pre", "task_test"): '#BF360C',
        ("task_learn", "rest_post"): '#AD1457',
        ("task_test", "rest_post"): '#880E4F',
    }

    for p_idx, (pat, phases) in enumerate(PATIENT_PHASES.items()):
        ax = axes3[p_idx // 3, p_idx % 3]
        pairs = list(combinations(phases, 2))
        for pair in pairs:
            p1, p2 = pair
            cat = PCAT.get(pair, "?")
            color = pair_colors.get(pair, 'gray')
            ls = '-' if cat == "W" else '--'
            lw = 2.5 if cat == "W" else 1.0
            ms = 7 if cat == "W" else 4
            mk = 's' if cat == "W" else 'o'
            vals = [vi_sim.get((pat, band, p1, p2), np.nan) for band in BANDS]
            ax.plot(range(len(BANDS)), vals, f'{mk}{ls}', color=color,
                    linewidth=lw, markersize=ms,
                    label=PS.get(pair, f"{p1}↔{p2}"), alpha=0.8)

        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([bl(b) for b in BANDS], fontsize=8, rotation=20)
        ax.set_ylabel("Similarity (−VI)")
        ax.set_title(f"{pat} ({len(phases)} phases)", fontsize=11)
        ax.legend(fontsize=6, ncol=2, loc='lower right')

    fig3.tight_layout(rect=[0, 0, 1, 0.93])
    fig3.savefig(OUT_DIR / "all_patients_mean_VI.pdf", bbox_inches='tight', dpi=150)
    plt.close(fig3)
    print("Saved all_patients_mean_VI.pdf")

    # ================================================================
    # Print final summary
    # ================================================================
    print("\n" + "=" * 70)
    print("FINAL RESULTS")
    print("=" * 70)

    for mname, sim_vals, contrasts in metrics:
        print(f"\n--- {mname} ---")

        # H3
        print("  H3 (Within > Cross):")
        for band in BANDS:
            vals = [contrasts[(pat, band)]["within_cross"] for pat in PATIENTS_4PH
                    if (pat, band) in contrasts]
            na, N, d, u, r = unanimity(vals)
            flag = "★" if u else ("◆" if r else " ")
            print(f"    {flag} {bl(band):<14} {d} ({na}/{N}) "
                  f"  {' '.join(f'{v:+.4f}' for v in vals)}")

        # H1
        print("  H1 (Task stability TL↔TT − Cross):")
        for band in BANDS:
            vals = [contrasts[(pat, band)]["task_stability"] for pat in PATIENTS_4PH
                    if (pat, band) in contrasts]
            na, N, d, u, r = unanimity(vals)
            flag = "★" if u else ("◆" if r else " ")
            print(f"    {flag} {bl(band):<14} {d} ({na}/{N}) "
                  f"  {' '.join(f'{v:+.4f}' for v in vals)}")

        # H2 variants
        for h2_key, h2_label in [
            ("tt_post_vs_pre_post", "H2a (TT↔Post − Pre↔Post)"),
            ("tt_post_vs_pre_tt", "H2b (TT↔Post − Pre↔TT)"),
        ]:
            print(f"  {h2_label}:")
            for band in BANDS:
                vals = [contrasts[(pat, band)][h2_key] for pat in PATIENTS_4PH
                        if (pat, band) in contrasts]
                na, N, d, u, r = unanimity(vals)
                flag = "★" if u else ("◆" if r else " ")
                meaning = ""
                if u:
                    meaning = "→ PERSISTS" if d == "+" else "→ RECOVERS"
                print(f"    {flag} {bl(band):<14} {d} ({na}/{N}) "
                      f"  {' '.join(f'{v:+.4f}' for v in vals)}  {meaning}")

    print(f"\nAll figures saved to: {OUT_DIR}")
    print("Done!")


if __name__ == "__main__":
    main()
