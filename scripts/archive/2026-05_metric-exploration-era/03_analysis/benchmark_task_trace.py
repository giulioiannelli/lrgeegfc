#!/usr/bin/env python3
"""Benchmark ALL methods for the band-dependent task trace hypothesis.

For each method and each band, compute the contrast:
  C = sim(task_test, rest_post) - sim(rest_pre, rest_post)

The hypothesis is NOT that C is always positive. Instead:
  - For SOME bands, C > 0 unanimously (task persists in rest_post)
  - For OTHER bands, C < 0 unanimously (rest recovers to rest_pre)
  - The PATTERN (which bands persist, which recover) is consistent across patients

A good method shows many bands with unanimous sign (either direction),
ideally with BOTH directions present.

Also tests the "directional task reorganization" contrast:
  C_dir = sim(task_test, rest_post) - sim(rest_pre, task_learn)
  → Is the post-task rest closer to the task than the pre-task rest was?
  If positive: task pulled rest_post toward task structure.

Output: data/figures/metric_exploration/synthesis/
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.config.paths import FIGURES_ROOT

BANDS = BRAIN_BANDS_NAMES
PATIENTS_4PH = PATIENTS_4PHASE
BASE = FIGURES_ROOT / "metric_exploration"
OUT_DIR = BASE / "synthesis"

ALL_PAIRS = [
    ("rest_pre", "rest_post"), ("task_learn", "task_test"),
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"),
    ("task_learn", "rest_post"), ("task_test", "rest_post"),
]
PAIR_CAT = {
    ("rest_pre", "rest_post"): "within", ("task_learn", "task_test"): "within",
    ("rest_pre", "task_learn"): "cross", ("rest_pre", "task_test"): "cross",
    ("task_learn", "rest_post"): "cross", ("task_test", "rest_post"): "cross",
}

def bl(b):
    return BRAIN_BAND_TEX_DICT.get(b, b)


def load_wide(csv_path, col):
    """Load wide-format CSV, return (pat, band, p1, p2) → value."""
    df = pd.read_csv(csv_path)
    if col not in df.columns:
        return None
    vals = {}
    for _, r in df.iterrows():
        vals[(r["patient"], r["band"], r["phase1"], r["phase2"])] = r[col]
    return vals


def load_long(csv_path, metric_name):
    """Load long-format CSV (D(τ) distances), return (pat, band, p1, p2) → value."""
    df = pd.read_csv(csv_path)
    if "metric" not in df.columns:
        return None
    sub = df[df["metric"] == metric_name]
    vals = {}
    for _, r in sub.iterrows():
        vals[(r["patient"], r["band"], r["phase1"], r["phase2"])] = r["sim_value"]
    return vals


def test_band_trace(vals, is_similarity=True):
    """For each band, compute multiple contrasts and check sign consistency.

    Returns dict: band → {contrast_name → {vals, sign, n_agree, N, unanimous}}
    """
    contrasts = {
        "TT↔Post − Pre↔Post": lambda pv: pv.get(("task_test", "rest_post"), np.nan) -
                                           pv.get(("rest_pre", "rest_post"), np.nan),
        "TL↔Post − Pre↔Post": lambda pv: pv.get(("task_learn", "rest_post"), np.nan) -
                                           pv.get(("rest_pre", "rest_post"), np.nan),
        "TT↔Post − Pre↔TL": lambda pv: pv.get(("task_test", "rest_post"), np.nan) -
                                         pv.get(("rest_pre", "task_learn"), np.nan),
        "TT↔Post − Pre↔TT": lambda pv: pv.get(("task_test", "rest_post"), np.nan) -
                                         pv.get(("rest_pre", "task_test"), np.nan),
        "Within − Cross": lambda pv: np.mean([pv[p] for p in pv if PAIR_CAT[p] == "within"]) -
                                     np.mean([pv[p] for p in pv if PAIR_CAT[p] == "cross"]),
        "TL↔TT − Cross": lambda pv: pv.get(("task_learn", "task_test"), np.nan) -
                                     np.mean([pv[p] for p in pv if PAIR_CAT[p] == "cross"]),
    }

    results = {}
    for band in BANDS:
        band_res = {}
        for cname, cfunc in contrasts.items():
            c_vals = []
            for pat in PATIENTS_4PH:
                pv = {}
                for pair in ALL_PAIRS:
                    p1, p2 = pair
                    k = (pat, band, p1, p2)
                    if k in vals:
                        v = vals[k] if is_similarity else -vals[k]
                        pv[pair] = v
                if len(pv) == 6:
                    try:
                        c = cfunc(pv)
                        if not np.isnan(c):
                            c_vals.append(c)
                    except Exception:
                        pass

            if c_vals:
                n_pos = sum(v > 0 for v in c_vals)
                N = len(c_vals)
                n_agree = max(n_pos, N - n_pos)
                band_res[cname] = {
                    "vals": c_vals,
                    "mean": np.mean(c_vals),
                    "n_pos": n_pos,
                    "N": N,
                    "n_agree": n_agree,
                    "unanimous": n_agree == N and N >= 3,
                    "relaxed": n_agree >= N - 1 and N >= 3,
                    "direction": "+" if n_pos > N - n_pos else "−",
                }

        results[band] = band_res
    return results


def score_method(results, contrast_name):
    """Score a method for a given contrast:
    - n_unanimous: bands with unanimous sign (either direction)
    - n_relaxed: bands with N-1 agreement
    - has_both: does it show BOTH positive and negative unanimous bands?
    """
    n_unan = 0
    n_relax = 0
    dirs = []
    for band in BANDS:
        info = results[band].get(contrast_name)
        if info:
            if info["unanimous"]:
                n_unan += 1
                dirs.append(info["direction"])
            if info["relaxed"]:
                n_relax += 1

    has_both = "+" in dirs and "−" in dirs
    return n_unan, n_relax, has_both, dirs


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Build list of all methods
    methods = []  # (name, vals_dict, is_similarity)

    # logCosine D(τ)
    lc = BASE / "logcosine_dtau" / "logcosine_dtau_results.csv"
    if lc.exists():
        v = load_wide(lc, "logCosine")
        if v: methods.append(("logCosine D(τ)", v, True))

    # D(τ) distance variants (long format)
    dtau = BASE / "dtau_distances" / "results.csv"
    if dtau.exists():
        df_d = pd.read_csv(dtau)
        for m in df_d["metric"].unique():
            v = load_long(dtau, m)
            if v: methods.append((f"D(τ) {m}", v, True))

    # Partition metrics
    part = BASE / "partition_multiscale" / "results.csv"
    if part.exists():
        for col in ["mean_ARI", "max_ARI", "weighted_ARI", "ARI_k3", "ARI_k5",
                     "ARI_k10", "ARI_k15", "mean_NMI", "mean_VI",
                     "ARI_natural_k", "NMI_natural_k"]:
            v = load_wide(part, col)
            if v:
                is_sim = "VI" not in col
                methods.append((f"Part. {col}", v, is_sim))

    # Community tracking
    comm = BASE / "community_tracking" / "results.csv"
    if comm.exists():
        for col in ["mean_fraction_stable", "min_fraction_stable",
                     "fraction_stable_natural", "fraction_stable_best_k",
                     "jaccard_mean_natural", "jaccard_mean_best_k",
                     "fraction_stable_k3", "fraction_stable_k5"]:
            v = load_wide(comm, col)
            if v: methods.append((f"Comm. {col}", v, True))

    # Tree structure
    tree = BASE / "tree_structure" / "results.csv"
    if tree.exists():
        for col in ["RF_sim", "weighted_RF", "branch_score", "merge_order", "height_dist"]:
            v = load_wide(tree, col)
            if v: methods.append((f"Tree {col}", v, True))

    print(f"Loaded {len(methods)} methods\n")

    # Test all methods
    key_contrasts = [
        "TT↔Post − Pre↔Post",   # H2a: task trace
        "TL↔Post − Pre↔Post",   # H2b: task influence
        "TT↔Post − Pre↔TL",     # H2c: directional reorganization
        "TT↔Post − Pre↔TT",     # H2d: alternative directional
        "Within − Cross",        # H3: category effect
        "TL↔TT − Cross",        # H1-like: task stability
    ]

    all_scores = {}
    all_results = {}
    for name, vals, is_sim in methods:
        results = test_band_trace(vals, is_sim)
        all_results[name] = results
        scores = {}
        for cname in key_contrasts:
            n_u, n_r, has_both, dirs = score_method(results, cname)
            scores[cname] = {"unan": n_u, "relax": n_r, "both": has_both, "dirs": dirs}
        all_scores[name] = scores

    # ================================================================
    # Print comprehensive ranking for each contrast
    # ================================================================
    for cname in key_contrasts:
        print(f"\n{'='*80}")
        print(f"CONTRAST: {cname}")
        print(f"{'='*80}")

        ranked = sorted(all_scores.items(),
                         key=lambda x: (-x[1][cname]["unan"], -x[1][cname]["relax"]))
        print(f"{'Rank':<5} {'Method':<32} {'Unan':>5} {'Relax':>6} {'Both?':>6} {'Per-band dirs':>20}")
        for rank, (name, scores) in enumerate(ranked[:15], 1):
            s = scores[cname]
            dirs_str = " ".join(s["dirs"]) if s["dirs"] else "—"
            both = "YES" if s["both"] else ("no" if s["dirs"] else "—")
            print(f"{rank:<5} {name:<32} {s['unan']:>5} {s['relax']:>6} {both:>6} {dirs_str:>20}")

    # ================================================================
    # Composite score: balance all hypotheses
    # ================================================================
    print(f"\n{'='*80}")
    print("COMPOSITE SCORE: Balancing all hypotheses")
    print(f"{'='*80}")
    print("Score = H3_unan + H1_unan + H2a_unan + H2a_both_bonus")
    print("        (within>cross, task persist, task trace, band-dependent trace)\n")

    composite = {}
    for name, scores in all_scores.items():
        h3 = scores["Within − Cross"]["unan"]
        h1 = scores["TL↔TT − Cross"]["unan"]
        # For H2, pick the best task-trace contrast
        h2_variants = [scores[c]["unan"] for c in
                       ["TT↔Post − Pre↔Post", "TL↔Post − Pre↔Post",
                        "TT↔Post − Pre↔TL", "TT↔Post − Pre↔TT"]]
        h2_best = max(h2_variants)
        h2_both = any(scores[c]["both"] for c in
                      ["TT↔Post − Pre↔Post", "TL↔Post − Pre↔Post",
                       "TT↔Post − Pre↔TL", "TT↔Post − Pre↔TT"])

        # Relaxed counts
        h3r = scores["Within − Cross"]["relax"]
        h1r = scores["TL↔TT − Cross"]["relax"]
        h2r_variants = [scores[c]["relax"] for c in
                        ["TT↔Post − Pre↔Post", "TL↔Post − Pre↔Post",
                         "TT↔Post − Pre↔TL", "TT↔Post − Pre↔TT"]]
        h2r_best = max(h2r_variants)

        score = h3 * 3 + h1 * 2 + h2_best * 2 + (2 if h2_both else 0)
        score_r = h3r * 3 + h1r * 2 + h2r_best * 2 + (2 if h2_both else 0)

        composite[name] = {
            "score": score, "score_r": score_r,
            "H3_u": h3, "H3_r": h3r,
            "H1_u": h1, "H1_r": h1r,
            "H2_u": h2_best, "H2_r": h2r_best,
            "H2_both": h2_both,
        }

    ranked_composite = sorted(composite.items(), key=lambda x: (-x[1]["score"], -x[1]["score_r"]))

    print(f"{'Rank':<5} {'Method':<32} {'Score':>6} {'Sc_r':>5} "
          f"{'H3u':>4} {'H3r':>4} {'H1u':>4} {'H1r':>4} "
          f"{'H2u':>4} {'H2r':>4} {'Both':>5}")
    print("-" * 95)
    for rank, (name, c) in enumerate(ranked_composite[:20], 1):
        print(f"{rank:<5} {name:<32} {c['score']:>6} {c['score_r']:>5} "
              f"{c['H3_u']:>4} {c['H3_r']:>4} {c['H1_u']:>4} {c['H1_r']:>4} "
              f"{c['H2_u']:>4} {c['H2_r']:>4} {'YES' if c['H2_both'] else 'no':>5}")

    # ================================================================
    # Detail for top 3 composite methods
    # ================================================================
    top3 = ranked_composite[:3]
    for name, c in top3:
        print(f"\n--- {name} ---")
        results = all_results[name]
        for cname in key_contrasts:
            print(f"  {cname}:")
            for band in BANDS:
                info = results[band].get(cname)
                if info:
                    flag = "★" if info["unanimous"] else ("◆" if info["relaxed"] else " ")
                    vals_str = " ".join(f"{v:+.4f}" for v in info["vals"])
                    print(f"    {flag} {bl(band):<14} {info['direction']} "
                          f"({info['n_pos']}/{info['N']} pos)  [{vals_str}]")
                else:
                    print(f"      {bl(band):<14} N/A")

    # ================================================================
    # Figure: Composite ranking
    # ================================================================
    fig, ax = plt.subplots(figsize=(14, 10))
    top_n = min(20, len(ranked_composite))
    methods_plot = ranked_composite[:top_n]

    cols = ["H3_u", "H1_u", "H2_u"]
    col_labels = ["H3\nWithin>Cross", "H1\nTask persist", "H2\nTask trace"]
    mat = np.zeros((top_n, len(cols)))
    for i, (name, c) in enumerate(methods_plot):
        for j, col in enumerate(cols):
            mat[i, j] = c[col]

    im = ax.imshow(mat, cmap='YlGn', aspect='auto', vmin=0, vmax=6)
    for i in range(top_n):
        for j in range(len(cols)):
            v = int(mat[i, j])
            both = methods_plot[i][1]["H2_both"]
            txt = f"{v}/6"
            if j == 2 and both:
                txt += "\n(±)"
            color = 'white' if v >= 4 else 'black'
            ax.text(j, i, txt, ha='center', va='center', fontsize=9, color=color)

    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(col_labels, fontsize=10)
    ax.set_yticks(range(top_n))
    labels = [f"{name} [S={c['score']}]" for name, c in methods_plot]
    ax.set_yticklabels(labels, fontsize=8)
    plt.colorbar(im, ax=ax, shrink=0.7, label="# bands unanimous (/6)")
    ax.set_title("Composite ranking: balancing all reorganization hypotheses\n"
                 "H3=within>cross, H1=task persist, H2=band-dependent task trace in rest_post\n"
                 "(±) = has BOTH persistent and recovering bands",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "composite_ranking.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved composite_ranking.pdf")

    # Save ranking CSV
    rows = []
    for rank, (name, c) in enumerate(ranked_composite, 1):
        rows.append({"rank": rank, "method": name, **c})
    pd.DataFrame(rows).to_csv(OUT_DIR / "composite_ranking.csv", index=False)
    print("Saved composite_ranking.csv")
    print("\nDone!")


if __name__ == "__main__":
    main()
