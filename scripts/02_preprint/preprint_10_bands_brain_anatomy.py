#!/usr/bin/env python3
"""Six-band cohort brain map of per-contact trace contribution.

For each band, every cohort contact is placed at its MNI coordinate
on a glass-brain projection and colored by its **per-contact trace
contribution**:

    s(c) = mean over pairs (c, j) of sign(Δ_task(c,j)) × sign(Δ_rest(c,j))

s(c) ∈ [−1, +1]; green = positive (this contact's pairs tend to be
preserved across the task→rest transition), red = negative (anti).

The cross-band visual: bands with cohort-wide trace (β, α) show many
green contacts spread across the implanted cortex; null bands (θ,
γ_h) show a mixed-color scatter with no spatial coherence.

Inputs
------
data/raw/stereoeeg_patients/<Pat>/implant_pat_<NN>.csv     (MNI coords)
data/reports/imcoh_continuous_trace/per_pair_split/<Pat>_<band>.npz

Output (PDF only)
-----------------
data/preprint/figures/all_bands/per_pair_trace/fig_bands_brain_anatomy.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from nilearn import plotting

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()


COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]

LRG_CTM_DIR = ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"
PAIR_SPLIT_DIR = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "per_pair_split"
SEEG_RAW = ROOT / "data" / "raw" / "stereoeeg_patients"

TRACE_CMAP = LinearSegmentedColormap.from_list(
    "trace_rwg",
    ["#5b1212", "#a02525", "#d8584c", "#f0a89e", "#fdf2ef",
     "#bfddc7", "#5fac7e", "#1a7c3e", "#0a3a1d"],
    N=512,
)

SHARED_VABS = 0.18  # ±0.18 on per-contact trace score s(c)


def load_pair_data(pat: str, band: str) -> dict:
    npz_path = PAIR_SPLIT_DIR / f"{pat}_{band}.npz"
    d = np.load(npz_path)
    return dict(
        dD_task=np.asarray(d["dD_task"]),
        dD_rest=np.asarray(d["dD_rest"]),
        iu_i=np.asarray(d["iu_i"]).astype(int),
        iu_j=np.asarray(d["iu_j"]).astype(int),
    )


def load_implant(pat: str) -> pd.DataFrame:
    """Load (label, x, y, z) for each contact. MNI coords in mm."""
    nn = pat.replace("Pat_", "")
    df = pd.read_csv(SEEG_RAW / pat / f"implant_pat_{nn}.csv")
    return df[["label", "x", "y", "z"]].copy().assign(
        x=df["x"] / 1000.0,
        y=df["y"] / 1000.0,
        z=df["z"] / 1000.0,
    )


def per_contact_trace_score(pat: str, band: str) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-contact trace score s(c) for one patient × band.

    Returns
    -------
    scores : (N_contacts,) array of s(c) ∈ [−1, +1].
    coords : (N_contacts, 3) array of MNI coordinates (mm).
    """
    pair = load_pair_data(pat, band)
    sign_task = np.sign(pair["dD_task"])
    sign_rest = np.sign(pair["dD_rest"])
    concordance = sign_task * sign_rest

    iu_i = pair["iu_i"]
    iu_j = pair["iu_j"]
    n_contacts = int(max(iu_i.max(), iu_j.max())) + 1

    sums = np.zeros(n_contacts)
    counts = np.zeros(n_contacts)
    np.add.at(sums, iu_i, concordance)
    np.add.at(sums, iu_j, concordance)
    np.add.at(counts, iu_i, 1.0)
    np.add.at(counts, iu_j, 1.0)
    scores = sums / np.maximum(counts, 1.0)

    implant = load_implant(pat)
    coords = implant[["x", "y", "z"]].values[:n_contacts]

    return scores, coords


def aggregate_cohort(band: str) -> tuple[np.ndarray, np.ndarray]:
    all_scores, all_coords = [], []
    for pat in COHORT:
        s, c = per_contact_trace_score(pat, band)
        n = min(len(s), len(c))
        all_scores.append(s[:n])
        all_coords.append(c[:n])
    return (np.concatenate(all_scores),
            np.concatenate(all_coords, axis=0))


def plot_band_brain(ax, band: str, scores: np.ndarray, coords: np.ndarray,
                    cohort_row: pd.Series) -> None:
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
    rho = float(cohort_row["obs_median_rho"])
    p_val = float(cohort_row["paired_wilcoxon_p"])
    is_sig = p_val < 0.05

    plotting.plot_markers(
        node_values=scores,
        node_coords=coords,
        node_size=22,
        node_cmap=TRACE_CMAP,
        node_vmin=-SHARED_VABS,
        node_vmax=SHARED_VABS,
        display_mode="lzr",
        axes=ax,
        annotate=False,
        colorbar=False,
        alpha=0.88,
        node_threshold=None,
    )

    title_color = "#0a3a1d" if is_sig else "0.25"
    pos_pct = 100.0 * float(np.sum(scores > 0)) / max(len(scores), 1)
    ax.text(0.5, 0.97,
            band_tex, transform=ax.transAxes,
            ha="center", va="top",
            fontsize=18, fontweight="bold", color=title_color)
    ax.text(0.5, 0.04,
            rf"$\rho_{{\mathrm{{split}}}}^{{\mathrm{{coph}}}} = {rho:+.3f}$  |  "
            rf"$p = {p_val:.3f}$  |  "
            rf"{pos_pct:.0f}% green contacts ($n={len(scores)}$)",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=9.0, color="0.15")


def main() -> Path:
    cohort = pd.read_csv(LRG_CTM_DIR / "cohort_summary.csv")

    band_data = {}
    for band in BAND_ORDER:
        s, c = aggregate_cohort(band)
        band_data[band] = (s, c)
        n_pos = int(np.sum(s > 0))
        print(f"  {band}: {len(s)} contacts, {n_pos} green ({100*n_pos/len(s):.1f}%), "
              f"s range=[{s.min():+.3f}, {s.max():+.3f}]")

    fig = plt.figure(figsize=(18.0, 9.5))
    gs = fig.add_gridspec(
        2, 3, wspace=0.05, hspace=0.20,
        left=0.03, right=0.95, top=0.92, bottom=0.05,
    )

    for i, band in enumerate(BAND_ORDER):
        r, c = i // 3, i % 3
        ax = fig.add_subplot(gs[r, c])
        s, coords = band_data[band]
        cohort_row = cohort[cohort.band == band].iloc[0]
        plot_band_brain(ax, band, s, coords, cohort_row)

    cb_ax = fig.add_axes([0.96, 0.20, 0.012, 0.60])
    import matplotlib as mpl
    norm = mpl.colors.Normalize(vmin=-SHARED_VABS, vmax=SHARED_VABS)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=TRACE_CMAP)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cb_ax, orientation="vertical", extend="both")
    cb.set_label(
        r"per-contact trace score $s(c)$  "
        r"(green = pairs of this contact concord; red = discord)",
        fontsize=9, rotation=270, labelpad=18,
    )
    cb.ax.tick_params(labelsize=8)

    fig.text(0.5, 0.965,
             "Cohort-pooled per-contact trace contribution on glass-brain projections  "
             r"($n = 10$, ImCoh$|\cdot|$, lzry views)",
             ha="center", va="bottom", fontsize=11)
    fig.text(0.5, 0.945,
             r"each marker is a sEEG contact at its MNI coordinate; $s(c) = $ mean over "
             r"pairs $(c,j)$ of sign$(\Delta_t)\times$sign$(\Delta_r)$; "
             r"green = trace, red = anti",
             ha="center", va="bottom", fontsize=9, color="0.30")

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands" / "per_pair_trace"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_bands_brain_anatomy.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
