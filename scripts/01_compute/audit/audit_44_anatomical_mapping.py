#!/usr/bin/env python3
"""Audit 44 -- anatomical mapping for calibrated trace-leaves.

Maps the leaf indices flagged is_trace_p95 in measure 12 to:
  - sEEG contact label (via channel_labels.csv row index = leaf_id)
  - 3D MNI coordinates (via implant_pat_NN.csv)
  - Desikan-Killany dominant region (parsed from implant CSV's
    "Desikan-Killany" column: format
    "region1,pct1,region2,pct2,...,PTD,confidence" where region1 is
    the highest-weight parcellation region).

Outputs
-------
data/audit/lrg_localization_anatomy/per_trace_leaf.csv
    per-(patient, band, leaf) trace-leaf row with channel + coordinates
    + dominant region + pct + parcellation confidence.
data/audit/lrg_localization_anatomy/cohort_band_region_frequency.csv
    per-(band, region) count and patient-coverage summary.
data/outputs/figures/section_5_lrg_trace/anatomy/region_frequency_grid.pdf
    grid of per-band region-frequency bar charts (top regions).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

LEAF_CAL = ROOT / "data" / "audit" / "per_leaf_rho_null" / "calibrated_trace_leaves.csv"
RAW_ROOT = ROOT / "data" / "raw" / "stereoeeg_patients"
OUT_DIR = ROOT / "data" / "audit" / "lrg_localization_anatomy"
FIG_DIR = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "anatomy"

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}

NON_PARCEL = {"Wm", "Unk", "PTD"}


def normalize_label(raw: str) -> str:
    """Channel label cleanup: 'A 1,G2' -> 'A1' (strip space + drop group)."""
    main = raw.split(",")[0]
    return main.replace(" ", "").upper()


def parse_dk_dominant(dk: str) -> tuple[str, float, float]:
    """Parse 'region1,pct1,region2,pct2,...,PTD,confidence' -> (region1, pct1, conf).

    Returns (NaN-coded) ('Unk', np.nan, np.nan) on parse failure.
    """
    if pd.isna(dk):
        return "Unk", float("nan"), float("nan")
    parts = [p.strip() for p in str(dk).split(",")]
    region1 = parts[0] if parts else "Unk"
    pct1 = float("nan")
    if len(parts) > 1:
        try:
            pct1 = float(parts[1])
        except ValueError:
            pct1 = float("nan")
    # PTD confidence is the trailing pair "PTD, <num>"
    conf = float("nan")
    for k in range(len(parts) - 1):
        if parts[k] == "PTD":
            try:
                conf = float(parts[k + 1])
            except ValueError:
                conf = float("nan")
            break
    return region1, pct1, conf


def load_channel_to_anatomy(pat: str) -> pd.DataFrame:
    """Build leaf_id -> {channel, x, y, z, region, pct, conf} table for one patient."""
    labels_csv = RAW_ROOT / pat / "channel_labels.csv"
    impl_csv = RAW_ROOT / pat / f"implant_pat_{pat[-2:]}.csv"
    labels = pd.read_csv(labels_csv)
    impl = pd.read_csv(impl_csv)

    # Detect the Desikan column (varies "Desikan-Killany" / "Desikan-K")
    dk_col = next((c for c in impl.columns if c.startswith("Desikan")), None)
    if dk_col is None:
        raise RuntimeError(f"{pat}: no Desikan column in {impl_csv}")

    # Build a normalized-label -> implant row map
    impl["norm_label"] = impl["label"].astype(str).str.replace(" ", "").str.upper()
    impl_map = impl.set_index("norm_label")

    rows = []
    for leaf_id, raw in enumerate(labels["label"].astype(str).tolist()):
        ch = normalize_label(raw)
        rec = {"patient": pat, "leaf_id": leaf_id, "channel_raw": raw,
               "channel": ch, "x": np.nan, "y": np.nan, "z": np.nan,
               "region": "Unk", "region_pct": np.nan, "ptd_conf": np.nan}
        if ch in impl_map.index:
            r = impl_map.loc[ch]
            if isinstance(r, pd.DataFrame):  # duplicate key
                r = r.iloc[0]
            rec["x"] = float(r["x"])
            rec["y"] = float(r["y"])
            rec["z"] = float(r["z"])
            region, pct, conf = parse_dk_dominant(r[dk_col])
            rec["region"] = region
            rec["region_pct"] = pct
            rec["ptd_conf"] = conf
        rows.append(rec)
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    leaf_cal = pd.read_csv(LEAF_CAL)
    # patient-level anatomy maps
    per_patient = {pat: load_channel_to_anatomy(pat) for pat in PATIENTS}

    # Trace-leaf rows only
    trace_rows = []
    for pat in PATIENTS:
        anat = per_patient[pat].set_index("leaf_id")
        sub = leaf_cal[(leaf_cal.patient == pat) & (leaf_cal.is_trace_p95 == 1)]
        for _, r in sub.iterrows():
            ell = int(r["leaf_id"])
            if ell not in anat.index:
                continue
            a = anat.loc[ell]
            trace_rows.append({
                "patient": pat, "band": r["band"], "leaf_id": ell,
                "channel": a["channel"], "channel_raw": a["channel_raw"],
                "x": a["x"], "y": a["y"], "z": a["z"],
                "region": a["region"], "region_pct": a["region_pct"],
                "ptd_conf": a["ptd_conf"],
                "rho_demeaned": float(r["rho_demeaned"]),
                "tau_95": float(r["tau_95"]),
            })
    trace_df = pd.DataFrame(trace_rows)
    trace_df.to_csv(OUT_DIR / "per_trace_leaf.csv", index=False)
    print(f"[44] wrote {OUT_DIR / 'per_trace_leaf.csv'} ({len(trace_df)} rows)")

    # Cohort band-region frequency: count per (band, region) and number of patients
    band_freq = (
        trace_df.groupby(["band", "region"])
        .agg(
            n_leaves=("leaf_id", "count"),
            n_patients=("patient", "nunique"),
            mean_rho_demeaned=("rho_demeaned", "mean"),
            mean_pct=("region_pct", "mean"),
            mean_conf=("ptd_conf", "mean"),
        )
        .reset_index()
        .sort_values(["band", "n_leaves"], ascending=[True, False])
    )
    band_freq.to_csv(OUT_DIR / "cohort_band_region_frequency.csv", index=False)
    print(f"[44] wrote {OUT_DIR / 'cohort_band_region_frequency.csv'}")

    # Print summary: top 10 regions per band
    print("\n=== Top regions per band (calibrated trace-leaves cohort-pooled) ===")
    for band in BANDS:
        sub = band_freq[band_freq.band == band].head(10)
        if sub.empty:
            print(f"\n[{band}] no trace-leaves")
            continue
        print(f"\n[{band}] (total trace-leaves={int(sub.n_leaves.sum())} "
              f"across {sub.n_patients.max()} patients)")
        for _, r in sub.iterrows():
            print(f"  {r.region:32s}  n_leaves={int(r.n_leaves):3d}  "
                  f"n_pat={int(r.n_patients)}  mean_pct={r.mean_pct:.0f}%  "
                  f"mean_conf={r.mean_conf:.2f}")

    # ---------------- figure: per-band top-region bar chart -----------------
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    for ax, band in zip(axes.ravel(), BANDS):
        sub = band_freq[band_freq.band == band].head(8).iloc[::-1]
        if sub.empty:
            ax.text(0.5, 0.5, "no trace-leaves", ha="center", va="center",
                    transform=ax.transAxes, fontsize=10)
            ax.set_title(BAND_TEX[band], fontsize=12)
            ax.set_xticks([]); ax.set_yticks([])
            continue
        # Color by hemisphere/structure type
        colors = []
        for r in sub["region"]:
            if r.startswith("ctx-lh"):
                colors.append("#1f77b4")  # blue = left cortex
            elif r.startswith("ctx-rh"):
                colors.append("#d62728")  # red = right cortex
            elif r in ("Wm", "Unk"):
                colors.append("#999999")
            elif "Cerebellum" in r:
                colors.append("#2ca02c")
            elif r in ("Left-Hippocampus", "Right-Hippocampus",
                       "Left-Amygdala", "Right-Amygdala",
                       "Left-Thalamus-Proper", "Right-Thalamus-Proper"):
                colors.append("#9467bd")  # purple = subcortical
            else:
                colors.append("#7f7f7f")
        y = np.arange(len(sub))
        ax.barh(y, sub["n_leaves"].to_numpy(), color=colors, edgecolor="0.3", lw=0.5)
        ax.set_yticks(y)
        ax.set_yticklabels(sub["region"].tolist(), fontsize=8)
        for i, (n, npat) in enumerate(zip(sub["n_leaves"], sub["n_patients"])):
            ax.text(n + 0.1, i, f"{int(n)} ({int(npat)} pat)", va="center", fontsize=7)
        ax.set_xlabel("trace-leaves", fontsize=9)
        ax.set_title(BAND_TEX[band], fontsize=12)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "region_frequency_grid.pdf")
    plt.close(fig)
    print(f"\n[44] wrote {FIG_DIR / 'region_frequency_grid.pdf'}")


if __name__ == "__main__":
    main()
