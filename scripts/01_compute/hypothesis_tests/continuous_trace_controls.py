#!/usr/bin/env python3
"""Continuous-trace matrix — three-control consolidation.

Combines the per-pair vectors already saved by ``continuous_trace_matrix.py``
(in both ``shared`` and ``split-baseline`` modes) with:

* H2e drift-only ρ (Run C, surfaced from
  ``data/reports/imcoh_vi/h2e_split_half_rho_raw.csv``):
  ``ρ_null_drift = Spearman(D^pre_B − D^pre_A, D^post_B − D^post_A)``.
  Same noise regime as Run A — provides the within-session-only floor
  for the split-baseline ρ.

* Probe-bias split (Run B): re-evaluates Run A's per-(patient, band) ρ
  on (i) cross-probe pairs only and (ii) same-probe pairs only.
  Uses ``lrg_eegfc.utils.probe.build_probe_mask`` on the canonical
  ``channel_labels.csv``. Same-probe pairs carry trivial volume-conduction-style
  inflation; cross-probe pairs carry the load-bearing cohort signal.

Outputs:
    data/reports/imcoh_continuous_trace/controls_summary.csv
        per (patient, band): rho_shared, rho_split, rho_null_drift,
        rho_split_same_probe, rho_split_cross_probe,
        n_pairs_total, n_pairs_same, n_pairs_cross
    data/reports/imcoh_continuous_trace/controls_band_stats.md
        per band: cohort medians, IQR, Wilcoxon paired tests vs null_drift
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE,
    PATIENT_CHANNEL_DROP,
)
from lrg_eegfc.config.paths import REPORTS_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.probe import build_probe_mask
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, rank_biserial


CT_DIR = REPORTS_ROOT / "imcoh_continuous_trace"
H2E_RAW = REPORTS_ROOT / "imcoh_vi" / "h2e_split_half_rho_raw.csv"


def load_channel_labels(pat: str) -> list[str] | None:
    """Read patient's channel_labels.csv and apply canonical row drops."""
    fp = SEEG_DATAPATH / pat / "channel_labels.csv"
    if not fp.exists():
        return None
    df = pd.read_csv(fp)
    label_drops = PATIENT_CHANNEL_DROP.get(pat, {}).get("__labels__", [])
    if label_drops:
        df = df.drop(index=list(label_drops)).reset_index(drop=True)
    return df["label"].astype(str).tolist()


def probe_split_rho(npz_path: Path,
                    same_probe_mask: np.ndarray) -> tuple[float, float, int, int]:
    """Recompute ρ on same-probe vs cross-probe upper-triangle pairs.

    Returns (rho_same, rho_cross, n_same, n_cross). NaN where degenerate.
    """
    pp = np.load(npz_path)
    iu_i = pp["iu_i"]; iu_j = pp["iu_j"]
    dT = pp["dD_task"]; dR = pp["dD_rest"]
    sp = same_probe_mask[iu_i, iu_j]
    cp = ~sp

    def _rho(mask):
        if mask.sum() < 5:
            return np.nan
        x = dT[mask]; y = dR[mask]
        if x.std() == 0 or y.std() == 0:
            return np.nan
        r = spearmanr(x, y).statistic
        return float(r) if np.isfinite(r) else np.nan

    return _rho(sp), _rho(cp), int(sp.sum()), int(cp.sum())


def collect(patients: list[str], bands: list[str]) -> pd.DataFrame:
    summary_shared = pd.read_csv(CT_DIR / "per_cell_summary.csv")
    summary_split = pd.read_csv(CT_DIR / "per_cell_summary_split.csv")
    h2e = pd.read_csv(H2E_RAW)[["patient", "band", "rho_null_drift"]]

    pair_dir_shared = CT_DIR / "per_pair"
    pair_dir_split = CT_DIR / "per_pair_split"

    rows: list[dict] = []
    for pat in patients:
        labels = load_channel_labels(pat)
        if labels is None:
            print(f"[ctrl] WARN {pat}: missing channel_labels.csv; skipping probe split")
            sp_mask = None
        else:
            sp_mask = build_probe_mask(labels)
        for band in bands:
            row_shared = summary_shared[
                (summary_shared["patient"] == pat) & (summary_shared["band"] == band)
            ]
            row_split = summary_split[
                (summary_split["patient"] == pat) & (summary_split["band"] == band)
            ]
            row_h2e = h2e[(h2e["patient"] == pat) & (h2e["band"] == band)]
            if row_shared.empty or row_split.empty:
                continue
            rho_sh = float(row_shared["rho"].iloc[0])
            rho_sp = float(row_split["rho"].iloc[0])
            rho_nd = (float(row_h2e["rho_null_drift"].iloc[0])
                      if not row_h2e.empty else np.nan)

            rho_same, rho_cross, n_same, n_cross = (np.nan, np.nan, 0, 0)
            if sp_mask is not None:
                npz_path = pair_dir_split / f"{pat}_{band}.npz"
                if npz_path.exists():
                    N = sp_mask.shape[0]
                    if int(row_split["N_p"].iloc[0]) == N:
                        rho_same, rho_cross, n_same, n_cross = probe_split_rho(
                            npz_path, sp_mask
                        )

            rows.append({
                "patient": pat, "band": band,
                "rho_shared": rho_sh,
                "rho_split": rho_sp,
                "rho_null_drift": rho_nd,
                "rho_split_same_probe": rho_same,
                "rho_split_cross_probe": rho_cross,
                "n_pairs_same": n_same,
                "n_pairs_cross": n_cross,
            })
    return pd.DataFrame(rows)


def band_stats(df: pd.DataFrame, bands: list[str]) -> list[dict]:
    """Per-band cohort summary and paired Wilcoxon ρ_split vs ρ_null_drift."""
    out = []
    for band in bands:
        sub = df[df["band"] == band].copy()
        n = len(sub)
        d = sub.dropna(subset=["rho_split", "rho_null_drift"])
        if len(d) >= 3:
            diff = (d["rho_split"] - d["rho_null_drift"]).to_numpy()
            z, p = wilcoxon_z(diff)
            r_rb = rank_biserial(diff)
        else:
            z, p, r_rb = np.nan, np.nan, np.nan
        cross = sub["rho_split_cross_probe"].dropna()
        same = sub["rho_split_same_probe"].dropna()
        out.append({
            "band": band, "n": n,
            "med_shared": sub["rho_shared"].median(),
            "med_split":  sub["rho_split"].median(),
            "med_null_drift": sub["rho_null_drift"].median(),
            "med_cross":  cross.median() if len(cross) else np.nan,
            "med_same":   same.median()  if len(same)  else np.nan,
            "n_pos_shared": int((sub["rho_shared"] > 0).sum()),
            "n_pos_split":  int((sub["rho_split"]  > 0).sum()),
            "n_pos_cross":  int((cross > 0).sum()) if len(cross) else 0,
            "n_pos_same":   int((same  > 0).sum()) if len(same)  else 0,
            "n_split_above_drift": int((d["rho_split"] > d["rho_null_drift"]).sum())
                                      if len(d) else 0,
            "n_paired_drift": len(d),
            "wilcoxon_z": z, "wilcoxon_p": p, "r_rb": r_rb,
        })
    return out


def write_md(stats_rows: list[dict], outpath: Path) -> None:
    tex = BRAIN_BAND_TEX_DICT
    lines = [
        "# Continuous-trace controls — three-test consolidation",
        "",
        "Cohort N=10 (N=9 for `null_drift` — Pat_14 not in H2e).",
        "",
        "## Per-band cohort summary (median ρ)",
        "",
        "| band | shared | split (Run A) | null_drift (Run C) | cross-probe (Run B) | same-probe (Run B) |",
        "|------|-------:|--------------:|-------------------:|--------------------:|-------------------:|",
    ]
    for r in stats_rows:
        lines.append(
            f"| {tex[r['band']]} | {r['med_shared']:+.3f} | {r['med_split']:+.3f} | "
            f"{r['med_null_drift']:+.3f} | {r['med_cross']:+.3f} | {r['med_same']:+.3f} |"
        )
    lines += [
        "",
        "## Cohort sign-test counts (patients with ρ > 0)",
        "",
        "| band | shared | split | cross-probe | same-probe | split>drift |",
        "|------|-------:|------:|------------:|-----------:|------------:|",
    ]
    for r in stats_rows:
        lines.append(
            f"| {tex[r['band']]} | {r['n_pos_shared']}/{r['n']} | "
            f"{r['n_pos_split']}/{r['n']} | "
            f"{r['n_pos_cross']}/{r['n']} | {r['n_pos_same']}/{r['n']} | "
            f"{r['n_split_above_drift']}/{r['n_paired_drift']} |"
        )
    lines += [
        "",
        "## Paired Wilcoxon: ρ_split (Run A) > ρ_null_drift (Run C)?",
        "",
        "Same noise regime → fair comparison. Tests whether the residual "
        "split-baseline ρ exceeds the within-session drift floor.",
        "",
        "| band | n_paired | z | p | r_rb |",
        "|------|---------:|--:|--:|-----:|",
    ]
    for r in stats_rows:
        lines.append(
            f"| {tex[r['band']]} | {r['n_paired_drift']} | "
            f"{r['wilcoxon_z']:+.2f} | {r['wilcoxon_p']:.4f} | "
            f"{r['r_rb']:+.3f} |"
        )
    lines += [
        "",
        "## Reading",
        "",
        "* **shared** is the original `H2c`-style ρ — sharing D_pre across "
        "  Δ_task and Δ_rest. Vulnerable to spurious correlation.",
        "* **split (Run A)** uses independent half-baselines D_pre_A ⊥ D_pre_B. "
        "  Removes the shared-baseline-noise contribution. Lives in halved-"
        "  data noise regime.",
        "* **null_drift (Run C, from H2e)** uses only halves of rest_pre and "
        "  rest_post (no task), `Spearman(D^pre_B − D^pre_A, D^post_B − D^post_A)`. "
        "  Same halved-data noise regime as Run A. Captures session-drift only.",
        "* **cross-probe / same-probe (Run B)** restrict the Run A pair set. "
        "  Same-probe pairs are anatomy-dominated; cross-probe carry the "
        "  load-bearing volume-conduction-immune signal.",
        "",
        "Pass criteria after controls:",
        "1. `split` cohort median > 0 (residual genuine signal).",
        "2. `split > null_drift` in ≥ 7 of 10 patients (above noise floor).",
        "3. `cross-probe` cohort median > 0 (not anatomy-bias artefact).",
        "",
    ]
    outpath.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", default=None)
    ap.add_argument("--bands", default=None)
    args = ap.parse_args()

    patients = (args.patients.split(",") if args.patients
                else list(PATIENTS_4PHASE))
    bands = args.bands.split(",") if args.bands else list(BRAIN_BANDS_NAMES)

    df = collect(patients, bands)
    csv_out = CT_DIR / "controls_summary.csv"
    df.to_csv(csv_out, index=False)
    print(f"[ctrl] wrote {csv_out}  ({len(df)} rows)")

    stats = band_stats(df, bands)
    md_out = CT_DIR / "controls_band_stats.md"
    write_md(stats, md_out)
    print(f"[ctrl] wrote {md_out}")

    print()
    print(f"{'band':<10s}  {'shared':>7s}  {'split':>7s}  {'drift':>7s}  "
          f"{'cross':>7s}  {'same':>7s}  {'split>drift':>11s}")
    for r in stats:
        print(f"{r['band']:<10s}  {r['med_shared']:>+7.3f}  "
              f"{r['med_split']:>+7.3f}  {r['med_null_drift']:>+7.3f}  "
              f"{r['med_cross']:>+7.3f}  {r['med_same']:>+7.3f}  "
              f"{r['n_split_above_drift']:>3d}/{r['n_paired_drift']:<3d}")


if __name__ == "__main__":
    main()
