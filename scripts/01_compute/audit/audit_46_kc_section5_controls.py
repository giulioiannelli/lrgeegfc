#!/usr/bin/env python3
"""Audit 46 — Section 5.2 KC controls: cross-probe restriction + Pat_03 dropout
extensions + within-probe BH-FDR. Mirrors the CTM control regime of section 5.3.

Three control layers per (patient, band, lambda):

    1. Within-baseline null (DONE, audit_41) — read from
       ``data/audit/lrg_global_probe_controls/cohort_controls_summary.csv``;
       this script only ingests the existing numbers, does not recompute.

    2. Cross-probe KC — for each phase pair, compute
       ``d_KC_xprobe(lam; Z_a, Z_b) = || ((1-lam)*m_norm + lam*M_norm)|_xprobe ||_2``
       where the restriction picks only leaf-pairs (i, j) with different
       sEEG probes; m, M are the standard Kendall-Colijn (m, M) per leaf
       pair, normalized to [0, 1] using the pooled-max convention of
       ``lrg_eegfc.utils.metrics.tree_distance.kc_distance``. Triangle:
       ``T_KC_xprobe(lam; p, b) = d_KC_xprobe(lam; Z^TT, Z^RPost)
                                 - d_KC_xprobe(lam; Z^RPre, Z^TT)``.
       Cohort one-sided Wilcoxon T_KC_xprobe > 0 on n=10 patients.

    3. Pat_03 dropout — extend the existing partial table to all (band,
       lambda) cells under three flavors: absolute Wilcoxon T_KC > 0,
       within-baseline-null Wilcoxon T_KC > T_KC_null, cross-probe
       Wilcoxon T_KC_xprobe > 0.

Within-probe BH-FDR at m=12 (6 bands × 2 lambda values, λ=0 and λ=1) on
each control layer.

Outputs
-------
- ``data/audit/section5_v2_kc_controls/per_patient_table.csv`` — 60 rows
  (10 patients × 6 bands) × {lam=0.0, lam=0.5, lam=1.0} with real T_KC,
  null T_KC (from audit_41), xprobe T_KC, and per-row indicator booleans
  for each control's "below-floor" event.
- ``data/audit/section5_v2_kc_controls/cohort_summary.csv`` — 18 rows
  (6 bands × 3 lambda) with absolute / within-baseline / cross-probe
  cohort Wilcoxon p, n_above, rank-biserial, BH-corrected q (m=12 for
  λ=0/λ=1, λ=0.5 reported uncorrected).
- ``data/audit/section5_v2_kc_controls/pat03_dropout.csv`` — 18 rows
  with p_full and p_drop for all three control layers.
- ``data/outputs/figures/section_5_lrg_trace/lrg_controls/kc_controls_summary.pdf``
  — companion figure: per-band cohort medians for absolute / null /
  xprobe alongside per-band rank-biserial bars; β / γ_l highlighted.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.patient import parse_seeg_label
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr, rank_biserial
from lrg_eegfc.utils.metrics.tree_distance import kc_vectors
from lrg_eegfc.workflow.lrg import load_lrg_result

PATIENTS = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
BANDS = list(BRAIN_BANDS_NAMES)
KC_LAMBDAS = [0.0, 0.5, 1.0]
PHASES = ("rest_pre", "task_test", "rest_post")
SEEG_ROOT = ROOT / "data" / "raw" / "stereoeeg_patients"

OUT = ROOT / "data" / "audit" / "section5_v2_kc_controls"
FIG_DIR = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "lrg_controls"

EXISTING_NULL = ROOT / "data" / "audit" / "lrg_global_probe_controls" / "kc_null.csv"
EXISTING_REAL = (
    ROOT
    / "data"
    / "reports"
    / "section_5_lrg_trace"
    / "03_kc_lambda_triangle"
    / "tables"
    / "Td_per_patient_per_band_lambda.csv"
)


# ---------------------------------------------------------------------------
# Probe assignment per leaf
# ---------------------------------------------------------------------------

def load_probe_ids(patient: str) -> np.ndarray:
    """Return per-leaf probe identifier (string) of length n_leaves.

    Uses ``channel_labels.csv`` row order, which matches the LRG leaf
    indices (verified n_nodes == channel_labels rows for all 10 patients
    of the current cohort).
    """
    csv = SEEG_ROOT / patient / "channel_labels.csv"
    df = pd.read_csv(csv)
    labels = df.iloc[:, 0].astype(str).tolist()
    probes = []
    for lab in labels:
        probe, _ = parse_seeg_label(lab)
        probes.append(probe if probe is not None else "_unknown")
    return np.asarray(probes, dtype=object)


def cross_probe_pair_mask(probe_ids: np.ndarray) -> np.ndarray:
    """Boolean mask over scipy condensed pair order, True for cross-probe."""
    n = len(probe_ids)
    mask = np.zeros(n * (n - 1) // 2, dtype=bool)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            mask[idx] = probe_ids[i] != probe_ids[j]
            idx += 1
    return mask


# ---------------------------------------------------------------------------
# Cross-probe KC distance
# ---------------------------------------------------------------------------

def kc_xprobe_distance(
    Z1: np.ndarray, Z2: np.ndarray, lam: float, mask: np.ndarray
) -> float:
    """Cross-probe-restricted KC distance with the standard pooled-max norm.

    Mirrors ``lrg_eegfc.utils.metrics.tree_distance.kc_distance`` exactly,
    except the L2 norm is taken over the cross-probe subset only.
    """
    m1, M1 = kc_vectors(Z1)
    m2, M2 = kc_vectors(Z2)
    mmax = max(int(m1.max()), int(m2.max()), 1)
    Mmax = max(float(M1.max()), float(M2.max()), 1e-12)
    m1n, m2n = m1 / mmax, m2 / mmax
    M1n, M2n = M1 / Mmax, M2 / Mmax
    v1 = (1 - lam) * m1n + lam * M1n
    v2 = (1 - lam) * m2n + lam * M2n
    diff = (v1 - v2)[mask]
    return float(np.linalg.norm(diff))


def xprobe_triangle_for(patient: str, band: str, mask: np.ndarray) -> dict:
    """Compute cross-probe T_KC for the three lambdas; return None on any miss."""
    Z = {}
    for phi in PHASES:
        res = load_lrg_result(patient, phi, band, fc_method="imcoh_abs")
        if res is None:
            return None
        Z[phi] = res.linkage_matrix
    out = {}
    for lam in KC_LAMBDAS:
        d_pre_tt = kc_xprobe_distance(Z["rest_pre"], Z["task_test"], lam, mask)
        d_tt_post = kc_xprobe_distance(Z["task_test"], Z["rest_post"], lam, mask)
        # T_d > 0 = trace (rsPost closer to task than rsPre).
        out[f"xprobe_T_lam{lam:.1f}"] = d_pre_tt - d_tt_post
        out[f"xprobe_d_pre_tt_lam{lam:.1f}"] = d_pre_tt
        out[f"xprobe_d_tt_post_lam{lam:.1f}"] = d_tt_post
    return out


# ---------------------------------------------------------------------------
# Wilcoxon one-sided p
# ---------------------------------------------------------------------------

def wilcoxon_one_sided_greater(x: np.ndarray, y: np.ndarray | None = None) -> float:
    """Paired one-sided Wilcoxon (alternative='greater'); returns nan on degenerate input.

    Tests x > 0 (or x > y) in the trace direction under the project-wide
    sign convention T_d > 0 = trace.
    """
    try:
        if y is None:
            return float(wilcoxon(x, alternative="greater", zero_method="wilcox").pvalue)
        return float(wilcoxon(x, y, alternative="greater", zero_method="wilcox").pvalue)
    except ValueError:
        return float("nan")


def rb_paired_greater(x: np.ndarray, y: np.ndarray | None = None) -> float:
    """Rank-biserial in [-1, +1] for a one-sided 'x > y' (or 'x > 0') Wilcoxon."""
    try:
        diff = x if y is None else (np.asarray(x) - np.asarray(y))
        diff = diff[diff != 0]
        if len(diff) == 0:
            return float("nan")
        return float(rank_biserial(diff))
    except Exception:
        return float("nan")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # ---- 1. cross-probe T_KC per (patient, band, lambda) ----
    print("[audit_46] computing cross-probe T_KC ...")
    rows_xpr = []
    for pat in PATIENTS:
        try:
            probe_ids = load_probe_ids(pat)
            mask = cross_probe_pair_mask(probe_ids)
            n_pairs = len(mask)
            n_xpr = int(mask.sum())
            print(f"  {pat}: n_leaves={len(probe_ids)}, n_pairs={n_pairs}, "
                  f"n_xprobe_pairs={n_xpr} ({100 * n_xpr / n_pairs:.1f}%)")
        except Exception as e:
            print(f"  {pat}: WARN failed probe assignment: {e}")
            continue
        for band in BANDS:
            try:
                tri = xprobe_triangle_for(pat, band, mask)
                if tri is None:
                    print(f"    {band}: missing LRG cache, skip")
                    continue
                row = dict(patient=pat, band=band,
                           n_leaves=len(probe_ids),
                           n_xprobe_pairs=n_xpr)
                row.update(tri)
                rows_xpr.append(row)
            except Exception as e:
                print(f"    {band}: WARN {e}")
    df_xpr = pd.DataFrame(rows_xpr)

    # ---- 2. read existing real T_KC and within-baseline null ----
    real_df = pd.read_csv(EXISTING_REAL)
    real_df = real_df[real_df["band"].isin(BANDS)].copy()
    null_df = pd.read_csv(EXISTING_NULL)

    # ---- 3. assemble per-patient long table ----
    pieces = []
    for lam in KC_LAMBDAS:
        sub_real = real_df[real_df["lam"] == lam][["patient", "band", "T_KC"]].rename(
            columns={"T_KC": "real_T_KC"}
        )
        null_col = f"null_T_lam{lam:.1f}"
        sub_null = null_df[["patient", "band", null_col]].rename(
            columns={null_col: "null_T_KC"}
        )
        sub_xpr = df_xpr[["patient", "band", f"xprobe_T_lam{lam:.1f}"]].rename(
            columns={f"xprobe_T_lam{lam:.1f}": "xprobe_T_KC"}
        )
        merged = sub_real.merge(sub_null, on=["patient", "band"]).merge(
            sub_xpr, on=["patient", "band"]
        )
        merged["lam"] = lam
        merged["abs_above_zero"] = merged["real_T_KC"] > 0
        merged["real_above_null"] = merged["real_T_KC"] > merged["null_T_KC"]
        merged["xprobe_above_zero"] = merged["xprobe_T_KC"] > 0
        pieces.append(merged)
    per_patient = pd.concat(pieces, ignore_index=True)
    per_patient = per_patient[
        ["patient", "band", "lam", "real_T_KC", "null_T_KC", "xprobe_T_KC",
         "abs_above_zero", "real_above_null", "xprobe_above_zero"]
    ]
    per_patient.to_csv(OUT / "per_patient_table.csv", index=False)

    # ---- 4. cohort summary per (band, lambda, control) ----
    coh_rows = []
    for band in BANDS:
        for lam in KC_LAMBDAS:
            sub = per_patient[(per_patient["band"] == band) & (per_patient["lam"] == lam)]
            real = sub["real_T_KC"].values.astype(float)
            null = sub["null_T_KC"].values.astype(float)
            xpr = sub["xprobe_T_KC"].values.astype(float)
            n = len(sub)
            row = dict(
                band=band, lam=lam, n_pat=n,
                # absolute Wilcoxon (T_KC > 0)
                n_abs_above_zero=int(np.sum(real > 0)),
                median_real=float(np.median(real)),
                p_abs=wilcoxon_one_sided_greater(real),
                rb_abs=rb_paired_greater(real),
                # within-baseline null (T_KC > T_KC_null)
                n_real_above_null=int(np.sum(real > null)),
                median_null=float(np.median(null)),
                p_null=wilcoxon_one_sided_greater(real, null),
                rb_null=rb_paired_greater(real, null),
                # cross-probe (T_KC_xprobe > 0)
                n_xprobe_above_zero=int(np.sum(xpr > 0)),
                median_xprobe=float(np.median(xpr)),
                p_xprobe=wilcoxon_one_sided_greater(xpr),
                rb_xprobe=rb_paired_greater(xpr),
            )
            coh_rows.append(row)
    coh = pd.DataFrame(coh_rows)

    # within-probe BH-FDR at m=12 over (6 bands × 2 lambdas), separately per layer
    coh["q_abs_within_m12"] = np.nan
    coh["q_null_within_m12"] = np.nan
    coh["q_xprobe_within_m12"] = np.nan
    mask12 = coh["lam"].isin([0.0, 1.0])
    for col, qcol in [
        ("p_abs", "q_abs_within_m12"),
        ("p_null", "q_null_within_m12"),
        ("p_xprobe", "q_xprobe_within_m12"),
    ]:
        p = coh.loc[mask12, col].values
        q = bh_fdr(p)
        coh.loc[mask12, qcol] = q
    coh.to_csv(OUT / "cohort_summary.csv", index=False)

    # ---- 5. Pat_03 dropout under all three layers ----
    drop_rows = []
    full = per_patient.copy()
    drop = per_patient[per_patient["patient"] != "Pat_03"].copy()
    for band in BANDS:
        for lam in KC_LAMBDAS:
            f = full[(full["band"] == band) & (full["lam"] == lam)]
            d = drop[(drop["band"] == band) & (drop["lam"] == lam)]
            row = dict(
                band=band, lam=lam,
                n_full=len(f), n_drop=len(d),
                # absolute
                p_abs_full=wilcoxon_one_sided_greater(f["real_T_KC"].values),
                p_abs_drop=wilcoxon_one_sided_greater(d["real_T_KC"].values),
                n_abs_above_zero_full=int((f["real_T_KC"] > 0).sum()),
                n_abs_above_zero_drop=int((d["real_T_KC"] > 0).sum()),
                # within-baseline null
                p_null_full=wilcoxon_one_sided_greater(f["real_T_KC"].values,
                                                    f["null_T_KC"].values),
                p_null_drop=wilcoxon_one_sided_greater(d["real_T_KC"].values,
                                                    d["null_T_KC"].values),
                n_real_above_null_full=int((f["real_T_KC"] > f["null_T_KC"]).sum()),
                n_real_above_null_drop=int((d["real_T_KC"] > d["null_T_KC"]).sum()),
                # cross-probe
                p_xprobe_full=wilcoxon_one_sided_greater(f["xprobe_T_KC"].values),
                p_xprobe_drop=wilcoxon_one_sided_greater(d["xprobe_T_KC"].values),
                n_xprobe_above_zero_full=int((f["xprobe_T_KC"] > 0).sum()),
                n_xprobe_above_zero_drop=int((d["xprobe_T_KC"] > 0).sum()),
            )
            drop_rows.append(row)
    drop_df = pd.DataFrame(drop_rows)
    drop_df.to_csv(OUT / "pat03_dropout.csv", index=False)

    # ---- 6. companion figure (compact summary) ----
    layers = [("p_abs", "absolute"), ("p_null", "within-baseline null"), ("p_xprobe", "cross-probe")]
    fig, axes = plt.subplots(2, 1, figsize=(8.4, 6.2), sharex=True)
    band_x = np.arange(len(BANDS))
    width = 0.21
    colors = {"absolute": "#888", "within-baseline null": "#1f77b4", "cross-probe": "#e07b00"}

    # top: -log10 p at lambda=0
    ax = axes[0]
    for k, (pcol, lab) in enumerate(layers):
        sub = coh[coh["lam"] == 0.0]
        vals = -np.log10(np.maximum(sub[pcol].values, 1e-4))
        ax.bar(band_x + (k - 1) * width, vals, width=width, color=colors[lab],
               edgecolor="black", linewidth=0.5, label=lab)
    ax.set_xticks(band_x)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS])
    ax.axhline(-np.log10(0.05), color="0.4", lw=0.7, ls=":")
    ax.set_ylabel(r"$-\log_{10} p$  ($\lambda=0$, topology)")
    ax.legend(loc="upper left", frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    # bottom: -log10 p at lambda=1
    ax = axes[1]
    for k, (pcol, lab) in enumerate(layers):
        sub = coh[coh["lam"] == 1.0]
        vals = -np.log10(np.maximum(sub[pcol].values, 1e-4))
        ax.bar(band_x + (k - 1) * width, vals, width=width, color=colors[lab],
               edgecolor="black", linewidth=0.5, label=lab)
    ax.set_xticks(band_x)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS])
    ax.axhline(-np.log10(0.05), color="0.4", lw=0.7, ls=":")
    ax.set_ylabel(r"$-\log_{10} p$  ($\lambda=1$, heights)")
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "kc_controls_summary.pdf")
    plt.close(fig)

    # ---- 7. console summary ----
    print("\n=== cohort summary (sample) ===")
    print(coh[coh["lam"].isin([0.0, 1.0])][[
        "band", "lam", "n_pat",
        "n_abs_above_zero", "p_abs", "q_abs_within_m12",
        "n_real_above_null", "p_null", "q_null_within_m12",
        "n_xprobe_above_zero", "p_xprobe", "q_xprobe_within_m12",
    ]].to_string(index=False))
    print(f"\n[audit_46] outputs at {OUT}\n         figure at {FIG_DIR / 'kc_controls_summary.pdf'}")


if __name__ == "__main__":
    main()
