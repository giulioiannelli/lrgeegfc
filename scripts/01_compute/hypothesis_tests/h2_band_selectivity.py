#!/usr/bin/env python3
"""Band-selectivity test for H2 sub-hypotheses.

For each of the three H2 operationalisations (H2a VI magnitude, H2c ultrametric
drift ρ, H2d block-persistence Δρ) we ask a single question:

    Are the 6 bands behaving the same, or are some bands different?

Test: Friedman χ²(5) on per-patient vectors of 6 band values.
Post-hoc: paired Wilcoxon for every band pair, Bonferroni-corrected (m=15).
Focused contrast: θ vs each other band (m=5), one-sided (θ smaller).

This formalises "θ is the ergodic exception" as a statistical claim and gives
a principled non-ergodicity test for every sub-hypothesis.

Outputs:
    data/reports/imcoh_vi/h2_band_selectivity.md
    data/reports/imcoh_vi/h2_band_selectivity.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT

OUT_DIR = REPORTS_ROOT / "imcoh_vi"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT


# ─────────────────────────── loaders ───────────────────────────

def load_h2a_contrasts() -> pd.DataFrame:
    """k-averaged VI contrast per (patient, band) for H2a."""
    df = pd.read_csv(OUT_DIR / "hypothesis_contrasts.csv")
    h2a = df[df["hypothesis"] == "H2a"].copy()
    return (
        h2a.groupby(["patient", "band"])["contrast"].mean()
        .unstack("band").reindex(columns=BANDS)
    )


def load_h2c_rho(target: str = "task_test") -> pd.DataFrame:
    """Per-patient ρ per band for H2c. target ∈ {task_test, task_learn}."""
    df = pd.read_csv(OUT_DIR / "h2c_ultrametric_drift_raw.csv")
    col = "rho_task" if target == "task_test" else "rho_learn"
    have = "has_ttest" if target == "task_test" else "has_tlearn"
    df = df[df["has_rpre"] & df["has_rpost"] & df[have]].copy()
    wide = df.pivot(index="patient", columns="band", values=col).reindex(columns=BANDS)
    return wide.dropna()


def load_h2d_delta_rho() -> pd.DataFrame:
    """k-averaged Δρ per (patient, band) for H2d."""
    df = pd.read_csv(OUT_DIR / "h2d_persistence_raw.csv")
    agg = (
        df.groupby(["patient", "band"])["delta_rho"].mean()
        .unstack("band").reindex(columns=BANDS)
    )
    return agg.dropna()


# ─────────────────────────── stats ───────────────────────────

def friedman(table: pd.DataFrame) -> dict:
    """Friedman χ² on a (patients × bands) DataFrame."""
    data = table.to_numpy()
    n, k = data.shape
    chi2, p = stats.friedmanchisquare(*data.T)
    ranks = stats.rankdata(data, axis=1).mean(axis=0)
    kendall_w = chi2 / (n * (k - 1)) if n > 0 else np.nan
    return {
        "n": n,
        "chi2": chi2,
        "df": k - 1,
        "p": p,
        "kendall_W": kendall_w,
        "mean_rank": dict(zip(table.columns, ranks)),
    }


def pairwise_wilcoxon(table: pd.DataFrame, alternative: str = "two-sided") -> pd.DataFrame:
    """All pairs (i, j) with i < j. Returns DataFrame of (band_i, band_j, stat, p_raw, p_bonf)."""
    bands = list(table.columns)
    rows = []
    pairs = [(a, b) for i, a in enumerate(bands) for b in bands[i + 1:]]
    m = len(pairs)
    for a, b in pairs:
        x = table[a].to_numpy()
        y = table[b].to_numpy()
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]
        if len(x) < 3:
            rows.append((a, b, np.nan, np.nan, np.nan, len(x)))
            continue
        try:
            res = stats.wilcoxon(x, y, alternative=alternative, zero_method="wilcox")
            rows.append((a, b, float(res.statistic), float(res.pvalue),
                         min(1.0, float(res.pvalue) * m), len(x)))
        except ValueError:
            rows.append((a, b, np.nan, np.nan, np.nan, len(x)))
    return pd.DataFrame(rows, columns=["band_a", "band_b", "W", "p_raw", "p_bonf", "n"])


def theta_vs_others(table: pd.DataFrame, alt: str = "less") -> pd.DataFrame:
    """Focused θ-vs-others contrast, one-sided paired Wilcoxon (θ < other).

    alt='less' asks whether θ values are less than the comparison band values —
    correct for H2d Δρ (θ is the ergodic / lower band). Override for other
    sub-hypotheses if needed.
    """
    others = [b for b in table.columns if b != "theta"]
    m = len(others)
    rows = []
    for b in others:
        x = table["theta"].to_numpy()
        y = table[b].to_numpy()
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]
        if len(x) < 3:
            rows.append((b, np.nan, np.nan, np.nan, len(x)))
            continue
        res = stats.wilcoxon(x, y, alternative=alt, zero_method="wilcox")
        rows.append((b, float(res.statistic), float(res.pvalue),
                     min(1.0, float(res.pvalue) * m), len(x)))
    return pd.DataFrame(rows, columns=["other_band", "W", "p_raw", "p_bonf", "n"])


# ─────────────────────────── report ───────────────────────────

def band_rank_row(ranks: dict) -> str:
    return " | ".join(f"{TEX[b]}={ranks[b]:.2f}" for b in BANDS)


def fmt_friedman(label: str, fr: dict, note: str = "") -> str:
    sig = "★" if fr["p"] < 0.05 else " "
    out = [
        f"### {label}",
        "",
        f"n patients = {fr['n']}, k bands = {fr['df'] + 1}.",
        f"Friedman χ²({fr['df']}) = **{fr['chi2']:.3f}**, p = **{fr['p']:.4f}** {sig}",
        f"Kendall's W = {fr['kendall_W']:.3f}",
        "",
        "**Mean rank per band** (1 = lowest, 6 = highest):",
        "",
        "| " + " | ".join(TEX[b] for b in BANDS) + " |",
        "|" + "---|" * len(BANDS),
        "| " + " | ".join(f"{fr['mean_rank'][b]:.2f}" for b in BANDS) + " |",
        "",
    ]
    if note:
        out.append(note)
        out.append("")
    return "\n".join(out)


def fmt_pairs(df: pd.DataFrame, title: str = "") -> str:
    if df.empty:
        return ""
    lines = [title, "", "| band_a | band_b | W | p (raw) | p (Bonf) | n |",
             "|:-------|:-------|--:|--------:|---------:|--:|"]
    for _, r in df.iterrows():
        sig = "★" if (np.isfinite(r["p_bonf"]) and r["p_bonf"] < 0.05) else ""
        lines.append(
            f"| {TEX[r['band_a']]} | {TEX[r['band_b']]} | {r['W']:.1f} | "
            f"{r['p_raw']:.4f} | {r['p_bonf']:.4f}{sig} | {int(r['n'])} |"
        )
    return "\n".join(lines) + "\n"


def fmt_theta(df: pd.DataFrame, title: str) -> str:
    lines = [title, "",
             "One-sided paired Wilcoxon: H₀ θ ≥ other, H₁ θ < other. Bonferroni m=5.",
             "",
             "| other band | W | p (raw) | p (Bonf) | n |",
             "|:-----------|--:|--------:|---------:|--:|"]
    for _, r in df.iterrows():
        sig = "★" if (np.isfinite(r["p_bonf"]) and r["p_bonf"] < 0.05) else ""
        lines.append(
            f"| {TEX[r['other_band']]} | {r['W']:.1f} | {r['p_raw']:.4f} | "
            f"{r['p_bonf']:.4f}{sig} | {int(r['n'])} |"
        )
    return "\n".join(lines) + "\n"


# ─────────────────────────── main ───────────────────────────

def main():
    print("Loading per-patient×band tables...")
    tables = {
        "H2a (VI magnitude)": load_h2a_contrasts(),
        "H2c ρ (target=task_test)": load_h2c_rho("task_test"),
        "H2c ρ (target=task_learn)": load_h2c_rho("task_learn"),
        "H2d Δρ (block persistence)": load_h2d_delta_rho(),
    }

    # Theta direction: for H2a VI and H2d Δρ, smaller θ = more ergodic (alt="less").
    # For H2c ρ, "smaller ρ" = less direction-consistent; still alt="less".
    theta_alt = {
        "H2a (VI magnitude)": "less",
        "H2c ρ (target=task_test)": "less",
        "H2c ρ (target=task_learn)": "less",
        "H2d Δρ (block persistence)": "less",
    }

    md = ["# Band-selectivity: Friedman + post-hoc on H2 sub-hypotheses",
          "",
          "For each H2 sub-hypothesis we test whether the 6 bands behave the same.",
          "Rejection of Friedman H₀ = bands DIFFER = band-selectivity claim supported.",
          "Post-hoc: paired Wilcoxon θ vs each other band (one-sided, H₁: θ < other).",
          "Interpretation: θ significantly < other bands = θ is the ergodic exception.",
          ""]

    csv_rows = []

    for label, table in tables.items():
        print(f"\n{label}: n={len(table)} patients × {len(table.columns)} bands")
        fr = friedman(table)
        pairs = pairwise_wilcoxon(table, alternative="two-sided")
        theta = theta_vs_others(table, alt=theta_alt[label])

        md.append(f"## {label}")
        md.append("")
        md.append(fmt_friedman("Friedman omnibus", fr))
        md.append(fmt_theta(theta, "**Focused: θ vs each other band**"))
        md.append(fmt_pairs(pairs, "**All pairwise (two-sided)**"))
        md.append("")

        # Per-patient table
        md.append("**Per-patient values**")
        md.append("")
        md.append("| patient | " + " | ".join(TEX[b] for b in BANDS) + " |")
        md.append("|:--------|" + "-----:|" * len(BANDS))
        for pat, row in table.iterrows():
            md.append(f"| {pat} | " + " | ".join(
                f"{row[b]:+.3f}" if np.isfinite(row[b]) else "—" for b in BANDS
            ) + " |")
        md.append("")

        # CSV rows
        csv_rows.append({
            "measure": label,
            "n": fr["n"],
            "friedman_chi2": fr["chi2"],
            "friedman_p": fr["p"],
            "kendall_W": fr["kendall_W"],
            **{f"mean_rank_{b}": fr["mean_rank"][b] for b in BANDS},
            **{f"theta_vs_{r['other_band']}_p_bonf": r["p_bonf"]
               for _, r in theta.iterrows()},
        })

    # Write outputs
    md_path = OUT_DIR / "h2_band_selectivity.md"
    md_path.write_text("\n".join(md))
    print(f"\nWrote {md_path}")

    csv_path = OUT_DIR / "h2_band_selectivity.csv"
    pd.DataFrame(csv_rows).to_csv(csv_path, index=False)
    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()
