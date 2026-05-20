"""Regate audit_66 / audit_67 Grassmann cohort verdicts on Wilcoxon p<0.05
alone, dropping the hardcoded ``n_below >= 8`` patient-count filter and the
hardcoded ``|med_surr| < 0.05·|med_obs|`` magnitude-ratio filter that were
holdovers from §5.4 visualization conventions.

No new computation. Reads existing cohort_summary.csv from both audits and
recomputes:
    sig_cell(band, k) = (paired_wilcoxon_p < 0.05)   # cohort one-sided
plus longest contiguous run per band, side-by-side with the old gate.

Verifies the three manuscript windows that load-bear §7 errata:
    β       k = 27..55      (audit_66 cited 29 contiguous)
    γ_l     k = 12..23      (audit_66 cited 12 contiguous)
    γ_h     k = 19..27      (audit_66 cited 9 contiguous)
and the audit_67 persist counts within each window.

Outputs:
    data/audit/grassmann_regate_no_filter/per_cell_audit66.csv
    data/audit/grassmann_regate_no_filter/per_cell_audit67.csv
    data/audit/grassmann_regate_no_filter/contig_summary.csv
    .agents/reports/2026-05-15_grassmann-regating-no-7of10-filter.md
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
A66 = ROOT / "data" / "audit" / "grassmann_matched_strength_surrogate"
A67 = ROOT / "data" / "audit" / "grassmann_epi_exclusion"
OUT = ROOT / "data" / "audit" / "grassmann_regate_no_filter"
REPORT = ROOT / ".agents" / "reports" / "2026-05-15_grassmann-regating-no-7of10-filter.md"

BANDS_ORDER = ["delta", "theta", "alpha", "low_gamma", "beta", "high_gamma"]
WINDOWS = {
    "beta":       (27, 55),
    "low_gamma":  (12, 23),
    "high_gamma": (19, 27),
}


def _band_present(df: pd.DataFrame) -> list[str]:
    present = set(df.band.unique().tolist())
    return [b for b in BANDS_ORDER if b in present]


def regate_one(cohort: pd.DataFrame, label: str) -> pd.DataFrame:
    """Return per-cell verdict table with old + new gate side by side."""
    rows = []
    for _, r in cohort.iterrows():
        p = float(r["paired_wilcoxon_p"])
        med_obs = float(r["obs_median_T_G"])
        med_surr = float(r["surr_median_T_G_per_patient_median"])
        n_below = int(r["n_patients_below_own_surrogate"])
        # Old gate (3 conditions, all hardcoded)
        old_sep = (
            (p < 0.05)
            and (abs(med_surr) < 0.05 * max(1.0, abs(med_obs)))
            and (n_below >= 8)
        )
        # New gate (Wilcoxon only)
        new_sig = (p < 0.05)
        ratio = (med_obs / med_surr) if med_surr != 0 else float("nan")
        rows.append({
            "audit": label,
            "band": r["band"],
            "k": int(r["k"]),
            "T_G_obs_cohort_median": med_obs,
            "T_G_surr_cohort_median": med_surr,
            "ratio_obs_over_surr": ratio,
            "wilcoxon_p_one_sided_less": p,
            "n_trace_descriptive": n_below,
            "old_gate_separated": bool(old_sep),
            "new_gate_sig": bool(new_sig),
        })
    return pd.DataFrame(rows)


def contig_runs(mask: np.ndarray) -> tuple[int, int, int]:
    """Longest contiguous True run; returns (length, start_idx, end_idx)."""
    if not mask.any():
        return 0, -1, -1
    longest = 0; longest_start = -1; longest_end = -1
    cur = 0; cur_start = -1
    for i, v in enumerate(mask):
        if v:
            if cur == 0:
                cur_start = i
            cur += 1
            if cur > longest:
                longest = cur
                longest_start = cur_start
                longest_end = i
        else:
            cur = 0
    return longest, longest_start, longest_end


def contig_summary(per_cell: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = []
    for band in _band_present(per_cell):
        sub = per_cell[per_cell.band == band].sort_values("k").reset_index(drop=True)
        ks = sub.k.values
        old = sub.old_gate_separated.values.astype(bool)
        new = sub.new_gate_sig.values.astype(bool)
        n_old = int(old.sum())
        n_new = int(new.sum())
        L_old, s_old, e_old = contig_runs(old)
        L_new, s_new, e_new = contig_runs(new)
        # best k by ratio magnitude (only over new-gate-sig cells if any, else all)
        scope = sub[sub.new_gate_sig] if sub.new_gate_sig.any() else sub
        if not scope.empty:
            r = scope.assign(absratio=scope.ratio_obs_over_surr.abs()).sort_values(
                "absratio", ascending=False).iloc[0]
            best_k = int(r["k"])
            best_ratio = float(r["ratio_obs_over_surr"])
        else:
            best_k = -1; best_ratio = float("nan")
        rows.append({
            "audit": label,
            "band": band,
            "n_total_cells": int(len(ks)),
            "n_sig_old_gate": n_old,
            "n_sig_new_gate": n_new,
            "delta_n_sig": n_new - n_old,
            "longest_run_old": L_old,
            "longest_run_old_k_start": int(ks[s_old]) if s_old >= 0 else -1,
            "longest_run_old_k_end":   int(ks[e_old]) if e_old >= 0 else -1,
            "longest_run_new": L_new,
            "longest_run_new_k_start": int(ks[s_new]) if s_new >= 0 else -1,
            "longest_run_new_k_end":   int(ks[e_new]) if e_new >= 0 else -1,
            "best_k_ratio_magnitude": best_k,
            "best_k_ratio_value": best_ratio,
        })
    return pd.DataFrame(rows)


def window_check(per_cell: pd.DataFrame, band: str, k0: int, k1: int,
                 label: str) -> dict:
    sub = per_cell[(per_cell.band == band) & per_cell.k.between(k0, k1)].sort_values("k")
    n_window = len(sub)
    n_old = int(sub.old_gate_separated.sum())
    n_new = int(sub.new_gate_sig.sum())
    # Longest contiguous run within the window itself
    new_mask = sub.new_gate_sig.values.astype(bool)
    L_new, s_new, e_new = contig_runs(new_mask)
    return {
        "audit": label,
        "band": band,
        "window_k_start": k0,
        "window_k_end": k1,
        "window_n_cells": n_window,
        "n_sig_old_gate_in_window": n_old,
        "n_sig_new_gate_in_window": n_new,
        "longest_run_in_window_new": L_new,
        "longest_run_in_window_new_k_start":
            int(sub.k.values[s_new]) if s_new >= 0 else -1,
        "longest_run_in_window_new_k_end":
            int(sub.k.values[e_new]) if e_new >= 0 else -1,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)

    c66 = pd.read_csv(A66 / "cohort_summary.csv")
    c67 = pd.read_csv(A67 / "cohort_summary.csv")

    pc66 = regate_one(c66, "audit_66_full")
    pc67 = regate_one(c67, "audit_67_epiX")
    pc66.to_csv(OUT / "per_cell_audit66.csv", index=False)
    pc67.to_csv(OUT / "per_cell_audit67.csv", index=False)

    cs66 = contig_summary(pc66, "audit_66_full")
    cs67 = contig_summary(pc67, "audit_67_epiX")
    pd.concat([cs66, cs67], ignore_index=True).to_csv(
        OUT / "contig_summary.csv", index=False)

    win_rows = []
    for band, (k0, k1) in WINDOWS.items():
        win_rows.append(window_check(pc66, band, k0, k1, "audit_66_full"))
        win_rows.append(window_check(pc67, band, k0, k1, "audit_67_epiX"))
    win_df = pd.DataFrame(win_rows)
    win_df.to_csv(OUT / "window_check.csv", index=False)

    # ------------------------- Markdown report -----------------------------
    lines: list[str] = []
    lines.append("---")
    lines.append("name: 2026-05-15_grassmann-regating-no-7of10-filter")
    lines.append("era: IMCOH_ABS_COHORT_N10")
    lines.append("status: current")
    lines.append("kind: report")
    lines.append("audit: 69")
    lines.append("upstream: [audit_66, audit_67]")
    lines.append("---")
    lines.append("")
    lines.append("# Grassmann regating — drop hardcoded patient-count + magnitude filters")
    lines.append("")
    lines.append(
        "The §5.4 audit_66 / audit_67 cohort verdict tagged a (band, k) cell "
        "**separated** only when **three** criteria hit at once: cohort-paired "
        "one-sided Wilcoxon `p < 0.05` (the statistical gate), `|cohort median "
        "surr T_G| < 0.05 × max(1, |cohort median obs T_G|)` (a hardcoded "
        "magnitude-ratio filter), and `n_below ≥ 8/10` patients individually "
        "below their own surrogate at one-sided `p < 0.05` (a hardcoded "
        "patient-count direction-agreement filter). The two extra filters were "
        "holdovers from the §5.4 small-multiples visualization convention and "
        "were never disclosed as gating criteria in the §7 errata or in the "
        "§5.4 manuscript text. This regating run drops both and gates every "
        "(band, k) cell on the cohort Wilcoxon `p < 0.05` alone."
    )
    lines.append("")
    lines.append("## Per-band side-by-side: old gate vs Wilcoxon-only gate")
    lines.append("")
    lines.append(
        "| audit | band | n_k | n_sig old | n_sig new | Δ | "
        "old run (k_start–k_end, len) | new run (k_start–k_end, len) | "
        "best-k ratio |"
    )
    lines.append(
        "|---|---|---|---|---|---|---|---|---|"
    )
    for _, r in pd.concat([cs66, cs67], ignore_index=True).iterrows():
        flag = ""
        if abs(int(r["delta_n_sig"])) > 3:
            flag = " **⚠**"
        old_run_str = (
            f"k={r['longest_run_old_k_start']}–{r['longest_run_old_k_end']}, "
            f"L={r['longest_run_old']}" if r["longest_run_old"] > 0 else "—"
        )
        new_run_str = (
            f"k={r['longest_run_new_k_start']}–{r['longest_run_new_k_end']}, "
            f"L={r['longest_run_new']}" if r["longest_run_new"] > 0 else "—"
        )
        lines.append(
            f"| {r['audit']} | {r['band']} | {r['n_total_cells']} "
            f"| {r['n_sig_old_gate']} | {r['n_sig_new_gate']} "
            f"| {int(r['delta_n_sig']):+d}{flag} | {old_run_str} | {new_run_str} "
            f"| k={r['best_k_ratio_magnitude']} ({r['best_k_ratio_value']:+.2f}) |"
        )
    lines.append("")
    lines.append(
        "**Δ flag (⚠)** marks cells where the new gate changes the count by "
        "more than ±3 cells. Sign of Δ is `new − old` — positive means the "
        "Wilcoxon-only gate accepts cells the old composite gate rejected."
    )
    lines.append("")
    lines.append("## Manuscript-window load-bearing claims (§7 errata)")
    lines.append("")
    lines.append(
        "The §7 errata cites three Grassmann manuscript windows under the "
        "old composite gate. Recomputed under the Wilcoxon-only gate inside "
        "the same window:"
    )
    lines.append("")
    lines.append(
        "| audit | band | window k=k0..k1 | window cells | n_sig old in window "
        "| n_sig new in window | longest contig run (new) within window |"
    )
    lines.append("|---|---|---|---|---|---|---|")
    for _, r in win_df.iterrows():
        run_str = (
            f"L={r['longest_run_in_window_new']} "
            f"(k={r['longest_run_in_window_new_k_start']}–"
            f"{r['longest_run_in_window_new_k_end']})"
            if r["longest_run_in_window_new"] > 0 else "—"
        )
        lines.append(
            f"| {r['audit']} | {r['band']} "
            f"| {r['window_k_start']}..{r['window_k_end']} "
            f"| {r['window_n_cells']} "
            f"| {r['n_sig_old_gate_in_window']} "
            f"| {r['n_sig_new_gate_in_window']} "
            f"| {run_str} |"
        )
    lines.append("")
    lines.append("## Verdict per window (new vs old)")
    lines.append("")
    for band, (k0, k1) in WINDOWS.items():
        w66 = win_df[(win_df.audit == "audit_66_full") & (win_df.band == band)].iloc[0]
        w67 = win_df[(win_df.audit == "audit_67_epiX") & (win_df.band == band)].iloc[0]
        d66 = int(w66.n_sig_new_gate_in_window) - int(w66.n_sig_old_gate_in_window)
        d67 = int(w67.n_sig_new_gate_in_window) - int(w67.n_sig_old_gate_in_window)
        verdict_66 = (
            "unchanged" if d66 == 0 else
            "expanded" if d66 > 0 else "shrunk"
        )
        verdict_67 = (
            "unchanged" if d67 == 0 else
            "expanded" if d67 > 0 else "shrunk"
        )
        lines.append(
            f"- **{band} k={k0}..{k1}** — audit_66 window count "
            f"{w66.n_sig_old_gate_in_window} → {w66.n_sig_new_gate_in_window} "
            f"({verdict_66}, Δ={d66:+d}); "
            f"audit_67 window count {w67.n_sig_old_gate_in_window} → "
            f"{w67.n_sig_new_gate_in_window} ({verdict_67}, Δ={d67:+d})."
        )
    lines.append("")
    lines.append("## What the new gate means")
    lines.append("")
    lines.append(
        "Under the Wilcoxon-only gate, **sig_cell(band, k)** is the cohort-"
        "paired one-sided Wilcoxon signed-rank test of per-patient `T_G_obs` "
        "vs per-patient `T_G_surr_median` rejecting `H_0: T_G_obs ≥ T_G_surr` "
        "at α=0.05. `n_below_own_surrogate` is reported as a descriptive "
        "column but no longer enters the gate. The magnitude-ratio criterion "
        "is also dropped — Wilcoxon ranks already encode the magnitude "
        "evidence the ratio filter was duplicating."
    )
    lines.append("")
    lines.append("## Provenance")
    lines.append("")
    lines.append("- Inputs:")
    lines.append(f"  - `{A66.relative_to(ROOT)}/cohort_summary.csv`")
    lines.append(f"  - `{A67.relative_to(ROOT)}/cohort_summary.csv`")
    lines.append("- Outputs:")
    lines.append(f"  - `{(OUT / 'per_cell_audit66.csv').relative_to(ROOT)}`")
    lines.append(f"  - `{(OUT / 'per_cell_audit67.csv').relative_to(ROOT)}`")
    lines.append(f"  - `{(OUT / 'contig_summary.csv').relative_to(ROOT)}`")
    lines.append(f"  - `{(OUT / 'window_check.csv').relative_to(ROOT)}`")
    lines.append("- Build script:")
    lines.append("  - `scripts/01_compute/audit/audit_69_grassmann_regate_no_filter.py`")
    lines.append("")

    REPORT.write_text("\n".join(lines))
    print(f"[audit_69] wrote {REPORT}")
    print(f"[audit_69] per-cell CSVs and contig_summary in {OUT}")


if __name__ == "__main__":
    main()
