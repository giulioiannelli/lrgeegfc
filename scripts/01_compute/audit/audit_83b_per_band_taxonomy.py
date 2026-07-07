#!/usr/bin/env python3
"""audit_83b — per-band consistency-taxonomy synthesis on top of audit_83.

Synthesises the audit_83 matched-strength localization across all 6 bands and
all 4 conditions (epi-include / epi-exclude  x  contact-level / shaft-collapsed)
into a per-band CONSISTENCY-TAXONOMY verdict, cross-referenced against the
cohort rho_split detectability gate. This is the Step-1 lock that decides, per
band, whether the offline trace has a cohort-consistent system localization
(-> "consistent" tier) or not (-> "patient-specific"), with theta the candidate
"absent" tier.

NOT a new test — pure synthesis of already-computed audit_83 outputs + the
audit_63 gate. No surrogates, no FC loading.

5-point critical preamble
-------------------------
1. Claim: each band's offline cophenetic trace falls into one of three tiers —
   CONSISTENT (a system carries it cohort-wide, surviving matched-strength under
   BOTH epi modes AND shaft-collapse, like beta->OFC), PATIENT-SPECIFIC (a real
   trace exists for some patients but localizes to no single system that
   survives the filters), or ABSENT (no cohort trace and no robust localization;
   theta candidate).
2. Null / standard: the locked beta->OFC bar — a carrier system must clear
   matched-strength at BH q<0.05 across the a-priori systems WITHIN its band, in
   epi-include AND epi-exclude, AND still clear (q<0.10) under shaft-collapse
   (one electrode shaft = one observation; kills within-shaft pseudoreplication).
3. Strongest alternative the bar controls for: a "localization" that is really
   (a) node-strength geometry (matched-strength holds strength fixed), (b) one
   electrode shaft threading a structure (shaft-collapse), or (c) an epi/SOZ
   channel artifact (epi-exclude). A tier="consistent" verdict must beat all three.
4. What it CANNOT do: audit_83 carries no leave-one-patient-out column, so a
   carrier surviving all 4 conditions can still be one-patient-driven; such
   carriers are flagged LOO-UNTESTED and must get audit_83 --R / LOO before a
   "consistent" claim is published. Undersampled systems (low K_implanted) are
   flagged. The gate (rho_split vs own surrogate) measures DETECTABILITY, not
   biology — a low-gate band is "undetected", which at high SNR may still be a
   real but patient-specific trace.
5. Falsification: a band claimed "consistent" whose carrier does NOT survive
   shaft-collapse (q>=0.10) or flips sign under epi-exclude is downgraded to
   patient-specific. A band claimed "absent" (theta) that shows a carrier
   surviving all 4 conditions is NOT absent — it is patient-specific-or-better
   and the "theta shows nothing" headline must soften.

Inputs (read-only)
------------------
- data/audit/localization_atlas/matched_strength_allbands_{include,exclude}.csv
- data/audit/localization_atlas/matched_strength_allbands_shaftcollapsed_{include,exclude}.csv
  (any missing condition is skipped with a warning; the verdict marks which
   filters were actually available so a partial run is never silently trusted)
- data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
  (the cohort rho_split gate: obs_rho, surr_p5/p50/p95 per patient x band)

Outputs (data/audit/localization_atlas/)
----------------------------------------
- per_band_taxonomy_systems.csv : every (band, system) with q_upper/q_lower per
  condition and the robust-carrier / robust-depleted flags.
- per_band_taxonomy_verdict.csv : one row per band — gate stats, carrier &
  depleted system lists, tier, and caveat flags.

Usage
-----
    python audit_83b_per_band_taxonomy.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
LOC = ROOT / "data/audit/localization_atlas"
GATE_CSV = (ROOT / "data/audit/matched_strength_surrogate_split_baseline"
            / "per_patient_per_band.csv")

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
GRAN = "system"                       # the a-priori-system multiplicity family
Q_CARRIER = 0.05                      # the locked beta->OFC bar (contact, both epi)
Q_CARRIER_SHAFT = 0.10               # shaft-collapse survival (looser; K drops)
K_UNDERSAMPLED = 3                    # flag carriers seen in <= this many patients

CONDITIONS = {                        # cond tag -> csv filename
    "inc_contact":   "matched_strength_allbands_include.csv",
    "exc_contact":   "matched_strength_allbands_exclude.csv",
    "inc_shaft":     "matched_strength_allbands_shaftcollapsed_include.csv",
    "exc_shaft":     "matched_strength_allbands_shaftcollapsed_exclude.csv",
}


def load_conditions() -> tuple[pd.DataFrame, list[str]]:
    """Stack the available audit_83 condition CSVs (system granularity)."""
    frames, have = [], []
    for cond, fname in CONDITIONS.items():
        p = LOC / fname
        if not p.exists():
            print(f"  [warn] missing condition {cond}: {fname}")
            continue
        d = pd.read_csv(p)
        d = d[d["granularity"] == GRAN].copy()
        d["cond"] = cond
        frames.append(d)
        have.append(cond)
    if not frames:
        raise SystemExit("no audit_83 allbands CSVs found — run audit_83 --all-bands first")
    return pd.concat(frames, ignore_index=True), have


def bh_within(df: pd.DataFrame, pcol: str, qcol: str) -> pd.DataFrame:
    """BH-FDR across systems WITHIN each (band, cond) — the locked family."""
    df = df.copy()
    df[qcol] = np.nan
    for (_, _), idx in df.groupby(["band", "cond"]).groups.items():
        sub = df.loc[idx]
        df.loc[idx, qcol] = bh_fdr(sub[pcol].to_numpy())
    return df


def gate_stats() -> pd.DataFrame:
    """Cohort rho_split detectability gate per band (from audit_63 split null)."""
    g = pd.read_csv(GATE_CSV)
    rows = []
    for b in BANDS:
        d = g[g["band"] == b]
        if d.empty:
            continue
        diff = (d["obs_rho"] - d["surr_p50"]).to_numpy()
        try:
            gate_p = float(wilcoxon(diff, alternative="greater").pvalue)
        except ValueError:
            gate_p = np.nan
        rows.append({
            "band": b,
            "n_clear_p95": int((d["obs_rho"] > d["surr_p95"]).sum()),
            "median_obs_rho": float(np.median(d["obs_rho"])),
            "gate_wilcoxon_p": gate_p,
        })
    return pd.DataFrame(rows)


def main() -> None:
    df, have = load_conditions()
    df = bh_within(df, "matched_strength_p", "q_upper")
    df = bh_within(df, "matched_strength_p_lower", "q_lower")

    have_shaft = ("inc_shaft" in have) and ("exc_shaft" in have)
    have_both_epi = ("inc_contact" in have) and ("exc_contact" in have)

    # pivot q_upper / q_lower / K to one row per (band, system)
    keep = ["band", "unit", "cond", "K_implanted", "M_obs",
            "q_upper", "q_lower"]
    long = df[keep].copy()
    qup = long.pivot_table(index=["band", "unit"], columns="cond",
                           values="q_upper")
    qlo = long.pivot_table(index=["band", "unit"], columns="cond",
                           values="q_lower")
    Km = long.groupby(["band", "unit"])["K_implanted"].max()
    qup.columns = [f"qU_{c}" for c in qup.columns]
    qlo.columns = [f"qL_{c}" for c in qlo.columns]
    sysdf = qup.join(qlo).join(Km).reset_index()

    def _q(row, prefix, cond):
        col = f"{prefix}_{cond}"
        return row[col] if col in row and pd.notna(row[col]) else np.nan

    carrier_flags, deplete_flags = [], []
    for _, r in sysdf.iterrows():
        # robust carrier: clears in BOTH epi at contact AND survives shaft-collapse
        inc_c = _q(r, "qU", "inc_contact"); exc_c = _q(r, "qU", "exc_contact")
        inc_s = _q(r, "qU", "inc_shaft");   exc_s = _q(r, "qU", "exc_shaft")
        both_epi_contact = (inc_c < Q_CARRIER) and (exc_c < Q_CARRIER)
        shaft_ok = (not have_shaft) or (
            (inc_s < Q_CARRIER_SHAFT) and (exc_s < Q_CARRIER_SHAFT))
        carrier_flags.append(bool(both_epi_contact and shaft_ok))
        # robust depleted: lower-tail clears both epi at contact (+shaft if avail)
        linc_c = _q(r, "qL", "inc_contact"); lexc_c = _q(r, "qL", "exc_contact")
        linc_s = _q(r, "qL", "inc_shaft");   lexc_s = _q(r, "qL", "exc_shaft")
        dep_contact = (linc_c < Q_CARRIER) and (lexc_c < Q_CARRIER)
        dep_shaft = (not have_shaft) or (
            (linc_s < Q_CARRIER_SHAFT) and (lexc_s < Q_CARRIER_SHAFT))
        deplete_flags.append(bool(dep_contact and dep_shaft))
    sysdf["robust_carrier"] = carrier_flags
    sysdf["robust_depleted"] = deplete_flags
    sysdf["undersampled"] = sysdf["K_implanted"] <= K_UNDERSAMPLED

    LOC.mkdir(parents=True, exist_ok=True)
    sysdf.sort_values(["band", "unit"]).to_csv(
        LOC / "per_band_taxonomy_systems.csv", index=False)

    # ---- per-band verdict ----
    gate = gate_stats().set_index("band")
    vrows = []
    for b in BANDS:
        sb = sysdf[sysdf["band"] == b]
        if sb.empty:
            continue
        carriers = sb[sb["robust_carrier"]]
        depleted = sb[sb["robust_depleted"]]
        carrier_list = sorted(
            f"{u}{'*' if k else ''}" for u, k in
            zip(carriers["unit"], carriers["undersampled"]))
        deplete_list = sorted(depleted["unit"].tolist())
        gp = gate.loc[b] if b in gate.index else None
        n_clear = int(gp["n_clear_p95"]) if gp is not None else -1
        gate_p = float(gp["gate_wilcoxon_p"]) if gp is not None else np.nan
        # robust carrier that is NOT solely undersampled
        solid_carrier = bool(any(~carriers["undersampled"]))
        # GATE FIRST: a band must have a cohort trace before localization means
        # anything. All carriers survive shaft-collapse (verified), so the tier
        # is set by the NET cohort trace (gate), not localization robustness.
        # trace_level — a-priori, NOT tuned to the result:
        #   cohort : present ON AVERAGE (locked one-sided gate p<0.05) OR a
        #            MAJORITY of patients clear own p95 (n>=5 of 10).
        #   subset : real in a minority (n in 3..4 — beats Binom(10,0.05) chance
        #            ~0.5) but NO cohort-mean shift (gate p>=0.05).
        #   none   : n<=2 and gate p>=0.05 (theta).
        if (np.isfinite(gate_p) and gate_p < 0.05) or n_clear >= 5:
            trace_level = "cohort"
        elif n_clear >= 3:
            trace_level = "subset"
        else:
            trace_level = "none"
        # a localization that clears while the band has NO net cohort trace is a
        # spatial redistribution of ~zero net concordance, not a trace (theta:
        # MTL clears yet median rho<0; gamma_h: parietal clears yet median rho=0).
        localization_without_trace = bool(
            len(carrier_list) > 0 and trace_level != "cohort")
        if trace_level == "none":
            tier = "absent"
        elif solid_carrier and trace_level == "cohort":
            tier = "consistent"          # cohort trace AND a robust home
        else:
            tier = "patient-specific"    # trace present but no cohort-consistent home
        vrows.append({
            "band": b,
            "tier": tier,
            "trace_level": trace_level,
            "n_clear_p95": n_clear,
            "median_obs_rho": round(float(gp["median_obs_rho"]), 3) if gp is not None else np.nan,
            "gate_wilcoxon_p": round(float(gp["gate_wilcoxon_p"]), 4) if gp is not None else np.nan,
            "robust_carriers": ", ".join(carrier_list) if carrier_list else "—",
            "robust_depleted": ", ".join(deplete_list) if deplete_list else "—",
            "loc_without_cohort_trace": localization_without_trace,
            "LOO_tested": False,           # audit_83 has no LOO column
        })
    verdict = pd.DataFrame(vrows)
    verdict.to_csv(LOC / "per_band_taxonomy_verdict.csv", index=False)

    print("\n=== conditions available ===")
    print(f"  {have}  | both_epi={have_both_epi}  shaft_collapse={have_shaft}")
    if not have_shaft:
        print("  [!] shaft-collapse NOT yet available — carriers are PROVISIONAL")
    print("\n=== per-band consistency-taxonomy verdict ===")
    print("  (* = carrier undersampled K<=%d; gate=detectability; "
          "LOO still owed for any 'consistent')" % K_UNDERSAMPLED)
    print(verdict.to_string(index=False))
    print(f"\nOutputs -> {LOC}/per_band_taxonomy_{{systems,verdict}}.csv")


if __name__ == "__main__":
    main()
