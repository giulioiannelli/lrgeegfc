#!/usr/bin/env python3
"""diag — unpack the audit_83 cohort-median localization into its variability.

The matched-strength localization test (audit_83) reports a single cohort-median
of per-patient demeaned unit-means and a surrogate p. That summary hides exactly
the variability a referee will ask about:

  (1) BETWEEN-PATIENT: how many of the K implanted patients are actually positive
      in the unit? Is the median carried by 5/5 or by 3/5 with one big driver?
  (2) WITHIN-PATIENT: inside one patient's unit, what fraction of pair-incidences
      are trace (c_ij>0) vs anti-trace (c_ij<0)? Heavy cancellation = a near-zero
      net even when "trace pairs exist".
  (3) LEAVE-ONE-OUT on the OBSERVED cohort median: does dropping any single
      patient flip the sign of the median (a one-patient artifact)?

This is a DESCRIPTIVE decomposition of the already-scoped concordance measure
(no new measure, no surrogate): it reuses the exact concordance + global-demean +
per-unit-mean pipeline of audit_83 on the observed split-baseline trace, and just
reports the per-patient terms instead of their cohort median.

c_ij = (rank(dD_task)-mid)(rank(dD_rest)-mid)  -> per-pair rho_split contribution.
c_ij>0 = sign-consistent trace direction; <0 = anti-trace.

Outputs: prints tables; writes data/audit/localization_atlas/per_patient_decomp_{epi}.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.utils.surrogate import cophenetic_condensed_from_eigs
from lrg_eegfc.workflow.fc import load_fc_matrix

ROOT = setup_script_env()
OUT = ROOT / "data/audit/localization_atlas"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
PARACORE = {"MTL", "cingulate", "insula"}
BANDS = ["beta"]
# (granularity-column, unit) targets to unpack
TARGETS = [
    ("supersystem", "limbic"),
    ("system", "OFC"),
    ("paracore", "paralimbic_core"),
    ("system", "lateral_temporal"),
    ("system", "cingulate"),
    ("system", "insula"),
    ("system", "MTL"),
    ("region", "Hip"),
]


def concordance(dt, dr):
    n = dt.size
    mid = (n + 1) / 2.0
    return (rankdata(dt) - mid) * (rankdata(dr) - mid)


def _canon_cophenet(W):
    W = np.asarray(W, dtype=np.float64)
    ev, V = np.linalg.eigh(np.diag(W.sum(axis=1)) - W)
    return cophenetic_condensed_from_eigs(ev, V)


def _canon_trace(pat, band):
    """CANONICAL observed split-baseline trace (matches audit_83.obs_trace /
    the surrogate construction) — not the LRG-pipeline per_pair_split."""
    Wt = load_fc_matrix(pat, "task_test", band, "imcoh_abs")
    Wpost = load_fc_matrix(pat, "rest_post", band, "imcoh_abs")
    WA = np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_A_imcoh_abs.npy")
    WB = np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_B_imcoh_abs.npy")
    Dt, Dpost, DA, DB = (_canon_cophenet(W) for W in (Wt, Wpost, WA, WB))
    N = Wt.shape[0]
    iu = np.triu_indices(N, k=1)
    return (concordance(Dt - DA, Dpost - DB),
            iu[0].astype(int), iu[1].astype(int))


def patient_terms(pat, band, epi_excl):
    """Return dict (gran,unit) -> (mean_demeaned, n_incid, frac_raw_pos) for pat."""
    c, iu_i, iu_j = _canon_trace(pat, band)
    rdf = load_channel_regions(pat)
    rdf["paracore"] = rdf["system"].map(
        lambda s: "paralimbic_core" if s in PARACORE else "x")
    keep = epi_keep_mask(rdf, pat) if epi_excl else None

    vals = np.concatenate([c, c])
    nidx = np.concatenate([iu_i, iu_j])
    if keep is not None:
        m = keep[nidx]
        vals, nidx = vals[m], nidx[m]
    raw = vals.copy()
    dem = vals - vals.mean()

    res = {}
    for col, unit in TARGETS:
        uv = rdf[col].to_numpy().astype(str)[nidx]
        sel = uv == unit
        if not sel.any():
            continue
        res[(col, unit)] = (float(dem[sel].mean()),
                            int(sel.sum()),
                            float((raw[sel] > 0).mean()))
    return res


def run(band, epi_excl):
    etag = "exclude" if epi_excl else "include"
    rows = []
    bypat = {}
    for pat in COHORT:
        try:
            bypat[pat] = patient_terms(pat, band, epi_excl)
        except (FileNotFoundError, Exception) as e:           # noqa: BLE001
            print(f"  [skip] {pat}: {e}")

    print(f"\n{'='*78}\n  band={band}  epi-{etag}\n{'='*78}")
    hdr = (f"  {'unit':<22s} {'K':>2s} {'pos/K':>6s} {'median':>10s} "
           f"{'LOO-min':>10s} {'wMean+frac':>10s}")
    print(hdr)
    print("  " + "-" * 74)
    for col, unit in TARGETS:
        terms = [(p, *bypat[p][(col, unit)]) for p in bypat
                 if (col, unit) in bypat[p]]
        if not terms:
            continue
        means = np.array([t[1] for t in terms])
        fracs = np.array([t[3] for t in terms])
        K = len(means)
        npos = int((means > 0).sum())
        med = float(np.median(means))
        # leave-one-out on the cohort median (observed)
        loo = [float(np.median(np.delete(means, i))) for i in range(K)]
        loo_min = min(loo)
        wfrac = float(np.mean(fracs))      # typical within-unit trace fraction
        flag = ""
        if med > 0 and loo_min <= 0:
            flag = "  <- median flips sign under LOO"
        print(f"  {unit:<22s} {K:>2d} {npos:>3d}/{K:<2d} {med:>10.0f} "
              f"{loo_min:>10.0f} {wfrac:>9.2f}{flag}")
        # per-patient detail
        detail = "    " + "  ".join(
            f"{p.split('_')[1]}:{m:+.0f}({f:.0%},n{n})"
            for p, m, n, f in sorted(terms, key=lambda t: -t[1]))
        print(detail)
        for p, m, n, f in terms:
            rows.append({"band": band, "epi": etag, "granularity": col,
                         "unit": unit, "patient": p, "mean_demeaned": m,
                         "n_incid": n, "frac_trace_raw": f})
    return rows


def main():
    allrows = []
    for epi_excl in (False, True):
        for band in BANDS:
            allrows += run(band, epi_excl)
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(allrows)
    df.to_csv(OUT / "per_patient_decomp.csv", index=False)
    print(f"\nWrote {OUT}/per_patient_decomp.csv")
    print("\nLegend: pos/K = patients with positive demeaned unit-mean; "
          "median = cohort median (the audit_83 statistic, pre-surrogate); "
          "LOO-min = smallest median dropping one patient; "
          "wMean+frac = mean within-unit fraction of TRACE (c>0) incidences "
          "(0.50 = full cancellation).")


if __name__ == "__main__":
    main()
