#!/usr/bin/env python3
"""S1 -- is the incumbent trace diluted, and by what?

The question this stage answers, and the only one: at each diffusion scale, how
much of ``rho_sym`` is contributed by contact pairs the tree actually resolves at
that scale, versus pairs already merged below it or not yet separated above it?
If the contribution is flat across tree levels, the dilution hypothesis is wrong,
and the report must say so and stop pursuing scale-local band-passes -- that
condition is pre-registered in the scope report, §5.1, and is checked here
verbatim.

Five diagnostics, all observed-only (they describe the readout; they are not
tests and nothing here is compared to zero):

  A. cross-scale similarity of the readout's INPUT -- how much the cophenetic
     vector and the task reorganisation actually change along the scale axis,
     set against how much the hierarchy changes (N_eff runs 118 -> 4)
  B. how many distinct values a Spearman over ~7000 pairs really sees, given
     that a UPGMA tree offers at most N-1 cophenetic values
  C. how the pairs are distributed over tree levels
  D. the contribution decomposition: share of the trace by tree level, per scale
  E. the pre-registered §5.1 falsification check of the dilution hypothesis

Cross-scale correlations in (E) are computed on per-patient MARGINS against the
matched-strength null, the same object W0-C measured n_eff on, so the numbers are
comparable to its 1.2-1.9.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.metrics.cohort_gate import effective_tests
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from w0s_00_artifacts import load_cells                      # noqa: E402

OUT = Path(os.environ.get(
    "W0S_OUT_ANALYSIS",
    ROOT / "data" / "paper_final" / "lane_s_scale"))

NOCT, NQ = 8, 5


def main() -> None:
    cs = load_cells()
    OUT.mkdir(parents=True, exist_ok=True)
    s = cs.s
    print(f"[w0s-dilution] {len(cs.bands)} bands x {cs.n_scales} scales "
          f"| readouts {len(cs.readouts)}", flush=True)

    # ------------------------------------------------------------------ #
    # A/B/C -- what the readout sees at each scale
    # ------------------------------------------------------------------ #
    rows = []
    for b in cs.bands:
        neff = cs.diag(b, "n_eff")
        mcom = cs.diag(b, "m_comm")
        ndis = cs.diag(b, "n_distinct")
        pfr = cs.diag(b, "pairfrac")
        xs_c = cs.diag(b, "xs_coph")
        xs_r = cs.diag(b, "xs_reorg")
        for j in range(cs.n_scales):
            r = dict(band=b, s=float(s[j]),
                     n_eff=float(np.nanmedian(neff[:, j])),
                     m_comm=float(np.nanmedian(mcom[:, j])),
                     n_distinct=float(np.nanmedian(ndis[:, j])),
                     ostar=float(np.nanmedian(cs.diag(b, "ostar")[:, j])),
                     qstar=float(np.nanmedian(cs.diag(b, "qstar")[:, j])))
            # similarity of the input to the FINEST scale and to the NEXT scale
            r["xs_coph_vs_s0"] = float(np.nanmedian(xs_c[:, j, 0]))
            r["xs_reorg_vs_s0"] = float(np.nanmedian(xs_r[:, j, 0]))
            if j + 1 < cs.n_scales:
                r["xs_coph_adjacent"] = float(np.nanmedian(xs_c[:, j, j + 1]))
                r["xs_reorg_adjacent"] = float(np.nanmedian(xs_r[:, j, j + 1]))
            for o in range(NOCT):
                r[f"pairfrac_o{o+1}"] = float(np.nanmedian(pfr[:, j, o]))
            rows.append(r)
    df_scale = pd.DataFrame(rows)
    df_scale.to_csv(OUT / "dilution_by_scale.csv", index=False)

    # whole cross-scale matrices, cohort median, for the figure and the head
    xs_rows = []
    for b in cs.bands:
        for key in ("xs_coph", "xs_reorg"):
            M = np.nanmedian(cs.diag(b, key), axis=0)
            iu = np.triu_indices(cs.n_scales, 1)
            xs_rows.append(dict(band=b, quantity=key,
                                mean_offdiag=float(np.nanmean(M[iu])),
                                min_offdiag=float(np.nanmin(M[iu])),
                                corr_first_last=float(M[0, -1])))
    pd.DataFrame(xs_rows).to_csv(OUT / "dilution_input_similarity.csv", index=False)

    # ------------------------------------------------------------------ #
    # D -- contribution share by tree level, per scale
    # ------------------------------------------------------------------ #
    share_rows = []
    for b in cs.bands:
        for tag, nlab, pre in (("octave", NOCT, "Ccon_o"),
                               ("quintile", NQ, "Qcon_q")):
            con = np.stack([cs.cell(b, f"{pre}{i+1}")[0] for i in range(nlab)])
            mar = np.stack([cs.margins(b, f"{pre}{i+1}")[0] for i in range(nlab)])
            #  (nlab, K, nS): observed knob-integrated contributions, and the
            #  same net of each patient's own matched-strength null. The observed
            #  shares answer "what is the statistic made of"; the margin shares
            #  answer "which tree levels carry the excess over the null", and the
            #  two are different questions.
            tot = np.nansum(np.abs(con), axis=0)
            tot_m = np.nansum(np.abs(mar), axis=0)
            for i in range(nlab):
                for j in range(cs.n_scales):
                    with np.errstate(invalid="ignore", divide="ignore"):
                        sh = np.abs(con[i, :, j]) / tot[:, j]
                        shm = np.abs(mar[i, :, j]) / tot_m[:, j]
                    share_rows.append(dict(
                        band=b, stratification=tag, stratum=i + 1,
                        s=float(s[j]),
                        contribution=float(np.nanmedian(con[i, :, j])),
                        abs_share=float(np.nanmedian(sh)),
                        margin=float(np.nanmedian(mar[i, :, j])),
                        abs_share_margin=float(np.nanmedian(shm))))
    df_share = pd.DataFrame(share_rows)
    df_share.to_csv(OUT / "dilution_contribution_shares.csv", index=False)

    # ------------------------------------------------------------------ #
    # E -- the pre-registered §5.1 falsification check
    # ------------------------------------------------------------------ #
    ver_rows = []
    for b in cs.bands:
        inc_M, _ = cs.margins(b, "T")
        inc = effective_tests(inc_M)
        for tag, nlab, cpre, lpre in (("octave", NOCT, "Ccon_o", "Tloc_o"),
                                      ("quintile", NQ, "Qcon_q", "Qloc_q")):
            sub = df_share[(df_share.band == b)
                           & (df_share.stratification == tag)]
            piv = sub.pivot(index="s", columns="stratum", values="abs_share")
            ratio = (piv.max(axis=1) / piv.median(axis=1)).to_numpy()
            frac_flat = float(np.nanmean(ratio <= 2.0))
            offd = []
            for i in range(nlab):
                M, _ = cs.margins(b, f"{lpre}{i+1}")
                e = effective_tests(M)
                if np.isfinite(e["mean_offdiag"]) and e["n_finite"] >= 5:
                    offd.append(e["mean_offdiag"])
            med_off = float(np.median(offd)) if offd else np.nan
            ver_rows.append(dict(
                band=b, stratification=tag,
                incumbent_mean_offdiag=inc["mean_offdiag"],
                incumbent_n_eff_pr=inc["n_eff_pr"],
                share_ratio_median=float(np.nanmedian(ratio)),
                share_ratio_max=float(np.nanmax(ratio)),
                frac_scales_flat=frac_flat,
                stratum_mean_offdiag_median=med_off,
                n_strata_estimated=len(offd),
                cond1_shares_flat=bool(frac_flat > 0.5),
                cond2_strata_as_correlated=bool(
                    np.isfinite(med_off)
                    and abs(med_off - inc["mean_offdiag"]) <= 0.10),
            ))
    df_ver = pd.DataFrame(ver_rows)
    df_ver["dilution_falsified"] = (df_ver.cond1_shares_flat
                                    & df_ver.cond2_strata_as_correlated)
    df_ver.to_csv(OUT / "dilution_falsification_check.csv", index=False)

    # ------------------------------------------------------------------ #
    print("\n-- input similarity along the scale axis (cohort median) --")
    print(pd.DataFrame(xs_rows).to_string(index=False))
    print("\n-- §5.1 falsification check --")
    print(df_ver[["band", "stratification", "share_ratio_median",
                  "frac_scales_flat", "incumbent_mean_offdiag",
                  "stratum_mean_offdiag_median", "dilution_falsified"]]
          .to_string(index=False))
    print(f"\n[w0s-dilution] -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
