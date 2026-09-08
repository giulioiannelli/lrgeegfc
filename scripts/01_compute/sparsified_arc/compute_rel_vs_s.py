#!/usr/bin/env python3
r"""compute_rel_vs_s -- cache the per-scale split-half reliability rel(s) of the
cophenetic distance, cohort-median per band, on the mst@0.20 backbone.

rel(s) = Spearman(D_A(s), D_B(s)) with A,B the two rest_pre split-halves -- the
model-free reachable ceiling for any cophenetic concordance at scale s (you cannot
correlate a phase-shift with persistence better than the baseline correlates with
itself). This is the SAME reliability anchor the measure-slide bars are cut at
(reachable max = rho(A,B) ~ 0.5, NEVER 1); here it is resolved across the diffusion
scale so a scale-resolved concordance curve can be read as a fraction of what is
reliably measurable at that scale.

Reuses the eig / backbone / cophenetic machinery from 05_enc_inf_arc (no forks).

Writes: data/sparsified_arc/enc_inf_arc_mst020/rel_vs_s.csv   (band, s, rel_median, rel_iqr_lo, rel_iqr_hi, n)
"""
from __future__ import annotations
import importlib.util
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
_spec = importlib.util.spec_from_file_location(
    "enc05", ROOT / "scripts" / "01_compute" / "sparsified_arc" / "05_enc_inf_arc.py")
_m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_m)
from lrg_eegfc.utils.fc.heat_multiscale import scale_grid, cophenetic_at_scale  # noqa: E402

OUT = _m.OUT / "rel_vs_s.csv"
N_S = _m.N_S


def main():
    rows = []
    for band in _m.BANDS:
        RELS, s_axis = [], None
        for pat in _m.COHORT:
            try:
                WA, WB = _m.load_phase(pat, "A", band), _m.load_phase(pat, "B", band)
            except Exception:
                continue
            ea, eb = _m._eig(WA), _m._eig(WB)
            sg = scale_grid(ea[0], n=N_S)
            s_axis = sg
            RELS.append([spearmanr(cophenetic_at_scale(*ea, s),
                                   cophenetic_at_scale(*eb, s)).correlation for s in sg])
        if not RELS:
            continue
        RELS = np.asarray(RELS)
        med = np.nanmedian(RELS, axis=0)
        lo = np.nanquantile(RELS, 0.25, axis=0)
        hi = np.nanquantile(RELS, 0.75, axis=0)
        for j, s in enumerate(s_axis):
            rows.append(dict(band=band, s=float(s), rel_median=float(med[j]),
                             rel_iqr_lo=float(lo[j]), rel_iqr_hi=float(hi[j]),
                             n=int(np.isfinite(RELS[:, j]).sum())))
        print(f"  {band:9s} rel(s): {med.min():.2f}-{med.max():.2f} "
              f"(median {np.median(med):.2f}, n={RELS.shape[0]})", flush=True)
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"\nwrote {OUT}  ({len(df)} rows)", flush=True)


if __name__ == "__main__":
    main()
