#!/usr/bin/env python3
"""Acceptance check for the substrate entry point every downstream lane calls.

Asserts on real FC (hypothesis-level code never mocks FC):
  (a) canonical_graph(...) reproduces the incumbent inline pipeline
      (load_fc_matrix -> clip/symmetrise -> mst_union_top_fraction) bit-for-bit;
  (b) the five-phase path is a supported, tested path: canonical_phase_eigs
      returns all five phases and drives cross_phase_functionals_over_scales to
      all four functionals;
  (c) the four-phase path (task_learn absent) works through the SAME call and
      returns exactly {T_probe}, equal to the five-phase T_probe;
  (d) split-half phases resolve from the canonical halves cache, and a
      "<phase>_A" style name resolves too;
  (e) substrate overrides (transform, backbone, frac) take effect and the
      knob-integrated ensemble yields the declared plateau fractions.
"""
from __future__ import annotations
import sys
import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.paths import IMCOH_HALVES_CACHE
from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction, tmfg_backbone
from lrg_eegfc.utils.fc.heat_multiscale import cross_phase_functionals_over_scales
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.substrate import (
    CANONICAL, CANONICAL_PHASES, canonical_eig, canonical_graph,
    canonical_graph_ensemble, canonical_phase_eigs, canonical_scale_grid,
    canonical_structure,
)

fails = []


def check(name, ok, detail=""):
    print(f"  [{'OK ' if ok else 'FAIL'}] {name}{'  ' + detail if detail else ''}", flush=True)
    if not ok:
        fails.append(name)


def main():
    pat, band = "Pat_02", "beta"
    print(f"substrate = {CANONICAL.label()}  (parameter_free={CANONICAL.is_parameter_free})",
          flush=True)

    # (a) bit-for-bit vs the incumbent inline pipeline
    W = np.asarray(load_fc_matrix(pat, "task_test", band, fc_method="imcoh_abs"), float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    W = 0.5 * (W + W.T)
    ref = mst_union_top_fraction(W, 0.20)
    got = canonical_graph(pat, "task_test", band, backbone="mst", frac=0.20,
                          transform="abs")
    check("(a) reproduces incumbent inline pipeline",
          np.array_equal(ref, got), f"max|diff|={np.abs(ref-got).max():.2e}")

    # (b) five-phase path
    eig5 = canonical_phase_eigs(pat, band)
    s = canonical_scale_grid()
    f5 = cross_phase_functionals_over_scales(eig5, s)
    check("(b) five-phase path", sorted(eig5) == sorted(CANONICAL_PHASES)
          and sorted(f5) == ["T_encode", "T_probe", "T_probespec", "T_probespec_pe"],
          f"phases={sorted(eig5)} functionals={sorted(f5)}")

    # (c) four-phase path through the same call
    eig4 = canonical_phase_eigs(pat, band,
                                phases=("A", "B", "task_test", "rest_post"))
    f4 = cross_phase_functionals_over_scales(eig4, s)
    check("(c) four-phase path == five-phase T_probe",
          sorted(f4) == ["T_probe"]
          and np.nanmax(np.abs(f4["T_probe"] - f5["T_probe"])) < 1e-12,
          f"keys={sorted(f4)} max|diff|={np.nanmax(np.abs(f4['T_probe']-f5['T_probe'])):.2e}")

    # (d) split-half resolution, both spellings
    a1 = canonical_graph(pat, "A", band)
    a2 = canonical_graph(pat, "rest_pre_A", band)
    ref_half = np.load(IMCOH_HALVES_CACHE / pat / f"{band}_rest_pre_A_imcoh_abs.npy")
    check("(d) split-half phases resolve ('A' == 'rest_pre_A')",
          np.array_equal(a1, a2) and a1.shape == np.asarray(ref_half).shape)

    # (e) overrides + knob-integrated ensemble
    g_sq = canonical_graph(pat, "task_test", band, transform="sq")
    g_tmfg = canonical_graph(pat, "task_test", band, backbone="tmfg")
    ok_sq = not np.allclose(g_sq, got)
    ok_tmfg = np.array_equal(g_tmfg, tmfg_backbone(W))
    plateau = CANONICAL.with_(plateau=(0.10, 0.28), plateau_fracs=(0.10, 0.14, 0.20, 0.28))
    fr = [f for f, _ in canonical_graph_ensemble(pat, "task_test", band, substrate=plateau)]
    check("(e) overrides + ensemble", ok_sq and ok_tmfg and fr == [0.10, 0.14, 0.20, 0.28],
          f"transform override changes graph={ok_sq}, tmfg exact={ok_tmfg}, fracs={fr}")

    st = canonical_structure(pat, "task_test", band)
    print(f"  structure: density={st['density']:.4f} mean_deg={st['mean_degree']:.1f} "
          f"ncomp={st['n_components']} clustering={st['clustering']:.4f} "
          f"lam2/lam_max={st['spectral_gap']:.4g} wfrac={st['weight_fraction']:.3f}", flush=True)
    _ = canonical_eig(pat, "A", band)

    print("\nALL ACCEPTANCE CHECKS PASSED" if not fails else f"\nFAILED: {fails}", flush=True)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
