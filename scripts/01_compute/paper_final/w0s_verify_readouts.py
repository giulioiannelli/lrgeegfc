#!/usr/bin/env python3
"""Equivalence checks for lane W0-S, on real FC -- run before the grid.

Four things must be true before any verdict is worth reading, and each is
checked against the library rather than asserted in prose:

1. The graphs this lane builds are the graphs the locked contract builds.
   ``select_backbone(dense, CANONICAL.backbone, frac=f)`` must be bit-identical
   to ``canonical_graph_ensemble``'s yield at every plateau fraction.
2. The lane's ``T`` is the incumbent trace. It must agree with
   ``rho_sym_over_scales`` and with ``cross_phase_functionals_over_scales``'s
   ``T_probe`` to float precision at every scale.
3. The octave decomposition is exact: ``sum_o Ccon(s, o) == T(s)``.
4. The equal-count decomposition is exact: ``sum_q Qcon(s, q) == T(s)``.

Anything that fails here invalidates the grid, so this runs first and prints the
worst discrepancy rather than only a pass/fail.
"""
from __future__ import annotations

import numpy as np

from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import (
    cross_phase_functionals_over_scales,
    laplacian_eig,
    rho_sym_over_scales,
)
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.substrate import (
    CANONICAL,
    canonical_graph,
    canonical_graph_ensemble,
)

setup_script_env()

from w0s_01_scale_locality_grid import (            # noqa: E402
    FRACS, IX, NOCT, NQ, PHASES, SGRID, realization,
)

PAT, BAND = "Pat_02", "beta"


def main() -> None:
    print(f"[w0s-verify] {PAT}/{BAND} | substrate {CANONICAL.label()} "
          f"| fracs {FRACS}", flush=True)

    # 1 -- the graph is the locked graph -------------------------------------- #
    worst_g = 0.0
    for ph in PHASES:
        dense = canonical_graph(PAT, ph, BAND, dense=True)
        for (f_ens, G_ens), f in zip(
                canonical_graph_ensemble(PAT, ph, BAND), FRACS):
            assert abs(f_ens - f) < 1e-12, f"fraction order {f_ens} vs {f}"
            G_mine = select_backbone(dense, CANONICAL.backbone, frac=float(f))
            worst_g = max(worst_g, float(np.max(np.abs(G_mine - G_ens))))
    print(f"  1. graph == canonical_graph_ensemble       max|diff| = {worst_g:.3g}")
    assert worst_g == 0.0

    # 2/3/4 -- the readouts ---------------------------------------------------- #
    Ws = {ph: canonical_graph(PAT, ph, BAND, dense=True) for ph in PHASES}
    vals, _ = realization(Ws, want_diag=False)

    worst_rho = worst_probe = worst_o = worst_q = 0.0
    for fi, f in enumerate(FRACS):
        eig = {ph: laplacian_eig(select_backbone(Ws[ph], CANONICAL.backbone,
                                                 frac=float(f)))
               for ph in PHASES}
        ref_rho = rho_sym_over_scales(eig, SGRID)
        ref_probe = cross_phase_functionals_over_scales(eig, SGRID)["T_probe"]
        mine = vals[fi, :, IX["T"]]
        ok = np.isfinite(mine) & np.isfinite(ref_rho)
        worst_rho = max(worst_rho, float(np.max(np.abs(mine[ok] - ref_rho[ok]))))
        worst_probe = max(worst_probe,
                          float(np.max(np.abs(mine[ok] - ref_probe[ok]))))
        so = np.nansum([vals[fi, :, IX[f"Ccon_o{o}"]] for o in range(1, NOCT + 1)],
                       axis=0)
        sq = np.nansum([vals[fi, :, IX[f"Qcon_q{q}"]] for q in range(1, NQ + 1)],
                       axis=0)
        worst_o = max(worst_o, float(np.max(np.abs(so[ok] - mine[ok]))))
        worst_q = max(worst_q, float(np.max(np.abs(sq[ok] - mine[ok]))))

    print(f"  2. T == rho_sym_over_scales                max|diff| = {worst_rho:.3g}")
    print(f"     T == T_probe                            max|diff| = {worst_probe:.3g}")
    print(f"  3. sum_o Ccon == T                         max|diff| = {worst_o:.3g}")
    print(f"  4. sum_q Qcon == T                         max|diff| = {worst_q:.3g}")
    assert worst_rho < 1e-12 and worst_probe < 1e-12
    assert worst_o < 1e-10 and worst_q < 1e-10
    print("[w0s-verify] all equivalence checks passed", flush=True)


if __name__ == "__main__":
    main()
