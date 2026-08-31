"""Throwaway: synthesize a grid with the real shapes so the analysis stages can
be debugged before the real cells land. Not part of the pipeline."""
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0, "scripts/01_compute/paper_final")
from w0s_01_scale_locality_grid import READOUTS, SGRID, FRACS, NOCT, NQ  # noqa

OUT = Path("scratch/fake_grid/cells")
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(0)
nS, nF, nR = SGRID.size, len(FRACS), 20
bands = ["theta", "beta"]
pats = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
        "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
for b in bands:
    eff = 0.06 if b == "beta" else 0.0
    for p in pats:
        obs = rng.normal(0.1, 0.05, (nF, nS, len(READOUTS))) + eff
        surr = rng.normal(0.1, 0.05, (nR, nF, nS, len(READOUTS))).astype(np.float32)
        # starve the top octaves like the real thing
        obs[:, :, [READOUTS.index(f"Tloc_o{o}") for o in (7, 8)]] = np.nan
        surr[:, :, :, [READOUTS.index(f"Tloc_o{o}") for o in (7, 8)]] = np.nan
        d = dict(
            n_eff=np.tile(np.linspace(118, 4, nS), (nF, 1)),
            m_comm=np.tile(np.linspace(1, 113, nS), (nF, 1)),
            n_distinct=np.tile(np.linspace(20, 3, nS), (nF, 1)),
            pairfrac=rng.random((nF, nS, NOCT)),
            ostar=np.tile(np.linspace(6, 1, nS), (nF, 1)),
            qstar=np.tile(np.linspace(5, 1, nS), (nF, 1)),
            xs_coph=np.tile(np.eye(nS), (nF, 1, 1)),
            xs_reorg=np.tile(np.eye(nS), (nF, 1, 1)),
        )
        np.savez_compressed(OUT / f"{p}__{b}.npz", s=SGRID,
                            fracs=np.array(FRACS),
                            readouts=np.array(READOUTS), obs=obs, surr=surr,
                            N=np.array([118]), **d)
print("fake grid written", OUT)
