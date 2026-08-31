#!/usr/bin/env python3
"""Shared reader for the lane W0-S grid: cells on disk -> gate-ready arrays.

One place decides how a per-cell ``.npz`` becomes the ``(obs, surr)`` pair the
locked cohort gate consumes, so the diagnostics stage and the verdict stage
cannot disagree about the knob integration or the NaN policy.

Knob integration, stated once
-----------------------------
Every fraction of the locked plateau ``f in {0.07, 0.10, 0.14, 0.20}`` is a
legitimate reading of the same substrate, so no number is reported at one of
them. The knob-integrated statistic is the **median over fractions**, applied
identically to the observed graph and to every surrogate realization:

    O[k]    = median_f obs[k, f, s]
    S[k, r] = median_f surr[k, r, f, s]

The margin is then exactly ``patient_margin(O, S)`` from
:mod:`lrg_eegfc.utils.metrics.cohort_gate` -- the locked gate's own definition,
unmodified. Taking the median over the knob *after* the surrogate subtraction
would be an equally defensible order, but medians do not commute, and this order
is the one that leaves the gate's contract exact rather than approximated. The
across-fraction spread is carried alongside as the knob uncertainty.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, PATIENTS_4PHASE

__all__ = ["CellSet", "load_cells"]


class CellSet:
    """The whole grid in memory, addressed by (band, readout).

    Attributes
    ----------
    s : (nS,) scale grid
    fracs : (nF,) plateau fractions
    readouts : list[str]
    patients : list[str]
    bands : list[str]
    """

    def __init__(self, root: Path, patients: Sequence[str], bands: Sequence[str]):
        self.root = Path(root)
        self.patients = list(patients)
        self.bands = list(bands)
        self._obs: dict[tuple[str, str], np.ndarray] = {}
        self._surr: dict[tuple[str, str], np.ndarray] = {}
        self._diag: dict[tuple[str, str], dict] = {}
        self._cache: dict[tuple[str, str], tuple] = {}
        self.s = self.fracs = None
        self.readouts: list[str] = []
        self.N: dict[tuple[str, str], int] = {}
        for b in self.bands:
            for p in self.patients:
                f = self.root / "cells" / f"{p}__{b}.npz"
                if not f.exists():
                    continue
                z = np.load(f, allow_pickle=False)
                if self.s is None:
                    self.s = z["s"]
                    self.fracs = z["fracs"]
                    self.readouts = [str(x) for x in z["readouts"]]
                self._obs[(p, b)] = z["obs"]                 # (nF, nS, nRead)
                self._surr[(p, b)] = z["surr"]               # (R, nF, nS, nRead)
                self._diag[(p, b)] = {k: z[k] for k in
                                      ("n_eff", "m_comm", "n_distinct",
                                       "pairfrac", "ostar", "qstar",
                                       "xs_coph", "xs_reorg") if k in z}
                self.N[(p, b)] = int(z["N"][0])
        if self.s is None:
            raise FileNotFoundError(f"no cells under {self.root/'cells'}")
        self.ix = {k: i for i, k in enumerate(self.readouts)}
        self.n_scales = int(self.s.size)
        # Common surrogate depth. Cells may carry different R on disk (an early
        # pass ran deeper before the shared machine's contention forced the
        # ensemble down), and an unequal R would give some patients a
        # better-estimated null floor than others in the very statistic --
        # obs minus own-surrogate-median -- that exists to equalise them. Every
        # cell is therefore truncated to the smallest R present, so the extra
        # draws are retained on disk but never give one cell an advantage.
        self.R = min(v.shape[0] for v in self._surr.values())
        for k in list(self._surr):
            if self._surr[k].shape[0] > self.R:
                self._surr[k] = self._surr[k][: self.R]

    # ------------------------------------------------------------------ #
    def have(self, band: str) -> list[str]:
        """Patients with a cell on disk for this band, in cohort order."""
        return [p for p in self.patients if (p, band) in self._obs]

    def cell(self, band: str, readout: str) -> tuple[np.ndarray, np.ndarray, list[str]]:
        """Knob-integrated ``(O, S, labels)``: ``O`` is ``(K, nS)``, ``S`` is ``(K, R, nS)``.

        Memoised: the held-out-realization loops ask for the same cell hundreds
        of times, and re-reducing the ``(R, nF, nS)`` block each time dominates
        the analysis runtime for no benefit.
        """
        key = (band, readout)
        if key in self._cache:
            return self._cache[key]
        m = self.ix[readout]
        labs = self.have(band)
        O = np.stack([np.nanmedian(self._obs[(p, band)][:, :, m], axis=0)
                      for p in labs])
        S = np.stack([np.nanmedian(self._surr[(p, band)][:, :, :, m].astype(float),
                                   axis=1) for p in labs])
        self._cache[key] = (O, S, labs)
        return O, S, labs

    def knob_spread(self, band: str, readout: str) -> np.ndarray:
        """Across-fraction IQR of the observed statistic, ``(K, nS)`` -- knob uncertainty."""
        m = self.ix[readout]
        labs = self.have(band)
        return np.stack([
            np.nanpercentile(self._obs[(p, band)][:, :, m], 75, axis=0)
            - np.nanpercentile(self._obs[(p, band)][:, :, m], 25, axis=0)
            for p in labs])

    def margins(self, band: str, readout: str) -> tuple[np.ndarray, list[str]]:
        """Per-patient margin matrix ``(K, nS)`` = ``O - median_r S``."""
        O, S, labs = self.cell(band, readout)
        return O - np.nanmedian(S, axis=1), labs

    def heldout_margins(self, band: str, readout: str, draw: int
                        ) -> np.ndarray:
        """Margin matrix with surrogate realization ``draw`` promoted to observed.

        Under the null the true observation is exchangeable with its own
        surrogates, so this is a draw from the readout's null margin matrix --
        with the real cohort size, the real per-patient heteroscedasticity and,
        crucially, the readout's own noise level. It is what makes a
        cross-scale-independence statistic interpretable: a noisier readout has
        a higher effective number of independent scales for free, and only the
        comparison against these draws can tell that apart from scale-specific
        signal.
        """
        _, S, _ = self.cell(band, readout)
        keep = np.ones(S.shape[1], bool)
        keep[draw] = False
        return S[:, draw, :] - np.nanmedian(S[:, keep, :], axis=1)

    def diag(self, band: str, key: str) -> np.ndarray:
        """Cohort stack of an observed-only diagnostic, knob-median over fractions."""
        labs = self.have(band)
        return np.stack([np.nanmedian(self._diag[(p, band)][key], axis=0)
                         for p in labs])


def load_cells(root: Optional[Path] = None,
               patients: Optional[Sequence[str]] = None,
               bands: Optional[Sequence[str]] = None) -> CellSet:
    """Load the grid written by ``w0s_01_scale_locality_grid.py``."""
    from lrg_eegfc.utils.scripting import setup_script_env
    r = Path(root) if root is not None else (
        setup_script_env() / "data" / "paper_final" / "lane_s_scale" / "grid")
    return CellSet(r, patients or PATIENTS_4PHASE, bands or BRAIN_BANDS_NAMES)
