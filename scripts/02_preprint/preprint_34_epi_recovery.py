#!/usr/bin/env python3
"""Within-patient epi-node recovery — matched-strength CERTIFICATION (honest).

Supersedes the first version of this figure, which compared the propagator
subspace to an UNFAIR strength baseline (template-distance on a monotone
feature). This version reads the audit_86 certification, which scores strength
fairly (within-patient Mann-Whitney AUC) and removes strength from the
propagator features via the matched-strength surrogate.

Two panels, all six bands.

(a) **Recovery AUC, three honest scores** (chance = 0.5, dashed): fair node
    strength (Mann-Whitney) vs the raw propagator ρ(τ) subspace vs the
    strength-REMOVED (matched-strength surrogate-z) propagator subspace.
    Reading: *strength wins at β (hubness); the strength-removed propagator only
    exceeds strength at α — the one matched-strength-certified, strength-
    orthogonal signal.*

(b) **Strength-orthogonal margin** Δ = (strength-removed propagator recovery −
    fair strength recovery), per band (0 line). Reading: *α is the only band
    where the propagator carries epi information beyond hubness; β is negative
    (hubness dominates).*

No in-axes numeric text (project rule). Reads only the audit_86 certify CSV.

Input  : data/audit/epi_marker/propagator_ms_certify.csv
Output : data/preprint/figures/all_bands/fig_epi_recovery.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
SCORES = ["strength_mono", "diffusion_raw", "diffusion_msz"]
SC_C = {"strength_mono": "#2c5f8a", "diffusion_raw": "#8a8a8a",
        "diffusion_msz": "#d1495b"}
SC_L = {"strength_mono": "node strength (fair, Mann-Whitney)",
        "diffusion_raw": r"propagator $\rho(\tau)$ subspace (raw)",
        "diffusion_msz": r"propagator, strength-removed (matched-strength $z$)"}


def main() -> None:
    src = ROOT / "data" / "audit" / "epi_marker" / "propagator_ms_certify.csv"
    if not src.exists():
        raise SystemExit("[preprint_34] run audit_86 first")
    c = pd.read_csv(src)
    rec = c[c.kind == "recovery"].copy()
    rec["feat"] = rec.feature.str.replace("__recovery_", "", regex=False)

    def val(band, feat):
        r = rec[(rec.band == band) & (rec.feat == feat)]
        return float(r.median_auc.iloc[0]) if not r.empty else np.nan

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(13.0, 4.9))

    # (a) three-way recovery bars
    x = np.arange(len(BANDS)); w = 0.8 / len(SCORES)
    for si, s in enumerate(SCORES):
        vals = [val(b, s) for b in BANDS]
        axa.bar(x + (si - (len(SCORES) - 1) / 2) * w, vals, width=w,
                color=SC_C[s], edgecolor="none", zorder=3)
    axa.axhline(0.5, color="0.35", lw=1.1, ls="--", zorder=4)
    axa.text(-0.45, 0.5, "chance", color="0.35", va="bottom", ha="left",
             fontsize=9)
    axa.set_xticks(x)
    axa.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS])
    axa.set_ylabel("within-patient recovery AUC")
    axa.set_ylim(0.45, 0.78)
    axa.spines[["top", "right"]].set_visible(False)
    axa.set_title("(a) recovery: fair strength vs propagator (raw / strength-removed)",
                  fontsize=10.5, loc="left")

    # (b) strength-orthogonal margin Δ = msz − fair strength
    deltas = [val(b, "diffusion_msz") - val(b, "strength_mono") for b in BANDS]
    cols = ["#d1495b" if d > 0 else "#888888" for d in deltas]
    axb.bar(x, deltas, width=0.6, color=cols, edgecolor="none", zorder=3)
    axb.axhline(0.0, color="0.35", lw=1.1, zorder=4)
    axb.set_xticks(x)
    axb.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS])
    axb.set_ylabel(r"$\Delta$ AUC (strength-removed propagator $-$ fair strength)")
    axb.spines[["top", "right"]].set_visible(False)
    axb.set_title("(b) strength-orthogonal margin (positive = beats hubness)",
                  fontsize=10.5, loc="left")

    h = [Patch(facecolor=SC_C[s], edgecolor="none", label=SC_L[s])
         for s in SCORES]
    h.append(Line2D([], [], color="0.35", lw=1.1, ls="--", label="chance"))
    fig.legend(handles=h, loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, -0.10), fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 1))

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_epi_recovery.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[preprint_34] -> {out}")


if __name__ == "__main__":
    main()
