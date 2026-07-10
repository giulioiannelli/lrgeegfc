"""
M-1 explainer — how ImCoh is computed, and the band average.

For one representative contact pair, show the signed frequency-resolved
imaginary coherence ImCoh(f) (Nolte 2004, in [-1, 1], oscillating in sign),
its per-bin absolute value |ImCoh(f)|, and the in-band mean that becomes the
single number we feed the pipeline:

    imcoh_abs_ij = < |ImCoh_ij(f)| >_band     (|.| per bin FIRST, then average)

Vector PDF, transparent background, light ink → sits on a DARK slide.
Output: data/outputs/figures/talk/imcoh_band_average.pdf
QA (dark-matte PNG, opt-in): ... --qa /path/to/qa.png
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.config.paths import IMCOH_CACHE, FIGURES_ROOT
from lrg_eegfc.config.const import BRAIN_BANDS

use_lrg_style()

# ---- palette (dark-slide, light ink) --------------------------------------
INK = "#e6edf3"
MUTED = "#8b949e"
SIGNED = "#5aa9e6"    # signed ImCoh(f)
ABSC = "#e0a44a"     # |ImCoh(f)|
MEANC = "#e5006e"    # band average

PAT, BAND, PHASE, FS, NPERSEG = "Pat_05", "beta", "rest_pre", 2048, 4096

# ---- load freq-resolved SIGNED ImCoh (N, N, F_band) -----------------------
arr = np.load(IMCOH_CACHE / PAT / f"{BAND}_{PHASE}_imcoh_freqresolved_nperseg-{NPERSEG}.npy")
n, _, nf = arr.shape

# in-band frequency axis (rfft grid, restricted to the band)
lo, hi = BRAIN_BANDS[BAND]
freqs_full = np.fft.rfftfreq(NPERSEG, d=1.0 / FS)
band_freqs = freqs_full[(freqs_full >= lo) & (freqs_full <= hi)]
if band_freqs.size != nf:                       # fall back if the cache trimmed edges
    band_freqs = np.linspace(lo, hi, nf)

# pick a strong OFF-shaft pair whose ImCoh(f) CROSSES ZERO in-band, so the
# per-bin |.| visibly differs from the signed spectrum (the whole point).
absmat = np.abs(arr).mean(axis=2)               # <|ImCoh|>_f  (what we feed)
signedmean = arr.mean(axis=2)                    # <ImCoh>_f    (signed mean)
gap = absmat - np.abs(signedmean)                # Jensen gap; large ⇒ sign-crossing
iu = np.triu_indices(n, k=6)                     # off-shaft pairs only
strong = absmat[iu] >= np.percentile(absmat[iu], 85)
cand = np.where(strong)[0]
kbest = cand[np.argmax(gap[iu][cand])]           # strong AND most sign-crossing
i, j = int(iu[0][kbest]), int(iu[1][kbest])

signed = arr[i, j]                              # ImCoh(f), signed
absval = np.abs(signed)                         # |ImCoh(f)|
mean_abs = absval.mean()                         # the number we feed the pipeline

# ---- plot -----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.6, 4.3))
ax.axhline(0, color=MUTED, lw=0.9, alpha=0.7)
ax.plot(band_freqs, signed, color=SIGNED, lw=1.8, marker="o", ms=3.2,
        label="signed  ImCoh(f)")
ax.fill_between(band_freqs, 0, absval, color=ABSC, alpha=0.28)
ax.plot(band_freqs, absval, color=ABSC, lw=1.8, label="|ImCoh(f)|")
ax.axhline(mean_abs, color=MEANC, lw=2.0, ls="--",
           label=r"$\langle |\mathrm{ImCoh}| \rangle_{\beta}$" + f" = {mean_abs:.2f}")

ax.set_xlim(lo, hi)
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("ImCoh")
ax.set_title("")
leg = ax.legend(loc="lower right", frameon=False, fontsize=10)
for t in leg.get_texts():
    t.set_color(INK)

# light-ink styling for a dark slide
for spine in ax.spines.values():
    spine.set_color(MUTED)
ax.tick_params(colors=INK)
ax.xaxis.label.set_color(INK)
ax.yaxis.label.set_color(INK)

fig.tight_layout()

out = FIGURES_ROOT / "talk" / "imcoh_band_average.pdf"
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, transparent=True, bbox_inches="tight")
print(f"wrote {out}  (pair i={i}, j={j}, <|ImCoh|>={mean_abs:.3f}, "
      f"signed range [{signed.min():.2f}, {signed.max():.2f}], nf={nf})")

if "--qa" in sys.argv:
    qa = sys.argv[sys.argv.index("--qa") + 1]
    fig.savefig(qa, facecolor="#0b0e17", bbox_inches="tight", dpi=130)
    print("wrote QA", qa)

plt.close(fig)
