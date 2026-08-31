import sys
sys.path.insert(0, "scripts/01_compute/paper_final")
from lrg_eegfc.visuals import imshow_colorbar_caxdivider  # noqa
from lrg_eegfc.visuals.styles import band_color, use_lrg_style  # noqa
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
print("imports OK", band_color("beta"), BRAIN_BAND_TEX_DICT["beta"])
import w0s_05_figures  # noqa
print("figure module OK")
