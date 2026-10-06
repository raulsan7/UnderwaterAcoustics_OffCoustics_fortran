"""
Module: plot_utils/style.py
Description: Visual asthetics (fonts, colormaps, colors, line styles) shared by every plot. Text is rendered
             with LaTeX; if no working LaTeX installation is found (e.g. on a cluster) it falls back to
             Matplotlib's mathtext (text.usetex = False) with a warning.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 09/2026
"""

import io
import shutil
import matplotlib.pyplot as plt


class AcademicStyle:
    """
    Configures visual asthetics for scientific publications.
    """

    # Fonts
    FONT_FAMILY   = "serif"
    DEFAULT_CMAP  = "inferno"
    PRESSURE_CMAP = "seismic"
    CASE_COLORS   = ["darkblue", "darkorange", "teal"]    # Monopile, Floating (shallow), Floating (deep)
    BAND_COLORS   = ["steelblue", "tomato", "seagreen", "mediumpurple", "goldenrod", "darkcyan", "salmon", "olivedrab", "slategray", "peru"]
    LINE_STYLES   = ['-', '--', ':', '-.', '-o', '-v', '--o', ':^']
    WS_CMAP       = "viridis"                             # Wind speeds (--ws): colormap sampled on WS_CMAP_RANGE, low to high speed
    WS_CMAP_RANGE = (0.0, 0.8)                            # (the bright yellow end of viridis does not read on white)
    WS_LW_RANGE   = (1.2, 4.0)                            # Wind speeds (--ws): line width from low to high speed


    P_REF = 1e-6            # [Pa] Reference Pressure for water

    @classmethod
    def apply(cls):
        """
        Applies global configuration to Matplotlib
        """
        usetex = cls.latex_available()

        plt.rcParams.update({
            "font.family": cls.FONT_FAMILY,
            "text.usetex": usetex,
            "font.size": 11,
            "axes.labelsize": 16,
            "axes.titlesize": 20,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize":14,
            "figure.titlesize": 20,
            "savefig.dpi": 300,
            "savefig.bbox": 'tight'
        })

        # Without LaTeX, at least the math parts keep the Computer Modern look
        if not usetex:
            plt.rcParams["mathtext.fontset"] = "cm"

    @staticmethod
    def latex_available():
        """
        Checks whether Matplotlib can render text with LaTeX: latex and dvipng must be on the PATH and a test
        string with the symbols used by the plots must compile.

        Returns
        -------
        usetex : bool
            True if text.usetex can be enabled, False (with a warning) otherwise.
        """

        missing = [exe for exe in ("latex", "dvipng") if shutil.which(exe) is None]
        if missing:
            print(f"WARNING AcademicStyle.apply(): {', '.join(missing)} not found, falling back to text.usetex = False")
            return False

        fig = None
        try:
            with plt.rc_context({"text.usetex": True, "font.family": AcademicStyle.FONT_FAMILY}):
                fig = plt.figure(figsize=(1, 1))
                fig.text(0.5, 0.5, r"SPL [dB re 1 $\mu$Pa$^2$/Hz], $\theta$ = 10° $|$ $1/r^2$")
                fig.savefig(io.BytesIO(), format="png")
        except Exception as e:
            print(f"WARNING AcademicStyle.apply(): LaTeX test render failed ({type(e).__name__}), "
                  f"falling back to text.usetex = False")
            return False
        finally:
            if fig is not None:
                plt.close(fig)

        return True

