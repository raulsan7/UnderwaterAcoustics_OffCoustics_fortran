"""
Module: plot_utils/colors.py
Description: Colors and line styles of the plots: color of a case type, colors of the cases of a comparison,
             shades of a color and the palettes of AcademicStyle.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import matplotlib.colors as mcolors
from plot_utils.style import AcademicStyle


def get_type_color(tag: str = ""):     # [-] Case type (e.g. 'DTU10MWMonopile') or label containing it
    """
    Color of a case type: monopiles, shallow floating platforms and the rest (CASE_COLORS).

    Returns
    -------
    color : str
    """

    if "Monopile" in tag:
        return AcademicStyle.CASE_COLORS[0]
    if "shallow" in tag:
        return AcademicStyle.CASE_COLORS[1]

    return AcademicStyle.CASE_COLORS[2]

def get_case_colors(n: int = 1):    # [-] Number of cases
    """
    Colors of the cases of a comparison, by position: CASE_COLORS first, then BAND_COLORS.

    Returns
    -------
    colors : list of str
    """

    palette = AcademicStyle.CASE_COLORS + [c for c in AcademicStyle.BAND_COLORS if c not in AcademicStyle.CASE_COLORS]

    return [palette[i % len(palette)] for i in range(n)]

def get_color_shades(color,             # [-] Base color (any matplotlib color)
                     n: int = 1):       # [-] Number of shades
    """
    Shades of a single color, from darker to lighter, to tell apart several curves of the same case (e.g.
    the observers of a spectrum). The base color is mixed with black (darker shades) or with white (lighter
    tints), so every curve still reads as the same color. A single shade is the base color itself.

    Returns
    -------
    shades : list of tuple
        n RGB colors.
    """

    rgb = np.array(mcolors.to_rgb(color))
    if n <= 1:
        return [tuple(rgb)]

    shades = []
    for t in np.linspace(-0.35, 0.45, n):     # [-] < 0: fraction of black, > 0: fraction of white
        target = np.zeros(3) if t < 0 else np.ones(3)
        shades.append(tuple(rgb + abs(t) * (target - rgb)))

    return shades

def get_line_styles():
    return AcademicStyle.LINE_STYLES

def get_band_colors():
    return AcademicStyle.BAND_COLORS
