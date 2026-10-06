"""
Module: plot_utils/tex.py
Description: Text shown in figures when LaTeX renders it (text.usetex).

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import re
import matplotlib.pyplot as plt


def tex_escape(text: str = ""):     # [-] Text to show in a figure
    """
    Escapes the characters LaTeX treats as special (e.g. the '_' of file names) when text.usetex is on.

    Returns
    -------
    text : str
    """

    if not plt.rcParams["text.usetex"]:
        return text

    return re.sub(r"([_&%#$])", r"\\\1", text)
