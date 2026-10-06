"""
Module: plot_utils/saving.py
Description: Names of the figures and their saving with --save (one subfolder per case or comparison).

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import re
import argparse
import matplotlib.pyplot as plt
from plot_utils.files import rel


def set_figure_name(fig,        # [-] Figure created by a plot function
                    *parts):    # [-] Words of the name, e.g. "spectrum", "PSD", "octave"
    """
    Names a figure after the plot that created it (e.g. 'spectrum_PSD_octave', 'cylinder_SPL_0.88Hz').
    save_figures uses the name as file name.
    """

    words = []
    for part in parts:
        word = str(part).strip() if part is not None else ""
        word = re.sub(r"^at\s+", "", word)          # 'at 0.88 Hz' -> '0.88 Hz'
        word = re.sub(r"\s+", "", word)             # '0.88 Hz'    -> '0.88Hz'
        word = re.sub(r"[^\w.\-]", "", word)
        if word:
            words.append(word)

    fig.set_label("_".join(words))

def save_figures(args : argparse.Namespace = None,     # [-] Parsed command-line arguments
                 cases: list = None):                  # [-] Plotted case (dict) or cases (list of dict)
    """
    With --save, saves every open figure as <output_dir>/<case name>/<figure name>.<format>. A comparison of
    several cases goes to <output_dir>/comparison_<name 1>_vs_<name 2>...
    """

    if not args.save:
        return

    if isinstance(cases, dict):
        cases = [cases]

    names  = [case["Name"] for case in cases]
    folder = names[0] if len(names) == 1 else "comparison_" + "_vs_".join(names)
    if len(folder) > 200:
        folder = f"comparison_{len(names)}_cases_{names[0]}"

    out_dir = args.output_dir / folder
    out_dir.mkdir(parents=True, exist_ok=True)

    used = {}
    for num in plt.get_fignums():
        fig  = plt.figure(num)
        name = fig.get_label() or f"figure_{num}"

        # Repeated plots (same function and mode) get a counter
        used[name] = used.get(name, 0) + 1
        if used[name] > 1:
            name = f"{name}_{used[name]}"

        for fmt in args.format:
            fig.savefig(out_dir / f"{name}.{fmt}", format=fmt)

    print(f"\n{len(plt.get_fignums())} figures saved in {rel(out_dir)} ({', '.join(args.format)})")
