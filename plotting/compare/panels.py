"""
Module: compare/panels.py
Description: Panels, titles and legends of a comparison of two or more cases of any emitter (one panel per
             case, drawn by the functions of plots/ through ax_in).

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import re
import numpy as np
import matplotlib.pyplot as plt
from wrappers.case import emitter_label

MAX_COLS = 3    # [-] Panels per row: more cases wrap into new rows


def select_cases(cases: list = None,     # [-] List of Cases
                 flag : str  = None,     # [-] has_<step> flag the plot needs (or the name of what it needs)
                 name : str  = None,     # [-] Name of the calling function, for the messages
                 has         = None):    # [-] Function (case) -> bool telling whether a case has it (default: the flag)
    """
    Keeps the cases whose file contains the step (has_<step> flag True, or has(case) True). The others (older
    versions of the Fortran code or interrupted runs) are left out with a message.

    Returns
    -------
    cases : list of dict
    """

    if not cases:
        print(f"{name}(): no cases given, skipping")
        return []

    if has is None:
        has = lambda case: bool(case.get(flag, False))

    for case in cases:
        if not has(case):
            print(f"{name}(): {case['Name']} has no {flag.replace('has_', '')} data, left out of the comparison")

    return [case for case in cases if has(case)]

def make_panels(n      : int   = 1,         # [-] Number of cases
                figsize: tuple = (7, 5),    # [-] Size of a single panel
                stacked: bool  = False,     # [-] One column (True) or rows of up to MAX_COLS panels (False)
                **kwargs):                  # [-] Extra arguments of plt.subplots (sharey, subplot_kw, ...)
    """
    Creates one panel per case: a single row for up to MAX_COLS cases (more cases wrap into new rows), or a
    single column if stacked.

    Returns
    -------
    fig, axes, ncols : matplotlib.figure.Figure, np.ndarray of Axes (exactly n), int
    """

    if stacked:
        nrows, ncols = n, 1
    else:
        ncols = min(n, MAX_COLS)
        nrows = int(np.ceil(n / ncols))

    total_figsize = (figsize[0] * ncols, figsize[1] * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=total_figsize, squeeze=False, **kwargs)
    axes = axes.ravel()

    # Spare cells of the last row
    for ax in axes[n:]:
        fig.delaxes(ax)

    return fig, axes[:n], ncols

def remove_common_title_part(axes    : np.ndarray = None,     # [-] Axes of the comparison
                             pattern : str        = None):    # [-] Regex of the part to remove
    """
    Deletes from every title the part matching 'pattern' (e.g. ' - Depth = -15.0 [m]') if it is the same in
    all of them.
    """

    titles = [ax.get_title() for ax in axes]

    match = re.search(pattern, titles[0], re.IGNORECASE)

    if match:
        common_str = match.group(1)

        if all(common_str in t for t in titles):
            for ax in axes:
                new_title = ax.get_title().replace(common_str, '').strip()

                if new_title.endswith('-') or new_title.endswith(','):
                    new_title = new_title[:-1].strip()

                ax.set_title(new_title)

def wrap_case_titles(axes : np.ndarray = None,     # [-] Axes of the comparison
                     cases: list       = None):    # [-] Case plotted in each panel
    """
    Moves the case label of every panel title to its own line ('label - OASPL - ...' -> 'label\\nOASPL - ...'):
    comparison labels are long and the titles of neighbouring panels would overlap.
    """

    for ax, case in zip(axes, cases):
        title = ax.get_title()
        label = emitter_label(case)
        if title.startswith(label + " - "):
            ax.set_title(label + "\n" + title[len(label) + 3:])

def add_case_legend(fig,                    # [-] Figure of the comparison
                    handles: list = None,   # [-] One handle per case
                    labels : list = None,   # [-] One label per case
                    title  : str  = None):  # [-] Legend title
    """
    Adds one legend entry per case below the axes, one per row, for figures that overlay every case (their
    labels are long).

    Returns
    -------
    bottom : float
        Fraction of the figure height the legend takes, for fig.tight_layout(rect=[0, bottom, 1, 1]).
    """

    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=1, fontsize=11,
               title=title, title_fontsize=11, frameon=False)

    rows = len(labels) + (1 if title else 0)
    return min(0.5, (0.27 * rows + 0.15) / fig.get_figheight())
