"""
Module: plots/line.py
Description: Magnitude along the straight line of the 'line' step (P1 -> P2), for any emitter, with one curve
             per member of a family (--ws).

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import matplotlib.pyplot as plt
from wrappers.case import emitter_label, case_color, family_tag, get_family
from plots.common import no_absorption, psd_not_applicable
from plots.decay import get_family_curves
from plot_utils.drawing import add_family_legends
from plot_utils.saving import set_figure_name


def plot_line(case            : dict  = None,         # [-] Case
              mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
              target_frequency: float = 0.88,         # [Hz] Target frequency for single-frequency modes
              absorption      : bool  = True,         # [-] Whether to apply absorption attenuation
              filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
              filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
              figsize         : tuple = (10,5),       # [-] Figure size
              ax_in           : plt.Axes = None):     # [-] Optional external axes for subplots
    """
    Plots acoustic magnitude along an arbitrary straight line between two points. Points without data (NaN,
    inside the structure) are gaps in the curve. With a family (--ws) one curve per member (color = member,
    legend below the axes), not filled.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    if not case["has_line"]:
        print("plot_line(): no line data, skipping")
        return None, None

    absorption = no_absorption(case, absorption, "plot_line")
    if psd_not_applicable(get_family(case)[0], mode, "plot_line"):
        return None, None

    tag    = emitter_label(case)
    p1     = case["P1_line"]
    p2     = case["P2_line"]
    do_log = case["Logspace_line"]

    curves, ws_mode = get_family_curves(case, "line", "plot_line", mode=mode, target_frequency=target_frequency,
                                        absorption=absorption, filter_under=filter_under, filter_over=filter_over)
    if not curves:
        return None, None

    if ax_in is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = ax_in
        fig = ax.figure

    draw = ax.semilogx if do_log else ax.plot

    # Plot
    handles, labels = [], []
    for c, (r, vals, unit_label, title_extra) in curves:
        if ws_mode:
            color = c["family_color"]
            handles.append(draw(r, vals, lw=c["family_lw"], color=color, zorder=c["family_zorder"])[0])
            labels.append(c["family_label"])
        else:
            color = case_color(case)
            draw(r, vals, lw=2, color=color)

        if not ws_mode:
            ax.fill_between(r, vals, 0, alpha=0.1, color=color)

    ax.grid(True, which='both', ls='--', alpha=0.7)

    # Add text box with P1 and P2 coordinates for context
    coord_text = f"P1: {np.round(p1, 1)} m\nP2: {np.round(p2, 1)} m"
    ax.text(0.97, 0.95, coord_text, transform=ax.transAxes, ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'),
            fontsize=9)

    ax.set_xlabel("Distance from P1 [m]")
    ax.set_ylabel(unit_label)
    ax.set_title(f"{tag} - {mode.upper()} {title_extra}")
    ax.set_xlim(min(r_c.min() for _, (r_c, _, _, _) in curves), max(r_c.max() for _, (r_c, _, _, _) in curves))
    ax.set_ylim(min(np.nanmin(vals_c) for _, (_, vals_c, _, _) in curves),
                max(np.nanmax(vals_c) for _, (_, vals_c, _, _) in curves))

    if ws_mode:
        add_family_legends(ax, handles, labels, family_title=case["family_title"])

    if ax_in is None:
        fig.tight_layout()
        set_figure_name(fig, "line", mode.upper(), title_extra, family_tag([case]) if ws_mode else "")

    return fig, ax
