"""
Module: plots/specific/sources.py
Description: Point sources (emitters FrequencySource and TimeSource): their positions and the amplitude of
             every source at every frequency (the force of a dipole, the volume velocity Q of a monopole), and
             for a TimeSource also the series in time, read from case["specific"] (Source_xyz, Source_order,
             Source_Force, Source_time, Source_Force_time). Forces [N] and volume velocities [m^3/s] never
             share an axis: with both kinds of source, the monopoles go to a second axis on the right.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import matplotlib.pyplot as plt
from wrappers.case import emitter_label, seabed_depth, point
from plot_utils.colors import get_band_colors
from plot_utils.saving import set_figure_name

MAX_LEGEND_SOURCES = 10     # [-] More sources than this (e.g. read from a source_file) get no legend entry each
MAX_BAR_TONES      = 20     # [-] More tones than this (a grid, a file) are drawn as lines instead of grouped bars
MONOPOLE           = 0      # [-] Source_order of a monopole (SOURCE_MONOPOLE of the Fortran SourceTypes)
DIPOLE             = 1      # [-] Source_order of a dipole (SOURCE_DIPOLE)

AMPLITUDE_LABEL = {DIPOLE: r"$|F|$ [N]", MONOPOLE: r"$|Q|$ [m$^3$/s]"}       # [-] Axis of the amplitudes
SERIES_LABEL    = {DIPOLE: "F [N]", MONOPOLE: r"Q [m$^3$/s]"}                # [-] Axis of the series in time
COMPONENTS      = {DIPOLE: ("Fx", "Fy", "Fz"), MONOPOLE: ("Q", "", "")}       # [-] Names of the three components
MONOPOLE_STYLE  = dict(ls='--', hatch='//')                                  # [-] How a monopole is told apart


def quantity(order: np.ndarray) -> str:
    """
    Name of what the amplitude panels show: 'Force' for dipoles only (the titles of before), 'Volume velocity'
    for monopoles only, 'Source' for a mixed set.
    """

    present = set(int(o) for o in order)
    if present == {DIPOLE}:
        return "Force"
    if present == {MONOPOLE}:
        return "Volume velocity"
    return "Source"


def order_axes(ax    : plt.Axes   = None,    # [-] Axis of the panel
               order : np.ndarray = None,    # [-] Source_order (Nsources,)
               labels: dict       = None):   # [-] Axis label of each order (AMPLITUDE_LABEL or SERIES_LABEL)
    """
    Axes on which the sources of each multipole order are drawn, so that forces [N] and volume velocities
    [m^3/s] never share an axis: with a single order, ax itself; with both, the dipoles on ax and the monopoles
    on a twin axis on the right. Each axis gets the label of its quantity.

    Returns
    -------
    axes_by_order : dict {order: Axes}
    """

    present = set(int(o) for o in order)
    if present == {MONOPOLE}:
        ax.set_ylabel(labels[MONOPOLE])
        return {MONOPOLE: ax}

    ax.set_ylabel(labels[DIPOLE])
    if MONOPOLE not in present:
        return {DIPOLE: ax}

    twin = ax.twinx()
    twin.set_ylabel(labels[MONOPOLE])
    return {DIPOLE: ax, MONOPOLE: twin}


def joint_legend(axes_by_order: dict = None,    # [-] Result of order_axes
                 handles      : list = None,    # [-] Artists in the order of the sources (lines, bar containers)
                 **kwargs):                     # [-] Passed to Axes.legend
    """
    One legend with the artists of every axis of order_axes, in the order in which the sources were drawn
    (not grouped by axis), on the last axis so that it lies on top of the others.
    """

    handles = [h for h in handles if not h.get_label().startswith('_')]
    if handles:
        list(axes_by_order.values())[-1].legend(handles, [h.get_label() for h in handles], **kwargs)


def plot_sources(case   : dict  = None,       # [-] Case of point sources
                 figsize: tuple = (12, 5)):   # [-] Figure size (a TimeSource has a third panel: 1.5 times wider)
    """
    Plots the sources of the case: their positions in the vertical plane XZ (left; dipoles as circles,
    monopoles as squares, the surface and the seabed as dashed lines) and their amplitude at every tone
    (right; one bar per source, or one line per source with more than MAX_BAR_TONES tones, on a logarithmic
    axis if they span two decades or more): |F| = sqrt(|Fx|^2 + |Fy|^2 + |Fz|^2) [N] for a dipole, |Q| [m^3/s]
    for a monopole (its only component). With both kinds of source the monopoles are drawn on a second axis on
    the right, dashed (lines) or hatched (bars).
    A TimeSource has its series in time (middle; one color per source, a solid line for its first component
    that is not zero and dashed ones for the next, the components that are always zero left out; the volume
    velocity of a monopole dashed, on the right axis if there are dipoles too) and the amplitude of its spectrum
    as lines (right; the bins of the FFT, with the Hann window).
    With more than MAX_LEGEND_SOURCES sources the legends only say how many there are.

    Returns
    -------
    fig, axes : matplotlib.figure.Figure, np.ndarray of Axes
    """

    if case["Emitter"] != "source":
        print("plot_sources(): no point sources in the case, skipping")
        return None, None

    tag    = emitter_label(case)
    xyz    = case["specific"]["Source_xyz"]                 # [m] (Nsources, 3)
    order  = case["specific"]["Source_order"]               # [-] (Nsources,) 0 monopole, 1 dipole
    force  = case["specific"]["Source_Force"]               # [N] or Q [m^3/s] (Nfreqs, Nsources, 3) complex
    Freqs  = case["Freqs"]
    labels = case["group_labels"]
    colors = get_band_colors()
    Depth  = seabed_depth(case, *point(case))

    n_src   = len(xyz)
    F_abs   = np.sqrt(np.sum(np.abs(force)**2, axis=2))     # [N] or [m^3/s] (Nfreqs, Nsources)
    is_time = "Source_Force_time" in case["specific"]
    legend  = n_src <= MAX_LEGEND_SOURCES                   # [-] One legend entry per source
    name    = quantity(order)

    if is_time:
        fig, axes = plt.subplots(1, 3, figsize=(1.5 * figsize[0], figsize[1]))
    else:
        fig, axes = plt.subplots(1, 2, figsize=figsize)

    # 1. Positions (XZ plane)
    ax = axes[0]
    ax.axhline(0.0, color='deepskyblue', linestyle='--', linewidth=1.0, alpha=0.6, label='Surface')
    ax.axhline(-Depth, color='gray', linestyle='--', linewidth=1.0, alpha=0.6, label='Seabed')
    for i in range(n_src):
        marker = 's' if order[i] == MONOPOLE else 'o'
        if legend:
            label = f"{labels[i]}: ({xyz[i, 0]:.1f}, {xyz[i, 1]:.1f}, {xyz[i, 2]:.1f}) m"
        else:
            label = f"{n_src} sources" if i == 0 else None
        ax.plot(xyz[i, 0], xyz[i, 2], marker, color=colors[i % len(colors)], ms=8, label=label)
    xpad = max(1.0, 0.1 * np.ptp(xyz[:, 0]))
    ax.set_xlim(xyz[:, 0].min() - xpad, xyz[:, 0].max() + xpad)
    ax.set_ylim(-1.05 * Depth, 0.05 * Depth)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('z [m]')
    ax.set_title(f"{tag} - Positions")
    ax.grid(True, ls='--', alpha=0.4)
    ax.legend(fontsize=10)

    if is_time:
        # 2. Series (a color per source; the first component of a dipole that is not zero is a solid line,
        #    the next ones dashed and dash-dotted; the volume velocity of a monopole dashed)
        ax       = axes[1]
        time     = case["specific"]["Source_time"]            # [s] (Nt,)
        series   = case["specific"]["Source_Force_time"]      # [N] or Q [m^3/s] (Nt, Nsources, 3)
        by_order = order_axes(ax, order, SERIES_LABEL)
        styles   = ['-', '--', '-.']
        handles  = []
        for i in range(n_src):
            o, n_drawn = int(order[i]), 0
            for k, component in enumerate(COMPONENTS[o]):
                if np.any(series[:, i, k] != 0.0):
                    style = MONOPOLE_STYLE['ls'] if o == MONOPOLE else styles[n_drawn]
                    handles += by_order[o].plot(time, series[:, i, k], style, color=colors[i % len(colors)],
                                                lw=1.0, label=f"{labels[i]} {component}" if legend else None)
                    n_drawn += 1
        ax.set_xlabel("Time [s]")
        ax.set_title(f"{tag} - {name} series")
        ax.grid(True, ls='--', alpha=0.4)
        if legend: joint_legend(by_order, handles, fontsize=10)

        # 3. Amplitude of the spectrum, one line per source
        ax       = axes[2]
        by_order = order_axes(ax, order, AMPLITUDE_LABEL)
        handles  = []
        for i in range(n_src):
            o = int(order[i])
            handles += by_order[o].plot(Freqs, F_abs[:, i], color=colors[i % len(colors)], lw=1.5, label=labels[i],
                                        ls=MONOPOLE_STYLE['ls'] if o == MONOPOLE and len(by_order) > 1 else '-')
        ax.set_xlabel("Frequency [Hz]")
        ax.set_title(f"{tag} - {name} spectrum (Hann window)")
        ax.grid(True, ls='--', alpha=0.4)
        if legend: joint_legend(by_order, handles, fontsize=10)
    elif len(Freqs) <= MAX_BAR_TONES:
        # 2. Amplitude per source and tone (grouped bars)
        ax       = axes[1]
        by_order = order_axes(ax, order, AMPLITUDE_LABEL)
        x        = np.arange(len(Freqs))
        width    = 0.8 / n_src
        handles  = []
        for i in range(n_src):
            o     = int(order[i])
            hatch = MONOPOLE_STYLE['hatch'] if o == MONOPOLE and len(by_order) > 1 else None
            handles.append(by_order[o].bar(x + (i - (n_src - 1) / 2) * width, F_abs[:, i], width,
                                           color=colors[i % len(colors)], hatch=hatch,
                                           edgecolor='k' if hatch else None, linewidth=0.5 if hatch else None,
                                           label=labels[i]))
        ax.set_xticks(x)
        ax.set_xticklabels([f"{f:g}" for f in Freqs])
        ax.set_xlabel("Frequency [Hz]")
        ax.set_title(f"{tag} - {name} amplitude")
        ax.grid(True, axis='y', ls='--', alpha=0.4)
        if legend: joint_legend(by_order, handles, fontsize=10)
    else:
        # 2. Amplitude per source at many tones (one line per source, a point per tone)
        ax       = axes[1]
        by_order = order_axes(ax, order, AMPLITUDE_LABEL)
        handles  = []
        for i in range(n_src):
            o = int(order[i])
            handles += by_order[o].plot(Freqs, F_abs[:, i], '.', color=colors[i % len(colors)], lw=1.0, ms=3,
                                        ls=MONOPOLE_STYLE['ls'] if o == MONOPOLE and len(by_order) > 1 else '-',
                                        label=labels[i] if legend else None)
        if Freqs.max() / Freqs.min() >= 100.0:
            ax.set_xscale('log')
        ax.set_xlabel("Frequency [Hz]")
        ax.set_title(f"{tag} - {name} amplitude ({len(Freqs)} tones)")
        ax.grid(True, ls='--', alpha=0.4)
        if legend: joint_legend(by_order, handles, fontsize=10)

    fig.tight_layout()
    set_figure_name(fig, "sources")

    return fig, axes
