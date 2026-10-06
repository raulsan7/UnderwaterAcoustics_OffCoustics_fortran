"""
Module: plots/specific/sources.py
Description: Frequency-domain point sources (emitter FrequencySource): their positions and the amplitude of the
             force of every source at every tone, read from case["specific"] (Source_xyz, Source_Force).

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import matplotlib.pyplot as plt
from wrappers.case import emitter_label, seabed_depth, point
from plot_utils.colors import get_band_colors
from plot_utils.saving import set_figure_name


def plot_sources(case   : dict  = None,       # [-] Case of frequency-domain sources
                 figsize: tuple = (12, 5)):   # [-] Figure size
    """
    Plots the sources of the case: their positions in the vertical plane XZ (left; dipoles as circles,
    monopoles as squares, the surface and the seabed as dashed lines) and the amplitude of their force,
    |F| = sqrt(|Fx|^2 + |Fy|^2 + |Fz|^2), at every tone (right; one bar per source).

    Returns
    -------
    fig, axes : matplotlib.figure.Figure, np.ndarray of Axes
    """

    if case["Emitter"] != "source":
        print("plot_sources(): no frequency-domain sources in the case, skipping")
        return None, None

    tag    = emitter_label(case)
    xyz    = case["specific"]["Source_xyz"]                 # [m] (Nsources, 3)
    order  = case["specific"]["Source_order"]               # [-] (Nsources,) 0 monopole, 1 dipole
    force  = case["specific"]["Source_Force"]               # [N] (Nfreqs, Nsources, 3) complex
    Freqs  = case["Freqs"]
    labels = case["group_labels"]
    colors = get_band_colors()
    Depth  = seabed_depth(case, *point(case))

    n_src  = len(xyz)
    F_abs  = np.sqrt(np.sum(np.abs(force)**2, axis=2))      # [N] (Nfreqs, Nsources)

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # 1. Positions (XZ plane)
    ax = axes[0]
    ax.axhline(0.0, color='deepskyblue', linestyle='--', linewidth=1.0, alpha=0.6, label='Surface')
    ax.axhline(-Depth, color='gray', linestyle='--', linewidth=1.0, alpha=0.6, label='Seabed')
    for i in range(n_src):
        marker = 's' if order[i] == 0 else 'o'
        ax.plot(xyz[i, 0], xyz[i, 2], marker, color=colors[i % len(colors)], ms=8,
                label=f"{labels[i]}: ({xyz[i, 0]:.1f}, {xyz[i, 1]:.1f}, {xyz[i, 2]:.1f}) m")
    xpad = max(1.0, 0.1 * np.ptp(xyz[:, 0]))
    ax.set_xlim(xyz[:, 0].min() - xpad, xyz[:, 0].max() + xpad)
    ax.set_ylim(-1.05 * Depth, 0.05 * Depth)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('z [m]')
    ax.set_title(f"{tag} - Positions")
    ax.grid(True, ls='--', alpha=0.4)
    ax.legend(fontsize=10)

    # 2. Force amplitude per source and tone (grouped bars)
    ax    = axes[1]
    x     = np.arange(len(Freqs))
    width = 0.8 / n_src
    for i in range(n_src):
        ax.bar(x + (i - (n_src - 1) / 2) * width, F_abs[:, i], width, color=colors[i % len(colors)], label=labels[i])
    ax.set_xticks(x)
    ax.set_xticklabels([f"{f:g}" for f in Freqs])
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel(r"$|F|$ [N]")
    ax.set_title(f"{tag} - Force amplitude")
    ax.grid(True, axis='y', ls='--', alpha=0.4)
    ax.legend(fontsize=10)

    fig.tight_layout()
    set_figure_name(fig, "sources")

    return fig, axes
