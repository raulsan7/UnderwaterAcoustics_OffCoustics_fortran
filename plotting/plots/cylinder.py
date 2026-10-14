"""
Module: plots/cylinder.py
Description: Unwrapped map of the 'cylinder' step (azimuth x depth), for any emitter.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import matplotlib.pyplot as plt
import utils.AcousticUtils as au
from wrappers.case import emitter_label
from plots.common import no_absorption, psd_not_applicable
from plots.slices import get_contour_levels
from plot_utils.conversions import convert_to, is_valid
from plot_utils.saving import set_figure_name


def plot_cylinder(case            : dict = None,        # [-] Case
                  mode            : str = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                  target_frequency: float = 0.88,       # [Hz] Target frequency for single-frequency modes
                  absorption      : bool = True,        # [-] Whether to apply absorption attenuation
                  filter_under    : float = None,       # [Hz] Lower frequency cutoff (f >= filter_under)
                  filter_over     : float = None,       # [Hz] Upper frequency cutoff (f <= filter_over)
                  figsize         : tuple = (12,5),     # [-] Figure size
                  cmap            : str = 'inferno',    # [-] Contour colormap
                  ax_in           : plt.Axes = None):   # [-] Optional external axes for subplots
    """
    Plots 2D unwrapped cylindrical surface maps for acoustic magnitudes (OASPL, SPL, etc.). Observers without
    data (NaN, inside the structure) are left as holes in the contours. Mode 'SWL' is not implemented: the
    figure is left empty with a message.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    if not case["has_cylinder"]:
        print("plot_cylinder(): no cylinder data, skipping")
        return None, None

    absorption = no_absorption(case, absorption, "plot_cylinder")
    if psd_not_applicable([case], mode, "plot_cylinder"):
        return None, None

    tag = emitter_label(case)

    # Extract data
    theta     = case["Theta_deg_cylinder"]   # 1D array of angles [deg]
    z         = case["Z_cylinder"]           # 1D array of height/depth [m]
    p         = case["P_cylinder"]           # Pressure grid
    Freqs     = case["Freqs"]
    df        = case["bin_width"]            # Frequency bin width [Hz] for 'PSD'
    p_ref     = case["p_ref"]
    r         = case["R_cylinder"]

    if not is_valid(p, "plot_cylinder"):
            return None, None

    nf, ntheta, nz = len(Freqs), len(theta), len(z)

    # Reshape
    p         = p.reshape((nf, ntheta*nz))

    # Distance grid for absorption
    distances = np.full(ntheta*nz, r, dtype=float)

    # Filter frequency
    freq_mask = np.ones(nf, dtype=bool)
    if filter_under is not None:
        freq_mask &= (Freqs >= filter_under)
    if filter_over is not None:
        freq_mask &= (Freqs <= filter_over)
    if not np.any(freq_mask):
        print("plot_cylinder(): no frequencies match the specified filter bounds!")
        return None, None

    p, Freqs = p[freq_mask], Freqs[freq_mask]
    nf = len(Freqs)

    mode = mode.upper()

    if ax_in is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = ax_in
        fig = ax.figure

    title_extra = ""
    # Select mode
    if mode == 'OASPL':
        vals = au.pressure_to_OASPL(p, p_ref, absorption, Freqs, distances, enbw=case["ENBW"], water=case["water"])
        unit_label = r"OASPL [dB re 1$\mu$ Pa]"

    elif mode == 'SWL':
        print("plot_cylinder(): SWL not yet implemented")
        return fig, ax
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f = Freqs[idx]
        p = p[idx,:]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = convert_to(p, mode, df, enbw=case["ENBW"])

        if absorption and mode in ('SPL', 'PSD'): vals = au.add_absorption(f, vals, distances, water=case["water"])

    # Plot
    theta_plot = np.append(theta, 360.)
    vals = vals.reshape((ntheta, nz))
    vals_plot = np.append(vals, vals[0:1, :], axis=0)

    # Upper half of the range
    custom_levels = get_contour_levels(vals_plot)
    line_levels = custom_levels[::2]

    # Matplotlib contour expects Z shape (len(y), len(x)), so transpose vals_plot. NaN: holes
    contour_set = ax.contourf(theta_plot, z, vals_plot.T, levels=custom_levels, cmap=cmap, extend='min')
    ax.contour(theta_plot, z, vals_plot.T, levels=line_levels, colors='k', linewidths=0.4, alpha=0.5)
    cbar = fig.colorbar(contour_set, ax=ax, orientation='vertical', pad=0.02)
    cbar.set_label(unit_label)

    ax.set_xlabel(r'Azimuth Angle, $\theta$ [deg]')
    ax.set_ylabel('Depth, z [m]')
    ax.set_title(f"{tag} - {mode} {title_extra}")
    ax.set_xlim(theta_plot.min(), theta_plot.max())
    ax.set_ylim(z.min(), z.max())
    ax.set_xticks(np.arange(0, 361, 45))

    if ax_in is None:
        fig.tight_layout()
        set_figure_name(fig, "cylinder", mode, title_extra)

    return fig, ax
