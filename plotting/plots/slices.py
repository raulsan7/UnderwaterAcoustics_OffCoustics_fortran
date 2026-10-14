"""
Module: plots/slices.py
Description: Maps of the XY, XZ and vertical slices, for any emitter, with the nodes of the case drawn over them
             as points (monopoles as squares).

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import matplotlib.pyplot as plt
import utils.AcousticUtils as au
from wrappers.case import emitter_label
from plots.common import no_absorption, psd_not_applicable
from plot_utils.conversions import convert_to, is_valid
from plot_utils.saving import set_figure_name


# ---------- Helpers ---------- #
def get_contour_levels(vals,                # [-] Magnitude of the map (NaN: holes)
                       n   : int = 60):     # [-] Number of levels
    """
    Levels of a map that shows its upper half: from half of the maximum up to it (contourf with extend='min'
    paints what lies below). Without a positive maximum (non-dB modes) the levels span the whole range.

    Returns
    -------
    levels : np.ndarray
    """

    max_val = np.nanmax(vals)
    min_val = max_val / 2.0
    if not min_val < max_val: min_val = np.nanmin(vals)     # Non-dB modes may have a maximum <= 0

    return np.linspace(min_val, max_val, n)

def draw_nodes(ax,                          # [-] Axes of the map
               case: dict       = None,     # [-] Case whose nodes are drawn
               u   : np.ndarray = None,     # [m] Horizontal coordinate of every node on the map
               v   : np.ndarray = None):    # [m] Vertical coordinate of every node on the map
    """
    Draws the nodes of the case over a map: black points, and squares for the monopoles when the case gives
    the multipole order of its nodes (node_order 0). The artists carry gid 'nodes'.
    """

    order    = case.get("node_order")
    monopole = np.zeros(len(u), dtype=bool) if order is None else (order == 0)

    if np.any(~monopole):
        ax.plot(u[~monopole], v[~monopole], 'o', color='black', markersize=2, zorder=10, gid="nodes")
    if np.any(monopole):
        ax.plot(u[monopole], v[monopole], 's', color='black', markersize=3, zorder=10, gid="nodes")


# ---------- Slice Plots ---------- #
def plot_sliceXY(case            : dict  = None,         # [-] Case
                 mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                 target_frequency: float = 0.88,         # [Hz] Target frequency for single-frequency modes
                 absorption      : bool  = True,         # [-] Whether to apply absorption attenuation
                 filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                 filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                 structure       : bool = True,          # [-] Wheter to plot the nodes on the slice
                 figsize         : tuple = (7,7),        # [-] Figure size
                 cmap            : str   = 'inferno',    # [-] Contour colormap
                 ax_in           : plt.Axes = None):     # [-] Optional external axes for subplots
    """
    Plots a 2D spatial slice in the XY plane of acoustic magnitudes. Observers without data (NaN, inside the
    structure) are left as holes in the contours. Mode 'SWL' is not implemented: the figure is left empty
    with a message.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    if not case["has_slicexy"]:
        print("plot_sliceXY(): no XY slice data, skipping")
        return None, None

    absorption = no_absorption(case, absorption, "plot_sliceXY")
    if psd_not_applicable([case], mode, "plot_sliceXY"):
        return None, None

    tag = emitter_label(case)

    # Extract data
    x         = case["X_slicexy"]           # 1D array of X coordinates [m]
    y         = case["Y_slicexy"]           # 1D array of Y coordinates [m]
    z_slice   = case["Z_slicexy"]           # Z level of the slice [m]
    p         = case["P_slicexy"]           # Pressure grid
    Freqs     = case["Freqs"]
    df        = case["bin_width"]           # Frequency bin width [Hz] for 'PSD'
    observers = case["Obs_slicexy"]
    p_ref     = case["p_ref"]
    cx, cy    = case["Center_slicexy"]

    if not is_valid(p, "plot_sliceXY"):
        return None, None

    Nfreqs = len(Freqs)
    distances = np.sqrt((observers[..., 0] - cx)**2 + (observers[..., 1] - cy)**2)

    # Filter frequencies
    freq_mask = np.ones(Nfreqs, dtype=bool)
    if filter_under is not None:
        freq_mask &= (Freqs >= filter_under)
    if filter_over is not None:
        freq_mask &= (Freqs <= filter_over)
    if not np.any(freq_mask):
        print("plot_sliceXY(): no frequencies match the specified filter bounds!")
        return None, None

    p, Freqs = p[freq_mask], Freqs[freq_mask]

    mode = mode.upper()

    if ax_in is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = ax_in
        fig = ax.figure

    title_extra = ""
    X, Y = np.meshgrid(x, y, indexing='ij')

    # Select mode
    if mode == 'OASPL':
        vals = au.pressure_to_OASPL(p, p_ref, absorption, Freqs, distances, enbw=case["ENBW"], water=case["water"])
        unit_label = r"OASPL [dB re 1$\mu$ Pa]"
    elif mode == 'SWL':
        print("WARNING plot_sliceXY(): SWL not yet implemented")
        return fig, ax
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f   = Freqs[idx]
        p_f = p[idx, :]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = convert_to(p_f, mode, df, enbw=case["ENBW"])
        if absorption and mode in ('SPL', 'PSD'): vals = au.add_absorption(f, vals, distances, water=case["water"])

    # Plot (NaN: holes)
    contour_set = ax.contourf(X, Y, vals, levels=60, cmap=cmap)
    cbar = fig.colorbar(contour_set, ax=ax, orientation='horizontal', pad = 0.1)
    cbar.set_label(unit_label)
    if abs(x[-1]-x[0]) == abs(y[-1]-y[0]): ax.set_aspect('equal')

    # Nodes of the case
    if structure:
        nodes = case["nodes"]
        draw_nodes(ax, case, nodes[:, 0], nodes[:, 1])

    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title(f"{tag} - {mode} {title_extra} - Depth = {z_slice:.1f} [m]")
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())

    if ax_in is None:
        fig.tight_layout()
        set_figure_name(fig, "sliceXY", mode, title_extra)

    return fig, ax

def plot_sliceXZ(case            : dict  = None,         # [-] Case
                 mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                 target_frequency: float = 0.88,         # [Hz] Target frequency for single-frequency modes
                 absorption      : bool  = True,         # [-] Whether to apply absorption attenuation
                 filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                 filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                 structure       : bool  = True,         # [-] Wheter to plot the nodes on the slice
                 figsize         : tuple = (7,7),        # [-] Figure size
                 cmap            : str   = 'inferno',    # [-] Contour colormap
                 ax_in           : plt.Axes = None):     # [-] Optional external axes for subplots
    """
    Plots a 2D spatial slice in the XZ vertical plane of acoustic magnitudes. Observers without data (NaN,
    inside the structure) are left as holes in the contours. With the nodes drawn, the vertical axis reaches
    the highest of them. Mode 'SWL' is not implemented: the figure is left empty with a message.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    if not case["has_slicexz"]:
        print("plot_sliceXZ(): no XZ slice data, skipping")
        return None, None

    absorption = no_absorption(case, absorption, "plot_sliceXZ")
    if psd_not_applicable([case], mode, "plot_sliceXZ"):
        return None, None

    tag = emitter_label(case)

    # Extract data
    x         = case["X_slicexz"]           # 1D array of X coordinates [m]
    z         = case["Z_slicexz"]           # 1D array of Z coordinates [m]
    y_slice   = case["Y_slicexz"]           # Y level of the slice [m]
    p         = case["P_slicexz"]           # Pressure grid (Nf, nx, nz)
    Freqs     = case["Freqs"]
    df        = case["bin_width"]           # Frequency bin width [Hz] for 'PSD'
    observers = case["Obs_slicexz"]         # Observer grid (nx, nz, 3)
    p_ref     = case["p_ref"]

    if not is_valid(p, "plot_sliceXZ"):
        return None, None

    cx, cz = np.mean(x), np.mean(z)
    Nfreqs = len(Freqs)
    distances = np.sqrt((observers[...,0]-cx)**2 + (observers[...,2]-cz)**2)

    # Filter frequencies
    freq_mask = np.ones(Nfreqs, dtype=bool)
    if filter_under is not None:
        freq_mask &= (Freqs >= filter_under)
    if filter_over is not None:
        freq_mask &= (Freqs <= filter_over)
    if not np.any(freq_mask):
        print("plot_sliceXZ(): no frequencies match the specified filter bounds!")
        return None, None

    p, Freqs = p[freq_mask], Freqs[freq_mask]

    mode = mode.upper()

    if ax_in is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = ax_in
        fig = ax.figure

    title_extra = ""
    X, Z = np.meshgrid(x, z, indexing='ij')

    # Select mode
    if mode == 'OASPL':
        vals = au.pressure_to_OASPL(p, p_ref, absorption, Freqs, distances, enbw=case["ENBW"], water=case["water"])
        unit_label = r"OASPL [dB re 1$\mu$ Pa]"
    elif mode == 'SWL':
        print("WARNING plot_sliceXZ(): SWL not yet implemented")
        return fig, ax
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f   = Freqs[idx]
        p_f = p[idx, :, :]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = convert_to(p_f, mode, df, enbw=case["ENBW"])
        if absorption and mode in ('SPL', 'PSD'): vals = au.add_absorption(f, vals, distances, water=case["water"])

    # Plot (upper half of the range, NaN: holes)
    contour_set = ax.contourf(X, Z, vals, levels=get_contour_levels(vals), cmap=cmap, extend='min')
    cbar = fig.colorbar(contour_set, ax=ax, orientation='horizontal', pad=0.1)
    cbar.set_label(unit_label)
    if abs(x[-1]-x[0]) == abs(z[-1]-z[0]): ax.set_aspect('equal')

    # Nodes of the case: the vertical axis reaches the highest one
    z_top = z.max()
    if structure:
        nodes = case["nodes"]
        draw_nodes(ax, case, nodes[:, 0], nodes[:, 2])
        if len(nodes):
            z_top = max(z_top, np.nanmax(nodes[:, 2]))

    ax.set_xlabel('x [m]')
    ax.set_ylabel('z [m]')
    ax.set_title(f"{tag} - {mode} {title_extra} - Y = {y_slice:.1f} [m]")
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(z.min(), z_top)

    if ax_in is None:
        fig.tight_layout()
        set_figure_name(fig, "sliceXZ", mode, title_extra)

    return fig, ax

def plot_sliceVertical(case            : dict  = None,         # [-] Case
                       mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                       target_frequency: float = 0.88,         # [Hz] Target frequency for single-frequency modes
                       absorption      : bool  = True,         # [-] Whether to apply absorption attenuation
                       filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                       filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                       structure       : bool  = True,         # [-] Wheter to plot the nodes on the slice
                       figsize         : tuple = (7,7),        # [-] Figure size
                       cmap            : str   = 'inferno',    # [-] Contour colormap
                       ax_in           : plt.Axes = None):     # [-] Optional external axes for subplots
    """
    Plots a 2D spatial slice in a vertical plane of arbitrary azimuth of acoustic magnitudes. The plane normal
    is (cos az, sin az) and its horizontal axis u runs along (-sin az, cos az) from Center_sliceV
    (AcousticSolver.run_sliceVertical). The plane usually crosses the structure: observers without data (NaN,
    inside the structure) are left as holes in the contours. The nodes are projected onto the plane. Mode
    'SWL' is not implemented: the figure is left empty with a message.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    if not case["has_slicevertical"]:
        print("plot_sliceVertical(): no vertical slice data, skipping")
        return None, None

    absorption = no_absorption(case, absorption, "plot_sliceVertical")
    if psd_not_applicable([case], mode, "plot_sliceVertical"):
        return None, None

    tag = emitter_label(case)

    # Extract data
    u         = case["U_sliceV"]            # 1D array of along-slice coordinates [m]
    z         = case["Z_sliceV"]            # 1D array of Z coordinates [m]
    p         = case["P_sliceV"]            # Pressure grid (Nf, nu, nz)
    observers = case["Coords_sliceV"]       # Observer grid (nu, nz, 3)
    azimuth   = case["Azimuth_sliceV"]      # Azimuth angle of the plane normal [deg]
    center    = case["Center_sliceV"]       # Center point [cx, cy]
    Freqs     = case["Freqs"]
    df        = case["bin_width"]           # Frequency bin width [Hz] for 'PSD'
    p_ref     = case["p_ref"]

    if not is_valid(p, "plot_sliceVertical"):
        return None, None

    cx, cy = center[0], center[1]
    Nfreqs = len(Freqs)
    distances = np.sqrt((observers[...,0]-cx)**2 + (observers[...,1]-cy)**2)     # Horizontal distance to the centre

    # Filter frequencies
    freq_mask = np.ones(Nfreqs, dtype=bool)
    if filter_under is not None:
        freq_mask &= (Freqs >= filter_under)
    if filter_over is not None:
        freq_mask &= (Freqs <= filter_over)
    if not np.any(freq_mask):
        print("plot_sliceVertical(): no frequencies match the specified filter bounds!")
        return None, None

    p, Freqs = p[freq_mask], Freqs[freq_mask]

    mode = mode.upper()

    if ax_in is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = ax_in
        fig = ax.figure

    title_extra = ""
    U, Z = np.meshgrid(u, z, indexing='ij')

    # Select mode
    if mode == 'OASPL':
        vals = au.pressure_to_OASPL(p, p_ref, absorption, Freqs, distances, enbw=case["ENBW"], water=case["water"])
        unit_label = r"OASPL [dB re 1$\mu$ Pa]"
    elif mode == 'SWL':
        print("WARNING plot_sliceVertical(): SWL not yet implemented")
        return fig, ax
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f   = Freqs[idx]
        p_f = p[idx, :, :]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = convert_to(p_f, mode, df, enbw=case["ENBW"])
        if absorption and mode in ('SPL', 'PSD'): vals = au.add_absorption(f, vals, distances, water=case["water"])

    # Plot (upper half of the range, NaN: holes)
    contour_set = ax.contourf(U, Z, vals, levels=get_contour_levels(vals), cmap=cmap, extend='min')
    cbar = fig.colorbar(contour_set, ax=ax, orientation='horizontal', pad=0.1)
    cbar.set_label(unit_label)
    if abs(u[-1]-u[0]) == abs(z[-1]-z[0]): ax.set_aspect('equal')

    # Nodes projected onto the plane: u = (x - cx, y - cy) . (-sin az, cos az). The vertical axis reaches the highest one
    z_top = z.max()
    if structure:
        theta_rad = np.deg2rad(azimuth)
        cos_t = np.cos(theta_rad)
        sin_t = np.sin(theta_rad)

        nodes  = case["nodes"]
        u_proj = -(nodes[:,0] - cx) * sin_t + (nodes[:,1] - cy) * cos_t
        draw_nodes(ax, case, u_proj, nodes[:, 2])
        if len(nodes):
            z_top = max(z_top, np.nanmax(nodes[:, 2]))

    ax.set_xlabel('u (Distance along slice) [m]')
    ax.set_ylabel('z [m]')
    ax.set_title(f"{tag} - {mode} {title_extra} - Azimuth = {azimuth:.1f} [deg]")
    ax.set_xlim(u.min(), u.max())
    ax.set_ylim(z.min(), z_top)

    if ax_in is None:
        fig.tight_layout()
        set_figure_name(fig, "sliceVertical", mode, title_extra)

    return fig, ax
