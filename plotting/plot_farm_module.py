"""
Module: plot_farm_module.py
Description: Plots of a wind farm (.hdf5 files of FORTRAN/farm_acoustic_data): structure of every turbine,
             spectrums, polar directivity, cylinder, line and XY/XZ/vertical slices. Observers without data
             (NaN pressure, inside a structure) are shown as holes.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 09/2026
"""

import numpy as np
import matplotlib.pyplot as plt
import plot_utils as pu
import utils.AcousticUtils as au
import utils.EnvironmentalData as ed



# ---------- Structural Plots ---------- #
def plot_structure(case   : dict = None,        # [-] Dictionary containing simulation data
                   mode   : str = '3D',         # [-] Plot mode for the structure representation
                   WindDir: bool = True,        # [-] Wheter to plot or not wind direction
                   figsize: tuple = (8,6),      # [-] Figure size
                   ax_in  : plt.Axes = None):   # [-] optional external axes for subplots
    """
    Plots the structure nodes in shape(Nmembers, Nnodes, 3) per turbine in the farm

    Returns
    -------
    fig, ax: matplotlib.figure.Figure, matplotlib.axes.Axes
    The generated figure and axes objet for further editing.
    """

    if not all("Structure_nodes" in turb for turb in case["Turbine_parameters"]):
        print("plot_structure(): no structure data, skipping")
        return None, None

    if mode not in ("3D", "xy", "xz"):
        raise ValueError("plot_structure(): mode must be '3D', 'xy', 'xz'")

    # Extract nodes and depth. Turbines of different types have different (Nmembers, Nnodes), so the
    # nodes are kept as a list with one array shape(Nmembers, Nnodes, 3) per turbine
    Nturb = case["Nturbines"]
    nodes = [turb["Structure_nodes"] for turb in case["Turbine_parameters"]]
    Depth = case["Depth"]
    tags  = np.asarray([turb["Case_type"] for turb in case["Turbine_parameters"]])
    Wind_dir = np.asarray([turb["WindDir"] for turb in case["Turbine_parameters"]])

    # Every node of the farm, to compute extents. SAITEC's rigid members are NaN-padded: NaN-aware functions
    all_nodes = np.concatenate([turb.reshape(-1, 3) for turb in nodes])

    # Find node with maximum Z coordinate to place WindDir arrow
    max_z_nodes = np.zeros((Nturb, 3))
    min_z_nodes = np.zeros((Nturb, 3))
    for i in range(Nturb):
        xyz = nodes[i].reshape(-1,3)
        max_z_nodes[i] = xyz[np.nanargmax(xyz[:,2]),:]
        min_z_nodes[i] = xyz[np.nanargmin(xyz[:,2]),:]

    max_x, max_y, max_z = max_z_nodes[:,0], max_z_nodes[:,1], max_z_nodes[:,2]
    min_x, min_y, min_z = min_z_nodes[:,0], min_z_nodes[:,1], min_z_nodes[:,2]


    if mode == "3D":

        # Handle axes
        if ax_in is None:
                fig = plt.figure(figsize=figsize)
                ax  = fig.add_subplot(111, projection="3d")
        else:
            ax = ax_in
            fig = ax.figure

        # Compute extents for drawing planes
        xs = all_nodes[:,0]
        ys = all_nodes[:,1]
        zs = all_nodes[:,2]

        xmin, xmax = np.nanmin(xs), np.nanmax(xs)
        ymin, ymax = np.nanmin(ys), np.nanmax(ys)
        zmin, zmax = np.nanmin(zs), np.nanmax(zs)
        xpad = 0.1 * (xmax-xmin) if xmax>xmin else 1.0
        ypad = 0.1 * (ymax-ymin) if ymax>ymin else 1.0
        zpad = 0.1 * (zmax-zmin) if zmax>zmin else 1.0

        # Draw seabed and surface planes
        Xg = np.linspace(xmin-xpad, xmax+xpad, 2)
        Yg = np.linspace(ymin-ypad, ymax+ypad, 2)
        XX, YY = np.meshgrid(Xg, Yg, indexing='ij')

        ax.plot_surface(XX, YY, np.full_like(XX, -Depth), color='gray', alpha=0.5, shade=True, label="Seabed")
        ax.plot_surface(XX, YY, np.full_like(XX, 0.0), color="deepskyblue", alpha=0.25, shade=False, label="Surface")

        # Plot each member as connected nodes for each turbine
        for i in range(Nturb):
            turb = nodes[i]
            for member in turb:
                ax.plot(member[:,0], member[:,1], member[:,2], "-o", color="k", ms=3)

        # Draw wind direction arrow
        if WindDir:
            scale_length = np.mean([xmax - xmin, ymax - ymin, zmax - zmin])
            if scale_length <= 0:
                scale_length = 1.0
            arrow_length = scale_length * 0.5
            if arrow_length == 0:
                arrow_length = 0.75

            for i in range(Nturb):
                pu.draw_direction_arrow(ax, max_x[i], max_y[i], Wind_dir[i],
                                        length=arrow_length, z=max_z[i], color = "red",
                                        label=f"WindDir {Wind_dir[i]}")

        x_center = 0.5 * (xmin + xmax)
        y_center = 0.5 * (ymin + ymax)

        x_span = (xmax - xmin) + 2*xpad
        y_span = (ymax - ymin) + 2*ypad
        span_xy = max(x_span, y_span)

        x_limits = (x_center - 0.5*span_xy, x_center + 0.5*span_xy)
        y_limits = (y_center - 0.5*span_xy, y_center + 0.5*span_xy)

        z_span = (zmax - zmin) + 2*zpad
        z_limits = (zmin - zpad, zmax + zpad)

        ax.legend()
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        ax.set_title(f"3D View (XY)")
        ax.set_xlim(*x_limits)
        ax.set_ylim(*y_limits)
        ax.set_zlim(*z_limits)
        ax.set_box_aspect((span_xy, span_xy, z_span))

    else:

        # Handle axes
        if ax_in is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            ax = ax_in
            fig = ax.figure

        if mode == "xy":

            # Compute extents for drawing planes
            xs = all_nodes[:,0]
            ys = all_nodes[:,1]

            xmin, xmax = np.nanmin(xs), np.nanmax(xs)
            ymin, ymax = np.nanmin(ys), np.nanmax(ys)
            xpad = 0.1*(xmax-xmin) if xmax > xmin else 1.0
            ypad = 0.1*(ymax-ymin) if ymax > ymin else 1.0

            x_center = 0.5 * (xmin + xmax)
            y_center = 0.5 * (ymin + ymax)
            x_span = (xmax - xmin) + 2*xpad
            y_span = (ymax - ymin) + 2*ypad
            span_xy = max(x_span, y_span)

            x_limits = (x_center - 0.5*span_xy, x_center + 0.5*span_xy)
            y_limits = (y_center - 0.5*span_xy, y_center + 0.5*span_xy)

            # Project onto z = 0 plane
            for i in range(Nturb):
                turb = nodes[i]
                for member in turb:
                    ax.plot(member[:,0], member[:,1], '-o', color='k', ms=3)

            # Draw wind direction arrow
            if WindDir:
                for i in range(Nturb):
                    arrow_length = (xmax - xmin + 2*xpad) * 0.25
                    if arrow_length == 0: arrow_length = .75

                    pu.draw_direction_arrow(ax, max_x[i], max_y[i], Wind_dir[i],
                                            length=arrow_length, z=max_z[i], color = "red",
                                            label=f"WindDir {Wind_dir[i]}")

            ax.set_xlabel('x [m]')
            ax.set_ylabel('y [m]')
            ax.set_title(f"Top View (XY)")
            ax.set_xlim(*x_limits)
            ax.set_ylim(*y_limits)
            ax.set_aspect('equal', adjustable='box')
            if WindDir: ax.legend()     # The arrows are the only labelled artists of this view

        elif mode == 'xz':

            # Compute extents for drawing plane
            xs = all_nodes[:,0]
            zs = all_nodes[:,2]

            xmin, xmax = np.nanmin(xs), np.nanmax(xs)
            zmin, zmax = np.nanmin(zs), np.nanmax(zs)
            xpad = 0.1*(xmax-xmin) if xmax > xmin else 1.0
            zpad = 0.1*(zmax-zmin) if zmax > zmin else 1.0

            ax.axhline(0.0, color="deepskyblue", linestyle='--', linewidth=1.0, alpha=0.6, label='Surface')
            ax.axhline(-Depth, color='gray', linestyle='--', linewidth=1.0, alpha=0.6, label='Seabed')

            x_fill = np.array([xmin - xpad, xmax + xpad])
            ax.fill_between(x_fill, np.full_like(x_fill, -Depth), np.full_like(x_fill, 0.0), color='deepskyblue', alpha=0.25)

            # Project onto y = 0 plane
            for i in range(Nturb):
                turb = nodes[i]
                for member in turb:
                    ax.plot(member[:,0], member[:,2], '-o', color='k', ms=3)

            # Draw wind direction arrow (only visible in this plane when the wind blows along x)
            if WindDir:
                arrow_length = (xmax - xmin + 2*xpad) * 0.25
                if arrow_length == 0: arrow_length = .75
                for i in range(Nturb):
                    if Wind_dir[i] == 0.0 or Wind_dir[i] == 180.0:
                        pu.draw_direction_arrow(ax, max_x[i], max_z[i], Wind_dir[i],
                                            length=arrow_length, color = "red",
                                            label=f"WindDir {Wind_dir[i]}")

            ax.set_xlabel('x [m]')
            ax.set_ylabel('z [m]')
            ax.set_title(f"Side View (XZ)")
            ax.set_xlim(xmin - xpad, xmax + xpad)
            ax.set_ylim(zmin - zpad, zmax + zpad)
            ax.legend()

    if ax_in is None:
        pu.set_figure_name(fig, "structure", mode)

    return fig, ax


# ---------- Spectrum Plots ---------- #
def plot_spectrums(case         : dict     = None,      # [-] Dictionary containing simulation data
                   mode         : str      = 'SPL',     # [-] Selects Sound Pressure/Power Level 'SPL' or 'SWL'
                   absorption   : bool     = False,     # [-] Wheter to apply approximate absorption attenuation
                   octave       : bool     = True,      # [-] 1/3 octave or fine resolution
                   filter_under : float    = None,      # [Hz] Minimum frequency threshold
                   filter_over  : float    = None,      # [Hz] Maximum frequency threshold
                   do_thresholds: bool     = False,     # [-] Whether to plot regulatory/reference threholds
                   figsize      : tuple    = (10,5),    # [-] Figure size
                   ax_in        : plt.Axes = None):     # [-] Optional external axes for subplots
    """
    Plots the acoustic spectrum. Observers without data (NaN, inside a structure) are left out with a
    message.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The generated figure and axes objects for further editing.
    """

    if not case["has_spectrums"]:
        print("plot_spectrum(): no spectrum data, skipping")
        return None, None
    if absorption:
        print("plot_spectrums(): distances non defined since there are more than one turbine: Setting absorption = False")

    tags = np.asarray([turb["Case_type"] for turb in case["Turbine_parameters"]])
    colors = pu.get_case_color(tags=tags); color = colors[0]
    Nturb = case["Nturbines"]

    # Extract data
    Freqs     = case["Freqs"]
    Observers = case["Obs_spectrums"]
    p         = case["P_spectrums"]
    AxisPos   = np.asarray([turb["AxisPos"] for turb in case["Turbine_parameters"]])
    p_ref     = case["p_ref"]

    if not pu.is_valid(p, "plot_spectrums"):
        return None, None

    if mode == 'SPL':
        spl = au.pressure_to_SPL(p, p_ref)
        ylabel = r"SPL [dB re 1 $\mu$ Pa]"
    else:
        raise RuntimeError("plot_spectrum(): mode SWL not implemented yet")

    # Filter frequencies
    mask = np.ones_like(Freqs, dtype=bool)
    if filter_under is not None:
        mask &= (Freqs >= filter_under)
    if filter_over is not None:
        mask &= (Freqs <= filter_over)
    if not np.any(mask):
        print("plot_spectrums(): no frequencies match the specified filter bounds!")
        return None, None

    Freqs = Freqs[mask]
    spl   = spl[mask]
    oaspl = au.pressure_to_OASPL(p[mask], p_ref)

    if octave:
        Freqs, spl = au.to_third_octave(Freqs, spl, fmin=Freqs.min(), fmax=Freqs.max())

    # Plot
    if ax_in is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        ax = ax_in
        fig = ax.figure

    linestyles   = pu.get_line_styles()
    num_observers = spl.shape[1] if spl.ndim> 1 else 1
    if num_observers > len(linestyles): print("plot_spectrum(): too many observers, linestyles will be repeated")

    for i in range(num_observers):
        ls = linestyles[i % len(linestyles)]
        y_data = spl[:,i] if spl.ndim > 1 else spl
        if Observers.ndim > 1:
            coord = tuple(float(f"{val:.1f}") for val in Observers[i])
        else:
            coord = (float(f"{Observers[i]:.1f}"),)

        # Observers inside a structure carry no data (NaN at every frequency)
        if np.all(np.isnan(y_data)):
            print(f"plot_spectrums(): observer {coord} [m] has no data (NaN, inside a structure), skipping it")
            continue

        ax.semilogx(Freqs, y_data, ls, color=color, label=f"{coord} [m]")

        print(f"\nObserver coordinates: {coord} [m]")
        print(f"OASPL = {oaspl[i]:.1f}")

    # Animal thresholds
    if do_thresholds:
        thr = ed.Animal_threshold(do_plot=False)
        f_hz = thr['f'] * 1000 # kHz --> Hz

        THRESHOLD_STYLES = [
            (f_hz,                       thr['T_Lf'], ':',  [0.4, 0.4, 0.4], 'HT. Low-frequency'),
            (f_hz,                       thr['T_Mf'], '--', [0.7, 0.7, 0.7], 'HT. Mid-frequency'),
            (f_hz,                       thr['T_Hf'], '-.', [0.6, 0.6, 0.6], 'HT. High-frequency'),
            (f_hz,                       thr['T_PW'], '-',  [0.5, 0.5, 0.5], 'HT. Pinnipeds'),
            (thr['atlantic_herring'][:,0], thr['atlantic_herring'][:,1], '-+', [0.6, 0.6, 0.6], 'HT. Atlantic herring'),
            (thr['Atlantic_cod'][:,0],     thr['Atlantic_cod'][:,1],     '-.', [0.8, 0.8, 0.8], 'HT. Atlantic cod'),
        ]

        for (fx, fy, ls, color, label) in THRESHOLD_STYLES:
                ax.plot(fx, fy, ls, color=color, lw=1.2, alpha=0.6, label=label)

    ax.set_xlim(Freqs.min(), max(Freqs.max(), 1000))


    ax.set_ylabel(ylabel)
    ax.set_xlabel("Frequency [Hz]")
    ax.grid(True, which='both', ls='--', alpha=0.4)

    if do_thresholds:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.16), ncol=3, frameon=False)
        if ax_in is None:
            fig.tight_layout(rect=[0, 0.06, 1, 1])
    else:
        ax.legend()
        if ax_in is None:
            fig.tight_layout()

    if ax_in is None:
        pu.set_figure_name(fig, "spectrum", mode, "octave" if octave else "narrowband")

    return fig, ax


# ---------- Polar Plots ---------- #
def plot_polar(case             : dict  = None,        # [-] Dictionary containing simulation data
               mode             : str   = 'OASPL',     # [-] Magnitude to plot: 'OASPL', 'SPL', or 'ABS'
               target_frequency : float = 0.88,        # [Hz] Target frequency for 'SPL' or 'ABS' modes
               absorption       : bool  = False,        # [-] Whether to apply absorption attenuation
               bands            : list  = [10., 100.], # [Hz] Frequency band limits for OASPL
               figsize          : tuple = (7,7),       # [-] Figure size
               ax_in            : plt.Axes = None):    # [-] Optional external axes for subplots

    """
    Plots polar directivity diagram for the wind farm. Angles without data (NaN, inside a structure) are gaps
    in the curves and holes down to the centre in their fill.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The generated figure and axes objects for further editing.
    """

    if not case.get("has_polar", False):
        print("plot_polar(): no polar data, skipping")
        return None, None

    if absorption:
        print("plot_polar(): distances non defined since there are more than one turbine: Setting absorption = False")
        absorption = False

    colors = pu.get_band_colors()

    # Extract data
    theta      = case["Theta_deg_polar"]
    p          = case["P_polar"]
    Freqs      = case["Freqs"]
    p_ref      = case["p_ref"]
    r          = case["R_polar"]

    if not pu.is_valid(p, "plot_polar"):
        return None, None

    theta_rad  = np.deg2rad(theta)
    fmin, fmax = Freqs.min(), max(Freqs.max(), 1000)
    N_obs      = len(theta)

    # Filter frequencies within specified bands
    if not bands:
        band_labels = [f"All ({fmin:.2g}-{fmax:.2g} Hz)"]
        band_masks  = [np.ones(len(Freqs), dtype=bool)]
    else:
        cuts = sorted(bands)
        edges = [fmin] + cuts + [fmax]
        band_labels, band_masks = [], []
        for i in range(len(edges)-1):
            lo, hi = edges[i], edges[i+1]
            mask   = Freqs <= hi if i == 0 else (Freqs > lo) & (Freqs <= hi)
            band_masks.append(mask)
            band_labels.append(f"{lo:.4g}-{hi:.4g} Hz")

    n_bands = len(band_labels)
    if n_bands > len(colors):
        print("plot_polar(): too many bands, colors will be repeated")
    colors = [colors[i % len(colors)] for i in range(n_bands)]

    # Setup plot
    if ax_in is None:
        fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=figsize)
    else:
        ax = ax_in
        fig = ax.figure

    distances = np.full(N_obs, r)
    theta_plot = np.append(theta_rad, theta_rad[0])

    ax.set_theta_zero_location("W")
    ax.set_theta_direction(-1)

    mode = mode.upper()
    name_extra = ""

    if mode == "OASPL":
        all_oaspl_vals = []
        for label, mask, color in zip(band_labels, band_masks, colors):
            # A band without frequencies would give a fake -60 dB curve
            if not np.any(mask):
                print(f"plot_polar(): no frequencies in band {label}, skipping it")
                continue

            # p shape is (Nfreqs, Ntheta)
            oaspl = au.OASPL_band(p, mask, p_ref, absorption, Freqs, distances)
            all_oaspl_vals.append((oaspl, label, color))

        if all_oaspl_vals:
            vals = np.concatenate([oaspl for oaspl, _, _ in all_oaspl_vals])
            rmin, rmax = np.nanmin(vals), np.nanmax(vals)
            margin = abs(rmax-rmin) * 0.1 if rmax != rmin else abs(rmin) * 0.1
            if margin == 0: margin = 1.0

            for oaspl, label, color in all_oaspl_vals:
                oaspl_plot = np.append(oaspl, oaspl[0])

                # Plot (angles without data: gap in the curve, fill collapsed to the centre)
                ax.plot(theta_plot, oaspl_plot, lw=2, color=color, label=label)
                ax.fill(theta_plot, np.where(np.isfinite(oaspl_plot), oaspl_plot, rmin - margin), alpha=0.08, color=color)

            ax.set_rmin(rmin - margin)

    else:
        # Find closest frequency index for single-frequency modes
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f   = Freqs[idx]
        p_f = p[idx, :]
        name_extra = f"{f:.2f}Hz"

        vals, unit_label = pu.convert_to(p_f, mode=mode)

        if absorption and mode == "SPL":
            vals = au.add_absorption(f, vals, distances)

        vals_plot = np.append(vals, vals[0])
        color = colors[0]

        rmin, rmax = np.nanmin(vals), np.nanmax(vals)
        margin = abs(rmax-rmin) * 0.1 if rmax != rmin else abs(rmin) * 0.1
        if margin == 0: margin = 1.0 # Fallback

        # Plot (angles without data: gap in the curve, fill collapsed to the centre)
        ax.plot(theta_plot, vals_plot, lw=2, color=color, label= f"{unit_label} at {f:.2f} [Hz]")
        ax.fill(theta_plot, np.where(np.isfinite(vals_plot), vals_plot, rmin - margin), alpha=0.08, color=color)

        ax.set_rmin(rmin - margin)

    # Use a descriptive title since "tags" might involve multiple turb types in a farm
    ax.set_title(f"Wind Farm Polar Directivity (R = {r} m)")
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)

    if ax_in is None:
        fig.tight_layout(rect=[0, 0.06, 1, 1])
        pu.set_figure_name(fig, "polar", mode, name_extra)

    return fig, ax


# ---------- Cylindrical Plots ---------- #
def plot_cylinder(case            : dict = None,        # [-] Dictionary containing simulation data
                  mode            : str = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                  target_frequency: float = 0.88,       # [Hz] Target frequency for single-frequency modes
                  absorption      : bool = False,        # [-] Whether to apply absorption attenuation
                  filter_under    : float = None,       # [Hz] Lower frequency cutoff (f >= filter_under)
                  filter_over     : float = None,       # [Hz] Upper frequency cutoff (f <= filter_over)
                  figsize         : tuple = (12,5),     # [-] Figure size
                  cmap            : str = 'inferno',    # [-] Contour colormap
                  ax_in           : plt.Axes = None):   # [-] Optional external axes for subplots
    """
    Plots 2D unwrapped cylindrical surface maps for acoustic magnitudes (OASPL, SPL, SWL, etc.) for the wind
    farm. Observers without data (NaN, inside a structure) are left as holes in the contours.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The generated figure and axes objects for further editing.
    """

    if not case.get("has_cylinder", False):
        print("plot_cylinder(): no cylinder data, skipping")
        return None, None

    if absorption:
        print("plot_cylinder(): distances non defined since there are more than one turbine: Setting absorption = False")
        absorption = False

    # Extract data
    theta     = case["Theta_deg_cylinder"]   # 1D array of angles [deg]
    z         = case["Z_cylinder"]           # 1D array of height/depth [m]
    p         = case["P_cylinder"]           # Pressure grid
    Freqs     = case["Freqs"]
    p_ref     = case["p_ref"]
    r         = case["R_cylinder"]

    if not pu.is_valid(p, "plot_cylinder"):
        return None, None

    nf, ntheta, nz = len(Freqs), len(theta), len(z)

    # Reshape to flatten spatial dimensions for uniform processing
    p = p.reshape((nf, ntheta * nz))

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
        vals = au.pressure_to_OASPL(p, p_ref)
        unit_label = r"OASPL [dB re 1$\mu$ Pa]"

    elif mode == 'SWL':
        print("plot_cylinder(): SWL not yet implemented")
        return fig, ax
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f = Freqs[idx]
        p_f = p[idx, :]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = pu.convert_to(p_f, mode=mode)

    # Plot Setup: close the cylinder at 360 degrees
    theta_plot = np.append(theta, 360.)
    vals = vals.reshape((ntheta, nz))
    vals_plot = np.append(vals, vals[0:1, :], axis=0)

    # Matplotlib contour expects Z shape (len(y), len(x)), so transpose vals_plot. NaN: holes
    contour_set = ax.contourf(theta_plot, z, vals_plot.T, levels=60, cmap=cmap)
    ax.contour(theta_plot, z, vals_plot.T, levels=30, colors='k', linewidths=0.4, alpha=0.5)

    cbar = fig.colorbar(contour_set, ax=ax, orientation='vertical', pad=0.02)
    cbar.set_label(unit_label)

    ax.set_xlabel(r'Azimuth Angle, $\theta$ [deg]')
    ax.set_ylabel('Depth, z [m]')
    ax.set_title(f"Wind Farm Cylindrical Field - {mode} {title_extra}\n(R = {r} m)")
    ax.set_xlim(theta_plot.min(), theta_plot.max())
    ax.set_ylim(z.min(), z.max())
    ax.set_xticks(np.arange(0, 361, 45))

    if ax_in is None:
        fig.tight_layout()
        pu.set_figure_name(fig, "cylinder", mode, title_extra)

    return fig, ax


# ---------- Linear Plots ---------- #
def plot_line(case            : dict = None,        # [-] Dictionary containing simulation data
              mode            : str = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
              target_frequency: float = 0.88,       # [Hz] Target frequency for single-frequency modes
              absorption      : bool = False,        # [-] Whether to apply absorption attenuation
              filter_under    : float = None,       # [Hz] Lower frequency cutoff (f >= filter_under)
              filter_over     : float = None,       # [Hz] Upper frequency cutoff (f <= filter_over)
              figsize         : tuple = (10,5),     # [-] Figure size
              ax_in           : plt.Axes = None):   # [-] Optional external axes for subplots
    """
    Plots the acoustic pressure field along a 1D line for the wind farm. Points without data (NaN, inside a
    structure) are gaps in the curve.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The generated figure and axes objects for further editing.
    """

    if not case.get("has_line", False):
        print("plot_line(): no line data, skipping")
        return None, None

    if absorption:
        print("plot_line(): distances non defined since there are more than one turbine: Setting absorption = False")
        absorption = False

    # Extract data
    Freqs     = case["Freqs"]
    p         = case["P_line"]
    p_ref     = case["p_ref"]
    distances = case["Distances_line"]
    logspace  = case.get("Logspace_line", False)

    # Points info for the title (defaulting to 0s if something is missing)
    p1 = case.get("P1_line", [0, 0, 0])
    p2 = case.get("P2_line", [0, 0, 0])

    if not pu.is_valid(p, "plot_line"):
        return None, None

    nf = len(Freqs)

    # Filter frequency
    freq_mask = np.ones(nf, dtype=bool)
    if filter_under is not None:
        freq_mask &= (Freqs >= filter_under)
    if filter_over is not None:
        freq_mask &= (Freqs <= filter_over)

    if not np.any(freq_mask):
        print("plot_line(): no frequencies match the specified filter bounds!")
        return None, None

    p, Freqs = p[freq_mask], Freqs[freq_mask]

    mode = mode.upper()

    if ax_in is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = ax_in
        fig = ax.figure

    title_extra = ""

    # Select mode
    if mode == 'OASPL':
        vals = au.pressure_to_OASPL(p, p_ref)
        unit_label = r"OASPL [dB re 1$\mu$ Pa]"

    elif mode == 'SWL':
        print("plot_line(): SWL not yet implemented")
        return fig, ax
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f = Freqs[idx]
        p_f = p[idx, :]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = pu.convert_to(p_f, mode=mode)

    # Plot Setup
    ax.plot(distances, vals, '-', color='teal', linewidth=2.0, marker='o', markersize=4, markerfacecolor='white')

    ax.set_xlabel(r"Distance from along line [m]")
    ax.set_ylabel(unit_label)

    # Format points for title readability
    p1_str = f"({p1[0]:.1f}, {p1[1]:.1f}, {p1[2]:.1f})" if len(p1) == 3 else f"{p1}"
    p2_str = f"({p2[0]:.1f}, {p2[1]:.1f}, {p2[2]:.1f})" if len(p2) == 3 else f"{p2}"

    ax.set_title(f"Wind Farm Line Profile - {mode} {title_extra}\nFrom {p1_str} to {p2_str}")

    # Apply logarithmic scale to X axis if points were generated logarithmically
    if logspace:
        ax.set_xscale('log')

    ax.grid(True, which='both', linestyle='--', alpha=0.6)

    if ax_in is None:
        fig.tight_layout()
        pu.set_figure_name(fig, "line", mode, title_extra)

    return fig, ax


# ---------- Slice Plots ---------- #
def plot_sliceXY(case            : dict  = None,         # [-] Dictionary containing simulation data
                 mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                 target_frequency: float = 0.88,         # [Hz] Target frequency for single-frequency modes
                 absorption      : bool  = False,        # [-] Whether to apply absorption attenuation
                 filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                 filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                 structure       : bool  = True,         # [-] Wheter to plot structure on slice
                 figsize         : tuple = (7,7),        # [-] Figure size
                 cmap            : str   = 'inferno',    # [-] Contour colormap
                 ax_in           : plt.Axes = None):     # [-] Optional external axes for subplots
    """
    Plots a 2D spatial slice in the XY plane of acoustic magnitudes for the wind farm. Observers without data
    (NaN, inside a structure) are left as holes in the contours.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    if not case.get("has_slicexy", False):
        print("plot_sliceXY(): no XY slice data, skipping")
        return None, None

    if absorption:
        print("plot_sliceXY(): distances non defined since there are more than one turbine: Setting absorption = False")
        absorption = False

    # Extract data
    x         = case["X_slicexy"]           # 1D array of X coordinates [m]
    y         = case["Y_slicexy"]           # 1D array of Y coordinates [m]
    z_slice   = case["Z_slicexy"]           # Z level of the slice [m]
    p         = case["P_slicexy"]           # Pressure grid
    Freqs     = case["Freqs"]
    p_ref     = case["p_ref"]

    if not pu.is_valid(p, "plot_sliceXY"):
        return None, None

    Nfreqs, nx, ny = len(Freqs), len(x), len(y)

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
        vals = au.pressure_to_OASPL(p, p_ref)
        unit_label = r"OASPL [dB re 1$\mu$ Pa]"
    elif mode == 'SWL':
        print("WARNING plot_sliceXY(): SWL not yet implemented")
        return fig, ax
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f   = Freqs[idx]
        p_f = p[idx, :]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = pu.convert_to(p_f, mode=mode)

    # Plot (NaN: holes)
    contour_set = ax.contourf(X, Y, vals, levels=60, cmap=cmap)
    cbar = fig.colorbar(contour_set, ax=ax, orientation='horizontal', pad=0.1)
    cbar.set_label(unit_label)
    if abs(x[-1]-x[0]) == abs(y[-1]-y[0]): ax.set_aspect('equal')

    # Structure - Adapted for Wind Farm
    if structure:
        try:
            nodes_farm = [turb["Structure_nodes"] for turb in case.get("Turbine_parameters", []) if "Structure_nodes" in turb]
            for turb_nodes in nodes_farm:
                for member in turb_nodes:
                    ax.plot(member[:,0], member[:,1], '-o', color='black', markersize=2, linewidth=1.0, zorder=10)
        except Exception as e:
            print(f"plot_sliceXY() WARNING: Could not plot structure. ({type(e).__name__}: {e})")

    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title(f"Wind Farm XY Slice - {mode} {title_extra}\nDepth = {z_slice:.1f} [m]")
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())

    if ax_in is None:
        fig.tight_layout()
        pu.set_figure_name(fig, "sliceXY", mode, title_extra)

    return fig, ax

def plot_sliceXZ(case            : dict  = None,         # [-] Dictionary containing simulation data
                 mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                 target_frequency: float = 0.88,         # [Hz] Target frequency for single-frequency modes
                 absorption      : bool  = False,        # [-] Whether to apply absorption attenuation
                 filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                 filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                 structure       : bool  = True,         # [-] Wheter to plot structure on slice
                 figsize         : tuple = (7,7),        # [-] Figure size
                 cmap            : str   = 'inferno',    # [-] Contour colormap
                 ax_in           : plt.Axes = None):     # [-] Optional external axes for subplots
    """
    Plots a 2D spatial slice in the XZ vertical plane of acoustic magnitudes for the wind farm. Observers
    without data (NaN, inside a structure) are left as holes in the contours.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    if not case.get("has_slicexz", False):
        print("plot_sliceXZ(): no XZ slice data, skipping")
        return None, None

    if absorption:
        print("plot_sliceXZ(): distances non defined since there are more than one turbine: Setting absorption = False")
        absorption = False

    # Extract data
    x         = case["X_slicexz"]           # 1D array of X coordinates [m]
    z         = case["Z_slicexz"]           # 1D array of Z coordinates [m]
    y_slice   = case["Y_slicexz"]           # Y level of the slice [m]
    p         = case["P_slicexz"]           # Pressure grid (Nf, nx, nz)
    Freqs     = case["Freqs"]
    p_ref     = case["p_ref"]

    if not pu.is_valid(p, "plot_sliceXZ"):
        return None, None

    Nfreqs, nx, nz = len(Freqs), len(x), len(z)

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
        vals = au.pressure_to_OASPL(p, p_ref)
        unit_label = r"OASPL [dB re 1$\mu$ Pa]"
    elif mode == 'SWL':
        print("WARNING plot_sliceXZ(): SWL not yet implemented")
        return fig, ax
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f   = Freqs[idx]
        p_f = p[idx, :, :]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = pu.convert_to(p_f, mode=mode)

    # Plot (NaN: holes)
    contour_set = ax.contourf(X, Z, vals, levels=60, cmap=cmap)
    cbar = fig.colorbar(contour_set, ax=ax, orientation='horizontal', pad=0.1)
    cbar.set_label(unit_label)
    if abs(x[-1]-x[0]) == abs(z[-1]-z[0]): ax.set_aspect('equal')

    # Calculate max Z for plotting limits
    max_z_plot = z.max()

    # Structure - Adapted for Wind Farm
    if structure:
        try:
            nodes_farm = [turb["Structure_nodes"] for turb in case.get("Turbine_parameters", []) if "Structure_nodes" in turb]
            for turb_nodes in nodes_farm:
                for member in turb_nodes:
                    ax.plot(member[:,0], member[:,2], '-o', color='black', markersize=2, linewidth=1.0, zorder=10)

                    # Update max Z to properly scale the vertical axis (NaN-padded members: nanmax)
                    current_max_z = np.nanmax(member[:,2]) if np.any(np.isfinite(member[:,2])) else max_z_plot
                    if current_max_z > max_z_plot:
                        max_z_plot = current_max_z

        except Exception as e:
            print(f"plot_sliceXZ() WARNING: Could not plot structure. ({type(e).__name__}: {e})")

    ax.set_xlabel('x [m]')
    ax.set_ylabel('z [m]')
    ax.set_title(f"Wind Farm XZ Slice - {mode} {title_extra}\nY = {y_slice:.1f} [m]")
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(z.min(), max_z_plot)

    if ax_in is None:
        fig.tight_layout()
        pu.set_figure_name(fig, "sliceXZ", mode, title_extra)

    return fig, ax

def plot_sliceVertical(case            : dict  = None,         # [-] Dictionary containing simulation data
                       mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                       target_frequency: float = 0.88,         # [Hz] Target frequency for single-frequency modes
                       absorption      : bool  = False,        # [-] Whether to apply absorption attenuation
                       filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                       filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                       structure       : bool  = True,         # [-] Wheter to plot structure on slice
                       figsize         : tuple = (7,7),        # [-] Figure size
                       cmap            : str   = 'inferno',    # [-] Contour colormap
                       ax_in           : plt.Axes = None):     # [-] Optional external axes for subplots
    """
    Plots a 2D spatial slice in an arbitrary vertical plane of acoustic magnitudes for the wind farm. The
    plane normal is (cos az, sin az) and its horizontal axis u runs along (-sin az, cos az) from
    Center_sliceV (AcousticSolver.run_sliceVertical). The plane usually crosses the piles: observers without
    data (NaN, inside a structure) are left as holes in the contours.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    if not case.get("has_slicevertical", False):
        print("plot_sliceVertical(): no vertical slice data, skipping")
        return None, None

    if absorption:
        print("plot_sliceVertical(): distances non defined since there are more than one turbine: Setting absorption = False")
        absorption = False

    # Extract data
    u         = case["U_sliceV"]            # 1D array of along-slice coordinates [m]
    z         = case["Z_sliceV"]            # 1D array of Z coordinates [m]
    p         = case["P_sliceV"]            # Pressure grid (Nf, nu, nz)
    azimuth   = case["Azimuth_sliceV"]      # Azimuth angle of the slice [deg]
    center    = case["Center_sliceV"]       # Center point [cx, cy]
    Freqs     = case["Freqs"]
    p_ref     = case["p_ref"]

    if not pu.is_valid(p, "plot_sliceVertical"):
        return None, None

    Nfreqs, nu, nz = len(Freqs), len(u), len(z)

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
        vals = au.pressure_to_OASPL(p, p_ref)
        unit_label = r"OASPL [dB re 1$\mu$ Pa]"
    elif mode == 'SWL':
        print("WARNING plot_sliceVertical(): SWL not yet implemented")
        return fig, ax
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f   = Freqs[idx]
        p_f = p[idx, :, :]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = pu.convert_to(p_f, mode=mode)

    # Plot (NaN: holes)
    contour_set = ax.contourf(U, Z, vals, levels=60, cmap=cmap)
    cbar = fig.colorbar(contour_set, ax=ax, orientation='horizontal', pad=0.1)
    cbar.set_label(unit_label)
    if abs(u[-1]-u[0]) == abs(z[-1]-z[0]): ax.set_aspect('equal')

    # Calculate max Z for plotting limits
    max_z_plot = z.max()

    # Structure - Projected onto the arbitrary vertical slice
    if structure:
        try:
            cx, cy = center[0], center[1]
            theta_rad = np.deg2rad(azimuth)
            cos_t = np.cos(theta_rad)
            sin_t = np.sin(theta_rad)

            nodes_farm = [turb["Structure_nodes"] for turb in case.get("Turbine_parameters", []) if "Structure_nodes" in turb]
            for turb_nodes in nodes_farm:
                for member in turb_nodes:
                    # Project X, Y nodes onto the U axis (along the slice direction): the plane normal is
                    # (cos az, sin az), so the axis U_sliceV runs along is (-sin az, cos az)
                    x_nodes = member[:,0]
                    y_nodes = member[:,1]
                    z_nodes = member[:,2]

                    u_proj = -(x_nodes - cx) * sin_t + (y_nodes - cy) * cos_t

                    ax.plot(u_proj, z_nodes, '-o', color='black', markersize=2, linewidth=1.0, zorder=10)

                    # Update max Z to properly scale the vertical axis (NaN-padded members: nanmax)
                    current_max_z = np.nanmax(z_nodes) if np.any(np.isfinite(z_nodes)) else max_z_plot
                    if current_max_z > max_z_plot:
                        max_z_plot = current_max_z

        except Exception as e:
            print(f"plot_sliceVertical() WARNING: Could not plot structure. ({type(e).__name__}: {e})")

    ax.set_xlabel('u (Distance along slice) [m]')
    ax.set_ylabel('z [m]')
    ax.set_title(f"Wind Farm Vertical Slice - {mode} {title_extra}\nAzimuth = {azimuth:.1f}° $|$ Center = ({center[0]:.1f}, {center[1]:.1f})")
    ax.set_xlim(u.min(), u.max())
    ax.set_ylim(z.min(), max_z_plot)

    if ax_in is None:
        fig.tight_layout()
        pu.set_figure_name(fig, "sliceVertical", mode, title_extra)

    return fig, ax

