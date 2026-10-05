"""
Module: plot_single_case.py
Description: Plots of a single wind turbine (.hdf5 files of FORTRAN/turbine_acoustic_data): structure,
             spectrum, polar directivity, cylinder, distance decay, line, XY/XZ/vertical slices, and the
             sphere/cylinder radiation metrics. Observers without data (NaN pressure, inside the structure)
             are shown as holes and left out of the integral metrics.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 09/2026
"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
import plot_utils as pu
import utils.AcousticUtils as au
import utils.EnvironmentalData as ed
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# NumPy >= 2.0 renamed trapz to trapezoid
trapezoid = getattr(np, "trapezoid", None) or np.trapz


# ---------- Helpers ---------- #
def get_curve(case            : dict  = None,     # [-] Dictionary containing simulation data
              step            : str   = None,     # [-] Step whose observers form the curve: 'decay' or 'line'
              mode            : str   = 'OASPL',  # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
              target_frequency: float = 0.88,     # [Hz] Target frequency for single-frequency modes
              absorption      : bool  = True,     # [-] Whether to apply absorption attenuation
              filter_under    : float = None,     # [Hz] Lower frequency cutoff (f >= filter_under)
              filter_over     : float = None,     # [Hz] Upper frequency cutoff (f <= filter_over)
              name            : str   = None):    # [-] Name of the calling plot, for the messages
    """
    Magnitude of the pressure of one case along the distances of its 'decay' or 'line' step: the curve plot_distance_decay
    and plot_line draw (one per wind speed with --ws). Distances without data (NaN) stay NaN.

    Returns
    -------
    r, vals, unit_label, title_extra : distances, magnitude, label of the y axis and note of the title
        ('at 0.88 Hz' for the single-frequency modes, empty for OASPL). None if the case cannot be plotted.
    """

    p     = case[f"P_{step}"]
    Freqs = case["Freqs"]
    df    = pu.get_bin_width(Freqs)         # Frequency bin width [Hz] for 'PSD'
    r     = case[f"Distances_{step}"]
    p_ref = case["p_ref"]

    if not pu.is_valid(p, name):
        return None

    # Filter frequencies
    freq_mask = np.ones(len(Freqs), dtype=bool)
    if filter_under is not None:
        freq_mask &= (Freqs >= filter_under)
    if filter_over is not None:
        freq_mask &= (Freqs <= filter_over)
    if not np.any(freq_mask):
        print(f"{name}(): no frequencies match the specified filter bounds!")
        return None

    p, Freqs = p[freq_mask], Freqs[freq_mask]

    mode = mode.upper()

    # Select mode
    title_extra = ""
    if mode == 'OASPL':
        vals = au.pressure_to_OASPL(p, p_ref, absorption, Freqs, r, enbw=case["ENBW"])
        unit_label = r"OASPL [dB re $\mu$ Pa]"
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f   = Freqs[idx]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = pu.convert_to(p[idx, :], mode, df, enbw=case["ENBW"])
        if absorption and mode in ('SPL', 'PSD'): vals = au.add_absorption(f, vals, r)

    return r, vals, unit_label, title_extra

def get_ws_curves(case            : dict  = None,     # [-] Dictionary containing simulation data
                  step            : str   = None,     # [-] 'decay' or 'line'
                  name            : str   = None,     # [-] Name of the calling plot, for the messages
                  **kwargs):                          # [-] Arguments of get_curve
    """
    Curves of a plot of a 'decay' or 'line' step: the case alone or, with --ws (case["WS_cases"]), one per wind
    speed. Wind speeds whose file lacks the step, or that cannot be plotted, are left out with a message.

    Returns
    -------
    curves, ws_mode : list of (dict, tuple), bool
        Case and its get_curve() result, by increasing wind speed; whether the curves are wind speeds.
    """

    cases, ws_mode = pu.get_ws_cases(case)

    curves = []
    for c in cases:
        if not c[f"has_{step}"]:
            print(f"{name}(): {c['Name']} has no {step} data, left out")
            continue

        curve = get_curve(c, step, name=name, **kwargs)
        if curve is not None:
            curves.append((c, curve))

    return curves, ws_mode


# ---------- Structural Plots ---------- #
def plot_structure(case   : dict  = None,       # [-] Dictionary containing simulation data
                   mode   : str   = '3D',       # [-] Plot mode for the structure representation
                   WindDir: bool  = True,       # [-] Wheter to plot or not wind direction
                   figsize: tuple = (8,6),      # [-] Figure size
                   ax_in  : plt.Axes = None):   # [-] Optional external axes for subplots
    """
    Plots the structure nodes in shape (Nmembers, Nnodes, 3)

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The generated figure and axes objects for further editing.
    """

    case["has_nodes"] = "Structure_nodes" in case
    if not case["has_nodes"]:
        print("plot_structure(): no structure data, skipping")
        return None, None

    if mode not in ('3D', 'xy', 'xz'):
        raise ValueError("mode must be '3D', 'xy', or 'xz'")

    # Extract nodes and depth
    nodes = case["Structure_nodes"]
    Depth = case["Depth"]
    tag   = case["Case_type"]

    # Find node with maximum Z coordinate to place WindDir arrow. Some
    # members (e.g. SAITEC's rigid links, which have fewer real nodes than
    # Nnodes) are NaN-padded, so use the NaN-aware argmax.
    flat_nodes = nodes.reshape(-1, 3)
    max_z_node = flat_nodes[np.nanargmax(flat_nodes[:,2])]
    max_x, max_y, max_z = max_z_node

    if mode == '3D':

        if ax_in is None:
            fig = plt.figure(figsize=figsize)
            ax  = fig.add_subplot(111, projection='3d')
        else:
            ax = ax_in
            fig = ax.figure

        # Compute extents for drawing planes (NaN-safe: some members may be
        # NaN-padded, e.g. SAITEC's rigid-link joints)
        xs = nodes[..., 0].ravel()
        ys = nodes[..., 1].ravel()
        xmin, xmax = np.nanmin(xs), np.nanmax(xs)
        ymin, ymax = np.nanmin(ys), np.nanmax(ys)
        xpad = 0.1 * (xmax - xmin) if xmax > xmin else 1.0
        ypad = 0.1 * (ymax - ymin) if ymax > ymin else 1.0

        # Draw seabed and surface planes
        Xg = np.linspace(xmin - xpad, xmax + xpad, 2)
        Yg = np.linspace(ymin - ypad, ymax + ypad, 2)
        XX, YY = np.meshgrid(Xg, Yg)

        ax.plot_surface(XX, YY, np.full_like(XX, -Depth), color = 'gray', alpha = 0.5, shade = True, label = "Seabed")
        ax.plot_surface(XX, YY, np.full_like(XX, 0.0), color='deepskyblue', alpha=0.25, shade=False, label = "Surface")

        # Plot each member as connected nodes
        for member in nodes:
            ax.plot(member[:, 0], member[:, 1], member[:, 2], '-o', color='k', ms=3)

        # Draw Wind Direction Arrow
        if WindDir:
            arrow_length = (xmax - xmin + 2*xpad) * 0.25
            if arrow_length == 0: arrow_length = .75

            pu.draw_direction_arrow(ax, max_x, max_y, case["WindDir"],
                                  length=arrow_length, z=max_z,
                                  color='red',
                                  label=f'WindDir: {case["WindDir"]}°')
        ax.legend()

        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_zlabel('z [m]')
        ax.set_title(tag)

        ax.set_xlim(xmin - xpad, xmax + xpad)
        ax.set_ylim(ymin - ypad, ymax + ypad)
        ax.set_zlim(-Depth, max_z + max(10, abs(max_z)*0.1))

    else:

        if ax_in is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            ax = ax_in
            fig = ax.figure

        if mode == "xy":
            xs = nodes[..., 0].ravel()
            ys = nodes[..., 1].ravel()
            xmin, xmax = np.nanmin(xs), np.nanmax(xs)
            ymin, ymax = np.nanmin(ys), np.nanmax(ys)
            xpad = 0.1 * (xmax - xmin) if xmax > xmin else 1.0
            ypad = 0.1 * (ymax - ymin) if ymax > ymin else 1.0

            # Project onto z = 0 plane
            for member in nodes:
                ax.plot(member[:, 0], member[:, 1], '-o', color='k', ms=3)

            if WindDir:
                arrow_length = (xmax - xmin + 2*xpad) * 0.25
                if arrow_length == 0: arrow_length = .75

                pu.draw_direction_arrow(ax, max_x, max_y, case["WindDir"],
                      length=arrow_length, color='red',
                      label=f'WindDir: {case["WindDir"]}°')

            ax.set_xlabel('x [m]')
            ax.set_ylabel('y [m]')
            ax.set_title(f"{tag} - Top View (XY)")
            ax.set_xlim(xmin - xpad, xmax + xpad)
            ax.set_ylim(ymin - ypad, ymax + ypad)
            ax.set_aspect('equal', adjustable='box')
            if WindDir: ax.legend()     # The arrow is the only labelled artist of this view

        elif mode == "xz":
            xs = nodes[..., 0].ravel()
            zs = nodes[..., 2].ravel()
            xmin, xmax = np.nanmin(xs), np.nanmax(xs)
            zmin, zmax = np.nanmin(zs), np.nanmax(zs)
            xpad = 0.1 * (xmax - xmin) if xmax > xmin else 1.0
            zpad = 0.1 * (zmax - zmin) if zmax > zmin else 1.0

            ax.axhline(0.0, color='deepskyblue', linestyle='--', linewidth=1.0, alpha=0.6, label='Surface')
            ax.axhline(-Depth, color='gray', linestyle='--', linewidth=1.0, alpha=0.6, label="Seabed")

            for member in nodes:
                ax.plot(member[:, 0], member[:, 2], '-o', color='k', ms=3)

            if WindDir:
                arrow_length = (xmax - xmin + 2*xpad) * 0.25
                if arrow_length == 0: arrow_length = .75

                pu.draw_direction_arrow(ax, max_x, max_z, case["WindDir"],
                      length=arrow_length, color='red',
                      label=f'WindDir: {case["WindDir"]}°')

            ax.set_xlabel('x [m]')
            ax.set_ylabel('z [m]')
            ax.set_title(f"{tag} - Side View (XZ)")
            ax.set_xlim(xmin - xpad, xmax + xpad)
            ax.set_ylim(zmin - zpad, zmax + zpad)
            ax.legend()

    if ax_in is None:
        fig.tight_layout()
        pu.set_figure_name(fig, "structure", mode)

    return fig, ax


# ---------- Spectrum Plots ---------- #
def get_spectrum(case         : dict  = None,      # [-] Dictionary containing simulation data
                 mode         : str   = 'SPL',     # [-] 'SPL' or 'PSD'
                 absorption   : bool  = True,      # [-] Wheter to apply approximate absorption attenuation
                 octave       : bool  = True,      # [-] 1/3 octave or fine resolution
                 filter_under : float = None,      # [Hz] Minimum frequency threshold
                 filter_over  : float = None):     # [Hz] Maximum frequency threshold
    """
    Spectrum of every observer of one case, as plot_spectrum draws it (once per wind speed with --ws): SPL or
    PSD in dB, narrowband or 1/3 octave, within the frequency filter, and the OASPL of every observer.

    Returns
    -------
    spectrum : dict
        Freqs, spl (Nf, Nobs), oaspl (Nobs,), distances (Nobs,) to the turbine axis and ylabel. None if the
        case cannot be plotted.
    """

    Freqs     = case["Freqs"]
    Observers = case["Obs_spectrums"]
    p         = case["P_spectrums"]
    AxisPos   = case["AxisPos"]
    p_ref     = case["p_ref"]

    if not pu.is_valid(p, "plot_spectrum"):
        return None

    # Distance to turbine axis
    distances = np.sqrt((Observers[:,0] - AxisPos[0])**2 + (Observers[:,1] - AxisPos[1])**2)

    mode = mode.upper()

    if mode == "SPL":
        spl = au.pressure_to_SPL(p, p_ref, absorption, Freqs, distances)
        ylabel = r"SPL [dB re 1 $\mu$ Pa]"
    elif mode == "PSD":
        df = pu.get_bin_width(Freqs)
        if df is None: raise ValueError("plot_spectrum(): mode 'PSD' needs at least two frequencies to get the bin width")
        spl = au.pressure_to_PSD(p, df, p_ref, absorption, Freqs, distances, enbw=case["ENBW"])
        ylabel = r"PSD [dB re 1 $\mu$Pa$^2$/Hz]"
    else:
        raise RuntimeError(f"plot_spectrum(): mode {mode} not implemented yet")


    # Filter frequencies
    mask = np.ones_like(Freqs, dtype=bool)
    if filter_under is not None:
        mask &= (Freqs >= filter_under)
    if filter_over is not None:
        mask &= (Freqs <= filter_over)
    if not np.any(mask):
        print("plot_spectrum(): no frequencies match the specified filter bounds!")
        return None

    Freqs = Freqs[mask]
    spl   = spl[mask]
    oaspl = au.pressure_to_OASPL(p[mask], p_ref, absorption, Freqs, distances, enbw=case["ENBW"])

    if octave:
        Freqs, spl = au.to_third_octave(Freqs, spl, fmin=Freqs.min(), fmax=Freqs.max())

        # PSD is a density: the energetic sum of its bins is band power / df, so rescale
        # it to band power / band width (spectrum level), independent of the resolution
        if mode == "PSD":
            bw    = Freqs * (2**(1/6) - 2**(-1/6))      # [Hz] 1/3 octave band widths
            shape = [-1] + [1] * (spl.ndim - 1)
            spl   = spl - 10.0 * np.log10(bw / df).reshape(shape)

    return dict(Freqs=Freqs, spl=spl, oaspl=oaspl, distances=distances, ylabel=ylabel)

def plot_spectrum(case         : dict  = None,      # [-] Dictionary containing simulation data
                  mode         : str   = 'SPL',     # [-] Selects Sound Pressure/Power Level or PSD. 'SPL', 'PSD' or 'SWL'
                  absorption   : bool  = True,      # [-] Wheter to apply approximate absorption attenuation
                  octave       : bool  = True,      # [-] 1/3 octave or fine resolution
                  filter_under : float = None,      # [Hz] Minimum frequency threshold
                  filter_over  : float = None,      # [Hz] Maximum frequency threshold
                  do_thresholds: bool  = False,     # [-] Whether to plot regulatory/reference thresholds
                  smooth       : bool  = True,      # [-] Narrowband only: smooth every observer curve too (False: only the means)
                  smooth_fraction: int = 12,        # [-] Narrowband only: smoothing window, 1/N octave (None: no smoothing at all)
                  figsize      : tuple = (10,5),    # [-] Figure size
                  ax_in        : plt.Axes = None):  # [-] Optional external axes for subplots
    """
    Plots the acoustic spectrum. Observers without data (NaN, inside the structure) are left out with a
    message. If every observer with data lies at the same distance to the axis, their energetic mean
    (solid black) and arithmetic dB mean (dashed black) are added. In narrowband the means are always
    smoothed over 1/smooth_fraction octave, so they read as a trend instead of a saw, and with smooth=True
    so is every observer curve. The means are computed from the raw curves either way, and OASPL never
    sees the smoothing.

    With --ws the spectra of every wind speed are drawn: color = wind speed and line style = observer, with a
    legend of each below the axes. Observers at the same distance to the axis are not drawn one by one but
    as their energetic mean and the range they span (shaded), so the plot stays readable.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The generated figure and axes objects for further editing.
    """

    if not case["has_spectrums"]:
        print("plot_spectrum(): no spectrum data, skipping")
        return None, None

    tag   = case["Case_type"]
    color = pu.get_case_color(case)

    cases, ws_mode = pu.get_ws_cases(case)

    spectra = []
    for c in cases:
        if not c["has_spectrums"]:
            print(f"plot_spectrum(): {c['Name']} has no spectrum data, left out")
            continue

        spectrum = get_spectrum(c, mode, absorption, octave, filter_under, filter_over)
        if spectrum is not None:
            spectra.append((c, spectrum))

    if not spectra:
        return None, None

    ylabel = spectra[0][1]["ylabel"]
    Freqs  = spectra[0][1]["Freqs"]


    # Plot
    if ax_in is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        ax = ax_in
        fig = ax.figure

    linestyles = pu.get_line_styles()

    # Narrowband smoothing: of every curve (smooth=True) or only of the means; the raw spl feeds the means
    smooth_means  = (not octave) and bool(smooth_fraction)
    smooth_curves = smooth_means and smooth
    smooth_note   = f"1/{smooth_fraction} oct. smoothed" if smooth_means else ""

    ws_handles, ws_labels = [], []                        # --ws: one entry per wind speed (color)
    style_handles, style_labels = [], []                  # --ws: one entry per observer (line style) or mean

    for c, spectrum in spectra:
        Freqs, spl, oaspl, distances = spectrum["Freqs"], spectrum["spl"], spectrum["oaspl"], spectrum["distances"]
        Observers = c["Obs_spectrums"]

        num_observers = spl.shape[1] if spl.ndim> 1 else 1
        if num_observers > len(linestyles): print("plot_spectrum(): too many observers, linestyles will be repeated")

        spl_plot = au.smooth_fractional_octave(Freqs, spl, smooth_fraction, domain='db') if smooth_curves else spl

        # Observers with data: those inside the structure carry NaN at every frequency
        valid = [i for i in range(num_observers) if not np.all(np.isnan(spl_plot[:,i] if spl_plot.ndim > 1 else spl_plot))]

        # Observers all at the same distance to the axis (e.g. around the platform): energetic and dB means
        same_distance = len(valid) > 1 and np.allclose(distances[valid], distances[valid[0]], rtol=0.0, atol=0.1)

        if ws_mode:
            color  = c["WS_color"]
            kwargs = dict(lw=c["WS_lw"], zorder=c["WS_zorder"])
            ws_handles.append(Line2D([], [], color=color, lw=c["WS_lw"]))
            ws_labels.append(pu.get_ws_label(c))
        else:
            kwargs = {}

        shades = [color] * num_observers if ws_mode else pu.get_color_shades(color, num_observers)    # One shade per observer

        for i in range(num_observers):
            ls = linestyles[i % len(linestyles)]
            y_data = spl_plot[:,i] if spl_plot.ndim > 1 else spl_plot

            # Observers are told apart by their coordinates: several of them may share the same distance
            x_obs, y_obs, z_obs = Observers[i]
            coords = f"({x_obs:.1f}, {y_obs:.1f}, {z_obs:.1f}) m"

            # Observers inside the structure carry no data (NaN at every frequency)
            if i not in valid:
                print(f"plot_spectrum(): observer {coords} has no data (NaN, inside the structure), skipping it")
                continue

            if not (ws_mode and same_distance):
                ax.semilogx(Freqs, y_data, ls, color=shades[i], label=None if ws_mode else coords, **kwargs)

            # The observers are the same at every wind speed: the line styles are explained once
            if ws_mode and c is spectra[0][0] and not same_distance:
                style_handles.append(ax.plot([], [], ls, color='k')[0])
                style_labels.append(coords)

            print(f"\nObserver: {coords} | r = {distances[i]:.1f} m")
            print(f"{tag}{', ' + pu.get_ws_label(c) if ws_mode else ''}: OASPL = {oaspl[i]:.1f}")

        if same_distance:
            spl_valid = spl[:, valid]
            mean_energetic = 10.0 * np.log10(np.nanmean(10.0**(spl_valid / 10.0), axis=1))    # Mean of p^2, in dB
            mean_db        = np.nanmean(spl_valid, axis=1)                                     # Mean of the levels

            # Narrowband: each mean smoothed in its own domain (power for the energetic one, dB for the other)
            if smooth_means:
                mean_energetic = au.smooth_fractional_octave(Freqs, mean_energetic, smooth_fraction, domain='db')
                mean_db        = au.smooth_fractional_octave(Freqs, mean_db, smooth_fraction, domain='mean')

            # Only the means smoothed: the note goes on their labels (otherwise it goes on the title)
            label_extra = f", {smooth_note}" if smooth_means and not smooth_curves else ""

            if ws_mode:
                with warnings.catch_warnings():     # A frequency where every observer is NaN stays NaN
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    spl_min, spl_max = np.nanmin(spl_valid, axis=1), np.nanmax(spl_valid, axis=1)
                if smooth_means:
                    spl_min = au.smooth_fractional_octave(Freqs, spl_min, smooth_fraction, domain='db')
                    spl_max = au.smooth_fractional_octave(Freqs, spl_max, smooth_fraction, domain='db')

                ax.semilogx(Freqs, mean_energetic, '-', color=color, **kwargs)
                ax.fill_between(Freqs, spl_min, spl_max, color=color, alpha=0.15, lw=0, zorder=kwargs["zorder"])

                if c is spectra[0][0]:
                    style_handles += [Line2D([], [], color='k', lw=2), Patch(facecolor='k', alpha=0.15, lw=0)]
                    style_labels  += [f"Energetic mean{label_extra}", "Range of the observers"]
            else:
                ax.semilogx(Freqs, mean_energetic, '-' , color='k', label=f"Energetic mean{label_extra}")
                ax.semilogx(Freqs, mean_db       , '--', color='k', label=f"Arithmetic mean (dB){label_extra}")

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
    ax.set_title(f"{tag} ({smooth_note})" if smooth_curves else tag)
    ax.grid(True, which='both', ls='--', alpha=0.4)
    ax.set_ylim(bottom=0.)

    if ws_mode:
        # Thresholds join the line styles in the second legend
        if do_thresholds:
            handles, labels = ax.get_legend_handles_labels()
            style_handles += handles
            style_labels  += labels

        pu.add_ws_legends(ax, ws_handles, ws_labels, style_handles, style_labels,
                          style_title=None if do_thresholds else ("Observers" if style_handles else None),
                          ncol_style=2 if do_thresholds else 1)
        if ax_in is None:
            fig.tight_layout()
    elif do_thresholds:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.16), ncol=3, frameon=False)
        if ax_in is None:
            fig.tight_layout(rect=[0, 0.06, 1, 1])
    else:
        ax.legend()
        if ax_in is None:
            fig.tight_layout()

    if ax_in is None:
        pu.set_figure_name(fig, "spectrum", mode.upper(), "octave" if octave else "narrowband", "allws" if ws_mode else "")

    return fig, ax


# ---------- Polar Plots ---------- #
def get_bands(Freqs: np.ndarray = None,     # [Hz] Frequencies of the case
              bands: list       = None):    # [Hz] Frequency band limits for OASPL
    """
    Frequency bands of the OASPL polar plot, from the lowest frequency to the highest (at least 1000 Hz).

    Returns
    -------
    band_labels, band_masks : list of str, list of np.ndarray
        Label of every band and boolean mask of its frequencies (all of them if there are no limits).
    """

    fmin, fmax = Freqs.min(), max(Freqs.max(), 1000)

    if not bands:
        return [f"All ({fmin:.2g}-{fmax:.2g} Hz)"], [np.ones(len(Freqs), dtype=bool)]

    cuts = sorted(bands)
    edges = [fmin] + cuts + [fmax]
    band_labels, band_masks = [], []
    for i in range(len(edges)-1):
        lo, hi = edges[i], edges[i+1]
        mask   = Freqs <= hi if i == 0 else (Freqs > lo) & (Freqs <= hi)
        band_masks.append(mask)
        band_labels.append(f"{lo:.4g}-{hi:.4g} Hz")

    return band_labels, band_masks

def plot_polar(case             : dict  = None,        # [-] Dictionary containing simulation data
               mode             : str   = 'OASPL',     # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', or 'ABS'
               target_frequency : float = 0.88,        # [Hz] Target frequency for 'SPL', 'PSD' or 'ABS' modes
               absorption       : bool  = True,        # [-] Whether to apply absorption attenuation
               bands            : list  = [10., 100.], # [Hz] Frequency band limits for OASPL
               figsize          : tuple = (7,7),       # [-] Figure size
               ax_in            : plt.Axes = None):    # [-] Optional external axes for subplots
    """
    Plots polar directivity diagram. Angles without data (NaN, inside the structure) are gaps in the curves
    and holes down to the centre in their fill.

    With --ws every wind speed is drawn: color = wind speed and, in OASPL, line style = frequency band, with
    a legend of each below the axes. The curves are not filled.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The generated figure and axes objects for further editing.
    """

    if not case["has_polar"]:
        print("plot_polar(): no polar data, skipping")
        return None, None

    tag    = case["Case_type"]
    colors = pu.get_band_colors()

    cases, ws_mode = pu.get_ws_cases(case)

    polars = []
    for c in cases:
        if not c["has_polar"]:
            print(f"plot_polar(): {c['Name']} has no polar data, left out")
        elif pu.is_valid(c["P_polar"], "plot_polar"):
            polars.append(c)

    if not polars:
        return None, None

    band_labels, _ = get_bands(polars[0]["Freqs"], bands)
    n_bands = len(band_labels)
    if n_bands > len(colors): print("plot_polar(): too many bands, colors will be repeated")
    colors = [colors[i % len(colors)] for i in range(n_bands)]


    # Setup plot
    if ax_in is None:
        fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=figsize)
    else:
        ax = ax_in
        fig = ax.figure

    ax.set_theta_zero_location("W")
    ax.set_theta_direction(-1)

    mode = mode.upper()
    name_extra = ""

    # Curves of every case: (case, angles, magnitude, label, band)
    series = []
    for c in polars:
        p          = c["P_polar"]
        Freqs      = c["Freqs"]
        df         = pu.get_bin_width(Freqs)   # Frequency bin width [Hz] for 'PSD'
        p_ref      = c["p_ref"]
        r          = c["R_polar"]
        N_obs      = len(c["Obs_polar"])

        theta_rad  = np.deg2rad(c["Theta_deg_polar"])
        distances  = np.full(N_obs, r)
        theta_plot = np.append(theta_rad, theta_rad[0])

        if mode == "OASPL":
            labels, masks = get_bands(Freqs, bands)
            for b, (label, mask) in enumerate(zip(labels, masks)):
                # A band without frequencies would give a fake -60 dB curve
                if not np.any(mask):
                    print(f"plot_polar(): no frequencies in band {label}, skipping it")
                    continue

                oaspl = au.OASPL_band(p, mask, p_ref, absorption, Freqs, distances, enbw=c["ENBW"])
                series.append((c, theta_plot, np.append(oaspl, oaspl[0]), label, b))

        else:
            # Find closest frequency index for single-frequency modes
            idx = np.argmin(np.abs(Freqs - target_frequency))
            f   = Freqs[idx]
            name_extra = f"{f:.2f}Hz"

            vals, unit_label = pu.convert_to(p[idx, :], mode, df, enbw=c["ENBW"])

            if absorption and mode in ("SPL", "PSD"): vals = au.add_absorption(f, vals, distances)

            series.append((c, theta_plot, np.append(vals, vals[0]), f"{unit_label} at {f:.2f} [Hz]", 0))

    if series:
        vals = np.concatenate([vals_plot for _, _, vals_plot, _, _ in series])
        rmin, rmax = np.nanmin(vals), np.nanmax(vals)
        if mode == "OASPL" and rmax > 0: rmin = max(0., rmin)
        margin = abs(rmax-rmin) * 0.1 if rmax != rmin else abs(rmin) * 0.1
        if margin == 0: margin = 1.0 # Fallback

        # Plot (angles without data: gap in the curve, fill collapsed to the centre)
        linestyles = pu.get_line_styles()
        for c, theta_plot, vals_plot, label, b in series:
            if ws_mode:
                ax.plot(theta_plot, vals_plot, linestyles[b % len(linestyles)] if mode == "OASPL" else '-',
                        lw=c["WS_lw"], color=c["WS_color"], zorder=c["WS_zorder"])
            else:
                ax.plot(theta_plot, vals_plot, lw=2, color=colors[b], label=label)
                ax.fill(theta_plot, np.where(np.isfinite(vals_plot), vals_plot, rmin - margin), alpha=0.08, color=colors[b])

        ax.set_rmin(rmin - margin)

    ax.set_title(f"{tag}")

    if ws_mode:
        ws_handles = [Line2D([], [], color=c["WS_color"], lw=c["WS_lw"]) for c in polars]
        ws_labels  = [pu.get_ws_label(c) for c in polars]

        # OASPL: the line styles of the bands (with a single band there is nothing to tell apart)
        band_handles, band_labels = [], []
        if mode == "OASPL":
            labels = [label for c, _, _, label, b in series if c is polars[0]]
            bs     = [b for c, _, _, _, b in series if c is polars[0]]
            if len(labels) > 1:
                band_handles = [ax.plot([], [], linestyles[b % len(linestyles)], color='k')[0] for b in bs]
                band_labels  = labels
        elif series:
            ax.set_title(f"{tag}\n{series[0][3]}")

        pu.add_ws_legends(ax, ws_handles, ws_labels, band_handles, band_labels, style_title="Frequency band", offset=0.5)
    else:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)

    if ax_in is None:
        if ws_mode:
            fig.tight_layout()
        else:
            fig.tight_layout(rect=[0, 0.06, 1, 1])
        pu.set_figure_name(fig, "polar", mode, name_extra, "allws" if ws_mode else "")

    return fig, ax


# ---------- Cylindrical Plots ---------- #
def plot_cylinder(case            : dict = None,        # [-] Dictionary containing simulation data
                  mode            : str = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                  target_frequency: float = 0.88,       # [Hz] Target frequency for single-frequency modes
                  absorption      : bool = True,        # [-] Whether to apply absorption attenuation
                  filter_under    : float = None,       # [Hz] Lower frequency cutoff (f >= filter_under)
                  filter_over     : float = None,       # [Hz] Upper frequency cutoff (f <= filter_over)
                  figsize         : tuple = (12,5),     # [-] Figure size
                  cmap            : str = 'inferno',    # [-] Contour colormap
                  ax_in           : plt.Axes = None):   # [-] Optional external axes for subplots
    """
    Plots 2D unwrapped cylindrical surface maps for acoustic magnitudes (OASPL, SPL, SWL, etc.). Observers
    without data (NaN, inside the structure) are left as holes in the contours.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes (or array of Axes for multi-band OASPL)
    """

    if not case["has_cylinder"]:
        print("plot_cylinder(): no cylinder data, skipping")
        return None, None

    tag = case["Case_type"]

    # Extract data
    theta     = case["Theta_deg_cylinder"]   # 1D array of angles [deg]
    z         = case["Z_cylinder"]           # 1D array of height/depth [m]
    p         = case["P_cylinder"]           # Pressure grid
    Freqs     = case["Freqs"]
    df        = pu.get_bin_width(Freqs)      # Frequency bin width [Hz] for 'PSD'
    p_ref     = case["p_ref"]
    r         = case["R_cylinder"]
    observers = case["Obs_cylinder"]
    cx, cy    = case["Center_cylinder"]

    if not pu.is_valid(p, "plot_cylinder"):
            return None, None

    nf, ntheta, nz = len(Freqs), len(theta), len(z)

    # Reshape
    observers = observers.reshape((ntheta*nz,3))
    p         = p.reshape((nf, ntheta*nz))

    # Distance grid for absorption
    distances = np.full(ntheta*nz, r, dtype=float)

    # Meshgrid and differential of area
    if "dA_cylinder" in case:
        dA = case["dA_cylinder"]
    else:
        dtheta = 2*np.pi / ntheta
        dz     = (z.max() - z.min())/(nz-1) if nz > 1 else 1.0
        dA     = r * dz * dtheta

    THETA, Z = np.meshgrid(theta, z)

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
        vals = au.pressure_to_OASPL(p, p_ref, absorption, Freqs, distances, enbw=case["ENBW"])
        unit_label = r"OASPL [dB re 1$\mu$ Pa]"

    elif mode == 'SWL':
        print("plot_cylinder(): SWL not yet implemented")
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f = Freqs[idx]
        p = p[idx,:]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = pu.convert_to(p, mode, df, enbw=case["ENBW"])

        if absorption and mode in ('SPL', 'PSD'): vals = au.add_absorption(f, vals, distances)

    # Plot
    theta_plot = np.append(theta, 360.)
    vals = vals.reshape((ntheta, nz))
    vals_plot = np.append(vals, vals[0:1, :], axis=0)

    # Min range and their
    max_val = np.nanmax(vals_plot)
    min_val = max_val / 2.0
    if not min_val < max_val: min_val = np.nanmin(vals_plot)     # Non-dB modes may have a maximum <= 0
    custom_levels = np.linspace(min_val, max_val, 60)
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
        pu.set_figure_name(fig, "cylinder", mode, title_extra)

    return fig, ax


# ---------- Linear Plots ---------- #
def plot_distance_decay(case            : dict  = None,         # [-] Dictionary containing simulation data
                        mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                        target_frequency: float = 0.88,         # [Hz] Target frequency for single-frequency modes
                        absorption      : bool  = True,         # [-] Whether to apply absorption attenuation
                        filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                        filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                        turbine_data    : bool  = False,        # [-] Wheter to plot additional turbine data from GE 2025
                        figsize         : tuple = (10,5),       # [-] Figure size
                        ax_in           : plt.Axes = None):     # [-] Optional external axes for subplots
    """
    Plots acoustic magnitude decay over distance (1D line plot). Distances without data (NaN, inside the
    structure) are gaps in the curve; the 1/r and 1/r^2 references start at the first distance with data.
    With --ws one curve per wind speed (color = wind speed, legend below the axes); the references (black)
    start at the curve of the case itself and the curves are not filled.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    if not case["has_decay"]:
        print("plot_distance_decay(): no decay data, skipping")
        return None, None

    tag    = case["Case_type"]
    do_log = case["Logspace_decay"]

    curves, ws_mode = get_ws_curves(case, "decay", "plot_distance_decay", mode=mode, target_frequency=target_frequency,
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
            color = c["WS_color"]
            handles.append(draw(r, vals, lw=c["WS_lw"], color=color, zorder=c["WS_zorder"])[0])
        else:
            color = pu.get_case_color(case)
            draw(r, vals, lw=2, color=color)
        labels.append(pu.get_ws_label(c))

    # References from the first distance with data of the case itself (the first curve if it was left out)
    c, (r, vals, unit_label, title_extra) = next((curve for curve in curves if curve[0] is case), curves[0])

    i0 = np.flatnonzero(np.isfinite(vals))[0]      # First distance with data: origin of the 1/r references
    ref_color  = 'k' if ws_mode else color
    ref_vals_1 = vals[i0] - 10 * np.log10(np.maximum(r, np.finfo(float).eps) / r[i0])
    ref_vals_2 = vals[i0] - 20 * np.log10(np.maximum(r, np.finfo(float).eps) / r[i0])
    draw(r, ref_vals_1, '--', color=ref_color, alpha=0.4, lw=1.0, label=r"$1/r$")
    draw(r, ref_vals_2, ':', color=ref_color, alpha=0.4, lw=1.0, label=r"$1/r^2$")

    min_val = min(ref_vals_2.min(), *(np.nanmin(vals_c) for _, (_, vals_c, _, _) in curves))

    # Plot data from Ge2025
    if turbine_data:
        TURBINE_STYLES = pu.get_turbine_styles()
        all_turbine_x, all_turbine_y = [], []
        for (data, label, marker, color) in TURBINE_STYLES:
            ax.scatter(data[:, 0], data[:, 1],
                    marker=marker, color=color, s=40,
                    alpha=1.0, zorder=5, label=label, edgecolors='k', linewidths=0.3)

            all_turbine_x.append(data[:,0])
            all_turbine_y.append(data[:,1])
            min_val = min(min_val, data[:,1].min()*0.95)

        # Regresion
        tx, ty = np.concatenate(all_turbine_x), np.concatenate(all_turbine_y)

        log_tx = np.log10(tx)
        a, b   = np.polyfit(log_tx, ty, 1)
        x_fit  = np.logspace(np.log10(tx.min()), np.log10(tx.max()), len(r))
        y_fit  = a * np.log10(x_fit) + b

        ax.semilogx(x_fit, y_fit, lw=2.0, color="lightcoral", ls="--", alpha=0.5, label="Regresion -20.4 dB/decade")

    ax.grid(True, which='both', ls='--', alpha=0.7)
    if not ws_mode:
        ax.fill_between(r, vals, min_val, alpha=0.1, color=color)

    ax.set_xlabel("Distance, r [m]")
    ax.set_ylabel(unit_label)
    ax.set_title(f"{tag} - {mode.upper()} {title_extra}")
    ax.set_xlim(min(r_c.min() for _, (r_c, _, _, _) in curves), max(r_c.max() for _, (r_c, _, _, _) in curves))
    ax.set_ylim(bottom=min_val)

    if ws_mode:
        pu.add_ws_legends(ax, handles, labels)

    if ax_in is None:
        fig.tight_layout()
        pu.set_figure_name(fig, "decay", mode.upper(), title_extra, "allws" if ws_mode else "")

    return fig, ax

def plot_line(case            : dict  = None,         # [-] Dictionary containing simulation data
              mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
              target_frequency: float = 0.88,         # [Hz] Target frequency for single-frequency modes
              absorption      : bool  = True,         # [-] Whether to apply absorption attenuation
              filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
              filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
              figsize         : tuple = (10,5),       # [-] Figure size
              ax_in           : plt.Axes = None):     # [-] Optional external axes for subplots
    """
    Plots acoustic magnitude along an arbitrary straight line between two points. Points without data (NaN,
    inside the structure) are gaps in the curve. With --ws one curve per wind speed (color = wind speed,
    legend below the axes), not filled.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    if not case["has_line"]:
        print("plot_line(): no line data, skipping")
        return None, None

    tag    = case["Case_type"]
    p1     = case["P1_line"]
    p2     = case["P2_line"]
    do_log = case["Logspace_line"]

    curves, ws_mode = get_ws_curves(case, "line", "plot_line", mode=mode, target_frequency=target_frequency,
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
            color = c["WS_color"]
            handles.append(draw(r, vals, lw=c["WS_lw"], color=color, zorder=c["WS_zorder"])[0])
        else:
            color = pu.get_case_color(case)
            draw(r, vals, lw=2, color=color)
        labels.append(pu.get_ws_label(c))

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
        pu.add_ws_legends(ax, handles, labels)

    if ax_in is None:
        fig.tight_layout()
        pu.set_figure_name(fig, "line", mode.upper(), title_extra, "allws" if ws_mode else "")

    return fig, ax


# ---------- Slice Plots ---------- #
def plot_sliceXY(case            : dict  = None,         # [-] Dictionary containing simulation data
                 mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                 target_frequency: float = 0.88,         # [Hz] Target frequency for single-frequency modes
                 absorption      : bool  = True,         # [-] Whether to apply absorption attenuation
                 filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                 filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                 structure       : bool = True,          # [-] Wheter to plot structure on slice
                 figsize         : tuple = (7,7),        # [-] Figure size
                 cmap            : str   = 'inferno',    # [-] Contour colormap
                 ax_in           : plt.Axes = None):     # [-] Optional external axes for subplots
    """
    Plots a 2D spatial slice in the XY plane of acoustic magnitudes. Observers without data (NaN, inside the
    structure) are left as holes in the contours.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    if not case["has_slicexy"]:
        print("plot_sliceXY(): no XY slice data, skipping")
        return None, None

    tag = case["Case_type"]

    # Extract data
    x         = case["X_slicexy"]           # 1D array of X coordinates [m]
    y         = case["Y_slicexy"]           # 1D array of Y coordinates [m]
    z_slice   = case["Z_slicexy"]           # Z level of the slice [m]
    p         = case["P_slicexy"]           # Pressure grid
    Freqs     = case["Freqs"]
    df        = pu.get_bin_width(Freqs)     # Frequency bin width [Hz] for 'PSD'
    observers = case["Obs_slicexy"]
    p_ref     = case["p_ref"]
    cx, cy    = case["Center_slicexy"]

    if not pu.is_valid(p, "plot_sliceXY"):
        return None, None

    Nfreqs, nx, ny = len(Freqs), len(x), len(y)
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
    NFreqs = len(Freqs)

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
        vals = au.pressure_to_OASPL(p, p_ref, absorption, Freqs, distances, enbw=case["ENBW"])
        unit_label = r"OASPL [dB re 1$\mu$ Pa]"
    elif mode == 'SWL':
        print("WARNING PLot_sliceXY(): SWL not yet implemented")
        pass
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f   = Freqs[idx]
        p_f = p[idx, :]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = pu.convert_to(p_f, mode, df, enbw=case["ENBW"])
        if absorption and mode in ('SPL', 'PSD'): vals = au.add_absorption(f, vals, distances)

    # Plot (NaN: holes)
    contour_set = ax.contourf(X, Y, vals, levels=60, cmap=cmap)
    cbar = fig.colorbar(contour_set, ax=ax, orientation='horizontal', pad = 0.1)
    cbar.set_label(unit_label)
    if abs(x[-1]-x[0]) == abs(y[-1]-y[0]): ax.set_aspect('equal')

    # Structure
    if structure:
        try:
            nodes = case["Structure_nodes"]
            for member in nodes:
                ax.plot(member[:,0], member[:,1], '-o', color='black', markersize=2, linewidth=1.0, zorder=10)
        except Exception as e:
            print(f"plot_sliceXY() WARNING: Could not plot structure. ({type(e).__name__}: {e})")


    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title(f"{tag} - {mode} {title_extra} - Depth = {z_slice:.1f} [m]")
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())

    if ax_in is None:
        fig.tight_layout()
        pu.set_figure_name(fig, "sliceXY", mode, title_extra)

    return fig, ax

def plot_sliceXZ(case            : dict  = None,         # [-] Dictionary containing simulation data
                 mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                 target_frequency: float = 0.88,         # [Hz] Target frequency for single-frequency modes
                 absorption      : bool  = True,         # [-] Whether to apply absorption attenuation
                 filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                 filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                 structure       : bool  = True,         # [-] Wheter to plot structure on slice
                 figsize         : tuple = (7,7),        # [-] Figure size
                 cmap            : str   = 'inferno',    # [-] Contour colormap
                 ax_in           : plt.Axes = None):     # [-] Optional external axes for subplots
    """
    Plots a 2D spatial slice in the XZ vertical plane of acoustic magnitudes. Observers without data (NaN,
    inside the structure) are left as holes in the contours.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    if not case["has_slicexz"]:
        print("plot_sliceXZ(): no XZ slice data, skipping")
        return None, None

    tag = case["Case_type"]

    # Extract data
    x         = case["X_slicexz"]           # 1D array of X coordinates [m]
    z         = case["Z_slicexz"]           # 1D array of Z coordinates [m]
    y_slice   = case["Y_slicexz"]           # Y level of the slice [m]
    p         = case["P_slicexz"]           # Pressure grid (Nf, nx, nz)
    Freqs     = case["Freqs"]
    df        = pu.get_bin_width(Freqs)     # Frequency bin width [Hz] for 'PSD'
    observers = case["Obs_slicexz"]         # Observer grid (nx, nz, 3)
    p_ref     = case["p_ref"]

    if not pu.is_valid(p, "plot_sliceXZ"):
        return None, None

    cx, cz = np.mean(x), np.mean(z)
    Nfreqs, nx, nz = len(Freqs), len(x), len(z)
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
    NFreqs = len(Freqs)

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
        vals = au.pressure_to_OASPL(p, p_ref, absorption, Freqs, distances, enbw=case["ENBW"])
        unit_label = r"OASPL [dB re 1$\mu$ Pa]"
    elif mode == 'SWL':
        print("WARNING plot_sliceXZ(): SWL not yet implemented")
        pass
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f   = Freqs[idx]
        p_f = p[idx, :, :]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = pu.convert_to(p_f, mode, df, enbw=case["ENBW"])
        if absorption and mode in ('SPL', 'PSD'): vals = au.add_absorption(f, vals, distances)

    # Plot
    # Min range and their
    max_val = np.nanmax(vals)
    min_val = max_val / 2.0
    if not min_val < max_val: min_val = np.nanmin(vals)     # Non-dB modes may have a maximum <= 0
    custom_levels = np.linspace(min_val, max_val, 60)
    line_levels = custom_levels[::2]

    # NaN: holes
    contour_set = ax.contourf(X, Z, vals, levels=custom_levels, cmap=cmap, extend='min')
    cbar = fig.colorbar(contour_set, ax=ax, orientation='horizontal', pad=0.1)
    cbar.set_label(unit_label)
    if abs(x[-1]-x[0]) == abs(z[-1]-z[0]): ax.set_aspect('equal')

    z_top = z.max()
    if structure:
            try:
                nodes = case["Structure_nodes"]
                for member in nodes:
                    ax.plot(member[:,0], member[:,2], '-o', color='black', markersize=2, linewidth=1.0, zorder=10)
                z_top = np.nanmax(nodes[:,:,2])
            except Exception as e:
                print(f"plot_sliceXZ() WARNING: Could not plot structure. ({type(e).__name__}: {e})")

    ax.set_xlabel('x [m]')
    ax.set_ylabel('z [m]')
    ax.set_title(f"{tag} - {mode} {title_extra} - Y = {y_slice:.1f} [m]")
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(z.min(), z_top)

    if ax_in is None:
        fig.tight_layout()
        pu.set_figure_name(fig, "sliceXZ", mode, title_extra)

    return fig, ax

def plot_sliceVertical(case            : dict  = None,         # [-] Dictionary containing simulation data
                       mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                       target_frequency: float = 0.88,         # [Hz] Target frequency for single-frequency modes
                       absorption      : bool  = True,         # [-] Whether to apply absorption attenuation
                       filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                       filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                       structure       : bool  = True,         # [-] Wheter to plot structure on slice
                       figsize         : tuple = (7,7),        # [-] Figure size
                       cmap            : str   = 'inferno',    # [-] Contour colormap
                       ax_in           : plt.Axes = None):     # [-] Optional external axes for subplots
    """
    Plots a 2D spatial slice in a vertical plane of arbitrary azimuth of acoustic magnitudes. The plane normal
    is (cos az, sin az) and its horizontal axis u runs along (-sin az, cos az) from Center_sliceV
    (AcousticSolver.run_sliceVertical). The plane usually crosses the structure: observers without data (NaN,
    inside the structure) are left as holes in the contours.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    if not case["has_slicevertical"]:
        print("plot_sliceVertical(): no vertical slice data, skipping")
        return None, None

    tag = case["Case_type"]

    # Extract data
    u         = case["U_sliceV"]            # 1D array of along-slice coordinates [m]
    z         = case["Z_sliceV"]            # 1D array of Z coordinates [m]
    p         = case["P_sliceV"]            # Pressure grid (Nf, nu, nz)
    observers = case["Coords_sliceV"]       # Observer grid (nu, nz, 3)
    azimuth   = case["Azimuth_sliceV"]      # Azimuth angle of the plane normal [deg]
    center    = case["Center_sliceV"]       # Center point [cx, cy]
    Freqs     = case["Freqs"]
    df        = pu.get_bin_width(Freqs)     # Frequency bin width [Hz] for 'PSD'
    p_ref     = case["p_ref"]

    if not pu.is_valid(p, "plot_sliceVertical"):
        return None, None

    cx, cy = center[0], center[1]
    Nfreqs, nu, nz = len(Freqs), len(u), len(z)
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
    NFreqs = len(Freqs)

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
        vals = au.pressure_to_OASPL(p, p_ref, absorption, Freqs, distances, enbw=case["ENBW"])
        unit_label = r"OASPL [dB re 1$\mu$ Pa]"
    elif mode == 'SWL':
        print("WARNING plot_sliceVertical(): SWL not yet implemented")
        pass
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f   = Freqs[idx]
        p_f = p[idx, :, :]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = pu.convert_to(p_f, mode, df, enbw=case["ENBW"])
        if absorption and mode in ('SPL', 'PSD'): vals = au.add_absorption(f, vals, distances)

    # Plot
    # Min range and their
    max_val = np.nanmax(vals)
    min_val = max_val / 2.0
    if not min_val < max_val: min_val = np.nanmin(vals)     # Non-dB modes may have a maximum <= 0
    custom_levels = np.linspace(min_val, max_val, 60)

    # NaN: holes
    contour_set = ax.contourf(U, Z, vals, levels=custom_levels, cmap=cmap, extend='min')
    cbar = fig.colorbar(contour_set, ax=ax, orientation='horizontal', pad=0.1)
    cbar.set_label(unit_label)
    if abs(u[-1]-u[0]) == abs(z[-1]-z[0]): ax.set_aspect('equal')

    # Structure projected onto the plane: u = (x - cx, y - cy) . (-sin az, cos az)
    z_top = z.max()
    if structure:
        try:
            theta_rad = np.deg2rad(azimuth)
            cos_t = np.cos(theta_rad)
            sin_t = np.sin(theta_rad)

            nodes = case["Structure_nodes"]
            for member in nodes:
                u_proj = -(member[:,0] - cx) * sin_t + (member[:,1] - cy) * cos_t
                ax.plot(u_proj, member[:,2], '-o', color='black', markersize=2, linewidth=1.0, zorder=10)
            z_top = np.nanmax(nodes[:,:,2])
        except Exception as e:
            print(f"plot_sliceVertical() WARNING: Could not plot structure. ({type(e).__name__}: {e})")

    ax.set_xlabel('u (Distance along slice) [m]')
    ax.set_ylabel('z [m]')
    ax.set_title(f"{tag} - {mode} {title_extra} - Azimuth = {azimuth:.1f} [deg]")
    ax.set_xlim(u.min(), u.max())
    ax.set_ylim(z.min(), z_top)

    if ax_in is None:
        fig.tight_layout()
        pu.set_figure_name(fig, "sliceVertical", mode, title_extra)

    return fig, ax


# ---------- Compute metrics ---------- #
def compute_sphere_metrics(case            : dict  = None,         # [-] Dictionary containing simulation data
                           absorption      : bool  = True,         # [-] Whether to apply absorption attenuation for OASPL
                           filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                           filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                           rho_water       : float = 1025.0,       # [kg/m³] Water density
                           c_water         : float = 1500.0,       # [m/s] Speed of sound in water
                           verbose         : bool  = True):        # [-] Print metrics to console
    """
    Computes and optionally prints radiation metrics from a spherical pressure field. Observers without data
    (NaN, inside the structure) are left out and the excluded fraction of the sphere area is reported.

    Returns
    -------
    metrics : dict
        Dictionary containing all computed acoustic metrics (Power, SL, etc.)
    """

    if not case["has_sphere"]:
        print("compute_sphere_metrics(): no sphere data found in case, skipping.")
        return None

    rho_wat = case.get("rho_wat", rho_water)
    c_wat   = case.get("c_wat", c_water)
    tag     = case["Case_type"]

    # Extract data
    p       = case["P_sphere"]
    Freqs   = case["Freqs"]
    r       = case["R_sphere"]
    dA      = case["dA_sphere"]
    p_ref   = case["p_ref"]
    n_theta = case["N_theta"]
    nz      = case["Nz_sphere"]

    if not pu.is_valid(p, "compute_sphere_metrics"):
        return None

    # Wetted area saved by the Fortran code (save_specific_params of every turbine type)
    A_wet = case.get("wet_area")
    if A_wet is None:
        A_wet = 848.12 if "Monopile" in tag else 9359.07
        print(f"WARNING compute_sphere_metrics(): wet_area not found in the file, using the hardcoded {A_wet} m²")

    # Filter frequencies
    Nfreqs = len(Freqs)
    freq_mask = np.ones(Nfreqs, dtype=bool)
    if filter_under is not None:
        freq_mask &= (Freqs >= filter_under)
    if filter_over is not None:
        freq_mask &= (Freqs <= filter_over)
    if not np.any(freq_mask):
        print("compute_sphere_metrics(): no frequencies match the specified filter bounds!")
        return None

    p, Freqs  = p[freq_mask], Freqs[freq_mask]
    distances = np.full((n_theta, nz), r)

    # --- 1. Compute OASPL and Directional SL --- #
    oaspl = au.pressure_to_OASPL(p, p_ref, absorption, Freqs, distances, enbw=case["ENBW"])
    sl    = oaspl + 20 * np.log10(r)

    # --- 2. Compute Radiated Power (W) --- #
    # Observers without data (NaN) are left out of the surface integral (equal dA for every observer)
    p_sq  = np.sum(np.abs(p)**2, axis=0)
    valid = np.isfinite(p_sq)
    excluded = 1.0 - np.count_nonzero(valid) / valid.size
    I = p_sq / (2.*rho_wat*c_wat)
    W = np.sum(I[valid]*dA)
    W_dB = 10* np.log10(W/1e-12)

    # --- 3. Compute Power per Wet Area --- #
    W_per_wet = W / A_wet
    W_per_wet_dB = 10.0 * np.log10(W_per_wet / 1e-12)

    # --- 4. Equivalent Omnidirectional SL --- #
    W_ref_pressure = (4.0 * np.pi * p_ref**2) / (rho_wat * c_wat)
    SL_omni = 10.0 * np.log10(W / W_ref_pressure)

    # --- 5. Directional SL Statistics --- #
    SL_max = np.nanmax(sl)
    SL_min = np.nanmin(sl)
    SL_mean_power = 10.0 * np.log10(np.nanmean(10.0**(sl / 10.0))) # Energetic mean

    # --- 6. Mean SPL on the sphere --- #
    # p is a peak amplitude (p_rms^2 = |p|^2 / 2), as in the OASPL of pressure_to_SPL
    SPL_sphere = 10.0 * np.log10(np.mean(p_sq[valid]) / (2.0 * p_ref**2))

    # --- Package results --- #
    metrics = {
        "OASPL_sphere": oaspl,
        "SL_sphere": sl,
        "W": W,
        "W_dB": W_dB,
        "W_per_wet": W_per_wet,
        "W_per_wet_dB": W_per_wet_dB,
        "SL_omni": SL_omni,
        "SL_max": SL_max,
        "SL_min": SL_min,
        "SL_mean_power": SL_mean_power,
        "SPL_sphere_mean": SPL_sphere,
        "Excluded_area_fraction": excluded
    }

    # --- Print formatting --- #
    if verbose:
        print("\n" + "=" * 80)
        print(f"  RADIATION METRICS (SPHERE r = {r:.0f} m, free-field)")

        if filter_under is not None or filter_over is not None:
            lo = f"{filter_under:.1f}" if filter_under is not None else "0"
            hi = f"{filter_over:.1f}"  if filter_over is not None else "∞"
            print(f"  Frequency band: {lo} - {hi} Hz")
        print("=" * 80)

        label_w = 42
        print(f"\n[Case: {tag}]")
        print(f"  {'Radius of sphere':<{label_w}}: {r:.1f} m")
        print(f"  {'Wet area of structure':<{label_w}}: {A_wet:.1f} m²")
        if excluded > 0:
            print(f"  {'Excluded area (observers without data)':<{label_w}}: {100*excluded:.2f} % "
                  f"({valid.size - np.count_nonzero(valid)} of {valid.size} observers)")
        print(f"  {'Total radiated power (W)':<{label_w}}: {W:.4e} W")
        print(f"  {'Total radiated power (dB)':<{label_w}}: {W_dB:.1f} dB re 1 pW")
        print(f"  {'Power / wet area (W)':<{label_w}}: {W_per_wet:.4e} W/m²")
        print(f"  {'Power / wet area (dB)':<{label_w}}: {W_per_wet_dB:.1f} dB re 1 pW/m²")
        print(f"  {'Equivalent omnidirectional SL':<{label_w}}: {SL_omni:.1f} dB re 1 µPa @ 1 m")
        print(f"  {'Directional SL MAX':<{label_w}}: {SL_max:.1f} dB re 1 µPa @ 1 m")
        print(f"  {'Directional SL MIN':<{label_w}}: {SL_min:.1f} dB re 1 µPa @ 1 m")
        print(f"  {'Directional SL mean (energ.)':<{label_w}}: {SL_mean_power:.1f} dB re 1 µPa @ 1 m")
        print(f"  {'Mean SPL on sphere':<{label_w}}: {SPL_sphere:.1f} dB re 1 µPa")

    return metrics

def compute_cylinder_metrics(case            : dict  = None,         # [-] Dictionary containing simulation data
                             filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                             filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                             rho_water       : float = 1025.0,       # [kg/m³] Water density
                             c_water         : float = 1500.0,       # [m/s] Speed of sound in water
                             verbose         : bool  = True):        # [-] Print metrics to console
    """
    Computes and optionally prints propagation metrics from a cylindrical pressure field. Observers without
    data (NaN, inside the structure) are left out and the excluded fraction of the lateral area is reported.

    Returns
    -------
    metrics : dict
        Dictionary containing all computed acoustic metrics (Power, SPL, TL, etc.)
    """

    if not case["has_cylinder"]:
        print("compute_cylinder_metrics(): no cylinder data, skipping")
        return None

    rho_wat = case.get("rho_wat", rho_water)
    c_wat   = case.get("c_wat", c_water)
    tag     = case["Case_type"]

    # Extract data
    p       = case["P_cylinder"]
    Freqs   = case["Freqs"]
    r       = case["R_cylinder"]
    z       = case["Z_cylinder"]
    theta   = case["Theta_deg_cylinder"]
    p_ref   = case["p_ref"]

    if not pu.is_valid(p, "compute_cylinder_metrics"):
            return None

    n_theta, nz = len(theta), len(z)

    # Ensure theta is in radians for integration
    theta_rad = theta if np.max(theta) <= 2*np.pi + 0.1 else np.deg2rad(theta)

    # Filter frequencies
    Nfreqs = len(Freqs)
    freq_mask = np.ones(Nfreqs, dtype=bool)
    if filter_under is not None:
        freq_mask &= (Freqs >= filter_under)
    if filter_over is not None:
        freq_mask &= (Freqs <= filter_over)
    if not np.any(freq_mask):
        print("compute_cylinder_metrics(): no frequencies match the specified filter bounds!")
        return None

    p, Freqs = p[freq_mask], Freqs[freq_mask]

    # --- 1. Intensity and Integration Setup --- #
    p_sq = np.sum(np.abs(p)**2, axis=0)      # Sum energy over frequencies: shape (n_theta, nz)
    I = p_sq / (2.0 * rho_wat * c_wat)

    A = 2 * np.pi * r * np.max(np.abs(z))    # Lateral area of the cylinder

    # Observers without data (NaN) are left out of the surface integral
    valid = np.isfinite(I)
    I     = np.where(valid, I, 0.0)

    # Close the cylinder: the azimuthal integral covers the full 360 deg (last point = first point)
    theta_closed = np.append(theta_rad, theta_rad[0] + 2*np.pi)
    I_closed     = np.append(I, I[:1, :], axis=0)
    w_closed     = np.append(valid, valid[:1, :], axis=0).astype(float)

    # --- 2. Compute Radiated Power (W) --- #
    I_vs_z = trapezoid(I_closed * r, theta_closed, axis=0)
    W = np.abs(trapezoid(I_vs_z, z))          # Absolute to handle z-axis integration direction safely
    W_dB = 10 * np.log10(W / 1e-12)

    # Excluded fraction of the lateral area, with the same integration weights
    A_valid  = np.abs(trapezoid(trapezoid(w_closed * r, theta_closed, axis=0), z))
    A_total  = np.abs(trapezoid(trapezoid(np.ones_like(w_closed) * r, theta_closed, axis=0), z))
    excluded = 1.0 - A_valid / A_total

    # --- 3. Compute Power per Cylinder Area --- #
    W_per_area = W / A
    W_per_area_dB = 10 * np.log10(W_per_area / 1e-12)

    # --- 4. Mean SPL on the cylinder --- #
    # p is a peak amplitude (p_rms^2 = |p|^2 / 2), as in the OASPL of pressure_to_SPL
    p2_avg = np.mean(p_sq[valid])
    SPL_avg = 10 * np.log10(p2_avg / (2.0 * p_ref**2))

    # --- 5. Apparent Omnidirectional SL (from Far Field) --- #
    W_ref_press = (4.0 * np.pi * p_ref**2) / (rho_wat * c_wat)
    SL_omni_far = 10.0 * np.log10(W / W_ref_press)

    # --- 6. Comparative Metrics (if Sphere data exists) --- #
    TL = None
    power_loss = None
    TL_str = ""
    power_loss_str = ""

    # --- Package results ---
    metrics = {
        "W_cyl": W,
        "W_cyl_dB": W_dB,
        "W_cyl_per_area": W_per_area,
        "W_cyl_per_area_dB": W_per_area_dB,
        "SPL_cyl_mean": SPL_avg,
        "SL_omni_far": SL_omni_far,
        "TL_cyl": TL,
        "Power_loss_cyl": power_loss,
        "Excluded_area_fraction": excluded
    }

    # --- Print formatting ---
    if verbose:
        print("\n" + "=" * 80)
        print(f"  PROPAGATION METRICS (CYLINDER r = {r:.0f} m)")

        if filter_under is not None or filter_over is not None:
            lo = f"{filter_under:.1f}" if filter_under is not None else "0"
            hi = f"{filter_over:.1f}"  if filter_over is not None else "∞"
            print(f"  Frequency band: {lo} - {hi} Hz")
        print("=" * 80)

        label_w = 44
        print(f"\n[Case: {tag}]")
        print(f"  {'Radius of cylinder':<{label_w}}: {r:.1f} m")
        print(f"  {'Height of cylinder (H)':<{label_w}}: {np.max(np.abs(z)):.1f} m")
        print(f"  {'Lateral area of cylinder':<{label_w}}: {A:.1f} m²")
        if excluded > 0:
            print(f"  {'Excluded area (observers without data)':<{label_w}}: {100*excluded:.2f} % "
                  f"({valid.size - np.count_nonzero(valid)} of {valid.size} observers)")
        print(f"  {'Total power flow (W)':<{label_w}}: {W:.4e} W")
        print(f"  {'Total power flow (dB)':<{label_w}}: {W_dB:.1f} dB re 1 pW")
        print(f"  {'Power / cylinder area (W)':<{label_w}}: {W_per_area:.4e} W/m²")
        print(f"  {'Power / cylinder area (dB)':<{label_w}}: {W_per_area_dB:.1f} dB re 1 pW/m²")
        print(f"  {'Mean SPL on cylinder':<{label_w}}: {SPL_avg:.1f} dB re 1 µPa")
        print(f"  {'Apparent omnid. SL (from far-field)':<{label_w}}: {SL_omni_far:.1f} dB re 1 µPa @ 1 m")

        if TL_str:
            print(f"  {'Transmission Loss (TL)':<{label_w}}: {TL_str}")
        if power_loss_str:
            print(f"  {'Power loss (W_sphere / W_cyl)':<{label_w}}: {power_loss_str}")

    return metrics

def compute_point_source(case            : dict  = None,         # [-] Dictionary containing simulation data
                         filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                         filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                         enbw            : float = None,         # [bins] Equivalent noise bandwidth of the FFT window (None: case["ENBW"], from the file)
                         smooth_fraction : int   = 12,           # [-] Plots only: narrowband curves smoothed over 1/N octave (None: raw)
                         rho_water       : float = 1025.0,       # [kg/m³] Water density (if not in the file)
                         c_water         : float = 1500.0,       # [m/s] Speed of sound in water (if not in the file)
                         verbose         : bool  = True,         # [-] Print totals and the 1/3 octave table
                         do_plot         : bool  = True,         # [-] Plot SL spectra, radiation diagnostics and PL
                         figsize         : tuple = (10, 5)):     # [-] Size of the SL figure
    """
    Reduces the turbine to an equivalent point source: the free-field monopole that radiates the same power,
    SL = 10 log10(rho c W / (4 pi p_ref^2 r_ref^2)) [dB re 1 uPa^2 m^2], with r_ref = 1 m. Two powers are used:

    - Free field (sphere, Image Method with N_images = 0): power of the source itself -> SL_free.
    - In situ (cylinder, waveguide): power that the turbine injects in the actual water column -> SL_insitu.
      With lossless boundaries every watt radiated by the nodes crosses the cylinder, whatever its radius.

    Both powers come from W = sum |p|^2 / (2 rho c) dA (far-field plane-wave intensity). Known biases, quantified
    with exact dipole power on the SAITEC geometry (10/2026): the sphere overestimates when it hugs the
    structure (geometric factor, reported as an estimate) and at low kr (near field); the cylinder
    overestimates by ~0.5 to 4.7 dB because the modes cross it obliquely (only the radial intensity carries
    power). Exact values need the active intensity 1/2 Re{p v_r*} (e.g. cylinders at r +- 5 cm).

    The FFT of the Fortran code (compute_rfft) is a Hann window scaled for amplitude, so |p|^2/2 per bin is the
    power of a tone on that bin while sums over bins are ENBW times too large: densities and band/total
    powers are divided by enbw, the per-bin tone level is not. By default enbw is the one of the file
    (case["ENBW"], see plot_utils.read_enbw): 1.5 for forces from compute_rfft, 1 for forces without a window.

    Also computed:
    - Directivity index per 1/3 octave band on the sphere, DI = 10 log10(max I / mean I).
    - Multipole fit on the sphere around its centre, p = A h0(kr) + (D . n) h1(kr): monopole and dipole parts
      of the power (the dipole one with its exact near-field factor 1 + 1/(kr)^2), the fraction of |p|^2 they
      explain (how well the turbine reduces to a single point) and the dipole axis per band.
    - Propagation loss PL = SL - RL per 1/3 octave band at the spectrum observers (and their energetic mean),
      to turn a measured band level into a source level: SL_measured = L_measured + PL.

    Returns
    -------
    metrics : dict
        Narrowband arrays on "Freqs" (W_*_bin [W], SL_*_bin [dB re 1 uPa^2 m^2], SL_*_density
        [dB re 1 uPa^2 m^2/Hz], kr, ratio_insitu_free [dB], multipole fractions), 1/3 octave arrays on "fc"
        (SL_*_band, DI_band, PL_band, dipole axis, ...) and broadband totals (SL_*_total, W_*_total).
    """

    has_sphere, has_cyl = case.get("has_sphere", False), case.get("has_cylinder", False)
    if not (has_sphere or has_cyl):
        print("compute_point_source(): no sphere nor cylinder data, skipping")
        return None

    rho_wat = case.get("rho_wat", rho_water)
    c_wat   = case.get("c_wat", c_water)
    p_ref   = case["p_ref"]
    tag     = case["Case_type"]
    Freqs   = case["Freqs"]
    enbw    = case["ENBW"] if enbw is None else enbw
    df      = pu.get_bin_width(Freqs)     # [Hz] Bin width, from the full grid
    if df is None:
        print("compute_point_source(): needs at least two frequencies, skipping")
        return None

    # --- 0. Frequency band --- #
    mask = np.ones(len(Freqs), dtype=bool)
    if filter_under is not None:
        mask &= (Freqs >= filter_under)
    if filter_over is not None:
        mask &= (Freqs <= filter_over)
    if not np.any(mask):
        print("compute_point_source(): no frequencies match the specified filter bounds!")
        return None

    f  = Freqs[mask]
    k  = 2.0 * np.pi * f / c_wat
    fc = au.generate_iso_third_octave(f.min(), f.max())[1]

    def to_SL(W):                       # [W] -> [dB re 1 uPa^2 m^2]: free-field monopole of the same power at 1 m
        with np.errstate(divide='ignore'):
            return 10.0 * np.log10(rho_wat * c_wat * W / (4.0 * np.pi * p_ref**2))

    def band_sum(values):               # Sum over every 1/3 octave band (same bands and edges as to_third_octave)
        return au.to_third_octave(f, values, fmin=f.min(), fmax=f.max(), domain='linear', linear_mode='power')[1]

    bw      = fc * (2**(1/6) - 2**(-1/6))     # [Hz] 1/3 octave band widths
    metrics = {"Freqs": f, "fc": fc, "df": df, "enbw": enbw}

    # --- 1. Free field: sphere --- #
    if has_sphere:
        p_s    = case["P_sphere"][mask].reshape(len(f), -1)              # (Nf, Nobs), peak amplitudes
        obs    = case["Obs_sphere"].reshape(-1, 3)
        center = np.asarray(case["Center_sphere"], dtype=float)
        r_s    = case["R_sphere"]
        dA_s   = case["dA_sphere"]
        valid  = np.all(np.isfinite(p_s), axis=0)
        p_s    = p_s[:, valid]

        # Radiated power of every bin [W], tone convention (not divided by the ENBW)
        W_free_bin = np.sum(np.abs(p_s)**2, axis=1) / (2.0 * rho_wat * c_wat) * dA_s

        # Validity: near field (kr) and how tightly the sphere hugs the wet nodes (geometric bias estimate: monopole at a)
        nodes = case["Structure_nodes"].reshape(-1, 3)
        nodes = nodes[np.all(np.isfinite(nodes), axis=1) & (nodes[:, 2] < 0.0)]
        a_max = np.max(np.linalg.norm(nodes - center, axis=1)) if len(nodes) else np.nan
        ratio = a_max / r_s
        geo_bias_dB = 10*np.log10((1/(2*ratio)) * np.log((1+ratio)/(1-ratio))) if 0 < ratio < 1 else np.nan

        # Multipole fit p = c0 + c . n (angular basis, constant over frequency since r is fixed)
        n_hat  = (obs[valid] - center) / r_s
        B      = np.column_stack([np.ones(len(n_hat)), n_hat])            # (Nobs, 4)
        coef   = np.linalg.lstsq(B, p_s.T, rcond=None)[0]                 # (4, Nf)
        G      = B.T @ B
        p2_fit = np.real(np.einsum('ik,ij,jk->k', coef.conj(), G, coef))  # ||B c||^2 per bin
        p2_all = np.sum(np.abs(p_s)**2, axis=1)
        kr     = k * r_s
        W_mono_bin = 4.0*np.pi * np.abs(coef[0])**2 * r_s**2 / (2.0*rho_wat*c_wat)
        W_dip_bin  = (4.0*np.pi/3.0) * np.sum(np.abs(coef[1:])**2, axis=0) * r_s**2 / (2.0*rho_wat*c_wat) / (1.0 + 1.0/kr**2)

        # Directivity per band (band-summed intensity pattern on the sphere)
        I_band  = au.to_third_octave(f, np.abs(p_s)**2, fmin=f.min(), fmax=f.max(), domain='linear', linear_mode='power')[1]
        DI_band = 10.0 * np.log10(np.max(I_band, axis=1) / np.mean(I_band, axis=1))

        # Dipole axis per band: principal direction of Re(sum D D^H), D = c / h1 (only its direction matters)
        axis_elev, axis_azim = np.full(len(fc), np.nan), np.full(len(fc), np.nan)
        for b, f_c in enumerate(fc):
            in_b = (f >= f_c / 2**(1/6)) & (f < f_c * 2**(1/6))
            if np.any(in_b):
                M = np.real(coef[1:, in_b] @ coef[1:, in_b].conj().T)
                u = np.linalg.eigh(M)[1][:, -1]
                u = u if u[2] >= 0 else -u
                axis_elev[b] = np.degrees(np.arcsin(np.clip(u[2], -1, 1)))     # [deg] from the horizontal plane
                axis_azim[b] = np.degrees(np.arctan2(u[1], u[0])) % 180.0       # [deg] from +x (axis, so modulo 180)

        W_free_band = band_sum(W_free_bin) / enbw
        metrics.update({
            "W_free_bin": W_free_bin, "SL_free_bin": to_SL(W_free_bin), "SL_free_density": to_SL(W_free_bin / (enbw * df)),
            "W_free_band": W_free_band, "SL_free_band": to_SL(W_free_band), "W_free_total": np.sum(W_free_bin) / enbw,
            "SL_free_total": to_SL(np.sum(W_free_bin) / enbw),
            "kr": kr, "R_sphere": r_s, "a_over_r": ratio, "geo_bias_dB": geo_bias_dB,
            "multipole_explained": p2_fit / p2_all, "W_mono_bin": W_mono_bin, "W_dip_bin": W_dip_bin,
            "multipole_explained_band": band_sum(p2_fit) / band_sum(p2_all),
            "mono_fraction_band": band_sum(W_mono_bin) / band_sum(W_free_bin),
            "dip_fraction_band" : band_sum(W_dip_bin) / band_sum(W_free_bin),
            "SL_multipole_band" : to_SL(band_sum(W_mono_bin + W_dip_bin) / enbw),
            "DI_band": DI_band, "dipole_elevation_band": axis_elev, "dipole_azimuth_band": axis_azim,
            "Excluded_sphere_fraction": 1.0 - np.count_nonzero(valid) / valid.size,
        })

    # --- 2. In situ: cylinder (waveguide) --- #
    if has_cyl:
        p_c   = case["P_cylinder"][mask]                                   # (Nf, Ntheta, Nz)
        r_c   = case["R_cylinder"]
        z_c   = case["Z_cylinder"]
        theta = np.deg2rad(case["Theta_deg_cylinder"])
        I_c   = np.abs(p_c)**2 / (2.0 * rho_wat * c_wat)
        I_c   = np.where(np.isfinite(I_c), I_c, 0.0)                       # Masked observers left out

        # Closed azimuthal integral (last point = first point), then trapezoid along z, for every bin
        theta_closed = np.append(theta, theta[0] + 2*np.pi)
        I_closed     = np.concatenate([I_c, I_c[:, :1, :]], axis=1)
        W_insitu_bin = np.abs(trapezoid(trapezoid(I_closed * r_c, theta_closed, axis=1), z_c, axis=1))

        W_insitu_band = band_sum(W_insitu_bin) / enbw
        metrics.update({
            "W_insitu_bin": W_insitu_bin, "SL_insitu_bin": to_SL(W_insitu_bin),
            "SL_insitu_density": to_SL(W_insitu_bin / (enbw * df)), "W_insitu_band": W_insitu_band,
            "SL_insitu_band": to_SL(W_insitu_band), "W_insitu_total": np.sum(W_insitu_bin) / enbw,
            "SL_insitu_total": to_SL(np.sum(W_insitu_bin) / enbw), "R_cylinder": r_c,
        })

    if has_sphere and has_cyl:
        with np.errstate(divide='ignore', invalid='ignore'):
            metrics["ratio_insitu_free"]      = 10.0 * np.log10(metrics["W_insitu_bin"] / metrics["W_free_bin"])
            metrics["ratio_insitu_free_band"] = 10.0 * np.log10(metrics["W_insitu_band"] / metrics["W_free_band"])

    # --- 3. Propagation loss at the spectrum observers: PL = SL - RL per band --- #
    if case.get("has_spectrums", False):
        SL_src  = metrics["SL_free_band"] if has_sphere else metrics["SL_insitu_band"]
        p_o     = case["P_spectrums"][mask]                                # (Nf, Nobs)
        RL_band = 10.0 * np.log10(band_sum(np.abs(p_o)**2 / 2.0) / enbw / p_ref**2)    # Band SPL [dB re 1 uPa]
        obs_o   = case["Obs_spectrums"]
        with np.errstate(invalid='ignore'):
            RL_mean = 10.0 * np.log10(np.nanmean(10.0**(RL_band / 10.0), axis=1))
        metrics.update({
            "PL_band": SL_src[:, None] - RL_band, "PL_mean_band": SL_src - RL_mean, "RL_band": RL_band,
            "PL_observers": obs_o, "PL_source": "free" if has_sphere else "insitu",
        })

    # --- 4. Print --- #
    if verbose:
        lo = f"{f.min():.1f}"; hi = f"{f.max():.1f}"
        print("\n" + "=" * 96)
        print(f"  EQUIVALENT POINT SOURCE [{tag}]  |  {lo} - {hi} Hz  |  df = {df:.3f} Hz, ENBW = {enbw:g} bins")
        print("=" * 96)
        label_w = 46
        if has_sphere:
            print(f"  {'SL free field, broadband':<{label_w}}: {metrics['SL_free_total']:.1f} dB re 1 µPa² m²  "
                  f"(W = {metrics['W_free_total']:.3e} W)")
            print(f"  {'Sphere r / farthest wet node (a/r)':<{label_w}}: {r_s:.1f} m / {a_max:.1f} m ({ratio:.2f})")
            if np.isfinite(geo_bias_dB) and geo_bias_dB > 0.5:
                print(f"  WARNING sphere hugs the structure: |p|² power ~+{geo_bias_dB:.1f} dB too high (geometric estimate) "
                      f"(r_sphere ≳ 3·a = {3*a_max:.0f} m keeps it below ~0.2 dB)")
            if np.any(kr < 3):
                print(f"  WARNING near field (kr < 3) below {f[kr < 3].max():.1f} Hz: |p|² power too high "
                      f"(dipole: up to +{10*np.log10(1 + 1/kr.min()**2):.1f} dB at {f.min():.1f} Hz)")
        if has_cyl:
            print(f"  {'SL in situ (waveguide), broadband':<{label_w}}: {metrics['SL_insitu_total']:.1f} dB re 1 µPa² m²  "
                  f"(W = {metrics['W_insitu_total']:.3e} W, r_cyl = {r_c:.0f} m)")
            print(f"  WARNING cylinder |p|² power ignores the oblique modal propagation: ~+0.5 to +4.7 dB too high "
                  f"(exact: 1/2 Re{{p v_r*}})")
        if has_sphere and has_cyl:
            print(f"  {'In situ / free power, broadband':<{label_w}}: "
                  f"{10*np.log10(metrics['W_insitu_total'] / metrics['W_free_total']):+.1f} dB")

        head = f"\n  {'fc [Hz]':>8} {'SL free':>8} {'SL situ':>8} {'situ-free':>9} {'DI':>5} {'mono%':>6} {'dip%':>6} " \
               f"{'expl%':>6} {'dip.elev':>8} {'PL mean':>8}"
        print(head + "\n  " + " " * 8 + f" {'[dB re 1 µPa² m²]':>17} {'[dB]':>9} {'[dB]':>5}" + " " * 21 + f"{'[deg]':>8} {'[dB]':>8}")
        for b, f_c in enumerate(fc):
            g = lambda key, fmt: format(metrics[key][b], fmt) if key in metrics and np.isfinite(metrics[key][b]) else "-"
            pct = lambda key: format(100*metrics[key][b], ".0f") if key in metrics and np.isfinite(metrics[key][b]) else "-"
            print(f"  {f_c:8.1f} {g('SL_free_band', '.1f'):>8} {g('SL_insitu_band', '.1f'):>8} "
                  f"{g('ratio_insitu_free_band', '+.1f'):>9} {g('DI_band', '.1f'):>5} {pct('mono_fraction_band'):>6} "
                  f"{pct('dip_fraction_band'):>6} {pct('multipole_explained_band'):>6} {g('dipole_elevation_band', '.0f'):>8} "
                  f"{g('PL_mean_band', '.1f'):>8}")
        print("  SL: power-equivalent free-field monopole, band levels (power in the band). mono%/dip%: share of the sphere "
              "power,\n  expl%: share of |p|² explained by monopole + dipole, PL: SL - band SPL at the spectrum observers "
              "(energetic mean).")

    # --- 5. Plots --- #
    if do_plot:
        color = pu.get_case_color(case)
        smooth = lambda y: au.smooth_fractional_octave(f, y, smooth_fraction, domain='db') if smooth_fraction else y
        note   = f" (1/{smooth_fraction} oct. smoothed)" if smooth_fraction else ""

        # 5.1 Source level spectra: narrowband densities and band-averaged densities (same units)
        fig, ax = plt.subplots(figsize=figsize)
        if has_sphere:
            ax.semilogx(f, smooth(metrics["SL_free_density"]), '-', color=color, label=f"Free field{note}")
            ax.semilogx(fc, metrics["SL_free_band"] - 10*np.log10(bw), 'o', color=color, label="Free field, 1/3 oct.")
            if np.any(kr < 3):
                ax.axvspan(f.min(), f[kr < 3].max(), color='0.85', zorder=0, label=r"Near field ($kr < 3$)")
        if has_cyl:
            ax.semilogx(f, smooth(metrics["SL_insitu_density"]), '--', color='k', label=f"In situ{note}")
            ax.semilogx(fc, metrics["SL_insitu_band"] - 10*np.log10(bw), 's', color='k', mfc='none', label="In situ, 1/3 oct.")
        totals = r" $|$ ".join(f"{name} {metrics[key]:.1f} dB" for name, key in (("free", "SL_free_total"), ("in situ", "SL_insitu_total"))
                               if key in metrics)
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel(r"SL [dB re 1 $\mu$Pa$^2$ m$^2$/Hz]")
        ax.set_title(f"{tag} - Equivalent point source\nBroadband SL: {totals} re 1 " + r"$\mu$Pa$^2$ m$^2$")
        ax.grid(True, which='both', ls='--', alpha=0.4)
        ax.legend(fontsize=11)
        fig.tight_layout()
        pu.set_figure_name(fig, "point_source", "SL")

        # 5.2 Radiation diagnostics per 1/3 octave band
        if has_sphere:
            n_rows = 3 if has_cyl else 2
            fig2, axs = plt.subplots(n_rows, 1, figsize=(figsize[0], 2.6 * n_rows), sharex=True)
            row = 0
            if has_cyl:
                axs[row].semilogx(fc, metrics["ratio_insitu_free_band"], '-o', color='k')
                axs[row].axhline(0.0, color='0.5', lw=0.8)
                axs[row].set_ylabel("In situ − free\n[dB]")
                row += 1
            axs[row].semilogx(fc, DI_band, '-o', color=color)
            axs[row].set_ylabel("DI [dB]")
            row += 1
            axs[row].semilogx(fc, 100*metrics["dip_fraction_band"], '-o', color=color, label="Dipole power")
            axs[row].semilogx(fc, 100*metrics["mono_fraction_band"], '--s', color=color, label="Monopole power")
            axs[row].semilogx(fc, 100*metrics["multipole_explained_band"], ':^', color='k', label=r"$|p|^2$ explained")
            axs[row].set_ylabel(r"Share [\%]" if plt.rcParams["text.usetex"] else "Share [%]")    # % starts a LaTeX comment
            axs[row].legend(fontsize=10)
            axs[-1].set_xlabel("1/3 octave band centre frequency [Hz]")
            for a in axs:
                a.grid(True, which='both', ls='--', alpha=0.4)
                if np.any(kr < 3):
                    a.axvspan(fc.min() / 2**(1/6), f[kr < 3].max(), color='0.85', zorder=0)
            axs[0].set_title(f"{tag} - Radiation diagnostics (sphere r = {r_s:.0f} m)")
            fig2.tight_layout()
            pu.set_figure_name(fig2, "point_source", "diagnostics")

        # 5.3 Propagation loss at the spectrum observers
        if "PL_band" in metrics:
            fig3, ax3 = plt.subplots(figsize=figsize)
            shades = pu.get_color_shades(color, metrics["PL_band"].shape[1])
            styles = pu.get_line_styles()
            for i, o in enumerate(metrics["PL_observers"]):
                ax3.semilogx(fc, metrics["PL_band"][:, i], styles[i % len(styles)], marker='o', color=shades[i],
                             label=f"({o[0]:.1f}, {o[1]:.1f}, {o[2]:.1f}) m")
            ax3.semilogx(fc, metrics["PL_mean_band"], '-', color='k', lw=2, label="Energetic mean")
            ax3.set_xlabel("1/3 octave band centre frequency [Hz]")
            ax3.set_ylabel(r"PL [dB re 1 m$^2$]")
            ax3.set_title(f"{tag} - Propagation loss, PL = SL$_{{{metrics['PL_source']}}}$ − SPL")
            ax3.grid(True, which='both', ls='--', alpha=0.4)
            ax3.legend(fontsize=11)
            fig3.tight_layout()
            pu.set_figure_name(fig3, "point_source", "PL")

    return metrics

