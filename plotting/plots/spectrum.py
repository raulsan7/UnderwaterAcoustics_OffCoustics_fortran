"""
Module: plots/spectrum.py
Description: Spectrum of the observers of the 'spectrums' step, for any emitter: SPL or PSD, 1/3 octave or
             narrowband (smoothed), the means of observers at the same distance and, with a family (--ws), one
             set of curves per member.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
import utils.AcousticUtils as au
import utils.EnvironmentalData as ed
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from wrappers.case import emitter_label, case_color, get_family, family_tag
from plots.common import no_absorption, psd_not_applicable, octave_not_applicable
from plot_utils.colors import get_line_styles, get_color_shades
from plot_utils.conversions import is_valid
from plot_utils.drawing import add_family_legends
from plot_utils.saving import set_figure_name


def get_spectrum(case         : dict  = None,      # [-] Case
                 mode         : str   = 'SPL',     # [-] 'SPL' or 'PSD'
                 absorption   : bool  = True,      # [-] Wheter to apply approximate absorption attenuation
                 octave       : bool  = True,      # [-] 1/3 octave or fine resolution
                 filter_under : float = None,      # [Hz] Minimum frequency threshold
                 filter_over  : float = None):     # [Hz] Maximum frequency threshold
    """
    Spectrum of every observer of one case, as plot_spectrum draws it (once per member of a family): SPL or
    PSD in dB, narrowband or 1/3 octave, within the frequency filter, and the OASPL of every observer.

    Returns
    -------
    spectrum : dict
        Freqs, spl (Nf, Nobs), oaspl (Nobs,), distances (Nobs,) to the origin (or the centre) of the case and
        ylabel. None if the case cannot be plotted.
    """

    Freqs     = case["Freqs"]
    Observers = case["Obs_spectrums"]
    p         = case["P_spectrums"]
    p_ref     = case["p_ref"]

    if not is_valid(p, "plot_spectrum"):
        return None

    # Distance to the reference of the case (turbine axis; centre of a farm): absorption and equal distances
    ref = case["origin"] if case["origin"] is not None else case["centre"]
    if ref is not None:
        distances = np.sqrt((Observers[:,0] - ref[0])**2 + (Observers[:,1] - ref[1])**2)
    else:
        distances = np.full(len(Observers), np.nan)

    mode = mode.upper()

    if mode == "SPL":
        spl = au.pressure_to_SPL(p, p_ref, absorption, Freqs, distances, water=case["water"])
        ylabel = r"SPL [dB re 1 $\mu$ Pa]"
    elif mode == "PSD":
        df = case["bin_width"]
        if df is None: raise ValueError("plot_spectrum(): mode 'PSD' needs at least two frequencies to get the bin width")
        spl = au.pressure_to_PSD(p, df, p_ref, absorption, Freqs, distances, enbw=case["ENBW"], water=case["water"])
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
    oaspl = au.pressure_to_OASPL(p[mask], p_ref, absorption, Freqs, distances, enbw=case["ENBW"], water=case["water"])

    if octave:
        Freqs, spl = au.to_third_octave(Freqs, spl, fmin=Freqs.min(), fmax=Freqs.max())

        # PSD is a density: the energetic sum of its bins is band power / df, so rescale
        # it to band power / band width (spectrum level), independent of the resolution
        if mode == "PSD":
            bw    = Freqs * (2**(1/6) - 2**(-1/6))      # [Hz] 1/3 octave band widths
            shape = [-1] + [1] * (spl.ndim - 1)
            spl   = spl - 10.0 * np.log10(bw / df).reshape(shape)

    return dict(Freqs=Freqs, spl=spl, oaspl=oaspl, distances=distances, ylabel=ylabel)

def plot_spectrum(case         : dict  = None,      # [-] Case
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
    message. If every observer with data lies at the same distance to the reference of the case (turbine
    axis, centre of a farm), their energetic mean (solid black) and arithmetic dB mean (dashed black) are
    added. In narrowband the means are always smoothed over 1/smooth_fraction octave, so they read as a trend
    instead of a saw, and with smooth=True so is every observer curve. The means are computed from the raw
    curves either way, and OASPL never sees the smoothing.

    With a family (--ws) the spectra of every member are drawn: color = member and line style = observer,
    with a legend of each below the axes. Observers at the same distance are not drawn one by one but as their
    energetic mean and the range they span (shaded), so the plot stays readable.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The generated figure and axes objects for further editing.
    """

    if not case["has_spectrums"]:
        print("plot_spectrum(): no spectrum data, skipping")
        return None, None

    absorption = no_absorption(case, absorption, "plot_spectrum")

    tag   = emitter_label(case)
    color = case_color(case)

    cases, ws_mode = get_family(case)
    if psd_not_applicable(cases, mode, "plot_spectrum") or octave_not_applicable(cases, octave, "plot_spectrum"):
        return None, None

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

    linestyles = get_line_styles()

    # Narrowband smoothing: of every curve (smooth=True) or only of the means; the raw spl feeds the means
    smooth_means  = (not octave) and bool(smooth_fraction)
    smooth_curves = smooth_means and smooth
    smooth_note   = f"1/{smooth_fraction} oct. smoothed" if smooth_means else ""

    ws_handles, ws_labels = [], []                        # Family: one entry per member (color)
    style_handles, style_labels = [], []                  # Family: one entry per observer (line style) or mean

    for c, spectrum in spectra:
        Freqs, spl, oaspl, distances = spectrum["Freqs"], spectrum["spl"], spectrum["oaspl"], spectrum["distances"]
        Observers = c["Obs_spectrums"]

        num_observers = spl.shape[1] if spl.ndim> 1 else 1
        if num_observers > len(linestyles): print("plot_spectrum(): too many observers, linestyles will be repeated")

        spl_plot = au.smooth_fractional_octave(Freqs, spl, smooth_fraction, domain='db') if smooth_curves else spl

        # Observers with data: those inside the structure carry NaN at every frequency
        valid = [i for i in range(num_observers) if not np.all(np.isnan(spl_plot[:,i] if spl_plot.ndim > 1 else spl_plot))]

        # Observers all at the same distance to the reference (e.g. around the platform): energetic and dB means
        same_distance = len(valid) > 1 and np.allclose(distances[valid], distances[valid[0]], rtol=0.0, atol=0.1)

        if ws_mode:
            color  = c["family_color"]
            kwargs = dict(lw=c["family_lw"], zorder=c["family_zorder"])
            ws_handles.append(Line2D([], [], color=color, lw=c["family_lw"]))
            ws_labels.append(c["family_label"])
        else:
            kwargs = {}

        shades = [color] * num_observers if ws_mode else get_color_shades(color, num_observers)    # One shade per observer

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

            # The observers are the same for every member: the line styles are explained once
            if ws_mode and c is spectra[0][0] and not same_distance:
                style_handles.append(ax.plot([], [], ls, color='k')[0])
                style_labels.append(coords)

            print(f"\nObserver: {coords}" + (f" | r = {distances[i]:.1f} m" if np.isfinite(distances[i]) else ""))
            print(f"{tag}{', ' + c['family_label'] if ws_mode else ''}: OASPL = {oaspl[i]:.1f}")

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

        add_family_legends(ax, ws_handles, ws_labels, style_handles, style_labels,
                           style_title=None if do_thresholds else ("Observers" if style_handles else None),
                           ncol_style=2 if do_thresholds else 1, family_title=case["family_title"])
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
        set_figure_name(fig, "spectrum", mode.upper(), "octave" if octave else "narrowband", family_tag([case]) if ws_mode else "")

    return fig, ax
