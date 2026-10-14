"""
Module: plots/polar.py
Description: Polar directivity of the 'polar' step, for any emitter: OASPL per frequency band or a single
             frequency and, with a family (--ws), one curve per member.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import matplotlib.pyplot as plt
import utils.AcousticUtils as au
from matplotlib.lines import Line2D
from wrappers.case import emitter_label, get_family, family_tag
from plots.common import no_absorption, psd_not_applicable
from plot_utils.colors import get_band_colors, get_line_styles
from plot_utils.conversions import convert_to, is_valid
from plot_utils.drawing import add_family_legends
from plot_utils.saving import set_figure_name


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

def plot_polar(case             : dict  = None,        # [-] Case
               mode             : str   = 'OASPL',     # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', or 'ABS'
               target_frequency : float = 0.88,        # [Hz] Target frequency for 'SPL', 'PSD' or 'ABS' modes
               absorption       : bool  = True,        # [-] Whether to apply absorption attenuation
               bands            : list  = [10., 100.], # [Hz] Frequency band limits for OASPL
               figsize          : tuple = (7,7),       # [-] Figure size
               ax_in            : plt.Axes = None):    # [-] Optional external axes for subplots
    """
    Plots polar directivity diagram. Angles without data (NaN, inside the structure) are gaps in the curves
    and holes down to the centre in their fill.

    With a family (--ws) every member is drawn: color = member and, in OASPL, line style = frequency band,
    with a legend of each below the axes. The curves are not filled.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The generated figure and axes objects for further editing.
    """

    if not case["has_polar"]:
        print("plot_polar(): no polar data, skipping")
        return None, None

    absorption = no_absorption(case, absorption, "plot_polar")

    tag    = emitter_label(case)
    colors = get_band_colors()

    cases, ws_mode = get_family(case)
    if psd_not_applicable(cases, mode, "plot_polar"):
        return None, None

    polars = []
    for c in cases:
        if not c["has_polar"]:
            print(f"plot_polar(): {c['Name']} has no polar data, left out")
        elif is_valid(c["P_polar"], "plot_polar"):
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
        df         = c["bin_width"]         # Frequency bin width [Hz] for 'PSD'
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

                oaspl = au.OASPL_band(p, mask, p_ref, absorption, Freqs, distances, enbw=c["ENBW"], water=c["water"])
                series.append((c, theta_plot, np.append(oaspl, oaspl[0]), label, b))

        else:
            # Find closest frequency index for single-frequency modes
            idx = np.argmin(np.abs(Freqs - target_frequency))
            f   = Freqs[idx]
            name_extra = f"{f:.2f}Hz"

            vals, unit_label = convert_to(p[idx, :], mode, df, enbw=c["ENBW"])

            if absorption and mode in ("SPL", "PSD"): vals = au.add_absorption(f, vals, distances, water=c["water"])

            series.append((c, theta_plot, np.append(vals, vals[0]), f"{unit_label} at {f:.2f} [Hz]", 0))

    if series:
        vals = np.concatenate([vals_plot for _, _, vals_plot, _, _ in series])
        rmin, rmax = np.nanmin(vals), np.nanmax(vals)
        if mode == "OASPL" and rmax > 0: rmin = max(0., rmin)
        margin = abs(rmax-rmin) * 0.1 if rmax != rmin else abs(rmin) * 0.1
        if margin == 0: margin = 1.0 # Fallback

        # Plot (angles without data: gap in the curve, fill collapsed to the centre)
        linestyles = get_line_styles()
        for c, theta_plot, vals_plot, label, b in series:
            if ws_mode:
                ax.plot(theta_plot, vals_plot, linestyles[b % len(linestyles)] if mode == "OASPL" else '-',
                        lw=c["family_lw"], color=c["family_color"], zorder=c["family_zorder"])
            else:
                ax.plot(theta_plot, vals_plot, lw=2, color=colors[b], label=label)
                ax.fill(theta_plot, np.where(np.isfinite(vals_plot), vals_plot, rmin - margin), alpha=0.08, color=colors[b])

        ax.set_rmin(rmin - margin)

    ax.set_title(f"{tag}")

    if ws_mode:
        ws_handles = [Line2D([], [], color=c["family_color"], lw=c["family_lw"]) for c in polars]
        ws_labels  = [c["family_label"] for c in polars]

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

        add_family_legends(ax, ws_handles, ws_labels, band_handles, band_labels, style_title="Frequency band", offset=0.5,
                           family_title=case["family_title"])
    else:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)

    if ax_in is None:
        if ws_mode:
            fig.tight_layout()
        else:
            fig.tight_layout(rect=[0, 0.06, 1, 1])
        set_figure_name(fig, "polar", mode, name_extra, family_tag([case]) if ws_mode else "")

    return fig, ax
