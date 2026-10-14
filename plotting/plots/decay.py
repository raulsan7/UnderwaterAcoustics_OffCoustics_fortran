"""
Module: plots/decay.py
Description: Curves along the distances of the 'decay' and 'line' steps (get_curve, shared with plots/line.py)
             and the decay with distance of the 'decay' step, for any emitter: one curve per member of a family
             (--ws), 1/r and 1/r^2 references and the external reference data of the case (measured turbine
             noise, Ge 2025).

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import matplotlib.pyplot as plt
import utils.AcousticUtils as au
from wrappers.case import emitter_label, case_color, get_family, family_tag
from plots.common import no_absorption, psd_not_applicable
from plot_utils.conversions import convert_to, is_valid
from plot_utils.drawing import add_family_legends
from plot_utils.saving import set_figure_name


# ---------- Helpers ---------- #
def get_curve(case            : dict  = None,     # [-] Case
              step            : str   = None,     # [-] Step whose observers form the curve: 'decay' or 'line'
              mode            : str   = 'OASPL',  # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
              target_frequency: float = 0.88,     # [Hz] Target frequency for single-frequency modes
              absorption      : bool  = True,     # [-] Whether to apply absorption attenuation
              filter_under    : float = None,     # [Hz] Lower frequency cutoff (f >= filter_under)
              filter_over     : float = None,     # [Hz] Upper frequency cutoff (f <= filter_over)
              name            : str   = None):    # [-] Name of the calling plot, for the messages
    """
    Magnitude of the pressure of one case along the distances of its 'decay' or 'line' step: the curve plot_distance_decay
    and plot_line draw (one per member of a family). Distances without data (NaN) stay NaN.

    Returns
    -------
    r, vals, unit_label, title_extra : distances, magnitude, label of the y axis and note of the title
        ('at 0.88 Hz' for the single-frequency modes, empty for OASPL). None if the case cannot be plotted.
    """

    p     = case[f"P_{step}"]
    Freqs = case["Freqs"]
    df    = case["bin_width"]               # Frequency bin width [Hz] for 'PSD'
    r     = case[f"Distances_{step}"]
    p_ref = case["p_ref"]

    if not is_valid(p, name):
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
        vals = au.pressure_to_OASPL(p, p_ref, absorption, Freqs, r, enbw=case["ENBW"], water=case["water"])
        unit_label = r"OASPL [dB re $\mu$ Pa]"
    else:
        idx = np.argmin(np.abs(Freqs - target_frequency))
        f   = Freqs[idx]
        title_extra = f"at {f:.2f} Hz"

        vals, unit_label = convert_to(p[idx, :], mode, df, enbw=case["ENBW"])
        if absorption and mode in ('SPL', 'PSD'): vals = au.add_absorption(f, vals, r, water=case["water"])

    return r, vals, unit_label, title_extra

def get_family_curves(case            : dict  = None,     # [-] Case
                      step            : str   = None,     # [-] 'decay' or 'line'
                      name            : str   = None,     # [-] Name of the calling plot, for the messages
                      **kwargs):                          # [-] Arguments of get_curve
    """
    Curves of a plot of a 'decay' or 'line' step: the case alone or, with a family (--ws), one per member.
    Members whose file lacks the step, or that cannot be plotted, are left out with a message.

    Returns
    -------
    curves, family_mode : list of (dict, tuple), bool
        Case and its get_curve() result, in the order of the family; whether the curves are a family.
    """

    cases, ws_mode = get_family(case)

    curves = []
    for c in cases:
        if not c[f"has_{step}"]:
            print(f"{name}(): {c['Name']} has no {step} data, left out")
            continue

        curve = get_curve(c, step, name=name, **kwargs)
        if curve is not None:
            curves.append((c, curve))

    return curves, ws_mode


# ---------- Decay ---------- #
def plot_distance_decay(case            : dict  = None,         # [-] Case
                        mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                        target_frequency: float = 0.88,         # [Hz] Target frequency for single-frequency modes
                        absorption      : bool  = True,         # [-] Whether to apply absorption attenuation
                        filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                        filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                        turbine_data    : bool  = False,        # [-] Wheter to plot the reference data of the case (turbines: Ge 2025)
                        figsize         : tuple = (10,5),       # [-] Figure size
                        ax_in           : plt.Axes = None):     # [-] Optional external axes for subplots
    """
    Plots acoustic magnitude decay over distance (1D line plot). Distances without data (NaN, inside the
    structure) are gaps in the curve; the 1/r and 1/r^2 references start at the first distance with data.
    With a family (--ws) one curve per member (color = member, legend below the axes); the references (black)
    start at the curve of the case itself and the curves are not filled. turbine_data draws the external
    reference data of the case (case["references"]: measured turbine noise of Ge 2025 and its regression).

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    if not case["has_decay"]:
        print("plot_distance_decay(): no decay data, skipping")
        return None, None

    absorption = no_absorption(case, absorption, "plot_distance_decay")
    if psd_not_applicable(get_family(case)[0], mode, "plot_distance_decay"):
        return None, None

    tag    = emitter_label(case)
    do_log = case["Logspace_decay"]

    curves, ws_mode = get_family_curves(case, "decay", "plot_distance_decay", mode=mode, target_frequency=target_frequency,
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

    # References from the first distance with data of the case itself (the first curve if it was left out)
    c, (r, vals, unit_label, title_extra) = next((curve for curve in curves if curve[0] is case), curves[0])

    i0 = np.flatnonzero(np.isfinite(vals))[0]      # First distance with data: origin of the 1/r references
    ref_color  = 'k' if ws_mode else color
    ref_vals_1 = vals[i0] - 10 * np.log10(np.maximum(r, np.finfo(float).eps) / r[i0])
    ref_vals_2 = vals[i0] - 20 * np.log10(np.maximum(r, np.finfo(float).eps) / r[i0])
    draw(r, ref_vals_1, '--', color=ref_color, alpha=0.4, lw=1.0, label=r"$1/r$")
    draw(r, ref_vals_2, ':', color=ref_color, alpha=0.4, lw=1.0, label=r"$1/r^2$")

    min_val = min(ref_vals_2.min(), *(np.nanmin(vals_c) for _, (_, vals_c, _, _) in curves))

    # External reference data of the case (turbines: Ge 2025). The curve below is filled with the color of
    # the last data set, as so far
    references = case["references"] if turbine_data else None
    if turbine_data and references is None:
        print(f"plot_distance_decay(): {case['Name']} has no reference data, turbine_data ignored")
    if references:
        all_turbine_x, all_turbine_y = [], []
        for (data, label, marker, color) in references:
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
        add_family_legends(ax, handles, labels, family_title=case["family_title"])

    if ax_in is None:
        fig.tight_layout()
        set_figure_name(fig, "decay", mode.upper(), title_extra, family_tag([case]) if ws_mode else "")

    return fig, ax
