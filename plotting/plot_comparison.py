"""
Module: plot_comparison.py
Description: Side-by-side comparison of two or more single-turbine cases (one panel per case, drawn by the
             functions of plot_single_case). Cases whose file lacks a step are left out of that plot.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 09/2026
"""

import re
import numpy as np
import matplotlib.pyplot as plt
import plot_utils as pu
import plot_single_case as psc

MAX_COLS = 3    # [-] Panels per row: more cases wrap into new rows


# ---------- Helpers ---------- #
def select_cases(cases: list = None,     # [-] List of dictionaries containing simulation data
                 flag : str  = None,     # [-] has_<step> flag the plot needs
                 name : str  = None):    # [-] Name of the calling function, for the messages
    """
    Keeps the cases whose file contains the step (has_<step> flag True, or the given key present). The others
    (older versions of the Fortran code or interrupted runs) are left out with a message.

    Returns
    -------
    cases : list of dict
    """

    if not cases:
        print(f"{name}(): no cases given, skipping")
        return []

    def has(case):
        return bool(case.get(flag, False)) if flag.startswith("has_") else flag in case

    for case in cases:
        if not has(case):
            print(f"{name}(): {case['Name']} has no {flag.replace('has_', '')} data, left out of the comparison")

    return [case for case in cases if has(case)]

def make_panels(n      : int   = 1,         # [-] Number of cases
                figsize: tuple = (7, 5),    # [-] Size of a single panel
                stacked: bool  = False,     # [-] One column (True) or rows of up to MAX_COLS panels (False)
                **kwargs):                  # [-] Extra arguments of plt.subplots (sharey, subplot_kw, ...)
    """
    Creates one panel per case: a single row for up to MAX_COLS cases (more cases wrap into new rows), or a
    single column if stacked.

    Returns
    -------
    fig, axes, ncols : matplotlib.figure.Figure, np.ndarray of Axes (exactly n), int
    """

    if stacked:
        nrows, ncols = n, 1
    else:
        ncols = min(n, MAX_COLS)
        nrows = int(np.ceil(n / ncols))

    total_figsize = (figsize[0] * ncols, figsize[1] * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=total_figsize, squeeze=False, **kwargs)
    axes = axes.ravel()

    # Spare cells of the last row
    for ax in axes[n:]:
        fig.delaxes(ax)

    return fig, axes[:n], ncols

def remove_common_title_part(axes    : np.ndarray = None,     # [-] Axes of the comparison
                             pattern : str        = None):    # [-] Regex of the part to remove
    """
    Deletes from every title the part matching 'pattern' (e.g. ' - Depth = -15.0 [m]') if it is the same in
    all of them.
    """

    titles = [ax.get_title() for ax in axes]

    match = re.search(pattern, titles[0], re.IGNORECASE)

    if match:
        common_str = match.group(1)

        if all(common_str in t for t in titles):
            for ax in axes:
                new_title = ax.get_title().replace(common_str, '').strip()

                if new_title.endswith('-') or new_title.endswith(','):
                    new_title = new_title[:-1].strip()

                ax.set_title(new_title)

def wrap_case_titles(axes : np.ndarray = None,     # [-] Axes of the comparison
                     cases: list       = None):    # [-] Case plotted in each panel
    """
    Moves the case label of every panel title to its own line ('label - OASPL - ...' -> 'label\\nOASPL - ...'):
    comparison labels are long and the titles of neighbouring panels would overlap.
    """

    for ax, case in zip(axes, cases):
        title = ax.get_title()
        label = case["Case_type"]
        if title.startswith(label + " - "):
            ax.set_title(label + "\n" + title[len(label) + 3:])


# ---------- Structural Plots ---------- #
def plot_structure(cases  : list  = None,      # [-] List of dictionaries containing simulation data
                   mode   : str   = '3D',       # [-] Plot mode for the structure representation ('3D', 'xy', 'xz')
                   WindDir: bool  = True,       # [-] Whether to plot or not wind direction
                   figsize: tuple = (8, 6)):    # [-] Size of a single panel
    """
    Plots structure nodes of every case side by side (one panel per case).
    """

    cases = select_cases(cases, "Structure_nodes", "plot_structure")
    if not cases:
        return None, None

    subplot_kw = {'projection': '3d'} if mode == '3D' else None
    sharey = False if mode == '3D' else True

    fig, axes, ncols = make_panels(len(cases), figsize, subplot_kw=subplot_kw, sharey=sharey)

    for i, case in enumerate(cases):
        psc.plot_structure(case=case,
                           mode=mode,
                           WindDir=WindDir,
                           figsize=figsize,
                           ax_in=axes[i])
        if i % ncols > 0 and mode != '3D':
            axes[i].set_ylabel("")

    wrap_case_titles(axes, cases)

    fig.tight_layout()
    pu.set_figure_name(fig, "structure", mode)
    return fig, axes


# ---------- Spectrum Plots ---------- #
def plot_spectrum(cases        : list  = None,      # [-] List of dictionaries containing simulation data
                  mode         : str   = 'SPL',     # [-] Selects Sound Pressure/Power Level. 'SPL' or 'SWL'
                  absorption   : bool  = True,      # [-] Whether to apply approximate absorption attenuation
                  octave       : bool  = True,      # [-] 1/3 octave or fine resolution
                  filter_under : float = None,      # [Hz] Minimum frequency threshold
                  filter_over  : float = None,      # [Hz] Maximum frequency threshold
                  do_thresholds: bool  = False,     # [-] Whether to plot regulatory/reference thresholds
                  figsize      : tuple = (7, 5)):  # [-] Size of a single panel
    """
    Plots the acoustic spectrum of every case side by side (one panel per case).
    """

    cases = select_cases(cases, "has_spectrums", "plot_spectrum")
    if not cases:
        return None, None

    fig, axes, ncols = make_panels(len(cases), figsize, sharey=True)

    for i, case in enumerate(cases):
        psc.plot_spectrum(case=case,
                          mode=mode,
                          absorption=absorption,
                          octave=octave,
                          filter_under=filter_under,
                          filter_over=filter_over,
                          do_thresholds=do_thresholds,
                          figsize=figsize,
                          ax_in=axes[i])
        if i % ncols > 0:
            axes[i].set_ylabel("")

    fig.tight_layout()
    pu.set_figure_name(fig, "spectrum", mode, "octave" if octave else "narrowband")
    return fig, axes


# ---------- Polar Plots ---------- #
def plot_polar(cases           : list  = None,         # [-] List of dictionaries containing simulation data
               mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', or 'ABS'
               target_frequency: float = 0.88,         # [Hz] Target frequency for 'SPL' or 'ABS' modes
               absorption      : bool  = True,         # [-] Whether to apply absorption attenuation
               bands           : list  = [10., 100.],  # [Hz] Frequency band limits for OASPL
               figsize         : tuple = (7, 7)):      # [-] Size of a single panel
    """
    Plots polar directivity diagram of every case side by side (one panel per case).
    """

    cases = select_cases(cases, "has_polar", "plot_polar")
    if not cases:
        return None, None

    fig, axes, ncols = make_panels(len(cases), figsize, subplot_kw={"projection": "polar"})

    for i, case in enumerate(cases):
        psc.plot_polar(case=case,
                       mode=mode,
                       target_frequency=target_frequency,
                       absorption=absorption,
                       bands=bands,
                       figsize=figsize,
                       ax_in=axes[i])

    # Legend handling
    handles, labels = axes[0].get_legend_handles_labels()

    for ax in axes:
        if ax.get_legend() is not None:
            ax.get_legend().remove()

    fig.legend(handles, labels, loc='lower center', bbox_to_anchor = (0.5, 0.02), ncol=len(handles), frameon=False)

    fig.tight_layout(rect=[0, 0.1, 1, 0.95])

    pu.set_figure_name(fig, "polar", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "")
    return fig, axes


# ---------- Cylindrical Plots ---------- #
def plot_cylinder(cases           : list  = None,      # [-] List of dictionaries containing simulation data
                  mode            : str   = 'OASPL',   # [-] Magnitude to plot: 'OASPL', 'SPL', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                  target_frequency: float = 0.88,      # [Hz] Target frequency for single-frequency modes
                  absorption      : bool  = True,      # [-] Whether to apply absorption attenuation
                  filter_under    : float = None,      # [Hz] Lower frequency cutoff (f >= filter_under)
                  filter_over     : float = None,      # [Hz] Upper frequency cutoff (f <= filter_over)
                  figsize         : tuple = (10, 5),   # [-] Size of a single panel
                  cmap            : str   = 'inferno'):# [-] Contour colormap
    """
    Plots 2D unwrapped cylindrical surface maps of every case, stacked in a single column.
    """

    cases = select_cases(cases, "has_cylinder", "plot_cylinder")
    if not cases:
        return None, None

    fig, axes, ncols = make_panels(len(cases), figsize, stacked=True, sharex=True)

    case_names = []
    mode_suffix = ""

    for i, case in enumerate(cases):
        psc.plot_cylinder(case=case,
                          mode=mode,
                          target_frequency=target_frequency,
                          absorption=absorption,
                          filter_under=filter_under,
                          filter_over=filter_over,
                          figsize=figsize,
                          cmap=cmap,
                          ax_in=axes[i])

        # Hide x lables
        if i < len(cases) - 1:
            axes[i].set_xlabel("")

        # suptitle handling
        title_text = axes[i].get_title()
        axes[i].set_title("")

        # Separate name from magnitude label
        if " - " in title_text:
            name, suffix = title_text.rsplit(" - ", 1)
            case_names.append(name)
            # Save first suffix
            if i == 0:
                mode_suffix = suffix
        else:
            case_names.append(title_text)

        # Every panel keeps its case name
        axes[i].set_title(case_names[-1], fontsize=14)

    # Reconstruct gloabl title
    if mode_suffix:
        global_title = mode_suffix
    else:
        global_title = ", ".join(case_names)

    fig.suptitle(global_title, fontsize=16, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    pu.set_figure_name(fig, "cylinder", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "")
    return fig, axes


# ---------- Linear Plots ---------- #
def plot_distance_decay(cases           : list  = None,      # [-] List of dictionaries containing simulation data
                        mode            : str   = 'OASPL',   # [-] Magnitude to plot: 'OASPL', 'SPL', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                        target_frequency: float = 0.88,      # [Hz] Target frequency for single-frequency modes
                        absorption      : bool  = True,      # [-] Whether to apply absorption attenuation
                        filter_under    : float = None,      # [Hz] Lower frequency cutoff (f >= filter_under)
                        filter_over     : float = None,      # [Hz] Upper frequency cutoff (f <= filter_over)
                        turbine_data    : bool  = False,     # [-] Whether to plot additional turbine data from GE 2025
                        figsize         : tuple = (7, 5)):  # [-] Size of a single panel
    """
    Plots acoustic magnitude decay over distance of every case side by side (one panel per case).
    """

    cases = select_cases(cases, "has_decay", "plot_distance_decay")
    if not cases:
        return None, None

    fig, axes, ncols = make_panels(len(cases), figsize, sharey=True)

    for i, case in enumerate(cases):
        psc.plot_distance_decay(case=case,
                                mode=mode,
                                target_frequency=target_frequency,
                                absorption=absorption,
                                filter_under=filter_under,
                                filter_over=filter_over,
                                turbine_data=turbine_data,
                                figsize=figsize,
                                ax_in=axes[i])
        if i % ncols > 0:
            axes[i].set_ylabel("")

    wrap_case_titles(axes, cases)

    handles, labels = axes[0].get_legend_handles_labels()

    for ax in axes:
        if ax.get_legend() is not None:
            ax.get_legend().remove()

    if handles:
        fig.legend(handles, labels,
                   loc='lower center',
                   bbox_to_anchor=(0.5, 0.02),
                   ncol=len(handles),
                   fontsize=12,
                   frameon=False)

        fig.tight_layout(rect=[0, 0.12, 1, 1])
    else:
        fig.tight_layout()

    pu.set_figure_name(fig, "decay", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "")
    return fig, axes

def plot_line(cases           : list  = None,      # [-] List of dictionaries containing simulation data
              mode            : str   = 'OASPL',   # [-] Magnitude to plot: 'OASPL', 'SPL', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
              target_frequency: float = 0.88,      # [Hz] Target frequency for single-frequency modes
              absorption      : bool  = True,      # [-] Whether to apply absorption attenuation
              filter_under    : float = None,      # [Hz] Lower frequency cutoff (f >= filter_under)
              filter_over     : float = None,      # [Hz] Upper frequency cutoff (f <= filter_over)
              figsize         : tuple = (7, 5)):  # [-] Size of a single panel
    """
    Plots acoustic magnitude along a straight line of every case side by side (one panel per case).
    """

    cases = select_cases(cases, "has_line", "plot_line")
    if not cases:
        return None, None

    fig, axes, ncols = make_panels(len(cases), figsize, sharey=True)

    ymin, ymax = float('inf'), float('-inf')
    for i, case in enumerate(cases):
        psc.plot_line(case=case,
                      mode=mode,
                      target_frequency=target_frequency,
                      absorption=absorption,
                      filter_under=filter_under,
                      filter_over=filter_over,
                      figsize=figsize,
                      ax_in=axes[i])
        if i % ncols > 0:
            axes[i].set_ylabel("")

        ymin_local, ymax_local = axes[i].get_ylim()

        ymin = min(ymin, ymin_local)
        ymax = max(ymax, ymax_local)

        axes[i].set_ylim(ymin, ymax)

    wrap_case_titles(axes, cases)

    pu.set_figure_name(fig, "line", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "")
    return fig, axes


# ---------- Slice Plots ---------- #
def plot_sliceXY(cases           : list  = None,      # [-] List of dictionaries containing simulation data
                 mode            : str   = 'OASPL',   # [-] Magnitude to plot: 'OASPL', 'SPL', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                 target_frequency: float = 0.88,      # [Hz] Target frequency for single-frequency modes
                 absorption      : bool  = True,      # [-] Whether to apply absorption attenuation
                 filter_under    : float = None,      # [Hz] Lower frequency cutoff (f >= filter_under)
                 filter_over     : float = None,      # [Hz] Upper frequency cutoff (f <= filter_over)
                 structure       : bool  = True,      # [-] Whether to plot structure on slice
                 figsize         : tuple = (7, 7),    # [-] Size of a single panel
                 cmap            : str   = 'inferno'):# [-] Contour colormap
    """
    Plots a 2D spatial XY slice of every case side by side (one panel per case).
    """

    cases = select_cases(cases, "has_slicexy", "plot_sliceXY")
    if not cases:
        return None, None

    fig, axes, ncols = make_panels(len(cases), figsize)

    for i, case in enumerate(cases):
        psc.plot_sliceXY(case=case,
                         mode=mode,
                         target_frequency=target_frequency,
                         absorption=absorption,
                         filter_under=filter_under,
                         filter_over=filter_over,
                         structure=structure,
                         figsize=figsize,
                         cmap=cmap,
                         ax_in=axes[i])

        if i % ncols > 0:
            axes[i].set_ylabel("")

    # Delete Depth label if it is the same for all plots
    remove_common_title_part(axes, r'([-\(,\s]*(?:Depth|Z)\s*[:=]\s*[-0-9\.]+\s*(?:\[m\]|m)?[\)]*)')
    wrap_case_titles(axes, cases)

    fig.tight_layout()
    pu.set_figure_name(fig, "sliceXY", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "")
    return fig, axes

def plot_sliceXZ(cases           : list  = None,      # [-] List of dictionaries containing simulation data
                 mode            : str   = 'OASPL',   # [-] Magnitude to plot: 'OASPL', 'SPL', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                 target_frequency: float = 0.88,      # [Hz] Target frequency for single-frequency modes
                 absorption      : bool  = True,      # [-] Whether to apply absorption attenuation
                 filter_under    : float = None,      # [Hz] Lower frequency cutoff (f >= filter_under)
                 filter_over     : float = None,      # [Hz] Upper frequency cutoff (f <= filter_over)
                 structure       : bool  = True,      # [-] Whether to plot structure on slice
                 figsize         : tuple = (7, 7),    # [-] Size of a single panel
                 cmap            : str   = 'inferno'):# [-] Contour colormap
    """
    Plots a 2D spatial XZ slice of every case side by side (one panel per case).
    """

    cases = select_cases(cases, "has_slicexz", "plot_sliceXZ")
    if not cases:
        return None, None

    fig, axes, ncols = make_panels(len(cases), figsize)

    for i, case in enumerate(cases):
        psc.plot_sliceXZ(case=case,
                         mode=mode,
                         target_frequency=target_frequency,
                         absorption=absorption,
                         filter_under=filter_under,
                         filter_over=filter_over,
                         structure=structure,
                         figsize=figsize,
                         cmap=cmap,
                         ax_in=axes[i])

        if i % ncols > 0: axes[i].set_ylabel("")

    # Delete the plane position (Y = ...) if it is the same for all plots
    remove_common_title_part(axes, r'([-\(,\s]*(?:Angle|Azimuth|Y)\s*[:=]\s*[-0-9\.]+\s*(?:\[.*?\]|deg|°|m)?[\)]*)')
    wrap_case_titles(axes, cases)

    fig.tight_layout()
    pu.set_figure_name(fig, "sliceXZ", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "")
    return fig, axes

def plot_sliceVertical(cases           : list  = None,      # [-] List of dictionaries containing simulation data
                       mode            : str   = 'OASPL',   # [-] Magnitude to plot: 'OASPL', 'SPL', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                       target_frequency: float = 0.88,      # [Hz] Target frequency for single-frequency modes
                       absorption      : bool  = True,      # [-] Whether to apply absorption attenuation
                       filter_under    : float = None,      # [Hz] Lower frequency cutoff (f >= filter_under)
                       filter_over     : float = None,      # [Hz] Upper frequency cutoff (f <= filter_over)
                       structure       : bool  = True,      # [-] Whether to plot structure on slice
                       figsize         : tuple = (7, 7),    # [-] Size of a single panel
                       cmap            : str   = 'inferno'):# [-] Contour colormap
    """
    Plots a 2D spatial vertical slice (any azimuth) of every case side by side (one panel per case).
    """

    cases = select_cases(cases, "has_slicevertical", "plot_sliceVertical")
    if not cases:
        return None, None

    fig, axes, ncols = make_panels(len(cases), figsize)

    for i, case in enumerate(cases):
        psc.plot_sliceVertical(case=case,
                               mode=mode,
                               target_frequency=target_frequency,
                               absorption=absorption,
                               filter_under=filter_under,
                               filter_over=filter_over,
                               structure=structure,
                               figsize=figsize,
                               cmap=cmap,
                               ax_in=axes[i])

        if i % ncols > 0: axes[i].set_ylabel("")

    # Delete the plane azimuth if it is the same for all plots
    remove_common_title_part(axes, r'([-\(,\s]*(?:Angle|Azimuth)\s*[:=]\s*[-0-9\.]+\s*(?:\[.*?\]|deg|°)?[\)]*)')
    wrap_case_titles(axes, cases)

    fig.tight_layout()
    pu.set_figure_name(fig, "sliceVertical", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "")
    return fig, axes


# ---------- Compute metrics ---------- #
def compute_sphere_metrics(cases       : list  = None,     # [-] List of dictionaries containing simulation data
                           **kwargs):                      # [-] Arguments of plot_single_case.compute_sphere_metrics
    """
    Computes and prints the sphere radiation metrics of every case.

    Returns
    -------
    metrics : list of dict (None for the cases without sphere data)
    """

    return [psc.compute_sphere_metrics(case, **kwargs) for case in cases]

def compute_cylinder_metrics(cases       : list  = None,     # [-] List of dictionaries containing simulation data
                             **kwargs):                      # [-] Arguments of plot_single_case.compute_cylinder_metrics
    """
    Computes and prints the cylinder propagation metrics of every case.

    Returns
    -------
    metrics : list of dict (None for the cases without cylinder data)
    """

    return [psc.compute_cylinder_metrics(case, **kwargs) for case in cases]

