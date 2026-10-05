"""
Module: plot_comparison.py
Description: Side-by-side comparison of two or more single-turbine cases (one panel per case, drawn by the
             functions of plot_single_case). Cases whose file lacks a step are left out of that plot. The
             equivalent point source overlays the cases instead (same units, the point is their levels).

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 09/2026
"""

import re
import warnings
import numpy as np
import matplotlib.pyplot as plt
import plot_utils as pu
import plot_single_case as psc
import utils.AcousticUtils as au
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

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

def ws_suffix(cases: list = None):     # [-] List of dictionaries containing simulation data
    """
    Word that tells the figures drawn with every wind speed (--ws) from the plain ones in the name they are
    saved with.
    """

    return "allws" if any("WS_cases" in case for case in cases) else ""

def add_case_legend(fig,                    # [-] Figure of the comparison
                    handles: list = None,   # [-] One handle per case
                    labels : list = None,   # [-] One label per case
                    title  : str  = None):  # [-] Legend title
    """
    Adds one legend entry per case below the axes, one per row, for figures that overlay every case (their
    labels are long).

    Returns
    -------
    bottom : float
        Fraction of the figure height the legend takes, for fig.tight_layout(rect=[0, bottom, 1, 1]).
    """

    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=1, fontsize=11,
               title=title, title_fontsize=11, frameon=False)

    rows = len(labels) + (1 if title else 0)
    return min(0.5, (0.27 * rows + 0.15) / fig.get_figheight())


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
                  smooth       : bool  = True,      # [-] Narrowband only: smooth every observer curve too (False: only the means)
                  smooth_fraction: int = 12,        # [-] Narrowband only: smoothing window, 1/N octave (None: no smoothing at all)
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
                          smooth=smooth,
                          smooth_fraction=smooth_fraction,
                          figsize=figsize,
                          ax_in=axes[i])
        if i % ncols > 0:
            axes[i].set_ylabel("")

    fig.tight_layout()
    pu.set_figure_name(fig, "spectrum", mode, "octave" if octave else "narrowband", ws_suffix(cases))
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

    if ws_suffix(cases):
        # Every panel keeps the legends of its wind speeds (and bands) below it
        fig.tight_layout(rect=[0, 0, 1, 0.95])
    else:
        # Legend handling
        handles, labels = axes[0].get_legend_handles_labels()

        for ax in axes:
            if ax.get_legend() is not None:
                ax.get_legend().remove()

        fig.legend(handles, labels, loc='lower center', bbox_to_anchor = (0.5, 0.02), ncol=len(handles), frameon=False)

        fig.tight_layout(rect=[0, 0.1, 1, 0.95])

    pu.set_figure_name(fig, "polar", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "", ws_suffix(cases))
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

    # References and Ge 2025 data (the wind speeds of --ws have a legend of their own below every panel)
    handles, labels = axes[0].get_legend_handles_labels()

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

    pu.set_figure_name(fig, "decay", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "", ws_suffix(cases))
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
    if ws_suffix(cases):
        fig.tight_layout()      # Room for the wind speed legends below the panels

    pu.set_figure_name(fig, "line", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "", ws_suffix(cases))
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

def compute_point_source(cases           : list  = None,      # [-] List of dictionaries containing simulation data
                         filter_under    : float = None,      # [Hz] Lower frequency cutoff (f >= filter_under)
                         filter_over     : float = None,      # [Hz] Upper frequency cutoff (f <= filter_over)
                         enbw            : float = None,      # [bins] ENBW of the FFT window (None: case["ENBW"] of every file)
                         smooth_fraction : int   = 12,        # [-] Plots only: narrowband curves smoothed over 1/N octave (None: raw)
                         rho_water       : float = 1025.0,    # [kg/m³] Water density (if not in the file)
                         c_water         : float = 1500.0,    # [m/s] Speed of sound in water (if not in the file)
                         verbose         : bool  = True,      # [-] Print every case and the table of broadband totals
                         do_plot         : bool  = True,      # [-] Plot SL spectra, radiation diagnostics and PL of every case
                         figsize         : tuple = (10, 5)):  # [-] Size of the SL and PL figures
    """
    Equivalent point source of every case (plot_single_case.compute_point_source) and their comparison. Every
    case prints its own output (totals and 1/3 octave table), followed by a table with the broadband totals of
    all of them. The figures are those of a single case with every case overlaid in its own color, instead of
    one panel per case, since they share units and the point is to compare their levels:

    - SL spectra: narrowband densities (free field solid, in situ dashed) and 1/3 octave band levels.
    - Radiation diagnostics per 1/3 octave band (cases with sphere): in situ - free, DI, and the dipole and
      monopole shares of the sphere power with the share of |p|^2 they explain.
    - Propagation loss at the spectrum observers: energetic mean of every case, shaded between its observers.

    The near-field limit (kr = 3) of every sphere is a dotted vertical line of its color. Cases without sphere
    nor cylinder data are left out with a message.

    Returns
    -------
    metrics : list of dict (None for the cases left out)
    """

    if not cases:
        print("compute_point_source(): no cases given, skipping")
        return []

    # --- 1. Metrics of every case, printed one after another (the figures below replace theirs) --- #
    metrics = []
    for case in cases:
        if not (case.get("has_sphere", False) or case.get("has_cylinder", False)):
            print(f"compute_point_source(): {case['Name']} has no sphere nor cylinder data, left out of the comparison")
            metrics.append(None)
            continue

        metrics.append(psc.compute_point_source(case=case,
                                                filter_under=filter_under,
                                                filter_over=filter_over,
                                                enbw=enbw,
                                                smooth_fraction=smooth_fraction,
                                                rho_water=rho_water,
                                                c_water=c_water,
                                                verbose=verbose,
                                                do_plot=False))

    pairs = [(case, m) for case, m in zip(cases, metrics) if m is not None]
    if not pairs:
        return metrics

    def near_field_limit(m):     # [Hz] Highest frequency with kr < 3 on the sphere, None without near field
        if "kr" not in m or not np.any(m["kr"] < 3):
            return None
        return m["Freqs"][m["kr"] < 3].max()

    # --- 2. Broadband totals of every case --- #
    if verbose:
        label_w = max(len("Case"), *(len(case["Case_type"]) for case, _ in pairs))
        width   = label_w + 70
        print("\n" + "=" * width)
        print("  EQUIVALENT POINT SOURCE - COMPARISON  |  broadband SL [dB re 1 µPa² m²]")
        print("=" * width)
        print(f"  {'Case':<{label_w}} {'Band [Hz]':>13} {'ENBW':>5} {'SL free':>8} {'SL situ':>8} {'situ-free':>9} "
              f"{'r sphere':>9} {'r cyl':>7}")
        for case, m in pairs:
            band      = f"{m['Freqs'].min():.1f}-{m['Freqs'].max():.1f}"
            sl_free   = f"{m['SL_free_total']:.1f}" if "SL_free_total" in m else "-"
            sl_situ   = f"{m['SL_insitu_total']:.1f}" if "SL_insitu_total" in m else "-"
            situ_free = (f"{10*np.log10(m['W_insitu_total'] / m['W_free_total']):+.1f}"
                         if "W_insitu_total" in m and "W_free_total" in m else "-")
            r_sphere  = f"{m['R_sphere']:.0f} m" if "R_sphere" in m else "-"
            r_cyl     = f"{m['R_cylinder']:.0f} m" if "R_cylinder" in m else "-"
            print(f"  {case['Case_type']:<{label_w}} {band:>13} {m['enbw']:>5g} {sl_free:>8} {sl_situ:>8} "
                  f"{situ_free:>9} {r_sphere:>9} {r_cyl:>7}")

    if not do_plot:
        return metrics

    note = f" (1/{smooth_fraction} oct. smoothed)" if smooth_fraction else ""

    def smooth(f, y):            # Narrowband curve smoothed as in the single-case figure
        return au.smooth_fractional_octave(f, y, smooth_fraction, domain='db') if smooth_fraction else y

    # --- 3. Source level spectra: narrowband densities and band-averaged densities (same units) --- #
    fig, ax = plt.subplots(figsize=figsize)
    handles, labels = [], []
    has_free = has_situ = has_near = False

    for case, m in pairs:
        color  = pu.get_case_color(case)
        f, fc  = m["Freqs"], m["fc"]
        bw     = fc * (2**(1/6) - 2**(-1/6))      # [Hz] 1/3 octave band widths
        totals = []

        if "SL_free_density" in m:
            has_free = True
            ax.semilogx(f, smooth(f, m["SL_free_density"]), '-', color=color)
            ax.semilogx(fc, m["SL_free_band"] - 10*np.log10(bw), 'o', color=color)
            totals.append(f"free {m['SL_free_total']:.1f}")
        if "SL_insitu_density" in m:
            has_situ = True
            ax.semilogx(f, smooth(f, m["SL_insitu_density"]), '--', color=color)
            ax.semilogx(fc, m["SL_insitu_band"] - 10*np.log10(bw), 's', color=color, mfc='none')
            totals.append(f"in situ {m['SL_insitu_total']:.1f}")
        if near_field_limit(m) is not None:
            has_near = True
            ax.axvline(near_field_limit(m), color=color, ls=':', lw=1.2)

        handles.append(Line2D([], [], color=color, lw=2))
        labels.append(f"{case['Case_type']}: " + ", ".join(totals))

    styles = []
    if has_free:
        styles += [Line2D([], [], color='k', ls='-', label=f"Free field{note}"),
                   Line2D([], [], color='k', ls='none', marker='o', label="Free field, 1/3 oct.")]
    if has_situ:
        styles += [Line2D([], [], color='k', ls='--', label=f"In situ{note}"),
                   Line2D([], [], color='k', ls='none', marker='s', mfc='none', label="In situ, 1/3 oct.")]
    if has_near:
        styles += [Line2D([], [], color='k', ls=':', lw=1.2, label=r"Near field limit ($kr = 3$)")]

    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel(r"SL [dB re 1 $\mu$Pa$^2$ m$^2$/Hz]")
    ax.set_title("Equivalent point source")
    ax.grid(True, which='both', ls='--', alpha=0.4)
    ax.legend(handles=styles, fontsize=10)
    bottom = add_case_legend(fig, handles, labels, title=r"Broadband SL [dB re 1 $\mu$Pa$^2$ m$^2$]")
    fig.tight_layout(rect=[0, bottom, 1, 1])
    pu.set_figure_name(fig, "point_source", "SL")

    # --- 4. Radiation diagnostics per 1/3 octave band (cases with sphere) --- #
    sphere_pairs = [(case, m) for case, m in pairs if "DI_band" in m]
    if sphere_pairs:
        show_ratio = any("ratio_insitu_free_band" in m for _, m in sphere_pairs)
        n_rows     = 3 if show_ratio else 2
        fig2, axs  = plt.subplots(n_rows, 1, figsize=(figsize[0], 2.6 * n_rows), sharex=True)
        handles, labels = [], []

        for case, m in sphere_pairs:
            color = pu.get_case_color(case)
            fc    = m["fc"]
            row   = 0
            if show_ratio:
                if "ratio_insitu_free_band" in m:
                    axs[row].semilogx(fc, m["ratio_insitu_free_band"], '-o', color=color)
                row += 1
            axs[row].semilogx(fc, m["DI_band"], '-o', color=color)
            row += 1
            axs[row].semilogx(fc, 100*m["dip_fraction_band"], '-o', color=color)
            axs[row].semilogx(fc, 100*m["mono_fraction_band"], '--s', color=color)
            axs[row].semilogx(fc, 100*m["multipole_explained_band"], ':^', color=color)
            if near_field_limit(m) is not None:
                for a in axs:
                    a.axvline(near_field_limit(m), color=color, ls=':', lw=1.2)

            handles.append(Line2D([], [], color=color, lw=2))
            labels.append(f"{case['Case_type']} (sphere r = {m['R_sphere']:.0f} m)")

        if show_ratio:
            axs[0].axhline(0.0, color='0.5', lw=0.8)
            axs[0].set_ylabel("In situ − free\n[dB]")
        axs[-2].set_ylabel("DI [dB]")
        axs[-1].set_ylabel(r"Share [\%]" if plt.rcParams["text.usetex"] else "Share [%]")    # % starts a LaTeX comment
        axs[-1].legend(handles=[Line2D([], [], color='k', ls='-', marker='o', label="Dipole power"),
                                Line2D([], [], color='k', ls='--', marker='s', label="Monopole power"),
                                Line2D([], [], color='k', ls=':', marker='^', label=r"$|p|^2$ explained")],
                       fontsize=10)
        axs[-1].set_xlabel("1/3 octave band centre frequency [Hz]")
        for a in axs:
            a.grid(True, which='both', ls='--', alpha=0.4)
        axs[0].set_title("Radiation diagnostics")
        bottom = add_case_legend(fig2, handles, labels)
        fig2.tight_layout(rect=[0, bottom, 1, 1])
        pu.set_figure_name(fig2, "point_source", "diagnostics")

    # --- 5. Propagation loss at the spectrum observers --- #
    pl_pairs = [(case, m) for case, m in pairs if "PL_band" in m]
    if pl_pairs:
        fig3, ax3 = plt.subplots(figsize=figsize)
        handles, labels = [], []

        for case, m in pl_pairs:
            color = pu.get_case_color(case)
            fc    = m["fc"]
            with warnings.catch_warnings():     # A band where every observer is NaN stays NaN
                warnings.simplefilter("ignore", category=RuntimeWarning)
                pl_min, pl_max = np.nanmin(m["PL_band"], axis=1), np.nanmax(m["PL_band"], axis=1)
            ax3.semilogx(fc, m["PL_mean_band"], '-o', color=color, lw=2)
            ax3.fill_between(fc, pl_min, pl_max, color=color, alpha=0.15, lw=0)

            source = "free" if m["PL_source"] == "free" else "in situ"
            handles.append(Line2D([], [], color=color, lw=2))
            labels.append(f"{case['Case_type']} (SL {source}, {m['PL_band'].shape[1]} observers)")

        ax3.set_xlabel("1/3 octave band centre frequency [Hz]")
        ax3.set_ylabel(r"PL [dB re 1 m$^2$]")
        ax3.set_title("Propagation loss, PL = SL − SPL")
        ax3.grid(True, which='both', ls='--', alpha=0.4)
        ax3.legend(handles=[Line2D([], [], color='k', lw=2, marker='o', label="Energetic mean of the observers"),
                            Patch(facecolor='k', alpha=0.15, label="Range of the observers")], fontsize=10)
        bottom = add_case_legend(fig3, handles, labels)
        fig3.tight_layout(rect=[0, bottom, 1, 1])
        pu.set_figure_name(fig3, "point_source", "PL")

    return metrics

