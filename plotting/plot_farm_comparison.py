"""
Module: plot_farm_comparison.py
Description: Side-by-side comparison of two or more wind farms (one panel per farm, drawn by the functions of
             plot_farm_module, with the helpers of plot_comparison). Farms whose file lacks a step are left out
             of that plot. The equivalent point source overlays the farms instead (same units, the point is
             their levels).

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import plot_utils as pu
import plot_farm_module as pf
import plot_comparison as pc


# ---------- Helpers ---------- #
def label_titles(axes : list = None,     # [-] Axes of the comparison
                 cases: list = None):    # [-] Farm plotted in each panel
    """
    Puts the label of every farm on its own line above the title of its panel ('Wind Farm XY Slice - ...' ->
    'label\\nXY Slice - ...'): comparison labels are long and the titles of neighbouring panels would overlap.
    """

    for ax, case in zip(axes, cases):
        title = ax.get_title().removeprefix("Wind Farm ")
        ax.set_title(case["Case_type"] + ("\n" + title if title else ""))


# ---------- Structural Plots ---------- #
def plot_structure(cases  : list  = None,      # [-] List of dictionaries containing simulation data
                   mode   : str   = '3D',       # [-] Plot mode for the structure representation ('3D', 'xy', 'xz')
                   WindDir: bool  = True,       # [-] Whether to plot or not wind direction
                   figsize: tuple = (8, 6)):    # [-] Size of a single panel
    """
    Plots the structure of every farm side by side (one panel per farm).
    """

    cases = pc.select_cases(cases, "has_nodes", "plot_structure")
    if not cases:
        return None, None

    subplot_kw = {'projection': '3d'} if mode == '3D' else None
    sharey = False if mode == '3D' else True

    fig, axes, ncols = pc.make_panels(len(cases), figsize, subplot_kw=subplot_kw, sharey=sharey)

    for i, case in enumerate(cases):
        pf.plot_structure(case=case,
                          mode=mode,
                          WindDir=WindDir,
                          figsize=figsize,
                          ax_in=axes[i])
        if i % ncols > 0 and mode != '3D':
            axes[i].set_ylabel("")

    label_titles(axes, cases)

    fig.tight_layout()
    pu.set_figure_name(fig, "structure", mode)
    return fig, axes


# ---------- Spectrum Plots ---------- #
def plot_spectrums(cases        : list  = None,      # [-] List of dictionaries containing simulation data
                   mode         : str   = 'SPL',     # [-] Selects Sound Pressure Level or PSD. 'SPL' or 'PSD'
                   absorption   : bool  = False,     # [-] Whether to apply approximate absorption attenuation
                   octave       : bool  = True,      # [-] 1/3 octave or fine resolution
                   filter_under : float = None,      # [Hz] Minimum frequency threshold
                   filter_over  : float = None,      # [Hz] Maximum frequency threshold
                   do_thresholds: bool  = False,     # [-] Whether to plot regulatory/reference thresholds
                   smooth       : bool  = True,      # [-] Narrowband only: smooth every observer curve too (False: only the means)
                   smooth_fraction: int = 12,        # [-] Narrowband only: smoothing window, 1/N octave (None: no smoothing at all)
                   figsize      : tuple = (7, 5)):   # [-] Size of a single panel
    """
    Plots the acoustic spectrum of every farm side by side (one panel per farm).
    """

    cases = pc.select_cases(cases, "has_spectrums", "plot_spectrums")
    if not cases:
        return None, None

    fig, axes, ncols = pc.make_panels(len(cases), figsize, sharey=True)

    for i, case in enumerate(cases):
        pf.plot_spectrums(case=case,
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

    # The lower limit of every panel switched off the autoscale of the shared axis: the top fits every farm
    axes[0].autoscale(axis='y')
    axes[0].set_ylim(bottom=0.)

    label_titles(axes, cases)

    fig.tight_layout()
    pu.set_figure_name(fig, "spectrum", mode, "octave" if octave else "narrowband")
    return fig, axes


# ---------- Polar Plots ---------- #
def plot_polar(cases           : list  = None,         # [-] List of dictionaries containing simulation data
               mode            : str   = 'OASPL',      # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', or 'ABS'
               target_frequency: float = 0.88,         # [Hz] Target frequency for 'SPL', 'PSD' or 'ABS' modes
               absorption      : bool  = False,        # [-] Whether to apply absorption attenuation
               bands           : list  = [10., 100.],  # [Hz] Frequency band limits for OASPL
               figsize         : tuple = (7, 7)):      # [-] Size of a single panel
    """
    Plots polar directivity diagram of every farm side by side (one panel per farm).
    """

    cases = pc.select_cases(cases, "has_polar", "plot_polar")
    if not cases:
        return None, None

    fig, axes, ncols = pc.make_panels(len(cases), figsize, subplot_kw={"projection": "polar"})

    for i, case in enumerate(cases):
        pf.plot_polar(case=case,
                      mode=mode,
                      target_frequency=target_frequency,
                      absorption=absorption,
                      bands=bands,
                      figsize=figsize,
                      ax_in=axes[i])

    label_titles(axes, cases)

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
                  mode            : str   = 'OASPL',   # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                  target_frequency: float = 0.88,      # [Hz] Target frequency for single-frequency modes
                  absorption      : bool  = False,     # [-] Whether to apply absorption attenuation
                  filter_under    : float = None,      # [Hz] Lower frequency cutoff (f >= filter_under)
                  filter_over     : float = None,      # [Hz] Upper frequency cutoff (f <= filter_over)
                  figsize         : tuple = (10, 5),   # [-] Size of a single panel
                  cmap            : str   = 'inferno'):# [-] Contour colormap
    """
    Plots 2D unwrapped cylindrical surface maps of every farm, stacked in a single column.
    """

    cases = pc.select_cases(cases, "has_cylinder", "plot_cylinder")
    if not cases:
        return None, None

    fig, axes, ncols = pc.make_panels(len(cases), figsize, stacked=True, sharex=True)

    global_title = ""

    for i, case in enumerate(cases):
        pf.plot_cylinder(case=case,
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

        # The magnitude is the title of the figure, the radius stays with every panel
        title_lines = axes[i].get_title().removeprefix("Wind Farm ").split("\n")
        if i == 0:
            global_title = title_lines[0]
        axes[i].set_title("\n".join([case["Case_type"]] + title_lines[1:]), fontsize=14)

    fig.suptitle(global_title, fontsize=16, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    pu.set_figure_name(fig, "cylinder", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "")
    return fig, axes


# ---------- Linear Plots ---------- #
def plot_line(cases           : list  = None,      # [-] List of dictionaries containing simulation data
              mode            : str   = 'OASPL',   # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
              target_frequency: float = 0.88,      # [Hz] Target frequency for single-frequency modes
              absorption      : bool  = False,     # [-] Whether to apply absorption attenuation
              filter_under    : float = None,      # [Hz] Lower frequency cutoff (f >= filter_under)
              filter_over     : float = None,      # [Hz] Upper frequency cutoff (f <= filter_over)
              figsize         : tuple = (7, 5)):   # [-] Size of a single panel
    """
    Plots acoustic magnitude along a straight line of every farm side by side (one panel per farm).
    """

    cases = pc.select_cases(cases, "has_line", "plot_line")
    if not cases:
        return None, None

    fig, axes, ncols = pc.make_panels(len(cases), figsize, sharey=True)

    ymin, ymax = float('inf'), float('-inf')
    for i, case in enumerate(cases):
        pf.plot_line(case=case,
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

    label_titles(axes, cases)

    pu.set_figure_name(fig, "line", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "")
    return fig, axes


# ---------- Slice Plots ---------- #
def plot_sliceXY(cases           : list  = None,      # [-] List of dictionaries containing simulation data
                 mode            : str   = 'OASPL',   # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                 target_frequency: float = 0.88,      # [Hz] Target frequency for single-frequency modes
                 absorption      : bool  = False,     # [-] Whether to apply absorption attenuation
                 filter_under    : float = None,      # [Hz] Lower frequency cutoff (f >= filter_under)
                 filter_over     : float = None,      # [Hz] Upper frequency cutoff (f <= filter_over)
                 structure       : bool  = True,      # [-] Whether to plot structure on slice
                 figsize         : tuple = (7, 7),    # [-] Size of a single panel
                 cmap            : str   = 'inferno'):# [-] Contour colormap
    """
    Plots a 2D spatial XY slice of every farm side by side (one panel per farm).
    """

    cases = pc.select_cases(cases, "has_slicexy", "plot_sliceXY")
    if not cases:
        return None, None

    fig, axes, ncols = pc.make_panels(len(cases), figsize)

    for i, case in enumerate(cases):
        pf.plot_sliceXY(case=case,
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
    pc.remove_common_title_part(axes, r'([-\(,\s]*(?:Depth|Z)\s*[:=]\s*[-0-9\.]+\s*(?:\[m\]|m)?[\)]*)')
    label_titles(axes, cases)

    fig.tight_layout()
    pu.set_figure_name(fig, "sliceXY", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "")
    return fig, axes

def plot_sliceXZ(cases           : list  = None,      # [-] List of dictionaries containing simulation data
                 mode            : str   = 'OASPL',   # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                 target_frequency: float = 0.88,      # [Hz] Target frequency for single-frequency modes
                 absorption      : bool  = False,     # [-] Whether to apply absorption attenuation
                 filter_under    : float = None,      # [Hz] Lower frequency cutoff (f >= filter_under)
                 filter_over     : float = None,      # [Hz] Upper frequency cutoff (f <= filter_over)
                 structure       : bool  = True,      # [-] Whether to plot structure on slice
                 figsize         : tuple = (7, 7),    # [-] Size of a single panel
                 cmap            : str   = 'inferno'):# [-] Contour colormap
    """
    Plots a 2D spatial XZ slice of every farm side by side (one panel per farm).
    """

    cases = pc.select_cases(cases, "has_slicexz", "plot_sliceXZ")
    if not cases:
        return None, None

    fig, axes, ncols = pc.make_panels(len(cases), figsize)

    for i, case in enumerate(cases):
        pf.plot_sliceXZ(case=case,
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
    pc.remove_common_title_part(axes, r'([-\(,\s]*(?:Angle|Azimuth|Y)\s*[:=]\s*[-0-9\.]+\s*(?:\[.*?\]|deg|°|m)?[\)]*)')
    label_titles(axes, cases)

    fig.tight_layout()
    pu.set_figure_name(fig, "sliceXZ", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "")
    return fig, axes

def plot_sliceVertical(cases           : list  = None,      # [-] List of dictionaries containing simulation data
                       mode            : str   = 'OASPL',   # [-] Magnitude to plot: 'OASPL', 'SPL', 'PSD', 'SWL', 'ABS', 'REAL', 'IMAG', 'PHASE'
                       target_frequency: float = 0.88,      # [Hz] Target frequency for single-frequency modes
                       absorption      : bool  = False,     # [-] Whether to apply absorption attenuation
                       filter_under    : float = None,      # [Hz] Lower frequency cutoff (f >= filter_under)
                       filter_over     : float = None,      # [Hz] Upper frequency cutoff (f <= filter_over)
                       structure       : bool  = True,      # [-] Whether to plot structure on slice
                       figsize         : tuple = (7, 7),    # [-] Size of a single panel
                       cmap            : str   = 'inferno'):# [-] Contour colormap
    """
    Plots a 2D spatial vertical slice (any azimuth) of every farm side by side (one panel per farm).
    """

    cases = pc.select_cases(cases, "has_slicevertical", "plot_sliceVertical")
    if not cases:
        return None, None

    fig, axes, ncols = pc.make_panels(len(cases), figsize)

    for i, case in enumerate(cases):
        pf.plot_sliceVertical(case=case,
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

    # Delete the plane azimuth and centre if they are the same for all plots
    pc.remove_common_title_part(axes, r'([-\(,\s]*(?:Angle|Azimuth)\s*[:=]\s*[-0-9\.]+\s*(?:\[.*?\]|deg|°)?\s*'
                                      r'(?:\$\|\$\s*Center\s*=\s*\([^)]*\))?[\)]*)')
    label_titles(axes, cases)

    fig.tight_layout()
    pu.set_figure_name(fig, "sliceVertical", mode, f"{target_frequency:.2f}Hz" if mode.upper() != "OASPL" else "")
    return fig, axes


# ---------- Compute metrics ---------- #
def compute_cylinder_metrics(cases       : list  = None,     # [-] List of dictionaries containing simulation data
                             **kwargs):                      # [-] Arguments of plot_farm_module.compute_cylinder_metrics
    """
    Computes and prints the cylinder propagation metrics of every farm.

    Returns
    -------
    metrics : list of dict (None for the farms without cylinder data)
    """

    return [pf.compute_cylinder_metrics(case, **kwargs) for case in cases]

def compute_point_source(cases           : list  = None,      # [-] List of dictionaries containing simulation data
                         **kwargs):                           # [-] Arguments of plot_comparison.compute_point_source
    """
    Equivalent point source of every farm (plot_farm_module.compute_point_source) and their comparison: every
    farm prints its own output, followed by a table with the broadband totals of all of them, and the figures
    overlay the farms in their own color (see plot_comparison.compute_point_source). Farms without cylinder
    data are left out with a message.

    Returns
    -------
    metrics : list of dict (None for the farms left out)
    """

    return pc.compute_point_source([pf.as_single_case(case) for case in cases], **kwargs)
