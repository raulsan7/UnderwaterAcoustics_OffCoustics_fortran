"""
Module: plot_utils/drawing.py
Description: Drawing helpers: direction arrow (e.g. wind direction) on 2D or 3D axes, and the legends of a plot
             that draws a family of cases (e.g. every wind speed of --ws) below its axes.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def draw_direction_arrow(ax,                                # [-] Target axis. 3D vs 2D is auto-detected.
                          x              : float,           # [m] Tip location (in the horizontal plane).
                          y              : float,           # [m] Tip location (in the horizontal plane).
                          direction_deg  : float,           # [deg] Arrow direction
                          length         : float,           # [m] Arrow length. Must be > 0.
                          z              : float = 0.0,     # [m] Tip Z coordinate. Ignored for 2D axes.
                          color          : str   = 'red',   # [-] Arrow color
                          label          : str   = None,    # [.] Arrow label
                          head_frac      : float = 0.25,    # [-] Fraction of 'length' used for the arrowhead ignored in 2D
                          head_width_frac: float = 0.25,    # [-] Arrowhead half-width as a fraction of head length ignored in 2D
                          linewidth      : float = 1.8):    # [-] Arrow line width
    """
    Draws a direction arrow (e.g. wind direction) on a 2D or 3D matplotlib axis.
    The arrow tip is placed at (x, y[, z]); it points backward from there
    along 'direction_deg' (degrees, standard math convention: 0=+x, 90=+y).

    Returns
    -------
    The created artist (FancyArrow for 2D, Poly3DCollection for 3D).
    """

    if length is None or length <= 0:
        raise ValueError(f"draw_direction_arrow(): length must be > 0, got {length!r}")

    theta = np.deg2rad(direction_deg)
    u, v = np.cos(theta), np.sin(theta)

    is_3d = isinstance(ax, Axes3D)

    if not is_3d:
        # 2D case: matplotlib's arrow/annotate already draws a proper filled head.
        tail_x = x - length * u
        tail_y = y - length * v
        arrow = ax.annotate(
            '', xy=(x, y), xytext=(tail_x, tail_y),
            arrowprops=dict(arrowstyle='-|>', color=color, linewidth=linewidth,
                             mutation_scale=15),
        )
        if label is not None:
            # dummy line for legend, since annotate() doesn't register one
            ax.plot([], [], color=color, linewidth=linewidth, label=label)
        return arrow

    # 3D case: draw shaft + solid triangular head manually.
    tip = np.array([x, y, z])
    tail = tip - length * np.array([u, v, 0])

    ax.plot([tail[0], tip[0]], [tail[1], tip[1]], [tail[2], tip[2]],
            color=color, linewidth=linewidth)

    head_len = length * head_frac
    head_width = head_len * head_width_frac
    perp = np.array([-v, u, 0])  # in-plane perpendicular to arrow direction
    base_center = tip - head_len * np.array([u, v, 0])
    p1 = tip
    p2 = base_center + head_width * perp
    p3 = base_center - head_width * perp

    head = Poly3DCollection([[p1, p2, p3]], color=color)
    ax.add_collection3d(head)

    if label is not None:
        ax.plot([], [], color=color, linewidth=linewidth, label=label)

    return head

def add_family_legends(ax,                             # [-] Axes of the plot
                       family_handles: list = None,    # [-] One handle per member of the family (color)
                       family_labels : list = None,    # [-] Their labels (e.g. '11.4 m/s')
                       style_handles : list = None,    # [-] Handles of the second dimension of the plot (observers, bands...)
                       style_labels  : list = None,    # [-] Their labels
                       style_title   : str  = None,    # [-] Title of the second legend
                       ncol_style    : int  = 1,       # [-] Columns of the second legend
                       offset        : float = 0.68,   # [in] Distance of the legends below the axes
                       family_title  : str  = None):   # [-] Title of the family legend (e.g. 'Wind speed')
    """
    Legends of a plot that draws a family of cases (e.g. every wind speed of --ws), below its axes: the colors
    are the members of the family and, when the plot needs a second dimension, its line styles (black) are
    told apart by a second legend. Both sit side by side, centred under the axes, the family on the left. The
    offset is absolute, so it clears the tick labels and the label of the x axis whatever the height of the
    axes. Axes legends are part of the layout of fig.tight_layout(), which leaves them room.
    """

    kw = dict(frameon=False, title_fontsize=plt.rcParams["legend.fontsize"], bbox_to_anchor=(0.5, 0.0),
              bbox_transform=ax.transAxes + ScaledTranslation(0.0, -offset, ax.figure.dpi_scale_trans))

    if not style_handles:
        ax.legend(family_handles, family_labels, title=family_title, ncol=min(len(family_handles), 5), loc="upper center", **kw)
        return

    first = ax.legend(family_handles, family_labels, title=family_title, ncol=1 if len(family_handles) <= 2 else 2,
                      loc="upper right", **kw)
    ax.add_artist(first)        # The second ax.legend() would replace it
    ax.legend(style_handles, style_labels, title=style_title, ncol=ncol_style, loc="upper left", **kw)
