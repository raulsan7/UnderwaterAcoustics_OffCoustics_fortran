"""
Module: plots/specific/structure.py
Description: Structure of a single wind turbine or of every turbine of a farm, drawn with its members (lines
             between the nodes of every member, case["specific"]) and the wind direction arrows: 3D view, top
             view (XY) or side view (XZ). Other emitters have no structure.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import matplotlib.pyplot as plt
from wrappers.case import emitter_label, seabed_depth, point
from plot_utils.drawing import draw_direction_arrow
from plot_utils.saving import set_figure_name


def has_structure(case: dict = None):     # [-] Case
    """
    Whether the case has the members of a structure to draw (every turbine of a farm needs them).

    Returns
    -------
    has : bool
    """

    if case["Emitter"] == "turbine":
        return "Structure_nodes" in case["specific"]
    if case["Emitter"] == "farm":
        return case["specific"]["has_nodes"]

    return False

def plot_structure(case   : dict  = None,       # [-] Case of a turbine or of a farm
                   mode   : str   = '3D',       # [-] Plot mode for the structure representation
                   WindDir: bool  = True,       # [-] Wheter to plot or not wind direction
                   figsize: tuple = (8,6),      # [-] Figure size
                   ax_in  : plt.Axes = None):   # [-] Optional external axes for subplots
    """
    Plots the structure of a turbine (plot_turbine_structure) or of every turbine of a farm
    (plot_farm_structure).

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The generated figure and axes objects for further editing.
    """

    if case["Emitter"] == "farm":
        return plot_farm_structure(case, mode, WindDir, figsize, ax_in)

    return plot_turbine_structure(case, mode, WindDir, figsize, ax_in)

def plot_turbine_structure(case   : dict  = None,       # [-] Case of a turbine
                           mode   : str   = '3D',       # [-] Plot mode for the structure representation
                           WindDir: bool  = True,       # [-] Wheter to plot or not wind direction
                           figsize: tuple = (8,6),      # [-] Figure size
                           ax_in  : plt.Axes = None):   # [-] Optional external axes for subplots
    """
    Plots the structure nodes in shape (Nmembers, Nnodes, 3), every member as a line between its nodes.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The generated figure and axes objects for further editing.
    """

    if not has_structure(case):
        print("plot_structure(): no structure data, skipping")
        return None, None

    if mode not in ('3D', 'xy', 'xz'):
        raise ValueError("mode must be '3D', 'xy', or 'xz'")

    # Extract nodes, depth and wind direction (arrow only when the case has a direction)
    nodes     = case["specific"]["Structure_nodes"]
    Depth     = seabed_depth(case, *point(case))
    tag       = emitter_label(case)
    direction = case["direction"]
    WindDir   = WindDir and direction is not None
    arrow     = f'{case.get("direction_name")}: {direction}°'

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

            draw_direction_arrow(ax, max_x, max_y, direction,
                                 length=arrow_length, z=max_z,
                                 color='red',
                                 label=arrow)
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

                draw_direction_arrow(ax, max_x, max_y, direction,
                      length=arrow_length, color='red',
                      label=arrow)

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

                draw_direction_arrow(ax, max_x, max_z, direction,
                      length=arrow_length, color='red',
                      label=arrow)

            ax.set_xlabel('x [m]')
            ax.set_ylabel('z [m]')
            ax.set_title(f"{tag} - Side View (XZ)")
            ax.set_xlim(xmin - xpad, xmax + xpad)
            ax.set_ylim(zmin - zpad, zmax + zpad)
            ax.legend()

    if ax_in is None:
        fig.tight_layout()
        set_figure_name(fig, "structure", mode)

    return fig, ax

def plot_farm_structure(case   : dict = None,        # [-] Case of a wind farm
                        mode   : str = '3D',         # [-] Plot mode for the structure representation
                        WindDir: bool = True,        # [-] Wheter to plot or not wind direction
                        figsize: tuple = (8,6),      # [-] Figure size
                        ax_in  : plt.Axes = None):   # [-] optional external axes for subplots
    """
    Plots the structure nodes in shape(Nmembers, Nnodes, 3) of every turbine of the farm, with the wind
    direction arrow of each turbine.

    Returns
    -------
    fig, ax: matplotlib.figure.Figure, matplotlib.axes.Axes
    The generated figure and axes objet for further editing.
    """

    if not has_structure(case):
        print("plot_structure(): no structure data, skipping")
        return None, None

    if mode not in ("3D", "xy", "xz"):
        raise ValueError("plot_structure(): mode must be '3D', 'xy', 'xz'")

    # Extract nodes and depth. Turbines of different types have different (Nmembers, Nnodes), so the
    # nodes are kept as a list with one array shape(Nmembers, Nnodes, 3) per turbine
    turbines = case["specific"]["Turbine_parameters"]
    Nturb    = case["specific"]["Nturbines"]
    nodes    = [turb["Structure_nodes"] for turb in turbines]
    Depth    = seabed_depth(case, *point(case))
    tag      = emitter_label(case)
    Wind_dir = np.asarray([turb["WindDir"] for turb in turbines])

    # Every node of the farm, to compute extents. SAITEC's rigid members are NaN-padded: NaN-aware functions
    all_nodes = np.concatenate([turb.reshape(-1, 3) for turb in nodes])

    # Find node with maximum Z coordinate to place WindDir arrow
    max_z_nodes = np.zeros((Nturb, 3))
    for i in range(Nturb):
        xyz = nodes[i].reshape(-1,3)
        max_z_nodes[i] = xyz[np.nanargmax(xyz[:,2]),:]

    max_x, max_y, max_z = max_z_nodes[:,0], max_z_nodes[:,1], max_z_nodes[:,2]


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
                draw_direction_arrow(ax, max_x[i], max_y[i], Wind_dir[i],
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
        ax.set_title(tag)
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

                    draw_direction_arrow(ax, max_x[i], max_y[i], Wind_dir[i],
                                         length=arrow_length, z=max_z[i], color = "red",
                                         label=f"WindDir {Wind_dir[i]}")

            ax.set_xlabel('x [m]')
            ax.set_ylabel('y [m]')
            ax.set_title(f"{tag} - Top View (XY)")
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
                        draw_direction_arrow(ax, max_x[i], max_z[i], Wind_dir[i],
                                             length=arrow_length, color = "red",
                                             label=f"WindDir {Wind_dir[i]}")

            ax.set_xlabel('x [m]')
            ax.set_ylabel('z [m]')
            ax.set_title(f"{tag} - Side View (XZ)")
            ax.set_xlim(xmin - xpad, xmax + xpad)
            ax.set_ylim(zmin - zpad, zmax + zpad)
            ax.legend()

    if ax_in is None:
        set_figure_name(fig, "structure", mode)

    return fig, ax
