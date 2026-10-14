"""
Module: wrappers/case.py
Description: The Case contract between the wrappers (one per emitter) and the generic plots. A Case is a plain
             dict built by a wrapper and checked by validate_case at its end:

             Path, Name           file and stem
             Emitter              'turbine' | 'farm' | 'source'
             label                text of the titles (replaced by a distinct label in comparisons)
             Method, Method_name, Method_tag, c_wat, rho_wat, p_ref, Lower_HBC
                                  solver parameters, as readers.solver.read_solver_params
             Depth                water depth [m], read only through seabed_depth(case, x, y)
             water, Alpha_dBkm    water of the seawater absorption (T, S, pH, z) and alpha(f) [dB/km] when the
                                  solver already applied it to the pressures (None otherwise)
             Freqs, ENBW          frequencies [Hz] and window bandwidth [bins] (1.0 without window)
             bin_width            spectral bin width [Hz] for PSD; None for tonal sources
             nodes                (N, 3) points drawn over the maps, NaN rows removed (may be (0, 3))
             node_group, group_labels, node_order
                                  optional: group of every node, label of every group, multipole order
             origin               (2,) reference point of the observer distances; None: no single reference,
                                  absorption is switched off
             centre               (2,) point the spectrum groups its observers by distance to (None if unknown)
             direction, direction_name
                                  direction [deg] of an arrow (e.g. wind) and its label; None: no arrow
             references           external reference curves of the decay plot (None: none)
             family               the case at several values of a parameter (--ws), each member with
                                  family_label, family_title, family_color, family_lw, family_zorder, family_tag
             Color                color of the case (by type, or by position in a comparison)
             specific             emitter-specific data, read only by the exclusive functions
             <steps>              arrays of readers.solver.read_steps (P_polar, Obs_polar, ..., has_<step>)

             The plots read the Case through these keys and the accessors of this module.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import matplotlib.colors as mcolors
from pathlib import Path
from plot_utils.tex import tex_escape
from plot_utils.colors import get_case_colors


EMITTERS = ("turbine", "farm", "source")

# Datasets of every post-processing step (readers.solver.read_steps): has_<step> flag -> keys it carries
STEP_KEYS = {
    "has_spectrums"    : ["P_spectrums", "Obs_spectrums"],
    "has_polar"        : ["P_polar", "R_polar", "Z_polar", "Obs_polar", "Theta_deg_polar", "Center_polar"],
    "has_cylinder"     : ["P_cylinder", "R_cylinder", "Z_cylinder", "Obs_cylinder", "Theta_deg_cylinder", "Center_cylinder"],
    "has_decay"        : ["P_decay", "Z_decay", "Distance_decay", "Distances_decay", "Obs_decay", "Logspace_decay"],
    "has_line"         : ["P_line", "P1_line", "P2_line", "Distance_line", "Distances_line", "Obs_line", "Logspace_line"],
    "has_slicexy"      : ["P_slicexy", "Obs_slicexy", "X_slicexy", "Y_slicexy", "Z_slicexy", "Center_slicexy",
                          "Nx_slicexy", "Ny_slicexy"],
    "has_slicexz"      : ["P_slicexz", "Obs_slicexz", "X_slicexz", "Z_slicexz", "Y_slicexz", "Nx_slicexz", "Nz_slicexz"],
    "has_slicevertical": ["P_sliceV", "Coords_sliceV", "U_sliceV", "Z_sliceV", "Azimuth_sliceV", "Center_sliceV"],
    "has_sphere"       : ["P_sphere", "Obs_sphere", "Theta_sphere", "Z_sphere", "R_sphere", "Center_sphere", "dA_sphere",
                          "N_theta", "Nz_sphere"],
}
OPTIONAL_STEP_KEYS = {"has_cylinder": ["dA_cylinder"]}     # Older files lack it: plot_cylinder computes it

# Solver parameters, the same for every emitter (the method-specific ones only when the file has them)
SOLVER_KEYS   = ["Method", "Method_name", "Method_tag", "ENBW", "c_wat", "rho_wat", "p_ref", "Lower_HBC", "Depth", "Freqs",
                 "water", "Alpha_dBkm"]
METHOD_KEYS   = ["Nmodes", "H", "N_images"]

REQUIRED_KEYS = ["Path", "Name", "Emitter", "label", "bin_width", "nodes", "origin", "centre", "direction", "specific"] \
                + SOLVER_KEYS + list(STEP_KEYS)
FAMILY_KEYS   = ["family_label", "family_title", "family_color", "family_lw", "family_zorder", "family_tag"]


# ---------- CONSTRUCTION ---------- #
def solver_part(raw: dict = None):     # [-] Raw dictionary of a reader
    """
    The part of a Case every emitter shares: solver parameters and post-processing steps, copied from the raw
    dictionary of a reader as they are (same names, same arrays: nothing is duplicated).

    Returns
    -------
    case : dict
    """

    case = {key: raw[key] for key in SOLVER_KEYS + METHOD_KEYS if key in raw}

    for flag, keys in STEP_KEYS.items():
        case[flag] = raw[flag]
        if raw[flag]:
            case.update({key: raw[key] for key in keys + OPTIONAL_STEP_KEYS.get(flag, []) if key in raw})

    return case

def flatten_nodes(nodes_list: list = None):     # [-] Arrays of nodes (..., 3), one per group (member, turbine, ...)
    """
    Flat list of the nodes drawn over the maps: every array reshaped to (n, 3), the rows with NaN dropped
    (NaN-padded members, e.g. SAITEC's rigid links), and the group of every node.

    Returns
    -------
    nodes, node_group : np.ndarray (N, 3), np.ndarray (N,) of int
    """

    nodes, groups = [np.zeros((0, 3))], [np.zeros(0, dtype=int)]
    for g, array in enumerate(nodes_list):
        xyz  = np.asarray(array, dtype=float).reshape(-1, 3)
        keep = np.all(np.isfinite(xyz), axis=1)
        nodes.append(xyz[keep])
        groups.append(np.full(np.count_nonzero(keep), g, dtype=int))

    return np.concatenate(nodes), np.concatenate(groups)

def validate_case(case: dict = None):     # [-] Case built by a wrapper
    """
    Checks a Case against the contract (see the description of this module) and raises a single ValueError
    listing every missing or malformed key.

    Returns
    -------
    case : dict
        The same Case, so a wrapper can end with 'return validate_case(case)'.
    """

    problems = []

    def check(condition: bool, message: str):
        if not condition:
            problems.append(message)

    def is_float(value):
        return isinstance(value, (float, np.floating)) and not isinstance(value, bool)

    def is_point(value):
        return isinstance(value, np.ndarray) and value.shape == (2,) and np.all(np.isfinite(value))

    for key in REQUIRED_KEYS:
        check(key in case, f"missing key '{key}'")
    if problems:
        raise ValueError(f"wrappers.case.validate_case(): {case.get('Name', '?')} breaks the Case contract:\n  - "
                         + "\n  - ".join(problems))

    check(isinstance(case["Path"], Path), "Path must be a pathlib.Path")
    check(isinstance(case["Name"], str), "Name must be a str")
    check(case["Emitter"] in EMITTERS, f"Emitter must be one of {EMITTERS}, not {case['Emitter']!r}")
    for key in ("label", "Method", "Method_name", "Method_tag"):
        check(isinstance(case[key], str), f"{key} must be a str")
    for key in ("c_wat", "rho_wat", "p_ref", "Lower_HBC", "Depth", "ENBW"):
        check(is_float(case[key]) and np.isfinite(case[key]), f"{key} must be a finite float")

    freqs = case["Freqs"]
    check(isinstance(freqs, np.ndarray) and freqs.ndim == 1 and freqs.size > 0, "Freqs must be a non-empty 1D array")
    water = case["water"]
    check(isinstance(water, dict) and all(is_float(water.get(key)) for key in ("T", "S", "pH", "z")),
          "water must be a dict of floats T, S, pH, z")
    alpha = case["Alpha_dBkm"]
    check(alpha is None or (isinstance(alpha, np.ndarray) and isinstance(freqs, np.ndarray) and alpha.shape == freqs.shape),
          "Alpha_dBkm must be None or an array shaped as Freqs")
    check(case["bin_width"] is None or (is_float(case["bin_width"]) and case["bin_width"] > 0),
          "bin_width must be a positive float or None")

    # Nodes and their optional descriptions
    nodes = case["nodes"]
    n     = len(nodes) if isinstance(nodes, np.ndarray) else 0
    check(isinstance(nodes, np.ndarray) and nodes.ndim == 2 and nodes.shape[1] == 3 and np.all(np.isfinite(nodes)),
          "nodes must be a (N, 3) array without NaN")
    for key in ("node_group", "node_order"):
        value = case.get(key)
        check(value is None or (isinstance(value, np.ndarray) and value.shape == (n,) and value.dtype.kind in "iu"),
              f"{key} must be None or a (N,) int array")
    labels = case.get("group_labels")
    check(labels is None or (isinstance(labels, list) and all(isinstance(label, str) for label in labels)),
          "group_labels must be None or a list of str")

    # Reference points and direction
    check(case["origin"] is None or is_point(case["origin"]), "origin must be None or a finite (2,) array")
    check(case["centre"] is None or is_point(case["centre"]), "centre must be None or a finite (2,) array")
    if case["direction"] is not None:
        check(is_float(case["direction"]), "direction must be None or a float [deg]")
        check(isinstance(case.get("direction_name"), str), "direction_name must be a str when direction is given")
    check(case.get("references") is None or isinstance(case["references"], list), "references must be None or a list")

    # Family (--ws): the case itself is one of its members
    if "family" in case:
        family = case["family"]
        check(isinstance(family, list) and any(member is case for member in family), "family must be a list holding the case")
        for member in family if isinstance(family, list) else []:
            missing = [key for key in FAMILY_KEYS if key not in member]
            check(not missing, f"family member {member.get('Name', '?')} lacks {missing}")

    if "Color" in case:
        check(mcolors.is_color_like(case["Color"]), f"Color {case['Color']!r} is not a color")
    check(isinstance(case["specific"], dict), "specific must be a dict")

    # Steps: every flag a bool, every array of a present step there
    for flag, keys in STEP_KEYS.items():
        check(isinstance(case[flag], bool), f"{flag} must be a bool")
        if case[flag]:
            missing = [key for key in keys if key not in case]
            check(not missing, f"{flag} is True but {missing} are missing")

    if problems:
        raise ValueError(f"wrappers.case.validate_case(): {case['Name']} breaks the Case contract:\n  - "
                         + "\n  - ".join(problems))

    return case


# ---------- ACCESSORS ---------- #
def seabed_depth(case: dict  = None,      # [-] Case
                 x   : float = None,      # [m] Horizontal position (None: anywhere)
                 y   : float = None):     # [m] Horizontal position (None: anywhere)
    """
    Water depth at (x, y). Constant today (-Lower_HBC of the solver); with bathymetry it becomes a function of
    the position, so every reader of the depth goes through here.

    Returns
    -------
    depth : float
        [m]
    """

    return case["Depth"]

def emitter_label(case: dict = None):     # [-] Case
    """
    Label of the case in titles and prints (its type, or its distinct label in a comparison).

    Returns
    -------
    label : str
    """

    return case["label"]

def case_color(case: dict = None):     # [-] Case
    """
    Color of the case: that of its type, or that of its position in a comparison.

    Returns
    -------
    color : matplotlib color
    """

    return case["Color"]

def get_family(case: dict = None):     # [-] Case
    """
    Cases drawn by a line plot: with a family (e.g. --ws) one per member, otherwise the case alone.

    Returns
    -------
    cases, family_mode : list of dict, bool
        The cases to draw and whether they are the members of a family.
    """

    if "family" in case:
        return case["family"], True

    return [case], False

def family_tag(cases: list = None):     # [-] Cases of a plot
    """
    Word that tells the figures drawn with a family (e.g. 'allws', every wind speed of --ws) from the plain
    ones in the name they are saved with.

    Returns
    -------
    tag : str
        Empty if no case has a family.
    """

    for case in cases:
        if "family" in case:
            return case["family_tag"]

    return ""

def apply_comparison_labels(cases : list = None,     # [-] Cases of a comparison
                            labels: list = None):    # [-] One label per case
    """
    Stores the labels of a comparison in the cases (label) with a color each, by position. The water depth is
    added to the labels that match, and the file name to those that still do.
    """

    if len(set(labels)) < len(labels):
        labels = [f"{label}, H = {seabed_depth(case, *point(case)):.0f} m" for label, case in zip(labels, cases)]
    if len(set(labels)) < len(labels):
        labels = [f"{label} ({tex_escape(case['Name'])})" for label, case in zip(labels, cases)]

    for case, label, color in zip(cases, labels, get_case_colors(len(cases))):
        case["label"] = label
        case["Color"] = color

def point(case: dict = None):     # [-] Case
    """
    Horizontal position that stands for the case (its centre), where its depth is read.

    Returns
    -------
    x, y : float or None
    """

    return (None, None) if case["centre"] is None else tuple(case["centre"])
