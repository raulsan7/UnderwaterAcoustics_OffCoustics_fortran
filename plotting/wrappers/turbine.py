"""
Module: wrappers/turbine.py
Description: Single wind turbine: command line of plotting/turbine.py, raw dictionary of readers/turbine.py ->
             Case, labels of a comparison of turbines and the wind speeds of --ws (the family of a case: the
             files of the same folder whose name only differs in the number after 'ws').

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from readers.turbine import read_turbine
from readers.solver import WS_STEPS
from readers.describe import describe_hdf5
from wrappers.case import solver_part, flatten_nodes, validate_case, apply_comparison_labels
from plot_utils.cli import add_common_args, finish_args
from plot_utils.files import TURBINE_DATA_DIR, resolve_paths
from plot_utils.style import AcademicStyle
from plot_utils.colors import get_type_color
from plot_utils.tex import tex_escape
from plot_utils.conversions import get_bin_width

# --ws: the wind speeds of a case are found by the number after "ws" in the file name
WS_REGEX = r"ws(\d+(?:\.\d+)?)"

# Turbine data kept apart from the contract, for the exclusive functions (structure, sphere metrics)
SPECIFIC = ["case_type", "WindSpeed", "WindDir", "AxisPos", "BariPos", "Structure_nodes", "Nm", "Nn",
            "D", "rho_mat", "wet_area", "col_D", "XsecA", "XsecB", "t", "col_members", "pon_members",
            "Diam_y", "Diam_z", "float_members", "joint_members"]


# ---------- COMMAND LINE ---------- #
def parse_turbine_args():
    """
    Parses arguments from the terminal.
    """

    parser = argparse.ArgumentParser(description="Plots a single wind turbine (one .hdf5 file) or compares "
                                                 "several of them (two or more files).")
    parser.add_argument(
        "files", type = str, nargs = "*",
        help = "Turbine .hdf5 files written by Fortran: one file plots a single case, several files are compared. "
               "A bare file name is also searched in --data_dir (default: most recent .hdf5 in --data_dir)"
    )
    parser.add_argument(
        "--ws", action = "store_true",
        help = "Line plots (spectrum, polar, decay, line) show every wind speed of each case: the files that only "
               "differ in the ws<speed> of their name are searched next to the given one. Maps and metrics "
               "keep using the given file"
    )
    add_common_args(parser, TURBINE_DATA_DIR)

    return finish_args(parser.parse_args(), describe_hdf5)


# ---------- CASES ---------- #
def load_turbine(path   : Path = None,      # [-] Path to the single-turbine .hdf5 file
                 verbose: bool = False,     # [-] if True, prints a list of keys available in the file
                 steps  : list = None):     # [-] Post-processing steps to read (keys of STEPS). None: all of them
    """
    Loads a single-turbine .hdf5 file coming from Fortran as a Case.

    Returns
    -------
    case : dict
        Case of the contract (wrappers/case.py). Every post-processing step has a has_<step> flag, since any of
        them may be missing (files from older versions of the Fortran code or from interrupted runs).
    """

    return to_case(read_turbine(path, verbose=verbose, steps=steps))

def to_case(raw: dict = None):     # [-] Raw dictionary of readers.turbine.read_turbine
    """
    Translates the datasets of a turbine into the Case contract: label = case type, nodes = every member of
    Structure_nodes (one group per member), origin and centre = turbine axis, direction = wind direction,
    references = measured turbine noise (Ge 2025) of the decay plot.

    Returns
    -------
    case : dict
    """

    case = {"Path"   : raw["Path"],
            "Name"   : raw["Name"],
            "Emitter": "turbine",
            "label"  : tex_escape(raw["case_type"])}
    case.update(solver_part(raw))

    case["bin_width"] = get_bin_width(case["Freqs"])

    # Structure nodes shape(Nm, Nn, 3): one group per member
    members = raw.get("Structure_nodes")
    members = [] if members is None else (list(members) if members.ndim == 3 else [members])
    case["nodes"], case["node_group"] = flatten_nodes(members)
    case["group_labels"] = None
    case["node_order"]   = None

    # Reference points: the turbine axis
    axis = raw["AxisPos"]
    case["origin"] = case["centre"] = axis[:2] if axis is not None and np.all(np.isfinite(axis[:2])) else None

    case["direction"]      = raw["WindDir"]
    case["direction_name"] = "WindDir"
    case["references"]     = get_turbine_references()
    case["Color"]          = get_type_color(case["label"])
    case["specific"]       = {key: raw[key] for key in SPECIFIC if key in raw}

    return validate_case(case)

def load_turbines(args   : argparse.Namespace = None,   # [-] Parsed command-line arguments
                  verbose: bool = False):               # [-] if True, prints a list of keys available in the file
    """
    Loads the case given in the terminal, or all of them to compare them. With --ws every case also gets the
    other wind speeds of its file name (its family, see load_ws_cases).

    Returns
    -------
    cases : list of dict
        One Case per file. With several files every case gets a distinct label and color.
    """

    paths = resolve_paths(args.files, args.data_dir)
    if args.ws:
        paths = drop_ws_twins(paths)

    cases = [load_turbine(path, verbose=verbose or args.verbose) for path in paths]

    if args.ws:
        for case in cases:
            load_ws_cases(case, verbose=verbose or args.verbose)
        set_ws_styles(cases)

    if len(cases) > 1:
        set_comparison_labels(cases)

    return [validate_case(case) for case in cases]

def get_turbine_references():
    """
    Measured noise of operating offshore wind turbines (Ge 2025), drawn by plot_distance_decay(turbine_data=True).

    Returns
    -------
    references : list of (data, label, marker, color)
        data shape(N, 2): distance [m] and level [dB].
    """

    from utils.EnvironmentalData import Turbine_noise
    Tripile, Gravity_based, Monopile, Jacket, Suction_bucket = Turbine_noise()
    TURBINE_STYLES = [
        (Tripile,        "Tripile (Ge 2025)",        "s", "gold"      ),
        (Gravity_based,  "Gravity based (Ge 2025)",  "o", "limegreen" ),
        (Monopile,       "Monopile (Ge 2025)",       "^", "darkred"   ),
        (Jacket,         "Jacket (Ge 2025)",         "v", "steelblue" ),
        (Suction_bucket, "Suction bucket (Ge 2025)", "D", "silver"    ),
    ]

    return TURBINE_STYLES


# ---------- COMPARISON ---------- #
def set_comparison_labels(cases: list = None):     # [-] List of Cases
    """
    Gives every case of a comparison a distinct label, used by titles, built from its case type, wind speed
    and method (e.g. 'DTU10MWFloating, 11.4 m/s, ANM231'; without the wind speed with --ws, where every panel
    shows all of them). The water depth is added if two labels match, and the file name if they still do.
    Colors are assigned by position.
    """

    labels = [", ".join([tex_escape(case["specific"]["case_type"])] +
                        ([] if "family" in case else [get_ws_label(case)]) +
                        [case['Method_tag']]) for case in cases]

    apply_comparison_labels(cases, labels)


# ---------- WIND SPEEDS (--ws) ---------- #
def get_ws_label(case: dict = None):     # [-] Case of a turbine
    return f"{case['specific']['WindSpeed']:.1f} m/s"

def find_ws_paths(path: Path = None):     # [-] Path to a single-turbine .hdf5 file
    """
    Finds the files of the same case at other wind speeds: those of the folder whose name only differs from
    the given one in the number after 'ws' (plot_mn_IM30_ws11.4.hdf5 -> ..._ws4.5, ..._ws8.0, ...). The name is
    matched whole, so plot_fl_IM30_ws11.4.hdf5 does not collect the ..._ws<speed>_shallow.hdf5 files.

    Returns
    -------
    found : list of (float, Path)
        Wind speed of the name and path of every file, the given one included, by increasing speed. Empty if
        the name has no 'ws<speed>'.
    """

    path    = Path(path)
    matches = list(re.finditer(WS_REGEX, path.name))
    if not matches:
        return []

    last    = matches[-1]
    pattern = re.compile(re.escape(path.name[:last.start()]) + WS_REGEX + re.escape(path.name[last.end():]))

    found = []
    for other in path.parent.glob("*" + path.suffix):
        match = pattern.fullmatch(other.name)
        if match:
            found.append((float(match.group(1)), other))

    return sorted(found)

def drop_ws_twins(paths: list = None):     # [-] Paths given in the terminal
    """
    With --ws every file stands for all its wind speeds: a file whose wind speeds were already requested by
    another one is left out (with a message), so its panel is not repeated.

    Returns
    -------
    kept : list of Path
    """

    kept, families = [], []
    for path in paths:
        twin = next((other for other, family in zip(kept, families) if path.resolve() in family), None)
        if twin is not None:
            print(f"wrappers.turbine.drop_ws_twins(): {path.name} is another wind speed of {twin.name}, loaded only once")
            continue

        kept.append(path)
        families.append({other.resolve() for _, other in find_ws_paths(path)})

    return kept

def load_ws_cases(case   : dict = None,     # [-] Case of a turbine
                  verbose: bool = False):   # [-] if True, prints a list of keys available in the files
    """
    Adds the family of the case (case["family"]): the case at every wind speed found next to its file, by
    increasing speed. The case itself is one of them, so it keeps every step; the other wind speeds only read
    the steps drawn as lines (WS_STEPS: a few tens of MB instead of the GB of a whole file). The line plots
    draw one curve per member. A name without 'ws<speed>' leaves the case as it is.
    """

    found = find_ws_paths(case["Path"])
    if not found:
        print(f"wrappers.turbine.load_ws_cases(): no ws<speed> in the name of {case['Name']}, --ws does not apply to it")
        return

    members = []
    for ws, path in found:
        if path.resolve() == case["Path"].resolve():
            other = case
        else:
            other = load_turbine(path, verbose=verbose, steps=WS_STEPS)

        if not np.isfinite(other["specific"]["WindSpeed"]):     # Files without the WindSpeed dataset: the one of the name
            other["specific"]["WindSpeed"] = ws
        members.append(other)

    case["family"] = sorted(members, key=lambda other: other["specific"]["WindSpeed"])
    for other in case["family"]:
        other["family_label"] = get_ws_label(other)
        other["family_title"] = "Wind speed"
        other["family_tag"]   = "allws"

    print(f"{case['Name']}: wind speeds {', '.join(other['family_label'] for other in case['family'])}")

def set_ws_styles(cases: list = None):     # [-] List of Cases
    """
    Gives every wind speed of the plot its own color (family_color, from the colormap AcademicStyle.WS_CMAP in
    increasing speed), line width (family_lw) and drawing order (family_zorder): the faster the wind, the
    thicker the line and the further down it is drawn, so curves that coincide (e.g. every speed above rated)
    stay visible as concentric lines. The speeds of every case are ranked together, so the same speed has the
    same style in every panel.
    """

    speeds = sorted({other["specific"]["WindSpeed"] for case in cases for other in case.get("family", [])})
    if not speeds:
        return

    n       = len(speeds)
    colors  = plt.get_cmap(AcademicStyle.WS_CMAP)(np.linspace(*AcademicStyle.WS_CMAP_RANGE, n))
    widths  = np.linspace(*AcademicStyle.WS_LW_RANGE, n)

    for case in cases:
        for other in case.get("family", []):
            rank = speeds.index(other["specific"]["WindSpeed"])
            other["family_color"]  = tuple(colors[rank])
            other["family_lw"]     = float(widths[rank])
            other["family_zorder"] = 2 + (n - 1 - rank) / n
