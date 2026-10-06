"""
Module: wrappers/farm.py
Description: Wind farm: command line of plotting/farm.py, raw dictionary of readers/farm.py -> Case, and the
             labels of a comparison of farms. A farm has no single axis (origin None: no absorption) and no
             --ws (declined for farms).

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import argparse
import numpy as np
from pathlib import Path
from readers.farm import read_farm
from readers.describe import describe_hdf5
from wrappers.case import solver_part, flatten_nodes, validate_case, apply_comparison_labels
from plot_utils.cli import add_common_args, finish_args
from plot_utils.files import FARM_DATA_DIR, resolve_paths
from plot_utils.colors import get_type_color
from plot_utils.tex import tex_escape
from plot_utils.conversions import get_bin_width


# ---------- COMMAND LINE ---------- #
def parse_farm_args():
    """
    Parses arguments from the terminal.
    """

    parser = argparse.ArgumentParser(description="Plots a wind farm (one .hdf5 file) or compares several of them "
                                                 "(two or more files).")
    parser.add_argument(
        "files", type = str, nargs = "*",
        help = "Wind farm .hdf5 files written by Fortran: one file plots a single farm, several files are compared. "
               "A bare file name is also searched in --data_dir (default: most recent .hdf5 in --data_dir)"
    )
    add_common_args(parser, FARM_DATA_DIR)

    return finish_args(parser.parse_args(), describe_hdf5)


# ---------- CASES ---------- #
def load_farm(path   : Path = None,      # [-] Path to the wind farm .hdf5 file
              verbose: bool = False):    # [-] if True, prints a list of keys available in the file
    """
    Loads a wind farm .hdf5 file coming from Fortran as a Case.

    Returns
    -------
    case : dict
        Case of the contract (wrappers/case.py). The parameters of each turbine are in
        case["specific"]["Turbine_parameters"] (one dict per turbine, as readers.turbine.read_turbine_params).
    """

    return to_case(read_farm(path, verbose=verbose))

def to_case(raw: dict = None):     # [-] Raw dictionary of readers.farm.read_farm
    """
    Translates the datasets of a farm into the Case contract: label = number of turbines and their types,
    nodes = the members of every turbine (one group per turbine), no origin (a farm has no single axis, so
    absorption is off), centre = mean turbine axis, no direction (every turbine has its own, drawn by the
    exclusive plot_structure), color = that of the type of the first turbine.

    Returns
    -------
    case : dict
    """

    turbines = raw["Turbine_parameters"]

    case = {"Path"   : raw["Path"],
            "Name"   : raw["Name"],
            "Emitter": "farm",
            "label"  : get_farm_label(raw)}
    case.update(solver_part(raw))

    case["bin_width"] = get_bin_width(case["Freqs"])

    # Nodes of every turbine: one group per turbine
    case["nodes"], case["node_group"] = flatten_nodes([turb.get("Structure_nodes", np.zeros((0, 3))) for turb in turbines])
    case["group_labels"] = [f"Turbine {i}" for i in range(1, raw["Nturbines"] + 1)]
    case["node_order"]   = None

    # Reference points: no single axis; the centre of the farm is the mean turbine axis
    axes = [turb["AxisPos"] for turb in turbines]
    case["origin"] = None
    case["centre"] = np.asarray(axes)[:, :2].mean(axis=0) if all(axis is not None for axis in axes) else None

    case["direction"]  = None
    case["references"] = None
    case["Color"]      = get_type_color(tex_escape(turbines[0]["case_type"]))
    case["specific"]   = {"Turbine_parameters": turbines,
                          "Nturbines"         : raw["Nturbines"],
                          "Source_offset"     : raw["Source_offset"],
                          "has_nodes"         : raw["has_nodes"]}

    return validate_case(case)

def load_farms(args   : argparse.Namespace = None,   # [-] Parsed command-line arguments
               verbose: bool = False):               # [-] if True, prints a list of keys available in the file
    """
    Loads the farm given in the terminal, or all of them to compare them.

    Returns
    -------
    cases : list of dict
        One Case per file. With several files every case gets a distinct label and color.
    """

    paths = resolve_paths(args.files, args.data_dir)
    cases = [load_farm(path, verbose=verbose or args.verbose) for path in paths]

    if len(cases) > 1:
        set_farm_comparison_labels(cases)

    return [validate_case(case) for case in cases]


# ---------- LABELS ---------- #
def get_farm_label(raw: dict = None):     # [-] Raw dictionary or Case of a wind farm
    """
    Label of a farm: number of turbines and their types, e.g. '5 x DTU10MWFloating' or
    '2 x DTU10MWMonopile + DTU10MWFloating'.

    Returns
    -------
    label : str
    """

    turbines = raw["Turbine_parameters"] if "Turbine_parameters" in raw else raw["specific"]["Turbine_parameters"]
    types    = list(dict.fromkeys(tex_escape(turb["case_type"]) for turb in turbines))
    count    = raw["Nturbines"] if "Nturbines" in raw else raw["specific"]["Nturbines"]

    return f"{count} x {' + '.join(types)}"

def get_farm_ws_label(case: dict = None):     # [-] Case of a wind farm
    """
    Wind speed of a farm: that of its turbines ('11.4 m/s') or the range they span ('8.0-11.4 m/s').

    Returns
    -------
    label : str
        Empty if the file has no wind speed.
    """

    speeds = [turb["WindSpeed"] for turb in case["specific"]["Turbine_parameters"] if np.isfinite(turb["WindSpeed"])]
    if not speeds:
        return ""

    low, high = min(speeds), max(speeds)

    return f"{low:.1f} m/s" if f"{low:.1f}" == f"{high:.1f}" else f"{low:.1f}-{high:.1f} m/s"

def set_farm_comparison_labels(cases: list = None):     # [-] List of Cases of wind farms
    """
    Gives every farm of a comparison a distinct label, used by titles, built from its turbines, wind speed and
    method (e.g. '5 x DTU10MWFloating, 11.4 m/s, ANM24'). The water depth is added if two labels match, and the
    file name if they still do. Colors are assigned by position.
    """

    labels = [", ".join(part for part in (get_farm_label(case), get_farm_ws_label(case), case['Method_tag']) if part)
              for case in cases]

    apply_comparison_labels(cases, labels)
