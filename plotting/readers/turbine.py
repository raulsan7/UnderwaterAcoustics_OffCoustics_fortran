"""
Module: readers/turbine.py
Description: Reader of the single wind turbine files (FORTRAN/turbine_acoustic_data): generic turbine datasets
             (WindTurbine.save_turbine_params), those of its type (save_specific_params of
             FORTRAN/src/emitters/TurbineTypes/*.f90) and the steps of the solver. The per-turbine datasets of
             a farm are read with the same function (suffix "_i").

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import h5py
from pathlib import Path
from readers.hdf5_io import read_ds, read_array, read_float, print_keys, detect_emitter, wrong_emitter
from readers.solver import read_solver_params, read_steps


# Type-specific turbine parameters (save_specific_params of FORTRAN/src/emitters/TurbineTypes/*.f90)
SPECIFIC_KEYS = ["D", "rho_mat", "wet_area",                                      # DTU10MWMonopile (+ the rest)
                 "col_D", "XsecA", "XsecB", "t", "col_members", "pon_members",    # DTU10MWFloating
                 "Diam_y", "Diam_z", "float_members", "joint_members"]            # SAITEC2MWFloating


def read_turbine(path   : Path = None,      # [-] Path to the single-turbine .hdf5 file
                 verbose: bool = False,     # [-] if True, prints a list of keys available in the file
                 steps  : list = None):     # [-] Post-processing steps to read (keys of STEPS). None: all of them
    """
    Reads a single-turbine .hdf5 file coming from Fortran.

    Returns
    -------
    raw : dict
        Datasets with their Fortran names: turbine parameters (read_turbine_params), solver parameters
        (read_solver_params), Freqs and the post-processing steps, each with a has_<step> flag since any of
        them may be missing (files from older versions of the Fortran code or from interrupted runs).
    """

    path = Path(path)
    print(f"Loading HDF5 file: {path}")

    if not path.exists():
        raise FileNotFoundError(f"readers.turbine.read_turbine(): File not found {path}")

    with h5py.File(path, 'r') as f:
        if verbose:
            print_keys(f)

        emitter = detect_emitter(f)
        if emitter != "turbine":
            raise wrong_emitter(path, emitter, "readers.turbine.read_turbine")

        raw = {}
        raw["Path"] = path
        raw["Name"] = path.stem

        # Turbine parameters: generic ones and those of its type (D, wet_area, ...)
        raw.update(read_turbine_params(f))

        # Files saved before "case_type" was persisted: guess it from the file name
        if raw["case_type"] == "Unknown":
            raw["case_type"] = guess_case_type(path)

        # Solver parameters: method autodetected, water depth from the boundaries
        read_solver_params(f, raw)

        # Common variables
        raw["Freqs"] = read_array(f, "Freqs")

        # Post-processing steps
        read_steps(f, raw, steps)

    return raw

def read_turbine_params(f,                  # [-] Open h5py.File
                        suffix: str = ""):  # [-] "" for a single turbine, "_i" for turbine i of a farm
    """
    Reads the generic parameters of a turbine (WindTurbine.save_turbine_params) and those of its type
    (save_specific_params: D, wet_area, col_D, ...).

    Returns
    -------
    turb : dict
        case_type is "Unknown" when the file has neither case_type nor Case_type.
    """

    turb = {}

    # Literal case type, e.g. "DTU10MWMonopile" (saved both as case_type and Case_type)
    case_type = read_ds(f, "case_type" + suffix)
    if case_type is None:
        case_type = read_ds(f, "Case_type" + suffix)

    turb["case_type"] = case_type if case_type is not None else "Unknown"

    # Common variables
    turb["WindSpeed"] = read_float(f, "WindSpeed" + suffix)
    turb["WindDir"]   = read_float(f, "WindDir" + suffix)
    turb["Depth"]     = read_float(f, "Depth" + suffix)
    turb["AxisPos"]   = read_array(f, "AxisPos" + suffix)
    turb["BariPos"]   = read_array(f, "BariPos" + suffix)

    # Structure nodes shape(Nm, Nn, 3). SAITEC's rigid members are NaN-padded: plots skip those nodes
    nodes = read_ds(f, "Structure_nodes" + suffix)
    if nodes is not None:
        turb["Structure_nodes"] = nodes
        turb["Nm"], turb["Nn"]  = (nodes.shape[0], nodes.shape[1]) if nodes.ndim == 3 else (None, None)

    # Type-specific parameters
    for key in SPECIFIC_KEYS:
        val = read_ds(f, key + suffix)
        if val is not None:
            turb[key] = val

    return turb

def guess_case_type(path: Path = None):     # [-] Path to a single-turbine .hdf5 file
    """
    Fallback for files saved before "case_type" was persisted: plot_mn_* files are monopiles, the rest
    floating platforms.

    Returns
    -------
    case_type : str
    """

    if "mn" in path.name or "monopile" in str(path).lower():
        return "Monopile"

    return "Floating"
