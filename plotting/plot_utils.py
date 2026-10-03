"""
Module: plot_utils.py
Description: Command line interface, discovery and loading of the .hdf5 files written by the Fortran code
             (single turbine and wind farm), figure saving and small helpers shared by every plot.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 09/2026
"""

import re
import sys
import h5py
import argparse
import numpy as np
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import utils.AcousticUtils as au
import plot_style as ps
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# ---------- DEFAULT PATHS ---------- #
FORTRAN_DIR      = Path(__file__).resolve().parent.parent     # [-] FORTRAN/ (this file lives in FORTRAN/plotting/)
TURBINE_DATA_DIR = FORTRAN_DIR / "turbine_acoustic_data"      # [-] Default folder of the single-turbine files
FARM_DATA_DIR    = FORTRAN_DIR / "farm_acoustic_data"         # [-] Default folder of the wind farm files
FIGURES_DIR      = FORTRAN_DIR / "figures"                    # [-] Default folder of the saved figures

# Post-processing steps (AcousticSolver.run_*) and the dataset holding the pressure of each one
STEPS = {
    "spectrums"    : "P_spectrums",
    "polar"        : "P_polar",
    "cylinder"     : "P_cylinder",
    "decay"        : "P_decay",
    "line"         : "P_line",
    "sliceXY"      : "P_slicexy",
    "sliceXZ"      : "P_slicexz",
    "sliceVertical": "P_sliceVertical",
    "sphere"       : "P_sphere",
}

# Type-specific turbine parameters (save_specific_params of FORTRAN/src/emitters/TurbineTypes/*.f90)
SPECIFIC_KEYS = ["D", "rho_mat", "wet_area",                                      # DTU10MWMonopile (+ the rest)
                 "col_D", "XsecA", "XsecB", "t", "col_members", "pon_members",    # DTU10MWFloating
                 "Diam_y", "Diam_z", "float_members", "joint_members"]            # SAITEC2MWFloating

# Name written by the Fortran get_name() of every method ("Method" dataset), for the files that lack it
METHOD_NAMES = {"ImageMethod"          : "Images Method",
                "AnalyticalNormalModes": "Analytical Normal Modes"}


# ---------- TURBINE UTILITIES ---------- #
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
    add_common_args(parser, TURBINE_DATA_DIR)

    return finish_args(parser.parse_args())

def load_turbine(path   : Path = None,      # [-] Path to the single-turbine .hdf5 file
                 verbose: bool = False):    # [-] if True, prints a list of keys available in the file
    """
    Loads all available data from a single-turbine .hdf5 file coming from Fortran.

    Returns
    -------
    case : dict
        Simulation data. Every post-processing step has a has_<step> flag, since any of them may be missing
        (files from older versions of the Fortran code or from interrupted runs).
    """

    path = Path(path)
    print(f"Loading HDF5 file: {path}")

    if not path.exists():
        raise FileNotFoundError(f"plot_utils.load_turbine(): File not found {path}")

    with h5py.File(path, 'r') as f:
        if verbose:
            print_keys(f)

        if "Nturbines" in f:
            raise ValueError(f"plot_utils.load_turbine(): {path.name} is a wind farm file, "
                             f"plot it with plotting/plot_farm.py")

        case = {}
        case["Path"] = path
        case["Name"] = path.stem

        # Turbine parameters: generic ones and those of its type (D, wet_area, ...)
        case.update(read_turbine_params(f))

        # Files saved before "case_type" was persisted: guess it from the file name
        if case["case_type"] == "Unknown":
            case["case_type"] = case["Case_type"] = guess_case_type(path)

        # Solver parameters: method autodetected, water depth from the boundaries
        read_solver_params(f, case)

        # Common variables
        case["Freqs"] = read_array(f, "Freqs")

        # Post-processing steps
        read_steps(f, case)

    return case

def load_turbines(args   : argparse.Namespace = None,   # [-] Parsed command-line arguments
                  verbose: bool = False):               # [-] if True, prints a list of keys available in the file
    """
    Loads the case given in the terminal, or all of them to compare them.

    Returns
    -------
    cases : list of dict
        One dictionary per file. With several files every case gets a distinct label (Case_type) and color.
    """

    paths = resolve_paths(args.files, args.data_dir)
    cases = [load_turbine(path, verbose=verbose or args.verbose) for path in paths]

    if len(cases) > 1:
        set_comparison_labels(cases)

    return cases

def set_comparison_labels(cases: list = None):     # [-] List of dictionaries containing simulation data
    """
    Gives every case of a comparison a distinct label, used by titles (Case_type), built from its case type,
    wind speed and method (e.g. 'DTU10MWFloating, 11.4 m/s, ANM231'). The water depth is added if two labels
    match, and the file name if they still do. Colors are assigned by position.
    """

    labels = [f"{tex_escape(case['case_type'])}, {case['WindSpeed']:.1f} m/s, {case['Method_tag']}" for case in cases]

    if len(set(labels)) < len(labels):
        labels = [f"{label}, H = {case['Depth']:.0f} m" for label, case in zip(labels, cases)]
    if len(set(labels)) < len(labels):
        labels = [f"{label} ({tex_escape(case['Name'])})" for label, case in zip(labels, cases)]

    for case, label, color in zip(cases, labels, get_case_colors(len(cases))):
        case["Case_type"] = label
        case["Color"]     = color


# ---------- FARM UTILITIES ---------- #
def parse_farm_args():
    """
    Parses arguments from the terminal.
    """

    parser = argparse.ArgumentParser(description="Plots a wind farm (one .hdf5 file).")
    parser.add_argument(
        "file", type = str, nargs = "?", default = None,
        help = "Wind farm .hdf5 file written by Fortran. A bare file name is also searched in --data_dir "
               "(default: most recent .hdf5 in --data_dir)"
    )
    add_common_args(parser, FARM_DATA_DIR)

    args = parser.parse_args()
    args.files = [args.file] if args.file else []

    return finish_args(args)

def load_farm(args   : argparse.Namespace = None,   # [-] Parsed command-line arguments
              verbose: bool = False):               # [-] if True, prints a list of keys available in the file
    """
    Loads files from farm simulation

    Returns
    -------
    case : dict
        Simulation data. The parameters of each turbine are in case["Turbine_parameters"] (one dict per
        turbine, same keys as a single-turbine case) and every post-processing step has a has_<step> flag.
    """

    path = resolve_paths(args.files, args.data_dir)[0]
    print(f"Loading HDF5 file: {path}")

    with h5py.File(path, 'r') as f:
        if verbose or args.verbose:
            print_keys(f)

        if "Nturbines" not in f:
            raise ValueError(f"plot_utils.load_farm(): {path.name} is a single-turbine file, "
                             f"plot it with plotting/plot_turbine.py")

        # ---------- Create case dictionary ---------- #
        case = {}
        case["Path"] = path
        case["Name"] = path.stem

        # Farm parameters: turbine i is saved with the suffix "_i" (WindFarm.save_params)
        case["Nturbines"]          = int(read_ds(f, "Nturbines"))
        case["Source_offset"]      = read_array(f, "Source_offset")
        case["Turbine_parameters"] = [read_turbine_params(f, f"_{i}") for i in range(1, case["Nturbines"] + 1)]

        # Solver parameters: method autodetected, water depth from the boundaries (no "Depth" in farm files)
        read_solver_params(f, case)

        # Common variables
        case["Freqs"] = read_array(f, "Freqs")

        # Post-processing steps
        read_steps(f, case)

    return case


# ---------- FILE UTILITIES ---------- #
def add_common_args(parser  : argparse.ArgumentParser = None,   # [-] Parser of the turbine or farm script
                    data_dir: Path = None):                     # [-] Default folder of the .hdf5 files
    """
    Adds the arguments shared by the turbine and the farm scripts.
    """

    parser.add_argument(
        "--data_dir", type = str, default = str(data_dir),
        help = f"Directory where the .hdf5 files are stored (default: {rel(data_dir)})"
    )
    parser.add_argument(
        "--list", action = "store_true",
        help = "Lists the .hdf5 files in --data_dir (date, case, method and saved steps) and exits."
    )
    parser.add_argument(
        "--save", action = "store_true",
        help = "Saves generated plots in --output_dir, one subfolder per case."
    )
    parser.add_argument(
        "--output_dir", type = str, default = str(FIGURES_DIR),
        help = f"Path where plots will be saved (default: {rel(FIGURES_DIR)})"
    )
    parser.add_argument(
        "--format", type = str, nargs = "+", default = ["png"], choices = ["png", "pdf", "svg", "eps"],
        help = "File format(s) of the saved plots (default: png)"
    )
    parser.add_argument(
        "--no-show", action = "store_true",
        help = "Does not open the plot windows (e.g. together with --save on a cluster)."
    )
    parser.add_argument(
        "--verbose", action = "store_true",
        help = "Prints the datasets available in every loaded file."
    )

def finish_args(args: argparse.Namespace = None):     # [-] Parsed command-line arguments
    """
    Converts the path arguments into Path objects and serves --list (prints the available files and exits).
    """

    args.data_dir   = Path(args.data_dir)
    args.output_dir = Path(args.output_dir)

    if args.list:
        list_hdf5(args.data_dir)
        sys.exit(0)

    return args

def find_hdf5(data_dir: Path = None):     # [-] Folder of the .hdf5 files
    """
    Returns
    -------
    files : list of Path
        The .hdf5 files of data_dir, newest first.
    """

    if not data_dir.is_dir():
        return []

    return sorted(data_dir.glob("*.hdf5"), key=lambda path: path.stat().st_mtime, reverse=True)

def resolve_paths(files   : list = None,    # [-] File names given in the terminal
                  data_dir: Path = None):   # [-] Folder of the .hdf5 files
    """
    Turns the file names given in the terminal into existing paths. Each name is tried as given and inside
    data_dir, with and without the .hdf5 extension. Without names, the most recent .hdf5 of data_dir is used.

    Returns
    -------
    paths : list of Path
    """

    if not files:
        found = find_hdf5(data_dir)
        if not found:
            raise FileNotFoundError(f"plot_utils.resolve_paths(): no .hdf5 files found in {data_dir}")
        print(f"No file given: using the most recent one in {rel(data_dir)}")
        return [found[0]]

    paths = []
    for name in files:
        candidates = [Path(name), data_dir / name]
        if not name.endswith(".hdf5"):
            candidates += [Path(name + ".hdf5"), data_dir / (name + ".hdf5")]

        path = next((candidate for candidate in candidates if candidate.is_file()), None)
        if path is None:
            raise FileNotFoundError(f"plot_utils.resolve_paths(): file not found {name} "
                                    f"(use --list to see the files in {rel(data_dir)})")
        paths.append(path)

    return paths

def list_hdf5(data_dir: Path = None):     # [-] Folder of the .hdf5 files
    """
    Prints the .hdf5 files of data_dir, newest first, with the case they contain and the steps saved in them.
    """

    files = find_hdf5(data_dir)
    if not files:
        print(f"No .hdf5 files found in {rel(data_dir)}")
        return

    print(f"\n.hdf5 files in {rel(data_dir)} (newest first):\n")
    for idx, path in enumerate(files):
        try:
            with h5py.File(path, 'r') as f:
                description, steps = describe_hdf5(f)
        except Exception as e:
            description, steps = f"unreadable file ({type(e).__name__})", "-"

        stat   = path.stat()
        date   = datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M")
        marker = "*" if idx == 0 else " "

        print(f"  {marker} {date}  {stat.st_size / 1e6:6.0f} MB  {path.name}")
        print(f"        {description}")
        print(f"        steps: {steps}\n")

    print("  * = file plotted when no file is given")

def describe_hdf5(f):     # [-] Open h5py.File
    """
    Returns
    -------
    description, steps : str, str
        Case, method, wind speed and depth of the file, and the post-processing steps it contains.
    """

    if "Nmodes" in f:
        method = f"ANM{read_ds(f, 'Nmodes')}"
    elif "N_images" in f:
        method = f"IM{read_ds(f, 'N_images')}"
    else:
        method = read_ds(f, "Method") or "unknown method"

    if "Nturbines" in f:
        n_turb = int(read_ds(f, "Nturbines"))
        types  = [read_ds(f, f"case_type_{i}") or read_ds(f, f"Case_type_{i}") for i in range(1, n_turb + 1)]
        speeds = [read_ds(f, f"WindSpeed_{i}") for i in range(1, n_turb + 1)]
        name   = "Wind farm: " + ", ".join(f"{types.count(t)} x {t}" for t in dict.fromkeys(types))
        ws     = f"{min(speeds):.1f}-{max(speeds):.1f}" if min(speeds) != max(speeds) else f"{speeds[0]:.1f}"
    else:
        name = read_ds(f, "case_type") or read_ds(f, "Case_type") or f"{guess_case_type(Path(f.filename))} (from file name)"
        ws   = f"{read_ds(f, 'WindSpeed'):.1f}" if "WindSpeed" in f else "?"

    depth = ""
    if "Upper_HBC" in f and "Lower_HBC" in f:
        depth = f" | H = {read_ds(f, 'Upper_HBC') - read_ds(f, 'Lower_HBC'):.0f} m"

    steps = ", ".join(step for step, key in STEPS.items() if key in f)

    return f"{name} | {method} | {ws} m/s{depth}", (steps if steps else "none")

def print_keys(f):     # [-] Open h5py.File
    """
    Prints an elegant list of the datasets available in the file.
    """

    keys = sorted(f.keys())
    print("--------------------")
    for idx, key in enumerate(keys, start=1):
        print(f"{idx:2d}. {key}")
    print("--------------------")

def read_ds(f,          # [-] Open h5py.File
            key: str):  # [-] Dataset name
    """
    Auxiliar function to read datasets written by the Fortran code and process them safely.

    - Fortran is column-major, so h5py sees its dimensions reversed: every array with ndim >= 2 is
      transposed back to the Fortran order, e.g. P_cylinder -> (Nf, Ntheta, Nz).
    - Complex arrays (every P_* dataset) are saved as reals with shape (2, N1, N2, ...), i.e. (..., N2, N1, 2)
      in h5py: the real/imaginary pair is ALWAYS the last h5py axis, whatever the size of the other axes.
    - Scalars are shape-(1,) datasets and fixed-length strings may carry padding spaces.

    Returns
    -------
    None if the dataset does not exist; str for strings; int/float for scalars; np.ndarray otherwise.
    """

    if key not in f:
        return None

    val = f[key][()]

    # Decode text strings (spare bytes for fortran)
    if isinstance(val, bytes):
        return val.decode("utf-8").strip()

    # Decode string arrays: unwrap to a plain string when there is only one
    if isinstance(val, np.ndarray) and val.dtype.kind in ('S', 'U'):
        decoded = [x.decode('utf-8').strip() if isinstance(x, bytes) else x.strip() for x in val.ravel()]
        return decoded[0] if len(decoded) == 1 else decoded

    # Complex pressure: combine the (real, imaginary) pair of the last axis before transposing
    if key.startswith("P_") and isinstance(val, np.ndarray) and val.ndim >= 2:
        if val.shape[-1] != 2:
            raise ValueError(f"plot_utils.read_ds(): complex dataset {key} with shape {val.shape} "
                             f"has no (real, imaginary) pair on its last axis")
        val = np.ascontiguousarray(val, dtype=np.float64).view(np.complex128)[..., 0]

    # Transpose multidimensional arrays for Column-Major
    if isinstance(val, np.ndarray) and val.ndim >= 2:
        val = val.T

    # Collect scalars: NumPy 2.x no longer allows int()/float() on a shape-(1,) array
    if isinstance(val, np.ndarray) and val.shape == (1,):
        return val.item()

    return val

def read_array(f,          # [-] Open h5py.File
               key: str):  # [-] Dataset name
    """
    Reads a dataset that must stay an array even with a single element (read_ds returns a scalar then).

    Returns
    -------
    None if the dataset does not exist, np.ndarray (ndim >= 1) otherwise.
    """

    val = read_ds(f, key)

    return None if val is None else np.atleast_1d(np.asarray(val))

def read_float(f,          # [-] Open h5py.File
               key: str):  # [-] Dataset name
    """
    Reads a scalar dataset that older files may lack.

    Returns
    -------
    val : float
        The value, or NaN if the dataset does not exist.
    """

    val = read_ds(f, key)

    return np.nan if val is None else float(val)

def read_solver_params(f,                   # [-] Open h5py.File
                       case: dict = None):  # [-] Dictionary where the parameters are stored
    """
    Reads the solver and environment parameters, valid for turbines and farms:

    - Method_name: the "Method" dataset (Fortran get_name(), e.g. "Images Method"), written in every file
      since 10/2026. Older files lack it: their method is told by its parameters (Nmodes: Analytical Normal
      Modes, N_images: Image Method). Method ("ImageMethod", "AnalyticalNormalModes") and Method_tag
      ("IM30", "ANM24") keep the names used so far; any other method gets its name without blanks and its
      initials (e.g. "Parabolic Equation" -> "ParabolicEquation", "PE").
    - ENBW: equivalent noise bandwidth of the FFT window of the forces, in bins (see read_enbw).
    - Water depth: Upper_HBC - Lower_HBC.
    """

    method_name = read_ds(f, "Method")

    if "Nmodes" in f:
        case["Method"]     = "AnalyticalNormalModes"
        case["Nmodes"]     = int(read_ds(f, "Nmodes"))
        case["H"]          = read_ds(f, "H")
        case["Method_tag"] = f"ANM{case['Nmodes']}"
    elif "N_images" in f:
        case["Method"]     = "ImageMethod"
        case["N_images"]   = int(read_ds(f, "N_images"))
        case["Method_tag"] = f"IM{case['N_images']}"
    elif method_name is not None:
        case["Method"]     = method_name.replace(" ", "")
        case["Method_tag"] = "".join(word[0] for word in method_name.split()).upper()
    else:
        raise ValueError(f"plot_utils.read_solver_params(): unknown acoustic method in {case['Path']} "
                         f"(none of Method, Nmodes or N_images found)")

    case["Method_name"] = method_name if method_name is not None else METHOD_NAMES[case["Method"]]
    case["ENBW"]        = read_enbw(f)

    case["c_wat"]     = float(read_ds(f, "c_wat"))
    case["rho_wat"]   = float(read_ds(f, "rho_wat"))
    case["p_ref"]     = float(read_ds(f, "p_ref"))
    case["Upper_HBC"] = float(read_ds(f, "Upper_HBC"))
    case["Lower_HBC"] = float(read_ds(f, "Lower_HBC"))
    case["Depth"]     = case["Upper_HBC"] - case["Lower_HBC"]     # [m] Water depth

def read_enbw(f):     # [-] Open h5py.File
    """
    Equivalent noise bandwidth, in bins, of the FFT window applied to the forces. Every bin holds the
    amplitude of a tone, so the sums over bins (OASPL, PSD, band powers) are divided by it. The Fortran
    save_shared_parameters writes "ENBW" only when the forces come from a windowed FFT:

    - "ENBW" in the file: that value (1.5, Hann window of compute_rfft);
    - neither "ENBW" nor "Method": file written before 10/2026, always a turbine or a farm whose forces come
      from compute_rfft, so au.HANN_ENBW;
    - "Method" without "ENBW": forces without a window (e.g. manual sources given in frequency), so 1.0, i.e.
      no correction.

    Returns
    -------
    enbw : float
    """

    enbw = read_ds(f, "ENBW")
    if enbw is not None:
        return float(enbw)

    return au.HANN_ENBW if "Method" not in f else 1.0

def read_turbine_params(f,                  # [-] Open h5py.File
                        suffix: str = ""):  # [-] "" for a single turbine, "_i" for turbine i of a farm
    """
    Reads the generic parameters of a turbine (WindTurbine.save_turbine_params) and those of its type
    (save_specific_params: D, wet_area, col_D, ...).

    Returns
    -------
    turb : dict
    """

    turb = {}

    # Literal case type, e.g. "DTU10MWMonopile" (saved both as case_type and Case_type)
    case_type = read_ds(f, "case_type" + suffix)
    if case_type is None:
        case_type = read_ds(f, "Case_type" + suffix)

    turb["case_type"] = case_type if case_type is not None else "Unknown"
    turb["Case_type"] = tex_escape(turb["case_type"])     # Label used in titles, replaced in comparisons

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

def read_steps(f,                   # [-] Open h5py.File
               case: dict = None):  # [-] Dictionary where the data are stored
    """
    Reads every post-processing step present in the file (AcousticSolver.run_*). Any of them may be missing,
    so each one gets a has_<step> flag. Shapes follow the Fortran order: pressures are (Nf, ...), observers
    (..., 3).
    """

    # 1. Spectrum variables
    case["has_spectrums"] = "P_spectrums" in f
    if case["has_spectrums"]:
        case["P_spectrums"]   = read_ds(f, "P_spectrums")         # (Nf, Nobs)
        case["Obs_spectrums"] = read_ds(f, "Obs_spectrums")       # (Nobs, 3)

    # 2. Polar
    case["has_polar"] = "P_polar" in f
    if case["has_polar"]:
        case["P_polar"]         = read_ds(f, "P_polar")           # (Nf, Ntheta)
        case["R_polar"]         = float(read_ds(f, "R_polar"))
        case["Z_polar"]         = float(read_ds(f, "Z_polar"))
        case["Obs_polar"]       = read_ds(f, "Obs_polar")         # (Ntheta, 3)
        case["Theta_deg_polar"] = read_array(f, "Theta_deg_polar")
        case["Center_polar"]    = read_array(f, "Center_polar")

    # 3. Cylinder
    case["has_cylinder"] = "P_cylinder" in f
    if case["has_cylinder"]:
        case["P_cylinder"]         = read_ds(f, "P_cylinder")     # (Nf, Ntheta, Nz)
        case["R_cylinder"]         = float(read_ds(f, "R_cylinder"))
        case["Z_cylinder"]         = read_array(f, "Z_cylinder")
        case["Obs_cylinder"]       = read_ds(f, "Obs_cylinder")   # (Ntheta, Nz, 3)
        case["Theta_deg_cylinder"] = read_array(f, "Theta_deg_cylinder")
        case["Center_cylinder"]    = read_array(f, "Center_cylinder")
        if "dA_cylinder" in f:                                    # Older files: plot_cylinder computes it
            case["dA_cylinder"]    = float(read_ds(f, "dA_cylinder"))

    # 4. Decay (single turbine only)
    case["has_decay"] = "P_decay" in f
    if case["has_decay"]:
        case["P_decay"]         = read_ds(f, "P_decay")           # (Nf, N)
        case["Z_decay"]         = float(read_ds(f, "Z_decay"))
        case["Distance_decay"]  = read_array(f, "Distance_decay")
        case["Distances_decay"] = read_array(f, "Distances_decay")
        case["Obs_decay"]       = read_ds(f, "Obs_decay")         # (N, 3)
        case["Logspace_decay"]  = bool(read_ds(f, "Logspace_decay"))

    # 5. Line
    case["has_line"] = "P_line" in f
    if case["has_line"]:
        case["P_line"]         = read_ds(f, "P_line")             # (Nf, N)
        case["P1_line"]        = read_array(f, "P1_line")
        case["P2_line"]        = read_array(f, "P2_line")
        case["Distance_line"]  = float(read_ds(f, "Distance_line"))
        case["Distances_line"] = read_array(f, "Distances_line")
        case["Obs_line"]       = read_ds(f, "Obs_line")           # (N, 3)
        case["Logspace_line"]  = bool(read_ds(f, "Logspace_line"))

        # The plots use Distances_line: warn if the observers were not placed at those distances from P1
        obs_distances = np.linalg.norm(case["Obs_line"] - case["P1_line"], axis=1)
        mismatch      = np.max(np.abs(obs_distances - case["Distances_line"]))
        if mismatch > 1e-3 * max(case["Distance_line"], 1.0):
            print(f"WARNING plot_utils.read_steps(): the line observers (Obs_line) lie up to {obs_distances.max():.3g} m "
                  f"from P1, but Distances_line goes up to {case['Distances_line'].max():.3g} m. The line plot uses "
                  f"Distances_line (check AcousticSolver.run_line)")

    # 6. SliceXY
    case["has_slicexy"] = "P_slicexy" in f
    if case["has_slicexy"]:
        case["P_slicexy"]      = read_ds(f, "P_slicexy")          # (Nf, nx, ny)
        case["Obs_slicexy"]    = read_ds(f, "Obs_slicexy")        # (nx, ny, 3)
        case["X_slicexy"]      = read_array(f, "X_slicexy")
        case["Y_slicexy"]      = read_array(f, "Y_slicexy")
        case["Z_slicexy"]      = float(read_ds(f, "Z_slicexy"))
        case["Center_slicexy"] = read_array(f, "Center_slicexy")
        case["Nx_slicexy"]     = int(read_ds(f, "Nx_slicexy"))
        case["Ny_slicexy"]     = int(read_ds(f, "Ny_slicexy"))

    # 7. SliceXZ
    case["has_slicexz"] = "P_slicexz" in f
    if case["has_slicexz"]:
        case["P_slicexz"]   = read_ds(f, "P_slicexz")             # (Nf, nx, nz)
        case["Obs_slicexz"] = read_ds(f, "Obs_slicexz")           # (nx, nz, 3)
        case["X_slicexz"]   = read_array(f, "X_slicexz")
        case["Z_slicexz"]   = read_array(f, "Z_slicexz")
        case["Y_slicexz"]   = float(read_ds(f, "Y_slicexz"))
        case["Nx_slicexz"]  = int(read_ds(f, "Nx_slicexz"))
        case["Nz_slicexz"]  = int(read_ds(f, "Nz_slicexz"))

    # 8. Vertical slice: normal (cos az, sin az), horizontal axis (-sin az, cos az), observer = centre + u * axis
    case["has_slicevertical"] = "P_sliceVertical" in f
    if case["has_slicevertical"]:
        case["P_sliceV"]       = read_ds(f, "P_sliceVertical")    # (Nf, nu, nz)
        case["Coords_sliceV"]  = read_ds(f, "Coords_sliceV")      # (nu, nz, 3)
        case["U_sliceV"]       = read_array(f, "U_sliceV")
        case["Z_sliceV"]       = read_array(f, "Z_sliceV")
        case["Azimuth_sliceV"] = float(read_ds(f, "Azimuth_sliceV"))
        case["Center_sliceV"]  = read_array(f, "Center_sliceV")

    # 9. Sphere (single turbine and Image Method only)
    case["has_sphere"] = "P_sphere" in f
    if case["has_sphere"]:
        case["P_sphere"]      = read_ds(f, "P_sphere")            # (Nf, Ntheta, Nz)
        case["Obs_sphere"]    = read_ds(f, "Obs_sphere")          # (Ntheta, Nz, 3)
        case["Theta_sphere"]  = read_array(f, "Theta_sphere")
        case["Z_sphere"]      = read_array(f, "Z_sphere")
        case["R_sphere"]      = float(read_ds(f, "R_sphere"))
        case["Center_sphere"] = read_array(f, "Center_sphere")
        case["dA_sphere"]     = float(read_ds(f, "dA_sphere"))

        # Extra dimensions reconstruction
        case["N_theta"]   = len(case["Theta_sphere"])
        case["Nz_sphere"] = len(case["Z_sphere"])

def rel(path: Path = None):     # [-] Any path
    """
    Returns
    -------
    path : str
        The path relative to the working directory when it is below it, absolute otherwise.
    """

    try:
        return str(Path(path).resolve().relative_to(Path.cwd()))
    except ValueError:
        return str(path)


# ---------- FIGURE UTILITIES ---------- #
def set_figure_name(fig,        # [-] Figure created by a plot function
                    *parts):    # [-] Words of the name, e.g. "spectrum", "PSD", "octave"
    """
    Names a figure after the plot that created it (e.g. 'spectrum_PSD_octave', 'cylinder_SPL_0.88Hz').
    save_figures uses the name as file name.
    """

    words = []
    for part in parts:
        word = str(part).strip() if part is not None else ""
        word = re.sub(r"^at\s+", "", word)          # 'at 0.88 Hz' -> '0.88 Hz'
        word = re.sub(r"\s+", "", word)             # '0.88 Hz'    -> '0.88Hz'
        word = re.sub(r"[^\w.\-]", "", word)
        if word:
            words.append(word)

    fig.set_label("_".join(words))

def save_figures(args : argparse.Namespace = None,     # [-] Parsed command-line arguments
                 cases: list = None):                  # [-] Plotted case (dict) or cases (list of dict)
    """
    With --save, saves every open figure as <output_dir>/<case name>/<figure name>.<format>. A comparison of
    several turbines goes to <output_dir>/comparison_<name 1>_vs_<name 2>...
    """

    if not args.save:
        return

    if isinstance(cases, dict):
        cases = [cases]

    names  = [case["Name"] for case in cases]
    folder = names[0] if len(names) == 1 else "comparison_" + "_vs_".join(names)
    if len(folder) > 200:
        folder = f"comparison_{len(names)}_cases_{names[0]}"

    out_dir = args.output_dir / folder
    out_dir.mkdir(parents=True, exist_ok=True)

    used = {}
    for num in plt.get_fignums():
        fig  = plt.figure(num)
        name = fig.get_label() or f"figure_{num}"

        # Repeated plots (same function and mode) get a counter
        used[name] = used.get(name, 0) + 1
        if used[name] > 1:
            name = f"{name}_{used[name]}"

        for fmt in args.format:
            fig.savefig(out_dir / f"{name}.{fmt}", format=fmt)

    print(f"\n{len(plt.get_fignums())} figures saved in {rel(out_dir)} ({', '.join(args.format)})")

def tex_escape(text: str = ""):     # [-] Text to show in a figure
    """
    Escapes the characters LaTeX treats as special (e.g. the '_' of file names) when text.usetex is on.

    Returns
    -------
    text : str
    """

    if not plt.rcParams["text.usetex"]:
        return text

    return re.sub(r"([_&%#$])", r"\\\1", text)


# ---------- COMMON UTILITIES ---------- #
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

def get_case_color(case: dict = None,       # [-] Dictionary containing simulation data
                   tags: list = None):      # [-] List with tags
    """
    Gets color based on case type. In a comparison every case carries its own color (by position).
    """

    def color_from_tag(tag: str):
        if "Monopile" in tag:
            return ps.AcademicStyle.CASE_COLORS[0]
        if "shallow" in tag:
            return ps.AcademicStyle.CASE_COLORS[1]
        return ps.AcademicStyle.CASE_COLORS[2]

    if case is not None:
        if "Color" in case:
            return case["Color"]
        return color_from_tag(case["Case_type"])

    if tags is not None:
        return np.array([color_from_tag(tag) for tag in tags])

    raise ValueError("get_case_color(): either 'case' or 'tags' must be provided")

def get_case_colors(n: int = 1):    # [-] Number of cases
    """
    Colors of the cases of a comparison, by position: CASE_COLORS first, then BAND_COLORS.

    Returns
    -------
    colors : list of str
    """

    palette = ps.AcademicStyle.CASE_COLORS + [c for c in ps.AcademicStyle.BAND_COLORS if c not in ps.AcademicStyle.CASE_COLORS]

    return [palette[i % len(palette)] for i in range(n)]

def get_color_shades(color,             # [-] Base color (any matplotlib color)
                     n: int = 1):       # [-] Number of shades
    """
    Shades of a single color, from darker to lighter, to tell apart several curves of the same case (e.g.
    the observers of a spectrum). The base color is mixed with black (darker shades) or with white (lighter
    tints), so every curve still reads as the same color. A single shade is the base color itself.

    Returns
    -------
    shades : list of tuple
        n RGB colors.
    """

    rgb = np.array(mcolors.to_rgb(color))
    if n <= 1:
        return [tuple(rgb)]

    shades = []
    for t in np.linspace(-0.35, 0.45, n):     # [-] < 0: fraction of black, > 0: fraction of white
        target = np.zeros(3) if t < 0 else np.ones(3)
        shades.append(tuple(rgb + abs(t) * (target - rgb)))

    return shades

def get_line_styles():
    return ps.AcademicStyle.LINE_STYLES

def get_band_colors():
    return ps.AcademicStyle.BAND_COLORS

def get_bin_width(freqs):               # [Hz] Frequency array shape(NFreqs,)
    """
    Gets the frequency bin width df [Hz] needed to convert pressure into PSD.

    Freqs lies on the uniform FFT grid but may have gaps (compute_force(filter_freqs=True)
    drops non-excited bins), so df is the smallest spacing, not the local one. If the
    spectrum was decimated (skipf > 1) the true df cannot be recovered from Freqs.
    Returns None if it cannot be determined (less than two frequencies).
    """

    freqs = np.unique(np.asarray(freqs, dtype=float))
    if freqs.size < 2:
        return None

    return np.min(np.diff(freqs))

def convert_to(p,                               # [Pa] Pressure array shape(NFreqs, Nobs)
               mode: str   = "SPL",             # [-] What to convert to
               df  : float = None,              # [Hz] Frequency bin width, only needed for 'PSD'
               enbw: float = au.HANN_ENBW):     # [bins] ENBW of the FFT window (case["ENBW"]), only used by 'PSD'
    """
    Converts pressure into desired units.
    """

    mode = mode.upper()

    if mode == "SPL":
        vals = au.pressure_to_SPL(p)
        unit_label = r"SPL [dB re 1 $\mu$Pa]"
    elif mode == "PSD":
        if df is None: raise ValueError("convert_to(): mode 'PSD' needs the frequency bin width 'df'")
        vals = au.pressure_to_PSD(p, df, enbw=enbw)
        unit_label = r"PSD [dB re 1 $\mu$Pa$^2$/Hz]"
    elif mode == "ABS":
        vals = np.abs(p)
        unit_label = r"Magnitude [Pa]"
    elif mode == "REAL":
        vals = np.real(p)
        unit_label = r"Real Part [Pa]"
    elif mode == "IMAG":
        vals = np.imag(p)
        unit_label = r"Imaginary Part [Pa]"
    elif mode == "PHASE":
        vals = np.angle(p, deg=True)
        if np.any(vals < -90) and np.any(vals > 90):
            vals = np.mod(vals, 360)
        unit_label = r"Phase [deg]"
    else:
        print(f"convert_to(): Unrecognized mode '{mode}', defaulting to 'SPL'")
        vals = au.pressure_to_SPL(p)
        unit_label = r"SPL [dB re 1 $\mu$Pa]"

    return vals, unit_label

def get_turbine_styles():
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

def is_valid(array,                 # [Pa] Pressure array shape(NFreqs, ...)
             name: str = None):     # [-] Name of the calling function, for the messages
    """
    Returns True if the pressure array can be plotted: it exists, has no +-inf entries and at least one
    finite value.

    NaN entries are valid on purpose: the Fortran code writes NaN (real and imaginary parts, at every
    frequency) at the observers lying closer than 5 m to a source node, i.e. inside the structure. They
    carry no data and every plot shows them as holes. When there are some, a message tells how many.
    """

    if array is None:
        return False

    a = np.asarray(array)
    if a.size == 0:
        if name is not None: print(f"{name}(): empty pressure array, skipping")
        return False

    if np.any(np.isinf(a)):
        if name is not None: print(f"{name}(): inf found in pressure array, skipping")
        return False

    nan = np.isnan(a)
    if np.all(nan):
        if name is not None: print(f"{name}(): no observer has data (NaN everywhere), skipping")
        return False

    if name is not None and np.any(nan):
        masked = np.all(nan.reshape(a.shape[0], -1), axis=0) if a.ndim > 1 else nan
        print(f"{name}(): {np.count_nonzero(masked)} of {masked.size} observers without data "
              f"(NaN, inside the structure)")

    return True

