"""
Module: readers/solver.py
Description: Datasets the solver writes in every file, whatever the emitter: acoustic method, window
             correction (ENBW), environment and the post-processing steps (AcousticSolver.run_*).

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import utils.AcousticUtils as au
from readers.hdf5_io import read_ds, read_array


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

# --ws: only the steps drawn as lines are read from the files other than the given one (their sphere,
# cylinder and slices are GB)
WS_STEPS = ["spectrums", "polar", "decay", "line"]

# Name written by the Fortran get_name() of every method ("Method" dataset), for the files that lack it
METHOD_NAMES = {"ImageMethod"          : "Images Method",
                "AnalyticalNormalModes": "Analytical Normal Modes"}


def read_solver_params(f,                   # [-] Open h5py.File
                       case: dict = None):  # [-] Dictionary where the parameters are stored
    """
    Reads the solver and environment parameters, valid for every emitter:

    - Method_name: the "Method" dataset (Fortran get_name(), e.g. "Images Method"), written in every file
      since 10/2026. Older files lack it: their method is told by its parameters (Nmodes: Analytical Normal
      Modes, N_images: Image Method). Method ("ImageMethod", "AnalyticalNormalModes") and Method_tag
      ("IM30", "ANM24") keep the names used so far; any other method gets its name without blanks and its
      initials (e.g. "Parabolic Equation" -> "ParabolicEquation", "PE").
    - ENBW: equivalent noise bandwidth of the FFT window of the forces, in bins (see read_enbw).
    - Water depth: - Lower_HBC.
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
        raise ValueError(f"readers.solver.read_solver_params(): unknown acoustic method in {case['Path']} "
                         f"(none of Method, Nmodes or N_images found)")

    case["Method_name"] = method_name if method_name is not None else METHOD_NAMES[case["Method"]]
    case["ENBW"]        = read_enbw(f)

    case["c_wat"]     = float(read_ds(f, "c_wat"))
    case["rho_wat"]   = float(read_ds(f, "rho_wat"))
    case["p_ref"]     = float(read_ds(f, "p_ref"))
    case["Lower_HBC"] = float(read_ds(f, "Lower_HBC"))
    case["Depth"]     = -case["Lower_HBC"]     # [m] Water depth

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

def read_steps(f,                   # [-] Open h5py.File
               case : dict = None,  # [-] Dictionary where the data are stored
               steps: list = None): # [-] Steps to read (keys of STEPS, e.g. "decay"). None: all of them
    """
    Reads every post-processing step present in the file (AcousticSolver.run_*), or only those of 'steps'.
    Any of them may be missing, so each one gets a has_<step> flag (False for a step left unread). Shapes
    follow the Fortran order: pressures are (Nf, ...), observers (..., 3).
    """

    def has(step):      # Step requested and present
        return (steps is None or step in steps) and STEPS[step] in f

    # 1. Spectrum variables
    case["has_spectrums"] = has("spectrums")
    if case["has_spectrums"]:
        case["P_spectrums"]   = read_ds(f, "P_spectrums")         # (Nf, Nobs)
        case["Obs_spectrums"] = read_ds(f, "Obs_spectrums")       # (Nobs, 3)

    # 2. Polar
    case["has_polar"] = has("polar")
    if case["has_polar"]:
        case["P_polar"]         = read_ds(f, "P_polar")           # (Nf, Ntheta)
        case["R_polar"]         = float(read_ds(f, "R_polar"))
        case["Z_polar"]         = float(read_ds(f, "Z_polar"))
        case["Obs_polar"]       = read_ds(f, "Obs_polar")         # (Ntheta, 3)
        case["Theta_deg_polar"] = read_array(f, "Theta_deg_polar")
        case["Center_polar"]    = read_array(f, "Center_polar")

    # 3. Cylinder
    case["has_cylinder"] = has("cylinder")
    if case["has_cylinder"]:
        case["P_cylinder"]         = read_ds(f, "P_cylinder")     # (Nf, Ntheta, Nz)
        case["R_cylinder"]         = float(read_ds(f, "R_cylinder"))
        case["Z_cylinder"]         = read_array(f, "Z_cylinder")
        case["Obs_cylinder"]       = read_ds(f, "Obs_cylinder")   # (Ntheta, Nz, 3)
        case["Theta_deg_cylinder"] = read_array(f, "Theta_deg_cylinder")
        case["Center_cylinder"]    = read_array(f, "Center_cylinder")
        if "dA_cylinder" in f:                                    # Older files: plot_cylinder computes it
            case["dA_cylinder"]    = float(read_ds(f, "dA_cylinder"))

    # 4. Decay (single emitters only)
    case["has_decay"] = has("decay")
    if case["has_decay"]:
        case["P_decay"]         = read_ds(f, "P_decay")           # (Nf, N)
        case["Z_decay"]         = float(read_ds(f, "Z_decay"))
        case["Distance_decay"]  = read_array(f, "Distance_decay")
        case["Distances_decay"] = read_array(f, "Distances_decay")
        case["Obs_decay"]       = read_ds(f, "Obs_decay")         # (N, 3)
        case["Logspace_decay"]  = bool(read_ds(f, "Logspace_decay"))

    # 5. Line
    case["has_line"] = has("line")
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
            print(f"WARNING readers.solver.read_steps(): the line observers (Obs_line) lie up to {obs_distances.max():.3g} m "
                  f"from P1, but Distances_line goes up to {case['Distances_line'].max():.3g} m. The line plot uses "
                  f"Distances_line (check AcousticSolver.run_line)")

    # 6. SliceXY
    case["has_slicexy"] = has("sliceXY")
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
    case["has_slicexz"] = has("sliceXZ")
    if case["has_slicexz"]:
        case["P_slicexz"]   = read_ds(f, "P_slicexz")             # (Nf, nx, nz)
        case["Obs_slicexz"] = read_ds(f, "Obs_slicexz")           # (nx, nz, 3)
        case["X_slicexz"]   = read_array(f, "X_slicexz")
        case["Z_slicexz"]   = read_array(f, "Z_slicexz")
        case["Y_slicexz"]   = float(read_ds(f, "Y_slicexz"))
        case["Nx_slicexz"]  = int(read_ds(f, "Nx_slicexz"))
        case["Nz_slicexz"]  = int(read_ds(f, "Nz_slicexz"))

    # 8. Vertical slice: normal (cos az, sin az), horizontal axis (-sin az, cos az), observer = centre + u * axis
    case["has_slicevertical"] = has("sliceVertical")
    if case["has_slicevertical"]:
        case["P_sliceV"]       = read_ds(f, "P_sliceVertical")    # (Nf, nu, nz)
        case["Coords_sliceV"]  = read_ds(f, "Coords_sliceV")      # (nu, nz, 3)
        case["U_sliceV"]       = read_array(f, "U_sliceV")
        case["Z_sliceV"]       = read_array(f, "Z_sliceV")
        case["Azimuth_sliceV"] = float(read_ds(f, "Azimuth_sliceV"))
        case["Center_sliceV"]  = read_array(f, "Center_sliceV")

    # 9. Sphere (Image Method only)
    case["has_sphere"] = has("sphere")
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
