"""
Module: readers/farm.py
Description: Reader of the wind farm files (FORTRAN/farm_acoustic_data): number of turbines, offsets of their
             sources, the parameters of every turbine (WindFarm.save_params, suffix "_i") and the steps of the
             solver.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import h5py
from pathlib import Path
from readers.hdf5_io import read_ds, read_array, print_keys, detect_emitter, wrong_emitter
from readers.solver import read_solver_params, read_steps
from readers.turbine import read_turbine_params


def read_farm(path   : Path = None,      # [-] Path to the wind farm .hdf5 file
              verbose: bool = False):    # [-] if True, prints a list of keys available in the file
    """
    Reads a wind farm .hdf5 file coming from Fortran.

    Returns
    -------
    raw : dict
        Datasets with their Fortran names. The parameters of each turbine are in raw["Turbine_parameters"]
        (one dict per turbine, as read_turbine_params) and every post-processing step has a has_<step> flag.
    """

    path = Path(path)
    print(f"Loading HDF5 file: {path}")

    if not path.exists():
        raise FileNotFoundError(f"readers.farm.read_farm(): File not found {path}")

    with h5py.File(path, 'r') as f:
        if verbose:
            print_keys(f)

        emitter = detect_emitter(f)
        if emitter != "farm":
            raise wrong_emitter(path, emitter, "readers.farm.read_farm")

        raw = {}
        raw["Path"] = path
        raw["Name"] = path.stem

        # Farm parameters: turbine i is saved with the suffix "_i" (WindFarm.save_params)
        raw["Nturbines"]          = int(read_ds(f, "Nturbines"))
        raw["Source_offset"]      = read_array(f, "Source_offset")
        raw["Turbine_parameters"] = [read_turbine_params(f, f"_{i}") for i in range(1, raw["Nturbines"] + 1)]
        raw["has_nodes"]          = all("Structure_nodes" in turb for turb in raw["Turbine_parameters"])

        # Solver parameters: method autodetected, water depth from the boundaries (no "Depth" in farm files)
        read_solver_params(f, raw)

        # Common variables
        raw["Freqs"] = read_array(f, "Freqs")

        # Post-processing steps
        read_steps(f, raw)

    return raw
