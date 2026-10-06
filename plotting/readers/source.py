"""
Module: readers/source.py
Description: Reader of the frequency-domain source files (FORTRAN/source_acoustic_data, emitter
             FrequencySource, FrequencySource.save_params): positions, multipole order and complex forces of
             the point sources, and the steps of the solver. Its frequencies are the tones given by the user,
             not the bins of an FFT.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import h5py
import numpy as np
from pathlib import Path
from readers.hdf5_io import read_ds, read_array, read_complex, print_keys, detect_emitter, wrong_emitter
from readers.solver import read_solver_params, read_steps


def read_source(path   : Path = None,      # [-] Path to the source .hdf5 file
                verbose: bool = False):    # [-] if True, prints a list of keys available in the file
    """
    Reads a frequency-domain source .hdf5 file coming from Fortran.

    Returns
    -------
    raw : dict
        Datasets with their Fortran names: Emitter, Nsources, NFreqs, Source_xyz (Nsources, 3) [m],
        Source_order (Nsources,) (0 monopole, 1 dipole), Source_Force (NFreqs, Nsources, 3) [N] complex (peak
        amplitude, time factor exp(-i w t)), the solver parameters, Freqs and the post-processing steps.
    """

    path = Path(path)
    print(f"Loading HDF5 file: {path}")

    if not path.exists():
        raise FileNotFoundError(f"readers.source.read_source(): File not found {path}")

    with h5py.File(path, 'r') as f:
        if verbose:
            print_keys(f)

        # Source files written before the Emitter dataset (06/10/2026) saved the positions as Structure_nodes
        if "Emitter" not in f and "Nsources" in f:
            raise ValueError(f"readers.source.read_source(): {path.name} was written before the Emitter dataset and "
                             f"is not supported: regenerate it with the current Fortran code")

        emitter = detect_emitter(f)
        if emitter != "source":
            raise wrong_emitter(path, emitter, "readers.source.read_source")

        raw = {}
        raw["Path"] = path
        raw["Name"] = path.stem

        # Sources: positions (Nsources, 3), multipole order (Nsources,) and forces (NFreqs, Nsources, 3)
        raw["Emitter"]      = read_ds(f, "Emitter")
        raw["Nsources"]     = int(read_ds(f, "Nsources"))
        raw["NFreqs"]       = int(read_ds(f, "NFreqs"))
        raw["Source_xyz"]   = np.atleast_2d(read_ds(f, "Source_xyz"))
        raw["Source_order"] = read_array(f, "Source_order").astype(int)
        raw["Source_Force"] = read_complex(f, "Source_Force").reshape(raw["NFreqs"], raw["Nsources"], 3)

        # Solver parameters: method autodetected, water depth from the boundaries
        read_solver_params(f, raw)

        # Common variables
        raw["Freqs"] = read_array(f, "Freqs")

        # Post-processing steps
        read_steps(f, raw)

    return raw
