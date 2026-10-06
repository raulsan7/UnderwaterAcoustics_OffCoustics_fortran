"""
Module: readers/hdf5_io.py
Description: Reading of the datasets written by the Fortran code (column-major arrays, complex pairs, padded
             strings, shape-(1,) scalars) and detection of the emitter that wrote a file.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
from pathlib import Path

# Value of the "Emitter" dataset -> emitter type of the plotting package
EMITTERS = {"FrequencySource": "source"}

# Entry script of every emitter type, for the messages of a file given to the wrong script
SCRIPTS = {"turbine": "plotting/turbine.py",
           "farm"   : "plotting/farm.py",
           "source" : "plotting/source.py"}


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
      Other complex datasets are combined by read_complex.
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
        val = combine_complex(val, key)

    # Transpose multidimensional arrays for Column-Major
    if isinstance(val, np.ndarray) and val.ndim >= 2:
        val = val.T

    # Collect scalars: NumPy 2.x no longer allows int()/float() on a shape-(1,) array
    if isinstance(val, np.ndarray) and val.shape == (1,):
        return val.item()

    return val

def combine_complex(val: np.ndarray = None,     # [-] Real dataset as h5py reads it, (real, imaginary) pair on its last axis
                    key: str        = ""):      # [-] Dataset name, for the messages
    """
    Combines the (real, imaginary) pair of the last h5py axis of a complex dataset written by Fortran.

    Returns
    -------
    val : np.ndarray
        Complex array, h5py order without the last axis (not transposed).
    """

    if val.shape[-1] != 2:
        raise ValueError(f"readers.hdf5_io.read_ds(): complex dataset {key} with shape {val.shape} "
                         f"has no (real, imaginary) pair on its last axis")

    return np.ascontiguousarray(val, dtype=np.float64).view(np.complex128)[..., 0]

def read_complex(f,          # [-] Open h5py.File
                 key: str):  # [-] Dataset name
    """
    Reads a complex dataset whose name does not start with "P_" (read_ds only combines those), e.g. the
    forces Source_Force: Fortran complex (N1, N2, ...) is (..., N2, N1, 2) in h5py; the pair is combined and the
    array transposed back to the Fortran order (N1, N2, ...).

    Returns
    -------
    None if the dataset does not exist, complex np.ndarray otherwise.
    """

    if key not in f:
        return None

    val = combine_complex(f[key][()], key)

    return val.T if val.ndim >= 2 else val

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

def detect_emitter(f):     # [-] Open h5py.File
    """
    Emitter type of a file: the "Emitter" dataset if present (e.g. "FrequencySource" -> 'source'); otherwise
    "Nturbines" means a wind farm; otherwise a single turbine (every file written before the Emitter dataset).

    Returns
    -------
    emitter : str
        'turbine', 'farm', 'source', or the value of an unknown Emitter dataset as it is.
    """

    emitter = read_ds(f, "Emitter")
    if emitter is not None:
        return EMITTERS.get(emitter, emitter)

    return "farm" if "Nturbines" in f else "turbine"

def wrong_emitter(path   : Path = None,     # [-] Path of the file
                  emitter: str  = None,     # [-] Emitter type of the file (detect_emitter)
                  caller : str  = None):    # [-] Name of the calling function, for the message
    """
    Error for a file of another emitter type than the one of the script.

    Returns
    -------
    error : ValueError
    """

    names  = {"turbine": "a single-turbine file", "farm": "a wind farm file", "source": "a FrequencySource file"}
    what   = names.get(emitter, f"a file of emitter '{emitter}'")
    script = SCRIPTS.get(emitter)
    hint   = f", plot it with {script}" if script else ", no plotting script reads it yet"

    return ValueError(f"{caller}(): {path.name} is {what}{hint}")
