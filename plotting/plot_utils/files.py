"""
Module: plot_utils/files.py
Description: Default folders and discovery of the .hdf5 files written by the Fortran code: newest file of a
             folder, file names given in the terminal turned into paths, and the --list of a folder.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import h5py
from pathlib import Path
from datetime import datetime


# ---------- DEFAULT PATHS ---------- #
FORTRAN_DIR      = Path(__file__).resolve().parent.parent.parent   # [-] FORTRAN/ (this file lives in FORTRAN/plotting/plot_utils/)
TURBINE_DATA_DIR = FORTRAN_DIR / "turbine_acoustic_data"          # [-] Default folder of the single-turbine files
FARM_DATA_DIR    = FORTRAN_DIR / "farm_acoustic_data"             # [-] Default folder of the wind farm files
SOURCE_DATA_DIR  = FORTRAN_DIR / "source_acoustic_data"           # [-] Default folder of the frequency-domain source files
FIGURES_DIR      = FORTRAN_DIR / "figures"                        # [-] Default folder of the saved figures


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

def list_hdf5(data_dir: Path = None,     # [-] Folder of the .hdf5 files
              describe       = None):    # [-] Function (open h5py.File) -> (description, steps) of a file
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
                description, steps = describe(f)
        except Exception as e:
            description, steps = f"unreadable file ({type(e).__name__})", "-"

        stat   = path.stat()
        date   = datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M")
        marker = "*" if idx == 0 else " "

        print(f"  {marker} {date}  {stat.st_size / 1e6:6.0f} MB  {path.name}")
        print(f"        {description}")
        print(f"        steps: {steps}\n")

    print("  * = file plotted when no file is given")

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
