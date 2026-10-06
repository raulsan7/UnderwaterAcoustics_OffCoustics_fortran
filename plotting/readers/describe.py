"""
Module: readers/describe.py
Description: One-line description of any .hdf5 file written by the Fortran code (single turbine, wind farm or
             frequency-domain sources), for the --list of every entry script.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

from pathlib import Path
from readers.hdf5_io import read_ds, detect_emitter
from readers.solver import STEPS
from readers.turbine import guess_case_type


def describe_hdf5(f):     # [-] Open h5py.File
    """
    Returns
    -------
    description, steps : str, str
        Case, method, wind speed and depth of the file (sources: number of sources and tones instead of the
        wind speed), and the post-processing steps it contains.
    """

    if "Nmodes" in f:
        method = f"ANM{read_ds(f, 'Nmodes')}"
    elif "N_images" in f:
        method = f"IM{read_ds(f, 'N_images')}"
    else:
        method = read_ds(f, "Method") or "unknown method"

    steps   = ", ".join(step for step, key in STEPS.items() if key in f)
    steps   = steps if steps else "none"
    emitter = detect_emitter(f)

    # Frequency-domain sources: no wind, the depth is the one of the solver
    if emitter == "source":
        name  = f"{read_ds(f, 'Emitter')}: {read_ds(f, 'Nsources')} sources, {read_ds(f, 'NFreqs')} tones"
        depth = f" | H = {-read_ds(f, 'Lower_HBC'):.0f} m" if "Lower_HBC" in f else ""
        return f"{name} | {method}{depth}", steps

    if emitter == "farm":
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

    return f"{name} | {method} | {ws} m/s{depth}", steps
