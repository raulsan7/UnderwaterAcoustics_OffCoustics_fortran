"""
Module: wrappers/source.py
Description: Point sources (emitters FrequencySource and TimeSource): command line of plotting/source.py, raw
             dictionary of readers/source.py -> Case, and the labels of a comparison of sources. The frequencies of
             a FrequencySource are tones, so its Case has a bin width only if they are uniformly spaced (otherwise
             no PSD); those of a TimeSource are the bins of an FFT, with the bin width 1 / (Nt dt) and a Hann
             window (ENBW). A single source is the origin of the distances, several sources have none (no
             absorption). No --ws.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import argparse
from pathlib import Path
from readers.source import read_source
from readers.describe import describe_hdf5
from wrappers.case import solver_part, flatten_nodes, validate_case, apply_comparison_labels
from plot_utils.cli import add_common_args, finish_args
from plot_utils.files import SOURCE_DATA_DIR, resolve_paths
from plot_utils.colors import get_type_color
from plot_utils.conversions import get_uniform_bin_width


# ---------- COMMAND LINE ---------- #
def parse_source_args():
    """
    Parses arguments from the terminal.
    """

    parser = argparse.ArgumentParser(description="Plots point sources, FrequencySource or TimeSource (one .hdf5 "
                                                 "file), or compares several of them (two or more files).")
    parser.add_argument(
        "files", type = str, nargs = "*",
        help = "Source .hdf5 files written by Fortran: one file plots a single case, several files are compared. "
               "A bare file name is also searched in --data_dir (default: most recent .hdf5 in --data_dir)"
    )
    add_common_args(parser, SOURCE_DATA_DIR)

    return finish_args(parser.parse_args(), describe_hdf5)


# ---------- CASES ---------- #
def load_source(path   : Path = None,      # [-] Path to the source .hdf5 file
                verbose: bool = False):    # [-] if True, prints a list of keys available in the file
    """
    Loads a point-source .hdf5 file (FrequencySource or TimeSource) coming from Fortran as a Case.

    Returns
    -------
    case : dict
        Case of the contract (wrappers/case.py). The forces of the sources are in case["specific"]["Source_Force"]
        and, for a TimeSource, their series in case["specific"]["Source_Force_time"] with the times Source_time.
    """

    return to_case(read_source(path, verbose=verbose))

def to_case(raw: dict = None):     # [-] Raw dictionary of readers.source.read_source
    """
    Translates the datasets of the sources into the Case contract: label 'FrequencySource' or 'TimeSource',
    nodes = the sources (one group each, multipole order), bin width = the spacing of the tones of a
    FrequencySource if it is uniform (None otherwise: no PSD nor 1/3 octave bands) or the frequency resolution
    1 / (Nt dt) of a TimeSource, origin = the source if there is only one (several sources: no origin, so no
    absorption), centre = mean position, no direction nor references.

    Returns
    -------
    case : dict
    """

    xyz = raw["Source_xyz"]

    case = {"Path"   : raw["Path"],
            "Name"   : raw["Name"],
            "Emitter": "source",
            "label"  : raw["Emitter"]}
    case.update(solver_part(raw))

    # Tones given one by one, or the bins of an FFT (the true resolution, whatever frequencies were kept)
    if raw["Emitter"] == "TimeSource":
        case["bin_width"] = 1.0 / (raw["Nt"] * raw["Time_step"])
    else:
        case["bin_width"] = get_uniform_bin_width(case["Freqs"])

    # One group per source
    case["nodes"], case["node_group"] = flatten_nodes([xyz[i] for i in range(raw["Nsources"])])
    case["group_labels"] = [f"Source {i}" for i in range(1, raw["Nsources"] + 1)]
    case["node_order"]   = raw["Source_order"][case["node_group"]]

    # Reference points: a single source is the origin of the distances; the centre is the mean position
    case["origin"] = xyz[0, :2].copy() if raw["Nsources"] == 1 else None
    case["centre"] = xyz[:, :2].mean(axis=0)

    case["direction"]  = None
    case["references"] = None
    case["Color"]      = get_type_color(case["label"])
    case["specific"]   = {key: raw[key] for key in ("Emitter", "Nsources", "NFreqs", "Source_xyz", "Source_order",
                                                    "Source_Force", "Nt", "Time_step", "Source_time",
                                                    "Source_Force_time") if key in raw}

    return validate_case(case)

def load_sources(args   : argparse.Namespace = None,   # [-] Parsed command-line arguments
                 verbose: bool = False):               # [-] if True, prints a list of keys available in the file
    """
    Loads the sources given in the terminal, or all of them to compare them.

    Returns
    -------
    cases : list of dict
        One Case per file. With several files every case gets a distinct label and color.
    """

    paths = resolve_paths(args.files, args.data_dir)
    cases = [load_source(path, verbose=verbose or args.verbose) for path in paths]

    if len(cases) > 1:
        set_source_comparison_labels(cases)

    return [validate_case(case) for case in cases]


# ---------- LABELS ---------- #
def set_source_comparison_labels(cases: list = None):     # [-] List of Cases of sources
    """
    Gives every case of a comparison a distinct label, used by titles, built from its number of sources, its
    tones or the samples of its force series, and its method (e.g. '2 sources, 3 tones, ANM9' or
    '1 sources, 2000 samples, ANM8'). The water depth is added if two labels match, and the file name if they
    still do. Colors are assigned by position.
    """

    labels = []
    for case in cases:
        spec  = case["specific"]
        count = f"{spec['Nt']} samples" if "Nt" in spec else f"{spec['NFreqs']} tones"
        labels.append(f"{spec['Nsources']} sources, {count}, {case['Method_tag']}")

    apply_comparison_labels(cases, labels)
