"""
Module: plot_utils/cli.py
Description: Arguments shared by the entry scripts (turbine.py, farm.py, source.py). Every wrapper builds its
             own parser (description, files, flags of its emitter) on top of them.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import sys
import argparse
from pathlib import Path
from plot_utils.files import FIGURES_DIR, list_hdf5, rel


def add_common_args(parser  : argparse.ArgumentParser = None,   # [-] Parser of an entry script
                    data_dir: Path = None):                     # [-] Default folder of the .hdf5 files
    """
    Adds the arguments shared by every entry script.
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

def finish_args(args    : argparse.Namespace = None,   # [-] Parsed command-line arguments
                describe       = None):                # [-] Function (open h5py.File) -> (description, steps), for --list
    """
    Converts the path arguments into Path objects and serves --list (prints the available files and exits).
    """

    args.data_dir   = Path(args.data_dir)
    args.output_dir = Path(args.output_dir)

    if args.list:
        list_hdf5(args.data_dir, describe)
        sys.exit(0)

    return args
