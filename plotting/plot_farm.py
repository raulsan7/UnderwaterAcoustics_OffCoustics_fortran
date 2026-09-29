"""
Plots wind farms (.hdf5 files written by Fortran in FORTRAN/farm_acoustic_data).

    python plotting/plot_farm.py                                   # most recent file
    python plotting/plot_farm.py Farm_2_DTU10MN_ANM24_ws8.0-11.4   # given file
    python plotting/plot_farm.py --list                            # available files
    python plotting/plot_farm.py ... --save [--no-show]            # figures saved in FORTRAN/figures/<case>/

The plots shown are selected below.
"""

import matplotlib.pyplot as plt
import plot_utils as pu
import plot_farm_module as p
from plot_utils import parse_farm_args, load_farm
from plot_style import AcademicStyle


def main():

    # 1. Visual asthetics
    AcademicStyle.apply()

    args = parse_farm_args()
    case = load_farm(args, verbose=False)

    # Select plot
    p.plot_structure(case, mode='xy')
    p.plot_spectrums(case, do_thresholds=True)
    p.plot_polar(case)
    p.plot_cylinder(case)
    p.plot_line(case)
    p.plot_sliceXY(case, mode='SPL', target_frequency=12.)
    p.plot_sliceXZ(case, mode='SPL', target_frequency=12.)
    p.plot_sliceVertical(case, mode='SPL', target_frequency=12.)

    # 2. Save (only with --save) and show
    pu.save_figures(args, case)
    if not args.no_show:
        plt.show()



if __name__ == "__main__":
    main()
