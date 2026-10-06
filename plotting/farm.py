"""
Plots wind farms (.hdf5 files written by Fortran in FORTRAN/farm_acoustic_data).

    python plotting/farm.py                                   # most recent file
    python plotting/farm.py Farm_2_DTU10MN_ANM24_ws8.0-11.4   # one file: single farm
    python plotting/farm.py a.hdf5 b.hdf5 ...                 # several files: comparison
    python plotting/farm.py --list                            # available files
    python plotting/farm.py ... --save [--no-show]            # figures saved in FORTRAN/figures/<case>/

The plots shown are selected below.
"""

import matplotlib.pyplot as plt
import plot_utils as pu
import plot_farm_module as p
import plot_farm_comparison as pc
from plot_utils import parse_farm_args, load_farms
from plot_style import AcademicStyle


def main():

    # 1. Visual asthetics
    AcademicStyle.apply()

    args  = parse_farm_args()
    cases = load_farms(args, verbose=False)

    # Select plot
    if len(cases) == 1:
        case = cases[0]
        p.plot_structure(case, mode='xy')
        p.plot_spectrums(case, do_thresholds=False, octave=True)
        p.plot_spectrums(case, do_thresholds=False, octave=False, smooth=True, mode='PSD')
        p.plot_polar(case)
        p.plot_cylinder(case)
        p.plot_line(case)
        p.plot_sliceXY(case, mode='SPL', target_frequency=12.)
        p.plot_sliceXY(case)
        p.plot_sliceXZ(case, mode='SPL', target_frequency=12.)
        p.plot_sliceVertical(case, mode='SPL', target_frequency=12.)
        p.compute_cylinder_metrics(case)
        p.compute_point_source(case)
    else:
        pc.plot_spectrums(cases, do_thresholds=False, octave=True, mode='SPL')
        pc.plot_spectrums(cases, do_thresholds=False, octave=False, mode='PSD', smooth=True)
        # pc.plot_structure(cases, mode='xy')
        # pc.plot_polar(cases)
        # pc.plot_cylinder(cases)
        # pc.plot_line(cases)
        # pc.plot_sliceXY(cases, mode='SPL', target_frequency=12.)
        # pc.plot_sliceXZ(cases, mode='SPL', target_frequency=12.)
        # pc.plot_sliceVertical(cases, mode='SPL', target_frequency=12.)
        # pc.compute_cylinder_metrics(cases)
        # pc.compute_point_source(cases)

    # 2. Save (only with --save) and show
    pu.save_figures(args, cases)
    if not args.no_show:
        plt.show()



if __name__ == "__main__":
    main()
