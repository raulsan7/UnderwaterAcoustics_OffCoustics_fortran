"""
Plots single wind turbines (.hdf5 files written by Fortran in FORTRAN/turbine_acoustic_data).

    python plotting/plot_turbine.py                          # most recent file
    python plotting/plot_turbine.py plot_mn_ANM24_ws11.4     # one file: single case
    python plotting/plot_turbine.py a.hdf5 b.hdf5 ...        # several files: comparison
    python plotting/plot_turbine.py --list                   # available files
    python plotting/plot_turbine.py ... --save [--no-show]   # figures saved in FORTRAN/figures/<case>/

The plots shown are selected below.
"""

import matplotlib.pyplot as plt
import plot_utils as pu
import plot_single_case as p
import plot_comparison as pc
from plot_style import AcademicStyle
from plot_utils import parse_turbine_args, load_turbines

def main():
    # 1. Visual asthetics
    AcademicStyle.apply()

    args  = parse_turbine_args()
    cases = load_turbines(args, verbose=False)

    if len(cases) == 1:
        case = cases[0]
        # p.plot_structure(case, mode="xz")
        p.plot_spectrum(case, do_thresholds=False, absorption=False, filter_under=10.0, mode='PSD')
        p.plot_spectrum(case, do_thresholds=False, absorption=False, filter_under=10.0, octave=False)
        p.plot_polar(case)
        p.plot_cylinder(case)
        p.plot_distance_decay(case, turbine_data=True)
        p.plot_line(case)
        p.plot_sliceXY(case)
        p.plot_sliceXZ(case)
        p.plot_sliceVertical(case)
        p.compute_sphere_metrics(case)
        p.compute_cylinder_metrics(case)
    else:
        pc.plot_spectrum(cases, do_thresholds=False, absorption=False, filter_under=10.0, mode='PSD')
        pc.plot_spectrum(cases, do_thresholds=False, absorption=False, filter_under=10.0, octave=False)
        pc.plot_polar(cases)
        pc.plot_cylinder(cases)
        pc.plot_distance_decay(cases, turbine_data=True)
        pc.plot_line(cases)
        pc.plot_sliceXY(cases)
        pc.plot_sliceXZ(cases)
        pc.plot_sliceVertical(cases)
        # pc.compute_sphere_metrics(cases, filter_under=1.0)
        # pc.compute_cylinder_metrics(cases, filter_under=1.0)

    # 2. Save (only with --save) and show
    pu.save_figures(args, cases)
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
