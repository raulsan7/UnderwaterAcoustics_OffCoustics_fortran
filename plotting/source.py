"""
Plots point sources (.hdf5 files written by Fortran in FORTRAN/source_acoustic_data: emitter FrequencySource,
dipoles or monopoles at given positions forced at given tones, and emitter TimeSource, forced by a force series
in time).

    python plotting/source.py                            # most recent file
    python plotting/source.py Source_Ns2_Nf3_ANM9        # one file: single case
    python plotting/source.py Source_Ns1_Nt2000_ANM8     # a TimeSource file
    python plotting/source.py a.hdf5 b.hdf5 ...          # several files: comparison
    python plotting/source.py --list                     # available files
    python plotting/source.py ... --save [--no-show]     # figures saved in FORTRAN/figures/<case>/

The frequencies of a FrequencySource are tones, not the bins of an FFT: there is a PSD only if they are uniformly
spaced (otherwise its plots are skipped with a message). Those of a TimeSource are the bins of an FFT with a
Hann window, so it has a PSD, as a turbine has.
The plots shown are selected below.
"""

import matplotlib.pyplot as plt
import plots as p
import compare as pc
from plot_utils.style import AcademicStyle
from plot_utils.saving import save_figures
from wrappers.source import parse_source_args, load_sources

# Single-frequency plots (line and XY slice) at the frequencies of Figs. 5 and 6 of Deng et al. (the nearest
# frequency of the file is taken), with the colormap of its Fig. 5. 100 Hz, also in those figures, is left out:
# it is not a line of its source spectrum.
TARGET_FREQUENCIES = [36.0, 48.0, 72.0]     # [Hz]
SLICE_CMAP         = 'jet'                  # [-] Colormap of the XY slices

def main():
    # 1. Visual asthetics
    AcademicStyle.apply()

    args  = parse_source_args()
    cases = load_sources(args, verbose=False)

    if len(cases) == 1:
        case = cases[0]
        p.plot_sources(case)
        p.plot_spectrum(case, do_thresholds=False, absorption=False,)
        p.plot_spectrum(case, do_thresholds=False, absorption=False, octave=False, smooth=True)
        p.plot_spectrum(case, do_thresholds=False, absorption=False, octave=False, smooth=False)
        p.plot_polar(case, bands = [])
        p.plot_cylinder(case)
        p.plot_distance_decay(case)
        for f in TARGET_FREQUENCIES:
            p.plot_line(case, mode='SPL', target_frequency=f)
            p.plot_sliceXY(case, mode='SPL', target_frequency=f, cmap=SLICE_CMAP)
        p.plot_sliceXZ(case)
        p.plot_sliceVertical(case)
        p.compute_sphere_metrics(case)
        p.compute_cylinder_metrics(case)
        p.compute_point_source(case)
    else:
        pc.plot_spectrum(cases, do_thresholds=False, absorption=False, filter_under=10.0, mode='SPL')
        pc.plot_spectrum(cases, do_thresholds=False, absorption=False, filter_over=10.0, octave=False, mode='PSD', smooth=False)
        # pc.plot_polar(cases)
        # pc.plot_cylinder(cases)
        # pc.plot_distance_decay(cases)
        # pc.plot_line(cases)
        # pc.plot_sliceXY(cases, mode='SPL', target_frequency=50.)
        # pc.plot_sliceXZ(cases, mode='SPL', target_frequency=50.)
        # pc.plot_sliceVertical(cases)
        # pc.compute_sphere_metrics(cases)
        # pc.compute_cylinder_metrics(cases)
        # pc.compute_point_source(cases)


    # 2. Save (only with --save) and show
    save_figures(args, cases)
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
