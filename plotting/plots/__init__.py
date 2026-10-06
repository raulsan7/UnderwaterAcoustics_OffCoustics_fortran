"""
Module: plots/__init__.py
Description: Plots of one case, for any emitter: one generic function per post-processing step (spectrum,
             polar, cylinder, decay, line, slices, cylinder metrics and equivalent point source) and the
             exclusive functions of plots/specific (structure, sphere metrics, sources). They read a Case
             (contract of wrappers/case.py) and never a reader or a wrapper. Imported as 'p' by the entry
             scripts:

             import plots as p
             p.plot_spectrum(case, mode='PSD', octave=False)

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

from plots.spectrum import get_spectrum, plot_spectrum
from plots.polar import get_bands, plot_polar
from plots.cylinder import plot_cylinder
from plots.decay import get_curve, get_family_curves, plot_distance_decay
from plots.line import plot_line
from plots.slices import plot_sliceXY, plot_sliceXZ, plot_sliceVertical
from plots.metrics import compute_cylinder_metrics, compute_point_source
from plots.specific.structure import plot_structure
from plots.specific.sphere import compute_sphere_metrics
from plots.specific.sources import plot_sources
