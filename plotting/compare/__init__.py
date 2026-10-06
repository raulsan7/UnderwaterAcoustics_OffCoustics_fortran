"""
Module: compare/__init__.py
Description: Comparison of two or more cases of any emitter: one panel per case, drawn by the functions of
             plots/ through ax_in (compare/plots.py), with the panels, titles and legends of compare/panels.py.
             Imported as 'pc' by the entry scripts:

             import compare as pc
             pc.plot_spectrum(cases, mode='PSD', octave=False)

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

from compare.panels import MAX_COLS, select_cases, make_panels, remove_common_title_part, wrap_case_titles, add_case_legend
from compare.plots import (plot_structure, plot_spectrum, plot_polar, plot_cylinder, plot_distance_decay, plot_line,
                           plot_sliceXY, plot_sliceXZ, plot_sliceVertical, compute_sphere_metrics,
                           compute_cylinder_metrics, compute_point_source)
