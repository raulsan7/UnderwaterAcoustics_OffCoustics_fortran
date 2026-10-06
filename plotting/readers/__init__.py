"""
Module: readers/__init__.py
Description: Readers of the .hdf5 files written by the Fortran code. They return a raw dictionary in the
             vocabulary of the emitter (dataset names as Fortran wrote them, Fortran order restored) and do not
             interpret it: the wrappers turn it into the Case the plots work with.

             hdf5_io   datasets (read_ds, read_complex, ...) and the emitter of a file
             solver    datasets written by the solver in every file (method, environment, post-processing steps)
             turbine   single wind turbine (and the per-turbine datasets of a farm)
             farm      wind farm
             source    frequency-domain point sources (FrequencySource)
             describe  one-line description of any file for --list

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""
