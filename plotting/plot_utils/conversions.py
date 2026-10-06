"""
Module: plot_utils/conversions.py
Description: Conversions of the complex pressure into the magnitude of a plot (SPL, PSD, ABS, REAL, IMAG,
             PHASE), frequency bin width of an FFT grid, and the check of the pressure arrays before plotting.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import utils.AcousticUtils as au


def get_bin_width(freqs):               # [Hz] Frequency array shape(NFreqs,)
    """
    Gets the frequency bin width df [Hz] needed to convert pressure into PSD.

    Freqs lies on the uniform FFT grid but may have gaps (compute_force(filter_freqs=True)
    drops non-excited bins), so df is the smallest spacing, not the local one. If the
    spectrum was decimated (skipf > 1) the true df cannot be recovered from Freqs.
    Returns None if it cannot be determined (less than two frequencies).
    """

    freqs = np.unique(np.asarray(freqs, dtype=float))
    if freqs.size < 2:
        return None

    return np.min(np.diff(freqs))

def get_uniform_bin_width(freqs,                # [Hz] Frequency array shape(NFreqs,)
                          rtol: float = 1e-6):  # [-] Relative tolerance of the spacing
    """
    Gets the bin width df [Hz] of frequencies given one by one (e.g. the tones of frequency-domain sources):
    their spacing if it is uniform, so that every frequency stands for a band of that width. Frequencies that
    are not uniformly spaced have no bin width: the PSD and the 1/3 octave bands are not applicable to them.
    FFT grids, whose gaps (filter_freqs=True) keep the bin width, use get_bin_width.

    Returns
    -------
    df : float or None
        None with less than two frequencies or a non-uniform spacing.
    """

    freqs = np.unique(np.asarray(freqs, dtype=float))
    if freqs.size < 2:
        return None

    spacing = np.diff(freqs)
    if not np.allclose(spacing, spacing[0], rtol=rtol, atol=0.0):
        return None

    return np.min(spacing)

def convert_to(p,                               # [Pa] Pressure array shape(NFreqs, Nobs)
               mode: str   = "SPL",             # [-] What to convert to
               df  : float = None,              # [Hz] Frequency bin width, only needed for 'PSD'
               enbw: float = au.HANN_ENBW):     # [bins] ENBW of the FFT window (case["ENBW"]), only used by 'PSD'
    """
    Converts pressure into desired units.
    """

    mode = mode.upper()

    if mode == "SPL":
        vals = au.pressure_to_SPL(p)
        unit_label = r"SPL [dB re 1 $\mu$Pa]"
    elif mode == "PSD":
        if df is None: raise ValueError("convert_to(): mode 'PSD' needs the frequency bin width 'df'")
        vals = au.pressure_to_PSD(p, df, enbw=enbw)
        unit_label = r"PSD [dB re 1 $\mu$Pa$^2$/Hz]"
    elif mode == "ABS":
        vals = np.abs(p)
        unit_label = r"Magnitude [Pa]"
    elif mode == "REAL":
        vals = np.real(p)
        unit_label = r"Real Part [Pa]"
    elif mode == "IMAG":
        vals = np.imag(p)
        unit_label = r"Imaginary Part [Pa]"
    elif mode == "PHASE":
        vals = np.angle(p, deg=True)
        if np.any(vals < -90) and np.any(vals > 90):
            vals = np.mod(vals, 360)
        unit_label = r"Phase [deg]"
    else:
        print(f"convert_to(): Unrecognized mode '{mode}', defaulting to 'SPL'")
        vals = au.pressure_to_SPL(p)
        unit_label = r"SPL [dB re 1 $\mu$Pa]"

    return vals, unit_label

def is_valid(array,                 # [Pa] Pressure array shape(NFreqs, ...)
             name: str = None):     # [-] Name of the calling function, for the messages
    """
    Returns True if the pressure array can be plotted: it exists, has no +-inf entries and at least one
    finite value.

    NaN entries are valid on purpose: the Fortran code writes NaN (real and imaginary parts, at every
    frequency) at the observers lying closer than 5 m to a source node, i.e. inside the structure. They
    carry no data and every plot shows them as holes. When there are some, a message tells how many.
    """

    if array is None:
        return False

    a = np.asarray(array)
    if a.size == 0:
        if name is not None: print(f"{name}(): empty pressure array, skipping")
        return False

    if np.any(np.isinf(a)):
        if name is not None: print(f"{name}(): inf found in pressure array, skipping")
        return False

    nan = np.isnan(a)
    if np.all(nan):
        if name is not None: print(f"{name}(): no observer has data (NaN everywhere), skipping")
        return False

    if name is not None and np.any(nan):
        masked = np.all(nan.reshape(a.shape[0], -1), axis=0) if a.ndim > 1 else nan
        print(f"{name}(): {np.count_nonzero(masked)} of {masked.size} observers without data "
              f"(NaN, inside the structure)")

    return True
