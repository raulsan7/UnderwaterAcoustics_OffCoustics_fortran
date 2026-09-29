"""
Module: AcousticUtils.py
Description: Acoustic conversions used by the plots (SPL, PSD, OASPL, band OASPL, seawater absorption and
             1/3 octave bands). Copy of PYTHON/utils/AcousticUtils.py, so the plots do not depend on it.
             Masked observers (NaN pressure at every frequency, see plot_utils.is_valid) stay NaN through
             every conversion: frequency sums use np.sum, never np.nansum (which would turn them into a
             fake -60 dB OASPL).

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 09/2026
"""

import numpy as np


def pressure_to_SPL(p, p_ref=1e-6, absorption=False, freq=None, distance=None):
    """
    Convert complex pressure amplitudes (Helmholtz) to SPL [dB re p_ref].

    Parameters
    ----------
    p : array_like
        Complex pressure amplitudes, arbitrary shape. (Nf, -, - ...)
    p_ref : float, optional
        Reference pressure in Pa. Default is 1e-6 Pa (underwater acoustics).
    absorption : bool, optional
        Whether to apply seawater absorption correction. Default is False.
    freq : ndarray, shape (Nf,), optional
        Frequencies in Hz. Required if ``absorption=True``.
    distance : float, optional
        Wave travel distance in m. Required if ``absorption=True``.

    Returns
    -------
    spl : ndarray
        SPL in dB, same shape as ``p``.
    """
    p = np.asarray(p, dtype=complex)
    p_rms = np.abs(p)/np.sqrt(2)            # p_rms(fi)
    p_rms = np.maximum(p_rms, p_ref)

    spl = 20.0 * np.log10(p_rms / p_ref)

    if absorption:
        if freq is None: raise ValueError("frequency array needed if absorption is enabled")
        if distance is None: raise ValueError("Distance is needed if absorption is enabled")
        spl = add_absorption(freq, spl, distance)

    return spl

def SPL_to_OASPL(spl, axis=0):
    """
    Convert narrowband SPL [dB] to OASPL [dB] by summing over frequency.

    Parameters
    ----------
    spl : array_like
        SPL values in dB. Frequency axis must match ``axis``.
    axis : int, optional
        Axis to sum over (frequency axis). Default is 0.

    Returns
    -------
    oaspl : ndarray
        OASPL in dB, with the frequency axis removed.
    """
    spl = np.asarray(spl, dtype=float)

    return 10.0 * np.log10(np.sum(10 ** (spl / 10), axis=axis)+1e-6)

def pressure_to_OASPL(p, p_ref=1e-6, absorption=False, freq=None, distance=None, axis=0):
    """
    Convert complex pressure amplitudes (Helmholtz) to OASPL [dB re p_ref].

    Convenience wrapper around :func:`pressure_to_SPL` and :func:`SPL_to_OASPL`.

    Parameters
    ----------
    p : array_like
        Complex pressure amplitudes. Frequency axis must match ``axis``.
    p_ref : float, optional
        Reference pressure in Pa. Default is 1e-6 Pa.
    axis : int, optional
        Frequency axis. Default is 0.
    absorption : bool, optional
        Whether to apply seawater absorption correction. Default is False.
    freq : ndarray, shape (Nf,), optional
        Frequencies in Hz. Required if ``absorption=True``.
    distance : float, optional
        Wave travel distance in m. Required if ``absorption=True``.

    Returns
    -------
    oaspl : ndarray
        OASPL in dB, with the frequency axis removed.
    """
    p = np.asarray(p)
    spl = pressure_to_SPL(p, p_ref=p_ref, absorption=absorption, freq=freq, distance=distance)

    return SPL_to_OASPL(spl, axis=axis)

def OASPL_band(p, mask, p_ref=1e-6, absorption=False, freq=None, distance=None, axis=0):
    """
    OASPL restricted to a frequency-band boolean mask.

    Parameters
    ----------
    p : array_like
        Complex pressure amplitudes. Frequency axis must match ``axis``.
    mask : array_like of bool, shape (Nf,)
        Boolean mask selecting the frequency band of interest.
    p_ref : float, optional
        Reference pressure in Pa. Default is 1e-6 Pa.
    axis : int, optional
        Frequency axis. Default is 0.
    absorption : bool, optional
        Whether to apply seawater absorption correction. Default is False.
    freq : ndarray, shape (Nf,), optional
        Frequencies in Hz. Required if ``absorption=True``.
    distance : float, optional
        Wave travel distance in m. Required if ``absorption=True``.

    Returns
    -------
    oaspl : ndarray
        Band-limited OASPL in dB, with the frequency axis removed.
    """

    p_band = np.asarray(p)[mask]
    if freq is not None:
        freq = freq[mask]
    else:
        freq = None

    return pressure_to_OASPL(p_band, p_ref=p_ref, axis=axis, absorption=absorption, freq=freq, distance=distance)

def add_absorption(freq, SPL, distance):
    """
    Apply seawater absorption attenuation to SPL values.

    Uses the empirical formula for frequency-dependent absorption
    coefficient alpha(f) [dB/km] and subtracts the total attenuation
    over ``distance``.

    Parameters
    ----------
    freq : array_like
        Frequencies in Hz.
    SPL : array_like
        SPL values in dB, same shape as ``freq``.
    distance : float
        Wave travel distance in m.

    Returns
    -------
    spl_corrected : ndarray
        Absorption-corrected SPL in dB.
    """
    freq = np.asarray(freq)
    distance = np.asarray(distance)
    SPL = np.asarray(SPL)

    if freq.size > 1:
        if distance.ndim > 0 and distance.shape != SPL.shape[1:]:
            raise ValueError(f"Distance shape {distance.shape} must match SPL spatial shape {SPL.shape[1:]}")
    else:
        if distance.shape != SPL.shape:
            raise ValueError(f"Distance shape {distance.shape} must match SPL spatial shape {SPL.shape[1:]}")

    def alpha(f):
        fx = (f/1000)**2
        return 0.0033 + 0.11*fx/(1+fx) + 44*fx/(4100+fx) + 0.0003*fx
    
    # Broadcast
    alpha_vals = alpha(freq)
    for _ in range(distance.ndim):
        alpha_vals = alpha_vals[..., np.newaxis]
    
    return SPL - alpha_vals*distance

def to_third_octave(freqs, spectrum, fmin=20.0, fmax=20000.0, domain='db', linear_mode='power'):
    """
    Compute a 1/3 octave ISO band spectrum from any input spectrum.

    Band limits follow IEC 61260-1:
        f_low(n)  = fc(n) / 2^(1/6)
        f_high(n) = fc(n) * 2^(1/6)

    Bands are contiguous and non-overlapping: the upper edge of band n is
    exactly the lower edge of band n+1.

    Parameters
    ----------
    freqs : array_like
        Frequency vector of the input spectrum (Hz). Any resolution accepted
        (e.g. 1 Hz FFT bins, irregular measurement grid, etc.).
    spectrum : array_like
        Spectrum values corresponding to `freqs`.
    fmin : float
        Lower frequency limit for output bands (Hz). Default: 20 Hz.
    fmax : float
        Upper frequency limit for output bands (Hz). Default: 20000 Hz.
    domain : str
        'db'     → spectrum in dB; energy summation is used.
                   Suitable for: SPL (dB), Lw (dB), PSD (dB/Hz), etc.
        'linear' → spectrum in linear units; summation mode set by
                   `linear_mode`.
                   Suitable for: Pa, Pa²/Hz, m/s², V²/Hz, etc.
    linear_mode : str  (only relevant when domain='linear')
        'power'     → direct sum:        band = Σ xᵢ
                      Use for power/energy quantities: Pa², W/Hz, V²/Hz.
        'amplitude' → quadratic sum:     band = √(Σ xᵢ²)
                      Use for amplitude quantities: Pa, m/s², g.
        'mean'      → arithmetic mean:   band = mean(xᵢ)
                      Use for non-energetic quantities: temperature, etc.

    Returns
    -------
    fc : np.ndarray
        ISO 1/3 octave center frequencies (Hz).
    bands : np.ndarray
        Band values in the same units as the input `spectrum`.

    Raises
    ------
    ValueError
        If `domain` or `linear_mode` are not recognised strings.
    ValueError
        If `freqs` and `spectrum` have different lengths.

    Notes
    -----
    If a band contains no input bins, the value of the nearest bin is used
    as a fallback instead of returning NaN. This can happen when the input
    frequency resolution is coarser than a 1/3 octave band width at low
    frequencies.

    Quick reference
    ---------------
    Input spectrum            domain      linear_mode
    ─────────────────────── ─────────── ─────────────
    SPL in dB, Lw in dB      'db'        —
    PSD in dB/Hz              'db'        —
    Pa² (mean square press.)  'linear'    'power'
    Pa²/Hz (linear PSD)       'linear'    'power'
    Pa (pressure amplitude)   'linear'    'amplitude'
    m/s², g (acceleration)    'linear'    'amplitude'
    Temperature, non-energy   'linear'    'mean'
    """
    if np.any(np.iscomplex(spectrum)): raise ValueError("Spectrum has to be real")

    freqs    = np.asarray(freqs,    dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)

    if freqs.shape[0] != spectrum.shape[0]:
        raise ValueError("`freqs` length must match spectrum's first dimension.")

    if domain not in ('db', 'linear'):
        raise ValueError(f"Unknown domain '{domain}'. Use 'db' or 'linear'.")

    if domain == 'linear' and linear_mode not in ('power', 'amplitude', 'mean'):
        raise ValueError(f"Unknown linear_mode '{linear_mode}'.")

    # Guardar shape original y aplanar todo excepto el eje de frecuencias
    # (nfreqs, d1, d2, ...) → (nfreqs, d1*d2*...)
    original_shape = spectrum.shape          # e.g. (nfreqs, Nnodes, 3)
    trailing_shape  = original_shape[1:]     # e.g. (Nnodes, 3)
    n_cols = int(np.prod(trailing_shape)) if trailing_shape else 1

    spectrum_2d = spectrum.reshape(original_shape[0], n_cols)  # (nfreqs, n_cols)

    _, fc   = generate_iso_third_octave(fmin, fmax)
    half_bw = 2 ** (1 / 6)
    n_bands = len(fc)
    bands   = np.empty((n_bands, n_cols))

    for i, f_center in enumerate(fc):
        f_low  = f_center / half_bw
        f_high = f_center * half_bw
        mask   = (freqs >= f_low) & (freqs < f_high)

        if np.any(mask):
            vals = spectrum_2d[mask, :]           # (bins_in_band, n_cols)
        else:
            idx  = np.argmin(np.abs(freqs - f_center))
            vals = spectrum_2d[[idx], :]          # (1, n_cols)

        if domain == 'db':
            bands[i, :] = 10 * np.log10(np.sum(10 ** (vals / 10), axis=0))
        else:
            if linear_mode == 'power':
                bands[i, :] = np.sum(vals, axis=0)
            elif linear_mode == 'amplitude':
                bands[i, :] = np.sqrt(np.sum(vals ** 2, axis=0))
            else:
                bands[i, :] = np.mean(vals, axis=0)

    # Restaurar shape original: (n_bands, d1, d2, ...)
    # Si el input era 1D (solo frecuencias), devolver 1D
    if trailing_shape:
        out = bands.reshape((n_bands,) + trailing_shape)
    else:
        out = bands[:, 0]   # caso 1D: (nfreqs,) → (n_bands,)

    return fc, out

def generate_iso_third_octave(fmin=0.0, fmax=20000.0):
    """
    Generates central frequencies ofISO 1/3 octave between fmin and fmax
    
    IEC 61260-1: fc(n) = 1000*2^(n/3)

    - fmin: minimum value [Hz]
    - fmax: maximum value [Hz]

    Returns:
        - n_vals: band indices
        - fc    : central frequencies [h]
    """

    n_vals = np.arange(-40, 40)
    fc_all = 1000.0 * 2 ** (n_vals / 3)
    
    mask = (fc_all >= fmin) & (fc_all <= fmax)
    return n_vals[mask], np.asarray(fc_all[mask])

def pressure_to_PSD(p, df, p_ref=1e-6, absorption=False, freq=None, distance=None):
    """
    Convert complex pressure amplitudes to Power Spectral Density (PSD) in dB
    re (p_ref²/Hz).

    The PSD is obtained from the SPL of each frequency bin by subtracting
    ``10·log10(df)``, which accounts for the bin bandwidth. This assumes that
    the input complex pressure ``p`` represents the peak amplitude at each
    frequency (the same convention used by :func:`pressure_to_SPL`).

    Parameters
    ----------
    p : array_like
        Complex pressure amplitudes, arbitrary shape.
        The first axis is assumed to be the frequency axis.
    df : float or array_like
        Frequency bin width in Hz. Must be broadcastable with the frequency
        dimension of ``p`` (i.e., ``df.shape == p.shape[0]`` or scalar).
    p_ref : float, optional
        Reference pressure in Pa. Default 1e-6 Pa (underwater acoustics).
    absorption : bool, optional
        Whether to apply seawater absorption correction. Default False.
    freq : ndarray, shape (Nf,), optional
        Frequencies in Hz. Required if ``absorption=True``.
    distance : float, optional
        Wave travel distance in m. Required if ``absorption=True``.

    Returns
    -------
    psd : ndarray
        PSD in dB re (p_ref²/Hz), same shape as ``p``.

    Notes
    -----
    The conversion is simply ``PSD [dB] = SPL [dB] - 10·log10(df)``, where
    SPL is computed by :func:`pressure_to_SPL`. For broadband signals, this
    yields a density that is independent of the frequency resolution of the
    original FFT.
    """
    # Ensure df is array-like and broadcastable
    df = np.asarray(df)
    p = np.asarray(p, dtype=complex)

    # Compute SPL in dB re p_ref
    spl = pressure_to_SPL(p, p_ref=p_ref, absorption=absorption,
                          freq=freq, distance=distance)

    # Convert to PSD: PSD = SPL - 10*log10(df)
    # Use np.subtract with broadcasting along frequency axis
    # If df is scalar, broadcasting is trivial; if vector, align with first axis
    # Ensure the shape of df matches the frequency axis length of p
    if df.ndim == 1 and df.size == p.shape[0]:
        # Reshape df to allow broadcasting over trailing dimensions
        shape = [-1] + [1] * (p.ndim - 1)
        df = df.reshape(shape)
    # else if scalar, it will broadcast automatically

    psd = spl - 10.0 * np.log10(df)

    return psd
