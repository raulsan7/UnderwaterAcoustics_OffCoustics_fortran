"""
Module: plots/specific/sphere.py
Description: Radiation metrics of the 'sphere' step (free field, Image Method with N_images = 0). The power per
             wetted area needs the wetted area of the structure (case["specific"]["wet_area"], saved by every
             turbine type); turbine files without it fall back to the values of the DTU 10 MW turbines.

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import utils.AcousticUtils as au
from wrappers.case import emitter_label
from plots.common import no_absorption
from plot_utils.conversions import is_valid


def compute_sphere_metrics(case            : dict  = None,         # [-] Case
                           absorption      : bool  = True,         # [-] Whether to apply absorption attenuation for OASPL
                           filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                           filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                           rho_water       : float = 1025.0,       # [kg/m³] Water density
                           c_water         : float = 1500.0,       # [m/s] Speed of sound in water
                           verbose         : bool  = True):        # [-] Print metrics to console
    """
    Computes and optionally prints radiation metrics from a spherical pressure field. Observers without data
    (NaN, inside the structure) are left out and the excluded fraction of the sphere area is reported. The
    power per wetted area is left out for emitters without a wetted area.

    Returns
    -------
    metrics : dict
        Dictionary containing all computed acoustic metrics (Power, SL, etc.)
    """

    if not case["has_sphere"]:
        print("compute_sphere_metrics(): no sphere data found in case, skipping.")
        return None

    absorption = no_absorption(case, absorption, "compute_sphere_metrics")

    rho_wat = case.get("rho_wat", rho_water)
    c_wat   = case.get("c_wat", c_water)
    tag     = emitter_label(case)

    # Extract data
    p       = case["P_sphere"]
    Freqs   = case["Freqs"]
    r       = case["R_sphere"]
    dA      = case["dA_sphere"]
    p_ref   = case["p_ref"]
    n_theta = case["N_theta"]
    nz      = case["Nz_sphere"]

    if not is_valid(p, "compute_sphere_metrics"):
        return None

    # Wetted area saved by the Fortran code (save_specific_params of every turbine type)
    A_wet = case["specific"].get("wet_area")
    if A_wet is None and case["Emitter"] == "turbine":
        A_wet = 848.12 if "Monopile" in tag else 9359.07
        print(f"WARNING compute_sphere_metrics(): wet_area not found in the file, using the hardcoded {A_wet} m²")

    # Filter frequencies
    Nfreqs = len(Freqs)
    freq_mask = np.ones(Nfreqs, dtype=bool)
    if filter_under is not None:
        freq_mask &= (Freqs >= filter_under)
    if filter_over is not None:
        freq_mask &= (Freqs <= filter_over)
    if not np.any(freq_mask):
        print("compute_sphere_metrics(): no frequencies match the specified filter bounds!")
        return None

    p, Freqs  = p[freq_mask], Freqs[freq_mask]
    distances = np.full((n_theta, nz), r)

    # --- 1. Compute OASPL and Directional SL --- #
    oaspl = au.pressure_to_OASPL(p, p_ref, absorption, Freqs, distances, enbw=case["ENBW"])
    sl    = oaspl + 20 * np.log10(r)

    # --- 2. Compute Radiated Power (W) --- #
    # Observers without data (NaN) are left out of the surface integral (equal dA for every observer)
    p_sq  = np.sum(np.abs(p)**2, axis=0)
    valid = np.isfinite(p_sq)
    excluded = 1.0 - np.count_nonzero(valid) / valid.size
    I = p_sq / (2.*rho_wat*c_wat)
    W = np.sum(I[valid]*dA)
    W_dB = 10* np.log10(W/1e-12)

    # --- 3. Compute Power per Wet Area (emitters with a wetted area) --- #
    if A_wet is not None:
        W_per_wet = W / A_wet
        W_per_wet_dB = 10.0 * np.log10(W_per_wet / 1e-12)
    else:
        W_per_wet = W_per_wet_dB = None

    # --- 4. Equivalent Omnidirectional SL --- #
    W_ref_pressure = (4.0 * np.pi * p_ref**2) / (rho_wat * c_wat)
    SL_omni = 10.0 * np.log10(W / W_ref_pressure)

    # --- 5. Directional SL Statistics --- #
    SL_max = np.nanmax(sl)
    SL_min = np.nanmin(sl)
    SL_mean_power = 10.0 * np.log10(np.nanmean(10.0**(sl / 10.0))) # Energetic mean

    # --- 6. Mean SPL on the sphere --- #
    # p is a peak amplitude (p_rms^2 = |p|^2 / 2), as in the OASPL of pressure_to_SPL
    SPL_sphere = 10.0 * np.log10(np.mean(p_sq[valid]) / (2.0 * p_ref**2))

    # --- Package results --- #
    metrics = {
        "OASPL_sphere": oaspl,
        "SL_sphere": sl,
        "W": W,
        "W_dB": W_dB,
        "W_per_wet": W_per_wet,
        "W_per_wet_dB": W_per_wet_dB,
        "SL_omni": SL_omni,
        "SL_max": SL_max,
        "SL_min": SL_min,
        "SL_mean_power": SL_mean_power,
        "SPL_sphere_mean": SPL_sphere,
        "Excluded_area_fraction": excluded
    }

    # --- Print formatting --- #
    if verbose:
        print("\n" + "=" * 80)
        print(f"  RADIATION METRICS (SPHERE r = {r:.0f} m, free-field)")

        if filter_under is not None or filter_over is not None:
            lo = f"{filter_under:.1f}" if filter_under is not None else "0"
            hi = f"{filter_over:.1f}"  if filter_over is not None else "∞"
            print(f"  Frequency band: {lo} - {hi} Hz")
        print("=" * 80)

        label_w = 42
        print(f"\n[Case: {tag}]")
        print(f"  {'Radius of sphere':<{label_w}}: {r:.1f} m")
        if A_wet is not None:
            print(f"  {'Wet area of structure':<{label_w}}: {A_wet:.1f} m²")
        if excluded > 0:
            print(f"  {'Excluded area (observers without data)':<{label_w}}: {100*excluded:.2f} % "
                  f"({valid.size - np.count_nonzero(valid)} of {valid.size} observers)")
        print(f"  {'Total radiated power (W)':<{label_w}}: {W:.4e} W")
        print(f"  {'Total radiated power (dB)':<{label_w}}: {W_dB:.1f} dB re 1 pW")
        if A_wet is not None:
            print(f"  {'Power / wet area (W)':<{label_w}}: {W_per_wet:.4e} W/m²")
            print(f"  {'Power / wet area (dB)':<{label_w}}: {W_per_wet_dB:.1f} dB re 1 pW/m²")
        print(f"  {'Equivalent omnidirectional SL':<{label_w}}: {SL_omni:.1f} dB re 1 µPa @ 1 m")
        print(f"  {'Directional SL MAX':<{label_w}}: {SL_max:.1f} dB re 1 µPa @ 1 m")
        print(f"  {'Directional SL MIN':<{label_w}}: {SL_min:.1f} dB re 1 µPa @ 1 m")
        print(f"  {'Directional SL mean (energ.)':<{label_w}}: {SL_mean_power:.1f} dB re 1 µPa @ 1 m")
        print(f"  {'Mean SPL on sphere':<{label_w}}: {SPL_sphere:.1f} dB re 1 µPa")

    return metrics
