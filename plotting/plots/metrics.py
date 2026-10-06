"""
Module: plots/metrics.py
Description: Radiation metrics of any emitter from the 'cylinder' step (power flow through the cylinder) and its
             reduction to an equivalent point source (sphere and/or cylinder).

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""

import numpy as np
import matplotlib.pyplot as plt
import utils.AcousticUtils as au
from wrappers.case import emitter_label, case_color
from plot_utils.colors import get_color_shades, get_line_styles
from plot_utils.conversions import is_valid
from plot_utils.saving import set_figure_name

# NumPy >= 2.0 renamed trapz to trapezoid
trapezoid = getattr(np, "trapezoid", None) or np.trapz


def compute_cylinder_metrics(case            : dict  = None,         # [-] Case
                             filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                             filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                             rho_water       : float = 1025.0,       # [kg/m³] Water density
                             c_water         : float = 1500.0,       # [m/s] Speed of sound in water
                             verbose         : bool  = True):        # [-] Print metrics to console
    """
    Computes and optionally prints propagation metrics from a cylindrical pressure field. Observers without
    data (NaN, inside the structure) are left out and the excluded fraction of the lateral area is reported.
    If the cylinder encloses the whole emitter (e.g. every turbine of a farm), the power flow is that of the
    emitter as a whole.

    Returns
    -------
    metrics : dict
        Dictionary containing all computed acoustic metrics (Power, SPL, TL, etc.)
    """

    if not case["has_cylinder"]:
        print("compute_cylinder_metrics(): no cylinder data, skipping")
        return None

    rho_wat = case.get("rho_wat", rho_water)
    c_wat   = case.get("c_wat", c_water)
    tag     = emitter_label(case)

    # Extract data
    p       = case["P_cylinder"]
    Freqs   = case["Freqs"]
    r       = case["R_cylinder"]
    z       = case["Z_cylinder"]
    theta   = case["Theta_deg_cylinder"]
    p_ref   = case["p_ref"]

    if not is_valid(p, "compute_cylinder_metrics"):
            return None

    # Ensure theta is in radians for integration
    theta_rad = theta if np.max(theta) <= 2*np.pi + 0.1 else np.deg2rad(theta)

    # Filter frequencies
    Nfreqs = len(Freqs)
    freq_mask = np.ones(Nfreqs, dtype=bool)
    if filter_under is not None:
        freq_mask &= (Freqs >= filter_under)
    if filter_over is not None:
        freq_mask &= (Freqs <= filter_over)
    if not np.any(freq_mask):
        print("compute_cylinder_metrics(): no frequencies match the specified filter bounds!")
        return None

    p, Freqs = p[freq_mask], Freqs[freq_mask]

    # --- 1. Intensity and Integration Setup --- #
    p_sq = np.sum(np.abs(p)**2, axis=0)      # Sum energy over frequencies: shape (n_theta, nz)
    I = p_sq / (2.0 * rho_wat * c_wat)

    A = 2 * np.pi * r * np.max(np.abs(z))    # Lateral area of the cylinder

    # Observers without data (NaN) are left out of the surface integral
    valid = np.isfinite(I)
    I     = np.where(valid, I, 0.0)

    # Close the cylinder: the azimuthal integral covers the full 360 deg (last point = first point)
    theta_closed = np.append(theta_rad, theta_rad[0] + 2*np.pi)
    I_closed     = np.append(I, I[:1, :], axis=0)
    w_closed     = np.append(valid, valid[:1, :], axis=0).astype(float)

    # --- 2. Compute Radiated Power (W) --- #
    I_vs_z = trapezoid(I_closed * r, theta_closed, axis=0)
    W = np.abs(trapezoid(I_vs_z, z))          # Absolute to handle z-axis integration direction safely
    W_dB = 10 * np.log10(W / 1e-12)

    # Excluded fraction of the lateral area, with the same integration weights
    A_valid  = np.abs(trapezoid(trapezoid(w_closed * r, theta_closed, axis=0), z))
    A_total  = np.abs(trapezoid(trapezoid(np.ones_like(w_closed) * r, theta_closed, axis=0), z))
    excluded = 1.0 - A_valid / A_total

    # --- 3. Compute Power per Cylinder Area --- #
    W_per_area = W / A
    W_per_area_dB = 10 * np.log10(W_per_area / 1e-12)

    # --- 4. Mean SPL on the cylinder --- #
    # p is a peak amplitude (p_rms^2 = |p|^2 / 2), as in the OASPL of pressure_to_SPL
    p2_avg = np.mean(p_sq[valid])
    SPL_avg = 10 * np.log10(p2_avg / (2.0 * p_ref**2))

    # --- 5. Apparent Omnidirectional SL (from Far Field) --- #
    W_ref_press = (4.0 * np.pi * p_ref**2) / (rho_wat * c_wat)
    SL_omni_far = 10.0 * np.log10(W / W_ref_press)

    # --- 6. Comparative Metrics (if Sphere data exists) --- #
    TL = None
    power_loss = None
    TL_str = ""
    power_loss_str = ""

    # --- Package results ---
    metrics = {
        "W_cyl": W,
        "W_cyl_dB": W_dB,
        "W_cyl_per_area": W_per_area,
        "W_cyl_per_area_dB": W_per_area_dB,
        "SPL_cyl_mean": SPL_avg,
        "SL_omni_far": SL_omni_far,
        "TL_cyl": TL,
        "Power_loss_cyl": power_loss,
        "Excluded_area_fraction": excluded
    }

    # --- Print formatting ---
    if verbose:
        print("\n" + "=" * 80)
        print(f"  PROPAGATION METRICS (CYLINDER r = {r:.0f} m)")

        if filter_under is not None or filter_over is not None:
            lo = f"{filter_under:.1f}" if filter_under is not None else "0"
            hi = f"{filter_over:.1f}"  if filter_over is not None else "∞"
            print(f"  Frequency band: {lo} - {hi} Hz")
        print("=" * 80)

        label_w = 44
        print(f"\n[Case: {tag}]")
        print(f"  {'Radius of cylinder':<{label_w}}: {r:.1f} m")
        print(f"  {'Height of cylinder (H)':<{label_w}}: {np.max(np.abs(z)):.1f} m")
        print(f"  {'Lateral area of cylinder':<{label_w}}: {A:.1f} m²")
        if excluded > 0:
            print(f"  {'Excluded area (observers without data)':<{label_w}}: {100*excluded:.2f} % "
                  f"({valid.size - np.count_nonzero(valid)} of {valid.size} observers)")
        print(f"  {'Total power flow (W)':<{label_w}}: {W:.4e} W")
        print(f"  {'Total power flow (dB)':<{label_w}}: {W_dB:.1f} dB re 1 pW")
        print(f"  {'Power / cylinder area (W)':<{label_w}}: {W_per_area:.4e} W/m²")
        print(f"  {'Power / cylinder area (dB)':<{label_w}}: {W_per_area_dB:.1f} dB re 1 pW/m²")
        print(f"  {'Mean SPL on cylinder':<{label_w}}: {SPL_avg:.1f} dB re 1 µPa")
        print(f"  {'Apparent omnid. SL (from far-field)':<{label_w}}: {SL_omni_far:.1f} dB re 1 µPa @ 1 m")

        if TL_str:
            print(f"  {'Transmission Loss (TL)':<{label_w}}: {TL_str}")
        if power_loss_str:
            print(f"  {'Power loss (W_sphere / W_cyl)':<{label_w}}: {power_loss_str}")

    return metrics

def compute_point_source(case            : dict  = None,         # [-] Case
                         filter_under    : float = None,         # [Hz] Lower frequency cutoff (f >= filter_under)
                         filter_over     : float = None,         # [Hz] Upper frequency cutoff (f <= filter_over)
                         enbw            : float = None,         # [bins] Equivalent noise bandwidth of the FFT window (None: case["ENBW"], from the file)
                         smooth_fraction : int   = 12,           # [-] Plots only: narrowband curves smoothed over 1/N octave (None: raw)
                         rho_water       : float = 1025.0,       # [kg/m³] Water density (if not in the file)
                         c_water         : float = 1500.0,       # [m/s] Speed of sound in water (if not in the file)
                         verbose         : bool  = True,         # [-] Print totals and the 1/3 octave table
                         do_plot         : bool  = True,         # [-] Plot SL spectra, radiation diagnostics and PL
                         figsize         : tuple = (10, 5)):     # [-] Size of the SL figure
    """
    Reduces the emitter to an equivalent point source: the free-field monopole that radiates the same power,
    SL = 10 log10(rho c W / (4 pi p_ref^2 r_ref^2)) [dB re 1 uPa^2 m^2], with r_ref = 1 m. Two powers are used:

    - Free field (sphere, Image Method with N_images = 0): power of the source itself -> SL_free.
    - In situ (cylinder, waveguide): power that the emitter injects in the actual water column -> SL_insitu.
      With lossless boundaries every watt radiated by the nodes crosses the cylinder, whatever its radius.
      A wind farm has no sphere: only the in situ level, from the cylinder that must enclose every turbine.

    Both powers come from W = sum |p|^2 / (2 rho c) dA (far-field plane-wave intensity). Known biases, quantified
    with exact dipole power on the SAITEC geometry (10/2026): the sphere overestimates when it hugs the
    structure (geometric factor, reported as an estimate) and at low kr (near field); the cylinder
    overestimates by ~0.5 to 4.7 dB because the modes cross it obliquely (only the radial intensity carries
    power). Exact values need the active intensity 1/2 Re{p v_r*} (e.g. cylinders at r +- 5 cm).

    The FFT of the Fortran code (compute_rfft) is a Hann window scaled for amplitude, so |p|^2/2 per bin is the
    power of a tone on that bin while sums over bins are ENBW times too large: densities and band/total
    powers are divided by enbw, the per-bin tone level is not. By default enbw is the one of the file
    (case["ENBW"], see readers.solver.read_enbw): 1.5 for forces from compute_rfft, 1 for forces without a window.

    Also computed:
    - Directivity index per 1/3 octave band on the sphere, DI = 10 log10(max I / mean I).
    - Multipole fit on the sphere around its centre, p = A h0(kr) + (D . n) h1(kr): monopole and dipole parts
      of the power (the dipole one with its exact near-field factor 1 + 1/(kr)^2), the fraction of |p|^2 they
      explain (how well the emitter reduces to a single point) and the dipole axis per band.
    - Propagation loss PL = SL - RL per 1/3 octave band at the spectrum observers (and their energetic mean),
      to turn a measured band level into a source level: SL_measured = L_measured + PL.

    Returns
    -------
    metrics : dict
        Narrowband arrays on "Freqs" (W_*_bin [W], SL_*_bin [dB re 1 uPa^2 m^2], SL_*_density
        [dB re 1 uPa^2 m^2/Hz], kr, ratio_insitu_free [dB], multipole fractions), 1/3 octave arrays on "fc"
        (SL_*_band, DI_band, PL_band, dipole axis, ...) and broadband totals (SL_*_total, W_*_total).
    """

    has_sphere, has_cyl = case.get("has_sphere", False), case.get("has_cylinder", False)
    if not (has_sphere or has_cyl):
        print("compute_point_source(): no sphere nor cylinder data, skipping")
        return None

    rho_wat = case.get("rho_wat", rho_water)
    c_wat   = case.get("c_wat", c_water)
    p_ref   = case["p_ref"]
    tag     = emitter_label(case)
    Freqs   = case["Freqs"]
    enbw    = case["ENBW"] if enbw is None else enbw
    df      = case["bin_width"]           # [Hz] Bin width, from the full grid
    if df is None:
        print(f"compute_point_source(): its densities and 1/3 octave bands need uniformly spaced frequencies "
              f"(a bin width) and those of {case['Name']} are not: not applicable, skipping")
        return None

    # --- 0. Frequency band --- #
    mask = np.ones(len(Freqs), dtype=bool)
    if filter_under is not None:
        mask &= (Freqs >= filter_under)
    if filter_over is not None:
        mask &= (Freqs <= filter_over)
    if not np.any(mask):
        print("compute_point_source(): no frequencies match the specified filter bounds!")
        return None

    f  = Freqs[mask]
    k  = 2.0 * np.pi * f / c_wat
    fc = au.generate_iso_third_octave(f.min(), f.max())[1]

    def to_SL(W):                       # [W] -> [dB re 1 uPa^2 m^2]: free-field monopole of the same power at 1 m
        with np.errstate(divide='ignore'):
            return 10.0 * np.log10(rho_wat * c_wat * W / (4.0 * np.pi * p_ref**2))

    def band_sum(values):               # Sum over every 1/3 octave band (same bands and edges as to_third_octave)
        return au.to_third_octave(f, values, fmin=f.min(), fmax=f.max(), domain='linear', linear_mode='power')[1]

    bw      = fc * (2**(1/6) - 2**(-1/6))     # [Hz] 1/3 octave band widths
    metrics = {"Freqs": f, "fc": fc, "df": df, "enbw": enbw}

    # --- 1. Free field: sphere --- #
    if has_sphere:
        p_s    = case["P_sphere"][mask].reshape(len(f), -1)              # (Nf, Nobs), peak amplitudes
        obs    = case["Obs_sphere"].reshape(-1, 3)
        center = np.asarray(case["Center_sphere"], dtype=float)
        r_s    = case["R_sphere"]
        dA_s   = case["dA_sphere"]
        valid  = np.all(np.isfinite(p_s), axis=0)
        p_s    = p_s[:, valid]

        # Radiated power of every bin [W], tone convention (not divided by the ENBW)
        W_free_bin = np.sum(np.abs(p_s)**2, axis=1) / (2.0 * rho_wat * c_wat) * dA_s

        # Validity: near field (kr) and how tightly the sphere hugs the wet nodes (geometric bias estimate: monopole at a)
        nodes = case["nodes"]
        nodes = nodes[nodes[:, 2] < 0.0]
        a_max = np.max(np.linalg.norm(nodes - center, axis=1)) if len(nodes) else np.nan
        ratio = a_max / r_s
        geo_bias_dB = 10*np.log10((1/(2*ratio)) * np.log((1+ratio)/(1-ratio))) if 0 < ratio < 1 else np.nan

        # Multipole fit p = c0 + c . n (angular basis, constant over frequency since r is fixed)
        n_hat  = (obs[valid] - center) / r_s
        B      = np.column_stack([np.ones(len(n_hat)), n_hat])            # (Nobs, 4)
        coef   = np.linalg.lstsq(B, p_s.T, rcond=None)[0]                 # (4, Nf)
        G      = B.T @ B
        p2_fit = np.real(np.einsum('ik,ij,jk->k', coef.conj(), G, coef))  # ||B c||^2 per bin
        p2_all = np.sum(np.abs(p_s)**2, axis=1)
        kr     = k * r_s
        W_mono_bin = 4.0*np.pi * np.abs(coef[0])**2 * r_s**2 / (2.0*rho_wat*c_wat)
        W_dip_bin  = (4.0*np.pi/3.0) * np.sum(np.abs(coef[1:])**2, axis=0) * r_s**2 / (2.0*rho_wat*c_wat) / (1.0 + 1.0/kr**2)

        # Directivity per band (band-summed intensity pattern on the sphere)
        I_band  = au.to_third_octave(f, np.abs(p_s)**2, fmin=f.min(), fmax=f.max(), domain='linear', linear_mode='power')[1]
        DI_band = 10.0 * np.log10(np.max(I_band, axis=1) / np.mean(I_band, axis=1))

        # Dipole axis per band: principal direction of Re(sum D D^H), D = c / h1 (only its direction matters)
        axis_elev, axis_azim = np.full(len(fc), np.nan), np.full(len(fc), np.nan)
        for b, f_c in enumerate(fc):
            in_b = (f >= f_c / 2**(1/6)) & (f < f_c * 2**(1/6))
            if np.any(in_b):
                M = np.real(coef[1:, in_b] @ coef[1:, in_b].conj().T)
                u = np.linalg.eigh(M)[1][:, -1]
                u = u if u[2] >= 0 else -u
                axis_elev[b] = np.degrees(np.arcsin(np.clip(u[2], -1, 1)))     # [deg] from the horizontal plane
                axis_azim[b] = np.degrees(np.arctan2(u[1], u[0])) % 180.0       # [deg] from +x (axis, so modulo 180)

        W_free_band = band_sum(W_free_bin) / enbw
        metrics.update({
            "W_free_bin": W_free_bin, "SL_free_bin": to_SL(W_free_bin), "SL_free_density": to_SL(W_free_bin / (enbw * df)),
            "W_free_band": W_free_band, "SL_free_band": to_SL(W_free_band), "W_free_total": np.sum(W_free_bin) / enbw,
            "SL_free_total": to_SL(np.sum(W_free_bin) / enbw),
            "kr": kr, "R_sphere": r_s, "a_over_r": ratio, "geo_bias_dB": geo_bias_dB,
            "multipole_explained": p2_fit / p2_all, "W_mono_bin": W_mono_bin, "W_dip_bin": W_dip_bin,
            "multipole_explained_band": band_sum(p2_fit) / band_sum(p2_all),
            "mono_fraction_band": band_sum(W_mono_bin) / band_sum(W_free_bin),
            "dip_fraction_band" : band_sum(W_dip_bin) / band_sum(W_free_bin),
            "SL_multipole_band" : to_SL(band_sum(W_mono_bin + W_dip_bin) / enbw),
            "DI_band": DI_band, "dipole_elevation_band": axis_elev, "dipole_azimuth_band": axis_azim,
            "Excluded_sphere_fraction": 1.0 - np.count_nonzero(valid) / valid.size,
        })

    # --- 2. In situ: cylinder (waveguide) --- #
    if has_cyl:
        p_c   = case["P_cylinder"][mask]                                   # (Nf, Ntheta, Nz)
        r_c   = case["R_cylinder"]
        z_c   = case["Z_cylinder"]
        theta = np.deg2rad(case["Theta_deg_cylinder"])
        I_c   = np.abs(p_c)**2 / (2.0 * rho_wat * c_wat)
        I_c   = np.where(np.isfinite(I_c), I_c, 0.0)                       # Masked observers left out

        # Closed azimuthal integral (last point = first point), then trapezoid along z, for every bin
        theta_closed = np.append(theta, theta[0] + 2*np.pi)
        I_closed     = np.concatenate([I_c, I_c[:, :1, :]], axis=1)
        W_insitu_bin = np.abs(trapezoid(trapezoid(I_closed * r_c, theta_closed, axis=1), z_c, axis=1))

        W_insitu_band = band_sum(W_insitu_bin) / enbw
        metrics.update({
            "W_insitu_bin": W_insitu_bin, "SL_insitu_bin": to_SL(W_insitu_bin),
            "SL_insitu_density": to_SL(W_insitu_bin / (enbw * df)), "W_insitu_band": W_insitu_band,
            "SL_insitu_band": to_SL(W_insitu_band), "W_insitu_total": np.sum(W_insitu_bin) / enbw,
            "SL_insitu_total": to_SL(np.sum(W_insitu_bin) / enbw), "R_cylinder": r_c,
        })

    if has_sphere and has_cyl:
        with np.errstate(divide='ignore', invalid='ignore'):
            metrics["ratio_insitu_free"]      = 10.0 * np.log10(metrics["W_insitu_bin"] / metrics["W_free_bin"])
            metrics["ratio_insitu_free_band"] = 10.0 * np.log10(metrics["W_insitu_band"] / metrics["W_free_band"])

    # --- 3. Propagation loss at the spectrum observers: PL = SL - RL per band --- #
    if case.get("has_spectrums", False):
        SL_src  = metrics["SL_free_band"] if has_sphere else metrics["SL_insitu_band"]
        p_o     = case["P_spectrums"][mask]                                # (Nf, Nobs)
        RL_band = 10.0 * np.log10(band_sum(np.abs(p_o)**2 / 2.0) / enbw / p_ref**2)    # Band SPL [dB re 1 uPa]
        obs_o   = case["Obs_spectrums"]
        with np.errstate(invalid='ignore'):
            RL_mean = 10.0 * np.log10(np.nanmean(10.0**(RL_band / 10.0), axis=1))
        metrics.update({
            "PL_band": SL_src[:, None] - RL_band, "PL_mean_band": SL_src - RL_mean, "RL_band": RL_band,
            "PL_observers": obs_o, "PL_source": "free" if has_sphere else "insitu",
        })

    # --- 4. Print --- #
    if verbose:
        lo = f"{f.min():.1f}"; hi = f"{f.max():.1f}"
        print("\n" + "=" * 96)
        print(f"  EQUIVALENT POINT SOURCE [{tag}]  |  {lo} - {hi} Hz  |  df = {df:.3f} Hz, ENBW = {enbw:g} bins")
        print("=" * 96)
        label_w = 46
        if has_sphere:
            print(f"  {'SL free field, broadband':<{label_w}}: {metrics['SL_free_total']:.1f} dB re 1 µPa² m²  "
                  f"(W = {metrics['W_free_total']:.3e} W)")
            print(f"  {'Sphere r / farthest wet node (a/r)':<{label_w}}: {r_s:.1f} m / {a_max:.1f} m ({ratio:.2f})")
            if np.isfinite(geo_bias_dB) and geo_bias_dB > 0.5:
                print(f"  WARNING sphere hugs the structure: |p|² power ~+{geo_bias_dB:.1f} dB too high (geometric estimate) "
                      f"(r_sphere ≳ 3·a = {3*a_max:.0f} m keeps it below ~0.2 dB)")
            if np.any(kr < 3):
                print(f"  WARNING near field (kr < 3) below {f[kr < 3].max():.1f} Hz: |p|² power too high "
                      f"(dipole: up to +{10*np.log10(1 + 1/kr.min()**2):.1f} dB at {f.min():.1f} Hz)")
        if has_cyl:
            print(f"  {'SL in situ (waveguide), broadband':<{label_w}}: {metrics['SL_insitu_total']:.1f} dB re 1 µPa² m²  "
                  f"(W = {metrics['W_insitu_total']:.3e} W, r_cyl = {r_c:.0f} m)")
            print(f"  WARNING cylinder |p|² power ignores the oblique modal propagation: ~+0.5 to +4.7 dB too high "
                  f"(exact: 1/2 Re{{p v_r*}})")
        if has_sphere and has_cyl:
            print(f"  {'In situ / free power, broadband':<{label_w}}: "
                  f"{10*np.log10(metrics['W_insitu_total'] / metrics['W_free_total']):+.1f} dB")

        head = f"\n  {'fc [Hz]':>8} {'SL free':>8} {'SL situ':>8} {'situ-free':>9} {'DI':>5} {'mono%':>6} {'dip%':>6} " \
               f"{'expl%':>6} {'dip.elev':>8} {'PL mean':>8}"
        print(head + "\n  " + " " * 8 + f" {'[dB re 1 µPa² m²]':>17} {'[dB]':>9} {'[dB]':>5}" + " " * 21 + f"{'[deg]':>8} {'[dB]':>8}")
        for b, f_c in enumerate(fc):
            g = lambda key, fmt: format(metrics[key][b], fmt) if key in metrics and np.isfinite(metrics[key][b]) else "-"
            pct = lambda key: format(100*metrics[key][b], ".0f") if key in metrics and np.isfinite(metrics[key][b]) else "-"
            print(f"  {f_c:8.1f} {g('SL_free_band', '.1f'):>8} {g('SL_insitu_band', '.1f'):>8} "
                  f"{g('ratio_insitu_free_band', '+.1f'):>9} {g('DI_band', '.1f'):>5} {pct('mono_fraction_band'):>6} "
                  f"{pct('dip_fraction_band'):>6} {pct('multipole_explained_band'):>6} {g('dipole_elevation_band', '.0f'):>8} "
                  f"{g('PL_mean_band', '.1f'):>8}")
        print("  SL: power-equivalent free-field monopole, band levels (power in the band). mono%/dip%: share of the sphere "
              "power,\n  expl%: share of |p|² explained by monopole + dipole, PL: SL - band SPL at the spectrum observers "
              "(energetic mean).")

    # --- 5. Plots --- #
    if do_plot:
        color = case_color(case)
        smooth = lambda y: au.smooth_fractional_octave(f, y, smooth_fraction, domain='db') if smooth_fraction else y
        note   = f" (1/{smooth_fraction} oct. smoothed)" if smooth_fraction else ""

        # 5.1 Source level spectra: narrowband densities and band-averaged densities (same units)
        fig, ax = plt.subplots(figsize=figsize)
        if has_sphere:
            ax.semilogx(f, smooth(metrics["SL_free_density"]), '-', color=color, label=f"Free field{note}")
            ax.semilogx(fc, metrics["SL_free_band"] - 10*np.log10(bw), 'o', color=color, label="Free field, 1/3 oct.")
            if np.any(kr < 3):
                ax.axvspan(f.min(), f[kr < 3].max(), color='0.85', zorder=0, label=r"Near field ($kr < 3$)")
        if has_cyl:
            ax.semilogx(f, smooth(metrics["SL_insitu_density"]), '--', color='k', label=f"In situ{note}")
            ax.semilogx(fc, metrics["SL_insitu_band"] - 10*np.log10(bw), 's', color='k', mfc='none', label="In situ, 1/3 oct.")
        totals = r" $|$ ".join(f"{name} {metrics[key]:.1f} dB" for name, key in (("free", "SL_free_total"), ("in situ", "SL_insitu_total"))
                               if key in metrics)
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel(r"SL [dB re 1 $\mu$Pa$^2$ m$^2$/Hz]")
        ax.set_title(f"{tag} - Equivalent point source\nBroadband SL: {totals} re 1 " + r"$\mu$Pa$^2$ m$^2$")
        ax.grid(True, which='both', ls='--', alpha=0.4)
        ax.legend(fontsize=11)
        fig.tight_layout()
        set_figure_name(fig, "point_source", "SL")

        # 5.2 Radiation diagnostics per 1/3 octave band
        if has_sphere:
            n_rows = 3 if has_cyl else 2
            fig2, axs = plt.subplots(n_rows, 1, figsize=(figsize[0], 2.6 * n_rows), sharex=True)
            row = 0
            if has_cyl:
                axs[row].semilogx(fc, metrics["ratio_insitu_free_band"], '-o', color='k')
                axs[row].axhline(0.0, color='0.5', lw=0.8)
                axs[row].set_ylabel("In situ − free\n[dB]")
                row += 1
            axs[row].semilogx(fc, DI_band, '-o', color=color)
            axs[row].set_ylabel("DI [dB]")
            row += 1
            axs[row].semilogx(fc, 100*metrics["dip_fraction_band"], '-o', color=color, label="Dipole power")
            axs[row].semilogx(fc, 100*metrics["mono_fraction_band"], '--s', color=color, label="Monopole power")
            axs[row].semilogx(fc, 100*metrics["multipole_explained_band"], ':^', color='k', label=r"$|p|^2$ explained")
            axs[row].set_ylabel(r"Share [\%]" if plt.rcParams["text.usetex"] else "Share [%]")    # % starts a LaTeX comment
            axs[row].legend(fontsize=10)
            axs[-1].set_xlabel("1/3 octave band centre frequency [Hz]")
            for a in axs:
                a.grid(True, which='both', ls='--', alpha=0.4)
                if np.any(kr < 3):
                    a.axvspan(fc.min() / 2**(1/6), f[kr < 3].max(), color='0.85', zorder=0)
            axs[0].set_title(f"{tag} - Radiation diagnostics (sphere r = {r_s:.0f} m)")
            fig2.tight_layout()
            set_figure_name(fig2, "point_source", "diagnostics")

        # 5.3 Propagation loss at the spectrum observers
        if "PL_band" in metrics:
            fig3, ax3 = plt.subplots(figsize=figsize)
            shades = get_color_shades(color, metrics["PL_band"].shape[1])
            styles = get_line_styles()
            for i, o in enumerate(metrics["PL_observers"]):
                ax3.semilogx(fc, metrics["PL_band"][:, i], styles[i % len(styles)], marker='o', color=shades[i],
                             label=f"({o[0]:.1f}, {o[1]:.1f}, {o[2]:.1f}) m")
            ax3.semilogx(fc, metrics["PL_mean_band"], '-', color='k', lw=2, label="Energetic mean")
            ax3.set_xlabel("1/3 octave band centre frequency [Hz]")
            ax3.set_ylabel(r"PL [dB re 1 m$^2$]")
            ax3.set_title(f"{tag} - Propagation loss, PL = SL$_{{{metrics['PL_source']}}}$ − SPL")
            ax3.grid(True, which='both', ls='--', alpha=0.4)
            ax3.legend(fontsize=11)
            fig3.tight_layout()
            set_figure_name(fig3, "point_source", "PL")

    return metrics
