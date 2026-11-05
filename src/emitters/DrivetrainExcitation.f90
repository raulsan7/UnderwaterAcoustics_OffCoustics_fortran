! ===============================================================================
! MODULE: DrivetrainExcitation
! ===============================================================================
! Defines the drivetrain excitation of the turbine types: the excitation
! spectra of the DTU 10 MW (drivetrain10MW_excitation_spectrum), SAITEC
! 2 MW (drivetrain2MW_excitation_spectrum) and MySE 3.0-135
! (drivetrain3MW_excitation_spectrum, with the reader of its force files
! read_drivetrain3MW_force_file) drivetrains, the synthetic time series
! built from them (generate_timeseries_banded_sines), and the mask of the
! frequencies worth keeping (filter_non_useful_freqs).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE DrivetrainExcitation

USE Kinds, ONLY: I32, WP, PI
USE Sorting, ONLY: sort_real_c
USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_nan

IMPLICIT NONE

PRIVATE
PUBLIC :: generate_timeseries_banded_sines, filter_non_useful_freqs, drivetrain10MW_excitation_spectrum, &
          drivetrain2MW_excitation_spectrum, drivetrain3MW_excitation_spectrum, read_drivetrain3MW_force_file


CONTAINS


! ===============================================================================
! SYNTHETIC EXCITATION AND FREQUENCY FILTERING
! ===============================================================================
SUBROUTINE generate_timeseries_banded_sines(peaks, keys, t, zeta, nfreq, seed, &
                                            used_freqs, fcut, a, freqs_out)
    ! -------------------------------------------------------------------------
    ! Reconstructs a synthetic time-series signal from spectral peak
    ! distributions. Optimized for HPC: aggregates all frequency components
    ! first, then performs a single OpenMP+SIMD pass over the time domain to
    ! maximize cache locality. The phases are reproducible for a given seed,
    ! and nfreq = 1 gives a single line at f0 for every banded peak.
    ! -------------------------------------------------------------------------

    REAL(WP)         , INTENT(IN)                         :: peaks(:,:)      ! [-] Peaks (Npeaks, 2): freq [Hz], RMS
    CHARACTER(len=20), INTENT(IN)                         :: keys(:)         ! [-] Peak labels, 'mesh' = gear mesh
    REAL(WP)         , INTENT(IN)                         :: t(:)            ! [s] Time vector (nt)
    REAL(WP)         , INTENT(IN), OPTIONAL               :: zeta            ! [-] Bandwidth ratio (sigma = zeta * f0)
    INTEGER(I32)     , INTENT(IN), OPTIONAL               :: nfreq           ! [-] Sine components per banded peak
    INTEGER(I32)     , INTENT(IN), OPTIONAL               :: seed            ! [-] Random seed
    LOGICAL          , INTENT(IN), OPTIONAL               :: used_freqs      ! [-] Whether to return freqs_out
    REAL(WP)         , INTENT(IN), OPTIONAL               :: fcut            ! [Hz] Peaks above fcut are always banded
    REAL(WP)         , INTENT(OUT), ALLOCATABLE           :: a(:)            ! [-] Synthesized signal (nt)
    REAL(WP)         , INTENT(OUT), OPTIONAL, ALLOCATABLE :: freqs_out(:)    ! [Hz] Sorted unique frequencies used

    ! Local variables
    INTEGER(I32)              :: i, j, seed_size, nt, n_peaks, n_active, nfreq_, seed_
    INTEGER(I32)              :: max_freq_list, n_kept
    REAL(WP)                  :: f0, Arms, sigma, phase, rnd, zeta_, fcut_, sqrt2
    REAL(WP)    , ALLOCATABLE :: fk(:), g(:), Arms_k(:), phases(:)

    ! Global accumulation arrays for all sine wave components
    REAL(WP)         , ALLOCATABLE :: all_omega(:), all_amp(:), all_phase(:), freq_list(:)
    INTEGER(I32)     , ALLOCATABLE :: seed_arr(:)
    LOGICAL                        :: used_freqs_
    CHARACTER(len=20)              :: key_label

    ! Time loop variables
    REAL(WP) :: current_t, current_sum


    ! Defaults
    sqrt2       = sqrt(2.0_WP)
    zeta_       = 0.02_WP ; if (present(zeta))       zeta_       = zeta
    nfreq_      = 50_I32  ; if (present(nfreq))      nfreq_      = nfreq
    seed_       = 42_I32  ; if (present(seed))       seed_       = seed
    fcut_       = 10.0_WP ; if (present(fcut))       fcut_       = fcut
    used_freqs_ = .false. ; if (present(used_freqs)) used_freqs_ = used_freqs

    if (nfreq_ < 1_I32) error stop "generate_timeseries_banded_sines: nfreq must be >= 1"

    nt = size(t)
    ALLOCATE(a(nt))
    n_peaks = size(peaks, 1)

    if (size(keys) < n_peaks) error stop "generate_timeseries_banded_sines: keys shorter than peaks"
    if (size(peaks, 2) < 2) error stop "generate_timeseries_banded_sines: peaks must have at least 2 columns"

    ! Pre-allocate global lists based on the maximum possible components
    max_freq_list = max(1_I32, n_peaks * nfreq_)
    ALLOCATE(freq_list(max_freq_list))
    ALLOCATE(all_omega(max_freq_list), all_amp(max_freq_list), all_phase(max_freq_list))
    n_active = 0_I32

    ! Initialize random seed
    CALL random_seed(size=seed_size)
    ALLOCATE(seed_arr(seed_size))
    do i = 1, seed_size
        seed_arr(i) = seed_ + 37_I32 * (i - 1_I32)
    end do
    CALL random_seed(put=seed_arr)

    ! --- 1. Aggregate all frequency components (no time-domain computation here) --- !
    do i = 1, n_peaks
        f0        = peaks(i, 1)
        Arms      = peaks(i, 2)
        key_label = trim(keys(i))

        ! Gear mesh or high frequency -> banded Gaussian
        if (index(key_label, "mesh") > 0 .or. f0 > fcut_) then
            sigma = zeta_ * f0
            if (sigma <= 0.0_WP) cycle

            ALLOCATE(fk(nfreq_), g(nfreq_), Arms_k(nfreq_), phases(nfreq_))

            if (nfreq_ == 1_I32) then
                fk(1) = f0
            else
                do j = 1, nfreq_
                    fk(j) = f0 - 4.0_WP * sigma + (j - 1_I32) * (8.0_WP * sigma) / real(nfreq_ - 1_I32, WP)
                end do
            end if

            g      = exp(-0.5_WP * ((fk - f0) / sigma)**2)
            g      = g / sqrt(sum(g**2))
            Arms_k = Arms * g

            do j = 1, nfreq_
                CALL random_number(rnd)
                phases(j) = 2.0_WP * PI * rnd

                ! Store properties globally
                n_active            = n_active + 1_I32
                freq_list(n_active) = fk(j)
                all_omega(n_active) = 2.0_WP * PI * fk(j)
                all_amp(n_active)   = Arms_k(j)
                all_phase(n_active) = phases(j)
            end do

            DEALLOCATE(fk, g, Arms_k, phases)

        else
            ! Pure tone (shaft)
            CALL random_number(rnd)
            phase = 2.0_WP * PI * rnd

            ! Store properties globally
            n_active            = n_active + 1_I32
            freq_list(n_active) = f0
            all_omega(n_active) = 2.0_WP * PI * f0
            all_amp(n_active)   = sqrt2 * Arms
            all_phase(n_active) = phase
        end if
    end do

    ! --- 2. Time-domain synthesis --- !
    ! Time is the outer loop: the main memory array 'a' is written exactly once per time step
    !$OMP PARALLEL DO PRIVATE(i, j, current_t, current_sum)
    do i = 1, nt
        current_t   = t(i)
        current_sum = 0.0_WP

        ! Inner loop over frequencies: sum into a CPU register via SIMD
        !$OMP SIMD REDUCTION(+:current_sum)
        do j = 1, n_active
            current_sum = current_sum + all_amp(j) * sin(all_omega(j) * current_t + all_phase(j))
        end do

        a(i) = current_sum
    end do
    !$OMP END PARALLEL DO

    ! --- 3. Frequency sorting and cleanup --- !
    if (used_freqs_) then
        if (n_active > 0_I32) then
            CALL sort_real_c(freq_list(1:n_active))
            ALLOCATE(freqs_out(n_active))
            n_kept = 0_I32
            do i = 1, n_active
                if (i == 1_I32) then
                    n_kept            = n_kept + 1_I32
                    freqs_out(n_kept) = freq_list(i)
                else if (abs(freq_list(i) - freq_list(i-1)) > 1.0e-6_WP) then
                    n_kept            = n_kept + 1_I32
                    freqs_out(n_kept) = freq_list(i)
                end if
            end do
            if (n_kept < n_active) freqs_out = freqs_out(1:n_kept)
        else
            ALLOCATE(freqs_out(0))
        end if
    else
        if (present(freqs_out)) ALLOCATE(freqs_out(0))
    end if

    DEALLOCATE(freq_list, seed_arr, all_omega, all_amp, all_phase)

END SUBROUTINE generate_timeseries_banded_sines


FUNCTION filter_non_useful_freqs(freqs, maintain_freqs, freqs_over, freqs_under, width_bins) RESULT(mask)
    ! -------------------------------------------------------------------------
    ! Builds a logical mask selecting which frequencies to keep. Frequencies
    ! outside the [freqs_over, freqs_under] zone are always kept; frequencies
    ! inside that zone are kept only if they lie within one frequency bin (df,
    ! or width_bins bins when given) of one of the maintain_freqs values, so
    ! that physically relevant peaks are preserved while the rest of the band
    ! is discarded.
    ! -------------------------------------------------------------------------

    REAL(WP), INTENT(IN)           :: freqs(:)             ! [Hz] Frequency vector (Nfreqs)
    REAL(WP), INTENT(IN)           :: maintain_freqs(:)    ! [Hz] Frequencies to preserve inside the filter zone
    REAL(WP), INTENT(IN), OPTIONAL :: freqs_over           ! [Hz] Lower limit of the filter zone
    REAL(WP), INTENT(IN), OPTIONAL :: freqs_under          ! [Hz] Upper limit of the filter zone
    REAL(WP), INTENT(IN), OPTIONAL :: width_bins           ! [-] Half-width of the kept band around each frequency
    LOGICAL , ALLOCATABLE          :: mask(:)              ! [-] Frequency mask (Nfreqs), .true. to keep

    ! Local variables
    INTEGER(I32)              :: n, j
    REAL(WP)                  :: df, half_width, freqs_over_, freqs_under_
    LOGICAL     , ALLOCATABLE :: filter_zone(:)


    n = size(freqs)
    ALLOCATE(mask(n))
    ALLOCATE(filter_zone(n))

    ! Determine frequency step
    if (n > 1_I32) then
        df = freqs(2) - freqs(1)
        ! Optionally check for equal spacing (omitted for performance)
    else
        df = 0.0_WP
    end if

    half_width = df ; if (present(width_bins)) half_width = width_bins * df

    ! Define filter zone (frequencies to consider for removal)
    filter_zone  = .true.
    freqs_over_  = minval(freqs)
    freqs_under_ = maxval(freqs)
    if (present(freqs_over))  freqs_over_  = freqs_over
    if (present(freqs_under)) freqs_under_ = freqs_under

    if (present(freqs_over))  filter_zone = filter_zone .and. (freqs >= freqs_over_)
    if (present(freqs_under)) filter_zone = filter_zone .and. (freqs <= freqs_under_)

    ! Outside filter zone: keep all frequencies
    mask = .not. filter_zone

    ! Inside filter zone: keep only if within half_width of any maintain_freqs
    do j = 1, size(maintain_freqs)
        WHERE (filter_zone .and. abs(freqs - maintain_freqs(j)) <= half_width)
            mask = .true.
        END WHERE
    end do

    DEALLOCATE(filter_zone)

END FUNCTION filter_non_useful_freqs


! ===============================================================================
! DTU 10 MW DRIVETRAIN EXCITATIONS
! ===============================================================================
SUBROUTINE drivetrain10MW_excitation_spectrum(rpm, damping, p_shaft, alpha_mesh, freqs_amp, keys)
    ! -------------------------------------------------------------------------
    ! Constructs the physically-consistent excitation spectrum (frequencies +
    ! amplitudes) for the DTU 10 MW medium-speed drivetrain, combining shaft
    ! and gear-mesh harmonics, RPM-dependent interpolation and a torsional
    ! modal transfer function, normalized to unit RMS.
    ! -------------------------------------------------------------------------

    REAL(WP)         , INTENT(IN), OPTIONAL     :: rpm               ! [rpm] Rotor speed
    REAL(WP)         , INTENT(IN), OPTIONAL     :: damping           ! [-] Modal damping ratio of the torsional modes
    REAL(WP)         , INTENT(IN), OPTIONAL     :: p_shaft           ! [-] Decay exponent of shaft harmonic amplitudes
    REAL(WP)         , INTENT(IN), OPTIONAL     :: alpha_mesh        ! [-] Decay rate of gear-mesh harmonic amplitudes
    REAL(WP)         , INTENT(OUT), ALLOCATABLE :: freqs_amp(:,:)    ! [-] Peaks (N, 2): frequency [Hz], amplitude
    CHARACTER(len=20), INTENT(OUT), ALLOCATABLE :: keys(:)           ! [-] Peak labels (N)

    ! Local variables
    REAL(WP), PARAMETER :: min_rpm = 6.0_WP, max_rpm = 9.6_WP
    REAL(WP)            :: factor, rpm_clipped, modal_sum

    ! Torsional eigenfrequencies [Hz]
    REAL(WP), PARAMETER :: torsional_modes(14) = [ &
        3.942_WP, 16.683_WP, 33.811_WP, 61.736_WP, 72.015_WP, &
        108.447_WP, 175.792_WP, 184.351_WP, 214.403_WP, &
        245.476_WP, 274.443_WP, 336.874_WP, 346.436_WP, 391.730_WP]

    ! Shaft frequencies at min and max rpm
    REAL(WP), PARAMETER :: shaft_min(13) = [ &
        0.1_WP, 0.2_WP, 0.3_WP, 0.6_WP, &
        0.322_WP, 0.644_WP, 0.966_WP, &
        1.383_WP, 2.766_WP, 4.149_WP, &
        5.017_WP, 10.034_WP, 15.051_WP]
    REAL(WP), PARAMETER :: shaft_max(13) = [ &
        0.160_WP, 0.320_WP, 0.480_WP, 0.640_WP, &
        0.515_WP, 1.030_WP, 1.545_WP, &
        2.213_WP, 4.426_WP, 6.639_WP, &
        8.028_WP, 16.056_WP, 24.024_WP]
    CHARACTER(len=6), PARAMETER :: shaft_keys(13) = [ &
        "lss_1p", "lss_2p", "lss_3p", "lss_6p", &
        "is1_1p", "is1_2p", "is1_3p", &
        "is2_1p", "is2_2p", "is2_3p", &
        "hss_1p", "hss_2p", "hss_3p"]

    ! Gear mesh frequencies at min and max rpm
    REAL(WP), PARAMETER :: mesh_min(9) = [ &
        10.30_WP, 20.6_WP, 30.9_WP, &
        48.952_WP, 97.904_WP, 146.856_WP, &
        93.423_WP, 186.846_WP, 280.269_WP]
    REAL(WP), PARAMETER :: mesh_max(9) = [ &
        16.48_WP, 32.96_WP, 49.44_WP, &
        78.3_WP, 156.6_WP, 234.9_WP, &
        149.492_WP, 298.984_WP, 448.476_WP]
    CHARACTER(len=11), PARAMETER :: mesh_keys(9) = [ &
        "lss_mesh_1p", "lss_mesh_2p", "lss_mesh_3p", &
        "ims_mesh_1p", "ims_mesh_2p", "ims_mesh_3p", &
        "hss_mesh_1p", "hss_mesh_2p", "hss_mesh_3p"]

    ! Base amplitudes (hardcoded)
    REAL(WP), PARAMETER :: A_base_shaft(13) = [ &
        0.6_WP, 0.3_WP, 1.0_WP, 0.2_WP, &
        0.4_WP, 0.2_WP, 0.6_WP, &
        0.3_WP, 0.15_WP, 0.4_WP, &
        0.2_WP, 0.1_WP, 0.3_WP]
    REAL(WP), PARAMETER :: A_base_mesh(9) = [ &
        0.8_WP, 0.8_WP, 0.8_WP, &
        0.5_WP, 0.5_WP, 0.5_WP, &
        0.3_WP, 0.3_WP, 0.3_WP]

    INTEGER(I32)                   :: n_total, i, m, pos_p, pos_underscore, n
    REAL(WP)         , ALLOCATABLE :: freqs(:), amps_source(:)
    REAL(WP)                       :: H, rms, fn, rpm_, damping_, alpha_mesh_, p_shaft_
    CHARACTER(len=20), ALLOCATABLE :: keys_work(:)
    CHARACTER(len=20)              :: key_part


    ! Defaults
    rpm_        = 9.6_WP  ; if (present(rpm))        rpm_        = rpm
    damping_    = 0.02_WP ; if (present(damping))    damping_    = damping
    p_shaft_    = 1.0_WP  ; if (present(p_shaft))    p_shaft_    = p_shaft
    alpha_mesh_ = 0.8_WP  ; if (present(alpha_mesh)) alpha_mesh_ = alpha_mesh

    ! Clip rpm
    rpm_clipped = min(max(rpm_, min_rpm), max_rpm)
    factor      = (rpm_clipped - min_rpm) / (max_rpm - min_rpm)

    ! Interpolate frequencies
    n_total = size(shaft_min) + size(mesh_min)
    ALLOCATE(freqs(n_total), amps_source(n_total), keys_work(n_total))

    ! Shaft frequencies and amplitudes
    do i = 1, size(shaft_min)
        freqs(i)     = shaft_min(i) + factor * (shaft_max(i) - shaft_min(i))
        keys_work(i) = shaft_keys(i)

        pos_underscore = index(trim(shaft_keys(i)), "_", back=.true.)
        pos_p          = index(trim(shaft_keys(i)), "p", back=.true.)

        if (pos_underscore > 0_I32 .and. pos_p > pos_underscore) then
            key_part = trim(shaft_keys(i)(pos_underscore+1:pos_p-1))
            read(key_part, *) n
        else
            n = 1_I32
        end if

        amps_source(i) = A_base_shaft(i) * (1.0_WP / real(n, WP)**p_shaft_)
    end do

    ! Mesh frequencies and amplitudes
    do i = 1, size(mesh_min)
        pos_underscore = index(trim(mesh_keys(i)), "_", back=.true.)
        pos_p          = index(trim(mesh_keys(i)), "p", back=.true.)

        if (pos_underscore > 0_I32 .and. pos_p > pos_underscore) then
            key_part = trim(mesh_keys(i)(pos_underscore+1:pos_p-1))
            read(key_part, *) m
        else
            m = 1_I32
        end if

        freqs(size(shaft_min)+i)       = mesh_min(i) + factor * (mesh_max(i) - mesh_min(i))
        keys_work(size(shaft_min)+i)   = mesh_keys(i)
        amps_source(size(shaft_min)+i) = A_base_mesh(i) * exp(-alpha_mesh_ * (real(m, WP) - 1.0_WP))
    end do

    ! Modal transfer function H(f)
    do i = 1, n_total
        modal_sum = 0.0_WP
        do m = 1, size(torsional_modes)
            fn        = torsional_modes(m)
            modal_sum = modal_sum + 1.0_WP / ((1.0_WP - (freqs(i) / fn)**2)**2 &
                                              + (2.0_WP * damping_ * freqs(i) / fn)**2)
        end do
        H              = sqrt(modal_sum)
        amps_source(i) = amps_source(i) * H
    end do

    ! Normalize to unit RMS: sqrt(sum(A^2/2)) = 1
    rms         = sqrt(sum(amps_source**2) / 2.0_WP)
    amps_source = amps_source / rms

    ! Build output array
    ALLOCATE(freqs_amp(n_total, 2))
    freqs_amp(:, 1) = freqs
    freqs_amp(:, 2) = amps_source
    ALLOCATE(keys(n_total))
    keys = keys_work

    DEALLOCATE(freqs, amps_source, keys_work)

END SUBROUTINE drivetrain10MW_excitation_spectrum


! ===============================================================================
! SAITEC 2 MW DRIVETRAIN EXCITATION
! ===============================================================================
SUBROUTINE drivetrain2MW_excitation_spectrum(rpm, path_barRB, path_barGB, fmin, fmax, tol, Normalised, &
                                             freqs_amp, keys)
    ! -------------------------------------------------------------------------
    ! Constructs the excitation spectrum (frequencies + amplitudes) for the
    ! SAITEC 2 MW variable-speed drivetrain described in Escaler & Mebarki
    ! (2018), combining rotor-bearing and gearbox kinematic harmonics with
    ! externally-supplied bar-chart amplitude data, merging near-duplicate
    ! frequencies (keeping the larger amplitude) and optionally normalizing to
    ! unit RMS.
    ! -------------------------------------------------------------------------

    REAL(WP)         , INTENT(IN)               :: rpm               ! [rpm] Rotor speed
    CHARACTER(len=*) , INTENT(IN)               :: path_barRB        ! [-] Path to the rotor-bearing bar-chart CSV file
    CHARACTER(len=*) , INTENT(IN)               :: path_barGB        ! [-] Path to the gearbox bar-chart CSV file
    REAL(WP)         , INTENT(IN), OPTIONAL     :: fmin              ! [Hz] Lower frequency limit (inclusive)
    REAL(WP)         , INTENT(IN), OPTIONAL     :: fmax              ! [Hz] Upper frequency limit (exclusive)
    REAL(WP)         , INTENT(IN), OPTIONAL     :: tol               ! [Hz] Tolerance to merge close frequencies
    LOGICAL          , INTENT(IN), OPTIONAL     :: Normalised        ! [-] Whether to normalize to unit RMS
    REAL(WP)         , INTENT(OUT), ALLOCATABLE :: freqs_amp(:,:)    ! [-] Peaks (N, 2): frequency [Hz], amplitude
    CHARACTER(len=20), INTENT(OUT), ALLOCATABLE :: keys(:)           ! [-] Peak labels (N)

    ! Local variables
    ! Kinematic peak frequencies (rotor-bearing RB1/RB2, gearbox GB1/GB2/GB3)
    REAL(WP) :: RB1f(12), RB2f(11), GB1f(9), GB2f(9), GB3f(11)
    REAL(WP) :: f3, f_em
    REAL(WP) :: fmin_, fmax_, tol_
    LOGICAL  :: normalised_

    REAL(WP), ALLOCATABLE :: RB_data(:,:), GB_data(:,:)
    REAL(WP), ALLOCATABLE :: RB1_amp(:), RB2_amp(:), GB1_amp(:), GB2_amp(:), GB3_amp(:)

    REAL(WP)         , ALLOCATABLE :: freqs_raw(:), amps_raw(:)
    CHARACTER(len=20), ALLOCATABLE :: keys_raw(:)

    REAL(WP)         , ALLOCATABLE :: merged_freqs(:), merged_amps(:)
    CHARACTER(len=20), ALLOCATABLE :: merged_keys(:)
    LOGICAL          , ALLOCATABLE :: mask(:)

    REAL(WP)          :: rms, fk, ak
    CHARACTER(len=20) :: key_tmp
    INTEGER(I32)      :: n_total, i, j, kk, best_idx, n_merged, n_kept


    ! Defaults
    fmin_       = 0.0_WP       ; if (present(fmin))       fmin_       = fmin
    fmax_       = huge(1.0_WP) ; if (present(fmax))       fmax_       = fmax
    tol_        = 1.0e-3_WP    ; if (present(tol))        tol_        = tol
    normalised_ = .true.       ; if (present(Normalised)) normalised_ = Normalised

    CALL compute_freqs2MW(rpm, RB1f, RB2f, GB1f, GB2f, GB3f, f3)

    f_em = 84.0_WP * f3

    ! --- 1. Read bar-chart amplitude data (CSV columns, "Nan" = missing) ---
    CALL read_bar_csv(path_barRB, 2_I32, RB_data)
    CALL read_bar_csv(path_barGB, 3_I32, GB_data)

    RB1_amp = pack(RB_data(:, 1), .not. ieee_is_nan(RB_data(:, 1)))
    RB2_amp = pack(RB_data(:, 2), .not. ieee_is_nan(RB_data(:, 2)))
    GB1_amp = pack(GB_data(:, 1), .not. ieee_is_nan(GB_data(:, 1)))
    GB2_amp = pack(GB_data(:, 2), .not. ieee_is_nan(GB_data(:, 2)))
    GB3_amp = pack(GB_data(:, 3), .not. ieee_is_nan(GB_data(:, 3)))

    if (size(RB1_amp) /= size(RB1f) .or. size(RB2_amp) /= size(RB2f) .or. &
        size(GB1_amp) /= size(GB1f) .or. size(GB2_amp) /= size(GB2f) .or. &
        size(GB3_amp) /= size(GB3f)) then
        error stop "drivetrain2MW_excitation_spectrum: bar-chart CSV row counts do not match " // &
                   "the expected kinematic peak counts"
    end if

    ! --- 2. Assemble raw (freq, amp, key) lists ---
    ! (kinematic frequency order preserved; amplitude re-paired positionally
    ! from each family's own frequency-sorted bar-chart column)
    n_total = size(RB1f) + size(RB2f) + size(GB1f) + size(GB2f) + size(GB3f) + 1_I32
    ALLOCATE(freqs_raw(n_total), amps_raw(n_total), keys_raw(n_total))

    j = 0_I32
    CALL append_family(RB1f, sorted_amps_by_freq(RB1f, RB1_amp), 'RB1_f_', freqs_raw, amps_raw, keys_raw, j)
    CALL append_family(RB2f, sorted_amps_by_freq(RB2f, RB2_amp), 'RB2_f_', freqs_raw, amps_raw, keys_raw, j)
    CALL append_family(GB1f, sorted_amps_by_freq(GB1f, GB1_amp), 'GB1_f_', freqs_raw, amps_raw, keys_raw, j)
    CALL append_family(GB2f, sorted_amps_by_freq(GB2f, GB2_amp), 'GB2_f_', freqs_raw, amps_raw, keys_raw, j)
    CALL append_family(GB3f, sorted_amps_by_freq(GB3f, GB3_amp), 'GB3_f_', freqs_raw, amps_raw, keys_raw, j)

    ! Electromagnetic peak at 84xf3 (observed in generator)
    j            = j + 1_I32
    freqs_raw(j) = f_em
    amps_raw(j)  = 0.2_WP
    keys_raw(j)  = 'generator_84xf3'

    DEALLOCATE(RB_data, GB_data, RB1_amp, RB2_amp, GB1_amp, GB2_amp, GB3_amp)

    ! --- 3. Sort everything by frequency (carrying amps/keys along, insertion sort) ---
    do i = 2, n_total
        fk      = freqs_raw(i)
        ak      = amps_raw(i)
        key_tmp = keys_raw(i)
        kk      = i - 1_I32
        do while (kk >= 1_I32)
            if (freqs_raw(kk) <= fk) exit
            freqs_raw(kk+1) = freqs_raw(kk)
            amps_raw(kk+1)  = amps_raw(kk)
            keys_raw(kk+1)  = keys_raw(kk)
            kk              = kk - 1_I32
        end do
        freqs_raw(kk+1) = fk
        amps_raw(kk+1)  = ak
        keys_raw(kk+1)  = key_tmp
    end do

    ! --- 4. Merge near-duplicate frequencies (within tol), keep larger amplitude ---
    ALLOCATE(merged_freqs(n_total), merged_amps(n_total), merged_keys(n_total))
    n_merged = 0_I32
    i        = 1_I32
    do while (i <= n_total)
        kk = i
        do while (kk + 1_I32 <= n_total)
            if ((freqs_raw(kk+1) - freqs_raw(i)) < tol_) then
                kk = kk + 1_I32
            else
                exit
            end if
        end do
        best_idx               = i - 1_I32 + maxloc(amps_raw(i:kk), dim=1)
        n_merged               = n_merged + 1_I32
        merged_freqs(n_merged) = freqs_raw(best_idx)
        merged_amps(n_merged)  = amps_raw(best_idx)
        merged_keys(n_merged)  = keys_raw(best_idx)
        i                      = kk + 1_I32
    end do
    DEALLOCATE(freqs_raw, amps_raw, keys_raw)

    ! --- 5. Apply frequency range filter (fmin, fmax) ---
    ALLOCATE(mask(n_merged))
    mask   = (merged_freqs(1:n_merged) >= fmin_) .and. (merged_freqs(1:n_merged) < fmax_)
    n_kept = count(mask)

    ALLOCATE(freqs_amp(n_kept, 2), keys(n_kept))
    j = 0_I32
    do i = 1, n_merged
        if (mask(i)) then
            j               = j + 1_I32
            freqs_amp(j, 1) = merged_freqs(i)
            freqs_amp(j, 2) = merged_amps(i)
            keys(j)         = merged_keys(i)
        end if
    end do
    DEALLOCATE(merged_freqs, merged_amps, merged_keys, mask)

    ! --- 6. Normalize to unit RMS: sqrt(sum(A^2/2)) = 1 ---
    if (normalised_ .and. n_kept > 0_I32) then
        rms             = sqrt(sum(freqs_amp(:, 2)**2) / 2.0_WP)
        freqs_amp(:, 2) = freqs_amp(:, 2) / rms
    end if

END SUBROUTINE drivetrain2MW_excitation_spectrum


SUBROUTINE compute_freqs2MW(rpm, RB1f, RB2f, GB1f, GB2f, GB3f, f3_out)
    ! -------------------------------------------------------------------------
    ! Computes the SAITEC 2 MW drivetrain kinematic peak frequencies (rotor
    ! bearings RB1/RB2 and gearbox stages GB1/GB2/GB3) at the given rotor
    ! speed, following the planetary+parallel gear-train kinematics of the
    ! Escaler & Mebarki (2018) 2 MW drivetrain.
    ! -------------------------------------------------------------------------

    REAL(WP), INTENT(IN)  :: rpm         ! [rpm] Rotor speed
    REAL(WP), INTENT(OUT) :: RB1f(12)    ! [Hz] Rotor bearing 1 kinematic peak frequencies
    REAL(WP), INTENT(OUT) :: RB2f(11)    ! [Hz] Rotor bearing 2 kinematic peak frequencies
    REAL(WP), INTENT(OUT) :: GB1f(9)     ! [Hz] Gearbox stage 1 kinematic peak frequencies
    REAL(WP), INTENT(OUT) :: GB2f(9)     ! [Hz] Gearbox stage 2 kinematic peak frequencies
    REAL(WP), INTENT(OUT) :: GB3f(11)    ! [Hz] Gearbox stage 3 kinematic peak frequencies
    REAL(WP), INTENT(OUT) :: f3_out      ! [Hz] High-speed shaft rotation frequency

    ! Local variables
    REAL(WP), PARAMETER :: zs = 18.0_WP, zp = 34.0_WP, zr = 87.0_WP
    REAL(WP), PARAMETER :: z12 = 70.0_WP, z21 = 16.0_WP, z23 = 84.0_WP, z32 = 19.0_WP
    REAL(WP)            :: f0, f1, fp, f2, f3, gmfp, gmf12, gmf23, fb


    f0    = rpm / 60.0_WP
    f1    = f0 * (zs + zr) / zs
    fp    = zs / zp * (f1 - f0)
    f2    = z12 / z21 * f1
    f3    = z23 / z32 * f2
    gmfp  = zp * fp
    gmf12 = z12 * f1
    gmf23 = z23 * f2
    fb    = 3.0_WP * f0

    RB1f = [f0, fb, 2.0_WP * fb, 3.0_WP * fb, 4.0_WP * fb, gmfp, f3, 2.0_WP * gmfp, 2.0_WP * f3, &
            3.0_WP * gmfp, 3.0_WP * f3, 4.0_WP * gmfp]
    RB2f = [f0, fb,              3.0_WP * fb, 4.0_WP * fb, gmfp, f3, 2.0_WP * gmfp, 2.0_WP * f3, &
            3.0_WP * gmfp, 3.0_WP * f3, 4.0_WP * gmfp]
    GB1f = [gmfp, f3, 2.0_WP * gmfp, 2.0_WP * f3, 3.0_WP * gmfp, 3.0_WP * f3, 4.0_WP * gmfp, 5.0_WP * gmfp, &
            6.0_WP * gmfp]
    GB2f = [gmfp, f3, gmf12, 2.0_WP * gmf12, 3.0_WP * gmf12, 4.0_WP * gmf12, 5.0_WP * gmf12, 6.0_WP * gmf12, &
            gmf23]
    GB3f = [gmf12, 2.0_WP * gmf12, 3.0_WP * gmf12, 4.0_WP * gmf12, 5.0_WP * gmf12, 6.0_WP * gmf12, &
            gmf23, 2.0_WP * gmf23, 3.0_WP * gmf23, 4.0_WP * gmf23, 5.0_WP * gmf23]

    f3_out = f3

END SUBROUTINE compute_freqs2MW


FUNCTION sorted_amps_by_freq(freqs_kin, amp_col) RESULT(amps_sorted)
    ! -------------------------------------------------------------------------
    ! Pairs a kinematic frequency family with its bar-chart amplitude column
    ! and sorts the pair by frequency, returning only the resulting amplitude
    ! order (matching the validated reference implementation, which re-pairs
    ! this frequency-sorted amplitude array positionally with the ORIGINAL,
    ! unsorted kinematic frequency array).
    ! -------------------------------------------------------------------------

    REAL(WP), INTENT(IN)  :: freqs_kin(:)      ! [Hz] Kinematic frequency family
    REAL(WP), INTENT(IN)  :: amp_col(:)        ! [-] Bar-chart amplitude column of the family
    REAL(WP), ALLOCATABLE :: amps_sorted(:)    ! [-] Amplitudes sorted by their paired frequency

    ! Local variables
    REAL(WP)    , ALLOCATABLE :: f_sort(:)
    INTEGER(I32)              :: n, i, k
    REAL(WP)                  :: fk, ak


    n = size(freqs_kin)
    ALLOCATE(f_sort(n), amps_sorted(n))
    f_sort      = freqs_kin
    amps_sorted = amp_col

    do i = 2, n
        fk = f_sort(i)
        ak = amps_sorted(i)
        k  = i - 1_I32
        do while (k >= 1_I32)
            if (f_sort(k) <= fk) exit
            f_sort(k+1)      = f_sort(k)
            amps_sorted(k+1) = amps_sorted(k)
            k                = k - 1_I32
        end do
        f_sort(k+1)      = fk
        amps_sorted(k+1) = ak
    end do

    DEALLOCATE(f_sort)

END FUNCTION sorted_amps_by_freq


SUBROUTINE append_family(freqs_kin, amps_kin, prefix, freqs_raw, amps_raw, keys_raw, j)
    ! -------------------------------------------------------------------------
    ! Appends one kinematic-frequency family (with its already frequency-
    ! sorted amplitudes) into the raw (freq, amp, key) accumulation arrays,
    ! advancing the shared cursor j and building a "<prefix>f.ffHz" label key.
    ! -------------------------------------------------------------------------

    REAL(WP)         , INTENT(IN)    :: freqs_kin(:)    ! [Hz] Kinematic frequency family
    REAL(WP)         , INTENT(IN)    :: amps_kin(:)     ! [-] Amplitudes paired with freqs_kin
    CHARACTER(len=*) , INTENT(IN)    :: prefix          ! [-] Label prefix of the family (e.g. 'RB1_f_')
    REAL(WP)         , INTENT(INOUT) :: freqs_raw(:)    ! [Hz] Accumulated frequencies
    REAL(WP)         , INTENT(INOUT) :: amps_raw(:)     ! [-] Accumulated amplitudes
    CHARACTER(len=20), INTENT(INOUT) :: keys_raw(:)     ! [-] Accumulated labels
    INTEGER(I32)     , INTENT(INOUT) :: j               ! [-] Shared cursor: last filled position

    ! Local variables
    INTEGER(I32) :: i


    do i = 1, size(freqs_kin)
        j            = j + 1_I32
        freqs_raw(j) = freqs_kin(i)
        amps_raw(j)  = amps_kin(i)
        write(keys_raw(j), '(A,F0.2,A)') trim(prefix), freqs_kin(i), 'Hz'
    end do

END SUBROUTINE append_family


SUBROUTINE read_bar_csv(filename, ncols, data)
    ! -------------------------------------------------------------------------
    ! Reads a comma-separated bar-chart amplitude file (one header line
    ! followed by ncols numeric columns, "Nan" marking a missing entry) into a
    ! (nrows, ncols) array, preserving IEEE NaN for missing entries so callers
    ! can strip them with ieee_is_nan. Blank lines are skipped in both passes,
    ! and a data row that cannot be read stops the program.
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, IEEE_QUIET_NAN

    CHARACTER(len=*), INTENT(IN)               :: filename     ! [-] Path to the bar-chart CSV file
    INTEGER(I32)    , INTENT(IN)               :: ncols        ! [-] Number of numeric columns to read
    REAL(WP)        , INTENT(OUT), ALLOCATABLE :: data(:,:)    ! [-] Bar-chart amplitudes (nrows, ncols), NaN if missing

    ! Local variables
    INTEGER(I32)        :: fu, ios, i, nrows
    CHARACTER(len=2048) :: line


    open(newunit=fu, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0_I32) error stop "read_bar_csv: could not open file " // trim(filename)

    read(fu, '(A)', iostat=ios) line    ! Header line

    nrows = 0_I32
    do
        read(fu, '(A)', iostat=ios) line
        if (ios /= 0_I32) exit
        if (len_trim(line) == 0) cycle
        nrows = nrows + 1_I32
    end do

    ALLOCATE(data(nrows, ncols))
    data = ieee_value(1.0_WP, IEEE_QUIET_NAN)

    rewind(fu)
    read(fu, '(A)') line    ! Skip header again

    i = 0_I32
    do while (i < nrows)
        read(fu, '(A)', iostat=ios) line
        if (ios /= 0_I32) exit
        if (len_trim(line) == 0) cycle
        i = i + 1_I32
        read(line, *, iostat=ios) data(i, 1:ncols)
        if (ios /= 0_I32) error stop "read_bar_csv: could not read a data row of " // trim(filename)
    end do

    close(fu)

END SUBROUTINE read_bar_csv


! ===============================================================================
! MINGYANG SMART ENERGY 3 MW DRIVETRAIN EXCITATION
! ===============================================================================
SUBROUTINE drivetrain3MW_excitation_spectrum(rpm, T_win, freqs_amp, keys)
    ! -------------------------------------------------------------------------
    ! Constructs the excitation spectrum (frequencies + amplitudes) of the
    ! MySE 3.0-135 drivetrain (two-stage planetary gearbox), the prescribed
    ! nacelle force of the Beleolico OpenFAST runs, reproducing the CONFIG of
    ! scripts/07_drivetrain_timeseries.py (Beleolico_offshore_turbine): pure
    ! tones at the rotor-locked orders N (f = N rpm / 60) of the stage-2 mesh
    ! (A: N = i1 Zr2 k = 540 k, k = 1-8), the stage-1 mesh (B: N = Zr1 k =
    ! 100 k, k = 1-30), the high-speed-shaft sidebands (C: N = 540 k +- i1 i2,
    ! k = 1-2) and the stage intermodulation (D: N = 540 k +- Zr1, k = 1-2).
    ! The amplitudes are relative to the first A line: the torque scaling of
    ! the force is common to every line and cancels in the normalization. As
    ! in the force files, only lines with exact frequency in [f_min, f_max]
    ! are kept, the frequencies are rounded to the 1/T_win grid of the
    ! analysis window of the runs, and lines closer than min_sep_bins grid
    ! steps are merged (amplitudes added, label of the largest one first, e.g.
    ! 'A5+B27'), so the frequencies are those of the force file of that rpm.
    ! Normalized to unit RMS.
    ! -------------------------------------------------------------------------

    REAL(WP)         , INTENT(IN)               :: rpm               ! [rpm] Rotor speed
    REAL(WP)         , INTENT(IN)               :: T_win             ! [s] Analysis window of the force file
    REAL(WP)         , INTENT(OUT), ALLOCATABLE :: freqs_amp(:,:)    ! [-] Peaks (N, 2): frequency [Hz], amplitude
    CHARACTER(len=20), INTENT(OUT), ALLOCATABLE :: keys(:)           ! [-] Peak labels (N)

    ! Local variables
    ! Gearbox tooth numbers (stage-1 ring and sun, stage-2 ring and sun)
    REAL(WP), PARAMETER :: zr1 = 100.0_WP, zs1 = 25.0_WP, zr2 = 108.0_WP, zs2 = 28.0_WP

    ! Harmonics k of the families, amplitude ratios to the A line of the same k and harmonic decay
    ! (S_k tabulated for k = 1-4, S_4 (4 / k)**s_exp beyond)
    INTEGER(I32)    , PARAMETER :: nk_A = 8_I32, nk_B = 30_I32, nk_sb = 2_I32
    REAL(WP)        , PARAMETER :: ratio_B = 3.2_WP
    REAL(WP)        , PARAMETER :: ratio_sb(2) = [0.5_WP, 0.25_WP]
    CHARACTER(len=1), PARAMETER :: name_sb(2) = ['C', 'D']
    REAL(WP)        , PARAMETER :: s_tab(4) = [1.0_WP, 0.29_WP, 0.15_WP, 0.10_WP], s_exp = 1.7_WP

    ! Band of the lines and merging distance [grid steps]
    REAL(WP)    , PARAMETER :: f_min = 1.0_WP, f_max = 1000.0_WP
    INTEGER(I32), PARAMETER :: min_sep_bins = 3_I32

    REAL(WP)                       :: i1, i2, N_A, N_HSS, f_exact, fk, ak, rms
    REAL(WP)                       :: N_sb(2)
    INTEGER(I32)                   :: n_raw, n_band, n_lines, i, j, k, kk, sb, best
    REAL(WP)         , ALLOCATABLE :: S_k(:), orders(:), amps_raw(:), freqs(:), amps(:)
    CHARACTER(len=20), ALLOCATABLE :: labels_raw(:), labels(:)
    CHARACTER(len=20)              :: key_tmp


    if (T_win <= 0.0_WP) error stop "drivetrain3MW_excitation_spectrum: T_win must be positive"

    ! --- 1. Gear orders and harmonic decay --- !
    i1    = 1.0_WP + zr1 / zs1
    i2    = 1.0_WP + zr2 / zs2
    N_A   = i1 * zr2
    N_HSS = i1 * i2
    N_sb  = [N_HSS, zr1]

    ALLOCATE(S_k(max(nk_A, nk_B, nk_sb)))
    do k = 1, size(S_k)
        if (k <= size(s_tab)) then
            S_k(k) = s_tab(k)
        else
            S_k(k) = s_tab(size(s_tab)) * (real(size(s_tab), WP) / real(k, WP))**s_exp
        end if
    end do

    ! --- 2. Orders, amplitudes and labels of every line (A, B, then C and D) --- !
    n_raw = nk_A + nk_B + 4_I32 * nk_sb
    ALLOCATE(orders(n_raw), amps_raw(n_raw), labels_raw(n_raw))

    j = 0_I32
    do k = 1, nk_A
        j           = j + 1_I32
        orders(j)   = N_A * real(k, WP)
        amps_raw(j) = S_k(k)
        write(labels_raw(j), '(A,I0)') 'A', k
    end do

    do k = 1, nk_B
        j           = j + 1_I32
        orders(j)   = zr1 * real(k, WP)
        amps_raw(j) = ratio_B * S_k(k)
        write(labels_raw(j), '(A,I0)') 'B', k
    end do

    do i = 1, 2
        do k = 1, nk_sb
            do sb = -1, 1, 2
                j           = j + 1_I32
                orders(j)   = N_A * real(k, WP) + real(sb, WP) * N_sb(i)
                amps_raw(j) = ratio_sb(i) * S_k(k)
                write(labels_raw(j), '(A,I0,A)') name_sb(i), k, merge('+', '-', sb > 0_I32)
            end do
        end do
    end do

    ! --- 3. Frequencies on the 1/T_win grid, only the lines in [f_min, f_max] --- !
    ALLOCATE(freqs(n_raw), amps(n_raw), labels(n_raw))
    n_band = 0_I32
    do j = 1, n_raw
        f_exact = orders(j) * rpm / 60.0_WP
        if (f_exact < f_min .or. f_exact > f_max) cycle

        n_band         = n_band + 1_I32
        freqs(n_band)  = anint(f_exact * T_win) / T_win
        amps(n_band)   = amps_raw(j)
        labels(n_band) = labels_raw(j)
    end do
    DEALLOCATE(S_k, orders, amps_raw, labels_raw)

    ! --- 4. Sort by frequency (carrying amps/labels along, stable insertion sort) --- !
    do i = 2, n_band
        fk      = freqs(i)
        ak      = amps(i)
        key_tmp = labels(i)
        kk      = i - 1_I32
        do while (kk >= 1_I32)
            if (freqs(kk) <= fk) exit
            freqs(kk+1)  = freqs(kk)
            amps(kk+1)   = amps(kk)
            labels(kk+1) = labels(kk)
            kk           = kk - 1_I32
        end do
        freqs(kk+1)  = fk
        amps(kk+1)   = ak
        labels(kk+1) = key_tmp
    end do

    ! --- 5. Merge lines closer than min_sep_bins grid steps to the previous one of their group --- !
    ALLOCATE(freqs_amp(n_band, 2), keys(n_band))
    n_lines = 0_I32
    i       = 1_I32
    do while (i <= n_band)
        kk = i
        do while (kk + 1_I32 <= n_band)
            if (nint((freqs(kk+1) - freqs(kk)) * T_win) < min_sep_bins) then
                kk = kk + 1_I32
            else
                exit
            end if
        end do

        best    = i - 1_I32 + maxloc(amps(i:kk), dim=1)
        key_tmp = labels(best)
        do j = i, kk
            if (j /= best) key_tmp = trim(key_tmp) // '+' // trim(labels(j))
        end do

        n_lines               = n_lines + 1_I32
        freqs_amp(n_lines, 1) = freqs(best)
        freqs_amp(n_lines, 2) = sum(amps(i:kk))
        keys(n_lines)         = key_tmp
        i                     = kk + 1_I32
    end do
    DEALLOCATE(freqs, amps, labels)

    if (n_lines < n_band) then
        freqs_amp = freqs_amp(1:n_lines, :)
        keys      = keys(1:n_lines)
    end if

    ! --- 6. Normalize to unit RMS: sqrt(sum(A^2/2)) = 1 --- !
    if (n_lines > 0_I32) then
        rms             = sqrt(sum(freqs_amp(:, 2)**2) / 2.0_WP)
        freqs_amp(:, 2) = freqs_amp(:, 2) / rms
    end if

END SUBROUTINE drivetrain3MW_excitation_spectrum


SUBROUTINE read_drivetrain3MW_force_file(path, rpm, T_win)
    ! -------------------------------------------------------------------------
    ! Reads the rotor speed and the analysis window of a drivetrain force
    ! file of the Beleolico runs from its description <rootname>.force.json
    ! (written by 07_drivetrain_timeseries.py and copied next to the .outb by
    ! the run scripts): its first "rpm" (the top-level one, the rotor speed of
    ! the force lines) and "T_win" (of its configuration). A missing file or
    ! key stops the program.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)  :: path     ! [-] Path to the .force.json file
    REAL(WP)        , INTENT(OUT) :: rpm      ! [rpm] Rotor speed of the force lines
    REAL(WP)        , INTENT(OUT) :: T_win    ! [s] Analysis window of the force lines (grid 1/T_win)


    rpm   = read_json_number(path, 'rpm')
    T_win = read_json_number(path, 'T_win')

END SUBROUTINE read_drivetrain3MW_force_file


FUNCTION read_json_number(path, key) RESULT(val)
    ! -------------------------------------------------------------------------
    ! Returns the number of the first '"key": value' pair of a JSON file
    ! written with one pair per line (Python json.dump with indent). A file
    ! that cannot be opened, a missing key or a value that is not a number
    ! stops the program.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: path    ! [-] Path to the JSON file
    CHARACTER(len=*), INTENT(IN) :: key     ! [-] Key, without quotes
    REAL(WP)                     :: val     ! [-] Value of the first pair with that key

    ! Local variables
    INTEGER(I32)        :: fu, ios, pos
    CHARACTER(len=2048) :: line


    open(newunit=fu, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0_I32) error stop "read_json_number: could not open file " // trim(path)

    do
        read(fu, '(A)', iostat=ios) line
        if (ios /= 0_I32) exit

        pos = index(line, '"' // key // '":')
        if (pos == 0_I32) cycle

        read(line(pos+len(key)+3:), *, iostat=ios) val
        close(fu)
        if (ios /= 0_I32) error stop "read_json_number: value of " // key // " is not a number in " // trim(path)
        return
    end do

    close(fu)
    error stop "read_json_number: key " // key // " not found in " // trim(path)

END FUNCTION read_json_number

END MODULE DrivetrainExcitation
