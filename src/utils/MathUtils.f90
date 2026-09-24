! ===============================================================================
! MODULE: MathUtils
! ===============================================================================
! Grab-bag of numerical and signal-processing helpers shared across the code-
! base: elapsed wall-clock time formatting, duplicate structural-node removal,
! FFTW-based real FFT of time series, synthetic drivetrain excitation signal
! generation, frequency-mask filtering, a C qsort wrapper for sorting reals,
! and Bessel/Hankel-function acoustic scattering helpers (K0/K1 modified
! Bessel approximations) used to build the impedance/shielding correction.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE MathUtils

USE FFTW3
USE omp_lib
USE Kinds, ONLY: WP, I32, PI
USE, INTRINSIC :: iso_c_binding

IMPLICIT NONE
PUBLIC

! -------------------------------------------------------------------------
! C interface: qsort
! Binds to the C standard library's qsort(3), used by sort_real_c to sort
! a REAL(WP) array in ascending order via the compare_real callback.
! -------------------------------------------------------------------------
INTERFACE
    SUBROUTINE qsort(base, nmemb, size, compar) BIND(c, name='qsort')
        IMPORT :: C_PTR, C_SIZE_T, C_FUNPTR
        TYPE(C_PTR), VALUE :: base
        INTEGER(C_SIZE_T), VALUE :: nmemb
        INTEGER(C_SIZE_T), VALUE :: size
        TYPE(C_FUNPTR), VALUE :: compar
    END SUBROUTINE qsort
END INTERFACE


CONTAINS


! -------------------------------------------------------------------------
! Subroutine: format_elapsed
! Converts a wall-clock elapsed time (in seconds, measured via OpenMP's
! timer) into a human-readable value/unit pair, automatically switching
! among seconds, minutes, hours and days depending on magnitude.
! -------------------------------------------------------------------------
SUBROUTINE format_elapsed(start_time, elapsed_display, tag)
    REAL(WP), INTENT(IN) :: start_time
    REAL(WP), INTENT(OUT) :: elapsed_display
    CHARACTER(len=*), INTENT(OUT) :: tag
    REAL(WP) :: end_time, elapsed_time
    
    ! CORRECTED: Use OpenMP wall-clock timer instead of processor CPU_TIME
    end_time = real(omp_get_wtime(), WP)
    elapsed_time = end_time - start_time
    elapsed_display = elapsed_time
    
    if (elapsed_display > 60.0_WP .and. elapsed_display <= 3600.0_WP) then
        elapsed_display = elapsed_display / 60.0_WP
        tag = "min"
    elseif (elapsed_display > 3600.0_WP .and. elapsed_display <= 86400.0_WP) then
        elapsed_display = elapsed_display / 3600.0_WP
        tag = "h"
    elseif (elapsed_display > 86400.0_WP) then
        elapsed_display = elapsed_display / 86400.0_WP
        tag = "days"
    else
        tag = "s"
    end if
END SUBROUTINE format_elapsed


! -------------------------------------------------------------------------
! Function: remove_duplicate_nodes
! Removes duplicate nodes shared between consecutive structural members.
! Given an array of nodes of shape (Nm, Nn, 3), where Nm is the number of
! members and Nn the number of nodes per member, identifies nodes that are
! duplicated at the junction between member i and member i+1 (the last
! node of member i is assumed to coincide with the first node of member
! i+1) and marks one of the two for removal, based on delete_last.
!
! Arguments:
!   x_all      : REAL(WP), dimension(:,:,:)  (Nm, Nn, 3)
!   delete_last: LOGICAL  If .true., the first node of the following member
!                          is removed; if .false., the last node of the
!                          previous member is removed.
! Returns:
!   mask       : LOGICAL, dimension(Nm*Nn)  Flat mask with .true. for nodes to keep
! -------------------------------------------------------------------------
FUNCTION remove_duplicate_nodes(x_all, delete_last) RESULT(mask)
    REAL(WP), INTENT(IN) :: x_all(:,:,:)   ! (Nm, Nn, 3)
    LOGICAL,  INTENT(IN) :: delete_last
    LOGICAL, ALLOCATABLE :: mask(:)

    integer(i32) :: Nm, Nn, i, idx_prev, idx_next
    REAL(WP) :: tol
    LOGICAL :: are_equal

    Nm = size(x_all, 1)
    Nn = size(x_all, 2)

    allocate(mask(Nm*Nn))
    mask = .true.   ! start keeping all nodes

    ! Tolerance for comparing coordinates (relative)
    tol = 1.0e-6_WP

    do i = 1, Nm-1
        ! Last node of current member
        idx_prev = (i-1)*Nn + Nn
        ! First node of next member
        idx_next = i*Nn + 1

        ! Compare coordinates (simple Euclidean distance, but we can check
        ! each coordinate difference)
        are_equal = all(abs(x_all(i, Nn, :) - x_all(i+1, 1, :)) < tol)

        if (are_equal) then
            if (delete_last) then
                ! Remove first node of next member
                mask(idx_next) = .false.
            else
                ! Remove last node of current member
                mask(idx_prev) = .false.
            end if
        end if
    end do

END FUNCTION remove_duplicate_nodes


! -------------------------------------------------------------------------
! Function: divide_span
! Computes representative span intervals (dx) for an ordered spatial vector.
! Fully vectorized to eliminate conditionals and temporary arrays.
! -------------------------------------------------------------------------
FUNCTION divide_span(x) RESULT(dx)
    REAL(WP), INTENT(IN) :: x(:)
    REAL(WP), ALLOCATABLE :: dx(:)

    INTEGER(I32) :: n
    REAL(WP), ALLOCATABLE :: x_work(:)
    LOGICAL :: flipped

    n = size(x)
    ALLOCATE(dx(n))

    if (n == 0) return

    ALLOCATE(x_work(n))

    ! Check direction and possibly reverse
    flipped = .false.
    if (n > 1 .and. x(1) > x(n)) then
        x_work = x(n:1:-1)
        flipped = .true.
    else
        x_work = x
    end if

    ! Vectorized span calculation via algebraic simplification
    if (n == 1) then
        dx(1) = 0.0_WP
    else
        dx(1)     = 0.5_WP * (x_work(2) - x_work(1))
        dx(2:n-1) = 0.5_WP * (x_work(3:n) - x_work(1:n-2))
        dx(n)     = 0.5_WP * (x_work(n) - x_work(n-1))
    end if

    ! If we had reversed, reverse dx back
    if (flipped) then
        dx = dx(n:1:-1)
    end if

    ! Check consistency: all dx must have the same sign
    if (.not. (all(dx > 0.0_WP) .or. all(dx < 0.0_WP))) then
        error stop "divide_span: check input, non-monotonic grid"
    end if

END FUNCTION divide_span


! -------------------------------------------------------------------------
! Function: integral_elipse_z
! Area under the horizontal line z = k inside the ellipse
! (2y/Dy)^2 + (2z/Dz)^2 <= 1, i.e. the submerged cross-sectional area of an
! elliptical member (semi-axes Dy/2, Dz/2) up to elevation k. Used to build
! the outer/inner submerged shell area of the SAITEC elliptical floaters.
! -------------------------------------------------------------------------
PURE FUNCTION integral_elipse_z(Dy, Dz, k) RESULT(area)
    REAL(WP), INTENT(IN) :: Dy   ! [m] Diameter in y-direction
    REAL(WP), INTENT(IN) :: Dz   ! [m] Diameter in z-direction
    REAL(WP), INTENT(IN) :: k    ! [m] z(y) = k line
    REAL(WP) :: area

    REAL(WP) :: a, b, elipse_area, k_clip, u

    a = Dy / 2.0_WP
    b = Dz / 2.0_WP

    elipse_area = PI * a * b

    k_clip = min(max(k, -b), b)
    u = k_clip / b

    area = 0.5_WP * elipse_area + a * b * (asin(u) + u * sqrt(1.0_WP - u*u))

END FUNCTION integral_elipse_z


! -------------------------------------------------------------------------
! Subroutine: compute_rfft
! Computes the Real Fast Fourier Transform (RFFT) of a 3D time-series array
! via FFTW. Optimized for HPC: eliminates all intermediate arrays, fuses
! the mean-removal and windowing passes, and performs direct-to-final-size
! allocation to avoid memory thrashing.
! -------------------------------------------------------------------------
SUBROUTINE compute_rfft(array, nt, dt, skipf, remove_zero, array_out, freqs)
    REAL(WP)                , INTENT(IN)           :: array(:,:,:)
    INTEGER(I32)            , INTENT(IN)           :: nt
    REAL(WP)                , INTENT(IN)           :: dt
    INTEGER(I32)            , INTENT(IN), OPTIONAL :: skipf
    LOGICAL                 , INTENT(IN), OPTIONAL :: remove_zero
    COMPLEX(WP), ALLOCATABLE, INTENT(OUT)          :: array_out(:,:,:)
    REAL(WP), ALLOCATABLE   , INTENT(OUT)          :: freqs(:)

    ! Local variables
    INTEGER(I32) :: nfreq_raw, nfreq_final, n_signals, i, j, k, idx, pos, start_idx
    INTEGER(I32) :: skipf_, n2, n3
    LOGICAL :: remove_zero_
    REAL(WP) :: win_norm, mean_val, scale_fact
    REAL(WP), ALLOCATABLE :: window(:), data_flat(:,:)
    COMPLEX(WP), ALLOCATABLE :: out_flat(:,:)

    ! FFTW variables
    TYPE(C_PTR) :: plan
    INTEGER(C_INT) :: rank, howmany, istride, idist, ostride, odist, flags
    INTEGER(C_INT), ALLOCATABLE :: n_(:), inembed(:), onembed(:)

    ! Defaults
    skipf_ = 1           ; if (PRESENT(skipf))       skipf_       = skipf
    remove_zero_ = .true.; if (PRESENT(remove_zero)) remove_zero_ = remove_zero

    n2 = size(array,2)
    n3 = size(array,3)
    n_signals = n2 * n3
    nfreq_raw = nt / 2 + 1

    ! 1. Compute window and its norm once
    ALLOCATE(window(nt))
    win_norm = 0.0_WP
    !$omp parallel do reduction(+:win_norm)
    do i = 1, nt
        window(i) = 0.5_WP - 0.5_WP*cos(2.0_WP*PI * REAL(i-1,WP) / REAL(nt-1, WP))
        win_norm = win_norm + window(i)
    end do
    !$omp end parallel do
    win_norm = win_norm / REAL(nt, WP)

    ! 2. Fused Mean, Windowing, and Flattening (No temporary 3D arrays)
    ALLOCATE(data_flat(nt, n_signals))
    
    !$omp parallel do private(idx, j, k, i, mean_val)
    do idx = 1, n_signals
        ! Map flat index back to j, k to read directly from 3D array
        k = (idx - 1) / n2 + 1
        j = mod(idx - 1, n2) + 1
        
        ! Pass A: Compute mean
        mean_val = 0.0_WP
        do i = 1, nt
            mean_val = mean_val + array(i, j, k)
        end do
        mean_val = mean_val / REAL(nt, WP)
        
        ! Pass B: Subtract mean, apply window, store flat (column-major)
        do i = 1, nt
            data_flat(i, idx) = (array(i, j, k) - mean_val) * window(i)
        end do
    end do
    !$omp end parallel do
    DEALLOCATE(window)

    ! 3. FFTW Execution
    ALLOCATE(out_flat(nfreq_raw, n_signals))
    rank = 1
    ALLOCATE(n_(1), inembed(1), onembed(1))
    n_(1)      = nt
    inembed(1) = nt
    onembed(1) = nfreq_raw
    howmany    = n_signals
    istride    = 1
    idist      = nt
    ostride    = 1
    odist      = nfreq_raw
    flags      = FFTW_ESTIMATE

    plan = fftw_plan_many_dft_r2c(rank, n_, howmany, &
                                  data_flat, inembed, istride, idist, &
                                  out_flat, onembed, ostride, odist, flags)

    if (.not. c_associated(plan)) error stop "compute_rfft: FFTW plan creation failed"
    CALL fftw_execute_dft_r2c(plan, data_flat, out_flat)
    CALL fftw_destroy_plan(plan)
    
    DEALLOCATE(data_flat, n_, inembed, onembed)

    ! 4. Direct-to-Final Allocation (Eliminates array re-allocations)
    start_idx = 1
    if (remove_zero_) start_idx = 2

    ! Pre-calculate final exact size
    nfreq_final = 0
    do i = start_idx, nfreq_raw, skipf_
        nfreq_final = nfreq_final + 1
    end do

    ALLOCATE(array_out(nfreq_final, n2, n3))
    ALLOCATE(freqs(nfreq_final))

    ! Fill frequencies
    pos = 1
    do i = start_idx, nfreq_raw, skipf_
        freqs(pos) = REAL(i-1, WP) / (REAL(nt, WP) * dt)
        pos = pos + 1
    end do

    ! Scale and map directly to 3D array_out
    scale_fact = 2.0_WP / (REAL(nt, WP) * win_norm)

    !$omp parallel do private(idx, j, k, pos, i)
    do idx = 1, n_signals
        k = (idx - 1) / n2 + 1
        j = mod(idx - 1, n2) + 1
        
        pos = 1
        ! Inner loop over frequencies guarantees contiguous writes to array_out
        do i = start_idx, nfreq_raw, skipf_
            array_out(pos, j, k) = out_flat(i, idx) * scale_fact
            pos = pos + 1
        end do
    end do
    !$omp end parallel do

    DEALLOCATE(out_flat)

END SUBROUTINE compute_rfft


! -------------------------------------------------------------------------
! Subroutine: generate_timeseries_banded_sines
! Reconstructs a synthetic time-series signal from spectral peak
! distributions. Optimized for HPC: aggregates all frequency components
! first, then performs a single OpenMP+SIMD pass over the time domain to
! maximize cache locality.
! -------------------------------------------------------------------------
SUBROUTINE generate_timeseries_banded_sines(peaks, keys, t , zeta, nfreq, seed, &
                                            used_freqs, fcut, a, freqs_out)
    REAL(WP)         , INTENT(IN) :: peaks(:,:)
    CHARACTER(len=20), INTENT(IN) :: keys(:)
    REAL(WP)         , INTENT(IN) :: t(:)
    REAL(WP)         , INTENT(IN), OPTIONAL :: zeta
    INTEGER(I32)     , INTENT(IN), OPTIONAL :: nfreq
    INTEGER(I32)     , INTENT(IN), OPTIONAL :: seed
    LOGICAL          , INTENT(IN), OPTIONAL :: used_freqs
    REAL(WP)         , INTENT(IN), OPTIONAL :: fcut
    REAL(WP), ALLOCATABLE, INTENT(OUT) :: a(:)
    REAL(WP), ALLOCATABLE, INTENT(OUT), OPTIONAL :: freqs_out(:)

    ! Local variables
    INTEGER(I32) :: i, j, seed_size, nt, n_peaks, n_active, nfreq_, seed_
    INTEGER(I32) :: max_freq_list, n_kept
    REAL(WP) :: f0, Arms, sigma, phase, rnd, zeta_, fcut_, sqrt2
    REAL(WP), ALLOCATABLE :: fk(:), g(:), Arms_k(:), phases(:)
    
    ! Global accumulation arrays for all sine wave components
    REAL(WP), ALLOCATABLE :: all_omega(:), all_amp(:), all_phase(:), freq_list(:)
    INTEGER(I32), ALLOCATABLE :: seed_arr(:)
    LOGICAL :: used_freqs_
    CHARACTER(len=20) :: key_label
    
    ! Time loop variables
    REAL(WP) :: current_t, current_sum

    ! Defaults
    sqrt2       = sqrt(2.0_WP)
    zeta_       = 0.02_WP; if (PRESENT(zeta))       zeta_       = zeta
    nfreq_      = 50_I32 ; if (PRESENT(nfreq))      nfreq_      = nfreq
    seed_       = 42_I32 ; if (PRESENT(seed))       seed_       = seed
    fcut_       = 10.0_WP; if (PRESENT(fcut))       fcut_       = fcut
    used_freqs_ = .false.; if (PRESENT(used_freqs)) used_freqs_ = used_freqs
    
    nt = size(t)
    ALLOCATE(a(nt))
    n_peaks = size(peaks,1)
    
    if (size(keys) < n_peaks) error stop "generate_timeseries_banded_sines: keys shorter than peaks"
    if (size(peaks,2) < 2) error stop "generate_timeseries_banded_sines: peaks must have at least 2 columns"

    ! Pre-allocate global lists based on the maximum possible components
    max_freq_list = max(1_I32, n_peaks * nfreq_)
    ALLOCATE(freq_list(max_freq_list))
    ALLOCATE(all_omega(max_freq_list), all_amp(max_freq_list), all_phase(max_freq_list))
    n_active = 0_I32
    
    ! Initialize random seed
    CALL random_seed(size=seed_size)
    ALLOCATE(seed_arr(seed_size))
    CALL random_seed(put=seed_arr)
    
    ! Phase 1: Aggregate all frequency components (No time-domain computation here)
    do i = 1, n_peaks
        f0   = peaks(i,1)
        Arms = peaks(i,2)
        key_label = trim(keys(i))

        ! Gear mesh or high frequency -> banded Gaussian
        if (index(key_label, "mesh") > 0 .or. f0 > fcut_) then
            sigma = zeta_ * f0
            if (sigma <= 0.0_WP) cycle

            ALLOCATE(fk(nfreq_), g(nfreq_), Arms_k(nfreq_), phases(nfreq_))
            
            do j = 1, nfreq_
                fk(j) = f0 - 4.0_WP * sigma + (j-1) * (8.0_WP * sigma) / REAL(nfreq_-1, WP)
            end do

            g = exp(-0.5_WP * ((fk-f0)/sigma)**2)
            g = g / sqrt(sum(g**2))
            Arms_k = Arms * g

            do j = 1, nfreq_
                CALL random_number(rnd)
                phases(j) = 2.0_WP * PI * rnd
                
                ! Store properties globally
                n_active = n_active + 1
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
            n_active = n_active + 1
            freq_list(n_active) = f0
            all_omega(n_active) = 2.0_WP * PI * f0
            all_amp(n_active)   = sqrt2 * Arms
            all_phase(n_active) = phase
        end if
    end do

    ! Phase 2: HPC Time-domain Synthesis
    ! Time is the outer loop. We write to main memory array 'a' exactly once per time step.
    !$omp parallel do private(i, j, current_t, current_sum)
    do i = 1, nt
        current_t = t(i)
        current_sum = 0.0_WP
        
        ! Inner loop over frequencies: Sum into an ultra-fast CPU register via SIMD
        !$omp simd reduction(+:current_sum)
        do j = 1, n_active
            current_sum = current_sum + all_amp(j) * sin(all_omega(j) * current_t + all_phase(j))
        end do
        
        a(i) = current_sum
    end do
    !$omp end parallel do

    ! Phase 3: Frequency Sorting and Cleanup
    if (used_freqs_) then
        if (n_active > 0) then
            CALL sort_real_c(freq_list(1:n_active))
            ALLOCATE(freqs_out(n_active))
            n_kept = 0_I32
            do i = 1, n_active
                if (i == 1) then
                    n_kept            = n_kept + 1
                    freqs_out(n_kept) = freq_list(i)
                else if ( abs(freq_list(i) - freq_list(i-1)) > 1.0e-6_WP ) then
                    n_kept            = n_kept + 1
                    freqs_out(n_kept) = freq_list(i)
                end if
            end do
            if (n_kept < n_active) freqs_out = freqs_out(1:n_kept)
        else
             ALLOCATE(freqs_out(0))
        end if
    else
        if (PRESENT(freqs_out)) ALLOCATE (freqs_out(0))
    end if

    DEALLOCATE(freq_list, seed_arr, all_omega, all_amp, all_phase)

END SUBROUTINE generate_timeseries_banded_sines


! -------------------------------------------------------------------------
! Function: filter_non_useful_freqs
! Builds a logical mask selecting which frequencies to keep. Frequencies
! outside the [freqs_over, freqs_under] zone are always kept; frequencies
! inside that zone are kept only if they lie within one frequency bin (df)
! of one of the maintain_freqs values, so that physically relevant peaks
! are preserved while the rest of the band is discarded.
! -------------------------------------------------------------------------
FUNCTION filter_non_useful_freqs(freqs, maintain_freqs, freqs_over, freqs_under) RESULT(mask)
    REAL(WP), INTENT(IN) :: freqs(:)
    REAL(WP), INTENT(IN) :: maintain_freqs(:)
    REAL(WP), INTENT(IN), OPTIONAL :: freqs_over, freqs_under
    LOGICAL, ALLOCATABLE :: mask(:)

    INTEGER(I32) :: n, j
    REAL(WP) :: df, freqs_over_, freqs_under_
    LOGICAL, ALLOCATABLE :: filter_zone(:)

    n = SIZE(freqs)
    ALLOCATE(mask(n))
    ALLOCATE(filter_zone(n))

    ! Determine frequency step
    IF (n > 1) THEN
        df = freqs(2) - freqs(1)
        ! Optionally check for equal spacing (omit for performance)
    ELSE
        df = 0.0_WP
    END IF

    ! Define filter zone (frequencies to consider for removal)
    filter_zone = .TRUE.
    freqs_over_  = MINVAL(freqs)
    freqs_under_ = MAXVAL(freqs)
    IF (PRESENT(freqs_over))  freqs_over_  = freqs_over
    IF (PRESENT(freqs_under)) freqs_under_ = freqs_under

    IF (PRESENT(freqs_over))  filter_zone = filter_zone .AND. (freqs >= freqs_over_)
    IF (PRESENT(freqs_under)) filter_zone = filter_zone .AND. (freqs <= freqs_under_)

    ! Outside filter zone: keep all frequencies
    mask = .NOT. filter_zone

    ! Inside filter zone: keep only if within df of any maintain_freqs
    DO j = 1, SIZE(maintain_freqs)
        WHERE (filter_zone .AND. ABS(freqs - maintain_freqs(j)) <= df)
            mask = .TRUE.
        END WHERE
    END DO

    DEALLOCATE(filter_zone)
END FUNCTION filter_non_useful_freqs


! -------------------------------------------------------------------------
! Function: compare_real
! Comparison callback for qsort: ascending order for real(c_double) values.
! -------------------------------------------------------------------------
FUNCTION compare_real(a,b) BIND(c) RESULT(res)
    TYPE(C_PTR), VALUE :: a, b
    INTEGER(C_INT) :: res
    REAL(C_DOUBLE), POINTER :: fa, fb

    CALL c_f_pointer(a, fa)
    CALL c_F_pointer(b, fb)

    if (fa < fb) then
        res = -1_C_INT
    else if (fa > fb) then
        res = 1_C_INT
    else
        res = 0_C_INT
    end if

END FUNCTION compare_real


! -------------------------------------------------------------------------
! Subroutine: sort_real_c
! Sorts a REAL(WP) array in ascending order using the C standard library's
! qsort, via the compare_real comparator callback.
! -------------------------------------------------------------------------
SUBROUTINE sort_real_c(arr)
    REAL(WP), INTENT(INOUT), TARGET :: arr(:)
    REAL(WP), POINTER :: parr(:)
    INTEGER(C_SIZE_T) :: nmemb, size_elem
    TYPE(C_PTR) :: ptr

    parr => arr
    nmemb = size(parr, kind=C_SIZE_T)
    size_elem = c_sizeof(parr(1))
    ptr = c_loc(parr(1))

    CALL qsort(ptr, nmemb, size_elem, c_funloc(compare_real))
END SUBROUTINE sort_real_c


! -------------------------------------------------------------------------
! Function: alpha_hankel
! Analytical scattering-correction coefficient alpha, based on cylindrical
! Hankel functions (orders 0 and 2, combined into the derivative of the
! order-1 function). Used to correct the radiated acoustic force for the
! finite diameter of a cylindrical structural member.
! -------------------------------------------------------------------------
FUNCTION alpha_hankel(k, D) RESULT(alpha)
    USE Kinds, ONLY: I1

    REAL(WP), INTENT(in) :: k(:)            ! [1/m] Wavenumber array shape(nfreqs)
    REAL(WP), INTENT(in) :: D               ! [m] Cylinder diameter

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: H0(:), H2(:), H1p(:), alpha(:)
    REAL(WP)   , ALLOCATABLE :: ka(:)
    INTEGER(I32)             :: Nfreqs

    Nfreqs = size(k)
    ka = k * D / 2.0_WP

    ALLOCATE(H0(Nfreqs), H2(Nfreqs), H1p(Nfreqs), alpha(Nfreqs))
        
    H0 = CMPLX(bessel_jn(0, ka), -bessel_yn(0, ka), kind=WP)      ! H_0^(2) order 0 second specie (H_n^(2) = Jn - iYn)
    H2 = CMPLX(bessel_jn(2, ka), -bessel_yn(2, ka), kind=WP)      ! H_2^(2) order 2 second specie (H_n^(2) = Jn - iYn)
    H1p = 0.5_WP * (H0 - H2)                                      ! dH_1^(1) order 1 first specie (H_n^(1) = Jn + iYn)

    alpha = -3.0_WP * I1 /(PI * ka**2 * H1p)

END FUNCTION alpha_hankel


! -------------------------------------------------------------------------
! Subroutine: compute_hankel_complex
! Efficiently computes the Hankel functions of the first kind (H0 and H1)
! for a wavenumber that is either purely real (Propagative Mode, using the
! intrinsic Bessel functions) or purely imaginary (Evanescent Mode, using
! the bessel_k0/bessel_k1 polynomial approximations below).
! -------------------------------------------------------------------------
SUBROUTINE compute_hankel_complex(k_rm, r, H0, H1)
    COMPLEX(WP), INTENT(IN)  :: k_rm
    REAL(WP), INTENT(IN)     :: r
    COMPLEX(WP), INTENT(OUT) :: H0, H1

    ! Local variables
    REAL(WP) :: Z_real, Y

    if (abs(aimag(k_rm)) <= 1.0e-12_WP) then
        ! Propagative Mode: k_rm is Real
        Z_real = real(k_rm, WP) * r
        H0 = CMPLX(bessel_jn(0, Z_real), bessel_yn(0, Z_real), kind=WP)
        H1 = CMPLX(bessel_jn(1, Z_real), bessel_yn(1, Z_real), kind=WP)
    else
        ! Evanescent Mode: k_rm is Pure Imaginary
        Y = abs(aimag(k_rm)) * r
        H0 = CMPLX(0.0_WP, -2.0_WP / PI * bessel_k0(Y), kind=WP)
        H1 = CMPLX(-2.0_WP / PI * bessel_k1(Y), 0.0_WP, kind=WP)   
    end if

END SUBROUTINE compute_hankel_complex


! -------------------------------------------------------------------------
! Function: bessel_k0
! Abramowitz and Stegun polynomial approximation (9.8.5 - 9.8.6) of the
! modified Bessel function of the second kind, order 0: K0(x) for x > 0.
! Used to build the evanescent-mode Hankel function for the cylindrical
! scattering correction.
! -------------------------------------------------------------------------
PURE FUNCTION bessel_k0(x) RESULT(k0)
    REAL(WP), INTENT(IN) :: x
    REAL(WP) :: k0, y, t, i0

    if (x <= 0.0_WP) then
        k0 = 0.0_WP
        return
    end if

    if (x <= 2.0_WP) then
        y = (x / 3.75_WP)**2
        i0 = 1.0_WP + y * (3.5156229_WP + y * (3.0899424_WP + y * (1.2067492_WP + &
             y * (0.2659732_WP + y * (0.0360768_WP + y * 0.0045813_WP)))))

        t = x / 2.0_WP
        y = t * t
        k0 = -log(t) * i0 + (-0.57721566_WP + y * (0.42278420_WP + y * (0.23069756_WP + &
             y * (0.03488590_WP + y * (0.00262698_WP + y * (0.00010750_WP + &
             y * 0.00000740_WP))))))
    else
        y = 2.0_WP / x
        k0 = (exp(-x) / sqrt(x)) * (1.25331414_WP + y * (-0.07832358_WP + &
             y * (0.02189568_WP + y * (-0.01062449_WP + y * (0.00587872_WP + &
             y * (-0.00251540_WP + y * 0.00053208_WP))))))
    end if
END FUNCTION bessel_k0


! -------------------------------------------------------------------------
! Function: bessel_k1
! Abramowitz and Stegun polynomial approximation (9.8.9 - 9.8.10) of the
! modified Bessel function of the second kind, order 1: K1(x) for x > 0.
! Used to build the evanescent-mode Hankel function for the cylindrical
! scattering correction.
! -------------------------------------------------------------------------
PURE FUNCTION bessel_k1(x) RESULT(k1)
    REAL(WP), INTENT(IN) :: x
    REAL(WP) :: k1, y, t, i1

    if (x <= 0.0_WP) then
        k1 = 0.0_WP
        return
    end if

    if (x <= 2.0_WP) then
        y = (x / 3.75_WP)**2
        i1 = x * (0.5_WP + y * (0.87890594_WP + y * (0.51498869_WP + &
             y * (0.15084934_WP + y * (0.02658733_WP + y * (0.00301532_WP + &
             y * 0.00032411_WP))))))

        t = x / 2.0_WP
        y = t * t
        k1 = log(t) * i1 + (1.0_WP / x) * (1.0_WP + y * (0.15443144_WP + &
             y * (-0.67278579_WP + y * (-0.18156897_WP + y * (-0.01919402_WP + &
             y * (-0.00110404_WP + y * (-0.00004686_WP)))))))
    else
        y = 2.0_WP / x
        k1 = (exp(-x) / sqrt(x)) * (1.25331414_WP + y * (0.23498619_WP + &
             y * (-0.03655620_WP + y * (0.01504268_WP + y * (-0.00780353_WP + &
             y * (0.00325614_WP + y * (-0.00068245_WP)))))))
    end if
END FUNCTION bessel_k1

END MODULE MathUtils