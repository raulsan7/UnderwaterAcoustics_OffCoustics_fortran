! ===============================================================================
! MODULE: SignalProcessing
! ===============================================================================
! Defines the signal-processing helpers: the FFTW-based real FFT of the
! nodal time series, with mean removal and Hann windowing (compute_rfft).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE SignalProcessing

USE FFTW3
USE Kinds, ONLY: I32, WP, PI
USE, INTRINSIC :: iso_c_binding, ONLY: C_INT, C_PTR, c_associated

IMPLICIT NONE

PRIVATE
PUBLIC :: compute_rfft


CONTAINS


SUBROUTINE compute_rfft(array, nt, dt, skipf, remove_zero, array_out, freqs)
    ! -------------------------------------------------------------------------
    ! Computes the Real Fast Fourier Transform (RFFT) of a 3D time-series array
    ! via FFTW. Optimized for HPC: eliminates all intermediate arrays, fuses
    ! the mean-removal and windowing passes, and performs direct-to-final-size
    ! allocation to avoid memory thrashing. The FFTW forward transform uses
    ! exp(-i*w*t) (time factor exp(+i*w*t)), while the acoustic kernels of the
    ! solvers (exp(+i*k*r), Hankel functions of the first kind) assume the
    ! time factor exp(-i*w*t). The spectrum is therefore returned conjugated,
    ! X(f) = sum x(t) exp(+i*w*t), so that every propagated wave is outgoing
    ! (causal) and the relative phases between sources are kept.
    ! -------------------------------------------------------------------------

    REAL(WP)    , INTENT(IN)               :: array(:,:,:)        ! [-] Time series (nt, N2, N3)
    INTEGER(I32), INTENT(IN)               :: nt                  ! [-] Number of time samples
    REAL(WP)    , INTENT(IN)               :: dt                  ! [s] Time step
    INTEGER(I32), INTENT(IN), OPTIONAL     :: skipf               ! [-] Keep one out of every skipf frequency bins
    LOGICAL     , INTENT(IN), OPTIONAL     :: remove_zero         ! [-] Whether to drop the zero-frequency bin
    COMPLEX(WP) , INTENT(OUT), ALLOCATABLE :: array_out(:,:,:)    ! [-] Scaled one-sided spectrum (Nfreqs, N2, N3)
    REAL(WP)    , INTENT(OUT), ALLOCATABLE :: freqs(:)            ! [Hz] Frequency vector (Nfreqs)

    ! Local variables
    INTEGER(I32)              :: nfreq_raw, nfreq_final, n_signals, i, j, k, idx, pos, start_idx
    INTEGER(I32)              :: skipf_, n2, n3
    LOGICAL                   :: remove_zero_
    REAL(WP)                  :: win_norm, mean_val, scale_fact
    REAL(WP)    , ALLOCATABLE :: window(:), data_flat(:,:)
    COMPLEX(WP) , ALLOCATABLE :: out_flat(:,:)

    ! FFTW variables
    TYPE(C_PTR)                 :: plan
    INTEGER(C_INT)              :: rank, howmany, istride, idist, ostride, odist, flags
    INTEGER(C_INT), ALLOCATABLE :: n_(:), inembed(:), onembed(:)


    ! Defaults
    skipf_       = 1_I32  ; if (present(skipf))       skipf_       = skipf
    remove_zero_ = .true. ; if (present(remove_zero)) remove_zero_ = remove_zero

    n2        = size(array, 2)
    n3        = size(array, 3)
    n_signals = n2 * n3
    nfreq_raw = nt / 2_I32 + 1_I32

    ! --- 1. Compute window and its norm once --- !
    ALLOCATE(window(nt))
    win_norm = 0.0_WP
    !$OMP PARALLEL DO REDUCTION(+:win_norm)
    do i = 1, nt
        window(i) = 0.5_WP - 0.5_WP * cos(2.0_WP * PI * real(i - 1_I32, WP) / real(nt - 1_I32, WP))
        win_norm  = win_norm + window(i)
    end do
    !$OMP END PARALLEL DO
    win_norm = win_norm / real(nt, WP)

    ! --- 2. Fused mean, windowing and flattening (no temporary 3D arrays) --- !
    ALLOCATE(data_flat(nt, n_signals))

    !$OMP PARALLEL DO PRIVATE(idx, j, k, i, mean_val)
    do idx = 1, n_signals
        ! Map flat index back to j, k to read directly from 3D array
        k = (idx - 1_I32) / n2 + 1_I32
        j = mod(idx - 1_I32, n2) + 1_I32

        ! Pass A: compute mean
        mean_val = 0.0_WP
        do i = 1, nt
            mean_val = mean_val + array(i, j, k)
        end do
        mean_val = mean_val / real(nt, WP)

        ! Pass B: subtract mean, apply window, store flat (column-major)
        do i = 1, nt
            data_flat(i, idx) = (array(i, j, k) - mean_val) * window(i)
        end do
    end do
    !$OMP END PARALLEL DO
    DEALLOCATE(window)

    ! --- 3. FFTW execution ---
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

    ! --- 4. Direct-to-final allocation (eliminates array re-allocations) --- !
    start_idx = 1_I32
    if (remove_zero_) start_idx = 2_I32

    ! Pre-calculate final exact size
    nfreq_final = 0_I32
    do i = start_idx, nfreq_raw, skipf_
        nfreq_final = nfreq_final + 1_I32
    end do

    ALLOCATE(array_out(nfreq_final, n2, n3))
    ALLOCATE(freqs(nfreq_final))

    ! Fill frequencies
    pos = 1_I32
    do i = start_idx, nfreq_raw, skipf_
        freqs(pos) = real(i - 1_I32, WP) / (real(nt, WP) * dt)
        pos        = pos + 1_I32
    end do

    ! Scale, conjugate to the exp(-i*w*t) convention and map directly to 3D array_out
    scale_fact = 2.0_WP / (real(nt, WP) * win_norm)

    !$OMP PARALLEL DO PRIVATE(idx, j, k, pos, i)
    do idx = 1, n_signals
        k = (idx - 1_I32) / n2 + 1_I32
        j = mod(idx - 1_I32, n2) + 1_I32

        pos = 1_I32
        ! Inner loop over frequencies guarantees contiguous writes to array_out
        do i = start_idx, nfreq_raw, skipf_
            array_out(pos, j, k) = conjg(out_flat(i, idx)) * scale_fact
            pos                  = pos + 1_I32
        end do
    end do
    !$OMP END PARALLEL DO

    DEALLOCATE(out_flat)

END SUBROUTINE compute_rfft

END MODULE SignalProcessing
