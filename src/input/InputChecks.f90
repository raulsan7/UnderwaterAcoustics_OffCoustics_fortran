! ===============================================================================
! MODULE: InputChecks
! ===============================================================================
! Defines the consistency checks shared by the namelist readers: the
! OpenFAST reading window and decimation factors (check_time_window), the
! length of a string against its Config_t component (check_fits), vectors
! given in full or not at all (check_all_or_none), and observer counts
! (check_count). Each one stops the run when its check fails.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE InputChecks

USE Kinds, ONLY: I32, WP, DEFAULT_R, DEFAULT_I

IMPLICIT NONE

PRIVATE
PUBLIC :: check_time_window, check_fits, check_all_or_none, check_count, check_source_positions


CONTAINS


SUBROUTINE check_time_window(From, Upto, skip, skipf, block_name)
    ! -------------------------------------------------------------------------
    ! Stops the run when the OpenFAST reading window or the decimation
    ! factors cannot be used: 0 <= From < Upto <= 1 and skip, skipf >= 1 are
    ! required. Otherwise the reader divides by zero (skip, skipf = 0), reads
    ! an empty window (From >= Upto) or rows that do not exist (From < 0,
    ! Upto > 1), or the FFT keeps no frequency at all (skipf < 0).
    ! -------------------------------------------------------------------------

    REAL(WP)        , INTENT(IN) :: From          ! [-] Fraction of the time series to start reading from
    REAL(WP)        , INTENT(IN) :: Upto          ! [-] Fraction of the time series to stop reading at
    INTEGER(I32)    , INTENT(IN) :: skip          ! [-] Keep 1 out of every skip time samples
    INTEGER(I32)    , INTENT(IN) :: skipf         ! [-] Keep 1 out of every skipf frequency bins
    CHARACTER(len=*), INTENT(IN) :: block_name    ! [-] Namelist group shown in the error message


    ! Written with .not. so that a NaN value is rejected as well
    if (.not. (From >= 0.0_WP .and. From < Upto .and. Upto <= 1.0_WP)) then
        write(*, '(A,F0.4,A,F0.4,A)') " ERROR: &" // block_name // ": From = ", From, ", Upto = ", Upto, &
                                      " (0 <= From < Upto <= 1 is required)"
        error stop
    end if

    if (skip < 1_I32 .or. skipf < 1_I32) then
        write(*, '(A,I0,A,I0,A)') " ERROR: &" // block_name // ": skip = ", skip, ", skipf = ", skipf, &
                                  " (both must be >= 1)"
        error stop
    end if

END SUBROUTINE check_time_window


SUBROUTINE check_fits(value, component, name)
    ! -------------------------------------------------------------------------
    ! Stops the run when a string read from a namelist is longer than the
    ! Config_t component that will store it, which would silently truncate
    ! it (e.g. an output path cut at 200 characters). Only the length of
    ! component is used, never its value.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: value        ! [-] String as read from the namelist
    CHARACTER(len=*), INTENT(IN) :: component    ! [-] Config_t component that will store it
    CHARACTER(len=*), INTENT(IN) :: name         ! [-] Namelist variable shown in the error message


    if (len_trim(value) > len(component)) then
        write(*, '(A,I0,A,I0,A)') " ERROR: " // name // " has ", len_trim(value), &
                                  " characters, the maximum is ", len(component), ": " // trim(value)
        error stop
    end if

END SUBROUTINE check_fits


SUBROUTINE check_all_or_none(values, name)
    ! -------------------------------------------------------------------------
    ! Stops the run when a vector parameter whose -999.0 sentinel applies to
    ! the whole vector (a centre, x/y limits, decay distances) is only partly
    ! given. The run_* routines would otherwise either use the missing
    ! component literally as -999.0 m or silently drop the given ones.
    ! -------------------------------------------------------------------------

    REAL(WP)        , INTENT(IN) :: values(:)    ! [m] Components as read from the namelist
    CHARACTER(len=*), INTENT(IN) :: name         ! [-] Namelist variable shown in the error message


    if (any(values == DEFAULT_R) .and. .not. all(values == DEFAULT_R)) then
        write(*, '(A)') " ERROR: " // name // " is only partly given: set every component, " // &
                        "or leave all of them at -999.0 (auto)"
        error stop
    end if

END SUBROUTINE check_all_or_none


SUBROUTINE check_count(n, name)
    ! -------------------------------------------------------------------------
    ! Stops the run when an observer count is neither -999 (auto) nor at
    ! least 1. A zero or negative count gives empty datasets (or a NaN
    ! spacing) and a run that ends "successfully" without any observer.
    ! -------------------------------------------------------------------------

    INTEGER(I32)    , INTENT(IN) :: n       ! [-] Count as read from the namelist
    CHARACTER(len=*), INTENT(IN) :: name    ! [-] Namelist variable shown in the error message


    if (n /= DEFAULT_I .and. n < 1_I32) then
        write(*, '(A,I0,A)') " ERROR: " // name // " = ", n, " (use -999 for the automatic value, or at least 1)"
        error stop
    end if

END SUBROUTINE check_count


SUBROUTINE check_source_positions(source_xyz, Nsources, block_name)
    ! -------------------------------------------------------------------------
    ! Stops the run when the positions of the point sources of a namelist
    ! cannot be used: rows given beyond Nsources, a source with a missing or
    ! non-finite coordinate, or a source above the free surface (z > 0). The
    ! seabed is not checked here, since it can depend on (x, y) and on the
    ! method: the solver checks that every source is inside the water column
    ! (check_sources_in_water). The buffer is initialised to NaN, so NaN marks
    ! an entry not given.
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite

    REAL(WP)        , INTENT(IN) :: source_xyz(:,:)    ! [m] Positions as read (rows of the buffer, 3)
    INTEGER(I32)    , INTENT(IN) :: Nsources           ! [-] Number of sources
    CHARACTER(len=*), INTENT(IN) :: block_name         ! [-] Namelist group shown in the error message

    ! Local variables
    INTEGER(I32) :: i_s


    if (any(ieee_is_finite(source_xyz(Nsources+1:, :)))) then
        write(*, '(A,I0,A)') " ERROR: &" // block_name // ": source_xyz has more than Nsources = ", Nsources, " rows"
        error stop
    end if

    do i_s = 1, Nsources
        if (.not. all(ieee_is_finite(source_xyz(i_s, :)))) then
            write(*, '(A,I0,A)') " ERROR: &" // block_name // ": source_xyz(", i_s, ", 1:3) is missing or not finite"
            error stop
        end if

        if (source_xyz(i_s, 3) > 0.0_WP) then
            write(*, '(A,I0,A,F0.2,A)') " ERROR: &" // block_name // ": source ", i_s, " is at z = ", &
                                        source_xyz(i_s, 3), " m, above the free surface (z = 0)"
            error stop
        end if
    end do

END SUBROUTINE check_source_positions

END MODULE InputChecks
