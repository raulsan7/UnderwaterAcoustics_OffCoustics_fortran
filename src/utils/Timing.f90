! ===============================================================================
! MODULE: Timing
! ===============================================================================
! Defines the formatting of a wall-clock elapsed time into a human-readable
! value/unit pair (format_elapsed).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE Timing

USE omp_lib
USE Kinds, ONLY: WP

IMPLICIT NONE

PRIVATE
PUBLIC :: format_elapsed


CONTAINS


SUBROUTINE format_elapsed(start_time, elapsed_display, tag)
    ! -------------------------------------------------------------------------
    ! Converts a wall-clock elapsed time (in seconds, measured via OpenMP's
    ! timer) into a human-readable value/unit pair, automatically switching
    ! among seconds, minutes, hours and days depending on magnitude.
    ! -------------------------------------------------------------------------

    REAL(WP)        , INTENT(IN)  :: start_time         ! [s] Wall-clock start time (from omp_get_wtime)
    REAL(WP)        , INTENT(OUT) :: elapsed_display    ! [-] Elapsed time, expressed in the unit given by tag
    CHARACTER(len=*), INTENT(OUT) :: tag                ! [-] Unit of elapsed_display ("s", "min", "h" or "days")

    ! Local variables
    REAL(WP) :: end_time, elapsed_time


    ! Use the OpenMP wall-clock timer instead of the processor CPU_TIME
    end_time        = real(omp_get_wtime(), WP)
    elapsed_time    = end_time - start_time
    elapsed_display = elapsed_time

    if (elapsed_display > 60.0_WP .and. elapsed_display <= 3600.0_WP) then
        elapsed_display = elapsed_display / 60.0_WP
        tag = "min"
    else if (elapsed_display > 3600.0_WP .and. elapsed_display <= 86400.0_WP) then
        elapsed_display = elapsed_display / 3600.0_WP
        tag = "h"
    else if (elapsed_display > 86400.0_WP) then
        elapsed_display = elapsed_display / 86400.0_WP
        tag = "days"
    else
        tag = "s"
    end if

END SUBROUTINE format_elapsed

END MODULE Timing
