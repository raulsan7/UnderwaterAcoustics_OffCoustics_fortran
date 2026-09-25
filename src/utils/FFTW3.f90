! ===============================================================================
! MODULE: FFTW3
! ===============================================================================
! Thin ISO_C_BINDING interface to the FFTW3 C library, exposing only the
! subset needed by this project: planning and executing many-array real-to-
! complex DFTs, destroying the resulting plan, and the FFTW planning-strategy
! constants (ESTIMATE/MEASURE/PATIENT/EXHAUSTIVE). No Fortran-side logic
! lives here; it is a direct pass-through to the underlying C API.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE FFTW3

USE, INTRINSIC :: iso_c_binding, ONLY: C_INT, C_PTR, C_DOUBLE, C_DOUBLE_COMPLEX

IMPLICIT NONE

PRIVATE
PUBLIC :: FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE
PUBLIC :: fftw_plan_many_dft_r2c, fftw_execute_dft_r2c, fftw_destroy_plan


! ===============================================================================
! FFTW PLANNING-STRATEGY CONSTANTS
! ===============================================================================
! Mirror the values defined in fftw3.h
INTEGER(C_INT), PARAMETER :: FFTW_ESTIMATE   = 64    ! [-] Quick heuristic plan
INTEGER(C_INT), PARAMETER :: FFTW_MEASURE    = 0     ! [-] Plan chosen by timing several candidates
INTEGER(C_INT), PARAMETER :: FFTW_PATIENT    = 32    ! [-] Wider candidate search than FFTW_MEASURE
INTEGER(C_INT), PARAMETER :: FFTW_EXHAUSTIVE = 8     ! [-] Exhaustive candidate search


! ===============================================================================
! FFTW C BINDINGS
! ===============================================================================
! See http://fftw.org for the underlying C API semantics
INTERFACE

    TYPE(C_PTR) FUNCTION fftw_plan_many_dft_r2c(rank, n, howmany, &
                                                in, inembed, istride, idist, &
                                                out, onembed, ostride, odist, &
                                                flags) BIND(C, name='fftw_plan_many_dft_r2c')
        ! -------------------------------------------------------------------------
        ! Creates (plans) a many-array real-to-complex DFT transform.
        ! -------------------------------------------------------------------------

        IMPORT :: C_INT, C_PTR, C_DOUBLE, C_DOUBLE_COMPLEX

        INTEGER(C_INT)           , VALUE        :: rank       ! [-] Rank (dimensionality) of each transform
        INTEGER(C_INT)           , DIMENSION(*) :: n          ! [-] Size of each transform dimension
        INTEGER(C_INT)           , VALUE        :: howmany    ! [-] Number of transforms to compute
        REAL(C_DOUBLE)           , DIMENSION(*) :: in         ! [-] Real input data
        INTEGER(C_INT)           , DIMENSION(*) :: inembed    ! [-] Embedding dimensions of the input array
        INTEGER(C_INT)           , VALUE        :: istride    ! [-] Stride between consecutive input elements
        INTEGER(C_INT)           , VALUE        :: idist      ! [-] Distance between consecutive input arrays
        COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(*) :: out        ! [-] Complex output data
        INTEGER(C_INT)           , DIMENSION(*) :: onembed    ! [-] Embedding dimensions of the output array
        INTEGER(C_INT)           , VALUE        :: ostride    ! [-] Stride between consecutive output elements
        INTEGER(C_INT)           , VALUE        :: odist      ! [-] Distance between consecutive output arrays
        INTEGER(C_INT)           , VALUE        :: flags      ! [-] Planning-strategy flags (FFTW_ESTIMATE, ...)

    END FUNCTION fftw_plan_many_dft_r2c

    SUBROUTINE fftw_execute_dft_r2c(plan, in, out) BIND(C, name='fftw_execute_dft_r2c')
        ! -------------------------------------------------------------------------
        ! Executes a previously created real-to-complex DFT plan on given data.
        ! -------------------------------------------------------------------------

        IMPORT :: C_PTR, C_DOUBLE, C_DOUBLE_COMPLEX

        TYPE(C_PTR)              , VALUE        :: plan    ! [-] Plan created by fftw_plan_many_dft_r2c
        REAL(C_DOUBLE)           , DIMENSION(*) :: in      ! [-] Real input data
        COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(*) :: out     ! [-] Complex output data

    END SUBROUTINE fftw_execute_dft_r2c

    SUBROUTINE fftw_destroy_plan(plan) BIND(C, name='fftw_destroy_plan')
        ! -------------------------------------------------------------------------
        ! Frees the memory associated with a plan, once it is no longer needed.
        ! -------------------------------------------------------------------------

        IMPORT :: C_PTR

        TYPE(C_PTR), VALUE :: plan    ! [-] Plan to destroy

    END SUBROUTINE fftw_destroy_plan

END INTERFACE

END MODULE FFTW3
