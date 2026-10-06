! ===============================================================================
! MODULE: Kinds
! ===============================================================================
! Defines the numeric kind parameters (integer and real precisions) and the
! fundamental physical/mathematical constants used consistently across the
! whole off-coustics codebase. The working precision (WP) must stay DP: the
! FFTW3 bindings (C_DOUBLE) and the HDF5 writers (H5T_NATIVE_DOUBLE) assume
! double precision.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE Kinds

IMPLICIT NONE

PUBLIC


! ===============================================================================
! INTEGER PRECISION KINDS
! ===============================================================================
INTEGER, PARAMETER :: I32 = selected_int_kind(9)     ! [-] 32-bit integer kind
INTEGER, PARAMETER :: I64 = selected_int_kind(18)    ! [-] 64-bit integer kind


! ===============================================================================
! REAL PRECISION KINDS
! ===============================================================================
INTEGER, PARAMETER :: SP = selected_real_kind(6, 37)      ! [-] Single precision kind
INTEGER, PARAMETER :: DP = selected_real_kind(15, 307)    ! [-] Double precision kind
INTEGER, PARAMETER :: WP = DP                             ! [-] Working precision (must stay DP, see header)


! ===============================================================================
! PHYSICAL AND MATHEMATICAL CONSTANTS
! ===============================================================================
REAL(WP)   , PARAMETER :: PI = 3.14159265358979323846_WP    ! [-] Pi constant
COMPLEX(WP), PARAMETER :: I1 = (0.0_WP, 1.0_WP)             ! [-] Imaginary unit


! ===============================================================================
! ARGUMENTS PARSER PARAMETERS
! ===============================================================================
REAL(WP)    , PARAMETER :: DEFAULT_R      = -999.0_WP    ! [-] Sentinel for an unset REAL parameter
INTEGER(I32), PARAMETER :: DEFAULT_I      = -999_I32     ! [-] Sentinel for an unset INTEGER parameter
INTEGER(I32), PARAMETER :: MAX_PARSE_OBS  = 200_I32      ! [-] Internal buffer size for the input observers
INTEGER(I32), PARAMETER :: MAX_PARSE_SRC  = 20_I32       ! [-] Internal buffer size for the input sources
INTEGER(I32), PARAMETER :: MAX_PARSE_FREQ = 200_I32      ! [-] Internal buffer size for the input frequencies

END MODULE Kinds
