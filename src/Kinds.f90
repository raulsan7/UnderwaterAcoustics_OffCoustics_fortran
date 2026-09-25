! ===============================================================================
! MODULE: Kinds
! ===============================================================================
! Defines the numeric kind parameters (integer and real precisions) and the
! fundamental physical/mathematical constants used consistently across the
! whole off-coustics codebase. Centralizing them here means the working
! precision (WP) can be changed for the entire project from a single place.
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
INTEGER, PARAMETER :: WP = DP                             ! [-] Working precision, change when desired


! ===============================================================================
! PHYSICAL AND MATHEMATICAL CONSTANTS
! ===============================================================================
REAL(WP)   , PARAMETER :: PI = 3.14159265358979323846_WP    ! [-] Pi constant
COMPLEX(WP), PARAMETER :: I1 = (0.0_WP, 1.0_WP)             ! [-] Imaginary unit


! ===============================================================================
! ARGUMENTS PARSER PARAMETERS
! ===============================================================================
REAL(WP)    , PARAMETER :: DEFAULT_R     = -999.0_WP    ! [-] Sentinel for an unset REAL parameter
INTEGER(I32), PARAMETER :: DEFAULT_I     = -999_I32     ! [-] Sentinel for an unset INTEGER parameter
INTEGER(I32), PARAMETER :: MAX_PARSE_OBS = 200_I32      ! [-] Internal buffer size for the input observers

END MODULE Kinds
