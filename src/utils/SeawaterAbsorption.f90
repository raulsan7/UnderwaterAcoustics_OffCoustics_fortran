! ===============================================================================
! MODULE: SeawaterAbsorption
! ===============================================================================
! Defines the volume absorption coefficient of sea water: the simplified
! formula of Ainslie and McColm (1998, J. Acoust. Soc. Am. 103, 1671-1672,
! Eqs. 1-3) for the viscous and chemical (boric acid and magnesium sulphate
! relaxation) absorption, within 10 % of the Francois-Garrison (1982)
! equation between 100 Hz and 1 MHz, and the decibel/neper conversion used
! to turn it into an attenuation rate of the pressure amplitude.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE SeawaterAbsorption

USE Kinds, ONLY: WP

IMPLICIT NONE

PRIVATE
PUBLIC :: absorption_ainslie_mccolm, DB_PER_NEPER


REAL(WP), PARAMETER :: DB_PER_NEPER = 20.0_WP / log(10.0_WP)    ! [-] Decibels per neper of pressure amplitude


CONTAINS


ELEMENTAL FUNCTION absorption_ainslie_mccolm(f, T, S, pH, z) RESULT(alpha)
    ! -------------------------------------------------------------------------
    ! Volume absorption coefficient of sea water [dB/km] from Ainslie and
    ! McColm (1998), Eqs. (1)-(3), with f in kHz and z in km as in the paper:
    !     f1    = 0.78 * sqrt(S/35) * exp(T/26)               (boric acid)
    !     f2    = 42 * exp(T/17)                               (magnesium sulphate)
    !     alpha = 0.106 * f1*f^2/(f^2 + f1^2) * exp((pH - 8)/0.56)
    !           + 0.52 * (1 + T/43) * (S/35) * f2*f^2/(f^2 + f2^2) * exp(-z/6)
    !           + 0.00049 * f^2 * exp(-(T/27 + z/17))
    ! The paper checks it against Francois-Garrison for 100 Hz - 1 MHz,
    ! -6 < T < 35 degC, 7.7 < pH < 8.3, 5 < S < 50 ppt and 0 < z < 7 km; below
    ! 100 Hz it is extrapolated (absorption below 1e-3 dB/km there).
    ! -------------------------------------------------------------------------

    REAL(WP), INTENT(IN) :: f        ! [Hz] Frequency
    REAL(WP), INTENT(IN) :: T        ! [degC] Water temperature
    REAL(WP), INTENT(IN) :: S        ! [ppt] Water salinity
    REAL(WP), INTENT(IN) :: pH       ! [-] Water acidity
    REAL(WP), INTENT(IN) :: z        ! [m] Depth (positive downwards)
    REAL(WP)             :: alpha    ! [dB/km] Absorption coefficient

    ! Local variables
    REAL(WP) :: f_sq, z_km, f1, f2


    f_sq = (f * 1.0e-3_WP)**2      ! [kHz^2]
    z_km = z * 1.0e-3_WP           ! [km]
    f1   = 0.78_WP * sqrt(S / 35.0_WP) * exp(T / 26.0_WP)
    f2   = 42.0_WP * exp(T / 17.0_WP)

    alpha = 0.106_WP * f1 * f_sq / (f_sq + f1**2) * exp((pH - 8.0_WP) / 0.56_WP)                     &
          + 0.52_WP * (1.0_WP + T / 43.0_WP) * (S / 35.0_WP) * f2 * f_sq / (f_sq + f2**2) * exp(-z_km / 6.0_WP) &
          + 0.00049_WP * f_sq * exp(-(T / 27.0_WP + z_km / 17.0_WP))

END FUNCTION absorption_ainslie_mccolm

END MODULE SeawaterAbsorption
