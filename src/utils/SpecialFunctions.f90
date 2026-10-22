! ===============================================================================
! MODULE: SpecialFunctions
! ===============================================================================
! Defines the special functions of the normal-mode propagator: the Hankel
! functions (compute_hankel_complex), whose evanescent branch uses the K0/K1
! modified Bessel approximations (bessel_k0, bessel_k1).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE SpecialFunctions

USE Kinds, ONLY: I32, WP, PI, I1

IMPLICIT NONE

PRIVATE
PUBLIC :: compute_hankel_complex


CONTAINS


SUBROUTINE compute_hankel_complex(k_rm, r, H0, H1)
    ! -------------------------------------------------------------------------
    ! Efficiently computes the Hankel functions of the first kind (H0 and H1)
    ! for a wavenumber that is either purely real (propagative mode, using the
    ! intrinsic Bessel functions) or purely imaginary (evanescent mode, using
    ! the bessel_k0/bessel_k1 polynomial approximations below). H0 and H1 are
    ! optional and only the ones present are evaluated: the monopole kernel
    ! asks for H0 alone and saves the Bessel functions of order 1, the dipole
    ! kernel asks for both and gets the same values as before.
    ! -------------------------------------------------------------------------

    COMPLEX(WP), INTENT(IN)            :: k_rm    ! [1/m] Horizontal wavenumber (purely real or purely imaginary)
    REAL(WP)   , INTENT(IN)            :: r       ! [m] Horizontal distance
    COMPLEX(WP), INTENT(OUT), OPTIONAL :: H0      ! [-] Hankel function of the first kind, order 0, at k_rm * r
    COMPLEX(WP), INTENT(OUT), OPTIONAL :: H1      ! [-] Hankel function of the first kind, order 1, at k_rm * r

    ! Local variables
    REAL(WP) :: Z_real, Y


    if (abs(aimag(k_rm)) <= 1.0e-12_WP) then
        ! Propagative mode: k_rm is real
        Z_real = real(k_rm, WP) * r
        if (present(H0)) H0 = cmplx(bessel_jn(0, Z_real), bessel_yn(0, Z_real), kind=WP)
        if (present(H1)) H1 = cmplx(bessel_jn(1, Z_real), bessel_yn(1, Z_real), kind=WP)
    else
        ! Evanescent mode: k_rm is purely imaginary
        Y = abs(aimag(k_rm)) * r
        if (present(H0)) H0 = cmplx(0.0_WP, -2.0_WP / PI * bessel_k0(Y), kind=WP)
        if (present(H1)) H1 = cmplx(-2.0_WP / PI * bessel_k1(Y), 0.0_WP, kind=WP)
    end if

END SUBROUTINE compute_hankel_complex


PURE FUNCTION bessel_k0(x) RESULT(k0)
    ! -------------------------------------------------------------------------
    ! Abramowitz and Stegun polynomial approximation (9.8.5 - 9.8.6) of the
    ! modified Bessel function of the second kind, order 0: K0(x) for x > 0
    ! (relative error below 2e-7). Used by compute_hankel_complex to build the
    ! evanescent-mode Hankel function of the normal-mode propagator.
    ! -------------------------------------------------------------------------

    REAL(WP), INTENT(IN) :: x     ! [-] Argument (x > 0)
    REAL(WP)             :: k0    ! [-] K0(x)

    ! Local variables
    REAL(WP) :: y, t, i0


    if (x <= 0.0_WP) then
        k0 = 0.0_WP
        return
    end if

    if (x <= 2.0_WP) then
        y  = (x / 3.75_WP)**2
        i0 = 1.0_WP + y * (3.5156229_WP + y * (3.0899424_WP + y * (1.2067492_WP + &
             y * (0.2659732_WP + y * (0.0360768_WP + y * 0.0045813_WP)))))

        t  = x / 2.0_WP
        y  = t * t
        k0 = -log(t) * i0 + (-0.57721566_WP + y * (0.42278420_WP + y * (0.23069756_WP + &
             y * (0.03488590_WP + y * (0.00262698_WP + y * (0.00010750_WP + &
             y * 0.00000740_WP))))))
    else
        y  = 2.0_WP / x
        k0 = (exp(-x) / sqrt(x)) * (1.25331414_WP + y * (-0.07832358_WP + &
             y * (0.02189568_WP + y * (-0.01062449_WP + y * (0.00587872_WP + &
             y * (-0.00251540_WP + y * 0.00053208_WP))))))
    end if

END FUNCTION bessel_k0


PURE FUNCTION bessel_k1(x) RESULT(k1)
    ! -------------------------------------------------------------------------
    ! Abramowitz and Stegun polynomial approximation (9.8.9 - 9.8.10) of the
    ! modified Bessel function of the second kind, order 1: K1(x) for x > 0
    ! (relative error below 2e-7). Used by compute_hankel_complex to build the
    ! evanescent-mode Hankel function of the normal-mode propagator.
    ! -------------------------------------------------------------------------

    REAL(WP), INTENT(IN) :: x     ! [-] Argument (x > 0)
    REAL(WP)             :: k1    ! [-] K1(x)

    ! Local variables
    REAL(WP) :: y, t, i11


    if (x <= 0.0_WP) then
        k1 = 0.0_WP
        return
    end if

    if (x <= 2.0_WP) then
        y  = (x / 3.75_WP)**2
        i11 = x * (0.5_WP + y * (0.87890594_WP + y * (0.51498869_WP + &
             y * (0.15084934_WP + y * (0.02658733_WP + y * (0.00301532_WP + &
             y * 0.00032411_WP))))))

        t  = x / 2.0_WP
        y  = t * t
        k1 = log(t) * i11 + (1.0_WP / x) * (1.0_WP + y * (0.15443144_WP + &
             y * (-0.67278579_WP + y * (-0.18156897_WP + y * (-0.01919402_WP + &
             y * (-0.00110404_WP + y * (-0.00004686_WP)))))))
    else
        y  = 2.0_WP / x
        k1 = (exp(-x) / sqrt(x)) * (1.25331414_WP + y * (0.23498619_WP + &
             y * (-0.03655620_WP + y * (0.01504268_WP + y * (-0.00780353_WP + &
             y * (0.00325614_WP + y * (-0.00068245_WP)))))))
    end if

END FUNCTION bessel_k1

END MODULE SpecialFunctions
