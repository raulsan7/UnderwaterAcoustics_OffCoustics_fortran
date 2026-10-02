! ===============================================================================
! MODULE: StructureModel
! ===============================================================================
! Defines the geometry and hydrodynamic corrections of the structural
! members: the removal of the nodes shared by consecutive members
! (remove_duplicate_nodes), the span represented by every node
! (divide_span), the submerged area of an elliptical cross-section
! (integral_elipse_z), and the cylindrical scattering-correction
! coefficient of the radiated force (alpha_hankel).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE StructureModel

USE Kinds, ONLY: I32, WP, PI

IMPLICIT NONE

PRIVATE
PUBLIC :: remove_duplicate_nodes, divide_span, integral_elipse_z, alpha_hankel


CONTAINS


FUNCTION remove_duplicate_nodes(x_all, delete_last) RESULT(mask)
    ! -------------------------------------------------------------------------
    ! Removes duplicate nodes shared between consecutive structural members.
    ! Given an array of nodes of shape (Nm, Nn, 3), where Nm is the number of
    ! members and Nn the number of nodes per member, identifies nodes that are
    ! duplicated at the junction between member i and member i+1 (the last
    ! node of member i is assumed to coincide with the first node of member
    ! i+1) and marks one of the two for removal, based on delete_last.
    ! -------------------------------------------------------------------------

    REAL(WP), INTENT(IN)  :: x_all(:,:,:)    ! [m] Nodes coordinates (Nm, Nn, 3)
    LOGICAL , INTENT(IN)  :: delete_last     ! [-] Drop the next member's first node (.true.) or the last one (.false.)
    LOGICAL , ALLOCATABLE :: mask(:)         ! [-] Flat mask (Nm*Nn), .true. for nodes to keep

    ! Local variables
    INTEGER(I32) :: Nm, Nn, i, idx_prev, idx_next
    REAL(WP)     :: tol
    LOGICAL      :: are_equal


    Nm = size(x_all, 1)
    Nn = size(x_all, 2)

    ALLOCATE(mask(Nm*Nn))
    mask = .true.    ! Start keeping all nodes

    ! Tolerance for comparing coordinates
    tol = 1.0e-6_WP

    do i = 1, Nm - 1
        ! Last node of current member
        idx_prev = (i - 1_I32) * Nn + Nn
        ! First node of next member
        idx_next = i * Nn + 1_I32

        ! Compare coordinates component-wise within the tolerance
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


FUNCTION divide_span(x) RESULT(dx)
    ! -------------------------------------------------------------------------
    ! Computes representative span intervals (dx) for an ordered spatial vector.
    ! Fully vectorized to eliminate conditionals and temporary arrays.
    ! -------------------------------------------------------------------------

    REAL(WP), INTENT(IN)  :: x(:)     ! [m] Ordered (monotonic) spatial coordinates
    REAL(WP), ALLOCATABLE :: dx(:)    ! [m] Span interval represented by each point

    ! Local variables
    INTEGER(I32)              :: n
    REAL(WP)    , ALLOCATABLE :: x_work(:)
    LOGICAL                   :: flipped


    n = size(x)
    ALLOCATE(dx(n))

    if (n == 0_I32) return

    ALLOCATE(x_work(n))

    ! Check direction and possibly reverse
    flipped = .false.
    if (n > 1_I32 .and. x(1) > x(n)) then
        x_work  = x(n:1:-1)
        flipped = .true.
    else
        x_work = x
    end if

    ! Vectorized span calculation via algebraic simplification
    if (n == 1_I32) then
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


PURE FUNCTION integral_elipse_z(Dy, Dz, k) RESULT(area)
    ! -------------------------------------------------------------------------
    ! Area under the horizontal line z = k inside the ellipse
    ! (2y/Dy)^2 + (2z/Dz)^2 <= 1, i.e. the submerged cross-sectional area of an
    ! elliptical member (semi-axes Dy/2, Dz/2) up to elevation k. Used to build
    ! the outer/inner submerged shell area of the SAITEC elliptical floaters.
    ! -------------------------------------------------------------------------

    REAL(WP), INTENT(IN) :: Dy      ! [m] Diameter in y-direction
    REAL(WP), INTENT(IN) :: Dz      ! [m] Diameter in z-direction
    REAL(WP), INTENT(IN) :: k       ! [m] z(y) = k line
    REAL(WP)             :: area    ! [m^2] Area of the ellipse below the line z = k

    ! Local variables
    REAL(WP) :: a, b, elipse_area, k_clip, u


    a = Dy / 2.0_WP
    b = Dz / 2.0_WP

    elipse_area = PI * a * b

    k_clip = min(max(k, -b), b)
    u      = k_clip / b

    area = 0.5_WP * elipse_area + a * b * (asin(u) + u * sqrt(1.0_WP - u * u))

END FUNCTION integral_elipse_z


FUNCTION alpha_hankel(k, D) RESULT(alpha)
    ! -------------------------------------------------------------------------
    ! Analytical scattering-correction coefficient alpha, based on cylindrical
    ! Hankel functions (orders 0 and 2, combined into the derivative of the
    ! order-1 function). Used to correct the radiated acoustic force for the
    ! finite diameter of a cylindrical structural member.
    ! -------------------------------------------------------------------------

    USE Kinds, ONLY: I1

    REAL(WP)   , INTENT(IN)  :: k(:)        ! [1/m] Wavenumber array (Nfreqs)
    REAL(WP)   , INTENT(IN)  :: D           ! [m] Cylinder diameter
    COMPLEX(WP), ALLOCATABLE :: alpha(:)    ! [-] Scattering-correction coefficient (Nfreqs)

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: H0(:), H2(:), H1p(:)
    REAL(WP)    , ALLOCATABLE :: ka(:)
    INTEGER(I32)              :: Nfreqs


    Nfreqs = size(k)
    ka     = k * D / 2.0_WP

    ALLOCATE(H0(Nfreqs), H2(Nfreqs), H1p(Nfreqs), alpha(Nfreqs))

    H0  = cmplx(bessel_jn(0, ka), -bessel_yn(0, ka), kind=WP)    ! H_0^(2), order 0, second kind (H_n^(2) = Jn - iYn)
    H2  = cmplx(bessel_jn(2, ka), -bessel_yn(2, ka), kind=WP)    ! H_2^(2), order 2, second kind (H_n^(2) = Jn - iYn)
    H1p = 0.5_WP * (H0 - H2)                                     ! Derivative of H_1^(2): (H_0^(2) - H_2^(2)) / 2

    alpha = -3.0_WP * I1 / (PI * ka**2 * H1p)

END FUNCTION alpha_hankel

END MODULE StructureModel
