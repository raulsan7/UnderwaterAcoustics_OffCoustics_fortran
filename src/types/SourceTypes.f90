! ===============================================================================
! MODULE: SourceTypes
! ===============================================================================
! Defines the source representation handed from an emitter (a wind turbine,
! a wind farm or hand-made sources) to an acoustic solver: the point sources
! with their coordinates, complex forces, frequency vector and multipole
! order (AcousticSource_t, labelled with SOURCE_MONOPOLE, SOURCE_DIPOLE, ...),
! their grouping into vertical axes (AxisGroup_t, built by group_by_axes),
! and the reference geometry and output labels of the emitter
! (GeometricData_t).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE SourceTypes

USE Kinds, ONLY: I32, WP

IMPLICIT NONE

PRIVATE
PUBLIC :: GeometricData_t, AcousticSource_t, AxisGroup_t, SOURCE_MONOPOLE, SOURCE_DIPOLE


! ===============================================================================
! MULTIPOLE ORDER LABELS
! A label per source, so that higher orders can be added later without new arrays.
! The strength of a source is in Force(:, s, :) of AcousticSource_t:
! - dipole: the force applied on the fluid, components (Fx, Fy, Fz) [N];
! - monopole: the volume velocity Q [m^3/s] in component 1, components 2 and 3
!   are 0. Its pressure is p = -i*w*rho*Q * exp(i*k*r) / (4*pi*r).
! All of them are peak complex amplitudes with time factor exp(-i*w*t).
! ===============================================================================
INTEGER(I32), PARAMETER :: SOURCE_MONOPOLE = 0_I32    ! [-] Monopole source (order 0)
INTEGER(I32), PARAMETER :: SOURCE_DIPOLE   = 1_I32    ! [-] Dipole source (order 1)


TYPE :: GeometricData_t
    ! -------------------------------------------------------------------------
    ! Reference geometry and labels of whatever produced the acoustic sources
    ! (a turbine, a wind farm or hand-made sources). Used by the automatic
    ! defaults of the run_* routines and by the default output file, built as
    ! out_dir // name_tag // <solver method tag> // name_suffix // ".hdf5".
    ! -------------------------------------------------------------------------

    REAL(WP)           :: center(2)   = 0.0_WP                ! [m] Reference centre (auto: polar, cylinder, sphere)
    REAL(WP)           :: origin(2)   = 0.0_WP                ! [m] Reference origin (auto start of decay, y of sliceXZ)
    REAL(WP)           :: direction   = 0.0_WP                ! [deg] Reference direction, anticlockwise from +x
    CHARACTER(len=200) :: out_dir     = './acoustic_data/'    ! [-] Default output folder ('./farm_acoustic_data/')
    CHARACTER(len=128) :: name_tag    = 'plot_sources'        ! [-] File name head ('plot_mn', 'Farm_2_DTU10MN')
    CHARACTER(len=64)  :: name_suffix = ''                    ! [-] File name tail ('_ws11.4', '_ws8.0-11.4')

END TYPE GeometricData_t


TYPE :: AxisGroup_t
    ! -------------------------------------------------------------------------
    ! A vertical axis of sources: the sources that share the same horizontal
    ! position (x, y), within the tolerance of group_by_axes, at any depth
    ! and of any multipole order. Every source belongs to exactly one axis.
    ! Built by AcousticSource_t % group_by_axes; the normal-mode kernels
    ! evaluate the terms that depend only on the horizontal distance once
    ! per axis and reuse them for every source of it.
    ! -------------------------------------------------------------------------

    REAL(WP)                  :: x             ! [m] x of the axis (that of its first source)
    REAL(WP)                  :: y             ! [m] y of the axis (that of its first source)
    INTEGER(I32)              :: Nsources      ! [-] Number of sources on the axis
    INTEGER(I32)              :: Nmonopoles    ! [-] Number of monopoles on the axis
    INTEGER(I32)              :: Ndipoles      ! [-] Number of dipoles on the axis
    INTEGER(I32), ALLOCATABLE :: idx(:)        ! [-] Indices of its sources in xyz, Force, order; ascending (Nsources)

END TYPE AxisGroup_t


TYPE :: AcousticSource_t
    ! -------------------------------------------------------------------------
    ! Abstract representation of the acoustic excitation produced by an
    ! emitter (a wind turbine, a wind farm or hand-made sources): a set of
    ! point sources with their coordinates, complex strengths per frequency
    ! (forces of the dipoles, volume velocities of the monopoles; see the
    ! labels above), the associated frequency vector and the multipole order
    ! of every source. This is the hand-off object between
    ! Emitter_t % to_source and the acoustic solvers. The grouping of the
    ! sources into vertical axes (axis) is not filled in by the emitters: a
    ! solver that uses it builds it with group_by_axes.
    ! -------------------------------------------------------------------------

    REAL(WP)         , ALLOCATABLE :: xyz(:,:)        ! [m] Acoustic source coordinates (Nsources, 3)
    COMPLEX(WP)      , ALLOCATABLE :: Force(:,:,:)    ! [N] Forces; monopole: Q [m^3/s] in 1 (Nfreqs, Nsources, 3)
    REAL(WP)         , ALLOCATABLE :: Freqs(:)        ! [Hz] Frequency vector (Nfreqs)
    INTEGER(I32)     , ALLOCATABLE :: order(:)        ! [-] Multipole order per source (SOURCE_MONOPOLE, ...)
    TYPE(AxisGroup_t), ALLOCATABLE :: axis(:)         ! [-] Vertical axes of the sources, from group_by_axes (Naxes)
    TYPE(GeometricData_t)          :: geometry        ! [-] Reference geometry and labels of the emitter

    CONTAINS

    ! --- Public type-bound procedures --- !
    PROCEDURE :: group_by_axes

END TYPE AcousticSource_t


CONTAINS


SUBROUTINE group_by_axes(self, tolerance)
    ! -------------------------------------------------------------------------
    ! Groups the sources into vertical axes (self % axis): a source belongs
    ! to the first axis whose first source lies at a horizontal distance
    ! <= tolerance (1 cm by default; 0 asks for exactly the same x and y),
    ! or else opens a new axis. The axes are numbered in order of first
    ! appearance, the indices of each axis are ascending and each axis counts
    ! its monopoles and dipoles. The horizontal position of an axis is that
    ! of its first source, so a source within the tolerance but not exactly
    ! on it is evaluated at that position (an error of at most the tolerance
    ! in its horizontal distance). Any previous grouping is replaced. Used by
    ! the normal-mode kernels, whose Hankel functions depend only on the
    ! horizontal source-observer distance.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSource_t), INTENT(INOUT)        :: self         ! [-] Acoustic sources
    REAL(WP)               , INTENT(IN), OPTIONAL :: tolerance    ! [m] Maximum horizontal distance within an axis

    ! Local variables
    REAL(WP)                  :: tol
    INTEGER(I32)              :: i, a, Nsources, Naxes
    INTEGER(I32), ALLOCATABLE :: axis_of(:), first_of(:), src_ids(:)


    tol = 1.0e-2_WP ; if (present(tolerance)) tol = tolerance

    if (.not. ALLOCATED(self % xyz) .or. .not. ALLOCATED(self % order)) &
        error stop "AcousticSource_t % group_by_axes(): xyz and order must be allocated"

    Nsources = size(self % xyz, 1)
    ALLOCATE(axis_of(Nsources), first_of(Nsources))
    src_ids = [(i, i = 1, Nsources)]

    ! Axis of every source: the first axis within tol, or a new one
    Naxes = 0_I32
    do i = 1, Nsources
        axis_of(i) = 0_I32
        do a = 1, Naxes
            if (hypot(self % xyz(i, 1) - self % xyz(first_of(a), 1), &
                      self % xyz(i, 2) - self % xyz(first_of(a), 2)) <= tol) then
                axis_of(i) = a
                exit
            end if
        end do

        if (axis_of(i) == 0_I32) then
            Naxes           = Naxes + 1_I32
            first_of(Naxes) = i
            axis_of(i)      = Naxes
        end if
    end do

    ! Position, sources (ascending) and counts of every axis
    if (ALLOCATED(self % axis)) DEALLOCATE(self % axis)
    ALLOCATE(self % axis(Naxes))

    do a = 1, Naxes
        self % axis(a) % x          = self % xyz(first_of(a), 1)
        self % axis(a) % y          = self % xyz(first_of(a), 2)
        self % axis(a) % idx        = pack(src_ids, axis_of == a)
        self % axis(a) % Nsources   = size(self % axis(a) % idx)
        self % axis(a) % Nmonopoles = count(self % order(self % axis(a) % idx) == SOURCE_MONOPOLE)
        self % axis(a) % Ndipoles   = count(self % order(self % axis(a) % idx) == SOURCE_DIPOLE)
    end do

    DEALLOCATE(axis_of, first_of, src_ids)

END SUBROUTINE group_by_axes

END MODULE SourceTypes
