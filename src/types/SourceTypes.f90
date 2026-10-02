! ===============================================================================
! MODULE: SourceTypes
! ===============================================================================
! Defines the source representation handed from an emitter (a wind turbine,
! a wind farm or hand-made sources) to an acoustic solver: the point sources
! with their coordinates, complex forces, frequency vector and multipole
! order (AcousticSource_t, labelled with SOURCE_MONOPOLE, SOURCE_DIPOLE, ...),
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
PUBLIC :: GeometricData_t, AcousticSource_t, SOURCE_MONOPOLE, SOURCE_DIPOLE


! ===============================================================================
! MULTIPOLE ORDER LABELS
! A label per source, so that higher orders can be added later without new arrays
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


TYPE :: AcousticSource_t
    ! -------------------------------------------------------------------------
    ! Abstract representation of the acoustic excitation produced by a wind
    ! turbine: a set of point sources with their coordinates, complex nodal
    ! forces per frequency, the associated frequency vector and the multipole
    ! order of every source. This is the hand-off object between
    ! Emitter_t % to_source and the acoustic solvers.
    ! -------------------------------------------------------------------------

    REAL(WP)    , ALLOCATABLE :: xyz(:,:)        ! [m] Acoustic source coordinates (Nsources, 3)
    COMPLEX(WP) , ALLOCATABLE :: Force(:,:,:)    ! [N] Complex source forces (Nfreqs, Nsources, 3)
    REAL(WP)    , ALLOCATABLE :: Freqs(:)        ! [Hz] Frequency vector (Nfreqs)
    INTEGER(I32), ALLOCATABLE :: order(:)        ! [-] Multipole order per source (SOURCE_MONOPOLE, SOURCE_DIPOLE, ...)
    TYPE(GeometricData_t)     :: geometry        ! [-] Reference geometry and labels of the emitter

END TYPE AcousticSource_t

END MODULE SourceTypes
