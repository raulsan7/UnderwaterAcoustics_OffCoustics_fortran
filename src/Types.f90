! ===============================================================================
! MODULE: Types
! ===============================================================================
! Defines the plain derived types shared across the whole off-coustics codebase:
! the configuration types filled in from the input namelist files (turbine,
! environment, solver, and the per-method solver settings), the abstract
! source representation passed from a wind turbine model into an acoustic
! solver, and the small helper types used internally by the Method of Images.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 2026-09-22
! ===============================================================================
MODULE TYPEs

USE Kinds

IMPLICIT NONE

PUBLIC


! ===============================================================================
! TYPE: Config_WindTurbine_t
! ===============================================================================
! Holds every wind-turbine-related setting read from the input namelist file.
! It is used to build the concrete WindTurbine_t object (Monopile, Floating,
! ...): where to find the OpenFAST outputs, how to read/filter/reduce the
! nodal time series, and where the turbine sits in the global reference frame.
! ===============================================================================
TYPE :: Config_WindTurbine_t

    ! --- Input parameters --- !
    CHARACTER(len=200) :: case_type             ! [-] Tag to distinguish turbine
    CHARACTER(len=200) :: rootname              ! [-] Name without extensions of the OpenFAST output files
    CHARACTER(len=200) :: OP_output_dir         ! [-] Directory with OpenFast output files
    CHARACTER(len=200) :: OP_output_path        ! [-] Openfast output file path
    CHARACTER(len=200) :: OP_SubDyn_path        ! [-] Openfast SubDyn file path
    LOGICAL            :: Binary                ! [-] Whether OpenFast files are binary (.outb) or not (.out)
    REAL(WP)           :: WindSpeed             ! [m/s] Wind Speed in norm
    REAL(WP)           :: WindDir               ! [deg] Wind direction 0 deg points to +x axis (anticlockwise from +x)
    REAL(WP)           :: Depth                 ! [m] Water depth
    INTEGER(I32)       :: Nmembers              ! [-] Number of structural members in the OpenFAST model
    INTEGER(I32)       :: Nnodes                ! [-] Number of structural nodes in the OpenFAST model
    REAL(WP)           :: AxisPos(2)            ! [m] Position of the turbine axis in xy plane
    REAL(WP)           :: BariPos(2)            ! [m] Position of the turbine baricenter in xy plane
    REAL(WP)           :: From                  ! [-] Fraction to start reading input
    REAL(WP)           :: Upto                  ! [-] Fraction to stop  reading input
    INTEGER(I32)       :: skip                  ! [-] Skips time snapshots as time(::skip)
    INTEGER(I32)       :: skipf                 ! [-] Skips time snapshots as freq(::skipf)
    LOGICAL            :: filter_freqs          ! [-] Whether to delete non useful frequencies
    REAL(WP)           :: freqs_over            ! [-] Frequencies smaller than freqs_over wont be filtered
    LOGICAL            :: is_monopile           ! [-] Whether the turbine is a monopile or a floating platform

END TYPE Config_WindTurbine_t


! ===============================================================================
! TYPE: Config_Environment_t
! ===============================================================================
! Holds the base environmental/fluid properties read from the input namelist
! file (site depth, speed of sound, density, reference pressure). These act
! as the default medium properties used by every acoustic solver, unless a
! solver-method-specific configuration overrides them.
! ===============================================================================
TYPE :: Config_Environment_t

    ! --- Input parameters --- !
    REAL(WP)           :: depth_base            ! [m] Site depth
    REAL(WP)           :: speed_sound_base      ! [m/s] Base speed of sound in fluid
    REAL(WP)           :: density_base          ! [kg/m3] Base fluid density
    REAL(WP)           :: reference_pressure    ! [Pa] Fluid reference pressure

END TYPE Config_Environment_t


! ===============================================================================
! TYPE: Config_Solver_t
! ===============================================================================
! Holds the generic acoustic-solver settings read from the input namelist
! file: which method to run, where to save the results, and the parallel
! execution parameters. Method-specific settings live in Config_Method_t
! and its extensions instead.
! ===============================================================================
TYPE :: Config_Solver_t

    ! --- Input parameters --- !
    LOGICAL            :: debug                 ! [-] More info for debugin
    LOGICAL            :: in_cluster            ! [-] Whether we are working on a cluster
    CHARACTER(len=200) :: save_path             ! [-] Where to save results
    CHARACTER(len=50)  :: method                ! [-] Acoustic method to use
    INTEGER(I32)       :: block_size            ! [-] Block size for each CPU in parallel computations

END TYPE Config_Solver_t


! ===============================================================================
! TYPE: Config_t
! ===============================================================================
! Top-level configuration object built by Args_parser from the input
! namelist files. It aggregates the turbine, environment and solver
! configurations together with the polymorphic, method-specific
! configuration (e.g. Config_ImageMethod_t) selected at run time.
! ===============================================================================
TYPE :: Config_t

    TYPE(Config_WindTurbine_t)          :: turbine    ! [-] Wind turbine configuration
    TYPE(Config_Environment_t)          :: env         ! [-] Environment/fluid configuration
    TYPE(Config_Solver_t)               :: solver      ! [-] Acoustic solver configuration
    CLASS(Config_Method_t), ALLOCATABLE :: method      ! [-] Polymorphic, solver-method-specific configuration

END TYPE Config_t


! ===============================================================================
! TYPE: AcousticSource_t
! ===============================================================================
! Abstract representation of the acoustic excitation produced by a wind
! turbine: a set of point sources with their coordinates, complex nodal
! forces per frequency, and the associated frequency vector. This is the
! hand-off object between WindTurbine_t % turbine_to_source and the
! acoustic solvers.
! ===============================================================================
TYPE:: AcousticSource_t

    REAL(WP)   , ALLOCATABLE :: xyz(:,:)            ! [m] Acoustic source coordinates (Nsources, 3)
    COMPLEX(WP), ALLOCATABLE :: Force(:,:,:)        ! [N] Complex source forces (Nfreqs, Nsources, 3)
    REAL(WP)   , ALLOCATABLE :: Freqs(:)            ! [Hz] Frequency verctor (Nfreqs)

END TYPE AcousticSource_t


! ===============================================================================
! ABSTRACT: Config_Method_t
! ===============================================================================
! Empty abstract base type used purely as a common label for every
! solver-method-specific configuration type (e.g. Config_ImageMethod_t).
! It lets Config_t % method be declared polymorphic and allocated at run
! time to whichever concrete configuration matches the selected method.
! ===============================================================================
TYPE, ABSTRACT :: Config_Method_t
END TYPE Config_Method_t


! ===============================================================================
! TYPE: Config_ImageMethod_t
! ===============================================================================
! Method-of-Images solver configuration, read from the ImagesMethod.nml
! namelist. Sets how many image sources are used to enforce the boundary
! conditions at the free surface and the seabed, the boundary types and
! heights, and the per-boundary reflection attenuation coefficients.
! ===============================================================================
TYPE, EXTENDS(Config_Method_t) :: Config_ImageMethod_t

    INTEGER(I32) :: N_images               ! [-] Number of image source levels used to enforce boundary conditions
    REAL(WP)     :: speed_sound            ! [m/s] Speed of sound override for the image method (reserved, currently unused)
    REAL(WP)     :: density                ! [kg/m3] Fluid density override for the image method (reserved, currently unused)
    INTEGER(I32) :: UpperBC                ! [-] Upper boundary condition flag (-1: pressure-release surface, p = 0)
    INTEGER(I32) :: LowerBC                ! [-] Lower boundary condition flag (+1: rigid boundary, dp/dn = 0)
    REAL(WP)     :: UpperHBC               ! [m] z-coordinate of the upper reflecting plane (e.g. free surface)
    REAL(WP)     :: LowerHBC               ! [m] z-coordinate of the lower reflecting plane (e.g. seabed)
    REAL(WP)     :: attenuation_surface    ! [-] Attenuation coefficient applied to surface reflections
    REAL(WP)     :: attenuation_bottom     ! [-] Attenuation coefficient applied to bottom reflections
    REAL(WP)     :: epsilon                ! [-] Small numerical tolerance used to avoid singularities in image computations

END TYPE Config_ImageMethod_t


! ===============================================================================
! TYPE: ImageSystem_t
! ===============================================================================
! Helper type used internally by the Method of Images to hold the full set
! of real and image sources built from the original AcousticSource_t: their
! coordinates, complex forces, and which boundary condition each image was
! generated from.
! ===============================================================================
TYPE :: ImageSystem_t

    REAL(WP)    , ALLOCATABLE :: xyz(:,:)      ! [m] Nodes coordinates (Nnodes, 3)
    COMPLEX(WP) , ALLOCATABLE :: Force(:,:,:)  ! [N] Force array (Nfreqs, Nnodes, 3)
    INTEGER(I32), ALLOCATABLE :: BC_all(:)     ! [-] Boundary conditions array

END TYPE ImageSystem_t


END MODULE TYPEs
