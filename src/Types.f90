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
! Date    : 09/2026
! ===============================================================================
MODULE Types

USE Kinds

IMPLICIT NONE

PUBLIC


TYPE :: Config_WindTurbine_t
    ! -------------------------------------------------------------------------
    ! Holds every wind-turbine-related setting read from the input namelist file.
    ! It is used to build the concrete WindTurbine_t object (Monopile, Floating,
    ! ...): where to find the OpenFAST outputs, how to read/filter/reduce the
    ! nodal time series, and where the turbine sits in the global reference frame.
    ! -------------------------------------------------------------------------

    ! --- Input parameters --- !
    CHARACTER(len=200) :: case_type         ! [-] Tag to distinguish turbine
    CHARACTER(len=200) :: rootname          ! [-] Name without extensions of the OpenFAST output files
    CHARACTER(len=200) :: OP_output_dir     ! [-] Directory with OpenFAST output files
    CHARACTER(len=200) :: OP_output_path    ! [-] OpenFAST output file path
    CHARACTER(len=200) :: OP_SubDyn_path    ! [-] OpenFAST SubDyn file path
    LOGICAL            :: Binary            ! [-] Whether OpenFAST files are binary (.outb) or not (.out)
    REAL(WP)           :: WindSpeed         ! [m/s] Wind speed in norm
    REAL(WP)           :: WindDir           ! [deg] Wind direction, 0 deg points to +x axis (anticlockwise from +x)
    REAL(WP)           :: Depth             ! [m] Water depth copied from Environment.nml
    INTEGER(I32)       :: Nmembers          ! [-] Number of structural members in the OpenFAST model
    INTEGER(I32)       :: Nnodes            ! [-] Number of structural nodes in the OpenFAST model
    REAL(WP)           :: AxisPos(2)        ! [m] Position of the turbine axis in xy plane
    REAL(WP)           :: BariPos(2)        ! [m] Position of the turbine barycenter in xy plane
    REAL(WP)           :: From              ! [-] Fraction to start reading input
    REAL(WP)           :: Upto              ! [-] Fraction to stop reading input
    INTEGER(I32)       :: skip              ! [-] Skips time snapshots as time(::skip)
    INTEGER(I32)       :: skipf             ! [-] Skips frequency bins as freq(::skipf)
    LOGICAL            :: filter_freqs      ! [-] Whether to delete non useful frequencies
    REAL(WP)           :: freqs_over        ! [Hz] Frequencies smaller than freqs_over won't be filtered
    LOGICAL            :: is_monopile       ! [-] Whether the turbine is a monopile or a floating platform

END TYPE Config_WindTurbine_t


TYPE :: Config_Environment_t
    ! -------------------------------------------------------------------------
    ! Holds the base environmental/fluid properties read from the input namelist
    ! file (site depth, speed of sound, density, reference pressure). These act
    ! as the default medium properties used by every acoustic solver, unless a
    ! solver-method-specific configuration overrides them.
    ! -------------------------------------------------------------------------

    ! --- Input parameters --- !
    REAL(WP) :: depth_base            ! [m] Site depth
    REAL(WP) :: speed_sound_base      ! [m/s] Base speed of sound in fluid
    REAL(WP) :: density_base          ! [kg/m^3] Base fluid density
    REAL(WP) :: reference_pressure    ! [Pa] Fluid reference pressure

END TYPE Config_Environment_t


TYPE :: Config_Solver_t
    ! -------------------------------------------------------------------------
    ! Holds the generic acoustic-solver settings read from the input namelist
    ! file: which method to run, where to save the results, and the parallel
    ! execution parameters. Method-specific settings live in Config_Method_t
    ! and its extensions instead.
    ! -------------------------------------------------------------------------

    ! --- Input parameters --- !
    LOGICAL            :: debug         ! [-] More info for debugging
    LOGICAL            :: in_cluster    ! [-] Whether we are working on a cluster
    CHARACTER(len=200) :: save_path     ! [-] Where to save results
    CHARACTER(len=50)  :: method        ! [-] Acoustic method to use
    INTEGER(I32)       :: block_size    ! [-] Block size for each CPU in parallel computations

END TYPE Config_Solver_t


TYPE :: GeometricData_t
    ! -------------------------------------------------------------------------
    ! Reference geometry and labels of whatever produced the acoustic sources
    ! (a turbine, a wind farm or hand-made sources). Used by the automatic
    ! defaults of the run_* routines and by the default output file name.
    ! -------------------------------------------------------------------------

    REAL(WP)          :: center(2)   = 0.0_WP       ! [m] Reference centre (auto centre of polar, cylinder, sphere)
    REAL(WP)          :: origin(2)   = 0.0_WP       ! [m] Reference origin (auto start of decay, auto y of sliceXZ)
    REAL(WP)          :: direction   = 0.0_WP       ! [deg] Reference direction, anticlockwise from +x (decay)
    CHARACTER(len=64) :: name_tag    = 'sources'    ! [-] Default output file name tag ('mn', 'fl', 'farm', ...)
    CHARACTER(len=64) :: name_suffix = ''           ! [-] Default output file name suffix ('_ws11.4', ...)

END TYPE GeometricData_t


TYPE :: AcousticSource_t
    ! -------------------------------------------------------------------------
    ! Abstract representation of the acoustic excitation produced by a wind
    ! turbine: a set of point sources with their coordinates, complex nodal
    ! forces per frequency, and the associated frequency vector. This is the
    ! hand-off object between WindTurbine_t % turbine_to_source and the
    ! acoustic solvers.
    ! -------------------------------------------------------------------------

    REAL(WP)   , ALLOCATABLE :: xyz(:,:)        ! [m] Acoustic source coordinates (Nsources, 3)
    COMPLEX(WP), ALLOCATABLE :: Force(:,:,:)    ! [N] Complex source forces (Nfreqs, Nsources, 3)
    REAL(WP)   , ALLOCATABLE :: Freqs(:)        ! [Hz] Frequency vector (Nfreqs)
    TYPE(GeometricData_t)    :: geometry        ! [-] Reference geometry and labels of the emitter

END TYPE AcousticSource_t


TYPE, ABSTRACT :: Config_Method_t
    ! -------------------------------------------------------------------------
    ! Empty abstract base type used purely as a common label for every
    ! solver-method-specific configuration type (e.g. Config_ImageMethod_t).
    ! It lets Config_t % method be declared polymorphic and allocated at run
    ! time to whichever concrete configuration matches the selected method.
    ! -------------------------------------------------------------------------

END TYPE Config_Method_t


TYPE, EXTENDS(Config_Method_t) :: Config_ImageMethod_t
    ! -------------------------------------------------------------------------
    ! Method-of-Images solver configuration, read from the ImagesMethod.nml
    ! namelist. Sets how many image sources are used to enforce the boundary
    ! conditions at the free surface and the seabed, the boundary types and
    ! heights, and the per-boundary reflection attenuation coefficients.
    ! -------------------------------------------------------------------------

    INTEGER(I32) :: N_images               ! [-] Number of image source levels used to enforce boundary conditions
    INTEGER(I32) :: UpperBC                ! [-] Upper boundary condition flag (-1: pressure-release surface, p = 0)
    INTEGER(I32) :: LowerBC                ! [-] Lower boundary condition flag (+1: rigid boundary, dp/dn = 0)
    REAL(WP)     :: UpperHBC               ! [m] z-coordinate of the upper plane; seabed = UpperHBC - depth_base
    REAL(WP)     :: attenuation_surface    ! [-] Attenuation coefficient applied to surface reflections
    REAL(WP)     :: attenuation_bottom     ! [-] Attenuation coefficient applied to bottom reflections
    REAL(WP)     :: epsilon                ! [-] Small numerical tolerance used to avoid singularities

END TYPE Config_ImageMethod_t


TYPE :: ImageSystem_t
    ! -------------------------------------------------------------------------
    ! Helper type used internally by the Method of Images to hold the full set
    ! of real and image sources built from the original AcousticSource_t: their
    ! coordinates, complex forces, and which boundary condition each image was
    ! generated from.
    ! -------------------------------------------------------------------------

    REAL(WP)    , ALLOCATABLE :: xyz(:,:)        ! [m] Nodes coordinates (Nnodes, 3)
    COMPLEX(WP) , ALLOCATABLE :: Force(:,:,:)    ! [N] Force array (Nfreqs, Nnodes, 3)
    INTEGER(I32), ALLOCATABLE :: BC_all(:)       ! [-] Boundary conditions array

END TYPE ImageSystem_t


TYPE, EXTENDS(Config_Method_t) :: Config_AnalyticalNormalModes_t
    ! -------------------------------------------------------------------------
    ! Analytical Normal Modes solver configuration, read from the
    ! AnalyticalNormalModes.nml namelist. Sets how many normal modes are
    ! retained in the modal expansion, the heights of the upper and lower
    ! waveguide boundaries, and the singularity tolerance.
    ! -------------------------------------------------------------------------

    INTEGER(I32) :: Nmodes      ! [-] Number of retained modes
    REAL(WP)     :: UpperHBC    ! [m] z-coordinate of the upper plane; seabed = UpperHBC - depth_base
    REAL(WP)     :: epsilon     ! [-] Small numerical tolerance used to avoid singularities

END TYPE Config_AnalyticalNormalModes_t


TYPE :: Config_DoRun_t
    ! -------------------------------------------------------------------------
    ! Which post-processing steps run_all() executes. do_all is a master
    ! switch: when .true. every step runs regardless of its individual flag;
    ! when .false., only the steps whose own do_<name> is .true. run.
    ! -------------------------------------------------------------------------

    LOGICAL :: do_all          ! [-] Master switch: run every step regardless of the flags below
    LOGICAL :: do_spectrums    ! [-] Spectrums at observer points (run_spectrums)
    LOGICAL :: do_polar        ! [-] Polar contour (run_polar)
    LOGICAL :: do_cylinder     ! [-] Cylinder surface (run_cylinder)
    LOGICAL :: do_decay        ! [-] Distance decay along wind direction (run_decay)
    LOGICAL :: do_line         ! [-] Line between two points (run_line)
    LOGICAL :: do_sliceXY      ! [-] Horizontal slice XY (run_sliceXY)
    LOGICAL :: do_sliceXZ      ! [-] Vertical slice XZ (run_sliceXZ)
    LOGICAL :: do_sphere       ! [-] Spherical radiation, free-field (run_sphere)

END TYPE Config_DoRun_t


TYPE :: Config_Spectrums_t
    ! -------------------------------------------------------------------------
    ! Configuration for run_spectrums(): pressure at a handful of fixed
    ! observer points read from the namelist. If n_observers is left unset,
    ! 3 hardcoded default points are used instead.
    ! -------------------------------------------------------------------------

    INTEGER(I32)              :: n_observers       ! [-] Number of observers
    REAL(WP)    , ALLOCATABLE :: observers_x(:)    ! [m] x-coordinates for the observers
    REAL(WP)    , ALLOCATABLE :: observers_y(:)    ! [m] y-coordinates for the observers
    REAL(WP)    , ALLOCATABLE :: observers_z(:)    ! [m] z-coordinates for the observers

END TYPE Config_Spectrums_t


TYPE :: Config_Polar_t
    ! -------------------------------------------------------------------------
    ! Configuration for run_polar(): pressure on a horizontal circle of
    ! observers around the turbine.
    ! -------------------------------------------------------------------------

    REAL(WP)     :: r_polar            ! [m] Circle radius
    REAL(WP)     :: z_polar            ! [m] Circle depth
    INTEGER(I32) :: n_polar            ! [-] Number of azimuthal observer points
    REAL(WP)     :: center_polar(2)    ! [m] Circle centre (x, y)

END TYPE Config_Polar_t


TYPE :: Config_Cylinder_t
    ! -------------------------------------------------------------------------
    ! Configuration for run_cylinder(): pressure on a vertical cylindrical
    ! surface of observers around the turbine.
    ! -------------------------------------------------------------------------

    REAL(WP)     :: r_cyl            ! [m] Cylinder radius
    INTEGER(I32) :: n_cyl            ! [-] Number of azimuthal observer points
    INTEGER(I32) :: nz_cyl           ! [-] Number of vertical observer levels
    REAL(WP)     :: center_cyl(2)    ! [m] Cylinder axis centre (x, y)

END TYPE Config_Cylinder_t


TYPE :: Config_Decay_t
    ! -------------------------------------------------------------------------
    ! Configuration for run_decay(): pressure decay along a line in the wind
    ! direction, starting from the turbine.
    ! -------------------------------------------------------------------------

    REAL(WP)     :: distance_decay(2)    ! [m] Start and end distance from the turbine
    INTEGER(I32) :: n_decay              ! [-] Number of observer points along the line
    REAL(WP)     :: z_decay              ! [m] Line depth
    LOGICAL      :: logspace_decay       ! [-] Whether to space points logarithmically instead of linearly

END TYPE Config_Decay_t


TYPE :: Config_Line_t
    ! -------------------------------------------------------------------------
    ! Configuration for run_line(): pressure along an arbitrary straight line
    ! between two points.
    ! -------------------------------------------------------------------------

    REAL(WP)     :: p1_line(3)       ! [m] Line start point (x, y, z)
    REAL(WP)     :: p2_line(3)       ! [m] Line end point (x, y, z)
    INTEGER(I32) :: n_line           ! [-] Number of observer points along the line
    LOGICAL      :: logspace_line    ! [-] Whether to space points logarithmically instead of linearly

END TYPE Config_Line_t


TYPE :: Config_SliceXY_t
    ! -------------------------------------------------------------------------
    ! Configuration for run_sliceXY(): pressure on a horizontal (x, y) grid of
    ! observers at a fixed depth.
    ! -------------------------------------------------------------------------

    REAL(WP)     :: z_xy            ! [m] Slice depth
    INTEGER(I32) :: nx_xy           ! [-] Number of observer points along x
    INTEGER(I32) :: ny_xy           ! [-] Number of observer points along y
    REAL(WP)     :: xlim_xy(2)      ! [m] Slice extent along x (min, max)
    REAL(WP)     :: ylim_xy(2)      ! [m] Slice extent along y (min, max)
    REAL(WP)     :: center_xy(2)    ! [m] Slice centre (x, y)

END TYPE Config_SliceXY_t


TYPE :: Config_SliceXZ_t
    ! -------------------------------------------------------------------------
    ! Configuration for run_sliceXZ(): pressure on a vertical (x, z) grid of
    ! observers at a fixed y-position.
    ! -------------------------------------------------------------------------

    REAL(WP)     :: y_xz          ! [m] Slice y-position
    INTEGER(I32) :: nx_xz         ! [-] Number of observer points along x
    INTEGER(I32) :: nz_xz         ! [-] Number of observer points along z
    REAL(WP)     :: xlim_xz(2)    ! [m] Slice extent along x (min, max)
    REAL(WP)     :: zlim_xz(2)    ! [m] Slice extent along z (min, max)

END TYPE Config_SliceXZ_t


TYPE :: Config_Sphere_t
    ! -------------------------------------------------------------------------
    ! Configuration for run_sphere(): free-field pressure on a spherical
    ! surface of observers around the turbine (auto-enlarged if too small to
    ! enclose every turbine node).
    ! -------------------------------------------------------------------------

    REAL(WP)     :: r_sphere            ! [m] Sphere radius
    INTEGER(I32) :: n_theta_sphere      ! [-] Number of azimuthal observer points
    INTEGER(I32) :: nz_sphere           ! [-] Number of vertical observer levels
    REAL(WP)     :: center_sphere(3)    ! [m] Sphere centre (x, y, z)

END TYPE Config_Sphere_t


TYPE :: Config_Computation_t
    ! -------------------------------------------------------------------------
    ! Aggregates the per-step configuration for every run_* post-processing
    ! routine that run_all() can execute, read from ComputationSetupTurbine.nml.
    ! -------------------------------------------------------------------------

    TYPE(Config_DoRun_t)     :: do_run              ! [-] Which steps to run
    TYPE(Config_Spectrums_t) :: params_spectrums    ! [-] run_spectrums() configuration
    TYPE(Config_Polar_t)     :: params_polar        ! [-] run_polar() configuration
    TYPE(Config_Cylinder_t)  :: params_cylinder     ! [-] run_cylinder() configuration
    TYPE(Config_Decay_t)     :: params_decay        ! [-] run_decay() configuration
    TYPE(Config_Line_t)      :: params_line         ! [-] run_line() configuration
    TYPE(Config_SliceXY_t)   :: params_sliceXY      ! [-] run_sliceXY() configuration
    TYPE(Config_SliceXZ_t)   :: params_sliceXZ      ! [-] run_sliceXZ() configuration
    TYPE(Config_Sphere_t)    :: params_sphere       ! [-] run_sphere() configuration

END TYPE Config_Computation_t


TYPE :: Config_t
    ! -------------------------------------------------------------------------
    ! Top-level configuration object built by Args_parser from the input
    ! namelist files. It aggregates the turbine, environment and solver
    ! configurations together with the polymorphic, method-specific
    ! configuration (e.g. Config_ImageMethod_t) selected at run time.
    ! -------------------------------------------------------------------------

    TYPE(Config_WindTurbine_t)              :: turbine       ! [-] Wind turbine configuration
    TYPE(Config_Environment_t)              :: env           ! [-] Environment/fluid configuration
    TYPE(Config_Solver_t)                   :: solver        ! [-] Acoustic solver configuration
    TYPE(Config_Computation_t)              :: run_inputs    ! [-] Per-run_* post-processing step configuration
    CLASS(Config_Method_t)    , ALLOCATABLE :: method        ! [-] Polymorphic, solver-method-specific configuration

END TYPE Config_t

END MODULE Types
