! ===============================================================================
! MODULE: ConfigTypes
! ===============================================================================
! Defines the configuration types filled in from the input namelist files:
! turbine, wind farm, environment, solver, the per-method solver settings
! (Config_Method_t and its extensions), the per-step post-processing
! settings, and the top-level Config_t that aggregates all of them.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE ConfigTypes

USE Kinds, ONLY: I32, WP

IMPLICIT NONE

PRIVATE
PUBLIC :: Config_WindTurbine_t, Config_WindFarm_t, Config_Environment_t, Config_Solver_t, Config_Method_t,          &
          Config_ImageMethod_t, Config_AnalyticalNormalModes_t, Config_DoRun_t, Config_Spectrums_t, Config_Polar_t, &
          Config_Cylinder_t, Config_Decay_t, Config_Line_t, Config_SliceXY_t, Config_SliceXZ_t,                     &
          Config_SliceVertical_t, Config_Sphere_t, Config_Computation_t, Config_t


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

END TYPE Config_WindTurbine_t


TYPE :: Config_WindFarm_t
    ! -------------------------------------------------------------------------
    ! Holds the wind farm configuration read from WindFarm.nml: one
    ! Config_WindTurbine_t per &turbine block, each one already completed
    ! with the settings shared by the whole farm (&farm_config) and with the
    ! water depth of Environment.nml.
    ! -------------------------------------------------------------------------

    INTEGER(I32)                            :: Nturbines      ! [-] Number of turbines in the farm
    CHARACTER(len=100)                      :: farm_name      ! [-] File name label ('None': automatic)
    TYPE(Config_WindTurbine_t), ALLOCATABLE :: turbines(:)    ! [-] Per-turbine configuration (Nturbines)

END TYPE Config_WindFarm_t


TYPE :: Config_FrequencySource_t
    ! -------------------------------------------------------------------------
    ! Holds the frequency-domain point sources read from FrequencySource.nml:
    ! their positions and the complex force of every source at every frequency.
    ! -------------------------------------------------------------------------

    ! --- Input parameters --- !
    INTEGER(I32)             :: Nsources         ! [-] Number of sources
    INTEGER(I32)             :: Nfreqs           ! [-] Number of frequencies
    REAL(WP)   , ALLOCATABLE :: Freqs(:)         ! [Hz] Frequency vector (Nfreqs)
    REAL(WP)   , ALLOCATABLE :: source_xyz(:,:)  ! [m] Source coordinates (Nsources, 3)
    COMPLEX(WP), ALLOCATABLE :: Force(:,:,:)     ! [N] Complex force on the fluid (Nfreqs, Nsources, 3)

END TYPE Config_FrequencySource_t


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
    LOGICAL            :: debug             ! [-] More info for debugging
    LOGICAL            :: in_cluster        ! [-] Whether we are working on a cluster
    LOGICAL            :: update_existing   ! [-] Update an existing output file of the same case instead of rebuilding it
    CHARACTER(len=200) :: save_path         ! [-] Where to save results
    CHARACTER(len=100) :: extra_label       ! [-] Label appended to the output file name ('None' or '': none)
    CHARACTER(len=50)  :: method            ! [-] Acoustic method to use
    CHARACTER(len=20)  :: simulation_type   ! [-] WindTurbine, WindFarm or FrequencySource
    INTEGER(I32)       :: block_size        ! [-] Block size for each CPU in parallel computations
    REAL(WP)           :: mask_distance     ! [m] Observers closer than this to a source node get NaN (0: none)

END TYPE Config_Solver_t


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
    ! conditions at the free surface and the seabed, the boundary types, 
    ! and the per-boundary reflection attenuation coefficients.
    ! -------------------------------------------------------------------------

    INTEGER(I32) :: N_images               ! [-] Number of image source levels used to enforce boundary conditions
    INTEGER(I32) :: UpperBC                ! [-] Upper boundary condition flag (-1: pressure-release surface, p = 0)
    INTEGER(I32) :: LowerBC                ! [-] Lower boundary condition flag (+1: rigid boundary, dp/dn = 0)
    REAL(WP)     :: attenuation_surface    ! [-] Attenuation coefficient applied to surface reflections
    REAL(WP)     :: attenuation_bottom     ! [-] Attenuation coefficient applied to bottom reflections
    REAL(WP)     :: epsilon                ! [m] Minimum source-observer distance used to floor r

END TYPE Config_ImageMethod_t


TYPE, EXTENDS(Config_Method_t) :: Config_AnalyticalNormalModes_t
    ! -------------------------------------------------------------------------
    ! Analytical Normal Modes solver configuration, read from the
    ! AnalyticalNormalModes.nml namelist. Sets how many normal modes are
    ! retained in the modal expansion and the singularity tolerance.
    ! -------------------------------------------------------------------------

    INTEGER(I32) :: Nmodes      ! [-] Number of retained modes
    REAL(WP)     :: epsilon     ! [m] Minimum horizontal source-observer distance used to floor r

END TYPE Config_AnalyticalNormalModes_t


TYPE :: Config_DoRun_t
    ! -------------------------------------------------------------------------
    ! Which post-processing steps run() executes. do_all is a master switch
    ! resolved by the parser: when .true. every step that makes sense for the
    ! emitter (turbine or farm) is switched on; when .false., only the steps
    ! whose own do_<name> is .true. run. Steps without meaning for a wind farm
    ! (sphere, decay) are always switched off by the parser in farm mode.
    ! -------------------------------------------------------------------------

    LOGICAL :: do_all              ! [-] Master switch: run every step valid for the emitter
    LOGICAL :: do_spectrums        ! [-] Spectrums at observer points (run_spectrums)
    LOGICAL :: do_polar            ! [-] Polar contour (run_polar)
    LOGICAL :: do_cylinder         ! [-] Cylinder surface (run_cylinder)
    LOGICAL :: do_decay            ! [-] Distance decay along wind direction (run_decay), turbine only
    LOGICAL :: do_line             ! [-] Line between two points (run_line)
    LOGICAL :: do_sliceXY          ! [-] Horizontal slice XY (run_sliceXY)
    LOGICAL :: do_sliceXZ          ! [-] Vertical slice XZ (run_sliceXZ)
    LOGICAL :: do_sliceVertical    ! [-] Vertical slice at any azimuth (run_sliceVertical)
    LOGICAL :: do_sphere           ! [-] Spherical radiation, free-field (run_sphere), turbine only

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


TYPE :: Config_SliceVertical_t
    ! -------------------------------------------------------------------------
    ! Configuration for run_sliceVertical(): pressure on a vertical plane of
    ! observers whose horizontal normal points along an arbitrary azimuth.
    ! The plane spans the horizontal axis (-sin, cos)(azimuth) and the depth.
    ! Valid for a single turbine and for a wind farm alike.
    ! -------------------------------------------------------------------------

    REAL(WP)     :: azimuth_sv      ! [deg] Azimuth of the plane normal (0 deg: normal along +x, YZ plane; anticlockwise)
    REAL(WP)     :: center_sv(2)    ! [m] Plane centre (x, y)
    REAL(WP)     :: width_sv        ! [m] Full horizontal span of the plane
    REAL(WP)     :: zlim_sv(2)      ! [m] Depth extent (min, max)
    INTEGER(I32) :: nu_sv           ! [-] Number of observer points along the span
    INTEGER(I32) :: nz_sv           ! [-] Number of observer points along z
    REAL(WP)     :: margin_sv       ! [m] Margin added on both sides to the automatic width
    REAL(WP)     :: ulim_sv(2)      ! [m] Span limits actually used (u_min, u_max), resolved at run time

END TYPE Config_SliceVertical_t


TYPE :: Config_Sphere_t
    ! -------------------------------------------------------------------------
    ! Configuration for run_sphere(): free-field pressure on a spherical
    ! surface of observers around the turbine (auto-enlarged if too small to
    ! enclose every turbine node).
    ! -------------------------------------------------------------------------

    REAL(WP)     :: r_sphere            ! [m] Sphere radius
    INTEGER(I32) :: n_theta_sphere      ! [-] Number of azimuthal observer points
    INTEGER(I32) :: nz_sphere           ! [-] Number of vertical observer levels
    REAL(WP)     :: f_farfield_sphere   ! [-] Frequency from which all f> f_farfield are considered to be in far-field
    REAL(WP)     :: center_sphere(3)    ! [m] Sphere centre (x, y, z)

END TYPE Config_Sphere_t


TYPE :: Config_Computation_t
    ! -------------------------------------------------------------------------
    ! Aggregates the per-step configuration for every run_* post-processing
    ! routine that run() can execute, read from ComputationSetupTurbine.nml
    ! (single turbine) or ComputationSetupFarm.nml (wind farm).
    ! -------------------------------------------------------------------------

    TYPE(Config_DoRun_t)         :: do_run                  ! [-] Which steps to run
    TYPE(Config_Spectrums_t)     :: params_spectrums        ! [-] run_spectrums() configuration
    TYPE(Config_Polar_t)         :: params_polar            ! [-] run_polar() configuration
    TYPE(Config_Cylinder_t)      :: params_cylinder         ! [-] run_cylinder() configuration
    TYPE(Config_Decay_t)         :: params_decay            ! [-] run_decay() configuration
    TYPE(Config_Line_t)          :: params_line             ! [-] run_line() configuration
    TYPE(Config_SliceXY_t)       :: params_sliceXY          ! [-] run_sliceXY() configuration
    TYPE(Config_SliceXZ_t)       :: params_sliceXZ          ! [-] run_sliceXZ() configuration
    TYPE(Config_SliceVertical_t) :: params_sliceVertical    ! [-] run_sliceVertical() configuration
    TYPE(Config_Sphere_t)        :: params_sphere           ! [-] run_sphere() configuration

END TYPE Config_Computation_t


TYPE :: Config_t
    ! -------------------------------------------------------------------------
    ! Top-level configuration object built by Args_parser from the input
    ! namelist files. It aggregates the turbine, environment and solver
    ! configurations together with the polymorphic, method-specific
    ! configuration (e.g. Config_ImageMethod_t) selected at run time.
    ! -------------------------------------------------------------------------

    TYPE(Config_WindTurbine_t)              :: turbine          ! [-] Wind turbine configuration
    TYPE(Config_WindFarm_t)                 :: farm             ! [-] Wind farm configuration
    TYPE(Config_Environment_t)              :: env              ! [-] Environment/fluid configuration
    TYPE(Config_Solver_t)                   :: solver           ! [-] Acoustic solver configuration
    TYPE(Config_Computation_t)              :: run_inputs       ! [-] Per-run_* post-processing step configuration
    TYPE(Config_FrequencySource_t)          :: freq_source      ! [-] Frequency-domain source configuration
    CLASS(Config_Method_t)    , ALLOCATABLE :: method           ! [-] Polymorphic, solver-method-specific configuration

END TYPE Config_t

END MODULE ConfigTypes
