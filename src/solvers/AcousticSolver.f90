! ===============================================================================
! MODULE: AcousticSolver
! ===============================================================================
! Defines AcousticSolver_t, the abstract parent type for all acoustic
! propagation solvers (ImageMethod_t and AnalyticalNormalModes_t). It holds
! the shared solver state (sources and their reference geometry,
! environment, save path, block size, debug/cluster flags), the observer
! masking and the water-column queries (z_surface, z_seabed, in_water), and
! declares the generic observer-geometry post-processing
! routines (spectrums, polar, cylinder, decay, line, sliceXY, sliceXZ,
! sliceVertical, sphere) shared by every concrete solver, implemented in the
! submodule AcousticSolver_steps, and the output-file management,
! implemented in the submodule AcousticSolver_output. The solver never
! knows what produced the sources (a turbine, a wind farm or hand-made
! sources): which steps make sense for each emitter is decided by the
! parser through the do_run flags.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE AcousticSolver

USE Kinds
USE Emitter, ONLY: Emitter_t
USE SourceTypes, ONLY: AcousticSource_t
USE ConfigTypes, ONLY: Config_Environment_t, Config_Computation_t, Config_t

IMPLICIT NONE

PRIVATE
PUBLIC :: AcousticSolver_t, pressure_norm


TYPE, ABSTRACT :: AcousticSolver_t
    ! -------------------------------------------------------------------------
    ! Parent class for all acoustic propagation solvers (ImageMethod_t and
    ! AnalyticalNormalModes_t). Holds the shared sources/environment/turbine
    ! state and implements the generic observer-geometry post-processing
    ! routines; init, compute_pressure and the parameter/path/summary
    ! routines are deferred.
    ! -------------------------------------------------------------------------

    TYPE(AcousticSource_t)     :: sources               ! [-] Acoustic point sources and their reference geometry
    TYPE(Config_Environment_t) :: env                   ! [-] Environment configuration (depth, sound speed, ...)
    TYPE(Config_Computation_t) :: run_inputs            ! [-] Per-run_* post-processing step configuration
    INTEGER(I32)               :: block_size            ! [-] Number of observers processed simultaneously
    LOGICAL                    :: debug                 ! [-] Whether to print extra debug info (e.g. pressure norms)
    LOGICAL                    :: in_cluster            ! [-] Whether the solver runs on a cluster
    LOGICAL                    :: update_existing       ! [-] Whether an existing output file of the same case is updated
    REAL(WP)                   :: mask_distance         ! [m] Observers closer than this to a source node get NaN
    CHARACTER(len=50)          :: solver_name = 'None'  ! [-] Solver name
    CHARACTER(len=512)         :: save_path   = 'None'  ! [-] Path to save acoustics
    CHARACTER(len=100)         :: extra_label = 'None'  ! [-] Label appended to the file name ('None' or '': none)

    CONTAINS

    ! --- Public type-bound procedures --- !
    PROCEDURE :: check_observers_distances              ! Observer mask: .false. closer than mask_distance to a source
    PROCEDURE :: compute_pressure_masked                ! compute_pressure with NaN at the masked observers
    PROCEDURE :: resolve_save_path                      ! Output path shared by every solver (build_save_path helper)
    PROCEDURE :: check_update                           ! Keeps an existing output only if it holds the same case
    PROCEDURE :: save_shared_parameters                 ! Datasets written by every solver (Method, ENBW)
    PROCEDURE :: update_parameters                      ! Appends run_* parameters actually used to the .dat
    PROCEDURE :: run                                    ! Config-driven entry point: every enabled do_* step
    PROCEDURE :: run_spectrums                          ! Subroutines to compute acoustic pressure
    PROCEDURE :: run_polar
    PROCEDURE :: run_cylinder
    PROCEDURE :: run_decay
    PROCEDURE :: run_line
    PROCEDURE :: run_sliceXY
    PROCEDURE :: run_sliceXZ
    PROCEDURE :: run_sliceVertical
    PROCEDURE :: run_sphere => dummy_run_sphere

    ! --- Water-column queries (a method overrides them to publish its boundaries) --- !
    PROCEDURE :: z_surface => default_z_surface         ! z of the water surface at (x, y)
    PROCEDURE :: z_seabed  => default_z_seabed          ! z of the seabed at (x, y)
    PROCEDURE :: in_water  => default_in_water          ! Observer mask: .true. inside the water column

    ! --- Deferred (abstract) procedures --- !
    PROCEDURE(init_i)             , DEFERRED :: init
    PROCEDURE(get_name_i)         , DEFERRED :: get_name
    PROCEDURE(compute_pressure_i) , DEFERRED :: compute_pressure
    PROCEDURE(save_parameters_i)  , DEFERRED :: save_parameters
    PROCEDURE(build_save_path_i)  , DEFERRED :: build_save_path
    PROCEDURE(write_txt_summary_i), DEFERRED :: write_txt_summary

END TYPE AcousticSolver_t


! ===============================================================================
! ABSTRACT INTERFACES FOR DEFERRED PROCEDURES
! ===============================================================================
ABSTRACT INTERFACE

    SUBROUTINE init_i(self, sources, config)
        ! -------------------------------------------------------------------------
        ! Constructor: copies the shared and method-specific configuration into
        ! the solver and prepares it for the pressure computations.
        ! -------------------------------------------------------------------------

        IMPORT :: AcousticSolver_t, AcousticSource_t, Config_t

        CLASS(AcousticSolver_t), INTENT(INOUT)              :: self       ! [-] Solver instance
        TYPE(AcousticSource_t) , INTENT(INOUT), ALLOCATABLE :: sources    ! [-] Acoustic point-dipole sources
        TYPE(Config_t)         , INTENT(INOUT), ALLOCATABLE :: config     ! [-] Global configuration

    END SUBROUTINE init_i

    FUNCTION get_name_i(self) RESULT(name)
        ! -------------------------------------------------------------------------
        ! Returns the human-readable solver name.
        ! -------------------------------------------------------------------------

        IMPORT :: AcousticSolver_t

        CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance
        CHARACTER(len=:)       , ALLOCATABLE   :: name    ! [-] Human-readable solver name

    END FUNCTION get_name_i

    SUBROUTINE compute_pressure_i(self, observers, block_size, total_pressure)
        ! -------------------------------------------------------------------------
        ! Evaluates the total radiated pressure field at a set of observer points
        ! for every source frequency.
        ! -------------------------------------------------------------------------

        IMPORT :: AcousticSolver_t, AcousticSource_t, WP, I32

        CLASS(AcousticSolver_t), INTENT(INOUT)            :: self                   ! [-] Solver instance
        REAL(WP)               , INTENT(IN)               :: observers(:,:)         ! [m] Observers coords (Nobs, 3)
        INTEGER(I32)           , INTENT(IN), OPTIONAL     :: block_size             ! [-] Observers per block
        COMPLEX(WP)            , INTENT(OUT), ALLOCATABLE :: total_pressure(:,:)    ! [Pa] Pressure (Nfreqs, Nobs)

    END SUBROUTINE compute_pressure_i

    SUBROUTINE save_parameters_i(self, emitter_obj)
        ! -------------------------------------------------------------------------
        ! Saves the solver, environment and emitter parameters to the HDF5
        ! output file and writes the .dat text summary.
        ! -------------------------------------------------------------------------

        IMPORT :: AcousticSolver_t, Emitter_t

        CLASS(AcousticSolver_t), INTENT(IN) :: self       ! [-] Solver instance
        CLASS(Emitter_t)       , INTENT(IN) :: emitter_obj    ! [-] Emitter that produced the sources

    END SUBROUTINE save_parameters_i

    SUBROUTINE build_save_path_i(self)
        ! -------------------------------------------------------------------------
        ! Builds the name of the HDF5 output file.
        ! -------------------------------------------------------------------------

        IMPORT :: AcousticSolver_t

        CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    END SUBROUTINE build_save_path_i

    SUBROUTINE write_txt_summary_i(self, emitter_obj)
        ! -------------------------------------------------------------------------
        ! Writes the human-readable .dat summary next to the HDF5 output file.
        ! -------------------------------------------------------------------------

        IMPORT :: AcousticSolver_t, Emitter_t

        CLASS(AcousticSolver_t), INTENT(IN) :: self       ! [-] Solver instance
        CLASS(Emitter_t)       , INTENT(IN) :: emitter_obj    ! [-] Emitter that produced the sources

    END SUBROUTINE write_txt_summary_i

END INTERFACE


! ===============================================================================
! INTERFACES OF THE PROCEDURES IMPLEMENTED IN SUBMODULES
! AcousticSolver_output: output path, update of an existing output, shared datasets
! AcousticSolver_steps  : run() and the run_* post-processing steps
! ===============================================================================
INTERFACE

    MODULE SUBROUTINE resolve_save_path(self, method_tag)
        ! -------------------------------------------------------------------------
        ! Resolves the final HDF5 output path, shared by every solver's
        ! build_save_path. If the configured save_path is "None"/"Default", the
        ! name is built from the emitter labels carried by the sources around the
        ! solver's method tag: out_dir // name_tag // method_tag // name_suffix,
        ! e.g. ./turbine_acoustic_data/plot_mn_ANM24_ws11.4.hdf5 or
        ! ./farm_acoustic_data/Farm_2_DTU10MN_ANM24_ws8.0-11.4.hdf5. Otherwise the
        ! .hdf5 extension is appended to the given path. In both cases, a set
        ! extra_label (neither 'None' nor empty) is appended before the extension
        ! as "_" // extra_label. The destination folder is created if needed. A
        ! previous file with the same name is deleted when update_existing is
        ! .false. (save_variable only overwrites the datasets it writes, so those
        ! of an older run would otherwise survive). With update_existing it is
        ! kept and its .dat moved to .dat.prev: check_update then decides whether
        ! it holds the same case and update_parameters merges the step blocks.
        ! -------------------------------------------------------------------------

        CLASS(AcousticSolver_t), INTENT(INOUT) :: self          ! [-] Solver instance
        CHARACTER(len=*)       , INTENT(IN)    :: method_tag    ! [-] Method and resolution tag, e.g. "_ANM24"

    END SUBROUTINE resolve_save_path

    MODULE SUBROUTINE check_update(self, emitter_obj)
        ! -------------------------------------------------------------------------
        ! Called by save_parameters before anything is written. If resolve_save_path
        ! kept an existing output file (update_existing), it is updated only if it
        ! holds the same case: the same frequency vector (OpenFAST run, From/Upto,
        ! skipf, frequency filtering) and the same .dat header (emitter,
        ! environment and solver settings; block size and cluster flag excluded).
        ! The new header is written here to compare it. Otherwise the file and its
        ! .dat.prev are deleted with a warning and the output is rebuilt.
        ! -------------------------------------------------------------------------

        CLASS(AcousticSolver_t), INTENT(IN) :: self       ! [-] Solver instance
        CLASS(Emitter_t)       , INTENT(IN) :: emitter_obj    ! [-] Emitter that produced the sources

    END SUBROUTINE check_update

    MODULE SUBROUTINE update_parameters(self)
        ! -------------------------------------------------------------------------
        ! Appends to the .dat summary already written by save_parameters() one
        ! block per post-processing step that actually ran (per do_run), with the
        ! resolved parameters run() used for it, then closes the file with its
        ! final separator. In an update of an existing output (.dat.prev kept by
        ! resolve_save_path and check_update), the steps that did not run now
        ! keep their block from the previous .dat, tagged "(previous run)", and
        ! .dat.prev is deleted. Must be called after run(), since every step only
        ! resolves its own DEFAULT_R/DEFAULT_I sentinels when it actually
        ! executes.
        ! -------------------------------------------------------------------------

        CLASS(AcousticSolver_t), INTENT(IN) :: self    ! [-] Solver instance

    END SUBROUTINE update_parameters

    MODULE SUBROUTINE save_shared_parameters(self, emitter_obj)
        ! -------------------------------------------------------------------------
        ! Writes the datasets shared by every solver, called by each solver's
        ! save_parameters right after its own datasets: the solver name (Method,
        ! from get_name) and, only when the emitter forces come from a windowed
        ! FFT (get_enbw > 0), the equivalent noise bandwidth of that window in
        ! bins (ENBW). Neither is written in the .dat summary.
        ! -------------------------------------------------------------------------

        CLASS(AcousticSolver_t), INTENT(IN) :: self       ! [-] Solver instance
        CLASS(Emitter_t)       , INTENT(IN) :: emitter_obj    ! [-] Emitter that produced the sources

    END SUBROUTINE save_shared_parameters

    MODULE SUBROUTINE run_spectrums(self)
        ! -------------------------------------------------------------------------
        ! Computes the pressure spectrum at a small set of fixed observer points
        ! (default: 3 points along +x at mid-depth). Saves the pressure, observer
        ! coordinates and frequency vector to the solver's HDF5 file.
        ! -------------------------------------------------------------------------

        CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    END SUBROUTINE run_spectrums

    MODULE SUBROUTINE run_polar(self)
        ! -------------------------------------------------------------------------
        ! Computes pressure on a horizontal circle of radius r_polar around the
        ! turbine barycenter, at fixed depth z_polar. Saves the pressure plus the
        ! polar geometry (radius, depth, azimuth angles, observer coordinates) to
        ! the solver's HDF5 file.
        ! -------------------------------------------------------------------------

        CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    END SUBROUTINE run_polar

    MODULE SUBROUTINE run_cylinder(self)
        ! -------------------------------------------------------------------------
        ! Computes pressure on a vertical cylindrical surface (radius r_cyl,
        ! n_cyl azimuth points x nz_cyl depth levels) around the turbine
        ! barycenter. Saves the reshaped pressure/observer grids and the cylinder
        ! geometry metadata to the solver's HDF5 file. The depth levels include
        ! both ends (z = -Depth and z = 0), so dA_cylinder is the area of an
        ! interior cell: a surface integral must weight the first and last
        ! levels by 1/2 (trapezoid rule).
        ! -------------------------------------------------------------------------

        CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    END SUBROUTINE run_cylinder

    MODULE SUBROUTINE run_decay(self)
        ! -------------------------------------------------------------------------
        ! Computes pressure along a straight line following the wind direction
        ! from the turbine axis, at fixed depth z_decay, with linear or log-spaced
        ! distances. Saves the pressure and the observer/distance geometry to the
        ! solver's HDF5 file. n_decay = 1 places a single observer at the start
        ! distance; n_decay < 1, and an end distance <= 0 with log spacing, stop
        ! the run. A start distance <= 0 is moved to 1 m and persisted, so the
        ! .dat summary reports the distances actually used.
        ! -------------------------------------------------------------------------

        CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    END SUBROUTINE run_decay

    MODULE SUBROUTINE run_line(self)
        ! -------------------------------------------------------------------------
        ! Computes pressure along an arbitrary straight line between two points
        ! p1_line and p2_line, with linear or log-spaced distances. Saves the
        ! pressure and the observer/distance geometry to the solver's HDF5 file.
        ! The automatic n_line always keeps at least the two end points;
        ! n_line = 1 places a single observer at the first spacing point, and
        ! n_line < 1 stops the run.
        ! -------------------------------------------------------------------------

        CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    END SUBROUTINE run_line

    MODULE SUBROUTINE run_sliceXY(self)
        ! -------------------------------------------------------------------------
        ! Computes pressure on a horizontal (x,y) grid at fixed depth z_xy, over
        ! the given x/y limits. Saves the reshaped pressure/observer grids and the
        ! slice geometry metadata to the solver's HDF5 file.
        ! -------------------------------------------------------------------------

        CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    END SUBROUTINE run_sliceXY

    MODULE SUBROUTINE run_sliceXZ(self)
        ! -------------------------------------------------------------------------
        ! Computes pressure on a vertical (x,z) grid at fixed y_xz, over the given
        ! x/z limits. Saves the reshaped pressure/observer grids and the slice
        ! geometry metadata to the solver's HDF5 file.
        ! -------------------------------------------------------------------------

        CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    END SUBROUTINE run_sliceXZ

    MODULE SUBROUTINE run_sliceVertical(self)
        ! -------------------------------------------------------------------------
        ! Computes pressure on a vertical plane whose horizontal normal points
        ! along azimuth_sv (0 deg: normal along +x, i.e. the YZ plane; measured
        ! anticlockwise, as the wind direction and the polar angle). The plane
        ! spans the horizontal axis span = (-sin, cos)(azimuth) and the depth.
        ! With the automatic width every source is projected onto the span axis
        ! and the plane covers all of them plus margin_sv on both sides, so the
        ! same routine is valid for a single turbine and for a wind farm. If
        ! run_sliceXZ already covers the same plane (normal along +-y through
        ! y = y_xz) the step is skipped and switched off in do_run, so the .dat
        ! summary does not report it. Observers are ordered with z varying
        ! fastest, as in the Python implementation. Saves the reshaped
        ! pressure/observer grids and the plane geometry metadata to the
        ! solver's HDF5 file.
        ! -------------------------------------------------------------------------

        CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    END SUBROUTINE run_sliceVertical

    MODULE SUBROUTINE dummy_run_sphere(self)
        ! -------------------------------------------------------------------------
        ! Default (no-op) implementation of run_sphere, used by solvers that do
        ! not support a spherical free-field observer geometry. Prints a
        ! "skipped" notice and switches do_sphere off, so the .dat summary does
        ! not report a step that never ran.
        ! -------------------------------------------------------------------------

        CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    END SUBROUTINE dummy_run_sphere

    MODULE SUBROUTINE run(self)
        ! -------------------------------------------------------------------------
        ! Config-driven entry point for the post-processing pipeline: runs, in a
        ! fixed order, every step whose do_<name> flag is enabled. The parser has
        ! already resolved do_all and switched off the steps that make no sense
        ! for the emitter (sphere and decay for a wind farm), so the same routine
        ! serves a single turbine, a wind farm or hand-made sources.
        ! -------------------------------------------------------------------------

        CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    END SUBROUTINE run

END INTERFACE


CONTAINS


! ===============================================================================
! HELPERS
! ===============================================================================
FUNCTION check_observers_distances(self, observers) RESULT(valid)
    ! -------------------------------------------------------------------------
    ! Flags the observers lying closer than mask_distance to any source node.
    ! At such distances the observer sits inside (or right on) the structure,
    ! where the point-dipole field is meaningless and nearly singular. Returns
    ! .true. for every valid observer; the run continues either way.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(IN)  :: self              ! [-] Solver instance
    REAL(WP)               , INTENT(IN)  :: observers(:,:)    ! [m] Observers coords (Nobs, 3)
    LOGICAL                , ALLOCATABLE :: valid(:)          ! [-] Observer mask (Nobs), .true. if far enough

    ! Local variables
    INTEGER(I32) :: i, j, n_nodes, n_obs
    REAL(WP)     :: dx, dy, dz, d2_mask


    n_nodes = size(self % sources % xyz, 1)
    n_obs   = size(observers, 1)
    d2_mask = self % mask_distance**2

    ALLOCATE(valid(n_obs))
    valid = .true.

    !$OMP PARALLEL DO PRIVATE(i, j, dx, dy, dz)
    do i = 1, n_obs
        do j = 1, n_nodes
            dx = self % sources % xyz(j, 1) - observers(i, 1)
            dy = self % sources % xyz(j, 2) - observers(i, 2)
            dz = self % sources % xyz(j, 3) - observers(i, 3)

            if (dx * dx + dy * dy + dz * dz < d2_mask) then
                valid(i) = .false.
                exit
            end if
        end do
    end do
    !$OMP END PARALLEL DO

END FUNCTION check_observers_distances


SUBROUTINE compute_pressure_masked(self, label, observers, block_size, total_pressure)
    ! -------------------------------------------------------------------------
    ! Wrapper of compute_pressure used by every run_* routine. Observers lying
    ! closer than mask_distance to any source node, or outside the water
    ! column (in_water), are not computed: their pressure is set to NaN (real
    ! and imaginary parts) and one warning per criterion reports how many were
    ! masked. The observer grid keeps its regular shape, so the masked points
    ! show up as holes when plotted. Without masked observers it is a plain
    ! call to compute_pressure (identical results).
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, IEEE_QUIET_NAN

    CLASS(AcousticSolver_t), INTENT(INOUT)            :: self                   ! [-] Solver instance
    CHARACTER(len=*)       , INTENT(IN)               :: label                  ! [-] Step name shown in the warning
    REAL(WP)               , INTENT(IN)               :: observers(:,:)         ! [m] Observers coords (Nobs, 3)
    INTEGER(I32)           , INTENT(IN)               :: block_size             ! [-] Observers per block
    COMPLEX(WP)            , INTENT(OUT), ALLOCATABLE :: total_pressure(:,:)    ! [Pa] Pressure (Nf, Nobs); masked: NaN

    ! Local variables
    LOGICAL     , ALLOCATABLE :: valid(:), wet(:)
    REAL(WP)    , ALLOCATABLE :: obs_valid(:,:)
    COMPLEX(WP) , ALLOCATABLE :: p_valid(:,:)
    INTEGER(I32)              :: n_obs, n_valid, n_close, n_dry, i, j
    REAL(WP)                  :: nan


    n_obs   = size(observers, 1)
    valid   = self % check_observers_distances(observers)
    wet     = self % in_water(observers)
    n_close = int(count(.not. valid), I32)    ! Masked by distance
    n_dry   = int(count(.not. wet), I32)      ! Masked by in_water
    valid   = valid .and. wet
    n_valid = int(count(valid), I32)

    ! No masked observer: plain call
    if (n_valid == n_obs) then
        CALL self % compute_pressure(observers, block_size, total_pressure)
        return
    end if

    if (n_close > 0_I32) &
        write(*, '(A,I0,A,I0,A,F4.1,A)') ' WARNING [' // label // '] ', n_close, ' of ', n_obs, &
                                         ' observers closer than ', self % mask_distance, &
                                         ' m to a source node: masked (NaN)'
    if (n_dry > 0_I32) &
        write(*, '(A,I0,A,I0,A)') ' WARNING [' // label // '] ', n_dry, ' of ', n_obs, &
                                  ' observers outside the water column: masked (NaN)'

    nan = ieee_value(1.0_WP, IEEE_QUIET_NAN)
    ALLOCATE(total_pressure(size(self % sources % Freqs), n_obs))
    total_pressure = cmplx(nan, nan, kind=WP)

    if (n_valid == 0_I32) return

    ! Compute only the valid observers, then scatter them back into the full grid
    ALLOCATE(obs_valid(n_valid, 3))
    j = 0_I32
    do i = 1, n_obs
        if (valid(i)) then
            j               = j + 1_I32
            obs_valid(j, :) = observers(i, :)
        end if
    end do

    CALL self % compute_pressure(obs_valid, block_size, p_valid)

    j = 0_I32
    do i = 1, n_obs
        if (valid(i)) then
            j                    = j + 1_I32
            total_pressure(:, i) = p_valid(:, j)
        end if
    end do

END SUBROUTINE compute_pressure_masked


FUNCTION pressure_norm(p) RESULT(pnorm)
    ! -------------------------------------------------------------------------
    ! Norm of |p| over the computed entries only (masked NaN observers are
    ! left out), used by the debug output of the run_* routines. Without NaN
    ! it is exactly norm2(abs(p)).
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_nan

    COMPLEX(WP), INTENT(IN) :: p(:,:)    ! [Pa] Pressure (Nfreqs, Nobs), possibly with NaN
    REAL(WP)                :: pnorm     ! [Pa] Norm of the computed pressures


    if (any(ieee_is_nan(real(p, WP)))) then
        pnorm = norm2(abs(pack(p, .not. ieee_is_nan(real(p, WP)))))
    else
        pnorm = norm2(abs(p))
    end if

END FUNCTION pressure_norm


! ===============================================================================
! WATER-COLUMN QUERIES
! Defaults: flat surface at z = 0 and flat seabed at z = -depth_base
! ===============================================================================
FUNCTION default_z_surface(self, x, y) RESULT(z)
    ! -------------------------------------------------------------------------
    ! Default water-column query: z of the water surface at (x, y), the plane
    ! z = 0.
    ! A method overrides z_surface, z_seabed and in_water to publish its own
    ! boundaries (a surface not at z = 0, a variable bathymetry, ...). The
    ! run_* steps never read depth_base directly: they only ask these
    ! queries, so they stay generic for every method.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(IN) :: self    ! [-] Solver instance
    REAL(WP)               , INTENT(IN) :: x       ! [m] x-coordinate
    REAL(WP)               , INTENT(IN) :: y       ! [m] y-coordinate
    REAL(WP)                            :: z       ! [m] z-coordinate of the water surface


    z = 0.0_WP

END FUNCTION default_z_surface


FUNCTION default_z_seabed(self, x, y) RESULT(z)
    ! -------------------------------------------------------------------------
    ! Default water-column query: z of the seabed at (x, y), the plane
    ! z = -depth_base.
    ! A method overrides z_surface, z_seabed and in_water to publish its own
    ! boundaries (a surface not at z = 0, a variable bathymetry, ...). The
    ! run_* steps never read depth_base directly: they only ask these
    ! queries, so they stay generic for every method.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(IN) :: self    ! [-] Solver instance
    REAL(WP)               , INTENT(IN) :: x       ! [m] x-coordinate
    REAL(WP)               , INTENT(IN) :: y       ! [m] y-coordinate
    REAL(WP)                            :: z       ! [m] z-coordinate of the seabed


    z = -self % env % depth_base

END FUNCTION default_z_seabed


FUNCTION default_in_water(self, observers) RESULT(valid)
    ! -------------------------------------------------------------------------
    ! Default water-column query: flags the observers inside the water
    ! column. Every observer is valid by default, so compute_pressure_masked
    ! never masks an observer for this reason.
    ! A method overrides z_surface, z_seabed and in_water to publish its own
    ! boundaries (a surface not at z = 0, a variable bathymetry, ...). The
    ! run_* steps never read depth_base directly: they only ask these
    ! queries, so they stay generic for every method.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(IN)  :: self              ! [-] Solver instance
    REAL(WP)               , INTENT(IN)  :: observers(:,:)    ! [m] Observers coords (Nobs, 3)
    LOGICAL                , ALLOCATABLE :: valid(:)          ! [-] Observer mask (Nobs), .true. inside the water


    ALLOCATE(valid(size(observers, 1)))
    valid = .true.

END FUNCTION default_in_water

END MODULE AcousticSolver
