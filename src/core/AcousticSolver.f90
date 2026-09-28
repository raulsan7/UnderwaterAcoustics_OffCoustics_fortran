! ===============================================================================
! MODULE: AcousticSolver
! ===============================================================================
! Defines AcousticSolver_t, the abstract parent type for all acoustic
! propagation solvers (ImageMethod_t and AnalyticalNormalModes_t). It holds
! the shared solver state (sources and their reference geometry,
! environment, save path, block size, debug/cluster flags) and implements
! the generic observer-geometry post-processing routines (spectrums, polar,
! cylinder, decay, line, sliceXY, sliceXZ, sliceVertical, sphere) shared by
! every concrete solver. The solver never knows what produced the sources
! (a turbine, a wind farm or hand-made sources): which steps make sense for
! each emitter is decided by the parser through the do_run flags.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE AcousticSolver

USE Kinds
USE Types
USE IOUtils
USE WindFarm, ONLY: WindFarm_t
USE WindTurbine, ONLY: WindTurbine_t, save_turbine_params, write_turbine_txt

IMPLICIT NONE

PRIVATE
PUBLIC :: AcousticSolver_t, save_variable, save_emitter_params, write_emitter_txt, pressure_norm


! ===============================================================================
! OBSERVER MASKING
! ===============================================================================
REAL(WP), PARAMETER :: MASK_DISTANCE = 5.0_WP    ! [m] Observers closer than this to a source node get NaN pressure


TYPE, ABSTRACT :: AcousticSolver_t
    ! -------------------------------------------------------------------------
    ! Parent class for all acoustic propagation solvers (ImageMethod_t and
    ! AnalyticalNormalModes_t). Holds the shared sources/environment/turbine
    ! state and implements the generic observer-geometry post-processing
    ! routines; init, compute_pressure and the parameter/path/summary
    ! routines are deferred.
    ! -------------------------------------------------------------------------

    TYPE(AcousticSource_t)     :: sources                 ! [-] Acoustic point sources and their reference geometry
    TYPE(Config_Environment_t) :: env                     ! [-] Environment configuration (depth, sound speed, ...)
    TYPE(Config_Computation_t) :: run_inputs              ! [-] Per-run_* post-processing step configuration
    INTEGER(I32)               :: block_size              ! [-] Number of observers processed simultaneously
    LOGICAL                    :: debug                   ! [-] Whether to print extra debug info (e.g. pressure norms)
    LOGICAL                    :: in_cluster              ! [-] Whether the solver runs on a cluster
    CHARACTER(len=50)          :: solver_name = 'None'    ! [-] Solver name
    CHARACTER(len=512)         :: save_path   = 'None'    ! [-] Path to save acoustics

    CONTAINS

    ! --- Public type-bound procedures --- !
    PROCEDURE :: check_observers_distances    ! Observer mask: .false. closer than MASK_DISTANCE to a source
    PROCEDURE :: compute_pressure_masked      ! compute_pressure with NaN at the masked observers
    PROCEDURE :: resolve_save_path            ! Output path shared by every solver (build_save_path helper)
    PROCEDURE :: update_parameters            ! Appends run_* parameters actually used to the .dat
    PROCEDURE :: run                          ! Config-driven entry point: every enabled do_* step
    PROCEDURE :: run_spectrums                ! Subroutines to compute acoustic pressure
    PROCEDURE :: run_polar
    PROCEDURE :: run_cylinder
    PROCEDURE :: run_decay
    PROCEDURE :: run_line
    PROCEDURE :: run_sliceXY
    PROCEDURE :: run_sliceXZ
    PROCEDURE :: run_sliceVertical
    PROCEDURE :: run_sphere => dummy_run_sphere

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

    SUBROUTINE save_parameters_i(self, turbine, farm)
        ! -------------------------------------------------------------------------
        ! Saves the solver, environment and turbine parameters to the HDF5
        ! output file and writes the .dat text summary.
        ! -------------------------------------------------------------------------

        IMPORT :: AcousticSolver_t, WindTurbine_t, WindFarm_t

        CLASS(AcousticSolver_t), INTENT(IN)           :: self       ! [-] Solver instance
        CLASS(WindTurbine_t)   , INTENT(IN), OPTIONAL :: turbine    ! [-] Turbine instance (optional)
        TYPE(WindFarm_t)       , INTENT(IN), OPTIONAL :: farm       ! [-] Wind farm instance (optional)

    END SUBROUTINE save_parameters_i

    SUBROUTINE build_save_path_i(self)
        ! -------------------------------------------------------------------------
        ! Builds the name of the HDF5 output file.
        ! -------------------------------------------------------------------------

        IMPORT :: AcousticSolver_t

        CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    END SUBROUTINE build_save_path_i

    SUBROUTINE write_txt_summary_i(self, turbine, farm)
        ! -------------------------------------------------------------------------
        ! Writes the human-readable .dat summary next to the HDF5 output file.
        ! -------------------------------------------------------------------------

        IMPORT :: AcousticSolver_t, WindTurbine_t, WindFarm_t

        CLASS(AcousticSolver_t), INTENT(IN)           :: self       ! [-] Solver instance
        CLASS(WindTurbine_t)   , INTENT(IN), OPTIONAL :: turbine    ! [-] Turbine instance (optional)
        TYPE(WindFarm_t)       , INTENT(IN), OPTIONAL :: farm       ! [-] Wind farm instance (optional)

    END SUBROUTINE write_txt_summary_i

END INTERFACE


CONTAINS


! ===============================================================================
! HELPERS
! ===============================================================================
FUNCTION check_observers_distances(self, observers) RESULT(valid)
    ! -------------------------------------------------------------------------
    ! Flags the observers lying closer than MASK_DISTANCE to any source node.
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
    d2_mask = MASK_DISTANCE**2

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
    ! closer than MASK_DISTANCE to any source node are not computed: their
    ! pressure is set to NaN (real and imaginary parts) and a warning reports
    ! how many were masked. The observer grid keeps its regular shape, so the
    ! masked points show up as holes when plotted. Without masked observers
    ! it is a plain call to compute_pressure (identical results).
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, IEEE_QUIET_NAN

    CLASS(AcousticSolver_t), INTENT(INOUT)            :: self                   ! [-] Solver instance
    CHARACTER(len=*)       , INTENT(IN)               :: label                  ! [-] Step name shown in the warning
    REAL(WP)               , INTENT(IN)               :: observers(:,:)         ! [m] Observers coords (Nobs, 3)
    INTEGER(I32)           , INTENT(IN)               :: block_size             ! [-] Observers per block
    COMPLEX(WP)            , INTENT(OUT), ALLOCATABLE :: total_pressure(:,:)    ! [Pa] Pressure (Nf, Nobs); masked: NaN

    ! Local variables
    LOGICAL     , ALLOCATABLE :: valid(:)
    REAL(WP)    , ALLOCATABLE :: obs_valid(:,:)
    COMPLEX(WP) , ALLOCATABLE :: p_valid(:,:)
    INTEGER(I32)              :: n_obs, n_valid, i, j
    REAL(WP)                  :: nan


    n_obs   = size(observers, 1)
    valid   = self % check_observers_distances(observers)
    n_valid = int(count(valid), I32)

    ! No masked observer: plain call
    if (n_valid == n_obs) then
        CALL self % compute_pressure(observers, block_size, total_pressure)
        return
    end if

    write(*, '(A,I0,A,I0,A,F4.1,A)') ' WARNING [' // label // '] ', n_obs - n_valid, ' of ', n_obs, &
                                     ' observers closer than ', MASK_DISTANCE, ' m to a source node: masked (NaN)'

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


SUBROUTINE resolve_save_path(self, method_tag)
    ! -------------------------------------------------------------------------
    ! Resolves the final HDF5 output path, shared by every solver's
    ! build_save_path. If the configured save_path is "None"/"Default", the
    ! name is built from the emitter labels carried by the sources around the
    ! solver's method tag: out_dir // name_tag // method_tag // name_suffix,
    ! e.g. ./turbine_acoustic_data/plot_mn_ANM24_ws11.4.hdf5 or
    ! ./farm_acoustic_data/Farm_2_DTU10MN_ANM24_ws8.0-11.4.hdf5. Otherwise the
    ! .hdf5 extension is appended to the given path. The destination folder
    ! is created if needed and a previous file with the same name is deleted,
    ! since save_variable only overwrites the datasets it writes and datasets
    ! of an older run would otherwise survive in the new file.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self          ! [-] Solver instance
    CHARACTER(len=*)       , INTENT(IN)    :: method_tag    ! [-] Method and resolution tag, e.g. "_ANM24"

    ! Local variables
    CHARACTER(len=512) :: path_upper
    INTEGER(I32)       :: u, ierr


    path_upper = to_upper(trim(self % save_path))

    if (trim(path_upper) == 'NONE' .or. trim(path_upper) == 'DEFAULT') then
        associate(g => self % sources % geometry)
            self % save_path = trim(g % out_dir) // trim(g % name_tag) // method_tag // trim(g % name_suffix) // ".hdf5"
        end associate
    else
        self % save_path = trim(self % save_path) // ".hdf5"
    end if

    ! Make sure the destination folder exists before anything tries to write into it
    CALL ensure_parent_dir(trim(self % save_path))

    ! Start from an empty file, as the .dat summary does (status='replace')
    open(newunit=u, file=trim(self % save_path), status='old', iostat=ierr)
    if (ierr == 0_I32) then
        close(u, status='delete')
        write(*,*) '   ---> Previous output file replaced: ', trim(self % save_path)
    end if

END SUBROUTINE resolve_save_path


SUBROUTINE update_parameters(self)
    ! -------------------------------------------------------------------------
    ! Appends to the .dat summary already written by save_parameters() one
    ! block per post-processing step that actually ran (per do_run), with the
    ! resolved parameters run() used for it, then closes the file with its
    ! final separator. Must be called after run(), since every step only
    ! resolves its own DEFAULT_R/DEFAULT_I sentinels when it actually
    ! executes.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(IN) :: self    ! [-] Solver instance

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: dat_path
    INTEGER(I32)                  :: u, i


    dat_path = swap_extension(trim(self % save_path), ".dat")

    open(newunit=u, file=dat_path, status='old', position='append', action='write')

    write(u, '(A)') repeat('-', 70)
    write(u, '(A)') 'POST-PROCESSING STEPS EXECUTED'
    write(u, '(A)') repeat('-', 70)

    associate(d => self % run_inputs % do_run, pp => self % run_inputs)

        if (d % do_spectrums) then
            write(u, '(A)') ''
            write(u, '(A)') 'Spectrums:'
            if (.not. ALLOCATED(pp % params_spectrums % observers_x)) then
                write(u, '(A)') '   observers : 3 default points (auto, x=10/250/500 m, z=-Depth/2)'
            else
                write(u, '(A,I0)') '   n_observers : ', pp % params_spectrums % n_observers
                do i = 1, pp % params_spectrums % n_observers
                    write(u, '(A,I0,A,F8.2,A,F8.2,A,F8.2,A)') '   obs ', i, ' (x,y,z) : ', &
                                                              pp % params_spectrums % observers_x(i), ', ', &
                                                              pp % params_spectrums % observers_y(i), ', ', &
                                                              pp % params_spectrums % observers_z(i), ' m'
                end do
            end if
        end if

        if (d % do_polar) then
            write(u, '(A)') ''
            write(u, '(A)')        'Polar:'
            write(u, '(A,F8.2,A)') '   r_polar      : ', pp % params_polar % r_polar, ' m'
            write(u, '(A,F8.2,A)') '   z_polar      : ', pp % params_polar % z_polar, ' m'
            write(u, '(A,I0)')     '   n_polar      : ', pp % params_polar % n_polar
            write(u, '(A,2F10.3)') '   center_polar : ', pp % params_polar % center_polar
        end if

        if (d % do_cylinder) then
            write(u, '(A)') ''
            write(u, '(A)')        'Cylinder:'
            write(u, '(A,F8.2,A)') '   r_cyl      : ', pp % params_cylinder % r_cyl, ' m'
            write(u, '(A,I0)')     '   n_cyl      : ', pp % params_cylinder % n_cyl
            write(u, '(A,I0)')     '   nz_cyl     : ', pp % params_cylinder % nz_cyl
            write(u, '(A,2F10.3)') '   center_cyl : ', pp % params_cylinder % center_cyl
        end if

        if (d % do_decay) then
            write(u, '(A)') ''
            write(u, '(A)')        'Decay:'
            write(u, '(A,2F10.3)') '   distance_decay : ', pp % params_decay % distance_decay
            write(u, '(A,I0)')     '   n_decay        : ', pp % params_decay % n_decay
            write(u, '(A,F8.2,A)') '   z_decay        : ', pp % params_decay % z_decay, ' m'
            write(u, '(A,L1)')     '   logspace_decay : ', pp % params_decay % logspace_decay
        end if

        if (d % do_line) then
            write(u, '(A)') ''
            write(u, '(A)')        'Line:'
            write(u, '(A,3F10.3)') '   p1_line       : ', pp % params_line % p1_line
            write(u, '(A,3F10.3)') '   p2_line       : ', pp % params_line % p2_line
            write(u, '(A,I0)')     '   n_line        : ', pp % params_line % n_line
            write(u, '(A,L1)')     '   logspace_line : ', pp % params_line % logspace_line
        end if

        if (d % do_sliceXY) then
            write(u, '(A)') ''
            write(u, '(A)')        'SliceXY:'
            write(u, '(A,F8.2,A)') '   z_xy       : ', pp % params_sliceXY % z_xy, ' m'
            write(u, '(A,I0)')     '   nx_xy      : ', pp % params_sliceXY % nx_xy
            write(u, '(A,I0)')     '   ny_xy      : ', pp % params_sliceXY % ny_xy
            write(u, '(A,2F10.3)') '   xlim_xy    : ', pp % params_sliceXY % xlim_xy
            write(u, '(A,2F10.3)') '   ylim_xy    : ', pp % params_sliceXY % ylim_xy
            write(u, '(A,2F10.3)') '   center_xy  : ', pp % params_sliceXY % center_xy
        end if

        if (d % do_sliceXZ) then
            write(u, '(A)') ''
            write(u, '(A)')        'SliceXZ:'
            write(u, '(A,F8.2,A)') '   y_xz       : ', pp % params_sliceXZ % y_xz, ' m'
            write(u, '(A,I0)')     '   nx_xz      : ', pp % params_sliceXZ % nx_xz
            write(u, '(A,I0)')     '   nz_xz      : ', pp % params_sliceXZ % nz_xz
            write(u, '(A,2F10.3)') '   xlim_xz    : ', pp % params_sliceXZ % xlim_xz
            write(u, '(A,2F10.3)') '   zlim_xz    : ', pp % params_sliceXZ % zlim_xz
        end if

        if (d % do_sliceVertical) then
            write(u, '(A)') ''
            write(u, '(A)')         'SliceVertical:'
            write(u, '(A,F8.2,A)')  '   azimuth_sv : ', pp % params_sliceVertical % azimuth_sv, ' deg'
            write(u, '(A,2F10.3)')  '   center_sv  : ', pp % params_sliceVertical % center_sv
            write(u, '(A,F10.3,A)') '   width_sv   : ', pp % params_sliceVertical % width_sv, ' m'
            write(u, '(A,2F10.3)')  '   ulim_sv    : ', pp % params_sliceVertical % ulim_sv
            write(u, '(A,2F10.3)')  '   zlim_sv    : ', pp % params_sliceVertical % zlim_sv
            write(u, '(A,I0)')      '   nu_sv      : ', pp % params_sliceVertical % nu_sv
            write(u, '(A,I0)')      '   nz_sv      : ', pp % params_sliceVertical % nz_sv
            write(u, '(A,F8.2,A)')  '   margin_sv  : ', pp % params_sliceVertical % margin_sv, ' m'
        end if

        if (d % do_sphere) then
            write(u, '(A)') ''
            write(u, '(A)')        'Sphere:'
            write(u, '(A,F8.2,A)') '   r_sphere       : ', pp % params_sphere % r_sphere, ' m'
            write(u, '(A,I0)')     '   n_theta_sphere : ', pp % params_sphere % n_theta_sphere
            write(u, '(A,I0)')     '   nz_sphere      : ', pp % params_sphere % nz_sphere
            write(u, '(A,3F10.3)') '   center_sphere  : ', pp % params_sphere % center_sphere
        end if

    end associate

    write(u, *)
    write(u, '(A)') repeat('=', 70)

    close(u)

    write(*,*) '   ---> Parameters updated at ', dat_path

END SUBROUTINE update_parameters


SUBROUTINE save_emitter_params(save_path, turbine, farm)
    ! -------------------------------------------------------------------------
    ! Saves the parameters of whatever produced the sources into the HDF5
    ! output file. Nothing is saved for hand-made sources.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*)    , INTENT(IN)           :: save_path    ! [-] Path to the HDF5 output file
    CLASS(WindTurbine_t), INTENT(IN), OPTIONAL :: turbine      ! [-] Turbine instance (optional)
    TYPE(WindFarm_t)    , INTENT(IN), OPTIONAL :: farm         ! [-] Wind farm instance (optional)


    if (present(turbine)) then
        CALL save_turbine_params(turbine, save_path)
        CALL turbine % save_specific_params(save_path)
    end if

    if (present(farm)) CALL farm % save_params(save_path)

END SUBROUTINE save_emitter_params


SUBROUTINE write_emitter_txt(unit, turbine, farm)
    ! -------------------------------------------------------------------------
    ! Writes the block of whatever produced the sources into the .dat
    ! summary file. Nothing is written for hand-made sources.
    ! -------------------------------------------------------------------------

    INTEGER(I32)        , INTENT(IN)           :: unit       ! [-] Unit of the open .dat file
    CLASS(WindTurbine_t), INTENT(IN), OPTIONAL :: turbine    ! [-] Turbine instance (optional)
    TYPE(WindFarm_t)    , INTENT(IN), OPTIONAL :: farm       ! [-] Wind farm instance (optional)


    if (present(turbine)) CALL write_turbine_txt(turbine, unit)
    if (present(farm)) CALL farm % write_txt(unit)

END SUBROUTINE write_emitter_txt


! ===============================================================================
! RUN SUBROUTINES
! ===============================================================================
SUBROUTINE run_spectrums(self)
    ! -------------------------------------------------------------------------
    ! Computes the pressure spectrum at a small set of fixed observer points
    ! (default: 3 points along +x at mid-depth). Saves the pressure, observer
    ! coordinates and frequency vector to the solver's HDF5 file.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: p(:,:)
    REAL(WP)    , ALLOCATABLE :: observers_(:,:)
    INTEGER(I32)              :: bs


    write(*, '(A)') '>> [Spectrums] computing pressure at observer points...'

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_spectrums)
        if (pp % n_observers == DEFAULT_I) then
            ALLOCATE(observers_(3, 3))
            observers_(:, 1) = [10.0_WP, 250.0_WP, 500.0_WP]
            observers_(:, 2) = 0.0_WP
            observers_(:, 3) = -self % env % depth_base  / 2.0_WP
        else
            ALLOCATE(observers_(pp % n_observers, 3))
            observers_(:, 1) = pp % observers_x
            observers_(:, 2) = pp % observers_y
            observers_(:, 3) = pp % observers_z
        end if
    end associate

    bs = self % block_size

    if (any(observers_(:, 3) > 0.0_WP)) &
        write(*,*) "WARNING AcousticSolver % run_spectrums: there are observers placed over z = 0 m"

    ! Compute pressure
    CALL self % compute_pressure_masked('Spectrums', observers_, bs, p)

    if (self % debug) write(*, '(A,F10.4)') '   Pressure norm : ', pressure_norm(p)

    ! Save data
    CALL save_variable(self % save_path, "P_spectrums", p)
    CALL save_variable(self % save_path, "Obs_spectrums", observers_)
    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)

    write(*, '(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_spectrums


SUBROUTINE run_polar(self)
    ! -------------------------------------------------------------------------
    ! Computes pressure on a horizontal circle of radius r_polar around the
    ! turbine barycenter, at fixed depth z_polar. Saves the pressure plus the
    ! polar geometry (radius, depth, azimuth angles, observer coordinates) to
    ! the solver's HDF5 file.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: p(:,:)
    REAL(WP)    , ALLOCATABLE :: observers(:,:), theta_deg_polar(:)
    REAL(WP)                  :: r, z, center(2), theta_deg, theta_rad, d_theta, Depth
    INTEGER(I32)              :: n_theta, bs, i


    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_polar)
        if (pp % r_polar == DEFAULT_R) pp % r_polar = 500.0_WP
        if (pp % z_polar == DEFAULT_R) pp % z_polar = -Depth / 2.0_WP
        if (pp % n_polar == DEFAULT_I) pp % n_polar = 72_I32
        if (all(pp % center_polar == DEFAULT_R)) &    ! Not a logical .and.
            pp % center_polar = self % sources % geometry % center

        r       = pp % r_polar
        z       = pp % z_polar
        n_theta = pp % n_polar
        center  = pp % center_polar
    end associate

    ! Build observers array
    ALLOCATE(observers(n_theta, 3), theta_deg_polar(n_theta))
    d_theta = 360.0_WP / real(n_theta, WP)

    do i = 1, n_theta
        theta_deg          = real(i - 1_I32, WP) * d_theta
        theta_rad          = theta_deg * PI / 180.0_WP
        observers(i, 1)    = r * cos(theta_rad) + center(1)
        observers(i, 2)    = r * sin(theta_rad) + center(2)
        observers(i, 3)    = z
        theta_deg_polar(i) = theta_deg
    end do

    write(*, '(A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,I0)') &
        '>> [Polar]    r=', r, ' m | center=(', center(1), ',', center(2), &
        ') m | z=', z, ' m | observers=', n_theta

    ! Compute pressure
    CALL self % compute_pressure_masked('Polar', observers, bs, p)
    if (self % debug) write(*, '(A,F10.4)') '   Pressure norm : ', pressure_norm(p)

    ! Save data
    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_polar", p)
    CALL save_variable(self % save_path, "R_polar", r)
    CALL save_variable(self % save_path, "Z_polar", z)
    CALL save_variable(self % save_path, "Theta_deg_polar", theta_deg_polar)
    CALL save_variable(self % save_path, "Obs_polar", observers)
    CALL save_variable(self % save_path, "Center_polar", center)

    write(*, '(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_polar


SUBROUTINE run_cylinder(self)
    ! -------------------------------------------------------------------------
    ! Computes pressure on a vertical cylindrical surface (radius r_cyl,
    ! n_cyl azimuth points x nz_cyl depth levels) around the turbine
    ! barycenter. Saves the reshaped pressure/observer grids and the cylinder
    ! geometry metadata to the solver's HDF5 file.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: p(:,:)
    REAL(WP)    , ALLOCATABLE :: observers(:,:)
    REAL(WP)                  :: r, center(2), theta_rad, d_theta, z_val, dz, Depth
    INTEGER(I32)              :: n_theta, nz, bs, i, j, idx, Nobs

    REAL(WP)   , ALLOCATABLE :: theta_deg_arr(:)
    REAL(WP)   , ALLOCATABLE :: z_cylinder(:)
    REAL(WP)   , ALLOCATABLE :: obs_reshaped(:,:,:)
    COMPLEX(WP), ALLOCATABLE :: p_reshaped(:,:,:)
    REAL(WP)                 :: dA


    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_cylinder)
        if (pp % r_cyl == DEFAULT_R)  pp % r_cyl  = 500.0_WP
        if (pp % n_cyl == DEFAULT_I)  pp % n_cyl  = 72_I32
        if (pp % nz_cyl == DEFAULT_I) pp % nz_cyl = 20_I32
        if (any(pp % center_cyl == DEFAULT_R)) &
            pp % center_cyl = self % sources % geometry % center

        r       = pp % r_cyl
        n_theta = pp % n_cyl
        nz      = pp % nz_cyl
        center  = pp % center_cyl
    end associate

    Nobs = n_theta * nz
    ALLOCATE(observers(Nobs, 3))

    d_theta = 360.0_WP / real(n_theta, WP)
    dz      = (0.0_WP - (-Depth)) / max(real(nz - 1_I32, WP), 1.0_WP)

    ALLOCATE(theta_deg_arr(n_theta))
    ALLOCATE(z_cylinder(nz))

    do j = 1, nz
        z_cylinder(j) = -Depth + real(j - 1_I32, WP) * dz
    end do

    idx = 1_I32
    do i = 1, n_theta
        theta_deg_arr(i) = real(i - 1_I32, WP) * d_theta
        theta_rad        = theta_deg_arr(i) * PI / 180.0_WP
        do j = 1, nz
            z_val             = z_cylinder(j)
            observers(idx, 1) = r * cos(theta_rad) + center(1)
            observers(idx, 2) = r * sin(theta_rad) + center(2)
            observers(idx, 3) = z_val
            idx               = idx + 1_I32
        end do
    end do

    write(*, '(A,F8.2,A,F8.2,A,F8.2,A,I0,A,I0,A,I0)') &
        '>> [Cylinder]   r=', r, ' m | center=(', center(1), ',', center(2), &
        ') m | grid=', n_theta, ' x ', nz, ' = ', Nobs

    CALL self % compute_pressure_masked('Cylinder', observers, bs, p)
    if (self % debug) write(*, '(A,F10.4)') '   Pressure norm : ', pressure_norm(p)

    ALLOCATE(p_reshaped(size(p, 1), n_theta, nz))
    ALLOCATE(obs_reshaped(n_theta, nz, 3))

    idx = 1_I32
    do i = 1, n_theta
        do j = 1, nz
            obs_reshaped(i, j, :) = observers(idx, :)
            p_reshaped(:, i, j)   = p(:, idx)
            idx                   = idx + 1_I32
        end do
    end do

    dA = r * dz * (2.0_WP * PI / real(n_theta, WP))

    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_cylinder", p_reshaped)
    CALL save_variable(self % save_path, "R_cylinder", r)
    CALL save_variable(self % save_path, "Z_cylinder", z_cylinder)
    CALL save_variable(self % save_path, "Theta_deg_cylinder", theta_deg_arr)
    CALL save_variable(self % save_path, "Obs_cylinder", obs_reshaped)
    CALL save_variable(self % save_path, "Center_cylinder", center)
    CALL save_variable(self % save_path, "dA_cylinder", dA)

    write(*, '(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_cylinder


SUBROUTINE run_decay(self)
    ! -------------------------------------------------------------------------
    ! Computes pressure along a straight line following the wind direction
    ! from the turbine axis, at fixed depth z_decay, with linear or log-spaced
    ! distances. Saves the pressure and the observer/distance geometry to the
    ! solver's HDF5 file.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: p(:,:)
    REAL(WP)    , ALLOCATABLE :: observers(:,:)
    REAL(WP)    , ALLOCATABLE :: distances(:)
    REAL(WP)                  :: distance(2), z, wind_unit(2), d_val, WindDir_rad
    REAL(WP)                  :: d_log1, d_log2, d_step, linear_step, Depth
    INTEGER(I32)              :: n_points, bs, i
    LOGICAL                   :: logspace
    INTEGER(I32)              :: logspace_int


    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_decay)
        if (pp % z_decay == DEFAULT_R)             pp % z_decay        = -Depth / 2.0_WP
        if (any(pp % distance_decay == DEFAULT_R)) pp % distance_decay = [10.0_WP, 500.0_WP]
        if (pp % n_decay == DEFAULT_I)             pp % n_decay        = 200_I32
        ! logspace_decay has no sentinel: it is read directly as a LOGICAL from the .nml

        z        = pp % z_decay
        distance = pp % distance_decay
        n_points = pp % n_decay
        logspace = pp % logspace_decay
    end associate

    WindDir_rad = self % sources % geometry % direction * PI / 180.0_WP
    wind_unit   = [cos(WindDir_rad), sin(WindDir_rad)]

    ALLOCATE(observers(n_points, 3))
    ALLOCATE(distances(n_points))

    if (distance(1) <= 0.0_WP) then
        write(*,*) "WARNING AcousticSolver % run_decay(): first point can't be placed at r = 0 m"
        distance(1) = 1.0_WP
    end if

    if (logspace) then
        d_log1 = log10(distance(1))
        d_log2 = log10(distance(2))
        d_step = (d_log2 - d_log1) / real(n_points - 1_I32, WP)
        do i = 1, n_points
            d_val             = 10.0_WP**(d_log1 + real(i - 1_I32, WP) * d_step)
            distances(i)      = d_val
            observers(i, 1:2) = self % sources % geometry % origin + d_val * wind_unit
            observers(i, 3)   = z
        end do
    else
        linear_step = (distance(2) - distance(1)) / real(n_points - 1_I32, WP)
        do i = 1, n_points
            d_val             = distance(1) + real(i - 1_I32, WP) * linear_step
            distances(i)      = d_val
            observers(i, 1:2) = self % sources % geometry % origin + d_val * wind_unit
            observers(i, 3)   = z
        end do
    end if

    write(*, '(A,I0,A,F8.2,A,F8.2,A)') &
        '>> [Decay]      points=', n_points, ' | direction=', &
        self % sources % geometry % direction, ' deg | z=', z, ' m'

    CALL self % compute_pressure_masked('Decay', observers, bs, p)
    if (self % debug) write(*, '(A,F10.4)') '   Pressure norm : ', pressure_norm(p)

    logspace_int = merge(1_I32, 0_I32, logspace)

    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_decay", p)
    CALL save_variable(self % save_path, "Z_decay", z)
    CALL save_variable(self % save_path, "Distance_decay", distance)
    CALL save_variable(self % save_path, "Distances_decay", distances)
    CALL save_variable(self % save_path, "Obs_decay", observers)
    CALL save_variable(self % save_path, "Logspace_decay", logspace_int)

    write(*, '(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_decay


SUBROUTINE run_line(self)
    ! -------------------------------------------------------------------------
    ! Computes pressure along an arbitrary straight line between two points
    ! p1_line and p2_line, with linear or log-spaced distances. Saves the
    ! pressure and the observer/distance geometry to the solver's HDF5 file.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: p(:,:)
    REAL(WP)    , ALLOCATABLE :: observers(:,:)
    REAL(WP)    , ALLOCATABLE :: distances(:)
    REAL(WP)                  :: p1(3), p2(3), dir_vec(3), dist_total, d_val
    REAL(WP)                  :: d_log1, d_log2, d_step, linear_step, Depth
    INTEGER(I32)              :: n_points, i, bs
    LOGICAL                   :: logspace
    INTEGER(I32)              :: logspace_int


    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_line)
        ! Per coordinate: only the components left at -999 take the default point's value
        WHERE (pp % p1_line == DEFAULT_R) pp % p1_line = [100.0_WP, 0.0_WP, 0.0_WP]
        WHERE (pp % p2_line == DEFAULT_R) pp % p2_line = [100.0_WP, 0.0_WP, -Depth]
        ! logspace_line has no sentinel: it is read directly as a LOGICAL from the .nml

        p1       = pp % p1_line
        p2       = pp % p2_line
        logspace = pp % logspace_line
    end associate

    dist_total = norm2(p2 - p1)
    if (dist_total == 0.0_WP) error stop "AcousticSolver.run_line(): points p1 and p2 are coincident."

    ! n_line depends on dist_total, so it is resolved separately
    if (self % run_inputs % params_line % n_line == DEFAULT_I) then
        if (dist_total > 200.0_WP) then
            self % run_inputs % params_line % n_line = int(dist_total / 5.0_WP, I32)
        else
            self % run_inputs % params_line % n_line = int(dist_total, I32)
        end if
    end if
    n_points = self % run_inputs % params_line % n_line

    ALLOCATE(observers(n_points, 3))
    ALLOCATE(distances(n_points))
    dir_vec = (p2 - p1) / dist_total

    if (logspace) then
        d_log1 = log10(dist_total * 1.0e-6_WP)
        d_log2 = log10(dist_total)
        d_step = (d_log2 - d_log1) / real(n_points - 1_I32, WP)
        do i = 1, n_points
            d_val           = 10.0_WP**(d_log1 + real(i - 1_I32, WP) * d_step)
            distances(i)    = d_val
            observers(i, :) = p1 + (d_val / dist_total) * dir_vec
        end do
    else
        linear_step = dist_total / real(n_points - 1_I32, WP)
        do i = 1, n_points
            d_val           = real(i - 1_I32, WP) * linear_step
            distances(i)    = d_val
            observers(i, :) = p1 + (d_val / dist_total) * dir_vec
        end do
    end if

    write(*, '(A,I0,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A)') &
        '>> [Line]       points=', n_points, ' | p1=(', p1(1), ',', p1(2), ',', p1(3), &
        ') m | p2=(', p2(1), ',', p2(2), ',', p2(3), ') m'

    CALL self % compute_pressure_masked('Line', observers, bs, p)
    if (self % debug) write(*, '(A,F10.4)') '   Pressure norm : ', pressure_norm(p)

    logspace_int = merge(1_I32, 0_I32, logspace)

    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_line", p)
    CALL save_variable(self % save_path, "P1_line", p1)
    CALL save_variable(self % save_path, "P2_line", p2)
    CALL save_variable(self % save_path, "Obs_line", observers)
    CALL save_variable(self % save_path, "Distance_line", dist_total)
    CALL save_variable(self % save_path, "Distances_line", distances)
    CALL save_variable(self % save_path, "Logspace_line", logspace_int)

    write(*, '(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_line


SUBROUTINE run_sliceXY(self)
    ! -------------------------------------------------------------------------
    ! Computes pressure on a horizontal (x,y) grid at fixed depth z_xy, over
    ! the given x/y limits. Saves the reshaped pressure/observer grids and the
    ! slice geometry metadata to the solver's HDF5 file.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: p(:,:)
    COMPLEX(WP) , ALLOCATABLE :: p_reshaped(:,:,:)
    REAL(WP)    , ALLOCATABLE :: observers(:,:)
    REAL(WP)    , ALLOCATABLE :: obs_reshaped(:,:,:)
    REAL(WP)    , ALLOCATABLE :: xs(:), ys(:)
    REAL(WP)                  :: z, xlim(2), ylim(2), center(2), Depth
    REAL(WP)                  :: mid_x, mid_y, dx, dy, x_val, y_val
    INTEGER(I32)              :: nx, ny, i, j, idx, bs
    REAL(WP)    , PARAMETER   :: tol = 1.0e-5_WP


    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_sliceXY)
        if (pp % z_xy == DEFAULT_R) pp % z_xy = -Depth / 2.0_WP
        if (pp % z_xy > 0.0_WP) then
            write(*,*) "AcousticSolver.run_sliceXY(): z must be <= 0. Switching to -z"
            pp % z_xy = -pp % z_xy
        end if

        if (pp % nx_xy == DEFAULT_I)        pp % nx_xy   = 26_I32
        if (pp % ny_xy == DEFAULT_I)        pp % ny_xy   = 26_I32
        if (any(pp % xlim_xy == DEFAULT_R)) pp % xlim_xy = [-500.0_WP, 500.0_WP]
        if (any(pp % ylim_xy == DEFAULT_R)) pp % ylim_xy = [-500.0_WP, 500.0_WP]

        z    = pp % z_xy
        nx   = pp % nx_xy
        ny   = pp % ny_xy
        xlim = pp % xlim_xy
        ylim = pp % ylim_xy
    end associate

    mid_x = 0.5_WP * (xlim(1) + xlim(2))
    mid_y = 0.5_WP * (ylim(1) + ylim(2))

    associate(pp => self % run_inputs % params_sliceXY)
        if (any(pp % center_xy == DEFAULT_R)) then
            pp % center_xy = [mid_x, mid_y]
        else if (abs(pp % center_xy(1) - mid_x) > tol .or. abs(pp % center_xy(2) - mid_y) > tol) then
            write(*,*) "AcousticSolver.run_sliceXY(): provided center does not match midpoint. " // &
                       "Switching to proper center."
            pp % center_xy = [mid_x, mid_y]
        end if
        center = pp % center_xy
    end associate

    ALLOCATE(observers(nx*ny, 3))
    ALLOCATE(xs(nx), ys(ny))
    dx = (xlim(2) - xlim(1)) / max(real(nx - 1_I32, WP), 1.0_WP)
    dy = (ylim(2) - ylim(1)) / max(real(ny - 1_I32, WP), 1.0_WP)

    do j = 1, nx
        xs(j) = xlim(1) + real(j - 1_I32, WP) * dx
    end do
    do i = 1, ny
        ys(i) = ylim(1) + real(i - 1_I32, WP) * dy
    end do

    idx = 1_I32
    do i = 1, ny
        y_val = ys(i)
        do j = 1, nx
            x_val             = xs(j)
            observers(idx, 1) = x_val
            observers(idx, 2) = y_val
            observers(idx, 3) = z
            idx               = idx + 1_I32
        end do
    end do

    write(*, '(A,I0,A,I0,A,F8.2,A,F8.2,A,F8.2,A)') &
        '>> [SliceXY]    nx=', nx, ' | ny=', ny, ' | z=', z, &
        ' m | center=(', center(1), ',', center(2), ') m'

    CALL self % compute_pressure_masked('SliceXY', observers, bs, p)
    if (self % debug) write(*, '(A,F10.4)') '   Pressure norm : ', pressure_norm(p)

    ALLOCATE(p_reshaped(size(p, 1), nx, ny))
    ALLOCATE(obs_reshaped(nx, ny, 3))

    idx = 1_I32
    do i = 1, ny
        do j = 1, nx
            obs_reshaped(j, i, :) = observers(idx, :)
            p_reshaped(:, j, i)   = p(:, idx)
            idx                   = idx + 1_I32
        end do
    end do

    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_slicexy", p_reshaped)
    CALL save_variable(self % save_path, "Obs_slicexy", obs_reshaped)
    CALL save_variable(self % save_path, "X_slicexy", xs)
    CALL save_variable(self % save_path, "Y_slicexy", ys)
    CALL save_variable(self % save_path, "Z_slicexy", z)
    CALL save_variable(self % save_path, "Center_slicexy", center)
    CALL save_variable(self % save_path, "Nx_slicexy", nx)
    CALL save_variable(self % save_path, "Ny_slicexy", ny)

    write(*, '(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_sliceXY


SUBROUTINE run_sliceXZ(self)
    ! -------------------------------------------------------------------------
    ! Computes pressure on a vertical (x,z) grid at fixed y_xz, over the given
    ! x/z limits. Saves the reshaped pressure/observer grids and the slice
    ! geometry metadata to the solver's HDF5 file.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: p(:,:)
    COMPLEX(WP) , ALLOCATABLE :: p_reshaped(:,:,:)
    REAL(WP)    , ALLOCATABLE :: observers(:,:)
    REAL(WP)    , ALLOCATABLE :: obs_reshaped(:,:,:)
    REAL(WP)    , ALLOCATABLE :: xs(:), zs(:)
    REAL(WP)                  :: y, xlim(2), zlim(2), dx, dz, x_val, z_val, Depth
    INTEGER(I32)              :: nx, nz, i, j, idx, bs


    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_sliceXZ)
        if (pp % y_xz  == DEFAULT_R)        pp % y_xz    = self % sources % geometry % origin(2)
        if (pp % nx_xz == DEFAULT_I)        pp % nx_xz   = 26_I32
        if (pp % nz_xz == DEFAULT_I)        pp % nz_xz   = 26_I32
        if (any(pp % xlim_xz == DEFAULT_R)) pp % xlim_xz = [-500.0_WP, 500.0_WP]
        WHERE (pp % zlim_xz == DEFAULT_R)   pp % zlim_xz = [0.0_WP, -Depth]    ! Per limit, as documented

        y    = pp % y_xz
        nx   = pp % nx_xz
        nz   = pp % nz_xz
        xlim = pp % xlim_xz
        zlim = pp % zlim_xz
    end associate

    ALLOCATE(observers(nx*nz, 3))
    ALLOCATE(xs(nx), zs(nz))
    dx = (xlim(2) - xlim(1)) / max(real(nx - 1_I32, WP), 1.0_WP)
    dz = (zlim(2) - zlim(1)) / max(real(nz - 1_I32, WP), 1.0_WP)

    do j = 1, nx
        xs(j) = xlim(1) + real(j - 1_I32, WP) * dx
    end do
    do i = 1, nz
        zs(i) = zlim(1) + real(i - 1_I32, WP) * dz
    end do

    idx = 1_I32
    do i = 1, nz
        z_val = zs(i)
        do j = 1, nx
            x_val             = xs(j)
            observers(idx, 1) = x_val
            observers(idx, 2) = y
            observers(idx, 3) = z_val
            idx               = idx + 1_I32
        end do
    end do

    write(*, '(A,I0,A,I0,A,F8.2,A,F8.2,A,F8.2,A)') &
        '>> [SliceXZ]    nx=', nx, ' | nz=', nz, ' | y=', y, &
        ' m | z from ', zlim(1), ' to ', zlim(2), ' m'

    CALL self % compute_pressure_masked('SliceXZ', observers, bs, p)
    if (self % debug) write(*, '(A,F10.4)') '   Pressure norm : ', pressure_norm(p)

    ALLOCATE(p_reshaped(size(p, 1), nx, nz))
    ALLOCATE(obs_reshaped(nx, nz, 3))

    idx = 1_I32
    do i = 1, nz
        do j = 1, nx
            obs_reshaped(j, i, :) = observers(idx, :)
            p_reshaped(:, j, i)   = p(:, idx)
            idx                   = idx + 1_I32
        end do
    end do

    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_slicexz", p_reshaped)
    CALL save_variable(self % save_path, "Obs_slicexz", obs_reshaped)
    CALL save_variable(self % save_path, "X_slicexz", xs)
    CALL save_variable(self % save_path, "Z_slicexz", zs)
    CALL save_variable(self % save_path, "Y_slicexz", y)
    CALL save_variable(self % save_path, "Nx_slicexz", nx)
    CALL save_variable(self % save_path, "Nz_slicexz", nz)

    write(*, '(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_sliceXZ


SUBROUTINE run_sliceVertical(self)
    ! -------------------------------------------------------------------------
    ! Computes pressure on a vertical plane whose horizontal normal points
    ! along azimuth_sv (0 deg: normal to +x, anticlockwise). The plane spans
    ! the horizontal axis span = (-sin, cos)(azimuth) and the depth. With the
    ! automatic width every source is projected onto the span axis and the
    ! plane covers all of them plus margin_sv on both sides, so the same
    ! routine is valid for a single turbine and for a wind farm. Observers
    ! are ordered with z varying fastest, as in the Python implementation.
    ! Saves the reshaped pressure/observer grids and the plane geometry
    ! metadata to the solver's HDF5 file.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: p(:,:)
    COMPLEX(WP) , ALLOCATABLE :: p_reshaped(:,:,:)
    REAL(WP)    , ALLOCATABLE :: observers(:,:)
    REAL(WP)    , ALLOCATABLE :: obs_reshaped(:,:,:)
    REAL(WP)    , ALLOCATABLE :: us(:), zs(:), proj(:)
    REAL(WP)                  :: azimuth, az_rad, normal(2), span(2), center(2), ulim(2), zlim(2)
    REAL(WP)                  :: du, dz, Depth
    INTEGER(I32)              :: nu, nz, i, j, idx, bs


    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_sliceVertical)

        ! --- 1. Plane orientation: normal and in-plane horizontal (span) axis --- !
        azimuth = pp % azimuth_sv
        az_rad  = azimuth * PI / 180.0_WP
        normal  = [ cos(az_rad), sin(az_rad)]
        span    = [-sin(az_rad), cos(az_rad)]

        ! --- 2. Centre: turbine axis, or mean turbine axis for a wind farm --- !
        if (any(pp % center_sv == DEFAULT_R)) pp % center_sv = self % sources % geometry % origin

        ! --- 3. Depth extent, clipped to the water column --- !
        if (pp % zlim_sv(1) == DEFAULT_R) pp % zlim_sv(1) = -Depth
        if (pp % zlim_sv(2) == DEFAULT_R) pp % zlim_sv(2) = 0.0_WP
        pp % zlim_sv(1) = max(pp % zlim_sv(1), -Depth)
        pp % zlim_sv(2) = min(pp % zlim_sv(2), 0.0_WP)

        ! --- 4. Span: every source projected onto the span axis (auto), or symmetric width --- !
        if (pp % width_sv == DEFAULT_R) then
            proj          = (self % sources % xyz(:, 1) - pp % center_sv(1)) * span(1) + &
                            (self % sources % xyz(:, 2) - pp % center_sv(2)) * span(2)
            pp % ulim_sv  = [minval(proj) - pp % margin_sv, maxval(proj) + pp % margin_sv]
            pp % width_sv = pp % ulim_sv(2) - pp % ulim_sv(1)
        else
            pp % ulim_sv = [-0.5_WP * pp % width_sv, 0.5_WP * pp % width_sv]
        end if

        ! --- 5. Grid sizing: ~5 m along the span, ~2 m along z, at least 30 points each --- !
        if (pp % nu_sv == DEFAULT_I) &
            pp % nu_sv = max(30_I32, ceiling((pp % ulim_sv(2) - pp % ulim_sv(1)) / 5.0_WP, I32))
        if (pp % nz_sv == DEFAULT_I) &
            pp % nz_sv = max(30_I32, ceiling((pp % zlim_sv(2) - pp % zlim_sv(1)) / 2.0_WP, I32))

        center = pp % center_sv
        ulim   = pp % ulim_sv
        zlim   = pp % zlim_sv
        nu     = pp % nu_sv
        nz     = pp % nz_sv
    end associate

    ! --- 6. Observer grid: span outer, z inner (z varies fastest) --- !
    ALLOCATE(us(nu), zs(nz), observers(nu*nz, 3))
    du = (ulim(2) - ulim(1)) / max(real(nu - 1_I32, WP), 1.0_WP)
    dz = (zlim(2) - zlim(1)) / max(real(nz - 1_I32, WP), 1.0_WP)

    do i = 1, nu
        us(i) = ulim(1) + real(i - 1_I32, WP) * du
    end do
    do j = 1, nz
        zs(j) = zlim(1) + real(j - 1_I32, WP) * dz
    end do

    idx = 1_I32
    do i = 1, nu
        do j = 1, nz
            observers(idx, 1) = center(1) + us(i) * span(1)
            observers(idx, 2) = center(2) + us(i) * span(2)
            observers(idx, 3) = zs(j)
            idx               = idx + 1_I32
        end do
    end do

    write(*, '(A,F7.2,A,I0,A,I0,A,I0,A,F9.2,A,F9.2,A)') &
        '>> [SliceVert]  azimuth=', azimuth, ' deg | grid=', nu, ' x ', nz, ' = ', nu * nz, &
        ' | center=(', center(1), ',', center(2), ') m'
    write(*, '(A,F6.2,A,F6.2,A,F6.2,A,F6.2,A)') &
        '   normal=(', normal(1), ', ', normal(2), ') | span=(', span(1), ', ', span(2), ')'

    ! --- 7. Compute pressure (observers closer than MASK_DISTANCE to a node get NaN) --- !
    CALL self % compute_pressure_masked('SliceVert', observers, bs, p)
    if (self % debug) write(*, '(A,F10.4)') '   Pressure norm : ', pressure_norm(p)

    ! --- 8. Reshape to (Nfreqs, nu, nz) and save --- !
    ALLOCATE(p_reshaped(size(p, 1), nu, nz))
    ALLOCATE(obs_reshaped(nu, nz, 3))

    idx = 1_I32
    do i = 1, nu
        do j = 1, nz
            obs_reshaped(i, j, :) = observers(idx, :)
            p_reshaped(:, i, j)   = p(:, idx)
            idx                   = idx + 1_I32
        end do
    end do

    CALL save_variable(self % save_path, "Freqs"          , self % sources % Freqs)
    CALL save_variable(self % save_path, "P_sliceVertical", p_reshaped)
    CALL save_variable(self % save_path, "Coords_sliceV"  , obs_reshaped)
    CALL save_variable(self % save_path, "U_sliceV"       , us)
    CALL save_variable(self % save_path, "Z_sliceV"       , zs)
    CALL save_variable(self % save_path, "Azimuth_sliceV" , azimuth)
    CALL save_variable(self % save_path, "Center_sliceV"  , center)

    write(*, '(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_sliceVertical


SUBROUTINE dummy_run_sphere(self)
    ! -------------------------------------------------------------------------
    ! Default (no-op) implementation of run_sphere, used by solvers that do
    ! not support a spherical free-field observer geometry. Only consumes its
    ! arguments and prints a "skipped" notice.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance


    write(*, '(A)') '>> [Sphere]     skipped, not available for: ' // trim(self % solver_name)

END SUBROUTINE dummy_run_sphere


SUBROUTINE run(self)
    ! -------------------------------------------------------------------------
    ! Config-driven entry point for the post-processing pipeline: runs, in a
    ! fixed order, every step whose do_<name> flag is enabled. The parser has
    ! already resolved do_all and switched off the steps that make no sense
    ! for the emitter (sphere and decay for a wind farm), so the same routine
    ! serves a single turbine, a wind farm or hand-made sources.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    INTEGER(I32) :: k, n_steps


    k = 0_I32

    associate(d => self % run_inputs % do_run)
        n_steps = int(count([d % do_spectrums, d % do_polar, d % do_decay, d % do_line, d % do_sphere, &
                             d % do_sliceXY, d % do_sliceXZ, d % do_sliceVertical, d % do_cylinder]), I32)

        write(*, '(A)') ''
        write(*, '(A)') repeat('=', 70)
        write(*, '(A,I0,A)') ' RUNNING ACOUSTIC POST-PROCESSING (', n_steps, ' steps)'
        write(*, '(A)') repeat('=', 70)

        if (d % do_spectrums) then
            CALL print_step_header(k, n_steps, 'Spectrums at observer points')
            CALL self % run_spectrums()
        end if

        if (d % do_polar) then
            CALL print_step_header(k, n_steps, 'Polar contour')
            CALL self % run_polar()
        end if

        if (d % do_decay) then
            CALL print_step_header(k, n_steps, 'Distance decay along wind')
            CALL self % run_decay()
        end if

        if (d % do_line) then
            CALL print_step_header(k, n_steps, 'Line between two points')
            CALL self % run_line()
        end if

        if (d % do_sphere) then
            CALL print_step_header(k, n_steps, 'Spherical radiation (free-field)')
            CALL self % run_sphere()
        end if

        if (d % do_sliceXY) then
            CALL print_step_header(k, n_steps, 'Horizontal slice XY')
            CALL self % run_sliceXY()
        end if

        if (d % do_sliceXZ) then
            CALL print_step_header(k, n_steps, 'Vertical slice XZ')
            CALL self % run_sliceXZ()
        end if

        if (d % do_sliceVertical) then
            CALL print_step_header(k, n_steps, 'Vertical slice at any azimuth')
            CALL self % run_sliceVertical()
        end if

        if (d % do_cylinder) then
            CALL print_step_header(k, n_steps, 'Cylinder surface')
            CALL self % run_cylinder()
        end if
    end associate

    write(*, '(A)') ''
    write(*, '(A)') repeat('=', 70)
    write(*, '(A)') ' ALL ACOUSTIC COMPUTATIONS FINISHED'
    write(*, '(A)') repeat('=', 70)

END SUBROUTINE run


SUBROUTINE print_step_header(k, n_steps, title)
    ! -------------------------------------------------------------------------
    ! Advances the post-processing step counter and prints its "[k/N] title"
    ! header line.
    ! -------------------------------------------------------------------------

    INTEGER(I32)    , INTENT(INOUT) :: k          ! [-] Step counter, incremented on exit
    INTEGER(I32)    , INTENT(IN)    :: n_steps    ! [-] Number of enabled steps
    CHARACTER(len=*), INTENT(IN)    :: title      ! [-] Step title


    k = k + 1_I32

    write(*, '(A)') ''
    write(*, '(A,I0,A,I0,A)') '[', k, '/', n_steps, '] ' // title

END SUBROUTINE print_step_header

END MODULE AcousticSolver
