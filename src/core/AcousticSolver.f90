! ===============================================================================
! MODULE: AcousticSolver
! ===============================================================================
! Defines AcousticSolver_t, the abstract parent type for all acoustic
! propagation solvers (currently only ImageMethod_t). It holds the shared
! solver state (sources, environment, generic turbine configuration, save
! path, block size, debug/cluster flags) and implements the generic
! observer-geometry post-processing routines (spectrums, polar, cylinder,
! decay, line, sliceXY, sliceXZ, sphere) shared by every concrete solver.
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
USE WindTurbine, ONLY: WindTurbine_t

IMPLICIT NONE

PRIVATE
PUBLIC :: AcousticSolver_t, save_variable


! ===============================================================================
! ABSTRACT: AcousticSolver_t
! ===============================================================================
! Parent class for all acoustic propagation solvers (currently only
! ImageMethod_t). Holds the shared sources/environment/turbine state and
! implements the generic observer-geometry post-processing routines; init,
! compute_pressure and the parameter/path/summary routines are deferred.
! ===============================================================================
TYPE, ABSTRACT :: AcousticSolver_t

    TYPE(AcousticSource_t)     :: sources                       ! [-] Acoustic point-dipole sources (turbine nodes)
    TYPE(Config_Environment_t) :: env                           ! [-] Environment configuration (depth, sound speed, ...)
    TYPE(Config_WindTurbine_t) :: turbine_data                  ! [-] Generic turbine configuration (axis/baricenter, wind, depth, ...)
    TYPE(Config_Computation_t) :: run_inputs                    ! [-] Per-run_* post-processing step configuration
    INTEGER(I32)               :: block_size                    ! [-] Number of observers processed simultaneously
    LOGICAL                    :: debug                         ! [-] Whether to print extra debug info (e.g. pressure norms)
    LOGICAL                    :: in_cluster                    ! [-] Whether the solver runs on a cluster (affects I/O behavior)
    CHARACTER(len=50)          :: solver_name = 'None'          ! [-] Solver name
    CHARACTER(len=512)         :: save_path   = 'None'          ! [-] Path to save acoustics


    CONTAINS

    ! --- Public type-bound procedures --- !
    PROCEDURE :: check_observers_distances      ! Helper for run subroutines
    PROCEDURE :: update_parameters              ! Appends run_* parameters actually used to the .dat
    PROCEDURE :: run                            ! Config-driven entry point (do_all / individual do_*)
    PROCEDURE :: run_spectrums                  ! Subroutines to compute acoustic pressure
    PROCEDURE :: run_polar
    PROCEDURE :: run_cylinder
    PROCEDURE :: run_decay
    PROCEDURE :: run_line
    PROCEDURE :: run_sliceXY
    PROCEDURE :: run_sliceXZ
    PROCEDURE :: run_all
    PROCEDURE :: run_sphere => dummy_run_sphere

    ! --- Deferred(abstract) procedures --- !
    PROCEDURE(init_i)            , DEFERRED :: init
    PROCEDURE(get_name_i)        , DEFERRED :: get_name
    PROCEDURE(compute_pressure_i), DEFERRED :: compute_pressure
    PROCEDURE(save_parameters_i)   , DEFERRED :: save_parameters
    PROCEDURE(build_save_path_i)   , DEFERRED :: build_save_path
    PROCEDURE(write_txt_summary_i) , DEFERRED :: write_txt_summary

END TYPE AcousticSolver_t


! =========================================================================
! Abstract Interfaces for Deferred Procedures
! =========================================================================
! Initialization
ABSTRACT INTERFACE
    SUBROUTINE init_i(self, sources, config)
        IMPORT :: AcousticSolver_t, AcousticSource_t, Config_t
        CLASS(AcousticSolver_t), INTENT(INOUT)              :: self
        TYPE(AcousticSource_t) , INTENT(INOUT), ALLOCATABLE :: sources
        TYPE(Config_t)         , INTENT(INOUT), ALLOCATABLE :: config
    END SUBROUTINE init_i
END INTERFACE

! Get solver name
ABSTRACT INTERFACE
    FUNCTION get_name_i(self) RESULT(name)
        IMPORT :: AcousticSolver_t
        CLASS(AcousticSolver_t), INTENT(INOUT) :: self
        CHARACTER(len=:), ALLOCATABLE          :: name
    END FUNCTION get_name_i
END INTERFACE

! Compute pressure
ABSTRACT INTERFACE
    SUBROUTINE compute_pressure_i(self, observers, block_size, total_pressure)
        IMPORT AcousticSolver_t, AcousticSource_t, WP, I32
        CLASS(AcousticSolver_t) , INTENT(INOUT)        :: self
        REAL(WP)                , INTENT(IN)           :: observers(:,:)        ! [m] Observers coords (Nobs,3)
        INTEGER(I32)            , INTENT(IN), OPTIONAL :: block_size            ! [Pa] Pressure field (Nfreqs, Nobs)
        COMPLEX(WP), ALLOCATABLE, INTENT(OUT)          :: total_pressure(:,:)   ! [-] Number of observers to process simultaneously
    END SUBROUTINE compute_pressure_i
END INTERFACE

! Save parameters
ABSTRACT INTERFACE
    SUBROUTINE save_parameters_i(self, turbine)
        IMPORT:: AcousticSolver_t, WindTurbine_t
        CLASS(AcousticSolver_t), INTENT(IN) :: self
        CLASS(WindTurbine_t)   , INTENT(IN), OPTIONAL, ALLOCATABLE :: turbine
    END SUBROUTINE save_parameters_i
END INTERFACE

! Build name for the output file
ABSTRACT INTERFACE
    SUBROUTINE build_save_path_i(self)
        IMPORT:: AcousticSolver_t
        CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    END SUBROUTINE build_save_path_i
END INTERFACE

! Write the elegant .dat summary
ABSTRACT INTERFACE
    SUBROUTINE write_txt_summary_i(self, turbine)
        IMPORT:: AcousticSolver_t, WindTurbine_t
        CLASS(AcousticSolver_t), INTENT(IN)                        :: self
        CLASS(WindTurbine_t)   , INTENT(IN), OPTIONAL, ALLOCATABLE :: turbine
    END SUBROUTINE write_txt_summary_i
END INTERFACE


CONTAINS


! ---------- HELPERS ---------- !
! -------------------------------------------------------------------------
! Subroutine: check_observers_distances
! Safety check that aborts execution if any observer point coincides with
! (or sits closer than min_distance to) a turbine source node, which would
! make the pressure field computation singular.
! -------------------------------------------------------------------------
SUBROUTINE check_observers_distances(self, observers, min_distance)
    CLASS(AcousticSolver_t), INTENT(IN) :: self
    REAL(WP)               , INTENT(IN) :: observers(:,:)          ! (Nobs, 3)
    REAL(WP)               , INTENT(IN), OPTIONAL :: min_distance  ! [m]

    ! Local variables
    REAL(WP) :: min_dist
    INTEGER(I32) :: i, j, n_nodes, n_obs
    REAL(WP) :: dx, dy, dz, dist

    ! Default min value
    min_dist = 1.0_WP; if (PRESENT(min_distance)) min_dist = min_distance
    
    n_nodes = size(self % sources % xyz, 1)
    n_obs   = size(observers, 1)
    
    do i = 1, n_obs
        do j = 1, n_nodes
            dx = self % sources % xyz(j,1) - observers(i,1)
            dy = self % sources % xyz(j,2) - observers(i,2)
            dz = self % sources % xyz(j,3) - observers(i,3)
            dist = sqrt(dx*dx + dy*dy + dz*dz)

            if (abs(dist - min_dist)<1e-6_WP) then
                write(*, '(A,I0,A,3(F12.4,1X),A,F6.2,A)') &
                "Observer ", i, " at (", observers(i,:), ") is too close to a turbine node (min distance: ", min_dist, " m)."
                error stop
            end if
        end do
    end do

END SUBROUTINE check_observers_distances


! -------------------------------------------------------------------------
! Subroutine: update_parameters
! Appends to the .dat summary already written by save_parameters() one
! block per post-processing step that actually ran (per do_run), with the
! resolved parameters run()/run_all() used for it, then closes the file
! with its final separator. Must be called after run()/run_all(), since
! every step only resolves its own DEFAULT_R/DEFAULT_I sentinels when it
! actually executes.
! -------------------------------------------------------------------------
SUBROUTINE update_parameters(self)
    CLASS(AcousticSolver_t), INTENT(IN) :: self

    CHARACTER(len=:), ALLOCATABLE :: dat_path
    INTEGER(I32) :: u, i

    dat_path = swap_extension(trim(self % save_path), ".dat")

    open(NEWUNIT=u, FILE=dat_path, STATUS='old', POSITION='append', ACTION='write')

    write(u,'(A)') repeat('-', 70)
    write(u,'(A)') 'POST-PROCESSING STEPS EXECUTED'
    write(u,'(A)') repeat('-', 70)

    associate(d => self % run_inputs % do_run, pp => self % run_inputs)

                if (d % do_spectrums) then
            write(u,'(A)') ''
            write(u,'(A)') 'Spectrums:'
            if (.not. ALLOCATED(pp % params_spectrums % observers_x)) then
                write(u,'(A)') '   observers : 3 default points (auto, x=10/250/500 m, z=-Depth/2)'
            else
                write(u,'(A,I0)') '   n_observers : ', pp % params_spectrums % n_observers
                do i = 1, pp % params_spectrums % n_observers
                    write(u,'(A,I0,A,F8.2,A,F8.2,A,F8.2,A)') '   obs ', i, ' (x,y,z) : ', &
                        pp % params_spectrums % observers_x(i), ', ', &
                        pp % params_spectrums % observers_y(i), ', ', &
                        pp % params_spectrums % observers_z(i), ' m'
                end do
            end if
        end if

        if (d % do_polar) then
            write(u,'(A)') ''
            write(u,'(A)')         'Polar:'
            write(u,'(A,F8.2,A)')  '   r_polar      : ', pp % params_polar % r_polar, ' m'
            write(u,'(A,F8.2,A)')  '   z_polar      : ', pp % params_polar % z_polar, ' m'
            write(u,'(A,I0)')      '   n_polar      : ', pp % params_polar % n_polar
            write(u,'(A,2F10.3)')  '   center_polar : ', pp % params_polar % center_polar
        end if

        if (d % do_cylinder) then
            write(u,'(A)') ''
            write(u,'(A)')        'Cylinder:'
            write(u,'(A,F8.2,A)') '   r_cyl      : ', pp % params_cylinder % r_cyl, ' m'
            write(u,'(A,I0)')     '   n_cyl      : ', pp % params_cylinder % n_cyl
            write(u,'(A,I0)')     '   nz_cyl     : ', pp % params_cylinder % nz_cyl
            write(u,'(A,2F10.3)') '   center_cyl : ', pp % params_cylinder % center_cyl
        end if

        if (d % do_decay) then
            write(u,'(A)') ''
            write(u,'(A)')        'Decay:'
            write(u,'(A,2F10.3)') '   distance_decay : ', pp % params_decay % distance_decay
            write(u,'(A,I0)')     '   n_decay        : ', pp % params_decay % n_decay
            write(u,'(A,F8.2,A)') '   z_decay        : ', pp % params_decay % z_decay, ' m'
            write(u,'(A,L1)')     '   logspace_decay : ', pp % params_decay % logspace_decay
        end if

        if (d % do_line) then
            write(u,'(A)') ''
            write(u,'(A)')        'Line:'
            write(u,'(A,3F10.3)') '   p1_line       : ', pp % params_line % p1_line
            write(u,'(A,3F10.3)') '   p2_line       : ', pp % params_line % p2_line
            write(u,'(A,I0)')     '   n_line        : ', pp % params_line % n_line
            write(u,'(A,L1)')     '   logspace_line : ', pp % params_line % logspace_line
        end if

        if (d % do_sliceXY) then
            write(u,'(A)') ''
            write(u,'(A)')        'SliceXY:'
            write(u,'(A,F8.2,A)') '   z_xy       : ', pp % params_sliceXY % z_xy, ' m'
            write(u,'(A,I0)')     '   nx_xy      : ', pp % params_sliceXY % nx_xy
            write(u,'(A,I0)')     '   ny_xy      : ', pp % params_sliceXY % ny_xy
            write(u,'(A,2F10.3)') '   xlim_xy    : ', pp % params_sliceXY % xlim_xy
            write(u,'(A,2F10.3)') '   ylim_xy    : ', pp % params_sliceXY % ylim_xy
            write(u,'(A,2F10.3)') '   center_xy  : ', pp % params_sliceXY % center_xy
        end if

        if (d % do_sliceXZ) then
            write(u,'(A)') ''
            write(u,'(A)')        'SliceXZ:'
            write(u,'(A,F8.2,A)') '   y_xz       : ', pp % params_sliceXZ % y_xz, ' m'
            write(u,'(A,I0)')     '   nx_xz      : ', pp % params_sliceXZ % nx_xz
            write(u,'(A,I0)')     '   nz_xz      : ', pp % params_sliceXZ % nz_xz
            write(u,'(A,2F10.3)') '   xlim_xz    : ', pp % params_sliceXZ % xlim_xz
            write(u,'(A,2F10.3)') '   zlim_xz    : ', pp % params_sliceXZ % zlim_xz
        end if

        if (d % do_sphere) then
            write(u,'(A)') ''
            write(u,'(A)')        'Sphere:'
            write(u,'(A,F8.2,A)') '   r_sphere       : ', pp % params_sphere % r_sphere, ' m'
            write(u,'(A,I0)')     '   n_theta_sphere : ', pp % params_sphere % n_theta_sphere
            write(u,'(A,I0)')     '   nz_sphere      : ', pp % params_sphere % nz_sphere
            write(u,'(A,3F10.3)') '   center_sphere  : ', pp % params_sphere % center_sphere
        end if

    end associate

    write(u,*)
    write(u,'(A)') repeat('=', 70)

    close(u)

    print*, '   ---> Parameters updated at ', dat_path

END SUBROUTINE update_parameters


! ---------- RUN FUNCTIONS ---------- !
! -------------------------------------------------------------------------
! Subroutine: run_spectrums
! Computes the pressure spectrum at a small set of fixed observer points
! (default: 3 points along +x at mid-depth). Saves the pressure, observer
! coordinates and frequency vector to the solver's HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE run_spectrums(self)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    INTEGER(I32) :: bs

    write(*,'(A)') '>> [Spectrums] computing pressure at observer points...'

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_spectrums)
        if (pp % n_observers == DEFAULT_I) then
            ALLOCATE(observers_(3,3))
            observers_(:,1) = [10.0_WP, 250.0_WP, 500.0_WP]
            observers_(:,2) = 0.0_WP
            observers_(:,3) = - self % turbine_data % Depth / 2.0_WP
        else
            ALLOCATE(observers_(pp % n_observers, 3))
            observers_(:,1) = pp % observers_x
            observers_(:,2) = pp % observers_y
            observers_(:,3) = pp % observers_z
        end if
    end associate

    bs = self % block_size

    if (any(observers_(:,3) > 0.0_WP)) & 
    write(*,*) "WARNING AcousticSolver % run_spectrums: there are observers placed over z = 0 m"


    ! Compute pressure
    CALL self % check_observers_distances(observers_)
    CALL self % compute_pressure(observers_, bs, p)

    if (self % debug) write(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    ! Save data
    CALL save_variable(self % save_path, "P_spectrums", p)
    CALL save_variable(self % save_path, "Obs_spectrums", observers_)
    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)

    write(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_spectrums


! -------------------------------------------------------------------------
! Subroutine: run_polar
! Computes pressure on a horizontal circle of radius r_polar around the
! turbine barycenter, at fixed depth z_polar. Saves the pressure plus the
! polar geometry (radius, depth, azimuth angles, observer coordinates) to
! the solver's HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE run_polar(self)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers(:,:), theta_deg_polar(:)
    REAL(WP)                 :: r, z, center(2), theta_deg, theta_rad, d_theta, Depth
    INTEGER(I32)             :: n_theta, bs, i

    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_polar)
        if (pp % r_polar == DEFAULT_R) pp % r_polar = 500.0_WP
        if (pp % z_polar == DEFAULT_R) pp % z_polar = - Depth/2.0_WP
        if (pp % n_polar == DEFAULT_I) pp % n_polar = 72_I32
        if (all(pp % center_polar == DEFAULT_R)) &   ! not a logical .and.
            pp % center_polar = [self % turbine_data % BariPos(1), self % turbine_data % BariPos(2)]

        r       = pp % r_polar
        z       = pp % z_polar
        n_theta = pp % n_polar
        center  = pp % center_polar
    end associate

    ! Build observers array
    ALLOCATE(observers(n_theta, 3), theta_deg_polar(n_theta))
    d_theta = 360.0_WP / REAL(n_theta, WP)

    do i = 1, n_theta
        theta_deg = REAL(i - 1, WP) * d_theta
        theta_rad = theta_deg * PI / 180.0_WP
        observers(i, 1)  = r * cos(theta_rad) + center(1)
        observers(i, 2)  = r * sin(theta_rad) + center(2)
        observers(i, 3)  = z
        theta_deg_polar(i) = theta_deg
    end do

    CALL self % check_observers_distances(observers)

    write(*,'(A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,I0)') &
        '>> [Polar]    r=', r, ' m | center=(', center(1), ',', center(2), &
        ') m | z=', z, ' m | observers=', n_theta

    ! Compute_pressure
    CALL self % compute_pressure(observers, bs, p)
    if (self % debug) write(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    ! Save data
    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_polar", p)
    CALL save_variable(self % save_path, "R_polar", r)
    CALL save_variable(self % save_path, "Z_polar", z)
    CALL save_variable(self % save_path, "Theta_deg_polar", theta_deg_polar)
    CALL save_variable(self % save_path, "Obs_polar", observers)
    CALL save_variable(self % save_path, "Center_polar", center)

    write(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_polar


! -------------------------------------------------------------------------
! Subroutine: run_cylinder
! Computes pressure on a vertical cylindrical surface (radius r_cyl,
! n_cyl azimuth points x nz_cyl depth levels) around the turbine
! barycenter. Saves the reshaped pressure/observer grids and the cylinder
! geometry metadata to the solver's HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE run_cylinder(self)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers(:,:)
    REAL(WP)                 :: r, center(2), theta_rad, d_theta, z_val, dz, Depth
    INTEGER(I32)             :: n_theta, nz, bs, i, j, idx, Nobs

    REAL(WP), ALLOCATABLE    :: theta_deg_arr(:)
    REAL(WP), ALLOCATABLE    :: z_cylinder(:)
    REAL(WP), ALLOCATABLE    :: obs_reshaped(:,:,:)
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
            pp % center_cyl = [self % turbine_data % BariPos(1), self % turbine_data % BariPos(2)]

        r       = pp % r_cyl
        n_theta = pp % n_cyl
        nz      = pp % nz_cyl
        center  = pp % center_cyl
    end associate

    ! --- resto del cuerpo sin cambios ---
    Nobs = n_theta * nz
    ALLOCATE(observers(Nobs, 3))

    d_theta = 360.0_WP / REAL(n_theta, WP)
    dz = (0.0_WP - (-Depth)) / max(REAL(nz - 1, WP), 1.0_WP)

    ALLOCATE(theta_deg_arr(n_theta))
    ALLOCATE(z_cylinder(nz))

    do j = 1, nz
        z_cylinder(j) = - Depth + REAL(j - 1, WP) * dz
    end do

    idx = 1
    do i = 1, n_theta
        theta_deg_arr(i) = REAL(i - 1, WP) * d_theta
        theta_rad = theta_deg_arr(i) * PI / 180.0_WP
        do j = 1, nz
            z_val = z_cylinder(j)
            observers(idx, 1) = r * cos(theta_rad) + center(1)
            observers(idx, 2) = r * sin(theta_rad) + center(2)
            observers(idx, 3) = z_val
            idx = idx + 1
        end do
    end do

    CALL self % check_observers_distances(observers)

    write(*,'(A,F8.2,A,F8.2,A,F8.2,A,I0,A,I0,A,I0)') &
        '>> [Cylinder]   r=', r, ' m | center=(', center(1), ',', center(2), &
        ') m | grid=', n_theta, ' x ', nz, ' = ', Nobs

    CALL self % compute_pressure(observers, bs, p)
    if (self % debug) write(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    ALLOCATE(p_reshaped(SIZE(p,1), n_theta, nz))
    ALLOCATE(obs_reshaped(n_theta, nz, 3))

    idx = 1
    do i = 1, n_theta
        do j = 1, nz
            obs_reshaped(i,j,:) = observers(idx,:)
            p_reshaped(:,i,j)   = p(:,idx)
            idx = idx + 1
        end do
    end do

    dA = r * dz * (2.0_WP * PI / REAL(n_theta, WP))

    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_cylinder", p_reshaped)
    CALL save_variable(self % save_path, "R_cylinder", r)
    CALL save_variable(self % save_path, "Z_cylinder", z_cylinder)
    CALL save_variable(self % save_path, "Theta_deg_cylinder", theta_deg_arr)
    CALL save_variable(self % save_path, "Obs_cylinder", obs_reshaped)
    CALL save_variable(self % save_path, "Center_cylinder", center)
    CALL save_variable(self % save_path, "dA_cylinder", dA)

    write(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_cylinder


! -------------------------------------------------------------------------
! Subroutine: run_decay
! Computes pressure along a straight line following the wind direction
! from the turbine axis, at fixed depth z_decay, with linear or log-spaced
! distances. Saves the pressure and the observer/distance geometry to the
! solver's HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE run_decay(self)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers(:,:)
    REAL(WP)   , ALLOCATABLE :: distances(:)
    REAL(WP)                 :: distance(2), z, wind_unit(2), d_val, WindDir_rad
    REAL(WP)                 :: d_log1, d_log2, d_step, linear_step, Depth
    INTEGER(I32)             :: n_points, bs, i
    LOGICAL                  :: logspace
    INTEGER                  :: logspace_int

    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_decay)
        if (pp % z_decay == DEFAULT_R)             pp % z_decay        = -Depth/2.0_WP
        if (any(pp % distance_decay == DEFAULT_R)) pp % distance_decay = [10.0_WP, 500.0_WP]
        if (pp % n_decay == DEFAULT_I)             pp % n_decay        = 200_I32
        ! logspace_decay no usa sentinel: es un LOGICAL directo del .nml

        z        = pp % z_decay
        distance = pp % distance_decay
        n_points = pp % n_decay
        logspace = pp % logspace_decay
    end associate

    ! --- resto del cuerpo sin cambios ---
    WindDir_rad = self % turbine_data % WindDir * PI / 180.0_WP
    wind_unit = [cos(WindDir_rad), sin(WindDir_rad)]

    ALLOCATE(observers(n_points, 3))
    ALLOCATE(distances(n_points))

    if (logspace) then
        d_log1 = log10(distance(1))
        d_log2 = log10(distance(2))
        d_step = (d_log2 - d_log1) / REAL(n_points - 1, WP)
        do i = 1, n_points
            d_val = 10.0_WP ** (d_log1 + REAL(i - 1, WP) * d_step)
            distances(i)      = d_val
            observers(i, 1:2) = self % turbine_data % AxisPos(1:2) + d_val * wind_unit
            observers(i, 3)   = z
        end do
    else
        linear_step = (distance(2) - distance(1)) / REAL(n_points - 1, WP)
        do i = 1, n_points
            d_val = distance(1) + REAL(i - 1, WP) * linear_step
            distances(i)      = d_val
            observers(i, 1:2) = self % turbine_data % AxisPos(1:2) + d_val * wind_unit
            observers(i, 3)   = z
        end do
    end if

    CALL self % check_observers_distances(observers)

    write(*,'(A,I0,A,F8.2,A,F8.2,A)') &
        '>> [Decay]      points=', n_points, ' | direction=', &
        self % turbine_data % WindDir, ' deg | z=', z, ' m'

    CALL self % compute_pressure(observers, bs, p)
    if (self % debug) write(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    logspace_int = merge(1, 0, logspace)

    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_decay", p)
    CALL save_variable(self % save_path, "Z_decay", z)
    CALL save_variable(self % save_path, "Distance_decay", distance)
    CALL save_variable(self % save_path, "Distances_decay", distances)
    CALL save_variable(self % save_path, "Obs_decay", observers)
    CALL save_variable(self % save_path, "Logspace_decay", logspace_int)

    write(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_decay


! -------------------------------------------------------------------------
! Subroutine: run_line
! Computes pressure along an arbitrary straight line between two points
! p1_line and p2_line, with linear or log-spaced distances. Saves the
! pressure and the observer/distance geometry to the solver's HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE run_line(self)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers(:,:)
    REAL(WP)   , ALLOCATABLE :: distances(:)
    REAL(WP)                 :: p1(3), p2(3), dir_vec(3), dist_total, d_val
    REAL(WP)                 :: d_log1, d_log2, d_step, linear_step, Depth
    INTEGER(I32)             :: n_points, i, bs
    LOGICAL                  :: logspace
    INTEGER                  :: logspace_int

    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_line)
        if (any(pp % p1_line == DEFAULT_R)) pp % p1_line = [100.0_WP, 0.0_WP, 0.0_WP]
        if (any(pp % p2_line == DEFAULT_R)) pp % p2_line = [100.0_WP, 0.0_WP, - Depth]
        ! logspace_line no usa sentinel: es un LOGICAL directo del .nml

        p1       = pp % p1_line
        p2       = pp % p2_line
        logspace = pp % logspace_line
    end associate

    dist_total = norm2(p2 - p1)
    if (dist_total == 0.0_WP) error stop "AcousticSolver.run_line(): points p1 and p2 are coincident."

    ! n_line depende de dist_total, se resuelve aparte
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
        d_log1 = log10(dist_total * 1e-6_WP)
        d_log2 = log10(dist_total)
        d_step = (d_log2 - d_log1) / REAL(n_points - 1, WP)
        do i = 1, n_points
            d_val = 10.0_WP ** (d_log1 + REAL(i - 1, WP) * d_step)
            distances(i)    = d_val
            observers(i, :) = p1 + (d_val / dist_total) * dir_vec
        end do
    else
        linear_step = dist_total / REAL(n_points - 1, WP)
        do i = 1, n_points
            d_val = REAL(i - 1, WP) * linear_step
            distances(i) = d_val
            observers(i, :) = p1 + (d_val / dist_total) * dir_vec
        end do
    end if

    CALL self % check_observers_distances(observers)

    write(*,'(A,I0,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A)') &
        '>> [Line]       points=', n_points, ' | p1=(', p1(1), ',', p1(2), ',', p1(3), &
        ') m | p2=(', p2(1), ',', p2(2), ',', p2(3), ') m'

    CALL self % compute_pressure(observers, bs, p)
    if (self % debug) write(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    logspace_int = merge(1, 0, logspace)

    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_line", p)
    CALL save_variable(self % save_path, "P1_line", p1)
    CALL save_variable(self % save_path, "P2_line", p2)
    CALL save_variable(self % save_path, "Obs_line", observers)
    CALL save_variable(self % save_path, "Distance_line", dist_total)
    CALL save_variable(self % save_path, "Distances_line", distances)
    CALL save_variable(self % save_path, "Logspace_line", logspace_int)

    write(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_line


! -------------------------------------------------------------------------
! Subroutine: run_sliceXY
! Computes pressure on a horizontal (x,y) grid at fixed depth z_xy, over
! the given x/y limits. Saves the reshaped pressure/observer grids and the
! slice geometry metadata to the solver's HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE run_sliceXY(self)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    COMPLEX(WP), ALLOCATABLE :: p_reshaped(:,:,:)
    REAL(WP)   , ALLOCATABLE :: observers(:,:)
    REAL(WP)   , ALLOCATABLE :: obs_reshaped(:,:,:)
    REAL(WP)   , ALLOCATABLE :: xs(:), ys(:)
    REAL(WP)                 :: z, xlim(2), ylim(2), center(2), Depth
    REAL(WP)                 :: mid_x, mid_y, dx, dy, x_val, y_val
    INTEGER(I32)             :: nx, ny, i, j, idx, bs
    REAL(WP), PARAMETER      :: tol = 1e-5_WP

    Depth = self % env % depth_base
    bs   = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_sliceXY)
        if (pp % z_xy == DEFAULT_R) pp % z_xy = - Depth/2.0_WP
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
        else if (ABS(pp % center_xy(1) - mid_x) > tol .OR. ABS(pp % center_xy(2) - mid_y) > tol) then
            write(*,*) "AcousticSolver.run_sliceXY(): provided center does not match midpoint. Switching to proper center."
            pp % center_xy = [mid_x, mid_y]
        end if
        center = pp % center_xy
    end associate

    ! --- resto del cuerpo sin cambios ---
    ALLOCATE(observers(nx * ny, 3))
    ALLOCATE(xs(nx), ys(ny))
    dx = (xlim(2) - xlim(1)) / max(REAL(nx - 1, WP), 1.0_WP)
    dy = (ylim(2) - ylim(1)) / max(REAL(ny - 1, WP), 1.0_WP)

    do j = 1, nx
        xs(j) = xlim(1) + REAL(j - 1, WP) * dx
    end do
    do i = 1, ny
        ys(i) = ylim(1) + REAL(i - 1, WP) * dy
    end do

    idx = 1
    do i = 1, ny
        y_val = ys(i)
        do j = 1, nx
            x_val = xs(j)
            observers(idx, 1) = x_val
            observers(idx, 2) = y_val
            observers(idx, 3) = z
            idx = idx + 1
        end do
    end do

    CALL self % check_observers_distances(observers, 0.1_WP)

    write(*,'(A,I0,A,I0,A,F8.2,A,F8.2,A,F8.2,A)') &
        '>> [SliceXY]    nx=', nx, ' | ny=', ny, ' | z=', z, &
        ' m | center=(', center(1), ',', center(2), ') m'

    CALL self % compute_pressure(observers, bs, p)
    if (self % debug) write(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    ALLOCATE(p_reshaped(SIZE(p,1), nx, ny))
    ALLOCATE(obs_reshaped(nx, ny, 3))

    idx = 1
    do i = 1, ny
        do j = 1, nx
            obs_reshaped(j, i, :) = observers(idx, :)
            p_reshaped(:, j, i)   = p(:, idx)
            idx = idx + 1
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

    write(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_sliceXY


! -------------------------------------------------------------------------
! Subroutine: run_sliceXZ
! Computes pressure on a vertical (x,z) grid at fixed y_xz, over the given
! x/z limits. Saves the reshaped pressure/observer grids and the slice
! geometry metadata to the solver's HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE run_sliceXZ(self)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    COMPLEX(WP), ALLOCATABLE :: p_reshaped(:,:,:)
    REAL(WP)   , ALLOCATABLE :: observers(:,:)
    REAL(WP)   , ALLOCATABLE :: obs_reshaped(:,:,:)
    REAL(WP)   , ALLOCATABLE :: xs(:), zs(:)
    REAL(WP)                 :: y, xlim(2), zlim(2), dx, dz, x_val, z_val, Depth
    INTEGER(I32)             :: nx, nz, i, j, idx, bs

    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_sliceXZ)
        if (pp % y_xz == DEFAULT_R)         pp % y_xz    = self % turbine_data % AxisPos(2)
        if (pp % nx_xz == DEFAULT_I)        pp % nx_xz   = 26_I32
        if (pp % nz_xz == DEFAULT_I)        pp % nz_xz   = 26_I32
        if (any(pp % xlim_xz == DEFAULT_R)) pp % xlim_xz = [-500.0_WP, 500.0_WP]
        if (any(pp % zlim_xz == DEFAULT_R)) pp % zlim_xz = [0.0_WP, - Depth]

        y    = pp % y_xz
        nx   = pp % nx_xz
        nz   = pp % nz_xz
        xlim = pp % xlim_xz
        zlim = pp % zlim_xz
    end associate

    ! --- resto del cuerpo sin cambios ---
    ALLOCATE(observers(nx * nz, 3))
    ALLOCATE(xs(nx), zs(nz))
    dx = (xlim(2) - xlim(1)) / max(REAL(nx - 1, WP), 1.0_WP)
    dz = (zlim(2) - zlim(1)) / max(REAL(nz - 1, WP), 1.0_WP)

    do j = 1, nx
        xs(j) = xlim(1) + REAL(j - 1, WP) * dx
    end do
    do i = 1, nz
        zs(i) = zlim(1) + REAL(i - 1, WP) * dz
    end do

    idx = 1
    do i = 1, nz
        z_val = zs(i)
        do j = 1, nx
            x_val = xs(j)
            observers(idx, 1) = x_val
            observers(idx, 2) = y
            observers(idx, 3) = z_val
            idx = idx + 1
        end do
    end do

    CALL self % check_observers_distances(observers, 0.1_WP)

    write(*,'(A,I0,A,I0,A,F8.2,A,F8.2,A,F8.2,A)') &
        '>> [SliceXZ]    nx=', nx, ' | nz=', nz, ' | y=', y, &
        ' m | z from ', zlim(1), ' to ', zlim(2), ' m'

    CALL self % compute_pressure(observers, bs, p)
    if (self % debug) write(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    ALLOCATE(p_reshaped(SIZE(p,1), nx, nz))
    ALLOCATE(obs_reshaped(nx, nz, 3))

    idx = 1
    do i = 1, nz
        do j = 1, nx
            obs_reshaped(j, i, :) = observers(idx, :)
            p_reshaped(:, j, i)   = p(:, idx)
            idx = idx + 1
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

    write(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_sliceXZ


! -------------------------------------------------------------------------
! Subroutine: dummy_run_sphere
! Default (no-op) implementation of run_sphere, used by solvers that do
! not support a spherical free-field observer geometry. Only consumes its
! arguments and prints a "skipped" notice.
! -------------------------------------------------------------------------
SUBROUTINE dummy_run_sphere(self)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    write(*,'(A)') '>> [Sphere]     skipped, not available for: ' // trim(self % solver_name)

END SUBROUTINE dummy_run_sphere


! -------------------------------------------------------------------------
! Subroutine: run_all
! Orchestrates the full acoustic post-processing pipeline: runs, in order,
! spectrums, polar, decay, line, sphere, sliceXY, sliceXZ and cylinder,
! forwarding each routine's optional arguments and saving every result to
! the solver's HDF5 output file.
! -------------------------------------------------------------------------
SUBROUTINE run_all(self)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    write(*,'(A)') ''
    write(*,'(A)') repeat('=', 70)
    write(*,'(A)') ' RUNNING ALL ACOUSTIC POST-PROCESSING (8 steps)'
    write(*,'(A)') repeat('=', 70)

    write(*,'(A)') ''
    write(*,'(A)') '[1/8] Spectrums at observer points'
    CALL self % run_spectrums()

    write(*,'(A)') ''
    write(*,'(A)') '[2/8] Polar contour'
    CALL self % run_polar()

    write(*,'(A)') ''
    write(*,'(A)') '[3/8] Distance decay along wind'
    CALL self % run_decay()

    write(*,'(A)') ''
    write(*,'(A)') '[4/8] Line between two points'
    CALL self % run_line()

    write(*,'(A)') ''
    write(*,'(A)') '[5/8] Spherical radiation (free-field)'
    CALL self % run_sphere()

    write(*,'(A)') ''
    write(*,'(A)') '[6/8] Horizontal slice XY'
    CALL self % run_sliceXY()

    write(*,'(A)') ''
    write(*,'(A)') '[7/8] Vertical slice XZ'
    CALL self % run_sliceXZ()

    write(*,'(A)') ''
    write(*,'(A)') '[8/8] Cylinder surface'
    CALL self % run_cylinder()

    write(*,'(A)') ''
    write(*,'(A)') repeat('=', 70)
    write(*,'(A)') ' ALL ACOUSTIC COMPUTATIONS FINISHED'
    write(*,'(A)') repeat('=', 70)

END SUBROUTINE run_all


! -------------------------------------------------------------------------
! Subroutine: run
! Config-driven entry point for the post-processing pipeline. Calls
! run_all() when do_run % do_all is enabled; otherwise runs only the
! individual steps whose own do_<name> flag is enabled.
! -------------------------------------------------------------------------
SUBROUTINE run(self)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    associate(d => self % run_inputs % do_run)
        if (d % do_all) then
            CALL self % run_all()
        else
            if (d % do_spectrums) CALL self % run_spectrums()
            if (d % do_polar)     CALL self % run_polar()
            if (d % do_cylinder)  CALL self % run_cylinder()
            if (d % do_decay)     CALL self % run_decay()
            if (d % do_line)      CALL self % run_line()
            if (d % do_sliceXY)   CALL self % run_sliceXY()
            if (d % do_sliceXZ)   CALL self % run_sliceXZ()
            if (d % do_sphere)    CALL self % run_sphere()
        end if
    end associate

END SUBROUTINE run


END MODULE AcousticSolver