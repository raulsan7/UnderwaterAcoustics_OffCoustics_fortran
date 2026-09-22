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
! Date    : 2026-09-22
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

    TYPE(AcousticSource_t)     :: sources                        ! [-] Acoustic point-dipole sources (turbine nodes)
    TYPE(Config_Environment_t) :: env                            ! [-] Environment configuration (depth, sound speed, ...)
    TYPE(Config_WindTurbine_t) :: turbine_data                   ! [-] Generic turbine configuration (axis/baricenter, wind, depth, ...)
    INTEGER(I32)               :: block_size                     ! [-] Number of observers processed simultaneously
    LOGICAL                    :: debug                          ! [-] Whether to print extra debug info (e.g. pressure norms)
    LOGICAL                    :: in_cluster                     ! [-] Whether the solver runs on a cluster (affects I/O behavior)
    CHARACTER(len=50)          :: solver_name = 'None'           ! [-] Solver name
    CHARACTER(len=512)         :: save_path   = 'None'           ! [-] Path to save acoustics


    CONTAINS

    ! --- Public type-bound procedures --- !
    PROCEDURE :: check_observers_distances      ! Helper for run subroutines
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
        CLASS(AcousticSolver_t)  , INTENT(INOUT)           :: self
        TYPE(AcousticSource_t), INTENT(INOUT), ALLOCATABLE :: sources
        TYPE(Config_t)        , INTENT(IN)                 :: config
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


! ---------- RUN FUNCTIONS ---------- !
! -------------------------------------------------------------------------
! Subroutine: run_spectrums
! Computes the pressure spectrum at a small set of fixed observer points
! (default: 3 points along +x at mid-depth). Saves the pressure, observer
! coordinates and frequency vector to the solver's HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE run_spectrums(self, observers, z_spectrums)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    REAL(WP), INTENT(IN), OPTIONAL :: observers(:,:)
    REAL(WP), INTENT(IN), OPTIONAL :: z_spectrums

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)                 :: z_spectrums_
    INTEGER(I32) :: bs

    WRITE(*,'(A)') '>> [Spectrums] computing pressure at observer points...'

    ! Defaults
    z_spectrums_ = - self % turbine_data % Depth/2.0_WP; if (PRESENT(z_spectrums)) z_spectrums_ = z_spectrums
    
    if (PRESENT(observers)) then
        observers_ = observers
    else
        ALLOCATE(observers_(3,3))
        observers_(:,3) = z_spectrums_
        observers_(:,2) = 0.0_WP
        observers_(:,1) = [10.0_WP, 250.0_WP, 500.0_WP]
    end if

    bs = self % block_size

    ! Compute pressure
    CALL self % check_observers_distances(observers_)
    CALL self % compute_pressure(observers_, bs, p)

    if (self % debug) WRITE(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    ! Save data
    CALL save_variable(self % save_path, "P_spectrums", p)
    CALL save_variable(self % save_path, "Obs_spectrums", observers_)
    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)

    WRITE(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_spectrums


! -------------------------------------------------------------------------
! Subroutine: run_polar
! Computes pressure on a horizontal circle of radius r_polar around the
! turbine barycenter, at fixed depth z_polar. Saves the pressure plus the
! polar geometry (radius, depth, azimuth angles, observer coordinates) to
! the solver's HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE run_polar(self, r_polar, z_polar, n_polar, center_polar)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    REAL(WP)    , INTENT(IN), OPTIONAL  :: r_polar
    REAL(WP)    , INTENT(IN), OPTIONAL  :: z_polar
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_polar
    REAL(WP)    , INTENT(IN), OPTIONAL  :: center_polar(2)

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers(:,:), theta_deg_polar(:)
    REAL(WP)                 :: r, z, center(2), theta_deg, theta_rad, d_theta, Depth
    INTEGER(I32)             :: n_theta, bs, i

    
    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults
    r = 500.0_WP      ; if (PRESENT(r_polar)) r       = r_polar
    z = - Depth/2.0_WP; if (PRESENT(z_polar)) z       = z_polar
    n_theta = 72_I32  ; if (PRESENT(n_polar)) n_theta = n_polar

    if (PRESENT(center_polar)) then
        center = center_polar
    else
        center = [self % turbine_data % BariPos(1), self % turbine_data % BariPos(2)]
    end if

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

    WRITE(*,'(A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,I0)') &
        '>> [Polar]    r=', r, ' m | center=(', center(1), ',', center(2), &
        ') m | z=', z, ' m | observers=', n_theta

    ! Compute_pressure
    CALL self % compute_pressure(observers, bs, p)
    if (self % debug) WRITE(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    ! Save data
    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_polar", p)
    CALL save_variable(self % save_path, "R_polar", r)
    CALL save_variable(self % save_path, "Z_polar", z)
    CALL save_variable(self % save_path, "Theta_deg_polar", theta_deg_polar)
    CALL save_variable(self % save_path, "Obs_polar", observers)
    CALL save_variable(self % save_path, "Center_polar", center)

    WRITE(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_polar


! -------------------------------------------------------------------------
! Subroutine: run_cylinder
! Computes pressure on a vertical cylindrical surface (radius r_cyl,
! n_cyl azimuth points x nz_cyl depth levels) around the turbine
! barycenter. Saves the reshaped pressure/observer grids and the cylinder
! geometry metadata to the solver's HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE run_cylinder(self, r_cyl, n_cyl, nz_cyl, center_cyl)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    REAL(WP)    , INTENT(IN), OPTIONAL :: r_cyl
    INTEGER(I32), INTENT(IN), OPTIONAL :: n_cyl
    INTEGER(I32), INTENT(IN), OPTIONAL :: nz_cyl
    REAL(WP)    , INTENT(IN), OPTIONAL :: center_cyl(2)

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers(:,:)
    REAL(WP)                 :: r, center(2), theta_rad, d_theta, z_val, dz, Depth
    INTEGER(I32)             :: n_theta, nz, bs, i, j, idx, Nobs

    ! Additional for saving
    REAL(WP), ALLOCATABLE    :: theta_deg_arr(:)      ! Azimuth angles in degrees
    REAL(WP), ALLOCATABLE    :: z_cylinder(:)         ! Vertical coordinates
    REAL(WP), ALLOCATABLE    :: obs_reshaped(:,:,:)   ! Reshaped observers (n_theta, nz, 3)
    COMPLEX(WP), ALLOCATABLE :: p_reshaped(:,:,:)     ! Reshaped pressure (Nfreqs, n_theta, nz)
    REAL(WP)                 :: dA                    ! Differential area element


    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults
    r = 500.0_WP    ; if (PRESENT(r_cyl))  r       = r_cyl
    n_theta = 72_I32; if (PRESENT(n_cyl))  n_theta = n_cyl
    nz = 20_I32     ; if (PRESENT(nz_cyl)) nz      = nz_cyl

    if (PRESENT(center_cyl)) then
        center = center_cyl
    else
        center = [self % turbine_data % BariPos(1), self % turbine_data % BariPos(2)]
    end if    
    
    Nobs = n_theta * nz
    ALLOCATE(observers(Nobs, 3))

    d_theta = 360.0_WP / REAL(n_theta, WP)
    dz = (0.0_WP - (-Depth)) / max(REAL(nz - 1, WP), 1.0_WP)

    ! Allocate arrays for storing angles and vertical coordinates
    ALLOCATE(theta_deg_arr(n_theta))
    ALLOCATE(z_cylinder(nz))

    ! Fill vertical coordinate array
    do j = 1, nz
        z_cylinder(j) = - Depth + REAL(j - 1, WP) * dz
    end do

    ! Build observers array and angle array (meshgrid equivalent)
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

    WRITE(*,'(A,F8.2,A,F8.2,A,F8.2,A,I0,A,I0,A,I0)') &
        '>> [Cylinder]   r=', r, ' m | center=(', center(1), ',', center(2), &
        ') m | grid=', n_theta, ' x ', nz, ' = ', Nobs

    CALL self % compute_pressure(observers, bs, p)
    if (self % debug) WRITE(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    ! Reshape pressure and observers to 3D for storage
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

    ! Compute differential area element (as in Python)
    dA = r * dz * (2.0_WP * PI / REAL(n_theta, WP))

    ! Save results
    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_cylinder", p_reshaped)
    CALL save_variable(self % save_path, "R_cylinder", r)
    CALL save_variable(self % save_path, "Z_cylinder", z_cylinder)
    CALL save_variable(self % save_path, "Theta_deg_cylinder", theta_deg_arr)
    CALL save_variable(self % save_path, "Obs_cylinder", obs_reshaped)
    CALL save_variable(self % save_path, "Center_cylinder", center)
    CALL save_variable(self % save_path, "dA_cylinder", dA)

    WRITE(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_cylinder


! -------------------------------------------------------------------------
! Subroutine: run_decay
! Computes pressure along a straight line following the wind direction
! from the turbine axis, at fixed depth z_decay, with linear or log-spaced
! distances. Saves the pressure and the observer/distance geometry to the
! solver's HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE run_decay(self, distance_decay, n_decay, z_decay, logspace_decay)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    REAL(WP)    , INTENT(IN), OPTIONAL  :: distance_decay(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_decay
    REAL(WP)    , INTENT(IN), OPTIONAL  :: z_decay
    LOGICAL     , INTENT(IN), OPTIONAL  :: logspace_decay

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers(:,:)
    REAL(WP)   , ALLOCATABLE :: distances(:)        ! Actual distances along the line
    REAL(WP)                 :: distance(2), z, wind_unit(2), d_val, WindDir_rad
    REAL(WP)                 :: d_log1, d_log2, d_step, linear_step, Depth
    INTEGER(I32)             :: n_points, bs, i
    LOGICAL                  :: logspace
    INTEGER                  :: logspace_int        ! Integer representation for saving


    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults
    z = -Depth/2.0_WP             ; if (PRESENT(z_decay))        z        = z_decay
    distance = [10.0_WP, 500.0_WP]; if (PRESENT(distance_decay)) distance = distance_decay
    n_points = 200_I32            ; if (PRESENT(n_decay))        n_points = n_decay
    logspace = .true.             ; if (PRESENT(logspace_decay)) logspace = logspace_decay
    
    ! Wind direction
    WindDir_rad = self % turbine_data % WindDir * PI / 180.0_WP
    wind_unit = [cos(WindDir_rad), sin(WindDir_rad)]

    ALLOCATE(observers(n_points, 3))
    ALLOCATE(distances(n_points))

    ! Distances along the line
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

    WRITE(*,'(A,I0,A,F8.2,A,F8.2,A)') &
        '>> [Decay]      points=', n_points, ' | direction=', &
        self % turbine_data % WindDir, ' deg | z=', z, ' m'

    CALL self % compute_pressure(observers, bs, p)
    if (self % debug) WRITE(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    ! Convert logical to integer for saving
    logspace_int = merge(1, 0, logspace)

    ! Save results
    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_decay", p)
    CALL save_variable(self % save_path, "Z_decay", z)
    CALL save_variable(self % save_path, "Distance_decay", distance)
    CALL save_variable(self % save_path, "Distances_decay", distances)
    CALL save_variable(self % save_path, "Obs_decay", observers)
    CALL save_variable(self % save_path, "Logspace_decay", logspace_int)

    WRITE(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_decay


! -------------------------------------------------------------------------
! Subroutine: run_line
! Computes pressure along an arbitrary straight line between two points
! p1_line and p2_line, with linear or log-spaced distances. Saves the
! pressure and the observer/distance geometry to the solver's HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE run_line(self, p1_line, p2_line, n_line, logspace_line)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    REAL(WP)    , INTENT(IN), OPTIONAL  :: p1_line(3)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: p2_line(3)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_line
    LOGICAL     , INTENT(IN), OPTIONAL  :: logspace_line

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers(:,:)
    REAL(WP)   , ALLOCATABLE :: distances(:)     ! Distances from p1 along the line
    REAL(WP)                 :: p1(3), p2(3), dir_vec(3), dist_total, d_val
    REAL(WP)                 :: d_log1, d_log2, d_step, linear_step, Depth
    INTEGER(I32)             :: n_points, i, bs
    LOGICAL                  :: logspace
    INTEGER                  :: logspace_int     ! Integer representation for saving


    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults
    p1 = [100.0_WP, 0.0_WP, 0.0_WP] ; if (PRESENT(p1_line))       p1       = p1_line
    p2 = [100.0_WP, 0.0_WP, - Depth]; if (PRESENT(p2_line))       p2       = p2_line
    logspace = .false.              ; if (PRESENT(logspace_line)) logspace = logspace_line

    dist_total = norm2(p2 - p1)
    if (dist_total == 0.0_WP) error stop "AcousticSolver.run_line(): points p1 and p2 are coincident."

    if (PRESENT(n_line)) then
        n_points = n_line
    else
        if (dist_total > 200.0_WP) then
            n_points = int(dist_total / 5.0_WP, I32)
        else
            n_points = int(dist_total, I32)
        end if
    end if

    ALLOCATE(observers(n_points, 3))
    ALLOCATE(distances(n_points))
    dir_vec = (p2 - p1) / dist_total

    ! Distances along the line
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

    WRITE(*,'(A,I0,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A)') &
        '>> [Line]       points=', n_points, ' | p1=(', p1(1), ',', p1(2), ',', p1(3), &
        ') m | p2=(', p2(1), ',', p2(2), ',', p2(3), ') m'

    CALL self % compute_pressure(observers, bs, p)
    if (self % debug) WRITE(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    ! Convert logical to integer for saving
    logspace_int = merge(1, 0, logspace)

    ! Save results
    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_line", p)
    CALL save_variable(self % save_path, "P1_line", p1)
    CALL save_variable(self % save_path, "P2_line", p2)
    CALL save_variable(self % save_path, "Obs_line", observers)
    CALL save_variable(self % save_path, "Distance_line", dist_total)
    CALL save_variable(self % save_path, "Distances_line", distances)
    CALL save_variable(self % save_path, "Logspace_line", logspace_int)

    WRITE(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_line


! -------------------------------------------------------------------------
! Subroutine: run_sliceXY
! Computes pressure on a horizontal (x,y) grid at fixed depth z_xy, over
! the given x/y limits. Saves the reshaped pressure/observer grids and the
! slice geometry metadata to the solver's HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE run_sliceXY(self, z_xy, nx_xy, ny_xy, xlim_xy, ylim_xy, center_xy)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    REAL(WP)    , INTENT(IN), OPTIONAL  :: z_xy
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nx_xy
    INTEGER(I32), INTENT(IN), OPTIONAL  :: ny_xy
    REAL(WP)    , INTENT(IN), OPTIONAL  :: xlim_xy(2)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: ylim_xy(2)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: center_xy(2)

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    COMPLEX(WP), ALLOCATABLE :: p_reshaped(:,:,:)   ! (Nfreqs, nx, ny)
    REAL(WP)   , ALLOCATABLE :: observers(:,:)
    REAL(WP)   , ALLOCATABLE :: obs_reshaped(:,:,:) ! (nx, ny, 3)
    REAL(WP)   , ALLOCATABLE :: xs(:), ys(:)        ! Coordinate arrays
    REAL(WP)                 :: z, xlim(2), ylim(2), center(2), Depth
    REAL(WP)                 :: mid_x, mid_y, dx, dy, x_val, y_val
    INTEGER(I32)             :: nx, ny, i, j, idx, bs
    REAL(WP), PARAMETER      :: tol = 1e-5_WP


    Depth = self % env % depth_base
    bs   = self % block_size

    ! Defaults
    z = - Depth/2.0_WP; if (PRESENT(z_xy)) z = z_xy
    if (z > 0.0_WP) then
        write(*,*) "AcousticSolver.run_sliceXY(): z must be <= 0. Switching to -z"
        z = -z
    end if

    nx   = 26_I32               ; if (PRESENT(nx_xy))   nx         = nx_xy
    ny   = 26_I32               ; if (PRESENT(ny_xy))   ny         = ny_xy
    xlim = [-500.0_WP, 500.0_WP]; if (PRESENT(xlim_xy)) xlim       = xlim_xy
    ylim = [-500.0_WP, 500.0_WP]; if (PRESENT(ylim_xy)) ylim       = ylim_xy

    mid_x = 0.5_WP * (xlim(1) + xlim(2))
    mid_y = 0.5_WP * (ylim(1) + ylim(2))

    if (PRESENT(center_xy)) then
        if (ABS(center_xy(1) - mid_x) > tol .OR. ABS(center_xy(2) - mid_y) > tol) then
            write(*,*) "AcousticSolver.run_sliceXY(): provided center does not match midpoint. Switching to proper center."
            center = [mid_x, mid_y]
        else
            center = center_xy
        end if
    else
        center = [mid_x, mid_y]
    end if

    ALLOCATE(observers(nx * ny, 3))
    ALLOCATE(xs(nx), ys(ny))
    dx = (xlim(2) - xlim(1)) / max(REAL(nx - 1, WP), 1.0_WP)
    dy = (ylim(2) - ylim(1)) / max(REAL(ny - 1, WP), 1.0_WP)

    ! Fill coordinate arrays
    do j = 1, nx
        xs(j) = xlim(1) + REAL(j - 1, WP) * dx
    end do
    do i = 1, ny
        ys(i) = ylim(1) + REAL(i - 1, WP) * dy
    end do

    ! Build observers (meshgrid equivalent)
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

    WRITE(*,'(A,I0,A,I0,A,F8.2,A,F8.2,A,F8.2,A)') &
        '>> [SliceXY]    nx=', nx, ' | ny=', ny, ' | z=', z, &
        ' m | center=(', center(1), ',', center(2), ') m'

    CALL self % compute_pressure(observers, bs, p)
    if (self % debug) WRITE(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    ! Reshape pressure and observers to 3D (nx, ny, 3) and (Nfreqs, nx, ny)
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

    ! Save results
    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_slicexy", p_reshaped)
    CALL save_variable(self % save_path, "Obs_slicexy", obs_reshaped)
    CALL save_variable(self % save_path, "X_slicexy", xs)
    CALL save_variable(self % save_path, "Y_slicexy", ys)
    CALL save_variable(self % save_path, "Z_slicexy", z)
    CALL save_variable(self % save_path, "Center_slicexy", center)
    CALL save_variable(self % save_path, "Nx_slicexy", nx)
    CALL save_variable(self % save_path, "Ny_slicexy", ny)

    WRITE(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_sliceXY


! -------------------------------------------------------------------------
! Subroutine: run_sliceXZ
! Computes pressure on a vertical (x,z) grid at fixed y_xz, over the given
! x/z limits. Saves the reshaped pressure/observer grids and the slice
! geometry metadata to the solver's HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE run_sliceXZ(self, y_xz, nx_xz, nz_xz, xlim_xz, zlim_xz)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    REAL(WP)    , INTENT(IN), OPTIONAL  :: y_xz
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nx_xz
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nz_xz
    REAL(WP)    , INTENT(IN), OPTIONAL  :: xlim_xz(2)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: zlim_xz(2)

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    COMPLEX(WP), ALLOCATABLE :: p_reshaped(:,:,:)   ! (Nfreqs, nx, nz)
    REAL(WP)   , ALLOCATABLE :: observers(:,:)
    REAL(WP)   , ALLOCATABLE :: obs_reshaped(:,:,:) ! (nx, nz, 3)
    REAL(WP)   , ALLOCATABLE :: xs(:), zs(:)        ! Coordinate arrays
    REAL(WP)                 :: y, xlim(2), zlim(2), dx, dz, x_val, z_val, Depth
    INTEGER(I32)             :: nx, nz, i, j, idx, bs

    Depth = self % env % depth_base
    bs    = self % block_size

    ! Defaults
    y    = self % turbine_data % AxisPos(2); if (PRESENT(y_xz))    y    = y_xz
    nx   = 26_I32                       ; if (PRESENT(nx_xz))   nx   = nx_xz
    nz   = 26_I32                       ; if (PRESENT(nz_xz))   nz   = nz_xz
    xlim = [-500.0_WP, 500.0_WP]        ; if (PRESENT(xlim_xz)) xlim = xlim_xz
    zlim = [0.0_WP, - Depth]            ; if (PRESENT(zlim_xz)) zlim = zlim_xz

    ALLOCATE(observers(nx * nz, 3))
    ALLOCATE(xs(nx), zs(nz))
    dx = (xlim(2) - xlim(1)) / max(REAL(nx - 1, WP), 1.0_WP)
    dz = (zlim(2) - zlim(1)) / max(REAL(nz - 1, WP), 1.0_WP)

    ! Fill coordinate arrays
    do j = 1, nx
        xs(j) = xlim(1) + REAL(j - 1, WP) * dx
    end do
    do i = 1, nz
        zs(i) = zlim(1) + REAL(i - 1, WP) * dz
    end do

    ! Build observers
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

    WRITE(*,'(A,I0,A,I0,A,F8.2,A,F8.2,A,F8.2,A)') &
        '>> [SliceXZ]    nx=', nx, ' | nz=', nz, ' | y=', y, &
        ' m | z from ', zlim(1), ' to ', zlim(2), ' m'

    CALL self % compute_pressure(observers, bs, p)
    if (self % debug) WRITE(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    ! Reshape pressure and observers to 3D (nx, nz, 3) and (Nfreqs, nx, nz)
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

    ! Save results
    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_slicexz", p_reshaped)
    CALL save_variable(self % save_path, "Obs_slicexz", obs_reshaped)
    CALL save_variable(self % save_path, "X_slicexz", xs)
    CALL save_variable(self % save_path, "Z_slicexz", zs)
    CALL save_variable(self % save_path, "Y_slicexz", y)
    CALL save_variable(self % save_path, "Nx_slicexz", nx)
    CALL save_variable(self % save_path, "Nz_slicexz", nz)

    WRITE(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_sliceXZ


! -------------------------------------------------------------------------
! Subroutine: dummy_run_sphere
! Default (no-op) implementation of run_sphere, used by solvers that do
! not support a spherical free-field observer geometry. Only consumes its
! arguments and prints a "skipped" notice.
! -------------------------------------------------------------------------
SUBROUTINE dummy_run_sphere(self, r_sphere, n_theta_sphere, nz_sphere, center_sphere)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    REAL(WP)    , INTENT(IN), OPTIONAL  :: r_sphere
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_theta_sphere
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nz_sphere
    REAL(WP)    , INTENT(IN), OPTIONAL  :: center_sphere(3)

    ! local variables not to have interpreter warnings
    REAL(WP) :: r, center(3)
    INTEGER(I32) ::n_theta, nz

    
    r = 30.0_WP     ; if (PRESENT(r_sphere))       r       = r_sphere
    n_theta = 72_I32; if (PRESENT(n_theta_sphere)) n_theta = n_theta_sphere
    nz = 20_I32     ; if (PRESENT(nz_sphere))      nz      = nz_sphere
    center = 0.0_WP ; if (PRESENT(center_sphere))  center  = center_sphere

    ! Does nothing only skips
    WRITE(*,'(A)') '>> [Sphere]     skipped, not available for: ' // trim(self % solver_name)

END SUBROUTINE dummy_run_sphere


! -------------------------------------------------------------------------
! Subroutine: run_all
! Orchestrates the full acoustic post-processing pipeline: runs, in order,
! spectrums, polar, decay, line, sphere, sliceXY, sliceXZ and cylinder,
! forwarding each routine's optional arguments and saving every result to
! the solver's HDF5 output file.
! -------------------------------------------------------------------------
SUBROUTINE run_all(self, &
    observers_spectrum, z_spectrums, &
    r_polar, z_polar, n_polar, center_polar, &
    r_cyl, n_cyl, nz_cyl, center_cyl, &
    distance_decay, n_decay, z_decay, logspace_decay, &
    p1_line, p2_line, n_line, logspace_line, &
    z_xy, nx_xy, ny_xy, xlim_xy, ylim_xy, center_xy, &
    y_xz, nx_xz, nz_xz, xlim_xz, zlim_xz, &
    r_sphere, n_theta_sphere, nz_sphere, center_sphere)
        
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self

    ! --- spectrums ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: observers_spectrum(:,:)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: z_spectrums
    ! --- polar ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: r_polar, z_polar, center_polar(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_polar
    ! --- cylinder ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: r_cyl, center_cyl(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_cyl, nz_cyl
    ! --- decay ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: distance_decay(2), z_decay
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_decay
    LOGICAL     , INTENT(IN), OPTIONAL  :: logspace_decay
    ! --- line ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: p1_line(3), p2_line(3)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_line
    LOGICAL     , INTENT(IN), OPTIONAL  :: logspace_line
    ! --- sliceXY ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: z_xy, xlim_xy(2), ylim_xy(2), center_xy(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nx_xy, ny_xy
    ! --- sliceXZ ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: y_xz, xlim_xz(2), zlim_xz(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nx_xz, nz_xz
    ! --- sphere ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: r_sphere, center_sphere(3)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_theta_sphere, nz_sphere


    WRITE(*,'(A)') ''
    WRITE(*,'(A)') repeat('=', 70)
    WRITE(*,'(A)') ' RUNNING ALL ACOUSTIC POST-PROCESSING (8 steps)'
    WRITE(*,'(A)') repeat('=', 70)

    WRITE(*,'(A)') ''
    WRITE(*,'(A)') '[1/8] Spectrums at observer points'
    CALL self % run_spectrums(observers_spectrum, z_spectrums)

    WRITE(*,'(A)') ''
    WRITE(*,'(A)') '[2/8] Polar contour'
    CALL self % run_polar(r_polar, z_polar, n_polar, center_polar)

    WRITE(*,'(A)') ''
    WRITE(*,'(A)') '[3/8] Distance decay along wind'
    CALL self % run_decay(distance_decay, n_decay, z_decay, logspace_decay)

    WRITE(*,'(A)') ''
    WRITE(*,'(A)') '[4/8] Line between two points'
    CALL self % run_line(p1_line, p2_line, n_line, logspace_line)

    WRITE(*,'(A)') ''
    WRITE(*,'(A)') '[5/8] Spherical radiation (free-field)'
    CALL self % run_sphere(r_sphere, n_theta_sphere, nz_sphere, center_sphere)

    WRITE(*,'(A)') ''
    WRITE(*,'(A)') '[6/8] Horizontal slice XY'
    CALL self % run_sliceXY(z_xy, nx_xy, ny_xy, xlim_xy, ylim_xy, center_xy)

    WRITE(*,'(A)') ''
    WRITE(*,'(A)') '[7/8] Vertical slice XZ'
    CALL self % run_sliceXZ(y_xz, nx_xz, nz_xz, xlim_xz, zlim_xz)

    WRITE(*,'(A)') ''
    WRITE(*,'(A)') '[8/8] Cylinder surface'
    CALL self % run_cylinder(r_cyl, n_cyl, nz_cyl, center_cyl)

    WRITE(*,'(A)') ''
    WRITE(*,'(A)') repeat('=', 70)
    WRITE(*,'(A)') ' ALL ACOUSTIC COMPUTATIONS FINISHED'
    WRITE(*,'(A)') repeat('=', 70)

END SUBROUTINE run_all


END MODULE AcousticSolver