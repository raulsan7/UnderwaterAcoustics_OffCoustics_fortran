MODULE AcousticSolver

USE Kinds
USE Types

IMPLICIT NONE

PRIVATE
PUBLIC :: AcousticSolver_t


! ------------------
! Abstract Base Type
! ------------------
TYPE, ABSTRACT :: AcousticSolver_t

    TYPE(AcousticSource_t)     :: sources
    TYPE(Config_Environment_t) :: env
    
    CHARACTER(len=50)          :: solver_name = 'None'           ! [-] Solver name
    CHARACTER(len=512)         :: save_path   = 'None'           ! [-] Path to save acoustics


    CONTAINS

    ! --- Public type-bound procedures --- !
    PROCEDURE :: save_acoustics
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
    PROCEDURE(get_name_i)        , DEFERRED :: get_name
    PROCEDURE(compute_pressure_i), DEFERRED :: compute_pressure
    PROCEDURE(save_parameters_i) , DEFERRED :: save_parameters

END TYPE AcousticSolver_t


! --------------------------------------------------
! Abstract interfaces for deferred(abstract) methods
! --------------------------------------------------
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
    SUBROUTINE compute_pressure_i(self, sources, observers, block_size, total_pressure)
        IMPORT AcousticSolver_t, AcousticSource_t, WP, I32
        CLASS(AcousticSolver_t) , INTENT(INOUT)        :: self
        TYPE(AcousticSource_t) , INTENT(IN)           :: sources
        REAL(WP)                , INTENT(IN)           :: observers(:,:)        ! [m] Observers coords (Nobs,3)
        INTEGER(I32)            , INTENT(IN), OPTIONAL :: block_size            ! [Pa] Pressure field (Nfreqs, Nobs)
        COMPLEX(WP), ALLOCATABLE, INTENT(OUT)          :: total_pressure(:,:)   ! [-] Number of observers to process simultaneously
    END SUBROUTINE compute_pressure_i
END INTERFACE

! Save parameters
ABSTRACT INTERFACE
    SUBROUTINE save_parameters_i(self)
        IMPORT:: AcousticSolver_t
        CLASS(AcousticSolver_t), INTENT(IN) :: self
    END SUBROUTINE
END INTERFACE


CONTAINS


! ---------- HELPERS ---------- !
SUBROUTINE save_acoustics(Self, var_name, var)
    USE IOUtils, ONLY: save_to_hdf5

    CLASS(AcousticSolver_t), INTENT(IN) :: self
    CHARACTER(len=*)       , INTENT(IN) :: var_name
    CLASS(*)               , INTENT(IN) :: var(..)

    ! Local variables
    CHARACTER(len=512) :: save_path

    save_path = trim(self % save_path)
    
    SELECT RANK (v_rank => var)
    RANK(0) 
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP)) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP)); CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (LOGICAL) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        END SELECT
    RANK(1) 
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP)) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP)); CALL save_to_hdf5(save_path, trim(var_name), v)
        END SELECT
    RANK(2) 
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP)) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP)); CALL save_to_hdf5(save_path, trim(var_name), v)
        END SELECT
    RANK(3) 
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP)) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP)); CALL save_to_hdf5(save_path, trim(var_name), v)
        END SELECT
    END SELECT

END SUBROUTINE save_acoustics


SUBROUTINE check_observers_distances(self, sources, observers, min_distance)
    CLASS(AcousticSolver_t), INTENT(IN) :: self
    TYPE(AcousticSource_t) , INTENT(IN) :: sources
    REAL(WP)               , INTENT(IN) :: observers(:,:)          ! (Nobs, 3)
    REAL(WP)               , INTENT(IN), OPTIONAL :: min_distance  ! [m]

    ! Local variables
    REAL(WP) :: min_dist
    INTEGER(I32) :: i, j, n_nodes, n_obs
    REAL(WP) :: dx, dy, dz, dist

    ! Default min value
    min_dist = 1.0_WP; if (PRESENT(min_distance)) min_dist = min_distance
    
    n_nodes = size(sources % xyz, 1)
    n_obs   = size(observers, 1)
    
    do i = 1, n_obs
        do j = 1, n_nodes
            dx = sources % xyz(j,1) - observers(i,1)
            dy = sources % xyz(j,2) - observers(i,2)
            dz = sources % xyz(j,3) - observers(i,3)
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
SUBROUTINE run_spectrums(self, sources, config, observers, z_spectrums)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    TYPE(AcousticSource_t) , INTENT(IN)    :: sources
    TYPE(Config_t)         , INTENT(IN)    :: config

    REAL(WP), INTENT(IN), OPTIONAL :: observers(:,:)
    REAL(Wp), INTENT(IN), OPTIONAL :: z_spectrums

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)                 :: z_spectrums_
    INTEGER(I32) :: bs

    print*, ''; print*, 'Computing spectrums at observer points...'

    ! Defaults
    z_spectrums_ = - config % turbine % Depth/2.0_WP; if (PRESENT(z_spectrums)) z_spectrums_ = z_spectrums
    
    if (PRESENT(observers)) then
        observers_ = observers
    else
        ALLOCATE(observers_(3,3))
        observers_(:,3) = z_spectrums_
        observers_(:,2) = 0.0_WP
        observers_(:,1) = [10.0_WP, 250.0_WP, 500.0_WP]
    end if

    bs = config % solver % block_size

    ! Compute pressure
    CALL self % check_observers_distances(sources, observers_)
    CALL self % compute_pressure(sources, observers_, bs, p)

    if (config % solver % debug) write(*, '(A, F0.2)') 'Spectrum Pressure norm: ', norm2(abs(p))

    ! Save data
    CALL self % save_acoustics("P_spectrums", p)
    CALL self % save_acoustics("Obs_spectrums", observers_)
    CALL self % save_acoustics("Freqs", sources % Freqs)

    print*, ''
    print*, '   ---> Spectrum data saved at ', self % save_path

END SUBROUTINE run_spectrums


SUBROUTINE run_polar(self, sources, config, r_polar, z_polar, n_polar, center_polar)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    TYPE(AcousticSource_t) , INTENT(IN)    :: sources
    TYPE(Config_t)         , INTENT(IN)    :: config

    REAL(WP)    , INTENT(IN), OPTIONAL  :: r_polar
    REAL(WP)    , INTENT(IN), OPTIONAL  :: z_polar
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_polar
    REAL(WP)    , INTENT(IN), OPTIONAL  :: center_polar(2)

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers(:,:), theta_deg_polar(:)
    REAL(WP)                 :: r, z, center(2), theta_deg, theta_rad, d_theta, Depth
    INTEGER(I32)             :: n_theta, bs, i

    
    Depth = config % env % depth_base
    bs    = config % solver % block_size

    ! Defaults
    r = 500.0_WP      ; if (PRESENT(r_polar)) r       = r_polar
    z = - Depth/2.0_WP; if (PRESENT(z_polar)) z       = z_polar
    n_theta = 72_I32  ; if (PRESENT(n_polar)) n_theta = n_theta

    if (PRESENT(center_polar)) then
        center = center_polar
    else
        center = [config % turbine % BariPos(1), config % turbine % BariPos(2)]
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

    CALL self % check_observers_distances(sources, observers)

    write(*, '(A)') ''
    write(*, '(A, F0.2, A, F0.2, A, F0.2, A, F0.2, A, I0, A)') &
        'Computing polar pressure (r=', r, ' m, center=(', center(1), ', ', center(2), &
        ') m, z=', z, ' m), observers=', n_theta, '...'
        
    ! Compute_pressure
    CALL self % compute_pressure(sources, observers, bs, p)
    if (config % solver % debug) write(*, '(A, F0.2)') 'Polar Pressure Norm: ', norm2(abs(p))

    ! Save data
    CALL self % save_acoustics("Freqs", sources % Freqs)
    CALL self % save_acoustics("P_polar", p)
    CALL self % save_acoustics("R_polar", r)
    CALL self % save_acoustics("Z_polar", z)
    CALL self % save_acoustics("Theta_deg_polar", theta_deg_polar)
    CALL self % save_acoustics("Obs_polar", observers)
    CALL self % save_acoustics("Center_polar", center)

    print*, ''
    print*, '   ---> Polar data saved at ', self % save_path

END SUBROUTINE run_polar


SUBROUTINE run_cylinder(self, sources, config, r_cyl, n_cyl, nz_cyl, center_cyl)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    TYPE(AcousticSource_t) , INTENT(IN)    :: sources
    TYPE(Config_t)         , INTENT(IN)    :: config

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


    Depth = config % env % density_base
    bs    = config % solver % block_size

    ! Defaults
    r = 500.0_WP    ; if (PRESENT(r_cyl))  r       = r_cyl
    n_theta = 72_I32; if (PRESENT(n_cyl))  n_theta = n_cyl
    nz = 20_I32     ; if (PRESENT(nz_cyl)) nz      = nz_cyl

    if (PRESENT(center_cyl)) then
        center = center_cyl
    else
        center = [config % turbine % BariPos(1), config % turbine % BariPos(2)]
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

    CALL self % check_observers_distances(sources, observers)

    write(*, '(A)') ''
    write(*, '(A, F0.2, A, F0.2, A, F0.2, A, I0, A, I0, A, I0, A)') &
        'Computing cylindrical pressure (r=', r, ' m, center=(', center(1), ', ', center(2), &
        ') m, grid=', n_theta, ' x ', nz, ' = ', Nobs, ' obs)...'

    CALL self % compute_pressure(sources, observers, bs, p)
    if (config % solver % debug) write(*, '(A, F0.2)') 'Cylinder Pressure Norm: ', norm2(abs(p))

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
    CALL self % save_acoustics("Freqs", sources % Freqs)
    CALL self % save_acoustics("P_cylinder", p_reshaped)
    CALL self % save_acoustics("R_cylinder", r)
    CALL self % save_acoustics("Z_cylinder", z_cylinder)
    CALL self % save_acoustics("Theta_deg_cylinder", theta_deg_arr)
    CALL self % save_acoustics("Obs_cylinder", obs_reshaped)
    CALL self % save_acoustics("Center_cylinder", center)
    CALL self % save_acoustics("dA_cylinder", dA)

    print*, ''
    print*, '   ---> Cylinder data saved at ', self % save_path

END SUBROUTINE run_cylinder


SUBROUTINE run_decay(self, sources, config, distance_decay, n_decay, z_decay, logspace_decay)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    TYPE(AcousticSource_t) , INTENT(IN)    :: sources
    TYPE(Config_t)         , INTENT(IN)    :: config

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


    Depth = config % env % depth_base
    bs    = config % solver % block_size

    ! Defaults
    z = -Depth/2.0_WP             ; if (PRESENT(z_decay))        z        = z_decay
    distance = [10.0_WP, 500.0_WP]; if (PRESENT(distance_decay)) distance = distance_decay
    n_points = 200_I32            ; if (PRESENT(n_decay))        n_points = n_decay
    logspace = .true.             ; if (PRESENT(logspace_decay)) logspace = logspace_decay
    
    ! Wind direction
    WindDir_rad = config % turbine % WindDir * PI / 180.0_WP
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
            observers(i, 1:2) = config % turbine % AxisPos(1:2) + d_val * wind_unit
            observers(i, 3)   = z
        end do
    else
        linear_step = (distance(2) - distance(1)) / REAL(n_points - 1, WP)
        do i = 1, n_points
            d_val = distance(1) + REAL(i - 1, WP) * linear_step
            distances(i)      = d_val
            observers(i, 1:2) = config % turbine % AxisPos(1:2) + d_val * wind_unit
            observers(i, 3)   = z
        end do
    end if

    CALL self % check_observers_distances(sources, observers)

    write(*, '(A)') ''
    write(*, '(A, I0, A, F0.2, A, F0.2, A)') &
        'Computing distance decay (points=', n_points, ', direction=', &
        config % turbine % WindDir, ' deg, z=', z, ' m)...'

    CALL self % compute_pressure(sources, observers, bs, p)
    if (config % solver % debug) write(*, '(A, F0.2)') 'Decay Pressure Norm: ', norm2(abs(p))

    ! Convert logical to integer for saving
    logspace_int = merge(1, 0, logspace)

    ! Save results
    CALL self % save_acoustics("Freqs", sources % Freqs)
    CALL self % save_acoustics("P_decay", p)
    CALL self % save_acoustics("Z_decay", z)
    CALL self % save_acoustics("Distance_decay", distance)
    CALL self % save_acoustics("Distances_decay", distances)
    CALL self % save_acoustics("Obs_decay", observers)
    CALL self % save_acoustics("Logspace_decay", logspace_int)

    print*, ''
    print*, '   ---> Decay data saved at ', self % save_path

END SUBROUTINE run_decay


SUBROUTINE run_line(self, sources, config, p1_line, p2_line, n_line, logspace_line)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    TYPE(AcousticSource_t) , INTENT(IN)    :: sources
    TYPE(Config_t)         , INTENT(IN)    :: config

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


    Depth = config % env % depth_base
    bs    = config % solver % block_size

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

    bs = int(n_points / 10_I32, I32)
    if (bs < 1_I32) bs = 1_I32

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

    CALL self % check_observers_distances(sources, observers)

    write(*, '(A)') ''
    write(*, '(A, I0, A, F0.2, A, F0.2, A, F0.2, A, F0.2, A, F0.2, A, F0.2, A)') &
        'Computing line (points=', n_points, '), p1=(', p1(1), ', ', p1(2), ', ', p1(3), &
        ') m, p2=(', p2(1), ', ', p2(2), ', ', p2(3), ') m)...'

    CALL self % compute_pressure(sources, observers, bs, p)
    if (config % solver % debug) write(*, '(A, F0.2)') 'Line Pressure Norm: ', norm2(abs(p))

    ! Convert logical to integer for saving
    logspace_int = merge(1, 0, logspace)

    ! Save results
    CALL self % save_acoustics("Freqs", sources % Freqs)
    CALL self % save_acoustics("P_line", p)
    CALL self % save_acoustics("P1_line", p1)
    CALL self % save_acoustics("P2_line", p2)
    CALL self % save_acoustics("Obs_line", observers)
    CALL self % save_acoustics("Distance_line", dist_total)
    CALL self % save_acoustics("Distances_line", distances)
    CALL self % save_acoustics("Logspace_line", logspace_int)

    print*, ''
    print*, '   ---> Line data saved at ', self % save_path

END SUBROUTINE run_line


SUBROUTINE run_sliceXY(self, sources, config, z_xy, nx_xy, ny_xy, xlim_xy, ylim_xy, center_xy)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    TYPE(AcousticSource_t) , INTENT(IN)    :: sources
    TYPE(Config_t)         , INTENT(IN)    :: config

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


    Depth = config % env % depth_base
    bs   = config % solver % block_size

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

    CALL self % check_observers_distances(sources, observers, 0.1_WP)

    write(*, '(A)') ''
    write(*, '(A, I0, A, I0, A, F0.2, A, F0.2, A, F0.2, A)') &
        'Computing XY slice (nx=', nx, ', ny=', ny, ', z=', z, &
        ' m, center=(', center(1), ', ', center(2), ') m)...'

    CALL self % compute_pressure(sources, observers, bs, p)
    if (config % solver % debug) write(*, '(A, F0.2)') 'SliceXY Pressure Norm: ', norm2(abs(p))

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
    CALL self % save_acoustics("Freqs", sources % Freqs)
    CALL self % save_acoustics("P_slicexy", p_reshaped)
    CALL self % save_acoustics("Obs_slicexy", obs_reshaped)
    CALL self % save_acoustics("X_slicexy", xs)
    CALL self % save_acoustics("Y_slicexy", ys)
    CALL self % save_acoustics("Z_slicexy", z)
    CALL self % save_acoustics("Center_slicexy", center)
    CALL self % save_acoustics("Nx_slicexy", nx)
    CALL self % save_acoustics("Ny_slicexy", ny)

    print*, ''
    print*, '   ---> SliceXY data saved at ', self % save_path

END SUBROUTINE run_sliceXY


SUBROUTINE run_sliceXZ(self, sources, config, y_xz, nx_xz, nz_xz, xlim_xz, zlim_xz)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    TYPE(AcousticSource_t) , INTENT(IN)    :: sources
    TYPE(Config_t)         , INTENT(IN)    :: config

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

    Depth = config % env % depth_base
    bs    = config % solver % block_size

    ! Defaults
    y    = config % turbine % AxisPos(2); if (PRESENT(y_xz))    y    = y_xz
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

    CALL self % check_observers_distances(sources, observers, 0.1_WP)

    write(*, '(A)') ''
    write(*, '(A, I0, A, I0, A, F0.2, A, F0.2, A, F0.2, A)') &
        'Computing XZ slice (nx=', nx, ', nz=', nz, ', y=', y, &
        ' m, z from ', zlim(1), ' to ', zlim(2), ' m)...'

    CALL self % compute_pressure(sources, observers, bs, p)
    if (config % solver % debug) write(*, '(A, F0.2)') 'SliceXZ Pressure Norm: ', norm2(abs(p))

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
    CALL self % save_acoustics("Freqs", sources % Freqs)
    CALL self % save_acoustics("P_slicexz", p_reshaped)
    CALL self % save_acoustics("Obs_slicexz", obs_reshaped)
    CALL self % save_acoustics("X_slicexz", xs)
    CALL self % save_acoustics("Z_slicexz", zs)
    CALL self % save_acoustics("Y_slicexz", y)
    CALL self % save_acoustics("Nx_slicexz", nx)
    CALL self % save_acoustics("Nz_slicexz", nz)

    print*, ''
    print*, '   ---> SliceXZ data saved at ', self % save_path

END SUBROUTINE run_sliceXZ


SUBROUTINE dummy_run_sphere(self, sources, config, r_sphere, n_theta_sphere, nz_sphere, center_sphere)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    TYPE(AcousticSource_t) , INTENT(IN)    :: sources
    TYPE(Config_t)         , INTENT(IN)    :: config

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
    print*, '  [Skipped] Spherical radiation not available for: ', trim(self % solver_name)

END SUBROUTINE dummy_run_sphere


SUBROUTINE run_all(self, sources, config, &
    spectrums_observers, spectrums_z_obs, &
    polar_r, polar_z, polar_n_theta, polar_center, &
    cylinder_r, cylinder_n_theta, cylinder_nz, cylinder_center, &
    decay_distance, decay_n_points, decay_z, decay_logspace, &
    line_p1, line_p2, line_n_points, line_logspace, &
    sliceXY_z, sliceXY_nx, sliceXY_ny, sliceXY_xlim, sliceXY_ylim, sliceXY_center, &
    sliceXZ_y, sliceXZ_nx, sliceXZ_nz, sliceXZ_xlim, sliceXZ_zlim, &
    sphere_r, sphere_n_theta, sphere_nz, sphere_center)
        
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    TYPE(AcousticSource_t) , INTENT(IN)    :: sources
    TYPE(Config_t)         , INTENT(IN)    :: config

    ! --- spectrums ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: spectrums_observers(:,:)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: spectrums_z_obs
    ! --- polar ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: polar_r, polar_z, polar_center(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: polar_n_theta
    ! --- cylinder ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: cylinder_r, cylinder_center(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: cylinder_n_theta, cylinder_nz
    ! --- decay ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: decay_distance(2), decay_z
    INTEGER(I32), INTENT(IN), OPTIONAL  :: decay_n_points
    LOGICAL     , INTENT(IN), OPTIONAL  :: decay_logspace
    ! --- line ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: line_p1(3), line_p2(3)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: line_n_points
    LOGICAL     , INTENT(IN), OPTIONAL  :: line_logspace
    ! --- sliceXY ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: sliceXY_z, sliceXY_xlim(2), sliceXY_ylim(2), sliceXY_center(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: sliceXY_nx, sliceXY_ny
    ! --- sliceXZ ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: sliceXZ_y, sliceXZ_xlim(2), sliceXZ_zlim(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: sliceXZ_nx, sliceXZ_nz
    ! --- sphere ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: sphere_r, sphere_center(3)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: sphere_n_theta, sphere_nz


    print*, ''
    print*, '=================================================='
    print*, ' RUNNING ALL ACOUSTIC POST-PROCESSING'
    print*, ' Total steps: 8'
    print*, '=================================================='

    print*, ''
    print*, '[Step 1/8] Spectrums at observer points ...'
    CALL self % run_spectrums(sources, config, spectrums_observers, spectrums_z_obs)

    print*, ''
    print*, '[Step 2/8] Polar contour ...'
    CALL self % run_polar(sources, config, polar_r, polar_z, polar_n_theta, polar_center)

    print*, ''
    print*, '[Step 3/8] Distance decay along wind ...'
    CALL self % run_decay(sources, config, decay_distance, decay_n_points, decay_z, decay_logspace)

    print*, ''
    print*, '[Step 4/8] Line between two points ...'
    CALL self % run_line(sources, config, line_p1, line_p2, line_n_points, line_logspace)

    print*, ''
    print*, '[Step 5/8] Spherical radiation (free-field) ...'
    CALL self % run_sphere(sources, config, sphere_r, sphere_n_theta, sphere_nz, sphere_center)

    print*, ''
    print*, '[Step 6/8] Horizontal slice XY ...'
    CALL self % run_sliceXY(sources, config, sliceXY_z, sliceXY_nx, sliceXY_ny, sliceXY_xlim, sliceXY_ylim, sliceXY_center)

    print*, ''
    print*, '[Step 7/8] Vertical slice XZ ...'
    CALL self % run_sliceXZ(sources, config, sliceXZ_y, sliceXZ_nx, sliceXZ_nz, sliceXZ_xlim, sliceXZ_zlim)

    print*, ''
    print*, '[Step 8/8] Cylinder surface ...'
    CALL self % run_cylinder(sources, config, cylinder_r, cylinder_n_theta, cylinder_nz, cylinder_center)

    print*, ''
    print*, '=================================================='
    print*, ' ALL ACOUSTIC COMPUTATIONS FINISHED'
    print*, '=================================================='

END SUBROUTINE run_all


END MODULE AcousticSolver