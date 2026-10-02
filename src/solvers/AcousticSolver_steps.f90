! ===============================================================================
! SUBMODULE: AcousticSolver_steps
! ===============================================================================
! Defines the post-processing steps of AcousticSolver_t: the generic
! observer-geometry routines (run_spectrums, run_polar, run_cylinder,
! run_decay, run_line, run_sliceXY, run_sliceXZ, run_sliceVertical, and the
! no-op sphere of the solvers without one), and run(), the config-driven
! entry point that executes every enabled step in a fixed order.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
SUBMODULE (AcousticSolver) AcousticSolver_steps

USE HDF5Utils, ONLY: save_variable

IMPLICIT NONE


CONTAINS


! ===============================================================================
! RUN SUBROUTINES
! ===============================================================================
MODULE SUBROUTINE run_spectrums(self)
    ! -------------------------------------------------------------------------
    ! Computes the pressure spectrum at a small set of fixed observer points
    ! (default: 3 points along +x at mid-depth). Saves the pressure, observer
    ! coordinates and frequency vector to the solver's HDF5 file.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: p(:,:)
    REAL(WP)    , ALLOCATABLE :: observers_(:,:)
    REAL(WP)                  :: z_top, z_bot
    INTEGER(I32)              :: bs, i


    write(*, '(A)') '>> [Spectrums] computing pressure at observer points...'

    ! Defaults: the 3 automatic points are not persisted, update_parameters describes them itself
    associate(pp => self % run_inputs % params_spectrums)
        if (pp % n_observers == DEFAULT_I) then
            ALLOCATE(observers_(3, 3))
            observers_(:, 1) = [10.0_WP, 250.0_WP, 500.0_WP]
            observers_(:, 2) = 0.0_WP
            do i = 1, 3
                z_top            = self % z_surface(observers_(i, 1), observers_(i, 2))
                z_bot            = self % z_seabed(observers_(i, 1), observers_(i, 2))
                observers_(i, 3) = 0.5_WP * (z_top + z_bot)    ! Mid-depth at every default point
            end do
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


MODULE SUBROUTINE run_polar(self)
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
    REAL(WP)                  :: r, z, center(2), theta_deg, theta_rad, d_theta, z_top, z_bot
    INTEGER(I32)              :: n_theta, bs, i


    bs = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_polar)
        if (pp % r_polar == DEFAULT_R) pp % r_polar = 500.0_WP
        if (pp % n_polar == DEFAULT_I) pp % n_polar = 72_I32
        if (all(pp % center_polar == DEFAULT_R)) pp % center_polar = self % sources % geometry % center

        ! Water column at the circle centre
        z_top = self % z_surface(pp % center_polar(1), pp % center_polar(2))
        z_bot = self % z_seabed(pp % center_polar(1), pp % center_polar(2))
        if (pp % z_polar == DEFAULT_R) pp % z_polar = 0.5_WP * (z_top + z_bot)

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


MODULE SUBROUTINE run_cylinder(self)
    ! -------------------------------------------------------------------------
    ! Computes pressure on a vertical cylindrical surface (radius r_cyl,
    ! n_cyl azimuth points x nz_cyl depth levels) around the turbine
    ! barycenter. Saves the reshaped pressure/observer grids and the cylinder
    ! geometry metadata to the solver's HDF5 file. The depth levels include
    ! both ends (z_seabed and z_surface at the axis, z = -Depth and z = 0 by
    ! default), so dA_cylinder is the area of an interior cell: a surface
    ! integral must weight the first and last levels by 1/2 (trapezoid rule).
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: p(:,:)
    REAL(WP)    , ALLOCATABLE :: observers(:,:)
    REAL(WP)                  :: r, center(2), theta_rad, d_theta, z_val, dz, z_top, z_bot
    INTEGER(I32)              :: n_theta, nz, bs, i, j, idx, Nobs

    REAL(WP)   , ALLOCATABLE :: theta_deg_arr(:)
    REAL(WP)   , ALLOCATABLE :: z_cylinder(:)
    REAL(WP)   , ALLOCATABLE :: obs_reshaped(:,:,:)
    COMPLEX(WP), ALLOCATABLE :: p_reshaped(:,:,:)
    REAL(WP)                 :: dA


    bs = self % block_size

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

    ! Water column at the cylinder axis
    z_top = self % z_surface(center(1), center(2))
    z_bot = self % z_seabed(center(1), center(2))

    Nobs = n_theta * nz
    ALLOCATE(observers(Nobs, 3))

    d_theta = 360.0_WP / real(n_theta, WP)
    dz      = (z_top - z_bot) / max(real(nz - 1_I32, WP), 1.0_WP)

    ALLOCATE(theta_deg_arr(n_theta))
    ALLOCATE(z_cylinder(nz))

    do j = 1, nz
        z_cylinder(j) = z_bot + real(j - 1_I32, WP) * dz
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

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: p(:,:)
    REAL(WP)    , ALLOCATABLE :: observers(:,:)
    REAL(WP)    , ALLOCATABLE :: distances(:)
    REAL(WP)                  :: distance(2), z, wind_unit(2), d_val, WindDir_rad
    REAL(WP)                  :: d_log1, d_log2, d_step, linear_step, z_top, z_bot
    INTEGER(I32)              :: n_points, bs, i
    LOGICAL                   :: logspace
    INTEGER(I32)              :: logspace_int


    bs = self % block_size

    ! Water column at the start of the line (turbine axis or mean farm axis)
    z_top = self % z_surface(self % sources % geometry % origin(1), self % sources % geometry % origin(2))
    z_bot = self % z_seabed(self % sources % geometry % origin(1), self % sources % geometry % origin(2))

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_decay)
        if (pp % z_decay == DEFAULT_R)             pp % z_decay        = 0.5_WP * (z_top + z_bot)
        if (any(pp % distance_decay == DEFAULT_R)) pp % distance_decay = [10.0_WP, 500.0_WP]
        if (pp % n_decay == DEFAULT_I)             pp % n_decay        = 200_I32
        ! logspace_decay has no sentinel: it is read directly as a LOGICAL from the .nml

        if (pp % distance_decay(1) <= 0.0_WP) then
            write(*,*) "WARNING AcousticSolver % run_decay(): first point can't be placed at r = 0 m"
            pp % distance_decay(1) = 1.0_WP
        end if

        z        = pp % z_decay
        distance = pp % distance_decay
        n_points = pp % n_decay
        logspace = pp % logspace_decay
    end associate

    if (n_points < 1_I32) then
        write(*, '(A,I0)') " ERROR AcousticSolver % run_decay(): n_decay must be >= 1, got ", n_points
        error stop
    end if

    if (logspace .and. distance(2) <= 0.0_WP) then
        write(*, '(A,F0.3,A)') " ERROR AcousticSolver % run_decay(): log spacing needs distance_decay(2) > 0, got ", &
                               distance(2), " m"
        error stop
    end if

    WindDir_rad = self % sources % geometry % direction * PI / 180.0_WP
    wind_unit   = [cos(WindDir_rad), sin(WindDir_rad)]

    ALLOCATE(observers(n_points, 3))
    ALLOCATE(distances(n_points))

    ! max(): a single point sits at the start distance instead of 0/0 = NaN
    if (logspace) then
        d_log1 = log10(distance(1))
        d_log2 = log10(distance(2))
        d_step = (d_log2 - d_log1) / max(real(n_points - 1_I32, WP), 1.0_WP)
        do i = 1, n_points
            d_val             = 10.0_WP**(d_log1 + real(i - 1_I32, WP) * d_step)
            distances(i)      = d_val
            observers(i, 1:2) = self % sources % geometry % origin + d_val * wind_unit
            observers(i, 3)   = z
        end do
    else
        linear_step = (distance(2) - distance(1)) / max(real(n_points - 1_I32, WP), 1.0_WP)
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

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: p(:,:)
    REAL(WP)    , ALLOCATABLE :: observers(:,:)
    REAL(WP)    , ALLOCATABLE :: distances(:)
    REAL(WP)                  :: p1(3), p2(3), dir_vec(3), dist_total, d_val
    REAL(WP)                  :: d_log1, d_log2, d_step, linear_step, z_top, z_bot
    INTEGER(I32)              :: n_points, i, bs
    LOGICAL                   :: logspace
    INTEGER(I32)              :: logspace_int


    bs = self % block_size

    ! Water column at the default end points (x, y) = (100, 0) m
    z_top = self % z_surface(100.0_WP, 0.0_WP)
    z_bot = self % z_seabed(100.0_WP, 0.0_WP)

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_line)
        ! Per coordinate: only the components left at -999 take the default point's value
        WHERE (pp % p1_line == DEFAULT_R) pp % p1_line = [100.0_WP, 0.0_WP, z_top]
        WHERE (pp % p2_line == DEFAULT_R) pp % p2_line = [100.0_WP, 0.0_WP, z_bot]
        ! logspace_line has no sentinel: it is read directly as a LOGICAL from the .nml

        p1       = pp % p1_line
        p2       = pp % p2_line
        logspace = pp % logspace_line
    end associate

    dist_total = norm2(p2 - p1)
    if (dist_total == 0.0_WP) error stop "AcousticSolver.run_line(): points p1 and p2 are coincident."

    ! n_line depends on dist_total, so it is resolved separately (lines shorter than 2 m keep both end points)
    if (self % run_inputs % params_line % n_line == DEFAULT_I) then
        if (dist_total > 200.0_WP) then
            self % run_inputs % params_line % n_line = int(dist_total / 5.0_WP, I32)
        else
            self % run_inputs % params_line % n_line = max(2_I32, int(dist_total, I32))
        end if
    end if
    n_points = self % run_inputs % params_line % n_line

    if (n_points < 1_I32) then
        write(*, '(A,I0)') " ERROR AcousticSolver % run_line(): n_line must be >= 1, got ", n_points
        error stop
    end if

    ALLOCATE(observers(n_points, 3))
    ALLOCATE(distances(n_points))
    dir_vec = (p2 - p1) / dist_total    ! Unit vector from p1 to p2: observer i sits d_val metres from p1

    ! max(): a single point sits at the first spacing point instead of 0/0 = NaN
    if (logspace) then
        d_log1 = log10(dist_total * 1.0e-6_WP)
        d_log2 = log10(dist_total)
        d_step = (d_log2 - d_log1) / max(real(n_points - 1_I32, WP), 1.0_WP)
        do i = 1, n_points
            d_val           = 10.0_WP**(d_log1 + real(i - 1_I32, WP) * d_step)
            distances(i)    = d_val
            observers(i, :) = p1 + d_val * dir_vec
        end do
    else
        linear_step = dist_total / max(real(n_points - 1_I32, WP), 1.0_WP)
        do i = 1, n_points
            d_val           = real(i - 1_I32, WP) * linear_step
            distances(i)    = d_val
            observers(i, :) = p1 + d_val * dir_vec
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


MODULE SUBROUTINE run_sliceXY(self)
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
    REAL(WP)                  :: z, xlim(2), ylim(2), center(2), z_top, z_bot
    REAL(WP)                  :: mid_x, mid_y, dx, dy, x_val, y_val
    INTEGER(I32)              :: nx, ny, i, j, idx, bs
    REAL(WP)    , PARAMETER   :: tol = 1.0e-5_WP


    bs = self % block_size

    ! Water column at the emitter centre
    z_top = self % z_surface(self % sources % geometry % center(1), self % sources % geometry % center(2))
    z_bot = self % z_seabed(self % sources % geometry % center(1), self % sources % geometry % center(2))

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_sliceXY)
        if (pp % z_xy == DEFAULT_R) pp % z_xy = 0.5_WP * (z_top + z_bot)
        if (pp % z_xy > z_top) then
            write(*,*) "AcousticSolver.run_sliceXY(): z must be <= 0. Switching to -z"
            pp % z_xy = -pp % z_xy
        end if

        if (pp % nx_xy == DEFAULT_I)        pp % nx_xy   = 26_I32
        if (pp % ny_xy == DEFAULT_I)        pp % ny_xy   = 26_I32
        if (any(pp % xlim_xy == DEFAULT_R)) &
            pp % xlim_xy = self % sources % geometry % center(1) + [-500.0_WP, 500.0_WP]
        if (any(pp % ylim_xy == DEFAULT_R)) &
            pp % ylim_xy = self % sources % geometry % center(2) + [-500.0_WP, 500.0_WP]

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


MODULE SUBROUTINE run_sliceXZ(self)
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
    REAL(WP)                  :: y, xlim(2), zlim(2), dx, dz, x_val, z_val, z_top, z_bot
    INTEGER(I32)              :: nx, nz, i, j, idx, bs


    bs = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_sliceXZ)
        if (pp % y_xz  == DEFAULT_R)        pp % y_xz    = self % sources % geometry % origin(2)
        if (pp % nx_xz == DEFAULT_I)        pp % nx_xz   = 26_I32
        if (pp % nz_xz == DEFAULT_I)        pp % nz_xz   = 26_I32
        if (any(pp % xlim_xz == DEFAULT_R)) &
            pp % xlim_xz = self % sources % geometry % center(1) + [-500.0_WP, 500.0_WP]

        ! Water column at the middle of the slice
        z_top = self % z_surface(0.5_WP * (pp % xlim_xz(1) + pp % xlim_xz(2)), pp % y_xz)
        z_bot = self % z_seabed(0.5_WP * (pp % xlim_xz(1) + pp % xlim_xz(2)), pp % y_xz)
        WHERE (pp % zlim_xz == DEFAULT_R) pp % zlim_xz = [z_top, z_bot]    ! Per limit, as documented

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

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: p(:,:)
    COMPLEX(WP) , ALLOCATABLE :: p_reshaped(:,:,:)
    REAL(WP)    , ALLOCATABLE :: observers(:,:)
    REAL(WP)    , ALLOCATABLE :: obs_reshaped(:,:,:)
    REAL(WP)    , ALLOCATABLE :: us(:), zs(:), proj(:)
    REAL(WP)                  :: azimuth, az_rad, normal(2), span(2), center(2), ulim(2), zlim(2)
    REAL(WP)                  :: du, dz, z_top, z_bot, y_xz
    INTEGER(I32)              :: nu, nz, i, j, idx, bs


    bs = self % block_size

    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_sliceVertical)

        ! --- 1. Plane orientation: normal and in-plane horizontal (span) axis --- !
        azimuth = pp % azimuth_sv
        az_rad  = azimuth * PI / 180.0_WP
        normal  = [ cos(az_rad), sin(az_rad)]
        span    = [-sin(az_rad), cos(az_rad)]

        ! --- 2. Centre: turbine axis, or mean turbine axis for a wind farm --- !
        if (any(pp % center_sv == DEFAULT_R)) pp % center_sv = self % sources % geometry % origin

        ! Normal along +-y (azimuth = +-90 deg) through y = y_xz: same plane as run_sliceXZ, skipped
        if (self % run_inputs % do_run % do_sliceXZ .and. &
            abs(cos(az_rad)) < sin(1.0e-3_WP * PI / 180.0_WP)) then
            y_xz = self % run_inputs % params_sliceXZ % y_xz
            if (y_xz == DEFAULT_R) y_xz = self % sources % geometry % origin(2)    ! Same default as run_sliceXZ
            if (abs(y_xz - pp % center_sv(2)) < 1.0e-3_WP) then
                write(*, '(A,F8.2,A)') '>> [SliceVert]  skipped: same plane as SliceXZ (y = ', y_xz, ' m)'
                self % run_inputs % do_run % do_sliceVertical = .false.    ! Not run: left out of the .dat summary
                return
            end if
        end if

        ! --- 3. Depth extent, clipped to the water column at the plane centre --- !
        z_top = self % z_surface(pp % center_sv(1), pp % center_sv(2))
        z_bot = self % z_seabed(pp % center_sv(1), pp % center_sv(2))
        if (pp % zlim_sv(1) == DEFAULT_R) pp % zlim_sv(1) = z_bot
        if (pp % zlim_sv(2) == DEFAULT_R) pp % zlim_sv(2) = z_top
        pp % zlim_sv(1) = max(pp % zlim_sv(1), z_bot)
        pp % zlim_sv(2) = min(pp % zlim_sv(2), z_top)

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

    ! --- 7. Compute pressure (observers closer than mask_distance to a node get NaN) --- !
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


MODULE SUBROUTINE dummy_run_sphere(self)
    ! -------------------------------------------------------------------------
    ! Default (no-op) implementation of run_sphere, used by solvers that do
    ! not support a spherical free-field observer geometry. Prints a
    ! "skipped" notice and switches do_sphere off, so the .dat summary does
    ! not report a step that never ran.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self    ! [-] Solver instance


    write(*, '(A)') '>> [Sphere]     skipped, not available for: ' // trim(self % solver_name)
    self % run_inputs % do_run % do_sphere = .false.    ! Not run: left out of the .dat summary

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

END SUBMODULE AcousticSolver_steps
