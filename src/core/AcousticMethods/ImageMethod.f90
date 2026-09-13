MODULE ImageMethod

USE Kinds
USE Types
USE AcousticSolver, only: AcousticSolver_t
IMPLICIT NONE

PRIVATE
PUBLIC:: ImageMethod_t


! ------------------------
! Main Type Definition
! ------------------------
TYPE, EXTENDS(AcousticSolver_t):: ImageMethod_t

    INTEGER(I32) :: N_images         = 30           ! [-] Number of image levels to compute
    INTEGER(I32) :: default_N_images = 30
    INTEGER(I32) :: Upper_BC         = -1           ! [-] Upper Boundary condition -1 for  p = 0
    INTEGER(I32) :: Lower_BC         = +1           ! [-] Lower Boundary condition +1 for dp = 0
    REAL(WP)     :: Upper_HBC        = 0.0_WP       ! [m] z-coordinate of the upper reflecting plane (e.g. free surface)
    REAL(WP)     :: Lower_HBC        = -30.0_WP     ! [m] z-coordinate of the lower reflecting plane (e.g. seabed)
    REAL(WP)     :: up_R             = 1.0_WP       ! [-] Attenuation coefficient for upper boundary reflections
    REAL(WP)     :: lw_R             = 0.5_WP       ! [-] Attenuation coefficient for lower boundary reflections
 
    TYPE(ImageSystem_t) , ALLOCATABLE :: image_system(:)

    CONTAINS
    ! --- Public type-bound procedures --- !
    GENERIC   :: init        => init_array, init_single
    PROCEDURE :: init_array  => init_MethodImages
    PROCEDURE :: init_single => init_MethodImages_single
    PROCEDURE :: set_N_images
    PROCEDURE :: build_all_image_system
    PROCEDURE :: build_image_system
    PROCEDURE :: restore_default_images
    PROCEDURE :: method_of_images
    PROCEDURE :: dipole_pressure_images
    PROCEDURE :: run_sphere

    ! --- Deferred(abstract) procedures --- !
    PROCEDURE :: get_name           => get_name_MethodImages
    PROCEDURE :: compute_pressure   => compute_pressure_MethodImages
    PROCEDURE :: save_parameters    => save_parameters_MethodImages

END TYPE ImageMethod_t


CONTAINS


! ---------- DEFERRED ---------- !
SUBROUTINE init_MethodImages(self, sources, config, env, N_images, c_wat, rho_wat, &
                             Upper_BC, Lower_BC, Upper_HBC, Lower_HBC, &
                             attenuation_upper, attenuation_lower, eps, &
                             p_ref, cluster, debug, name)

    CLASS(MethodImages_t)     , INTENT(INOUT) :: self
    TYPE(AcousticSource_t)    , INTENT(IN)    :: sources
    TYPE(Config_Solver_t)     , INTENT(IN)    :: config
    TYPE(Config_Environment_t), INTENT(IN)    :: env

    ! Delete this in future
    INTEGER(I32), INTENT(IN), OPTIONAL        :: N_images
    INTEGER(I32), INTENT(IN), OPTIONAL        :: Upper_BC, Lower_BC
    REAL(WP)    , INTENT(IN), OPTIONAL        :: c_wat, rho_wat
    REAL(WP)    , INTENT(IN), OPTIONAL        :: Upper_HBC, Lower_HBC
    REAL(WP)    , INTENT(IN), OPTIONAL        :: attenuation_upper, attenuation_lower
    REAL(WP)    , INTENT(IN), OPTIONAL        :: eps, p_ref
    LOGICAL     , INTENT(IN), OPTIONAL        :: cluster, debug
    CHARACTER(len=*), INTENT(IN)              :: name

    ! Local variables
    CHARACTER(len=512) :: name_

    ! Defaults
    if (PRESENT(N_images))          self % N_images  = N_images
    if (PRESENT(c_wat))             self % c         = c_wat
    if (PRESENT(rho_wat))           self % rho       = rho_wat
    if (PRESENT(Upper_BC))          self % Upper_BC  = Upper_BC
    if (PRESENT(Lower_BC))          self % Lower_BC  = Lower_BC
    if (PRESENT(Upper_HBC))         self % Upper_HBC = Upper_HBC
    if (PRESENT(Lower_HBC))         self % Lower_HBC = Lower_HBC
    if (PRESENT(attenuation_upper)) self % up_R      = attenuation_upper
    if (PRESENT(attenuation_lower)) self % lw_R      = attenuation_lower
    if (PRESENT(eps))               self % eps       = eps
    if (PRESENT(p_ref))             self % p_ref     = p_ref
    if (PRESENT(cluster))           self % cluster   = cluster
    if (PRESENT(debug))             self % debug     = debug

    if (.not. PRESENT(Lower_HBC)) self % Lower_HBC = - env % depth_base
    
    ! Build save_path
    self % save_path = trim(self % turbines(1) % save_path) // trim(name) // ".hdf5"
    self % default_N_images = self % N_images
    
    if (self % Upper_HBC <= self % Lower_HBC) error stop "MethodImages % init: Upper_HBC must be strictly greater than Lower_HBC"
    
    print*, ''
    print*, 'Selected method: ', self % get_name()
    CALL self % build_all_image_system()
    CALL self % save_parameters()

END SUBROUTINE init_MethodImages


FUNCTION get_name_MethodImages(self) RESULT(name)
    CLASS(MethodImages_t), INTENT(INOUT) :: self
    CHARACTER(len=:)     , ALLOCATABLE   :: name
    name = "Images Method"; self % solver_name = trim(name)
END FUNCTION get_name_MethodImages


SUBROUTINE compute_pressure_MethodImages(self, observers, block_size, total_pressure)
    CLASS(MethodImages_t)   , INTENT(INOUT)        :: self
    REAL(WP)                , INTENT(IN)           :: observers(:,:)        ! [m] Observers coords (Nobs,3)
    INTEGER(I32)            , INTENT(IN), OPTIONAL :: block_size            ! [-] Number of observers to process simultaneously
    COMPLEX(WP), ALLOCATABLE, INTENT(OUT)          :: total_pressure(:,:)   ! [Hz] pressure array(Nfreqs, Nobs)

    ! Local variables
    INTEGER(I32) :: nobs, nf, nturb, nt, bs, n_chunks, ichunk
    INTEGER(I32) :: istart, iend, nb, i
    TYPE(ImageSystem_t) :: img
    REAL(WP), ALLOCATABLE :: freqs(:)
    COMPLEX(WP), ALLOCATABLE :: p_block(:,:)

    ! Defaults
    bs = 256_I32; if (PRESENT(block_size)) bs = block_size

    ! Basic checks
    nobs = SIZE(observers, 1)
    if (nobs <= 0) then
        allocate(total_pressure(0,0))
        return
    end if

    if (SIZE(observers, 2) /= 3) error stop "compute_pressure_MethodImages: observers must have shape (N,3)"

    ! Number of frequencies (take from first turbine)
    if (.not. ALLOCATED(self % turbines)) then
        error stop "compute_pressure_MethodImages: no turbines defined"
    end if

    nturb = SIZE(self % turbines)
    freqs = self % turbines(1) % Freqs
    nf = SIZE(freqs)

    allocate(total_pressure(nf, nobs))
    total_pressure = (0.0_WP, 0.0_WP)

    if (.not. self % cluster) bs = 1

    if (self % cluster) then
        do nt = 1, nturb
            img = self % image_system(nt)

            ! Print turbine header if more than one turbine
            if (nturb > 1) then
                print '(A,I0,A,I0,A)', 'Turbine ', nt, '/', nturb, ':'
                flush(6)
            end if

            n_chunks = (nobs + bs - 1) / bs

            !$OMP PARALLEL DO DEFAULT(SHARED) &
            !$OMP PRIVATE(ichunk, istart, iend, nb, p_block)&
            !$OMP SHARED(img, freqs, n_chunks, nobs, total_pressure, self)
            do ichunk = 0, n_chunks - 1
                istart = ichunk * bs + 1
                iend   = min(istart + bs - 1, nobs)
                nb     = iend - istart + 1

                ! Call kernel for this chunk (p_block will be allocated by the routine)
                call self % dipole_pressure_images(observers(istart:iend, :), freqs, img % nodes_pos, &
                                                   img % force, img % BC_all, p_block)

                ! Accumulate results into the global array (non-overlapping slices -> no race)
                total_pressure(:, istart:iend) = total_pressure(:, istart:iend) + p_block(:, 1:nb)

                if (ALLOCATED(p_block)) DEALLOCATE(p_block)
            end do
            !$OMP END PARALLEL DO
        end do
    else
        do nt = 1, nturb
            img = self % image_system(nt)

            ! Print turbine header if more than one turbine
            if (nturb > 1) then
                print '(A,I0,A,I0,A)', 'Turbine ', nt, '/', nturb, ':'
            end if

            ! Serial: loop observers (keep 2D slice for compatibility)
            do i = 1, nobs
                call self % dipole_pressure_images(observers(i:i, :), freqs, img % nodes_pos, &
                                                img % force, img % BC_all, p_block)
                total_pressure(:, i) = total_pressure(:, i) + p_block(:, 1)

                if (ALLOCATED(p_block)) DEALLOCATE(p_block)

                ! Print progress every 100 observers
                if (mod(i, 100) == 0) then
                    print '(A,I0,A,I0)', '  Progress: ', i, '/', nobs
                    flush(6)
                end if
            end do
        end do
    end if

END SUBROUTINE compute_pressure_MethodImages


SUBROUTINE save_parameters_MethodImages(self)
    CLASS(MethodImages_t), INTENT(IN) :: self

    ! Local variables
    CHARACTER(len=512) :: save_path

    save_path = trim(self % save_path)

    ! Solver
    CALL save_to_hdf5(save_path, "N_images" , self % N_images)
    CALL save_to_hdf5(save_path, "c_wat"    , self % c)
    CALL save_to_hdf5(save_path, "rho_wat"  , self % rho)
    CALL save_to_hdf5(save_path, "Upper_HBC", self % Upper_HBC)
    CALL save_to_hdf5(save_path, "Lower_HBC", self % Lower_HBC)
    CALL save_to_hdf5(save_path, "p_ref"    , self % p_ref)
    CALL save_to_hdf5(save_path, "Nturb"    , size(self % turbines))
    
    CALL save_turbine_params(self)

END SUBROUTINE save_parameters_MethodImages


! ---------- HELPERS ---------- !
SUBROUTINE build_all_image_system(self)

    CLASS(MethodImages_t), INTENT(INOUT) :: self
    INTEGER(I32) :: i, num_turbines
    COMPLEX(WP), ALLOCATABLE :: F_corr(:,:,:)

    if (.not. ALLOCATED(self % turbines)) return
    num_turbines = size(self % turbines)
    
    if (ALLOCATED(self % image_system)) deallocate(self % image_system)
    allocate(self % image_system(num_turbines))

    
    do i = 1, num_turbines
        F_corr = self % turbines(i) % get_impedance_corrected_force(self % c)
        CALL self % build_image_system(self % turbines(i) % x, F_corr, self % image_system(i))
    end do

END SUBROUTINE build_all_image_system


SUBROUTINE build_image_system(self, nodes_pos_real, force, img_sys)
    CLASS(MethodImages_t), INTENT(IN) :: self
    REAL(WP), INTENT(IN)              :: nodes_pos_real(:,:)   ! (Nnodes, 3)
    COMPLEX(WP), INTENT(IN)           :: force(:,:,:)          ! (Nfreqs, Nnodes, 3) or (Nfreqs, 3, Nnodes)
    TYPE(ImageSystem_t), INTENT(OUT)  :: img_sys

    REAL(WP), ALLOCATABLE :: zi(:), z_all(:)
    COMPLEX(WP), ALLOCATABLE :: force_prep(:,:,:), Force_out(:,:,:)
    INTEGER(I32), ALLOCATABLE :: parent(:), BC_all(:)
    INTEGER(I32) :: Nnodes, Nfreqs, total_nodes, j

    Nnodes = size(nodes_pos_real, 1)

    
    ! Ensure shape standard: (Nfreqs, Nnodes, 3)
    if (size(force, 3) == 3 .and. size(force, 2) == Nnodes) then
        force_prep = force
    else
        Nfreqs = size(force, 1)
        allocate(force_prep(Nfreqs, Nnodes, 3))
        do j = 1, Nnodes
            force_prep(:, j, 1) = force(:, 1, j)
            force_prep(:, j, 2) = force(:, 2, j)
            force_prep(:, j, 3) = force(:, 3, j)
        end do
    end if

    zi = nodes_pos_real(:, 3)
    
    CALL self % method_of_images(zi, force_prep, z_all, Force_out, parent, BC_all)

    total_nodes = size(z_all)
    allocate(img_sys % nodes_pos(total_nodes, 3))
    allocate(img_sys % force(size(Force_out, 1), total_nodes, 3))
    allocate(img_sys % BC_all(total_nodes))

    do j = 1, total_nodes
        img_sys % nodes_pos(j, 1) = nodes_pos_real(parent(j), 1)
        img_sys % nodes_pos(j, 2) = nodes_pos_real(parent(j), 2)
        img_sys % nodes_pos(j, 3) = z_all(j)
    end do
    
    img_sys % force  = Force_out
    img_sys % BC_all = BC_all

END SUBROUTINE build_image_system


SUBROUTINE set_N_images(self, new_N_images)
    CLASS(MethodImages_t), INTENT(INOUT) :: self
    INTEGER(I32), INTENT(IN)             :: new_N_images

    self % N_images = new_N_images
    CALL self % build_all_image_system()
END SUBROUTINE set_N_images


SUBROUTINE restore_default_images(self)
    CLASS(MethodImages_t), INTENT(INOUT) :: self

    if (self%N_images /= self%default_N_images) then
        self%N_images = self%default_N_images
        CALL self%build_all_image_system()
    end if
END SUBROUTINE restore_default_images


! ---------- METHOD ---------- !
SUBROUTINE method_of_images(self, zi, Force, z_all, Force_out, parent, BC_all)
    CLASS(MethodImages_t), INTENT(IN)         :: self
    REAL(WP), INTENT(IN)                      :: zi(:)
    COMPLEX(WP), INTENT(IN)                   :: Force(:,:,:)
    REAL(WP), ALLOCATABLE, INTENT(OUT)        :: z_all(:)
    COMPLEX(WP), ALLOCATABLE, INTENT(OUT)     :: Force_out(:,:,:)
    INTEGER(I32), ALLOCATABLE, INTENT(OUT)    :: parent(:)
    INTEGER(I32), ALLOCATABLE, INTENT(OUT)    :: BC_all(:)

    ! Local variables
    INTEGER(I32) :: Nnodes, Nfreqs, total_nodes, idx, i_node, img_lvl
    REAL(WP) :: z_upper, z_lower
    COMPLEX(WP), ALLOCATABLE :: F_upper(:,:), F_lower(:,:)
    CHARACTER(len=5) :: last_plane_upper, last_plane_lower
    INTEGER(I32) :: BC_u, BC_l
    
    Nnodes = size(zi); Nfreqs = size(Force, 1)
    total_nodes = Nnodes * (1 + 2 * self % N_images)
    
    ALLOCATE(z_all(total_nodes))
    ALLOCATE(Force_out(Nfreqs, total_nodes, 3))
    ALLOCATE(parent(total_nodes))
    ALLOCATE(BC_all(total_nodes))
    ALLOCATE(F_upper(Nfreqs, 3))
    ALLOCATE(F_lower(Nfreqs, 3))

    idx = 1
    do i_node = 1, Nnodes
        ! Real dipole nodes
        z_all(idx)         = zi(i_node)
        Force_out(:,idx,:) = Force(:,i_node,:)
        BC_all(idx)        = 1
        parent(idx)        = i_node
        idx                = idx + 1

        z_upper = zi(i_node)
        z_lower = zi(i_node)
        F_upper = Force(:,i_node,:)
        F_lower = Force(:,i_node,:)

        last_plane_upper = "upper"
        last_plane_lower = "lower"

        do img_lvl = 1, self % N_images
            ! ---------- Upper Reflection Chain ---------- !
            if (last_plane_upper == "upper") then
                z_upper = 2.0_WP * self % Lower_HBC - z_upper
                BC_u = self % Lower_BC
                last_plane_upper = "lower"
            else
                z_upper = 2.0_WP * self % Upper_HBC - z_upper
                BC_u = self % Upper_BC
                last_plane_upper = "upper"
            end if

            F_upper(:,3)       = - 1.0_WP * F_upper(:,3)
            F_upper            = F_upper * self % up_R

            z_all(idx)         = z_upper
            Force_out(:,idx,:) = F_upper
            BC_all(idx)        = BC_u
            parent(idx)        = i_node
            idx                = idx + 1

            ! ---------- Lower Reflection Chain ---------- !
            if (last_plane_lower == "lower") then
                z_lower = 2.0_WP * self % Upper_HBC - z_lower
                BC_l    = self % Upper_BC
                last_plane_lower = "upper"
            else 
                z_lower = 2.0_WP * self % Lower_HBC - z_lower
                BC_l    = self % Lower_BC
                last_plane_lower = "lower"
            end if

            F_lower(:,3) = -1.0_WP * F_lower(:,3)
            F_lower      = F_lower * self % lw_R

            z_all(idx)          = z_lower
            Force_out(:,idx, :) = F_lower
            BC_all(idx)         = BC_l
            parent(idx)         = i_node
            idx                 = idx + 1
        end do
    end do


END SUBROUTINE method_of_images


SUBROUTINE dipole_pressure_images(self, obs_pos, freqs, nodes_pos, force, BC_all, p_out)
    CLASS(MethodImages_t), INTENT(IN)       :: self
    REAL(WP), INTENT(IN)                    :: obs_pos(:,:)     ! (Nob, 3)  ; caller must pass a 2D slice (Nob,3)
    REAL(WP), INTENT(IN)                    :: freqs(:)         ! (Nfreqs)
    REAL(WP), INTENT(IN)                    :: nodes_pos(:,:)   ! (Ntotal,3)
    COMPLEX(WP), INTENT(IN)                 :: force(:,:,:)     ! (Nfreqs, Ntotal, 3)
    INTEGER(I32), INTENT(IN)                :: BC_all(:)        ! (Ntotal)
    COMPLEX(WP), ALLOCATABLE, INTENT(OUT)   :: p_out(:,:)       ! (Nfreqs, Nob)

    ! Local variables
    INTEGER(I32) :: Nfreqs, Ntotal, Nob
    INTEGER(I32) :: j, f, k
    REAL(WP), ALLOCATABLE :: dx(:), dy(:), dz(:), r(:), r_inv(:), k_wave(:), BC_real(:)
    REAL(WP) :: sx, sy, sz, epsv
    COMPLEX(WP) :: F_dot_r, term1, green, fx, fy, fz
    REAL(WP) :: fourpi

    ! Dimensions
    Nfreqs = SIZE(freqs)
    Ntotal = SIZE(nodes_pos, 1)
    Nob    = SIZE(obs_pos, 1)

    epsv = self % eps
    fourpi = 4.0_WP * PI

    ALLOCATE(p_out(Nfreqs, Nob))
    p_out = (0.0_WP, 0.0_WP)
    k_wave = 2.0_WP * PI * freqs / self % c
    BC_real = REAL(BC_all, WP)

    ! Temporary arrays to avoid recomputing distances for each frequency
    ALLOCATE(dx(Nob), dy(Nob), dz(Nob), r(Nob), r_inv(Nob))

    do j = 1, Ntotal
        sx = nodes_pos(j, 1)
        sy = nodes_pos(j, 2)
        sz = nodes_pos(j, 3)

        ! Compute relative vectors and distances for this source to all observers (vectorizable)
        !$OMP SIMD
        do k = 1, Nob
            dx(k) = obs_pos(k, 1) - sx
            dy(k) = obs_pos(k, 2) - sy
            dz(k) = obs_pos(k, 3) - sz
            r(k)  = sqrt(dx(k)*dx(k) + dy(k)*dy(k) + dz(k)*dz(k))
            if (r(k) < epsv) r(k) = epsv
            r_inv(k) = 1.0_WP / r(k)
        end do

        ! Loop over frequencies and accumulate contributions (inner loop vectorized over observers)
        do f = 1, Nfreqs
            fx = force(f, j, 1)
            fy = force(f, j, 2)
            fz = force(f, j, 3)
            !$OMP SIMD
            do k = 1, Nob
                F_dot_r = (fx * dx(k) + fy * dy(k) + fz * dz(k)) * r_inv(k)
                term1   = -I1 * k_wave(f) + r_inv(k)
                green   = exp(I1* k_wave(f) * r(k)) / (fourpi * r(k))
                p_out(f, k) = p_out(f, k) + F_dot_r * term1 * green * BC_real(j)
            end do
        end do

    end do

    DEALLOCATE(dx, dy, dz, r, r_inv)

END SUBROUTINE dipole_pressure_images


! ---------- RUN SPHERE ---------- !
SUBROUTINE run_sphere(self, r, n_theta, nz, center, block_size)
    CLASS(MethodImages_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: r
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_theta
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nz
    REAL(WP)    , INTENT(IN), OPTIONAL  :: center(3)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    COMPLEX(WP), ALLOCATABLE :: p_reshaped(:,:,:)   ! (Nfreqs, n_theta, nz)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)   , ALLOCATABLE :: obs_reshaped(:,:,:) ! (n_theta, nz, 3)
    REAL(WP)   , ALLOCATABLE :: theta_deg_arr(:)    ! Azimuth angles in degrees
    REAL(WP)   , ALLOCATABLE :: z_centers(:)        ! Vertical coordinates of layers
    REAL(WP)   , ALLOCATABLE :: z_edges(:), theta_edges(:)
    REAL(WP)                 :: r_, center_(3), max_dist, current_dist
    REAL(WP)                 :: dz, dtheta, zc, r_xy, theta_c
    REAL(WP)                 :: dA
    INTEGER(I32)             :: n_theta_, nz_, block_size_, i, j, idx, Nobs

    ! Defaults
    r_ = 30.0_WP;          if (PRESENT(r))          r_ = r
    n_theta_ = 72_I32;     if (PRESENT(n_theta))    n_theta_ = n_theta
    nz_ = 20_I32;          if (PRESENT(nz))         nz_ = nz
    block_size_ = 128_I32; if (PRESENT(block_size)) block_size_ = block_size

    if (PRESENT(center)) then
        center_ = center
    else
        center_ = [self % turbines(1) % BariPos(1), self % turbines(1) % BariPos(2), -self % turbines(1) % Depth / 2.0_WP]
    end if

    ! Basic enclosure check
    max_dist = 0.0_WP
    do i = 1, SIZE(self % turbines(1) % x, 1)
        current_dist = NORM2(self % turbines(1) % x(i,:) - center_)
        if (current_dist > max_dist) max_dist = current_dist
    end do
    max_dist = max_dist * 1.1_WP

    if (r_ < max_dist) then
        print*, "WindTurbine.run_sphere(): requested radius r=", r_, "m is smaller than max dist."
        print*, "Radius increased to ", max_dist * 1.05_WP, "m"
        r_ = max_dist * 1.05_WP
    end if

    Nobs = n_theta_ * nz_
    ALLOCATE(observers_(Nobs, 3))
    ALLOCATE(theta_deg_arr(n_theta_))
    ALLOCATE(z_centers(nz_), z_edges(nz_+1), theta_edges(n_theta_+1))

    ! Compute edges and centers for vertical and azimuthal discretization
    ! Vertical layers: centers between edges
    dz = (2.0_WP * r_) / REAL(nz_, WP)
    do i = 1, nz_ + 1
        z_edges(i) = (center_(3) - r_) + REAL(i - 1, WP) * dz
    end do
    do i = 1, nz_
        z_centers(i) = 0.5_WP * (z_edges(i) + z_edges(i+1))
    end do

    ! Azimuthal: centers between edges
    dtheta = (2.0_WP * PI) / REAL(n_theta_, WP)
    do j = 1, n_theta_ + 1
        theta_edges(j) = REAL(j - 1, WP) * dtheta
    end do
    do j = 1, n_theta_
        theta_deg_arr(j) = 0.5_WP * (theta_edges(j) + theta_edges(j+1)) * 180.0_WP / PI
    end do

    ! Build observers array
    idx = 1
    do i = 1, nz_
        zc = z_centers(i)
        r_xy = SQRT(MAX(0.0_WP, r_**2 - (zc - center_(3))**2))
        do j = 1, n_theta_
            theta_c = 0.5_WP * (theta_edges(j) + theta_edges(j+1))
            observers_(idx, 1) = center_(1) + r_xy * COS(theta_c)
            observers_(idx, 2) = center_(2) + r_xy * SIN(theta_c)
            observers_(idx, 3) = zc
            idx = idx + 1
        end do
    end do

    CALL self % check_observers_distances(observers_, 0.1_WP)

    write(*, '(A)') ''
    write(*, '(A)') 'Computing spherical pressure field ...'
    write(*, '(A, F0.2, A, F0.2, A, F0.2, A)') &
        '  Centre: (', center_(1), ', ', center_(2), ', ', center_(3), ') m'
    write(*, '(A, F0.2, A, I0, A, I0, A, I0, A)') &
        '  Radius: ', r_, ' m | grid: ', n_theta_, ' azim x ', nz_, ' z (', Nobs, ' observers)'

    CALL self % set_N_images(0_I32)
    CALL self % compute_pressure(observers_, block_size_, p)
    CALL self % restore_default_images()

    if (self % debug) write(*, '(A, F0.2)') 'Sphere Pressure Norm: ', norm2(abs(p))

    ! Reshape pressure and observers to 3D (n_theta, nz, 3) and (Nfreqs, n_theta, nz)
    ALLOCATE(p_reshaped(SIZE(p,1), n_theta_, nz_))
    ALLOCATE(obs_reshaped(n_theta_, nz_, 3))

    idx = 1
    do i = 1, nz_
        do j = 1, n_theta_
            obs_reshaped(j, i, :) = observers_(idx, :)
            p_reshaped(:, j, i)   = p(:, idx)
            idx = idx + 1
        end do
    end do

    ! Compute surface element (approximate)
    dA = (4.0_WP * PI * r_**2) / REAL(n_theta_ * nz_, WP)

    ! Save results
    CALL self % save_acoustics("Freqs", self % turbines(1) % Freqs)
    CALL self % save_acoustics("P_sphere", p_reshaped)
    CALL self % save_acoustics("Obs_sphere", obs_reshaped)
    CALL self % save_acoustics("Theta_sphere", theta_deg_arr)
    CALL self % save_acoustics("Z_sphere", z_centers)
    CALL self % save_acoustics("R_sphere", r_)
    CALL self % save_acoustics("Center_sphere", center_)
    CALL self % save_acoustics("dA_sphere", dA)

    print*, ''
    print*, '   ---> Sphere data saved at ', self % save_path

END SUBROUTINE run_sphere











END MODULE ImageMethod