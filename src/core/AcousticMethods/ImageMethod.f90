! ===============================================================================
! MODULE: ImageMethod
! ===============================================================================
! Defines ImageMethod_t, the concrete acoustic solver extending abstract
! AcousticSolver_t. It implements the Method of Images for a fluid layer
! bounded by two horizontal reflecting planes (e.g. free surface above,
! seabed below): each real dipole source is mirrored into an infinite chain
! of image dipoles above and below, with reflection coefficients and boundary
! conditions applied per bounce, then the pressure field is computed as the
! superposition of all real+image dipole contributions at each observer point
! (with optional multi-observer chunked/parallel evaluation).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 2026-09-22
! ===============================================================================
MODULE ImageMethod

USE Kinds
USE Types
USE IOUtils
USE WindTurbine, ONLY: WindTurbine_t, save_turbine_params
USE AcousticSolver, ONLY: AcousticSolver_t
IMPLICIT NONE

PRIVATE
PUBLIC:: ImageMethod_t


! ===============================================================================
! TYPE: ImageMethod_t
! ===============================================================================
! Concrete acoustic solver implementing the Method of Images: real dipole
! sources are reflected across an upper and a lower horizontal plane to build
! an equivalent free-space system (image_system) whose superposed dipole
! fields reproduce the effect of the bounded fluid layer on the pressure.
! ===============================================================================
TYPE, EXTENDS(AcousticSolver_t):: ImageMethod_t

    INTEGER(I32) :: N_images         = 30           ! [-] Number of image levels to compute
    INTEGER(I32) :: default_N_images = 30           ! [-] Original N_images from config, used to restore after temporary changes
    INTEGER(I32) :: Upper_BC         = -1           ! [-] Upper Boundary condition -1 for  p = 0
    INTEGER(I32) :: Lower_BC         = +1           ! [-] Lower Boundary condition +1 for dp = 0
    REAL(WP)     :: Upper_HBC        = 0.0_WP       ! [m] z-coordinate of the upper reflecting plane (e.g. free surface)
    REAL(WP)     :: Lower_HBC        = -30.0_WP     ! [m] z-coordinate of the lower reflecting plane (e.g. seabed)
    REAL(WP)     :: up_R             = 1.0_WP       ! [-] Attenuation coefficient for upper boundary reflections
    REAL(WP)     :: lw_R             = 0.5_WP       ! [-] Attenuation coefficient for lower boundary reflections
    REAL(WP)     :: eps              = 1e-6_WP      ! [m] Minimum source-observer distance, floors r to avoid singularities

    TYPE(ImageSystem_t) :: image_system             ! [-] Full real+image dipole system (coords, forces, boundary tags)

    CONTAINS
    ! --- Public type-bound procedures --- !
    PROCEDURE :: set_N_images
    PROCEDURE :: build_image_system
    PROCEDURE :: restore_default_images
    PROCEDURE :: method_of_images
    PROCEDURE :: dipole_pressure_images
    PROCEDURE :: run_sphere

    ! --- Deferred(abstract) procedures --- !
    PROCEDURE :: init               => init_ImageMethod
    PROCEDURE :: get_name           => get_name_ImageMethod
    PROCEDURE :: compute_pressure   => compute_pressure_ImageMethod
    PROCEDURE :: save_parameters    => save_parameters_ImageMethod
    PROCEDURE :: build_save_path    => build_save_path_ImageMethod
    PROCEDURE :: write_txt_summary  => write_txt_summary_ImageMethod

END TYPE ImageMethod_t


CONTAINS


! ===============================================================================
! DEFERRED (abstract AcousticSolver_t) PROCEDURES
! ===============================================================================
! -------------------------------------------------------------------------
! Subroutine: init_ImageMethod
! Constructor of the Image Method solver. Copies the shared solver/env/
! turbine data from the parent AcousticSolver_t setup, reads the method-
! specific parameters (N_images, boundary planes, reflection coefficients,
! singularity tolerance) from the configuration, resolves the output save
! path and builds the initial image system.
! -------------------------------------------------------------------------
SUBROUTINE init_ImageMethod(self, sources, config)

    CLASS(ImageMethod_t)  , INTENT(INOUT)              :: self
    TYPE(AcousticSource_t), INTENT(INOUT), ALLOCATABLE :: sources
    TYPE(Config_t)        , INTENT(IN)                 :: config
    CHARACTER(len=:), ALLOCATABLE :: name

    ! Parent variables
    self % sources      = sources
    self % save_path    = config % solver % save_path
    self % env          = config % env
    self % block_size   = config % solver % block_size
    self % debug        = config % solver % debug
    self % in_cluster   = config % solver % in_cluster
    self % turbine_data = config % turbine

    DEALLOCATE(sources)

    ! Son specific variables
    SELECT TYPE(m => config % method)
    TYPE IS(Config_ImageMethod_t)
        self % N_images         = m % N_images
        self % default_N_images = m % N_images
        self % Upper_BC         = m % UpperBC
        self % Lower_BC         = m % LowerBC
        self % Upper_HBC        = m % UpperHBC
        self % Lower_HBC        = m % LowerHBC
        self % up_R             = m % attenuation_surface
        self % lw_R             = m % attenuation_bottom
        self % eps              = m % epsilon
    END SELECT
    
    ! Build save_path
    CALL self % build_save_path()
    
    if (self % Upper_HBC <= self % Lower_HBC) error stop "ImageMethod % init: Upper_HBC must be strictly greater than Lower_HBC"
    
    print*, ''
    print*, 'Selected method: ', self % get_name()
    CALL self % build_image_system()

END SUBROUTINE init_ImageMethod


! -------------------------------------------------------------------------
! Function: get_name_ImageMethod
! Returns the human-readable solver name ("Images Method") and stores it
! in self % solver_name for later use (logs, summaries, HDF5 output).
! -------------------------------------------------------------------------
FUNCTION get_name_ImageMethod(self) RESULT(name)
    CLASS(ImageMethod_t), INTENT(INOUT) :: self
    CHARACTER(len=:)     , ALLOCATABLE   :: name
    name = "Images Method"; self % solver_name = trim(name)
END FUNCTION get_name_ImageMethod


! -------------------------------------------------------------------------
! Subroutine: compute_pressure_ImageMethod
! Evaluates the total radiated pressure field at an arbitrary set of
! observer points, for every source frequency. Observers are processed
! either one-by-one (serial) or split into chunks distributed over OpenMP
! threads (in_cluster = .true.), each chunk being handed to
! dipole_pressure_images together with the full real+image dipole system.
! -------------------------------------------------------------------------
SUBROUTINE compute_pressure_ImageMethod(self, observers, block_size, total_pressure)
    CLASS(ImageMethod_t)    , INTENT(INOUT)        :: self
    REAL(WP)                , INTENT(IN)           :: observers(:,:)        ! [m] Observers coords (Nobs,3)
    INTEGER(I32)            , INTENT(IN), OPTIONAL :: block_size            ! [-] Number of observers to process simultaneously
    COMPLEX(WP), ALLOCATABLE, INTENT(OUT)          :: total_pressure(:,:)   ! [Hz] pressure array(NFreqs, Nobs)

    ! Local variables
    INTEGER(I32)             :: nobs, nf, bs, n_chunks, ichunk
    INTEGER(I32)             :: istart, iend, nb, i
    TYPE(ImageSystem_t)      :: img
    REAL(WP), ALLOCATABLE    :: Freqs(:)
    COMPLEX(WP), ALLOCATABLE :: p_block(:,:)

    ! Defaults
    bs = 256_I32; if (PRESENT(block_size)) bs = block_size

    ! Basic checks
    nobs = SIZE(observers, 1)
    if (nobs <= 0) then
        allocate(total_pressure(0,0))
        return
    end if

    if (SIZE(observers, 2) /= 3) error stop "compute_pressure_ImageMethod: observers must have shape (N,3)"

    ! Number of frequencies (take from first turbine)
    if (.not. ALLOCATED(self % sources % xyz) .OR. .not. ALLOCATED(self % sources % Force) .OR. &
        .not. ALLOCATED(self % sources % Freqs)) then
        error stop "compute_pressure_ImageMethod: no sources defined"
    end if

    Freqs = self % sources % Freqs
    nf = SIZE(Freqs)

    allocate(total_pressure(nf, nobs))
    total_pressure = (0.0_WP, 0.0_WP)

    if (.not. self % in_cluster) bs = 1
    img = self % image_system

    if (self % in_cluster) then
            
        n_chunks = (nobs + bs - 1) / bs

        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP PRIVATE(ichunk, istart, iend, nb, p_block)&
        !$OMP SHARED(img, Freqs, n_chunks, nobs, total_pressure, self)
        do ichunk = 0, n_chunks - 1
            istart = ichunk * bs + 1
            iend   = min(istart + bs - 1, nobs)
            nb     = iend - istart + 1

            ! Call kernel for this chunk (p_block will be allocated by the routine)
            call self % dipole_pressure_images(observers(istart:iend, :), Freqs, img % xyz, &
                                                img % Force, img % BC_all, p_block)

            ! Accumulate results into the global array (non-overlapping slices -> no race)
            total_pressure(:, istart:iend) = total_pressure(:, istart:iend) + p_block(:, 1:nb)

            if (ALLOCATED(p_block)) DEALLOCATE(p_block)
        end do
        !$OMP END PARALLEL DO

    else !(.not. self % in_cluster)

        ! Serial: loop observers (keep 2D slice for compatibility)
        do i = 1, nobs
            call self % dipole_pressure_images(observers(i:i, :), Freqs, img % xyz, &
                                            img % Force, img % BC_all, p_block)
            total_pressure(:, i) = total_pressure(:, i) + p_block(:, 1)

            if (ALLOCATED(p_block)) DEALLOCATE(p_block)

            ! Print progress every 100 observers
            if (mod(i, 100) == 0) then
                print '(A,I0,A,I0)', '  Progress: ', i, '/', nobs
                flush(6)
            end if
        end do

    end if

END SUBROUTINE compute_pressure_ImageMethod


! -------------------------------------------------------------------------
! Subroutine: save_parameters_ImageMethod
! Saves the solver, environment and generic turbine parameters to the
! HDF5 output file, delegates the turbine-specific parameters to
! turbine % save_specific_params, and triggers the human-readable .dat
! text summary via write_txt_summary.
! -------------------------------------------------------------------------
SUBROUTINE save_parameters_ImageMethod(self, turbine)

    CLASS(ImageMethod_t), INTENT(IN)                        :: self
    CLASS(WindTurbine_t), INTENT(IN), OPTIONAL, ALLOCATABLE :: turbine

    ! Solver / environment / generic turbine data -> HDF5 (unchanged)
    CALL save_variable(self % save_path, "N_images" , self % N_images)
    CALL save_variable(self % save_path, "c_wat"    , self % env % speed_sound_base)
    CALL save_variable(self % save_path, "rho_wat"  , self % env % density_base)
    CALL save_variable(self % save_path, "Upper_HBC", self % Upper_HBC)
    CALL save_variable(self % save_path, "Lower_HBC", self % Lower_HBC)
    CALL save_variable(self % save_path, "p_ref"    , self % env % reference_pressure)

    CALL save_turbine_params(self % turbine_data, self % save_path)

    if (PRESENT(turbine)) THEN
        if (ALLOCATED(turbine)) THEN
            CALL turbine % save_specific_params(self % save_path)
        end if
    end if

    ! Human-readable .dat summary, delegated to the solver-specific implementation
    CALL self % write_txt_summary(turbine)

END SUBROUTINE save_parameters_ImageMethod


! -------------------------------------------------------------------------
! Subroutine: build_save_path_ImageMethod
! Resolves the final HDF5 output path. If the configured save_path is
! "None"/"Default", it auto-generates a descriptive filename (turbine
! type, number of images, wind speed) under a default output folder;
! otherwise it just appends the .hdf5 extension. Ensures the destination
! folder exists.
! -------------------------------------------------------------------------
SUBROUTINE build_save_path_ImageMethod(self)
    CLASS(ImageMethod_t), INTENT(INOUT) :: self
    
    ! Local variables
    CHARACTER(len=512) :: path_upper, filename, default_parent_path
    CHARACTER(len=16)  :: method_suffix, ws_str, imag_str
    CHARACTER(len=2)   :: turb_tag

    path_upper = to_upper(trim(self % save_path))
    path_upper = trim(path_upper)
    default_parent_path = "./turbine_acoustic_data/"

    method_suffix = "_IM"

    if (path_upper == 'NONE' .or. path_upper == 'DEFAULT') then

        if (self % turbine_data % is_monopile) then
            turb_tag = 'mn'
        else
            turb_tag = 'fl'
        end if

        write(imag_str,'(I0)') self % N_images
        write(ws_str, '(F6.1)') self % turbine_data % WindSpeed
        ws_str = adjustl(ws_str)

        filename = "plot_"             // trim(turb_tag) // &
                   trim(method_suffix) // trim(imag_str) // &
                   "_ws"               // trim(ws_str)

        self % save_path = trim(default_parent_path) // trim(filename) // ".hdf5"

    else
        self % save_path = trim(self % save_path) // ".hdf5"

    end if

    ! Make sure the destination folder exists before anything tries to write into it
    CALL ensure_parent_dir(trim(self % save_path))

END SUBROUTINE build_save_path_ImageMethod


! -------------------------------------------------------------------------
! Subroutine: write_txt_summary_ImageMethod
! Writes a plain-text (.dat) summary next to the HDF5 output, mirroring
! the saved turbine, environment and Image Method solver parameters in a
! human-readable form for quick inspection without opening the HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE write_txt_summary_ImageMethod(self, turbine)
    CLASS(ImageMethod_t), INTENT(IN)                        :: self
    CLASS(WindTurbine_t), INTENT(IN), OPTIONAL, ALLOCATABLE :: turbine

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: dat_path
    INTEGER(I32) :: u

    dat_path = swap_extension(trim(self % save_path), ".dat")

    OPEN(NEWUNIT=u, FILE=dat_path, STATUS='replace', ACTION='write')

    WRITE(u,'(A)') repeat('=', 70)
    WRITE(u,'(A)') 'OFF_COUSTICS SIMULATION SUMMARY'
    WRITE(u,'(A)') repeat('=', 70)
    WRITE(u,'(A)') 'HDF5 file: ' // trim(self % save_path)
    WRITE(u,*)

    WRITE(u,'(A)') repeat('-', 70)
    WRITE(u,'(A)') 'TURBINE'
    WRITE(u,'(A)') repeat('-', 70)
    WRITE(u,'(A,A)')      'Case type   : ', trim(self % turbine_data % case_type)
    WRITE(u,'(A,L1)')     'Is monopile : ', self % turbine_data % is_monopile
    WRITE(u,'(A,F8.2,A)') 'Wind speed  : ', self % turbine_data % WindSpeed, ' m/s'
    WRITE(u,'(A,F8.2,A)') 'Wind dir    : ', self % turbine_data % WindDir, ' deg'
    WRITE(u,'(A,F8.2,A)') 'Water depth : ', self % turbine_data % Depth, ' m'
    WRITE(u,'(A,2F10.3)') 'Axis pos    : ', self % turbine_data % AxisPos
    WRITE(u,'(A,2F10.3)') 'Bari pos    : ', self % turbine_data % BariPos
    WRITE(u,'(A,I0)')     'N members   : ', self % turbine_data % Nmembers
    WRITE(u,'(A,I0)')     'N nodes     : ', self % turbine_data % Nnodes
    WRITE(u,*)

    if (PRESENT(turbine)) then
        if (ALLOCATED(turbine)) then
            WRITE(u,'(A)') repeat('-', 70)
            WRITE(u,'(A)') 'TURBINE - SPECIFIC PARAMETERS'
            WRITE(u,'(A)') repeat('-', 70)
            CALL turbine % write_specific_params_txt(u)
            WRITE(u,*)
        end if
    end if

    WRITE(u,'(A)') repeat('-', 70)
    WRITE(u,'(A)') 'ENVIRONMENT'
    WRITE(u,'(A)') repeat('-', 70)
    WRITE(u,'(A,F8.2,A)')  'Speed of sound   : ', self % env % speed_sound_base, ' m/s'
    WRITE(u,'(A,F8.2,A)')  'Density          : ', self % env % density_base, ' kg/m^3'
    WRITE(u,'(A,ES10.2,A)') 'Reference press. : ', self % env % reference_pressure, ' Pa'
    WRITE(u,'(A,F8.2,A)')  'Depth            : ', self % env % depth_base, ' m'
    WRITE(u,*)

    WRITE(u,'(A)') repeat('-', 70)
    WRITE(u,'(A)') 'SOLVER: ' // trim(self % solver_name)
    WRITE(u,'(A)') repeat('-', 70)
    WRITE(u,'(A,I0)')     'N images     : ', self % N_images
    WRITE(u,'(A,F8.2,A)') 'Upper HBC    : ', self % Upper_HBC, ' m'
    WRITE(u,'(A,F8.2,A)') 'Lower HBC    : ', self % Lower_HBC, ' m'
    WRITE(u,'(A,F8.4)')   'Upper atten. : ', self % up_R
    WRITE(u,'(A,F8.4)')   'Lower atten. : ', self % lw_R
    WRITE(u,'(A,I0)')     'Block size   : ', self % block_size
    WRITE(u,'(A,L1)')     'In cluster   : ', self % in_cluster
    WRITE(u,*)
    WRITE(u,'(A)') repeat('=', 70)

    CLOSE(u)

    print*, '   ---> Text summary saved at ', dat_path

END SUBROUTINE write_txt_summary_ImageMethod


! ===============================================================================
! HELPERS
! ===============================================================================
! -------------------------------------------------------------------------
! Subroutine: build_image_system
! Builds (or rebuilds) self % image_system from the real sources: reshapes
! the per-source force array to the standard (NFreqs, Nnodes, 3) layout,
! calls method_of_images to generate the full chain of image z-positions,
! forces and boundary-condition tags, and assembles the resulting (x,y,z)
! coordinates by pairing each image with its parent source's horizontal
! position.
! -------------------------------------------------------------------------
SUBROUTINE build_image_system(self)
    CLASS(ImageMethod_t), INTENT(INOUT) :: self

    REAL(WP)    , ALLOCATABLE :: zi(:), z_all(:)
    COMPLEX(WP) , ALLOCATABLE :: Force(:,:,:), Force_prep(:,:,:), Force_out(:,:,:)
    INTEGER(I32), ALLOCATABLE :: parent(:), BC_all(:)
    INTEGER(I32)              :: Nnodes, NFreqs, total_nodes, j

    Nnodes = size(self % sources % xyz, 1)
    Force  = self % sources % Force 

    ! Ensure shape standard: (NFreqs, Nnodes, 3)
    if (size(Force, 3) == 3 .and. size(Force, 2) == Nnodes) then
        Force_prep = Force
    else
        NFreqs = size(self % sources % Freqs)
        allocate(Force_prep(NFreqs, Nnodes, 3))
        do j = 1, Nnodes
            Force_prep(:, j, 1) = Force(:, 1, j)
            Force_prep(:, j, 2) = Force(:, 2, j)
            Force_prep(:, j, 3) = Force(:, 3, j)
        end do
    end if

    zi = self % sources % xyz(:, 3)
    
    CALL self % method_of_images(zi, Force_prep, z_all, Force_out, parent, BC_all)

    total_nodes = size(z_all)
    ! Deallocate if rebuilding
    if (ALLOCATED(self % image_system % xyz)) DEALLOCATE(self % image_system % xyz)
    if (ALLOCATED(self % image_system % Force))     DEALLOCATE(self % image_system % Force)
    if (ALLOCATED(self % image_system % BC_all))    DEALLOCATE(self % image_system % BC_all)
    allocate(self % image_system % xyz(total_nodes, 3))
    allocate(self % image_system % Force(size(Force_out, 1), total_nodes, 3))
    allocate(self % image_system % BC_all(total_nodes))

    do j = 1, total_nodes
        self % image_system % xyz(j, 1) = self % sources % xyz(parent(j), 1)
        self % image_system % xyz(j, 2) = self % sources % xyz(parent(j), 2)
        self % image_system % xyz(j, 3) = z_all(j)
    end do
    
    self % image_system % Force  = Force_out
    self % image_system % BC_all = BC_all

END SUBROUTINE build_image_system


! -------------------------------------------------------------------------
! Subroutine: set_N_images
! Changes the number of image levels and rebuilds the image system to
! reflect it. Used e.g. to temporarily set N_images = 0 (free-field, no
! reflections) for calibration runs such as run_sphere.
! -------------------------------------------------------------------------
SUBROUTINE set_N_images(self, new_N_images)
    CLASS(ImageMethod_t), INTENT(INOUT) :: self
    INTEGER(I32), INTENT(IN)             :: new_N_images

    self % N_images = new_N_images
    CALL self % build_image_system()
END SUBROUTINE set_N_images


! -------------------------------------------------------------------------
! Subroutine: restore_default_images
! Restores N_images to its original (config-loaded) value and rebuilds
! the image system if it had been temporarily changed by set_N_images.
! -------------------------------------------------------------------------
SUBROUTINE restore_default_images(self)
    CLASS(ImageMethod_t), INTENT(INOUT) :: self

    if (self%N_images /= self%default_N_images) then
        self%N_images = self%default_N_images
        CALL self % build_image_system()
    end if
END SUBROUTINE restore_default_images


! ===============================================================================
! METHOD (core physics)
! ===============================================================================
! -------------------------------------------------------------------------
! Subroutine: method_of_images
! Builds the full chain of image dipoles for each real source by repeated
! reflection across the upper and lower horizontal planes, N_images times
! in each direction. Starting from the real source, the upper chain
! alternately reflects across the lower then upper plane (and vice versa
! for the lower chain), each bounce flipping the sign of the vertical
! force component, scaling the force by the corresponding reflection
! coefficient (up_R/lw_R) and tagging the image with the boundary
! condition (+-1) of the plane it was just reflected across. The result
! is a flat list of z-positions, forces and BC tags (real source + all
! images) later consumed by dipole_pressure_images.
! -------------------------------------------------------------------------
SUBROUTINE method_of_images(self, zi, Force, z_all, Force_out, parent, BC_all)
    CLASS(ImageMethod_t), INTENT(IN)         :: self
    REAL(WP), INTENT(IN)                      :: zi(:)
    COMPLEX(WP), INTENT(IN)                   :: Force(:,:,:)
    REAL(WP), ALLOCATABLE, INTENT(OUT)        :: z_all(:)
    COMPLEX(WP), ALLOCATABLE, INTENT(OUT)     :: Force_out(:,:,:)
    INTEGER(I32), ALLOCATABLE, INTENT(OUT)    :: parent(:)
    INTEGER(I32), ALLOCATABLE, INTENT(OUT)    :: BC_all(:)

    ! Local variables
    INTEGER(I32) :: Nnodes, NFreqs, total_nodes, idx, i_node, img_lvl
    REAL(WP) :: z_upper, z_lower
    COMPLEX(WP), ALLOCATABLE :: F_upper(:,:), F_lower(:,:)
    CHARACTER(len=5) :: last_plane_upper, last_plane_lower
    INTEGER(I32) :: BC_u, BC_l
    
    Nnodes = size(zi); NFreqs = size(Force, 1)
    total_nodes = Nnodes * (1 + 2 * self % N_images)
    
    ALLOCATE(z_all(total_nodes))
    ALLOCATE(Force_out(NFreqs, total_nodes, 3))
    ALLOCATE(parent(total_nodes))
    ALLOCATE(BC_all(total_nodes))
    ALLOCATE(F_upper(NFreqs, 3))
    ALLOCATE(F_lower(NFreqs, 3))

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


! -------------------------------------------------------------------------
! Subroutine: dipole_pressure_images
! Evaluates the acoustic dipole Green's function for a block of observers,
! summed over every real and image source, at every frequency. For each
! source-observer pair it computes the separation vector and distance
! (floored at eps to avoid the 1/r singularity), then the far/near-field
! dipole pressure term F.r_hat * (1/r - i*k) * exp(i*k*r) / (4*pi*r),
! weighted by the source's boundary-condition sign (BC_all), and
! accumulates it into p_out. Distances are precomputed once per source
! and reused across all frequencies.
! -------------------------------------------------------------------------
SUBROUTINE dipole_pressure_images(self, obs_pos, Freqs, nodes_pos, Force, BC_all, p_out)
    CLASS(ImageMethod_t), INTENT(IN)       :: self
    REAL(WP), INTENT(IN)                    :: obs_pos(:,:)     ! (Nob, 3)  ; caller must pass a 2D slice (Nob,3)
    REAL(WP), INTENT(IN)                    :: Freqs(:)         ! (NFreqs)
    REAL(WP), INTENT(IN)                    :: nodes_pos(:,:)   ! (Ntotal,3)
    COMPLEX(WP), INTENT(IN)                 :: Force(:,:,:)     ! (NFreqs, Ntotal, 3)
    INTEGER(I32), INTENT(IN)                :: BC_all(:)        ! (Ntotal)
    COMPLEX(WP), ALLOCATABLE, INTENT(OUT)   :: p_out(:,:)       ! (NFreqs, Nob)

    ! Local variables
    INTEGER(I32) :: NFreqs, Ntotal, Nob
    INTEGER(I32) :: j, f, k
    REAL(WP), ALLOCATABLE :: dx(:), dy(:), dz(:), r(:), r_inv(:), k_wave(:), BC_real(:)
    REAL(WP) :: sx, sy, sz, epsv
    COMPLEX(WP) :: F_dot_r, term1, green, fx, fy, fz
    REAL(WP) :: fourpi

    ! Dimensions
    NFreqs = SIZE(Freqs)
    Ntotal = SIZE(nodes_pos, 1)
    Nob    = SIZE(obs_pos, 1)

    epsv = self % eps
    fourpi = 4.0_WP * PI

    ALLOCATE(p_out(NFreqs, Nob))
    p_out = (0.0_WP, 0.0_WP)
    k_wave = 2.0_WP * PI * Freqs / self % env % speed_sound_base
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
        do f = 1, NFreqs
            fx = Force(f, j, 1)
            fy = Force(f, j, 2)
            fz = Force(f, j, 3)
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


! ===============================================================================
! RUN SPHERE
! ===============================================================================
! -------------------------------------------------------------------------
! Subroutine: run_sphere
! Builds a spherical (theta, z) grid of observer points enclosing the
! turbine, enlarging the requested radius if needed so every source lies
! inside it, then evaluates the free-field pressure on that surface
! (temporarily forcing N_images = 0, i.e. no boundary reflections) for
! later far-field/radiated-power post-processing. Reshapes and saves the
! pressure field, observer coordinates, angular/vertical grids and surface
! element area to the HDF5 output.
! -------------------------------------------------------------------------
SUBROUTINE run_sphere(self, r_sphere, n_theta_sphere, nz_sphere, center_sphere)
    CLASS(ImageMethod_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL :: r_sphere
    INTEGER(I32), INTENT(IN), OPTIONAL :: n_theta_sphere
    INTEGER(I32), INTENT(IN), OPTIONAL :: nz_sphere
    REAL(WP)    , INTENT(IN), OPTIONAL :: center_sphere(3)

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    COMPLEX(WP), ALLOCATABLE :: p_reshaped(:,:,:)   ! (NFreqs, n_theta, nz)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)   , ALLOCATABLE :: obs_reshaped(:,:,:) ! (n_theta, nz, 3)
    REAL(WP)   , ALLOCATABLE :: theta_deg_arr(:)    ! Azimuth angles in degrees
    REAL(WP)   , ALLOCATABLE :: z_centers(:)        ! Vertical coordinates of layers
    REAL(WP)   , ALLOCATABLE :: z_edges(:), theta_edges(:)
    REAL(WP)                 :: r_, center_(3), max_dist, current_dist
    REAL(WP)                 :: dz, dtheta, zc, r_xy, theta_c
    REAL(WP)                 :: dA
    INTEGER(I32)             :: n_theta_, nz_, i, j, idx, Nobs

    ! Defaults
    r_ = 30.0_WP;          if (PRESENT(r_sphere))          r_ = r_sphere
    n_theta_ = 72_I32;     if (PRESENT(n_theta_sphere))    n_theta_ = n_theta_sphere
    nz_ = 20_I32;          if (PRESENT(nz_sphere))         nz_ = nz_sphere
    if (PRESENT(center_sphere)) then
        center_ = center_sphere
    else
        center_ = [self % turbine_data % BariPos(1), self % turbine_data % BariPos(2), -self % turbine_data % Depth / 2.0_WP]
    end if

    ! Basic enclosure check
    max_dist = 0.0_WP
    do i = 1, SIZE(self % sources % xyz, 1)
        current_dist = norm2(self % sources % xyz(i,:) - center_)
        if (current_dist > max_dist) max_dist = current_dist
    end do
    max_dist = max_dist * 1.1_WP

    if (r_ < max_dist) then
        WRITE(*,'(A,F8.2,A)') '   requested radius r=', r_, ' m is smaller than max dist.'
        WRITE(*,'(A,F8.2,A)') '   radius increased to ', max_dist * 1.05_WP, ' m'
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
        r_xy = sqrt(max(0.0_WP, r_**2 - (zc - center_(3))**2))
        do j = 1, n_theta_
            theta_c = 0.5_WP * (theta_edges(j) + theta_edges(j+1))
            observers_(idx, 1) = center_(1) + r_xy * COS(theta_c)
            observers_(idx, 2) = center_(2) + r_xy * SIN(theta_c)
            observers_(idx, 3) = zc
            idx = idx + 1
        end do
    end do

    CALL self % check_observers_distances(observers_, 0.1_WP)

    WRITE(*,'(A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,I0,A,I0,A,I0)') &
        '>> [Sphere]     r=', r_, ' m | center=(', center_(1), ',', center_(2), ',', center_(3), &
        ') m | grid=', n_theta_, ' azim x ', nz_, ' z = ', Nobs

    CALL self % set_N_images(0_I32)
    CALL self % compute_pressure(observers_, self % block_size, p)
    CALL self % restore_default_images()

    if (self % debug) WRITE(*,'(A,F10.4)') '   Pressure norm : ', norm2(abs(p))

    ! Reshape pressure and observers to 3D (n_theta, nz, 3) and (NFreqs, n_theta, nz)
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
    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_sphere", p_reshaped)
    CALL save_variable(self % save_path, "Obs_sphere", obs_reshaped)
    CALL save_variable(self % save_path, "Theta_sphere", theta_deg_arr)
    CALL save_variable(self % save_path, "Z_sphere", z_centers)
    CALL save_variable(self % save_path, "R_sphere", r_)
    CALL save_variable(self % save_path, "Center_sphere", center_)
    CALL save_variable(self % save_path, "dA_sphere", dA)

    WRITE(*,'(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_sphere











END MODULE ImageMethod