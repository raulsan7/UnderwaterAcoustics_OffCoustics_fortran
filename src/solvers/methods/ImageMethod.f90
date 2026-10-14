! ===============================================================================
! MODULE: ImageMethod
! ===============================================================================
! Defines ImageMethod_t, the concrete acoustic solver extending abstract
! AcousticSolver_t. It implements the Method of Images for a fluid layer
! bounded by two horizontal reflecting planes (e.g. free surface above,
! seabed below): each real source, monopole or dipole, is mirrored into a
! chain of image sources above and below, with reflection coefficients and
! boundary conditions applied per bounce, then the pressure field is
! computed as the superposition of all real+image contributions at each
! observer point, with one free-field kernel per multipole order (with
! optional multi-observer chunked/parallel evaluation). With the volume
! absorption of the environment, the wavenumber of every kernel is the
! complex k + i*alpha, so each real or image contribution decays as
! exp(-alpha*r) along its own path.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE ImageMethod

USE Kinds
USE omp_lib
USE Emitter, ONLY: Emitter_t
USE HDF5Utils, ONLY: save_variable
USE FileUtils, ONLY: swap_extension
USE ConfigTypes, ONLY: Config_ImageMethod_t, Config_t
USE AcousticSolver, ONLY: AcousticSolver_t, pressure_norm
USE SourceTypes, ONLY: AcousticSource_t, SOURCE_DIPOLE, SOURCE_MONOPOLE

IMPLICIT NONE

PRIVATE
PUBLIC :: ImageMethod_t


TYPE :: ImageSystem_t
    ! -------------------------------------------------------------------------
    ! Helper type used internally by the Method of Images to hold the full set
    ! of real and image sources (the nodes) built from the original
    ! AcousticSource_t: their coordinates, complex strengths (vertical
    ! component flipped and attenuation applied once per reflection), the
    ! accumulated boundary-condition sign of each image (product of the BC
    ! signs of every plane it was reflected across; +1 for the real sources)
    ! and the multipole order of each node, that of its real source.
    ! -------------------------------------------------------------------------

    REAL(WP)    , ALLOCATABLE :: xyz(:,:)       ! [m] Nodes coordinates (Nnodes, 3)
    COMPLEX(WP) , ALLOCATABLE :: Force(:,:,:)   ! [N] Forces, or Q [m^3/s] in component 1 (Nfreqs, Nnodes, 3)
    INTEGER(I32), ALLOCATABLE :: BC_all(:)      ! [-] Accumulated boundary-condition sign of each node (Nnodes)
    INTEGER(I32), ALLOCATABLE :: order(:)       ! [-] Multipole order of each node, that of its real source (Nnodes)

END TYPE ImageSystem_t


TYPE, EXTENDS(AcousticSolver_t) :: ImageMethod_t
    ! -------------------------------------------------------------------------
    ! Concrete acoustic solver implementing the Method of Images: real
    ! sources (monopoles and dipoles) are reflected across an upper and a
    ! lower horizontal plane to build an equivalent free-space system
    ! (image_system) whose superposed monopole and dipole fields reproduce
    ! the effect of the bounded fluid layer on the pressure.
    ! -------------------------------------------------------------------------

    INTEGER(I32) :: N_images         = 30_I32       ! [-] Number of image levels to compute
    INTEGER(I32) :: default_N_images = 30_I32       ! [-] N_images read from config, restored after temporary changes
    INTEGER(I32) :: Upper_BC         = -1_I32       ! [-] Upper boundary condition, -1 for p = 0
    INTEGER(I32) :: Lower_BC         = +1_I32       ! [-] Lower boundary condition, +1 for dp = 0
    REAL(WP)     :: Lower_HBC        = -30.0_WP     ! [m] z-coordinate of the lower plane, - depth_base
    REAL(WP)     :: up_R             = 1.0_WP       ! [-] Attenuation coefficient for upper boundary reflections
    REAL(WP)     :: lw_R             = 1.0_WP       ! [-] Attenuation coefficient for lower boundary reflections
    REAL(WP)     :: eps              = 1.0e-6_WP    ! [m] Minimum source-observer distance used to floor r

    TYPE(ImageSystem_t) :: image_system             ! [-] Full real+image system (coords, forces, boundary tags, orders)

    CONTAINS

    ! --- Public type-bound procedures --- !
    PROCEDURE :: set_N_images
    PROCEDURE :: build_image_system
    PROCEDURE :: restore_default_images
    PROCEDURE :: method_of_images
    PROCEDURE :: dipole_pressure_images
    PROCEDURE :: monopole_pressure_images
    PROCEDURE :: run_sphere

    ! --- Deferred (abstract) procedures --- !
    PROCEDURE :: init              => init_ImageMethod
    PROCEDURE :: get_name          => get_name_ImageMethod
    PROCEDURE :: compute_pressure  => compute_pressure_ImageMethod
    PROCEDURE :: save_parameters   => save_parameters_ImageMethod
    PROCEDURE :: build_save_path   => build_save_path_ImageMethod
    PROCEDURE :: write_txt_summary => write_txt_summary_ImageMethod

END TYPE ImageMethod_t


CONTAINS


! ===============================================================================
! DEFERRED PROCEDURES
! ===============================================================================
SUBROUTINE init_ImageMethod(self, sources, config)
    ! -------------------------------------------------------------------------
    ! Constructor of the Image Method solver. Copies the shared solver/env/
    ! turbine data from the parent AcousticSolver_t setup, reads the method-
    ! specific parameters (N_images, boundary planes, reflection coefficients,
    ! singularity tolerance) from the configuration, tabulates the volume
    ! absorption alpha(f) if requested, resolves the output save path and
    ! builds the initial image system.
    ! -------------------------------------------------------------------------

    CLASS(ImageMethod_t)  , INTENT(INOUT)              :: self       ! [-] Solver instance
    TYPE(AcousticSource_t), INTENT(INOUT), ALLOCATABLE :: sources    ! [-] Acoustic sources (freed on exit)
    TYPE(Config_t)        , INTENT(INOUT), ALLOCATABLE :: config     ! [-] Global configuration (freed on exit)


    ! Parent variables (sources are moved, not copied)
    CALL MOVE_ALLOC(sources % xyz  , self % sources % xyz)
    CALL MOVE_ALLOC(sources % Force, self % sources % Force)
    CALL MOVE_ALLOC(sources % Freqs, self % sources % Freqs)
    CALL MOVE_ALLOC(sources % order, self % sources % order)

    ! Only the monopole and dipole kernels are implemented
    if (any(self % sources % order /= SOURCE_DIPOLE .and. self % sources % order /= SOURCE_MONOPOLE)) &
        error stop "ImageMethod % init(): only monopoles (order 0) and dipoles (order 1) are implemented"

    self % sources % geometry = sources % geometry
    self % env                = config % env
    self % run_inputs         = config % run_inputs
    self % save_path          = config % solver % save_path
    self % extra_label        = config % solver % extra_label
    self % block_size         = config % solver % block_size
    self % debug              = config % solver % debug
    self % in_cluster         = config % solver % in_cluster
    self % update_existing    = config % solver % update_existing
    self % mask_distance      = config % solver % mask_distance

    ! Child-specific variables
    SELECT TYPE (m => config % method)
    TYPE IS (Config_ImageMethod_t)
        self % N_images         = m % N_images
        self % default_N_images = m % N_images
        self % Upper_BC         = m % UpperBC
        self % Lower_BC         = m % LowerBC
        self % up_R             = m % attenuation_surface
        self % lw_R             = m % attenuation_bottom
        self % eps              = m % epsilon
    END SELECT

    ! Seabed plane built from the environment depth (single source of truth)
    if (self % env % depth_base <= 0.0_WP) error stop "ImageMethod % init(): depth_base must be strictly positive"
    self % Lower_HBC = - self % env % depth_base

    ! Sources inside the water column (the seabed is known from here on)
    CALL self % check_sources_in_water()

    ! Volume absorption of sea water at every frequency (only with env % absorption)
    CALL self % set_absorption()

    ! Build save_path
    CALL self % build_save_path()

    write(*,*) ''
    write(*,*) 'Selected method: ', self % get_name()
    CALL self % build_image_system()

    ! Cleanup memory
    DEALLOCATE(sources, config)

END SUBROUTINE init_ImageMethod


FUNCTION get_name_ImageMethod(self) RESULT(name)
    ! -------------------------------------------------------------------------
    ! Returns the human-readable solver name ("Images Method") and stores it
    ! in self % solver_name for later use (logs, summaries, HDF5 output).
    ! -------------------------------------------------------------------------

    CLASS(ImageMethod_t), INTENT(INOUT) :: self    ! [-] Solver instance
    CHARACTER(len=:)    , ALLOCATABLE   :: name    ! [-] Human-readable solver name


    name               = "Images Method"
    self % solver_name = trim(name)

END FUNCTION get_name_ImageMethod


SUBROUTINE compute_pressure_ImageMethod(self, observers, block_size, total_pressure)
    ! -------------------------------------------------------------------------
    ! Evaluates the total radiated pressure field at an arbitrary set of
    ! observer points, for every source frequency. Observers are processed
    ! either one-by-one (serial) or split into chunks distributed over OpenMP
    ! threads (in_cluster = .true.); the chunks split the observers, never the
    ! sources. Each chunk is handed to one kernel per multipole order present
    ! (dipole_pressure_images, monopole_pressure_images) together with the
    ! full real+image system: every kernel skips the nodes of the other
    ! orders, so the sources can be in any order, and the contributions are
    ! summed. With dipoles only, the operations are those of the dipole
    ! kernel alone.
    ! -------------------------------------------------------------------------

    CLASS(ImageMethod_t), INTENT(INOUT)            :: self                   ! [-] Solver instance
    REAL(WP)            , INTENT(IN)               :: observers(:,:)         ! [m] Observers coords (Nobs, 3)
    INTEGER(I32)        , INTENT(IN), OPTIONAL     :: block_size             ! [-] Observers per block
    COMPLEX(WP)         , INTENT(OUT), ALLOCATABLE :: total_pressure(:,:)    ! [Pa] Pressure (Nfreqs, Nobs)

    ! Local variables
    INTEGER(I32)              :: nobs, nf, bs, n_chunks, ichunk
    INTEGER(I32)              :: istart, iend, nb, i
    LOGICAL                   :: has_dipoles, has_monopoles
    COMPLEX(WP) , ALLOCATABLE :: p_block(:,:)


    ! Defaults
    bs = 256_I32 ; if (present(block_size)) bs = block_size

    ! Basic checks
    nobs = size(observers, 1)
    if (nobs <= 0_I32) then
        ALLOCATE(total_pressure(0, 0))
        return
    end if
    bs = max(1_I32, min(bs, (nobs + omp_get_max_threads() - 1_I32) / omp_get_max_threads()))

    if (size(observers, 2) /= 3) error stop "compute_pressure_ImageMethod: observers must have shape (N,3)"

    ! Sources must be defined
    if (.not. ALLOCATED(self % sources % xyz) .or. .not. ALLOCATED(self % sources % Force) .or. &
        .not. ALLOCATED(self % sources % Freqs)) then
        error stop "compute_pressure_ImageMethod: no sources defined"
    end if

    ! Number of frequencies
    nf = size(self % sources % Freqs)

    ALLOCATE(total_pressure(nf, nobs))
    total_pressure = (0.0_WP, 0.0_WP)

    if (.not. self % in_cluster) bs = 1_I32

    ! Multipole orders present in the image system: only their kernels are called
    has_dipoles   = any(self % image_system % order == SOURCE_DIPOLE)
    has_monopoles = any(self % image_system % order == SOURCE_MONOPOLE)

    if (self % in_cluster) then
        n_chunks = (nobs + bs - 1_I32) / bs

        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP PRIVATE(ichunk, istart, iend, nb, p_block) &
        !$OMP SHARED(n_chunks, nobs, total_pressure, self)
        do ichunk = 0, n_chunks - 1
            istart = ichunk * bs + 1_I32
            iend   = min(istart + bs - 1_I32, nobs)
            nb     = iend - istart + 1_I32

            ! Each kernel is called for this chunk only if there are sources of its order, and its
            ! result is accumulated before the next call reuses p_block (non-overlapping slices -> no race)
            if (has_dipoles) then
                CALL self % dipole_pressure_images(observers(istart:iend, :), self % sources % Freqs, &
                                                   self % image_system % xyz, self % image_system % Force, &
                                                   self % image_system % BC_all, self % image_system % order, p_block)
                total_pressure(:, istart:iend) = total_pressure(:, istart:iend) + p_block(:, 1:nb)
                if (ALLOCATED(p_block)) DEALLOCATE(p_block)
            end if

            if (has_monopoles) then
                CALL self % monopole_pressure_images(observers(istart:iend, :), self % sources % Freqs, &
                                                     self % image_system % xyz, self % image_system % Force, &
                                                     self % image_system % BC_all, self % image_system % order, p_block)
                total_pressure(:, istart:iend) = total_pressure(:, istart:iend) + p_block(:, 1:nb)
                if (ALLOCATED(p_block)) DEALLOCATE(p_block)
            end if
        end do
        !$OMP END PARALLEL DO

    else    ! .not. self % in_cluster

        ! Serial: loop observers (keep 2D slice for compatibility)
        do i = 1, nobs
            if (has_dipoles) then
                CALL self % dipole_pressure_images(observers(i:i, :), self % sources % Freqs, &
                                                   self % image_system % xyz, self % image_system % Force, &
                                                   self % image_system % BC_all, self % image_system % order, p_block)
                total_pressure(:, i) = total_pressure(:, i) + p_block(:, 1)
                if (ALLOCATED(p_block)) DEALLOCATE(p_block)
            end if

            if (has_monopoles) then
                CALL self % monopole_pressure_images(observers(i:i, :), self % sources % Freqs, &
                                                     self % image_system % xyz, self % image_system % Force, &
                                                     self % image_system % BC_all, self % image_system % order, p_block)
                total_pressure(:, i) = total_pressure(:, i) + p_block(:, 1)
                if (ALLOCATED(p_block)) DEALLOCATE(p_block)
            end if

            ! Print progress every 100 observers
            if (mod(i, 100_I32) == 0_I32) then
                write(*, '(A,I0,A,I0)') '  Progress: ', i, '/', nobs
                flush(6)
            end if
        end do

    end if

END SUBROUTINE compute_pressure_ImageMethod


SUBROUTINE save_parameters_ImageMethod(self, emitter_obj)
    ! -------------------------------------------------------------------------
    ! Saves the solver and environment parameters to the HDF5 output file,
    ! delegates the emitter parameters to emitter % save_params, and triggers
    ! the human-readable .dat text summary via write_txt_summary. An existing
    ! output kept for update is checked first (check_update).
    ! -------------------------------------------------------------------------

    CLASS(ImageMethod_t), INTENT(IN) :: self       ! [-] Solver instance
    CLASS(Emitter_t)    , INTENT(IN) :: emitter_obj    ! [-] Emitter that produced the sources


    ! Existing output of another case: deleted here, before anything is written
    CALL self % check_update(emitter_obj)

    ! Solver / environment / generic turbine data -> HDF5
    CALL save_variable(self % save_path, "N_images" , self % N_images)
    CALL save_variable(self % save_path, "c_wat"    , self % env % speed_sound_base)
    CALL save_variable(self % save_path, "rho_wat"  , self % env % density_base)
    CALL save_variable(self % save_path, "Lower_HBC", self % Lower_HBC)
    CALL save_variable(self % save_path, "p_ref"    , self % env % reference_pressure)

    ! Datasets shared by every solver (Method, ENBW) -> HDF5
    CALL self % save_shared_parameters(emitter_obj)

    ! Emitter data -> HDF5
    CALL emitter_obj % save_params(self % save_path)

    ! Human-readable .dat summary, delegated to the solver-specific implementation
    CALL self % write_txt_summary(emitter_obj)

END SUBROUTINE save_parameters_ImageMethod


SUBROUTINE build_save_path_ImageMethod(self)
    ! -------------------------------------------------------------------------
    ! Resolves the final HDF5 output path through the shared
    ! resolve_save_path, with this solver's method tag ("_IM" followed by the
    ! number of images). Folder and file name labels come from the emitter.
    ! -------------------------------------------------------------------------

    CLASS(ImageMethod_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    CHARACTER(len=16) :: imag_str


    write(imag_str, '(I0)') self % N_images
    CALL self % resolve_save_path("_IM" // trim(imag_str))

END SUBROUTINE build_save_path_ImageMethod


SUBROUTINE write_txt_summary_ImageMethod(self, emitter_obj)
    ! -------------------------------------------------------------------------
    ! Writes a plain-text (.dat) summary next to the HDF5 output, mirroring
    ! the saved emitter, environment and Image Method solver parameters in a
    ! human-readable form for quick inspection without opening the HDF5 file.
    ! -------------------------------------------------------------------------

    CLASS(ImageMethod_t), INTENT(IN) :: self       ! [-] Solver instance
    CLASS(Emitter_t)    , INTENT(IN) :: emitter_obj    ! [-] Emitter that produced the sources

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: dat_path
    INTEGER(I32)                  :: u


    dat_path = swap_extension(trim(self % save_path), ".dat")

    open(newunit=u, file=dat_path, status='replace', action='write')

    write(u, '(A)') repeat('=', 70)
    write(u, '(A)') 'OFF_COUSTICS SIMULATION SUMMARY'
    write(u, '(A)') repeat('=', 70)
    write(u, '(A)') 'HDF5 file: ' // trim(self % save_path)
    write(u, *)

    CALL emitter_obj % write_txt(u)

    write(u, '(A)') repeat('-', 70)
    write(u, '(A)') 'ENVIRONMENT'
    write(u, '(A)') repeat('-', 70)
    write(u, '(A,F8.2,A)')   'Speed of sound   : ', self % env % speed_sound_base, ' m/s'
    write(u, '(A,F8.2,A)')   'Density          : ', self % env % density_base, ' kg/m^3'
    write(u, '(A,ES10.2,A)') 'Reference press. : ', self % env % reference_pressure, ' Pa'
    write(u, '(A,F8.2,A)')   'Depth            : ', self % env % depth_base, ' m'
    CALL self % write_txt_water(u)
    write(u, *)

    write(u, '(A)') repeat('-', 70)
    write(u, '(A)') 'SOLVER: ' // trim(self % solver_name)
    write(u, '(A)') repeat('-', 70)
    write(u, '(A,I0)')     'N images     : ', self % N_images
    write(u, '(A,F8.2,A)') 'Lower HBC    : ', self % Lower_HBC, ' m'
    write(u, '(A,F8.4)')   'Upper atten. : ', self % up_R
    write(u, '(A,F8.4)')   'Lower atten. : ', self % lw_R
    write(u, '(A,I0)')     'Block size   : ', self % block_size
    write(u, '(A,L1)')     'In cluster   : ', self % in_cluster
    write(u, *)

    close(u)

    write(*,*) '   ---> Text summary saved at ', dat_path

END SUBROUTINE write_txt_summary_ImageMethod


! ===============================================================================
! HELPERS
! ===============================================================================
SUBROUTINE build_image_system(self)
    ! -------------------------------------------------------------------------
    ! Builds (or rebuilds) self % image_system from the real sources: calls
    ! method_of_images directly on the source arrays to generate the full
    ! chain of image z-positions, forces and boundary-condition tags, and
    ! assembles the resulting (x,y,z) coordinates by pairing each image with
    ! its parent source's horizontal position. Every node takes the multipole
    ! order of its parent source, which tells the kernels whether it is a
    ! monopole or a dipole. The image forces and tags are moved (not copied)
    ! into the image system, and any previous image system is freed first, so
    ! that old and new never coexist in memory.
    ! -------------------------------------------------------------------------

    CLASS(ImageMethod_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    REAL(WP)    , ALLOCATABLE :: z_all(:)
    COMPLEX(WP) , ALLOCATABLE :: Force_out(:,:,:)
    INTEGER(I32), ALLOCATABLE :: parent(:), BC_all(:)
    INTEGER(I32)              :: Nnodes, j


    ! Sources are always stored as (NFreqs, Nsources, 3)
    if (size(self % sources % Force, 2) /= size(self % sources % xyz, 1) .or. size(self % sources % Force, 3) /= 3) &
        error stop "ImageMethod % build_image_system(): sources % Force must have shape (NFreqs, Nsources, 3)"

    ! Rebuilding: free the previous image system first
    if (ALLOCATED(self % image_system % xyz))    DEALLOCATE(self % image_system % xyz)
    if (ALLOCATED(self % image_system % Force))  DEALLOCATE(self % image_system % Force)
    if (ALLOCATED(self % image_system % BC_all)) DEALLOCATE(self % image_system % BC_all)
    if (ALLOCATED(self % image_system % order))  DEALLOCATE(self % image_system % order)

    ! Source arrays passed directly (INTENT(IN) dummies), no local copies
    CALL self % method_of_images(self % sources % xyz(:, 3), self % sources % Force, &
                                 z_all, Force_out, parent, BC_all)

    ! Image coordinates: horizontal position of the parent source, reflected z
    Nnodes = size(z_all)
    ALLOCATE(self % image_system % xyz(Nnodes, 3))
    do j = 1, Nnodes
        self % image_system % xyz(j, 1) = self % sources % xyz(parent(j), 1)
        self % image_system % xyz(j, 2) = self % sources % xyz(parent(j), 2)
        self % image_system % xyz(j, 3) = z_all(j)
    end do

    ! Multipole order of every node: that of its parent source
    self % image_system % order = self % sources % order(parent)

    ! Image forces and boundary-condition tags are moved in, not copied
    CALL MOVE_ALLOC(Force_out, self % image_system % Force)
    CALL MOVE_ALLOC(BC_all   , self % image_system % BC_all)

END SUBROUTINE build_image_system


SUBROUTINE set_N_images(self, new_N_images)
    ! -------------------------------------------------------------------------
    ! Changes the number of image levels and rebuilds the image system to
    ! reflect it. Used e.g. to temporarily set N_images = 0 (free-field, no
    ! reflections) for calibration runs such as run_sphere.
    ! -------------------------------------------------------------------------

    CLASS(ImageMethod_t), INTENT(INOUT) :: self            ! [-] Solver instance
    INTEGER(I32)        , INTENT(IN)    :: new_N_images    ! [-] New number of image levels


    self % N_images = new_N_images
    CALL self % build_image_system()

END SUBROUTINE set_N_images


SUBROUTINE restore_default_images(self)
    ! -------------------------------------------------------------------------
    ! Restores N_images to its original (config-loaded) value and rebuilds
    ! the image system if it had been temporarily changed by set_N_images.
    ! -------------------------------------------------------------------------

    CLASS(ImageMethod_t), INTENT(INOUT) :: self    ! [-] Solver instance


    if (self % N_images /= self % default_N_images) then
        self % N_images = self % default_N_images
        CALL self % build_image_system()
    end if

END SUBROUTINE restore_default_images


! ===============================================================================
! METHOD (CORE PHYSICS)
! ===============================================================================
SUBROUTINE method_of_images(self, zi, Force, z_all, Force_out, parent, BC_all)
    ! -------------------------------------------------------------------------
    ! Builds the full chain of images for each real source by repeated
    ! reflection across the upper and lower horizontal planes, N_images times
    ! in each direction. Starting from the real source, the upper chain
    ! alternately reflects across the lower then upper plane (and vice versa
    ! for the lower chain). Each reflection flips the sign of the vertical
    ! force component and multiplies the image by the reflection coefficient
    ! of the plane just crossed: its boundary condition (+-1, Upper_BC or
    ! Lower_BC) times its attenuation (up_R or lw_R). Both factors accumulate
    ! along the chain, so an image reached after n_up upper and n_lw lower
    ! reflections carries (Upper_BC*up_R)**n_up * (Lower_BC*lw_R)**n_lw: the
    ! accumulated sign is stored in BC_all and the accumulated attenuation in
    ! the image force. The same chain serves monopoles: their volume velocity
    ! Q is in the first component (SourceTypes), which is only scaled by the
    ! attenuation, and with the BC sign this is the image of a monopole; the
    ! flip of the third component does not touch it (it is 0). The result is a
    ! flat list of z-positions, forces and BC signs (real source + all images)
    ! later consumed by dipole_pressure_images and monopole_pressure_images.
    ! -------------------------------------------------------------------------

    CLASS(ImageMethod_t), INTENT(IN)               :: self                ! [-] Solver instance
    REAL(WP)            , INTENT(IN)               :: zi(:)               ! [m] z of the real sources (Nsources)
    COMPLEX(WP)         , INTENT(IN)               :: Force(:,:,:)        ! [N] Real source forces (NFreqs, Nsources, 3)
    REAL(WP)            , INTENT(OUT), ALLOCATABLE :: z_all(:)            ! [m] z of the real + image sources (Nnodes)
    COMPLEX(WP)         , INTENT(OUT), ALLOCATABLE :: Force_out(:,:,:)    ! [N] Real + image forces (NFreqs, Nnodes, 3)
    INTEGER(I32)        , INTENT(OUT), ALLOCATABLE :: parent(:)           ! [-] Parent source of each node (Nnodes)
    INTEGER(I32)        , INTENT(OUT), ALLOCATABLE :: BC_all(:)           ! [-] Cumulative BC sign per node (Nnodes)

    ! Local variables
    INTEGER(I32)                  :: Nsources, NFreqs, Nnodes, idx, i_s, img_lvl
    REAL(WP)                      :: z_upper, z_lower, att_u, att_l, Upper_HBC
    COMPLEX(WP)     , ALLOCATABLE :: F_upper(:,:), F_lower(:,:)
    CHARACTER(len=5)              :: last_plane_upper, last_plane_lower
    INTEGER(I32)                  :: BC_u, BC_l

    Upper_HBC = 0.0_WP
    Nsources  = size(zi)
    NFreqs    = size(Force, 1)
    Nnodes    = Nsources * (1_I32 + 2_I32 * self % N_images)

    ALLOCATE(z_all(Nnodes))
    ALLOCATE(Force_out(NFreqs, Nnodes, 3))
    ALLOCATE(parent(Nnodes))
    ALLOCATE(BC_all(Nnodes))
    ALLOCATE(F_upper(NFreqs, 3))
    ALLOCATE(F_lower(NFreqs, 3))

    idx = 1_I32
    do i_s = 1, Nsources
        ! Real source node
        z_all(idx)           = zi(i_s)
        Force_out(:, idx, :) = Force(:, i_s, :)
        BC_all(idx)          = 1_I32
        parent(idx)          = i_s
        idx                  = idx + 1_I32

        z_upper = zi(i_s)
        z_lower = zi(i_s)
        F_upper = Force(:, i_s, :)
        F_lower = Force(:, i_s, :)
        BC_u    = 1_I32
        BC_l    = 1_I32

        last_plane_upper = "upper"
        last_plane_lower = "lower"

        do img_lvl = 1, self % N_images
            ! --- 1. Upper reflection chain ---
            if (last_plane_upper == "upper") then
                z_upper          = 2.0_WP * self % Lower_HBC - z_upper
                BC_u             = BC_u * self % Lower_BC
                att_u            = self % lw_R
                last_plane_upper = "lower"
            else
                z_upper          = 2.0_WP * Upper_HBC - z_upper
                BC_u             = BC_u * self % Upper_BC
                att_u            = self % up_R
                last_plane_upper = "upper"
            end if

            F_upper(:, 3) = -1.0_WP * F_upper(:, 3)
            F_upper       = F_upper * att_u

            z_all(idx)           = z_upper
            Force_out(:, idx, :) = F_upper
            BC_all(idx)          = BC_u
            parent(idx)          = i_s
            idx                  = idx + 1_I32

            ! --- 2. Lower reflection chain ---
            if (last_plane_lower == "lower") then
                z_lower          = 2.0_WP * Upper_HBC - z_lower
                BC_l             = BC_l * self % Upper_BC
                att_l            = self % up_R
                last_plane_lower = "upper"
            else
                z_lower          = 2.0_WP * self % Lower_HBC - z_lower
                BC_l             = BC_l * self % Lower_BC
                att_l            = self % lw_R
                last_plane_lower = "lower"
            end if

            F_lower(:, 3) = -1.0_WP * F_lower(:, 3)
            F_lower       = F_lower * att_l

            z_all(idx)           = z_lower
            Force_out(:, idx, :) = F_lower
            BC_all(idx)          = BC_l
            parent(idx)          = i_s
            idx                  = idx + 1_I32
        end do
    end do

END SUBROUTINE method_of_images


SUBROUTINE dipole_pressure_images(self, obs_pos, Freqs, nodes_pos, Force, BC_all, order, p_out)
    ! -------------------------------------------------------------------------
    ! Evaluates the acoustic dipole Green's function for a block of observers,
    ! summed over every real and image dipole, at every frequency; the nodes
    ! of another multipole order are skipped. For each source-observer pair
    ! it computes the separation vector and distance (floored at eps to avoid
    ! the 1/r singularity), then the far/near-field dipole pressure term
    ! F.r_hat * (1/r - i*k) * exp(i*k*r) / (4*pi*r), weighted by the source's
    ! boundary-condition sign (BC_all), and accumulates it into p_out. The
    ! frequency loop is the innermost one: p_out(:, k), Force(:, j, :) and the
    ! wavenumbers are accessed with unit stride, the frequency-independent
    ! factors of each pair are computed once, and exp(i*k*r) is evaluated as
    ! cmplx(cos(k*r), sin(k*r)), so the loop is vectorized with SIMD sin/cos.
    ! The sources are accumulated in the same order as before (j ascending).
    ! With the volume absorption of the environment the wavenumber is the
    ! complex k + i*alpha: the term gets exp(-alpha*r) and its near-field
    ! factor becomes 1/r + alpha - i*k (a separate loop, so the one without
    ! absorption is unchanged).
    ! -------------------------------------------------------------------------

    CLASS(ImageMethod_t), INTENT(IN)               :: self              ! [-] Solver instance
    REAL(WP)            , INTENT(IN)               :: obs_pos(:,:)      ! [m] Observers coords (Nob, 3)
    REAL(WP)            , INTENT(IN)               :: Freqs(:)          ! [Hz] Frequency vector (NFreqs)
    REAL(WP)            , INTENT(IN)               :: nodes_pos(:,:)    ! [m] Real + image coords (Nnodes, 3)
    COMPLEX(WP)         , INTENT(IN)               :: Force(:,:,:)      ! [N] Real + image forces (NFreqs, Nnodes, 3)
    INTEGER(I32)        , INTENT(IN)               :: BC_all(:)         ! [-] BC sign of each node (Nnodes)
    INTEGER(I32)        , INTENT(IN)               :: order(:)          ! [-] Multipole order of each node (Nnodes)
    COMPLEX(WP)         , INTENT(OUT), ALLOCATABLE :: p_out(:,:)        ! [Pa] Pressure (NFreqs, Nob)

    ! Local variables
    INTEGER(I32)              :: NFreqs, Nnodes, Nob
    INTEGER(I32)              :: j, f, k
    LOGICAL                   :: absorbing
    REAL(WP)    , ALLOCATABLE :: k_wave(:)
    REAL(WP)                  :: dx, dy, dz, r, r_inv, scal, kr, epsv
    COMPLEX(WP)               :: F_dot_r


    ! Dimensions
    NFreqs = size(Freqs)
    Nnodes = size(nodes_pos, 1)
    Nob    = size(obs_pos, 1)

    epsv   = self % eps

    ALLOCATE(p_out(NFreqs, Nob))
    p_out     = (0.0_WP, 0.0_WP)
    k_wave    = 2.0_WP * PI * Freqs / self % env % speed_sound_base
    absorbing = self % env % absorption

    do j = 1, Nnodes
        if (order(j) /= SOURCE_DIPOLE) cycle    ! Monopoles: monopole_pressure_images

        do k = 1, Nob
            ! Separation, distance (floored at eps) and frequency-independent factors of this pair
            dx = obs_pos(k, 1) - nodes_pos(j, 1)
            dy = obs_pos(k, 2) - nodes_pos(j, 2)
            dz = obs_pos(k, 3) - nodes_pos(j, 3)
            r  = sqrt(dx * dx + dy * dy + dz * dz)
            if (r < epsv) r = epsv
            r_inv = 1.0_WP / r
            scal  = real(BC_all(j), WP) * r_inv * r_inv / (4.0_WP * PI)

            ! Frequencies innermost: unit stride, vectorized (SIMD sin/cos)
            if (absorbing) then
                ! Volume absorption: complex wavenumber k + i*alpha
                !$OMP SIMD PRIVATE(F_dot_r, kr)
                do f = 1, NFreqs
                    kr          = k_wave(f) * r
                    F_dot_r     = Force(f, j, 1) * dx + Force(f, j, 2) * dy + Force(f, j, 3) * dz
                    p_out(f, k) = p_out(f, k) + F_dot_r * cmplx(r_inv + self % alpha(f), -k_wave(f), kind=WP) * &
                                  cmplx(cos(kr), sin(kr), kind=WP) * exp(-self % alpha(f) * r) * scal
                end do
            else
                !$OMP SIMD PRIVATE(F_dot_r, kr)
                do f = 1, NFreqs
                    kr          = k_wave(f) * r
                    F_dot_r     = Force(f, j, 1) * dx + Force(f, j, 2) * dy + Force(f, j, 3) * dz
                    p_out(f, k) = p_out(f, k) + F_dot_r * cmplx(r_inv, -k_wave(f), kind=WP) * &
                                  cmplx(cos(kr), sin(kr), kind=WP) * scal
                end do
            end if
        end do
    end do

END SUBROUTINE dipole_pressure_images


SUBROUTINE monopole_pressure_images(self, obs_pos, Freqs, nodes_pos, Force, BC_all, order, p_out)
    ! -------------------------------------------------------------------------
    ! Evaluates the acoustic monopole Green's function for a block of
    ! observers, summed over every real and image monopole, at every
    ! frequency; the nodes of another multipole order are skipped. A monopole
    ! of volume velocity Q (first component of Force, see SourceTypes)
    ! radiates p = -i*w*rho*Q * exp(i*k*r) / (4*pi*r), weighted by its
    ! boundary-condition sign (BC_all); the attenuation of the reflections is
    ! already in Force. The frequency factor -i*w*rho/(4*pi) = -i*f*rho/2 is
    ! computed once per frequency and the distance factor BC/r once per pair
    ! (r floored at eps). Same loop structure as dipole_pressure_images:
    ! frequencies innermost, with unit stride and SIMD sin/cos. With the
    ! volume absorption of the environment the wavenumber is the complex
    ! k + i*alpha, so the term gets exp(-alpha*r) (a separate loop, so the one
    ! without absorption is unchanged).
    ! -------------------------------------------------------------------------

    CLASS(ImageMethod_t), INTENT(IN)               :: self              ! [-] Solver instance
    REAL(WP)            , INTENT(IN)               :: obs_pos(:,:)      ! [m] Observers coords (Nob, 3)
    REAL(WP)            , INTENT(IN)               :: Freqs(:)          ! [Hz] Frequency vector (NFreqs)
    REAL(WP)            , INTENT(IN)               :: nodes_pos(:,:)    ! [m] Real + image coords (Nnodes, 3)
    COMPLEX(WP)         , INTENT(IN)               :: Force(:,:,:)      ! [m^3/s] Q in component 1 (NFreqs, Nnodes, 3)
    INTEGER(I32)        , INTENT(IN)               :: BC_all(:)         ! [-] BC sign of each node (Nnodes)
    INTEGER(I32)        , INTENT(IN)               :: order(:)          ! [-] Multipole order of each node (Nnodes)
    COMPLEX(WP)         , INTENT(OUT), ALLOCATABLE :: p_out(:,:)        ! [Pa] Pressure (NFreqs, Nob)

    ! Local variables
    INTEGER(I32)              :: NFreqs, Nnodes, Nob
    INTEGER(I32)              :: j, f, k
    LOGICAL                   :: absorbing
    REAL(WP)    , ALLOCATABLE :: k_wave(:)
    COMPLEX(WP) , ALLOCATABLE :: src_fac(:)
    REAL(WP)                  :: dx, dy, dz, r, scal, kr, epsv


    ! Dimensions
    NFreqs = size(Freqs)
    Nnodes = size(nodes_pos, 1)
    Nob    = size(obs_pos, 1)

    epsv   = self % eps

    ALLOCATE(p_out(NFreqs, Nob))
    p_out     = (0.0_WP, 0.0_WP)
    k_wave    = 2.0_WP * PI * Freqs / self % env % speed_sound_base
    src_fac   = cmplx(0.0_WP, -0.5_WP * Freqs * self % env % density_base, kind=WP)    ! -i*w*rho/(4*pi)
    absorbing = self % env % absorption

    do j = 1, Nnodes
        if (order(j) /= SOURCE_MONOPOLE) cycle    ! Dipoles: dipole_pressure_images

        do k = 1, Nob
            ! Distance (floored at eps) and frequency-independent factor of this pair
            dx = obs_pos(k, 1) - nodes_pos(j, 1)
            dy = obs_pos(k, 2) - nodes_pos(j, 2)
            dz = obs_pos(k, 3) - nodes_pos(j, 3)
            r  = sqrt(dx * dx + dy * dy + dz * dz)
            if (r < epsv) r = epsv
            scal = real(BC_all(j), WP) / r

            ! Frequencies innermost: unit stride, vectorized (SIMD sin/cos)
            if (absorbing) then
                ! Volume absorption: complex wavenumber k + i*alpha
                !$OMP SIMD PRIVATE(kr)
                do f = 1, NFreqs
                    kr          = k_wave(f) * r
                    p_out(f, k) = p_out(f, k) + Force(f, j, 1) * src_fac(f) * cmplx(cos(kr), sin(kr), kind=WP) * &
                                  exp(-self % alpha(f) * r) * scal
                end do
            else
                !$OMP SIMD PRIVATE(kr)
                do f = 1, NFreqs
                    kr          = k_wave(f) * r
                    p_out(f, k) = p_out(f, k) + Force(f, j, 1) * src_fac(f) * cmplx(cos(kr), sin(kr), kind=WP) * scal
                end do
            end if
        end do
    end do

END SUBROUTINE monopole_pressure_images


! ===============================================================================
! RUN SPHERE (overriden from AcousticSolver_t)
! ===============================================================================
SUBROUTINE run_sphere(self)
    ! -------------------------------------------------------------------------
    ! Builds a spherical (theta, z) grid of observer points around the
    ! sources, of the requested radius (used as it is) or, when it is 0 or
    ! -999, of the automatic one: the distance from the centre to the
    ! farthest source plus 3/k at f_farfield.
    ! Then evaluates the free-field pressure on that surface (temporarily
    ! forcing N_images = 0, i.e. no boundary reflections, and switching off
    ! the volume absorption, since the sphere measures the power radiated by
    ! the sources) for later far-field/radiated-power post-processing.
    ! Reshapes and saves the pressure field, observer coordinates,
    ! angular/vertical grids and surface element area to the HDF5 output.
    ! -------------------------------------------------------------------------

    CLASS(ImageMethod_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: p(:,:)
    COMPLEX(WP) , ALLOCATABLE :: p_reshaped(:,:,:)
    REAL(WP)    , ALLOCATABLE :: observers_(:,:)
    REAL(WP)    , ALLOCATABLE :: obs_reshaped(:,:,:)
    REAL(WP)    , ALLOCATABLE :: theta_deg_arr(:)
    REAL(WP)    , ALLOCATABLE :: z_centers(:)
    REAL(WP)    , ALLOCATABLE :: z_edges(:), theta_edges(:)
    REAL(WP)                  :: r_, center_(3)
    REAL(WP)                  :: dz, dtheta, zc, r_xy, theta_c, a_max
    REAL(WP)                  :: dA, f_farfield_, k_far
    INTEGER(I32)              :: n_theta_, nz_, i, j, idx, Nobs
    LOGICAL                   :: absorbing


    ! Defaults (resolved once, persisted back into run_inputs)
    associate(pp => self % run_inputs % params_sphere)
        if (pp % r_sphere == DEFAULT_R)          pp % r_sphere          = 0.0_WP
        if (pp % n_theta_sphere == DEFAULT_I)    pp % n_theta_sphere    = 144_I32
        if (pp % nz_sphere == DEFAULT_I)         pp % nz_sphere         = 40_I32
        if (pp % f_farfield_sphere == DEFAULT_R) pp % f_farfield_sphere = 10.0_WP
        if (any(pp % center_sphere == DEFAULT_R)) &
            pp % center_sphere = 0.5_WP * (minval(self % sources % xyz, dim=1) + maxval(self % sources % xyz, dim=1))

        r_          = pp % r_sphere
        n_theta_    = pp % n_theta_sphere
        nz_         = pp % nz_sphere
        f_farfield_ = pp % f_farfield_sphere
        center_     = pp % center_sphere
    end associate

    ! Basic enclosure check
    if (r_ == 0.0_WP) then
        a_max = 0.0_WP
        do i = 1, size(self % sources % xyz, 1)
            a_max = max(a_max, norm2(self % sources % xyz(i, :) - center_))
        end do

        k_far = 2.0_WP * PI * f_farfield_ / self % env % speed_sound_base
        r_ = a_max + 3.0_WP / k_far     ! [-] Farfield criterion: 3
    end if

    self % run_inputs % params_sphere % r_sphere = r_

    Nobs = n_theta_ * nz_
    ALLOCATE(observers_(Nobs, 3))
    ALLOCATE(theta_deg_arr(n_theta_))
    ALLOCATE(z_centers(nz_), z_edges(nz_+1), theta_edges(n_theta_+1))

    dz = (2.0_WP * r_) / real(nz_, WP)
    do i = 1, nz_ + 1
        z_edges(i) = (center_(3) - r_) + real(i - 1_I32, WP) * dz
    end do
    do i = 1, nz_
        z_centers(i) = 0.5_WP * (z_edges(i) + z_edges(i+1))
    end do

    dtheta = (2.0_WP * PI) / real(n_theta_, WP)
    do j = 1, n_theta_ + 1
        theta_edges(j) = real(j - 1_I32, WP) * dtheta
    end do
    do j = 1, n_theta_
        theta_deg_arr(j) = 0.5_WP * (theta_edges(j) + theta_edges(j+1)) * 180.0_WP / PI
    end do

    idx = 1_I32
    do i = 1, nz_
        zc   = z_centers(i)
        r_xy = sqrt(max(0.0_WP, r_**2 - (zc - center_(3))**2))
        do j = 1, n_theta_
            theta_c            = 0.5_WP * (theta_edges(j) + theta_edges(j+1))
            observers_(idx, 1) = center_(1) + r_xy * cos(theta_c)
            observers_(idx, 2) = center_(2) + r_xy * sin(theta_c)
            observers_(idx, 3) = zc
            idx                = idx + 1_I32
        end do
    end do

    write(*, '(A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,I0,A,I0,A,I0)') &
        '>> [Sphere]     r=', r_, ' m | center=(', center_(1), ',', center_(2), ',', center_(3), &
        ') m | grid=', n_theta_, ' azim x ', nz_, ' z = ', Nobs

    absorbing               = self % env % absorption
    self % env % absorption = .false.
    CALL self % set_N_images(0_I32)
    CALL self % compute_pressure_masked('Sphere', observers_, self % block_size, p)
    CALL self % restore_default_images()
    self % env % absorption = absorbing

    if (self % debug) write(*, '(A,F10.4)') '   Pressure norm : ', pressure_norm(p)

    ALLOCATE(p_reshaped(size(p, 1), n_theta_, nz_))
    ALLOCATE(obs_reshaped(n_theta_, nz_, 3))

    idx = 1_I32
    do i = 1, nz_
        do j = 1, n_theta_
            obs_reshaped(j, i, :) = observers_(idx, :)
            p_reshaped(:, j, i)   = p(:, idx)
            idx                   = idx + 1_I32
        end do
    end do

    dA = (4.0_WP * PI * r_**2) / real(n_theta_ * nz_, WP)

    CALL save_variable(self % save_path, "Freqs", self % sources % Freqs)
    CALL save_variable(self % save_path, "P_sphere", p_reshaped)
    CALL save_variable(self % save_path, "Obs_sphere", obs_reshaped)
    CALL save_variable(self % save_path, "Theta_sphere", theta_deg_arr)
    CALL save_variable(self % save_path, "Z_sphere", z_centers)
    CALL save_variable(self % save_path, "R_sphere", r_)
    CALL save_variable(self % save_path, "Center_sphere", center_)
    CALL save_variable(self % save_path, "dA_sphere", dA)

    write(*, '(A)') '   Saved to: ' // trim(self % save_path)

END SUBROUTINE run_sphere

END MODULE ImageMethod
