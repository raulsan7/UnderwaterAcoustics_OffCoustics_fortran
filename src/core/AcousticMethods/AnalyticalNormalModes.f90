! ===============================================================================
! MODULE: AnalyticalNormalModes
! ===============================================================================
! Defines AnalyticalNormalModes_t, the concrete acoustic solver extending
! abstract AcousticSolver_t. It models the shallow-water column as an ideal
! waveguide bounded by a pressure-release surface (p = 0) above and a rigid
! seabed (dp/dz = 0) below, and expands the field radiated by each real
! dipole source as a superposition of analytical normal modes: a vertical
! modal shape Psi_m(z) times a horizontal cylindrical wave given by Hankel
! functions of the first kind, either propagating (real k_rm) or evanescent
! (imaginary k_rm). The pressure is computed at each observer point with the
! same chunked/parallel evaluation scheme used by the Image Method.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE AnalyticalNormalModes

USE Kinds
USE Types
USE IOUtils
USE omp_lib
USE WindFarm, ONLY: WindFarm_t
USE WindTurbine, ONLY: WindTurbine_t
USE MathUtils, ONLY: compute_hankel_complex
USE AcousticSolver, ONLY: AcousticSolver_t, save_emitter_params, write_emitter_txt

IMPLICIT NONE

PRIVATE
PUBLIC :: AnalyticalNormalModes_t


TYPE, EXTENDS(AcousticSolver_t) :: AnalyticalNormalModes_t
    ! -------------------------------------------------------------------------
    ! Concrete acoustic solver implementing the analytical normal-mode expansion
    ! of the dipole field in an ideal waveguide of depth H = Upper_HBC - Lower_HBC
    ! (pressure-release top, rigid bottom). Each mode m has a fixed vertical
    ! wavenumber k_zm = (2m-1)*pi/(2H) and a frequency-dependent horizontal
    ! wavenumber k_rm = sqrt(k^2 - k_zm^2).
    ! -------------------------------------------------------------------------

    INTEGER(I32) :: Nmodes      = DEFAULT_I     ! [-] Number of normal modes retained in the expansion
    LOGICAL      :: auto_Nmodes = .false.       ! [-] Whether Nmodes was computed automatically from f_max
    REAL(WP)     :: H           = 30.0_WP       ! [m] Waveguide height (distance between both boundaries)
    REAL(WP)     :: Upper_HBC   = 0.0_WP        ! [m] z-coordinate of the upper boundary (free surface, p = 0)
    REAL(WP)     :: Lower_HBC   = -30.0_WP      ! [m] z-coordinate of the lower boundary (seabed, dp/dz = 0)
    REAL(WP)     :: eps         = 1.0e-12_WP    ! [m] Minimum horizontal source-observer distance, floors r

    CONTAINS

    ! --- Public type-bound procedures --- !
    PROCEDURE :: Psi_kzm
    PROCEDURE :: dPsi_kzm_dz
    PROCEDURE :: print_mode_summary
    PROCEDURE :: dipolar_pressure_NM

    ! --- Deferred (abstract) procedures --- !
    PROCEDURE :: init              => init_AnalyticalNormalModes
    PROCEDURE :: get_name          => get_name_AnalyticalNormalModes
    PROCEDURE :: compute_pressure  => compute_pressure_AnalyticalNormalModes
    PROCEDURE :: save_parameters   => save_parameters_AnalyticalNormalModes
    PROCEDURE :: build_save_path   => build_save_path_AnalyticalNormalModes
    PROCEDURE :: write_txt_summary => write_txt_summary_AnalyticalNormalModes

END TYPE AnalyticalNormalModes_t


CONTAINS


! ===============================================================================
! DEFERRED PROCEDURES
! ===============================================================================
SUBROUTINE init_AnalyticalNormalModes(self, sources, config)
    ! -------------------------------------------------------------------------
    ! Constructor of the Analytical Normal Modes solver. Copies the shared
    ! solver/env/turbine data from the parent AcousticSolver_t setup, reads the
    ! method-specific parameters (Nmodes, boundary planes, singularity
    ! tolerance) from the configuration, checks the waveguide depth against the
    ! turbine/environment depth, resolves the number of modes if requested
    ! automatically (propagating modes at f_max plus a near-field buffer),
    ! resolves the output save path and prints the mode summary table.
    ! -------------------------------------------------------------------------

    CLASS(AnalyticalNormalModes_t), INTENT(INOUT)              :: self       ! [-] Solver instance
    TYPE(AcousticSource_t)        , INTENT(INOUT), ALLOCATABLE :: sources    ! [-] Acoustic sources (freed on exit)
    TYPE(Config_t)                , INTENT(INOUT), ALLOCATABLE :: config     ! [-] Global configuration (freed on exit)

    ! Local variables
    REAL(WP)     :: f_max
    INTEGER(I32) :: m_prop, near_field_buffer


    ! Parent variables
    ! Parent variables (sources are moved, not copied)
    CALL MOVE_ALLOC(sources % xyz  , self % sources % xyz)
    CALL MOVE_ALLOC(sources % Force, self % sources % Force)
    CALL MOVE_ALLOC(sources % Freqs, self % sources % Freqs)
    self % sources % geometry = sources % geometry
    self % env                = config % env
    self % run_inputs         = config % run_inputs
    self % save_path          = config % solver % save_path
    self % block_size         = config % solver % block_size
    self % debug              = config % solver % debug
    self % in_cluster         = config % solver % in_cluster

    ! Child-specific variables
    SELECT TYPE (m => config % method)
    TYPE IS (Config_AnalyticalNormalModes_t)
        self % Nmodes    = m % Nmodes
        self % Upper_HBC = m % UpperHBC
        self % eps       = m % epsilon
    END SELECT

    ! Waveguide built from the environment depth (single source of truth)
    if (self % env % depth_base <= 0.0_WP) &
        error stop "AnalyticalNormalModes % init(): depth_base must be strictly positive"
    self % Lower_HBC = self % Upper_HBC - self % env % depth_base
    self % H         = self % Upper_HBC - self % Lower_HBC

    ! Number of modes
    if (self % Nmodes == DEFAULT_I) then
        f_max              = maxval(self % sources % Freqs)
        m_prop             = int(floor(0.5_WP * (4.0_WP * self % H * f_max / self % env % speed_sound_base &
                                                 + 1.0_WP)), I32)
        near_field_buffer  = 5_I32
        self % Nmodes      = max(1_I32, m_prop + near_field_buffer)
        self % auto_Nmodes = .true.
    else if (self % Nmodes <= 0_I32 .and. self % Nmodes /= DEFAULT_I) then
        error stop "AnalyticalNormalModes % init(): Nmodes must be > 0 (or -999 for automatic)"
    end if

    ! Build save_path (after Nmodes is resolved, it is part of the default filename)
    CALL self % build_save_path()

    write(*,*) ''
    write(*,*) 'Selected method: ', self % get_name()
    CALL self % print_mode_summary()

    ! Cleanup memory
    DEALLOCATE(sources, config)

END SUBROUTINE init_AnalyticalNormalModes


FUNCTION get_name_AnalyticalNormalModes(self) RESULT(name)
    ! -------------------------------------------------------------------------
    ! Returns the human-readable solver name ("Analytical Normal Modes") and
    ! stores it in self % solver_name for later use (logs, summaries, HDF5).
    ! -------------------------------------------------------------------------

    CLASS(AnalyticalNormalModes_t), INTENT(INOUT) :: self    ! [-] Solver instance
    CHARACTER(len=:)              , ALLOCATABLE   :: name    ! [-] Human-readable solver name


    name               = "Analytical Normal Modes"
    self % solver_name = trim(name)

END FUNCTION get_name_AnalyticalNormalModes


SUBROUTINE compute_pressure_AnalyticalNormalModes(self, observers, block_size, total_pressure)
    ! -------------------------------------------------------------------------
    ! Evaluates the total radiated pressure field at an arbitrary set of
    ! observer points, for every source frequency. Observers are processed
    ! either one-by-one (serial) or split into chunks distributed over OpenMP
    ! threads (in_cluster = .true.), each chunk being handed to
    ! dipolar_pressure_NM together with the real dipole sources, which are
    ! passed directly (no local copies).
    ! -------------------------------------------------------------------------

    CLASS(AnalyticalNormalModes_t), INTENT(INOUT)            :: self                   ! [-] Solver instance
    REAL(WP)                      , INTENT(IN)               :: observers(:,:)         ! [m] Observers coords (Nobs, 3)
    INTEGER(I32)                  , INTENT(IN), OPTIONAL     :: block_size             ! [-] Observers per block
    COMPLEX(WP)                   , INTENT(OUT), ALLOCATABLE :: total_pressure(:,:)    ! [Pa] Pressure (NFreqs, Nobs)

    ! Local variables
    INTEGER(I32)              :: nobs, nf, bs, n_chunks, ichunk
    INTEGER(I32)              :: istart, iend, nb, i
    COMPLEX(WP) , ALLOCATABLE :: p_block(:,:)


    ! Defaults
    bs = 256_I32 ; if (present(block_size)) bs = block_size

    ! Basic checks
    nobs = size(observers, 1)
    if (nobs <= 0_I32) then
        ALLOCATE(total_pressure(0, 0))
        return
    end if

    if (size(observers, 2) /= 3) &
        error stop "compute_pressure_AnalyticalNormalModes(): observers must have shape (N,3)"

    if (.not. ALLOCATED(self % sources % xyz) .or. .not. ALLOCATED(self % sources % Force) .or. &
        .not. ALLOCATED(self % sources % Freqs)) then
        error stop "compute_pressure_AnalyticalNormalModes(): no sources defined"
    end if

    nf = size(self % sources % Freqs)

    if (size(self % sources % Force, 1) /= nf .or. &
        size(self % sources % Force, 2) /= size(self % sources % xyz, 1) .or. &
        size(self % sources % Force, 3) /= 3) &
        error stop "compute_pressure_AnalyticalNormalModes(): sources % Force must have shape (NFreqs, Nsources, 3)"

    ! Never fewer blocks than threads, so small observer sets still use every thread
    bs = max(1_I32, min(bs, (nobs + omp_get_max_threads() - 1_I32) / omp_get_max_threads()))

    ALLOCATE(total_pressure(nf, nobs))
    total_pressure = (0.0_WP, 0.0_WP)

    if (.not. self % in_cluster) bs = 1_I32
    if (self % in_cluster) then
        n_chunks = (nobs + bs - 1_I32) / bs

        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP PRIVATE(ichunk, istart, iend, nb, p_block) &
        !$OMP SHARED(n_chunks, nobs, total_pressure, self)
        do ichunk = 0, n_chunks - 1
            istart = ichunk * bs + 1_I32
            iend   = min(istart + bs - 1_I32, nobs)
            nb     = iend - istart + 1_I32

            ! Call kernel for this chunk (p_block will be allocated in the routine)
            CALL self % dipolar_pressure_NM(observers(istart:iend, :), self % sources % Freqs, &
                                            self % sources % xyz, self % sources % Force, p_block)

            ! Accumulate results into the global array (non-overlapping slices)
            total_pressure(:, istart:iend) = total_pressure(:, istart:iend) + p_block(:, 1:nb)

            if (ALLOCATED(p_block)) DEALLOCATE(p_block)
        end do
        !$OMP END PARALLEL DO

    else    ! .not. self % in_cluster

        ! Serial: loop observers (keep 2D slice for compatibility)
        do i = 1, nobs
            CALL self % dipolar_pressure_NM(observers(i:i, :), self % sources % Freqs, &
                                            self % sources % xyz, self % sources % Force, p_block)

            total_pressure(:, i) = total_pressure(:, i) + p_block(:, 1)

            if (ALLOCATED(p_block)) DEALLOCATE(p_block)

            ! Print progress every 100 observers
            if (mod(i, 100_I32) == 0_I32) then
                write(*, '(A,I0,A,I0)') '  Progress: ', i, '/', nobs
                flush(6)
            end if
        end do

    end if

END SUBROUTINE compute_pressure_AnalyticalNormalModes


SUBROUTINE save_parameters_AnalyticalNormalModes(self, turbine, farm)
    ! -------------------------------------------------------------------------
    ! Saves the solver, environment and generic turbine parameters to the
    ! HDF5 output file, delegates the turbine-specific parameters to
    ! turbine % save_specific_params, and triggers the human-readable .dat
    ! text summary via write_txt_summary.
    ! -------------------------------------------------------------------------

    CLASS(AnalyticalNormalModes_t), INTENT(IN)           :: self       ! [-] Solver instance
    CLASS(WindTurbine_t)          , INTENT(IN), OPTIONAL :: turbine    ! [-] Turbine instance (optional)
    TYPE(WindFarm_t)              , INTENT(IN), OPTIONAL :: farm       ! [-] Wind farm instance (optional)


    ! Solver / environment / generic turbine data -> HDF5
    CALL save_variable(self % save_path, "Nmodes"   , self % Nmodes)
    CALL save_variable(self % save_path, "H"        , self % H)
    CALL save_variable(self % save_path, "c_wat"    , self % env % speed_sound_base)
    CALL save_variable(self % save_path, "rho_wat"  , self % env % density_base)
    CALL save_variable(self % save_path, "Upper_HBC", self % Upper_HBC)
    CALL save_variable(self % save_path, "Lower_HBC", self % Lower_HBC)
    CALL save_variable(self % save_path, "p_ref"    , self % env % reference_pressure)

    ! Emitter data -> HDF5
    CALL save_emitter_params(self % save_path, turbine, farm)

    ! Human-readable .dat summary, delegated to the solver-specific implementation
    CALL self % write_txt_summary(turbine, farm)

END SUBROUTINE save_parameters_AnalyticalNormalModes


SUBROUTINE build_save_path_AnalyticalNormalModes(self)
    ! -------------------------------------------------------------------------
    ! Resolves the final HDF5 output path through the shared
    ! resolve_save_path, with this solver's method tag ("_ANM" followed by the
    ! number of modes, so it must run after Nmodes is resolved). Folder and
    ! file name labels come from the emitter.
    ! -------------------------------------------------------------------------

    CLASS(AnalyticalNormalModes_t), INTENT(INOUT) :: self    ! [-] Solver instance

    ! Local variables
    CHARACTER(len=16) :: modes_str


    write(modes_str, '(I0)') self % Nmodes
    CALL self % resolve_save_path("_ANM" // trim(modes_str))

END SUBROUTINE build_save_path_AnalyticalNormalModes


SUBROUTINE write_txt_summary_AnalyticalNormalModes(self, turbine, farm)
    ! -------------------------------------------------------------------------
    ! Writes a plain-text (.dat) summary next to the HDF5 output, mirroring
    ! the saved turbine, environment and Analytical Normal Modes solver
    ! parameters in a human-readable form for quick inspection without
    ! opening the HDF5 file.
    ! -------------------------------------------------------------------------

    CLASS(AnalyticalNormalModes_t), INTENT(IN)           :: self       ! [-] Solver instance
    CLASS(WindTurbine_t)          , INTENT(IN), OPTIONAL :: turbine    ! [-] Turbine instance (optional)
    TYPE(WindFarm_t)              , INTENT(IN), OPTIONAL :: farm       ! [-] Wind farm instance (optional)

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

    CALL write_emitter_txt(u, turbine, farm)

    write(u, '(A)') repeat('-', 70)
    write(u, '(A)') 'ENVIRONMENT'
    write(u, '(A)') repeat('-', 70)
    write(u, '(A,F8.2,A)')   'Speed of sound   : ', self % env % speed_sound_base, ' m/s'
    write(u, '(A,F8.2,A)')   'Density          : ', self % env % density_base, ' kg/m^3'
    write(u, '(A,ES10.2,A)') 'Reference press. : ', self % env % reference_pressure, ' Pa'
    write(u, '(A,F8.2,A)')   'Depth            : ', self % env % depth_base, ' m'
    write(u, *)

    write(u, '(A)') repeat('-', 70)
    write(u, '(A)') 'SOLVER: ' // trim(self % solver_name)
    write(u, '(A)') repeat('-', 70)
    write(u, '(A,I0)')       'N modes      : ', self % Nmodes
    write(u, '(A,L1)')       'Auto N modes : ', self % auto_Nmodes
    write(u, '(A,F8.2,A)')   'Height H     : ', self % H, ' m'
    write(u, '(A,F8.2,A)')   'Upper HBC    : ', self % Upper_HBC, ' m'
    write(u, '(A,F8.2,A)')   'Lower HBC    : ', self % Lower_HBC, ' m'
    write(u, '(A,ES10.2,A)') 'Epsilon      : ', self % eps, ' m'
    write(u, '(A,I0)')       'Block size   : ', self % block_size
    write(u, '(A,L1)')       'In cluster   : ', self % in_cluster
    write(u, *)

    close(u)

    write(*,*) '   ---> Text summary saved at ', dat_path

END SUBROUTINE write_txt_summary_AnalyticalNormalModes


! ===============================================================================
! HELPERS
! ===============================================================================
SUBROUTINE print_mode_summary(self)
    ! -------------------------------------------------------------------------
    ! Prints a formatted table of the retained normal modes: vertical
    ! wavenumber and wavelength, cutoff frequency, and horizontal wavenumber
    ! and wavelength evaluated at f_max, flagging each mode as propagating or
    ! evanescent. For more than 20 modes only the first 10 and last 5 rows
    ! are shown.
    ! -------------------------------------------------------------------------

    CLASS(AnalyticalNormalModes_t), INTENT(IN) :: self    ! [-] Solver instance

    ! Local variables
    REAL(WP)                       :: fmin, fmax, k_medium, k_zm, lambda_zm, fc, kr_sq, krm_val
    INTEGER(I32)                   :: i, m_idx, line_width
    CHARACTER(len=20)              :: kzm_str, krm_str, lambda_str, status
    CHARACTER(len=16)              :: re_str, im_str
    CHARACTER(len=*) , PARAMETER   :: title = "ANALYTICAL NORMAL MODES SETUP"
    INTEGER(I32)     , ALLOCATABLE :: modes_to_show(:)


    fmin       = minval(self % sources % Freqs)
    fmax       = maxval(self % sources % Freqs)
    k_medium   = 2.0_WP * PI * fmax / self % env % speed_sound_base
    line_width = 115_I32

    write(*, '(A)') ''
    write(*, '(A)') repeat('=', line_width)
    write(*, '(A)') repeat(' ', (line_width - len(title)) / 2_I32) // title
    write(*, '(A)') repeat('=', line_width)
    write(*, '(A,F10.2,A)')         '  Water Depth (H):      ', self % H, ' m'
    write(*, '(A,F10.2,A)')         '  Speed of Sound (c):   ', self % env % speed_sound_base, ' m/s'
    write(*, '(A,F10.2,A,F10.2,A)') '  Frequency Band:       ', fmin, ' Hz - ', fmax, ' Hz'
    write(*, '(A,I10)')             '  Retained Modes (m):   ', self % Nmodes
    write(*, '(A,F10.2,A)')         '  * Note: k_rm and Horiz. Wavelength are evaluated at f_max (', fmax, ' Hz)'
    write(*, '(A)') repeat('-', line_width)
    write(*, '(2X,A5,A,A18,A,A12,A,A13,A,A20,A,A13,A,A)') &
        'Mode ', ' | ', 'k_zm [rad/m]      ', ' | ', 'Vert. WL [m]', ' | ', 'f_cutoff [Hz]', ' | ', &
        'k_rm [rad/m]        ', ' | ', 'Horiz. WL [m]', ' | ', 'Status'
    write(*, '(A)') repeat('-', line_width)

    if (self % Nmodes <= 20_I32) then
        modes_to_show = [(i, i = 1, self % Nmodes)]
    else
        modes_to_show = [(i, i = 1, 10), -1_I32, (i, i = self % Nmodes - 4, self % Nmodes)]
    end if

    do i = 1, size(modes_to_show)
        m_idx = modes_to_show(i)
        if (m_idx == -1_I32) then
            write(*, '(2X,A5,A,A18,A,A12,A,A13,A,A20,A,A13,A,A)') &
                '...  ', ' | ', '...', ' | ', '...', ' | ', '...', ' | ', '...', ' | ', '...', ' | ', '...'
            cycle
        end if

        ! Vertical properties (frequency independent for p = 0 top / dp/dz = 0 bottom)
        k_zm      = (2.0_WP * real(m_idx, WP) - 1.0_WP) * PI / (2.0_WP * self % H)
        lambda_zm = 4.0_WP * self % H / (2.0_WP * real(m_idx, WP) - 1.0_WP)
        fc        = (2.0_WP * real(m_idx, WP) - 1.0_WP) * self % env % speed_sound_base / (4.0_WP * self % H)

        write(re_str, '(F12.4)') k_zm
        kzm_str = trim(adjustl(re_str)) // '+0.0000j'

        ! Horizontal properties (frequency dependent, evaluated at f_max)
        kr_sq = k_medium**2 - k_zm**2
        if (kr_sq >= 0.0_WP) then
            krm_val = sqrt(kr_sq)
            write(re_str, '(F12.4)') krm_val
            krm_str = trim(adjustl(re_str)) // '+0.0000j'
            write(lambda_str, '(F12.2)') 2.0_WP * PI / krm_val
            lambda_str = adjustl(lambda_str)
            status     = 'Propagating'
        else
            krm_val = sqrt(abs(kr_sq))
            write(im_str, '(F12.4)') krm_val
            krm_str    = '0.0000+' // trim(adjustl(im_str)) // 'j'
            lambda_str = 'N/A (Decays)'
            status     = 'Evanescent'
        end if

        write(*, '(2X,I5,A,A18,A,F12.2,A,F13.2,A,A20,A,A13,A,A)') &
            m_idx, ' | ', kzm_str, ' | ', lambda_zm, ' | ', fc, ' | ', krm_str, ' | ', lambda_str, ' | ', trim(status)
    end do

    write(*, '(A)') repeat('=', line_width)
    write(*, '(A)') ''

END SUBROUTINE print_mode_summary


! ===============================================================================
! METHOD (CORE PHYSICS)
! ===============================================================================
FUNCTION Psi_kzm(self, z, m, kzm) RESULT(psi)
    ! -------------------------------------------------------------------------
    ! Evaluates the vertical modal shape Psi_m(z) = sqrt(2*rho/H) * sin(k_zm*(z - Upper_HBC))
    ! at the given depths. The vertical wavenumber is either given directly
    ! (kzm) or computed from the mode index m as k_zm = (2m-1)*pi/(2H). The
    ! shape satisfies Psi_m(Upper_HBC) = 0 = 0 (pressure-release surface).
    ! -------------------------------------------------------------------------

    CLASS(AnalyticalNormalModes_t), INTENT(IN)           :: self      ! [-] Solver instance
    REAL(WP)                      , INTENT(IN)           :: z(:)      ! [m] Depths at which the modal shape is evaluated
    INTEGER(I32)                  , INTENT(IN), OPTIONAL :: m         ! [-] Mode index (used if kzm is not given)
    REAL(WP)                      , INTENT(IN), OPTIONAL :: kzm       ! [1/m] Vertical wavenumber of the mode
    REAL(WP)                      , ALLOCATABLE          :: psi(:)    ! [-] Modal shape Psi_m(z)


    ALLOCATE(psi(size(z)))

    if (present(m) .and. .not. present(kzm)) then
        psi = sqrt(2.0_WP * self % env % density_base / self % H) * &
              sin((2.0_WP * real(m, WP) - 1.0_WP) * PI * (z - self % Upper_HBC) / (2.0_WP * self % H))
    else if (present(kzm)) then
        psi = sqrt(2.0_WP * self % env % density_base / self % H) * sin(kzm * (z - self % Upper_HBC))
    else
        error stop "Psi_kzm: Either m or kzm must be provided."
    end if

END FUNCTION Psi_kzm


FUNCTION dPsi_kzm_dz(self, z, m, kzm) RESULT(dpsi)
    ! -------------------------------------------------------------------------
    ! Evaluates the vertical derivative of the modal shape,
    ! dPsi_m/dz = sqrt(2*rho/H) * k_zm * cos(k_zm*z), at the given depths. The
    ! vertical wavenumber is either given directly (kzm) or computed from the
    ! mode index m as k_zm = (2m-1)*pi/(2H). The derivative satisfies
    ! dPsi_m/dz(Lower_HBC=-H) = 0 (rigid seabed).
    ! -------------------------------------------------------------------------

    CLASS(AnalyticalNormalModes_t), INTENT(IN)           :: self       ! [-] Solver instance
    REAL(WP)                      , INTENT(IN)           :: z(:)       ! [m] Depths at which the derivative is evaluated
    INTEGER(I32)                  , INTENT(IN), OPTIONAL :: m          ! [-] Mode index (used if kzm is not given)
    REAL(WP)                      , INTENT(IN), OPTIONAL :: kzm        ! [1/m] Vertical wavenumber of the mode
    REAL(WP)                      , ALLOCATABLE          :: dpsi(:)    ! [1/m] Modal shape derivative dPsi_m/dz

    ! Local variables
    REAL(WP) :: k_zm


    if (present(kzm)) then
        k_zm = kzm
    else if (present(m)) then
        k_zm = (2.0_WP * real(m, WP) - 1.0_WP) * PI / (2.0_WP * self % H)
    else
        error stop "dPsi_kzm_dz: Either m or kzm must be provided."
    end if

    ALLOCATE(dpsi(size(z)))
    dpsi = sqrt(2.0_WP * self % env % density_base / self % H) * k_zm * cos(k_zm * (z - self % Upper_HBC))

END FUNCTION dPsi_kzm_dz


SUBROUTINE dipolar_pressure_NM(self, obs_pos, Freqs, nodes_pos, Force, p_out)
    ! -------------------------------------------------------------------------
    ! Evaluates the normal-mode dipole Green's function for a block of
    ! observers, summed over every real source and every retained mode, at
    ! every frequency. The vertical modal shapes
    !     Psi_m(z)   = sqrt(2*rho/H) * sin(k_zm*z)
    !     dPsi_m/dz  = sqrt(2*rho/H) * k_zm * cos(k_zm*z)
    ! (Psi_m(0) = 0, dPsi_m/dz(-H) = 0) are frequency independent, so they
    ! are tabulated once per block for observers and sources. Horizontal
    ! distances are also precomputed (floored at eps). For each (f, m) the
    ! horizontal wavenumber k_rm = sqrt(k^2 - k_zm^2) is either real
    ! (propagating) or pure imaginary (evanescent), and each source-observer
    ! pair contributes
    !     Psi_m(z) * [ Fz*dPsi_m(zs)/rho * H0(k_rm*r)
    !                + Psi_m(zs)*k_rm/rho * H1(k_rm*r) * (Fx*dx + Fy*dy)/r ]
    ! with H0/H1 Hankel functions of the first kind. The global -i/4
    ! coefficient is applied at the end.
    ! -------------------------------------------------------------------------

    CLASS(AnalyticalNormalModes_t), INTENT(IN)               :: self              ! [-] Solver instance
    REAL(WP)                      , INTENT(IN)               :: obs_pos(:,:)      ! [m] Observers coords (Nob, 3)
    REAL(WP)                      , INTENT(IN)               :: Freqs(:)          ! [Hz] Frequency vector (NFreqs)
    REAL(WP)                      , INTENT(IN)               :: nodes_pos(:,:)    ! [m] Source coords (Nnodes, 3)
    COMPLEX(WP)                   , INTENT(IN)               :: Force(:,:,:)      ! [N] Forces (NFreqs, Nnodes, 3)
    COMPLEX(WP)                   , INTENT(OUT), ALLOCATABLE :: p_out(:,:)        ! [Pa] Pressure (NFreqs, Nob)

    ! Local variables
    INTEGER(I32)              :: NFreqs, Nnodes, Nob, Nm
    INTEGER(I32)              :: j, f, k, im
    REAL(WP)                  :: epsv, rho, A_norm, kappa
    COMPLEX(WP)               :: k_rm, H0, H1, fx, fy, fz, coef_z, coef_xy, F_dot_rhat
    REAL(WP)    , ALLOCATABLE :: dx(:,:), dy(:,:), r(:,:), r_inv(:,:)
    REAL(WP)    , ALLOCATABLE :: k_wave(:), k_zm(:)
    REAL(WP)    , ALLOCATABLE :: psi_z(:,:), psi_zs(:,:), dpsi_zs(:,:)


    ! Dimensions
    NFreqs = size(Freqs)
    Nnodes = size(nodes_pos, 1)
    Nob    = size(obs_pos, 1)
    Nm     = self % Nmodes

    epsv   = self % eps
    rho    = self % env % density_base
    A_norm = sqrt(2.0_WP * rho / self % H)

    ALLOCATE(p_out(NFreqs, Nob))
    p_out  = (0.0_WP, 0.0_WP)
    k_wave = 2.0_WP * PI * Freqs / self % env % speed_sound_base

    ! Horizontal source-observer geometry (frequency and node independent)
    ALLOCATE(dx(Nob, Nnodes), dy(Nob, Nnodes), r(Nob, Nnodes), r_inv(Nob, Nnodes))
    do j = 1, Nnodes
        !$OMP SIMD
        do k = 1, Nob
            dx(k, j) = obs_pos(k, 1) - nodes_pos(j, 1)
            dy(k, j) = obs_pos(k, 2) - nodes_pos(j, 2)
            r(k, j)  = sqrt(dx(k, j) * dx(k, j) + dy(k, j) * dy(k, j))

            if (r(k, j) < epsv) r(k, j) = epsv
            r_inv(k, j) = 1.0_WP / r(k, j)
        end do
    end do

    ! Vertical modal shapes (frequency independent)
    ALLOCATE(k_zm(Nm), psi_z(Nob, Nm), psi_zs(Nnodes, Nm), dpsi_zs(Nnodes, Nm))
    do im = 1, Nm
        k_zm(im)       = (2.0_WP * real(im, WP) - 1.0_WP) * PI / (2.0_WP * self % H)
        psi_z(:, im)   = self % Psi_kzm(z=obs_pos(:, 3), kzm=k_zm(im))
        psi_zs(:, im)  = self % Psi_kzm(z=nodes_pos(:, 3), kzm=k_zm(im))
        dpsi_zs(:, im) = self % dPsi_kzm_dz(z=nodes_pos(:, 3), kzm=k_zm(im))
    end do

    ! Loop over frequencies
    do f = 1, NFreqs
        ! Normal modes superposition
        do im = 1, Nm

            ! Horizontal wavenumber: real (propagating) or pure imaginary (evanescent)
            if (k_wave(f) >= k_zm(im)) then
                k_rm = cmplx(sqrt(k_wave(f)**2 - k_zm(im)**2), 0.0_WP, kind=WP)
            else
                kappa = sqrt(k_zm(im)**2 - k_wave(f)**2)
                k_rm  = I1 * kappa
            end if

            ! Loop over sources and observers
            do j = 1, Nnodes
                fx = Force(f, j, 1)
                fy = Force(f, j, 2)
                fz = Force(f, j, 3)

                coef_z  = fz * dpsi_zs(j, im) / rho     ! Vertical component (dipolar depth effect)
                coef_xy = psi_zs(j, im) * k_rm / rho    ! Horizontal component (dipolar directivity)

                do k = 1, Nob
                    CALL compute_hankel_complex(k_rm, r(k, j), H0, H1)

                    F_dot_rhat  = (fx * dx(k, j) + fy * dy(k, j)) * r_inv(k, j)
                    p_out(f, k) = p_out(f, k) + psi_z(k, im) * (coef_z * H0 + coef_xy * H1 * F_dot_rhat)
                end do
            end do

        end do
    end do

    ! Global coefficient
    p_out = -0.25_WP * I1 * p_out
    DEALLOCATE(dx, dy, r, r_inv, k_zm, psi_z, psi_zs, dpsi_zs)

END SUBROUTINE dipolar_pressure_NM

END MODULE AnalyticalNormalModes
