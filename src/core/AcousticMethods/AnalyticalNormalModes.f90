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
USE MathUtils, ONLY: compute_hankel_complex
USE WindTurbine, ONLY: WindTurbine_t, save_turbine_params
USE AcousticSolver, ONLY: AcousticSolver_t

IMPLICIT NONE

PRIVATE
PUBLIC:: AnalyticalNormalModes_t


! ===============================================================================
! TYPE: AnalyticalNormalModes_t
! ===============================================================================
! Concrete acoustic solver implementing the analytical normal-mode expansion
! of the dipole field in an ideal waveguide of depth H = Upper_HBC - Lower_HBC
! (pressure-release top, rigid bottom). Each mode m has a fixed vertical
! wavenumber k_zm = (2m-1)*pi/(2H) and a frequency-dependent horizontal
! wavenumber k_rm = sqrt(k^2 - k_zm^2).
! ===============================================================================
TYPE, EXTENDS(AcousticSolver_t) :: AnalyticalNormalModes_t

    INTEGER(I32) :: Nmodes           = DEFAULT_I    ! [-] Number of normal modes retained in the expansion
    LOGICAL      :: auto_Nmodes      = .false.      ! [-] Whether Nmodes was computed automatically from f_max
    REAL(WP)     :: H                = 30.0_WP      ! [m] Waveguide height (distance between both boundaries)
    REAL(WP)     :: Upper_HBC        = 0.0_WP       ! [m] z-coordinate of the upper boundary (free surface, p = 0)
    REAL(WP)     :: Lower_HBC        = -30.0_WP     ! [m] z-coordinate of the lower boundary (seabed, dp/dz = 0)
    REAL(WP)     :: eps              = 1e-12_WP     ! [m] Minimum horizontal source-observer distance, floors r

    CONTAINS
    ! --- Public type-bound procedures --- !
    PROCEDURE :: Psi_kzm
    PROCEDURE :: dPsi_kzm_dz
    PROCEDURE :: print_mode_summary
    PROCEDURE :: dipolar_pressure_NM

    ! --- Deferred(abstract) procedures --- !
    PROCEDURE :: init               => init_AnalyticalNormalModes
    PROCEDURE :: get_name           => get_name_AnalyticalNormalModes
    PROCEDURE :: compute_pressure   => compute_pressure_AnalyticalNormalModes
    PROCEDURE :: save_parameters    => save_parameters_AnalyticalNormalModes
    PROCEDURE :: build_save_path    => build_save_path_AnalyticalNormalModes
    PROCEDURE :: write_txt_summary  => write_txt_summary_AnalyticalNormalModes

END TYPE AnalyticalNormalModes_t


CONTAINS


! ===============================================================================
! DEFERRED (abstract AcousticSolver_t) PROCEDURES
! ===============================================================================
! -------------------------------------------------------------------------
! Subroutine: init_AnalyticalNormalModes
! Constructor of the Analytical Normal Modes solver. Copies the shared
! solver/env/turbine data from the parent AcousticSolver_t setup, reads the
! method-specific parameters (Nmodes, boundary planes, singularity
! tolerance) from the configuration, checks the waveguide depth against the
! turbine/environment depth, resolves the number of modes if requested
! automatically (propagating modes at f_max plus a near-field buffer),
! resolves the output save path and prints the mode summary table.
! -------------------------------------------------------------------------
SUBROUTINE init_AnalyticalNormalModes(self, sources, config)
    CLASS(AnalyticalNormalModes_t), INTENT(INOUT)              :: self
    TYPE(AcousticSource_t)        , INTENT(INOUT), ALLOCATABLE :: sources
    TYPE(Config_t)                , INTENT(INOUT), ALLOCATABLE :: config

    ! Local variables
    REAL(WP)     :: f_max
    INTEGER(I32) :: m_prop, near_field_buffer


    ! Parent variables
    self % sources      = sources
    self % env          = config % env
    self % run_inputs   = config % run_inputs
    self % save_path    = config % solver % save_path
    self % block_size   = config % solver % block_size
    self % debug        = config % solver % debug
    self % in_cluster   = config % solver % in_cluster
    self % turbine_data = config % turbine

    ! Son specific variables
    SELECT TYPE(m => config % method)
    TYPE IS(Config_AnalyticalNormalModes_t)
        self % Nmodes    = m % Nmodes
        self % Upper_HBC = m % UpperHBC
        self % Lower_HBC = m % LowerHBC
        self % eps       = m % epsilon
    END SELECT

    ! Waveguide height and depth consistency
    self % H = self % Upper_HBC - self % Lower_HBC
    if (self % H <= 0.0_WP) error stop "AnalyticalNormalModes % init(): Upper_HBC must be strictly greater than Lower_HBC"
    if ( abs(self % H - self % env % depth_base) > 1e-6 .or. &
    abs(self % H - self % turbine_data % Depth) > 1e-6 .or. &
    abs(self % env % depth_base - self % turbine_data % Depth) > 1e-6 ) then
        error stop "AnalyticalNormalModes % init(): turbine, environmental and solver depth must be equal"
    end if

    ! Number of modes
    if (self % Nmodes == DEFAULT_I) then
        f_max = maxval(self % sources % Freqs)
        m_prop = int(floor(0.5_WP * (4.0_WP * self % H * f_max / self % env % speed_sound_base + 1.0_WP)), I32)
        near_field_buffer = 5
        self % Nmodes = max(1_I32, m_prop + near_field_buffer)
        self % auto_Nmodes = .true.
    else if (self % Nmodes <= 0 .and. self % Nmodes /= DEFAULT_I) then
        error stop "AnalyticalNormalModes % init(): Nmodes must be > 0 (or -999 for automatic)"
    end if

    ! Build save_path (after Nmodes is resolved, it is part of the default filename)
    CALL self % build_save_path()

    print*, ''
    print*, 'Selected method: ', self % get_name()
    CALL self % print_mode_summary()

    ! Cleanup memory
    DEALLOCATE(sources, config)

END SUBROUTINE init_AnalyticalNormalModes


! -------------------------------------------------------------------------
! Function: get_name_AnalyticalNormalModes
! Returns the human-readable solver name ("Analytical Normal Modes") and
! stores it in self % solver_name for later use (logs, summaries, HDF5).
! -------------------------------------------------------------------------
FUNCTION get_name_AnalyticalNormalModes(self) RESULT(name)
    CLASS(AnalyticalNormalModes_t), INTENT(INOUT) :: self
    CHARACTER(len=:)              , ALLOCATABLE   :: name
    name = "Analytical Normal Modes"; self % solver_name = trim(name)
END FUNCTION get_name_AnalyticalNormalModes


! -------------------------------------------------------------------------
! Subroutine: compute_pressure_AnalyticalNormalModes
! Evaluates the total radiated pressure field at an arbitrary set of
! observer points, for every source frequency. Observers are processed
! either one-by-one (serial) or split into chunks distributed over OpenMP
! threads (in_cluster = .true.), each chunk being handed to
! dipolar_pressure_NM together with the real dipole sources.
! -------------------------------------------------------------------------
SUBROUTINE compute_pressure_AnalyticalNormalModes(self, observers, block_size, total_pressure)
    CLASS(AnalyticalNormalModes_t), INTENT(INOUT)        :: self
    REAL(WP)                      , INTENT(IN)           :: observers(:,:)        ! [m] Observers coords (Nobs,3)
    INTEGER(I32)                  , INTENT(IN), OPTIONAL :: block_size            ! [-] Number of observers to process simultaneously
    COMPLEX(WP), ALLOCATABLE      , INTENT(OUT)          :: total_pressure(:,:)   ! [Pa] pressure array(NFreqs, Nobs)

    ! Local variables
    INTEGER(I32)             :: nobs, nf, bs, n_chunks, ichunk
    INTEGER(I32)             :: istart, iend, nb, i
    REAL(WP), ALLOCATABLE    :: Freqs(:), xyz(:,:)
    COMPLEX(WP), ALLOCATABLE :: Force(:,:,:)
    COMPLEX(WP), ALLOCATABLE :: p_block(:,:)


    ! Defaults
    bs = 256_I32; if (PRESENT(block_size)) bs = block_size

    ! Basic checks
    nobs = size(observers,1)
    if (nobs <= 0) then
        ALLOCATE(total_pressure(0,0))
        return
    end if

    if (SIZE(observers, 2) /= 3) error stop "compute_pressure_AnalyticalNormalModes(): observers must have shape (N,3)"

    if (.not. ALLOCATED(self % sources % xyz) .OR. .not. ALLOCATED(self % sources % Force) .OR. &
        .not. ALLOCATED(self % sources % Freqs)) then
        error stop "compute_pressure_AnalyticalNormalModes(): no sources defined"
    end if

    Freqs = self % sources % Freqs
    Force = self % sources % Force
    xyz   = self % sources % xyz
    nf    = size(Freqs)

    if (SIZE(Force, 1) /= nf .or. SIZE(Force, 2) /= SIZE(xyz, 1) .or. SIZE(Force, 3) /= 3) &
        error stop "compute_pressure_AnalyticalNormalModes(): sources % Force must have shape (NFreqs, Nsources, 3)"

    ALLOCATE(total_pressure(nf, nobs))
    total_pressure = (0.0_WP, 0.0_WP)

    
    if (.not. self % in_cluster) bs = 1
    if (self % in_cluster) then

        n_chunks = (nobs + bs - 1) / bs

        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP PRIVATE(ichunk, istart, iend, nb, p_block) &
        !$OMP SHARED(xyz, FOrce, Freqs, n_chunks, nobs, total_pressure, self)
        do ichunk = 0, n_chunks -1

            istart = ichunk * bs + 1
            iend   = min(istart + bs -1, nobs)
            nb     = iend - istart + 1

            ! Call kernel for this chunk (p_block will be allocated in the routine)
            CALL self % dipolar_pressure_NM(observers(istart:iend, :), Freqs, xyz, Force, p_block)

            ! Acumulate results into the global array (non-overlapping slices)
            total_pressure(:, istart:iend) = total_pressure(:, istart:iend) + p_block(:, 1:nb)

            if (ALLOCATED(p_block)) DEALLOCATE(p_block)
        end do
        !$OMP END PARALLEL DO

    else ! (.not. self % in_cluster)

        ! Serial: loop observers (keep 2D slice for compatibility)
        do i =1, nobs

            CALL self % dipolar_pressure_NM(observers(i:i, :), Freqs, xyz, Force, p_block)
            total_pressure(:, i) = total_pressure(:, i) + p_block(:, 1)

            if (ALLOCATED(p_block)) DEALLOCATE(p_block)

            ! Print progress every 100 observers
            if (mod(i, 100) == 0) then
                print '(A,I0,A,I0)', '  Progress: ', i, '/', nobs
                flush(6)
            end if
        end do

    end if

END SUBROUTINE compute_pressure_AnalyticalNormalModes


! -------------------------------------------------------------------------
! Subroutine: save_parameters_AnalyticalNormalModes
! Saves the solver, environment and generic turbine parameters to the
! HDF5 output file, delegates the turbine-specific parameters to
! turbine % save_specific_params, and triggers the human-readable .dat
! text summary via write_txt_summary.
! -------------------------------------------------------------------------
SUBROUTINE save_parameters_AnalyticalNormalModes(self, turbine)

    CLASS(AnalyticalNormalModes_t), INTENT(IN)                        :: self
    CLASS(WindTurbine_t)          , INTENT(IN), OPTIONAL, ALLOCATABLE :: turbine

    ! Solver / environment / generic turbine data -> HDF5
    CALL save_variable(self % save_path, "Nmodes"   , self % Nmodes)
    CALL save_variable(self % save_path, "H"        , self % H)
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

END SUBROUTINE save_parameters_AnalyticalNormalModes


! -------------------------------------------------------------------------
! Subroutine: build_save_path_AnalyticalNormalModes
! Resolves the final HDF5 output path. If the configured save_path is
! "None"/"Default", it auto-generates a descriptive filename (turbine
! type, number of modes, wind speed) under a default output folder;
! otherwise it just appends the .hdf5 extension. Ensures the destination
! folder exists.
! -------------------------------------------------------------------------
SUBROUTINE build_save_path_AnalyticalNormalModes(self)
    CLASS(AnalyticalNormalModes_t), INTENT(INOUT) :: self

    ! Local variables
    CHARACTER(len=512) :: path_upper, filename, default_parent_path
    CHARACTER(len=16)  :: method_suffix, ws_str, modes_str
    CHARACTER(len=2)   :: turb_tag

    path_upper = to_upper(trim(self % save_path))
    path_upper = trim(path_upper)
    default_parent_path = "./turbine_acoustic_data/"

    method_suffix = "_ANM"

    if (path_upper == 'NONE' .or. path_upper == 'DEFAULT') then

        if (self % turbine_data % is_monopile) then
            turb_tag = 'mn'
        else
            turb_tag = 'fl'
        end if

        write(modes_str,'(I0)') self % Nmodes
        write(ws_str, '(F6.1)') self % turbine_data % WindSpeed
        ws_str = adjustl(ws_str)

        filename = "plot_"             // trim(turb_tag)  // &
                   trim(method_suffix) // trim(modes_str) // &
                   "_ws"               // trim(ws_str)

        self % save_path = trim(default_parent_path) // trim(filename) // ".hdf5"

    else
        self % save_path = trim(self % save_path) // ".hdf5"

    end if

    ! Make sure the destination folder exists before anything tries to write into it
    CALL ensure_parent_dir(trim(self % save_path))

END SUBROUTINE build_save_path_AnalyticalNormalModes


! -------------------------------------------------------------------------
! Subroutine: write_txt_summary_AnalyticalNormalModes
! Writes a plain-text (.dat) summary next to the HDF5 output, mirroring
! the saved turbine, environment and Analytical Normal Modes solver
! parameters in a human-readable form for quick inspection without
! opening the HDF5 file.
! -------------------------------------------------------------------------
SUBROUTINE write_txt_summary_AnalyticalNormalModes(self, turbine)
    CLASS(AnalyticalNormalModes_t), INTENT(IN)                        :: self
    CLASS(WindTurbine_t)          , INTENT(IN), OPTIONAL, ALLOCATABLE :: turbine

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: dat_path
    INTEGER(I32) :: u

    dat_path = swap_extension(trim(self % save_path), ".dat")

    open(NEWUNIT=u, FILE=dat_path, STATUS='replace', ACTION='write')

    write(u,'(A)') repeat('=', 70)
    write(u,'(A)') 'OFF_COUSTICS SIMULATION SUMMARY'
    write(u,'(A)') repeat('=', 70)
    write(u,'(A)') 'HDF5 file: ' // trim(self % save_path)
    write(u,*)

    write(u,'(A)') repeat('-', 70)
    write(u,'(A)') 'TURBINE'
    write(u,'(A)') repeat('-', 70)
    write(u,'(A,A)')      'Case type   : ', trim(self % turbine_data % case_type)
    write(u,'(A,L1)')     'Is monopile : ', self % turbine_data % is_monopile
    write(u,'(A,F8.2,A)') 'Wind speed  : ', self % turbine_data % WindSpeed, ' m/s'
    write(u,'(A,F8.2,A)') 'Wind dir    : ', self % turbine_data % WindDir, ' deg'
    write(u,'(A,F8.2,A)') 'Water depth : ', self % turbine_data % Depth, ' m'
    write(u,'(A,2F10.3)') 'Axis pos    : ', self % turbine_data % AxisPos
    write(u,'(A,2F10.3)') 'Bari pos    : ', self % turbine_data % BariPos
    write(u,'(A,I0)')     'N members   : ', self % turbine_data % Nmembers
    write(u,'(A,I0)')     'N nodes     : ', self % turbine_data % Nnodes
    write(u,*)

    if (PRESENT(turbine)) then
        if (ALLOCATED(turbine)) then
            write(u,'(A)') repeat('-', 70)
            write(u,'(A)') 'TURBINE - SPECIFIC PARAMETERS'
            write(u,'(A)') repeat('-', 70)
            CALL turbine % write_specific_params_txt(u)
            write(u,*)
        end if
    end if

    write(u,'(A)') repeat('-', 70)
    write(u,'(A)') 'ENVIRONMENT'
    write(u,'(A)') repeat('-', 70)
    write(u,'(A,F8.2,A)')  'Speed of sound   : ', self % env % speed_sound_base, ' m/s'
    write(u,'(A,F8.2,A)')  'Density          : ', self % env % density_base, ' kg/m^3'
    write(u,'(A,ES10.2,A)') 'Reference press. : ', self % env % reference_pressure, ' Pa'
    write(u,'(A,F8.2,A)')  'Depth            : ', self % env % depth_base, ' m'
    write(u,*)

    write(u,'(A)') repeat('-', 70)
    write(u,'(A)') 'SOLVER: ' // trim(self % solver_name)
    write(u,'(A)') repeat('-', 70)
    write(u,'(A,I0)')      'N modes      : ', self % Nmodes
    write(u,'(A,L1)')      'Auto N modes : ', self % auto_Nmodes
    write(u,'(A,F8.2,A)')  'Height H     : ', self % H, ' m'
    write(u,'(A,F8.2,A)')  'Upper HBC    : ', self % Upper_HBC, ' m'
    write(u,'(A,F8.2,A)')  'Lower HBC    : ', self % Lower_HBC, ' m'
    write(u,'(A,ES10.2,A)') 'Epsilon      : ', self % eps, ' m'
    write(u,'(A,I0)')      'Block size   : ', self % block_size
    write(u,'(A,L1)')      'In cluster   : ', self % in_cluster
    write(u,*)

    close(u)

    print*, '   ---> Text summary saved at ', dat_path

END SUBROUTINE write_txt_summary_AnalyticalNormalModes


! ===============================================================================
! HELPERS
! ===============================================================================
! -------------------------------------------------------------------------
! Subroutine: print_mode_summary
! Prints a formatted table of the retained normal modes: vertical
! wavenumber and wavelength, cutoff frequency, and horizontal wavenumber
! and wavelength evaluated at f_max, flagging each mode as propagating or
! evanescent. For more than 20 modes only the first 10 and last 5 rows
! are shown.
! -------------------------------------------------------------------------
SUBROUTINE print_mode_summary(self)
    CLASS(AnalyticalNormalModes_t), INTENT(IN) :: self

    ! Local variables
    REAL(WP)          :: fmin, fmax, k_medium, k_zm, lambda_zm, fc, kr_sq, krm_val
    INTEGER(I32)      :: i, m_idx, line_width
    CHARACTER(len=20) :: kzm_str, krm_str, lambda_str, status
    CHARACTER(len=16) :: re_str, im_str
    CHARACTER(len=*), PARAMETER :: title = "ANALYTICAL NORMAL MODES SETUP"
    CHARACTER(len=*), PARAMETER :: row_fmt = '(2X,A5,3(A,A18),A,A13,A,A20,A,A13,A,A)'
    INTEGER(I32), ALLOCATABLE :: modes_to_show(:)

    fmin       = minval(self % sources % Freqs)
    fmax       = maxval(self % sources % Freqs)
    k_medium   = 2.0_WP * PI * fmax / self % env % speed_sound_base
    line_width = 115

    print '(A)', ''
    print '(A)', repeat('=', line_width)
    print '(A)', repeat(' ', (line_width - len(title)) / 2) // title
    print '(A)', repeat('=', line_width)
    write(*, '(A,F10.2,A)')        '  Water Depth (H):      ', self % H, ' m'
    write(*, '(A,F10.2,A)')        '  Speed of Sound (c):   ', self % env % speed_sound_base, ' m/s'
    write(*, '(A,F10.2,A,F10.2,A)') '  Frequency Band:       ', fmin, ' Hz - ', fmax, ' Hz'
    write(*, '(A,I10)')            '  Retained Modes (m):   ', self % Nmodes
    write(*, '(A,F10.2,A)')        '  * Note: k_rm and Horiz. Wavelength are evaluated at f_max (', fmax, ' Hz)'
    print '(A)', repeat('-', line_width)
    write(*, '(2X,A5,A,A18,A,A12,A,A13,A,A20,A,A13,A,A)') &
        'Mode ', ' | ', 'k_zm [rad/m]      ', ' | ', 'Vert. WL [m]', ' | ', 'f_cutoff [Hz]', ' | ', &
        'k_rm [rad/m]        ', ' | ', 'Horiz. WL [m]', ' | ', 'Status'
    print '(A)', repeat('-', line_width)

    if (self % Nmodes <= 20) then
        modes_to_show = [(i, i = 1, self % Nmodes)]
    else
        modes_to_show = [(i, i = 1, 10), -1_I32, (i, i = self % Nmodes - 4, self % Nmodes)]
    end if

    do i = 1, size(modes_to_show)
        m_idx = modes_to_show(i)
        if (m_idx == -1) then
            write(*, '(2X,A5,A,A18,A,A12,A,A13,A,A20,A,A13,A,A)') &
                '...  ', ' | ', '...', ' | ', '...', ' | ', '...', ' | ', '...', ' | ', '...', ' | ', '...'
            cycle
        end if

        ! Vertical properties (frequency independent for p = 0 top / dp/dz = 0 bottom)
        k_zm      = (2.0_WP * REAL(m_idx, WP) - 1.0_WP) * PI / (2.0_WP * self % H)
        lambda_zm = 4.0_WP * self % H / (2.0_WP * REAL(m_idx, WP) - 1.0_WP)
        fc        = (2.0_WP * REAL(m_idx, WP) - 1.0_WP) * self % env % speed_sound_base / (4.0_WP * self % H)

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
            status = 'Propagating'
        else
            krm_val = sqrt(abs(kr_sq))
            write(im_str, '(F12.4)') krm_val
            krm_str = '0.0000+' // trim(adjustl(im_str)) // 'j'
            lambda_str = 'N/A (Decays)'
            status = 'Evanescent'
        end if

        write(*, '(2X,I5,A,A18,A,F12.2,A,F13.2,A,A20,A,A13,A,A)') &
            m_idx, ' | ', kzm_str, ' | ', lambda_zm, ' | ', fc, ' | ', krm_str, ' | ', lambda_str, ' | ', trim(status)
    end do

    print '(A)', repeat('=', line_width)
    print '(A)', ''

END SUBROUTINE print_mode_summary


! ===============================================================================
! METHOD (core physics)
! ===============================================================================
! -------------------------------------------------------------------------
! Subroutine: dipolar_pressure_NM
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
FUNCTION Psi_kzm(self, z, m, kzm) RESULT(psi)
    CLASS(ANalyticalNormalModes_t), INTENT(IN) :: self
    REAL(WP)    , INTENT(IN)           :: z(:)
    INTEGER(I32), INTENT(IN), OPTIONAL :: m
    REAL(WP)    , INTENT(IN), OPTIONAL :: kzm

    ! Local variables
    REAL(WP), ALLOCATABLE :: psi(:)

    ALLOCATE(psi(size(z)))

    if (PRESENT(m) .and. .not. PRESENT(kzm)) then
        psi = sqrt(2.0_WP * self % env % density_base / self % H) * &
              sin((2.0_WP * REAL(m, WP) - 1.0_WP) * PI * z / (2.0_WP * self % H))
    else if (PRESENT(kzm)) then
        psi = sqrt(2.0_WP * self % env % density_base / self % H) * sin(kzm * z)
    else
        error stop "Psi_kzm: Either m or kzm must be provided."
    end if
    
END FUNCTION Psi_kzm


FUNCTION dPsi_kzm_dz(self, z, m, kzm) RESULT(dpsi)
    CLASS(AnalyticalNormalModes_t), INTENT(IN) :: self
    REAL(WP)    , INTENT(IN)           :: z(:)
    INTEGER(I32), INTENT(IN), OPTIONAL :: m
    REAL(WP)    , INTENT(IN), OPTIONAL :: kzm

    ! Local variables
    REAL(WP), ALLOCATABLE :: dpsi(:)
    REAL(WP)              :: k_zm
    
    if (PRESENT(kzm)) then
        k_zm = kzm
    else if (PRESENT(m)) then
        k_zm = (2.0_WP * REAL(m, WP) - 1.0_WP) * PI / (2.0_WP * self % H)
    else
        error stop "dPsi_kzm_dz: Either m or kzm must be provided."
    end if

    allocate(dpsi(size(z)))
    dpsi = sqrt(2.0_WP * self % env % density_base / self % H) * k_zm * cos(k_zm * z)

END FUNCTION dPsi_kzm_dz


SUBROUTINE dipolar_pressure_NM(self, obs_pos, Freqs, nodes_pos, Force, p_out)
    CLASS(AnalyticalNormalModes_t), INTENT(IN) :: self
    REAL(WP), INTENT(IN)                        :: obs_pos(:,:)     ! (Nob, 3)  ; caller must pass a 2D slice (Nob,3)
    REAL(WP), INTENT(IN)                        :: Freqs(:)         ! (NFreqs)
    REAL(WP), INTENT(IN)                        :: nodes_pos(:,:)   ! (Nnodes, 3)
    COMPLEX(WP), INTENT(IN)                     :: Force(:,:,:)     ! (NFreqs, Nnodes, 3)
    COMPLEX(WP), ALLOCATABLE, INTENT(OUT)       :: p_out(:,:)       ! (NFreqs, Nob)

    ! Local variables
    INTEGER(I32) :: NFreqs, Nnodes, Nob, Nm
    INTEGER(I32) :: j, f, k, im
    REAL(WP)     :: epsv, rho, A_norm, kappa
    COMPLEX(WP)  :: k_rm, H0, H1, fx, fy, fz, coef_z, coef_xy, F_dot_rhat
    REAL(WP), ALLOCATABLE :: dx(:,:), dy(:,:), r(:,:), r_inv(:,:)
    REAL(WP), ALLOCATABLE :: k_wave(:), k_zm(:)
    REAL(WP), ALLOCATABLE :: psi_z(:,:), psi_zs(:,:), dpsi_zs(:,:)


    ! Dimensions
    Nfreqs = size(Freqs)
    Nnodes = size(nodes_pos, 1)
    Nob    = size(obs_pos, 1)
    Nm     = self % Nmodes

    epsv = self % eps
    rho  = self % env % density_base
    A_norm = sqrt(2.0_WP * rho / self % H)

    ALLOCATE(p_out(NFreqs, Nob))
    p_out = (0.0_WP, 0.0_WP)
    k_wave = 2.0_WP * PI * Freqs / self % env % speed_sound_base

    ! Horizontal source-observer geometry (frequency and node independent)
    ALLOCATE(dx(Nob, Nnodes), dy(Nob, Nnodes), r(Nob, Nnodes), r_inv(Nob, Nnodes))
    do j = 1, Nnodes
        !$OMP SIMD
        do k = 1, Nob
            dx(k,j) = obs_pos(k, 1) - nodes_pos(j, 1)
            dy(k,j) = obs_pos(k, 2) - nodes_pos(j, 2)
            r (k,j) = sqrt(dx(k,j)*dx(k,j) + dy(k,j)*dy(k,j))

            if (r(k, j) < epsv) r(k, j) = epsv
            r_inv(k, j) = 1.0_WP / r(k, j)
        end do
    end do

    ! Vertical modal shapes (Frequency independent)
    ALLOCATE(k_zm(Nm), psi_z(Nob, Nm), psi_zs(Nnodes, Nm), dpsi_zs(Nnodes, Nm))
    do im = 1, Nm
        k_zm(im)       = (2.0_WP * REAL(im, WP) - 1.0_WP) * PI / (2.0_WP * self % H)
        psi_z(:, im)   = self % Psi_kzm(z = obs_pos(:,3), kzm = k_zm(im))
        psi_zs(:, im)  = self % Psi_kzm(z = nodes_pos(:,3), kzm = k_zm(im))
        dpsi_zs(:, im) = self % dPsi_kzm_dz(z = nodes_pos(:,3), kzm = k_zm(im))
    end do

    ! Loop over frequencies
    do f = 1, NFreqs
        ! Normal modes superposition
        do im = 1, Nm

            ! Horizontal wavenumber: real (propagating) or pure imaginary (evanescent)
            if (k_wave(f) >= k_zm(im)) then
                k_rm = CMPLX(sqrt(k_wave(f)**2 - k_zm(im)**2), 0.0_WP, kind = WP)
            else
                kappa = sqrt(k_zm(im)**2 - k_wave(f)**2)
                k_rm = I1 * kappa
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

                    F_dot_rhat  = (fx * dx(k,j) + fy*dy(k,j)) * r_inv(k,j)
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
