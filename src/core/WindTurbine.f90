! ===============================================================================
! MODULE: WindTurbine
! ===============================================================================
! Defines WindTurbine_t, the abstract parent type for all wind turbine models
! (Monopile, Floating, etc.). It centralizes every shared parameter, the
! standard OpenFAST output/SubDyn reading pipeline, and the translation of
! structural nodes into acoustic point-dipole sources, so that each concrete
! turbine type only has to implement its own force/mass computation and
! platform-specific saved parameters.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE WindTurbine

USE Kinds
USE Types
USE IOUtils

IMPLICIT NONE

PRIVATE
PUBLIC :: WindTurbine_t, drivetrain10MW_excitation_spectrum, init, save_turbine_params, &
          write_turbine_txt


TYPE, ABSTRACT :: WindTurbine_t
    ! -------------------------------------------------------------------------
    ! This acts as the parent class for all turbine types (Monopile, Floating, etc.).
    ! It encapsulates all shared parameters, configuration, and standard OpenFAST
    ! input reading to completely eliminate duplicate code and variables.
    ! -------------------------------------------------------------------------

    ! --- Input parameters --- !
    CHARACTER(len=200) :: rootname         ! [-] Name without extensions of the OpenFAST output files
    CHARACTER(len=200) :: OP_output_dir    ! [-] Directory with OpenFAST output files
    LOGICAL            :: Binary           ! [-] Whether OpenFAST files are binary (.outb) or not (.out)
    REAL(WP)           :: WindSpeed        ! [m/s] Wind speed in norm
    REAL(WP)           :: WindDir          ! [deg] Wind direction, 0 deg points to +x axis (anticlockwise from +x)
    REAL(WP)           :: Depth            ! [m] Water depth
    REAL(WP)           :: AxisPos(2)       ! [m] Position of the turbine axis in xy plane
    REAL(WP)           :: BariPos(2)       ! [m] Position of the turbine barycenter in xy plane
    INTEGER(I32)       :: Nmembers         ! [-] Number of structural members in the OpenFAST model
    INTEGER(I32)       :: Nnodes           ! [-] Number of structural nodes in the OpenFAST model

    ! --- Dynamic physical arrays --- !
    REAL(WP), ALLOCATABLE :: Time(:)           ! [s] Simulation time vector
    REAL(WP), ALLOCATABLE :: acc(:,:,:,:)      ! [m/s^2] Acceleration array (nt, Nmembers, Nnodes, 3)
    REAL(WP), ALLOCATABLE :: xyz_all(:,:,:)    ! [m] Nodes coordinates (Nmembers, Nnodes, 3)
    REAL(WP), ALLOCATABLE :: xyz(:,:)          ! [m] Nodes coordinates list (Nnodes_total, 3)
    REAL(WP), ALLOCATABLE :: Freqs(:)          ! [Hz] Frequency array

    COMPLEX(WP), ALLOCATABLE :: Force(:,:,:)    ! [N] Force array (Nfreqs, Nnodes_total, 3)

    CHARACTER(len=:), ALLOCATABLE :: OP_output_path    ! [-] OpenFAST output file path
    CHARACTER(len=:), ALLOCATABLE :: OP_SubDyn_path    ! [-] OpenFAST SubDyn file path
    CHARACTER(len=:), ALLOCATABLE :: path_rpm          ! [-] Path where variables depending on rpm are located
    CHARACTER(len=:), ALLOCATABLE :: case_type         ! [-] Tag to distinguish turbine

    ! --- Derived variables --- !
    LOGICAL  :: is_monopile    ! [-] Whether turbine is monopile or floating
    REAL(WP) :: df             ! [Hz] FFT bin width 1/(nt*dt), defines the shared frequency grid of a wind farm

    CONTAINS

    ! --- Shared procedures --- !
    PROCEDURE :: init                            ! Constructor
    PROCEDURE :: read_input                      ! Reads OpenFAST outputs
    PROCEDURE :: translate_align_with_WindDir    ! Places nodes correctly
    PROCEDURE :: turbine_to_source               ! Translates turbine to source abstraction

    ! --- Deferred (abstract) procedures --- !
    PROCEDURE(compute_force_i)                , DEFERRED :: compute_force
    PROCEDURE(compute_mass_i)                 , DEFERRED :: compute_mass
    PROCEDURE(filter_frequencies_i)           , DEFERRED :: filter_frequencies
    PROCEDURE(save_specific_params_i)         , DEFERRED :: save_specific_params
    PROCEDURE(write_specific_params_txt_i)    , DEFERRED :: write_specific_params_txt
    PROCEDURE(get_impedance_corrected_force_i), DEFERRED :: get_impedance_corrected_force

END TYPE WindTurbine_t


! ===============================================================================
! ABSTRACT INTERFACES FOR DEFERRED PROCEDURES
! ===============================================================================
ABSTRACT INTERFACE

    SUBROUTINE compute_force_i(self, config, verbose)
        ! -------------------------------------------------------------------------
        ! Computes the nodal structural forces (Nfreqs, Nnodes_wet, 3) of the
        ! wetted structure from the OpenFAST accelerations.
        ! -------------------------------------------------------------------------

        IMPORT :: WindTurbine_t, Config_WindTurbine_t

        CLASS(WindTurbine_t)      , INTENT(INOUT) :: self       ! [-] Turbine instance
        TYPE(Config_WindTurbine_t), INTENT(IN)    :: config     ! [-] Turbine configuration
        LOGICAL                   , INTENT(IN)    :: verbose    ! [-] Whether to print a summary

    END SUBROUTINE compute_force_i

    SUBROUTINE compute_mass_i(self, mass, added_mass)
        ! -------------------------------------------------------------------------
        ! Computes the per-node dry (structural) and added (hydrodynamic) mass.
        ! -------------------------------------------------------------------------

        IMPORT :: WindTurbine_t, WP

        CLASS(WindTurbine_t), INTENT(INOUT)            :: self             ! [-] Turbine instance
        REAL(WP)            , INTENT(OUT), ALLOCATABLE :: mass(:)          ! [kg] Dry mass per node
        REAL(WP)            , INTENT(OUT), ALLOCATABLE :: added_mass(:)    ! [kg] Added mass per node

    END SUBROUTINE compute_mass_i

    FUNCTION filter_frequencies_i(self, freqs_over) RESULT(mask)
        ! -------------------------------------------------------------------------
        ! Builds the mask of the frequencies to keep, based on the turbine's
        ! drivetrain excitation spectrum.
        ! -------------------------------------------------------------------------

        IMPORT :: WindTurbine_t, WP

        CLASS(WindTurbine_t), INTENT(INOUT)        :: self          ! [-] Turbine instance
        REAL(WP)            , INTENT(IN), OPTIONAL :: freqs_over    ! [Hz] Frequencies below this are never filtered
        LOGICAL             , ALLOCATABLE          :: mask(:)       ! [-] Frequency mask (Nfreqs), .true. to keep

    END FUNCTION filter_frequencies_i

    SUBROUTINE save_specific_params_i(self, save_path, istr)
        ! -------------------------------------------------------------------------
        ! Saves the turbine-specific parameters into the HDF5 output file.
        ! -------------------------------------------------------------------------

        IMPORT :: WindTurbine_t

        CLASS(WindTurbine_t), INTENT(IN)           :: self         ! [-] Turbine instance
        CHARACTER(len=*)    , INTENT(IN)           :: save_path    ! [-] Path to the HDF5 output file
        CHARACTER(len=*)    , INTENT(IN), OPTIONAL :: istr         ! [-] Suffix appended to the variable names

    END SUBROUTINE save_specific_params_i

    SUBROUTINE write_specific_params_txt_i(self, unit)
        ! -------------------------------------------------------------------------
        ! Writes the turbine-specific parameters as human-readable lines into
        ! the .dat summary file.
        ! -------------------------------------------------------------------------

        IMPORT :: WindTurbine_t, I32

        CLASS(WindTurbine_t), INTENT(IN) :: self    ! [-] Turbine instance
        INTEGER(I32)        , INTENT(IN) :: unit    ! [-] Unit of the open .dat file

    END SUBROUTINE write_specific_params_txt_i

    FUNCTION get_impedance_corrected_force_i(self, c) RESULT(corrected_force)
        ! -------------------------------------------------------------------------
        ! Returns the nodal forces corrected by the cylindrical-scattering
        ! (impedance) factor of the turbine's structural members.
        ! -------------------------------------------------------------------------

        IMPORT :: WindTurbine_t, WP

        CLASS(WindTurbine_t), INTENT(IN)  :: self                      ! [-] Turbine instance
        REAL(WP)            , INTENT(IN)  :: c                         ! [m/s] Speed of sound in the fluid
        COMPLEX(WP)         , ALLOCATABLE :: corrected_force(:,:,:)    ! [N] Corrected force (Nfreqs, Nnodes_wet, 3)

    END FUNCTION get_impedance_corrected_force_i

END INTERFACE


CONTAINS


SUBROUTINE init(self, config)
    ! -------------------------------------------------------------------------
    ! Initializes the turbine properties using the configuration loaded from
    ! the input namelist files.
    ! -------------------------------------------------------------------------

    CLASS(WindTurbine_t)      , INTENT(INOUT) :: self      ! [-] Turbine instance
    TYPE(Config_WindTurbine_t), INTENT(IN)    :: config    ! [-] Turbine configuration

    ! Local variables
    CHARACTER(len=5) :: ext
    LOGICAL          :: exists


    self % rootname      = config % rootname
    self % case_type     = config % case_type
    self % OP_output_dir = config % OP_output_dir
    self % Binary        = config % Binary
    self % WindSpeed     = config % WindSpeed
    self % WindDir       = config % WindDir
    self % Depth         = config % Depth
    self % Nmembers      = config % Nmembers
    self % Nnodes        = config % Nnodes
    self % AxisPos       = config % AxisPos
    self % BariPos       = config % BariPos

    ! Set extension
    ext = ".out" ; if (self % Binary) ext = ".outb"

    ! Set paths
    if (index(config % OP_output_path, 'one') > 0 .or. &    ! "None" or "Default"
        index(config % OP_output_path, 'efault') > 0) then
        self % OP_output_path = trim(self % OP_output_dir) // trim(self % rootname) // trim(ext)
    else
        self % OP_output_path = trim(config % OP_output_path)
    end if

    if (index(config % OP_SubDyn_path, 'one') > 0 .or. &    ! "None" or "Default"
        index(config % OP_SubDyn_path, 'efault') > 0) then
        self % OP_SubDyn_path = trim(self % OP_output_dir) // trim(self % rootname) // ".SD.sum.yaml"
    else
        self % OP_SubDyn_path = trim(config % OP_SubDyn_path)
    end if

    ! Check if paths exist
    inquire(file=trim(self % OP_output_path), exist=exists)
    if (.not. exists) error stop "WindTurbine % init: output file not found " // trim(self % OP_output_path)
    inquire(file=trim(self % OP_SubDyn_path), exist=exists)
    if (.not. exists) error stop "WindTurbine % init: SubDyn file not found " // trim(self % OP_SubDyn_path)

    write(*,*) ">> Base turbine initialized for rootname: ", trim(self % rootname)

END SUBROUTINE init


SUBROUTINE read_input(self, config, verbose)
    ! -------------------------------------------------------------------------
    ! Reads either binary (.outb) or text (.out) OpenFAST output files.
    ! Since this formatting is standard for OpenFAST, having this in the base
    ! class avoids any code duplication between different turbine types.
    ! -------------------------------------------------------------------------

    USE IOUtils, ONLY: get_SDsum_variables, read_input_SD

    CLASS(WindTurbine_t)      , INTENT(INOUT) :: self       ! [-] Turbine instance
    TYPE(Config_WindTurbine_t), INTENT(INOUT) :: config     ! [-] Turbine configuration
    LOGICAL                   , INTENT(IN)    :: verbose    ! [-] Whether to print a summary

    ! Local variables
    REAL(WP)          :: dt
    INTEGER(I32)      :: status
    CHARACTER(len=10) :: unit


    if (verbose) then
        write(*,*) "---------------------------------------------------------"
        write(*,*) "LOADING OPENFAST DATA"
        write(*,*) "Filepath:  ", trim(self % OP_output_path)
        write(*,*) "---------------------------------------------------------"
    end if

    ! Load initial positions from SubDyn summary file
    if (self % Nmembers > 0_I32 .and. self % Nnodes > 0_I32) then
        CALL get_SDsum_variables(SD_path=self % OP_SubDyn_path, Nmembers=self % Nmembers, &
                                 Nnodes=self % Nnodes, verbose=.false., Nodes=self % xyz_all)
        CALL read_input_SD(filename=self % OP_output_path, what='acceleration', skip=config % skip, &
                           Nmembers=self % Nmembers, Nnodes=self % Nnodes, From=config % From, &
                           Upto=config % Upto, verbose=.false., Time_out=self % Time, &
                           Array_out=self % acc, unit_out=unit, status=status)
    else
        error stop "WindTurbine % read_input: Nmembers and Nnodes must be positive integers to read SubDyn nodes"
    end if

    ! FFT bin width of the time series, defines the shared frequency grid of a wind farm
    self % df = 1.0_WP / (real(size(self % Time), WP) * (self % Time(2) - self % Time(1)))

    ! Align with wind direction and translate to axis position
    CALL self % translate_align_with_WindDir()

    ! Compute barycenter if not provided (BariPos = (0, 0) is taken as not provided)
    if (all(self % BariPos == [0.0_WP, 0.0_WP])) then
        self % BariPos(1) = sum(self % xyz_all(:, :, 1)) / real(self % Nmembers * self % Nnodes, WP)
        self % BariPos(2) = sum(self % xyz_all(:, :, 2)) / real(self % Nmembers * self % Nnodes, WP)
    end if

    if (verbose) then
        dt = self % Time(2) - self % Time(1)
        write(*,*) ''
        write(*, '(A)') repeat('=', 70)
        write(*, '(A)') 'TURBINE: ' // trim(self % rootname)
        write(*, '(A)') repeat('-', 68)
        write(*, '(A, A)') 'Type: ', trim(self % case_type)
        write(*, '(A, F5.1, A, F5.1, A)') 'Wind Speed: ', self % WindSpeed, ' m/s   | Wind Direction: ', &
                                          self % WindDir, ' deg'
        write(*, '(A, F5.1, A, F4.1, A, F4.1, A)') 'Water Depth: ', self % Depth, &
                                                   ' m   | Position: (', self % AxisPos(1), ', ', &
                                                   self % AxisPos(2), ') m'
        write(*, '(A, I6, A, F6.4, A)') 'Time Series: ', size(self % Time), ' samples @ ', dt, ' s/sample'
        write(*, '(A)') repeat('=', 70)
        write(*,*) ">> OpenFAST data loaded successfully!"
        write(*,*)
    end if

END SUBROUTINE read_input


SUBROUTINE translate_align_with_WindDir(self)
    ! -------------------------------------------------------------------------
    ! Rotates the acceleration and node-coordinate arrays by the wind direction
    ! yaw angle and translates them to the turbine axis position, so that all
    ! downstream computations work in the global reference frame. If WindDir
    ! is (numerically) zero, only the translation is applied.
    ! -------------------------------------------------------------------------

    CLASS(WindTurbine_t), INTENT(INOUT) :: self    ! [-] Turbine instance

    ! Local variables
    REAL(WP)     :: yaw, c, s
    REAL(WP)     :: axis_xy(2)
    INTEGER(I32) :: i, j, m, n
    INTEGER(I32) :: nt, nm, nn
    REAL(WP)     :: x_old, y_old


    yaw = self % WindDir * PI / 180.0_WP
    c   = cos(yaw)
    s   = sin(yaw)

    if (abs(self % WindDir) < 1.0e-12_WP) then
        if (any(self % AxisPos /= 0.0_WP)) then
            self % xyz_all(:, :, 1) = self % xyz_all(:, :, 1) + self % AxisPos(1)
            self % xyz_all(:, :, 2) = self % xyz_all(:, :, 2) + self % AxisPos(2)
        end if
        return
    end if

    axis_xy = self % AxisPos

    if (ALLOCATED(self % acc)) then
        nt = size(self % acc, 1)
        nm = size(self % acc, 2)
        nn = size(self % acc, 3)
        !$OMP PARALLEL DO COLLAPSE(3)
        do n = 1, nn
            do m = 1, nm
                do i = 1, nt
                    x_old = self % acc(i, m, n, 1)
                    y_old = self % acc(i, m, n, 2)
                    self % acc(i, m, n, 1) = c * x_old - s * y_old
                    self % acc(i, m, n, 2) = s * x_old + c * y_old
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end if

    if (ALLOCATED(self % xyz_all)) then
        nm = size(self % xyz_all, 1)
        nn = size(self % xyz_all, 2)
        !$OMP PARALLEL DO COLLAPSE(2)
        do i = 1, nm
            do j = 1, nn
                x_old = self % xyz_all(i, j, 1)
                y_old = self % xyz_all(i, j, 2)
                self % xyz_all(i, j, 1) = c * x_old - s * y_old + axis_xy(1)
                self % xyz_all(i, j, 2) = s * x_old + c * y_old + axis_xy(2)
            end do
        end do
        !$OMP END PARALLEL DO
    end if

END SUBROUTINE translate_align_with_WindDir


SUBROUTINE turbine_to_source(self, config, sources)
    ! -------------------------------------------------------------------------
    ! Translates the physical turbine (nodes + nodal forces) into the abstract
    ! AcousticSource_t representation consumed by the acoustic solvers, applying
    ! the impedance/scattering force correction when the selected solver method
    ! requires it (currently the Image Method).
    ! -------------------------------------------------------------------------

    CLASS(WindTurbine_t)  , INTENT(INOUT)            :: self       ! [-] Turbine instance
    TYPE(Config_t)        , INTENT(IN)               :: config     ! [-] Global configuration
    TYPE(AcousticSource_t), INTENT(OUT), ALLOCATABLE :: sources    ! [-] Acoustic point-dipole sources

    ! Local variables
    CHARACTER(len=16) :: ws_str


    ALLOCATE(AcousticSource_t :: sources)

    if (.not. ALLOCATED(self % Force) .or. .not. ALLOCATED(self % xyz)) then
        error stop "turbine_to_source(): turbine forces not computed"
    end if

    ! Copy coordinates and active frequencies
    sources % xyz   = self % xyz
    sources % Freqs = self % Freqs

    ! Apply force correction
    if (config % solver % method == "ImageMethod") then
        sources % Force = self % get_impedance_corrected_force(config % env % speed_sound_base)
    else
        CALL MOVE_ALLOC(self % Force, sources % Force)
    end if

    ! Reference geometry and default file name labels read by the solver
    write(ws_str, '(F6.1)') self % WindSpeed

    sources % geometry % center      = self % BariPos
    sources % geometry % origin      = self % AxisPos
    sources % geometry % direction   = self % WindDir
    sources % geometry % name_tag    = merge('mn', 'fl', self % is_monopile)
    sources % geometry % name_suffix = "_ws" // trim(adjustl(ws_str))

    write(*,*) ">> [AcousticSource] Physical turbine forces translated into mathematical point dipoles"

END SUBROUTINE turbine_to_source


SUBROUTINE save_turbine_params(turbine, save_path, istr)
    ! -------------------------------------------------------------------------
    ! Saves the generic turbine parameters (shared by every turbine type) into
    ! the solver's HDF5 output file. Platform-specific parameters are handled
    ! separately by each concrete type's save_specific_params.
    ! -------------------------------------------------------------------------

    CLASS(WindTurbine_t), INTENT(IN)           :: turbine       ! [-] Turbine instance
    CHARACTER(len=*)    , INTENT(IN)           :: save_path     ! [-] Path to the HDF5 output file
    CHARACTER(len=*)    , INTENT(IN), OPTIONAL :: istr          ! [-] Suffix appended to the variable names

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: sfx

    
    sfx = ''; if (PRESENT(istr)) sfx = '_' // istr

    CALL save_variable(save_path, "case_type" // sfx , turbine % case_type)
    CALL save_variable(save_path, "WindSpeed" // sfx , turbine % WindSpeed)
    CALL save_variable(save_path, "WindDir"   // sfx , turbine % WindDir)
    CALL save_variable(save_path, "Depth"     // sfx , turbine % Depth)
    CALL save_variable(save_path, "AxisPos"   // sfx , turbine % AxisPos)
    CALL save_variable(save_path, "BariPos"   // sfx , turbine % BariPos)

END SUBROUTINE save_turbine_params


SUBROUTINE write_turbine_txt(turbine, unit)
    ! -------------------------------------------------------------------------
    ! Writes the generic and the type-specific turbine blocks of the .dat
    ! summary file.
    ! -------------------------------------------------------------------------

    CLASS(WindTurbine_t), INTENT(IN) :: turbine    ! [-] Turbine instance
    INTEGER(I32)        , INTENT(IN) :: unit       ! [-] Unit of the open .dat file


    write(unit, '(A)') repeat('-', 70)
    write(unit, '(A)') 'TURBINE'
    write(unit, '(A)') repeat('-', 70)
    write(unit, '(A,A)')      'Case type   : ', trim(turbine % case_type)
    write(unit, '(A,L1)')     'Is monopile : ', turbine % is_monopile
    write(unit, '(A,F8.2,A)') 'Wind speed  : ', turbine % WindSpeed, ' m/s'
    write(unit, '(A,F8.2,A)') 'Wind dir    : ', turbine % WindDir, ' deg'
    write(unit, '(A,F8.2,A)') 'Water depth : ', turbine % Depth, ' m'
    write(unit, '(A,2F10.3)') 'Axis pos    : ', turbine % AxisPos
    write(unit, '(A,2F10.3)') 'Bari pos    : ', turbine % BariPos
    write(unit, '(A,I0)')     'N members   : ', turbine % Nmembers
    write(unit, '(A,I0)')     'N nodes     : ', turbine % Nnodes
    write(unit, *)

    write(unit, '(A)') repeat('-', 70)
    write(unit, '(A)') 'TURBINE - SPECIFIC PARAMETERS'
    write(unit, '(A)') repeat('-', 70)
    CALL turbine % write_specific_params_txt(unit)
    write(unit, *)

END SUBROUTINE write_turbine_txt


! ===============================================================================
! DRIVETRAIN EXCITATIONS
! ===============================================================================
SUBROUTINE drivetrain10MW_excitation_spectrum(rpm, damping, p_shaft, alpha_mesh, freqs_amp, keys)
    ! -------------------------------------------------------------------------
    ! Constructs the physically-consistent excitation spectrum (frequencies +
    ! amplitudes) for the DTU 10 MW medium-speed drivetrain, combining shaft
    ! and gear-mesh harmonics, RPM-dependent interpolation and a torsional
    ! modal transfer function, normalized to unit RMS.
    ! -------------------------------------------------------------------------

    REAL(WP)         , INTENT(IN), OPTIONAL     :: rpm               ! [rpm] Rotor speed
    REAL(WP)         , INTENT(IN), OPTIONAL     :: damping           ! [-] Modal damping ratio of the torsional modes
    REAL(WP)         , INTENT(IN), OPTIONAL     :: p_shaft           ! [-] Decay exponent of shaft harmonic amplitudes
    REAL(WP)         , INTENT(IN), OPTIONAL     :: alpha_mesh        ! [-] Decay rate of gear-mesh harmonic amplitudes
    REAL(WP)         , INTENT(OUT), ALLOCATABLE :: freqs_amp(:,:)    ! [-] Peaks (N, 2): frequency [Hz], amplitude
    CHARACTER(len=20), INTENT(OUT), ALLOCATABLE :: keys(:)           ! [-] Peak labels (N)

    ! Local variables
    REAL(WP), PARAMETER :: min_rpm = 6.0_WP, max_rpm = 9.6_WP
    REAL(WP)            :: factor, rpm_clipped, modal_sum

    ! Torsional eigenfrequencies [Hz]
    REAL(WP), PARAMETER :: torsional_modes(14) = [ &
        3.942_WP, 16.683_WP, 33.811_WP, 61.736_WP, 72.015_WP, &
        108.447_WP, 175.792_WP, 184.351_WP, 214.403_WP, &
        245.476_WP, 274.443_WP, 336.874_WP, 346.436_WP, 391.730_WP]

    ! Shaft frequencies at min and max rpm
    REAL(WP), PARAMETER :: shaft_min(13) = [ &
        0.1_WP, 0.2_WP, 0.3_WP, 0.6_WP, &
        0.322_WP, 0.644_WP, 0.966_WP, &
        1.383_WP, 2.766_WP, 4.149_WP, &
        5.017_WP, 10.034_WP, 15.051_WP]
    REAL(WP), PARAMETER :: shaft_max(13) = [ &
        0.160_WP, 0.320_WP, 0.480_WP, 0.640_WP, &
        0.515_WP, 1.030_WP, 1.545_WP, &
        2.213_WP, 4.426_WP, 6.639_WP, &
        8.028_WP, 16.056_WP, 24.024_WP]
    CHARACTER(len=6), PARAMETER :: shaft_keys(13) = [ &
        "lss_1p", "lss_2p", "lss_3p", "lss_6p", &
        "is1_1p", "is1_2p", "is1_3p", &
        "is2_1p", "is2_2p", "is2_3p", &
        "hss_1p", "hss_2p", "hss_3p"]

    ! Gear mesh frequencies at min and max rpm
    REAL(WP), PARAMETER :: mesh_min(9) = [ &
        10.30_WP, 20.6_WP, 30.9_WP, &
        48.952_WP, 97.904_WP, 146.856_WP, &
        93.423_WP, 186.846_WP, 280.269_WP]
    REAL(WP), PARAMETER :: mesh_max(9) = [ &
        16.48_WP, 32.96_WP, 49.44_WP, &
        78.3_WP, 156.6_WP, 234.9_WP, &
        149.492_WP, 298.984_WP, 448.476_WP]
    CHARACTER(len=11), PARAMETER :: mesh_keys(9) = [ &
        "lss_mesh_1p", "lss_mesh_2p", "lss_mesh_3p", &
        "ims_mesh_1p", "ims_mesh_2p", "ims_mesh_3p", &
        "hss_mesh_1p", "hss_mesh_2p", "hss_mesh_3p"]

    ! Base amplitudes (hardcoded)
    REAL(WP), PARAMETER :: A_base_shaft(13) = [ &
        0.6_WP, 0.3_WP, 1.0_WP, 0.2_WP, &
        0.4_WP, 0.2_WP, 0.6_WP, &
        0.3_WP, 0.15_WP, 0.4_WP, &
        0.2_WP, 0.1_WP, 0.3_WP]
    REAL(WP), PARAMETER :: A_base_mesh(9) = [ &
        0.8_WP, 0.8_WP, 0.8_WP, &
        0.5_WP, 0.5_WP, 0.5_WP, &
        0.3_WP, 0.3_WP, 0.3_WP]

    INTEGER(I32)                   :: n_total, i, m, pos_p, pos_underscore, n
    REAL(WP)         , ALLOCATABLE :: freqs(:), amps_source(:)
    REAL(WP)                       :: H, rms, fn, rpm_, damping_, alpha_mesh_, p_shaft_
    CHARACTER(len=20), ALLOCATABLE :: keys_work(:)
    CHARACTER(len=20)              :: key_part


    ! Defaults
    rpm_        = 9.6_WP  ; if (present(rpm))        rpm_        = rpm
    damping_    = 0.02_WP ; if (present(damping))    damping_    = damping
    p_shaft_    = 1.0_WP  ; if (present(p_shaft))    p_shaft_    = p_shaft
    alpha_mesh_ = 0.8_WP  ; if (present(alpha_mesh)) alpha_mesh_ = alpha_mesh

    ! Clip rpm
    rpm_clipped = min(max(rpm_, min_rpm), max_rpm)
    factor      = (rpm_clipped - min_rpm) / (max_rpm - min_rpm)

    ! Interpolate frequencies
    n_total = size(shaft_min) + size(mesh_min)
    ALLOCATE(freqs(n_total), amps_source(n_total), keys_work(n_total))

    ! Shaft frequencies and amplitudes
    do i = 1, size(shaft_min)
        freqs(i)     = shaft_min(i) + factor * (shaft_max(i) - shaft_min(i))
        keys_work(i) = shaft_keys(i)

        pos_underscore = index(trim(shaft_keys(i)), "_", back=.true.)
        pos_p          = index(trim(shaft_keys(i)), "p", back=.true.)

        if (pos_underscore > 0_I32 .and. pos_p > pos_underscore) then
            key_part = trim(shaft_keys(i)(pos_underscore+1:pos_p-1))
            read(key_part, *) n
        else
            n = 1_I32
        end if

        amps_source(i) = A_base_shaft(i) * (1.0_WP / real(n, WP)**p_shaft_)
    end do

    ! Mesh frequencies and amplitudes
    do i = 1, size(mesh_min)
        pos_underscore = index(trim(mesh_keys(i)), "_", back=.true.)
        pos_p          = index(trim(mesh_keys(i)), "p", back=.true.)

        if (pos_underscore > 0_I32 .and. pos_p > pos_underscore) then
            key_part = trim(mesh_keys(i)(pos_underscore+1:pos_p-1))
            read(key_part, *) m
        else
            m = 1_I32
        end if

        freqs(size(shaft_min)+i)       = mesh_min(i) + factor * (mesh_max(i) - mesh_min(i))
        keys_work(size(shaft_min)+i)   = mesh_keys(i)
        amps_source(size(shaft_min)+i) = A_base_mesh(i) * exp(-alpha_mesh_ * (real(m, WP) - 1.0_WP))
    end do

    ! Modal transfer function H(f)
    do i = 1, n_total
        modal_sum = 0.0_WP
        do m = 1, size(torsional_modes)
            fn        = torsional_modes(m)
            modal_sum = modal_sum + 1.0_WP / ((1.0_WP - (freqs(i) / fn)**2)**2 &
                                              + (2.0_WP * damping_ * freqs(i) / fn)**2)
        end do
        H              = sqrt(modal_sum)
        amps_source(i) = amps_source(i) * H
    end do

    ! Normalize to unit RMS: sqrt(sum(A^2/2)) = 1
    rms         = sqrt(sum(amps_source**2) / 2.0_WP)
    amps_source = amps_source / rms

    ! Build output array
    ALLOCATE(freqs_amp(n_total, 2))
    freqs_amp(:, 1) = freqs
    freqs_amp(:, 2) = amps_source
    ALLOCATE(keys(n_total))
    keys = keys_work

    DEALLOCATE(freqs, amps_source, keys_work)

END SUBROUTINE drivetrain10MW_excitation_spectrum

END MODULE WindTurbine
