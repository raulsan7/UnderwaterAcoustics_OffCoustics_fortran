! ===============================================================================
! MODULE: WindTurbine
! ===============================================================================
! Defines WindTurbine_t, the abstract parent type for all wind turbine models
! (Monopile, Floating, etc.), an extension of Emitter_t. It centralizes every
! shared parameter, the standard OpenFAST output/SubDyn reading pipeline,
! the drivetrain frequency filtering and the Emitter_t procedures
! (translation of the structural nodes into acoustic point-dipole sources,
! HDF5 and .dat blocks), so that each concrete turbine type only has to
! implement its own force/mass computation, drivetrain excitation spectrum
! and platform-specific saved parameters.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE WindTurbine

USE Kinds
USE Emitter, ONLY: Emitter_t
USE StringUtils, ONLY: to_upper
USE HDF5Utils, ONLY: save_variable
USE SignalProcessing, ONLY: HANN_ENBW
USE ConfigTypes, ONLY: Config_WindTurbine_t, Config_t
USE SourceTypes, ONLY: AcousticSource_t, SOURCE_DIPOLE

IMPLICIT NONE

PRIVATE
PUBLIC :: WindTurbine_t, save_turbine_params, write_turbine_txt


TYPE, ABSTRACT, EXTENDS(Emitter_t) :: WindTurbine_t
    ! -------------------------------------------------------------------------
    ! This acts as the parent class for all turbine types (Monopile, Floating, etc.).
    ! It encapsulates all shared parameters, configuration, and standard OpenFAST
    ! input reading to completely eliminate duplicate code and variables.
    ! -------------------------------------------------------------------------

    ! --- Input parameters --- !
    CHARACTER(len=200) :: rootname                      ! [-] Name without extensions of the OpenFAST output files
    CHARACTER(len=200) :: OP_output_dir                 ! [-] Directory with OpenFAST output files
    LOGICAL            :: Binary                        ! [-] Whether OpenFAST files are binary (.outb) or not (.out)
    REAL(WP)           :: WindSpeed                     ! [m/s] Wind speed in norm
    REAL(WP)           :: WindDir                       ! [deg] Wind direction, 0 deg points to +x axis (anticlockwise from +x)
    REAL(WP)           :: Depth                         ! [m] Water depth
    REAL(WP)           :: AxisPos(2)                    ! [m] Position of the turbine axis in xy plane
    REAL(WP)           :: BariPos(2)                    ! [m] Position of the turbine barycenter in xy plane
    INTEGER(I32)       :: Nmembers                      ! [-] Number of structural members in the OpenFAST model
    INTEGER(I32)       :: Nnodes                        ! [-] Number of structural nodes in the OpenFAST model

    ! --- Dynamic physical arrays --- !
    REAL(WP), ALLOCATABLE :: Time(:)                    ! [s] Simulation time vector
    REAL(WP), ALLOCATABLE :: acc(:,:,:,:)               ! [m/s^2] Acceleration array (nt, Nmembers, Nnodes, 3)
    REAL(WP), ALLOCATABLE :: xyz_all(:,:,:)             ! [m] Nodes coordinates (Nmembers, Nnodes, 3)
    REAL(WP), ALLOCATABLE :: xyz(:,:)                   ! [m] Nodes coordinates list (Nnodes_total, 3)
    REAL(WP), ALLOCATABLE :: Freqs(:)                   ! [Hz] Frequency array

    COMPLEX(WP), ALLOCATABLE :: Force(:,:,:)            ! [N] Force array (Nfreqs, Nnodes_total, 3)

    CHARACTER(len=:), ALLOCATABLE :: OP_output_path     ! [-] OpenFAST output file path
    CHARACTER(len=:), ALLOCATABLE :: OP_SubDyn_path     ! [-] OpenFAST SubDyn file path
    CHARACTER(len=:), ALLOCATABLE :: OP_subdir          ! [-] Subfolder of OP_output_dir of the turbine type
    CHARACTER(len=:), ALLOCATABLE :: path_rpm           ! [-] Path where variables depending on rpm are located
    CHARACTER(len=:), ALLOCATABLE :: case_type          ! [-] Tag to distinguish turbine

    ! --- Derived variables --- !
    LOGICAL  :: is_monopile                             ! [-] Whether turbine is monopile or floating
    LOGICAL  :: pure_tones = .false.                    ! [-] Whether the drivetrain peaks are pure tones, not bands
    REAL(WP) :: df                                      ! [Hz] FFT bin width 1/(nt*dt), defines the shared frequency grid of a wind farm

    CONTAINS

    ! --- Shared procedures --- !
    PROCEDURE :: init                                   ! Polymorphic constructor; overridden by concrete turbine types.
    PROCEDURE :: init_base => init                      ! Alias for calling the shared parent initializer from a child override.
    PROCEDURE :: read_input                             ! Reads OpenFAST outputs
    PROCEDURE :: translate_align_with_WindDir           ! Places nodes correctly
    PROCEDURE :: filter_frequencies                     ! Keeps the frequencies near the drivetrain excitation peaks

    ! --- Emitter_t procedures --- !
    PROCEDURE :: to_source   => turbine_to_source       ! Translates turbine to source abstraction
    PROCEDURE :: save_params => save_params_turbine     ! Generic and type-specific parameters -> HDF5
    PROCEDURE :: write_txt   => write_txt_turbine       ! Generic and type-specific blocks of the .dat
    PROCEDURE :: get_enbw    => get_enbw_turbine        ! Hann window of compute_rfft (HANN_ENBW)

    ! --- Deferred (abstract) procedures --- !
    PROCEDURE(compute_force_i)                , DEFERRED :: compute_force
    PROCEDURE(compute_mass_i)                 , DEFERRED :: compute_mass
    PROCEDURE(excitation_spectrum_i)          , DEFERRED :: excitation_spectrum
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

    SUBROUTINE excitation_spectrum_i(self, rpm, freqs_amp, keys)
        ! -------------------------------------------------------------------------
        ! Builds the drivetrain excitation spectrum of the turbine type (peak
        ! frequencies and amplitudes, with their labels) at the given rotor
        ! speed. filter_frequencies keeps only the frequencies near those peaks.
        ! -------------------------------------------------------------------------

        IMPORT :: WindTurbine_t, WP

        CLASS(WindTurbine_t), INTENT(IN)               :: self              ! [-] Turbine instance
        REAL(WP)            , INTENT(IN)               :: rpm               ! [rpm] Rotor speed
        REAL(WP)            , INTENT(OUT), ALLOCATABLE :: freqs_amp(:,:)    ! [-] Peaks (N, 2): freq [Hz], amplitude
        CHARACTER(len=20)   , INTENT(OUT), ALLOCATABLE :: keys(:)           ! [-] Peak labels (N)

    END SUBROUTINE excitation_spectrum_i

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
    ! the input namelist files. The default OpenFAST paths are built inside
    ! the subfolder of OP_output_dir of the turbine type (OP_subdir), which
    ! every concrete type sets before calling this routine.
    ! -------------------------------------------------------------------------

    CLASS(WindTurbine_t)      , INTENT(INOUT) :: self      ! [-] Turbine instance
    TYPE(Config_WindTurbine_t), INTENT(IN)    :: config    ! [-] Turbine configuration

    ! Local variables
    CHARACTER(len=5) :: ext
    LOGICAL          :: exists


    if (.not. ALLOCATED(self % OP_subdir)) error stop "WindTurbine % init: OP_subdir not set by the turbine type " &
                                                      // trim(config % case_type)

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
    if (to_upper(trim(config % OP_output_path)) == 'NONE' .or. &    ! "None" or "Default", in any case
        to_upper(trim(config % OP_output_path)) == 'DEFAULT') then
        self % OP_output_path = trim(self % OP_output_dir) // self % OP_subdir // trim(self % rootname) // trim(ext)
    else
        self % OP_output_path = trim(config % OP_output_path)
    end if

    if (to_upper(trim(config % OP_SubDyn_path)) == 'NONE' .or. &    ! "None" or "Default", in any case
        to_upper(trim(config % OP_SubDyn_path)) == 'DEFAULT') then
        self % OP_SubDyn_path = trim(self % OP_output_dir) // self % OP_subdir // trim(self % rootname) &
                                // ".SD.sum.yaml"
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

    USE OpenFASTReader, ONLY: get_SDsum_variables, read_input_SD

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
        if (status /= 0_I32) error stop "WindTurbine % read_input: failed to read acceleration channels from " &
                                        // trim(self % OP_output_path)
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
        write(*, '(A, F5.1, A, F6.1, A)') 'Wind Speed: ', self % WindSpeed, ' m/s   | Wind Direction: ', &
                                          self % WindDir, ' deg'
        write(*, '(A, F5.1, A, F8.1, A, F8.1, A)') 'Water Depth: ', self % Depth, &
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
        !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(x_old, y_old)
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
        !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(x_old, y_old)
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


FUNCTION filter_frequencies(self, freqs_over) RESULT(mask)
    ! -------------------------------------------------------------------------
    ! Builds a logical mask selecting only the frequencies close to the
    ! drivetrain excitation peaks of the turbine type (its excitation_spectrum)
    ! at the current wind speed's rpm, so the force spectrum can be restricted
    ! to the physically relevant frequencies. Shared by every turbine type.
    ! A peak is a band of 50 sine components whose frequencies keep the bins
    ! within one bin of them. When the type flags its peaks as pure tones
    ! (pure_tones) a peak is a single line at its frequency, which keeps the
    ! bins within 1.5 bins of it: the three nearest bins, the main lobe of the
    ! Hann window of compute_rfft for a tone on a bin, whether or not the
    ! tone falls exactly on a bin.
    ! -------------------------------------------------------------------------

    USE CSVReader, ONLY: read_curve
    USE DrivetrainExcitation, ONLY: generate_timeseries_banded_sines, filter_non_useful_freqs

    CLASS(WindTurbine_t), INTENT(INOUT)        :: self          ! [-] Turbine instance
    REAL(WP)            , INTENT(IN), OPTIONAL :: freqs_over    ! [Hz] Frequencies below this are never filtered
    LOGICAL             , ALLOCATABLE          :: mask(:)       ! [-] Frequency mask (Nfreqs), .true. to keep

    ! Local variables
    CHARACTER(len=20), ALLOCATABLE :: keys(:)
    REAL(WP)         , ALLOCATABLE :: freqs_amp(:,:)
    REAL(WP)         , ALLOCATABLE :: freqs_to_use(:)
    REAL(WP)         , ALLOCATABLE :: tmp_signal(:)
    REAL(WP)                       :: rpm, freqs_over_, width_bins
    INTEGER(I32)                   :: nfreq


    freqs_over_ = 10.0_WP ; if (present(freqs_over)) freqs_over_ = freqs_over

    ! Sine components per peak and half-width [bins] of the band kept around each component
    nfreq      = merge(1_I32, 50_I32, self % pure_tones)
    width_bins = merge(1.5_WP, 1.0_WP, self % pure_tones)

    ! Read rpm at current wind speed
    rpm = read_curve(self % path_rpm, self % WindSpeed)

    ! Generate drivetrain excitation spectrum
    CALL self % excitation_spectrum(rpm=rpm, freqs_amp=freqs_amp, keys=keys)

    ! Get allowed frequencies
    CALL generate_timeseries_banded_sines(peaks=freqs_amp, keys=keys, t=self % Time, &
                                          zeta=0.02_WP, nfreq=nfreq, used_freqs=.true., a=tmp_signal, &
                                          freqs_out=freqs_to_use)

    ! Build mask
    mask = filter_non_useful_freqs(self % Freqs, freqs_to_use, freqs_over_, width_bins=width_bins)

    DEALLOCATE(tmp_signal)

    DEALLOCATE(freqs_amp, keys, freqs_to_use)

END FUNCTION filter_frequencies


SUBROUTINE turbine_to_source(self, config, sources)
    ! -------------------------------------------------------------------------
    ! Translates the physical turbine (nodes + nodal forces) into the abstract
    ! AcousticSource_t representation consumed by the acoustic solvers, applying
    ! the radiation-impedance correction of the forces, the same for every
    ! solver: IM and ANM propagate the same point dipoles, so both need it to
    ! account for the finite size of the members. The turbine forces are
    ! consumed: freed once the corrected copy exists, since nothing needs them
    ! afterwards (no second copy of Force stays in memory).
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

    ! Every node radiates as a point dipole
    ALLOCATE(sources % order(size(sources % xyz, 1)), source=SOURCE_DIPOLE)

    ! Radiation-impedance correction of the forces; the turbine forces are consumed
    sources % Force = self % get_impedance_corrected_force(config % env % speed_sound_base)
    DEALLOCATE(self % Force)

    ! Reference geometry and default output labels read by the solver (e.g. plot_mn_ANM24_ws11.4.hdf5)
    write(ws_str, '(F6.1)') self % WindSpeed

    sources % geometry % center      = self % BariPos
    sources % geometry % origin      = self % AxisPos
    sources % geometry % direction   = self % WindDir
    sources % geometry % out_dir     = "./turbine_acoustic_data/"
    sources % geometry % name_tag    = "plot_" // merge('mn', 'fl', self % is_monopile)
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


    sfx = '' ; if (present(istr)) sfx = '_' // istr

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


SUBROUTINE save_params_turbine(self, save_path)
    ! -------------------------------------------------------------------------
    ! Emitter_t % save_params of a single turbine: saves its generic
    ! parameters (save_turbine_params) and then its type-specific ones
    ! (save_specific_params) into the HDF5 output file.
    ! -------------------------------------------------------------------------

    CLASS(WindTurbine_t), INTENT(IN) :: self         ! [-] Turbine instance
    CHARACTER(len=*)    , INTENT(IN) :: save_path    ! [-] Path to the HDF5 output file


    CALL save_turbine_params(self, save_path)
    CALL self % save_specific_params(save_path)

END SUBROUTINE save_params_turbine


SUBROUTINE write_txt_turbine(self, unit)
    ! -------------------------------------------------------------------------
    ! Emitter_t % write_txt of a single turbine: writes its generic and
    ! type-specific blocks of the .dat summary file (write_turbine_txt).
    ! -------------------------------------------------------------------------

    CLASS(WindTurbine_t), INTENT(IN) :: self    ! [-] Turbine instance
    INTEGER(I32)        , INTENT(IN) :: unit    ! [-] Unit of the open .dat file


    CALL write_turbine_txt(self, unit)

END SUBROUTINE write_txt_turbine


FUNCTION get_enbw_turbine(self) RESULT(enbw)
    ! -------------------------------------------------------------------------
    ! Emitter_t % get_enbw of a single turbine: its forces come from
    ! compute_rfft, so their window is the Hann window (HANN_ENBW).
    ! -------------------------------------------------------------------------

    CLASS(WindTurbine_t), INTENT(IN) :: self    ! [-] Turbine instance
    REAL(WP)                         :: enbw    ! [-] Equivalent noise bandwidth in bins


    enbw = HANN_ENBW

END FUNCTION get_enbw_turbine

END MODULE WindTurbine
