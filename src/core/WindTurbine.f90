MODULE WindTurbine

USE Kinds
USE Types

IMPLICIT NONE

PRIVATE
PUBLIC :: WindTurbine_t, drivetrain10MW_excitation_spectrum, init



! ===============================================================================
! ABSTRACT: WindTurbine_t
! ===============================================================================
! This acts as the parent class for all turbine types (Monopile, Floating, etc.).
! It encapsulates all shared parameters, configuration, and standard OpenFAST
! input reading to completely eliminate duplicate code and variables.
! ===============================================================================
TYPE, ABSTRACT :: WindTurbine_t

    ! --- Input parameters --- !
    CHARACTER(len=200) :: rootname                      ! [-] Name without extensions of the OpenFAST output files
    CHARACTER(len=200) :: OP_output_dir                 ! [-] Directory with OpenFast output files
    LOGICAL            :: Binary                        ! [-] Wheter OpenFast files are binary (.outb) or not (.out)
    REAL(WP)           :: WindSpeed                     ! [m/s] Wind Speed in norm
    REAL(WP)           :: WindDir                       ! [deg] Wind direction 0 deg points to +x axis (anticlockwise from +x)
    REAL(WP)           :: Depth                         ! [m] Water depth
    REAL(WP)           :: AxisPos(2)                    ! [m] Position of the turbine axis in xy plane
    REAL(WP)           :: BariPos(2)                    ! [m] Position of the turbine baricenter in xy plane
    INTEGER(I32)       :: Nmembers                      ! [-] Number of structural members in the OpenFAST model
    INTEGER(I32)       :: Nnodes                        ! [-] Number of structural nodes in the OpenFAST model

    ! --- Dynamic Phisical arrays --- !
    REAL(WP), ALLOCATABLE :: Time(:)                    ! [s] Simulation time vector
    REAL(WP), ALLOCATABLE :: acc(:,:,:,:)               ! [m/s^2] Acceleration array (nt, Nm, Nn)
    REAL(WP), ALLOCATABLE :: x_all(:,:,:)               ! [m] Nodes coordinates (Nmembers, Nnodes,3)
    REAL(WP), ALLOCATABLE :: x(:,:)                     ! [m] Nodes cordinates list (Nnodes_total, 3)
    REAL(WP), ALLOCATABLE :: Freqs(:)                   ! [Hz] Frequency array

    COMPLEX(WP), ALLOCATABLE :: F(:,:,:)                ! [N] Force array (NFreqs, Nnodes_total, 3)

    CHARACTER(len=:), ALLOCATABLE :: OP_output_path     ! [-] Openfast output file path
    CHARACTER(len=:), ALLOCATABLE :: OP_SubDyn_path     ! [-] Openfast SubDyn file path
    CHARACTER(len=:), ALLOCATABLE :: path_rpm           ! [-] Path where variables depending on rpm are located
    CHARACTER(len=:), ALLOCATABLE :: case_type          ! [-] Tag to distinguish turbine

    ! --- Derived variables --- !
    LOGICAL :: is_monopile                              ! [-] Wheter turbine is monopile or floating

    CONTAINS
    ! --- Shared procedures --- !
    PROCEDURE :: init                           ! Constructor
    PROCEDURE :: read_input                     ! Reads OpenFast outputs
    PROCEDURE :: translate_align_with_WindDir   ! Place nodes correctly
    PROCEDURE :: turbine_to_source


    ! --- Deferred procedures --- !
    PROCEDURE(compute_force_i)                , DEFERRED :: compute_force
    PROCEDURE(compute_mass_i)                 , DEFERRED :: compute_mass
    PROCEDURE(filter_frequencies_i)           , DEFERRED :: filter_frequencies
    PROCEDURE(save_specific_params_i)         , DEFERRED :: save_specific_params
    PROCEDURE(get_impedance_corrected_force_i), DEFERRED :: get_impedance_corrected_force

END TYPE WindTurbine_t

! =========================================================================
! Abstract Interfaces for Deferred Procedures
! =========================================================================
! Compute force
ABSTRACT INTERFACE
    SUBROUTINE compute_force_i(self, config, verbose)
        IMPORT :: WindTurbine_t, Config_WindTurbine_t
        CLASS(WindTurbine_t)      , INTENT(INOUT) :: self
        TYPE(Config_WindTurbine_t), INTENT(IN)    :: config
        LOGICAL                   , INTENT(IN)    :: verbose
    END SUBROUTINE compute_force_i
END INTERFACE

! Compute mass
ABSTRACT INTERFACE
    SUBROUTINE compute_mass_i(self, mass, added_mass)
        IMPORT :: WindTurbine_t, WP
        CLASS(WindTurbine_t) , INTENT(INOUT):: self
        REAL(WP), ALLOCATABLE, INTENT(OUT) :: mass(:)
        REAL(WP), ALLOCATABLE, INTENT(OUT) :: added_mass(:)
    END SUBROUTINE compute_mass_i
END INTERFACE

! Filter frequencies
ABSTRACT INTERFACE
    FUNCTION filter_frequencies_i(self, freqs_over) RESULT(mask)
        IMPORT :: WindTurbine_t, WP
        CLASS(WindTurbine_t), INTENT(INOUT) :: self
        REAL(WP), INTENT(IN), OPTIONAL      :: freqs_over
        LOGICAL, ALLOCATABLE                :: mask(:)
    END FUNCTION filter_frequencies_i
END INTERFACE

! Save specific turbine parameters
ABSTRACT INTERFACE

    SUBROUTINE save_specific_params_i(self, save_path, ichar)
        IMPORT:: WindTurbine_t
        CLASS(WindTurbine_t), INTENT(IN) :: self
        CHARACTER(len=*)    , INTENT(IN) :: save_path
        CHARACTER(len=*)    , INTENT(IN), OPTIONAL :: ichar
    END SUBROUTINE save_specific_params_i
END INTERFACE

! Impedance correction
ABSTRACT INTERFACE
    FUNCTION get_impedance_corrected_force_i(self, c) RESULT(corrected_force)
        IMPORT :: WindTurbine_t, WP
        CLASS(WindTurbine_t), INTENT(IN) :: self
        REAL(WP), INTENT(IN) :: c
        COMPLEX(WP), ALLOCATABLE :: corrected_force(:,:,:)
    END FUNCTION get_impedance_corrected_force_i
END INTERFACE


CONTAINS

! --------------------------------------------------------------------------------
! Subroutine: init
! Initializes the turbine properties using the configuration loaded from input.nml
! --------------------------------------------------------------------------------
SUBROUTINE init(self, config)

    CLASS(WindTurbine_t)      , INTENT(INOUT) :: self
    TYPE(Config_WindTurbine_t), INTENT(INOUT) :: config

    ! Local variables
    CHARACTER(len=5) :: ext
    LOGICAL          :: exists = .false.

    self % rootname      = config % rootname
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
    ext = ".out"; if (self % Binary) ext = ".outb"
    
    ! Set paths
    if (config % OP_output_path == 'None') then
        self % OP_output_path = trim(self % OP_output_dir) // trim(self % rootname) // trim(ext)
    else
        self % OP_output_path = trim(config % OP_output_path)
    endif
    
    if (config % OP_SubDyn_path == 'None') then
        self % OP_SubDyn_path = trim(self % OP_output_dir) // trim(self % rootname) // ".SD.sum.yaml"
    else
        self % OP_SubDyn_path = trim(config % OP_SubDyn_path)
    end if
    
    ! Check if paths exists
    inquire(file=trim(self % OP_output_path), exist=exists)
    if (.not. exists) error stop "WindTurbine % init: output file not found " // trim(self % OP_output_path)
    inquire(file=trim(self % OP_SubDyn_path), exist=exists)
    if (.not. exists) error stop "WindTurbine % init: SubDyn file not found " // trim(self % OP_SubDyn_path)


    write(*,*) ">> Base turbine initialized for rootname: ", trim(self % rootname)

END SUBROUTINE init


! -------------------------------------------------------------------------
! Subroutine: read_input_base
! Reads either binary (.outb) or text (.out) OpenFAST output files.
! Since this formatting is standard for OpenFAST, having this in the base
! class avoids any code duplication between different turbine types!
! -------------------------------------------------------------------------
SUBROUTINE read_input(self, config, verbose)
    USE IOUtils, ONLY: get_SDsum_variables, read_input_SD

    CLASS(WindTurbine_t)      , INTENT(INOUT) :: self
    TYPE(Config_WindTurbine_t), INTENT(IN)    :: config
    LOGICAL                   , INTENT(IN)    :: verbose

    ! Local variables
    REAL(WP)          :: dt
    INTEGER(I32)      :: status
    CHARACTER(len=10) :: unit

    if (verbose) then
        write(*,*) "---------------------------------------------------------"
        write(*,*) "LOADING OPENFAST DATA"
        write(*,*) "Filepath:  ", TRIM(self % OP_output_path)
        write(*,*) "---------------------------------------------------------"
    end if
    ! Load initial positions from SubDyn summary file
    if (self % Nmembers > 0 .and. self % Nnodes >0) then
        CALL get_SDsum_variables(SD_path=self%OP_SubDyn_path, Nmembers=self%Nmembers, &
                                 Nnodes=self%Nnodes, verbose=.false., Nodes=self%x_all)
        CALL read_input_SD(filename=self%OP_output_path, what='acceleration', skip=config%skip,Nmembers=self%Nmembers, &
                           Nnodes=self%Nnodes, From=config%From, Upto=config%Upto, verbose=.false., Time_out=self%time, &
                           Array_out=self%acc, unit_out=unit, status=status)
    else 
        error stop "WindTurbine % read_input: Nmembers and Nnodes must be positive integers to read SubDyn nodes"
    end if

    ! Compute baricenter if not provided (Most probable is that if BariPos = 0 then is not provided)
    if (all(self % BariPos == (/ 0.0_WP, 0.0_WP /))) then
        self % BariPos(1) = sum(self%x_all(:,:,1)) / real(self%Nmembers * self%Nnodes, WP)
        self % BariPos(2) = sum(self%x_all(:,:,2)) / real(self%Nmembers * self%Nnodes, WP)
    end if

    ! Align with wind direction and translate to axis position
    CALL self % translate_align_with_WindDir()


    if (verbose) then
        dt = self%Time(2) - self%Time(1)
        write(*,*) ''
        write(*, '(A)') repeat('=', 70)
        write(*, '(A)') 'TURBINE: ' // trim(self%rootname)
        write(*, '(A)') repeat('-', 68)
        write(*, '(A, A)')  'Type: ', trim(self%case_type)
        write(*, '(A, F5.1, A, F5.1, A)') 'Wind Speed: ', self%WindSpeed, ' m/s   | Wind Direction: ', self%WindDir, ' deg'
        write(*, '(A, F5.1, A, F4.1, A, F4.1, A)') 'Water Depth: ', self%Depth,&
         ' m   | Position: (', self%AxisPos(1), ', ', self%AxisPos(2), ') m'
        write(*, '(A, I6, A, F6.4, A)') 'Time Series: ', size(self%Time), ' samples @ ', dt, ' s/sample'
        write(*, '(A)') repeat('=', 70)
        write(*,*) ">> OpenFAST data loaded successfully!"
        write(*,*)
    end if

END SUBROUTINE read_input


SUBROUTINE translate_align_with_WindDir(self)
    CLASS(WindTurbine_t), INTENT(INOUT) :: self

    REAL(WP) :: yaw, c, s
    REAL(WP) :: axis_xy(2)
    INTEGER(I32) :: i, j, m, n
    INTEGER(I32) :: nt, nm, nn
    REAL(WP) :: x_old, y_old

    yaw = self % WindDir * PI / 180.0_WP
    c   = cos(yaw)
    s   = sin(yaw)

    if (abs(self % WindDir) < 1.0e-12_WP) return

    axis_xy = self % AxisPos

    if (ALLOCATED(self % acc)) then
       nt = size(self % acc, 1)
       nm = size(self % acc, 2)
       nn = size(self % acc, 3)
       !$omp parallel do collapse(3)
       do n = 1, nn
          do m = 1, nm
             do i = 1, nt
                x_old = self % acc(i,m,n,1)
                y_old = self % acc(i,m,n,2)
                self % acc(i,m,n,1) = c * x_old - s * y_old
                self % acc(i,m,n,2) = s * x_old + c * y_old
             end do
          end do
       end do
       !$omp end parallel do
    end if

    if (ALLOCATED(self % x_all)) then
       nm = size(self % x_all, 1)
       nn = size(self % x_all, 2)
       !$omp parallel do collapse(2)
       do i = 1, nm
          do j = 1, nn
             x_old = self % x_all(i,j,1)
             y_old = self % x_all(i,j,2)
             self % x_all(i,j,1) = c * x_old - s * y_old + axis_xy(1)
             self % x_all(i,j,2) = s * x_old + c * y_old + axis_xy(2)
          end do
       end do
       !$omp end parallel do
    end if
END SUBROUTINE


SUBROUTINE turbine_to_source(self, config, sources)
    CLASS(WindTurbine_t)  , INTENT(IN)  :: self
    TYPE(Config_t)        , INTENT(IN)  :: config
    TYPE(AcousticSource_t), INTENT(OUT) :: sources


    if (.not. ALLOCATED(self % F) .or. .not. ALLOCATED(self % x)) then
        error stop "turbine_to_source(): turbine forces not computed"
    end if

    ! Copy coordinates and active frequencies
    sources % xyz  = self % x
    sources % Freqs = self % Freqs

    ! Apply force correction
    if (config % solver % method == "ImageMethod") then
        sources % force = self % get_impedance_corrected_force(config % env % speed_sound_base)
    else
        sources % force = self % F
    end if

    write(*,*) ">> [AcousticSource] Physical turbine forces translated into mathematical point dipoles"

END SUBROUTINE turbine_to_source




! ---------- Drivetrain excitations ---------- !
SUBROUTINE drivetrain10MW_excitation_spectrum(rpm, damping, p_shaft, alpha_mesh, freqs_amp, keys)
    !> Constructs the physically-consistent excitation spectrum for the
    !> DTU 10 MW medium-speed drivetrain.
    REAL(WP), INTENT(IN), OPTIONAL :: rpm
    REAL(WP), INTENT(IN), OPTIONAL :: damping
    REAL(WP), INTENT(IN), OPTIONAL :: p_shaft
    REAL(WP), INTENT(IN), OPTIONAL :: alpha_mesh
    REAL(WP), ALLOCATABLE, INTENT(OUT) :: freqs_amp(:,:)        ! shape(N,2)
    CHARACTER(len=20), ALLOCATABLE, INTENT(OUT) :: keys(:)      ! shape(N)

    ! Local variables
    REAL(WP), PARAMETER :: min_rpm = 6.0_WP, max_rpm = 9.6_WP
    REAL(WP) :: factor, rpm_clipped, modal_sum

    ! Torsional eigenfrequencies (Hz)
    REAL(WP), PARAMETER :: torsional_modes(14) = [ &
        3.942_WP, 16.683_WP, 33.811_WP, 61.736_WP, 72.015_WP, &
        108.447_WP, 175.792_WP, 184.351_WP, 214.403_WP, &
        245.476_WP, 274.443_WP, 336.874_WP, 346.436_WP, 391.730_WP ]

    ! Shaft frequenciess at min and max rpm
    REAL(WP), PARAMETER :: shaft_min(13) = [ &
        0.1_WP, 0.2_WP, 0.3_WP, 0.6_WP, &
        0.322_WP, 0.644_WP, 0.966_WP, &
        1.383_WP, 2.766_WP, 4.149_WP, &
        5.017_WP, 10.034_WP, 15.051_WP ]
    REAL(WP), PARAMETER :: shaft_max(13) = [ &
        0.160_WP, 0.320_WP, 0.480_WP, 0.640_WP, &
        0.515_WP, 1.030_WP, 1.545_WP, &
        2.213_WP, 4.426_WP, 6.639_WP, &
        8.028_WP, 16.056_WP, 24.024_WP ]
    CHARACTER(len=6), PARAMETER :: shaft_keys(13) = [ &
        "lss_1p", "lss_2p", "lss_3p", "lss_6p", &
        "is1_1p","is1_2p","is1_3p", &
        "is2_1p","is2_2p","is2_3p", &
        "hss_1p", "hss_2p", "hss_3p" ]

    ! Gear mesh frequencies at min and max rpm
    REAL(WP), PARAMETER :: mesh_min(9) = [ &
        10.30_WP, 20.6_WP, 30.9_WP, &
        48.952_WP, 97.904_WP, 146.856_WP, &
        93.423_WP, 186.846_WP, 280.269_WP ]
    REAL(WP), PARAMETER :: mesh_max(9) = [ &
        16.48_WP, 32.96_WP, 49.44_WP, &
        78.3_WP, 156.6_WP, 234.9_WP, &
        149.492_WP, 298.984_WP, 448.476_WP ]
    CHARACTER(len=11), PARAMETER :: mesh_keys(9) = [ &
        "lss_mesh_1p","lss_mesh_2p","lss_mesh_3p", &
        "ims_mesh_1p","ims_mesh_2p","ims_mesh_3p", &
        "hss_mesh_1p","hss_mesh_2p","hss_mesh_3p" ]

    ! Base amplitudes (hardcoded)
    REAL(WP), PARAMETER :: A_base_shaft(13) = [ &
        0.6_WP, 0.3_WP, 1.0_WP, 0.2_WP, &
        0.4_WP, 0.2_WP, 0.6_WP, &
        0.3_WP, 0.15_WP, 0.4_WP, &
        0.2_WP, 0.1_WP, 0.3_WP ]
    REAL(WP), PARAMETER :: A_base_mesh(9) = [ &
        0.8_WP, 0.8_WP, 0.8_WP, &
        0.5_WP, 0.5_WP, 0.5_WP, &
        0.3_WP, 0.3_WP, 0.3_WP ]

    INTEGER(I32) :: n_total, i, m, pos_p, pos_underscore, n
    REAL(WP), ALLOCATABLE :: freqs(:), amps_source(:)
    REAL(WP) :: H, rms, fn, rpm_, damping_, alpha_mesh_, p_shaft_
    CHARACTER(len=20), ALLOCATABLE :: keys_work(:)
    CHARACTER(len=20) :: key_part

    
    ! Defaults
    rpm_        = 9.6_WP ; if (PRESENT(rpm))        rpm_        = rpm
    damping_    = 0.02_WP; if (PRESENT(damping))    damping_    = damping
    p_shaft_    = 1.0_WP ; if (PRESENT(p_shaft))    p_shaft_    = p_shaft
    alpha_mesh_ = 0.8_WP ; if (PRESENT(alpha_mesh)) alpha_mesh_ = alpha_mesh

    ! Clip rpm
    rpm_clipped = min(max(rpm_, min_rpm), max_rpm)
    factor = (rpm_clipped - min_rpm) / (max_rpm - min_rpm)

    ! Interpolate frequencies
    n_total = size(shaft_min) + size(mesh_min)
    ALLOCATE(freqs(n_total), amps_source(n_total), keys_work(n_total))

    ! Shaft frequencies and amplitudes
    do i = 1, size(shaft_min)
        freqs(i) = shaft_min(i) + factor * (shaft_max(i) - shaft_min(i))
        keys_work(i) = shaft_keys(i)

        pos_underscore = index(trim(shaft_keys(i)), "_", back=.true.)
        pos_p         = index(trim(shaft_keys(i)), "p", back=.true.)

        if (pos_underscore > 0 .and. pos_p > pos_underscore) then
            key_part = trim(shaft_keys(i)(pos_underscore+1:pos_p-1))
            read(key_part, *) n
        else
            n = 1
        end if

        amps_source(i) = A_base_shaft(i) * (1.0_WP / real(n, WP) ** p_shaft_)
    end do

    ! Mesh frequencies and amplitudes
    do i = 1, size(mesh_min)
        pos_underscore = index(trim(mesh_keys(i)), "_", back=.true.)
        pos_p         = index(trim(mesh_keys(i)), "p", back=.true.)

        if (pos_underscore > 0 .and. pos_p > pos_underscore) then
            key_part = trim(mesh_keys(i)(pos_underscore+1:pos_p-1))
            read(key_part, *) m
        else
            m = 1
        end if

        freqs(size(shaft_min) + i) = mesh_min(i) + factor * (mesh_max(i) - mesh_min(i))
        keys_work(size(shaft_min) + i) = mesh_keys(i)
        amps_source(size(shaft_min) + i) = A_base_mesh(i) * exp(-alpha_mesh_ * (real(m, WP) - 1.0_WP))
    end do

    ! Modal transfer function H(f)
    do i = 1, n_total
        modal_sum = 0.0_WP
        do m = 1, size(torsional_modes)
            fn = torsional_modes(m)
            modal_sum = modal_sum + 1.0_WP /((1.0_WP - (freqs(i)/fn)**2)**2 + (2.0_WP * damping_ * freqs(i)/fn)**2)
        end do
        H = sqrt(modal_sum)
        amps_source(i) = amps_source(i) * H
    end do

    ! Normalize to unit RMS: sqrt(sum(A^2/2)) = 1
    rms = sqrt(sum(amps_source**2)/2.0_WP)
    amps_source = amps_source / rms

    ! Build output array
    ALLOCATE(freqs_amp(n_total, 2))
    freqs_amp(:,1) = freqs
    freqs_amp(:,2) = amps_source
    ALLOCATE(keys(n_total))
    keys = keys_work

    DEALLOCATE(freqs, amps_source, keys_work)

END SUBROUTINE drivetrain10MW_excitation_spectrum

END MODULE WindTurbine