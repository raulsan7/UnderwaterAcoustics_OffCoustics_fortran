! ===============================================================================
! MODULE: InputSolver
! ===============================================================================
! Defines the readers of the solver namelists: Solver.nml
! (read_solver_config), Environment.nml (read_environment_config), and the
! method-specific AcousticMethods/ImageMethod.nml (read_ImageMethod_config)
! and AcousticMethods/AnalyticalNormalModes.nml
! (read_AnalyticalNormalModes_config).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE InputSolver

USE InputChecks, ONLY: check_fits
USE Kinds, ONLY: I32, WP, DEFAULT_I
USE ConfigTypes, ONLY: Config_ImageMethod_t, Config_AnalyticalNormalModes_t, Config_t

IMPLICIT NONE

PRIVATE
PUBLIC :: read_solver_config, read_environment_config, read_ImageMethod_config, read_AnalyticalNormalModes_config


CONTAINS


SUBROUTINE read_solver_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "solver_config" namelist from Solver.nml (simulation type,
    ! which acoustic solver method to use, output save path and its optional
    ! extra label, block size, cluster/debug flags, whether an existing
    ! output file is updated and the observer masking distance) and packs the
    ! values into config % solver.
    ! Strings are read into long buffers and rejected when they do not fit in
    ! Config_Solver_t, so a long save_path is never cut.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to Solver.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    CHARACTER(len=1024) :: simulation_type
    CHARACTER(len=1024) :: method
    CHARACTER(len=1024) :: save_path
    CHARACTER(len=1024) :: extra_label
    INTEGER(I32)        :: block_size
    LOGICAL             :: in_cluster
    LOGICAL             :: debug
    LOGICAL             :: update_existing
    REAL(WP)            :: mask_distance

    ! --- 1. Configuration namelist --- !
    NAMELIST /solver_config/    &
        simulation_type,        &
        method,                 &
        save_path,              &
        extra_label,            &
        block_size,             &
        in_cluster,             &
        debug,                  &
        update_existing,        &
        mask_distance


    ! --- 2. Defaults --- !
    simulation_type = "WindTurbine"
    method          = "ImageMethod"
    save_path       = "None"
    extra_label     = "None"
    block_size      = 128_I32
    in_cluster      = .true.
    debug           = .false.
    update_existing = .true.
    mask_distance   = 5.0_WP

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=solver_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &solver_config: ", trim(iomsg)
        close(unit_num)
        error stop
    end if

    close(unit_num)

    ! --- 4. Validate what Config_t cannot hold --- !
    CALL check_fits(simulation_type, config % solver % simulation_type, 'simulation_type')
    CALL check_fits(method, config % solver % method, 'method')
    CALL check_fits(save_path, config % solver % save_path, 'save_path')
    CALL check_fits(extra_label, config % solver % extra_label, 'extra_label')
    if (mask_distance < 0.0_WP) error stop "read_solver_config(): mask_distance must be >= 0 (0 = no masking)"

    ! --- 5. Pack values into Config_t --- !
    config % solver % simulation_type = trim(simulation_type)
    config % solver % method          = trim(method)
    config % solver % save_path       = trim(save_path)
    config % solver % extra_label     = trim(extra_label)
    config % solver % block_size      = block_size
    config % solver % in_cluster      = in_cluster
    config % solver % debug           = debug
    config % solver % update_existing = update_existing
    config % solver % mask_distance   = mask_distance

END SUBROUTINE read_solver_config


SUBROUTINE read_environment_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "environment_config" namelist from Environment.nml (base
    ! water depth, speed of sound, density and reference pressure) and packs
    ! the values into config % env. Every value must be strictly positive:
    ! a zero or negative depth, sound speed or density gives NaN pressures
    ! or a late stop in the solvers, and p_ref <= 0 breaks every dB level.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to Environment.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    REAL(WP) :: depth_base
    REAL(WP) :: speed_sound_base
    REAL(WP) :: density_base
    REAL(WP) :: reference_pressure

    ! --- 1. Configuration namelist --- !
    NAMELIST /environment_config/ &
        depth_base,               &
        speed_sound_base,         &
        density_base,             &
        reference_pressure


    ! --- 2. Defaults --- !
    depth_base         = 30.0_WP
    speed_sound_base   = 1500.0_WP
    density_base       = 1025.0_WP
    reference_pressure = 1.0e-6_WP

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=environment_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &environment_config: ", trim(iomsg)
        close(unit_num)
        error stop
    end if

    close(unit_num)

    ! --- 4. Validate: every property must be strictly positive (written with .not. to reject NaN too) --- !
    if (.not. (depth_base > 0.0_WP .and. speed_sound_base > 0.0_WP .and. density_base > 0.0_WP .and. &
               reference_pressure > 0.0_WP)) then
        write(*, '(A)') " ERROR: &environment_config: depth_base, speed_sound_base, density_base and " // &
                        "reference_pressure must be > 0"
        write(*, '(A,4ES14.6)') "        values read:", depth_base, speed_sound_base, density_base, &
                                reference_pressure
        error stop
    end if

    ! --- 5. Pack values into Config_t --- !
    config % env % depth_base         = depth_base
    config % env % speed_sound_base   = speed_sound_base
    config % env % density_base       = density_base
    config % env % reference_pressure = reference_pressure

END SUBROUTINE read_environment_config


SUBROUTINE read_ImageMethod_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "ImagesMethod_config" namelist from ImageMethod.nml (number
    ! of images, surface/bottom boundary condition types,
    ! attenuation coefficients, epsilon) and packs the values into the
    ! Config_ImageMethod_t instance held by config % method. Only called when
    ! the Image Method solver was selected in Solver.nml.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ImageMethod.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    INTEGER(I32) :: N_images
    INTEGER(I32) :: UpperBC
    INTEGER(I32) :: LowerBC
    REAL(WP)     :: attenuation_surface
    REAL(WP)     :: attenuation_bottom
    REAL(WP)     :: epsilon

    ! --- 1. Configuration namelist --- !
    NAMELIST /ImagesMethod_config/ &
        N_images,                  &
        UpperBC,                   &
        LowerBC,                   &
        attenuation_surface,       &
        attenuation_bottom,        &
        epsilon


    ! --- 2. Defaults --- !
    N_images            = 30_I32
    UpperBC             = -1_I32
    LowerBC             = 1_I32
    attenuation_surface = 1.0_WP
    attenuation_bottom  = 0.5_WP
    epsilon             = 1.0e-6_WP

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=ImagesMethod_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &ImagesMethod_config: ", trim(iomsg)
        close(unit_num)
        error stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    SELECT TYPE (m => config % method)
    TYPE IS (Config_ImageMethod_t)
        m % N_images            = N_images
        m % UpperBC             = UpperBC
        m % LowerBC             = LowerBC
        m % attenuation_surface = attenuation_surface
        m % attenuation_bottom  = attenuation_bottom
        m % epsilon             = epsilon
    END SELECT

END SUBROUTINE read_ImageMethod_config


SUBROUTINE read_AnalyticalNormalModes_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "AnalyticalNormalModes_config" namelist from
    ! AnalyticalNormalModes.nml (number of retained modes, upper/lower
    ! epsilon) and packs the values into the
    ! Config_AnalyticalNormalModes_t instance held by config % method. Only
    ! called when the Analytical Normal Modes solver was selected in Solver.nml.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to AnalyticalNormalModes.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    INTEGER(I32) :: Nmodes
    REAL(WP)     :: epsilon

    ! --- 1. Configuration namelist --- !
    NAMELIST /AnalyticalNormalModes_config/ &
        Nmodes,                             &
        epsilon


    ! --- 2. Defaults --- !
    Nmodes   = DEFAULT_I
    epsilon  = 1.0e-12_WP

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=AnalyticalNormalModes_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &AnalyticalNormalModes_config: ", trim(iomsg)
        close(unit_num)
        error stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    SELECT TYPE (m => config % method)
    TYPE IS (Config_AnalyticalNormalModes_t)
        m % Nmodes   = Nmodes
        m % epsilon  = epsilon
    END SELECT

END SUBROUTINE read_AnalyticalNormalModes_config

END MODULE InputSolver
