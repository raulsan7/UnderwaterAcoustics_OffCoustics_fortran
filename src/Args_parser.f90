MODULE Args_parser

    USE Kinds
    USE Types
    
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: read_input_files

CONTAINS

SUBROUTINE read_input_files(files_path, config)
    CHARACTER(len=*), INTENT(IN)  :: files_path
    TYPE(Config_t)  , INTENT(OUT) :: config

    ! Local variables
    CHARACTER(len=200) :: turbine_config_path
    CHARACTER(len=200) :: solver_config_path
    CHARACTER(len=200) :: environment_config_path
    CHARACTER(len=200) :: ImagesMethod_config_path

    turbine_config_path      = files_path // "Turbine.nml"
    solver_config_path       = files_path // "Solver.nml"
    environment_config_path  = files_path // "Environment.nml"
    ImagesMethod_config_path = files_path // "ImagesMethod.nml"

    CALL read_turbine_config     (turbine_config_path     , config)
    CALL read_solver_config      (solver_config_path      , config)
    CALL read_environment_config (environment_config_path , config)

    SELECT CASE(trim(config % solver % method))
    CASE('ImageMethod')
        ALLOCATE(Config_ImageMethod_t :: config % method)
        CALL read_ImageMethod_config(ImagesMethod_config_path, config)
    
    CASE('AnalyticalNormalModes')
        error stop "Analytical normal modes: Not implemented yet"

    CASE DEFAULT
        error stop "read_input_files(): Unknown method "
    END SELECT


END SUBROUTINE read_input_files


SUBROUTINE read_turbine_config(path, config)
    CHARACTER(len=*), INTENT(IN)  :: path
    TYPE(Config_t)  , INTENT(INOUT) :: config

    ! Local variables
    INTEGER(I32) :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    CHARACTER(len=200) :: case_type
    CHARACTER(len=200) :: rootname
    CHARACTER(len=200) :: OP_output_dir
    CHARACTER(len=200) :: OP_output_filename
    CHARACTER(len=200) :: OP_SubDyn_filename
    INTEGER(I32)       :: Nmembers
    INTEGER(I32)       :: Nnodes
    INTEGER(I32)       :: skip
    INTEGER(I32)       :: skipf
    INTEGER(I32)       :: block_size
    REAL(WP)           :: WindSpeed
    REAL(WP)           :: WindDir
    REAL(WP)           :: Depth
    REAL(WP)           :: AxisPos(2)
    REAL(WP)           :: BariPos(2)
    REAL(WP)           :: From
    REAL(WP)           :: Upto
    REAL(WP)           :: freqs_over
    LOGICAL            :: Binary
    LOGICAL            :: filter_freqs

    
    ! ==============================
    ! Configuration namelist
    ! ==============================
    namelist                &
    /turbine_config/        &
        case_type,          &
        rootname,           &
        OP_output_dir,      &
        OP_output_filename, &
        OP_SubDyn_filename, &
        Binary,             &
        WindSpeed,          &
        WindDir,            &
        Depth,              &
        Nmembers,           &
        Nnodes,             &
        AxisPos,            &
        BariPos,            &
        From,               &
        Upto,               &
        skip,               &
        skipf,              &
        filter_freqs,       &
        freqs_over

    ! ==============================
    ! Defaults
    ! ==============================
    case_type          = 'DTU10MWMonopile'
    rootname           = 'DTU_DeltaWind_mn_ws11.4'
    OP_output_dir      = '../OP_output/'
    OP_output_filename = 'None'
    OP_SubDyn_filename = 'None'
    Binary             = .true.
    WindSpeed          = 11.4_WP
    WindDir            = 0.0_WP
    Depth              = 30.0_WP
    Nmembers           = 8_I32
    Nnodes             = 5_I32
    AxisPos            = (/ 0.0_WP, 0.0_WP /)
    BariPos            = (/ 0.0_WP, 0.0_WP /)
    From               = 0.0_WP
    Upto               = 1.0_WP
    skip               = 1_I32
    skipf              = 1_I32
    filter_freqs       = .true.
    freqs_over         = 10.0_WP

    ! ==============================
    ! Open and read namelist file
    ! ==============================
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=turbine_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0) then
        write(*,*) "ERROR: Failed to read &turbine_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! ==============================
    ! Pack values into Config_t 
    ! ==============================
    config % turbine % case_type       = trim(case_type)
    config % turbine % rootname        = trim(rootname)
    config % turbine % OP_output_dir   = trim(OP_output_dir)
    config % turbine % OP_output_path  = trim(OP_output_filename)
    config % turbine % OP_SubDyn_path  = trim(OP_SubDyn_filename)
    config % turbine % Binary          = Binary
    config % turbine % WindSpeed       = WindSpeed
    config % turbine % WindDir         = WindDir
    config % turbine % Depth           = Depth
    config % turbine % Nmembers        = Nmembers
    config % turbine % Nnodes          = Nnodes
    config % turbine % AxisPos         = AxisPos
    config % turbine % BariPos         = BariPos
    config % turbine % From            = From
    config % turbine % Upto            = Upto
    config % turbine % skip            = skip
    config % turbine % skipf           = skipf
    config % turbine % filter_freqs    = filter_freqs
    config % turbine % freqs_over      = freqs_over

END SUBROUTINE read_turbine_config


SUBROUTINE read_solver_config(path, config)
    CHARACTER(len=*), INTENT(IN)  :: path
    TYPE(Config_t)  , INTENT(INOUT) :: config

    ! Local variables
    INTEGER(I32) :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    CHARACTER(len=200) :: method
    CHARACTER(len=200) :: save_path
    INTEGER(I32)       :: block_size
    LOGICAL            :: in_cluster
    LOGICAL            :: debug

    
    ! ==============================
    ! Configuration namelist
    ! ==============================
    namelist            &
    /solver_config/     &
        method,         &
        save_path,      &
        block_size,     &
        in_cluster,     &
        debug

    ! ==============================
    ! Defaults
    ! ==============================
    method     = "ImageMethod"
    save_path  = "None"
    block_size = 128_I32
    in_cluster = .true.
    debug      = .false.

    ! ==============================
    ! Open and read namelist file
    ! ==============================
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=solver_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0) then
        write(*,*) "ERROR: Failed to read &solver_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! ==============================
    ! Pack values into Config_t 
    ! ==============================
    config % solver % method     = method
    config % solver % save_path  = save_path
    config % solver % block_size = block_size
    config % solver % in_cluster = in_cluster
    config % solver % debug      = debug

END SUBROUTINE read_solver_config


SUBROUTINE read_environment_config(path, config)
    CHARACTER(len=*), INTENT(IN)  :: path
    TYPE(Config_t)  , INTENT(INOUT) :: config

    ! Local variables
    INTEGER(I32) :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    REAL(WP) :: depth_base
    REAL(WP) :: speed_sound_base
    REAL(WP) :: density_base
    REAL(WP) :: reference_pressure

    
    ! ==============================
    ! Configuration namelist
    ! ==============================
    namelist                    &
    /environment_config/        &
        depth_base,             &
        speed_sound_base,       &
        density_base,           &
        reference_pressure

    ! ==============================
    ! Defaults
    ! ==============================
    depth_base         = 30.0_WP
    speed_sound_base   = 1500.0_WP
    density_base       = 1025.0_WP
    reference_pressure = 1e-6_WP

    ! ==============================
    ! Open and read namelist file
    ! ==============================
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=environment_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0) then
        write(*,*) "ERROR: Failed to read &environment_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! ==============================
    ! Pack values into Config_t 
    ! ==============================
    config % env % depth_base         = depth_base
    config % env % speed_sound_base   = speed_sound_base
    config % env % density_base       = density_base
    config % env % reference_pressure = reference_pressure

END SUBROUTINE read_environment_config


SUBROUTINE read_ImageMethod_config(path, config)
    CHARACTER(len=*), INTENT(IN)    :: path
    TYPE(Config_t)  , INTENT(INOUT) :: config

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    INTEGER(I32) :: N_images
    INTEGER(I32) :: UpperBC
    INTEGER(I32) :: LowerBC
    REAL(WP)     :: UpperHBC
    REAL(WP)     :: LowerHBC
    REAL(WP)     :: attenuation_surface
    REAL(WP)     :: attenuation_bottom
    REAL(WP)     :: epsilon

    namelist                    &
    /ImagesMethod_config/       &
        N_images,               &
        UpperBC,                &
        LowerBC,                &
        UpperHBC,                &
        LowerHBC,                &
        attenuation_surface,    &
        attenuation_bottom,     &
        epsilon

    ! Defaults
    N_images            = 30_I32
    UpperBC             = -1_I32
    LowerBC             = 1_I32
    UpperHBC            = 0.0_WP
    LowerHBC            = -30.0_WP
    attenuation_surface = 1.0_WP
    attenuation_bottom  = 0.5_WP
    epsilon             = 1.0e-6_WP

    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=ImagesMethod_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0) then
        write(*,*) "ERROR: Failed to read &ImagesMethod_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    SELECT TYPE(m => config % method)
    TYPE IS(Config_ImageMethod_t)
        m % N_images            = N_images
        m % UpperBC             = UpperBC
        m % LowerBC             = LowerBC
        m % UpperHBC            = UpperHBC
        m % LowerHBC            = LowerHBC
        m % attenuation_surface = attenuation_surface
        m % attenuation_bottom  = attenuation_bottom
        m % epsilon             = epsilon
    END SELECT


END SUBROUTINE read_ImageMethod_config


END MODULE Args_parser