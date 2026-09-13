MODULE Args_parser

    USE Kinds
    USE Types
    
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: read_input_file

CONTAINS

SUBROUTINE read_input_file(file_path, config)
    CHARACTER(len=*), INTENT(IN)  :: file_path
    TYPE(Config_t)  , INTENT(OUT) :: config

    ! Local variables
    INTEGER(I32) :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    LOGICAL            :: debug
    LOGICAL            :: in_cluster
    CHARACTER(len=200) :: save_path
    CHARACTER(len=50)  :: method

    REAL(WP)           :: depth_base
    REAL(WP)           :: speed_sound_base
    REAL(WP)           :: density_base
    REAL(WP)           :: reference_pressure
    REAL(WP)           :: attenuation_surface
    REAL(WP)           :: attenuation_bottom

    CHARACTER(len=200) :: rootname
    CHARACTER(len=200) :: OP_output_dir
    CHARACTER(len=200) :: OP_output_filename
    CHARACTER(len=200) :: OP_SubDyn_filename
    LOGICAL            :: Binary
    REAL(WP)           :: WindSpeed
    REAL(WP)           :: WindDir
    REAL(WP)           :: Depth
    INTEGER(I32)       :: Nmembers
    INTEGER(I32)       :: Nnodes
    INTEGER(I32)       :: skip
    INTEGER(I32)       :: block_size
    REAL(WP)           :: AxisPos(2)
    REAL(WP)           :: BariPos(2)
    REAL(WP)           :: From, Upto, freqs_over
    LOGICAL            :: filter_freqs

    ! Namelist definitions
    namelist /solver_config/ debug, in_cluster, save_path, method, block_size
    namelist /environment_config/ depth_base, speed_sound_base, density_base, &
                                    reference_pressure, attenuation_surface, &
                                    attenuation_bottom
    namelist /turbine_config/ rootname, OP_output_dir, OP_output_filename, OP_SubDyn_filename, &
                              Binary, WindSpeed, WindDir, Depth, Nmembers, Nnodes, AxisPos, BariPos, &
                              From, Upto, skip, filter_freqs

    ! ==============================
    ! Defaults
    ! ==============================
    ! Solver config
    debug      = .false.
    in_cluster = .false.
    save_path  = 'None'
    method     = 'ImageMethod'
    block_size = 128

    ! Environment config
    depth_base          = 30.0_WP
    speed_sound_base    = 1500.0_WP
    density_base        = 1025.0_WP
    reference_pressure  = 1E-06_WP
    attenuation_surface = 1.0_WP
    attenuation_bottom  = 0.5_WP

    ! Turbine config
    rootname           = 'DTU_DeltaWind_mn_ws11.4'
    OP_output_dir      = '../OP_output/'
    OP_output_filename = 'None'
    OP_SubDyn_filename = 'None'
    Binary             = .true.
    WindSpeed          = 11.4_WP
    WindDir            = 0.0_WP
    Depth              = 30.0_WP
    Nmembers           = 8
    Nnodes             = 5
    AxisPos            = (/ 0.0_WP, 0.0_WP /)
    BariPos            = (/ 0.0_WP, 0.0_WP /)
    From               = 0.0
    Upto               = 1.0
    skip               = 1
    filter_freqs       = .true.
    freqs_over         = 10.0_WP


    ! Open and read namelist file
    open(newunit=unit_num, file=file_path, status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
        write(*,*) "ERROR: Could not open input file: ", file_path
        stop
    end if
    
    REWIND(unit_num)
    read(unit_num, nml=solver_config, iostat=ierr)
    if (ierr /= 0) then
        write(*,*) "ERROR: Failed to read &solver_config"
        close(unit_num)
        stop
    end if

    REWIND(unit_num)
    read(unit_num, nml=environment_config, iostat=ierr)
    if (ierr /= 0) then
        write(*,*) "ERROR: Failed to read &environment_config"
        close(unit_num)
        stop
    end if

    REWIND(unit_num)
    read(unit_num, nml=turbine_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0) then
        write(*,*) "ERROR: Failed to read &turbine_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! Pack local variables back into Config_t structure
    config % solver % debug            = debug
    config % solver % in_cluster       = in_cluster
    config % solver % save_path        = save_path
    config % solver % method           = method
    config % solver %block_size        = blocK_size

    config % env % depth_base          = depth_base
    config % env % speed_sound_base    = speed_sound_base
    config % env % density_base        = density_base
    config % env % reference_pressure  = reference_pressure
    config % env % attenuation_surface = attenuation_surface
    config % env % attenuation_bottom  = attenuation_bottom

    config % turbine % rootname        = rootname
    config % turbine % OP_output_dir   = OP_output_dir
    config % turbine % OP_output_path  = OP_output_filename
    config % turbine % OP_SubDyn_path  = OP_SubDyn_filename
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
    config % turbine % filter_freqs    = filter_freqs
    config % turbine % freqs_over      = freqs_over

END SUBROUTINE read_input_file

END MODULE Args_parser