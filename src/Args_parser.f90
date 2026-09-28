! ===============================================================================
! MODULE: Args_parser
! ===============================================================================
! Reads the four .nml configuration files (Turbine.nml, Solver.nml,
! Environment.nml, ImagesMethod.nml) that describe a simulation case and
! packs their contents into the global Config_t structure defined in
! Types.f90. read_input_files is the single entry point: it locates the
! config files inside the given directory and dispatches to one reader
! subroutine per namelist, plus one extra reader chosen according to the
! selected acoustic solver method.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE Args_parser

USE Kinds
USE Types

IMPLICIT NONE

PRIVATE
PUBLIC :: read_input_files


CONTAINS


SUBROUTINE read_input_files(files_path, config)
    ! -------------------------------------------------------------------------
    ! Top-level entry point of the input-parsing stage. Builds the paths of
    ! the four .nml files from files_path, reads the turbine, solver and
    ! environment namelists, then dispatches to the method-specific config
    ! reader based on the selected acoustic solver (config % solver % method).
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)               :: files_path    ! [-] Directory of the .nml files (with trailing '/')
    TYPE(Config_t)  , INTENT(OUT), ALLOCATABLE :: config        ! [-] Global configuration filled in from the .nml files

    ! Local variables
    CHARACTER(len=200) :: turbine_config_path
    CHARACTER(len=200) :: farm_config_path
    CHARACTER(len=200) :: solver_config_path
    CHARACTER(len=200) :: environment_config_path
    CHARACTER(len=200) :: computation_config_path

    ! Method-specific configuration paths
    CHARACTER(len=200) :: ImagesMethod_config_path
    CHARACTER(len=200) :: AnalyticalNormalModes_config_path


    if (.not. ALLOCATED(config)) ALLOCATE(Config_t :: config)

    turbine_config_path     = files_path // "WindTurbine.nml"
    farm_config_path        = files_path // "WindFarm.nml"
    solver_config_path      = files_path // "Solver.nml"
    environment_config_path = files_path // "Environment.nml"

    ! Method-specific paths
    ImagesMethod_config_path          = files_path // "ImagesMethod.nml"
    AnalyticalNormalModes_config_path = files_path // "AnalyticalNormalModes.nml"

    ! Read principal files
    CALL read_solver_config(solver_config_path, config)
    CALL read_environment_config(environment_config_path, config)

    ! Emitter: a single turbine or a whole farm, each one with its own post-processing setup
    SELECT CASE (trim(config % solver % simulation_type))
    CASE ('turbine')
        CALL read_turbine_config(turbine_config_path, config)
        config % turbine % Depth = config % env % depth_base
        computation_config_path  = files_path // "ComputationSetupTurbine.nml"

    CASE ('farm')
        CALL read_farm_config(farm_config_path, config)
        computation_config_path = files_path // "ComputationSetupFarm.nml"

    CASE DEFAULT
        error stop "read_input_files(): Unknown simulation_type (use 'turbine' or 'farm')"
    END SELECT

    CALL read_computation_config(computation_config_path, config)

    ! Dispatch to the method-specific config reader based on which
    ! acoustic solver was selected in Solver.nml
    SELECT CASE (trim(config % solver % method))
    CASE ('ImageMethod')
        ALLOCATE(Config_ImageMethod_t :: config % method)
        CALL read_ImageMethod_config(ImagesMethod_config_path, config)

    CASE ('AnalyticalNormalModes')
        ALLOCATE(Config_AnalyticalNormalModes_t :: config % method)
        CALL read_AnalyticalNormalModes_config(AnalyticalNormalModes_config_path, config)

    CASE DEFAULT
        error stop "read_input_files(): Unknown method "
    END SELECT

END SUBROUTINE read_input_files


SUBROUTINE read_turbine_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "turbine_config" namelist from Turbine.nml (rootname, OpenFAST
    ! output paths, wind/depth conditions, node counts, time window, frequency
    ! filtering, etc.) and packs the values into config % turbine.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to Turbine.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
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
    REAL(WP)           :: WindSpeed
    REAL(WP)           :: WindDir
    REAL(WP)           :: AxisPos(2)
    REAL(WP)           :: BariPos(2)
    REAL(WP)           :: From
    REAL(WP)           :: Upto
    REAL(WP)           :: freqs_over
    LOGICAL            :: Binary
    LOGICAL            :: filter_freqs

    ! --- 1. Configuration namelist --- !
    NAMELIST /turbine_config/ &
        case_type,            &
        rootname,             &
        OP_output_dir,        &
        OP_output_filename,   &
        OP_SubDyn_filename,   &
        Binary,               &
        WindSpeed,            &
        WindDir,              &
        Nmembers,             &
        Nnodes,               &
        AxisPos,              &
        BariPos,              &
        From,                 &
        Upto,                 &
        skip,                 &
        skipf,                &
        filter_freqs,         &
        freqs_over


    ! --- 2. Defaults --- !
    case_type          = 'DTU10MWMonopile'
    rootname           = 'DTU_DeltaWind_mn_ws11.4'
    OP_output_dir      = '../OP_output/'
    OP_output_filename = 'None'
    OP_SubDyn_filename = 'None'
    Binary             = .true.
    WindSpeed          = 11.4_WP
    WindDir            = 0.0_WP
    Nmembers           = 8_I32
    Nnodes             = 5_I32
    AxisPos            = [0.0_WP, 0.0_WP]
    BariPos            = [0.0_WP, 0.0_WP]
    From               = 0.0_WP
    Upto               = 1.0_WP
    skip               = 1_I32
    skipf              = 1_I32
    filter_freqs       = .true.
    freqs_over         = 10.0_WP

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=turbine_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &turbine_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    config % turbine % case_type      = trim(case_type)
    config % turbine % rootname       = trim(rootname)
    config % turbine % OP_output_dir  = trim(OP_output_dir)
    config % turbine % OP_output_path = trim(OP_output_filename)
    config % turbine % OP_SubDyn_path = trim(OP_SubDyn_filename)
    config % turbine % Binary         = Binary
    config % turbine % WindSpeed      = WindSpeed
    config % turbine % WindDir        = WindDir
    config % turbine % Nmembers       = Nmembers
    config % turbine % Nnodes         = Nnodes
    config % turbine % AxisPos        = AxisPos
    config % turbine % BariPos        = BariPos
    config % turbine % From           = From
    config % turbine % Upto           = Upto
    config % turbine % skip           = skip
    config % turbine % skipf          = skipf
    config % turbine % filter_freqs   = filter_freqs
    config % turbine % freqs_over     = freqs_over

END SUBROUTINE read_turbine_config


SUBROUTINE read_farm_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads WindFarm.nml: the &farm_config block with the settings shared by
    ! every turbine (they define the common frequency grid) and one &turbine
    ! block per turbine, in order. Blocks are first counted and then read one
    ! by one, resetting the per-turbine defaults before each block so that no
    ! value leaks from the previous turbine. Every turbine takes the water
    ! depth from Environment.nml, which must therefore be read first.
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: iso_fortran_env, ONLY: IOSTAT_END

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to WindFarm.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr, n, i, j
    CHARACTER(len=200) :: iomsg

    ! Input variables: shared by the whole farm
    CHARACTER(len=100) :: farm_name
    REAL(WP)           :: From
    REAL(WP)           :: Upto
    INTEGER(I32)       :: skip
    INTEGER(I32)       :: skipf
    LOGICAL            :: filter_freqs
    REAL(WP)           :: freqs_over

    ! Input variables: one set per turbine
    CHARACTER(len=200) :: case_type
    CHARACTER(len=200) :: rootname
    CHARACTER(len=200) :: OP_output_dir
    CHARACTER(len=200) :: OP_output_filename
    CHARACTER(len=200) :: OP_SubDyn_filename
    LOGICAL            :: Binary
    REAL(WP)           :: WindSpeed
    REAL(WP)           :: WindDir
    REAL(WP)           :: AxisPos(2)
    REAL(WP)           :: BariPos(2)
    INTEGER(I32)       :: Nmembers
    INTEGER(I32)       :: Nnodes

    ! --- 1. Configuration namelists --- !
    NAMELIST /farm_config/ &
        farm_name,         &
        From,              &
        Upto,              &
        skip,              &
        skipf,             &
        filter_freqs,      &
        freqs_over

    NAMELIST /turbine/      &
        case_type,          &
        rootname,           &
        OP_output_dir,      &
        OP_output_filename, &
        OP_SubDyn_filename, &
        Binary,             &
        WindSpeed,          &
        WindDir,            &
        AxisPos,            &
        BariPos,            &
        Nmembers,           &
        Nnodes


    ! --- 2. Shared defaults and &farm_config block --- !
    farm_name    = 'None'
    From         = 0.0_WP
    Upto         = 1.0_WP
    skip         = 1_I32
    skipf        = 1_I32
    filter_freqs = .true.
    freqs_over   = 10.0_WP

    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=farm_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &farm_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    ! --- 3. Count the &turbine blocks --- !
    rewind(unit_num)
    n = 0_I32
    do
        read(unit_num, nml=turbine, iostat=ierr, iomsg=iomsg)
        if (ierr == IOSTAT_END) exit
        if (ierr /= 0_I32) then
            write(*, '(A,I0,A)') " ERROR: Failed to read &turbine block ", n + 1_I32, ": " // trim(iomsg)
            close(unit_num)
            stop
        end if
        n = n + 1_I32
    end do

    if (n == 0_I32) then
        write(*,*) "ERROR: No &turbine block found in ", path
        close(unit_num)
        stop
    end if

    config % farm % Nturbines = n
    config % farm % farm_name = trim(farm_name)
    ALLOCATE(config % farm % turbines(n))

    ! --- 4. Read every &turbine block, resetting the per-turbine defaults each time --- !
    rewind(unit_num)
    do i = 1, n
        case_type          = 'DTU10MWMonopile'
        rootname           = 'None'
        OP_output_dir      = '../OP_output/'
        OP_output_filename = 'None'
        OP_SubDyn_filename = 'None'
        Binary             = .true.
        WindSpeed          = 11.4_WP
        WindDir            = 0.0_WP
        AxisPos            = [0.0_WP, 0.0_WP]
        BariPos            = [0.0_WP, 0.0_WP]
        Nmembers           = 8_I32
        Nnodes             = 5_I32

        read(unit_num, nml=turbine)

        if (trim(rootname) == 'None') then
            write(*, '(A,I0)') " ERROR: rootname missing in &turbine block ", i
            close(unit_num)
            stop
        end if

        ! --- 5. Pack values into Config_t: per-turbine, shared and environment depth --- !
        associate(t => config % farm % turbines(i))
            t % case_type      = trim(case_type)
            t % rootname       = trim(rootname)
            t % OP_output_dir  = trim(OP_output_dir)
            t % OP_output_path = trim(OP_output_filename)
            t % OP_SubDyn_path = trim(OP_SubDyn_filename)
            t % Binary         = Binary
            t % WindSpeed      = WindSpeed
            t % WindDir        = WindDir
            t % AxisPos        = AxisPos
            t % BariPos        = BariPos
            t % Nmembers       = Nmembers
            t % Nnodes         = Nnodes
            t % From           = From
            t % Upto           = Upto
            t % skip           = skip
            t % skipf          = skipf
            t % filter_freqs   = filter_freqs
            t % freqs_over     = freqs_over
            t % Depth          = config % env % depth_base
        end associate
    end do

    close(unit_num)

    ! --- 6. Replicas of the same OpenFAST run radiate fully coherent sources --- !
    do i = 1, n
        do j = i + 1, n
            if (trim(config % farm % turbines(i) % rootname) == trim(config % farm % turbines(j) % rootname)) &
                write(*, '(A,I0,A,I0,A)') " WARNING read_farm_config: turbines ", i, " and ", j, &
                                          " share the same rootname (fully coherent sources)"
        end do
    end do

END SUBROUTINE read_farm_config


SUBROUTINE read_solver_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "solver_config" namelist from Solver.nml (which acoustic
    ! solver method to use, output save path, block size, cluster/debug
    ! flags) and packs the values into config % solver.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to Solver.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    CHARACTER(len=200) :: simulation_type
    CHARACTER(len=200) :: method
    CHARACTER(len=200) :: save_path
    INTEGER(I32)       :: block_size
    LOGICAL            :: in_cluster
    LOGICAL            :: debug

    ! --- 1. Configuration namelist --- !
    NAMELIST /solver_config/    &
        simulation_type,        &
        method,                 &
        save_path,              &
        block_size,             &
        in_cluster,             &
        debug


    ! --- 2. Defaults --- !
    simulation_type = "turbine"
    method          = "ImageMethod"
    save_path       = "None"
    block_size      = 128_I32
    in_cluster      = .true.
    debug           = .false.

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=solver_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &solver_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    config % solver % simulation_type = simulation_type
    config % solver % method          = method
    config % solver % save_path       = save_path
    config % solver % block_size      = block_size
    config % solver % in_cluster      = in_cluster
    config % solver % debug           = debug

END SUBROUTINE read_solver_config


SUBROUTINE read_environment_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "environment_config" namelist from Environment.nml (base
    ! water depth, speed of sound, density and reference pressure) and packs
    ! the values into config % env.
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
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=environment_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &environment_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    config % env % depth_base         = depth_base
    config % env % speed_sound_base   = speed_sound_base
    config % env % density_base       = density_base
    config % env % reference_pressure = reference_pressure

END SUBROUTINE read_environment_config


SUBROUTINE read_ImageMethod_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "ImagesMethod_config" namelist from ImagesMethod.nml (number
    ! of images, surface/bottom boundary condition types and heights,
    ! attenuation coefficients, epsilon) and packs the values into the
    ! Config_ImageMethod_t instance held by config % method. Only called when
    ! the Image Method solver was selected in Solver.nml.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ImagesMethod.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    INTEGER(I32) :: N_images
    INTEGER(I32) :: UpperBC
    INTEGER(I32) :: LowerBC
    REAL(WP)     :: UpperHBC
    REAL(WP)     :: attenuation_surface
    REAL(WP)     :: attenuation_bottom
    REAL(WP)     :: epsilon

    ! --- 1. Configuration namelist --- !
    NAMELIST /ImagesMethod_config/ &
        N_images,                  &
        UpperBC,                   &
        LowerBC,                   &
        UpperHBC,                  &
        attenuation_surface,       &
        attenuation_bottom,        &
        epsilon


    ! --- 2. Defaults --- !
    N_images            = 30_I32
    UpperBC             = -1_I32
    LowerBC             = 1_I32
    UpperHBC            = 0.0_WP
    attenuation_surface = 1.0_WP
    attenuation_bottom  = 0.5_WP
    epsilon             = 1.0e-6_WP

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=ImagesMethod_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &ImagesMethod_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    SELECT TYPE (m => config % method)
    TYPE IS (Config_ImageMethod_t)
        m % N_images            = N_images
        m % UpperBC             = UpperBC
        m % LowerBC             = LowerBC
        m % UpperHBC            = UpperHBC
        m % attenuation_surface = attenuation_surface
        m % attenuation_bottom  = attenuation_bottom
        m % epsilon             = epsilon
    END SELECT

END SUBROUTINE read_ImageMethod_config


SUBROUTINE read_AnalyticalNormalModes_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "AnalyticalNormalModes_config" namelist from
    ! AnalyticalNormalModes.nml (number of retained modes, upper/lower
    ! waveguide boundary heights, epsilon) and packs the values into the
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
    REAL(WP)     :: UpperHBC
    REAL(WP)     :: epsilon

    ! --- 1. Configuration namelist --- !
    NAMELIST /AnalyticalNormalModes_config/ &
        Nmodes,                             &
        UpperHBC,                           &
        epsilon


    ! --- 2. Defaults --- !
    Nmodes   = DEFAULT_I
    UpperHBC = 0.0_WP
    epsilon  = 1.0e-12_WP

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=AnalyticalNormalModes_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &AnalyticalNormalModes_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    SELECT TYPE (m => config % method)
    TYPE IS (Config_AnalyticalNormalModes_t)
        m % Nmodes   = Nmodes
        m % UpperHBC = UpperHBC
        m % epsilon  = epsilon
    END SELECT

END SUBROUTINE read_AnalyticalNormalModes_config


SUBROUTINE read_computation_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads every namelist group in ComputationSetupTurbine.nml or
    ! ComputationSetupFarm.nml (do_config plus one group per run_*
    ! post-processing step) and packs the values into config % run_inputs.
    ! Only the groups of the steps that will actually run are read.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml or ComputationSetupFarm.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration


    CALL read_do_config(path, config)
    if (config % run_inputs % do_run % do_spectrums)     CALL read_spectrums_config(path, config)
    if (config % run_inputs % do_run % do_polar)         CALL read_polar_config(path, config)
    if (config % run_inputs % do_run % do_cylinder)      CALL read_cylinder_config(path, config)
    if (config % run_inputs % do_run % do_decay)         CALL read_decay_config(path, config)
    if (config % run_inputs % do_run % do_line)          CALL read_line_config(path, config)
    if (config % run_inputs % do_run % do_sliceXY)       CALL read_sliceXY_config(path, config)
    if (config % run_inputs % do_run % do_sliceXZ)       CALL read_sliceXZ_config(path, config)
    if (config % run_inputs % do_run % do_sliceVertical) CALL read_sliceVertical_config(path, config)
    if (config % run_inputs % do_run % do_sphere)        CALL read_sphere_config(path, config)

END SUBROUTINE read_computation_config


SUBROUTINE read_do_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "do_config" namelist from ComputationSetupTurbine.nml or
    ! ComputationSetupFarm.nml and packs the post-processing execution
    ! switches into config % run_inputs. When do_all is enabled, every
    ! post-processing switch is forced to true. In farm mode the steps that
    ! make no sense for a wind farm (sphere, decay) are then switched off,
    ! so run() never needs to know what produced the sources.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml or ComputationSetupFarm.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration (solver settings already read)

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    LOGICAL :: do_all, do_spectrums, do_polar, do_cylinder, do_decay
    LOGICAL :: do_line, do_sliceXY, do_sliceXZ, do_sliceVertical, do_sphere

    ! --- 1. Configuration namelist --- !
    NAMELIST /do_config/  &
        do_all,           &
        do_spectrums,     &
        do_polar,         &
        do_cylinder,      &
        do_decay,         &
        do_line,          &
        do_sliceXY,       &
        do_sliceXZ,       &
        do_sliceVertical, &
        do_sphere


    ! --- 2. Defaults --- !
    do_all           = .true.
    do_spectrums     = .true.
    do_polar         = .true.
    do_cylinder      = .true.
    do_decay         = .true.
    do_line          = .true.
    do_sliceXY       = .true.
    do_sliceXZ       = .true.
    do_sliceVertical = .false.    ! New step: off unless requested, so older files keep their behaviour
    do_sphere        = .true.

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=do_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &do_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    config % run_inputs % do_run % do_all = do_all

    if (do_all) then
        config % run_inputs % do_run % do_spectrums     = .true.
        config % run_inputs % do_run % do_polar         = .true.
        config % run_inputs % do_run % do_cylinder      = .true.
        config % run_inputs % do_run % do_decay         = .true.
        config % run_inputs % do_run % do_line          = .true.
        config % run_inputs % do_run % do_sliceXY       = .true.
        config % run_inputs % do_run % do_sliceXZ       = .true.
        config % run_inputs % do_run % do_sliceVertical = .true.
        config % run_inputs % do_run % do_sphere        = .true.
    else
        config % run_inputs % do_run % do_spectrums     = do_spectrums
        config % run_inputs % do_run % do_polar         = do_polar
        config % run_inputs % do_run % do_cylinder      = do_cylinder
        config % run_inputs % do_run % do_decay         = do_decay
        config % run_inputs % do_run % do_line          = do_line
        config % run_inputs % do_run % do_sliceXY       = do_sliceXY
        config % run_inputs % do_run % do_sliceXZ       = do_sliceXZ
        config % run_inputs % do_run % do_sliceVertical = do_sliceVertical
        config % run_inputs % do_run % do_sphere        = do_sphere
    end if

    ! --- 5. Steps without meaning for a wind farm: no single turbine to enclose or to follow downstream --- !
    if (trim(config % solver % simulation_type) == 'farm') then
        if (.not. do_all .and. (do_sphere .or. do_decay)) &
            write(*,*) "WARNING read_do_config: sphere and decay are not available for a wind farm, skipped"
        config % run_inputs % do_run % do_sphere = .false.
        config % run_inputs % do_run % do_decay  = .false.
    end if

END SUBROUTINE read_do_config


SUBROUTINE read_spectrums_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "spectrums_config" namelist from ComputationSetupTurbine.nml
    ! (number of spectrum observer points and their coordinates) and packs
    ! the values into config % run_inputs.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    INTEGER(I32) :: n_observers
    REAL(WP)     :: observers_x(MAX_PARSE_OBS)
    REAL(WP)     :: observers_y(MAX_PARSE_OBS)
    REAL(WP)     :: observers_z(MAX_PARSE_OBS)

    ! --- 1. Configuration namelist --- !
    NAMELIST /spectrums_config/ &
        n_observers,            &
        observers_x,            &
        observers_y,            &
        observers_z


    ! --- 2. Defaults --- !
    n_observers = DEFAULT_I
    observers_x = DEFAULT_R
    observers_y = DEFAULT_R
    observers_z = DEFAULT_R

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=spectrums_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &spectrums_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    if (n_observers > MAX_PARSE_OBS) then
        write(*,*) "ERROR read_spectrums_config: n_observers exceeds the internal parsing cap of ", MAX_PARSE_OBS
        stop
    end if

    ! --- 4. Pack values into Config_t --- !
    config % run_inputs % params_spectrums % n_observers = n_observers

    if (n_observers /= DEFAULT_I .and. n_observers > 0_I32) then
        config % run_inputs % params_spectrums % observers_x = observers_x(1:n_observers)
        config % run_inputs % params_spectrums % observers_y = observers_y(1:n_observers)
        config % run_inputs % params_spectrums % observers_z = observers_z(1:n_observers)
    end if

END SUBROUTINE read_spectrums_config


SUBROUTINE read_polar_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "polar_config" namelist from ComputationSetupTurbine.nml
    ! (circle radius, depth, observer count and centre) and packs the values
    ! into config % run_inputs.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    REAL(WP)     :: r_polar, z_polar, center_polar(2)
    INTEGER(I32) :: n_polar

    ! --- 1. Configuration namelist --- !
    NAMELIST /polar_config/ &
        r_polar,            &
        z_polar,            &
        n_polar,            &
        center_polar


    ! --- 2. Defaults --- !
    r_polar      = DEFAULT_R
    z_polar      = DEFAULT_R
    n_polar      = DEFAULT_I
    center_polar = DEFAULT_R

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=polar_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &polar_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    config % run_inputs % params_polar % r_polar      = r_polar
    config % run_inputs % params_polar % z_polar      = z_polar
    config % run_inputs % params_polar % n_polar      = n_polar
    config % run_inputs % params_polar % center_polar = center_polar

END SUBROUTINE read_polar_config


SUBROUTINE read_cylinder_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "cylinder_config" namelist from ComputationSetupTurbine.nml
    ! (cylinder radius, azimuthal and vertical observer counts, and centre)
    ! and packs the values into config % run_inputs.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    REAL(WP)     :: r_cyl, center_cyl(2)
    INTEGER(I32) :: n_cyl, nz_cyl

    ! --- 1. Configuration namelist --- !
    NAMELIST /cylinder_config/ &
        r_cyl,                 &
        n_cyl,                 &
        nz_cyl,                &
        center_cyl


    ! --- 2. Defaults --- !
    r_cyl      = DEFAULT_R
    n_cyl      = DEFAULT_I
    nz_cyl     = DEFAULT_I
    center_cyl = DEFAULT_R

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=cylinder_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &cylinder_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    config % run_inputs % params_cylinder % r_cyl      = r_cyl
    config % run_inputs % params_cylinder % n_cyl      = n_cyl
    config % run_inputs % params_cylinder % nz_cyl     = nz_cyl
    config % run_inputs % params_cylinder % center_cyl = center_cyl

END SUBROUTINE read_cylinder_config


SUBROUTINE read_decay_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "decay_config" namelist from ComputationSetupTurbine.nml
    ! (start/end distances, observer count, depth and spacing mode) and packs
    ! the values into config % run_inputs.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    REAL(WP)     :: distance_decay(2), z_decay
    INTEGER(I32) :: n_decay
    LOGICAL      :: logspace_decay

    ! --- 1. Configuration namelist --- !
    NAMELIST /decay_config/ &
        distance_decay,     &
        n_decay,            &
        z_decay,            &
        logspace_decay


    ! --- 2. Defaults --- !
    distance_decay = DEFAULT_R
    n_decay        = DEFAULT_I
    z_decay        = DEFAULT_R
    logspace_decay = .true.

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=decay_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &decay_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t ---
    config % run_inputs % params_decay % distance_decay = distance_decay
    config % run_inputs % params_decay % n_decay        = n_decay
    config % run_inputs % params_decay % z_decay        = z_decay
    config % run_inputs % params_decay % logspace_decay = logspace_decay

END SUBROUTINE read_decay_config


SUBROUTINE read_line_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "line_config" namelist from ComputationSetupTurbine.nml
    ! (line endpoints, observer count and spacing mode) and packs the values
    ! into config % run_inputs.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    REAL(WP)     :: p1_line(3), p2_line(3)
    INTEGER(I32) :: n_line
    LOGICAL      :: logspace_line

    ! --- 1. Configuration namelist --- !
    NAMELIST /line_config/ &
        p1_line,           &
        p2_line,           &
        n_line,            &
        logspace_line


    ! --- 2. Defaults --- !
    p1_line       = DEFAULT_R
    p2_line       = DEFAULT_R
    n_line        = DEFAULT_I
    logspace_line = .false.

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=line_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &line_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    config % run_inputs % params_line % p1_line       = p1_line
    config % run_inputs % params_line % p2_line       = p2_line
    config % run_inputs % params_line % n_line        = n_line
    config % run_inputs % params_line % logspace_line = logspace_line

END SUBROUTINE read_line_config


SUBROUTINE read_sliceXY_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "sliceXY_config" namelist from ComputationSetupTurbine.nml
    ! (slice depth, grid dimensions, extents and centre) and packs the values
    ! into config % run_inputs.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    REAL(WP)     :: z_xy, xlim_xy(2), ylim_xy(2), center_xy(2)
    INTEGER(I32) :: nx_xy, ny_xy

    ! --- 1. Configuration namelist --- !
    NAMELIST /sliceXY_config/ &
        z_xy,                 &
        nx_xy,                &
        ny_xy,                &
        xlim_xy,              &
        ylim_xy,              &
        center_xy


    ! --- 2. Defaults --- !
    z_xy      = DEFAULT_R
    nx_xy     = DEFAULT_I
    ny_xy     = DEFAULT_I
    xlim_xy   = DEFAULT_R
    ylim_xy   = DEFAULT_R
    center_xy = DEFAULT_R

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=sliceXY_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &sliceXY_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    config % run_inputs % params_sliceXY % z_xy      = z_xy
    config % run_inputs % params_sliceXY % nx_xy     = nx_xy
    config % run_inputs % params_sliceXY % ny_xy     = ny_xy
    config % run_inputs % params_sliceXY % xlim_xy   = xlim_xy
    config % run_inputs % params_sliceXY % ylim_xy   = ylim_xy
    config % run_inputs % params_sliceXY % center_xy = center_xy

END SUBROUTINE read_sliceXY_config


SUBROUTINE read_sliceXZ_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "sliceXZ_config" namelist from ComputationSetupTurbine.nml
    ! (slice position, grid dimensions and x/z extents) and packs the values
    ! into config % run_inputs.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    REAL(WP)     :: y_xz, xlim_xz(2), zlim_xz(2)
    INTEGER(I32) :: nx_xz, nz_xz

    ! --- 1. Configuration namelist --- !
    NAMELIST /sliceXZ_config/ &
        y_xz,                 &
        nx_xz,                &
        nz_xz,                &
        xlim_xz,              &
        zlim_xz


    ! --- 2. Defaults --- !
    y_xz    = DEFAULT_R
    nx_xz   = DEFAULT_I
    nz_xz   = DEFAULT_I
    xlim_xz = DEFAULT_R
    zlim_xz = DEFAULT_R

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=sliceXZ_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &sliceXZ_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    config % run_inputs % params_sliceXZ % y_xz    = y_xz
    config % run_inputs % params_sliceXZ % nx_xz   = nx_xz
    config % run_inputs % params_sliceXZ % nz_xz   = nz_xz
    config % run_inputs % params_sliceXZ % xlim_xz = xlim_xz
    config % run_inputs % params_sliceXZ % zlim_xz = zlim_xz

END SUBROUTINE read_sliceXZ_config


SUBROUTINE read_sliceVertical_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "sliceVertical_config" namelist from ComputationSetupTurbine.nml
    ! or ComputationSetupFarm.nml (plane azimuth, centre, width, depth extent,
    ! grid dimensions and automatic-width margin) and packs the values into
    ! config % run_inputs. Since do_all also switches this newer step on, a
    ! file without the group is accepted and every value stays automatic.
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: iso_fortran_env, ONLY: IOSTAT_END

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml or ComputationSetupFarm.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    REAL(WP)     :: azimuth_sv, center_sv(2), width_sv, zlim_sv(2), margin_sv
    INTEGER(I32) :: nu_sv, nz_sv

    ! --- 1. Configuration namelist --- !
    NAMELIST /sliceVertical_config/ &
        azimuth_sv,                 &
        center_sv,                  &
        width_sv,                   &
        zlim_sv,                    &
        nu_sv,                      &
        nz_sv,                      &
        margin_sv


    ! --- 2. Defaults --- !
    azimuth_sv = 0.0_WP
    center_sv  = DEFAULT_R
    width_sv   = DEFAULT_R
    zlim_sv    = DEFAULT_R
    nu_sv      = DEFAULT_I
    nz_sv      = DEFAULT_I
    margin_sv  = 10.0_WP

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=sliceVertical_config, iostat=ierr, iomsg=iomsg)
    if (ierr == IOSTAT_END) then
        write(*,*) "WARNING read_sliceVertical_config: no &sliceVertical_config group, automatic values used"
    else if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &sliceVertical_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    config % run_inputs % params_sliceVertical % azimuth_sv = azimuth_sv
    config % run_inputs % params_sliceVertical % center_sv  = center_sv
    config % run_inputs % params_sliceVertical % width_sv   = width_sv
    config % run_inputs % params_sliceVertical % zlim_sv    = zlim_sv
    config % run_inputs % params_sliceVertical % nu_sv      = nu_sv
    config % run_inputs % params_sliceVertical % nz_sv      = nz_sv
    config % run_inputs % params_sliceVertical % margin_sv  = margin_sv
    config % run_inputs % params_sliceVertical % ulim_sv    = DEFAULT_R    ! Resolved by run_sliceVertical

END SUBROUTINE read_sliceVertical_config


SUBROUTINE read_sphere_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "sphere_config" namelist from ComputationSetupTurbine.nml
    ! (sphere radius, angular and vertical observer counts, and centre) and
    ! packs the values into config % run_inputs.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    REAL(WP)     :: r_sphere, center_sphere(3)
    INTEGER(I32) :: n_theta_sphere, nz_sphere

    ! --- 1. Configuration namelist --- !
    NAMELIST /sphere_config/ &
        r_sphere,            &
        n_theta_sphere,      &
        nz_sphere,           &
        center_sphere


    ! --- 2. Defaults --- !
    r_sphere       = DEFAULT_R
    n_theta_sphere = DEFAULT_I
    nz_sphere      = DEFAULT_I
    center_sphere  = DEFAULT_R

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        stop
    end if

    rewind(unit_num)
    read(unit_num, nml=sphere_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &sphere_config: ", trim(iomsg)
        close(unit_num)
        stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    config % run_inputs % params_sphere % r_sphere       = r_sphere
    config % run_inputs % params_sphere % n_theta_sphere = n_theta_sphere
    config % run_inputs % params_sphere % nz_sphere      = nz_sphere
    config % run_inputs % params_sphere % center_sphere  = center_sphere

END SUBROUTINE read_sphere_config

END MODULE Args_parser
