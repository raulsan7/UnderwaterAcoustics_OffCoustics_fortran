! ===============================================================================
! MODULE: Args_parser
! ===============================================================================
! Reads the .nml configuration files that describe a simulation case
! (Solver.nml, Environment.nml, WindTurbine.nml or WindFarm.nml,
! ComputationSetupTurbine.nml or ComputationSetupFarm.nml, and
! ImagesMethod.nml or AnalyticalNormalModes.nml) and packs their
! contents into the global Config_t structure defined in
! ConfigTypes.f90. get_input_path tells which directory holds them (the
! command-line argument, ./input_files/ by default) and read_input_files is the
! entry point of the parsing: it locates the config files inside that
! directory and dispatches to one reader subroutine per namelist (InputSolver,
! InputEmitters and InputPostProcessing), plus one extra reader chosen
! according to the selected acoustic solver method.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE Args_parser

USE Kinds
USE InputPostProcessing, ONLY: read_computation_config
USE InputEmitters, ONLY: read_turbine_config, read_farm_config, read_frequency_source_config
USE ConfigTypes, ONLY: Config_ImageMethod_t, Config_AnalyticalNormalModes_t, Config_t
USE InputSolver, ONLY: read_solver_config, read_environment_config, read_ImageMethod_config, &
                       read_AnalyticalNormalModes_config

IMPLICIT NONE

PRIVATE
PUBLIC :: get_input_path, read_input_files


CONTAINS


SUBROUTINE get_input_path(files_path)
    ! -------------------------------------------------------------------------
    ! Directory of the .nml files: the command-line argument (./off_coustics
    ! ./path_to_input_files) or, without arguments, ./input_files/. The
    ! trailing '/' that read_input_files needs is added when it is missing.
    ! More than one argument, or an empty one, stops the program.
    ! -------------------------------------------------------------------------

    CHARACTER(len=:), INTENT(OUT), ALLOCATABLE :: files_path    ! [-] Directory of the .nml files (with trailing '/')

    ! Local variables
    INTEGER(I32) :: nargs, arg_length, ierr


    nargs = command_argument_count()

    if (nargs == 0_I32) then
        files_path = "./input_files/"
        return
    end if

    if (nargs > 1_I32) then
        error stop "get_input_path(): Too many arguments. Usage: ./off_coustics [directory_of_the_input_files]"
    end if

    CALL get_command_argument(1, length=arg_length)
    if (arg_length == 0_I32) error stop "get_input_path(): The directory of the input files is an empty argument"

    ALLOCATE(CHARACTER(len=arg_length) :: files_path)
    CALL get_command_argument(1, files_path, status=ierr)
    if (ierr /= 0_I32) error stop "get_input_path(): Could not read the directory of the input files"

    if (files_path(arg_length:arg_length) /= "/") files_path = files_path // "/"

    write(*,*) ">> Reading the input files from: ", files_path

END SUBROUTINE get_input_path


SUBROUTINE read_input_files(files_path, config)
    ! -------------------------------------------------------------------------
    ! Top-level entry point of the input-parsing stage. Builds the paths of
    ! the .nml files from files_path, reads the solver and environment
    ! namelists, the turbine or wind farm namelists and their
    ! post-processing setup, then dispatches to the method-specific config
    ! reader based on the selected acoustic solver (config % solver % method).
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)               :: files_path    ! [-] Directory of the .nml files (with trailing '/')
    TYPE(Config_t)  , INTENT(OUT), ALLOCATABLE :: config        ! [-] Global configuration filled in from the .nml files

    ! Local variables
    CHARACTER(len=1024) :: turbine_config_path
    CHARACTER(len=1024) :: farm_config_path
    CHARACTER(len=1024) :: solver_config_path
    CHARACTER(len=1024) :: environment_config_path
    CHARACTER(len=1024) :: computation_config_path
    CHARACTER(len=1024) :: method_config_path


    ALLOCATE(Config_t :: config)

    turbine_config_path     = files_path // "WindTurbine.nml"
    farm_config_path        = files_path // "WindFarm.nml"
    solver_config_path      = files_path // "Solver.nml"
    environment_config_path = files_path // "Environment.nml"

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

    CASE ('FrequencySource')
        CALL read_frequency_source_config(files_path // "FrequencySource.nml", config)
        computation_config_path = files_path // "ComputationSetupSource.nml"

    CASE DEFAULT
        error stop "read_input_files(): Unknown simulation_type (use 'turbine' or 'farm'): " // &
                   trim(config % solver % simulation_type)
    END SELECT

    CALL read_computation_config(computation_config_path, config)

    ! Dispatch to the method-specific config reader based on which
    ! acoustic solver was selected in Solver.nml
    SELECT CASE (trim(config % solver % method))
    CASE ('ImageMethod')
        method_config_path = files_path // "ImagesMethod.nml"
        ALLOCATE(Config_ImageMethod_t :: config % method)
        CALL read_ImageMethod_config(method_config_path, config)

    CASE ('AnalyticalNormalModes')
        method_config_path = files_path // "AnalyticalNormalModes.nml"
        ALLOCATE(Config_AnalyticalNormalModes_t :: config % method)
        CALL read_AnalyticalNormalModes_config(method_config_path, config)

    CASE DEFAULT
        error stop "read_input_files(): Unknown method " // trim(config % solver % method)
    END SELECT

END SUBROUTINE read_input_files

END MODULE Args_parser
