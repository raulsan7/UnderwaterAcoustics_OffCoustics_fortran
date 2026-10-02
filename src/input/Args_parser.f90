! ===============================================================================
! MODULE: Args_parser
! ===============================================================================
! Reads the .nml configuration files that describe a simulation case
! (Solver.nml, Environment.nml, WindTurbine.nml or WindFarm.nml,
! ComputationSetupTurbine.nml or ComputationSetupFarm.nml, and
! ImagesMethod.nml or AnalyticalNormalModes.nml) and packs their
! contents into the global Config_t structure defined in
! ConfigTypes.f90. read_input_files is the single entry point: it locates the
! config files inside the given directory and dispatches to one reader
! subroutine per namelist (InputSolver, InputEmitters and
! InputPostProcessing), plus one extra reader chosen according to the
! selected acoustic solver method.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE Args_parser

USE InputPostProcessing, ONLY: read_computation_config
USE InputEmitters, ONLY: read_turbine_config, read_farm_config
USE ConfigTypes, ONLY: Config_ImageMethod_t, Config_AnalyticalNormalModes_t, Config_t
USE InputSolver, ONLY: read_solver_config, read_environment_config, read_ImageMethod_config, &
                       read_AnalyticalNormalModes_config

IMPLICIT NONE

PRIVATE
PUBLIC :: read_input_files


CONTAINS


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
    CHARACTER(len=200) :: turbine_config_path
    CHARACTER(len=200) :: farm_config_path
    CHARACTER(len=200) :: solver_config_path
    CHARACTER(len=200) :: environment_config_path
    CHARACTER(len=200) :: computation_config_path

    ! Method-specific configuration paths
    CHARACTER(len=200) :: ImagesMethod_config_path
    CHARACTER(len=200) :: AnalyticalNormalModes_config_path


    ALLOCATE(Config_t :: config)

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
        error stop "read_input_files(): Unknown simulation_type (use 'turbine' or 'farm'): " // &
                   trim(config % solver % simulation_type)
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
        error stop "read_input_files(): Unknown method " // trim(config % solver % method)
    END SELECT

END SUBROUTINE read_input_files

END MODULE Args_parser
