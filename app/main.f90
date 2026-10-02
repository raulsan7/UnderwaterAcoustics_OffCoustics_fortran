! ===============================================================================
! PROGRAM: main
! ===============================================================================
! Entry point of the off_coustics_fortran application. Orchestrates the full
! simulation pipeline: load the .nml configuration, build the emitter (a single
! wind turbine or a whole wind farm) and solve its structural dynamics,
! translate the physical forces into acoustic point-dipole sources, then
! instantiate and run the requested acoustic propagation solver over all
! post-processing outputs.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! -------------------------------------------------------------------------------
! COMPILE AND RUN: fpm run --flag "-O3 -fopenmp"
! ===============================================================================
PROGRAM main

USE Kinds
USE omp_lib
USE WindFarm, ONLY: WindFarm_t
USE ConfigTypes, ONLY: Config_t
USE Timing, ONLY: format_elapsed
USE ImageMethod, ONLY: ImageMethod_t
USE WindTurbine, ONLY: WindTurbine_t
USE TurbineFactory, ONLY: new_turbine
USE SourceTypes, ONLY: AcousticSource_t
USE Args_parser, ONLY: read_input_files
USE AcousticSolver, ONLY: AcousticSolver_t
USE AnalyticalNormalModes, ONLY: AnalyticalNormalModes_t

IMPLICIT NONE

TYPE(Config_t)         , ALLOCATABLE :: config
TYPE(AcousticSource_t) , ALLOCATABLE :: active_sources
TYPE(WindFarm_t)       , ALLOCATABLE :: active_farm
CLASS(WindTurbine_t)   , ALLOCATABLE :: active_turbine
CLASS(AcousticSolver_t), ALLOCATABLE :: active_solver

REAL(WP)         :: start_time, elapsed_display
CHARACTER(len=8) :: tag


    write(*,*) "=== Off_coustics_fortran Initialization ==="
    start_time = real(omp_get_wtime(), WP)

    ! --- 1. Load namelists dynamically --- !
    CALL read_input_files("./input_files/", config)

    SELECT CASE (trim(config % solver % simulation_type))
    CASE ('turbine')
        ! --- 2. Polymorphic instance for wind turbine --- !
        CALL new_turbine(config % turbine % case_type, active_turbine)
        write(*,*) '>> Instantiating ', trim(config % turbine % case_type), ' Turbine...'

        ! --- 3. Solve structural dynamics --- !
        write(*,*) ">> Processing turbine dynamics and structural forces..."
        CALL active_turbine % init(config=config % turbine)
        CALL active_turbine % read_input(config=config % turbine, verbose=.true.)
        CALL active_turbine % compute_force(config=config % turbine, verbose=.true.)

        ! --- 4. Translate physical forces to mathematical acoustic sources --- !
        write(*,*) ">> Translating physical forces to mathematical acoustic sources..."
        CALL active_turbine % turbine_to_source(config, active_sources)

    CASE ('farm')
        ! --- 2. Wind farm instance: one polymorphic turbine per &turbine block --- !
        ALLOCATE(active_farm)
        write(*,*) ">> Instantiating wind farm..."
        CALL active_farm % init(config=config % farm)

        ! --- 3. Solve the structural dynamics of every turbine, one after another --- !
        write(*,*) ">> Processing the dynamics and structural forces of every turbine..."
        CALL active_farm % compute_forces(config=config % farm, verbose=.true.)

        ! --- 4. Merge every turbine into a single set of acoustic sources --- !
        write(*,*) ">> Merging every turbine into a single set of acoustic sources..."
        CALL active_farm % farm_to_source(config, active_sources)
    END SELECT

    ! --- 5. Execute generic solver --- !
    SELECT CASE (trim(config % solver % method))
    CASE ("ImageMethod")
        ALLOCATE(ImageMethod_t :: active_solver)
    CASE ("AnalyticalNormalModes")
        ALLOCATE(AnalyticalNormalModes_t :: active_solver)
    END SELECT
    write(*,*) ">> Launching acoustic propagation solver: ", trim(config % solver % method)

    ! --- 6. Set acoustic solver (init frees config: active_farm tells both paths apart) --- !
    CALL active_solver % init(active_sources, config)
    if (ALLOCATED(active_farm)) then
        CALL active_solver % save_parameters(farm=active_farm)
    else
        CALL active_solver % save_parameters(turbine=active_turbine)
    end if

    ! --- 7. Run acoustic propagation --- !
    if (ALLOCATED(active_turbine)) DEALLOCATE(active_turbine)    ! Memory cleanup
    if (ALLOCATED(active_farm))    DEALLOCATE(active_farm)
    CALL active_solver % run()
    CALL active_solver % update_parameters()

    ! --- 8. Finalize simulation --- !
    CALL format_elapsed(start_time, elapsed_display, tag)

    write(*,*)
    write(*,*)               "==================================================="
    write(*, '(A, F8.3, A)') '    Total elapsed time = ', elapsed_display, ' (' // trim(tag) // ')'
    write(*,*)               "==================================================="
    write(*,*)               "Simulation finished successfully!"

END PROGRAM main