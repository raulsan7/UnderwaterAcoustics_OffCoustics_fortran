! ===============================================================================
! PROGRAM: main
! ===============================================================================
! Entry point of the off_coustics_fortran application. Orchestrates the full
! simulation pipeline: load the .nml configuration, build the emitter (a single
! wind turbine, a whole wind farm or a few frequency-domain point sources) and
! solve its structural dynamics when it has any, translate the physical forces
! into acoustic point-dipole sources, then instantiate and run the requested
! acoustic propagation solver over all post-processing outputs.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! -------------------------------------------------------------------------------
! COMPILE AND RUN: fpm run --flag "-O3 -fopenmp"
! INPUT FILES    : ./input_files/ by default, or the directory given in the
!                  command line: ./off_coustics ./path_to_input_files
!                  (with fpm: fpm run --flag "-O3 -fopenmp" -- ./path_to_input_files)
! ===============================================================================
PROGRAM main

USE Kinds
USE omp_lib
USE Emitter, ONLY: Emitter_t
USE WindFarm, ONLY: WindFarm_t
USE ConfigTypes, ONLY: Config_t
USE Timing, ONLY: format_elapsed
USE SolverFactory, ONLY: new_solver
USE WindTurbine, ONLY: WindTurbine_t
USE TurbineFactory, ONLY: new_turbine
USE SourceTypes, ONLY: AcousticSource_t
USE AcousticSolver, ONLY: AcousticSolver_t
USE FrequencySource, ONLY: FrequencySource_t
USE Args_parser, ONLY: get_input_path, read_input_files

IMPLICIT NONE

TYPE(Config_t)         , ALLOCATABLE :: config
TYPE(AcousticSource_t) , ALLOCATABLE :: active_sources
TYPE(WindFarm_t)       , ALLOCATABLE :: active_farm
TYPE(FrequencySource_t), ALLOCATABLE :: active_freq_source
CLASS(WindTurbine_t)   , ALLOCATABLE :: active_turbine
CLASS(AcousticSolver_t), ALLOCATABLE :: active_solver
CLASS(Emitter_t)       , ALLOCATABLE :: active_emitter
CHARACTER(len=:)       , ALLOCATABLE :: files_path

REAL(WP)         :: start_time, elapsed_display
CHARACTER(len=8) :: tag


    write(*,*) "=== Off_coustics_fortran Initialization ==="
    start_time = real(omp_get_wtime(), WP)

    ! --- 1. Load namelists dynamically --- !
    CALL get_input_path(files_path)
    CALL read_input_files(files_path, config)

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
        CALL active_turbine % to_source(config, active_sources)
        CALL MOVE_ALLOC(active_turbine, active_emitter)

    CASE ('farm')
        ! --- 2. Wind farm instance: one polymorphic turbine per &turbine block --- !
        ALLOCATE(active_farm)
        write(*,*) ">> Instantiating wind farm..."
        CALL active_farm % init(config=config % farm)

        ! --- 3. Solve the structural dynamics of every turbine, one after another --- !
        write(*,*) ">> Processing the dynamics and structural forces of every turbine..."
        CALL active_farm % compute_forces(config=config % farm, verbose=.true.)     ! COntains read_input()

        ! --- 4. Merge every turbine into a single set of acoustic sources --- !
        write(*,*) ">> Merging every turbine into a single set of acoustic sources..."
        CALL active_farm % to_source(config, active_sources)
        CALL MOVE_ALLOC(active_farm, active_emitter)

    CASE ('FrequencySource')
        ! --- 2. FrequencySource instance: a single emitter with user-defined sources --- !
        ALLOCATE(active_freq_source)
        write(*,*) ">> Instantiating FrequencySource..."
        
        ! --- 3. Nothing to solve: copy the forces given in FrequencySource.nml --- !
        CALL active_freq_source % init(config)

        ! --- 4. Translate the forces to mathematical acoustic sources --- !
        write(*,*) ">> Translating frequency-domain forces to mathematical acoustic sources..."
        CALL active_freq_source % to_source(config, active_sources)
        CALL MOVE_ALLOC(active_freq_source, active_emitter)

    END SELECT

    ! --- 5. Execute generic solver --- !
    CALL new_solver(trim(config % solver % method), active_solver)
    write(*,*) ">> Launching acoustic propagation solver: ", trim(config % solver % method)

    ! --- 6. Set acoustic solver (init frees config: the allocated emitter tells the paths apart) --- !
    CALL active_solver % init(active_sources, config)
    CALL active_solver % save_parameters(active_emitter)

    ! --- 7. Run acoustic propagation --- !
    DEALLOCATE(active_emitter)    ! Memory cleanup
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