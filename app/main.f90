! ===============================================================================
! PROGRAM: main
! ===============================================================================
! Entry point of the off_coustics_fortran application. Orchestrates the full
! simulation pipeline: load the .nml configuration, instantiate the requested
! wind turbine type and solve its structural dynamics, translate the physical
! forces into acoustic point-dipole sources, then instantiate and run the
! requested acoustic propagation solver over all post-processing outputs.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
!-------------------------------------------------------------------------------
! COMPILE AND RUN: fpm run --flag "-O3 -fopenmp"
! ===============================================================================
PROGRAM main

USE Kinds
USE Types
USE omp_lib
USE MathUtils, ONLY: format_elapsed
USE ImageMethod, ONLY: ImageMethod_t
USE WindTurbine, ONLY: WindTurbine_t
USE Args_parser, ONLY: read_input_files
USE AcousticSolver, ONLY: AcousticSolver_t
USE DTU10MWMonopile, ONLY: DTU10MWMonopile_t
USE DTU10MWFloating, ONLY: DTU10MWFloating_t

IMPLICIT NONE

TYPE(Config_t)         , ALLOCATABLE :: config
TYPE(AcousticSource_t) , ALLOCATABLE :: active_sources
CLASS(WindTurbine_t)   , ALLOCATABLE :: active_turbine
CLASS(AcousticSolver_t), ALLOCATABLE :: active_solver

REAL(WP)         :: start_time, elapsed_display
CHARACTER(len=8) :: tag

    WRITE(*,*) "=== Off_coustics_fortran Initialization ==="
    start_time = real(omp_get_wtime(), WP)

    ! 1. Load namelist dinamically
    CALL read_input_files("./input_files/", config)


    ! 2. Polymorphic instance for wind turbine
    SELECT CASE (trim(config % turbine % case_type))
    CASE ('DTU10MWMonopile')
        ALLOCATE(DTU10MWMonopile_t :: active_turbine)
    CASE ('DTU10MWFloating')
        ALLOCATE(DTU10MWFloating_t :: active_turbine)
    CASE DEFAULT
        write(*,*) 'Unsupported turbine type: ' // trim(config % turbine % case_type)
        stop
    END SELECT
    write(*,*) '>> Instantiating ', trim(config % turbine % case_type), ' Turbine...'


    ! 3. Solve structual dynamics
    write(*,*) ">> Processing turbine dynamics and structural forces..."
    CALL active_turbine % init(config = config % turbine)
    CALL active_turbine % read_input(config = config % turbine, verbose = .true.)
    CALL active_turbine % compute_force(config = config % turbine, verbose=.true.)


    ! 4. Translate physical forces to mathematical acoustic sources
    write(*,*) ">> Translating physical forces to mathematical acoustic sources..."
    CALL active_turbine % turbine_to_source(config, active_sources)

    
    ! 5. Execute generic solver
    SELECT CASE(trim(config % solver % method))
    CASE("ImageMethod")
        ALLOCATE(ImageMethod_t :: active_solver)
    CASE("AnalyticalNormalModes")
        ! ALLOCATE(AnalyticalNormalModes_t :: active_solver)
        error stop "Not implemented yet"
    END SELECT
    write(*,*) ">> Launching acoustic propagation solver: ", trim(config % solver % method)
    

    ! 6. Set acoustic solver
    CALL active_solver % init(active_sources, config)
    CALL active_solver % save_parameters(active_turbine)


    ! 7. Run acoustic propagation
    DEALLOCATE(active_turbine)      ! Memory cleanup
    CALL active_solver % run_all()


    ! Finalize simulation
    CALL format_elapsed(start_time, elapsed_display, tag)

    write(*,*)
    write(*,*)               "==================================================="
    write(*, '(A, F8.3, A)') '    Total elapsed time = ', elapsed_display, ' (' // trim(tag) // ')'
    write(*,*)               "==================================================="
    write(*,*)               "Simulation finished successfully!"

END PROGRAM main