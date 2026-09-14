PROGRAM main
    
USE Kinds
USE Types
USE omp_lib
! USE AcousticSource
USE MathUtils, ONLY: format_elapsed
USE WindTurbine, ONLY: WindTurbine_t
USE Args_parser, ONLY: read_input_files
USE AcousticSolver, ONLY: AcousticSolver_t
USE DTU10MWMonopile, ONLY: DTU10MWMonopile_t
USE DTU10MWFloating, ONLY: DTU10MWFloating_t

IMPLICIT NONE

TYPE(Config_t)                       :: config
TYPE(AcousticSource_t)               :: active_sources
CLASS(WindTurbine_t)   , ALLOCATABLE :: active_turbine
CLASS(AcousticSolver_t), ALLOCATABLE :: active_solver

REAL(WP)         :: start_time, elapsed_display
CHARACTER(len=8) :: tag

    WRITE(*,*) "=== Off_coustics_fortran Initialization ==="
    start_time = real(omp_get_wtime(), WP)

    ! 1. Load namelist dinamically
    CALL read_input_files("./input_files/", config)

    ! (Optional: print Debug if config % solver % debug = .true.)
    if (config % solver % debug) then
        write(*,*) "--- Configuration Loaded ---"
        write(*,*) "Acoustic Method:  ", trim(config%solver%method)
        write(*,*) "Data Save Path:   ", trim(config%solver%save_path)
        write(*,*) "Is Running in Cluster: ", config%solver%in_cluster
        write(*,*) "--- Environment base parameters ---"
        write(*,*) "Water depth:      ", config%env%depth_base
        write(*,*) "Speed of sound:   ", config%env%speed_sound_base
        write(*,*) "Density base:     ", config%env%density_base
        write(*,*) "Ref pressure:     ", config%env%reference_pressure
        write(*,*) "--- Turbine physical properties ---"
        write(*,*) "OpenFAST Rootname: ", trim(config%turbine%rootname)
        write(*,*) "OpenFAST Out Dir:  ", trim(config%turbine%OP_output_dir)
        write(*,*) "Is Binary Out:     ", config%turbine%Binary
        write(*,*) "Wind speed:        ", config%turbine%WindSpeed
        write(*,*) "Wind direction:    ", config%turbine%WindDir
        write(*,*) "Turbine depth:     ", config%turbine%Depth
        write(*,*) "Nmembers (max 9):  ", config%turbine%Nmembers
        write(*,*) "Nnodes (max 9):    ", config%turbine%Nnodes
        write(*,*) "Axis Position X,Y: ", config%turbine%AxisPos
        write(*,*) "Baricenter X,Y:    ", config%turbine%BariPos
        write(*,*) "------------------------------------"
    end if

    ! 2. Polimorfic instance for wind turbine
    if (INDEX(config % turbine % rootname, "_mn") > 0) then
        ALLOCATE(DTU10MWMonopile_t :: active_turbine)
        write(*,*) "Instantiating Monopile Turbine..."
    elseif (INDEX(config%turbine%rootname, "_fl") > 0) then
        ALLOCATE(DTU10MWFloating_t :: active_turbine)
        write(*,*) "Instantiating Floating Turbine..."
    end if

    ! 3. Solve structual dynamics
    write(*,*) "Processing turbine dynamics and structural forces..."
    CALL active_turbine % init(config = config % turbine)
    CALL active_turbine % read_input(config = config % turbine, verbose = .true.)
    CALL active_turbine % compute_force(config = config % turbine, verbose=.true.)

    ! 4. Translate physical forces to mathematical acoustic sources
    write(*,*) "Translating physical forces to mathematical acoustic sources..."
    CALL active_turbine % turbine_to_source(config, active_sources)

    DEALLOCATE(active_turbine)      ! Memory cleanup

    ! ! 5. Execute generic solver
    ! write(*,*) "Launching acoustic propagation solver..."
    ! select case(trim(config % solver % method))
    ! case("ImageMethod")
    !     allocate(MethodImages_t :: active_solver)
    ! case("AnalyticalNormalModes")
    !     allocate(AnalyticalNormalModes_t :: active_solver)
    ! end select
    ! !
    ! CALL active_solver % init(active_sources, config % solver, config % solver % env)
    ! CALL active_solver % run_all()

    ! Finalize simulation
    CALL format_elapsed(start_time, elapsed_display, tag)

    write(*,*)
    write(*,*)               "==================================================="
    write(*, '(A, F8.3, A)') '    Total elapsed time = ', elapsed_display, ' (' // trim(tag) // ')'
    write(*,*)               "==================================================="
    write(*,*)               "Simulation finished successfully!"

END PROGRAM main