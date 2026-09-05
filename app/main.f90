! fpm run --flag "-O3 -fopenmp"


PROGRAM main

USE Kinds
USE Types
USE Args_parser
USE DTU10MWMonopile
USE DTU10MWFloating

IMPLICIT NONE

TYPE(Config_t) :: config
TYPE(DTU10MWMonopile_t) :: monopile
TYPE(DTU10MWFloating_t) :: floating
LOGICAL :: verbose = .true.

write(*,*) "=== Off_coustics_fortran Initialization ==="

! 1. Parse configuration in runtime without recompilation
CALL read_input_file("input.nml", config)
CALL monopile % init(config % turbine)
CALL monopile % read_input(config % turbine, verbose)
CALL monopile % compute_force(config % turbine, verbose)

! 2. Validation / Debug block
! if (config % solver % debug) then
!     write(*,*) "--- Configuration Loaded ---"
!     write(*,*) "Acoustic Method:  ", trim(config%solver%method)
!     write(*,*) "Data Save Path:   ", trim(config%solver%save_path)
!     write(*,*) "Is Running in Cluster: ", config%solver%in_cluster
!     write(*,*) "--- Environment base parameters ---"
!     write(*,*) "Water depth:      ", config%env%depth_base
!     write(*,*) "Speed of sound:   ", config%env%speed_sound_base
!     write(*,*) "Density base:     ", config%env%density_base
!     write(*,*) "Ref pressure:     ", config%env%reference_pressure
!     write(*,*) "Surf attenuation: ", config%env%attenuation_surface
!     write(*,*) "Bot attenuation:  ", config%env%attenuation_bottom
!     write(*,*) "--- Turbine physical properties ---"
!     write(*,*) "OpenFAST Rootname: ", trim(config%turbine%rootname)
!     write(*,*) "OpenFAST Out Dir:  ", trim(config%turbine%OP_output_dir)
!     write(*,*) "Is Binary Out:     ", config%turbine%Binary
!     write(*,*) "Wind speed:        ", config%turbine%WindSpeed
!     write(*,*) "Wind direction:    ", config%turbine%WindDir
!     write(*,*) "Turbine depth:     ", config%turbine%Depth
!     write(*,*) "Nmembers (max 9):  ", config%turbine%Nmembers
!     write(*,*) "Nnodes (max 9):    ", config%turbine%Nnodes
!     write(*,*) "Axis Position X,Y: ", config%turbine%AxisPos
!     write(*,*) "Baricenter X,Y:    ", config%turbine%BariPos
!     write(*,*) "------------------------------------"
! end if

! ! 3. Core calculation steps to be integrated
! write(*,*) "Processing turbine dynamics and structural forces..."
! ! [To Do]: call init_turbine(...) and compute_forces(...)

! write(*,*) "Translating forces to acoustic source representations..."
! ! [To Do]: call convert_turbine_to_sources(...)

! write(*,*) "Launching acoustic propagation solver..."
! ! [To Do]: call run_acoustic_simulation(...)

! write(*,*) "Simulation finished successfully!"

END PROGRAM main