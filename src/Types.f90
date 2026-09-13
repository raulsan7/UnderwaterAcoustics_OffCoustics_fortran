MODULE TYPEs

USE Kinds

IMPLICIT NONE

PUBLIC


! ========================================
! WIND TURBINE (For input)
! ========================================
TYPE :: Config_WindTurbine_t

    ! --- Input parameters --- !
    CHARACTER(len=200) :: rootname              ! [-] Name without extensions of the OpenFAST output files
    CHARACTER(len=200) :: OP_output_dir         ! [-] Directory with OpenFast output files
    CHARACTER(len=200) :: OP_output_path        ! [-] Openfast output file path
    CHARACTER(len=200) :: OP_SubDyn_path        ! [-] Openfast SubDyn file path
    LOGICAL            :: Binary                ! [-] Wheter OpenFast files are binary (.outb) or not (.out)
    REAL(WP)           :: WindSpeed             ! [m/s] Wind Speed in norm    
    REAL(WP)           :: WindDir               ! [deg] Wind direction 0 deg points to +x axis (anticlockwise from +x)
    REAL(WP)           :: Depth                 ! [m] Water depth
    INTEGER(I32)       :: Nmembers              ! [-] Number of structural members in the OpenFAST model
    INTEGER(I32)       :: Nnodes                ! [-] Number of structural nodes in the OpenFAST model
    REAL(WP)           :: AxisPos(2)            ! [m] Position of the turbine axis in xy plane
    REAL(WP)           :: BariPos(2)            ! [m] Position of the turbine baricenter in xy plane
    REAL(WP)           :: From                  ! [-] Fraction to start reading input
    REAL(WP)           :: Upto                  ! [-] Fraction to stop  reading input
    INTEGER(I32)       :: skip                  ! [-] skips time snapshots as time(::skip)
    LOGICAL            :: filter_freqs          ! [-] Wheter to delete non usefull frequencies
    REAL(WP)           :: freqs_over            ! [-] Frequencies smaller than freqs_over wont be filtered
    LOGICAL            :: is_monopile           ! [-] Wheter turbine is monopile or floating

END TYPE Config_WindTurbine_t


! ========================================
! ENVIRONMENT (For input)
! ========================================
TYPE :: Config_Environment_t

    ! --- Input parameters --- !
    REAL(WP)           :: depth_base            ! [m] Site depth
    REAL(WP)           :: speed_sound_base      ! [m/s] Base speed of sound in fluid
    REAL(WP)           :: density_base          ! [kg/m3] Base fluid density
    REAL(WP)           :: reference_pressure    ! [Pa] Fluid reference pressure
    REAL(WP)           :: attenuation_surface   ! [-] Surface attenuation coefficient
    REAL(WP)           :: attenuation_bottom    ! [-] Seabed  attenuation coefficient

END TYPE Config_Environment_t


! ========================================
! ACOUSTIC SOLVER (For input)
! ========================================
TYPE :: Config_Solver_t

    ! --- Input parameters --- !
    LOGICAL            :: debug                 ! [-] More info for debugin
    LOGICAL            :: in_cluster            ! [-] Wheter we are working on a cluster
    CHARACTER(len=200) :: save_path             ! [-] Where to save results
    CHARACTER(len=50)  :: method                ! [-] Acoustic method to use
    INTEGER(I32)       :: block_size

END TYPE Config_Solver_t


! ========================================
! CONFIGURATION (For input)
! ========================================
TYPE :: Config_t

    TYPE(Config_WindTurbine_t) :: turbine
    TYPE(Config_Environment_t) :: env
    TYPE(Config_Solver_t)      :: solver

END TYPE Config_t


! ========================================
! ACOUSTIC SORCE ABSTRACTION
! ========================================
TYPE:: AcousticSource_t

    REAL(WP)   , ALLOCATABLE :: xyz(:,:)            ! [m] Acoustic source coordinates (Nsources, 3)
    COMPLEX(WP), ALLOCATABLE :: force(:,:,:)        ! [N] Complex dource forces (Nfreqs, Nsources, 3)
    REAL(WP)   , ALLOCATABLE :: Freqs(:)            ! [Hz] Frequency verctor (Nfreqs)

END TYPE AcousticSource_t


! ========================================
! HELPER FOR METHOD OF IMAGES
! ========================================
TYPE :: ImageSystem_t

    REAL(WP)    , ALLOCATABLE :: nodes_pos(:,:)     ! [m] Nodes coordinates (Nnodes, 3)
    COMPLEX(WP) , ALLOCATABLE :: force(:,:,:)       ! [N] Force array (Nfreqs, Nnodes, 3)
    INTEGER(I32), ALLOCATABLE :: BC_all(:)          ! [-] Boundary conditions array

END TYPE ImageSystem_t





END MODULE TYPEs