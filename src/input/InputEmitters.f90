! ===============================================================================
! MODULE: InputEmitters
! ===============================================================================
! Defines the readers of the emitter namelists: WindTurbine.nml for a single
! wind turbine (read_turbine_config) and WindFarm.nml for a wind farm
! (read_farm_config).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE InputEmitters

USE Kinds, ONLY: I32, WP
USE ConfigTypes, ONLY: Config_t
USE InputChecks, ONLY: check_time_window, check_fits

IMPLICIT NONE

PRIVATE
PUBLIC :: read_turbine_config, read_farm_config


CONTAINS


SUBROUTINE read_turbine_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "turbine_config" namelist from WindTurbine.nml (rootname,
    ! OpenFAST output paths, wind conditions, placement, node counts, time
    ! window, frequency filtering, etc.) and packs the values into
    ! config % turbine. Strings are read into long buffers and rejected when
    ! they do not fit in Config_WindTurbine_t (no silent truncation), and the
    ! time window and decimation factors are checked before any file is read.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to WindTurbine.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg

    ! Input variables
    CHARACTER(len=1024) :: case_type
    CHARACTER(len=1024) :: rootname
    CHARACTER(len=1024) :: OP_output_dir
    CHARACTER(len=1024) :: OP_output_filename
    CHARACTER(len=1024) :: OP_SubDyn_filename
    INTEGER(I32)        :: Nmembers
    INTEGER(I32)        :: Nnodes
    INTEGER(I32)        :: skip
    INTEGER(I32)        :: skipf
    REAL(WP)            :: WindSpeed
    REAL(WP)            :: WindDir
    REAL(WP)            :: AxisPos(2)
    REAL(WP)            :: BariPos(2)
    REAL(WP)            :: From
    REAL(WP)            :: Upto
    REAL(WP)            :: freqs_over
    LOGICAL             :: Binary
    LOGICAL             :: filter_freqs

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
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=turbine_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &turbine_config: ", trim(iomsg)
        close(unit_num)
        error stop
    end if

    close(unit_num)

    ! --- 4. Validate what Config_t cannot hold or the OpenFAST reader cannot use --- !
    CALL check_fits(case_type, config % turbine % case_type, 'case_type')
    CALL check_fits(rootname, config % turbine % rootname, 'rootname')
    CALL check_fits(OP_output_dir, config % turbine % OP_output_dir, 'OP_output_dir')
    CALL check_fits(OP_output_filename, config % turbine % OP_output_path, 'OP_output_filename')
    CALL check_fits(OP_SubDyn_filename, config % turbine % OP_SubDyn_path, 'OP_SubDyn_filename')
    CALL check_time_window(From, Upto, skip, skipf, 'turbine_config')

    ! --- 5. Pack values into Config_t --- !
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
    ! depth from Environment.nml, which must therefore be read first. Strings
    ! that do not fit in Config_t and invalid time windows or decimation
    ! factors stop the run before any OpenFAST file is read.
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: iso_fortran_env, ONLY: IOSTAT_END

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to WindFarm.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr, n, i, j
    CHARACTER(len=200) :: iomsg

    ! Input variables: shared by the whole farm
    CHARACTER(len=1024) :: farm_name
    REAL(WP)            :: From
    REAL(WP)            :: Upto
    INTEGER(I32)        :: skip
    INTEGER(I32)        :: skipf
    LOGICAL             :: filter_freqs
    REAL(WP)            :: freqs_over

    ! Input variables: one set per turbine
    CHARACTER(len=1024) :: case_type
    CHARACTER(len=1024) :: rootname
    CHARACTER(len=1024) :: OP_output_dir
    CHARACTER(len=1024) :: OP_output_filename
    CHARACTER(len=1024) :: OP_SubDyn_filename
    LOGICAL             :: Binary
    REAL(WP)            :: WindSpeed
    REAL(WP)            :: WindDir
    REAL(WP)            :: AxisPos(2)
    REAL(WP)            :: BariPos(2)
    INTEGER(I32)        :: Nmembers
    INTEGER(I32)        :: Nnodes

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
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=farm_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &farm_config: ", trim(iomsg)
        close(unit_num)
        error stop
    end if

    CALL check_fits(farm_name, config % farm % farm_name, 'farm_name')
    CALL check_time_window(From, Upto, skip, skipf, 'farm_config')

    ! --- 3. Count the &turbine blocks --- !
    rewind(unit_num)
    n = 0_I32
    do
        read(unit_num, nml=turbine, iostat=ierr, iomsg=iomsg)
        if (ierr == IOSTAT_END) exit
        if (ierr /= 0_I32) then
            write(*, '(A,I0,A)') " ERROR: Failed to read &turbine block ", n + 1_I32, ": " // trim(iomsg)
            close(unit_num)
            error stop
        end if
        n = n + 1_I32
    end do

    if (n == 0_I32) then
        write(*,*) "ERROR: No &turbine block found in ", path
        close(unit_num)
        error stop
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
            error stop
        end if

        CALL check_fits(case_type, config % farm % turbines(i) % case_type, 'case_type')
        CALL check_fits(rootname, config % farm % turbines(i) % rootname, 'rootname')
        CALL check_fits(OP_output_dir, config % farm % turbines(i) % OP_output_dir, 'OP_output_dir')
        CALL check_fits(OP_output_filename, config % farm % turbines(i) % OP_output_path, 'OP_output_filename')
        CALL check_fits(OP_SubDyn_filename, config % farm % turbines(i) % OP_SubDyn_path, 'OP_SubDyn_filename')

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

END MODULE InputEmitters
