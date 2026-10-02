! ===============================================================================
! MODULE: InputPostProcessing
! ===============================================================================
! Defines the readers of the post-processing namelists of
! ComputationSetupTurbine.nml or ComputationSetupFarm.nml: the entry point
! read_computation_config, the do_config switches (read_do_config), and one
! reader per run_* step (read_spectrums_config ... read_sphere_config).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE InputPostProcessing

USE ConfigTypes, ONLY: Config_t
USE InputChecks, ONLY: check_all_or_none, check_count
USE Kinds, ONLY: I32, WP, DEFAULT_R, DEFAULT_I, MAX_PARSE_OBS

IMPLICIT NONE

PRIVATE
PUBLIC :: read_computation_config


CONTAINS


SUBROUTINE read_computation_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads every namelist group in ComputationSetupTurbine.nml or
    ! ComputationSetupFarm.nml (do_config plus one group per run_*
    ! post-processing step) and packs the values into config % run_inputs.
    ! Only the groups of the steps that will actually run are read. Every
    ! observer count must be -999 (auto) or at least 1, and the vectors whose
    ! -999.0 sentinel applies to the whole vector (centres, x/y limits, decay
    ! distances) must be given in full or left at -999.0; p1_line, p2_line,
    ! zlim_xz and zlim_sv resolve each component on its own.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml or ComputationSetupFarm.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration


    CALL read_do_config(path, config)

    if (config % run_inputs % do_run % do_spectrums) CALL read_spectrums_config(path, config)

    if (config % run_inputs % do_run % do_polar) then
        CALL read_polar_config(path, config)
        CALL check_count(config % run_inputs % params_polar % n_polar, 'n_polar')
        CALL check_all_or_none(config % run_inputs % params_polar % center_polar, 'center_polar')
    end if

    if (config % run_inputs % do_run % do_cylinder) then
        CALL read_cylinder_config(path, config)
        CALL check_count(config % run_inputs % params_cylinder % n_cyl, 'n_cyl')
        CALL check_count(config % run_inputs % params_cylinder % nz_cyl, 'nz_cyl')
        CALL check_all_or_none(config % run_inputs % params_cylinder % center_cyl, 'center_cyl')
    end if

    if (config % run_inputs % do_run % do_decay) then
        CALL read_decay_config(path, config)
        CALL check_count(config % run_inputs % params_decay % n_decay, 'n_decay')
        CALL check_all_or_none(config % run_inputs % params_decay % distance_decay, 'distance_decay')
    end if

    if (config % run_inputs % do_run % do_line) then
        CALL read_line_config(path, config)
        CALL check_count(config % run_inputs % params_line % n_line, 'n_line')
    end if

    if (config % run_inputs % do_run % do_sliceXY) then
        CALL read_sliceXY_config(path, config)
        CALL check_count(config % run_inputs % params_sliceXY % nx_xy, 'nx_xy')
        CALL check_count(config % run_inputs % params_sliceXY % ny_xy, 'ny_xy')
        CALL check_all_or_none(config % run_inputs % params_sliceXY % xlim_xy, 'xlim_xy')
        CALL check_all_or_none(config % run_inputs % params_sliceXY % ylim_xy, 'ylim_xy')
        CALL check_all_or_none(config % run_inputs % params_sliceXY % center_xy, 'center_xy')
    end if

    if (config % run_inputs % do_run % do_sliceXZ) then
        CALL read_sliceXZ_config(path, config)
        CALL check_count(config % run_inputs % params_sliceXZ % nx_xz, 'nx_xz')
        CALL check_count(config % run_inputs % params_sliceXZ % nz_xz, 'nz_xz')
        CALL check_all_or_none(config % run_inputs % params_sliceXZ % xlim_xz, 'xlim_xz')
    end if

    if (config % run_inputs % do_run % do_sliceVertical) then
        CALL read_sliceVertical_config(path, config)
        CALL check_count(config % run_inputs % params_sliceVertical % nu_sv, 'nu_sv')
        CALL check_count(config % run_inputs % params_sliceVertical % nz_sv, 'nz_sv')
        CALL check_all_or_none(config % run_inputs % params_sliceVertical % center_sv, 'center_sv')
    end if

    if (config % run_inputs % do_run % do_sphere) then
        CALL read_sphere_config(path, config)
        CALL check_count(config % run_inputs % params_sphere % n_theta_sphere, 'n_theta_sphere')
        CALL check_count(config % run_inputs % params_sphere % nz_sphere, 'nz_sphere')
        CALL check_all_or_none(config % run_inputs % params_sphere % center_sphere, 'center_sphere')
    end if

END SUBROUTINE read_computation_config


SUBROUTINE read_do_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "do_config" namelist from ComputationSetupTurbine.nml or
    ! ComputationSetupFarm.nml and packs the post-processing execution
    ! switches into config % run_inputs. When do_all is enabled, every
    ! post-processing switch is forced to true. In farm mode the steps that
    ! make no sense for a wind farm (sphere, decay) default to false and are
    ! then always switched off, so run() never needs to know what produced
    ! the sources; the warning is only printed when the file requests them.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml or ComputationSetupFarm.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration (solver settings already read)

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr
    CHARACTER(len=200) :: iomsg
    LOGICAL            :: is_farm

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


    ! --- 2. Defaults (sphere and decay are only requested by default for a single turbine) --- !
    is_farm          = (trim(config % solver % simulation_type) == 'farm')
    do_all           = .true.
    do_spectrums     = .true.
    do_polar         = .true.
    do_cylinder      = .true.
    do_decay         = .not. is_farm
    do_line          = .true.
    do_sliceXY       = .true.
    do_sliceXZ       = .true.
    do_sliceVertical = .false.    ! Off unless requested: files without this switch do not run the step
    do_sphere        = .not. is_farm

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=do_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &do_config: ", trim(iomsg)
        close(unit_num)
        error stop
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
    if (is_farm) then
        if (.not. do_all .and. (do_sphere .or. do_decay)) &
            write(*,*) "WARNING read_do_config: sphere and decay are not available for a wind farm, skipped"
        config % run_inputs % do_run % do_sphere = .false.
        config % run_inputs % do_run % do_decay  = .false.
    end if

END SUBROUTINE read_do_config


SUBROUTINE read_spectrums_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "spectrums_config" namelist from ComputationSetupTurbine.nml
    ! or ComputationSetupFarm.nml (number of spectrum observer points and
    ! their coordinates) and packs the values into config % run_inputs.
    ! n_observers = -999 selects the 3 automatic points; any other value must
    ! lie in 1..MAX_PARSE_OBS and come with its n_observers (x, y, z)
    ! triplets. The coordinate buffers start as NaN, so a missing coordinate
    ! is detected without mistaking a real -999.0 m for the sentinel.
    ! Coordinates listed beyond n_observers are ignored with a warning.
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, ieee_is_nan, IEEE_QUIET_NAN

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml or ComputationSetupFarm.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)       :: unit_num, ierr, n
    CHARACTER(len=200) :: iomsg
    REAL(WP)           :: nan

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


    ! --- 2. Defaults (NaN marks the coordinates not given in the file) --- !
    nan         = ieee_value(1.0_WP, IEEE_QUIET_NAN)
    n_observers = DEFAULT_I
    observers_x = nan
    observers_y = nan
    observers_z = nan

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=spectrums_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &spectrums_config: ", trim(iomsg)
        close(unit_num)
        error stop
    end if

    close(unit_num)

    ! --- 4. Validate: n_observers in range and one (x, y, z) triplet per observer --- !
    if (n_observers > MAX_PARSE_OBS) then
        write(*,*) "ERROR read_spectrums_config: n_observers exceeds the internal parsing cap of ", MAX_PARSE_OBS
        error stop
    end if

    if (n_observers /= DEFAULT_I) then
        if (n_observers < 1_I32) then
            write(*, '(A,I0,A)') " ERROR read_spectrums_config: n_observers = ", n_observers, &
                                 " (use -999 for the automatic points, or list at least one observer)"
            error stop
        end if

        n = n_observers
        if (any(ieee_is_nan(observers_x(1:n))) .or. any(ieee_is_nan(observers_y(1:n))) .or. &
            any(ieee_is_nan(observers_z(1:n)))) then
            write(*, '(A,I0,A)') " ERROR read_spectrums_config: n_observers = ", n, &
                                 " but fewer x, y or z coordinates are listed"
            error stop
        end if

        if (.not. (all(ieee_is_nan(observers_x(n+1:))) .and. all(ieee_is_nan(observers_y(n+1:))) .and. &
                   all(ieee_is_nan(observers_z(n+1:))))) &
            write(*, '(A,I0,A)') " WARNING read_spectrums_config: coordinates listed beyond n_observers = ", n, &
                                 " are ignored"
    end if

    ! --- 5. Pack values into Config_t --- !
    config % run_inputs % params_spectrums % n_observers = n_observers

    if (n_observers /= DEFAULT_I) then
        config % run_inputs % params_spectrums % observers_x = observers_x(1:n_observers)
        config % run_inputs % params_spectrums % observers_y = observers_y(1:n_observers)
        config % run_inputs % params_spectrums % observers_z = observers_z(1:n_observers)
    end if

END SUBROUTINE read_spectrums_config


SUBROUTINE read_polar_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "polar_config" namelist from ComputationSetupTurbine.nml or
    ! ComputationSetupFarm.nml (circle radius, depth, observer count and
    ! centre) and packs the values into config % run_inputs.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml or ComputationSetupFarm.nml
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
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=polar_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &polar_config: ", trim(iomsg)
        close(unit_num)
        error stop
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
    ! or ComputationSetupFarm.nml (cylinder radius, azimuthal and vertical
    ! observer counts, and centre) and packs the values into
    ! config % run_inputs.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml or ComputationSetupFarm.nml
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
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=cylinder_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &cylinder_config: ", trim(iomsg)
        close(unit_num)
        error stop
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
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=decay_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &decay_config: ", trim(iomsg)
        close(unit_num)
        error stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    config % run_inputs % params_decay % distance_decay = distance_decay
    config % run_inputs % params_decay % n_decay        = n_decay
    config % run_inputs % params_decay % z_decay        = z_decay
    config % run_inputs % params_decay % logspace_decay = logspace_decay

END SUBROUTINE read_decay_config


SUBROUTINE read_line_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads the "line_config" namelist from ComputationSetupTurbine.nml or
    ! ComputationSetupFarm.nml (line endpoints, observer count and spacing
    ! mode) and packs the values into config % run_inputs.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml or ComputationSetupFarm.nml
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
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=line_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &line_config: ", trim(iomsg)
        close(unit_num)
        error stop
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
    ! or ComputationSetupFarm.nml (slice depth, grid dimensions, extents and
    ! centre) and packs the values into config % run_inputs.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml or ComputationSetupFarm.nml
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
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=sliceXY_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &sliceXY_config: ", trim(iomsg)
        close(unit_num)
        error stop
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
    ! or ComputationSetupFarm.nml (slice position, grid dimensions and x/z
    ! extents) and packs the values into config % run_inputs.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to ComputationSetupTurbine.nml or ComputationSetupFarm.nml
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
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=sliceXZ_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &sliceXZ_config: ", trim(iomsg)
        close(unit_num)
        error stop
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
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=sliceVertical_config, iostat=ierr, iomsg=iomsg)
    if (ierr == IOSTAT_END) then
        write(*,*) "WARNING read_sliceVertical_config: no &sliceVertical_config group, automatic values used"
    else if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &sliceVertical_config: ", trim(iomsg)
        close(unit_num)
        error stop
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
    REAL(WP)     :: r_sphere, center_sphere(3), f_farfield_sphere
    INTEGER(I32) :: n_theta_sphere, nz_sphere

    ! --- 1. Configuration namelist --- !
    NAMELIST /sphere_config/ &
        r_sphere,            &
        n_theta_sphere,      &
        nz_sphere,           &
        f_farfield_sphere,   &
        center_sphere


    ! --- 2. Defaults --- !
    r_sphere          = DEFAULT_R
    n_theta_sphere    = DEFAULT_I
    nz_sphere         = DEFAULT_I
    f_farfield_sphere = DEFAULT_R
    center_sphere     = DEFAULT_R

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=sphere_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &sphere_config: ", trim(iomsg)
        close(unit_num)
        error stop
    end if

    close(unit_num)

    ! --- 4. Pack values into Config_t --- !
    config % run_inputs % params_sphere % r_sphere          = r_sphere
    config % run_inputs % params_sphere % n_theta_sphere    = n_theta_sphere
    config % run_inputs % params_sphere % nz_sphere         = nz_sphere
    config % run_inputs % params_sphere % f_farfield_sphere = f_farfield_sphere
    config % run_inputs % params_sphere % center_sphere     = center_sphere

END SUBROUTINE read_sphere_config

END MODULE InputPostProcessing
