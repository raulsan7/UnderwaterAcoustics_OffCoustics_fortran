! ===============================================================================
! MODULE: InputEmitters
! ===============================================================================
! Defines the readers of the emitter namelists:
! - WindTurbine.nml for a single wind turbine (read_turbine_config)
! - WindFarm.nml for a wind farm (read_farm_config)
! - FrequencySource.nml for frequency-domain point sources (read_frequency_source_config).
! - TimeSource.nml for time-domain point sources (read_time_source_config).
! The last two share read_source_file, which reads the positions of any number
! of sources from a text file instead of the namelist.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE InputEmitters

USE Kinds
USE CSVReader, ONLY: read_table
USE ConfigTypes, ONLY: Config_t
USE StringUtils, ONLY: to_upper
USE ExpressionParser, ONLY: Expression_t, compile_expression
USE InputChecks, ONLY: check_time_window, check_fits, check_source_positions, check_file_exists

IMPLICIT NONE

PRIVATE
PUBLIC :: read_turbine_config, read_farm_config, read_frequency_source_config, read_time_source_config


INTEGER(I32), PARAMETER :: MIN_NT = 4_I32    ! [-] Fewest samples of a time axis that compute_rfft can use


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


SUBROUTINE read_frequency_source_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads FrequencySource.nml: the &frequency_source_config block with the
    ! sources, the frequency axis and the complex forces, each of them given in
    ! one of these ways:
    ! - positions: source_xyz written in the namelist (at most MAX_PARSE_SRC
    !   sources), or source_file, one row x, y, z per source (any number);
    ! - frequency axis: the list freqs (at most MAX_PARSE_FREQ), a grid of
    !   Nfreqs frequencies from fmin to fmax (freq_spacing 'linear' or 'log'),
    !   or freq_file, one frequency per row (any number);
    ! - forces: force_re (mandatory) and force_im (optional, 0 by default)
    !   written in the namelist (at most MAX_PARSE_FREQ frequencies and
    !   MAX_PARSE_SRC sources), text expressions of the frequency f [Hz]
    !   (force_re_expr and force_im_expr, '0' by default, one per source and
    !   component), or the files force_re_file and force_im_file (Nfreqs rows,
    !   3 * Nsources columns; no force_im_file means an imaginary part 0).
    ! A way is chosen by giving its values (an entry of a list, a value other
    ! than -999.0, a file name other than 'None' or an expression other than
    ! '0'). Two ways of the same item stop the run; with none, the written
    ! list and forces are taken, as before. The values written in the
    ! namelist are read into buffers of fixed size initialised to NaN, so that
    ! a missing entry, or one given beyond Nsources or Nfreqs, stops the run
    ! instead of being used as zero. The expressions are compiled here to
    ! report a syntax error before any work is done. The grid, the files and
    ! the expressions are turned into the frequency vector and the forces by
    ! FrequencySource % init. The result is packed into config % freq_source.
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, ieee_is_finite, IEEE_QUIET_NAN

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to FrequencySource.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)                  :: unit_num, ierr, i_s, i_f, k, n_expr
    CHARACTER(len=200)            :: iomsg
    CHARACTER(len=:), ALLOCATABLE :: message
    CHARACTER(len=2)              :: component(3) = ['Fx', 'Fy', 'Fz']
    CHARACTER(len=10)             :: freq_mode, force_mode
    LOGICAL                       :: source_from_file, ways(3)
    REAL(WP)                      :: nan
    REAL(WP)        , ALLOCATABLE :: xyz(:,:)
    TYPE(Expression_t)            :: expr

    ! Input variables
    INTEGER(I32)        :: Nsources
    INTEGER(I32)        :: NFreqs
    REAL(WP)            :: Freqs(MAX_PARSE_FREQ)
    REAL(WP)            :: source_xyz(MAX_PARSE_SRC, 3)
    REAL(WP)            :: force_re(MAX_PARSE_FREQ, MAX_PARSE_SRC, 3)
    REAL(WP)            :: force_im(MAX_PARSE_FREQ, MAX_PARSE_SRC, 3)
    REAL(WP)            :: fmin
    REAL(WP)            :: fmax
    CHARACTER(len=1024) :: source_file
    CHARACTER(len=1024) :: freq_spacing
    CHARACTER(len=1024) :: freq_file
    CHARACTER(len=1024) :: force_re_file
    CHARACTER(len=1024) :: force_im_file
    CHARACTER(len=256)  :: force_re_expr(MAX_PARSE_SRC, 3)
    CHARACTER(len=256)  :: force_im_expr(MAX_PARSE_SRC, 3)

    ! --- 1. Configuration namelist --- !
    NAMELIST /frequency_source_config/ &
        Nsources,                      &
        source_xyz,                    &
        source_file,                   &
        NFreqs,                        &
        Freqs,                         &
        fmin,                          &
        fmax,                          &
        freq_spacing,                  &
        freq_file,                     &
        force_re,                      &
        force_im,                      &
        force_re_expr,                 &
        force_im_expr,                 &
        force_re_file,                 &
        force_im_file


    ! --- 2. Defaults: NaN marks an entry that was not given --- !
    nan = ieee_value(0.0_WP, IEEE_QUIET_NAN)

    Nsources      = DEFAULT_I
    Nfreqs        = DEFAULT_I
    Freqs         = nan
    source_xyz    = nan
    force_re      = nan
    force_im      = 0.0_WP
    fmin          = DEFAULT_R
    fmax          = DEFAULT_R
    source_file   = 'None'
    freq_spacing  = 'linear'
    freq_file     = 'None'
    force_re_file = 'None'
    force_im_file = 'None'
    force_re_expr = '0'
    force_im_expr = '0'

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=frequency_source_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &frequency_source_config: ", trim(iomsg)
        close(unit_num)
        error stop
    end if

    close(unit_num)

    ! --- 4. Number of sources: written in the namelist, or the rows of source_file --- !
    CALL check_fits(source_file, config % freq_source % source_file, 'source_file')
    source_from_file = to_upper(trim(source_file)) /= 'NONE'

    if (source_from_file) then
        CALL read_source_file(source_file, source_xyz, Nsources, 'frequency_source_config', xyz)
    else if (Nsources < 1_I32 .or. Nsources > MAX_PARSE_SRC) then
        write(*, '(A,I0,A,I0,A)') "ERROR: &frequency_source_config: Nsources = ", Nsources, &
                                  " is out of bounds [1, ", MAX_PARSE_SRC, "]"
        error stop
    end if

    ! --- 5. Frequency axis: the list freqs, a grid from fmin to fmax, or freq_file --- !
    CALL check_fits(freq_file, config % freq_source % freq_file, 'freq_file')

    ways = [any(ieee_is_finite(Freqs)), fmin /= DEFAULT_R .or. fmax /= DEFAULT_R, to_upper(trim(freq_file)) /= 'NONE']
    if (count(ways) > 1) then
        write(*, '(A)') " ERROR: &frequency_source_config: the frequency axis is given in more than one way (" // &
                        join_given([CHARACTER(len=9) :: 'freqs', 'fmin/fmax', 'freq_file'], ways) // &
                        "): give only one of them"
        error stop
    end if

    freq_mode = 'list'
    if (ways(2)) freq_mode = 'grid'
    if (ways(3)) freq_mode = 'file'

    SELECT CASE (trim(freq_mode))
    CASE ('list')
        if (Nfreqs < 1_I32 .or. Nfreqs > MAX_PARSE_FREQ) then
            write(*, '(A,I0,A,I0,A)') "ERROR: &frequency_source_config: Nfreqs = ", Nfreqs, &
                                      " is out of bounds [1, ", MAX_PARSE_FREQ, "]"
            error stop
        end if

        if (.not. all(ieee_is_finite(Freqs(1:Nfreqs)))) then
            write(*, '(A,I0,A)') " ERROR: &frequency_source_config: Freqs needs ", Nfreqs, " finite values"
            error stop
        end if

        if (any(ieee_is_finite(Freqs(Nfreqs+1:)))) then
            write(*, '(A,I0,A)') " ERROR: &frequency_source_config: Freqs has more than Nfreqs = ", Nfreqs, " values"
            error stop
        end if

        if (any(Freqs(1:Nfreqs) <= 0.0_WP) .or. any(Freqs(2:Nfreqs) <= Freqs(1:Nfreqs-1))) then
            write(*, '(A)') " ERROR: &frequency_source_config: Freqs must be positive and strictly increasing"
            error stop
        end if

    CASE ('grid')
        if (Nfreqs < 2_I32) then
            write(*, '(A,I0,A)') " ERROR: &frequency_source_config: Nfreqs = ", Nfreqs, &
                                 " (a grid from fmin to fmax needs at least 2 frequencies)"
            error stop
        end if

        if (fmin == DEFAULT_R .or. fmax == DEFAULT_R .or. .not. (ieee_is_finite(fmin) .and. ieee_is_finite(fmax))) then
            write(*, '(A)') " ERROR: &frequency_source_config: a grid of frequencies needs both fmin and fmax"
            error stop
        end if

        if (.not. fmin > 0.0_WP) then
            write(*, '(A,ES12.5,A)') " ERROR: &frequency_source_config: fmin = ", fmin, &
                                     " Hz (the frequencies must be positive)"
            error stop
        end if

        if (.not. fmax > fmin) then
            write(*, '(A,ES12.5,A,ES12.5)') " ERROR: &frequency_source_config: fmax must be greater than fmin: ", &
                                            fmax, " <= ", fmin
            error stop
        end if

        SELECT CASE (to_upper(trim(freq_spacing)))
        CASE ('LINEAR')
            freq_spacing = 'linear'
        CASE ('LOG')
            freq_spacing = 'log'
        CASE DEFAULT
            write(*, '(A)') " ERROR: &frequency_source_config: freq_spacing = '" // trim(freq_spacing) // &
                            "' (it must be 'linear' or 'log')"
            error stop
        END SELECT

    CASE ('file')
        CALL check_file_exists(freq_file, 'freq_file', 'frequency_source_config')

        if (Nfreqs /= DEFAULT_I .and. Nfreqs < 1_I32) then
            write(*, '(A,I0,A)') " ERROR: &frequency_source_config: Nfreqs = ", Nfreqs, &
                                 " (with freq_file it is the number of rows of the file, or -999 to take it from it)"
            error stop
        end if
    END SELECT

    if (trim(freq_mode) /= 'grid' .and. to_upper(trim(freq_spacing)) /= 'LINEAR') &
        write(*,*) " WARNING &frequency_source_config: freq_spacing is ignored because fmin and fmax are not given"

    ! --- 6. Source positions: given for every source, inside the water column --- !
    if (source_from_file) then
        CALL check_source_positions(xyz, Nsources, 'frequency_source_config')
    else
        CALL check_source_positions(source_xyz, Nsources, 'frequency_source_config')
        xyz = source_xyz(1:Nsources, :)
    end if

    ! --- 7. Forces: written in the namelist, expressions of f, or files --- !
    CALL check_fits(force_re_file, config % freq_source % force_re_file, 'force_re_file')
    CALL check_fits(force_im_file, config % freq_source % force_im_file, 'force_im_file')

    ways = [any(ieee_is_finite(force_re)) .or. any(force_im /= 0.0_WP),                          &
            any(force_re_expr /= '0') .or. any(force_im_expr /= '0'),                            &
            to_upper(trim(force_re_file)) /= 'NONE' .or. to_upper(trim(force_im_file)) /= 'NONE']
    if (count(ways) > 1) then
        write(*, '(A)') " ERROR: &frequency_source_config: the forces are given in more than one way (" // &
                        join_given([CHARACTER(len=27) :: 'force_re/force_im', 'force_re_expr/force_im_expr', &
                                    'force_re_file/force_im_file'], ways) // "): give only one of them"
        error stop
    end if

    force_mode = 'inline'
    if (ways(2)) force_mode = 'expression'
    if (ways(3)) force_mode = 'file'

    SELECT CASE (trim(force_mode))
    CASE ('inline')
        ! The buffers of the namelist hold MAX_PARSE_FREQ frequencies and MAX_PARSE_SRC sources
        if (Nsources > MAX_PARSE_SRC) then
            write(*, '(A,I0,A,I0,A)') " ERROR: &frequency_source_config: force_re and force_im hold at most ", &
                                      MAX_PARSE_SRC, " sources and source_file has ", Nsources, &
                                      ": give the forces in files"
            error stop
        end if

        if (trim(freq_mode) == 'file' .and. Nfreqs == DEFAULT_I) then
            write(*, '(A)') " ERROR: &frequency_source_config: force_re written in the namelist needs Nfreqs, " // &
                            "the number of rows of freq_file"
            error stop
        end if

        if (Nfreqs > MAX_PARSE_FREQ) then
            write(*, '(A,I0,A,I0,A)') " ERROR: &frequency_source_config: force_re and force_im hold at most ", &
                                      MAX_PARSE_FREQ, " frequencies and Nfreqs = ", Nfreqs, &
                                      ": give the forces in files or as expressions of f"
            error stop
        end if

        ! Real part given for every (frequency, source), nothing beyond Nfreqs or Nsources
        if (any(ieee_is_finite(force_re(Nfreqs+1:, :, :))) .or. any(ieee_is_finite(force_re(:, Nsources+1:, :))) .or. &
            any(force_im(Nfreqs+1:, :, :) /= 0.0_WP) .or. any(force_im(:, Nsources+1:, :) /= 0.0_WP)) then
            write(*, '(A,I0,A,I0)') " ERROR: &frequency_source_config: force_re or force_im have entries beyond " // &
                                    "Nfreqs = ", Nfreqs, " or Nsources = ", Nsources
            error stop
        end if

        do i_s = 1, Nsources
            do i_f = 1, Nfreqs
                if (.not. all(ieee_is_finite(force_re(i_f, i_s, :)))) then
                    write(*, '(A,I0,A,I0,A)') " ERROR: &frequency_source_config: force_re(", i_f, ", ", i_s, &
                                              ", 1:3) is missing or not finite"
                    error stop
                end if

                if (.not. all(ieee_is_finite(force_im(i_f, i_s, :)))) then
                    write(*, '(A,I0,A,I0,A)') " ERROR: &frequency_source_config: force_im(", i_f, ", ", i_s, &
                                              ", 1:3) is not finite"
                    error stop
                end if
            end do
        end do

    CASE ('expression')
        ! The buffers of the namelist hold MAX_PARSE_SRC sources
        if (Nsources > MAX_PARSE_SRC) then
            write(*, '(A,I0,A,I0,A)') " ERROR: &frequency_source_config: force_re_expr and force_im_expr hold at " // &
                                      "most ", MAX_PARSE_SRC, " sources and source_file has ", Nsources, &
                                      ": give the forces in files"
            error stop
        end if

        if (any(force_re_expr(Nsources+1:, :) /= '0') .or. any(force_im_expr(Nsources+1:, :) /= '0')) then
            write(*, '(A,I0,A)') " ERROR: &frequency_source_config: force_re_expr or force_im_expr have entries " // &
                                 "beyond Nsources = ", Nsources, " rows"
            error stop
        end if

        do i_s = 1, Nsources
            do k = 1, 3
                CALL compile_expression(force_re_expr(i_s, k), expr, message, variable='f')
                if (len(message) > 0) then
                    write(*, '(A,I0,A)') " ERROR: &frequency_source_config: force_re_expr(", i_s, ", " // &
                                         component(k) // "): " // message
                    error stop
                end if

                CALL compile_expression(force_im_expr(i_s, k), expr, message, variable='f')
                if (len(message) > 0) then
                    write(*, '(A,I0,A)') " ERROR: &frequency_source_config: force_im_expr(", i_s, ", " // &
                                         component(k) // "): " // message
                    error stop
                end if
            end do
        end do

    CASE ('file')
        if (to_upper(trim(force_re_file)) == 'NONE') then
            write(*, '(A)') " ERROR: &frequency_source_config: force_im_file is given without force_re_file " // &
                            "(the real part)"
            error stop
        end if

        CALL check_file_exists(force_re_file, 'force_re_file', 'frequency_source_config')
        if (to_upper(trim(force_im_file)) /= 'NONE') &
            CALL check_file_exists(force_im_file, 'force_im_file', 'frequency_source_config')
    END SELECT

    ! --- 8. Pack values into Config_t --- !
    n_expr = merge(Nsources, 0_I32, trim(force_mode) == 'expression')

    associate(fs => config % freq_source)
        fs % Nsources         = Nsources
        fs % Nfreqs           = Nfreqs
        fs % source_xyz       = xyz
        fs % source_from_file = source_from_file
        fs % source_file      = trim(source_file)
        fs % freq_mode        = freq_mode
        fs % fmin             = fmin
        fs % fmax             = fmax
        fs % freq_spacing     = trim(freq_spacing)
        fs % freq_file        = trim(freq_file)
        fs % force_mode       = force_mode
        fs % force_re_file    = trim(force_re_file)
        fs % force_im_file    = trim(force_im_file)
        fs % force_re_expr    = force_re_expr(1:n_expr, :)
        fs % force_im_expr    = force_im_expr(1:n_expr, :)

        if (trim(freq_mode) == 'list') fs % Freqs = Freqs(1:Nfreqs)

        if (trim(force_mode) == 'inline') &
            fs % Force = cmplx(force_re(1:Nfreqs, 1:Nsources, :), force_im(1:Nfreqs, 1:Nsources, :), kind=WP)
    end associate

END SUBROUTINE read_frequency_source_config


SUBROUTINE read_time_source_config(path, config)
    ! -------------------------------------------------------------------------
    ! Reads TimeSource.nml: the &time_source_config block with the source
    ! positions (source_xyz written in the namelist, at most MAX_PARSE_SRC
    ! sources, or source_file, one row x, y, z per source, any number), the
    ! time axis (Nt samples of [t_start, t_end), or the times of a file), the
    ! force of every source and component (a text expression of t, default
    ! '0', at most MAX_PARSE_SRC sources, or a file with Nt rows and
    ! 3 * Nsources columns) and the frequencies kept after the FFT (skipf, fmin,
    ! fmax). The values written in the namelist are read into buffers of fixed
    ! size initialised to NaN, so that a position that is missing, or given
    ! beyond Nsources, stops the run. Every expression is compiled here to
    ! report a syntax error before any work is done; the series themselves are
    ! built by TimeSource % init. The result is packed into config % time_source.
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, ieee_is_finite, IEEE_QUIET_NAN

    CHARACTER(len=*), INTENT(IN)    :: path      ! [-] Path to TimeSource.nml
    TYPE(Config_t)  , INTENT(INOUT) :: config    ! [-] Global configuration

    ! Local variables
    INTEGER(I32)                  :: unit_num, ierr, i_s, k, n_expr
    CHARACTER(len=200)            :: iomsg
    CHARACTER(len=:), ALLOCATABLE :: message
    CHARACTER(len=2)              :: component(3) = ['Fx', 'Fy', 'Fz']
    LOGICAL                       :: source_from_file
    REAL(WP)                      :: nan
    REAL(WP)        , ALLOCATABLE :: xyz(:,:)
    TYPE(Expression_t)            :: expr

    ! Input variables
    INTEGER(I32)        :: Nsources
    INTEGER(I32)        :: Nt
    INTEGER(I32)        :: skipf
    REAL(WP)            :: source_xyz(MAX_PARSE_SRC, 3)
    REAL(WP)            :: t_start
    REAL(WP)            :: t_end
    REAL(WP)            :: fmin
    REAL(WP)            :: fmax
    CHARACTER(len=1024) :: source_file
    CHARACTER(len=1024) :: time_file
    CHARACTER(len=1024) :: force_file
    CHARACTER(len=256)  :: force_expr(MAX_PARSE_SRC, 3)

    ! --- 1. Configuration namelist --- !
    NAMELIST /time_source_config/ &
        Nsources,                 &
        source_xyz,               &
        source_file,              &
        time_file,                &
        t_start,                  &
        t_end,                    &
        Nt,                       &
        force_file,               &
        force_expr,               &
        skipf,                    &
        fmin,                     &
        fmax


    ! --- 2. Defaults: NaN marks a position that was not given --- !
    nan = ieee_value(0.0_WP, IEEE_QUIET_NAN)

    Nsources    = DEFAULT_I
    source_xyz  = nan
    source_file = 'None'
    time_file   = 'None'
    t_start    = DEFAULT_R
    t_end      = DEFAULT_R
    Nt         = DEFAULT_I
    force_file = 'None'
    force_expr = '0'
    skipf      = 1_I32
    fmin       = DEFAULT_R
    fmax       = DEFAULT_R

    ! --- 3. Open and read namelist file --- !
    open(newunit=unit_num, file=path, status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Could not open input file: ", path
        error stop
    end if

    rewind(unit_num)
    read(unit_num, nml=time_source_config, iostat=ierr, iomsg=iomsg)
    if (ierr /= 0_I32) then
        write(*,*) "ERROR: Failed to read &time_source_config: ", trim(iomsg)
        close(unit_num)
        error stop
    end if

    close(unit_num)

    ! --- 4. Sources: count and positions, written in the namelist or the rows of source_file --- !
    CALL check_fits(source_file, config % time_source % source_file, 'source_file')
    source_from_file = to_upper(trim(source_file)) /= 'NONE'

    if (source_from_file) then
        CALL read_source_file(source_file, source_xyz, Nsources, 'time_source_config', xyz)
        CALL check_source_positions(xyz, Nsources, 'time_source_config')
    else
        if (Nsources < 1_I32 .or. Nsources > MAX_PARSE_SRC) then
            write(*, '(A,I0,A,I0,A)') " ERROR: &time_source_config: Nsources = ", Nsources, &
                                      " is out of bounds [1, ", MAX_PARSE_SRC, "]"
            error stop
        end if

        CALL check_source_positions(source_xyz, Nsources, 'time_source_config')
        xyz = source_xyz(1:Nsources, :)
    end if

    ! --- 5. Time axis: a uniform one, or the times of a file --- !
    CALL check_fits(time_file, config % time_source % time_file, 'time_file')
    config % time_source % time_from_file = to_upper(trim(time_file)) /= 'NONE'

    if (config % time_source % time_from_file) then
        CALL check_file_exists(time_file, 'time_file', 'time_source_config')

        if (t_start /= DEFAULT_R .or. t_end /= DEFAULT_R .or. Nt /= DEFAULT_I) &
            write(*,*) " WARNING &time_source_config: t_start, t_end and Nt are ignored because time_file is given"
    else
        if (.not. (ieee_is_finite(t_start) .and. ieee_is_finite(t_end)) .or. &
            t_start == DEFAULT_R .or. t_end == DEFAULT_R) then
            write(*,*) " ERROR: &time_source_config: t_start and t_end are needed when time_file is 'None'"
            error stop
        end if

        if (.not. (t_end > t_start)) then
            write(*, '(A,ES12.5,A,ES12.5)') " ERROR: &time_source_config: t_end must be greater than t_start: ", &
                                            t_end, " <= ", t_start
            error stop
        end if

        if (Nt < MIN_NT) then
            write(*, '(A,I0,A,I0,A)') " ERROR: &time_source_config: Nt = ", Nt, " (at least ", MIN_NT, &
                                      " samples are required)"
            error stop
        end if
    end if

    ! --- 6. Forces: expressions of t, or the columns of a file --- !
    CALL check_fits(force_file, config % time_source % force_file, 'force_file')
    config % time_source % force_from_file = to_upper(trim(force_file)) /= 'NONE'

    if (any(force_expr(Nsources+1:, :) /= '0')) then
        write(*, '(A,I0,A)') " ERROR: &time_source_config: force_expr has entries beyond Nsources = ", Nsources, " rows"
        error stop
    end if

    if (config % time_source % force_from_file) then
        CALL check_file_exists(force_file, 'force_file', 'time_source_config')

        ! The rows beyond Nsources are '0' (checked above)
        if (any(force_expr /= '0')) &
            write(*,*) " WARNING &time_source_config: force_expr is ignored because force_file is given"
    else
        ! The buffer of the namelist holds MAX_PARSE_SRC sources
        if (Nsources > MAX_PARSE_SRC) then
            write(*, '(A,I0,A,I0,A)') " ERROR: &time_source_config: force_expr holds at most ", MAX_PARSE_SRC, &
                                      " sources and source_file has ", Nsources, ": give the forces in force_file"
            error stop
        end if

        do i_s = 1, Nsources
            do k = 1, 3
                CALL compile_expression(force_expr(i_s, k), expr, message)
                if (len(message) > 0) then
                    write(*, '(A,I0,A)') " ERROR: &time_source_config: force_expr(", i_s, ", " // component(k) // &
                                         "): " // message
                    error stop
                end if
            end do
        end do

        if (all(force_expr(1:Nsources, :) == '0')) &
            write(*,*) " WARNING &time_source_config: every force_expr is '0', the sources do not radiate"
    end if

    ! --- 7. Frequencies kept after the FFT --- !
    if (skipf < 1_I32) then
        write(*, '(A,I0,A)') " ERROR: &time_source_config: skipf = ", skipf, " (it must be at least 1)"
        error stop
    end if

    if ((fmin /= DEFAULT_R .and. .not. fmin >= 0.0_WP) .or. (fmax /= DEFAULT_R .and. .not. fmax > 0.0_WP)) then
        write(*, '(A,ES12.5,A,ES12.5,A)') " ERROR: &time_source_config: fmin = ", fmin, " and fmax = ", fmax, &
                                          " (they must be positive, or -999.0 for all)"
        error stop
    end if

    if (fmin /= DEFAULT_R .and. fmax /= DEFAULT_R .and. .not. fmax > fmin) then
        write(*, '(A,ES12.5,A,ES12.5)') " ERROR: &time_source_config: fmax must be greater than fmin: ", &
                                        fmax, " <= ", fmin
        error stop
    end if

    ! --- 8. Pack values into Config_t --- !
    n_expr = merge(0_I32, Nsources, config % time_source % force_from_file)

    associate(ts => config % time_source)
        ts % Nsources         = Nsources
        ts % source_from_file = source_from_file
        ts % source_file      = trim(source_file)
        ts % time_file        = trim(time_file)
        ts % t_start          = t_start
        ts % t_end            = t_end
        ts % Nt               = Nt
        ts % force_file       = trim(force_file)
        ts % skipf            = skipf
        ts % fmin             = fmin
        ts % fmax             = fmax
        ts % source_xyz       = xyz
        ts % force_expr       = force_expr(1:n_expr, :)
    end associate

END SUBROUTINE read_time_source_config


SUBROUTINE read_source_file(source_file, source_xyz, Nsources, block_name, xyz)
    ! -------------------------------------------------------------------------
    ! Reads the positions of the sources from source_file, one row x, y, z per
    ! source (blanks, tabs or commas; lines starting with # or ! are
    ! comments), with no limit on their number. Nsources becomes the number of
    ! rows; if it was also given in the namelist, it must be that number. The
    ! positions cannot be written in the namelist (source_xyz) as well. A
    ! coordinate that is not finite stops the run; the rest of the checks are
    ! those of check_source_positions.
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite

    CHARACTER(len=*), INTENT(IN)               :: source_file        ! [-] Path to the file of positions
    REAL(WP)        , INTENT(IN)               :: source_xyz(:,:)    ! [m] Positions of the namelist (NaN: not given)
    INTEGER(I32)    , INTENT(INOUT)            :: Nsources           ! [-] Number of sources: as read, then the rows
    CHARACTER(len=*), INTENT(IN)               :: block_name         ! [-] Namelist group shown in the error messages
    REAL(WP)        , INTENT(OUT), ALLOCATABLE :: xyz(:,:)           ! [m] Positions of the file (Nsources, 3)

    ! Local variables
    INTEGER(I32) :: i_s


    CALL check_file_exists(source_file, 'source_file', block_name)

    if (any(ieee_is_finite(source_xyz))) then
        write(*, '(A)') " ERROR: &" // block_name // ": the positions are given twice, in source_xyz and in source_file"
        error stop
    end if

    CALL read_table(trim(source_file), 3_I32, xyz)

    if (Nsources /= DEFAULT_I .and. Nsources /= size(xyz, 1)) then
        write(*, '(A,I0,A,I0,A)') " ERROR: &" // block_name // ": Nsources = ", Nsources, " and source_file has ", &
                                  size(xyz, 1), " rows"
        error stop
    end if

    Nsources = size(xyz, 1)

    do i_s = 1, Nsources
        if (.not. all(ieee_is_finite(xyz(i_s, :)))) then
            write(*, '(A,I0,A)') " ERROR: &" // block_name // ": the position of source ", i_s, &
                                 " in source_file is not finite"
            error stop
        end if
    end do

END SUBROUTINE read_source_file


FUNCTION join_given(names, given) RESULT(list)
    ! -------------------------------------------------------------------------
    ! Names of the inputs that are given, joined by ' and ', for the message
    ! of an item given in more than one way.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)  :: names(:)    ! [-] Names of the inputs
    LOGICAL         , INTENT(IN)  :: given(:)    ! [-] Whether each one is given (same size)
    CHARACTER(len=:), ALLOCATABLE :: list        ! [-] The names given, joined by ' and '

    ! Local variables
    INTEGER(I32) :: k


    list = ''
    do k = 1, size(names)
        if (.not. given(k)) cycle
        if (len(list) > 0) list = list // ' and '
        list = list // trim(names(k))
    end do

END FUNCTION join_given

END MODULE InputEmitters
