! ===============================================================================
! MODULE: FrequencySource
! ===============================================================================
! Defines FrequencySource_t, an extension of ManualSources_t whose forces are
! given directly in the frequency domain (FrequencySource.nml): a complex
! force per frequency, source and component. The frequency axis is a list, a
! grid (linear or logarithmic) or the rows of a file, and the forces are
! written in the namelist, read from files or evaluated from text expressions
! of the frequency f [Hz]; init builds both. There is no FFT window (get_enbw
! keeps -1, so no ENBW dataset is written) and no radiation-impedance
! correction, so the forces handed to the solver are the ones of the input.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE FrequencySource

USE CSVReader, ONLY: read_table
USE ConfigTypes, ONLY: Config_t
USE StringUtils, ONLY: to_upper
USE SourceTypes, ONLY: SOURCE_DIPOLE
USE ManualSources, ONLY: ManualSources_t
USE Kinds, ONLY: I32, WP, DEFAULT_I, MAX_PARSE_FREQ, MAX_PARSE_SRC
USE ExpressionParser, ONLY: Expression_t, compile_expression, evaluate_expression

IMPLICIT NONE

PRIVATE
PUBLIC :: FrequencySource_t


CHARACTER(len=2), PARAMETER :: COMPONENT(3) = ['Fx', 'Fy', 'Fz']    ! [-] Names of the force components


TYPE, EXTENDS(ManualSources_t) :: FrequencySource_t
    ! -------------------------------------------------------------------------
    ! Point dipoles with a complex force per frequency. What the solver
    ! receives is in ManualSources_t; this type builds it and keeps how the
    ! input was given, to list it in the .dat summary with the forces.
    ! -------------------------------------------------------------------------

    LOGICAL                         :: source_from_file      ! [-] Whether the positions were read from source_file
    CHARACTER(len=200)              :: source_file           ! [-] File of the positions
    CHARACTER(len=10)               :: freq_mode             ! [-] Frequency axis: 'list', 'grid' or 'file'
    REAL(WP)                        :: fmin                  ! [Hz] First frequency of the grid
    REAL(WP)                        :: fmax                  ! [Hz] Last frequency of the grid
    CHARACTER(len=10)               :: freq_spacing          ! [-] Spacing of the grid: 'linear' or 'log'
    CHARACTER(len=200)              :: freq_file             ! [-] File of the frequencies
    CHARACTER(len=10)               :: force_mode            ! [-] Forces: 'inline', 'expression' or 'file'
    CHARACTER(len=200)              :: force_re_file         ! [-] File of the real part of the forces
    CHARACTER(len=200)              :: force_im_file         ! [-] File of the imaginary part ('None': 0)
    CHARACTER(len=256), ALLOCATABLE :: force_re_expr(:,:)    ! [-] Real part as text expressions of f (Nsources, 3)
    CHARACTER(len=256), ALLOCATABLE :: force_im_expr(:,:)    ! [-] Imaginary part as text expressions of f (Nsources, 3)

    CONTAINS

    ! --- Public type-bound procedures --- !
    PROCEDURE :: init      => init_FrequencySource         ! Builds the frequency axis and the forces of the input
    PROCEDURE :: write_txt => write_txt_FrequencySource    ! Adds the input modes and the force table to the .dat

END TYPE FrequencySource_t


CONTAINS


SUBROUTINE init_FrequencySource(self, config)
    ! -------------------------------------------------------------------------
    ! Builds the emitter from the configuration: the positions, the frequency
    ! axis and the forces, each one from the mode in which it was given.
    ! -------------------------------------------------------------------------

    CLASS(FrequencySource_t), INTENT(INOUT) :: self      ! [-] Emitter instance
    TYPE(Config_t)          , INTENT(IN)    :: config    ! [-] Global configuration

    ! Local variables
    CHARACTER(len=16) :: nf_str


    self % kind_name = 'FrequencySource'
    self % title     = 'FREQUENCY-DOMAIN SOURCES'

    self % Nsources         = config % freq_source % Nsources
    self % xyz              = config % freq_source % source_xyz
    self % source_from_file = config % freq_source % source_from_file
    self % source_file      = config % freq_source % source_file
    self % freq_mode        = config % freq_source % freq_mode
    self % fmin             = config % freq_source % fmin
    self % fmax             = config % freq_source % fmax
    self % freq_spacing     = config % freq_source % freq_spacing
    self % freq_file        = config % freq_source % freq_file
    self % force_mode       = config % freq_source % force_mode
    self % force_re_file    = config % freq_source % force_re_file
    self % force_im_file    = config % freq_source % force_im_file
    self % force_re_expr    = config % freq_source % force_re_expr
    self % force_im_expr    = config % freq_source % force_im_expr

    CALL build_frequency_axis(self, config)
    CALL build_forces(self, config)

    ! Every source radiates as a point dipole
    ALLOCATE(self % order(self % Nsources), source=SOURCE_DIPOLE)

    write(nf_str, '(I0)') size(self % Freqs)
    self % name_suffix = '_Nf' // trim(nf_str)

    write(*, '(A,I0,A,I0,A)') ' >> [FrequencySource] Initialized with ', self % Nsources, ' sources at ', &
                              size(self % Freqs), ' frequencies'

END SUBROUTINE init_FrequencySource


SUBROUTINE build_frequency_axis(self, config)
    ! -------------------------------------------------------------------------
    ! Builds the frequency vector Freqs: the list of the namelist; Nfreqs
    ! frequencies from fmin to fmax with a constant step (linear) or a constant
    ! ratio (log), spaced as run_decay spaces its points; or the rows of
    ! freq_file, which must be Nfreqs if Nfreqs was given. In every mode the
    ! frequencies must be positive and strictly increasing (the reader has
    ! already checked the list and the limits of the grid; rounding or the
    ! values of a file can still break it). Except for the list, the axis is
    ! printed.
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite

    CLASS(FrequencySource_t), INTENT(INOUT) :: self      ! [-] Emitter instance
    TYPE(Config_t)          , INTENT(IN)    :: config    ! [-] Global configuration

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: origin
    REAL(WP)        , ALLOCATABLE :: table(:,:)
    REAL(WP)                      :: step, f_log1
    INTEGER(I32)                  :: i, Nf, first_bad


    SELECT CASE (trim(self % freq_mode))
    CASE ('list')
        self % Freqs = config % freq_source % Freqs
        origin       = 'freqs'

    CASE ('grid')
        Nf = config % freq_source % Nfreqs
        ALLOCATE(self % Freqs(Nf))

        if (trim(self % freq_spacing) == 'log') then
            f_log1 = log10(self % fmin)
            step   = (log10(self % fmax) - f_log1) / real(Nf - 1_I32, WP)
            do i = 1, Nf
                self % Freqs(i) = 10.0_WP**(f_log1 + real(i - 1_I32, WP) * step)
            end do
        else
            step = (self % fmax - self % fmin) / real(Nf - 1_I32, WP)
            do i = 1, Nf
                self % Freqs(i) = self % fmin + real(i - 1_I32, WP) * step
            end do
        end if
        origin = 'the ' // trim(self % freq_spacing) // ' grid from fmin to fmax'

    CASE ('file')
        CALL read_table(self % freq_file, 1_I32, table)

        if (config % freq_source % Nfreqs /= DEFAULT_I .and. size(table, 1) /= config % freq_source % Nfreqs) then
            write(*, '(A,I0,A,I0)') " ERROR FrequencySource: freq_file " // trim(self % freq_file) // " has ", &
                                    size(table, 1), " rows and Nfreqs = ", config % freq_source % Nfreqs
            error stop
        end if

        self % Freqs = table(:, 1)
        origin       = 'freq_file ' // trim(self % freq_file)

    CASE DEFAULT
        error stop "build_frequency_axis: unknown freq_mode " // trim(self % freq_mode)
    END SELECT

    ! Positive and strictly increasing, whatever the mode
    Nf        = size(self % Freqs)
    first_bad = 0_I32
    do i = 1, Nf
        if (.not. (ieee_is_finite(self % Freqs(i)) .and. self % Freqs(i) > 0.0_WP)) first_bad = i
        if (i > 1) then
            if (.not. self % Freqs(i) > self % Freqs(i-1)) first_bad = i
        end if
        if (first_bad > 0_I32) exit
    end do

    if (first_bad > 0_I32) then
        write(*, '(A,I0,A,ES12.5,A)') " ERROR FrequencySource: the frequencies of " // origin // " must be " // &
                                      "positive and strictly increasing, and frequency ", first_bad, " (", &
                                      self % Freqs(first_bad), " Hz) is not"
        error stop
    end if

    if (trim(self % freq_mode) /= 'list') &
        write(*, '(A,I0,A,ES12.5,A,ES12.5,A)') ' >> [FrequencySource] Frequency axis (' // origin // '): ', Nf, &
                                                ' frequencies from ', self % Freqs(1), ' to ', self % Freqs(Nf), ' Hz'

END SUBROUTINE build_frequency_axis


SUBROUTINE build_forces(self, config)
    ! -------------------------------------------------------------------------
    ! Builds the forces Force (Nfreqs, Nsources, 3) of the mode in which they
    ! were given: the values of the namelist; the text expressions of f
    ! evaluated at the frequency vector (real part and imaginary part); or the
    ! columns of force_re_file and force_im_file (Fx, Fy, Fz of source 1, then
    ! those of source 2, ...), one row per frequency, with an imaginary part 0
    ! if there is no force_im_file. A value that is not finite (a division by
    ! zero at one frequency, a logarithm of a negative number) stops the run,
    ! naming the input, the source, the component and the first frequency at
    ! which it happens.
    ! -------------------------------------------------------------------------

    CLASS(FrequencySource_t), INTENT(INOUT) :: self      ! [-] Emitter instance
    TYPE(Config_t)          , INTENT(IN)    :: config    ! [-] Global configuration

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: message
    REAL(WP)        , ALLOCATABLE :: re(:), im(:), re_table(:,:), im_table(:,:)
    TYPE(Expression_t)            :: expr
    INTEGER(I32)                  :: i_s, k, Nf


    Nf = size(self % Freqs)

    SELECT CASE (trim(self % force_mode))
    CASE ('inline')
        self % Force = config % freq_source % Force

    CASE ('expression')
        ALLOCATE(self % Force(Nf, self % Nsources, 3), re(Nf), im(Nf))

        do i_s = 1, self % Nsources
            do k = 1, 3
                CALL compile_expression(self % force_re_expr(i_s, k), expr, message, variable='f')
                if (len(message) > 0) then
                    write(*, '(A,I0,A)') " ERROR FrequencySource: force_re_expr(", i_s, ", " // COMPONENT(k) // &
                                         "): " // message
                    error stop
                end if
                CALL evaluate_expression(expr, self % Freqs, re)
                CALL check_finite_force(re, 'real part (force_re_expr)', i_s, k, self % Freqs)

                CALL compile_expression(self % force_im_expr(i_s, k), expr, message, variable='f')
                if (len(message) > 0) then
                    write(*, '(A,I0,A)') " ERROR FrequencySource: force_im_expr(", i_s, ", " // COMPONENT(k) // &
                                         "): " // message
                    error stop
                end if
                CALL evaluate_expression(expr, self % Freqs, im)
                CALL check_finite_force(im, 'imaginary part (force_im_expr)', i_s, k, self % Freqs)

                self % Force(:, i_s, k) = cmplx(re, im, kind=WP)
            end do
        end do

        write(*, '(A)') ' >> [FrequencySource] Forces evaluated from the expressions of f'

    CASE ('file')
        CALL read_force_table(self % force_re_file, 'force_re_file', self % Nsources, Nf, re_table)

        if (to_upper(trim(self % force_im_file)) /= 'NONE') then
            CALL read_force_table(self % force_im_file, 'force_im_file', self % Nsources, Nf, im_table)
        else
            ALLOCATE(im_table(Nf, 3*self % Nsources), source=0.0_WP)
        end if

        ALLOCATE(self % Force(Nf, self % Nsources, 3))
        do i_s = 1, self % Nsources
            do k = 1, 3
                CALL check_finite_force(re_table(:, 3*(i_s-1)+k), 'real part (force_re_file)', i_s, k, self % Freqs)
                CALL check_finite_force(im_table(:, 3*(i_s-1)+k), 'imaginary part (force_im_file)', i_s, k, &
                                        self % Freqs)
            end do

            self % Force(:, i_s, :) = cmplx(re_table(:, 3*(i_s-1)+1:3*i_s), im_table(:, 3*(i_s-1)+1:3*i_s), kind=WP)
        end do

        if (to_upper(trim(self % force_im_file)) /= 'NONE') then
            write(*, '(A)') ' >> [FrequencySource] Forces read from ' // trim(self % force_re_file) // &
                            ' (real part) and ' // trim(self % force_im_file) // ' (imaginary part)'
        else
            write(*, '(A)') ' >> [FrequencySource] Forces read from ' // trim(self % force_re_file) // &
                            ' (real part; imaginary part 0)'
        end if

    CASE DEFAULT
        error stop "build_forces: unknown force_mode " // trim(self % force_mode)
    END SELECT

END SUBROUTINE build_forces


SUBROUTINE read_force_table(path, name, Nsources, Nf, table)
    ! -------------------------------------------------------------------------
    ! Reads a file of forces, one row per frequency and 3 * Nsources columns
    ! (Fx1 Fy1 Fz1 Fx2 ...), and stops the run if its rows are not the
    ! frequencies of the axis.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)               :: path          ! [-] Path to the file
    CHARACTER(len=*), INTENT(IN)               :: name          ! [-] Namelist variable shown in the error message
    INTEGER(I32)    , INTENT(IN)               :: Nsources      ! [-] Number of sources
    INTEGER(I32)    , INTENT(IN)               :: Nf            ! [-] Number of frequencies of the axis
    REAL(WP)        , INTENT(OUT), ALLOCATABLE :: table(:,:)    ! [N] Forces (Nf, 3 * Nsources)


    CALL read_table(path, 3_I32 * Nsources, table)

    if (size(table, 1) /= Nf) then
        write(*, '(A,I0,A,I0,A)') " ERROR FrequencySource: " // name // " " // trim(path) // " has ", size(table, 1), &
                                  " rows and the frequency axis has ", Nf, " frequencies"
        error stop
    end if

END SUBROUTINE read_force_table


SUBROUTINE check_finite_force(values, part, i_s, k, freqs)
    ! -------------------------------------------------------------------------
    ! Stops the run if a part of the force of one source and component is not
    ! finite at some frequency, naming the first one.
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite

    REAL(WP)        , INTENT(IN) :: values(:)    ! [N] Part of the force at every frequency
    CHARACTER(len=*), INTENT(IN) :: part         ! [-] Which part and input, shown in the error message
    INTEGER(I32)    , INTENT(IN) :: i_s          ! [-] Source
    INTEGER(I32)    , INTENT(IN) :: k            ! [-] Component (1, 2, 3: Fx, Fy, Fz)
    REAL(WP)        , INTENT(IN) :: freqs(:)     ! [Hz] Frequency vector (same size)

    ! Local variables
    INTEGER(I32) :: first_bad


    if (all(ieee_is_finite(values))) return

    first_bad = minloc(merge(1, 0, ieee_is_finite(values)), dim=1)
    write(*, '(A,I0,A,ES12.5,A)') " ERROR FrequencySource: the " // part // " of the force " // COMPONENT(k) // &
                                  " of source ", i_s, " is not finite at f = ", freqs(first_bad), " Hz"
    error stop

END SUBROUTINE check_finite_force


SUBROUTINE write_txt_FrequencySource(self, unit)
    ! -------------------------------------------------------------------------
    ! Writes the sources block of the .dat summary: the title and positions of
    ! ManualSources_t, how the input was given when it was not written in the
    ! namelist (files, grid, expressions), and the forces, so that a change in
    ! any of them is seen by check_update:
    ! - the complex force of every source at every frequency, with nine
    !   significant digits, if the table fits in what the namelist can hold
    !   by hand (MAX_PARSE_FREQ frequencies and MAX_PARSE_SRC sources). When
    !   the frequencies and the forces were written in the namelist, this
    !   table is the input itself, as it always was;
    ! - otherwise, or when they come from a grid, files or expressions (whose
    !   values have any number of digits, and whose table can have
    !   Nfreqs * Nsources lines), checksums at full precision, as TimeSource
    !   does: the sums of the frequencies and of their squares and, for every
    !   source and component, the sums of the real and imaginary parts, of
    !   |F|**2 and of both parts weighted by the frequency index.
    ! -------------------------------------------------------------------------

    CLASS(FrequencySource_t), INTENT(IN) :: self    ! [-] Emitter instance
    INTEGER(I32)            , INTENT(IN) :: unit    ! [-] Unit of the open .dat file

    ! Local variables
    INTEGER(I32)          :: i, i_f, i_s, k, Nf
    COMPLEX(WP)           :: F_is(3)
    REAL(WP), ALLOCATABLE :: weight(:)
    LOGICAL               :: described, fits, by_hand


    CALL self % write_txt_base(unit)

    ! --- 1. How the input was given, if not written in the namelist --- !
    described = .false.

    if (self % source_from_file) then
        write(unit, '(A)') 'Source file : ' // trim(self % source_file)
        described = .true.
    end if

    SELECT CASE (trim(self % freq_mode))
    CASE ('grid')
        write(unit, '(A,ES23.15,A,ES23.15,A)') 'Freq. grid  : ' // trim(self % freq_spacing) // ', from ', &
                                                self % fmin, ' to ', self % fmax, ' Hz'
        described = .true.
    CASE ('file')
        write(unit, '(A)') 'Freq. file  : ' // trim(self % freq_file)
        described = .true.
    END SELECT

    SELECT CASE (trim(self % force_mode))
    CASE ('expression')
        write(unit, '(A)') 'Force [N] as expressions of f [Hz] (those left out are 0)'
        do i_s = 1, self % Nsources
            do k = 1, 3
                if (trim(self % force_re_expr(i_s, k)) /= '0') &
                    write(unit, '(I4,2X,A,2X,A,2X,A)') i_s, COMPONENT(k), 'Re', trim(self % force_re_expr(i_s, k))
                if (trim(self % force_im_expr(i_s, k)) /= '0') &
                    write(unit, '(I4,2X,A,2X,A,2X,A)') i_s, COMPONENT(k), 'Im', trim(self % force_im_expr(i_s, k))
            end do
        end do
        described = .true.
    CASE ('file')
        write(unit, '(A)') 'Force file  : ' // trim(self % force_re_file) // ' (real part)'
        write(unit, '(A)') 'Force file  : ' // trim(self % force_im_file) // ' (imaginary part)'
        described = .true.
    END SELECT

    if (described) write(unit, *)

    ! --- 2. Forces: the full table if it fits, checksums if it does not or if it was not written by hand --- !
    Nf      = size(self % Freqs)
    fits    = Nf <= MAX_PARSE_FREQ .and. self % Nsources <= MAX_PARSE_SRC
    by_hand = trim(self % freq_mode) == 'list' .and. trim(self % force_mode) == 'inline'

    write(unit, '(A)') 'Force [N]: peak complex amplitude, time factor exp(-i*w*t)'
    if (fits) then
        write(unit, '(A)') '   f   s   Freq [Hz]          Re(Fx)          Im(Fx)          Re(Fy)' // &
                           '          Im(Fy)          Re(Fz)          Im(Fz)'
        do i_f = 1, Nf
            do i_s = 1, self % Nsources
                F_is = self % Force(i_f, i_s, :)
                write(unit, '(I4,I4,F12.4,6ES16.8)') i_f, i_s, self % Freqs(i_f), &
                                                     (real(F_is(k), WP), aimag(F_is(k)), k = 1, 3)
            end do
        end do
    end if

    if (.not. (fits .and. by_hand)) then
        if (fits) write(unit, *)
        write(unit, '(A)') 'Checksums of the force table at full precision (i is the frequency index)'
        write(unit, '(A,ES23.15)') 'sum(f)      : ', sum(self % Freqs)
        write(unit, '(A,ES23.15)') 'sum(f**2)   : ', sum(self % Freqs**2)
        write(unit, '(A4,2X,A4,5A23)') '   s', 'comp', 'sum(Re F)', 'sum(Im F)', 'sum(|F|**2)', 'sum(i*Re F)', &
                                       'sum(i*Im F)'

        weight = [(real(i, WP), i = 1_I32, Nf)]
        do i_s = 1, self % Nsources
            do k = 1, 3
                associate(re => real(self % Force(:, i_s, k), WP), im => aimag(self % Force(:, i_s, k)))
                    write(unit, '(I4,2X,A4,5ES23.15)') i_s, COMPONENT(k), sum(re), sum(im), sum(re**2 + im**2), &
                                                       sum(weight * re), sum(weight * im)
                end associate
            end do
        end do
    end if
    write(unit, *)

END SUBROUTINE write_txt_FrequencySource

END MODULE FrequencySource
