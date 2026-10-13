! ===============================================================================
! MODULE: TimeSource
! ===============================================================================
! Defines TimeSource_t, an extension of ManualSources_t whose sources are
! given in the time domain (TimeSource.nml): a real series of every source
! and component, the force of a dipole or the volume velocity Q of a
! monopole (first component only), read from a file or evaluated from a text
! expression of t. The solver works in the frequency domain, so init takes
! the series there with compute_rfft, the same routine that does it for the
! accelerations of a turbine: mean removed, Hann window, one-sided spectrum.
! The amplitudes are therefore windowed and get_enbw returns HANN_ENBW, and a
! tone is read with its full amplitude only if its frequency is a multiple of
! the frequency resolution 1 / (record length).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE TimeSource

USE Kinds, ONLY: I32, WP, DEFAULT_R
USE CSVReader, ONLY: read_table
USE ConfigTypes, ONLY: Config_t
USE HDF5Utils, ONLY: save_variable
USE SourceTypes, ONLY: SOURCE_MONOPOLE
USE ManualSources, ONLY: ManualSources_t
USE SignalProcessing, ONLY: compute_rfft, HANN_ENBW
USE ExpressionParser, ONLY: Expression_t, compile_expression, evaluate_expression

IMPLICIT NONE

PRIVATE
PUBLIC :: TimeSource_t


TYPE, EXTENDS(ManualSources_t) :: TimeSource_t
    ! -------------------------------------------------------------------------
    ! Point monopoles and dipoles with a series in time (force, or volume
    ! velocity Q of a monopole). The series and the input that produced them
    ! are kept, besides the spectrum of ManualSources_t, to be saved and
    ! listed in the .dat summary.
    ! -------------------------------------------------------------------------

    INTEGER(I32)                    :: Nt                 ! [-] Number of time samples
    INTEGER(I32)                    :: skipf              ! [-] Keep 1 out of every skipf frequency bins
    REAL(WP)                        :: dt                 ! [s] Time step
    REAL(WP)                        :: fmin               ! [Hz] Lowest frequency kept (DEFAULT_R: all)
    REAL(WP)                        :: fmax               ! [Hz] Highest frequency kept (DEFAULT_R: Nyquist)
    REAL(WP)          , ALLOCATABLE :: time(:)            ! [s] Time axis (Nt)
    REAL(WP)          , ALLOCATABLE :: Ft(:,:,:)          ! [N] Force series (Nt, Nsources, 3)
    LOGICAL                         :: source_from_file   ! [-] Whether the positions were read from source_file
    LOGICAL                         :: time_from_file     ! [-] Whether the time axis was read from time_file
    LOGICAL                         :: force_from_file    ! [-] Whether the forces were read from force_file
    CHARACTER(len=200)              :: source_file        ! [-] File of the positions
    CHARACTER(len=200)              :: time_file          ! [-] File of the time axis
    CHARACTER(len=200)              :: force_file         ! [-] File of the forces
    CHARACTER(len=256), ALLOCATABLE :: force_expr(:,:)    ! [-] Text expressions of the forces (Nsources, 3)

    CONTAINS

    ! --- Public type-bound procedures --- !
    PROCEDURE :: init        => init_TimeSource           ! Builds the series and takes them to the frequency domain
    PROCEDURE :: save_params => save_params_TimeSource    ! Adds the series to the HDF5 output
    PROCEDURE :: write_txt   => write_txt_TimeSource      ! Adds the time axis and the series to the .dat
    PROCEDURE :: get_enbw    => get_enbw_TimeSource       ! Hann window of compute_rfft (HANN_ENBW)

END TYPE TimeSource_t


CONTAINS


SUBROUTINE init_TimeSource(self, config)
    ! -------------------------------------------------------------------------
    ! Builds the emitter from the configuration: the multipole orders, the
    ! time axis (printing its time step), the series and their spectrum, which
    ! is what the solver receives. A monopole with something in the second or
    ! third component of its series stops the run; this is checked on the
    ! series, before the FFT removes their mean.
    ! -------------------------------------------------------------------------

    CLASS(TimeSource_t), INTENT(INOUT) :: self      ! [-] Emitter instance
    TYPE(Config_t)     , INTENT(IN)    :: config    ! [-] Global configuration

    ! Local variables
    CHARACTER(len=16) :: nt_str


    self % kind_name = 'TimeSource'
    self % title     = 'TIME-DOMAIN SOURCES'

    self % Nsources         = config % time_source % Nsources
    self % order            = config % time_source % source_order
    self % xyz              = config % time_source % source_xyz
    self % skipf            = config % time_source % skipf
    self % fmin             = config % time_source % fmin
    self % fmax             = config % time_source % fmax
    self % source_from_file = config % time_source % source_from_file
    self % source_file      = config % time_source % source_file
    self % time_from_file   = config % time_source % time_from_file
    self % force_from_file  = config % time_source % force_from_file
    self % time_file        = config % time_source % time_file
    self % force_file       = config % time_source % force_file
    self % force_expr       = config % time_source % force_expr

    CALL build_time_axis(self, config)
    CALL build_force_series(self)
    CALL self % check_monopole_components(any(self % Ft /= 0.0_WP, dim=1), 'force series')
    CALL take_to_frequency_domain(self)

    write(nt_str, '(I0)') self % Nt
    self % name_suffix = '_Nt' // trim(nt_str)

END SUBROUTINE init_TimeSource


SUBROUTINE build_time_axis(self, config)
    ! -------------------------------------------------------------------------
    ! Builds the time axis: Nt samples of [t_start, t_end), or the times of a
    ! file, which must be uniformly spaced because the FFT needs a constant
    ! time step. It prints the time step, the record length, the frequency
    ! resolution and the Nyquist frequency, which the user needs to choose
    ! the tones and the sampling.
    ! -------------------------------------------------------------------------

    CLASS(TimeSource_t), INTENT(INOUT) :: self      ! [-] Emitter instance
    TYPE(Config_t)     , INTENT(IN)    :: config    ! [-] Global configuration

    ! Local variables
    REAL(WP), ALLOCATABLE :: table(:,:)
    REAL(WP)              :: deviation
    INTEGER(I32)          :: i


    if (self % time_from_file) then
        CALL read_table(self % time_file, 1_I32, table)

        self % Nt = size(table, 1)
        if (self % Nt < 4_I32) then
            write(*, '(A,I0,A)') " ERROR TimeSource: the time file has ", self % Nt, &
                                 " samples (at least 4 are required)"
            error stop
        end if

        self % time = table(:, 1)
        self % dt   = (self % time(self % Nt) - self % time(1)) / real(self % Nt - 1_I32, WP)

        if (.not. self % dt > 0.0_WP) then
            write(*,*) " ERROR TimeSource: the times of the time file must increase"
            error stop
        end if

        deviation = maxval(abs(self % time(2:) - self % time(:self % Nt-1) - self % dt)) / self % dt
        if (deviation > 1.0e-6_WP) then
            write(*, '(A,ES9.2,A)') " ERROR TimeSource: the time file is not uniformly spaced (its time step " // &
                                    "varies by ", deviation, " of its mean value); the FFT needs a constant time step"
            error stop
        end if
    else
        self % Nt = config % time_source % Nt
        self % dt = (config % time_source % t_end - config % time_source % t_start) / real(self % Nt, WP)

        ALLOCATE(self % time(self % Nt))
        do i = 1, self % Nt
            self % time(i) = config % time_source % t_start + real(i - 1_I32, WP) * self % dt
        end do
    end if

    write(*, '(A,I0,A,ES12.5,A,ES12.5,A)') ' >> [TimeSource] Time axis: ', self % Nt, ' samples, time step dt = ', &
                                            self % dt, ' s (sampling frequency ', 1.0_WP / self % dt, ' Hz)'
    write(*, '(A,ES12.5,A,ES12.5,A,ES12.5,A)') ' >> [TimeSource] Record length = ', real(self % Nt, WP) * self % dt, &
                                                ' s, frequency resolution df = ', &
                                                1.0_WP / (real(self % Nt, WP) * self % dt), &
                                                ' Hz, Nyquist frequency = ', 0.5_WP / self % dt, ' Hz'

END SUBROUTINE build_time_axis


SUBROUTINE build_force_series(self)
    ! -------------------------------------------------------------------------
    ! Builds the force series Ft (Nt, Nsources, 3): the columns of a file (Fx,
    ! Fy, Fz of source 1, then those of source 2, ...) or the text expressions
    ! evaluated at the time axis. A value that is not finite (a logarithm of a
    ! negative number, a division by zero) stops the run, naming the force and
    ! the first time at which it happens.
    ! -------------------------------------------------------------------------

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite

    CLASS(TimeSource_t), INTENT(INOUT) :: self    ! [-] Emitter instance

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: message
    CHARACTER(len=2)              :: component(3) = ['Fx', 'Fy', 'Fz']
    REAL(WP)        , ALLOCATABLE :: table(:,:)
    TYPE(Expression_t)            :: expr
    INTEGER(I32)                  :: i_s, k, first_bad


    ALLOCATE(self % Ft(self % Nt, self % Nsources, 3))

    if (self % force_from_file) then
        CALL read_table(self % force_file, 3_I32 * self % Nsources, table)

        if (size(table, 1) /= self % Nt) then
            write(*, '(A,I0,A,I0,A)') " ERROR TimeSource: the force file has ", size(table, 1), &
                                      " rows and the time axis has ", self % Nt, " samples"
            error stop
        end if

        do i_s = 1, self % Nsources
            self % Ft(:, i_s, :) = table(:, 3*(i_s-1)+1:3*i_s)
        end do
    else
        do i_s = 1, self % Nsources
            do k = 1, 3
                CALL compile_expression(self % force_expr(i_s, k), expr, message)
                if (len(message) > 0) then
                    write(*, '(A,I0,A)') " ERROR TimeSource: force_expr(", i_s, ", " // component(k) // "): " // message
                    error stop
                end if

                CALL evaluate_expression(expr, self % time, self % Ft(:, i_s, k))
            end do
        end do
    end if

    do i_s = 1, self % Nsources
        do k = 1, 3
            if (.not. all(ieee_is_finite(self % Ft(:, i_s, k)))) then
                first_bad = minloc(merge(1, 0, ieee_is_finite(self % Ft(:, i_s, k))), dim=1)
                write(*, '(A,I0,A,ES12.5,A)') " ERROR TimeSource: the force " // component(k) // " of source ", i_s, &
                                              " is not finite at t = ", self % time(first_bad), " s"
                error stop
            end if
        end do
    end do

END SUBROUTINE build_force_series


SUBROUTINE take_to_frequency_domain(self)
    ! -------------------------------------------------------------------------
    ! Takes the force series to the frequency domain with compute_rfft and keeps
    ! the frequencies of the band fmin to fmax (and one out of every skipf
    ! bins). They are the Freqs and Force that the solver receives. The zero
    ! frequency is never kept: a constant force does not radiate.
    ! -------------------------------------------------------------------------

    CLASS(TimeSource_t), INTENT(INOUT) :: self    ! [-] Emitter instance

    ! Local variables
    COMPLEX(WP) , ALLOCATABLE :: spectrum(:,:,:)
    REAL(WP)    , ALLOCATABLE :: f(:)
    INTEGER(I32), ALLOCATABLE :: kept(:)
    LOGICAL     , ALLOCATABLE :: in_band(:)
    REAL(WP)                  :: f_low, f_high
    INTEGER(I32)              :: i


    CALL compute_rfft(self % Ft, self % Nt, self % dt, skipf=self % skipf, remove_zero=.true., &
                      array_out=spectrum, freqs=f)

    f_low  = 0.0_WP
    f_high = huge(1.0_WP)
    if (self % fmin /= DEFAULT_R) f_low  = self % fmin
    if (self % fmax /= DEFAULT_R) f_high = self % fmax

    if (f_high > 0.5_WP / self % dt .and. self % fmax /= DEFAULT_R) &
        write(*, '(A,ES10.3,A,ES10.3,A)') " WARNING TimeSource: fmax = ", self % fmax, &
                                          " Hz is above the Nyquist frequency ", 0.5_WP / self % dt, " Hz"

    in_band = f >= f_low .and. f <= f_high
    if (.not. any(in_band)) then
        write(*, '(A,ES10.3,A,ES10.3,A)') " ERROR TimeSource: no frequency of the spectrum (", f(1), " to ", &
                                          f(size(f)), " Hz) is inside fmin and fmax"
        error stop
    end if

    kept         = pack([(i, i = 1_I32, size(f))], in_band)
    self % Freqs = f(kept)
    self % Force = spectrum(kept, :, :)

    write(*, '(A,I0,A,ES10.3,A,ES10.3,A)') ' >> [TimeSource] Forces taken to the frequency domain (Hann window): ', &
                                           size(kept), ' frequencies from ', self % Freqs(1), ' to ', &
                                           self % Freqs(size(kept)), ' Hz'

END SUBROUTINE take_to_frequency_domain


SUBROUTINE save_params_TimeSource(self, save_path)
    ! -------------------------------------------------------------------------
    ! Saves the datasets of ManualSources_t and the force series: the number of
    ! samples, the time step, the time axis (Nt) and the forces (Nt, Nsources, 3).
    ! -------------------------------------------------------------------------

    CLASS(TimeSource_t), INTENT(IN) :: self         ! [-] Emitter instance
    CHARACTER(len=*)   , INTENT(IN) :: save_path    ! [-] Path to the HDF5 output file


    CALL self % save_params_base(save_path)

    CALL save_variable(save_path, "Nt"               , self % Nt)
    CALL save_variable(save_path, "Time_step"        , self % dt)
    CALL save_variable(save_path, "Source_time"      , self % time)
    CALL save_variable(save_path, "Source_Force_time", self % Ft)

END SUBROUTINE save_params_TimeSource


SUBROUTINE write_txt_TimeSource(self, unit)
    ! -------------------------------------------------------------------------
    ! Writes the sources block of the .dat summary: the title, orders and
    ! positions of ManualSources_t, the file of the positions if they were read
    ! from one, the time axis, the frequencies kept and the origin of the
    ! forces (their expressions, or the file), with a line saying that the Fx
    ! series of a monopole is its volume velocity Q [m^3/s] if there is one. A
    ! series is too long to be listed, so for every force it writes the sum of
    ! its values, of their squares and of the values weighted by the sample
    ! index: any change in the series changes them, and check_update sees it.
    ! -------------------------------------------------------------------------

    CLASS(TimeSource_t), INTENT(IN) :: self    ! [-] Emitter instance
    INTEGER(I32)       , INTENT(IN) :: unit    ! [-] Unit of the open .dat file

    ! Local variables
    CHARACTER(len=2) :: component(3) = ['Fx', 'Fy', 'Fz']
    REAL(WP)         :: weight(self % Nt)
    INTEGER(I32)     :: i, i_s, k


    CALL self % write_txt_base(unit)

    if (self % source_from_file) then
        write(unit, '(A)') 'Source file : ' // trim(self % source_file)
        write(unit, *)
    end if

    write(unit, '(A)') 'Time axis'
    write(unit, '(A,I0)')      'N samples   : ', self % Nt
    write(unit, '(A,ES23.15)') 'First time  : ', self % time(1)
    write(unit, '(A,ES23.15)') 'Time step   : ', self % dt
    write(unit, '(A,I0)')      'Skipf       : ', self % skipf
    write(unit, '(A,ES23.15)') 'Fmin        : ', self % fmin
    write(unit, '(A,ES23.15)') 'Fmax        : ', self % fmax
    if (self % time_from_file) write(unit, '(A)') 'Time file   : ' // trim(self % time_file)
    write(unit, *)

    write(unit, '(A)') 'Force [N] as a function of time (the Hann-windowed FFT of every series goes to the solver)'
    if (any(self % order == SOURCE_MONOPOLE)) &
        write(unit, '(A)') 'Monopoles   : volume velocity Q [m^3/s] in the Fx columns (Fy and Fz are 0)'
    if (self % force_from_file) then
        write(unit, '(A)') 'Force file  : ' // trim(self % force_file)
    else
        do i_s = 1, self % Nsources
            do k = 1, 3
                if (trim(self % force_expr(i_s, k)) /= '0') &
                    write(unit, '(I4,2X,A,2X,A)') i_s, component(k), trim(self % force_expr(i_s, k))
            end do
        end do
    end if

    write(unit, '(A)') '   s  comp               sum(F)             sum(F**2)               sum(i*F)'
    weight = [(real(i, WP), i = 1_I32, self % Nt)]
    do i_s = 1, self % Nsources
        do k = 1, 3
            write(unit, '(I4,2X,A,3ES23.15)') i_s, component(k), sum(self % Ft(:, i_s, k)), &
                                              sum(self % Ft(:, i_s, k)**2), sum(weight * self % Ft(:, i_s, k))
        end do
    end do
    write(unit, *)

END SUBROUTINE write_txt_TimeSource


FUNCTION get_enbw_TimeSource(self) RESULT(enbw)
    ! -------------------------------------------------------------------------
    ! Emitter_t % get_enbw of a time-domain source: the forces come from
    ! compute_rfft, so their window is the Hann window (HANN_ENBW).
    ! -------------------------------------------------------------------------

    CLASS(TimeSource_t), INTENT(IN) :: self    ! [-] Emitter instance
    REAL(WP)                        :: enbw    ! [-] Equivalent noise bandwidth in bins


    enbw = HANN_ENBW

END FUNCTION get_enbw_TimeSource

END MODULE TimeSource
