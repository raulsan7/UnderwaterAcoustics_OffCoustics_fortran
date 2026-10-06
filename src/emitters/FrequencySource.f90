! ===============================================================================
! MODULE: FrequencySource
! ===============================================================================
! Defines FrequencySource_t, an extension of Emitter_t made of a few point
! dipoles placed by hand, whose complex forces are given directly in the
! frequency domain (FrequencySource.nml). It is the simplest emitter: there is
! no OpenFAST run, no FFT window (get_enbw keeps -1, so no ENBW dataset is
! written) and no radiation-impedance correction, so the forces handed to the
! solver are the ones of the input file.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================

MODULE FrequencySource

USE Kinds, ONLY: I32, WP
USE Emitter, ONLY: Emitter_t
USE ConfigTypes, ONLY: Config_t
USE HDF5Utils, ONLY: save_variable
USE SourceTypes, ONLY: AcousticSource_t, SOURCE_DIPOLE

IMPLICIT NONE

PRIVATE
PUBLIC :: FrequencySource_t


TYPE, EXTENDS(Emitter_t) :: FrequencySource_t
    ! -------------------------------------------------------------------------
    ! Point dipoles with a complex force per frequency. The emitter keeps its
    ! own copy of the input (init) because it outlives to_source: the solver
    ! moves the arrays of the AcousticSource_t away, while save_params and
    ! write_txt still need them.
    ! -------------------------------------------------------------------------

    INTEGER(I32)              :: Nsources       ! [-] Number of sources
    INTEGER(I32), ALLOCATABLE :: order(:)       ! [-] Multipole order per source
    REAL(WP)    , ALLOCATABLE :: Freqs(:)       ! [Hz] Frequency vector (Nfreqs)
    REAL(WP)    , ALLOCATABLE :: xyz(:,:)       ! [m] Source coordinates (Nsources, 3)
    COMPLEX(WP) , ALLOCATABLE :: Force(:,:,:)   ! [N] Complex force on the fluid (Nfreqs, Nsources, 3)

    CONTAINS

    ! --- Public type-bound procedures --- !
    PROCEDURE :: init => init_FrequencySource                   ! Copies the source of the configuration
    PROCEDURE :: to_source => to_source_FrequencySource         ! Converts the emitter to an AcousticSource_t
    PROCEDURE :: save_params => save_params_FrequencySource     ! Saves the emitter parameters to HDF5
    PROCEDURE :: write_txt => write_txt_FrequencySource         ! Writes the emitter parameters to a text file

END TYPE FrequencySource_t


CONTAINS

SUBROUTINE init_FrequencySource(self, config)
    ! -------------------------------------------------------------------------
    ! Copies the frequency-domain sources of the configuration into the
    ! emitter.
    ! -------------------------------------------------------------------------

    CLASS(FrequencySource_t), INTENT(INOUT) :: self     ! [-] Emitter instance
    TYPE(Config_t)          , INTENT(IN)    :: config   ! [-] Global configuration


    self % Nsources = config % freq_source % Nsources
    self % Freqs    = config % freq_source % Freqs
    self % xyz      = config % freq_source % source_xyz
    self % Force    = config % freq_source % Force

    ! Every source radiates as a point dipole
    ALLOCATE(self % order(self % Nsources), source=SOURCE_DIPOLE)

    write(*, '(A,I0,A,I0,A)') ' >> [FrequencySource] Initialized with ', self % Nsources, ' sources at ', &
                              size(self % Freqs), ' frequencies'

END SUBROUTINE init_FrequencySource


SUBROUTINE to_source_FrequencySource(self, config, sources)
    ! -------------------------------------------------------------------------
    ! Translates the input sources into the AcousticSource_t consumed by the
    ! acoustic solvers. No radiation-impedance correction is applied: the
    ! forces are the ones of FrequencySource.nml. The arrays are copied, not
    ! moved, because save_params and write_txt need the emitter data after the
    ! solver has taken the sources. Every source is a point dipole; the
    ! reference centre and origin are the mean of the source positions (x, y)
    ! and the reference direction is 0 deg (+x).
    ! -------------------------------------------------------------------------

    CLASS(FrequencySource_t), INTENT(INOUT)            :: self      ! [-] Emitter instance
    TYPE(Config_t)          , INTENT(IN)               :: config    ! [-] Global configuration
    TYPE(AcousticSource_t)  , INTENT(OUT), ALLOCATABLE :: sources   ! [-] Acoustic point-dipole sources

    ! Local variables
    CHARACTER(len=16) :: ns_str, nf_str


    ALLOCATE(AcousticSource_t :: sources)

    sources % xyz   = self % xyz
    sources % Force = self % Force
    sources % Freqs = self % Freqs
    sources % order = self % order

    ! Reference geometry and default output lables read by the solver
    write(ns_str, '(I0)') self % Nsources
    write(nf_str, '(I0)') size(self % Freqs)

    sources % geometry % center = sum(self % xyz(:,1:2), dim=1) / real(self % Nsources, kind=WP)
    sources % geometry % origin = sources % geometry % center
    sources % geometry % direction = 0.0_WP
    sources % geometry % out_dir = './source_acoustic_data/'
    sources % geometry % name_tag = 'Source_Ns' // trim(ns_str) // '_Nf' // trim(nf_str)
    sources % geometry % name_suffix = ''

    write(*,*) ">> [AcousticSource] Frequency-domain forces handed to the solver as point dipoles"

END SUBROUTINE to_source_FrequencySource


SUBROUTINE save_params_FrequencySource(self, save_path)
    ! -------------------------------------------------------------------------
    ! Saves the emitter type, the number of sources and frequencies, and the
    ! sources themselves into the HDF5 output file: positions (Nsources, 3),
    ! multipole order (Nsources) and complex forces (Nfreqs, Nsources, 3).
    ! -------------------------------------------------------------------------

    CLASS(FrequencySource_t), INTENT(IN) :: self         ! [-] Emitter instance
    CHARACTER(len=*)        , INTENT(IN) :: save_path    ! [-] Path to the HDF5 output file


    CALL save_variable(save_path, "Emitter"     , "FrequencySource")
    CALL save_variable(save_path, "Nsources"    , self % Nsources)
    CALL save_variable(save_path, "NFreqs"      , size(self % Freqs))
    CALL save_variable(save_path, "Source_xyz"  , self % xyz)
    CALL save_variable(save_path, "Source_order", self % order)
    CALL save_variable(save_path, "Source_Force", self % Force)

END SUBROUTINE save_params_FrequencySource


SUBROUTINE write_txt_FrequencySource(self, unit)
    ! -------------------------------------------------------------------------
    ! Writes the sources block of the .dat summary: the positions and the
    ! complex force of every source at every frequency. check_update compares
    ! this block to decide whether an existing output holds the same case, so
    ! it lists every input value (nine significant digits), not only the counts.
    ! -------------------------------------------------------------------------

    CLASS(FrequencySource_t), INTENT(IN) :: self    ! [-] Emitter instance
    INTEGER(I32)            , INTENT(IN) :: unit    ! [-] Unit of the open .dat file

    ! Local variables
    INTEGER(I32) :: i_f, i_s, k
    COMPLEX(WP)  :: F_is(3)

    write(unit, '(A)') repeat('-', 70)
    write(unit, '(A)') 'FREQUENCY-DOMAIN SOURCES'
    write(unit, '(A)') repeat('-', 70)
    write(unit, '(A,I0)')     'N sources   : ', self % Nsources
    write(unit, '(A,I0)')     'N freqs     : ', size(self % Freqs)
    write(unit, '(A)')        'Multipole   : dipole'
    write(unit, *)

    write(unit, '(A)') '   #         x [m]         y [m]         z [m]'
    do i_s = 1, self % Nsources
        write(unit, '(I4,3F14.6)') i_s, self % xyz(i_s, :)
    end do
    write(unit, *)

    write(unit, '(A)') 'Force [N]: peak complex amplitude, time factor exp(-i*w*t)'
    write(unit, '(A)') '   f   s   Freq [Hz]          Re(Fx)          Im(Fx)          Re(Fy)' // &
                       '          Im(Fy)          Re(Fz)          Im(Fz)'
    do i_f = 1, size(self % Freqs)
        do i_s = 1, self % Nsources
            F_is = self % Force(i_f, i_s, :)
            write(unit, '(I4,I4,F12.4,6ES16.8)') i_f, i_s, self % Freqs(i_f), &
                                                 (real(F_is(k), WP), aimag(F_is(k)), k = 1, 3)
        end do
    end do
    write(unit, *)

END SUBROUTINE write_txt_FrequencySource

END MODULE FrequencySource