! ===============================================================================
! MODULE: FrequencySource
! ===============================================================================
! Defines FrequencySource_t, an extension of ManualSources_t whose forces are
! given directly in the frequency domain (FrequencySource.nml): a complex
! force per frequency, source and component. There is no FFT window (get_enbw
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

USE Kinds, ONLY: I32, WP
USE ConfigTypes, ONLY: Config_t
USE SourceTypes, ONLY: SOURCE_DIPOLE
USE ManualSources, ONLY: ManualSources_t

IMPLICIT NONE

PRIVATE
PUBLIC :: FrequencySource_t


TYPE, EXTENDS(ManualSources_t) :: FrequencySource_t
    ! -------------------------------------------------------------------------
    ! Point dipoles with a complex force per frequency, as read from the
    ! namelist. Everything is in ManualSources_t; this type only fills it and
    ! lists the forces in the .dat summary.
    ! -------------------------------------------------------------------------

    CONTAINS

    ! --- Public type-bound procedures --- !
    PROCEDURE :: init      => init_FrequencySource         ! Copies the sources of the configuration
    PROCEDURE :: write_txt => write_txt_FrequencySource    ! Adds the force table to the .dat

END TYPE FrequencySource_t


CONTAINS


SUBROUTINE init_FrequencySource(self, config)
    ! -------------------------------------------------------------------------
    ! Copies the frequency-domain sources of the configuration into the
    ! emitter.
    ! -------------------------------------------------------------------------

    CLASS(FrequencySource_t), INTENT(INOUT) :: self      ! [-] Emitter instance
    TYPE(Config_t)          , INTENT(IN)    :: config    ! [-] Global configuration

    ! Local variables
    CHARACTER(len=16) :: nf_str


    self % kind_name = 'FrequencySource'
    self % title     = 'FREQUENCY-DOMAIN SOURCES'

    self % Nsources = config % freq_source % Nsources
    self % Freqs    = config % freq_source % Freqs
    self % xyz      = config % freq_source % source_xyz
    self % Force    = config % freq_source % Force

    ! Every source radiates as a point dipole
    ALLOCATE(self % order(self % Nsources), source=SOURCE_DIPOLE)

    write(nf_str, '(I0)') size(self % Freqs)
    self % name_suffix = '_Nf' // trim(nf_str)

    write(*, '(A,I0,A,I0,A)') ' >> [FrequencySource] Initialized with ', self % Nsources, ' sources at ', &
                              size(self % Freqs), ' frequencies'

END SUBROUTINE init_FrequencySource


SUBROUTINE write_txt_FrequencySource(self, unit)
    ! -------------------------------------------------------------------------
    ! Writes the sources block of the .dat summary: the title and positions of
    ! ManualSources_t followed by the complex force of every source at every
    ! frequency, with nine significant digits, so that a change in any force is
    ! seen by check_update.
    ! -------------------------------------------------------------------------

    CLASS(FrequencySource_t), INTENT(IN) :: self    ! [-] Emitter instance
    INTEGER(I32)            , INTENT(IN) :: unit    ! [-] Unit of the open .dat file

    ! Local variables
    INTEGER(I32) :: i_f, i_s, k
    COMPLEX(WP)  :: F_is(3)


    CALL self % write_txt_base(unit)

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
