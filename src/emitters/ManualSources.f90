! ===============================================================================
! MODULE: ManualSources
! ===============================================================================
! Defines ManualSources_t, the parent of the emitters made of a few
! point dipoles placed by hand (FrequencySource, TimeSource). Whatever way
! their forces are given, they reach the solver as complex forces per
! frequency, so everything from there on is shared: the sources handed to the
! solver (to_source), the datasets of the HDF5 output (save_params) and the
! block of the .dat summary with the positions (write_txt). A child builds
! its data in init, sets the three labels below, and adds what is its own;
! the parent version of save_params and write_txt stays reachable from a
! child override through save_params_base and write_txt_base.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE ManualSources

USE Kinds, ONLY: I32, WP
USE Emitter, ONLY: Emitter_t
USE ConfigTypes, ONLY: Config_t
USE HDF5Utils, ONLY: save_variable
USE SourceTypes, ONLY: AcousticSource_t

IMPLICIT NONE

PRIVATE
PUBLIC :: ManualSources_t


TYPE, ABSTRACT, EXTENDS(Emitter_t) :: ManualSources_t
    ! -------------------------------------------------------------------------
    ! Point dipoles with a complex force per frequency. The emitter keeps its
    ! own copy of the data because it outlives to_source: the solver moves the
    ! arrays of the AcousticSource_t away, while save_params and write_txt
    ! still need them.
    ! -------------------------------------------------------------------------

    ! --- Labels, set by every child --- !
    CHARACTER(len=32) :: kind_name      ! [-] Name of the emitter, saved as the Emitter dataset
    CHARACTER(len=40) :: title          ! [-] Title of its block of the .dat summary
    CHARACTER(len=16) :: name_suffix    ! [-] End of the default output name (e.g. '_Nf3')

    ! --- Sources, as the solver receives them --- !
    INTEGER(I32)              :: Nsources       ! [-] Number of sources
    INTEGER(I32), ALLOCATABLE :: order(:)       ! [-] Multipole order per source
    REAL(WP)    , ALLOCATABLE :: Freqs(:)       ! [Hz] Frequency vector (Nfreqs)
    REAL(WP)    , ALLOCATABLE :: xyz(:,:)       ! [m] Source coordinates (Nsources, 3)
    COMPLEX(WP) , ALLOCATABLE :: Force(:,:,:)   ! [N] Complex force on the fluid (Nfreqs, Nsources, 3)

    CONTAINS

    ! --- Shared procedures --- !
    PROCEDURE :: to_source   => to_source_ManualSources      ! Converts the emitter to an AcousticSource_t
    PROCEDURE :: save_params => save_params_ManualSources    ! Saves the sources to HDF5
    PROCEDURE :: write_txt   => write_txt_ManualSources      ! Writes the title and the positions to the .dat

    ! --- Aliases to call the parent version from a child override --- !
    PROCEDURE :: save_params_base => save_params_ManualSources    ! Parent save_params
    PROCEDURE :: write_txt_base   => write_txt_ManualSources      ! Parent write_txt

    ! --- Deferred (abstract) procedures --- !
    PROCEDURE(init_i), DEFERRED :: init                           ! Builds the sources from the configuration, sets the labels

END TYPE ManualSources_t


! ===============================================================================
! ABSTRACT INTERFACES FOR DEFERRED PROCEDURES
! ===============================================================================
ABSTRACT INTERFACE

    SUBROUTINE init_i(self, config)
        ! -------------------------------------------------------------------------
        ! Builds the emitter from the global configuration: fills Nsources, order,
        ! Freqs, xyz and Force and sets kind_name, title and name_suffix.
        ! -------------------------------------------------------------------------

        IMPORT :: ManualSources_t, Config_t

        CLASS(ManualSources_t), INTENT(INOUT) :: self      ! [-] Emitter instance
        TYPE(Config_t)        , INTENT(IN)    :: config    ! [-] Global configuration

    END SUBROUTINE init_i

END INTERFACE


CONTAINS


SUBROUTINE to_source_ManualSources(self, config, sources)
    ! -------------------------------------------------------------------------
    ! Translates the sources into the AcousticSource_t consumed by the acoustic
    ! solvers. No radiation-impedance correction is applied: the forces are the
    ! ones of the input. The arrays are copied, not moved, because save_params
    ! and write_txt need the emitter data after the solver has taken the
    ! sources. The reference centre and origin are the mean of the source
    ! positions (x, y) and the reference direction is 0 deg (+x).
    ! -------------------------------------------------------------------------

    CLASS(ManualSources_t), INTENT(INOUT)            :: self      ! [-] Emitter instance
    TYPE(Config_t)        , INTENT(IN)               :: config    ! [-] Global configuration
    TYPE(AcousticSource_t), INTENT(OUT), ALLOCATABLE :: sources   ! [-] Acoustic point-dipole sources

    ! Local variables
    CHARACTER(len=16) :: ns_str


    ALLOCATE(AcousticSource_t :: sources)

    sources % xyz   = self % xyz
    sources % Force = self % Force
    sources % Freqs = self % Freqs
    sources % order = self % order

    ! Reference geometry and default output labels read by the solver
    write(ns_str, '(I0)') self % Nsources

    sources % geometry % center      = sum(self % xyz(:, 1:2), dim=1) / real(self % Nsources, WP)
    sources % geometry % origin      = sources % geometry % center
    sources % geometry % direction   = 0.0_WP
    sources % geometry % out_dir     = './source_acoustic_data/'
    sources % geometry % name_tag    = 'Source_Ns' // trim(ns_str) // trim(self % name_suffix)
    sources % geometry % name_suffix = ''

    write(*,*) ">> [AcousticSource] Frequency-domain forces handed to the solver as point dipoles"

END SUBROUTINE to_source_ManualSources


SUBROUTINE save_params_ManualSources(self, save_path)
    ! -------------------------------------------------------------------------
    ! Saves the emitter type, the number of sources and frequencies, and the
    ! sources themselves into the HDF5 output file: positions (Nsources, 3),
    ! multipole order (Nsources) and complex forces (Nfreqs, Nsources, 3).
    ! -------------------------------------------------------------------------

    CLASS(ManualSources_t), INTENT(IN) :: self         ! [-] Emitter instance
    CHARACTER(len=*)      , INTENT(IN) :: save_path    ! [-] Path to the HDF5 output file


    CALL save_variable(save_path, "Emitter"     , trim(self % kind_name))
    CALL save_variable(save_path, "Nsources"    , self % Nsources)
    CALL save_variable(save_path, "NFreqs"      , size(self % Freqs))
    CALL save_variable(save_path, "Source_xyz"  , self % xyz)
    CALL save_variable(save_path, "Source_order", self % order)
    CALL save_variable(save_path, "Source_Force", self % Force)

END SUBROUTINE save_params_ManualSources


SUBROUTINE write_txt_ManualSources(self, unit)
    ! -------------------------------------------------------------------------
    ! Writes the first part of the sources block of the .dat summary: the title,
    ! the counts and the positions. A child adds the rest of the input after it.
    ! check_update compares this block to decide whether an existing output
    ! holds the same case, so the block must list every input value.
    ! -------------------------------------------------------------------------

    CLASS(ManualSources_t), INTENT(IN) :: self    ! [-] Emitter instance
    INTEGER(I32)          , INTENT(IN) :: unit    ! [-] Unit of the open .dat file

    ! Local variables
    INTEGER(I32) :: i_s


    write(unit, '(A)') repeat('-', 70)
    write(unit, '(A)') trim(self % title)
    write(unit, '(A)') repeat('-', 70)
    write(unit, '(A,I0)') 'N sources   : ', self % Nsources
    write(unit, '(A,I0)') 'N freqs     : ', size(self % Freqs)
    write(unit, '(A)')    'Multipole   : dipole'
    write(unit, *)

    write(unit, '(A)') '   #         x [m]         y [m]         z [m]'
    do i_s = 1, self % Nsources
        write(unit, '(I4,3F14.6)') i_s, self % xyz(i_s, :)
    end do
    write(unit, *)

END SUBROUTINE write_txt_ManualSources

END MODULE ManualSources
