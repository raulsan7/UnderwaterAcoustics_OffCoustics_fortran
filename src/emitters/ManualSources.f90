! ===============================================================================
! MODULE: ManualSources
! ===============================================================================
! Defines ManualSources_t, the parent of the emitters made of point sources
! placed by hand (FrequencySource, TimeSource), each one a monopole or a
! dipole (its multipole order). Whatever way their input is given, they reach
! the solver as complex amplitudes per frequency (the force of a dipole, the
! volume velocity Q of a monopole in the first component; see SourceTypes),
! so everything from there on is shared: the sources handed to the solver
! (to_source), the datasets of the HDF5 output (save_params), the block of
! the .dat summary with the orders and positions (write_txt), and the check
! that a monopole has nothing in its second and third components
! (check_monopole_components). A child builds its data in init, sets the
! three labels below, and adds what is its own; the parent version of
! save_params and write_txt stays reachable from a child override through
! save_params_base and write_txt_base.
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
USE SourceTypes, ONLY: AcousticSource_t, SOURCE_MONOPOLE, SOURCE_DIPOLE

IMPLICIT NONE

PRIVATE
PUBLIC :: ManualSources_t


TYPE, ABSTRACT, EXTENDS(Emitter_t) :: ManualSources_t
    ! -------------------------------------------------------------------------
    ! Point monopoles and dipoles with a complex amplitude per frequency. The
    ! emitter keeps its own copy of the data because it outlives to_source:
    ! the solver moves the arrays of the AcousticSource_t away, while
    ! save_params and write_txt still need them.
    ! -------------------------------------------------------------------------

    ! --- Labels, set by every child --- !
    CHARACTER(len=32) :: kind_name      ! [-] Name of the emitter, saved as the Emitter dataset
    CHARACTER(len=40) :: title          ! [-] Title of its block of the .dat summary
    CHARACTER(len=16) :: name_suffix    ! [-] End of the default output name (e.g. '_Nf3')

    ! --- Sources, as the solver receives them --- !
    INTEGER(I32)              :: Nsources       ! [-] Number of sources
    INTEGER(I32), ALLOCATABLE :: order(:)       ! [-] Multipole order per source (SOURCE_MONOPOLE, SOURCE_DIPOLE)
    REAL(WP)    , ALLOCATABLE :: Freqs(:)       ! [Hz] Frequency vector (Nfreqs)
    REAL(WP)    , ALLOCATABLE :: xyz(:,:)       ! [m] Source coordinates (Nsources, 3)
    COMPLEX(WP) , ALLOCATABLE :: Force(:,:,:)   ! [N] Force on the fluid; monopole: Q [m^3/s] in 1 (Nfreqs, Nsources, 3)

    CONTAINS

    ! --- Shared procedures --- !
    PROCEDURE :: to_source   => to_source_ManualSources      ! Converts the emitter to an AcousticSource_t
    PROCEDURE :: save_params => save_params_ManualSources    ! Saves the sources to HDF5
    PROCEDURE :: write_txt   => write_txt_ManualSources      ! Writes the title, orders and positions to the .dat

    ! --- Shared check, called by every child in init: a monopole has nothing in components 2 and 3 --- !
    PROCEDURE :: check_monopole_components => check_monopole_components_ManualSources

    ! --- Aliases to call the parent version from a child override --- !
    PROCEDURE :: save_params_base => save_params_ManualSources    ! Parent save_params
    PROCEDURE :: write_txt_base   => write_txt_ManualSources      ! Parent write_txt

    ! --- Deferred (abstract) procedures --- !
    PROCEDURE(init_i), DEFERRED :: init                           ! Builds the sources and sets the labels

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
    ! solvers, with the multipole order of every source. No
    ! radiation-impedance correction is applied: the forces and volume
    ! velocities are the ones of the input. The arrays are copied, not moved,
    ! because save_params and write_txt need the emitter data after the solver
    ! has taken the sources. The reference centre and origin are the mean of
    ! the source positions (x, y) and the reference direction is 0 deg (+x).
    ! The message names the orders handed over (it is the one of before when
    ! every source is a dipole).
    ! -------------------------------------------------------------------------

    CLASS(ManualSources_t), INTENT(INOUT)            :: self      ! [-] Emitter instance
    TYPE(Config_t)        , INTENT(IN)               :: config    ! [-] Global configuration
    TYPE(AcousticSource_t), INTENT(OUT), ALLOCATABLE :: sources   ! [-] Acoustic point sources

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

    if (all(self % order == SOURCE_DIPOLE)) then
        write(*,*) ">> [AcousticSource] Frequency-domain forces handed to the solver as point dipoles"
    else if (all(self % order == SOURCE_MONOPOLE)) then
        write(*,*) ">> [AcousticSource] Frequency-domain forces handed to the solver as point monopoles"
    else
        write(*,*) ">> [AcousticSource] Frequency-domain forces handed to the solver as point monopoles and dipoles"
    end if

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
    ! the counts, the multipole orders and the positions. The line Multipole
    ! says 'dipole' or 'monopole' when every source has that order (so a set
    ! of dipoles keeps the block it always had) and 'mixed' otherwise; only
    ! then the table of positions gets the column order. A child adds the rest
    ! of the input after it. check_update compares this block to decide
    ! whether an existing output holds the same case, so the block must list
    ! every input value, the orders included.
    ! -------------------------------------------------------------------------

    CLASS(ManualSources_t), INTENT(IN) :: self    ! [-] Emitter instance
    INTEGER(I32)          , INTENT(IN) :: unit    ! [-] Unit of the open .dat file

    ! Local variables
    CHARACTER(len=8) :: multipole
    INTEGER(I32)     :: i_s


    if (all(self % order == SOURCE_DIPOLE)) then
        multipole = 'dipole'
    else if (all(self % order == SOURCE_MONOPOLE)) then
        multipole = 'monopole'
    else
        multipole = 'mixed'
    end if

    write(unit, '(A)') repeat('-', 70)
    write(unit, '(A)') trim(self % title)
    write(unit, '(A)') repeat('-', 70)
    write(unit, '(A,I0)') 'N sources   : ', self % Nsources
    write(unit, '(A,I0)') 'N freqs     : ', size(self % Freqs)
    write(unit, '(A)')    'Multipole   : ' // trim(multipole)
    write(unit, *)

    if (trim(multipole) == 'mixed') then
        write(unit, '(A)') '   #         x [m]         y [m]         z [m]  order [-]'
        do i_s = 1, self % Nsources
            write(unit, '(I4,3F14.6,I7)') i_s, self % xyz(i_s, :), self % order(i_s)
        end do
    else
        write(unit, '(A)') '   #         x [m]         y [m]         z [m]'
        do i_s = 1, self % Nsources
            write(unit, '(I4,3F14.6)') i_s, self % xyz(i_s, :)
        end do
    end if
    write(unit, *)

END SUBROUTINE write_txt_ManualSources


SUBROUTINE check_monopole_components_ManualSources(self, nonzero, what)
    ! -------------------------------------------------------------------------
    ! Stops the run if a monopole has something in the second or third
    ! component of its input: its volume velocity Q goes in the first one
    ! (SourceTypes) and the other two must be 0. Each child calls it in init
    ! on what it reads: FrequencySource on the forces, TimeSource on the
    ! series before the FFT (which removes the mean and would hide a constant
    ! component).
    ! -------------------------------------------------------------------------

    CLASS(ManualSources_t), INTENT(IN) :: self            ! [-] Emitter instance
    LOGICAL               , INTENT(IN) :: nonzero(:,:)    ! [-] Whether each component is not 0 (Nsources, 3)
    CHARACTER(len=*)      , INTENT(IN) :: what            ! [-] Input checked, shown in the error message

    ! Local variables
    CHARACTER(len=2), PARAMETER :: COMPONENT(3) = ['Fx', 'Fy', 'Fz']
    INTEGER(I32)                :: i_s, k


    do i_s = 1, self % Nsources
        if (self % order(i_s) /= SOURCE_MONOPOLE) cycle
        do k = 2, 3
            if (nonzero(i_s, k)) then
                write(*, '(A,I0,A)') " ERROR " // trim(self % kind_name) // ": source ", i_s, &
                                     " is a monopole (its volume velocity Q goes in the first component) " // &
                                     "and its " // what // " " // COMPONENT(k) // " is not 0"
                error stop
            end if
        end do
    end do

END SUBROUTINE check_monopole_components_ManualSources

END MODULE ManualSources
