! ===============================================================================
! MODULE: Emitter
! ===============================================================================
! Defines Emitter_t, the abstract parent type of everything that produces
! acoustic sources (a wind turbine, a wind farm, and in the future manually
! defined sources). The acoustic solvers only see this interface: the
! AcousticSource_t it fills in (to_source), the blocks it writes in the HDF5
! (save_params) and .dat (write_txt) outputs, and the equivalent noise
! bandwidth of the window applied to its forces (get_enbw).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE Emitter

USE Kinds, ONLY: I32, WP
USE ConfigTypes, ONLY: Config_t
USE SourceTypes, ONLY: AcousticSource_t

IMPLICIT NONE

PRIVATE
PUBLIC :: Emitter_t


TYPE, ABSTRACT :: Emitter_t
    ! -------------------------------------------------------------------------
    ! Parent of everything that produces acoustic sources (a wind turbine, a
    ! wind farm, and in the future manually defined sources). The solvers only
    ! see this interface: the sources it fills in and the blocks it writes in
    ! the HDF5 and .dat outputs.
    ! -------------------------------------------------------------------------

    CONTAINS

    ! --- Shared procedures --- !
    PROCEDURE :: get_enbw                                   ! ENBW of the window applied to the forces (-1: no window)

    ! --- Deferred (abstract) procedures --- !
    PROCEDURE(to_source_i)  , DEFERRED :: to_source         ! Fills the AcousticSource_t handed to the solver
    PROCEDURE(save_params_i), DEFERRED :: save_params       ! Emitter parameters -> HDF5 output
    PROCEDURE(write_txt_i)  , DEFERRED :: write_txt         ! Emitter block of the .dat summary

END TYPE Emitter_t


! ===============================================================================
! ABSTRACT INTERFACES FOR DEFERRED PROCEDURES
! ===============================================================================
ABSTRACT INTERFACE

    SUBROUTINE to_source_i(self, config, sources)
        ! -------------------------------------------------------------------------
        ! Translates the emitter into the AcousticSource_t consumed by the
        ! acoustic solvers: source coordinates, complex forces, frequency vector
        ! and reference geometry and labels of the emitter.
        ! -------------------------------------------------------------------------

        IMPORT :: Emitter_t, Config_t, AcousticSource_t

        CLASS(Emitter_t)      , INTENT(INOUT)            :: self       ! [-] Emitter instance
        TYPE(Config_t)        , INTENT(IN)               :: config     ! [-] Global configuration
        TYPE(AcousticSource_t), INTENT(OUT), ALLOCATABLE :: sources    ! [-] Acoustic point sources

    END SUBROUTINE to_source_i

    SUBROUTINE save_params_i(self, save_path)
        ! -------------------------------------------------------------------------
        ! Saves the emitter parameters into the HDF5 output file.
        ! -------------------------------------------------------------------------

        IMPORT :: Emitter_t

        CLASS(Emitter_t), INTENT(IN) :: self         ! [-] Emitter instance
        CHARACTER(len=*), INTENT(IN) :: save_path    ! [-] Path to the HDF5 output file

    END SUBROUTINE save_params_i

    SUBROUTINE write_txt_i(self, unit)
        ! -------------------------------------------------------------------------
        ! Writes the emitter block of the .dat summary file.
        ! -------------------------------------------------------------------------

        IMPORT :: Emitter_t, I32

        CLASS(Emitter_t), INTENT(IN) :: self    ! [-] Emitter instance
        INTEGER(I32)    , INTENT(IN) :: unit    ! [-] Unit of the open .dat file

    END SUBROUTINE write_txt_i

END INTERFACE


CONTAINS


FUNCTION get_enbw(self) RESULT(enbw)
    ! -------------------------------------------------------------------------
    ! Equivalent noise bandwidth (ENBW), in bins, of the window applied to the
    ! emitter forces. The default -1 means "forces not obtained from a windowed
    ! FFT". Emitters whose forces come from compute_rfft override it and
    ! return HANN_ENBW; the solvers write the ENBW dataset only when it is > 0.
    ! -------------------------------------------------------------------------

    CLASS(Emitter_t), INTENT(IN) :: self    ! [-] Emitter instance
    REAL(WP)                     :: enbw    ! [-] Equivalent noise bandwidth in bins (-1: no window)


    enbw = -1.0_WP

END FUNCTION get_enbw

END MODULE Emitter
