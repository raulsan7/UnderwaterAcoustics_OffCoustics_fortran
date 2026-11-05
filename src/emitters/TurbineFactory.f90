! ===============================================================================
! MODULE: TurbineFactory
! ===============================================================================
! Defines new_turbine, the factory that allocates a wind turbine of the
! concrete type named by its case_type (DTU10MWMonopile, DTU10MWFloating,
! SAITEC2MWFloating, MYSE3MWMonopile). Shared by the single-turbine and the
! wind-farm paths.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE TurbineFactory

USE WindTurbine, ONLY: WindTurbine_t

IMPLICIT NONE

PRIVATE
PUBLIC :: new_turbine


CONTAINS


SUBROUTINE new_turbine(case_type, wt)
    ! -------------------------------------------------------------------------
    ! Allocates a turbine of the concrete type named by case_type. Shared by
    ! the single-turbine and the wind-farm paths.
    ! -------------------------------------------------------------------------

    USE DTU10MWMonopile  , ONLY: DTU10MWMonopile_t
    USE DTU10MWFloating  , ONLY: DTU10MWFloating_t
    USE MYSE3MWMonopile  , ONLY: MYSE3MWMonopile_t
    USE SAITEC2MWFloating, ONLY: SAITEC2MWFloating_t

    CHARACTER(len=*)    , INTENT(IN)               :: case_type    ! [-] Turbine model tag
    CLASS(WindTurbine_t), INTENT(OUT), ALLOCATABLE :: wt           ! [-] New turbine instance


    SELECT CASE (trim(case_type))
    CASE ('DTU10MWMonopile')
        ALLOCATE(DTU10MWMonopile_t :: wt)
    CASE ('DTU10MWFloating')
        ALLOCATE(DTU10MWFloating_t :: wt)
    CASE ('SAITEC2MWFloating')
        ALLOCATE(SAITEC2MWFloating_t :: wt)
    CASE ('MYSE3MWMonopile')
        ALLOCATE(MYSE3MWMonopile_t :: wt)
    CASE DEFAULT
        write(*,*) 'Unsupported turbine type: ' // trim(case_type)
        error stop
    END SELECT

END SUBROUTINE new_turbine

END MODULE TurbineFactory
