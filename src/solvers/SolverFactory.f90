! ===============================================================================
! MODULE: SolverFactory
! ===============================================================================
! Defines new_solver, the factory that allocates the acoustic solver of the
! concrete type named by the method of Solver.nml (ImageMethod,
! AnalyticalNormalModes).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE SolverFactory

USE AcousticSolver, ONLY: AcousticSolver_t

IMPLICIT NONE

PRIVATE
PUBLIC :: new_solver


CONTAINS


SUBROUTINE new_solver(method, solver)
    ! -------------------------------------------------------------------------
    ! Allocates an acoustic solver of the concrete type named by method
    ! (config % solver % method). The parser already rejects unknown methods,
    ! so the CASE DEFAULT never fires.
    ! -------------------------------------------------------------------------

    USE ImageMethod          , ONLY: ImageMethod_t
    USE AnalyticalNormalModes, ONLY: AnalyticalNormalModes_t

    CHARACTER(len=*)       , INTENT(IN)               :: method    ! [-] Acoustic method tag
    CLASS(AcousticSolver_t), INTENT(OUT), ALLOCATABLE :: solver    ! [-] New solver instance


    SELECT CASE (trim(method))
    CASE ("ImageMethod")
        ALLOCATE(ImageMethod_t :: solver)
    CASE ("AnalyticalNormalModes")
        ALLOCATE(AnalyticalNormalModes_t :: solver)
    CASE DEFAULT
        error stop "new_solver(): Unknown method " // trim(method)
    END SELECT

END SUBROUTINE new_solver

END MODULE SolverFactory
