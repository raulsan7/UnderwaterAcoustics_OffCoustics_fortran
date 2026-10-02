! ===============================================================================
! MODULE: Sorting
! ===============================================================================
! Defines the ascending sort of a REAL(WP) array (sort_real_c) through a
! binding to the C standard library's qsort(3) and its comparator callback
! (compare_real).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE Sorting

USE Kinds, ONLY: WP
USE, INTRINSIC :: iso_c_binding, ONLY: C_INT, C_PTR, C_DOUBLE, C_SIZE_T, C_FUNPTR, c_f_pointer, c_loc, c_funloc, &
                                       c_sizeof

IMPLICIT NONE

PRIVATE
PUBLIC :: sort_real_c


! ===============================================================================
! C INTERFACES
! ===============================================================================
INTERFACE

    SUBROUTINE qsort(base, nmemb, size, compar) BIND(C, name='qsort')
        ! -------------------------------------------------------------------------
        ! Binds to the C standard library's qsort(3), used by sort_real_c to sort
        ! a REAL(WP) array in ascending order via the compare_real callback.
        ! -------------------------------------------------------------------------

        IMPORT :: C_PTR, C_SIZE_T, C_FUNPTR

        TYPE(C_PTR)      , VALUE :: base      ! [-] Pointer to the first element of the array
        INTEGER(C_SIZE_T), VALUE :: nmemb     ! [-] Number of elements in the array
        INTEGER(C_SIZE_T), VALUE :: size      ! [B] Size of each element
        TYPE(C_FUNPTR)   , VALUE :: compar    ! [-] Comparison callback

    END SUBROUTINE qsort

END INTERFACE


CONTAINS


FUNCTION compare_real(a, b) BIND(C) RESULT(res)
    ! -------------------------------------------------------------------------
    ! Comparison callback for qsort: ascending order for REAL(C_DOUBLE) values.
    ! -------------------------------------------------------------------------

    TYPE(C_PTR)   , VALUE :: a      ! [-] Pointer to the first value
    TYPE(C_PTR)   , VALUE :: b      ! [-] Pointer to the second value
    INTEGER(C_INT)        :: res    ! [-] -1, 0 or 1 if a is smaller than, equal to or greater than b

    ! Local variables
    REAL(C_DOUBLE), POINTER :: fa, fb


    CALL c_f_pointer(a, fa)
    CALL c_f_pointer(b, fb)

    if (fa < fb) then
        res = -1_C_INT
    else if (fa > fb) then
        res = 1_C_INT
    else
        res = 0_C_INT
    end if

END FUNCTION compare_real


SUBROUTINE sort_real_c(arr)
    ! -------------------------------------------------------------------------
    ! Sorts a REAL(WP) array in ascending order using the C standard library's
    ! qsort, via the compare_real comparator callback.
    ! -------------------------------------------------------------------------

    REAL(WP), INTENT(INOUT), TARGET :: arr(:)    ! [-] Array to sort in place

    ! Local variables
    REAL(WP)         , POINTER :: parr(:)
    INTEGER(C_SIZE_T)          :: nmemb, size_elem
    TYPE(C_PTR)                :: ptr


    parr      => arr
    nmemb     = size(parr, kind=C_SIZE_T)
    size_elem = c_sizeof(parr(1))
    ptr       = c_loc(parr(1))

    CALL qsort(ptr, nmemb, size_elem, c_funloc(compare_real))

END SUBROUTINE sort_real_c

END MODULE Sorting
