! ===============================================================================
! MODULE: StringUtils
! ===============================================================================
! Defines small string helpers: upper-case conversion (to_upper), the
! splitting of a line into its blank-separated tokens (split_tokens), and
! the detection of the data lines of an OpenFAST .out file (is_data_line).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE StringUtils

USE Kinds, ONLY: I32

IMPLICIT NONE

PRIVATE
PUBLIC :: to_upper, split_tokens, is_data_line


CONTAINS


PURE SUBROUTINE split_tokens(text_line, tokens)
    ! -------------------------------------------------------------------------
    ! Splits "text_line" into its tokens separated by blanks, TABs or carriage
    ! returns (each token truncated to 32 characters).
    ! -------------------------------------------------------------------------

    CHARACTER(len=*) , INTENT(IN)               :: text_line    ! [-] Line to split
    CHARACTER(len=32), INTENT(OUT), ALLOCATABLE :: tokens(:)    ! [-] Tokens found, in order

    ! Local variables
    INTEGER(I32) :: i, n, istart, pass
    LOGICAL      :: in_token, is_sep


    do pass = 1, 2
        n        = 0_I32
        in_token = .false.
        istart   = 0_I32
        do i = 1, len(text_line) + 1
            if (i > len(text_line)) then
                is_sep = .true.
            else
                is_sep = text_line(i:i) == ' ' .or. text_line(i:i) == achar(9) .or. text_line(i:i) == achar(13)
            end if
            if (.not. is_sep .and. .not. in_token) then
                in_token = .true.
                istart   = i
            else if (is_sep .and. in_token) then
                in_token = .false.
                n        = n + 1_I32
                if (pass == 2_I32) tokens(n) = text_line(istart:i-1)
            end if
        end do
        if (pass == 1_I32) ALLOCATE(tokens(n))
    end do

END SUBROUTINE split_tokens


PURE FUNCTION is_data_line(text_line) RESULT(is_data)
    ! -------------------------------------------------------------------------
    ! True when the first non-blank character among the first 256 characters
    ! of "text_line" exists and is not '#': blank lines and comment lines such
    ! as "#Restarting here" are not data rows of an OpenFAST .out file.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: text_line    ! [-] Line of an OpenFAST .out file
    LOGICAL                      :: is_data      ! [-] Whether the line holds a data row

    ! Local variables
    CHARACTER(len=256) :: head
    INTEGER(I32)       :: i


    head = text_line(1:min(len(text_line), 256))
    do i = 1, len_trim(head)
        if (head(i:i) == achar(9) .or. head(i:i) == achar(13)) head(i:i) = ' '
    end do
    head    = adjustl(head)
    is_data = len_trim(head) > 0 .and. head(1:1) /= '#'

END FUNCTION is_data_line


PURE FUNCTION to_upper(str) RESULT(res)
    ! -------------------------------------------------------------------------
    ! Returns "str" with every lowercase ASCII letter converted to uppercase.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*)       , INTENT(IN) :: str    ! [-] Input string
    CHARACTER(len=len(str))             :: res    ! [-] Uppercase copy of str

    ! Local variables
    INTEGER(I32) :: i, c


    res = str
    do i = 1, len(str)
        c = iachar(str(i:i))
        if (c >= 97_I32 .and. c <= 122_I32) res(i:i) = achar(c - 32_I32)
    end do

END FUNCTION to_upper

END MODULE StringUtils
