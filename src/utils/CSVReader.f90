! ===============================================================================
! MODULE: CSVReader
! ===============================================================================
! Defines the reader of comma-separated lookup tables: the linear
! interpolation of one column of a CSV table at a given wind speed
! (read_curve), and the parser of one of its rows (parse_csv_row).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE CSVReader

USE Kinds, ONLY: I32, WP
USE FileUtils, ONLY: read_full_line
USE StringUtils, ONLY: split_tokens

IMPLICIT NONE

PRIVATE
PUBLIC :: read_curve, read_table


CONTAINS


FUNCTION read_curve(filename, ws, col_idx) RESULT(val)
    ! -------------------------------------------------------------------------
    ! Reads a CSV lookup table (windspeed vs. one or more dependent columns)
    ! and linearly interpolates (with linear extrapolation outside the table's
    ! range) the requested column at the given wind speed. A line is a data
    ! row when its first col_idx comma-separated fields are all numeric; any
    ! other line (headers such as "rpm," and "X,Y" or "ws,rpm", blank lines,
    ! short rows) is skipped by both passes, so tables with one or two header
    ! lines are both read. A missing file, an invalid column or fewer than 2
    ! data rows stops the program.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)           :: filename    ! [-] Path to the CSV curve file
    REAL(WP)        , INTENT(IN)           :: ws          ! [m/s] Wind speed
    INTEGER(I32)    , INTENT(IN), OPTIONAL :: col_idx     ! [-] Column to evaluate
    REAL(WP)                               :: val         ! [-] Interpolated scalar

    ! Local variables
    INTEGER(I32)                     :: fu, status, ios, n_rows, i, k, target_col, pass
    CHARACTER(len=2048)              :: line_buffer
    REAL(WP)           , ALLOCATABLE :: ws_array(:), var_array(:), row_buff(:)
    REAL(WP)                         :: dx, t
    LOGICAL                          :: is_row


    val        = 0.0_WP
    target_col = 2_I32
    if (present(col_idx)) target_col = col_idx

    if (target_col < 1_I32) error stop "read_curve: col_idx must be a positive column index"

    open(newunit=fu, file=trim(filename), status='old', action='read', iostat=status)
    if (status /= 0_I32) error stop "read_curve: could not open curve CSV file " // trim(filename)

    ! Pass 1 counts the data rows, pass 2 stores them (same row test in both)
    ALLOCATE(row_buff(target_col))
    n_rows = 0_I32
    do pass = 1, 2
        if (pass == 2_I32) then
            if (n_rows < 2_I32) then
                close(fu)
                error stop "read_curve: fewer than 2 numeric rows for interpolation in " // trim(filename)
            end if
            ALLOCATE(ws_array(n_rows), var_array(n_rows))
            rewind(fu)
            n_rows = 0_I32
        end if

        do
            read(fu, '(A)', iostat=ios) line_buffer
            if (ios /= 0_I32) exit

            CALL parse_csv_row(line_buffer, row_buff, is_row)
            if (.not. is_row) cycle

            n_rows = n_rows + 1_I32
            if (pass == 2_I32) then
                ws_array(n_rows)  = row_buff(1)
                var_array(n_rows) = row_buff(target_col)
            end if
        end do
    end do

    close(fu)

    ! Linear interpolation with extrapolation
    if (ws <= ws_array(1)) then
        k = 1_I32
    else if (ws >= ws_array(n_rows)) then
        k = n_rows - 1_I32
    else
        k = 1_I32
        do i = 1, n_rows - 1
            if (ws >= ws_array(i) .and. ws <= ws_array(i+1)) then
                k = i
                exit
            end if
        end do
    end if

    dx = ws_array(k+1) - ws_array(k)
    if (abs(dx) < 1.0e-12_WP) then
        t = 0.0_WP
    else
        t = (ws - ws_array(k)) / dx
    end if

    val = var_array(k) + t * (var_array(k+1) - var_array(k))

END FUNCTION read_curve


SUBROUTINE parse_csv_row(text_line, values, is_row)
    ! -------------------------------------------------------------------------
    ! Reads the first size(values) comma-separated fields of "text_line" as
    ! reals. is_row is .true. only if all of them exist and are numeric, so
    ! header lines, blank lines and short rows are all rejected the same way.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)  :: text_line    ! [-] Line of the CSV file
    REAL(WP)        , INTENT(OUT) :: values(:)    ! [-] Values of the first size(values) fields
    LOGICAL         , INTENT(OUT) :: is_row       ! [-] Whether all the fields were read

    ! Local variables
    INTEGER(I32) :: i, i_start, i_comma, ios


    is_row  = .false.
    values  = 0.0_WP
    i_start = 1_I32
    do i = 1, size(values)
        if (i_start > len_trim(text_line)) return

        i_comma = index(text_line(i_start:), ',')
        if (i_comma == 0_I32) then
            i_comma = len_trim(text_line) + 1_I32
        else
            i_comma = i_start + i_comma - 1_I32
        end if

        if (len_trim(text_line(i_start:i_comma-1)) == 0) return
        read(text_line(i_start:i_comma-1), *, iostat=ios) values(i)
        if (ios /= 0_I32) return

        i_start = i_comma + 1_I32
    end do
    is_row = .true.

END SUBROUTINE parse_csv_row


SUBROUTINE read_table(path, ncols, table)
    ! -------------------------------------------------------------------------
    ! Reads a text table of ncols numeric columns separated by blanks, TABs or
    ! commas. Blank lines and lines starting with '#' or '!' are skipped. A data
    ! row with another number of values, or with a value that is not a number,
    ! stops the run with the number of its line.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)               :: path          ! [-] Path to the text file
    INTEGER(I32)    , INTENT(IN)               :: ncols         ! [-] Number of columns of every data row
    REAL(WP)        , INTENT(OUT), ALLOCATABLE :: table(:,:)    ! [-] Values (number of data rows, ncols)

    ! Local variables
    CHARACTER(len=:) , ALLOCATABLE :: line
    CHARACTER(len=32), ALLOCATABLE :: tokens(:)
    INTEGER(I32)                   :: unit_num, ios, pass, n_rows, i_line, i, j


    open(newunit=unit_num, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0_I32) then
        write(*,*) "ERROR: Could not open the table file: ", path
        error stop
    end if

    ! Two passes: the first counts the data rows, the second reads them
    do pass = 1, 2
        n_rows = 0_I32
        i_line = 0_I32
        if (pass == 2_I32) rewind(unit_num)

        do
            CALL read_full_line(unit_num, line, ios)
            if (ios /= 0_I32) exit
            i_line = i_line + 1_I32

            ! Commas and TABs separate values as blanks do
            do i = 1, len(line)
                if (line(i:i) == ',' .or. line(i:i) == achar(9)) line(i:i) = ' '
            end do

            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#' .or. line(1:1) == '!') cycle

            n_rows = n_rows + 1_I32
            if (pass == 1_I32) cycle

            CALL split_tokens(line, tokens)
            if (size(tokens) /= ncols) then
                write(*, '(A,I0,A,I0,A,I0,A)') " ERROR: line ", i_line, " of " // trim(path) // " has ", size(tokens), &
                                               " values, ", ncols, " are required"
                close(unit_num)
                error stop
            end if

            do j = 1, ncols
                read(tokens(j), *, iostat=ios) table(n_rows, j)
                if (ios /= 0_I32) then
                    write(*, '(A,I0,A,I0,A)') " ERROR: line ", i_line, " of " // trim(path) // ", column ", j, &
                                              " is not a number"
                    close(unit_num)
                    error stop
                end if
            end do
        end do

        if (pass == 1_I32) then
            if (n_rows == 0_I32) then
                write(*,*) "ERROR: No data rows in the table file: ", path
                close(unit_num)
                error stop
            end if
            ALLOCATE(table(n_rows, ncols))
        end if
    end do

    close(unit_num)

END SUBROUTINE read_table

END MODULE CSVReader
