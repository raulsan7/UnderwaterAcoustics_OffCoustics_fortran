! ===============================================================================
! MODULE: FileUtils
! ===============================================================================
! Defines small file helpers: reading a whole record of any length
! (read_full_line), swapping the extension of a path (swap_extension),
! creating the parent folder of a path (ensure_parent_dir), moving a text
! file (move_text_file) and deleting a file if it exists (delete_file).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE FileUtils

USE Kinds, ONLY: I32

IMPLICIT NONE

PRIVATE
PUBLIC :: read_full_line, swap_extension, ensure_parent_dir, move_text_file, delete_file


CONTAINS


SUBROUTINE read_full_line(unit, line, ios)
    ! -------------------------------------------------------------------------
    ! Reads the next record of a formatted sequential unit into "line",
    ! whatever its length, with non-advancing reads of 4096-character pieces.
    ! -------------------------------------------------------------------------

    INTEGER(I32)    , INTENT(IN)               :: unit    ! [-] Open formatted sequential unit
    CHARACTER(len=:), INTENT(OUT), ALLOCATABLE :: line    ! [-] Whole record
    INTEGER(I32)    , INTENT(OUT)              :: ios     ! [-] 0 on success, nonzero at end of file or on error

    ! Local variables
    CHARACTER(len=4096) :: piece
    INTEGER(I32)        :: n_read


    line = ''
    do
        read(unit, '(A)', advance='no', iostat=ios, size=n_read) piece
        line = line // piece(1:n_read)
        if (is_iostat_eor(ios)) then
            ios = 0_I32
            exit
        end if
        if (ios /= 0_I32) exit
    end do

END SUBROUTINE read_full_line


PURE FUNCTION swap_extension(path, new_ext) RESULT(new_path)
    ! -------------------------------------------------------------------------
    ! Replaces the extension of "path" (everything after the last '.' of the
    ! file name, i.e. of the part after the last '/') with "new_ext". If the
    ! file name has no extension, or its only '.' is the leading one of a
    ! hidden file (".name"), "new_ext" is simply appended.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)  :: path        ! [-] Original file path
    CHARACTER(len=*), INTENT(IN)  :: new_ext     ! [-] New extension, including the leading '.'
    CHARACTER(len=:), ALLOCATABLE :: new_path    ! [-] Path with the new extension

    ! Local variables
    INTEGER(I32) :: dot_pos, slash_pos


    slash_pos = index(trim(path), '/', back=.true.)
    dot_pos   = index(trim(path), '.', back=.true.)

    if (dot_pos > slash_pos + 1_I32) then
        new_path = trim(path(1:dot_pos-1)) // trim(new_ext)
    else
        new_path = trim(path) // trim(new_ext)
    end if

END FUNCTION swap_extension


SUBROUTINE ensure_parent_dir(path)
    ! -------------------------------------------------------------------------
    ! Given a file path, creates its parent directory (and any missing
    ! intermediate directories) if it doesn't already exist. No-op if "path"
    ! has no directory component (i.e. it's a bare filename) or its parent is
    ! the root directory. The directory reaches the shell single-quoted (with
    ! any embedded quote escaped), so no character of the path is expanded,
    ! and a failed mkdir stops the program.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: path    ! [-] File path whose parent directory must exist

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: dir, quoted
    INTEGER(I32)                  :: slash_pos, i, exit_stat, cmd_stat


    slash_pos = index(trim(path), '/', back=.true.)
    if (slash_pos > 1_I32) then
        dir = trim(path(1:slash_pos-1))

        ! Single-quote the directory for the shell: each ' becomes '\''
        quoted = "'"
        do i = 1, len(dir)
            if (dir(i:i) == "'") then
                quoted = quoted // "'\''"
            else
                quoted = quoted // dir(i:i)
            end if
        end do
        quoted = quoted // "'"

        CALL execute_command_line('mkdir -p -- ' // quoted, exitstat=exit_stat, cmdstat=cmd_stat)
        if (cmd_stat /= 0_I32 .or. exit_stat /= 0_I32) then
            error stop "ensure_parent_dir: could not create directory " // dir
        end if
    end if

END SUBROUTINE ensure_parent_dir


SUBROUTINE move_text_file(src, dst)
    ! -------------------------------------------------------------------------
    ! Moves a text file: copies it line by line to "dst" (replacing any
    ! previous file there) and deletes "src". Portable (no RENAME extension).
    ! No-op if "src" does not exist.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: src    ! [-] Existing text file
    CHARACTER(len=*), INTENT(IN) :: dst    ! [-] Destination path

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: line
    INTEGER(I32)                  :: u_in, u_out, ios


    open(newunit=u_in, file=src, status='old', action='read', iostat=ios)
    if (ios /= 0_I32) return

    open(newunit=u_out, file=dst, status='replace', action='write')
    do
        CALL read_full_line(u_in, line, ios)
        if (ios /= 0_I32) exit
        write(u_out, '(A)') line
    end do
    close(u_out)
    close(u_in, status='delete')

END SUBROUTINE move_text_file


SUBROUTINE delete_file(path)
    ! -------------------------------------------------------------------------
    ! Deletes a file if it exists.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: path    ! [-] File to delete

    ! Local variables
    INTEGER(I32) :: u, ios


    open(newunit=u, file=path, status='old', iostat=ios)
    if (ios == 0_I32) close(u, status='delete')

END SUBROUTINE delete_file

END MODULE FileUtils
