! ===============================================================================
! MODULE: IOUtils
! ===============================================================================
! Handles all file I/O for the codebase: parsing OpenFAST SubDyn summary files
! (node/member topology), parsing OpenFAST binary/ASCII output channel files,
! a generic HDF5 variable-saving facility (the save_variable dispatcher, whose
! CLASS(*) assumed-rank argument is forwarded to the save_to_hdf5 generic
! interface: real/integer/complex data of rank 0 to 3, plus logical and
! character scalars; complex arrays are stored as (2, N...) in Fortran order,
! i.e. (..., N, 2) when read with h5py), small string helpers (to_upper,
! swap_extension, ensure_parent_dir), and a curve-lookup helper (read_curve).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE IOUtils

USE HDF5
USE Kinds, ONLY: I32, WP
USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, IEEE_QUIET_NAN

IMPLICIT NONE

PRIVATE
PUBLIC :: get_SDsum_variables, read_input_SD, read_curve, save_variable, save_to_hdf5, to_upper, swap_extension, &
          ensure_parent_dir, move_text_file, read_hdf5_real_1d, read_full_line


! ===============================================================================
! HDF5 GENERIC INTERFACE
! ===============================================================================
INTERFACE save_to_hdf5
    ! -------------------------------------------------------------------------
    ! Dispatches to the correct save_hdf5_<type>_<rank> specific procedure based
    ! on the actual type/rank of the variable passed in. This is the low-level
    ! layer used internally by the save_variable dispatcher below.
    ! -------------------------------------------------------------------------

    MODULE PROCEDURE save_hdf5_real_scalar
    MODULE PROCEDURE save_hdf5_int_scalar
    MODULE PROCEDURE save_hdf5_cplx_scalar
    MODULE PROCEDURE save_hdf5_logical_scalar
    MODULE PROCEDURE save_hdf5_char_scalar
    MODULE PROCEDURE save_hdf5_real_1d
    MODULE PROCEDURE save_hdf5_int_1d
    MODULE PROCEDURE save_hdf5_cplx_1d
    MODULE PROCEDURE save_hdf5_real_2d
    MODULE PROCEDURE save_hdf5_int_2d
    MODULE PROCEDURE save_hdf5_cplx_2d
    MODULE PROCEDURE save_hdf5_real_3d
    MODULE PROCEDURE save_hdf5_int_3d
    MODULE PROCEDURE save_hdf5_cplx_3d

END INTERFACE save_to_hdf5


CONTAINS


! ===============================================================================
! OPENFAST READERS
! ===============================================================================
SUBROUTINE get_SDsum_variables(SD_path, Nmembers, Nnodes, verbose, Nodes, status)
    ! -------------------------------------------------------------------------
    ! Parses an OpenFAST SubDyn summary file (SD.sum.yaml) to recover the
    ! structural topology: it locates the "#Member I Joint1_ID Joint2_ID"
    ! table to map each member to its node IDs, then locates the node
    ! coordinate table ("Node_[#] X_[m] Y_[m] Z_[m]") to recover each node's
    ! (x,y,z) position, and finally combines both into a (Nmembers, Nnodes, 3)
    ! tensor, padding any missing node with NaN. Any read or format error
    ! (including a node ID outside the node table) sets status = -1 when
    ! status is present and stops the program otherwise, so that a failed
    ! read never returns a geometry collapsed to the origin.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)               :: SD_path         ! [-] Path to OpenFAST SubDyn sum file
    INTEGER(I32)    , INTENT(IN), OPTIONAL     :: Nmembers        ! [-] Number of OpenFAST members
    INTEGER(I32)    , INTENT(IN), OPTIONAL     :: Nnodes          ! [-] Number of OpenFAST nodes
    LOGICAL         , INTENT(IN), OPTIONAL     :: verbose         ! [-] Flag to print more info
    REAL(WP)        , INTENT(OUT), ALLOCATABLE :: Nodes(:,:,:)    ! [m] Nodes array (Nmembers, Nnodes, 3)
    INTEGER(I32)    , INTENT(OUT), OPTIONAL    :: status          ! [-] 0 on success, -1 on error (stops if absent)

    ! Local variables
    INTEGER(I32) :: Nmembers_, Nnodes_
    LOGICAL      :: verbose_

    INTEGER(I32)       , ALLOCATABLE :: member_nodes(:,:)
    REAL(WP)           , ALLOCATABLE :: Nodes_flat(:,:)
    INTEGER(I32)                     :: ierr, file_unit, ios, i, j, k, idx, num_nodes
    INTEGER(I32)                     :: pos_open, pos_close
    REAL(WP)                         :: nan_val, x, y, z, idx_real
    LOGICAL                          :: found
    CHARACTER(len=1024)              :: line_buf
    CHARACTER(len=1024)              :: data_str


    ! Defaults
    Nmembers_ = 8_I32   ; if (present(Nmembers)) Nmembers_ = Nmembers
    Nnodes_   = 5_I32   ; if (present(Nnodes))   Nnodes_   = Nnodes
    verbose_  = .false. ; if (present(verbose))  verbose_  = verbose

    ! Initialize IEEE NaN for padding missing nodes
    nan_val = ieee_value(1.0_WP, IEEE_QUIET_NAN)
    if (present(status)) status = 0_I32

    if (verbose_) write(*,*) "Reading SD.sum nodes in ", trim(SD_path)

    ALLOCATE(Nodes(Nmembers_, Nnodes_, 3), source=0.0_WP)
    ALLOCATE(member_nodes(Nmembers_, Nnodes_), source=0_I32)

    open(newunit=file_unit, file=trim(SD_path), status='old', action='read', iostat=ierr)
    if (ierr /= 0_I32) then
        write(*,*) "Error: Could not open ", trim(SD_path)
        CALL set_error_status("get_SDsum_variables", SD_path, status)
        return
    end if

    ! --- 1. Extract connection map: member -> node IDs --- !
    found = .false.
    do while (.true.)
        read(file_unit, '(A)', iostat=ios) line_buf
        if (ios /= 0_I32) exit
        if (index(line_buf, "#Member I Joint1_ID Joint2_ID") > 0) then
            found = .true.
            exit
        end if
    end do

    if (.not. found) then
        write(*,*) "Error: Member mapping header not found."
        close(file_unit)
        CALL set_error_status("get_SDsum_variables", SD_path, status)
        return
    end if

    ! Parse node IDs for each member
    do i = 1, Nmembers_
        read(file_unit, '(A)', iostat=ios) line_buf
        if (ios /= 0_I32) then
            write(*,*) "Error: Member table ends before member ", i
            close(file_unit)
            CALL set_error_status("get_SDsum_variables", SD_path, status)
            return
        end if
        CALL parse_member_nodes(line_buf, Nnodes_, member_nodes(i, :))
    end do

    ! --- 2. Extract global node table: ID -> (x, y, z) --- !
    rewind(file_unit)
    found = .false.
    do while (.true.)
        read(file_unit, '(A)', iostat=ios) line_buf
        if (ios /= 0_I32) exit
        if (index(line_buf, '#     Node_[#]          X_[m]           Y_[m]           Z_[m]') > 0) then
            found = .true.
            exit
        end if
    end do

    if (.not. found) then
        write(*,*) "Error: Node coordinate header not found."
        close(file_unit)
        CALL set_error_status("get_SDsum_variables", SD_path, status)
        return
    end if

    ! Read the number of nodes (e.g., "Nodes: # 33 x 9")
    read(file_unit, '(A)', iostat=ios) line_buf
    if (ios /= 0_I32) then
        write(*,*) "Error: Node table truncated in ", trim(SD_path)
        close(file_unit)
        CALL set_error_status("get_SDsum_variables", SD_path, status)
        return
    end if
    num_nodes = extract_number_nodes(line_buf)

    ALLOCATE(Nodes_flat(num_nodes, 3), source=0.0_WP)

    ! Read coordinates: each line looks like
    !   - [      1.,          0.000,          0.000,        -30.000, ... ]
    ! Extract only the content between '[' and ']' before the numeric read,
    ! so the leading '-' YAML list marker never gets fed to the parser
    ! (removing it globally would corrupt negative coordinate values)
    do i = 1, num_nodes
        read(file_unit, '(A)', iostat=ios) line_buf
        if (ios /= 0_I32) then
            write(*,*) "Error: Node table truncated in ", trim(SD_path)
            close(file_unit)
            CALL set_error_status("get_SDsum_variables", SD_path, status)
            return
        end if

        pos_open  = index(line_buf, '[')
        pos_close = index(line_buf, ']')
        if (pos_open == 0_I32 .or. pos_close == 0_I32 .or. pos_close <= pos_open) then
            write(*,*) "Error: malformed node coordinate line: ", trim(line_buf)
            close(file_unit)
            CALL set_error_status("get_SDsum_variables", SD_path, status)
            return
        end if

        data_str = line_buf(pos_open+1:pos_close-1)

        ! Clear commas so list-directed read treats fields as blank-separated
        do k = 1, len_trim(data_str)
            if (data_str(k:k) == ',') data_str(k:k) = ' '
        end do

        ! Node ID is stored as a real literal (e.g. "1."), read as real then
        ! convert, mirroring the Python int(float(...)) pattern
        read(data_str, *, iostat=ios) idx_real, x, y, z
        if (ios /= 0_I32) then
            write(*,*) "Error: could not parse node coordinate line: ", trim(line_buf)
            close(file_unit)
            CALL set_error_status("get_SDsum_variables", SD_path, status)
            return
        end if

        idx = nint(idx_real)
        if (idx >= 1_I32 .and. idx <= num_nodes) then
            Nodes_flat(idx, :) = [x, y, z]
        else
            write(*,*) "Error: node index ", idx, " out of range [1, ", num_nodes, "]"
            close(file_unit)
            CALL set_error_status("get_SDsum_variables", SD_path, status)
            return
        end if
    end do

    close(file_unit)

    ! --- 3. Map flat nodes to output tensor (Nmembers, Nnodes, 3) --- !
    do i = 1, Nmembers_
        do j = 1, Nnodes_
            idx = member_nodes(i, j)
            if (idx < 0_I32 .or. idx > num_nodes) then
                write(*,*) "Error: member ", i, " refers to node ", idx, " outside the node table [1, ", &
                           num_nodes, "]"
                CALL set_error_status("get_SDsum_variables", SD_path, status)
                return
            else if (idx > 0_I32) then
                Nodes(i, j, :) = Nodes_flat(idx, :)
            else
                Nodes(i, j, :) = nan_val
            end if
        end do
    end do

END SUBROUTINE get_SDsum_variables


SUBROUTINE set_error_status(routine, path, status)
    ! -------------------------------------------------------------------------
    ! Reports a read failure: sets status = -1 when the caller passed it, and
    ! otherwise stops the program, so that a failed read never goes unnoticed.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)            :: routine    ! [-] Name of the calling routine
    CHARACTER(len=*), INTENT(IN)            :: path       ! [-] File that could not be read
    INTEGER(I32)    , INTENT(OUT), OPTIONAL :: status     ! [-] Set to -1 when present


    if (present(status)) then
        status = -1_I32
    else
        error stop routine // ": could not read " // trim(path)
    end if

END SUBROUTINE set_error_status


SUBROUTINE read_input_SD(filename, what, skip, Nmembers, Nnodes, From, Upto, verbose, &
                         Time_out, Array_out, unit_out, status)
    ! -------------------------------------------------------------------------
    ! Reads a time series of a requested structural quantity ("displacement",
    ! "acceleration", "force", "momentum", etc.) for every member/node from an
    ! OpenFAST output file. It builds the corresponding "MiNjXXX"-style channel
    ! names for all (member, node, xyz) combinations, delegates the actual file
    ! reading to parse_channels_auto, then decimates (skip) and reshapes the
    ! flat channel array into a (nt, Nmembers, Nnodes, 3) tensor.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)               :: filename              ! [-] Path to OpenFAST output file
    CHARACTER(len=*), INTENT(IN), OPTIONAL     :: what                  ! [-] Name of the output to read
    INTEGER(I32)    , INTENT(IN), OPTIONAL     :: skip                  ! [-] Decimation factor
    INTEGER(I32)    , INTENT(IN), OPTIONAL     :: Nmembers              ! [-] Total structural members
    INTEGER(I32)    , INTENT(IN), OPTIONAL     :: Nnodes                ! [-] Nodes per member
    REAL(WP)        , INTENT(IN), OPTIONAL     :: From                  ! [-] Start window fraction
    REAL(WP)        , INTENT(IN), OPTIONAL     :: Upto                  ! [-] End window fraction
    LOGICAL         , INTENT(IN), OPTIONAL     :: verbose               ! [-] Flag to print info
    REAL(WP)        , INTENT(OUT), ALLOCATABLE :: Time_out(:)           ! [s] Downsampled time array (nt)
    REAL(WP)        , INTENT(OUT), ALLOCATABLE :: Array_out(:,:,:,:)    ! [-] Data (nt, Nm, Nn, 3) in unit_out
    CHARACTER(len=*), INTENT(OUT)              :: unit_out              ! [-] Physical unit of Array_out
    INTEGER(I32)    , INTENT(OUT)              :: status                ! [-] 0 on success, nonzero on error

    ! Local variables
    CHARACTER(len=32) :: what_
    INTEGER(I32)      :: skip_, Nmembers_, Nnodes_
    REAL(WP)          :: From_, Upto_
    LOGICAL           :: verbose_

    ! Channel names, raw data and work indices
    CHARACTER(len=16)              :: strx, stry, strz
    CHARACTER(len=32), ALLOCATABLE :: out_channels(:), raw_units(:)
    REAL(WP)         , ALLOCATABLE :: raw_time(:), raw_array(:,:)
    INTEGER(I32)                   :: m, n, ch_idx, nt, nt_skip, t_skip, t_orig, n_chans, base, idx
    CHARACTER(len=32)              :: ch_name


    status = 0_I32

    ! Apply optional arguments
    what_     = "acceleration" ; if (present(what))     what_     = trim(what)
    skip_     = 1_I32          ; if (present(skip))     skip_     = skip
    Nmembers_ = 8_I32          ; if (present(Nmembers)) Nmembers_ = Nmembers
    Nnodes_   = 5_I32          ; if (present(Nnodes))   Nnodes_   = Nnodes
    From_     = 0.0_WP         ; if (present(From))     From_     = From
    Upto_     = 1.0_WP         ; if (present(Upto))     Upto_     = Upto
    verbose_  = .false.        ; if (present(verbose))  verbose_  = verbose

    ! --- 1. Determine channel suffixes based on 'what' --- !
    SELECT CASE (trim(what_))
    CASE ("displacement")
        strx = "TDxss"
        stry = "TDyss"
        strz = "TDzss"
    CASE ("acceleration")
        strx = "TAxe"
        stry = "TAye"
        strz = "TAze"
    CASE ("force")
        strx = "FKxe"
        stry = "FKye"
        strz = "FKze"
    CASE ("momentum")
        strx = "MKxe"
        stry = "MKye"
        strz = "MKze"
    CASE ("int_displacement")
        strx = "IntfTDXss"
        stry = "IntfTDYss"
        strz = "IntfTDZss"
    CASE ("int_acceleration")
        strx = "IntfTAXss"
        stry = "IntfTAYss"
        strz = "IntfTAZss"
    CASE DEFAULT
        write(*,*) "Error: Invalid option for 'what' = ", trim(what_)
        write(*,*) "Input 'displacement', 'acceleration', 'force', 'momentum', 'int_displacement' or 'int_acceleration'"
        status = -1_I32
        return
    END SELECT

    if (verbose_) then
        write(*,*) "Reading ", trim(filename), " --> ", trim(what_), &
                   ": [MiNj", trim(strx), ", MiNj", trim(stry), ", MiNj", trim(strz), "]"
    end if

    ! --- 2. Dynamically generate channel names --- !
    n_chans = Nmembers_ * Nnodes_ * 3_I32
    ALLOCATE(out_channels(n_chans))
    ch_idx = 1_I32

    do m = 1, Nmembers_
        do n = 1, Nnodes_
            ! X component
            write(ch_name, '("M",I0,"N",I0,A)') m, n, trim(strx)
            out_channels(ch_idx) = trim(ch_name)
            ch_idx               = ch_idx + 1_I32
            ! Y component
            write(ch_name, '("M",I0,"N",I0,A)') m, n, trim(stry)
            out_channels(ch_idx) = trim(ch_name)
            ch_idx               = ch_idx + 1_I32
            ! Z component
            write(ch_name, '("M",I0,"N",I0,A)') m, n, trim(strz)
            out_channels(ch_idx) = trim(ch_name)
            ch_idx               = ch_idx + 1_I32
        end do
    end do

    ! --- 3. Call the universal OpenFAST parser --- !
    CALL parse_channels_auto(full_path=filename, plotChannels=out_channels, &
                             From=From_, Upto=Upto_, &
                             available_channels=.false., verbose=verbose_, &
                             time_out=raw_time, array_out=raw_array, &
                             units_out=raw_units, status=status)
    if (status /= 0_I32) return

    nt = size(raw_time)
    if (nt == 0_I32) then
        write(*,*) "Warning: No data returned from parser."
        status = -1_I32
        return
    end if

    ! --- 4. Decimate (skip) and reshape into 4D tensor --- !
    nt_skip = (nt - 1_I32) / skip_ + 1_I32
    if (nt_skip < 2_I32) then
        write(*,*) "Error: At least 2 time samples are needed after decimation, got ", nt_skip
        status = -1_I32
        return
    end if

    ALLOCATE(Time_out(nt_skip))
    ALLOCATE(Array_out(nt_skip, Nmembers_, Nnodes_, 3))

    t_skip = 1_I32
    do t_orig = 1, nt, skip_
        Time_out(t_skip) = raw_time(t_orig)
        do m = 1, Nmembers_
            base = (m - 1_I32) * Nnodes_ * 3_I32 + 1_I32
            !$OMP SIMD PRIVATE(idx)
            do n = 1, Nnodes_
                idx = base + (n - 1_I32) * 3_I32
                Array_out(t_skip, m, n, 1) = raw_array(t_orig, idx)
                Array_out(t_skip, m, n, 2) = raw_array(t_orig, idx+1)
                Array_out(t_skip, m, n, 3) = raw_array(t_orig, idx+2)
            end do
        end do
        t_skip = t_skip + 1_I32
    end do

    ! Extract unit from the first valid channel
    unit_out = raw_units(1)

END SUBROUTINE read_input_SD


SUBROUTINE parse_channels_auto(full_path, plotChannels, From, Upto, &
                               available_channels, verbose, chunk_size, &
                               time_out, array_out, units_out, status)
    ! -------------------------------------------------------------------------
    ! Universal entry point for reading OpenFAST output channels: it inspects
    ! the file extension (".outb" vs ".out") and dispatches to the matching
    ! specialized parser, parse_channels_binary or parse_channels_ascii, so
    ! callers do not need to know which format a given run produced.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*) , INTENT(IN)               :: full_path             ! [-] Path to the OpenFAST output file
    CHARACTER(len=*) , INTENT(IN)               :: plotChannels(:)       ! [-] Names of the channels to read
    REAL(WP)         , INTENT(IN), OPTIONAL     :: From                  ! [-] Start window fraction
    REAL(WP)         , INTENT(IN), OPTIONAL     :: Upto                  ! [-] End window fraction
    LOGICAL          , INTENT(IN), OPTIONAL     :: available_channels    ! [-] Whether to print all channels in the file
    LOGICAL          , INTENT(IN), OPTIONAL     :: verbose               ! [-] Flag to print info
    INTEGER(I32)     , INTENT(IN), OPTIONAL     :: chunk_size            ! [-] Rows read per chunk (binary files only)
    REAL(WP)         , INTENT(OUT), ALLOCATABLE :: time_out(:)           ! [s] Time vector of the selected window
    REAL(WP)         , INTENT(OUT), ALLOCATABLE :: array_out(:,:)        ! [-] Channel data (nt, Nchannels)
    CHARACTER(len=32), INTENT(OUT), ALLOCATABLE :: units_out(:)          ! [-] Unit of each requested channel
    INTEGER(I32)     , INTENT(OUT)              :: status                ! [-] 0 on success, nonzero on error

    ! Local variables
    REAL(WP)     :: From_, Upto_
    LOGICAL      :: avail_ch, verb
    INTEGER(I32) :: c_size, ext_idx


    ! Defaults
    From_    = 0.0_WP   ; if (present(From))               From_    = From
    Upto_    = 1.0_WP   ; if (present(Upto))               Upto_    = Upto
    avail_ch = .false.  ; if (present(available_channels)) avail_ch = available_channels
    verb     = .true.   ; if (present(verbose))            verb     = verbose
    c_size   = 4096_I32 ; if (present(chunk_size))         c_size   = chunk_size
    status = 0_I32

    ! Detect extension
    ext_idx = index(full_path, '.outb', back=.true.)
    if (ext_idx > 0_I32 .and. ext_idx == len_trim(full_path) - 4_I32) then
        if (verb) write(*,*) "Detected binary OpenFAST output (.outb): ", trim(full_path)
        CALL parse_channels_binary(full_path, plotChannels, From_, Upto_, &
                                   avail_ch, verb, c_size, time_out, array_out, units_out, status)
    else
        ext_idx = index(full_path, '.out', back=.true.)
        if (ext_idx > 0_I32 .and. ext_idx == len_trim(full_path) - 3_I32) then
            if (verb) write(*,*) "Detected text OpenFAST output (.out): ", trim(full_path)
            CALL parse_channels_ascii(full_path, plotChannels, From_, Upto_, &
                                      avail_ch, verb, time_out, array_out, units_out, status)
        else
            write(*,*) "Error: Unsupported file extension. Expected .out or .outb."
            status = -1_I32
        end if
    end if

END SUBROUTINE parse_channels_auto


SUBROUTINE parse_channels_binary(full_path, plotChannels, From_val, Upto_val, &
                                 available_channels, verbose, chunk_size, &
                                 time_out, array_out, units_out, status)
    ! -------------------------------------------------------------------------
    ! Parses OpenFAST's packed binary output format (.outb), as written by
    ! WrBinFAST (NWTC_IO.f90): reads the fixed header (FileID, channel
    ! count/names/units, time and per-column scale and offset factors used to
    ! decompress the packed INTEGER(2) values), maps the requested channel
    ! names to column indices, then streams the data body in chunks (to bound
    ! memory use), decoding only the requested time window and channels back
    ! into REAL(WP) using the scale/offset factors. With FileID = 1 the packed
    ! INTEGER(4) times are one block of NT values stored before the INTEGER(2)
    ! channel matrix. The window must satisfy 0 <= From_val <= Upto_val <= 1,
    ! and every row of the window must be read, otherwise status is -1.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*) , INTENT(IN)               :: full_path             ! [-] Path to the .outb file
    CHARACTER(len=*) , INTENT(IN)               :: plotChannels(:)       ! [-] Names of the channels to read
    REAL(WP)         , INTENT(IN)               :: From_val              ! [-] Start window fraction
    REAL(WP)         , INTENT(IN)               :: Upto_val              ! [-] End window fraction
    LOGICAL          , INTENT(IN)               :: available_channels    ! [-] Whether to print all channels in the file
    LOGICAL          , INTENT(IN)               :: verbose               ! [-] Flag to print info
    INTEGER(I32)     , INTENT(IN)               :: chunk_size            ! [-] Rows read per chunk
    REAL(WP)         , INTENT(OUT), ALLOCATABLE :: time_out(:)           ! [s] Time vector of the selected window
    REAL(WP)         , INTENT(OUT), ALLOCATABLE :: array_out(:,:)        ! [-] Channel data (nt, Nchannels)
    CHARACTER(len=32), INTENT(OUT), ALLOCATABLE :: units_out(:)          ! [-] Unit of each requested channel
    INTEGER(I32)     , INTENT(OUT)              :: status                ! [-] 0 on success, nonzero on error

    ! Local variables
    ! OpenFAST binary format constants
    INTEGER(2), PARAMETER :: FMT_WITH_TIME    = 1
    INTEGER(2), PARAMETER :: FMT_WITHOUT_TIME = 2
    INTEGER(2), PARAMETER :: FMT_NO_COMPRESS  = 3
    INTEGER(2), PARAMETER :: FMT_CHAN_LEN     = 4

    ! File header variables
    INTEGER(I32)                   :: fu, i, j, k, req_chans, ios
    INTEGER(2)                     :: FileID, LenName
    INTEGER(I32)                   :: NumOutChans, NT
    REAL(8)                        :: TimeScl, TimeOff, TimeOut1, TimeIncr
    REAL(4)          , ALLOCATABLE :: ColScl(:), ColOff(:)
    INTEGER(I32)                   :: LenDesc
    CHARACTER(len=:) , ALLOCATABLE :: DescStr
    CHARACTER(len=32), ALLOCATABLE :: ChanName(:), ChanUnit(:)
    INTEGER(I32)     , ALLOCATABLE :: ReqIdx(:)
    INTEGER(4)       , ALLOCATABLE :: packed_time(:)

    ! Array streaming variables
    INTEGER(I32)              :: idx_0_start, idx_0_end, start_idx, end_idx, n_target_rows
    INTEGER(I32)              :: row_counter, written, chunk_rows, rows_to_read, current_row
    INTEGER(2)                :: val_i16
    INTEGER(2)  , ALLOCATABLE :: raw_data(:,:)


    status    = 0_I32
    req_chans = size(plotChannels)

    if (From_val < 0.0_WP .or. Upto_val > 1.0_WP .or. From_val > Upto_val) then
        write(*,*) "Error: Window fractions must satisfy 0 <= From <= Upto <= 1, got ", From_val, Upto_val
        status = -1_I32
        return
    end if

    open(newunit=fu, file=trim(full_path), access='stream', form='unformatted', status='old', iostat=status)
    if (status /= 0_I32) then
        write(*,*) "Error: Could not open binary file: ", trim(full_path)
        return
    end if

    ! --- 1. Read header --- !
    read(fu, iostat=ios) FileID

    if (ios == 0_I32 .and. FileID == FMT_CHAN_LEN) then
        read(fu, iostat=ios) LenName
    else
        LenName = 10_2
    end if

    if (ios == 0_I32) read(fu, iostat=ios) NumOutChans
    if (ios == 0_I32) read(fu, iostat=ios) NT
    if (ios /= 0_I32) then
        write(*,*) "Error: Could not read the header of binary file: ", trim(full_path)
        status = -1_I32
        close(fu)
        return
    end if

    if (FileID == FMT_NO_COMPRESS) then
        write(*,*) "Error: Uncompressed float64 binaries are not supported."
        status = -1_I32
        close(fu)
        return
    end if

    if (FileID /= FMT_WITH_TIME .and. FileID /= FMT_WITHOUT_TIME .and. FileID /= FMT_CHAN_LEN) then
        write(*,*) "Error: Unknown FileID ", FileID, " in binary file: ", trim(full_path)
        status = -1_I32
        close(fu)
        return
    end if

    if (LenName < 1_2 .or. LenName > 32_2 .or. NumOutChans < 1_I32 .or. NT < 0_I32) then
        write(*,*) "Error: Invalid header (LenName, NumOutChans, NT) = ", LenName, NumOutChans, NT, &
                   " in binary file: ", trim(full_path)
        status = -1_I32
        close(fu)
        return
    end if

    ! --- 2. Time info --- !
    if (FileID == FMT_WITH_TIME) then
        read(fu, iostat=ios) TimeScl
        if (ios == 0_I32) read(fu, iostat=ios) TimeOff
    else
        read(fu, iostat=ios) TimeOut1
        if (ios == 0_I32) read(fu, iostat=ios) TimeIncr
    end if

    ! --- 3. Channel scaling --- !
    ALLOCATE(ColScl(NumOutChans), ColOff(NumOutChans))
    if (ios == 0_I32) read(fu, iostat=ios) ColScl
    if (ios == 0_I32) read(fu, iostat=ios) ColOff

    ! --- 4. Description string --- !
    if (ios == 0_I32) read(fu, iostat=ios) LenDesc
    if (ios == 0_I32 .and. LenDesc > 0_I32) then
        ALLOCATE(CHARACTER(len=LenDesc) :: DescStr)
        read(fu, iostat=ios) DescStr
        DEALLOCATE(DescStr)
    end if

    ! --- 5. Channel names and units --- !
    ALLOCATE(ChanName(NumOutChans+1))
    ALLOCATE(ChanUnit(NumOutChans+1))
    ALLOCATE(ReqIdx(req_chans))

    ChanName = ""
    ChanUnit = ""
    do i = 1, NumOutChans + 1
        if (ios == 0_I32) read(fu, iostat=ios) ChanName(i)(1:LenName)
    end do
    do i = 1, NumOutChans + 1
        if (ios == 0_I32) read(fu, iostat=ios) ChanUnit(i)(1:LenName)
    end do

    ! Packed time column (FileID = 1 only): one block of NT values before the data matrix
    if (FileID == FMT_WITH_TIME) then
        ALLOCATE(packed_time(NT))
        if (ios == 0_I32) read(fu, iostat=ios) packed_time
    end if

    if (ios /= 0_I32) then
        write(*,*) "Error: Truncated header in binary file: ", trim(full_path)
        status = -1_I32
        close(fu)
        return
    end if

    if (available_channels) then
        write(*,*) "Available Channels in .outb: "
        do i = 1, NumOutChans + 1
            write(*, '(2A)') " - ", trim(ChanName(i))
        end do
    end if

    ! --- 6. Map requested channels --- !
    ALLOCATE(units_out(req_chans))
    do i = 1, req_chans
        ReqIdx(i) = 0_I32
        do j = 1, NumOutChans + 1
            if (trim(adjustl(plotChannels(i))) == trim(adjustl(ChanName(j)))) then
                if (j == 1_I32) then
                    ReqIdx(i) = -1_I32       ! Marker for the time column
                else
                    ReqIdx(i) = j - 1_I32    ! Data channel offset (1 to NumOutChans)
                end if
                units_out(i) = ChanUnit(j)
                exit
            end if
        end do
        if (ReqIdx(i) == 0_I32) then
            write(*,*) "Error: Channel '", trim(plotChannels(i)), "' not found in .outb."
            status = -1_I32
            close(fu)
            return
        end if
    end do

    ! --- 7. Window filtering logic --- !
    idx_0_start = int(real(NT, WP) * From_val)
    idx_0_end   = int(real(NT, WP) * Upto_val)

    start_idx     = idx_0_start + 1_I32
    end_idx       = idx_0_end
    n_target_rows = max(0_I32, end_idx - start_idx + 1_I32)

    ALLOCATE(time_out(n_target_rows))
    ALLOCATE(array_out(n_target_rows, req_chans))

    ! --- 8. Read data as massive contiguous arrays --- !
    chunk_rows = max(1_I32, chunk_size)
    ALLOCATE(raw_data(NumOutChans, chunk_rows))

    row_counter = 0_I32
    written     = 0_I32

    do while (row_counter < NT .and. written < n_target_rows)
        rows_to_read = min(chunk_rows, NT - row_counter)
        read(fu, iostat=ios) raw_data(:, 1:rows_to_read)
        if (ios /= 0_I32) exit

        do i = 1, rows_to_read
            current_row = row_counter + i
            if (current_row >= start_idx .and. current_row <= end_idx) then
                written = written + 1_I32

                ! Decode time
                if (FileID == FMT_WITH_TIME) then
                    time_out(written) = real((real(packed_time(current_row), 8) - TimeOff) / TimeScl, WP)
                else
                    time_out(written) = real(TimeOut1 + TimeIncr * real(current_row - 1_I32, 8), WP)
                end if

                ! Assign time column if any (assume at most one)
                do j = 1, req_chans
                    if (ReqIdx(j) == -1_I32) then
                        array_out(written, j) = time_out(written)
                        exit
                    end if
                end do

                ! Decode data channels with SIMD
                !$OMP SIMD PRIVATE(k, val_i16)
                do j = 1, req_chans
                    if (ReqIdx(j) /= -1_I32) then
                        k       = ReqIdx(j)
                        val_i16 = raw_data(k, i)
                        array_out(written, j) = real((real(val_i16, 4) - ColOff(k)) / ColScl(k), WP)
                    end if
                end do
            end if
        end do

        row_counter = row_counter + rows_to_read
    end do

    close(fu)

    if (written /= n_target_rows) then
        write(*,*) "Error: Only ", written, " of the ", n_target_rows, " requested rows could be read from ", &
                   trim(full_path)
        status = -1_I32
        return
    end if

    if (verbose) then
        write(*,*) "Output array shape: [", written, ", ", req_chans, "]"
    end if

END SUBROUTINE parse_channels_binary


SUBROUTINE parse_channels_ascii(full_path, plotChannels, From_val, Upto_val, &
                                available_channels, verbose, time_out, &
                                array_out, units_out, status)
    ! -------------------------------------------------------------------------
    ! Parses OpenFAST's plain-text output format (.out): skips the 6 metadata
    ! header lines, reads the channel-name and channel-unit header lines
    ! (line 7/8, of any length, blank- or TAB-delimited) to map requested
    ! channel names to column positions and to locate the "Time" column, does
    ! a fast pre-pass to count the data rows (blank lines and '#' comment
    ! lines such as "#Restarting here" are not data rows, as in numpy.loadtxt)
    ! to size the From/Upto time window (0 <= From_val <= Upto_val <= 1), then
    ! rewinds and parses only the rows of the window. An empty window, a row
    ! that cannot be parsed or a missing row returns status = -1.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*) , INTENT(IN)               :: full_path             ! [-] Path to the .out file
    CHARACTER(len=*) , INTENT(IN)               :: plotChannels(:)       ! [-] Names of the channels to read
    REAL(WP)         , INTENT(IN)               :: From_val              ! [-] Start window fraction
    REAL(WP)         , INTENT(IN)               :: Upto_val              ! [-] End window fraction
    LOGICAL          , INTENT(IN)               :: available_channels    ! [-] Whether to print all channels in the file
    LOGICAL          , INTENT(IN)               :: verbose               ! [-] Flag to print info
    REAL(WP)         , INTENT(OUT), ALLOCATABLE :: time_out(:)           ! [s] Time vector of the selected window
    REAL(WP)         , INTENT(OUT), ALLOCATABLE :: array_out(:,:)        ! [-] Channel data (nt, Nchannels)
    CHARACTER(len=32), INTENT(OUT), ALLOCATABLE :: units_out(:)          ! [-] Unit of each requested channel
    INTEGER(I32)     , INTENT(OUT)              :: status                ! [-] 0 on success, nonzero on error

    ! Local variables
    INTEGER(I32)                     :: fu, i, j, k, ios, n_req, n_all_chans, n_total, time_col
    INTEGER(I32)                     :: idx_0_start, idx_0_end, start_idx, end_idx, n_filtered, out_row
    INTEGER(I32)                     :: row_idx
    CHARACTER(len=256)               :: head_buffer
    CHARACTER(len=:)   , ALLOCATABLE :: line
    CHARACTER(len=32)  , ALLOCATABLE :: all_channels(:), all_units(:)
    INTEGER(I32)       , ALLOCATABLE :: col_map(:)
    REAL(WP)           , ALLOCATABLE :: row_buff(:)


    status = 0_I32
    n_req  = size(plotChannels)

    if (n_req == 0_I32) then
        write(*,*) "Error: plotChannels list cannot be empty."
        status = -1_I32
        return
    end if

    if (From_val < 0.0_WP .or. Upto_val > 1.0_WP .or. From_val > Upto_val) then
        write(*,*) "Error: Window fractions must satisfy 0 <= From <= Upto <= 1, got ", From_val, Upto_val
        status = -1_I32
        return
    end if

    open(newunit=fu, file=trim(full_path), status='old', action='read', iostat=status)
    if (status /= 0_I32) then
        write(*,*) "Error: Could not open ASCII file: ", trim(full_path)
        return
    end if

    ! --- 1. Read header lines --- !
    ! Skip lines 1 to 6 (description/metadata)
    do i = 1, 6
        read(fu, '(A)', iostat=status) head_buffer
        if (status /= 0_I32) then
            write(*,*) "Error reading header lines in file: ", trim(full_path)
            close(fu)
            status = -1_I32
            return
        end if
    end do

    ! Line 7: channel names (whole record, blank- or TAB-delimited)
    CALL read_full_line(fu, line, ios)
    if (ios /= 0_I32) then
        close(fu)
        status = -1_I32
        return
    end if

    CALL split_tokens(line, all_channels)
    n_all_chans = size(all_channels)

    if (n_all_chans == 0_I32) then
        write(*,*) "Error: No channel names found on line 7."
        close(fu)
        status = -1_I32
        return
    end if

    ! Line 8: channel units (split like the names, so a '/' inside a unit is kept)
    CALL read_full_line(fu, line, ios)
    if (ios /= 0_I32) then
        close(fu)
        status = -1_I32
        return
    end if

    CALL split_tokens(line, all_units)
    if (size(all_units) /= n_all_chans) then
        write(*,*) "Error: ", size(all_units), " units found on line 8 for ", n_all_chans, " channels."
        close(fu)
        status = -1_I32
        return
    end if

    if (available_channels) then
        write(*,*) "Available Channels in .out: "
        do i = 1, n_all_chans
            write(*, '(2A)') " - ", trim(all_channels(i))
        end do
    end if

    ! --- 2. Map requested channels and locate time column --- !
    ALLOCATE(col_map(n_req))
    ALLOCATE(units_out(n_req))
    time_col = 1_I32    ! Default assumption

    do j = 1, n_all_chans
        if (trim(adjustl(all_channels(j))) == "Time") time_col = j
    end do

    do i = 1, n_req
        col_map(i) = 0_I32
        do j = 1, n_all_chans
            if (trim(adjustl(plotChannels(i))) == trim(adjustl(all_channels(j)))) then
                col_map(i)   = j
                units_out(i) = all_units(j)
                exit
            end if
        end do
        if (col_map(i) == 0_I32) then
            write(*,*) "Error: Channel '", trim(plotChannels(i)), "' not found in ASCII file."
            close(fu)
            status = -1_I32
            return
        end if
    end do

    ! --- 3. Fast count of total data rows --- !
    ! Only the first 256 characters of each record are read to classify it
    n_total = 0_I32
    do
        read(fu, '(A)', iostat=ios) head_buffer
        if (ios /= 0_I32) exit
        if (is_data_line(head_buffer)) n_total = n_total + 1_I32
    end do

    if (n_total == 0_I32) then
        write(*,*) "Error: No data rows found in ASCII file."
        close(fu)
        status = -1_I32
        return
    end if

    ! --- 4. Compute filtering windows (matching Python and binary logic) --- !
    idx_0_start = int(real(n_total, WP) * From_val)
    idx_0_end   = int(real(n_total, WP) * Upto_val)

    start_idx = idx_0_start + 1_I32
    end_idx   = idx_0_end

    n_filtered = max(0_I32, end_idx - start_idx + 1_I32)
    if (n_filtered <= 0_I32) then
        write(*,*) "Warning: Filter range resulted in 0 output rows."
        close(fu)
        status = -1_I32
        return
    end if

    ! --- 5. Extract time series data --- !
    rewind(fu)
    do i = 1, 8
        read(fu, '(A)') head_buffer
    end do

    ALLOCATE(time_out(n_filtered))
    ALLOCATE(array_out(n_filtered, n_req))
    ALLOCATE(row_buff(n_all_chans))

    out_row = 0_I32
    row_idx = 0_I32
    do while (row_idx < end_idx)
        ! Rows before the window are only classified, never parsed
        if (row_idx + 1_I32 < start_idx) then
            read(fu, '(A)', iostat=ios) head_buffer
            if (ios /= 0_I32) exit
            if (is_data_line(head_buffer)) row_idx = row_idx + 1_I32
            cycle
        end if

        CALL read_full_line(fu, line, ios)
        if (ios /= 0_I32) exit
        if (.not. is_data_line(line)) cycle

        row_idx = row_idx + 1_I32
        read(line, *, iostat=ios) row_buff
        if (ios /= 0_I32) then
            write(*,*) "Error: Could not parse data row ", row_idx, " in file: ", trim(full_path)
            close(fu)
            status = -1_I32
            return
        end if

        out_row           = out_row + 1_I32
        time_out(out_row) = row_buff(time_col)
        do k = 1, n_req
            array_out(out_row, k) = row_buff(col_map(k))
        end do
    end do

    close(fu)

    if (out_row /= n_filtered) then
        write(*,*) "Error: Only ", out_row, " of the ", n_filtered, " requested rows could be read from ", &
                   trim(full_path)
        status = -1_I32
        return
    end if

    if (verbose) then
        write(*,*) "Output array shape: [", n_filtered, ", ", n_req, "]"
    end if

END SUBROUTINE parse_channels_ascii


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


! ===============================================================================
! HELPER FUNCTIONS
! ===============================================================================
SUBROUTINE parse_member_nodes(text_line, max_nodes, node_arr)
    ! -------------------------------------------------------------------------
    ! Tokenizes a whitespace-separated SD.sum.yaml member-mapping line (of the
    ! form "#Member I Joint1_ID Joint2_ID Prop_I Prop_J Mass Length Node
    ! IDs..."), skips the fixed 7-column prefix, and extracts the remaining
    ! tokens as this member's integer node IDs. That trailing node-ID list can
    ! be SHORTER than max_nodes (e.g. exactly 2 for a two-noded "Rigid link"
    ! member, versus Nnodes for a fully-discretized one); missing slots are
    ! left as 0, later padded to NaN coordinates by the caller. Any trailing
    ! "# <comment>" annotation (e.g. "# Rigid link") is truncated away rather
    ! than blanked, so it never turns into spurious extra tokens.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)  :: text_line              ! [-] Line read in SD.sum.yaml
    INTEGER(I32)    , INTENT(IN)  :: max_nodes              ! [-] Amount of nodes per member and line
    INTEGER(I32)    , INTENT(OUT) :: node_arr(max_nodes)    ! [-] Node indices found at the end of text_line

    ! Local variables
    INTEGER(I32), PARAMETER :: max_tokens = 200_I32    ! [-] Safety upper bound on the number of tokens
    INTEGER(I32), PARAMETER :: token_len  = 32_I32     ! [-] Max characters per token
    INTEGER(I32), PARAMETER :: n_prefix   = 7_I32      ! [-] Prefix: I, Joint1/2_ID, Prop_I/J, Mass, Length

    CHARACTER(len=token_len)      :: tokens(max_tokens)
    CHARACTER(len=len(text_line)) :: line_trim
    INTEGER(I32)                  :: n_tokens, i, istart, line_len, ierr, n_node_ids, pos_hash
    LOGICAL                       :: in_token


    line_trim = adjustl(text_line)
    line_len  = len_trim(line_trim)

    ! Strip the leading '#' comment marker
    if (line_len > 0_I32 .and. line_trim(1:1) == '#') line_trim(1:1) = ' '

    ! A second '#' marks a trailing annotation (e.g. "# Rigid link"):
    ! truncate the line there instead of blanking it
    pos_hash = index(line_trim(1:line_len), '#')
    if (pos_hash > 0_I32) line_len = pos_hash - 1_I32
    line_len = len_trim(line_trim(1:line_len))

    n_tokens = 0_I32
    in_token = .false.
    istart   = 0_I32
    do i = 1, line_len
        if (line_trim(i:i) /= ' ') then
            if (.not. in_token) then
                in_token = .true.
                istart   = i
            end if
        else
            if (in_token) then
                in_token = .false.
                n_tokens = n_tokens + 1_I32
                if (n_tokens > max_tokens) then
                    write(*,*) 'Error in parse_member_nodes: too many tokens in line (limit ', &
                               max_tokens, '): ', trim(text_line)
                    stop 1
                end if
                tokens(n_tokens) = line_trim(istart:i-1)
            end if
        end if
    end do

    ! Handle a trailing token if the line does not end with a blank space
    if (in_token) then
        n_tokens = n_tokens + 1_I32
        if (n_tokens > max_tokens) then
            write(*,*) 'Error in parse_member_nodes: too many tokens in line (limit ', &
                       max_tokens, '): ', trim(text_line)
            stop 1
        end if
        tokens(n_tokens) = line_trim(istart:line_len)
    end if

    ! Fast fail validation: need at least the fixed prefix columns
    if (n_tokens < n_prefix) then
        write(*,*) 'Error in parse_member_nodes: line has only ', n_tokens, &
                   ' tokens, expected at least ', n_prefix, ' prefix columns. Line: ', trim(text_line)
        stop 1
    end if

    ! Node IDs are whatever tokens remain after the fixed prefix columns
    ! (variable length: e.g. 2 for a rigid-link member, max_nodes for a
    ! fully-discretized one). Missing slots are left as 0
    node_arr   = 0_I32
    n_node_ids = min(n_tokens - n_prefix, max_nodes)
    do i = 1, n_node_ids
        read(tokens(n_prefix+i), *, iostat=ierr) node_arr(i)
        if (ierr /= 0_I32) then
            write(*,*) 'Error in parse_member_nodes: could not parse token "', &
                       trim(tokens(n_prefix+i)), '" as integer. Line: ', trim(text_line)
            stop 1
        end if
    end do

END SUBROUTINE parse_member_nodes


FUNCTION extract_number_nodes(text_line) RESULT(number_of_nodes)
    ! -------------------------------------------------------------------------
    ! Extracts the integer node count from an SD.sum.yaml line of the form
    ! "Nodes: # 33 x 9", i.e. the number found between the '#' and the 'x'.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: text_line          ! [-] Line read in SD.sum.yaml, e.g. "Nodes: # 33 x 9"
    INTEGER(I32)                 :: number_of_nodes    ! [-] Number of nodes extracted from between '#' and 'x'

    ! Local variables
    INTEGER(I32)                  :: pos_hash, pos_x, ierr
    CHARACTER(len=len(text_line)) :: between_str


    ! Locate the '#' and the 'x' that follow it
    pos_hash = index(text_line, '#')
    if (pos_hash == 0_I32) then
        write(*,*) 'ERROR in extract_number_nodes: no "#" found in line: ', trim(text_line)
        stop 1
    end if

    pos_x = index(text_line(pos_hash+1:), 'x')
    if (pos_x == 0_I32) then
        write(*,*) 'ERROR in extract_number_nodes: no "x" found after "#" in line: ', trim(text_line)
        stop 1
    end if
    pos_x = pos_x + pos_hash    ! Convert back to index within the full string

    ! Extract and convert the substring between '#' and 'x'
    between_str = text_line(pos_hash+1:pos_x-1)

    read(between_str, *, iostat=ierr) number_of_nodes
    if (ierr /= 0_I32) then
        write(*,*) 'ERROR in extract_number_nodes: could not parse "', trim(adjustl(between_str)), &
                   '" as integer. Line: ', trim(text_line)
        stop 1
    end if

END FUNCTION extract_number_nodes


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


SUBROUTINE read_hdf5_real_1d(filepath, var_name, data_array, found)
    ! -------------------------------------------------------------------------
    ! Reads a rank-1 REAL(WP) dataset (as written by save_hdf5_real_1d).
    ! found = .false. and data_array unallocated if the file cannot be opened
    ! or the dataset does not exist or is not rank-1.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)               :: filepath         ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN)               :: var_name         ! [-] Dataset name
    REAL(WP)        , INTENT(OUT), ALLOCATABLE :: data_array(:)    ! [-] Dataset values
    LOGICAL         , INTENT(OUT)              :: found            ! [-] Whether the dataset was read

    ! Local variables
    INTEGER(HID_T)   :: file_id, dspace_id, dset_id
    INTEGER          :: error, rank
    LOGICAL          :: file_exists
    INTEGER(HSIZE_T) :: dims(1), maxdims(1)


    found = .false.
    inquire(file=filepath, exist=file_exists)
    if (.not. file_exists) return

    CALL h5open_f(error)
    CALL h5fopen_f(filepath, H5F_ACC_RDONLY_F, file_id, error)
    if (error /= 0) then
        CALL h5close_f(error)
        return
    end if

    CALL h5lexists_f(file_id, var_name, found, error)
    if (found) then
        CALL h5dopen_f(file_id, var_name, dset_id, error)
        CALL h5dget_space_f(dset_id, dspace_id, error)
        CALL h5sget_simple_extent_ndims_f(dspace_id, rank, error)
        found = (rank == 1)
        if (found) then
            CALL h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)
            ALLOCATE(data_array(dims(1)))
            CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_array, dims, error)
        end if
        CALL h5sclose_f(dspace_id, error)
        CALL h5dclose_f(dset_id, error)
    end if

    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE read_hdf5_real_1d


! ===============================================================================
! HDF5 DISPATCHER
! ===============================================================================
SUBROUTINE save_variable(save_path, var_name, var)
    ! -------------------------------------------------------------------------
    ! Generic dispatcher: given a variable of assumed rank/type (CLASS(*)), it
    ! uses SELECT RANK / SELECT TYPE to identify the actual rank (0 to 3) and
    ! type (real/integer/complex, plus logical/character scalars) of "var" and
    ! forwards the call to the matching save_to_hdf5 (save_hdf5_<type>_<rank>)
    ! specific procedure. Any other type or rank stops the program.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: save_path    ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN) :: var_name     ! [-] Dataset name
    CLASS(*)        , INTENT(IN) :: var(..)      ! [-] Variable to save (any supported type, rank 0 to 3)


    SELECT RANK (v_rank => var)
    RANK (0)
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP))
            CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER)
            CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP))
            CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (LOGICAL)
            CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (CHARACTER(len=*))
            CALL save_to_hdf5(save_path, trim(var_name), v)
        CLASS DEFAULT
            error stop "save_variable: unsupported type for dataset " // trim(var_name)
        END SELECT
    RANK (1)
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP))
            CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER)
            CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP))
            CALL save_to_hdf5(save_path, trim(var_name), v)
        CLASS DEFAULT
            error stop "save_variable: unsupported type for dataset " // trim(var_name)
        END SELECT
    RANK (2)
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP))
            CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER)
            CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP))
            CALL save_to_hdf5(save_path, trim(var_name), v)
        CLASS DEFAULT
            error stop "save_variable: unsupported type for dataset " // trim(var_name)
        END SELECT
    RANK (3)
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP))
            CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER)
            CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP))
            CALL save_to_hdf5(save_path, trim(var_name), v)
        CLASS DEFAULT
            error stop "save_variable: unsupported type for dataset " // trim(var_name)
        END SELECT
    RANK DEFAULT
        error stop "save_variable: unsupported rank for dataset " // trim(var_name)
    END SELECT

END SUBROUTINE save_variable


! ===============================================================================
! HDF5 SAVE_TO_HDF5 SPECIFIC PROCEDURES: SCALARS
! ===============================================================================
! Each of the following procedures implements one (type, rank) overload of
! the save_to_hdf5 generic interface: it opens (or creates) the target
! HDF5 file, creates (or overwrites) a dataset named var_name, writes the
! given data into it, and closes the file. They are otherwise identical
! in structure; only the HDF5 datatype and dataspace rank differ.
SUBROUTINE save_hdf5_real_scalar(filepath, var_name, data_scalar)
    ! -------------------------------------------------------------------------
    ! save_to_hdf5 overload for a scalar REAL(WP), written as a 1-element
    ! H5T_NATIVE_DOUBLE dataset.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)         :: filepath       ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN)         :: var_name       ! [-] Dataset name
    REAL(WP)        , INTENT(IN), TARGET :: data_scalar    ! [-] Value to save

    ! Local variables
    INTEGER(HID_T)   :: file_id, dspace_id, dset_id
    INTEGER          :: error
    LOGICAL          :: file_exists, link_exists
    INTEGER(HSIZE_T) :: dims(1)


    CALL h5open_f(error)

    inquire(file=filepath, exist=file_exists)
    if (file_exists) then
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    else
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    end if

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    if (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = 1
    CALL h5screate_simple_f(1, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_scalar, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE save_hdf5_real_scalar


SUBROUTINE save_hdf5_int_scalar(filepath, var_name, data_scalar)
    ! -------------------------------------------------------------------------
    ! save_to_hdf5 overload for a scalar INTEGER, written as a 1-element
    ! H5T_NATIVE_INTEGER dataset.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)         :: filepath       ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN)         :: var_name       ! [-] Dataset name
    INTEGER         , INTENT(IN), TARGET :: data_scalar    ! [-] Value to save

    ! Local variables
    INTEGER(HID_T)   :: file_id, dspace_id, dset_id
    INTEGER          :: error
    LOGICAL          :: file_exists, link_exists
    INTEGER(HSIZE_T) :: dims(1)


    CALL h5open_f(error)

    inquire(file=filepath, exist=file_exists)
    if (file_exists) then
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    else
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    end if

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    if (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = 1
    CALL h5screate_simple_f(1, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_scalar, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE save_hdf5_int_scalar


SUBROUTINE save_hdf5_cplx_scalar(filepath, var_name, data_scalar)
    ! -------------------------------------------------------------------------
    ! save_to_hdf5 overload for a scalar COMPLEX(WP), written as a 2-element
    ! H5T_NATIVE_DOUBLE dataset holding (real part, imaginary part).
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: filepath       ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN) :: var_name       ! [-] Dataset name
    COMPLEX(WP)     , INTENT(IN) :: data_scalar    ! [-] Value to save

    ! Local variables
    INTEGER(HID_T)           :: file_id, dspace_id, dset_id
    INTEGER                  :: error
    LOGICAL                  :: file_exists, link_exists
    INTEGER(HSIZE_T)         :: dims(1)
    REAL(WP)        , TARGET :: tmp(2)


    CALL h5open_f(error)

    inquire(file=filepath, exist=file_exists)
    if (file_exists) then
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    else
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    end if

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    if (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    tmp(1) = real(data_scalar, WP)
    tmp(2) = aimag(data_scalar)

    dims(1) = 2
    CALL h5screate_simple_f(1, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tmp, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE save_hdf5_cplx_scalar


SUBROUTINE save_hdf5_logical_scalar(filepath, var_name, data_scalar)
    ! -------------------------------------------------------------------------
    ! save_to_hdf5 overload for a scalar LOGICAL, written as a 1-element
    ! H5T_NATIVE_INTEGER dataset (1 for .true., 0 for .false.).
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)         :: filepath       ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN)         :: var_name       ! [-] Dataset name
    LOGICAL         , INTENT(IN), TARGET :: data_scalar    ! [-] Value to save

    ! Local variables
    INTEGER(HID_T)   :: file_id, dspace_id, dset_id
    INTEGER          :: error
    LOGICAL          :: file_exists, link_exists
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER          :: int_val    ! Converted value


    int_val = merge(1, 0, data_scalar)

    CALL h5open_f(error)

    inquire(file=filepath, exist=file_exists)
    if (file_exists) then
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    else
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    end if

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    if (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = 1
    CALL h5screate_simple_f(1, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, int_val, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE save_hdf5_logical_scalar


SUBROUTINE save_hdf5_char_scalar(filepath, var_name, data_scalar)
    ! -------------------------------------------------------------------------
    ! save_to_hdf5 overload for a scalar CHARACTER string, written as a
    ! 1-element fixed-length string dataset sized to the trimmed text.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: filepath       ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN) :: var_name       ! [-] Dataset name
    CHARACTER(len=*), INTENT(IN) :: data_scalar    ! [-] String to save

    ! Local variables
    INTEGER(HID_T)   :: file_id, dspace_id, dset_id, dtype_id
    INTEGER          :: error
    LOGICAL          :: file_exists, link_exists
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER(SIZE_T)  :: str_len


    CALL h5open_f(error)

    inquire(file=filepath, exist=file_exists)
    if (file_exists) then
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    else
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    end if

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    if (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    ! Obtain trimmed text length
    str_len = int(len_trim(data_scalar), kind=SIZE_T)
    if (str_len == 0_SIZE_T) str_len = 1_SIZE_T

    ! Create datatype with adequate length
    CALL h5tcopy_f(H5T_FORTRAN_S1, dtype_id, error)
    CALL h5tset_size_f(dtype_id, str_len, error)

    dims(1) = 1
    CALL h5screate_simple_f(1, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, dtype_id, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, dtype_id, data_scalar, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5tclose_f(dtype_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE save_hdf5_char_scalar


! ===============================================================================
! HDF5 SAVE_TO_HDF5 SPECIFIC PROCEDURES: 1D ARRAYS
! ===============================================================================
SUBROUTINE save_hdf5_real_1d(filepath, var_name, data_array)
    ! -------------------------------------------------------------------------
    ! save_to_hdf5 overload for a rank-1 REAL(WP) array, written as a
    ! rank-1 H5T_NATIVE_DOUBLE dataset.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)         :: filepath         ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN)         :: var_name         ! [-] Dataset name
    REAL(WP)        , INTENT(IN), TARGET :: data_array(:)    ! [-] Array to save

    ! Local variables
    INTEGER(HID_T)   :: file_id, dspace_id, dset_id
    INTEGER          :: error
    LOGICAL          :: file_exists, link_exists
    INTEGER(HSIZE_T) :: dims(1)


    CALL h5open_f(error)

    inquire(file=filepath, exist=file_exists)
    if (file_exists) then
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    else
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    end if

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    if (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = size(data_array, 1)
    CALL h5screate_simple_f(1, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_array, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE save_hdf5_real_1d


SUBROUTINE save_hdf5_int_1d(filepath, var_name, data_array)
    ! -------------------------------------------------------------------------
    ! save_to_hdf5 overload for a rank-1 INTEGER array, written as a
    ! rank-1 H5T_NATIVE_INTEGER dataset.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)         :: filepath         ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN)         :: var_name         ! [-] Dataset name
    INTEGER         , INTENT(IN), TARGET :: data_array(:)    ! [-] Array to save

    ! Local variables
    INTEGER(HID_T)   :: file_id, dspace_id, dset_id
    INTEGER          :: error
    LOGICAL          :: file_exists, link_exists
    INTEGER(HSIZE_T) :: dims(1)


    CALL h5open_f(error)

    inquire(file=filepath, exist=file_exists)
    if (file_exists) then
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    else
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    end if

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    if (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = size(data_array, 1)
    CALL h5screate_simple_f(1, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_array, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE save_hdf5_int_1d


SUBROUTINE save_hdf5_cplx_1d(filepath, var_name, data_array)
    ! -------------------------------------------------------------------------
    ! save_to_hdf5 overload for a rank-1 COMPLEX(WP) array, written as a
    ! (2, N) H5T_NATIVE_DOUBLE dataset holding (real part, imaginary part),
    ! i.e. (N, 2) when read with h5py.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: filepath         ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN) :: var_name         ! [-] Dataset name
    COMPLEX(WP)     , INTENT(IN) :: data_array(:)    ! [-] Array to save

    ! Local variables
    INTEGER(HID_T)                        :: file_id, dspace_id, dset_id
    INTEGER                               :: error
    LOGICAL                               :: file_exists, link_exists
    INTEGER(HSIZE_T)                      :: dims(2)
    REAL(WP)        , ALLOCATABLE, TARGET :: tmp(:,:)


    CALL h5open_f(error)

    inquire(file=filepath, exist=file_exists)
    if (file_exists) then
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    else
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    end if

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    if (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = 2
    dims(2) = size(data_array, 1)

    ALLOCATE(tmp(2, dims(2)))
    tmp(1, :) = real(data_array, WP)
    tmp(2, :) = aimag(data_array)

    CALL h5screate_simple_f(2, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tmp, dims, error)

    DEALLOCATE(tmp)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE save_hdf5_cplx_1d


! ===============================================================================
! HDF5 SAVE_TO_HDF5 SPECIFIC PROCEDURES: 2D ARRAYS
! ===============================================================================
SUBROUTINE save_hdf5_real_2d(filepath, var_name, data_array)
    ! -------------------------------------------------------------------------
    ! save_to_hdf5 overload for a rank-2 REAL(WP) array, written as a
    ! rank-2 H5T_NATIVE_DOUBLE dataset.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)         :: filepath           ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN)         :: var_name           ! [-] Dataset name
    REAL(WP)        , INTENT(IN), TARGET :: data_array(:,:)    ! [-] Array to save

    ! Local variables
    INTEGER(HID_T)   :: file_id, dspace_id, dset_id
    INTEGER          :: error
    LOGICAL          :: file_exists, link_exists
    INTEGER(HSIZE_T) :: dims(2)


    CALL h5open_f(error)

    inquire(file=filepath, exist=file_exists)
    if (file_exists) then
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    else
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    end if

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    if (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = size(data_array, 1)
    dims(2) = size(data_array, 2)
    CALL h5screate_simple_f(2, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_array, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE save_hdf5_real_2d


SUBROUTINE save_hdf5_int_2d(filepath, var_name, data_array)
    ! -------------------------------------------------------------------------
    ! save_to_hdf5 overload for a rank-2 INTEGER array, written as a
    ! rank-2 H5T_NATIVE_INTEGER dataset.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)         :: filepath           ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN)         :: var_name           ! [-] Dataset name
    INTEGER         , INTENT(IN), TARGET :: data_array(:,:)    ! [-] Array to save

    ! Local variables
    INTEGER(HID_T)   :: file_id, dspace_id, dset_id
    INTEGER          :: error
    LOGICAL          :: file_exists, link_exists
    INTEGER(HSIZE_T) :: dims(2)


    CALL h5open_f(error)

    inquire(file=filepath, exist=file_exists)
    if (file_exists) then
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    else
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    end if

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    if (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = size(data_array, 1)
    dims(2) = size(data_array, 2)
    CALL h5screate_simple_f(2, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_array, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE save_hdf5_int_2d


SUBROUTINE save_hdf5_cplx_2d(filepath, var_name, data_array)
    ! -------------------------------------------------------------------------
    ! save_to_hdf5 overload for a rank-2 COMPLEX(WP) array, written as a
    ! (2, N1, N2) H5T_NATIVE_DOUBLE dataset holding (real, imaginary) parts,
    ! i.e. (N2, N1, 2) when read with h5py.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: filepath           ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN) :: var_name           ! [-] Dataset name
    COMPLEX(WP)     , INTENT(IN) :: data_array(:,:)    ! [-] Array to save

    ! Local variables
    INTEGER(HID_T)                        :: file_id, dspace_id, dset_id
    INTEGER                               :: error
    LOGICAL                               :: file_exists, link_exists
    INTEGER(HSIZE_T)                      :: dims(3)
    REAL(WP)        , ALLOCATABLE, TARGET :: tmp(:,:,:)


    CALL h5open_f(error)

    inquire(file=filepath, exist=file_exists)
    if (file_exists) then
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    else
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    end if

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    if (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = 2
    dims(2) = size(data_array, 1)
    dims(3) = size(data_array, 2)

    ALLOCATE(tmp(2, dims(2), dims(3)))
    tmp(1, :, :) = real(data_array, WP)
    tmp(2, :, :) = aimag(data_array)

    CALL h5screate_simple_f(3, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tmp, dims, error)

    DEALLOCATE(tmp)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE save_hdf5_cplx_2d


! ===============================================================================
! HDF5 SAVE_TO_HDF5 SPECIFIC PROCEDURES: 3D ARRAYS
! ===============================================================================
SUBROUTINE save_hdf5_real_3d(filepath, var_name, data_array)
    ! -------------------------------------------------------------------------
    ! save_to_hdf5 overload for a rank-3 REAL(WP) array, written as a
    ! rank-3 H5T_NATIVE_DOUBLE dataset.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)         :: filepath             ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN)         :: var_name             ! [-] Dataset name
    REAL(WP)        , INTENT(IN), TARGET :: data_array(:,:,:)    ! [-] Array to save

    ! Local variables
    INTEGER(HID_T)   :: file_id, dspace_id, dset_id
    INTEGER          :: error
    LOGICAL          :: file_exists, link_exists
    INTEGER(HSIZE_T) :: dims(3)


    CALL h5open_f(error)

    inquire(file=filepath, exist=file_exists)
    if (file_exists) then
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    else
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    end if

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    if (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = size(data_array, 1)
    dims(2) = size(data_array, 2)
    dims(3) = size(data_array, 3)
    CALL h5screate_simple_f(3, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_array, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE save_hdf5_real_3d


SUBROUTINE save_hdf5_int_3d(filepath, var_name, data_array)
    ! -------------------------------------------------------------------------
    ! save_to_hdf5 overload for a rank-3 INTEGER array, written as a
    ! rank-3 H5T_NATIVE_INTEGER dataset.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)         :: filepath             ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN)         :: var_name             ! [-] Dataset name
    INTEGER         , INTENT(IN), TARGET :: data_array(:,:,:)    ! [-] Array to save

    ! Local variables
    INTEGER(HID_T)   :: file_id, dspace_id, dset_id
    INTEGER          :: error
    LOGICAL          :: file_exists, link_exists
    INTEGER(HSIZE_T) :: dims(3)


    CALL h5open_f(error)

    inquire(file=filepath, exist=file_exists)
    if (file_exists) then
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    else
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    end if

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    if (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = size(data_array, 1)
    dims(2) = size(data_array, 2)
    dims(3) = size(data_array, 3)
    CALL h5screate_simple_f(3, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_array, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE save_hdf5_int_3d


SUBROUTINE save_hdf5_cplx_3d(filepath, var_name, data_array)
    ! -------------------------------------------------------------------------
    ! save_to_hdf5 overload for a rank-3 COMPLEX(WP) array, written as a
    ! (2, N1, N2, N3) H5T_NATIVE_DOUBLE dataset holding (real, imaginary) parts,
    ! i.e. (N3, N2, N1, 2) when read with h5py.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: filepath             ! [-] Path to the HDF5 file
    CHARACTER(len=*), INTENT(IN) :: var_name             ! [-] Dataset name
    COMPLEX(WP)     , INTENT(IN) :: data_array(:,:,:)    ! [-] Array to save

    ! Local variables
    INTEGER(HID_T)                        :: file_id, dspace_id, dset_id
    INTEGER                               :: error
    LOGICAL                               :: file_exists, link_exists
    INTEGER(HSIZE_T)                      :: dims(4)
    REAL(WP)        , ALLOCATABLE, TARGET :: tmp(:,:,:,:)


    CALL h5open_f(error)

    inquire(file=filepath, exist=file_exists)
    if (file_exists) then
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    else
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    end if

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    if (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = 2
    dims(2) = size(data_array, 1)
    dims(3) = size(data_array, 2)
    dims(4) = size(data_array, 3)

    ALLOCATE(tmp(2, dims(2), dims(3), dims(4)))
    tmp(1, :, :, :) = real(data_array, WP)
    tmp(2, :, :, :) = aimag(data_array)

    CALL h5screate_simple_f(4, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tmp, dims, error)

    DEALLOCATE(tmp)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

END SUBROUTINE save_hdf5_cplx_3d

END MODULE IOUtils
