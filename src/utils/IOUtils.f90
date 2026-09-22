! ===============================================================================
! MODULE: IOUtils
! ===============================================================================
! Handles all file I/O for the codebase: parsing OpenFAST SubDyn summary files
! (node/member topology), parsing OpenFAST binary/ASCII output channel files,
! a generic HDF5 variable-saving facility (via a save_variable generic
! interface overloaded for scalar/1D/2D/3D real/int/complex/logical/character
! data), small string helpers (to_upper, swap_extension, ensure_parent_dir),
! and a curve-lookup helper (read_curve).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 2026-09-22
! ===============================================================================
MODULE IOUtils

USE HDF5
USE KINDS, ONLY: I32, WP
USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN

IMPLICIT NONE

PRIVATE
PUBLIC :: get_SDsum_variables, read_input_SD, read_curve, save_variable, save_to_hdf5, to_upper, swap_extension, &
          ensure_parent_dir

! -------------------------------------------------------------------------
! Generic interface: save_to_hdf5
! Dispatches to the correct save_hdf5_<type>_<rank> specific procedure based
! on the actual type/rank of the variable passed in. This is the low-level
! layer used internally by the save_variable dispatcher below.
! -------------------------------------------------------------------------
INTERFACE save_to_hdf5
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
END INTERFACE

CONTAINS


! -------------------------------------------------------------------------
! Subroutine: get_SDsum_variables
! Parses an OpenFAST SubDyn summary file (SD.sum.yaml) to recover the
! structural topology: it locates the "#Member I Joint1_ID Joint2_ID"
! table to map each member to its node IDs, then locates the node
! coordinate table ("Node_[#] X_[m] Y_[m] Z_[m]") to recover each node's
! (x,y,z) position, and finally combines both into a (Nmembers, Nnodes, 3)
! tensor, padding any missing node with NaN.
! -------------------------------------------------------------------------
SUBROUTINE get_SDsum_variables(SD_path, Nmembers, Nnodes, verbose, Nodes)
    CHARACTER(len=*), INTENT(IN)           :: SD_path       ! [-] Path to OpenFAST SubDyn sum file
    INTEGER(I32)    , INTENT(IN), OPTIONAL :: Nmembers      ! [-] Number of OpenFAST members
    INTEGER(I32)    , INTENT(IN), OPTIONAL :: Nnodes        ! [-] Number of OpenFAST nodes
    LOGICAL         , INTENT(IN), OPTIONAL :: verbose       ! [-] Flag to print more info

    REAL(WP), ALLOCATABLE, INTENT(OUT) :: Nodes(:,:,:)      ! [m] Nodes array (Nmembers, Nnodes, 3)

    ! Local variables
    INTEGER(I32) :: Nmembers_ = 8, Nnodes_ = 5
    LOGICAL      :: verbose_ = .false.

    INTEGER(I32), ALLOCATABLE :: member_nodes(:,:)
    REAL(WP)    , ALLOCATABLE :: Nodes_flat(:,:)
    INTEGER(I32)              :: status, file_unit, ios, i, j, k, idx, num_nodes
    INTEGER(I32)              :: pos_open, pos_close
    REAL(WP)                  :: nan_val, x, y, z, idx_real
    LOGICAL                   :: found
    CHARACTER(len=1024)       :: line_buf
    CHARACTER(len=1024)       :: data_str


    ! Defaults
    if (present(Nmembers)) Nmembers_ = Nmembers
    if (present(Nnodes)) Nnodes_ = Nnodes
    if (present(verbose)) verbose_ = verbose

    ! Initialize IEEE Nan for padding missing nodes
    nan_val = ieee_value(1.0_wp, ieee_quiet_nan)
    status = 0

    if (verbose_) print *, "Reading SD.sum nodes in ", trim(SD_path)

    ALLOCATE(Nodes(Nmembers_, Nnodes_,3), source=0.0_WP)
    ALLOCATE(member_nodes(Nmembers_, Nnodes_), source=0_I32)

    open(newunit=file_unit, file= trim(SD_path), status='old', action='read', iostat=status)
    if (status /= 0) then
        print *, "Error: Could not open ", trim(SD_path)
        return
    end if

    ! ----------------------------------------------
    ! 1. Extract Connection Map: Member -> Node IDs
    ! ----------------------------------------------
    found = .false.
    do while (.true.)
        read(file_unit, '(A)', iostat=ios) line_buf
        if (ios /= 0) exit
        if(index(line_buf, "#Member I Joint1_ID Joint2_ID") > 0) then
            found = .true.
            exit
        end if
    end do

    if (.not. found) then
        print *, "Error: Member mapping header not found."
        status = -1; close(file_unit); return
    end if

    ! Parse node IDs for each member
    do i = 1, Nmembers_
        read(file_unit, '(A)') line_buf
        CALL parse_member_nodes(line_buf, Nnodes_, member_nodes(i,:))
    end do

    ! ----------------------------------------------
    ! 2. Extract Global Node Table: ID -> (x, y, z)
    ! ----------------------------------------------
    rewind(file_unit)
    found = .false.
    do while(.true.)
        read(file_unit, '(A)', iostat=ios) line_buf
        if (ios /= 0) exit
        if (index(line_buf, '#     Node_[#]          X_[m]           Y_[m]           Z_[m]') > 0) then
            found = .true.
            exit
        end if
    end do

    if (.not. found) then
        print *, "Error: Node coordinate header not found."
        status = -1; close(file_unit); return
    end if

    ! read the number of nodes (e.g., "Nodes: # 33 x 9")
    read(file_unit, '(A)') line_buf
    num_nodes = extract_number_nodes(line_buf)

    ALLOCATE(Nodes_flat(num_nodes,3), source=0.0_WP)

    ! read coordinates: each line looks like
    !   - [      1.,          0.000,          0.000,        -30.000, ... ]
    ! Extract only the content between '[' and ']' before the numeric READ,
    ! so the leading '-' YAML list marker never gets fed to the parser
    ! (removing it globally would corrupt negative coordinate values).
    do i = 1, num_nodes
        read(file_unit, '(A)') line_buf

        pos_open  = INDEX(line_buf, '[')
        pos_close = INDEX(line_buf, ']')
        if (pos_open == 0 .or. pos_close == 0 .or. pos_close <= pos_open) then
            print *, "Error: malformed node coordinate line: ", trim(line_buf)
            status = -1; close(file_unit); return
        end if

        data_str = line_buf(pos_open+1:pos_close-1)

        ! Clear commas so list-directed READ treats fields as blank-separated
        do k = 1, len_trim(data_str)
            if (data_str(k:k) == ',') data_str(k:k) = ' '
        end do

        ! Node ID is stored as a real literal (e.g. "1."), read as real then
        ! convert, mirroring the Python int(float(...)) pattern
        read(data_str, *, iostat=ios) idx_real, x, y, z
        if (ios /= 0) then
            print *, "Error: could not parse node coordinate line: ", trim(line_buf)
            status = -1; close(file_unit); return
        end if

        idx = NINT(idx_real)
        if (idx >= 1 .and. idx <= num_nodes) then
            Nodes_flat(idx, :) = [x, y, z]
        else
            print *, "Error: node index ", idx, " out of range [1, ", num_nodes, "]"
            status = -1; close(file_unit); return
        end if
    end do

    close(file_unit)

    !----------------------------------------------------------
    ! 3. Map Flat Nodes to Output Tensor (Nmembers, Nnodes, 3)
    ! ---------------------------------------------------------
    do i = 1, Nmembers_
        do j = 1, Nnodes_
            idx = member_nodes(i,j)
            if (idx > 0) then
                Nodes(i,j,:) = Nodes_flat(idx, :)
            else
                Nodes(i,j,:) = nan_val
            end if
        end do
    end do


END SUBROUTINE get_SDsum_variables


! -------------------------------------------------------------------------
! Subroutine: read_input_SD
! Reads a time series of a requested structural quantity ("displacement",
! "acceleration", "force", "momentum", etc.) for every member/node from an
! OpenFAST output file. It builds the corresponding "MiNjXXX"-style channel
! names for all (member, node, xyz) combinations, delegates the actual file
! reading to parse_channels_auto, then decimates (skip) and reshapes the
! flat channel array into a (nt, Nmembers, Nnodes, 3) tensor.
! -------------------------------------------------------------------------
SUBROUTINE read_input_SD(filename, what, skip, Nmembers, Nnodes, From, Upto, verbose, &
                         Time_out, Array_out, unit_out, status)
    CHARACTER(len=*), INTENT(IN)           :: filename          ! [-] Path to OpenFAST output file
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: what              ! [-] Name of the output to read
    INTEGER(I32)    , INTENT(IN), OPTIONAL :: skip              ! [-] Decimation factor
    INTEGER(I32)    , INTENT(IN), OPTIONAL :: Nmembers          ! [-] Total structural members
    INTEGER(I32)    , INTENT(IN), OPTIONAL :: Nnodes            ! [-] Nodes per member
    REAL(WP)        , INTENT(IN), OPTIONAL :: From              ! [-] Start window fraction
    REAL(WP)        , INTENT(IN), OPTIONAL :: Upto              ! [-] End window fraction
    LOGICAL         , INTENT(IN), OPTIONAL :: verbose           ! [-] Flag to print info

    REAL(WP), ALLOCATABLE, INTENT(OUT) :: Time_out(:)           ! [s] Downsampled time array (nt)
    REAL(WP), ALLOCATABLE, INTENT(OUT) :: Array_out(:,:,:,:)    ! [unit] Data tensor (nt, Nmembers, Nnodes, 3)
    CHARACTER(len=*)     , INTENT(OUT) :: unit_out              ! [-] Physical unit
    INTEGER(I32)         , INTENT(OUT) :: status

    ! Local Defaults Variables
    CHARACTER(len=32) :: what_ = "acceleration"
    INTEGER(I32)      :: skip_ = 1, Nmembers_ = 8, Nnodes_ = 5
    REAL(WP)          :: From_ = 0.0_WP, Upto_ = 1.0_WP
    LOGICAL           :: verbose_ = .false.

    ! Internal Variables
    CHARACTER(len=16) :: strx, stry, strz
    CHARACTER(len=32), ALLOCATABLE :: out_channels(:), raw_units(:)
    REAL(WP), ALLOCATABLE :: raw_time(:), raw_array(:,:)
    INTEGER(I32) :: m, n, ch_idx, nt, nt_skip, t_skip, t_orig, n_chans, base, idx
    CHARACTER(len=32) :: ch_name

    status = 0

    ! Apply optional arguments
    if (present(what))     what_     = trim(what)
    if (present(skip))     skip_     = skip
    if (present(Nmembers)) Nmembers_ = Nmembers
    if (present(Nnodes))   Nnodes_   = Nnodes
    if (present(From))     From_     = From
    if (present(Upto))     Upto_     = Upto
    if (present(verbose))  verbose_  = verbose

    ! -------------------------------------------------------------
    ! 1. Determine Channel Suffixes Based on 'what'
    ! -------------------------------------------------------------
    select case (trim(what_))
        case ("displacement")
            strx = "TDxss"; stry = "TDyss"; strz = "TDzss"
        case ("acceleration")
            strx = "TAxe"; stry = "TAye"; strz = "TAze"
        case ("force")
            strx = "FKxe"; stry = "FKye"; strz = "FKze"
        case ("momentum")
            strx = "MKxe"; stry = "MKye"; strz = "MKze"
        case ("int_displacement")
            strx = "IntfTDXss"; stry = "IntfTDYss"; strz = "IntfTDZss"
        case ("int_acceleration")
            strx = "IntfTAXss"; stry = "IntfTAYss"; strz = "IntfTAZss"
        case default
            print *, "Error: Invalid option for 'what' = ", trim(what_)
            print *, "Input 'displacement', 'acceleration', 'force', 'momentum', 'int_displacement' or 'int_acceleration'"
            status = -1; return
    end select

    if (verbose_) then
        print *, "Reading ", trim(filename), " --> ", trim(what_), &
                 ": [MiNj", trim(strx), ", MiNj", trim(stry), ", MiNj", trim(strz), "]"
    end if

    ! -------------------------------------------------------------
    ! 2. Dynamically Generate Channel Names
    ! -------------------------------------------------------------
    n_chans = Nmembers_ * Nnodes_ * 3
    allocate(out_channels(n_chans))
    ch_idx = 1

    do m = 1, Nmembers_
        do n = 1, Nnodes_
            ! X component
            write(ch_name, '("M",I0,"N",I0,A)') m, n, trim(strx)
            out_channels(ch_idx) = trim(ch_name); ch_idx = ch_idx + 1
            ! Y component
            write(ch_name, '("M",I0,"N",I0,A)') m, n, trim(stry)
            out_channels(ch_idx) = trim(ch_name); ch_idx = ch_idx + 1
            ! Z component
            write(ch_name, '("M",I0,"N",I0,A)') m, n, trim(strz)
            out_channels(ch_idx) = trim(ch_name); ch_idx = ch_idx + 1
        end do
    end do

    ! -------------------------------------------------------------
    ! 3. Call the Universal OpenFAST Parser
    ! -------------------------------------------------------------
    CALL parse_channels_auto(full_path=filename, plotChannels=out_channels, &
                             From=From_, Upto=Upto_,              &
                             available_channels=.false., verbose=verbose_,     &
                             time_out=raw_time, array_out=raw_array,        &
                             units_out=raw_units, status=status)
    if (status /= 0) return

    nt = size(raw_time)
    if (nt == 0) then
        print *, "Warning: No data returned from parser."
        return
    end if

    ! -------------------------------------------------------------
    ! 4. Decimate (Skip) and Reshape into 4D Tensor
    ! -------------------------------------------------------------
    nt_skip = (nt - 1) / skip_ + 1
    allocate(Time_out(nt_skip))
    allocate(Array_out(nt_skip, Nmembers_, Nnodes_, 3))

    t_skip = 1
    do t_orig = 1, nt, skip_
        Time_out(t_skip) = raw_time(t_orig)
        do m = 1, Nmembers_
            base = (m - 1) * Nnodes_ * 3 + 1
            !$OMP SIMD
            do n = 1, Nnodes_
                idx = base + (n - 1) * 3
                Array_out(t_skip, m, n, 1) = raw_array(t_orig, idx)
                Array_out(t_skip, m, n, 2) = raw_array(t_orig, idx + 1)
                Array_out(t_skip, m, n, 3) = raw_array(t_orig, idx + 2)
            end do
        end do
        t_skip = t_skip + 1
    end do

    ! Extract unit from the first valid channel
    unit_out = raw_units(1)

END SUBROUTINE read_input_SD


! -------------------------------------------------------------------------
! Subroutine: parse_channels_auto
! Universal entry point for reading OpenFAST output channels: it inspects
! the file extension (".outb" vs ".out") and dispatches to the matching
! specialized parser, parse_channels_binary or parse_channels_ascii, so
! callers do not need to know which format a given run produced.
! -------------------------------------------------------------------------
SUBROUTINE parse_channels_auto(full_path, plotChannels, From, Upto, &
                               available_channels, verbose, chunk_size,     &
                               time_out, array_out, units_out, status)
    CHARACTER(len=*), INTENT(IN)           :: full_path
    CHARACTER(len=*), INTENT(IN)           :: plotChannels(:)
    REAL(WP)        , INTENT(IN), OPTIONAL :: From, Upto
    LOGICAL         , INTENT(IN), OPTIONAL :: available_channels, verbose
    INTEGER(I32)    , INTENT(IN), OPTIONAL :: chunk_size

    REAL(WP)         , ALLOCATABLE, INTENT(OUT) :: time_out(:)
    REAL(WP)         , ALLOCATABLE, INTENT(OUT) :: array_out(:,:)
    CHARACTER(len=32), ALLOCATABLE, INTENT(OUT) :: units_out(:)
    INTEGER(I32)                  , INTENT(OUT) :: status

    ! Local Variables
    REAL(WP)     :: From_ = 0.0_WP, Upto_ = 1.0_WP
    LOGICAL      :: avail_ch = .false., verb = .true.
    INTEGER(I32) :: c_size = 4096_I32, ext_idx

    ! Defaults
    if (present(From))               From_    = From
    if (present(Upto))               Upto_    = Upto
    if (present(available_channels)) avail_ch = available_channels
    if (present(verbose))            verb     = verbose
    if (present(chunk_size))        c_size    = chunk_size
    status = 0

    ! Detect Extension
    ext_idx = index(full_path, '.outb', back=.true.)
    if (ext_idx > 0 .and. ext_idx == len_trim(full_path) - 4) then
        if (verb) print *, "Detected binary OpenFAST output (.outb): ", trim(full_path)
        CALL parse_channels_binary(full_path, plotChannels, From_, Upto_, &
                                   avail_ch, verb, c_size, time_out, array_out, units_out, status)
    else
        ext_idx = index(full_path, '.out', back=.true.)
        if (ext_idx > 0 .and. ext_idx == len_trim(full_path) - 3) then
            if (verb) print *, "Detected text OpenFAST output (.out): ", trim(full_path)
            CALL PARSE_CHANNELS_ASCII(full_path, plotChannels, From_, Upto_, &
                                      avail_ch, verb, time_out, array_out, units_out, status)
        else
            print *, "Error: Unsupported file extension. Expected .out or .outb."
            status = -1
        end if
    end if

END SUBROUTINE parse_channels_auto


! -------------------------------------------------------------------------
! Subroutine: parse_channels_binary
! Parses OpenFAST's packed binary output format (.outb): reads the fixed
! header (FileID, channel count/names/units, time and per-column scale
! and offset factors used to decompress the packed INTEGER(2) values),
! maps the requested channel names to column indices, then streams the
! data body in chunks (to bound memory use), decoding only the requested
! time window and channels back into REAL(WP) using the scale/offset
! factors.
! -------------------------------------------------------------------------
SUBROUTINE parse_channels_binary(full_path, plotChannels, From_val, Upto_val, &
                                 available_channels, verbose, chunk_size,     &
                                 time_out, array_out, units_out, status)
    CHARACTER(len=*), INTENT(IN)               :: full_path
    CHARACTER(len=*), INTENT(IN)               :: plotChannels(:)
    REAL(WP), INTENT(IN)                       :: From_val, Upto_val
    LOGICAL, INTENT(IN)                        :: available_channels, verbose
    INTEGER(I32), INTENT(IN)                   :: chunk_size

    REAL(WP), ALLOCATABLE, INTENT(OUT)         :: time_out(:)
    REAL(WP), ALLOCATABLE, INTENT(OUT)         :: array_out(:,:)
    CHARACTER(len=32), ALLOCATABLE, INTENT(OUT):: units_out(:)
    INTEGER(I32), INTENT(OUT)                  :: status

    ! OpenFAST Binary Format Constants
    INTEGER(2), PARAMETER :: FMT_WITH_TIME    = 1
    INTEGER(2), PARAMETER :: FMT_WITHOUT_TIME = 2
    INTEGER(2), PARAMETER :: FMT_NO_COMPRESS  = 3
    INTEGER(2), PARAMETER :: FMT_CHAN_LEN     = 4

    ! File Header Variables
    INTEGER(I32) :: fu, i, j, k, req_chans
    INTEGER(2)   :: FileID, LenName
    INTEGER(I32) :: NumOutChans, NT
    REAL(8)      :: TimeScl, TimeOff, TimeOut1, TimeIncr
    REAL(4), ALLOCATABLE :: ColScl(:), ColOff(:)
    INTEGER(I32) :: LenDesc
    CHARACTER(len=:), ALLOCATABLE  :: DescStr
    CHARACTER(len=32), ALLOCATABLE :: ChanName(:), ChanUnit(:)
    INTEGER(I32), ALLOCATABLE      :: ReqIdx(:)

    ! Array Streaming Variables
    INTEGER(I32) :: time_bytes, N_cols_raw
    INTEGER(I32) :: idx_0_start, idx_0_end, start_idx, end_idx, n_target_rows
    INTEGER(I32) :: row_counter, written, chunk_rows, rows_to_read, current_row
    INTEGER(4)   :: p_time
    INTEGER(2)   :: val_i16
    INTEGER(2), ALLOCATABLE :: raw_data(:,:)

    status = 0
    req_chans = size(plotChannels)

    open(newunit=fu, file=trim(full_path), access='stream', form='unformatted', status='old', iostat=status)
    if (status /= 0) then
        print *, "Error: Could not open binary file: ", trim(full_path)
        return
    end if

    ! --- 1. Read Header ---
    read(fu) FileID
    
    if (FileID == FMT_CHAN_LEN) then
        read(fu) LenName
    else
        LenName = 10_2
    end if

    read(fu) NumOutChans
    read(fu) NT

    ! --- 2. Time Info ---
    ! AQUÍ ESTABA EL BUG: FileID = 4 NO empaqueta el tiempo en la matriz
    if (FileID == FMT_WITH_TIME) then
        read(fu) TimeScl
        read(fu) TimeOff
        time_bytes = 4
    else
        read(fu) TimeOut1
        read(fu) TimeIncr
        time_bytes = 0
    end if

    if (FileID == FMT_NO_COMPRESS) then
        print *, "Error: Uncompressed float64 binaries are not supported."
        status = -1; close(fu); return
    end if

    ! --- 3. Channel Scaling ---
    allocate(ColScl(NumOutChans), ColOff(NumOutChans))
    read(fu) ColScl
    read(fu) ColOff

    ! --- 4. Description String ---
    read(fu) LenDesc
    if (LenDesc > 0) then
        allocate(CHARACTER(len=LenDesc) :: DescStr)
        read(fu) DescStr
        deallocate(DescStr)
    end if

    ! --- 5. Channel Names & Units ---
    allocate(ChanName(NumOutChans + 1))
    allocate(ChanUnit(NumOutChans + 1))
    allocate(ReqIdx(req_chans))
    
    ChanName = ""; ChanUnit = ""
    do i = 1, NumOutChans + 1; read(fu) ChanName(i)(1:LenName); end do
    do i = 1, NumOutChans + 1; read(fu) ChanUnit(i)(1:LenName); end do

    if (available_channels) then
        print *, "Available Channels in .outb: "
        do i = 1, NumOutChans + 1; print '(2A)', " - ", trim(ChanName(i)); end do
    end if

    ! --- 6. Map Requested Channels ---
    allocate(units_out(req_chans))
    do i = 1, req_chans
        ReqIdx(i) = 0
        do j = 1, NumOutChans + 1
            if (trim(adjustl(plotChannels(i))) == trim(adjustl(ChanName(j)))) then
                if (j == 1) then
                    ReqIdx(i) = -1 ! Marcador para columna de tiempo
                else
                    ReqIdx(i) = j - 1 ! Offset canales de datos (1 a NumOutChans)
                end if
                units_out(i) = ChanUnit(j)
                exit
            end if
        end do
        if (ReqIdx(i) == 0) then
            print *, "Error: Channel '", trim(plotChannels(i)), "' not found in .outb."
            status = -1; close(fu); return
        end if
    end do

    ! --- 7. Window Filtering Logic ---
    idx_0_start = int(real(NT, WP) * From_val)
    idx_0_end   = int(real(NT, WP) * Upto_val)
    
    start_idx = idx_0_start + 1
    end_idx   = idx_0_end
    n_target_rows = max(0, end_idx - start_idx + 1)

    allocate(time_out(n_target_rows))
    allocate(array_out(n_target_rows, req_chans))

    ! --- 8. Read Data as Massive Contiguous Arrays ---
    if (time_bytes == 4) then
        N_cols_raw = NumOutChans + 2
    else
        N_cols_raw = NumOutChans
    end if

    chunk_rows = chunk_size
    allocate(raw_data(N_cols_raw, chunk_rows))

    row_counter = 0
    written = 0

    do while (row_counter < NT)
        rows_to_read = min(chunk_rows, NT - row_counter)
        read(fu, iostat=status) raw_data(:, 1:rows_to_read)
        if (status /= 0) exit

        do i = 1, rows_to_read
            current_row = row_counter + i
            if (current_row >= start_idx .and. current_row <= end_idx) then
                written = written + 1

                ! Decode time
                if (time_bytes == 4) then
                    p_time = transfer(raw_data(1:2, i), p_time)
                    time_out(written) = real((real(p_time, 8) - TimeOff) / TimeScl, WP)
                else
                    time_out(written) = real(TimeOut1 + TimeIncr * real(current_row - 1, 8), WP)
                end if

                ! Assign time column if any (assume at most one)
                do j = 1, req_chans
                    if (ReqIdx(j) == -1) then
                        array_out(written, j) = time_out(written)
                        exit
                    end if
                end do

                ! Decode data channels with SIMD
                !$OMP SIMD
                do j = 1, req_chans
                    if (ReqIdx(j) /= -1) then
                        k = ReqIdx(j)
                        if (time_bytes == 4) then
                            val_i16 = raw_data(k + 2, i)
                        else
                            val_i16 = raw_data(k, i)
                        end if
                        array_out(written, j) = real((real(val_i16, 4) - ColOff(k)) / ColScl(k), WP)
                    end if
                end do
            end if
        end do

        row_counter = row_counter + rows_to_read
        if (written == n_target_rows) exit
    end do

    if (verbose) then
        print *, "Output array shape: [", written, ", ", req_chans, "]"
    end if


END SUBROUTINE parse_channels_binary


! -------------------------------------------------------------------------
! Subroutine: parse_channels_ascii
! Parses OpenFAST's plain-text output format (.out): skips the 6 metadata
! header lines, reads the channel-name and channel-unit header lines
! (line 7/8) to map requested channel names to column positions and to
! locate the "Time" column, does a fast pre-pass to count total data rows
! (needed to size the From/Upto time window), then rewinds and reads only
! the requested rows/columns into the output arrays.
! -------------------------------------------------------------------------
SUBROUTINE parse_channels_ascii(full_path, plotChannels, From_val, Upto_val, &
                                available_channels, verbose, time_out,       &
                                array_out, units_out, status)
    CHARACTER(len=*), INTENT(IN)               :: full_path
    CHARACTER(len=*), INTENT(IN)               :: plotChannels(:)
    REAL(WP), INTENT(IN)                       :: From_val, Upto_val
    LOGICAL, INTENT(IN)                        :: available_channels, verbose
    REAL(WP), ALLOCATABLE, INTENT(OUT)         :: time_out(:)
    REAL(WP), ALLOCATABLE, INTENT(OUT)         :: array_out(:,:)
    CHARACTER(len=32), ALLOCATABLE, INTENT(OUT):: units_out(:)
    INTEGER(I32), INTENT(OUT)                  :: status

    ! Local variables
    INTEGER(I32) :: fu, i, j, k, ios, n_req, n_all_chans, n_total, time_col
    INTEGER(I32) :: idx_0_start, idx_0_end, start_idx, end_idx, n_filtered, out_row
    INTEGER(I32) :: row_idx, pos, l_len
    CHARACTER(len=4096) :: line_buffer
    CHARACTER(len=1)    :: dummy_char
    CHARACTER(len=32), ALLOCATABLE :: all_channels(:), all_units(:)
    INTEGER(I32), ALLOCATABLE :: col_map(:)
    REAL(WP), ALLOCATABLE :: row_buff(:)

    status = 0
    n_req = size(plotChannels)

    if (n_req == 0) then
        print *, "Error: plotChannels list cannot be empty."
        status = -1; return
    end if

    open(newunit=fu, file=trim(full_path), status='old', action='read', iostat=status)
    if (status /= 0) then
        print *, "Error: Could not open ASCII file: ", trim(full_path)
        return
    end if

    ! --- 1. Read Header Lines ---
    ! Skip lines 1 to 6 (Description / Metadata)
    do i = 1, 6
        read(fu, '(A)', iostat=status) line_buffer
        if (status /= 0) then
            print *, "Error reading header lines in file: ", trim(full_path)
            close(fu); status = -1; return
        end if
    end do

    ! Line 7: Channel Names
    read(fu, '(A)', iostat=status) line_buffer
    if (status /= 0) then; close(fu); status = -1; return; end if

    ! Count whitespace-delimited tokens in Line 7
    n_all_chans = 0
    pos = 1
    l_len = len_trim(line_buffer)
    do while (pos <= l_len)
        do while (pos <= l_len .and. line_buffer(pos:pos) == ' ')
            pos = pos + 1
        end do
        if (pos > l_len) exit
        n_all_chans = n_all_chans + 1
        do while (pos <= l_len .and. line_buffer(pos:pos) /= ' ')
            pos = pos + 1
        end do
    end do

    if (n_all_chans == 0) then
        print *, "Error: No channel names found on line 7."
        close(fu); status = -1; return
    end if

    allocate(all_channels(n_all_chans))
    read(line_buffer, *) all_channels

    ! Line 8: Channel Units
    read(fu, '(A)', iostat=status) line_buffer
    if (status /= 0) then; close(fu); status = -1; return; end if

    allocate(all_units(n_all_chans))
    read(line_buffer, *) all_units

    if (available_channels) then
        print *, "Available Channels in .out: "
        do i = 1, n_all_chans
            print '(2A)', " - ", trim(all_channels(i))
        end do
    end if

    ! --- 2. Map Requested Channels & Locate Time Column ---
    allocate(col_map(n_req))
    allocate(units_out(n_req))
    time_col = 1 ! Default assumption

    do j = 1, n_all_chans
        if (trim(adjustl(all_channels(j))) == "Time") time_col = j
    end do

    do i = 1, n_req
        col_map(i) = 0
        do j = 1, n_all_chans
            if (trim(adjustl(plotChannels(i))) == trim(adjustl(all_channels(j)))) then
                col_map(i) = j
                units_out(i) = all_units(j)
                exit
            end if
        end do
        if (col_map(i) == 0) then
            print *, "Error: Channel '", trim(plotChannels(i)), "' not found in ASCII file."
            close(fu); status = -1; return
        end if
    end do

    ! --- 3. Fast Count Total Data Rows ---
    ! Al leer como texto '(A)' y un solo caracter, el puntero salta de linea de forma casi instantanea
    n_total = 0
    do
        read(fu, '(A)', iostat=ios) dummy_char
        if (ios /= 0) exit
        n_total = n_total + 1
    end do

    if (n_total == 0) then
        print *, "Error: No data rows found in ASCII file."
        close(fu); status = -1; return
    end if

    ! --- 4. Compute Filtering Windows (Matching Python & Binary Logic) ---
    idx_0_start = int(real(n_total, WP) * From_val)
    idx_0_end   = int(real(n_total, WP) * Upto_val)
    
    start_idx = idx_0_start + 1
    end_idx   = idx_0_end
    
    n_filtered = max(0, end_idx - start_idx + 1)
    if (n_filtered <= 0) then
        print *, "Warning: Filter range resulted in 0 output rows."
        close(fu); return
    end if

    ! --- 5. Extract Time Series Data ---
    rewind(fu)
    do i = 1, 8
        read(fu, '(A)') line_buffer
    end do

    allocate(time_out(n_filtered))
    allocate(array_out(n_filtered, n_req))
    allocate(row_buff(n_all_chans))

    out_row = 1
    do row_idx = 1, end_idx
        read(fu, *, iostat=ios) row_buff
        if (ios /= 0) exit

        if (row_idx >= start_idx) then
            time_out(out_row) = row_buff(time_col)
            do k = 1, n_req
                array_out(out_row, k) = row_buff(col_map(k))
            end do
            out_row = out_row + 1
        end if
    end do

    close(fu)

    if (verbose) then
        print *, "Output array shape: [", n_filtered, ", ", n_req, "]"
    end if

END SUBROUTINE PARSE_CHANNELS_ASCII


! -------------------------------------------------------------------------
! Function: read_curve
! Reads a two-line-header CSV lookup table (windspeed vs. one or more
! dependent columns) and linearly interpolates (with flat extrapolation
! outside the table's range) the requested column at the given wind speed.
! -------------------------------------------------------------------------
FUNCTION read_curve(filename, ws, col_idx) RESULT(val)
    CHARACTER(len=*), INTENT(IN)        :: filename
    REAL(WP), INTENT(IN)               :: ws        ! [m/s] WIndSpeed
    INTEGER(I32), INTENT(IN), OPTIONAL :: col_idx   ! [-] Column to evaluate
    REAL(WP)                           :: val       ! [-] Interpolated scalar

    INTEGER(I32) :: fu, status, ios, n_rows, i, j, k, target_col
    CHARACTER(len=2048) :: line_buffer
    REAL(WP), ALLOCATABLE :: ws_array(:), var_array(:), row_buff(:)
    REAL(WP) :: dx, t

    val = 0.0_WP
    target_col = 2
    if (present(col_idx)) target_col = col_idx

    open(newunit=fu, file=trim(filename), status='old', action='read', iostat=status)
    if (status /= 0) then
        print *, "Error: Could not open curve CSV file: ", trim(filename)
        return
    end if

    ! Skip CSV header: first line is the variable name, second line is the column names
    read(fu, '(A)', iostat=status) line_buffer
    if (status /= 0) then; close(fu); return; end if
    read(fu, '(A)', iostat=status) line_buffer
    if (status /= 0) then; close(fu); return; end if

    ! Count only numeric data rows
    allocate(row_buff(1))
    n_rows = 0
    do
        read(fu, '(A)', iostat=ios) line_buffer
        if (ios /= 0) exit
        if (len_trim(line_buffer) == 0) cycle
        ! Ignore header lines if they accidentally remain in the data block
        if (index(adjustl(line_buffer), 'X') == 1 .or. index(adjustl(line_buffer), 'Y') == 1) cycle

        ! Try to parse the first value as a real; if it fails, skip the line
        read(line_buffer, *, iostat=ios) row_buff(1)
        if (ios == 0) n_rows = n_rows + 1
    end do

    if (n_rows < 2) then
        print *, "Error: CSV requires at least 2 numeric points for interpolation."
        close(fu); deallocate(row_buff); return
    end if

    allocate(ws_array(n_rows))
    allocate(var_array(n_rows))
    deallocate(row_buff)
    allocate(row_buff(target_col))

    ! Read and parse data
    rewind(fu)
    read(fu, '(A)') line_buffer ! Skip first line
    read(fu, '(A)') line_buffer ! Skip second line

    i = 0
    do
        read(fu, '(A)', iostat=ios) line_buffer
        if (ios /= 0) exit
        if (len_trim(line_buffer) == 0) cycle

        ! Replace commas for blank spaces
        do j = 1, len_trim(line_buffer)
            if (line_buffer(j:j) == ',') line_buffer(j:j) = ' '
        end do

        read(line_buffer, *, iostat=ios) row_buff(1:target_col)
        if (ios /= 0) cycle

        i = i + 1
        ws_array(i)  = row_buff(1)
        var_array(i) = row_buff(target_col)
        if (i >= n_rows) exit
    end do

    close(fu)

    ! Linear interpolation with extrapolation
    if (ws <= ws_array(1)) then
        k = 1
    else if (ws >= ws_array(n_rows)) then
        k = n_rows - 1
    else
        k = 1
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


! -------------------- !
! HELPER FUNCTIONS
! -------------------- !

! -------------------------------------------------------------------------
! Subroutine: parse_member_nodes
! Tokenizes a whitespace-separated SD.sum.yaml member-mapping line (after
! stripping any trailing '#' comment) and extracts the last max_nodes
! tokens as integer node IDs for that member.
! -------------------------------------------------------------------------
SUBROUTINE parse_member_nodes(text_line, max_nodes, node_arr)
    CHARACTER(len=*), INTENT(IN) :: text_line           ! [-] Line read in SD.sum.yaml
    INTEGER(I32)    , INTENT(IN) :: max_nodes           ! [-] Amount of nodes per member and line

    INTEGER(I32), INTENT(OUT) :: node_arr(max_nodes)    ! [-] Array with node indices contained at the end of text_line

    ! Local variables
    INTEGER(I32), PARAMETER  :: max_tokens = 200        ! [-] Safety upper bound on number of whitespace-separated tokens
    INTEGER(I32), PARAMETER  :: token_len = 32          ! [-] Max characters per token
    CHARACTER(len=token_len) :: tokens(max_tokens)
    CHARACTER(len=len(text_line)) :: line_trim
    INTEGER(I32) :: n_tokens, i, istart, line_len, ierr
    LOGICAL :: in_token


    line_trim = adjustl(text_line)
    line_len  = len_trim(line_trim)

    do i = 1, line_len
        if (line_trim(i:i) == '#') line_trim(i:i) = ' '
    end do


    n_tokens = 0; in_token = .false.; istart = 0
    do i = 1, line_len
        if (line_trim(i:i) /= ' ') then
            if (.not. in_token) then
                in_token = .true.
                istart   = i
            end if
        else
            if (in_token) then
                in_token = .false.
                n_tokens = n_tokens + 1
                if (n_tokens > max_tokens) then
                    print*,'Error in parse_member_nodes: too many tokens in line (limit ', &
                                MAX_TOKENS, '): ', TRIM(text_line)
                    stop 1

                end if
                tokens(n_tokens) = line_trim(istart:i-1)
            end if
        end if
    end do


    ! Handle a trailing token if the line does not end with a blan space
    if (in_token) then
        n_tokens = n_tokens + 1
        if (n_tokens > max_tokens) then
            print *, 'Error in parse_member_nodes: too many tokens in line (limit ', &
                        MAX_TOKENS, '): ', TRIM(text_line)
            stop 1
        end if
        tokens(n_tokens) = line_trim(istart:line_len)
    end if

    ! Fast fall validation
    if (n_tokens < max_nodes) then
        print *, 'Error in parse_member_nodes: line has only ', n_tokens, &
                    ' tokens, but max_nodes = ', max_nodes, '. Line: ', TRIM(text_line)
        stop 1
    end if

    ! Extract and convert the LAST max_nodes tokens to integers
    do i = 1, max_nodes
        read(tokens(n_tokens - max_nodes + i), *, IOSTAT=ierr) node_arr(i)
        if (ierr /= 0) then
            print *, 'Error in parse_member_nodes: could not parse token "', &
                        TRIM(tokens(n_tokens - max_nodes + i)), '" as integer. Line: ', TRIM(text_line)
            stop 1
        end if
    end do

END SUBROUTINE parse_member_nodes


! -------------------------------------------------------------------------
! Function: extract_number_nodes
! Extracts the integer node count from an SD.sum.yaml line of the form
! "Nodes: # 33 x 9", i.e. the number found between the '#' and the 'x'.
! -------------------------------------------------------------------------
FUNCTION extract_number_nodes(text_line) RESULT(number_of_nodes)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: text_line          ! [-] Line read in SD.sum.yaml, e.g. "Nodes: # 33 x 9"

    INTEGER(I32) :: number_of_nodes                    ! [-] Number of nodes extracted from between '#' and 'x'

    ! Local variables
    INTEGER(I32) :: pos_hash, pos_x, ierr
    CHARACTER(len=LEN(text_line)) :: between_str


    ! Locate the '#' and the 'x' that follow it
    pos_hash = INDEX(text_line, '#')
    if (pos_hash == 0) then
        print *, 'ERROR in extract_number_nodes: no "#" found in line: ', TRIM(text_line)
        stop 1
    end if

    pos_x = INDEX(text_line(pos_hash+1:), 'x')
    if (pos_x == 0) then
        print *, 'ERROR in extract_number_nodes: no "x" found after "#" in line: ', TRIM(text_line)
        stop 1
    end if
    pos_x = pos_x + pos_hash   ! convert back to index within the full string

    ! Extract and convert the substring between '#' and 'x'
    between_str = text_line(pos_hash+1 : pos_x-1)

    read(between_str, *, IOSTAT=ierr) number_of_nodes
    if (ierr /= 0) then
        print *, 'ERROR in extract_number_nodes: could not parse "', TRIM(ADJUSTL(between_str)), &
                    '" as integer. Line: ', TRIM(text_line)
        stop 1
    end if

END FUNCTION extract_number_nodes


! -------------------------------------------------------------------------
! Function: to_upper
! Returns "str" with every lowercase ASCII letter converted to uppercase.
! -------------------------------------------------------------------------
PURE FUNCTION to_upper(str) RESULT(res)
    CHARACTER(len=*), INTENT(IN) :: str
    CHARACTER(len=LEN(str))      :: res
    INTEGER :: i, c
    res = str
    do i = 1, len(str)
        c = iachar(str(i:i))
        if (c >= 97 .and. c <= 122) res(i:i) = achar(c - 32)
    end do
END FUNCTION to_upper


! -----------------------------------------------------------------------
! Function: swap_extension
! Replaces the extension of "path" (everything after the last '.') with
! "new_ext". If "path" has no extension, "new_ext" is simply appended.
! -----------------------------------------------------------------------
PURE FUNCTION swap_extension(path, new_ext) RESULT(new_path)
    CHARACTER(len=*), INTENT(IN)  :: path
    CHARACTER(len=*), INTENT(IN)  :: new_ext
    CHARACTER(len=:), ALLOCATABLE :: new_path

    INTEGER :: dot_pos

    dot_pos = index(trim(path), '.', back=.true.)

    if (dot_pos > 0) then
        new_path = trim(path(1:dot_pos-1)) // trim(new_ext)
    else
        new_path = trim(path) // trim(new_ext)
    end if

END FUNCTION swap_extension


! -----------------------------------------------------------------------
! Subroutine: ensure_parent_dir
! Given a file path, creates its parent directory (and any missing
! intermediate directories) if it doesn't already exist. No-op if "path"
! has no directory component (i.e. it's a bare filename).
! -----------------------------------------------------------------------
SUBROUTINE ensure_parent_dir(path)
    CHARACTER(len=*), INTENT(IN) :: path

    CHARACTER(len=:), ALLOCATABLE :: dir
    INTEGER :: slash_pos

    slash_pos = index(trim(path), '/', back=.true.)
    if (slash_pos > 0) then
        dir = trim(path(1:slash_pos-1))
        CALL EXECUTE_COMMAND_LINE('mkdir -p "' // dir // '"')
    end if

END SUBROUTINE ensure_parent_dir


! -------------------- !
! HDF5 FUNCTIONS
! -------------------- !
! -------------------------------------------------------------------------
! Subroutine: save_variable
! Generic dispatcher: given a variable of assumed rank/type (CLASS(*)), it
! uses SELECT RANK / SELECT TYPE to identify the actual rank (0 to 3) and
! type (real/integer/complex/logical) of "var" and forwards the call to
! the matching save_to_hdf5 (save_hdf5_<type>_<rank>) specific procedure.
! -------------------------------------------------------------------------
SUBROUTINE save_variable(save_path, var_name, var)

    CHARACTER(len=*), INTENT(IN) :: save_path
    CHARACTER(len=*), INTENT(IN) :: var_name
    CLASS(*)        , INTENT(IN) :: var(..)
    
    SELECT RANK (v_rank => var)
    RANK(0) 
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP)) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP)); CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (LOGICAL) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        END SELECT
    RANK(1) 
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP)) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP)); CALL save_to_hdf5(save_path, trim(var_name), v)
        END SELECT
    RANK(2) 
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP)) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP)); CALL save_to_hdf5(save_path, trim(var_name), v)
        END SELECT
    RANK(3) 
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP)) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP)); CALL save_to_hdf5(save_path, trim(var_name), v)
        END SELECT
    END SELECT

END SUBROUTINE save_variable

! =========================================================================
! HDF5 save_variable specific procedures
! =========================================================================
! Each of the following procedures implements one (type, rank) overload of
! the save_variable generic interface: it opens (or creates) the target
! HDF5 file, creates (or overwrites) a dataset named var_name, writes the
! given data into it, and closes the file. They are otherwise identical
! in structure — only the HDF5 datatype and dataspace rank differ.
! =========================================================================
! ---------- SCALARS ---------- !
! save_variable overload: scalar REAL(WP)
SUBROUTINE save_hdf5_real_scalar(filepath, var_name, data_scalar)
    CHARACTER(len=*), INTENT(IN) :: filepath
    CHARACTER(len=*), INTENT(IN) :: var_name
    REAL(WP), INTENT(IN), TARGET :: data_scalar

    INTEGER(HID_T)                 :: file_id, dspace_id, dset_id
    INTEGER                        :: error
    LOGICAL                        :: file_exists, link_exists
    INTEGER(HSIZE_T), DIMENSION(1) :: dims

    CALL h5open_f(error)

    INQUIRE(FILE=filepath, EXIST=file_exists)
    IF (file_exists) THEN
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    ELSE
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    END IF

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    IF (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = 1
    CALL h5screate_simple_f(1, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_scalar, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
END SUBROUTINE save_hdf5_real_scalar

! save_variable overload: scalar INTEGER
SUBROUTINE save_hdf5_int_scalar(filepath, var_name, data_scalar)
    CHARACTER(len=*), INTENT(IN) :: filepath
    CHARACTER(len=*), INTENT(IN) :: var_name
    INTEGER, INTENT(IN), TARGET  :: data_scalar

    INTEGER(HID_T)                 :: file_id, dspace_id, dset_id
    INTEGER                        :: error
    LOGICAL                        :: file_exists, link_exists
    INTEGER(HSIZE_T), DIMENSION(1) :: dims

    CALL h5open_f(error)

    INQUIRE(FILE=filepath, EXIST=file_exists)
    IF (file_exists) THEN
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    ELSE
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    END IF

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    IF (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = 1
    CALL h5screate_simple_f(1, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_scalar, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
END SUBROUTINE save_hdf5_int_scalar

! save_variable overload: scalar COMPLEX(WP)
SUBROUTINE save_hdf5_cplx_scalar(filepath, var_name, data_scalar)
    CHARACTER(len=*), INTENT(IN) :: filepath
    CHARACTER(len=*), INTENT(IN) :: var_name
    COMPLEX(WP), INTENT(IN)      :: data_scalar

    INTEGER(HID_T)                 :: file_id, dspace_id, dset_id
    INTEGER                        :: error
    LOGICAL                        :: file_exists, link_exists
    INTEGER(HSIZE_T), DIMENSION(1) :: dims
    REAL(WP), TARGET               :: tmp(2)

    CALL h5open_f(error)

    INQUIRE(FILE=filepath, EXIST=file_exists)
    IF (file_exists) THEN
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    ELSE
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    END IF

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    IF (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    tmp(1) = REAL(data_scalar, WP)
    tmp(2) = AIMAG(data_scalar)

    dims(1) = 2
    CALL h5screate_simple_f(1, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tmp, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
END SUBROUTINE save_hdf5_cplx_scalar

! save_variable overload: scalar LOGICAL
SUBROUTINE save_hdf5_logical_scalar(filepath, var_name, data_scalar)
    CHARACTER(len=*), INTENT(IN) :: filepath
    CHARACTER(len=*), INTENT(IN) :: var_name
    LOGICAL, INTENT(IN), TARGET :: data_scalar

    INTEGER(HID_T)                 :: file_id, dspace_id, dset_id
    INTEGER                        :: error
    LOGICAL                        :: file_exists, link_exists
    INTEGER(HSIZE_T), DIMENSION(1) :: dims
    INTEGER                        :: int_val   ! Converted value

    int_val = MERGE(1, 0, data_scalar)

    CALL h5open_f(error)

    INQUIRE(FILE=filepath, EXIST=file_exists)
    IF (file_exists) THEN
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    ELSE
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    END IF

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    IF (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = 1
    CALL h5screate_simple_f(1, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, int_val, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
END SUBROUTINE save_hdf5_logical_scalar

! save_variable overload: scalar CHARACTER string
SUBROUTINE save_hdf5_char_scalar(filepath, var_name, data_scalar)
    CHARACTER(len=*), INTENT(IN) :: filepath
    CHARACTER(len=*), INTENT(IN) :: var_name
    CHARACTER(len=*), INTENT(IN) :: data_scalar

    INTEGER(HID_T)                 :: file_id, dspace_id, dset_id, dtype_id
    INTEGER                        :: error
    LOGICAL                        :: file_exists, link_exists
    INTEGER(HSIZE_T), DIMENSION(1) :: dims
    INTEGER(SIZE_T)                :: str_len

    CALL h5open_f(error)

    INQUIRE(FILE=filepath, EXIST=file_exists)
    IF (file_exists) THEN
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    ELSE
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    END IF

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    IF (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    ! Obtain trimed text length
    str_len = INT(LEN_TRIM(data_scalar), KIND=SIZE_T)
    IF (str_len == 0_SIZE_T) str_len = 1_SIZE_T

    ! Create datatype with adecuate lenght
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

! ---------- ARRAYS 1D ---------- !
! save_variable overload: rank-1 REAL(WP) array
SUBROUTINE save_hdf5_real_1d(filepath, var_name, data_array)
    CHARACTER(len=*), INTENT(IN) :: filepath
    CHARACTER(len=*), INTENT(IN) :: var_name
    REAL(WP), INTENT(IN), TARGET :: data_array(:)

    INTEGER(HID_T)                 :: file_id, dspace_id, dset_id
    INTEGER                        :: error
    LOGICAL                        :: file_exists, link_exists
    INTEGER(HSIZE_T), DIMENSION(1) :: dims

    CALL h5open_f(error)

    INQUIRE(FILE=filepath, EXIST=file_exists)
    IF (file_exists) THEN
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    ELSE
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    END IF

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    IF (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = SIZE(data_array, 1)
    CALL h5screate_simple_f(1, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_array, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
END SUBROUTINE save_hdf5_real_1d

! save_variable overload: rank-1 INTEGER array
SUBROUTINE save_hdf5_int_1d(filepath, var_name, data_array)
    CHARACTER(len=*), INTENT(IN) :: filepath
    CHARACTER(len=*), INTENT(IN) :: var_name
    INTEGER, INTENT(IN), TARGET  :: data_array(:)

    INTEGER(HID_T)                 :: file_id, dspace_id, dset_id
    INTEGER                        :: error
    LOGICAL                        :: file_exists, link_exists
    INTEGER(HSIZE_T), DIMENSION(1) :: dims

    CALL h5open_f(error)

    INQUIRE(FILE=filepath, EXIST=file_exists)
    IF (file_exists) THEN
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    ELSE
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    END IF

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    IF (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = SIZE(data_array, 1)
    CALL h5screate_simple_f(1, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_array, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
END SUBROUTINE save_hdf5_int_1d

! save_variable overload: rank-1 COMPLEX(WP) array
SUBROUTINE save_hdf5_cplx_1d(filepath, var_name, data_array)
    CHARACTER(len=*), INTENT(IN) :: filepath
    CHARACTER(len=*), INTENT(IN) :: var_name
    COMPLEX(WP), INTENT(IN)      :: data_array(:)

    INTEGER(HID_T)                 :: file_id, dspace_id, dset_id
    INTEGER                        :: error
    LOGICAL                        :: file_exists, link_exists
    INTEGER(HSIZE_T), DIMENSION(2) :: dims
    REAL(WP), ALLOCATABLE, TARGET  :: tmp(:,:)

    CALL h5open_f(error)

    INQUIRE(FILE=filepath, EXIST=file_exists)
    IF (file_exists) THEN
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    ELSE
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    END IF

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    IF (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = 2
    dims(2) = SIZE(data_array, 1)
    
    ALLOCATE(tmp(2, dims(2)))
    tmp(1,:) = REAL(data_array, WP)
    tmp(2,:) = AIMAG(data_array)

    CALL h5screate_simple_f(2, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tmp, dims, error)

    DEALLOCATE(tmp)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
END SUBROUTINE save_hdf5_cplx_1d

! ---------- ARRAYS 2D ---------- !
! save_variable overload: rank-2 REAL(WP) array
SUBROUTINE save_hdf5_real_2d(filepath, var_name, data_array)
    CHARACTER(len=*), INTENT(IN) :: filepath
    CHARACTER(len=*), INTENT(IN) :: var_name
    REAL(WP), INTENT(IN), TARGET :: data_array(:,:)

    INTEGER(HID_T)                 :: file_id, dspace_id, dset_id
    INTEGER                        :: error
    LOGICAL                        :: file_exists, link_exists
    INTEGER(HSIZE_T), DIMENSION(2) :: dims

    CALL h5open_f(error)

    INQUIRE(FILE=filepath, EXIST=file_exists)
    IF (file_exists) THEN
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    ELSE
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    END IF

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    IF (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = SIZE(data_array, 1)
    dims(2) = SIZE(data_array, 2)
    CALL h5screate_simple_f(2, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_array, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
END SUBROUTINE save_hdf5_real_2d

! save_variable overload: rank-2 INTEGER array
SUBROUTINE save_hdf5_int_2d(filepath, var_name, data_array)
    CHARACTER(len=*), INTENT(IN) :: filepath
    CHARACTER(len=*), INTENT(IN) :: var_name
    INTEGER, INTENT(IN), TARGET  :: data_array(:,:)

    INTEGER(HID_T)                 :: file_id, dspace_id, dset_id
    INTEGER                        :: error
    LOGICAL                        :: file_exists, link_exists
    INTEGER(HSIZE_T), DIMENSION(2) :: dims

    CALL h5open_f(error)

    INQUIRE(FILE=filepath, EXIST=file_exists)
    IF (file_exists) THEN
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    ELSE
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    END IF

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    IF (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = SIZE(data_array, 1)
    dims(2) = SIZE(data_array, 2)
    CALL h5screate_simple_f(2, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_array, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
END SUBROUTINE save_hdf5_int_2d

! save_variable overload: rank-2 COMPLEX(WP) array
SUBROUTINE save_hdf5_cplx_2d(filepath, var_name, data_array)
    CHARACTER(len=*), INTENT(IN) :: filepath
    CHARACTER(len=*), INTENT(IN) :: var_name
    COMPLEX(WP), INTENT(IN)      :: data_array(:,:)

    INTEGER(HID_T)                 :: file_id, dspace_id, dset_id
    INTEGER                        :: error
    LOGICAL                        :: file_exists, link_exists
    INTEGER(HSIZE_T), DIMENSION(3) :: dims
    REAL(WP), ALLOCATABLE, TARGET  :: tmp(:,:,:)

    CALL h5open_f(error)

    INQUIRE(FILE=filepath, EXIST=file_exists)
    IF (file_exists) THEN
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    ELSE
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    END IF

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    IF (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = 2
    dims(2) = SIZE(data_array, 1)
    dims(3) = SIZE(data_array, 2)

    ALLOCATE(tmp(2, dims(2), dims(3)))
    tmp(1,:,:) = REAL(data_array, WP)
    tmp(2,:,:) = AIMAG(data_array)

    CALL h5screate_simple_f(3, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tmp, dims, error)

    DEALLOCATE(tmp)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
END SUBROUTINE save_hdf5_cplx_2d

! ---------- ARRAYS 3D ---------- !
! save_variable overload: rank-3 REAL(WP) array
SUBROUTINE save_hdf5_real_3d(filepath, var_name, data_array)
    CHARACTER(len=*), INTENT(IN) :: filepath
    CHARACTER(len=*), INTENT(IN) :: var_name
    REAL(WP), INTENT(IN), TARGET :: data_array(:,:,:)

    INTEGER(HID_T)                 :: file_id, dspace_id, dset_id
    INTEGER                        :: error
    LOGICAL                        :: file_exists, link_exists
    INTEGER(HSIZE_T), DIMENSION(3) :: dims

    CALL h5open_f(error)

    INQUIRE(FILE=filepath, EXIST=file_exists)
    IF (file_exists) THEN
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    ELSE
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    END IF

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    IF (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = SIZE(data_array, 1)
    dims(2) = SIZE(data_array, 2)
    dims(3) = SIZE(data_array, 3)
    CALL h5screate_simple_f(3, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_array, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
END SUBROUTINE save_hdf5_real_3d

! save_variable overload: rank-3 INTEGER array
SUBROUTINE save_hdf5_int_3d(filepath, var_name, data_array)
    CHARACTER(len=*), INTENT(IN) :: filepath
    CHARACTER(len=*), INTENT(IN) :: var_name
    INTEGER, INTENT(IN), TARGET  :: data_array(:,:,:)

    INTEGER(HID_T)                 :: file_id, dspace_id, dset_id
    INTEGER                        :: error
    LOGICAL                        :: file_exists, link_exists
    INTEGER(HSIZE_T), DIMENSION(3) :: dims

    CALL h5open_f(error)

    INQUIRE(FILE=filepath, EXIST=file_exists)
    IF (file_exists) THEN
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    ELSE
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    END IF

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    IF (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = SIZE(data_array, 1)
    dims(2) = SIZE(data_array, 2)
    dims(3) = SIZE(data_array, 3)
    CALL h5screate_simple_f(3, dims, dspace_id, error)

    CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_array, dims, error)

    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
END SUBROUTINE save_hdf5_int_3d

! save_variable overload: rank-3 COMPLEX(WP) array
SUBROUTINE save_hdf5_cplx_3d(filepath, var_name, data_array)
    CHARACTER(len=*), INTENT(IN) :: filepath
    CHARACTER(len=*), INTENT(IN) :: var_name
    COMPLEX(WP), INTENT(IN)      :: data_array(:,:,:)

    INTEGER(HID_T)                 :: file_id, dspace_id, dset_id
    INTEGER                        :: error
    LOGICAL                        :: file_exists, link_exists
    INTEGER(HSIZE_T), DIMENSION(4) :: dims
    REAL(WP), ALLOCATABLE, TARGET  :: tmp(:,:,:,:)

    CALL h5open_f(error)

    INQUIRE(FILE=filepath, EXIST=file_exists)
    IF (file_exists) THEN
        CALL h5fopen_f(filepath, H5F_ACC_RDWR_F, file_id, error)
    ELSE
        CALL h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_id, error)
    END IF

    CALL h5lexists_f(file_id, var_name, link_exists, error)
    IF (link_exists) CALL h5ldelete_f(file_id, var_name, error)

    dims(1) = 2
    dims(2) = SIZE(data_array, 1)
    dims(3) = SIZE(data_array, 2)
    dims(4) = SIZE(data_array, 3)

    ALLOCATE(tmp(2, dims(2), dims(3), dims(4)))
    tmp(1,:,:,:) = REAL(data_array, WP)
    tmp(2,:,:,:) = AIMAG(data_array)

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