! ===============================================================================
! MODULE: HDF5Utils
! ===============================================================================
! Defines the HDF5 input/output layer, the only module that uses HDF5: a
! generic variable-saving facility (the save_variable dispatcher, whose
! CLASS(*) assumed-rank argument is forwarded to the save_to_hdf5 generic
! interface: real/integer/complex data of rank 0 to 3, plus logical and
! character scalars; complex arrays are stored as (2, N...) in Fortran order,
! i.e. (..., N, 2) when read with h5py), and the reader of a rank-1 real
! dataset (read_hdf5_real_1d).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE HDF5Utils

USE HDF5
USE Kinds, ONLY: WP

IMPLICIT NONE

PRIVATE
PUBLIC :: save_variable, read_hdf5_real_1d


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
! HDF5 READER
! ===============================================================================
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

END MODULE HDF5Utils
