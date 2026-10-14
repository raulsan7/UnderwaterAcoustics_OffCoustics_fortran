! ===============================================================================
! SUBMODULE: AcousticSolver_output
! ===============================================================================
! Defines the output-file management of AcousticSolver_t: the HDF5 output
! path (resolve_save_path), the update of an existing output of the same
! case (check_update, then update_parameters with the merge of the
! previous .dat), the .dat block of every post-processing step that ran
! (update_parameters), the datasets shared by every solver
! (save_shared_parameters) and the water lines of the .dat summary
! (write_txt_water).
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
SUBMODULE (AcousticSolver) AcousticSolver_output

USE StringUtils, ONLY: to_upper
USE HDF5Utils, ONLY: save_variable, read_hdf5_real_1d
USE FileUtils, ONLY: read_full_line, swap_extension, ensure_parent_dir, move_text_file, delete_file

IMPLICIT NONE


CONTAINS


! ===============================================================================
! OUTPUT FILE AND UPDATE OF AN EXISTING OUTPUT
! ===============================================================================
MODULE SUBROUTINE resolve_save_path(self, method_tag)
    ! -------------------------------------------------------------------------
    ! Resolves the final HDF5 output path, shared by every solver's
    ! build_save_path. If the configured save_path is "None"/"Default", the
    ! name is built from the emitter labels carried by the sources around the
    ! solver's method tag: out_dir // name_tag // method_tag // name_suffix,
    ! e.g. ./turbine_acoustic_data/plot_mn_ANM24_ws11.4.hdf5 or
    ! ./farm_acoustic_data/Farm_2_DTU10MN_ANM24_ws8.0-11.4.hdf5. Otherwise the
    ! .hdf5 extension is appended to the given path. In both cases, a set
    ! extra_label (neither 'None' nor empty) is appended before the extension
    ! as "_" // extra_label. The destination folder is created if needed. A
    ! previous file with the same name is deleted when update_existing is
    ! .false. (save_variable only overwrites the datasets it writes, so those
    ! of an older run would otherwise survive). With update_existing it is
    ! kept and its .dat moved to .dat.prev: check_update then decides whether
    ! it holds the same case and update_parameters merges the step blocks.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(INOUT) :: self          ! [-] Solver instance
    CHARACTER(len=*)       , INTENT(IN)    :: method_tag    ! [-] Method and resolution tag, e.g. "_ANM24"

    ! Local variables
    CHARACTER(len=512)            :: path_upper
    CHARACTER(len=:), ALLOCATABLE :: label_str, dat_path
    INTEGER(I32)                  :: u, ierr


    ! Optional user label: 'None' (any case) or an empty string adds nothing
    label_str = ''
    if (len_trim(self % extra_label) > 0 .and. to_upper(trim(self % extra_label)) /= 'NONE') &
        label_str = '_' // trim(self % extra_label)

    path_upper = to_upper(trim(self % save_path))

    if (trim(path_upper) == 'NONE' .or. trim(path_upper) == 'DEFAULT') then
        associate(g => self % sources % geometry)
            self % save_path = trim(g % out_dir) // trim(g % name_tag) // method_tag // trim(g % name_suffix) // &
                               label_str // ".hdf5"
        end associate
    else
        self % save_path = trim(self % save_path) // label_str // ".hdf5"
    end if

    ! Make sure the destination folder exists before anything tries to write into it
    CALL ensure_parent_dir(trim(self % save_path))

    ! A stale .dat.prev (e.g. from an interrupted run) must never be merged
    dat_path = swap_extension(trim(self % save_path), ".dat")
    CALL delete_file(dat_path // ".prev")

    open(newunit=u, file=trim(self % save_path), status='old', iostat=ierr)
    if (ierr == 0_I32) then
        if (self % update_existing) then
            ! Kept for an update: its .dat is set aside for check_update and update_parameters
            close(u)
            CALL move_text_file(dat_path, dat_path // ".prev")
            write(*,*) '   ---> Previous output file kept for update: ', trim(self % save_path)
        else
            ! Start from an empty file, as the .dat summary does (status='replace')
            close(u, status='delete')
            write(*,*) '   ---> Previous output file replaced: ', trim(self % save_path)
        end if
    end if

END SUBROUTINE resolve_save_path


MODULE SUBROUTINE check_update(self, emitter_obj)
    ! -------------------------------------------------------------------------
    ! Called by save_parameters before anything is written. If resolve_save_path
    ! kept an existing output file (update_existing), it is updated only if it
    ! holds the same case: the same frequency vector (OpenFAST run, From/Upto,
    ! skipf, frequency filtering) and the same .dat header (emitter,
    ! environment and solver settings; block size and cluster flag excluded).
    ! The new header is written here to compare it. Otherwise the file and its
    ! .dat.prev are deleted with a warning and the output is rebuilt.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(IN) :: self       ! [-] Solver instance
    CLASS(Emitter_t)       , INTENT(IN) :: emitter_obj    ! [-] Emitter that produced the sources

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: dat_path, prev_path, reason
    REAL(WP)        , ALLOCATABLE :: freqs_old(:)
    LOGICAL                       :: h5_exists, prev_exists, found


    inquire(file=trim(self % save_path), exist=h5_exists)
    if (.not. h5_exists) return                 ! Fresh output: nothing to update

    dat_path  = swap_extension(trim(self % save_path), ".dat")
    prev_path = dat_path // ".prev"
    reason    = ''

    ! --- 1. Same frequency vector --- !
    CALL read_hdf5_real_1d(trim(self % save_path), "Freqs", freqs_old, found)
    if (.not. found) then
        reason = 'no Freqs in the file'
    else if (size(freqs_old) /= size(self % sources % Freqs)) then
        reason = 'different number of frequencies'
    else if (maxval(abs(freqs_old - self % sources % Freqs)) > 1.0e-9_WP * maxval(abs(self % sources % Freqs))) then
        reason = 'different frequencies'
    end if

    ! --- 2. Same emitter, environment and solver (.dat header) --- !
    if (len(reason) == 0) then
        inquire(file=prev_path, exist=prev_exists)
        if (.not. prev_exists) then
            reason = 'no previous .dat to compare with'
        else
            CALL self % write_txt_summary(emitter_obj)
            CALL compare_dat_headers(dat_path, prev_path, reason)
        end if
    end if

    if (len(reason) == 0) then
        write(*,*) '   ---> Same case: only the steps run now are rewritten in ', trim(self % save_path)
    else
        write(*,*) 'WARNING check_update: ' // trim(self % save_path) // ' not updated (' // reason // &
                   '): deleted and rebuilt'
        CALL delete_file(trim(self % save_path))
        CALL delete_file(prev_path)
    end if

END SUBROUTINE check_update


SUBROUTINE compare_dat_headers(dat_path, prev_path, reason)
    ! -------------------------------------------------------------------------
    ! Compares the header just written in dat_path with the header of the
    ! previous .dat (everything before its post-processing block), ignoring
    ! the block size and cluster lines, which do not change the results.
    ! "reason" is empty if both match, otherwise it names the first difference.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*)             , INTENT(IN)  :: dat_path     ! [-] .dat holding only the new header
    CHARACTER(len=*)             , INTENT(IN)  :: prev_path    ! [-] Previous .dat
    CHARACTER(len=:), ALLOCATABLE, INTENT(OUT) :: reason       ! [-] '' if equal, else the first difference

    ! Local variables
    CHARACTER(len=256), ALLOCATABLE :: new_lines(:), old_lines(:)
    INTEGER(I32)                    :: i


    CALL read_header_lines(dat_path, new_lines)
    CALL read_header_lines(prev_path, old_lines)

    reason = ''
    do i = 1, min(size(new_lines), size(old_lines))
        if (new_lines(i) /= old_lines(i)) then
            reason = 'different setting: "' // trim(adjustl(old_lines(i))) // '" -> "' // trim(adjustl(new_lines(i))) // '"'
            return
        end if
    end do
    if (size(new_lines) /= size(old_lines)) reason = 'different .dat header length'

END SUBROUTINE compare_dat_headers


SUBROUTINE read_header_lines(path, lines)
    ! -------------------------------------------------------------------------
    ! Lines of a .dat summary before its post-processing block (the dashed
    ! separator ahead of "POST-PROCESSING STEPS EXECUTED" excluded), without
    ! trailing blanks and leaving out the block size and cluster lines.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*)               , INTENT(IN)  :: path        ! [-] .dat file
    CHARACTER(len=256), ALLOCATABLE, INTENT(OUT) :: lines(:)    ! [-] Header lines

    ! Local variables
    CHARACTER(len=256), ALLOCATABLE :: buffer(:)
    CHARACTER(len=256)              :: line
    INTEGER(I32)                    :: u, ios, n


    ALLOCATE(buffer(0))
    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios == 0_I32) then
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0_I32) exit
            if (trim(line) == 'POST-PROCESSING STEPS EXECUTED') exit
            if (index(line, 'Block size') == 1 .or. index(line, 'In cluster') == 1) cycle
            buffer = [CHARACTER(len=256) :: buffer, line]
        end do
        close(u)
    end if

    ! Drop the dashed separator written just before the post-processing title
    n = size(buffer)
    if (n > 0) then
        if (trim(buffer(n)) == repeat('-', 70)) n = n - 1_I32
    end if
    lines = buffer(1:n)

END SUBROUTINE read_header_lines


SUBROUTINE write_previous_block(u, prev_path, header)
    ! -------------------------------------------------------------------------
    ! Copies to the open .dat unit "u" the block of one post-processing step
    ! (header such as 'Polar:') from the previous .dat, if it has one, tagged
    ! "(previous run)": that step did not run now and its datasets were kept.
    ! A block ends at the first blank or non-indented line.
    ! -------------------------------------------------------------------------

    INTEGER(I32)    , INTENT(IN) :: u            ! [-] Unit of the open .dat file
    CHARACTER(len=*), INTENT(IN) :: prev_path    ! [-] Previous .dat ('' or missing: nothing to copy)
    CHARACTER(len=*), INTENT(IN) :: header       ! [-] Step header, e.g. 'Polar:'

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: line
    INTEGER(I32)                  :: u_prev, ios
    LOGICAL                       :: in_block


    open(newunit=u_prev, file=prev_path, status='old', action='read', iostat=ios)
    if (ios /= 0_I32) return

    in_block = .false.
    do
        CALL read_full_line(u_prev, line, ios)
        if (ios /= 0_I32) exit

        if (.not. in_block) then
            if (len(line) >= len(header)) then
                if (line(1:len(header)) == header) then
                    in_block = .true.
                    write(u, '(A)') ''
                    write(u, '(A)') header // ' (previous run)'
                end if
            end if
        else
            if (len_trim(line) == 0) exit
            if (line(1:1) /= ' ') exit
            write(u, '(A)') line
        end if
    end do

    close(u_prev)

END SUBROUTINE write_previous_block


MODULE SUBROUTINE update_parameters(self)
    ! -------------------------------------------------------------------------
    ! Appends to the .dat summary already written by save_parameters() one
    ! block per post-processing step that actually ran (per do_run), with the
    ! resolved parameters run() used for it, then closes the file with its
    ! final separator. In an update of an existing output (.dat.prev kept by
    ! resolve_save_path and check_update), the steps that did not run now
    ! keep their block from the previous .dat, tagged "(previous run)", and
    ! .dat.prev is deleted. Must be called after run(), since every step only
    ! resolves its own DEFAULT_R/DEFAULT_I sentinels when it actually
    ! executes.
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(IN) :: self    ! [-] Solver instance

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: dat_path, prev_path
    INTEGER(I32)                  :: u, i


    dat_path  = swap_extension(trim(self % save_path), ".dat")
    prev_path = dat_path // ".prev"

    open(newunit=u, file=dat_path, status='old', position='append', action='write')

    write(u, '(A)') repeat('-', 70)
    write(u, '(A)') 'POST-PROCESSING STEPS EXECUTED'
    write(u, '(A)') repeat('-', 70)

    associate(d => self % run_inputs % do_run, pp => self % run_inputs)

        if (d % do_spectrums) then
            write(u, '(A)') ''
            write(u, '(A)') 'Spectrums:'
            if (.not. ALLOCATED(pp % params_spectrums % observers_x)) then
                write(u, '(A)') '   observers : 3 default points (auto, x=10/250/500 m, z=-Depth/2)'
            else
                write(u, '(A,I0)') '   n_observers : ', pp % params_spectrums % n_observers
                do i = 1, pp % params_spectrums % n_observers
                    write(u, '(A,I0,A,F8.2,A,F8.2,A,F8.2,A)') '   obs ', i, ' (x,y,z) : ', &
                                                              pp % params_spectrums % observers_x(i), ', ', &
                                                              pp % params_spectrums % observers_y(i), ', ', &
                                                              pp % params_spectrums % observers_z(i), ' m'
                end do
            end if
        else
            CALL write_previous_block(u, prev_path, 'Spectrums:')
        end if

        if (d % do_polar) then
            write(u, '(A)') ''
            write(u, '(A)')        'Polar:'
            write(u, '(A,F8.2,A)') '   r_polar      : ', pp % params_polar % r_polar, ' m'
            write(u, '(A,F8.2,A)') '   z_polar      : ', pp % params_polar % z_polar, ' m'
            write(u, '(A,I0)')     '   n_polar      : ', pp % params_polar % n_polar
            write(u, '(A,2F10.3)') '   center_polar : ', pp % params_polar % center_polar
        else
            CALL write_previous_block(u, prev_path, 'Polar:')
        end if

        if (d % do_cylinder) then
            write(u, '(A)') ''
            write(u, '(A)')        'Cylinder:'
            write(u, '(A,F8.2,A)') '   r_cyl      : ', pp % params_cylinder % r_cyl, ' m'
            write(u, '(A,I0)')     '   n_cyl      : ', pp % params_cylinder % n_cyl
            write(u, '(A,I0)')     '   nz_cyl     : ', pp % params_cylinder % nz_cyl
            write(u, '(A,2F10.3)') '   center_cyl : ', pp % params_cylinder % center_cyl
        else
            CALL write_previous_block(u, prev_path, 'Cylinder:')
        end if

        if (d % do_decay) then
            write(u, '(A)') ''
            write(u, '(A)')        'Decay:'
            write(u, '(A,2F10.3)') '   distance_decay : ', pp % params_decay % distance_decay
            write(u, '(A,I0)')     '   n_decay        : ', pp % params_decay % n_decay
            write(u, '(A,F8.2,A)') '   z_decay        : ', pp % params_decay % z_decay, ' m'
            write(u, '(A,L1)')     '   logspace_decay : ', pp % params_decay % logspace_decay
        else
            CALL write_previous_block(u, prev_path, 'Decay:')
        end if

        if (d % do_line) then
            write(u, '(A)') ''
            write(u, '(A)')        'Line:'
            write(u, '(A,3F10.3)') '   p1_line       : ', pp % params_line % p1_line
            write(u, '(A,3F10.3)') '   p2_line       : ', pp % params_line % p2_line
            write(u, '(A,I0)')     '   n_line        : ', pp % params_line % n_line
            write(u, '(A,L1)')     '   logspace_line : ', pp % params_line % logspace_line
        else
            CALL write_previous_block(u, prev_path, 'Line:')
        end if

        if (d % do_sliceXY) then
            write(u, '(A)') ''
            write(u, '(A)')        'SliceXY:'
            write(u, '(A,F8.2,A)') '   z_xy       : ', pp % params_sliceXY % z_xy, ' m'
            write(u, '(A,I0)')     '   nx_xy      : ', pp % params_sliceXY % nx_xy
            write(u, '(A,I0)')     '   ny_xy      : ', pp % params_sliceXY % ny_xy
            write(u, '(A,2F10.3)') '   xlim_xy    : ', pp % params_sliceXY % xlim_xy
            write(u, '(A,2F10.3)') '   ylim_xy    : ', pp % params_sliceXY % ylim_xy
            write(u, '(A,2F10.3)') '   center_xy  : ', pp % params_sliceXY % center_xy
        else
            CALL write_previous_block(u, prev_path, 'SliceXY:')
        end if

        if (d % do_sliceXZ) then
            write(u, '(A)') ''
            write(u, '(A)')        'SliceXZ:'
            write(u, '(A,F8.2,A)') '   y_xz       : ', pp % params_sliceXZ % y_xz, ' m'
            write(u, '(A,I0)')     '   nx_xz      : ', pp % params_sliceXZ % nx_xz
            write(u, '(A,I0)')     '   nz_xz      : ', pp % params_sliceXZ % nz_xz
            write(u, '(A,2F10.3)') '   xlim_xz    : ', pp % params_sliceXZ % xlim_xz
            write(u, '(A,2F10.3)') '   zlim_xz    : ', pp % params_sliceXZ % zlim_xz
        else
            CALL write_previous_block(u, prev_path, 'SliceXZ:')
        end if

        if (d % do_sliceVertical) then
            write(u, '(A)') ''
            write(u, '(A)')         'SliceVertical:'
            write(u, '(A,F8.2,A)')  '   azimuth_sv : ', pp % params_sliceVertical % azimuth_sv, ' deg'
            write(u, '(A,2F10.3)')  '   center_sv  : ', pp % params_sliceVertical % center_sv
            write(u, '(A,F10.3,A)') '   width_sv   : ', pp % params_sliceVertical % width_sv, ' m'
            write(u, '(A,2F10.3)')  '   ulim_sv    : ', pp % params_sliceVertical % ulim_sv
            write(u, '(A,2F10.3)')  '   zlim_sv    : ', pp % params_sliceVertical % zlim_sv
            write(u, '(A,I0)')      '   nu_sv      : ', pp % params_sliceVertical % nu_sv
            write(u, '(A,I0)')      '   nz_sv      : ', pp % params_sliceVertical % nz_sv
            write(u, '(A,F8.2,A)')  '   margin_sv  : ', pp % params_sliceVertical % margin_sv, ' m'
        else
            CALL write_previous_block(u, prev_path, 'SliceVertical:')
        end if

        if (d % do_sphere) then
            write(u, '(A)') ''
            write(u, '(A)')        'Sphere:'
            write(u, '(A,F8.2,A)') '   r_sphere       : ', pp % params_sphere % r_sphere, ' m'
            write(u, '(A,I0)')     '   n_theta_sphere : ', pp % params_sphere % n_theta_sphere
            write(u, '(A,I0)')     '   nz_sphere      : ', pp % params_sphere % nz_sphere
            write(u, '(A,F8.2,A)') '   f_farfield     : ', pp % params_sphere % f_farfield_sphere, ' Hz'
            write(u, '(A,3F10.3)') '   center_sphere  : ', pp % params_sphere % center_sphere
        else
            CALL write_previous_block(u, prev_path, 'Sphere:')
        end if

    end associate

    write(u, *)
    write(u, '(A)') repeat('=', 70)

    close(u)
    CALL delete_file(prev_path)

    write(*,*) '   ---> Parameters updated at ', dat_path

END SUBROUTINE update_parameters


! ===============================================================================
! DATASETS SHARED BY EVERY SOLVER
! ===============================================================================
MODULE SUBROUTINE save_shared_parameters(self, emitter_obj)
    ! -------------------------------------------------------------------------
    ! Writes the datasets shared by every solver, called by each solver's
    ! save_parameters right after its own datasets: the solver name (Method,
    ! from get_name) and, only when the emitter forces come from a windowed
    ! FFT (get_enbw > 0), the equivalent noise bandwidth of that window in
    ! bins (ENBW). Neither is written in the .dat summary. Then the water of
    ! the volume absorption: temperature, salinity and pH when given (T_wat,
    ! S_wat, pH_wat) and, with env % absorption, the depth of the formula
    ! and alpha(f) in dB/km (z_alpha, Alpha_dBkm), whose presence tells that
    ! the pressures already include the absorption (write_txt_water writes
    ! the same in the .dat summary).
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(IN) :: self       ! [-] Solver instance
    CLASS(Emitter_t)       , INTENT(IN) :: emitter_obj    ! [-] Emitter that produced the sources

    ! Local variables
    REAL(WP) :: enbw


    CALL save_variable(self % save_path, "Method", trim(self % solver_name))

    ! Window of the forces: only when they come from a windowed FFT
    enbw = emitter_obj % get_enbw()
    if (enbw > 0.0_WP) CALL save_variable(self % save_path, "ENBW", enbw)

    ! Water of the volume absorption: the properties given, and alpha(f) when the pressures include it
    if (self % env % temperature /= DEFAULT_R) CALL save_variable(self % save_path, "T_wat", self % env % temperature)
    if (self % env % salinity /= DEFAULT_R)    CALL save_variable(self % save_path, "S_wat", self % env % salinity)
    if (self % env % pH /= DEFAULT_R)          CALL save_variable(self % save_path, "pH_wat", self % env % pH)

    if (self % env % absorption) then
        CALL save_variable(self % save_path, "z_alpha", self % z_alpha)
        CALL save_variable(self % save_path, "Alpha_dBkm", self % alpha_dB)
    end if

END SUBROUTINE save_shared_parameters


MODULE SUBROUTINE write_txt_water(self, u)
    ! -------------------------------------------------------------------------
    ! Appends to the ENVIRONMENT block of the .dat summary, open in unit u,
    ! the water of the volume absorption: one line per given temperature,
    ! salinity and pH and, with env % absorption, the formula, its depth and
    ! alpha at the lowest and highest frequency. Writes nothing when none
    ! of them is set, so the summaries of the cases without absorption
    ! keep their former content (check_update compares them).
    ! -------------------------------------------------------------------------

    CLASS(AcousticSolver_t), INTENT(IN) :: self    ! [-] Solver instance
    INTEGER(I32)           , INTENT(IN) :: u       ! [-] Unit of the open .dat file

    ! Local variables
    INTEGER(I32) :: i_min, i_max


    if (self % env % temperature /= DEFAULT_R) &
        write(u, '(A,F8.2,A)') 'Temperature      : ', self % env % temperature, ' degC'
    if (self % env % salinity /= DEFAULT_R) &
        write(u, '(A,F8.2,A)') 'Salinity         : ', self % env % salinity, ' ppt'
    if (self % env % pH /= DEFAULT_R) &
        write(u, '(A,F8.2)')   'pH               : ', self % env % pH

    if (.not. self % env % absorption) return

    i_min = int(minloc(self % sources % Freqs, dim=1), I32)
    i_max = int(maxloc(self % sources % Freqs, dim=1), I32)

    write(u, '(A)')                'Absorption       : Ainslie & McColm (1998)'
    write(u, '(A,F8.2,A)')         'Absorption depth : ', self % z_alpha, ' m'
    write(u, '(A,ES10.3,A,F0.2,A)') 'Alpha at f_min   : ', self % alpha_dB(i_min), ' dB/km at ', &
                                    self % sources % Freqs(i_min), ' Hz'
    write(u, '(A,ES10.3,A,F0.2,A)') 'Alpha at f_max   : ', self % alpha_dB(i_max), ' dB/km at ', &
                                    self % sources % Freqs(i_max), ' Hz'

END SUBROUTINE write_txt_water

END SUBMODULE AcousticSolver_output
