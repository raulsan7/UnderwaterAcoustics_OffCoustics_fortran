! ===============================================================================
! MODULE: DTU10MWFloating
! ===============================================================================
! Defines DTU10MWFloating_t, a concrete turbine type extending the abstract
! WindTurbine_t (src/core/WindTurbine.f90) for the DTU 10 MW floating platform
! (CENER DeltaWind version). It implements the deferred procedures:
! compute_force, compute_mass, filter_frequencies, save/write_specific_params
! and get_impedance_corrected_force, using a representative column/pontoon
! member to build the platform's distributed mass and scattering correction.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE DTU10MWFloating

USE Kinds
USE Types, ONLY: Config_WindTurbine_t
USE WindTurbine, ONLY: WindTurbine_t, drivetrain10MW_excitation_spectrum, init

IMPLICIT NONE

PRIVATE
PUBLIC :: DTU10MWFloating_t


TYPE, EXTENDS(WindTurbine_t) :: DTU10MWFloating_t
    ! -------------------------------------------------------------------------
    ! Concrete turbine type extending WindTurbine_t for the DTU 10 MW floating
    ! platform (CENER DeltaWind version). Stores the platform-specific geometry
    ! (column/pontoon cross-sections, wall thickness, densities, wetted area) and
    ! the member-ID maps used to build the representative 1D mass distributions.
    ! -------------------------------------------------------------------------

    ! Physical and structural parameters specific to the DTU 10 MW Floating CENER DeltaWind version
    REAL(WP) :: col_D    = 14.5_WP      ! [m] Columns diameter
    REAL(WP) :: XsecA    = 10.875_WP    ! [m] Pontoon width
    REAL(WP) :: XsecB    = 7.0_WP       ! [m] Pontoon height
    REAL(WP) :: rho_wat  = 1025.0_WP    ! [kg/m^3] Water density
    REAL(WP) :: rho_mat  = 7850.0_WP    ! [kg/m^3] Material density
    REAL(WP) :: t        = 0.023_WP     ! [m] Wall thickness
    REAL(WP) :: wet_area = 9659.07      ! [m^2] Wetted area

    ! [-] Column members ID list: [[0,1],[2,3],[4,5]]
    INTEGER(I32) :: col_members(3,2) = reshape([0_I32, 2_I32, 4_I32, 1_I32, 3_I32, 5_I32], [3, 2])
    ! [-] Pontoon members ID lists: [[6],[7],[8]]
    INTEGER(I32) :: pon_members(3,1) = reshape([6_I32, 7_I32, 8_I32], [3, 1])

    ! Built and stored by compute_force for later reuse by get_impedance_corrected_force
    LOGICAL, ALLOCATABLE :: keep_flat(:)    ! [-] Flat (Nmembers*Nnodes) wet and non-duplicate node keep-mask

    CONTAINS

    ! --- Shared procedures --- !
    PROCEDURE :: init => init_DTU10MWFloating

    ! --- Deferred (abstract) procedures --- !
    PROCEDURE :: compute_force                 => compute_force_DTU10MWFloating
    PROCEDURE :: compute_mass                  => compute_mass_DTU10MWFloating
    PROCEDURE :: filter_frequencies            => filter_frequencies_DTU10MWFloating
    PROCEDURE :: save_specific_params          => save_specific_params_DTU10MWFloating
    PROCEDURE :: write_specific_params_txt     => write_specific_params_txt_DTU10MWFloating
    PROCEDURE :: get_impedance_corrected_force => get_impedance_corrected_force_DTU10MWFloating

END TYPE DTU10MWFloating_t


CONTAINS


SUBROUTINE init_DTU10MWFloating(self, config)
    ! -------------------------------------------------------------------------
    ! Platform-specific constructor: resolves the (hardcoded) path to the DTU
    ! 10 MW rpm-vs-wind-speed curve, flags this turbine as non-monopile, then
    ! delegates to the shared base-class init for the common parameters.
    ! -------------------------------------------------------------------------

    CLASS(DTU10MWFloating_t)  , INTENT(INOUT) :: self      ! [-] Turbine instance
    TYPE(Config_WindTurbine_t), INTENT(IN)    :: config    ! [-] Turbine configuration

    ! Local variables
    LOGICAL :: exists


    ! Define path to rpm curve (hardcoded)
    self % path_rpm = "../wind_speed_curves_DTU_10MW/rpm_ws.csv"

    ! Verify that the rpm curve file exists
    inquire(file=self % path_rpm, exist=exists)
    if (.not. exists) error stop "DTU10MWFloating % init: RPM curve file not found: " // self % path_rpm

    ! Set the case type
    self % is_monopile   = .false.

    CALL init(self, config)

END SUBROUTINE init_DTU10MWFloating


SUBROUTINE compute_force_DTU10MWFloating(self, config, verbose)
    ! -------------------------------------------------------------------------
    ! Computes the structural nodal forces from the OpenFAST accelerations:
    ! flattens and RFFTs the acceleration time series, combines it with the
    ! per-node dry+added mass (F = -A*m_eff), optionally filters the result to
    ! the frequencies actually excited by the drivetrain, and finally strips
    ! out dry and duplicated nodes so only the wetted structure remains, ready
    ! to be treated as an array of point dipoles.
    ! -------------------------------------------------------------------------

    USE MathUtils, ONLY: compute_rfft

    CLASS(DTU10MWFloating_t)  , INTENT(INOUT) :: self       ! [-] Turbine instance
    TYPE(Config_WindTurbine_t), INTENT(IN)    :: config     ! [-] Turbine configuration
    LOGICAL                   , INTENT(IN)    :: verbose    ! [-] Whether to print a summary

    ! Local variables
    INTEGER(I32)              :: Nm, Nn, nt, nf, i, j, k, v, Nnodes_wet, Nfreqs_new
    INTEGER(I32)              :: m1, global_idx
    REAL(WP)                  :: dt, F_memory_mb
    LOGICAL                   :: verbose_, filter_freqs_
    INTEGER(I32)              :: skipf_
    REAL(WP)    , ALLOCATABLE :: acc_flat(:,:,:), mass(:), added_mass(:)
    REAL(WP)    , ALLOCATABLE :: freqs(:), mass_eff(:,:), freqs_new(:)
    COMPLEX(WP) , ALLOCATABLE :: A_flat(:,:,:), F_temp(:,:,:,:)
    COMPLEX(WP) , ALLOCATABLE :: A(:,:,:,:), F_filtered(:,:,:,:)
    LOGICAL     , ALLOCATABLE :: mask_freqs_to_use(:), keep(:,:), keep_flat(:)
    INTEGER(I32), ALLOCATABLE :: valid_node(:), valid_member(:)


    ! Use the parsed configuration for runtime options
    filter_freqs_ = config % filter_freqs
    verbose_      = verbose
    skipf_        = config % skipf

    Nm = self % Nmembers
    Nn = self % Nnodes
    nt = size(self % Time)
    dt = self % Time(2) - self % Time(1)

    ! Acceleration flatten
    ALLOCATE(acc_flat(nt, Nm*Nn, 3))
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i)
    do j = 1, Nm
        do k = 1, Nn
            global_idx = (j - 1_I32) * Nn + k
            do i = 1, nt
                acc_flat(i, global_idx, 1) = self % acc(i, j, k, 1)
                acc_flat(i, global_idx, 2) = self % acc(i, j, k, 2)
                acc_flat(i, global_idx, 3) = self % acc(i, j, k, 3)
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    ! Compute acceleration RFFT
    CALL compute_rfft(acc_flat, nt, dt, skipf=skipf_, remove_zero=.true., &
                      array_out=A_flat, freqs=freqs)
    self % Freqs = freqs
    nf           = size(freqs)

    ! Recover acceleration in shape (nf, Nm, Nn, 3)
    ALLOCATE(A(nf, Nm, Nn, 3))
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i)
    do j = 1, Nm
        do k = 1, Nn
            do i = 1, nf
                A(i, j, k, 1) = A_flat(i, (j-1)*Nn+k, 1)
                A(i, j, k, 2) = A_flat(i, (j-1)*Nn+k, 2)
                A(i, j, k, 3) = A_flat(i, (j-1)*Nn+k, 3)
            end do
        end do
    end do
    !$OMP END PARALLEL DO
    DEALLOCATE(acc_flat, A_flat)

    ! Compute mass properties
    CALL self % compute_mass(mass, added_mass)

    ! Convert the 1D arrays into a 2D matrix (Nm, Nn) matching C-order layout
    ALLOCATE(mass_eff(Nm, Nn))
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j)
    do i = 1, Nm
        do j = 1, Nn
            global_idx     = (i - 1_I32) * Nn + j
            mass_eff(i, j) = mass(global_idx) + added_mass(global_idx)
        end do
    end do
    !$OMP END PARALLEL DO

    DEALLOCATE(mass, added_mass)

    ! Compute force via F = -A * mass_eff
    ALLOCATE(F_temp(nf, Nm, Nn, 3))

    ! Parallelized over space and coordinate directions
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i, j, k)
    do k = 1, 3
        do j = 1, Nn
            do i = 1, Nm
                F_temp(:, i, j, k) = -A(:, i, j, k) * mass_eff(i, j)
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    DEALLOCATE(A, mass_eff)

    ! Filter frequencies
    if (filter_freqs_) then
        mask_freqs_to_use = self % filter_frequencies(config % freqs_over)
        Nfreqs_new        = count(mask_freqs_to_use)

        ALLOCATE(freqs_new(Nfreqs_new))
        ALLOCATE(F_filtered(Nfreqs_new, Nm, Nn, 3))

        v = 1_I32
        do i = 1, nf
            if (mask_freqs_to_use(i)) then
                freqs_new(v)           = self % Freqs(i)
                F_filtered(v, :, :, :) = F_temp(i, :, :, :)
                v                      = v + 1_I32
            end if
        end do

        CALL MOVE_ALLOC(freqs_new, self % Freqs)
        CALL MOVE_ALLOC(F_filtered, F_temp)
        nf = Nfreqs_new
        DEALLOCATE(mask_freqs_to_use)
    end if

    ! Deallocate time and raw accelerations
    if (ALLOCATED(self % Time)) DEALLOCATE(self % Time)
    if (ALLOCATED(self % acc))  DEALLOCATE(self % acc)

    ! Remove duplicate and dry nodes
    ALLOCATE(keep(Nm, Nn))
    keep = .true.

    ! Remove the first node of the second member of each column
    do i = 1, 3
        m1          = self % col_members(i, 2) + 1_I32    ! Convert to 1-based indexing
        keep(m1, 1) = .false.
    end do

    ! Remove dry nodes
    do j = 1, Nn
        do i = 1, Nm
            if (self % xyz_all(i, j, 3) > 0.0_WP) then
                keep(i, j) = .false.
            end if
        end do
    end do

    ! Flatten mask
    ALLOCATE(keep_flat(Nm*Nn))
    do i = 1, Nm
        do j = 1, Nn
            keep_flat((i-1)*Nn+j) = keep(i, j)
        end do
    end do
    DEALLOCATE(keep)

    Nnodes_wet = count(keep_flat)

    ! Build a thread-safe index map for retrieving valid wetted nodes
    ALLOCATE(valid_member(Nnodes_wet), valid_node(Nnodes_wet))
    v = 1_I32
    do i = 1, Nm
        do j = 1, Nn
            if (keep_flat((i-1)*Nn+j)) then
                valid_member(v) = i
                valid_node(v)   = j
                v               = v + 1_I32
            end if
        end do
    end do

    ! Extract coordinates of submerged/non-duplicate nodes: self % xyz (Nnodes_wet, 3)
    ALLOCATE(self % xyz(Nnodes_wet, 3))
    !$OMP PARALLEL DO PRIVATE(v, i, j)
    do v = 1, Nnodes_wet
        i                = valid_member(v)
        j                = valid_node(v)
        self % xyz(v, :) = self % xyz_all(i, j, :)
    end do
    !$OMP END PARALLEL DO

    ! Extract corresponding forces: self % Force (nf, Nnodes_wet, 3)
    ALLOCATE(self % Force(nf, Nnodes_wet, 3))
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(v, i, j, k)
    do k = 1, 3
        do v = 1, Nnodes_wet
            i                     = valid_member(v)
            j                     = valid_node(v)
            self % Force(:, v, k) = F_temp(:, i, j, k)
        end do
    end do
    !$OMP END PARALLEL DO

    ! Store the keep_flat array on the turbine instance for future impedance correction
    if (ALLOCATED(self % keep_flat)) DEALLOCATE(self % keep_flat)
    CALL MOVE_ALLOC(keep_flat, self % keep_flat)

    DEALLOCATE(F_temp, valid_member, valid_node)

    ! --- Verbose summary ---
    if (verbose_) then
        F_memory_mb = real(storage_size(self % Force) * size(self % Force) / 8, WP) / 1024.0_WP / 1024.0_WP
        write(*, '(A)') repeat('=', 70)
        write(*, '(A)') 'FORCE COMPUTATION COMPLETED (FLOATING PLATFORM)'
        write(*, '(A)') repeat('-', 68)
        write(*, '(A, I6)') 'Nodes (wet): ', Nnodes_wet
        write(*, '(A, I6)') 'Frequencies: ', nf
        write(*, '(A, F6.2, A)') 'Memory (F):  ', F_memory_mb, ' MB'
        write(*, '(A)') repeat('=', 70)
    end if

END SUBROUTINE compute_force_DTU10MWFloating


SUBROUTINE compute_mass_DTU10MWFloating(self, mass, added_mass)
    ! -------------------------------------------------------------------------
    ! Builds the per-node dry (structural) and added (hydrodynamic) mass for
    ! the whole platform. Rather than treating every member independently, it
    ! collapses each column/pontoon family into a single representative 1D
    ! line of nodes, computes the mass from its cylindrical/rectangular
    ! cross-section along that line, and then broadcasts the result back onto
    ! the full (Nmembers, Nnodes) grid, since all columns/pontoons share the
    ! same geometry.
    ! -------------------------------------------------------------------------

    USE MathUtils, ONLY: divide_span

    CLASS(DTU10MWFloating_t), INTENT(INOUT)            :: self             ! [-] Turbine instance
    REAL(WP)                , INTENT(OUT), ALLOCATABLE :: mass(:)          ! [kg] Dry mass per node (Nm*Nn)
    REAL(WP)                , INTENT(OUT), ALLOCATABLE :: added_mass(:)    ! [kg] Added mass per node (Nm*Nn)

    ! Local variables
    INTEGER(I32) :: Nm, Nn, Nnodes_total, i, j, col_m0, col_m1, pon_m, global_idx
    REAL(WP)     :: length, R_outer_col, R_inner_col
    REAL(WP)     :: A_struct_col, A_added_col, A_struct_pon, A_added_pon

    ! Representative geometry arrays
    REAL(WP), ALLOCATABLE :: col_pos(:,:)         ! Coordinates of full column (2*Nn - 1, 3)
    REAL(WP), ALLOCATABLE :: nodes_col_1d(:)      ! 1D coordinates along column axis (2*Nn - 1)
    REAL(WP), ALLOCATABLE :: ds_col(:)            ! Spacing intervals for columns (2*Nn - 1)
    REAL(WP), ALLOCATABLE :: mass_col(:)          ! Dry mass along representative column (2*Nn - 1)
    REAL(WP), ALLOCATABLE :: added_mass_col(:)    ! Fluid added mass along column (2*Nn - 1)

    REAL(WP), ALLOCATABLE :: arc_length(:)        ! Cumulative arc-length along pontoon (Nn)
    REAL(WP), ALLOCATABLE :: ds_pon(:)            ! Spacing intervals for pontoons (Nn)
    REAL(WP), ALLOCATABLE :: mass_pon(:)          ! Dry mass along representative pontoon (Nn)
    REAL(WP), ALLOCATABLE :: added_mass_pon(:)    ! Fluid added mass along pontoon (Nn)

    ! Intermediate 2D distribution arrays (Nm, Nn)
    REAL(WP), ALLOCATABLE :: mass_2D(:,:)
    REAL(WP), ALLOCATABLE :: added_mass_2D(:,:)


    Nm           = self % Nmembers
    Nn           = self % Nnodes
    Nnodes_total = Nm * Nn

    ! Retrieve representative member IDs (converting from 0-based to 1-based indexing)
    col_m0 = self % col_members(1, 1) + 1_I32
    col_m1 = self % col_members(1, 2) + 1_I32
    pon_m  = self % pon_members(1, 1) + 1_I32

    ! --- 1. Representative column calculations ---
    ! Concatenate column node coordinates dropping the duplicate shared node
    ALLOCATE(col_pos(2*Nn-1, 3))
    col_pos(1:Nn, :)        = self % xyz_all(col_m0, :, :)
    col_pos(Nn+1:2*Nn-1, :) = self % xyz_all(col_m1, 2:Nn, :)

    ! Compute physical column height
    length = norm2(col_pos(2*Nn-1, :) - col_pos(1, :))

    ! Generate linearly spaced 1D coordinates from 0 to total length
    ALLOCATE(nodes_col_1d(2*Nn-1))
    do i = 1, 2 * Nn - 1
        nodes_col_1d(i) = (real(i - 1_I32, WP) * length) / real(2_I32 * Nn - 2_I32, WP)
    end do

    ! Calculate spacing intervals
    ds_col = divide_span(nodes_col_1d)

    ! Cylindrical cross-section areas
    R_outer_col  = self % col_D / 2.0_WP
    R_inner_col  = R_outer_col - self % t
    A_struct_col = PI * (R_outer_col**2 - R_inner_col**2)
    A_added_col  = PI * R_outer_col**2

    ! Distributed structural mass and added mass for columns
    ALLOCATE(mass_col(2*Nn-1), added_mass_col(2*Nn-1))
    mass_col       = A_struct_col * self % rho_mat * ds_col
    added_mass_col = A_added_col * self % rho_wat * ds_col

    ! Add localized top/bottom steel closure lids mass (plate volume * density)
    mass_col(1)      = mass_col(1) + (PI * R_outer_col**2 * self % rho_mat * self % t)
    mass_col(2*Nn-1) = mass_col(2*Nn-1) + (PI * R_outer_col**2 * self % rho_mat * self % t)

    DEALLOCATE(col_pos, nodes_col_1d, ds_col)

    ! --- 2. Representative pontoon calculations ---
    ! Calculate cumulative horizontal arc-length along the pontoon axis
    ALLOCATE(arc_length(Nn))
    arc_length(1) = 0.0_WP
    do i = 2, Nn
        arc_length(i) = arc_length(i-1) + norm2(self % xyz_all(pon_m, i, :) - self % xyz_all(pon_m, i-1, :))
    end do

    ! Calculate spacing intervals
    ds_pon = divide_span(arc_length)

    ! Hollow square cross-section areas
    A_struct_pon = (self % XsecA * self % XsecB) - &
                   ((self % XsecA - 2.0_WP * self % t) * (self % XsecB - 2.0_WP * self % t))
    A_added_pon  = self % XsecA * self % XsecB

    ! Distributed structural mass and added mass for pontoons
    ALLOCATE(mass_pon(Nn), added_mass_pon(Nn))
    mass_pon       = A_struct_pon * self % rho_mat * ds_pon
    added_mass_pon = A_added_pon * self % rho_wat * ds_pon

    DEALLOCATE(arc_length, ds_pon)

    ! --- 3. Assemble structural mass grid and flatten ---
    ALLOCATE(mass_2D(Nm, Nn), added_mass_2D(Nm, Nn))
    mass_2D       = 0.0_WP
    added_mass_2D = 0.0_WP

    ! Broadcast columns (mapping split spans contiguously)
    do i = 1, 3
        col_m0 = self % col_members(i, 1) + 1_I32
        col_m1 = self % col_members(i, 2) + 1_I32

        mass_2D(col_m0, :) = mass_col(1:Nn)
        mass_2D(col_m1, :) = mass_col(Nn:2*Nn-1)

        added_mass_2D(col_m0, :) = added_mass_col(1:Nn)
        added_mass_2D(col_m1, :) = added_mass_col(Nn:2*Nn-1)
    end do

    ! Broadcast pontoons
    do i = 1, 3
        pon_m = self % pon_members(i, 1) + 1_I32

        mass_2D(pon_m, :)       = mass_pon(:)
        added_mass_2D(pon_m, :) = added_mass_pon(:)
    end do

    DEALLOCATE(mass_col, added_mass_col, mass_pon, added_mass_pon)

    ! Allocate outputs of size Nnodes_total (Nm * Nn)
    ALLOCATE(mass(Nnodes_total), added_mass(Nnodes_total))

    ! Flatten in C-order matching Python's flatten exactly
    do i = 1, Nm
        do j = 1, Nn
            global_idx             = (i - 1_I32) * Nn + j
            mass(global_idx)       = mass_2D(i, j)
            added_mass(global_idx) = added_mass_2D(i, j)
        end do
    end do

    DEALLOCATE(mass_2D, added_mass_2D)

END SUBROUTINE compute_mass_DTU10MWFloating


FUNCTION filter_frequencies_DTU10MWFloating(self, freqs_over) RESULT(mask)
    ! -------------------------------------------------------------------------
    ! Builds a mask that keeps only the computed frequencies lying close to the
    ! DTU 10 MW drivetrain excitation peaks (shaft/gear-mesh harmonics at the
    ! rpm corresponding to the current wind speed), so downstream force arrays
    ! can be trimmed to the physically relevant spectral content.
    ! -------------------------------------------------------------------------

    USE IOUtils, ONLY: read_curve
    USE MathUtils, ONLY: generate_timeseries_banded_sines, filter_non_useful_freqs

    CLASS(DTU10MWFloating_t), INTENT(INOUT)        :: self          ! [-] Turbine instance
    REAL(WP)                , INTENT(IN), OPTIONAL :: freqs_over    ! [Hz] Frequencies below this are never filtered
    LOGICAL                 , ALLOCATABLE          :: mask(:)       ! [-] Frequency mask (Nfreqs), .true. to keep

    ! Local variables
    CHARACTER(len=20), ALLOCATABLE :: keys(:)
    REAL(WP)         , ALLOCATABLE :: freqs_amp(:,:)
    REAL(WP)         , ALLOCATABLE :: freqs_to_use(:)
    REAL(WP)                       :: rpm, freqs_over_

    REAL(WP), ALLOCATABLE :: tmp_signal(:)


    freqs_over_ = 10.0_WP ; if (present(freqs_over)) freqs_over_ = freqs_over

    ! Read rpm at current wind speed
    rpm = read_curve(self % path_rpm, self % WindSpeed)

    ! Generate drivetrain excitation spectrum
    CALL drivetrain10MW_excitation_spectrum(rpm=rpm, freqs_amp=freqs_amp, keys=keys)

    ! Get allowed frequencies
    CALL generate_timeseries_banded_sines(peaks=freqs_amp, keys=keys, t=self % Time, &
                                          zeta=0.02_WP, used_freqs=.true., a=tmp_signal, &
                                          freqs_out=freqs_to_use)

    ! Build mask
    mask = filter_non_useful_freqs(self % Freqs, freqs_to_use, freqs_over_)

    DEALLOCATE(tmp_signal)

    DEALLOCATE(freqs_amp, keys, freqs_to_use)

END FUNCTION filter_frequencies_DTU10MWFloating


SUBROUTINE save_specific_params_DTU10MWFloating(self, save_path, istr)
    ! -------------------------------------------------------------------------
    ! Persists the platform-specific geometry and material parameters (column
    ! diameter, pontoon cross-section, wall thickness, densities, wetted area,
    ! member ID maps) to the solver's HDF5 output, optionally suffixing each
    ! variable name with an instance tag (istr) when several turbines are saved.
    ! -------------------------------------------------------------------------

    USE IOUtils, ONLY: save_variable

    CLASS(DTU10MWFloating_t), INTENT(IN)           :: self         ! [-] Turbine instance
    CHARACTER(len=*)        , INTENT(IN)           :: save_path    ! [-] Path to the HDF5 output file
    CHARACTER(len=*)        , INTENT(IN), OPTIONAL :: istr         ! [-] Suffix appended to the variable names


    if (present(istr)) then
        CALL save_variable(save_path, "Case_type_"       // istr, self % case_type)
        CALL save_variable(save_path, "Structure_nodes_" // istr, self % xyz_all)
        CALL save_variable(save_path, "col_D_"           // istr, self % col_D)
        CALL save_variable(save_path, "XsecA_"           // istr, self % XsecA)
        CALL save_variable(save_path, "XsecB_"           // istr, self % XsecB)
        CALL save_variable(save_path, "t_"               // istr, self % t)
        CALL save_variable(save_path, "rho_mat_"         // istr, self % rho_mat)
        CALL save_variable(save_path, "wet_area_"        // istr, self % wet_area)
        CALL save_variable(save_path, "col_members_"     // istr, self % col_members)
        CALL save_variable(save_path, "pon_members_"     // istr, self % pon_members)
    else
        CALL save_variable(save_path, "Case_type"      , self % case_type)
        CALL save_variable(save_path, "Structure_nodes", self % xyz_all)
        CALL save_variable(save_path, "col_D"          , self % col_D)
        CALL save_variable(save_path, "XsecA"          , self % XsecA)
        CALL save_variable(save_path, "XsecB"          , self % XsecB)
        CALL save_variable(save_path, "t"              , self % t)
        CALL save_variable(save_path, "rho_mat"        , self % rho_mat)
        CALL save_variable(save_path, "wet_area"       , self % wet_area)
        CALL save_variable(save_path, "col_members"    , self % col_members)
        CALL save_variable(save_path, "pon_members"    , self % pon_members)
    end if

END SUBROUTINE save_specific_params_DTU10MWFloating


SUBROUTINE write_specific_params_txt_DTU10MWFloating(self, unit)
    ! -------------------------------------------------------------------------
    ! Writes a short, human-readable summary of the platform-specific geometry
    ! and member maps to the given text unit, mirroring what
    ! save_specific_params_DTU10MWFloating stores in the HDF5 file.
    ! -------------------------------------------------------------------------

    CLASS(DTU10MWFloating_t), INTENT(IN) :: self    ! [-] Turbine instance
    INTEGER(I32)            , INTENT(IN) :: unit    ! [-] Unit of the open .dat file


    write(unit, '(A,F8.3,A)')   'Column diameter  : ', self % col_D, ' m'
    write(unit, '(A,F8.3,A)')   'Pontoon width    : ', self % XsecA, ' m'
    write(unit, '(A,F8.3,A)')   'Pontoon height   : ', self % XsecB, ' m'
    write(unit, '(A,F8.3,A)')   'Wall thickness   : ', self % t, ' m'
    write(unit, '(A,F8.1,A)')   'Water density    : ', self % rho_wat, ' kg/m^3'
    write(unit, '(A,F8.1,A)')   'Material density : ', self % rho_mat, ' kg/m^3'
    write(unit, '(A,F10.2,A)')  'Wetted area      : ', self % wet_area, ' m^2'
    write(unit, '(A,3(I0,1X))') 'Column members (top)   : ', self % col_members(:, 1)
    write(unit, '(A,3(I0,1X))') 'Column members (bottom): ', self % col_members(:, 2)
    write(unit, '(A,3(I0,1X))') 'Pontoon members        : ', self % pon_members(:, 1)

END SUBROUTINE write_specific_params_txt_DTU10MWFloating


FUNCTION get_impedance_corrected_force_DTU10MWFloating(self, c) RESULT(corrected_F)
    ! -------------------------------------------------------------------------
    ! Applies a cylindrical-scattering (Hankel-function based) shielding
    ! correction to the nodal forces, per member family: columns and pontoons
    ! are each treated as an infinite cylinder of their own (equivalent)
    ! diameter, and the resulting frequency-dependent correction factor scales
    ! the raw force so that a solver using the Image Method sees the acoustic
    ! shadowing effect of the structure itself.
    ! -------------------------------------------------------------------------

    USE MathUtils, ONLY: alpha_hankel

    CLASS(DTU10MWFloating_t), INTENT(IN)  :: self                  ! [-] Turbine instance
    REAL(WP)                , INTENT(IN)  :: c                     ! [m/s] Speed of sound in the fluid
    COMPLEX(WP)             , ALLOCATABLE :: corrected_F(:,:,:)    ! [N] Corrected force (Nfreqs, Nnodes_wet, 3)

    ! Local variables
    INTEGER(I32)              :: Nm, Nn, Nnodes_wet, nf, c_idx, p_idx, m, j, i, m0, m1, global_idx, idx
    INTEGER(I32), ALLOCATABLE :: original_to_wet(:)
    COMPLEX(WP) , ALLOCATABLE :: alpha_full(:,:)
    REAL(WP)    , ALLOCATABLE :: k(:)
    COMPLEX(WP) , ALLOCATABLE :: alpha_col(:), alpha_pon(:)
    REAL(WP)                  :: D_pon


    Nm         = self % Nmembers
    Nn         = self % Nnodes
    Nnodes_wet = size(self % xyz, 1)
    nf         = size(self % Freqs)

    ! --- 1. Map original flat node indices to wetted nodes (-1 for dry/removed nodes) ---
    ALLOCATE(original_to_wet(Nm*Nn))
    original_to_wet = -1_I32
    idx             = 1_I32
    do i = 1, Nm * Nn
        if (self % keep_flat(i)) then
            original_to_wet(i) = idx
            idx                = idx + 1_I32
        end if
    end do

    ! --- 2. Initialize the full scattering correction array (nf, Nnodes_wet) with 1.0 ---
    ALLOCATE(alpha_full(nf, Nnodes_wet))
    alpha_full = (1.0_WP, 0.0_WP)

    ! --- 3. Compute wavenumber vector k ---
    ALLOCATE(k(nf))
    k = (2.0_WP * PI * self % Freqs) / c

    ! --- 4. Columns shielding correction (cylindrical Hankel factors) ---
    alpha_col = alpha_hankel(k, self % col_D)    ! Shape (nf)

    do c_idx = 1, 3
        m0 = self % col_members(c_idx, 1) + 1_I32    ! Convert from 0-based to 1-based
        m1 = self % col_members(c_idx, 2) + 1_I32

        ! First column member (all nodes)
        do j = 1, Nn
            global_idx = (m0 - 1_I32) * Nn + j
            idx        = original_to_wet(global_idx)
            if (idx >= 1_I32) then
                alpha_full(:, idx) = alpha_col(:)
            end if
        end do

        ! Second column member (skip first node as it is duplicate/removed)
        do j = 2, Nn
            global_idx = (m1 - 1_I32) * Nn + j
            idx        = original_to_wet(global_idx)
            if (idx >= 1_I32) then
                alpha_full(:, idx) = alpha_col(:)
            end if
        end do
    end do
    DEALLOCATE(alpha_col)

    ! --- 5. Pontoons shielding correction (equivalent diameter Hankel factors) ---
    ! Equivalent circular diameter based on square cross-section area
    D_pon     = 2.0_WP * sqrt((self % XsecA * self % XsecB) / PI)
    alpha_pon = alpha_hankel(k, D_pon)    ! Shape (nf)

    do p_idx = 1, 3
        m = self % pon_members(p_idx, 1) + 1_I32    ! Convert from 0-based to 1-based

        do j = 1, Nn
            global_idx = (m - 1_I32) * Nn + j
            idx        = original_to_wet(global_idx)
            if (idx >= 1_I32) then
                alpha_full(:, idx) = alpha_pon(:)
            end if
        end do
    end do
    DEALLOCATE(alpha_pon)

    ! --- 6. Apply correction: F_corrected = F * abs(alpha_full) ---
    ALLOCATE(corrected_F(nf, Nnodes_wet, 3))

    ! High-performance vectorized multiplication over spatial components
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i, j, m)
    do m = 1, 3
        do j = 1, Nnodes_wet
            do i = 1, nf
                corrected_F(i, j, m) = self % Force(i, j, m) * abs(alpha_full(i, j))
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    ! Cleanup
    DEALLOCATE(original_to_wet, alpha_full, k)

END FUNCTION get_impedance_corrected_force_DTU10MWFloating

END MODULE DTU10MWFloating
