! ===============================================================================
! MODULE: SAITEC2MWFloating
! ===============================================================================
! Defines SAITEC2MWFloating_t, a concrete turbine type extending the abstract
! WindTurbine_t (src/core/WindTurbine.f90) for the SAITEC SENVION 2 MW floating
! platform (SATH dual-hull concrete technology). It implements the deferred
! procedures: compute_force, compute_mass, filter_frequencies,
! save/write_specific_params and get_impedance_corrected_force, treating the
! twin horizontal elliptical concrete hulls as the only structural members
! that generate acoustic point-dipole sources; the virtual rigid-link "joint"
! members that merely connect the two hulls to the tower base have no
! interior structural geometry of their own and are excluded entirely.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE SAITEC2MWFloating

USE Kinds
USE Types, ONLY: Config_WindTurbine_t
USE WindTurbine, ONLY: WindTurbine_t, init
USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_nan

IMPLICIT NONE

PRIVATE
PUBLIC :: SAITEC2MWFloating_t


TYPE, EXTENDS(WindTurbine_t) :: SAITEC2MWFloating_t
    ! -------------------------------------------------------------------------
    ! Concrete turbine type extending WindTurbine_t for the SAITEC SENVION 2 MW
    ! floating platform. Stores the platform-specific elliptical hull geometry and
    ! the member-ID maps used to build the representative 1D mass distribution and
    ! to identify the virtual joint members excluded from the acoustic source set.
    ! -------------------------------------------------------------------------

    ! Physical and structural parameters specific to the SAITEC SENVION 2 MW SATH platform
    REAL(WP) :: Diam_y   = 10.0_WP      ! [m] Elliptical floater diameter in y-direction
    REAL(WP) :: Diam_z   = 12.5_WP      ! [m] Elliptical floater diameter in z-direction
    REAL(WP) :: rho_wat  = 1025.0_WP    ! [kg/m^3] Water density
    REAL(WP) :: rho_mat  = 1850.0_WP    ! [kg/m^3] Material density (concrete)
    REAL(WP) :: t        = 1.0_WP       ! [m] Wall thickness
    REAL(WP) :: wet_area = 0.0_WP       ! [m^2] Wetted area (not computed for this platform)

    ! [-] Floater members ID list: [[0,1],[3,4]]
    INTEGER(I32) :: float_members(2,2) = reshape([0_I32, 3_I32, 1_I32, 4_I32], [2, 2])
    ! [-] Virtual rigid-link (joint) member IDs, excluded
    INTEGER(I32) :: joint_members(3) = [5_I32, 2_I32, 6_I32]

    ! Hardcoded paths to the SAITEC 2 MW drivetrain bar-chart amplitude data
    CHARACTER(len=:), ALLOCATABLE :: path_barRB    ! [-] Rotor-bearing bar-chart data file
    CHARACTER(len=:), ALLOCATABLE :: path_barGB    ! [-] Gearbox bar-chart data file

    CONTAINS

    ! --- Shared procedures --- !
    PROCEDURE :: init       => init_SAITEC2MWFloating
    PROCEDURE :: read_input => read_input_SAITEC2MWFloating

    ! --- Deferred (abstract) procedures --- !
    PROCEDURE :: compute_force                 => compute_force_SAITEC2MWFloating
    PROCEDURE :: compute_mass                  => compute_mass_SAITEC2MWFloating
    PROCEDURE :: filter_frequencies            => filter_frequencies_SAITEC2MWFloating
    PROCEDURE :: save_specific_params          => save_specific_params_SAITEC2MWFloating
    PROCEDURE :: write_specific_params_txt     => write_specific_params_txt_SAITEC2MWFloating
    PROCEDURE :: get_impedance_corrected_force => get_impedance_corrected_force_SAITEC2MWFloating

END TYPE SAITEC2MWFloating_t


CONTAINS


SUBROUTINE init_SAITEC2MWFloating(self, config)
    ! -------------------------------------------------------------------------
    ! Platform-specific constructor: resolves the (hardcoded) paths to the
    ! SAITEC 2 MW rpm-vs-wind-speed curve and drivetrain bar-chart amplitude
    ! data, flags this turbine as non-monopile, then delegates to the shared
    ! base-class init for the common parameters.
    ! -------------------------------------------------------------------------

    CLASS(SAITEC2MWFloating_t), INTENT(INOUT) :: self      ! [-] Turbine instance
    TYPE(Config_WindTurbine_t), INTENT(INOUT) :: config    ! [-] Turbine configuration

    ! Local variables
    LOGICAL :: exists = .false.


    ! Define paths to rpm curve and drivetrain bar-chart data (hardcoded)
    self % path_rpm   = "../wind_speed_curves_SAITEC_2MW/rpm_ws.csv"
    self % path_barRB = "../wind_speed_curves_SAITEC_2MW/Bar_RB12.csv"
    self % path_barGB = "../wind_speed_curves_SAITEC_2MW/Bar_GB123.csv"

    ! Verify that the required curve/data files exist
    inquire(file=self % path_rpm, exist=exists)
    if (.not. exists) error stop "SAITEC2MWFloating % init: RPM curve file not found: " // self % path_rpm
    inquire(file=self % path_barRB, exist=exists)
    if (.not. exists) error stop "SAITEC2MWFloating % init: Bar_RB12 data file not found: " // self % path_barRB
    inquire(file=self % path_barGB, exist=exists)
    if (.not. exists) error stop "SAITEC2MWFloating % init: Bar_GB123 data file not found: " // self % path_barGB

    ! Set the case type
    self % is_monopile   = .false.
    config % is_monopile = .false.

    CALL init(self, config)

    write(*,*) "SAITEC2MWFloating % init(): CHECK IF BARIPOS IS WELL COMPUTED, SHOULD BE (-8.9586, 0)"

END SUBROUTINE init_SAITEC2MWFloating


SUBROUTINE read_input_SAITEC2MWFloating(self, config, verbose)
    ! -------------------------------------------------------------------------
    ! Overrides the shared WindTurbine % read_input because, for this platform,
    ! the OpenFAST/SubDyn output does NOT follow the "one output channel group
    ! per structural member" assumption the base implementation relies on: the
    ! three virtual rigid-link "joint" members (float-to-float and float-to-
    ! tower connections) have no interior nodes and are never requested as
    ! SubDyn member outputs, so the .outb file only contains acceleration
    ! channels for the Nm_out=4 real members (SubDyn's MOutLst: MemberID
    ! 1, 2, 4, 5, output as sequential local groups M1..M4). The static
    ! geometry (SD.sum.yaml), by contrast, does list all self % Nmembers = 7
    ! members. This override reads the geometry at full member count, reads
    ! the acceleration at the reduced (compact) member count, then scatters
    ! the compact data back onto the full (Nmembers, Nnodes, 3) layout used by
    ! compute_mass/compute_force (joint members are left at exactly zero
    ! acceleration, which is harmless since they are excluded entirely by
    ! compute_force's keep mask). It also guards the automatic barycenter
    ! computation against the NaN padding that the rigid-link members (with
    ! fewer than Nnodes real nodes) introduce into self % xyz_all.
    ! -------------------------------------------------------------------------

    USE IOUtils, ONLY: get_SDsum_variables, read_input_SD

    CLASS(SAITEC2MWFloating_t), INTENT(INOUT) :: self       ! [-] Turbine instance
    TYPE(Config_WindTurbine_t), INTENT(INOUT) :: config     ! [-] Turbine configuration
    LOGICAL                   , INTENT(IN)    :: verbose    ! [-] Whether to print a summary

    ! Local variables
    INTEGER(I32), PARAMETER :: Nm_out = 4_I32    ! [-] Members actually output by SubDyn (the float_members)
    ! [-] SubDyn MOutLst MemberID (1-based) of each local output slot
    INTEGER(I32)     , PARAMETER   :: local_to_global(Nm_out) = [1_I32, 2_I32, 4_I32, 5_I32]
    REAL(WP)                       :: dt
    INTEGER(I32)                   :: status, i, gm, nt
    CHARACTER(len=10)              :: unit
    REAL(WP)         , ALLOCATABLE :: acc_compact(:,:,:,:)


    if (verbose) then
        write(*,*) "---------------------------------------------------------"
        write(*,*) "LOADING OPENFAST DATA"
        write(*,*) "Filepath:  ", trim(self % OP_output_path)
        write(*,*) "---------------------------------------------------------"
    end if

    if (self % Nmembers <= 0_I32 .or. self % Nnodes <= 0_I32) then
        error stop "SAITEC2MWFloating % read_input: Nmembers and Nnodes must be positive integers to read SubDyn nodes"
    end if

    ! Geometry: read the FULL structural model (all Nmembers, including the
    ! rigid links) from the SD.sum.yaml summary file
    CALL get_SDsum_variables(SD_path=self % OP_SubDyn_path, Nmembers=self % Nmembers, &
                             Nnodes=self % Nnodes, verbose=.false., Nodes=self % xyz_all)

    ! Acceleration: SubDyn only requested output for the Nm_out real members;
    ! read them compactly (local slots 1..Nm_out) and scatter into the full
    ! (nt, Nmembers, Nnodes, 3) tensor at their true global member position
    CALL read_input_SD(filename=self % OP_output_path, what='acceleration', skip=config % skip, &
                       Nmembers=Nm_out, Nnodes=self % Nnodes, From=config % From, Upto=config % Upto, &
                       verbose=.false., Time_out=self % Time, Array_out=acc_compact, &
                       unit_out=unit, status=status)
    if (status /= 0_I32) error stop "SAITEC2MWFloating % read_input: failed to read acceleration channels from " &
                                    // trim(self % OP_output_path)

    nt = size(self % Time)
    ALLOCATE(self % acc(nt, self % Nmembers, self % Nnodes, 3), source=0.0_WP)
    do i = 1, Nm_out
        gm                      = local_to_global(i)
        self % acc(:, gm, :, :) = acc_compact(:, i, :, :)
    end do
    DEALLOCATE(acc_compact)

    ! Compute barycenter if not provided. A plain node-count average (as the
    ! base class would do) is wrong here: node spacing along each floater is
    ! highly non-uniform (long vs. short segment) and the rigid-link members
    ! contribute duplicated/NaN-padded nodes clustered at the central pivot.
    ! Instead, use the actual structural (dry) mass distribution as weights,
    ! which correctly reduces to the true mass centroid of the wetted
    ! floaters and is exactly zero (by construction) at the excluded joints.
    if (all(self % BariPos == [0.0_WP, 0.0_WP])) then
        block
            REAL(WP)    , ALLOCATABLE :: mass(:), added_mass(:), x_flat(:), y_flat(:)
            INTEGER(I32)              :: Nm, Nn, ii, jj, idx, fm1


            CALL self % compute_mass(mass, added_mass)

            Nm = self % Nmembers
            Nn = self % Nnodes

            ! Zero out the shared/duplicate node (first node of the second
            ! member of each floater), exactly as compute_force's keep mask
            ! does, so its mass is not double-counted in the centroid
            do ii = 1, 2
                fm1                = self % float_members(ii, 2) + 1_I32
                mass((fm1-1)*Nn+1) = 0.0_WP
            end do

            ALLOCATE(x_flat(Nm*Nn), y_flat(Nm*Nn))
            do ii = 1, Nm
                do jj = 1, Nn
                    idx = (ii - 1_I32) * Nn + jj
                    ! mass(idx) is exactly 0 wherever xyz_all may be NaN
                    ! (joint/padding slots), so those entries never enter
                    ! the weighted sum; substitute 0 there to avoid 0*NaN
                    x_flat(idx) = merge(self % xyz_all(ii, jj, 1), 0.0_WP, mass(idx) > 0.0_WP)
                    y_flat(idx) = merge(self % xyz_all(ii, jj, 2), 0.0_WP, mass(idx) > 0.0_WP)
                end do
            end do

            self % BariPos(1) = sum(mass * x_flat) / sum(mass)
            self % BariPos(2) = sum(mass * y_flat) / sum(mass)

            DEALLOCATE(mass, added_mass, x_flat, y_flat)
        end block
    end if

    ! Align with wind direction and translate to axis position
    CALL self % translate_align_with_WindDir()
    config % BariPos = self % BariPos
    config % AxisPos = self % AxisPos

    if (verbose) then
        dt = self % Time(2) - self % Time(1)
        write(*,*) ''
        write(*, '(A)') repeat('=', 70)
        write(*, '(A)') 'TURBINE: ' // trim(self % rootname)
        write(*, '(A)') repeat('-', 68)
        write(*, '(A, A)') 'Type: ', trim(self % case_type)
        write(*, '(A, F5.1, A, F5.1, A)') 'Wind Speed: ', self % WindSpeed, ' m/s   | Wind Direction: ', &
                                          self % WindDir, ' deg'
        write(*, '(A, F5.1, A, F4.1, A, F4.1, A)') 'Water Depth: ', self % Depth, &
                                                   ' m   | Position: (', self % AxisPos(1), ', ', &
                                                   self % AxisPos(2), ') m'
        write(*, '(A, I6, A, F6.4, A)') 'Time Series: ', size(self % Time), ' samples @ ', dt, ' s/sample'
        write(*, '(A)') repeat('=', 70)
        write(*,*) ">> OpenFAST data loaded successfully!"
        write(*,*)
    end if

END SUBROUTINE read_input_SAITEC2MWFloating


SUBROUTINE compute_force_SAITEC2MWFloating(self, config, verbose)
    ! -------------------------------------------------------------------------
    ! Computes the structural nodal forces from the OpenFAST accelerations:
    ! flattens and RFFTs the acceleration time series, combines it with the
    ! per-node dry+added mass of the two elliptical hulls (F = -A*m_eff),
    ! optionally filters the result to the frequencies actually excited by the
    ! drivetrain, and finally strips out dry, duplicated and virtual-joint
    ! nodes so only the wetted floater structure remains, ready to be treated
    ! as an array of point dipoles.
    ! -------------------------------------------------------------------------

    USE MathUtils, ONLY: compute_rfft

    CLASS(SAITEC2MWFloating_t), INTENT(INOUT) :: self       ! [-] Turbine instance
    TYPE(Config_WindTurbine_t), INTENT(IN)    :: config     ! [-] Turbine configuration
    LOGICAL                   , INTENT(IN)    :: verbose    ! [-] Whether to print a summary

    ! Local variables
    INTEGER(I32)              :: Nm, Nn, nt, nf, i, j, k, v, Nnodes_wet, Nfreqs_new
    INTEGER(I32)              :: m1, m, global_idx, jm
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
    skipf_        = config % skip

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

    ! Remove duplicate, dry, and virtual-joint nodes
    ALLOCATE(keep(Nm, Nn))
    keep = .true.

    ! Remove the first node of the second member of each floater (shared/duplicate node)
    do i = 1, 2
        m1          = self % float_members(i, 2) + 1_I32    ! Convert to 1-based indexing
        keep(m1, 1) = .false.
    end do

    ! Exclude the virtual rigid-link (joint) members entirely: they have no
    ! independently-meaningful interior geometry and do not radiate
    do jm = 1, size(self % joint_members)
        m          = self % joint_members(jm) + 1_I32    ! Convert to 1-based indexing
        keep(m, :) = .false.
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

    DEALLOCATE(F_temp, valid_member, valid_node, keep_flat)

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

END SUBROUTINE compute_force_SAITEC2MWFloating


SUBROUTINE compute_mass_SAITEC2MWFloating(self, mass, added_mass)
    ! -------------------------------------------------------------------------
    ! Builds the per-node dry (structural) and added (hydrodynamic) mass for
    ! the two SATH concrete hulls. Collapses the representative floater
    ! (members float_members(1,1) and float_members(1,2), stitched together
    ! dropping the shared node) into a single 1D line, computes the submerged
    ! elliptical shell mass along that line at the hull's draft, and then
    ! broadcasts the result back onto both floaters (since they share the same
    ! geometry). All other members (tower/joints) receive zero mass here and
    ! are excluded later by compute_force's keep mask.
    ! -------------------------------------------------------------------------

    USE MathUtils, ONLY: divide_span, integral_elipse_z

    CLASS(SAITEC2MWFloating_t), INTENT(INOUT)            :: self             ! [-] Turbine instance
    REAL(WP)                  , INTENT(OUT), ALLOCATABLE :: mass(:)          ! [kg] Dry mass per node (Nm*Nn)
    REAL(WP)                  , INTENT(OUT), ALLOCATABLE :: added_mass(:)    ! [kg] Added mass per node (Nm*Nn)

    ! Local variables
    INTEGER(I32) :: Nm, Nn, Nnodes_total, i, j, float_m0, float_m1, global_idx
    REAL(WP)     :: z_draft, A_outer, A_inner, A_shell

    ! Representative geometry arrays
    REAL(WP), ALLOCATABLE :: float_pos(:,:)         ! Coordinates of full floater (2*Nn - 1, 3)
    REAL(WP), ALLOCATABLE :: s(:)                   ! Cumulative arc-length along the floater (2*Nn - 1)
    REAL(WP), ALLOCATABLE :: ds_float(:)            ! Spacing intervals (2*Nn - 1)
    REAL(WP), ALLOCATABLE :: mass_float(:)          ! Dry mass along representative floater (2*Nn - 1)
    REAL(WP), ALLOCATABLE :: added_mass_float(:)    ! Fluid added mass along floater (2*Nn - 1)

    ! Intermediate 2D distribution arrays (Nm, Nn)
    REAL(WP), ALLOCATABLE :: mass_2D(:,:)
    REAL(WP), ALLOCATABLE :: added_mass_2D(:,:)


    Nm           = self % Nmembers
    Nn           = self % Nnodes
    Nnodes_total = Nm * Nn

    ! Retrieve representative member IDs (converting from 0-based to 1-based indexing)
    float_m0 = self % float_members(1, 1) + 1_I32
    float_m1 = self % float_members(1, 2) + 1_I32

    ! --- 1. Representative floater calculations ---
    ! Concatenate floater node coordinates dropping the duplicate shared node
    ALLOCATE(float_pos(2*Nn-1, 3))
    float_pos(1:Nn, :)        = self % xyz_all(float_m0, :, :)
    float_pos(Nn+1:2*Nn-1, :) = self % xyz_all(float_m1, 2:Nn, :)

    ! Cumulative arc-length along the (possibly non-straight) floater axis
    ALLOCATE(s(2*Nn-1))
    s(1) = 0.0_WP
    do i = 2, 2 * Nn - 1
        s(i) = s(i-1) + norm2(float_pos(i, :) - float_pos(i-1, :))
    end do

    ! Calculate spacing intervals
    ds_float = divide_span(s)

    ! Submerged elliptical cross-section area at the floater's draft (z of its first node)
    z_draft = float_pos(1, 3)
    A_outer = integral_elipse_z(self % Diam_y, self % Diam_z, -z_draft)
    A_inner = integral_elipse_z(self % Diam_y - 2.0_WP * self % t, self % Diam_z - 2.0_WP * self % t, -z_draft)
    A_shell = A_outer - A_inner

    ! Distributed structural mass and added mass for the representative floater
    ALLOCATE(mass_float(2*Nn-1), added_mass_float(2*Nn-1))
    mass_float       = A_shell * self % rho_mat * ds_float
    added_mass_float = A_outer * self % rho_wat * ds_float

    DEALLOCATE(float_pos, s, ds_float)

    ! --- 2. Assemble structural mass grid and flatten ---
    ALLOCATE(mass_2D(Nm, Nn), added_mass_2D(Nm, Nn))
    mass_2D       = 0.0_WP
    added_mass_2D = 0.0_WP

    ! Broadcast both floaters (mapping split spans contiguously)
    do i = 1, 2
        float_m0 = self % float_members(i, 1) + 1_I32
        float_m1 = self % float_members(i, 2) + 1_I32

        mass_2D(float_m0, :) = mass_float(1:Nn)
        mass_2D(float_m1, :) = mass_float(Nn:2*Nn-1)

        added_mass_2D(float_m0, :) = added_mass_float(1:Nn)
        added_mass_2D(float_m1, :) = added_mass_float(Nn:2*Nn-1)
    end do

    DEALLOCATE(mass_float, added_mass_float)

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

END SUBROUTINE compute_mass_SAITEC2MWFloating


FUNCTION filter_frequencies_SAITEC2MWFloating(self, freqs_over) RESULT(mask)
    ! -------------------------------------------------------------------------
    ! Builds a mask that keeps only the computed frequencies lying close to the
    ! SAITEC 2 MW drivetrain excitation peaks (rotor-bearing and gearbox bar-
    ! chart harmonics at the rpm corresponding to the current wind speed), so
    ! downstream force arrays can be trimmed to the physically relevant
    ! spectral content.
    ! -------------------------------------------------------------------------

    USE IOUtils, ONLY: read_curve
    USE MathUtils, ONLY: generate_timeseries_banded_sines, filter_non_useful_freqs

    CLASS(SAITEC2MWFloating_t), INTENT(INOUT)        :: self          ! [-] Turbine instance
    REAL(WP)                  , INTENT(IN), OPTIONAL :: freqs_over    ! [Hz] Frequencies below this are never filtered
    LOGICAL                   , ALLOCATABLE          :: mask(:)       ! [-] Frequency mask (Nfreqs), .true. to keep

    ! Local variables
    CHARACTER(len=20), ALLOCATABLE :: keys(:)
    REAL(WP)         , ALLOCATABLE :: freqs_amp(:,:)
    REAL(WP)         , ALLOCATABLE :: freqs_to_use(:)
    REAL(WP)                       :: rpm, freqs_over_

    REAL(WP), ALLOCATABLE :: tmp_signal(:)


    freqs_over_ = 10.0_WP ; if (present(freqs_over)) freqs_over_ = freqs_over

    ! Read rpm at current wind speed
    rpm = read_curve(self % path_rpm, self % WindSpeed)

    ! Generate drivetrain excitation spectrum, capped at the max resolved frequency
    CALL drivetrain2MW_excitation_spectrum(rpm=rpm, path_barRB=self % path_barRB, path_barGB=self % path_barGB, &
                                           fmax=maxval(self % Freqs), Normalised=.false., &
                                           freqs_amp=freqs_amp, keys=keys)

    ! Get allowed frequencies
    CALL generate_timeseries_banded_sines(peaks=freqs_amp, keys=keys, t=self % Time, &
                                          zeta=0.02_WP, used_freqs=.true., a=tmp_signal, &
                                          freqs_out=freqs_to_use)

    ! Build mask
    mask = filter_non_useful_freqs(self % Freqs, freqs_to_use, freqs_over_)

    DEALLOCATE(tmp_signal)

    DEALLOCATE(freqs_amp, keys, freqs_to_use)

END FUNCTION filter_frequencies_SAITEC2MWFloating


SUBROUTINE save_specific_params_SAITEC2MWFloating(self, save_path, istr)
    ! -------------------------------------------------------------------------
    ! Persists the platform-specific geometry and material parameters
    ! (elliptical diameters, wall thickness, densities, wetted area, member ID
    ! maps) to the solver's HDF5 output, optionally suffixing each variable
    ! name with an instance tag (istr) when several turbines are saved.
    ! -------------------------------------------------------------------------

    USE IOUtils, ONLY: save_variable

    CLASS(SAITEC2MWFloating_t), INTENT(IN)           :: self         ! [-] Turbine instance
    CHARACTER(len=*)          , INTENT(IN)           :: save_path    ! [-] Path to the HDF5 output file
    CHARACTER(len=*)          , INTENT(IN), OPTIONAL :: istr         ! [-] Suffix appended to the variable names


    if (present(istr)) then
        CALL save_variable(save_path, "Case_type_"       // istr, self % case_type)
        CALL save_variable(save_path, "Structure_nodes_" // istr, self % xyz_all)
        CALL save_variable(save_path, "Diam_y_"          // istr, self % Diam_y)
        CALL save_variable(save_path, "Diam_z_"          // istr, self % Diam_z)
        CALL save_variable(save_path, "t_"               // istr, self % t)
        CALL save_variable(save_path, "rho_mat_"         // istr, self % rho_mat)
        CALL save_variable(save_path, "wet_area_"        // istr, self % wet_area)
        CALL save_variable(save_path, "float_members_"   // istr, self % float_members)
        CALL save_variable(save_path, "joint_members_"   // istr, self % joint_members)
    else
        CALL save_variable(save_path, "Case_type"      , self % case_type)
        CALL save_variable(save_path, "Structure_nodes", self % xyz_all)
        CALL save_variable(save_path, "Diam_y"         , self % Diam_y)
        CALL save_variable(save_path, "Diam_z"         , self % Diam_z)
        CALL save_variable(save_path, "t"              , self % t)
        CALL save_variable(save_path, "rho_mat"        , self % rho_mat)
        CALL save_variable(save_path, "wet_area"       , self % wet_area)
        CALL save_variable(save_path, "float_members"  , self % float_members)
        CALL save_variable(save_path, "joint_members"  , self % joint_members)
    end if

END SUBROUTINE save_specific_params_SAITEC2MWFloating


SUBROUTINE write_specific_params_txt_SAITEC2MWFloating(self, unit)
    ! -------------------------------------------------------------------------
    ! Writes a short, human-readable summary of the platform-specific geometry
    ! and member maps to the given text unit, mirroring what
    ! save_specific_params_SAITEC2MWFloating stores in the HDF5 file.
    ! -------------------------------------------------------------------------

    CLASS(SAITEC2MWFloating_t), INTENT(IN) :: self    ! [-] Turbine instance
    INTEGER(I32)              , INTENT(IN) :: unit    ! [-] Unit of the open .dat file


    write(unit, '(A,F8.3,A)')   'Floater diameter y   : ', self % Diam_y, ' m'
    write(unit, '(A,F8.3,A)')   'Floater diameter z   : ', self % Diam_z, ' m'
    write(unit, '(A,F8.3,A)')   'Wall thickness       : ', self % t, ' m'
    write(unit, '(A,F8.1,A)')   'Water density        : ', self % rho_wat, ' kg/m^3'
    write(unit, '(A,F8.1,A)')   'Material density     : ', self % rho_mat, ' kg/m^3'
    write(unit, '(A,F10.2,A)')  'Wetted area          : ', self % wet_area, ' m^2'
    write(unit, '(A,2(I0,1X))') 'Floater members (1st): ', self % float_members(:, 1)
    write(unit, '(A,2(I0,1X))') 'Floater members (2nd): ', self % float_members(:, 2)
    write(unit, '(A,3(I0,1X))') 'Joint members (excl.): ', self % joint_members

END SUBROUTINE write_specific_params_txt_SAITEC2MWFloating


FUNCTION get_impedance_corrected_force_SAITEC2MWFloating(self, c) RESULT(corrected_F)
    ! -------------------------------------------------------------------------
    ! Applies a cylindrical-scattering (Hankel-function based) shielding
    ! correction to the nodal forces. Both elliptical hulls share the same
    ! cross-section, so a single equivalent diameter D_eff = sqrt(Diam_y*Diam_z)/2
    ! is used uniformly for every wetted node, and the resulting
    ! frequency-dependent correction factor scales the raw force so that a
    ! solver using the Image Method sees the acoustic shadowing effect of the
    ! structure itself.
    ! -------------------------------------------------------------------------

    USE MathUtils, ONLY: alpha_hankel

    CLASS(SAITEC2MWFloating_t), INTENT(IN)  :: self                  ! [-] Turbine instance
    REAL(WP)                  , INTENT(IN)  :: c                     ! [m/s] Speed of sound in the fluid
    COMPLEX(WP)               , ALLOCATABLE :: corrected_F(:,:,:)    ! [N] Corrected force (Nfreqs, Nnodes_wet, 3)

    ! Local variables
    INTEGER(I32)              :: nf, Nnodes_wet, i, j, m
    REAL(WP)    , ALLOCATABLE :: k(:)
    COMPLEX(WP) , ALLOCATABLE :: alpha_floats(:)
    REAL(WP)                  :: D_eff


    nf         = size(self % Freqs)
    Nnodes_wet = size(self % xyz, 1)

    ! Wavenumber vector k
    ALLOCATE(k(nf))
    k = (2.0_WP * PI * self % Freqs) / c

    ! Equivalent diameter for the elliptical cross-section
    D_eff = sqrt(self % Diam_y * self % Diam_z) / 2.0_WP

    alpha_floats = alpha_hankel(k, D_eff)    ! Shape (nf)

    ! --- Apply correction: F_corrected = F * abs(alpha_floats) ---
    ALLOCATE(corrected_F(nf, Nnodes_wet, 3))

    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i, j, m)
    do m = 1, 3
        do j = 1, Nnodes_wet
            do i = 1, nf
                corrected_F(i, j, m) = self % Force(i, j, m) * abs(alpha_floats(i))
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    ! Cleanup
    DEALLOCATE(k, alpha_floats)

END FUNCTION get_impedance_corrected_force_SAITEC2MWFloating


! ===============================================================================
! SAITEC 2 MW DRIVETRAIN EXCITATION (PRIVATE HELPERS)
! ===============================================================================
SUBROUTINE drivetrain2MW_excitation_spectrum(rpm, path_barRB, path_barGB, fmin, fmax, tol, Normalised, &
                                             freqs_amp, keys)
    ! -------------------------------------------------------------------------
    ! Constructs the excitation spectrum (frequencies + amplitudes) for the
    ! SAITEC 2 MW variable-speed drivetrain described in Escaler & Mebarki
    ! (2018), combining rotor-bearing and gearbox kinematic harmonics with
    ! externally-supplied bar-chart amplitude data, merging near-duplicate
    ! frequencies (keeping the larger amplitude) and optionally normalizing to
    ! unit RMS.
    ! -------------------------------------------------------------------------

    REAL(WP)         , INTENT(IN)               :: rpm               ! [rpm] Rotor speed
    CHARACTER(len=*) , INTENT(IN)               :: path_barRB        ! [-] Path to the rotor-bearing bar-chart CSV file
    CHARACTER(len=*) , INTENT(IN)               :: path_barGB        ! [-] Path to the gearbox bar-chart CSV file
    REAL(WP)         , INTENT(IN), OPTIONAL     :: fmin              ! [Hz] Lower frequency limit (inclusive)
    REAL(WP)         , INTENT(IN), OPTIONAL     :: fmax              ! [Hz] Upper frequency limit (exclusive)
    REAL(WP)         , INTENT(IN), OPTIONAL     :: tol               ! [Hz] Tolerance to merge close frequencies
    LOGICAL          , INTENT(IN), OPTIONAL     :: Normalised        ! [-] Whether to normalize to unit RMS
    REAL(WP)         , INTENT(OUT), ALLOCATABLE :: freqs_amp(:,:)    ! [-] Peaks (N, 2): frequency [Hz], amplitude
    CHARACTER(len=20), INTENT(OUT), ALLOCATABLE :: keys(:)           ! [-] Peak labels (N)

    ! Local variables
    ! Kinematic peak frequencies (rotor-bearing RB1/RB2, gearbox GB1/GB2/GB3)
    REAL(WP) :: RB1f(12), RB2f(11), GB1f(9), GB2f(9), GB3f(11)
    REAL(WP) :: f3, f_em
    REAL(WP) :: fmin_, fmax_, tol_
    LOGICAL  :: normalised_

    REAL(WP), ALLOCATABLE :: RB_data(:,:), GB_data(:,:)
    REAL(WP), ALLOCATABLE :: RB1_amp(:), RB2_amp(:), GB1_amp(:), GB2_amp(:), GB3_amp(:)

    REAL(WP)         , ALLOCATABLE :: freqs_raw(:), amps_raw(:)
    CHARACTER(len=20), ALLOCATABLE :: keys_raw(:)

    REAL(WP)         , ALLOCATABLE :: merged_freqs(:), merged_amps(:)
    CHARACTER(len=20), ALLOCATABLE :: merged_keys(:)
    LOGICAL          , ALLOCATABLE :: mask(:)

    REAL(WP)          :: rms, fk, ak
    CHARACTER(len=20) :: key_tmp
    INTEGER(I32)      :: n_total, i, j, kk, best_idx, n_merged, n_kept


    ! Defaults
    fmin_       = 0.0_WP       ; if (present(fmin))       fmin_       = fmin
    fmax_       = huge(1.0_WP) ; if (present(fmax))       fmax_       = fmax
    tol_        = 1.0e-3_WP    ; if (present(tol))        tol_        = tol
    normalised_ = .true.       ; if (present(Normalised)) normalised_ = Normalised

    CALL compute_freqs2MW(rpm, RB1f, RB2f, GB1f, GB2f, GB3f, f3)

    f_em = 84.0_WP * f3

    ! --- 1. Read bar-chart amplitude data (CSV columns, "Nan" = missing) ---
    CALL read_bar_csv(path_barRB, 2_I32, RB_data)
    CALL read_bar_csv(path_barGB, 3_I32, GB_data)

    RB1_amp = pack(RB_data(:, 1), .not. ieee_is_nan(RB_data(:, 1)))
    RB2_amp = pack(RB_data(:, 2), .not. ieee_is_nan(RB_data(:, 2)))
    GB1_amp = pack(GB_data(:, 1), .not. ieee_is_nan(GB_data(:, 1)))
    GB2_amp = pack(GB_data(:, 2), .not. ieee_is_nan(GB_data(:, 2)))
    GB3_amp = pack(GB_data(:, 3), .not. ieee_is_nan(GB_data(:, 3)))

    if (size(RB1_amp) /= size(RB1f) .or. size(RB2_amp) /= size(RB2f) .or. &
        size(GB1_amp) /= size(GB1f) .or. size(GB2_amp) /= size(GB2f) .or. &
        size(GB3_amp) /= size(GB3f)) then
        error stop "drivetrain2MW_excitation_spectrum: bar-chart CSV row counts do not match " // &
                   "the expected kinematic peak counts"
    end if

    ! --- 2. Assemble raw (freq, amp, key) lists ---
    ! (kinematic frequency order preserved; amplitude re-paired positionally
    ! from each family's own frequency-sorted bar-chart column)
    n_total = size(RB1f) + size(RB2f) + size(GB1f) + size(GB2f) + size(GB3f) + 1_I32
    ALLOCATE(freqs_raw(n_total), amps_raw(n_total), keys_raw(n_total))

    j = 0_I32
    CALL append_family(RB1f, sorted_amps_by_freq(RB1f, RB1_amp), 'RB1_f_', freqs_raw, amps_raw, keys_raw, j)
    CALL append_family(RB2f, sorted_amps_by_freq(RB2f, RB2_amp), 'RB2_f_', freqs_raw, amps_raw, keys_raw, j)
    CALL append_family(GB1f, sorted_amps_by_freq(GB1f, GB1_amp), 'GB1_f_', freqs_raw, amps_raw, keys_raw, j)
    CALL append_family(GB2f, sorted_amps_by_freq(GB2f, GB2_amp), 'GB2_f_', freqs_raw, amps_raw, keys_raw, j)
    CALL append_family(GB3f, sorted_amps_by_freq(GB3f, GB3_amp), 'GB3_f_', freqs_raw, amps_raw, keys_raw, j)

    ! Electromagnetic peak at 84xf3 (observed in generator)
    j            = j + 1_I32
    freqs_raw(j) = f_em
    amps_raw(j)  = 0.2_WP
    keys_raw(j)  = 'generator_84xf3'

    DEALLOCATE(RB_data, GB_data, RB1_amp, RB2_amp, GB1_amp, GB2_amp, GB3_amp)

    ! --- 3. Sort everything by frequency (carrying amps/keys along, insertion sort) ---
    do i = 2, n_total
        fk      = freqs_raw(i)
        ak      = amps_raw(i)
        key_tmp = keys_raw(i)
        kk      = i - 1_I32
        do while (kk >= 1_I32)
            if (freqs_raw(kk) <= fk) exit
            freqs_raw(kk+1) = freqs_raw(kk)
            amps_raw(kk+1)  = amps_raw(kk)
            keys_raw(kk+1)  = keys_raw(kk)
            kk              = kk - 1_I32
        end do
        freqs_raw(kk+1) = fk
        amps_raw(kk+1)  = ak
        keys_raw(kk+1)  = key_tmp
    end do

    ! --- 4. Merge near-duplicate frequencies (within tol), keep larger amplitude ---
    ALLOCATE(merged_freqs(n_total), merged_amps(n_total), merged_keys(n_total))
    n_merged = 0_I32
    i        = 1_I32
    do while (i <= n_total)
        kk = i
        do while (kk + 1_I32 <= n_total)
            if ((freqs_raw(kk+1) - freqs_raw(i)) < tol_) then
                kk = kk + 1_I32
            else
                exit
            end if
        end do
        best_idx               = i - 1_I32 + maxloc(amps_raw(i:kk), dim=1)
        n_merged               = n_merged + 1_I32
        merged_freqs(n_merged) = freqs_raw(best_idx)
        merged_amps(n_merged)  = amps_raw(best_idx)
        merged_keys(n_merged)  = keys_raw(best_idx)
        i                      = kk + 1_I32
    end do
    DEALLOCATE(freqs_raw, amps_raw, keys_raw)

    ! --- 5. Apply frequency range filter (fmin, fmax) ---
    ALLOCATE(mask(n_merged))
    mask   = (merged_freqs(1:n_merged) >= fmin_) .and. (merged_freqs(1:n_merged) < fmax_)
    n_kept = count(mask)

    ALLOCATE(freqs_amp(n_kept, 2), keys(n_kept))
    j = 0_I32
    do i = 1, n_merged
        if (mask(i)) then
            j               = j + 1_I32
            freqs_amp(j, 1) = merged_freqs(i)
            freqs_amp(j, 2) = merged_amps(i)
            keys(j)         = merged_keys(i)
        end if
    end do
    DEALLOCATE(merged_freqs, merged_amps, merged_keys, mask)

    ! --- 6. Normalize to unit RMS: sqrt(sum(A^2/2)) = 1 ---
    if (normalised_ .and. n_kept > 0_I32) then
        rms             = sqrt(sum(freqs_amp(:, 2)**2) / 2.0_WP)
        freqs_amp(:, 2) = freqs_amp(:, 2) / rms
    end if

END SUBROUTINE drivetrain2MW_excitation_spectrum


SUBROUTINE compute_freqs2MW(rpm, RB1f, RB2f, GB1f, GB2f, GB3f, f3_out)
    ! -------------------------------------------------------------------------
    ! Computes the SAITEC 2 MW drivetrain kinematic peak frequencies (rotor
    ! bearings RB1/RB2 and gearbox stages GB1/GB2/GB3) at the given rotor
    ! speed, following the planetary+parallel gear-train kinematics of the
    ! Escaler & Mebarki (2018) 2 MW drivetrain.
    ! -------------------------------------------------------------------------

    REAL(WP), INTENT(IN)  :: rpm         ! [rpm] Rotor speed
    REAL(WP), INTENT(OUT) :: RB1f(12)    ! [Hz] Rotor bearing 1 kinematic peak frequencies
    REAL(WP), INTENT(OUT) :: RB2f(11)    ! [Hz] Rotor bearing 2 kinematic peak frequencies
    REAL(WP), INTENT(OUT) :: GB1f(9)     ! [Hz] Gearbox stage 1 kinematic peak frequencies
    REAL(WP), INTENT(OUT) :: GB2f(9)     ! [Hz] Gearbox stage 2 kinematic peak frequencies
    REAL(WP), INTENT(OUT) :: GB3f(11)    ! [Hz] Gearbox stage 3 kinematic peak frequencies
    REAL(WP), INTENT(OUT) :: f3_out      ! [Hz] High-speed shaft rotation frequency

    ! Local variables
    REAL(WP), PARAMETER :: zs = 18.0_WP, zp = 34.0_WP, zr = 87.0_WP
    REAL(WP), PARAMETER :: z12 = 70.0_WP, z21 = 16.0_WP, z23 = 84.0_WP, z32 = 19.0_WP
    REAL(WP)            :: f0, f1, fp, f2, f3, gmfp, gmf12, gmf23, fb


    f0    = rpm / 60.0_WP
    f1    = f0 * (zs + zr) / zs
    fp    = zs / zp * (f1 - f0)
    f2    = z12 / z21 * f1
    f3    = z23 / z32 * f2
    gmfp  = zp * fp
    gmf12 = z12 * f1
    gmf23 = z23 * f2
    fb    = 3.0_WP * f0

    RB1f = [f0, fb, 2.0_WP * fb, 3.0_WP * fb, 4.0_WP * fb, gmfp, f3, 2.0_WP * gmfp, 2.0_WP * f3, &
            3.0_WP * gmfp, 3.0_WP * f3, 4.0_WP * gmfp]
    RB2f = [f0, fb,              3.0_WP * fb, 4.0_WP * fb, gmfp, f3, 2.0_WP * gmfp, 2.0_WP * f3, &
            3.0_WP * gmfp, 3.0_WP * f3, 4.0_WP * gmfp]
    GB1f = [gmfp, f3, 2.0_WP * gmfp, 2.0_WP * f3, 3.0_WP * gmfp, 3.0_WP * f3, 4.0_WP * gmfp, 5.0_WP * gmfp, &
            6.0_WP * gmfp]
    GB2f = [gmfp, f3, gmf12, 2.0_WP * gmf12, 3.0_WP * gmf12, 4.0_WP * gmf12, 5.0_WP * gmf12, 6.0_WP * gmf12, &
            gmf23]
    GB3f = [gmf12, 2.0_WP * gmf12, 3.0_WP * gmf12, 4.0_WP * gmf12, 5.0_WP * gmf12, 6.0_WP * gmf12, &
            gmf23, 2.0_WP * gmf23, 3.0_WP * gmf23, 4.0_WP * gmf23, 5.0_WP * gmf23]

    f3_out = f3

END SUBROUTINE compute_freqs2MW


FUNCTION sorted_amps_by_freq(freqs_kin, amp_col) RESULT(amps_sorted)
    ! -------------------------------------------------------------------------
    ! Pairs a kinematic frequency family with its bar-chart amplitude column
    ! and sorts the pair by frequency, returning only the resulting amplitude
    ! order (matching the validated reference implementation, which re-pairs
    ! this frequency-sorted amplitude array positionally with the ORIGINAL,
    ! unsorted kinematic frequency array).
    ! -------------------------------------------------------------------------

    REAL(WP), INTENT(IN)  :: freqs_kin(:)      ! [Hz] Kinematic frequency family
    REAL(WP), INTENT(IN)  :: amp_col(:)        ! [-] Bar-chart amplitude column of the family
    REAL(WP), ALLOCATABLE :: amps_sorted(:)    ! [-] Amplitudes sorted by their paired frequency

    ! Local variables
    REAL(WP)    , ALLOCATABLE :: f_sort(:)
    INTEGER(I32)              :: n, i, k
    REAL(WP)                  :: fk, ak


    n = size(freqs_kin)
    ALLOCATE(f_sort(n), amps_sorted(n))
    f_sort      = freqs_kin
    amps_sorted = amp_col

    do i = 2, n
        fk = f_sort(i)
        ak = amps_sorted(i)
        k  = i - 1_I32
        do while (k >= 1_I32)
            if (f_sort(k) <= fk) exit
            f_sort(k+1)      = f_sort(k)
            amps_sorted(k+1) = amps_sorted(k)
            k                = k - 1_I32
        end do
        f_sort(k+1)      = fk
        amps_sorted(k+1) = ak
    end do

    DEALLOCATE(f_sort)

END FUNCTION sorted_amps_by_freq


SUBROUTINE append_family(freqs_kin, amps_kin, prefix, freqs_raw, amps_raw, keys_raw, j)
    ! -------------------------------------------------------------------------
    ! Appends one kinematic-frequency family (with its already frequency-
    ! sorted amplitudes) into the raw (freq, amp, key) accumulation arrays,
    ! advancing the shared cursor j and building a "<prefix>f.ffHz" label key.
    ! -------------------------------------------------------------------------

    REAL(WP)         , INTENT(IN)    :: freqs_kin(:)    ! [Hz] Kinematic frequency family
    REAL(WP)         , INTENT(IN)    :: amps_kin(:)     ! [-] Amplitudes paired with freqs_kin
    CHARACTER(len=*) , INTENT(IN)    :: prefix          ! [-] Label prefix of the family (e.g. 'RB1_f_')
    REAL(WP)         , INTENT(INOUT) :: freqs_raw(:)    ! [Hz] Accumulated frequencies
    REAL(WP)         , INTENT(INOUT) :: amps_raw(:)     ! [-] Accumulated amplitudes
    CHARACTER(len=20), INTENT(INOUT) :: keys_raw(:)     ! [-] Accumulated labels
    INTEGER(I32)     , INTENT(INOUT) :: j               ! [-] Shared cursor: last filled position

    ! Local variables
    INTEGER(I32) :: i


    do i = 1, size(freqs_kin)
        j            = j + 1_I32
        freqs_raw(j) = freqs_kin(i)
        amps_raw(j)  = amps_kin(i)
        write(keys_raw(j), '(A,F0.2,A)') trim(prefix), freqs_kin(i), 'Hz'
    end do

END SUBROUTINE append_family


SUBROUTINE read_bar_csv(filename, ncols, data)
    ! -------------------------------------------------------------------------
    ! Reads a comma-separated bar-chart amplitude file (one header line
    ! followed by ncols numeric columns, "Nan" marking a missing entry) into a
    ! (nrows, ncols) array, preserving IEEE NaN for missing entries so callers
    ! can strip them with ieee_is_nan.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)               :: filename     ! [-] Path to the bar-chart CSV file
    INTEGER(I32)    , INTENT(IN)               :: ncols        ! [-] Number of numeric columns to read
    REAL(WP)        , INTENT(OUT), ALLOCATABLE :: data(:,:)    ! [-] Bar-chart amplitudes (nrows, ncols), NaN if missing

    ! Local variables
    INTEGER(I32)        :: fu, ios, i, nrows
    CHARACTER(len=2048) :: line


    open(newunit=fu, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0_I32) error stop "read_bar_csv: could not open file " // trim(filename)

    read(fu, '(A)', iostat=ios) line    ! Header line

    nrows = 0_I32
    do
        read(fu, '(A)', iostat=ios) line
        if (ios /= 0_I32) exit
        if (len_trim(line) == 0) cycle
        nrows = nrows + 1_I32
    end do

    ALLOCATE(data(nrows, ncols))

    rewind(fu)
    read(fu, '(A)') line    ! Skip header again

    do i = 1, nrows
        read(fu, '(A)', iostat=ios) line
        read(line, *) data(i, 1:ncols)
    end do

    close(fu)

END SUBROUTINE read_bar_csv

END MODULE SAITEC2MWFloating
