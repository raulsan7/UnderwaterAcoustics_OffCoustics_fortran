! ===============================================================================
! MODULE: SAITEC2MWFloating
! ===============================================================================
! Defines SAITEC2MWFloating_t, a concrete turbine type extending the abstract
! WindTurbine_t (src/emitters/WindTurbine.f90) for the SAITEC SENVION 2 MW floating
! platform (SATH dual-hull concrete technology). It implements the deferred
! procedures: compute_force, compute_mass, excitation_spectrum,
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
USE WindTurbine, ONLY: WindTurbine_t
USE ConfigTypes, ONLY: Config_WindTurbine_t
USE DrivetrainExcitation, ONLY: drivetrain2MW_excitation_spectrum

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
    REAL(WP) :: Diam_y   = 10.0_WP                  ! [m] Elliptical floater diameter in y-direction
    REAL(WP) :: Diam_z   = 12.5_WP                  ! [m] Elliptical floater diameter in z-direction
    REAL(WP) :: rho_wat  = 1025.0_WP                ! [kg/m^3] Water density
    REAL(WP) :: rho_mat  = 1850.0_WP                ! [kg/m^3] Material density (concrete)
    REAL(WP) :: t        = 1.0_WP                   ! [m] Wall thickness
    REAL(WP) :: wet_area = 0.0_WP                   ! [m^2] Wetted area (not computed for this platform)

    ! [-] Floater members ID list: [[0,1],[3,4]]
    INTEGER(I32) :: float_members(2,2) = reshape([0_I32, 3_I32, 1_I32, 4_I32], [2, 2])
    ! [-] Virtual rigid-link (joint) member IDs, excluded
    INTEGER(I32) :: joint_members(3) = [5_I32, 2_I32, 6_I32]

    ! Hardcoded paths to the SAITEC 2 MW drivetrain bar-chart amplitude data
    CHARACTER(len=:), ALLOCATABLE :: path_barRB     ! [-] Rotor-bearing bar-chart data file
    CHARACTER(len=:), ALLOCATABLE :: path_barGB     ! [-] Gearbox bar-chart data file

    CONTAINS

    ! --- Shared procedures --- !
    PROCEDURE :: init       => init_SAITEC2MWFloating
    PROCEDURE :: read_input => read_input_SAITEC2MWFloating

    ! --- Deferred (abstract) procedures --- !
    PROCEDURE :: compute_force                 => compute_force_SAITEC2MWFloating
    PROCEDURE :: compute_mass                  => compute_mass_SAITEC2MWFloating
    PROCEDURE :: excitation_spectrum           => excitation_spectrum_SAITEC2MWFloating
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
    ! base-class init for the common parameters and checks that the hardcoded
    ! member IDs exist in the structural model.
    ! -------------------------------------------------------------------------

    CLASS(SAITEC2MWFloating_t), INTENT(INOUT) :: self      ! [-] Turbine instance
    TYPE(Config_WindTurbine_t), INTENT(IN)    :: config    ! [-] Turbine configuration

    ! Local variables
    LOGICAL      :: exists
    INTEGER(I32) :: n_req


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

    ! Flag the turbine as floating
    self % is_monopile = .false.

    CALL self % init_base(config)

    ! The 0-based float_members/joint_members IDs must exist in the structural model
    n_req = max(maxval(self % float_members), maxval(self % joint_members)) + 1_I32
    if (n_req > self % Nmembers) then
        write(*, '(A,I0,A,I0)') ' ERROR SAITEC2MWFloating % init: Nmembers = ', self % Nmembers, &
                                ' but float_members/joint_members need Nmembers >= ', n_req
        error stop
    end if

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
    ! compute_force's keep mask). Like the base implementation, it sets the
    ! FFT bin width df and computes the automatic barycenter in the global
    ! frame (after the rotation and translation), guarding it against the NaN
    ! padding that the rigid-link members (with fewer than Nnodes real nodes)
    ! introduce into self % xyz_all.
    ! -------------------------------------------------------------------------

    USE OpenFASTReader, ONLY: get_SDsum_variables, read_input_SD

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

    ! FFT bin width of the time series, defines the shared frequency grid of a wind farm
    self % df = 1.0_WP / (real(size(self % Time), WP) * (self % Time(2) - self % Time(1)))

    ! Align with wind direction and translate to axis position
    CALL self % translate_align_with_WindDir()

    ! Compute barycenter if not provided, in the global frame. A plain
    ! node-count average (as the base class would do) is wrong here: node
    ! spacing along each floater is highly non-uniform (long vs. short
    ! segment) and the rigid-link members contribute duplicated/NaN-padded
    ! nodes clustered at the central pivot. Instead, use the actual structural
    ! (dry) mass distribution as weights, which correctly reduces to the true
    ! mass centroid of the wetted floaters and is exactly zero (by
    ! construction) at the excluded joints. The masses only depend on arc
    ! lengths and on the draft, so they are not changed by the rotation and
    ! translation applied above.
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

    if (verbose) then
        dt = self % Time(2) - self % Time(1)
        write(*,*) ''
        write(*, '(A)') repeat('=', 70)
        write(*, '(A)') 'TURBINE: ' // trim(self % rootname)
        write(*, '(A)') repeat('-', 68)
        write(*, '(A, A)') 'Type: ', trim(self % case_type)
        write(*, '(A, F5.1, A, F6.1, A)') 'Wind Speed: ', self % WindSpeed, ' m/s   | Wind Direction: ', &
                                          self % WindDir, ' deg'
        write(*, '(A, F5.1, A, F8.1, A, F8.1, A)') 'Water Depth: ', self % Depth, &
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

    USE SignalProcessing, ONLY: compute_rfft

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
    skipf_        = config % skipf

    Nm = self % Nmembers
    Nn = self % Nnodes
    nt = size(self % Time)
    dt = self % Time(2) - self % Time(1)

    ! Acceleration flatten
    ALLOCATE(acc_flat(nt, Nm*Nn, 3))
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, global_idx)
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
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, global_idx)
    do i = 1, Nm
        do j = 1, Nn
            global_idx     = (i - 1_I32) * Nn + j
            mass_eff(i, j) = mass(global_idx) + added_mass(global_idx)
        end do
    end do
    !$OMP END PARALLEL DO

    DEALLOCATE(mass, added_mass)

    ! Compute force via F = A * mass_eff
    ALLOCATE(F_temp(nf, Nm, Nn, 3))

    ! Parallelized over space and coordinate directions
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i, j, k)
    do k = 1, 3
        do j = 1, Nn
            do i = 1, Nm
                F_temp(:, i, j, k) = A(:, i, j, k) * mass_eff(i, j)
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
        F_memory_mb = real(storage_size(self % Force), WP) * real(size(self % Force, kind=I64), WP) &
                      / 8.0_WP / 1024.0_WP / 1024.0_WP
        write(*, '(A)') repeat('=', 70)
        write(*, '(A)') 'FORCE COMPUTATION COMPLETED (FLOATING PLATFORM)'
        write(*, '(A)') repeat('-', 68)
        write(*, '(A, I6)') 'Nodes (wet): ', Nnodes_wet
        write(*, '(A, I6)') 'Frequencies: ', nf
        write(*, '(A, F10.2, A)') 'Memory (F):  ', F_memory_mb, ' MB'
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

    USE StructureModel, ONLY: divide_span, integral_elipse_z

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


SUBROUTINE excitation_spectrum_SAITEC2MWFloating(self, rpm, freqs_amp, keys)
    ! -------------------------------------------------------------------------
    ! Builds the SAITEC 2 MW drivetrain excitation spectrum (rotor-bearing and
    ! gearbox bar-chart harmonics) at the given rotor speed, capped at the
    ! maximum resolved frequency, used by filter_frequencies.
    ! -------------------------------------------------------------------------

    CLASS(SAITEC2MWFloating_t), INTENT(IN)               :: self              ! [-] Turbine instance
    REAL(WP)                  , INTENT(IN)               :: rpm               ! [rpm] Rotor speed
    REAL(WP)                  , INTENT(OUT), ALLOCATABLE :: freqs_amp(:,:)    ! [-] Peaks (N, 2): freq [Hz], amplitude
    CHARACTER(len=20)         , INTENT(OUT), ALLOCATABLE :: keys(:)           ! [-] Peak labels (N)


    CALL drivetrain2MW_excitation_spectrum(rpm=rpm, path_barRB=self % path_barRB, path_barGB=self % path_barGB, &
                                           fmax=maxval(self % Freqs), Normalised=.false., &
                                           freqs_amp=freqs_amp, keys=keys)

END SUBROUTINE excitation_spectrum_SAITEC2MWFloating


SUBROUTINE save_specific_params_SAITEC2MWFloating(self, save_path, istr)
    ! -------------------------------------------------------------------------
    ! Persists the platform-specific geometry and material parameters
    ! (elliptical diameters, wall thickness, material density, wetted area, member ID
    ! maps) to the solver's HDF5 output, optionally suffixing each variable
    ! name with an instance tag (istr) when several turbines are saved.
    ! -------------------------------------------------------------------------

    USE HDF5Utils, ONLY: save_variable

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

    USE StructureModel, ONLY: alpha_hankel

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
    D_eff = sqrt(self % Diam_y * self % Diam_z)

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

END MODULE SAITEC2MWFloating
