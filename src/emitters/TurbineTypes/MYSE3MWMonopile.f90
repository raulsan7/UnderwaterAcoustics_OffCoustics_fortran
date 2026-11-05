! ===============================================================================
! MODULE: MYSE3MWMonopile
! ===============================================================================
! Defines MYSE3MWMonopile_t, a concrete turbine type extending the abstract
! WindTurbine_t (src/emitters/WindTurbine.f90) for the MingYang MySE 3.0-135
! on a monopile of the Beleolico wind farm (Taranto), modelled in OpenFAST in
! Beleolico_offshore_turbine. It implements the deferred procedures: force
! computation from OpenFAST accelerations, distributed dry/added mass along the
! monopile, the drivetrain excitation spectrum (the pure tones of the
! prescribed nacelle force of the runs), saving/writing of platform-specific
! parameters, and the cylindrical-scattering impedance correction.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 11/2026
! ===============================================================================
MODULE MYSE3MWMonopile

USE Kinds
USE WindTurbine, ONLY: WindTurbine_t
USE ConfigTypes, ONLY: Config_WindTurbine_t
USE DrivetrainExcitation, ONLY: drivetrain3MW_excitation_spectrum

IMPLICIT NONE

PRIVATE
PUBLIC :: MYSE3MWMonopile_t


TYPE, EXTENDS(WindTurbine_t) :: MYSE3MWMonopile_t
    ! -------------------------------------------------------------------------
    ! Concrete turbine type extending WindTurbine_t for the MySE 3.0-135 of
    ! Beleolico on a monopile (14 m of water, SubDyn members M1-M7 from the
    ! seabed to z = 0). Stores the monopile diameters, material/water
    ! densities and wetted area, and the analysis window of the drivetrain
    ! force, and implements every deferred procedure of the abstract parent
    ! (force, mass, frequency filtering, parameter I/O and impedance
    ! correction).
    ! -------------------------------------------------------------------------

    ! Physical and structural parameters specific to the MySE 3.0-135 monopile of Beleolico
    REAL(WP) :: D        = 5.0_WP               ! [m] Monopile outer diameter
    REAL(WP) :: D_inner  = 4.88_WP              ! [m] Monopile inner diameter (wall thickness 60 mm)
    REAL(WP) :: rho_wat  = 1026.5_WP            ! [kg/m^3] Seawater density (WtrDens of the OpenFAST model)
    REAL(WP) :: rho_mat  = 7850.0_WP            ! [kg/m^3] Material density (steel)
    REAL(WP) :: wet_area = 219.91_WP            ! [m^2] Wetted area (pi D x 14 m of water)

    ! Read by init from the drivetrain force file of the run when the frequencies are filtered
    REAL(WP) :: T_win    = -1.0_WP              ! [s] Analysis window of the force lines (grid 1/T_win)

    CONTAINS

    ! --- Shared procedures --- !
    PROCEDURE :: init => init_MYSE3MWMonopile

    ! --- Deferred (abstract) procedures --- !
    PROCEDURE :: compute_force                 => compute_force_MYSE3MWMonopile
    PROCEDURE :: compute_mass                  => compute_mass_MYSE3MWMonopile
    PROCEDURE :: excitation_spectrum           => excitation_spectrum_MYSE3MWMonopile
    PROCEDURE :: save_specific_params          => save_specific_params_MYSE3MWMonopile
    PROCEDURE :: write_specific_params_txt     => write_specific_params_txt_MYSE3MWMonopile
    PROCEDURE :: get_impedance_corrected_force => get_impedance_corrected_force_MYSE3MWMonopile

END TYPE MYSE3MWMonopile_t


CONTAINS


SUBROUTINE init_MYSE3MWMonopile(self, config)
    ! -------------------------------------------------------------------------
    ! Sets the monopile-specific (hardcoded) rpm-curve path and subfolder of
    ! OP_output_dir, marks the turbine as monopile-founded with pure-tone
    ! drivetrain lines, and delegates to the base WindTurbine_t % init for the
    ! remaining shared initialization steps. When the frequencies are filtered
    ! it reads the rotor speed and the analysis window of the drivetrain force
    ! of the run (<rootname>.force.json, next to the OpenFAST output) and stops
    ! if the rpm curve at WindSpeed does not give the same lines.
    ! -------------------------------------------------------------------------

    USE CSVReader, ONLY: read_curve
    USE FileUtils, ONLY: swap_extension
    USE DrivetrainExcitation, ONLY: read_drivetrain3MW_force_file

    CLASS(MYSE3MWMonopile_t)  , INTENT(INOUT) :: self      ! [-] Turbine instance
    TYPE(Config_WindTurbine_t), INTENT(IN)    :: config    ! [-] Turbine configuration

    ! Local variables
    LOGICAL                        :: exists, same_lines
    REAL(WP)                       :: rpm, rpm_force
    REAL(WP)         , ALLOCATABLE :: freqs_amp(:,:), freqs_amp_force(:,:)
    CHARACTER(len=20), ALLOCATABLE :: keys(:)
    CHARACTER(len=:) , ALLOCATABLE :: path_force


    ! Path to the rpm curve file and subfolder of the OpenFAST files (hardcoded)
    self % path_rpm  = "../wind_speed_curves_MYSE_3MW/rpm_ws.csv"
    self % OP_subdir = "MYSE3MWMonopile/"

    ! Verify that the rpm curve file exists
    inquire(file=self % path_rpm, exist=exists)
    if (.not. exists) error stop "MYSE3MWMonopile % init: RPM curve file not found: " // self % path_rpm

    ! Flag the turbine as monopile; its drivetrain lines are pure tones
    self % is_monopile = .true.
    self % pure_tones  = .true.

    CALL self % init_base(config)

    if (.not. config % filter_freqs) return

    ! --- Drivetrain force of the run, needed by filter_frequencies --- !
    path_force = swap_extension(self % OP_output_path, ".force.json")
    inquire(file=path_force, exist=exists)
    if (.not. exists) error stop "MYSE3MWMonopile % init: drivetrain force file (needed by filter_freqs) " // &
                                 "not found: " // path_force

    CALL read_drivetrain3MW_force_file(path_force, rpm_force, self % T_win)

    ! The rpm that filter_frequencies takes from the curve must give the lines of the force file
    rpm = read_curve(self % path_rpm, self % WindSpeed)
    CALL drivetrain3MW_excitation_spectrum(rpm=rpm, T_win=self % T_win, freqs_amp=freqs_amp, keys=keys)
    CALL drivetrain3MW_excitation_spectrum(rpm=rpm_force, T_win=self % T_win, freqs_amp=freqs_amp_force, keys=keys)

    if (size(freqs_amp, 1) /= size(freqs_amp_force, 1)) then
        same_lines = .false.
    else
        same_lines = all(abs(freqs_amp(:, 1) - freqs_amp_force(:, 1)) < 0.5_WP / self % T_win)
    end if

    if (.not. same_lines) then
        write(*, '(A,F0.6,A,F0.4,A,F0.6,A)') ' ERROR MYSE3MWMonopile % init: rpm curve gives ', rpm, &
                                             ' rpm at WindSpeed = ', self % WindSpeed, &
                                             ' m/s, but the drivetrain force of the run is at ', rpm_force, ' rpm'
        write(*, '(A)') ' Set WindSpeed to the wind of the run (operating_points.csv of the OpenFAST case) ' // &
                        'or filter_freqs = .false.'
        error stop
    end if

    DEALLOCATE(freqs_amp, freqs_amp_force, keys)

END SUBROUTINE init_MYSE3MWMonopile


SUBROUTINE compute_force_MYSE3MWMonopile(self, config, verbose)
    ! -------------------------------------------------------------------------
    ! Computes the nodal structural forces from OpenFAST accelerations: removes
    ! duplicate/dry (above-water) nodes, computes the real FFT of the wet-node
    ! accelerations, and multiplies by the effective (dry + added) mass to get
    ! F = m*a in the frequency domain. Optionally keeps only the frequencies
    ! relevant to the drivetrain excitation.
    ! -------------------------------------------------------------------------

    USE SignalProcessing, ONLY: compute_rfft
    USE StructureModel, ONLY: remove_duplicate_nodes

    CLASS(MYSE3MWMonopile_t)  , INTENT(INOUT) :: self       ! [-] Turbine instance
    TYPE(Config_WindTurbine_t), INTENT(IN)    :: config     ! [-] Turbine configuration
    LOGICAL                   , INTENT(IN)    :: verbose    ! [-] Whether to print a summary

    ! Local variables
    INTEGER(I32) :: Nm, Nn, nt, Nnodes_wet, i, j, k, v, global_idx, Nfreqs_new
    REAL(WP)     :: dt, F_memory_mb

    LOGICAL     , ALLOCATABLE :: mask_duplicate(:), mask_freqs_to_use(:)
    REAL(WP)    , ALLOCATABLE :: mass(:), added_mass(:), mass_eff(:)
    REAL(WP)    , ALLOCATABLE :: x_flat(:,:), acc_flat(:,:,:)
    COMPLEX(WP) , ALLOCATABLE :: A(:,:,:)
    REAL(WP)    , ALLOCATABLE :: freqs(:), freqs_new(:)
    COMPLEX(WP) , ALLOCATABLE :: F_filtered(:,:,:)
    INTEGER(I32), ALLOCATABLE :: valid_j(:), valid_k(:)


    Nm = self % Nmembers
    Nn = self % Nnodes
    nt = size(self % Time)
    dt = self % Time(2) - self % Time(1)

    mask_duplicate = remove_duplicate_nodes(self % xyz_all, delete_last=.true.)

    ! Fuse masks: count all valid nodes at once
    Nnodes_wet = 0_I32
    do j = 1, Nm
        do k = 1, Nn
            global_idx = (j - 1_I32) * Nn + k
            if (mask_duplicate(global_idx) .and. self % xyz_all(j, k, 3) <= 0.0_WP) then
                Nnodes_wet = Nnodes_wet + 1_I32
            end if
        end do
    end do

    ALLOCATE(x_flat(Nnodes_wet, 3))
    ALLOCATE(acc_flat(nt, Nnodes_wet, 3))
    ALLOCATE(valid_j(Nnodes_wet), valid_k(Nnodes_wet))

    ! Build valid index maps and extract coordinates
    v = 1_I32
    do j = 1, Nm
        do k = 1, Nn
            global_idx = (j - 1_I32) * Nn + k
            if (mask_duplicate(global_idx) .and. self % xyz_all(j, k, 3) <= 0.0_WP) then
                valid_j(v)   = j
                valid_k(v)   = k
                x_flat(v, :) = self % xyz_all(j, k, :)
                v            = v + 1_I32
            end if
        end do
    end do

    self % xyz = x_flat

    ! Optimized flatten: temporal index i is the innermost to read RAM linearly
    !$OMP PARALLEL DO PRIVATE(v, j, k, i)
    do v = 1, Nnodes_wet
        j = valid_j(v)
        k = valid_k(v)
        !$OMP SIMD
        do i = 1, nt
            acc_flat(i, v, 1) = self % acc(i, j, k, 1)
            acc_flat(i, v, 2) = self % acc(i, j, k, 2)
            acc_flat(i, v, 3) = self % acc(i, j, k, 3)
        end do
    end do
    !$OMP END PARALLEL DO

    CALL self % compute_mass(mass, added_mass)
    ALLOCATE(mass_eff(Nnodes_wet))
    mass_eff = mass + added_mass
    DEALLOCATE(mass, added_mass)

    CALL compute_rfft(acc_flat, nt, dt, skipf=config % skipf, remove_zero=.true., &
                      array_out=A, freqs=freqs)
    self % Freqs = freqs

    ALLOCATE(self % Force(size(freqs), Nnodes_wet, 3))

    ! Vectorized force computation with OpenMP (implicit broadcast over frequencies)
    !$OMP PARALLEL DO PRIVATE(j)
    do j = 1, Nnodes_wet
        self % Force(:, j, 1) = A(:, j, 1) * mass_eff(j)
        self % Force(:, j, 2) = A(:, j, 2) * mass_eff(j)
        self % Force(:, j, 3) = A(:, j, 3) * mass_eff(j)
    end do
    !$OMP END PARALLEL DO

    ! Manual filtering using an explicit do-loop instead of the PACK intrinsic
    ! (avoids the extra temporary array reallocations PACK would need)
    if (config % filter_freqs) then
        mask_freqs_to_use = self % filter_frequencies(config % freqs_over)
        Nfreqs_new        = count(mask_freqs_to_use)

        ALLOCATE(freqs_new(Nfreqs_new))
        ALLOCATE(F_filtered(Nfreqs_new, Nnodes_wet, 3))

        v = 1_I32
        do i = 1, size(self % Freqs)
            if (mask_freqs_to_use(i)) then
                freqs_new(v)        = self % Freqs(i)
                F_filtered(v, :, :) = self % Force(i, :, :)
                v                   = v + 1_I32
            end if
        end do

        CALL MOVE_ALLOC(freqs_new, self % Freqs)
        CALL MOVE_ALLOC(F_filtered, self % Force)
    end if

    ! Cleanup
    DEALLOCATE(self % acc, self % Time, A, freqs, mass_eff)
    DEALLOCATE(mask_duplicate, valid_j, valid_k)
    if (ALLOCATED(mask_freqs_to_use)) DEALLOCATE(mask_freqs_to_use)

    ! Verbose summary
    if (verbose) then
        F_memory_mb = real(storage_size(self % Force), WP) * real(size(self % Force, kind=I64), WP) &
                      / 8.0_WP / 1024.0_WP / 1024.0_WP
        write(*, '(A)') repeat('=', 70)
        write(*, '(A)') 'FORCE COMPUTATION COMPLETED'
        write(*, '(A)') repeat('-', 68)
        write(*, '(A,I6)') 'Nodes (wet): ', Nnodes_wet
        write(*, '(A,I6)') 'Frequencies: ', size(self % Freqs)
        write(*, '(A,F10.2,A)') 'Memory (F):  ', F_memory_mb, ' MB'
        write(*, '(A)') repeat('=', 70)
    end if

END SUBROUTINE compute_force_MYSE3MWMonopile


SUBROUTINE compute_mass_MYSE3MWMonopile(self, mass, added_mass)
    ! -------------------------------------------------------------------------
    ! Computes the distributed dry (structural) and added (hydrodynamic) mass
    ! along the monopile nodes. Both diameters are constant with depth (wall
    ! thickness 60 mm, as the SubDyn members), the pile carries no internal
    ! water (as the OpenFAST model) and the added mass is that of the outer
    ! cross-section (Ca = 1, as the Morison members of HydroDyn). The nodal
    ! spans are evaluated at the real heights of the wetted nodes (self % xyz),
    ! so they do not depend on the water depth of Environment.nml.
    ! -------------------------------------------------------------------------

    USE StructureModel, ONLY: divide_span

    CLASS(MYSE3MWMonopile_t), INTENT(INOUT)            :: self             ! [-] Turbine instance
    REAL(WP)                , INTENT(OUT), ALLOCATABLE :: mass(:)          ! [kg] Dry mass per wetted node
    REAL(WP)                , INTENT(OUT), ALLOCATABLE :: added_mass(:)    ! [kg] Added mass per wetted node

    ! Local variables
    REAL(WP)                  :: A_mass, A_added_mass
    INTEGER(I32)              :: Nnodes
    REAL(WP)    , ALLOCATABLE :: ds(:)


    Nnodes = size(self % xyz, 1)

    ! Nodal spacings along the monopile axis
    ds = divide_span(self % xyz(:, 3))

    ! Cross-sectional areas
    A_mass       = PI * ((self % D / 2.0_WP)**2 - (self % D_inner / 2.0_WP)**2)
    A_added_mass = PI * (self % D / 2.0_WP)**2

    ! Distributed masses
    ALLOCATE(mass(Nnodes), added_mass(Nnodes))
    mass       = A_mass       * self % rho_mat * ds
    added_mass = A_added_mass * self % rho_wat * ds

    ! Cleanup
    DEALLOCATE(ds)

END SUBROUTINE compute_mass_MYSE3MWMonopile


SUBROUTINE excitation_spectrum_MYSE3MWMonopile(self, rpm, freqs_amp, keys)
    ! -------------------------------------------------------------------------
    ! Builds the MySE 3.0-135 drivetrain excitation spectrum (gear-mesh pure
    ! tones on the frequency grid of the force file of the run) at the given
    ! rotor speed, used by filter_frequencies.
    ! -------------------------------------------------------------------------

    CLASS(MYSE3MWMonopile_t), INTENT(IN)               :: self              ! [-] Turbine instance
    REAL(WP)                , INTENT(IN)               :: rpm               ! [rpm] Rotor speed
    REAL(WP)                , INTENT(OUT), ALLOCATABLE :: freqs_amp(:,:)    ! [-] Peaks (N, 2): freq [Hz], amplitude
    CHARACTER(len=20)       , INTENT(OUT), ALLOCATABLE :: keys(:)           ! [-] Peak labels (N)


    CALL drivetrain3MW_excitation_spectrum(rpm=rpm, T_win=self % T_win, freqs_amp=freqs_amp, keys=keys)

END SUBROUTINE excitation_spectrum_MYSE3MWMonopile


SUBROUTINE save_specific_params_MYSE3MWMonopile(self, save_path, istr)
    ! -------------------------------------------------------------------------
    ! Saves the monopile-specific geometric/material parameters (outer and
    ! inner diameters, material density, wetted area) into the solver's HDF5
    ! output file, tagging variable names with an optional instance suffix.
    ! -------------------------------------------------------------------------

    USE HDF5Utils, ONLY: save_variable

    CLASS(MYSE3MWMonopile_t), INTENT(IN)           :: self         ! [-] Turbine instance
    CHARACTER(len=*)        , INTENT(IN)           :: save_path    ! [-] Path to the HDF5 output file
    CHARACTER(len=*)        , INTENT(IN), OPTIONAL :: istr         ! [-] Suffix appended to the variable names


    if (present(istr)) then
        CALL save_variable(save_path, "Case_type_"       // istr, self % case_type)
        CALL save_variable(save_path, "Structure_nodes_" // istr, self % xyz_all)
        CALL save_variable(save_path, "D_"               // istr, self % D)
        CALL save_variable(save_path, "D_inner_"         // istr, self % D_inner)
        CALL save_variable(save_path, "rho_mat_"         // istr, self % rho_mat)
        CALL save_variable(save_path, "wet_area_"        // istr, self % wet_area)
    else
        CALL save_variable(save_path, "Case_type"      , self % case_type)
        CALL save_variable(save_path, "Structure_nodes", self % xyz_all)
        CALL save_variable(save_path, "D"              , self % D)
        CALL save_variable(save_path, "D_inner"        , self % D_inner)
        CALL save_variable(save_path, "rho_mat"        , self % rho_mat)
        CALL save_variable(save_path, "wet_area"       , self % wet_area)
    end if

END SUBROUTINE save_specific_params_MYSE3MWMonopile


SUBROUTINE write_specific_params_txt_MYSE3MWMonopile(self, unit)
    ! -------------------------------------------------------------------------
    ! Writes the monopile-specific parameters (outer and inner diameters,
    ! water/material density, wetted area) as human-readable lines into the
    ! given text unit.
    ! -------------------------------------------------------------------------

    CLASS(MYSE3MWMonopile_t), INTENT(IN) :: self    ! [-] Turbine instance
    INTEGER(I32)            , INTENT(IN) :: unit    ! [-] Unit of the open .dat file


    write(unit, '(A,F8.3,A)')  'Outer diameter   : ', self % D, ' m'
    write(unit, '(A,F8.3,A)')  'Inner diameter   : ', self % D_inner, ' m'
    write(unit, '(A,F8.1,A)')  'Water density    : ', self % rho_wat, ' kg/m^3'
    write(unit, '(A,F8.1,A)')  'Material density : ', self % rho_mat, ' kg/m^3'
    write(unit, '(A,F10.2,A)') 'Wetted area      : ', self % wet_area, ' m^2'

END SUBROUTINE write_specific_params_txt_MYSE3MWMonopile


FUNCTION get_impedance_corrected_force_MYSE3MWMonopile(self, c) RESULT(corrected_F)
    ! -------------------------------------------------------------------------
    ! Applies the cylindrical-scattering (Hankel-function) impedance correction
    ! for a monopile to the computed nodal forces, scaling each frequency by
    ! the magnitude of the scattering correction factor alpha(k, D).
    ! -------------------------------------------------------------------------

    USE StructureModel, ONLY: alpha_hankel

    CLASS(MYSE3MWMonopile_t), INTENT(IN)  :: self                  ! [-] Turbine instance
    REAL(WP)                , INTENT(IN)  :: c                     ! [m/s] Speed of sound in the fluid
    COMPLEX(WP)             , ALLOCATABLE :: corrected_F(:,:,:)    ! [N] Corrected force (Nfreqs, Nnodes_wet, 3)

    ! Local variables
    REAL(WP)    , ALLOCATABLE :: k(:)
    COMPLEX(WP) , ALLOCATABLE :: alpha(:)
    REAL(WP)    , ALLOCATABLE :: abs_alpha(:)
    INTEGER(I32)              :: dim1, dim2, dim3, m, j, i


    dim1 = size(self % Force, 1)
    dim2 = size(self % Force, 2)
    dim3 = size(self % Force, 3)

    ! Wave number vector
    ALLOCATE(k(dim1))
    k = (2.0_WP * PI * self % Freqs) / c

    ! Compute correction factor once
    alpha = alpha_hankel(k, self % D)

    ! Precompute absolute value of alpha (real), allocate temporary
    ALLOCATE(abs_alpha(dim1))
    abs_alpha = abs(alpha)

    ! Allocate corrected_F with explicit shape
    ALLOCATE(corrected_F(dim1, dim2, dim3))

    ! Parallelize over planes (j, m) and vectorize the inner loop over i
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i, j, m)
    do m = 1, dim3
        do j = 1, dim2
            !$OMP SIMD
            do i = 1, dim1
                corrected_F(i, j, m) = self % Force(i, j, m) * abs_alpha(i)
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    DEALLOCATE(abs_alpha)

END FUNCTION get_impedance_corrected_force_MYSE3MWMonopile

END MODULE MYSE3MWMonopile
