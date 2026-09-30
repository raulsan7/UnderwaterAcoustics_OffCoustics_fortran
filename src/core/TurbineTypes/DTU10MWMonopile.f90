! ===============================================================================
! MODULE: DTU10MWMonopile
! ===============================================================================
! Defines DTU10MWMonopile_t, a concrete turbine type extending the abstract
! WindTurbine_t (src/core/WindTurbine.f90) for the DTU 10 MW monopile-founded
! bottom-fixed platform. It implements the deferred procedures: force
! computation from OpenFAST accelerations, distributed dry/added mass along the
! monopile, drivetrain-based frequency filtering, saving/writing of platform-
! specific parameters, and the cylindrical-scattering impedance correction.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE DTU10MWMonopile

USE Kinds
USE Types, ONLY: Config_WindTurbine_t
USE WindTurbine, ONLY: WindTurbine_t, drivetrain10MW_excitation_spectrum, init

IMPLICIT NONE

PRIVATE
PUBLIC :: DTU10MWMonopile_t


TYPE, EXTENDS(WindTurbine_t) :: DTU10MWMonopile_t
    ! -------------------------------------------------------------------------
    ! Concrete turbine type extending WindTurbine_t for the DTU 10 MW monopile-
    ! founded bottom-fixed offshore wind turbine. Stores the monopile outer
    ! diameter, material/water densities and wetted area, and implements every
    ! deferred procedure of the abstract parent (force, mass, frequency filtering,
    ! parameter I/O and impedance correction).
    ! -------------------------------------------------------------------------

    ! Physical and structural parameters specific to the DTU 10 MW Monopile
    REAL(WP) :: D        = 9.0_WP       ! [m] Base monopile outer diameter
    REAL(WP) :: rho_wat  = 1025.0_WP    ! [kg/m^3] Seawater density
    REAL(WP) :: rho_mat  = 8500.0_WP    ! [kg/m^3] Material density (steel)
    REAL(WP) :: wet_area = 848.12_WP    ! [m^2] Wetted area

    CONTAINS

    ! --- Shared procedures --- !
    PROCEDURE :: init => init_DTU10MWMonopile

    ! --- Deferred (abstract) procedures --- !
    PROCEDURE :: compute_force                 => compute_force_DTU10MWMonopile
    PROCEDURE :: compute_mass                  => compute_mass_DTU10MWMonopile
    PROCEDURE :: filter_frequencies            => filter_frequencies_DTU10MWMonopile
    PROCEDURE :: save_specific_params          => save_specific_params_DTU10MWMonopile
    PROCEDURE :: write_specific_params_txt     => write_specific_params_txt_DTU10MWMonopile
    PROCEDURE :: get_impedance_corrected_force => get_impedance_corrected_force_DTU10MWMonopile

END TYPE DTU10MWMonopile_t


CONTAINS


SUBROUTINE init_DTU10MWMonopile(self, config)
    ! -------------------------------------------------------------------------
    ! Sets the monopile-specific (hardcoded) rpm-curve path, marks the turbine
    ! as monopile-founded, and delegates to the base WindTurbine_t % init for
    ! the remaining shared initialization steps.
    ! -------------------------------------------------------------------------

    CLASS(DTU10MWMonopile_t)  , INTENT(INOUT) :: self      ! [-] Turbine instance
    TYPE(Config_WindTurbine_t), INTENT(IN)    :: config    ! [-] Turbine configuration

    ! Local variables
    LOGICAL :: exists


    ! Path to the rpm curve file (hardcoded)
    self % path_rpm = "../wind_speed_curves_DTU_10MW/rpm_ws.csv"

    ! Verify that the rpm curve file exists
    inquire(file=self % path_rpm, exist=exists)
    if (.not. exists) error stop "DTU10MWMonopile % init: RPM curve file not found: " // self % path_rpm

    ! Flag the turbine as monopile
    self % is_monopile = .true.

    CALL init(self, config)

END SUBROUTINE init_DTU10MWMonopile


SUBROUTINE compute_force_DTU10MWMonopile(self, config, verbose)
    ! -------------------------------------------------------------------------
    ! Computes the nodal structural forces from OpenFAST accelerations: removes
    ! duplicate/dry (above-water) nodes, computes the real FFT of the wet-node
    ! accelerations, and multiplies by the effective (dry + added) mass to get
    ! F = -m*a in the frequency domain. Optionally keeps only the frequencies
    ! relevant to the drivetrain excitation.
    ! -------------------------------------------------------------------------

    USE MathUtils, ONLY: remove_duplicate_nodes, compute_rfft

    CLASS(DTU10MWMonopile_t)  , INTENT(INOUT) :: self       ! [-] Turbine instance
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
        self % Force(:, j, 1) = -A(:, j, 1) * mass_eff(j)
        self % Force(:, j, 2) = -A(:, j, 2) * mass_eff(j)
        self % Force(:, j, 3) = -A(:, j, 3) * mass_eff(j)
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
        F_memory_mb = real(storage_size(self % Force) * size(self % Force) / 8) / 1024.0_WP / 1024.0_WP
        write(*, '(A)') repeat('=', 70)
        write(*, '(A)') 'FORCE COMPUTATION COMPLETED'
        write(*, '(A)') repeat('-', 68)
        write(*, '(A,I6)') 'Nodes (wet): ', Nnodes_wet
        write(*, '(A,I6)') 'Frequencies: ', size(self % Freqs)
        write(*, '(A,F6.2,A)') 'Memory (F):  ', F_memory_mb, ' MB'
        write(*, '(A)') repeat('=', 70)
    end if

END SUBROUTINE compute_force_DTU10MWMonopile


SUBROUTINE compute_mass_DTU10MWMonopile(self, mass, added_mass)
    ! -------------------------------------------------------------------------
    ! Computes the distributed dry (structural) and added (hydrodynamic) mass
    ! along the monopile nodes. The outer diameter is constant with depth, but
    ! the inner diameter varies across three depth bands to represent the real
    ! wall-thickness profile: 8.70 m between +4 m and +10 m, 8.69 m between
    ! -10 m and +4 m, and 8.80 m below -10 m.
    ! -------------------------------------------------------------------------

    USE MathUtils, ONLY: divide_span

    CLASS(DTU10MWMonopile_t), INTENT(INOUT)            :: self             ! [-] Turbine instance
    REAL(WP)                , INTENT(OUT), ALLOCATABLE :: mass(:)          ! [kg] Dry mass per wetted node
    REAL(WP)                , INTENT(OUT), ALLOCATABLE :: added_mass(:)    ! [kg] Added mass per wetted node

    ! Local variables
    REAL(WP)                  :: D, Depth, rho_mat, rho_wat
    INTEGER(I32)              :: Nnodes, i
    REAL(WP)    , ALLOCATABLE :: z_nodes(:), ds(:), D_outer(:), D_inner(:)
    REAL(WP)    , ALLOCATABLE :: A_mass(:), A_added_mass(:)


    ! Retrieve parameters
    D       = self % D
    Depth   = self % Depth
    rho_mat = self % rho_mat
    rho_wat = self % rho_wat
    Nnodes  = size(self % xyz, 1)

    ! Generate vertical coordinate array linearly from -Depth to 0
    ALLOCATE(z_nodes(Nnodes))
    do i = 1, Nnodes
        z_nodes(i) = -Depth + real(i - 1_I32, WP) * Depth / real(Nnodes - 1_I32, WP)
    end do

    ! Nodal spacings
    ds = divide_span(z_nodes)

    ! Outer diameter constant
    ALLOCATE(D_outer(Nnodes))
    D_outer = D

    ! Inner diameter depending on elevation
    ALLOCATE(D_inner(Nnodes))
    do i = 1, Nnodes
        if (z_nodes(i) >= 4.0_WP .and. z_nodes(i) <= 10.0_WP) then
            D_inner(i) = 8.70_WP
        else if (z_nodes(i) >= -10.0_WP .and. z_nodes(i) < 4.0_WP) then
            D_inner(i) = 8.69_WP
        else if (z_nodes(i) < -10.0_WP) then
            D_inner(i) = 8.80_WP
        else
            D_inner(i) = 0.0_WP    ! Should not happen
        end if
    end do

    ! Cross-sectional areas
    ALLOCATE(A_mass(Nnodes), A_added_mass(Nnodes))
    A_mass       = PI * ((D_outer / 2.0_WP)**2 - (D_inner / 2.0_WP)**2)
    A_added_mass = PI * (D_outer / 2.0_WP)**2

    ! Distributed masses
    ALLOCATE(mass(Nnodes), added_mass(Nnodes))
    mass       = A_mass       * rho_mat * ds
    added_mass = A_added_mass * rho_wat * ds

    ! Cleanup
    DEALLOCATE(z_nodes, ds, D_outer, D_inner, A_mass, A_added_mass)

END SUBROUTINE compute_mass_DTU10MWMonopile


FUNCTION filter_frequencies_DTU10MWMonopile(self, freqs_over) RESULT(mask)
    ! -------------------------------------------------------------------------
    ! Builds a logical mask selecting only the frequencies close to the DTU
    ! 10 MW drivetrain excitation peaks (shaft/gear-mesh harmonics) at the
    ! current wind speed's rpm, so the force spectrum can be restricted to the
    ! physically relevant frequencies.
    ! -------------------------------------------------------------------------

    USE IOUtils, ONLY: read_curve
    USE MathUtils, ONLY: generate_timeseries_banded_sines, filter_non_useful_freqs

    CLASS(DTU10MWMonopile_t), INTENT(INOUT)        :: self          ! [-] Turbine instance
    REAL(WP)                , INTENT(IN), OPTIONAL :: freqs_over    ! [Hz] Frequencies below this are never filtered
    LOGICAL                 , ALLOCATABLE          :: mask(:)       ! [-] Frequency mask (Nfreqs), .true. to keep

    ! Local variables
    CHARACTER(len=20), ALLOCATABLE :: keys(:)
    REAL(WP)         , ALLOCATABLE :: freqs_amp(:,:)
    REAL(WP)         , ALLOCATABLE :: freqs_to_use(:)
    REAL(WP)         , ALLOCATABLE :: tmp_signal(:)
    REAL(WP)                       :: rpm, freqs_over_


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

END FUNCTION filter_frequencies_DTU10MWMonopile


SUBROUTINE save_specific_params_DTU10MWMonopile(self, save_path, istr)
    ! -------------------------------------------------------------------------
    ! Saves the monopile-specific geometric/material parameters (outer diameter,
    ! material density, wetted area) into the solver's HDF5 output file, tagging
    ! variable names with an optional instance suffix.
    ! -------------------------------------------------------------------------

    USE IOUtils, ONLY: save_variable

    CLASS(DTU10MWMonopile_t), INTENT(IN)           :: self         ! [-] Turbine instance
    CHARACTER(len=*)        , INTENT(IN)           :: save_path    ! [-] Path to the HDF5 output file
    CHARACTER(len=*)        , INTENT(IN), OPTIONAL :: istr         ! [-] Suffix appended to the variable names


    if (present(istr)) then
        CALL save_variable(save_path, "Case_type_"       // istr, self % case_type)
        CALL save_variable(save_path, "Structure_nodes_" // istr, self % xyz_all)
        CALL save_variable(save_path, "D_"               // istr, self % D)
        CALL save_variable(save_path, "rho_mat_"         // istr, self % rho_mat)
        CALL save_variable(save_path, "wet_area_"        // istr, self % wet_area)
    else
        CALL save_variable(save_path, "Case_type"      , self % case_type)
        CALL save_variable(save_path, "Structure_nodes", self % xyz_all)
        CALL save_variable(save_path, "D"              , self % D)
        CALL save_variable(save_path, "rho_mat"        , self % rho_mat)
        CALL save_variable(save_path, "wet_area"       , self % wet_area)
    end if

END SUBROUTINE save_specific_params_DTU10MWMonopile


SUBROUTINE write_specific_params_txt_DTU10MWMonopile(self, unit)
    ! -------------------------------------------------------------------------
    ! Writes the monopile-specific parameters (outer diameter, water/material
    ! density, wetted area) as human-readable lines into the given text unit.
    ! -------------------------------------------------------------------------

    CLASS(DTU10MWMonopile_t), INTENT(IN) :: self    ! [-] Turbine instance
    INTEGER(I32)            , INTENT(IN) :: unit    ! [-] Unit of the open .dat file


    write(unit, '(A,F8.3,A)')  'Outer diameter   : ', self % D, ' m'
    write(unit, '(A,F8.1,A)')  'Water density    : ', self % rho_wat, ' kg/m^3'
    write(unit, '(A,F8.1,A)')  'Material density : ', self % rho_mat, ' kg/m^3'
    write(unit, '(A,F10.2,A)') 'Wetted area      : ', self % wet_area, ' m^2'

END SUBROUTINE write_specific_params_txt_DTU10MWMonopile


FUNCTION get_impedance_corrected_force_DTU10MWMonopile(self, c) RESULT(corrected_F)
    ! -------------------------------------------------------------------------
    ! Applies the cylindrical-scattering (Hankel-function) impedance correction
    ! for a monopile to the computed nodal forces, scaling each frequency by
    ! the magnitude of the scattering correction factor alpha(k, D).
    ! -------------------------------------------------------------------------

    USE MathUtils, ONLY: alpha_hankel

    CLASS(DTU10MWMonopile_t), INTENT(IN)  :: self                  ! [-] Turbine instance
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

END FUNCTION get_impedance_corrected_force_DTU10MWMonopile

END MODULE DTU10MWMonopile
