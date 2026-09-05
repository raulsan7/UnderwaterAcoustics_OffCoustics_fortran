MODULE DTU10MWMonopile

USE Kinds
USE Types
USE WindTurbine, ONLY: WindTurbine_t, drivetrain10MW_excitation_spectrum, init

IMPLICIT NONE

PRIVATE
PUBLIC :: DTU10MWMonopile_t



! ------------------------------
! Derived type definition
! ------------------------------
TYPE, EXTENDS(WindTurbine_t) :: DTU10MWMonopile_t
    ! Physical and structural parameters specific to the DTU 10 MW Monopile
    REAL(WP) :: D        = 9.0_WP      ! [m] Base monopile outer diameter
    REAL(WP) :: rho_wat  = 1025.0_WP   ! [kg/m^3] Seawater density
    REAL(WP) :: rho_mat  = 8500.0_WP   ! [kg/m^3] Material density (steel)
    REAL(WP) :: wet_area = 848.12_WP   ! [m^2] Wetted area


    CONTAINS
    ! --- Shared procedures --- !
    PROCEDURE :: init                          => init_DTU10MWMonopile

    ! --- Deferred(abstract) procedures --- !
    PROCEDURE :: compute_force                 => compute_force_DTU10MWMonopile
    PROCEDURE :: compute_mass                  => compute_mass_DTU10MWMonopile
    PROCEDURE :: filter_frequencies            => filter_frequencies_DTU10MWMonopile
    PROCEDURE :: save_specific_params          => save_specific_params_DTU10MWMonopile
    PROCEDURE :: get_impedance_corrected_force => get_impedance_corrected_force_DTU10MWMonopile

END TYPE DTU10MWMonopile_t


CONTAINS


SUBROUTINE init_DTU10MWMonopile(self, config)

    CLASS(DTU10MWMonopile_t)  , INTENT(INOUT) :: self
    TYPE(Config_WindTurbine_t), INTENT(IN)    :: config

    ! Local variables
    LOGICAL :: exists = .false.


    ! Define path to rpm curve HARCODED
    self % path_rpm = "../wind_speed_curves_DTU_10MW/rpm_ws.csv"

    ! Verify that the rpm curve file exists
    inquire(file= self % path_rpm, exist=exists)
    if (.not. exists) error stop "DTU10MWMonopile % init: RPM curve file not found: " // self % path_rpm

    ! Set the case type
    self % case_type = "Monopile"

    CALL init(self, config)

END SUBROUTINE init_DTU10MWMonopile


SUBROUTINE compute_force_DTU10MWMonopile(self, config, verbose)
    USE MathUtils, ONLY: remove_duplicate_nodes, compute_rfft
    
    CLASS(DTU10MWMonopile_t)  , INTENT(INOUT) :: self
    TYPE(Config_WindTurbine_t), INTENT(IN)    :: config
    LOGICAL                   , INTENT(IN)    :: verbose

    ! Local variables
    INTEGER(I32) :: Nm, Nn, nt, Nnodes_wet, i, j, k, v, global_idx, Nfreqs_new
    REAL(WP)     :: dt, F_memory_mb

    LOGICAL     , ALLOCATABLE :: mask_duplicate(:), mask_freqs_to_use(:)
    REAL(WP)    , ALLOCATABLE :: mass(:), added_mass(:), mass_effective(:)
    REAL(WP)    , ALLOCATABLE :: x_flat(:,:), acc_flat(:,:,:)
    COMPLEX(WP) , ALLOCATABLE :: A(:,:,:)
    REAL(WP)    , ALLOCATABLE :: freqs(:), freqs_new(:)
    COMPLEX(WP) , ALLOCATABLE :: F_new(:,:,:)
    INTEGER(I32), ALLOCATABLE :: valid_j(:), valid_k(:)


    Nm = self % Nmembers
    Nn = self % Nnodes
    nt = size(self % Time)
    dt = self % Time(2) - self % Time(1)

    mask_duplicate = remove_duplicate_nodes(self % x_all, delete_last=.true.)

    ! Fuse masks: count all valid nodes at once
    Nnodes_wet = 0
    do j = 1, Nm
        do k = 1, Nn
            global_idx = (j - 1) * Nn + k
            if (mask_duplicate(global_idx) .and. self % x_all(j, k, 3) <= 0.0_WP) then
                Nnodes_wet = Nnodes_wet + 1
            end if
        end do
    end do

    ALLOCATE(x_flat(Nnodes_wet, 3))
    ALLOCATE(acc_flat(nt, Nnodes_wet, 3))
    ALLOCATE(valid_j(Nnodes_wet), valid_k(Nnodes_wet))

    ! Build valid index maps and extract coordinates
    v = 1
    do j = 1, Nm
        do k = 1, Nn
            global_idx = (j - 1) * Nn + k
            if (mask_duplicate(global_idx) .and. self % x_all(j, k, 3) <= 0.0_WP) then
                valid_j(v) = j
                valid_k(v) = k
                x_flat(v, :) = self % x_all(j, k, :)
                v = v + 1
            end if
        end do
    end do

    self % x = x_flat
    
    ! Optimized flatten: temporal index i is the most intern to read RAM linearly
    !$omp parallel do private(v, j, k, i)
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
    !$omp end parallel do

    CALL self % compute_mass(mass, added_mass)
    ALLOCATE(mass_effective(Nnodes_wet))
    mass_effective = mass + added_mass
    DEALLOCATE(mass, added_mass)

    CALL compute_rfft(acc_flat, nt, dt, skipf=config%skip, remove_zero=.true., &
                      array_out=A, freqs=freqs)
    self % Freqs = freqs

    ALLOCATE(self % F(size(freqs), Nnodes_wet, 3))

    ! Vectorized force computation with OpenMP (implicit broadcast over frequencies)
    !$omp parallel do private(j)
    do j = 1, Nnodes_wet
        self % F(:, j, 1) = - A(:, j, 1) * mass_effective(j)
        self % F(:, j, 2) = - A(:, j, 2) * mass_effective(j)
        self % F(:, j, 3) = - A(:, j, 3) * mass_effective(j)
    end do
    !$omp end parallel do

    ! Manual filtering without 'RACK' to avoid memory reallocations
    if (config % filter_freqs) then
        mask_freqs_to_use = self % filter_frequencies()
        Nfreqs_new = count(mask_freqs_to_use)
        
        ALLOCATE(freqs_new(Nfreqs_new))
        ALLOCATE(F_new(Nfreqs_new, Nnodes_wet, 3))

        v = 1
        do i = 1, size(self%Freqs)
            if (mask_freqs_to_use(i)) then
                freqs_new(v) = self%Freqs(i)
                F_new(v, :, :) = self%F(i, :, :)
                v = v + 1
            end if
        end do
        
        CALL move_alloc(freqs_new, self%Freqs)
        CALL move_alloc(F_new, self%F)
    end if

    ! Cleanup
    DEALLOCATE(self%acc, self%Time, A, freqs, mass_effective)
    DEALLOCATE(mask_duplicate, valid_j, valid_k)
    if (ALLOCATED(mask_freqs_to_use)) DEALLOCATE(mask_freqs_to_use)

    ! Verbose summary
    if (verbose) then
        F_memory_mb = real(storage_size(self%F) * size(self%F) / 8) / 1024.0_WP / 1024.0_WP
        print '(A)', repeat('=', 70)
        print '(A)', 'FORCE COMPUTATION COMPLETED'
        print '(A)', repeat('-', 68)
        print '(A,I6)',   'Nodes (wet): ', Nnodes_wet
        print '(A,I6)',   'Frequencies: ', size(self%Freqs)
        print '(A,F6.2,A)', 'Memory (F):  ', F_memory_mb, ' MB'
        print '(A)', repeat('=', 70)
    end if

END SUBROUTINE compute_force_DTU10MWMonopile


SUBROUTINE compute_mass_DTU10MWMonopile(self, mass, added_mass)
    USE MathUtils, ONLY: divide_span

    CLASS(DTU10MWMonopile_t), INTENT(INOUT) :: self
    REAL(WP), ALLOCATABLE   , INTENT(OUT)    :: mass(:)
    REAL(WP), ALLOCATABLE   , INTENT(OUT)    :: added_mass(:)

    ! Local variables
    REAL(WP) :: D, Depth, rho_mat, rho_wat
    INTEGER(I32) :: Nnodes, i
    REAL(WP), ALLOCATABLE :: z_nodes(:), ds(:), D_outer(:), D_inner(:)
    REAL(WP), ALLOCATABLE :: A_mass(:), A_added_mass(:)

    ! Retrieve parameters
    D       = self%D
    Depth   = self%Depth
    rho_mat = self%rho_mat
    rho_wat = self%rho_wat
    Nnodes  = size(self%x, 1)

    ! Generate vertical coordinate array linearly from -Depth to 0
    allocate(z_nodes(Nnodes))
    do i = 1, Nnodes
        z_nodes(i) = -Depth + real(i-1, WP) * Depth / real(Nnodes-1, WP)
    end do

    ! Nodal spacings
    ds = divide_span(z_nodes)

    ! Outer diameter constant
    allocate(D_outer(Nnodes))
    D_outer = D

    ! Inner diameter depending on elevation
    allocate(D_inner(Nnodes))
    do i = 1, Nnodes
        if (z_nodes(i) >= 4.0_WP .and. z_nodes(i) <= 10.0_WP) then
            D_inner(i) = 8.70_WP
        else if (z_nodes(i) >= -10.0_WP .and. z_nodes(i) < 4.0_WP) then
            D_inner(i) = 8.69_WP
        else if (z_nodes(i) < -10.0_WP) then
            D_inner(i) = 8.80_WP
        else
            D_inner(i) = 0.0_WP  ! should not happen
        end if
    end do

    ! Cross-sectional areas
    allocate(A_mass(Nnodes), A_added_mass(Nnodes))
    A_mass       = PI * ((D_outer/2.0_WP)**2 - (D_inner/2.0_WP)**2)
    A_added_mass = PI * (D_outer/2.0_WP)**2

    ! Distributed masses
    allocate(mass(Nnodes), added_mass(Nnodes))
    mass       = A_mass       * rho_mat * ds
    added_mass = A_added_mass * rho_wat * ds

    ! Cleanup local allocatables (optional, automatic on return)
    deallocate(z_nodes, ds, D_outer, D_inner, A_mass, A_added_mass)

END SUBROUTINE compute_mass_DTU10MWMonopile


FUNCTION filter_frequencies_DTU10MWMonopile(self, freqs_over) RESULT(mask)
    USE IOUtils, ONLY: read_curve
    USE MathUtils, ONLY: generate_timeseries_banded_sines, filter_non_usefull_freqs

    CLASS(DTU10MWMonopile_t), INTENT(INOUT)        :: self
    REAL(WP)                , INTENT(IN), OPTIONAL :: freqs_over 
    LOGICAL, ALLOCATABLE                           :: mask(:)

    ! Local variables
    CHARACTER(len=20), ALLOCATABLE :: keys(:)
    REAL(WP)         , ALLOCATABLE :: freqs_amp(:,:)
    REAL(WP)         , ALLOCATABLE :: freqs_to_use(:)
    REAL(WP)         , ALLOCATABLE :: tmp_signal(:)
    REAL(WP)                       :: rpm, freqs_over_

    freqs_over_ = 10.0_WP; if (PRESENT(freqs_over)) freqs_over_ = freqs_over

    ! Read rpm at current wind speed
    rpm = read_curve(self % path_rpm, self % WindSpeed)

    ! Generate drivetrain excitation spectrum
    CALL drivetrain10MW_excitation_spectrum(rpm=rpm, freqs_amp=freqs_amp, keys=keys)

    ! Get allowed frequencies
    CALL generate_timeseries_banded_sines(peaks=freqs_amp, keys=keys, t=self % Time, &
                                         zeta=0.02_WP, used_freqs=.true., a=tmp_signal, &
                                         freqs_out=freqs_to_use)

    ! Build mask
    mask = filter_non_usefull_freqs(self % Freqs, freqs_to_use, freqs_over_)

    DEALLOCATE(tmp_signal)

    DEALLOCATE(freqs_amp, keys, freqs_to_use)

END FUNCTION filter_frequencies_DTU10MWMonopile


SUBROUTINE save_specific_params_DTU10MWMonopile(self, save_path, ichar)
    USE IOUtils, ONLY: save_to_hdf5

    CLASS(DTU10MWMonopile_t), INTENT(IN) :: self
    CHARACTER(len=*)        , INTENT(IN) :: save_path
    CHARACTER(len=*)        , INTENT(IN), OPTIONAL :: ichar
    
    if (PRESENT(ichar)) then
        CALL save_to_hdf5(save_path, "D_"        // ichar, self % D)
        CALL save_to_hdf5(save_path, "rho_mat_"  // ichar, self % rho_mat)
        CALL save_to_hdf5(save_path, "wet_area_" // ichar, self % wet_area)
    else
        CALL save_to_hdf5(save_path, "D", self % D)
        CALL save_to_hdf5(save_path, "rho_mat", self % rho_mat)
        CALL save_to_hdf5(save_path, "wet_area", self % wet_area)
    end if

END SUBROUTINE save_specific_params_DTU10MWMonopile


FUNCTION get_impedance_corrected_force_DTU10MWMonopile(self, c) RESULT(corrected_F)
    USE MathUtils, ONLY: alpha_hankel

    CLASS(DTU10MWMonopile_t), INTENT(IN)  :: self
    REAL(WP)                , INTENT(IN)  :: c
    COMPLEX(WP)             , ALLOCATABLE :: corrected_F(:,:,:)

    ! Local variables
    REAL(WP)   , ALLOCATABLE :: k(:)
    COMPLEX(WP), ALLOCATABLE :: alpha(:)
    REAL(WP)   , ALLOCATABLE :: abs_alpha(:)
    INTEGER(I32)             :: dim1, dim2, dim3, m, j, i

    
    dim1 = size(self % F, 1)
    dim2 = size(self % F, 2)
    dim3 = size(self % F, 3)
    
    ! Wave number vector
    ALLOCATE(k(dim1))
    k = (2.0_WP * PI * self % Freqs) / c

    ! Compute correction factor once
    alpha = alpha_hankel(k, self % D)

    ! Precompute absolute value of alpha (real), allocate temporary
    ALLOCATE(abs_alpha(dim1))
    abs_alpha = ABS(alpha)

    ! Allocate corrected_F with explicit shape
    ALLOCATE(corrected_F(dim1, dim2, dim3))

    ! Parallelize over planes (j,m) and vectorize the inner loop over i
    !$omp parallel do collapse(2) default(shared) private(i,j,m)
    do m = 1, dim3
        do j = 1, dim2
            !$omp simd
            do i = 1, dim1
                corrected_F(i, j, m) = self%F(i, j, m) * abs_alpha(i)
            end do
        end do
    end do

    DEALLOCATE(abs_alpha)

END FUNCTION get_impedance_corrected_force_DTU10MWMonopile



END MODULE DTU10MWMonopile