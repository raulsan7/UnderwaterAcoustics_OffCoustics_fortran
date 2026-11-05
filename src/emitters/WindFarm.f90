! ===============================================================================
! MODULE: WindFarm
! ===============================================================================
! Defines WindFarm_t, an extension of Emitter_t that contains wind turbines
! of any type placed at arbitrary positions and operating points. It builds
! every turbine from its own configuration, computes their forces one after
! another and merges all of them into a single AcousticSource_t, relying on
! the linearity of the Helmholtz equation: the farm field is the coherent
! superposition of the fields radiated by every turbine.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 09/2026
! ===============================================================================
MODULE WindFarm

USE Kinds
USE Emitter, ONLY: Emitter_t
USE StringUtils, ONLY: to_upper
USE HDF5Utils, ONLY: save_variable
USE TurbineFactory, ONLY: new_turbine
USE SignalProcessing, ONLY: HANN_ENBW
USE ConfigTypes, ONLY: Config_WindFarm_t, Config_t
USE SourceTypes, ONLY: AcousticSource_t, SOURCE_DIPOLE
USE WindTurbine, ONLY: WindTurbine_t, save_turbine_params, write_turbine_txt

IMPLICIT NONE

PRIVATE
PUBLIC :: WindFarm_t


TYPE :: TurbineSlot_t
    ! -------------------------------------------------------------------------
    ! Wrapper that allows an array to hold turbines of different dynamic
    ! types, since Fortran has no arrays of polymorphic objects.
    ! -------------------------------------------------------------------------

    CLASS(WindTurbine_t), ALLOCATABLE :: wt    ! [-] Polymorphic turbine instance

END TYPE TurbineSlot_t


TYPE, EXTENDS(Emitter_t) :: WindFarm_t
    ! -------------------------------------------------------------------------
    ! Container of heterogeneous wind turbines. Turbine i owns the sources
    ! src_offset(i)+1 : src_offset(i+1) of the merged AcousticSource_t, so the
    ! contribution of every turbine can still be told apart after merging.
    ! -------------------------------------------------------------------------

    INTEGER(I32)                     :: Nturbines               ! [-] Number of turbines in the farm
    CHARACTER(len=100)               :: farm_name               ! [-] Label replacing the automatic file name ('None': auto)
    TYPE(TurbineSlot_t), ALLOCATABLE :: turbines(:)             ! [-] Turbines (Nturbines)
    INTEGER(I32)       , ALLOCATABLE :: src_offset(:)           ! [-] Source offset of every turbine (Nturbines + 1)

    CONTAINS

    ! --- Public type-bound procedures --- !
    PROCEDURE :: init           => init_WindFarm                ! Allocates and initializes every turbine
    PROCEDURE :: compute_forces => compute_forces_WindFarm      ! Reads OpenFAST outputs and computes forces
    PROCEDURE :: to_source      => farm_to_source_WindFarm      ! Merges every turbine into one source set
    PROCEDURE :: save_params    => save_params_WindFarm         ! Farm and per-turbine parameters -> HDF5
    PROCEDURE :: write_txt      => write_txt_WindFarm           ! Farm and per-turbine blocks of the .dat
    PROCEDURE :: get_enbw       => get_enbw_WindFarm            ! Hann window of compute_rfft (HANN_ENBW)

END TYPE WindFarm_t


CONTAINS


SUBROUTINE init_WindFarm(self, config)
    ! -------------------------------------------------------------------------
    ! Allocates every turbine with the concrete type given by its case_type
    ! and runs its constructor (paths, rpm curve, shared parameters).
    ! -------------------------------------------------------------------------

    CLASS(WindFarm_t)      , INTENT(INOUT) :: self      ! [-] Farm instance
    TYPE(Config_WindFarm_t), INTENT(INOUT) :: config    ! [-] Farm configuration

    ! Local variables
    INTEGER(I32) :: i


    self % Nturbines = config % Nturbines
    self % farm_name = config % farm_name
    ALLOCATE(self % turbines(self % Nturbines))

    do i = 1, self % Nturbines
        CALL new_turbine(config % turbines(i) % case_type, self % turbines(i) % wt)
        CALL self % turbines(i) % wt % init(config=config % turbines(i))
    end do

    write(*, '(A,I0,A)') ' >> [WindFarm] Initialized with ', self % Nturbines, ' turbines'

END SUBROUTINE init_WindFarm


SUBROUTINE compute_forces_WindFarm(self, config, verbose)
    ! -------------------------------------------------------------------------
    ! Reads the OpenFAST outputs and computes the nodal forces of every
    ! turbine, one after another, so that only one raw acceleration array
    ! lives in memory at a time (compute_force frees it).
    ! -------------------------------------------------------------------------

    CLASS(WindFarm_t)      , INTENT(INOUT) :: self       ! [-] Farm instance
    TYPE(Config_WindFarm_t), INTENT(INOUT) :: config     ! [-] Farm configuration
    LOGICAL                , INTENT(IN)    :: verbose    ! [-] Whether to print a summary per turbine

    ! Local variables
    INTEGER(I32) :: i


    do i = 1, self % Nturbines
        write(*, '(A,I0,A,I0,A)') ' >> [WindFarm] Turbine ', i, '/', self % Nturbines, &
                                  ': ' // trim(config % turbines(i) % rootname)
        CALL self % turbines(i) % wt % read_input(config=config % turbines(i), verbose=verbose)
        CALL self % turbines(i) % wt % compute_force(config=config % turbines(i), verbose=verbose)
    end do

END SUBROUTINE compute_forces_WindFarm


SUBROUTINE farm_to_source_WindFarm(self, config, sources)
    ! -------------------------------------------------------------------------
    ! Merges every turbine into one AcousticSource_t. Sources are stacked
    ! turbine after turbine (src_offset keeps track of them). Frequencies are
    ! merged on the union of the common FFT bins k = f/df: a turbine that has
    ! no bin k (e.g. removed by its own drivetrain filter) gets a zero force
    ! there. Impedance corrections are applied per turbine, through each
    ! turbine's own turbine_to_source.
    ! -------------------------------------------------------------------------

    CLASS(WindFarm_t)     , INTENT(INOUT)            :: self       ! [-] Farm instance
    TYPE(Config_t)        , INTENT(IN)               :: config     ! [-] Global configuration
    TYPE(AcousticSource_t), INTENT(OUT), ALLOCATABLE :: sources    ! [-] Farm acoustic point sources

    ! Local variables
    TYPE(AcousticSource_t), ALLOCATABLE :: src_i
    LOGICAL               , ALLOCATABLE :: bin_used(:)
    INTEGER(I32)          , ALLOCATABLE :: bin_pos(:), pos_i(:)
    INTEGER(I32)                        :: i, j, k, k_max, nf, Ns, n_i, off, m, s
    REAL(WP)                            :: df, sum_cos, sum_sin


    ! --- 1. Common FFT grid: every turbine must share df = 1/(nt*dt) --- !
    df = self % turbines(1) % wt % df
    do i = 2, self % Nturbines
        if (abs(self % turbines(i) % wt % df - df) > 1.0e-9_WP * df) then
            write(*, '(A,I0,A,ES14.6,A,ES14.6,A)') ' ERROR WindFarm % farm_to_source: turbine ', i, &
                                                   ' has df = ', self % turbines(i) % wt % df, &
                                                   ' Hz, turbine 1 has df = ', df, ' Hz'
            error stop
        end if
    end do

    ! --- 2. Union of frequency bins (integer indices, no float comparisons) --- !
    k_max = 0_I32
    do i = 1, self % Nturbines
        k_max = max(k_max, nint(maxval(self % turbines(i) % wt % Freqs) / df, I32))
    end do

    ALLOCATE(bin_used(0:k_max), bin_pos(0:k_max))
    bin_used = .false.
    do i = 1, self % Nturbines
        bin_used(nint(self % turbines(i) % wt % Freqs / df, I32)) = .true.
    end do

    nf      = count(bin_used)
    bin_pos = 0_I32
    j       = 0_I32
    do k = 0, k_max
        if (bin_used(k)) then
            j          = j + 1_I32
            bin_pos(k) = j
        end if
    end do

    ! --- 3. Source offsets and farm arrays --- !
    if (ALLOCATED(self % src_offset)) DEALLOCATE(self % src_offset)
    ALLOCATE(self % src_offset(self % Nturbines + 1))
    self % src_offset(1) = 0_I32
    do i = 1, self % Nturbines
        self % src_offset(i+1) = self % src_offset(i) + size(self % turbines(i) % wt % xyz, 1)
    end do
    Ns = self % src_offset(self % Nturbines + 1)

    ALLOCATE(sources)
    ALLOCATE(sources % xyz (Ns, 3), sources % Force(nf, Ns, 3), sources % Freqs(nf))
    sources % Force = (0.0_WP, 0.0_WP)
    ALLOCATE(sources % order(Ns), source=SOURCE_DIPOLE)    ! Every source radiates as a point dipole

    ! --- 4. Scatter every turbine into the farm arrays --- !
    do i = 1, self % Nturbines
        CALL self % turbines(i) % wt % to_source(config, src_i)

        off   = self % src_offset(i)
        n_i   = size(src_i % xyz, 1)
        pos_i = bin_pos(nint(src_i % Freqs / df, I32))

        sources % xyz(off+1:off+n_i, :) = src_i % xyz
        sources % Freqs(pos_i)          = src_i % Freqs    ! Copied, not recomputed as k*df

        !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(m, s, j)
        do m = 1, 3
            do s = 1, n_i
                do j = 1, size(pos_i)
                    sources % Force(pos_i(j), off+s, m) = src_i % Force(j, s, m)
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        if (size(pos_i) < nf) write(*, '(A,I0,A,I0,A)') ' WARNING WindFarm % farm_to_source: turbine ', i, &
                                                         ' has zero force in ', nf - size(pos_i), ' frequency bins'

        DEALLOCATE(src_i, pos_i)    ! The turbine forces were already consumed by turbine_to_source
    end do

    ! --- 5. Reference geometry and default file name labels --- !
    sum_cos                     = 0.0_WP
    sum_sin                     = 0.0_WP
    sources % geometry % center = 0.0_WP
    sources % geometry % origin = 0.0_WP

    do i = 1, self % Nturbines
        sources % geometry % center = sources % geometry % center + self % turbines(i) % wt % BariPos
        sources % geometry % origin = sources % geometry % origin + self % turbines(i) % wt % AxisPos

        sum_cos = sum_cos + cos(self % turbines(i) % wt % WindDir * PI / 180.0_WP)
        sum_sin = sum_sin + sin(self % turbines(i) % wt % WindDir * PI / 180.0_WP)
    end do

    sources % geometry % center    = sources % geometry % center / real(self % Nturbines, WP)
    sources % geometry % origin    = sources % geometry % origin / real(self % Nturbines, WP)
    sources % geometry % direction = atan2(sum_sin, sum_cos) * 180.0_WP / PI    ! Circular mean

    ! Default output folder and file name labels (the solver inserts its method tag in between)
    sources % geometry % out_dir = './farm_acoustic_data/'
    CALL build_default_name(self, sources % geometry % name_tag, sources % geometry % name_suffix)

    write(*, '(A,I0,A,I0,A)') ' >> [WindFarm] ', Ns, ' sources merged over ', nf, ' frequencies'

END SUBROUTINE farm_to_source_WindFarm


SUBROUTINE save_params_WindFarm(self, save_path)
    ! -------------------------------------------------------------------------
    ! Saves the farm layout and the generic and type-specific parameters of
    ! every turbine into the HDF5 output file. The datasets of turbine i
    ! carry the suffix "_i" (e.g. WindSpeed_2, Structure_nodes_2).
    ! -------------------------------------------------------------------------

    CLASS(WindFarm_t), INTENT(IN) :: self         ! [-] Farm instance
    CHARACTER(len=*) , INTENT(IN) :: save_path    ! [-] Path to the HDF5 output file

    ! Local variables
    INTEGER(I32)      :: i
    CHARACTER(len=16) :: istr


    CALL save_variable(save_path, "Nturbines"    , self % Nturbines)
    CALL save_variable(save_path, "Source_offset", self % src_offset)

    do i = 1, self % Nturbines
        write(istr, '(I0)') i
        CALL save_turbine_params(self % turbines(i) % wt, save_path, istr=trim(istr))
        CALL self % turbines(i) % wt % save_specific_params(save_path, istr=trim(istr))
    end do

END SUBROUTINE save_params_WindFarm


SUBROUTINE write_txt_WindFarm(self, unit)
    ! -------------------------------------------------------------------------
    ! Writes the farm block of the .dat summary (layout table) followed by
    ! the generic and type-specific block of every turbine.
    ! -------------------------------------------------------------------------

    CLASS(WindFarm_t), INTENT(IN) :: self    ! [-] Farm instance
    INTEGER(I32)     , INTENT(IN) :: unit    ! [-] Unit of the open .dat file

    ! Local variables
    INTEGER(I32) :: i


    write(unit, '(A)')    repeat('-', 70)
    write(unit, '(A)')    'WIND FARM'
    write(unit, '(A)')    repeat('-', 70)
    write(unit, '(A,I0)') 'N turbines  : ', self % Nturbines
    write(unit, '(A)')    '   #  Case type            WS [m/s]  Dir [deg]      x [m]      y [m]  Sources'
    do i = 1, self % Nturbines
        associate(wt => self % turbines(i) % wt)
            write(unit, '(I4,2X,A20,F9.2,F11.2,2F11.2,I9)') i, wt % case_type, wt % WindSpeed, wt % WindDir, &
                                                            wt % AxisPos, self % src_offset(i+1) - self % src_offset(i)
        end associate
    end do
    write(unit, *)

    do i = 1, self % Nturbines
        write(unit, '(A,I0)') '>>> Turbine ', i
        CALL write_turbine_txt(self % turbines(i) % wt, unit)
    end do

END SUBROUTINE write_txt_WindFarm


FUNCTION get_enbw_WindFarm(self) RESULT(enbw)
    ! -------------------------------------------------------------------------
    ! Emitter_t % get_enbw of a wind farm: the forces of every turbine come
    ! from compute_rfft, so their window is the Hann window (HANN_ENBW).
    ! -------------------------------------------------------------------------

    CLASS(WindFarm_t), INTENT(IN) :: self    ! [-] Farm instance
    REAL(WP)                      :: enbw    ! [-] Equivalent noise bandwidth in bins


    enbw = HANN_ENBW

END FUNCTION get_enbw_WindFarm


! ===============================================================================
! DEFAULT OUTPUT FILE NAME
! ===============================================================================
SUBROUTINE build_default_name(self, name_tag, name_suffix)
    ! -------------------------------------------------------------------------
    ! Builds the labels of the default output file name. The head counts the
    ! turbines of each type in order of first appearance, with the same
    ! aliases as the Python WindFarm (e.g. Farm_2_DTU10MN_1_SAITEC2FL), or is
    ! Farm_<farm_name> when a label is given in WindFarm.nml. The tail holds
    ! the wind speed: a single value when every turbine shares it (_ws11.4),
    ! otherwise its range (_ws8.0-11.4). The solver inserts its method tag in
    ! between, e.g. Farm_2_DTU10MN_ANM24_ws8.0-11.4.hdf5.
    ! -------------------------------------------------------------------------

    CLASS(WindFarm_t), INTENT(IN)  :: self           ! [-] Farm instance
    CHARACTER(len=*) , INTENT(OUT) :: name_tag       ! [-] File name head
    CHARACTER(len=*) , INTENT(OUT) :: name_suffix    ! [-] File name tail

    ! Local variables
    INTEGER(I32)      :: i, j, n_same
    LOGICAL           :: seen
    CHARACTER(len=16) :: n_str, ws_min_str, ws_max_str


    ! --- 1. Head: user label, or number of turbines of each type --- !
    if (to_upper(trim(self % farm_name)) /= 'NONE') then
        name_tag = 'Farm_' // trim(self % farm_name)
    else
        name_tag = 'Farm'
        do i = 1, self % Nturbines
            associate(type_i => self % turbines(i) % wt % case_type)
                ! Every type is counted once, at its first appearance
                seen = .false.
                do j = 1, i - 1
                    if (self % turbines(j) % wt % case_type == type_i) seen = .true.
                end do

                if (.not. seen) then
                    n_same = 0_I32
                    do j = i, self % Nturbines
                        if (self % turbines(j) % wt % case_type == type_i) n_same = n_same + 1_I32
                    end do

                    write(n_str, '(I0)') n_same
                    name_tag = trim(name_tag) // '_' // trim(n_str) // '_' // type_alias(type_i)
                end if
            end associate
        end do
    end if

    ! --- 2. Tail: common wind speed, or its range across the farm --- !
    write(ws_min_str, '(F6.1)') minval([(self % turbines(i) % wt % WindSpeed, i = 1, self % Nturbines)])
    write(ws_max_str, '(F6.1)') maxval([(self % turbines(i) % wt % WindSpeed, i = 1, self % Nturbines)])

    if (ws_min_str == ws_max_str) then
        name_suffix = '_ws' // trim(adjustl(ws_min_str))
    else
        name_suffix = '_ws' // trim(adjustl(ws_min_str)) // '-' // trim(adjustl(ws_max_str))
    end if

END SUBROUTINE build_default_name


FUNCTION type_alias(case_type) RESULT(alias)
    ! -------------------------------------------------------------------------
    ! Short turbine-type label used in the default farm file name, identical
    ! to the aliases of the Python WindFarm. Unknown types keep their name.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN)  :: case_type    ! [-] Turbine model tag
    CHARACTER(len=:), ALLOCATABLE :: alias        ! [-] Short label for the file name


    SELECT CASE (trim(case_type))
    CASE ('DTU10MWMonopile')
        alias = 'DTU10MN'
    CASE ('DTU10MWFloating')
        alias = 'DTU10FL'
    CASE ('SAITEC2MWFloating')
        alias = 'SAITEC2FL'
    CASE ('MYSE3MWMonopile')
        alias = 'MYSE3MN'
    CASE DEFAULT
        alias = trim(case_type)
    END SELECT

END FUNCTION type_alias

END MODULE WindFarm
