# Architecture of the off-coustics Fortran solver

This document describes how the solver in `FORTRAN/src` and `FORTRAN/app` is organised after the refactor of 10/2026: the data flow, the modules and their layers, the contracts between them, the output interface, and the decisions taken for future work. The code is the reference; this file tells where things are and why they are there.

All Fortran code follows `Programming_style_guidelines.txt`.

---

## 1. Pipeline

A run of `off_coustics_fortran` (`app/main.f90`) goes through these stages:

```
 input_files/  (Solver.nml, Environment.nml, Emitters/, ComputationSetups/, AcousticMethods/)
        |
        |  read_input_files                                   input/
        v
    Config_t ---------------------------------------------+
        |                                                 |
        |  new_turbine + init/read_input/compute_force,   |
        |  WindFarm_t % init/compute_forces,              |  emitters/
        |  or FrequencySource_t % init                    |
        v                                                 |
    Emitter_t  (WindTurbine_t, WindFarm_t, ...)           |
        |                                                 |
        |  emitter % to_source                            |
        v                                                 |
    AcousticSource_t  (xyz, Force, Freqs, order, geometry)|  types/
        |                                                 |
        |  new_solver + solver % init(sources, config) <--+  solvers/
        v
    AcousticSolver_t + method (ImageMethod_t, AnalyticalNormalModes_t, ...)
        |
        |  save_parameters(emitter)   -> HDF5 parameters + .dat header
        |  run()                      -> run_* steps -> compute_pressure_masked -> compute_pressure
        |  update_parameters()        -> .dat post-processing block
        v
    <out>.hdf5 + <out>.dat
        |
        v
    plotting/  (Python package)
```

In `main.f90`:

1. `get_input_path(files_path)` takes the directory of the namelists from the command line (`./off_coustics ./path_to_input_files`; `./input_files/` without arguments) and `read_input_files(files_path, config)` reads every namelist of it into `Config_t` (the files are described in §1.1). Paths written inside the namelists (e.g. `OP_output_dir`) stay relative to the working directory, not to that directory.
2. Depending on `simulation_type`, a single turbine (`new_turbine`, then `init`, `read_input`, `compute_force`), a wind farm (`WindFarm_t % init`, `compute_forces`) or a set of frequency-domain point sources (`FrequencySource_t % init`) is built.
3. `to_source` translates the emitter into an `AcousticSource_t`.
4. `new_solver` allocates the solver named by `method`; `init` moves the sources into it and frees `config`.
5. At the end of its branch the emitter is moved (`MOVE_ALLOC`) into a `CLASS(Emitter_t)` variable, which is handed to `save_parameters` and deallocated before `run()`.
6. `run()` executes every enabled post-processing step; `update_parameters()` closes the `.dat` summary.

The solvers never know what produced the sources. Which steps make sense for each emitter is decided by the parser through the `do_run` flags (sphere and decay are switched off for a wind farm).

### 1.1 Input files

`read_input_files` reads `Solver.nml` and `Environment.nml` from the input directory and then three files that it picks by name, so that every emitter and every method has its own file:

| File | Chosen by | Files |
|---|---|---|
| `Emitters/<simulation_type>.nml` | `simulation_type` of `Solver.nml` | `WindTurbine.nml`, `WindFarm.nml`, `FrequencySource.nml` |
| `ComputationSetups/<simulation_type>_CS.nml` | `simulation_type` | `WindTurbine_CS.nml`, `WindFarm_CS.nml`, `FrequencySource_CS.nml` |
| `AcousticMethods/<method>.nml` | `method` of `Solver.nml` | `ImageMethod.nml`, `AnalyticalNormalModes.nml` |

The values of `simulation_type` are `'WindTurbine'`, `'WindFarm'` and `'FrequencySource'` (the name of the emitter module, without `_t`), and those of `method` are `'ImageMethod'` and `'AnalyticalNormalModes'`. They are compared exactly (capitals included) and they are also the file names, so a new emitter or method needs no change in the path logic, only its `CASE` in `read_input_files` and in `main.f90`, and its files. The folder names are fixed in `read_input_files`; the directory itself is the command-line argument.

---

## 2. Directory tree

```
src/
├── Kinds.f90                    Kinds
├── types/
│   ├── SourceTypes.f90          SourceTypes
│   └── ConfigTypes.f90          ConfigTypes
├── utils/
│   ├── FFTW3.f90                FFTW3
│   ├── StringUtils.f90          StringUtils
│   ├── FileUtils.f90            FileUtils
│   ├── HDF5Utils.f90            HDF5Utils
│   ├── CSVReader.f90            CSVReader
│   ├── Sorting.f90              Sorting
│   ├── Timing.f90               Timing
│   ├── SignalProcessing.f90     SignalProcessing
│   └── SpecialFunctions.f90     SpecialFunctions
├── input/
│   ├── Args_parser.f90          Args_parser
│   ├── InputEmitters.f90        InputEmitters
│   ├── InputSolver.f90          InputSolver
│   ├── InputPostProcessing.f90  InputPostProcessing
│   └── InputChecks.f90          InputChecks
├── emitters/
│   ├── Emitter.f90              Emitter
│   ├── WindTurbine.f90          WindTurbine
│   ├── TurbineFactory.f90       TurbineFactory
│   ├── WindFarm.f90             WindFarm
│   ├── FrequencySource.f90      FrequencySource
│   ├── OpenFASTReader.f90       OpenFASTReader
│   ├── StructureModel.f90       StructureModel
│   ├── DrivetrainExcitation.f90 DrivetrainExcitation
│   └── TurbineTypes/
│       ├── DTU10MWMonopile.f90  DTU10MWMonopile
│       ├── DTU10MWFloating.f90  DTU10MWFloating
│       └── SAITEC2MWFloating.f90 SAITEC2MWFloating
└── solvers/
    ├── AcousticSolver.f90       AcousticSolver (parent module)
    ├── AcousticSolver_steps.f90 SUBMODULE (AcousticSolver) AcousticSolver_steps
    ├── AcousticSolver_output.f90 SUBMODULE (AcousticSolver) AcousticSolver_output
    ├── SolverFactory.f90        SolverFactory
    └── methods/
        ├── ImageMethod.f90      ImageMethod
        └── AnalyticalNormalModes.f90 AnalyticalNormalModes
app/main.f90                     PROGRAM main
```

What every module contains, and what must not go there:

| Module | Contains | Must not contain |
|---|---|---|
| `Kinds` | Kind parameters (`I32`, `WP`, ...), `PI`, `I1`, the parser sentinels `DEFAULT_R`, `DEFAULT_I` and the parser buffer sizes `MAX_PARSE_OBS`, `MAX_PARSE_SRC`, `MAX_PARSE_FREQ` | Types or procedures |
| `SourceTypes` | `AcousticSource_t`, `GeometricData_t`, the multipole labels `SOURCE_MONOPOLE`, `SOURCE_DIPOLE` | Procedures, configuration |
| `ConfigTypes` | Every `Config_*_t` filled in by the parser, including `Config_Method_t` and its extensions | Procedures, run-time state |
| `FFTW3` | C bindings of FFTW3 | Anything else |
| `StringUtils` | `to_upper`, `split_tokens`, `is_data_line` | File or HDF5 I/O |
| `FileUtils` | `read_full_line`, `swap_extension`, `ensure_parent_dir`, `move_text_file`, `delete_file` | HDF5, format-specific parsers |
| `HDF5Utils` | `save_variable` (dispatcher over the generic `save_to_hdf5`), `read_hdf5_real_1d`. The only module that does `USE HDF5` | Dataset names or knowledge of what is saved |
| `CSVReader` | `read_curve` (wind-speed lookup tables), `parse_csv_row` | Turbine-specific files |
| `Sorting` | `sort_real_c` (C `qsort`) | Anything else |
| `Timing` | `format_elapsed` | Anything else |
| `SignalProcessing` | `compute_rfft` (mean removal, Hann window, FFTW), `HANN_ENBW` | Turbine or solver knowledge |
| `SpecialFunctions` | `compute_hankel_complex`, `bessel_k0`, `bessel_k1` | Method-specific kernels |
| `Args_parser` | `get_input_path` (command-line directory of the namelists), `read_input_files` (finds the files of §1.1 and calls one reader per file), the entry of the parser | Readers of individual namelists |
| `InputEmitters` | `read_turbine_config`, `read_farm_config`, `read_frequency_source_config` | Solver or step settings |
| `InputSolver` | `read_solver_config`, `read_environment_config`, one `read_X_config` per method | Emitter or step settings |
| `InputPostProcessing` | `read_computation_config`, `read_do_config`, one reader per `run_*` step | Emitter or solver settings |
| `InputChecks` | `check_time_window`, `check_fits`, `check_all_or_none`, `check_count` | Namelist reading |
| `Emitter` | Abstract `Emitter_t` and its default `get_enbw` | Turbine or farm specifics |
| `WindTurbine` | Abstract `WindTurbine_t`: shared parameters, OpenFAST reading, `filter_frequencies`, the `Emitter_t` procedures of a turbine | Data of a concrete turbine |
| `TurbineTypes/*` | One concrete turbine each: its geometry, `compute_force`, `compute_mass`, `excitation_spectrum`, its saved parameters and impedance correction | Generic routines (they go to `StructureModel`, `DrivetrainExcitation`, `OpenFASTReader`) |
| `TurbineFactory` | `new_turbine`: one `CASE` per turbine type | Anything else |
| `WindFarm` | `WindFarm_t`: a container of turbines merged into one `AcousticSource_t` | Turbine-type specifics |
| `FrequencySource` | `FrequencySource_t`: point dipoles placed by hand, with forces given in the frequency domain | Turbine or farm specifics, time-domain forces |
| `OpenFASTReader` | SubDyn summary and `.out`/`.outb` channel readers | Physics |
| `StructureModel` | Geometry and hydrodynamic corrections of the members: `remove_duplicate_nodes`, `divide_span`, `integral_elipse_z`, `alpha_hankel` | Data of a concrete turbine |
| `DrivetrainExcitation` | Drivetrain spectra (`drivetrain10MW_excitation_spectrum`, `drivetrain2MW_excitation_spectrum` and its helpers), `generate_timeseries_banded_sines`, `filter_non_useful_freqs` | Turbine state |
| `AcousticSolver` | Abstract `AcousticSolver_t`, its deferred interfaces, the interfaces of the submodule procedures, `check_observers_distances`, `compute_pressure_masked`, `pressure_norm`, the water-column queries | Method kernels, emitter types |
| `AcousticSolver_steps` | `run()`, every `run_*` step, `dummy_run_sphere` | Reads of `depth_base` (use the queries), method kernels |
| `AcousticSolver_output` | `resolve_save_path`, `check_update`, `update_parameters`, `save_shared_parameters` and their private helpers | Pressure computations |
| `SolverFactory` | `new_solver`: one `CASE` per method | Anything else |
| `methods/*` | One acoustic method each, extending `AcousticSolver_t` | Emitter types (they only see `Emitter_t`) |

---

## 3. Dependency layers

```
Kinds  ->  types  ->  utils  ->  input / emitters  ->  solvers  ->  app
```

**Rule:** a module never `USE`s a module of a higher layer. `input` and `emitters` are on the same level and never use each other. Inside a layer the order is the one of the Makefile.

| Layer | Modules |
|---|---|
| Kinds | `Kinds` |
| types | `SourceTypes`, `ConfigTypes` |
| utils | `FFTW3`, `StringUtils`, `FileUtils`, `HDF5Utils`, `CSVReader`, `Sorting`, `Timing`, `SignalProcessing`, `SpecialFunctions` |
| input | `InputChecks`, `InputEmitters`, `InputSolver`, `InputPostProcessing`, `Args_parser` |
| emitters | `OpenFASTReader`, `StructureModel`, `DrivetrainExcitation`, `Emitter`, `WindTurbine`, `DTU10MWMonopile`, `DTU10MWFloating`, `SAITEC2MWFloating`, `TurbineFactory`, `WindFarm`, `FrequencySource` |
| solvers | `AcousticSolver`, `AcousticSolver_output`, `AcousticSolver_steps`, `ImageMethod`, `AnalyticalNormalModes`, `SolverFactory` |
| app | `main` |

Other rules:
- Every module is `PRIVATE` by default with an explicit `PUBLIC ::` list of what other modules use. There are no umbrella or compatibility modules.
- `USE Module, ONLY: ...` is preferred; `-Wall` reports unused imports.
- A submodule depends on its parent's object in the Makefile; gfortran writes its `.smod` into `build_make/`. fpm discovers `src/**` by itself.
- `AcousticSolver` cannot `USE` a method (that would be a cycle); this is why `run_sphere` is a no-op in the parent and overridden in `ImageMethod`.

---

## 4. Contracts

### 4.1 `Emitter_t` (`emitters/Emitter.f90`)

| Procedure | Kind | Meaning |
|---|---|---|
| `to_source(self, config, sources)` | deferred | Fills the `AcousticSource_t` handed to the solver: `xyz`, `Force`, `Freqs`, `order`, `geometry` |
| `save_params(self, save_path)` | deferred | Writes the emitter datasets into the HDF5 output |
| `write_txt(self, unit)` | deferred | Writes the emitter block of the `.dat` summary |
| `get_enbw(self)` | default `-1.0` | Equivalent noise bandwidth, in bins, of the window applied to the forces; `-1` means "forces not obtained from a windowed FFT". `WindTurbine_t` and `WindFarm_t` return `HANN_ENBW` (1.5); `FrequencySource_t` keeps `-1` |

`WindTurbine_t` implements them with `turbine_to_source`, `save_params_turbine` (`save_turbine_params` + `save_specific_params`), `write_txt_turbine` (`write_turbine_txt`) and `get_enbw_turbine`. `WindFarm_t` implements them with `farm_to_source_WindFarm`, `save_params_WindFarm`, `write_txt_WindFarm` and `get_enbw_WindFarm`. `FrequencySource_t` implements them with `to_source_FrequencySource`, `save_params_FrequencySource` and `write_txt_FrequencySource`; it also has a non-deferred `init(self, config)` that copies the sources of `config % freq_source`. Its `to_source` copies the arrays instead of moving them, because `save_params` and `write_txt` run after the solver has taken the `AcousticSource_t`.

`WindTurbine_t` is itself a template for the turbine types:
- **Shared:** `init`, `read_input`, `translate_align_with_WindDir`, `filter_frequencies` (rpm from `read_curve`, the type's `excitation_spectrum`, synthetic signal, frequency mask).
- **Deferred:** `compute_force`, `compute_mass`, `excitation_spectrum(self, rpm, freqs_amp, keys)`, `save_specific_params`, `write_specific_params_txt`, `get_impedance_corrected_force`.

### 4.2 `AcousticSolver_t` (`solvers/AcousticSolver.f90`)

Deferred procedures, implemented by every method:

| Procedure | Contract |
|---|---|
| `init(self, sources, config)` | Moves `xyz`, `Force`, `Freqs` and `order` out of `sources` (`MOVE_ALLOC`), stops if a source is not a dipole, copies `env`, `run_inputs` and the solver settings (`save_path`, `extra_label`, `block_size`, `debug`, `in_cluster`, `update_existing`, `mask_distance`), reads its own `Config_X_t`, calls `build_save_path`, prints `get_name()`, and frees `sources` and `config` |
| `get_name(self)` | Returns the human-readable name and stores it in `solver_name` (written as the `Method` dataset) |
| `compute_pressure(self, observers, block_size, total_pressure)` | Pressure `(Nfreqs, Nobs)` at `observers (Nobs, 3)` for every source frequency. Called only through `compute_pressure_masked` |
| `save_parameters(self, emitter)` | In this order: `check_update(emitter)`, its own datasets, `save_shared_parameters(emitter)`, `emitter % save_params`, `write_txt_summary(emitter)` |
| `build_save_path(self)` | Calls `resolve_save_path(method_tag)` with its tag, e.g. `_IM30`, `_ANM24` |
| `write_txt_summary(self, emitter)` | Writes the `.dat` header: title, `emitter % write_txt`, environment and solver blocks |

Water-column queries (non-deferred, overridable, implemented in the parent):

| Binding | Default | Use |
|---|---|---|
| `z_seabed(x, y)` | `-depth_base` | Bottom of the water column at `(x, y)` |
| `in_water(observers)` | all `.true.` | Mask of the observers inside the water column |

The `run_*` steps never read `depth_base`: at their reference point they ask `z_top = 0` and `z_bot = z_seabed(xr, yr)`, and use `z_top = 0`, `z_bot` and `0.5 * (z_top + z_bot)` for the surface, the seabed and mid-depth. With the defaults every expression gives the same bits as the former literals.

| Step | Reference point | What uses the queries |
|---|---|---|
| `run_spectrums` | each default observer (10/250/500, 0) | default z (mid-depth) |
| `run_polar` | `center_polar` | default `z_polar` |
| `run_cylinder` | `center_cyl` | z levels from `z_bot` to `z_top` |
| `run_decay` | `geometry % origin` | default `z_decay` |
| `run_line` | (100, 0) | default `p1_line(3)` (surface) and `p2_line(3)` (seabed) |
| `run_sliceXY` | `geometry % center` | default `z_xy`, check `z_xy > z_top` |
| `run_sliceXZ` | middle of `xlim_xz`, at `y_xz` | default `zlim_xz` |
| `run_sliceVertical` | `center_sv` | default and clipping of `zlim_sv` |

Other shared procedures:
- `compute_pressure_masked(label, observers, block_size, total_pressure)`: the wrapper every step uses. Observers closer than `mask_distance` to a source node (`check_observers_distances`) or outside the water column (`in_water`) are not computed and get NaN; one WARNING line per criterion reports how many. Without masked observers it is a plain `compute_pressure` call.
- `save_shared_parameters(emitter)`: writes `Method` (`solver_name`) and, if `emitter % get_enbw() > 0`, `ENBW`. Neither goes to the `.dat`.
- `run()`: executes the enabled steps in a fixed order (spectrums, polar, decay, line, sphere, sliceXY, sliceXZ, sliceVertical, cylinder).
- `run_sphere`: a no-op in the parent (`dummy_run_sphere`, which switches `do_sphere` off); `ImageMethod` overrides it with the free-field sphere.
- `resolve_save_path`, `check_update`, `update_parameters`: see §6.3.

### 4.3 `AcousticSource_t` and `GeometricData_t` (`types/SourceTypes.f90`)

| Field | Shape | Meaning |
|---|---|---|
| `xyz` | `(Nsources, 3)` [m] | Source coordinates |
| `Force` | `(Nfreqs, Nsources, 3)` [N] | Complex source forces (radiation-impedance corrected for turbines) |
| `Freqs` | `(Nfreqs)` [Hz] | Frequency vector |
| `order` | `(Nsources)` [-] | Multipole order label (`SOURCE_MONOPOLE = 0`, `SOURCE_DIPOLE = 1`); saved as `Source_order` by `FrequencySource_t` only |
| `geometry % center` | `(2)` [m] | Reference centre: turbine barycentre, or mean of the farm or of the sources |
| `geometry % origin` | `(2)` [m] | Reference origin: turbine axis, or mean of the farm or of the sources |
| `geometry % direction` | [deg] | Reference direction: wind direction, circular mean of the farm, or 0 for point sources |
| `geometry % out_dir`, `name_tag`, `name_suffix` | | Default output name `out_dir // name_tag // <method tag> // name_suffix // ".hdf5"` |

What each consumer takes from them:

| Consumer | Uses |
|---|---|
| `compute_pressure` (methods) | `xyz`, `Force`, `Freqs` |
| `check_observers_distances` | `xyz` |
| `check_update` | `Freqs` |
| `resolve_save_path` | `out_dir`, `name_tag`, `name_suffix` |
| `run_spectrums`, `run_line` | nothing (fixed default points) |
| `run_polar`, `run_cylinder` | `center` (default centre) |
| `run_decay` | `origin` (start), `direction` (line direction) |
| `run_sliceXY` | `center` (default limits, queries) |
| `run_sliceXZ` | `origin(2)` (default `y_xz`), `center(1)` (default `xlim_xz`) |
| `run_sliceVertical` | `origin` (default centre), `xyz` (automatic width) |
| `run_sphere` (IM) | `xyz` (default centre and radius) |
| every step | `Freqs` (saved with each step) |

---

## 5. How to add...

### 5.1 A turbine type
1. `src/emitters/TurbineTypes/X.f90`: `X_t EXTENDS(WindTurbine_t)` with its geometry and the deferred procedures, `excitation_spectrum` included. Its `init` sets `path_rpm` (and any data path), `is_monopile`, and calls the base `init`.
2. If its drivetrain is new, add the spectrum to `emitters/DrivetrainExcitation.f90` and call it from `excitation_spectrum`.
3. One `CASE` in `new_turbine` (`emitters/TurbineFactory.f90`).
4. Optionally its short label in `type_alias` (`emitters/WindFarm.f90`) for the default farm file name.
5. Makefile: the object in `OBJS` and its rule, before `TurbineFactory`.

### 5.2 An emitter (example: `FrequencySource`)
1. `src/emitters/X.f90`: `X_t EXTENDS(Emitter_t)` with:
   - `init(self, config)`, which copies its data from `Config_t`. It is not part of the contract because every emitter is built differently.
   - `to_source`, filling `xyz`, `Force`, `Freqs`, `order` and `geometry`.
   - `save_params`, which writes the marker `Emitter = '<X>'` and everything the plotting needs (§6.1).
   - `write_txt`, which must list every input value: `check_update` compares this block, so an input missing from it would reuse an old output after being changed.
   - `get_enbw`, if its forces come from `compute_rfft`.
2. Its configuration type in `types/ConfigTypes.f90` and its field in `Config_t`; its reader in `input/InputEmitters.f90`; its `CASE` in `read_input_files`; its files `Emitters/<X>.nml` and `ComputationSetups/<X>_CS.nml`, named exactly as its `simulation_type` (§1.1).
3. Its construction branch in `app/main.f90`: a typed variable (the abstract `active_emitter` can neither be allocated nor give access to `init`), `init`, `to_source`, and `MOVE_ALLOC` into `active_emitter` at the end of the branch.
4. Makefile: the object in `OBJS` and its rule, after `Emitter`.
5. Plotting: a reader and a wrapper in `plotting/` (see its README); the `Emitter` marker is how the scripts recognise the file.

The planned emitters are in §7.1.

### 5.3 An acoustic method
1. `src/solvers/methods/X.f90`: `X_t EXTENDS(AcousticSolver_t)` with the deferred procedures (§4.2); `save_parameters` calls `save_shared_parameters`.
2. `Config_X_t EXTENDS(Config_Method_t)` in `types/ConfigTypes.f90`.
3. `read_X_config` in `input/InputSolver.f90`, and its `CASE` in `read_input_files` (`input/Args_parser.f90`), with its own `.nml` in `input_files/AcousticMethods/<method>.nml`.
4. Its `CASE` in `new_solver` (`solvers/SolverFactory.f90`).
5. Makefile: the object before `SolverFactory`.
6. Nothing else: the `Method` dataset is written automatically.
7. Always surface is supposed to be at z = 0

If the method has its own boundaries, it overrides `z_seabed` and `in_water`, not the `run_*` steps.

### 5.4 A post-processing step
1. `run_X` in `solvers/AcousticSolver_steps.f90`, its interface body and type-bound binding in `solvers/AcousticSolver.f90`, and its call in `run()`.
2. `Config_X_t` and a `do_X` flag in `Config_DoRun_t` / `Config_Computation_t` (`types/ConfigTypes.f90`).
3. `read_X_config` and the flag in `read_do_config` (`input/InputPostProcessing.f90`), plus its namelist group in every `input_files/ComputationSetups/*_CS.nml`.
4. Its block in `update_parameters` (`solvers/AcousticSolver_output.f90`), with `write_previous_block` for updates.
5. Its plots in `plotting/`.

---

## 6. Output interface

> **Warning:** HDF5 dataset names, the `.dat` format and the namelist variable names are an interface. `plotting/` reads them, and `check_update` compares the `.dat` header. They are never renamed; new ones are added.

### 6.1 HDF5 datasets

Complex arrays are stored as `(2, N...)` in Fortran order, i.e. `(..., N, 2)` (real, imaginary) when read with h5py. Scalars are 1-element datasets; strings are fixed-length byte strings. Every step also writes `Freqs`.

| Writer | Datasets |
|---|---|
| `ImageMethod % save_parameters` | `N_images`, `c_wat`, `rho_wat`, `Lower_HBC`, `p_ref` |
| `AnalyticalNormalModes % save_parameters` | `Nmodes`, `H`, `c_wat`, `rho_wat`, `Lower_HBC`, `p_ref` |
| `save_shared_parameters` | `Method` (`"Images Method"`, `"Analytical Normal Modes"`), `ENBW` (1.5, only for forces from a windowed FFT) |
| turbine (`save_turbine_params`) | `case_type`, `WindSpeed`, `WindDir`, `Depth`, `AxisPos`, `BariPos` |
| turbine (`save_specific_params`) | `Case_type`, `Structure_nodes` plus the type's geometry (`D`, `rho_mat`, `wet_area`; `col_D`, `XsecA`, `XsecB`, `t`, `col_members`, `pon_members`; `Diam_y`, `Diam_z`, `float_members`, `joint_members`) |
| farm (`save_params_WindFarm`) | `Nturbines`, `Source_offset`, and every turbine dataset with the suffix `_i` (e.g. `WindSpeed_2`, `Structure_nodes_2`) |
| `FrequencySource` (`save_params_FrequencySource`) | `Emitter` (`"FrequencySource"`), `Nsources`, `NFreqs`, `Source_xyz` `(Nsources, 3)`, `Source_order` `(Nsources)`, `Source_Force` complex `(Nfreqs, Nsources, 3)` |
| `run_spectrums` | `P_spectrums`, `Obs_spectrums` |
| `run_polar` | `P_polar`, `R_polar`, `Z_polar`, `Theta_deg_polar`, `Obs_polar`, `Center_polar` |
| `run_cylinder` | `P_cylinder`, `R_cylinder`, `Z_cylinder`, `Theta_deg_cylinder`, `Obs_cylinder`, `Center_cylinder`, `dA_cylinder` |
| `run_decay` | `P_decay`, `Z_decay`, `Distance_decay`, `Distances_decay`, `Obs_decay`, `Logspace_decay` |
| `run_line` | `P_line`, `P1_line`, `P2_line`, `Obs_line`, `Distance_line`, `Distances_line`, `Logspace_line` |
| `run_sliceXY` | `P_slicexy`, `Obs_slicexy`, `X_slicexy`, `Y_slicexy`, `Z_slicexy`, `Center_slicexy`, `Nx_slicexy`, `Ny_slicexy` |
| `run_sliceXZ` | `P_slicexz`, `Obs_slicexz`, `X_slicexz`, `Z_slicexz`, `Y_slicexz`, `Nx_slicexz`, `Nz_slicexz` |
| `run_sliceVertical` | `P_sliceVertical`, `Coords_sliceV`, `U_sliceV`, `Z_sliceV`, `Azimuth_sliceV`, `Center_sliceV` |
| `run_sphere` (IM) | `P_sphere`, `Obs_sphere`, `Theta_sphere`, `Z_sphere`, `R_sphere`, `Center_sphere`, `dA_sphere` |

Masked observers (closer than `mask_distance` to a source, or outside the water column) have NaN pressure; the grids keep their shape.

`Emitter` identifies the emitter that wrote a file and is written only by the emitters added after the refactor (`FrequencySource`); the plotting recognises the older turbine and farm files by their datasets (`Nturbines`, `case_type`, ...). The water depth is not an emitter dataset: both solvers write `Lower_HBC` (`-depth_base`) and the plotting derives the depth from it.

### 6.2 The `.dat` summary

Written next to the HDF5 file (same name, `.dat`), with blocks between separators of 70 `=` or `-`:

1. Title `OFF_COUSTICS SIMULATION SUMMARY` and `HDF5 file: <path>`.
2. Emitter block (`emitter % write_txt`): `TURBINE` and `TURBINE - SPECIFIC PARAMETERS`, or `WIND FARM` (layout table) followed by `>>> Turbine i` and the blocks of every turbine, or `FREQUENCY-DOMAIN SOURCES` (positions and one line per frequency and source with the force, i.e. the whole input, because `check_update` compares this block).
3. `ENVIRONMENT` and `SOLVER: <name>` blocks (`write_txt_summary`).
4. `POST-PROCESSING STEPS EXECUTED` (`update_parameters`): one block per step (`Spectrums:`, `Polar:`, ...) with the parameters actually used, defaults resolved.
5. A final `=` separator.

Items 1-3 are the header compared by `check_update`.

### 6.3 `update_existing`

With `update_existing = .true.` in `Solver.nml`, an existing output of the same case is updated instead of rebuilt:
- `resolve_save_path` deletes any stale `.dat.prev`. If the HDF5 file exists, it is kept and its `.dat` is moved to `.dat.prev` (with `.false.` it is deleted).
- `check_update`, the first call of `save_parameters`, keeps it only if it holds the same case:
  - the same `Freqs` (same size, relative tolerance 1e-9);
  - the same `.dat` header (everything before the post-processing block, without the `Block size` and `In cluster` lines).

  Otherwise the file and its `.dat.prev` are deleted with a WARNING and the output is rebuilt.
- The steps run now overwrite their datasets; the other datasets stay.
- `update_parameters` writes the blocks of the steps run now and copies the blocks of the other steps from `.dat.prev`, tagged `(previous run)`, then deletes `.dat.prev`.

`Method` and `ENBW` are not in the `.dat`, so they do not take part in the comparison.

---

## 7. Planned extensions

Decided on 2026-10-02. Section 7.1 is partly implemented.

### 7.1 Manual sources
**Implemented (2026-10-06): `FrequencySource`.** `simulation_type = 'FrequencySource'`, `emitters/FrequencySource.f90` and `Emitters/FrequencySource.nml`. A few point dipoles with complex forces per frequency, with no window (`get_enbw` keeps `-1`, so no `ENBW` is written) and no radiation-impedance correction. `mask_distance` is set by the user in `Solver.nml`. The plotting reads `ENBW` (Hann correction only if present) and `Method`.

**Planned:**
- `TimeSource`: the same sources with forces given in the time domain (real), read from a file or computed from a function of time. They go through `compute_rfft`, so `get_enbw` returns `HANN_ENBW`.
- A parent `ManualSources_t EXTENDS(Emitter_t)` (`emitters/ManualSources.f90`) with what both manual emitters share.
- An input HDF5 file for many sources (`xyz`, `order`, `Freqs`, `Force`).
- One multipole-order label per source (already in place). Each method implements one kernel per order, and drops its dipole-only check in `init` when it supports monopoles.

### 7.2 Acoustic methods
- `Config_Environment_t` will grow optional arrays: c(z), ρ(z), attenuation, bathymetry H(x, y), sediment. It keeps the base scalars, so the simple methods keep working. Each method checks in `init` that the data it needs is present.
- The boundaries belong to each method: constants in IM and ANM, arrays in complex methods. A method publishes them by overriding `z_seabed` and `in_water`.
- Overriding a whole `run_*` is the exception (precedent: IM's `run_sphere`).

### 7.3 Sphere
The sphere stays IM-only: a no-op in the parent plus the IM override, because `AcousticSolver` cannot `USE ImageMethod` (cycle). A shared free-field kernel would make it generic for every method; that changes the outputs and is a future decision.

---

## 8. Verification protocol

Any change that is meant not to change the results follows `Programming_style_guidelines.txt` §8, plus the harness of `review/refactor_tools/` (git-ignored, kept locally):

1. Keep a copy of the original sources and build the reference binary from them. Run everything in an isolated copy, never in `turbine_acoustic_data/`, `farm_acoustic_data/` or `SAITEC_acoustic_data/`.
2. Builds: `make clean && make -j8`, `fpm build --flag "-O3 -fopenmp"`, and `make FFLAGS="-O3 -fopenmp -Wall"`, which must give no warning that the reference does not give. Since 10/2026 the only accepted `-Wall` warnings are the `-Wunused-dummy-argument` of the default procedures that ignore arguments by design: `default_z_seabed` (`x`, `y`), `default_in_water` (`self`), `excitation_spectrum` of the two DTU types (`self`), `get_enbw` of `Emitter_t`, `WindTurbine_t` and `WindFarm_t` (`self`), and `to_source` of `FrequencySource_t` (`config`, which the `Emitter_t` interface requires).
3. Cases, with `OMP_NUM_THREADS=1` (multi-threaded runs are not bit-reproducible):
   ```bash
   python refactor_cases.py setup   <reference_src_dir> <cases_root>
   python refactor_cases.py run     <binary> <cases_root> <tag>
   python refactor_cases.py compare <cases_root>/<tagA> <cases_root>/<tagB> [extra1,extra2,...]
   ```
   The cases are M1 (DTU 10 MW monopile, IM, all steps), M2 (the same with ANM), M3 (`update_existing`: polar + line, then sphere only, then `N_images = 20`), and L1, L2, L3 (DTU floating, SAITEC, 2-monopile farm, spectrums only). They take about 75 s.
4. Acceptance: `compare` prints `ALL IDENTICAL`. Every HDF5 dataset of the reference exists with the same dtype, shape and bytes; there is no extra dataset except the allowed ones; the `.dat` files are identical; stdout is identical except the `Total elapsed time` line. A difference is never accepted: its cause is found and fixed.
