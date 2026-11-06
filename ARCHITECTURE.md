# Architecture of the off-coustics Fortran solver

This document describes how the solver in `FORTRAN/src` and `FORTRAN/app` is organised after the refactor of 10/2026 and the manual sources added since (`FrequencySource`, `TimeSource`): the data flow, the modules and their layers, the contracts between them, the output interface, and the decisions taken for future work. The code is the reference; this file tells where things are and why they are there. It describes the code with `frequency_source_inputs.patch` applied (08/10/2026), the monopoles of the manual sources (13/10/2026), and the grouping of the sources into vertical axes used by the normal-mode kernels (06/11/2026).

All Fortran code follows `Programming_style_guidelines.txt`.

---

## 1. Pipeline

A run of `off_coustics` (`app/main.f90`; the executable of `make`, while `fpm` names it after the package, `Off_coustics`) goes through these stages:

```
 input_files/  (Solver.nml, Environment.nml, Emitters/, ComputationSetups/, AcousticMethods/)
        |
        |  get_input_path + read_input_files                   input/
        v
    Config_t ---------------------------------------------+
        |                                                 |
        |  new_turbine + init/read_input/compute_force,   |
        |  WindFarm_t % init/compute_forces,              |  emitters/
        |  FrequencySource_t % init or TimeSource_t % init|
        v                                                 |
    Emitter_t  (WindTurbine_t, WindFarm_t, ManualSources_t)
        |                                                 |
        |  emitter % to_source                            |
        v                                                 |
    AcousticSource_t  (xyz, Force, Freqs, order, geometry)|  types/
        |                                                 |
        |  new_solver + solver % init(sources, config) <--+  solvers/
        v
    AcousticSolver_t + method (ImageMethod_t, AnalyticalNormalModes_t)
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

1. `get_input_path(files_path)` takes the directory of the namelists from the command line (`./off_coustics ./path_to_input_files`; `./input_files/` without arguments) and `read_input_files(files_path, config)` reads every namelist of it into `Config_t` (the files are described in §1.1). Paths written inside the namelists (e.g. `OP_output_dir`, `freq_file`, `source_file`) stay relative to the working directory, not to that directory.
2. Depending on `simulation_type`, one emitter is built:
   - `'WindTurbine'`: `new_turbine`, then `init`, `read_input`, `compute_force`;
   - `'WindFarm'`: `WindFarm_t % init`, `compute_forces`;
   - `'FrequencySource'`: `FrequencySource_t % init` (point monopoles and dipoles with forces (volume velocities) given per frequency);
   - `'TimeSource'`: `TimeSource_t % init` (point monopoles and dipoles with forces (volume velocities) given in time, taken to the frequency domain with `compute_rfft`).
   Any other value stops the run (`CASE DEFAULT` of `read_input_files`, and again in `main.f90`).
3. `to_source` translates the emitter into an `AcousticSource_t`.
4. `new_solver` allocates the solver named by `method`; `init` moves the sources into it, checks that they lie in the water column (`check_sources_in_water`) and frees `config`.
5. At the end of its branch the emitter is moved (`MOVE_ALLOC`) into a `CLASS(Emitter_t)` variable, which is handed to `save_parameters` and deallocated before `run()`.
6. `run()` executes every enabled post-processing step; `update_parameters()` closes the `.dat` summary.

The solvers never know what produced the sources. Which steps run is decided by the parser through the `do_run` flags: `do_all = .true.` switches every step on, and for a wind farm `read_do_config` always switches sphere and decay off (there is no single structure to enclose or to follow downstream). The sphere is a no-op for every method but the image method (§7.3).

### 1.1 Input files

`read_input_files` reads `Solver.nml` and `Environment.nml` from the input directory and then three files that it picks by name, so that every emitter and every method has its own file:

| File | Chosen by | Files |
|---|---|---|
| `Emitters/<simulation_type>.nml` | `simulation_type` of `Solver.nml` | `WindTurbine.nml`, `WindFarm.nml`, `FrequencySource.nml`, `TimeSource.nml` |
| `ComputationSetups/<simulation_type>_CS.nml` | `simulation_type` | `WindTurbine_CS.nml`, `WindFarm_CS.nml`, `FrequencySource_CS.nml`, `TimeSource_CS.nml` |
| `AcousticMethods/<method>.nml` | `method` of `Solver.nml` | `ImageMethod.nml`, `AnalyticalNormalModes.nml` |

The values of `simulation_type` are `'WindTurbine'`, `'WindFarm'`, `'FrequencySource'` and `'TimeSource'` (the name of the emitter module, without `_t`), and those of `method` are `'ImageMethod'` and `'AnalyticalNormalModes'`. They are compared exactly (capitals included) and they are also the file names, so a new emitter or method needs no change in the path logic, only its `CASE` in `read_input_files` and in `main.f90`, and its files. The folder names are fixed in `read_input_files`; the directory itself is the command-line argument.

A namelist ends at the first `/` that is not inside a quoted string or a `!` comment, so comments in the `.nml` files start with `!` (a `//` "comment" silently ends the group: what follows is never read). The `/` must be followed by a newline: without it gfortran reads the values and then reports `End of file`, which stops the run; the shipped files keep a blank line and a comment after the group.

The shipped files (13/10/2026) run the `TimeSource` case, a monopole, with `method = 'AnalyticalNormalModes'` and `extra_label = 'None'` in `Solver.nml`, and `ImageMethod.nml` has `attenuation_bottom = 1.0`. The readers keep their own defaults: `attenuation_bottom` is also 1.0 when it is left out (`read_ImageMethod_config`, and the initial `lw_R` of `ImageMethod_t`), while `method` stays `'ImageMethod'` (`read_solver_config`).

### 1.2 Input modes of the manual sources

The two manual emitters accept their data written in the namelist or read from text files (`read_table`: blanks, tabs or commas; `#` and `!` lines are comments; a bad row stops the run with its line number). A way is chosen by giving its values (an entry of a list, a value other than `-999.0`, a file name other than `'None'`, an expression other than `'0'`); two ways of the same item stop the run with a message that names both, and with none the written list is taken, so an incomplete input of the old format gives the same error as before.

| Item | Written in the namelist | Read from files or generated |
|---|---|---|
| Positions (both emitters) | `source_xyz(i_s, 1:3)`, at most `MAX_PARSE_SRC` = 20 sources | `source_file`: one row `x y z` (or `x y z order`) per source, any number; `Nsources` is `-999` or the number of rows |
| Multipole orders (both emitters) | `source_order(i_s)`: 0 monopole, 1 dipole (an entry not given is 1), at most 20 sources | the fourth column of `source_file` (`x y z order`); with three columns every source is a dipole |
| Frequencies (`FrequencySource`) | `freqs` (`Nfreqs` values), at most `MAX_PARSE_FREQ` = 200 | a grid of `Nfreqs` (>= 2) frequencies from `fmin` (> 0) to `fmax`, `freq_spacing = 'linear'` (`fmin + (i-1) df`) or `'log'` (`10**(log10 fmin + (i-1) d)`, as `run_decay`); or `freq_file`, one frequency per row, any number |
| Forces (`FrequencySource`) | `force_re(i_f, i_s, 1:3)` (mandatory) and `force_im` (0 by default), at most 200 frequencies and 20 sources; a monopole uses only `force_re(i_f, i_s, 1)` and `force_im(i_f, i_s, 1)`, its components 2-3 may be left out (they become 0) | `force_re_expr(i_s, 1:3)` and `force_im_expr`, text expressions of `f` [Hz] ('0' by default, at most 20 sources); or `force_re_file` and `force_im_file` (`'None'`: imaginary part 0), one row per frequency and `3*Nsources` columns `Fx1 Fy1 Fz1 Fx2 ...`. A monopole uses only the first component (`force_re_expr(i_s, 1)`, `force_im_expr(i_s, 1)`; in the files its first column, 0 in the other two) |
| Time axis (`TimeSource`) | `t_start`, `t_end`, `Nt`: `Nt` samples of `[t_start, t_end)` | `time_file`, uniformly spaced times |
| Forces (`TimeSource`) | `force_expr(i_s, 1:3)`, text expressions of `t` [s], at most 20 sources; a monopole has `force_expr(i_s, 1)` = Q(t) and '0' in the other two | `force_file`, `Nt` rows and `3*Nsources` columns; a monopole has Q(t) in the first column of its three |

The 20/200 limits are the sizes of the parse buffers (`Kinds`), so they bind only what is written by hand. Frequencies must be positive and strictly increasing in every mode: the reader checks the list and the limits of the grid, `FrequencySource % init` checks the vector it builds (rounding, values of a file). Expressions are compiled by the reader (syntax errors and unknown names, e.g. `t` in an expression of `f`, before any work) and evaluated by `init` on the final axis, with a finite check that names the input, source, component and the first bad frequency or time. Positions of a file go through the same checks as the written ones (`check_source_positions`: finite, `z <= 0`; the solver then checks the seabed).

A monopole takes its volume velocity Q [m^3/s] from the first component of the force inputs (the convention of `SourceTypes`, §4.3). `source_order` given together with `source_file` stops the run (the orders go in its fourth column), and so do entries of `source_order` beyond `Nsources`, an order other than 0 or 1 (`check_source_orders`), a fourth column that is not exactly 0 or 1, and a `source_file` with a number of columns other than 3 or 4 (counted on its first data row by `count_table_columns`; `read_table` then stops at a row with another number). A monopole with a second or third component other than 0 stops the run in `init` (`check_monopole_components`); `TimeSource` checks it on the series, before the FFT, whose mean removal would hide a constant.

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
│   ├── ExpressionParser.f90     ExpressionParser
│   ├── Sorting.f90              Sorting
│   ├── Timing.f90               Timing
│   ├── SignalProcessing.f90     SignalProcessing
│   ├── SpecialFunctions.f90     SpecialFunctions
│   └── SeawaterAbsorption.f90   SeawaterAbsorption
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
│   ├── ManualSources.f90        ManualSources
│   ├── FrequencySource.f90      FrequencySource
│   ├── TimeSource.f90           TimeSource
│   ├── OpenFASTReader.f90       OpenFASTReader
│   ├── StructureModel.f90       StructureModel
│   ├── DrivetrainExcitation.f90 DrivetrainExcitation
│   └── TurbineTypes/
│       ├── DTU10MWMonopile.f90  DTU10MWMonopile
│       ├── DTU10MWFloating.f90  DTU10MWFloating
│       ├── SAITEC2MWFloating.f90 SAITEC2MWFloating
│       └── MYSE3MWMonopile.f90  MYSE3MWMonopile
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
| `SourceTypes` | `AcousticSource_t`, `GeometricData_t`, `AxisGroup_t`, the multipole labels `SOURCE_MONOPOLE`, `SOURCE_DIPOLE` and the strength convention of each order (force of a dipole, volume velocity of a monopole), and `group_by_axes`, the only procedure (it groups the sources into vertical axes from their own `xyz` and `order`) | Other procedures, configuration |
| `ConfigTypes` | Every `Config_*_t` filled in by the parser, including `Config_Method_t` and its extensions | Procedures, run-time state |
| `FFTW3` | C bindings of FFTW3 | Anything else |
| `StringUtils` | `to_upper`, `split_tokens`, `is_data_line` | File or HDF5 I/O |
| `FileUtils` | `read_full_line`, `swap_extension`, `ensure_parent_dir`, `move_text_file`, `delete_file` | HDF5, format-specific parsers |
| `HDF5Utils` | `save_variable` (dispatcher over the generic `save_to_hdf5`), `read_hdf5_real_1d`. The only module that does `USE HDF5` | Dataset names or knowledge of what is saved |
| `CSVReader` | `read_curve` (wind-speed lookup tables, with the private row parser `parse_csv_row`), `read_table` (numeric text tables of the manual sources), `count_table_columns` (values of the first data row, for a table with an optional column) | Turbine-specific files |
| `ExpressionParser` | `Expression_t`, `compile_expression(text, expr, message, variable)` (recursive descent to postfix; the variable is `t` unless named), `evaluate_expression(expr, x, values)` (vectorised), `function_list` | Knowledge of what the expression means |
| `Sorting` | `sort_real_c` (C `qsort`) | Anything else |
| `Timing` | `format_elapsed` | Anything else |
| `SignalProcessing` | `compute_rfft` (mean removal, Hann window, FFTW), `HANN_ENBW` | Turbine or solver knowledge |
| `SpecialFunctions` | `compute_hankel_complex` (public; `H0` and `H1` are optional and only the ones present are evaluated), the private `bessel_k0`, `bessel_k1` | Method-specific kernels |
| `SeawaterAbsorption` | `absorption_ainslie_mccolm(f, T, S, pH, z)` (elemental, dB/km: Ainslie & McColm 1998, Eqs. 1-3, f in Hz, z in m), `DB_PER_NEPER` | Solver state, the depth convention (the solver chooses z) |
| `Args_parser` | `get_input_path` (command-line directory of the namelists), `read_input_files` (finds the files of §1.1 and calls one reader per file), the entry of the parser | Readers of individual namelists |
| `InputEmitters` | `read_turbine_config`, `read_farm_config`, `read_frequency_source_config`, `read_time_source_config`, the private `read_source_file` (positions and orders from a file, both manual readers) and `join_given` | Solver or step settings |
| `InputSolver` | `read_solver_config`, `read_environment_config` (with the absorption: `absorption`, `temperature`, `salinity`, `pH`, their defaults and checks), one `read_X_config` per method | Emitter or step settings |
| `InputPostProcessing` | `read_computation_config`, `read_do_config`, one reader per `run_*` step | Emitter or solver settings |
| `InputChecks` | `check_time_window`, `check_fits`, `check_file_exists`, `check_all_or_none`, `check_count`, `check_source_positions`, `check_source_orders` | Namelist reading |
| `Emitter` | Abstract `Emitter_t` and its default `get_enbw` | Turbine or farm specifics |
| `WindTurbine` | Abstract `WindTurbine_t`: shared parameters, OpenFAST reading, `filter_frequencies`, the `Emitter_t` procedures of a turbine, `save_turbine_params`, `write_turbine_txt` | Data of a concrete turbine |
| `TurbineTypes/*` | One concrete turbine each: its geometry, `compute_force`, `compute_mass`, `excitation_spectrum`, its saved parameters and impedance correction | Generic routines (they go to `StructureModel`, `DrivetrainExcitation`, `OpenFASTReader`) |
| `TurbineFactory` | `new_turbine`: one `CASE` per turbine type | Anything else |
| `WindFarm` | `WindFarm_t`: a container of turbines (`TurbineSlot_t`) merged into one `AcousticSource_t` | Turbine-type specifics |
| `ManualSources` | Abstract `ManualSources_t`: the sources as the solver receives them (`Nsources`, `order`, `Freqs`, `xyz`, `Force`), the labels set by the children, and the shared `to_source`, `save_params`, `write_txt`, `check_monopole_components` | How the forces are given |
| `FrequencySource` | `FrequencySource_t`: the frequency axis (list, grid, file) and the forces (written, expressions of `f`, files) of `FrequencySource.nml`, and its `.dat` block | Time-domain forces, turbine or farm specifics |
| `TimeSource` | `TimeSource_t`: the time axis, the force series (expressions of `t` or a file), their spectrum with `compute_rfft`, its datasets and `.dat` block | Frequency-domain input, turbine or farm specifics |
| `OpenFASTReader` | SubDyn summary and `.out`/`.outb` channel readers | Physics |
| `StructureModel` | Geometry and hydrodynamic corrections of the members: `remove_duplicate_nodes`, `divide_span`, `integral_elipse_z`, `alpha_hankel` | Data of a concrete turbine |
| `DrivetrainExcitation` | Drivetrain spectra (`drivetrain10MW_excitation_spectrum`, `drivetrain2MW_excitation_spectrum` and its helpers, `drivetrain3MW_excitation_spectrum` and the reader of its force files `read_drivetrain3MW_force_file`), `generate_timeseries_banded_sines`, `filter_non_useful_freqs` | Turbine state |
| `AcousticSolver` | Abstract `AcousticSolver_t`, its deferred interfaces, the interfaces of the submodule procedures, `check_observers_distances`, `check_sources_in_water`, `set_absorption` (tabulates `alpha_dB`, `alpha` and `z_alpha`), `compute_pressure_masked`, `pressure_norm`, the water-column queries (`z_seabed_base`, `in_water_base`) | Method kernels, emitter types |
| `AcousticSolver_steps` | `run()`, every `run_*` step, `dummy_run_sphere` | Reads of `depth_base` (use the queries), method kernels |
| `AcousticSolver_output` | `resolve_save_path`, `check_update`, `update_parameters`, `save_shared_parameters`, `write_txt_water` and their private helpers | Pressure computations |
| `SolverFactory` | `new_solver`: one `CASE` per method | Anything else |
| `methods/*` | One acoustic method each, extending `AcousticSolver_t` | Emitter types (they only see `Emitter_t`) |

---

## 3. Dependency layers

```
Kinds  ->  types  ->  utils  ->  input / emitters  ->  solvers  ->  app
```

**Rule:** a module never `USE`s a module of a higher layer. `input` and `emitters` are on the same level and never use each other; both use `utils` (e.g. `read_table` and `ExpressionParser` serve `InputEmitters`, `FrequencySource` and `TimeSource`). Inside a layer the order is the one of the Makefile.

| Layer | Modules |
|---|---|
| Kinds | `Kinds` |
| types | `SourceTypes`, `ConfigTypes` |
| utils | `FFTW3`, `StringUtils`, `FileUtils`, `HDF5Utils`, `CSVReader`, `ExpressionParser`, `Sorting`, `Timing`, `SignalProcessing`, `SpecialFunctions`, `SeawaterAbsorption` |
| input | `InputChecks`, `InputEmitters`, `InputSolver`, `InputPostProcessing`, `Args_parser` |
| emitters | `OpenFASTReader`, `StructureModel`, `DrivetrainExcitation`, `Emitter`, `WindTurbine`, `DTU10MWMonopile`, `DTU10MWFloating`, `SAITEC2MWFloating`, `MYSE3MWMonopile`, `TurbineFactory`, `WindFarm`, `ManualSources`, `FrequencySource`, `TimeSource` |
| solvers | `AcousticSolver`, `AcousticSolver_output`, `AcousticSolver_steps`, `ImageMethod`, `AnalyticalNormalModes`, `SolverFactory` |
| app | `main` |

Other rules:
- Every module is `PRIVATE` by default with an explicit `PUBLIC ::` list of what other modules use (`Kinds` is the exception: everything in it is public). There are no umbrella or compatibility modules.
- `USE Module, ONLY: ...` is preferred; `-Wall` reports unused imports.
- A submodule depends on its parent's object in the Makefile; gfortran writes its `.smod` into `build_make/`. fpm discovers `src/**` by itself. The Makefile lists, for every object, the objects of the modules it uses: a missing one gives a build that works only in some orders (and a `make -j` that fails at random).
- `AcousticSolver` cannot `USE` a method (that would be a cycle); this is why `run_sphere` is a no-op in the parent and overridden in `ImageMethod`.

---

## 4. Contracts

### 4.1 The `Emitter_t` family (`emitters/`)

```
Emitter_t (ABSTRACT)                to_source, save_params, write_txt DEFERRED; get_enbw = -1
├── WindTurbine_t (ABSTRACT)        + 6 DEFERRED; init with the alias init_base
│   ├── DTU10MWMonopile_t
│   ├── DTU10MWFloating_t
│   ├── SAITEC2MWFloating_t
│   └── MYSE3MWMonopile_t
├── WindFarm_t                      turbines(:) of TurbineSlot_t (CLASS(WindTurbine_t))
└── ManualSources_t (ABSTRACT)      to_source, save_params, write_txt; aliases *_base; init DEFERRED
    ├── FrequencySource_t           init, write_txt
    └── TimeSource_t                init, save_params, write_txt, get_enbw
```

`Emitter_t` (`emitters/Emitter.f90`):

| Procedure | Kind | Meaning |
|---|---|---|
| `to_source(self, config, sources)` | deferred | Fills the `AcousticSource_t` handed to the solver: `xyz`, `Force`, `Freqs`, `order`, `geometry` |
| `save_params(self, save_path)` | deferred | Writes the emitter datasets into the HDF5 output |
| `write_txt(self, unit)` | deferred | Writes the emitter block of the `.dat` summary |
| `get_enbw(self)` | default `-1.0` | Equivalent noise bandwidth, in bins, of the window applied to the forces; `-1` means "forces not obtained from a windowed FFT". `WindTurbine_t`, `WindFarm_t` and `TimeSource_t` return `HANN_ENBW` (1.5); `FrequencySource_t` keeps `-1` |

`init` is not part of the `Emitter_t` contract because every emitter is built differently, so `main.f90` calls it through a variable of the concrete type (§5.2).

**`WindTurbine_t`** implements the contract with `turbine_to_source`, `save_params_turbine` (`save_turbine_params` + `save_specific_params`), `write_txt_turbine` (`write_turbine_txt`) and `get_enbw_turbine`. It is itself a template for the turbine types:
- **Shared:** `init`, `read_input`, `translate_align_with_WindDir`, `filter_frequencies` (rpm from `read_curve`, the type's `excitation_spectrum`, synthetic signal, frequency mask). By default every peak is a band of 50 sine components, each keeping the bins within one bin of it; a type whose drivetrain lines are pure tones sets `pure_tones = .true.`, and then each peak is a single line that keeps the bins within 1.5 bins (the three nearest, the main lobe of the Hann window for a tone on a bin, also when rounding or a window that is not a whole number of periods puts it slightly off the bin).
- **OpenFAST paths:** by default `init` builds `OP_output_dir // OP_subdir // rootname // '.outb'` (`.out`) and `// '.SD.sum.yaml'`. `OP_subdir` is the folder of the turbine type inside `OP_output_dir`, hardcoded by each type: `DTU10MWMonopile/`, `DTU10MWFloating/`, `SENVION2MWFloating/` (`SAITEC2MWFloating`), `MYSE3MWMonopile/`; `init` stops if a type leaves it unset. A full path in `OP_output_filename` or `OP_SubDyn_filename` is used as given.
- **Deferred:** `compute_force`, `compute_mass`, `excitation_spectrum(self, rpm, freqs_amp, keys)`, `save_specific_params`, `write_specific_params_txt`, `get_impedance_corrected_force`.
- A turbine type overrides `init` to set its own data (`path_rpm`, `OP_subdir`, `is_monopile`, ...) and calls the shared one as `CALL self % init_base(config)`: `init_base` is a second binding of the same procedure (`PROCEDURE :: init_base => init`), which the child does not override, so it always reaches the parent version.

**`MYSE3MWMonopile_t`** (MySE 3.0-135 of Beleolico, OpenFAST case `Beleolico_offshore_turbine`) follows `DTU10MWMonopile_t`: the pile is SubDyn members M1-M7 from the seabed (z = -14 m) to z = 0 (`Nmembers = 7`, `Nnodes = 5`, also in the `SDtwr` runs, whose tower members come after them), with constant diameters 5.0 / 4.88 m, steel and no internal water. Its drivetrain lines are the pure tones of the prescribed nacelle force of the runs (`pure_tones`): `drivetrain3MW_excitation_spectrum` reproduces the lines of `scripts/07_drivetrain_timeseries.py` (families A-D, frequencies rounded to the 1/T_win grid of the runs, close lines merged). With `filter_freqs`, `init` reads `rpm` and `T_win` from `<rootname>.force.json` (next to the `.outb`) and stops if the rpm of `wind_speed_curves_MYSE_3MW/rpm_ws.csv` at `WindSpeed` does not give the same lines: `WindSpeed` must be the wind of the run (`operating_points.csv` of the case).

**`WindFarm_t`** implements it with `farm_to_source_WindFarm`, `save_params_WindFarm`, `write_txt_WindFarm` and `get_enbw_WindFarm`; turbine `i` owns the sources `src_offset(i)+1 : src_offset(i+1)` of the merged set. It is not abstract.

**`ManualSources_t`** (`emitters/ManualSources.f90`, ABSTRACT) holds what both manual emitters share: point monopoles and dipoles that, whatever way their input is given, reach the solver as complex amplitudes per frequency (the force of a dipole, the volume velocity Q of a monopole in component 1; §4.3).
- Components: the labels `kind_name` (the `Emitter` dataset), `title` (of the `.dat` block) and `name_suffix` (end of the default file name), set by every child; and `Nsources`, `order`, `Freqs`, `xyz`, `Force (Nfreqs, Nsources, 3)`.
- Shared bindings: `to_source_ManualSources` (copies the arrays, because `save_params` and `write_txt` run after the solver has moved those of the `AcousticSource_t`; centre and origin = mean `(x, y)` of the sources, direction 0, folder `./source_acoustic_data/`, name `Source_Ns<N><name_suffix>`), `save_params_ManualSources`, `write_txt_ManualSources` (title, counts, the `Multipole` line and the positions, with an `order` column for a mixed set, §6.2), and `check_monopole_components` (stops the run if a monopole has something other than 0 in component 2 or 3; every child calls it in `init` on what it reads). The message of `to_source` names the orders handed over: point dipoles (as before), point monopoles, or point monopoles and dipoles.
- `init(self, config)` is DEFERRED: every child builds its data there and sets the three labels.
- **The `*_base` aliases.** A child that adds datasets or `.dat` lines overrides `save_params` or `write_txt` and calls the parent version through `save_params_base` and `write_txt_base`. The direct form `CALL self % ManualSources_t % write_txt(unit)` is rejected by gfortran 13 when the parent is ABSTRACT, so the alias (a second, never-overridden binding of the parent procedure, as `init_base` in `WindTurbine_t`) is the project's idiom.

**`FrequencySource_t`** (`init`, `write_txt`): `init` copies the positions and the orders (`source_order` of the configuration) and builds the frequency vector (`build_frequency_axis`: list, grid or the rows of `freq_file`) and the forces (`build_forces`: written, expressions of `f` evaluated on the final vector, or the columns of `force_re_file`/`force_im_file`) of §1.2, then calls `check_monopole_components` on `Force`. `name_suffix = '_Nf<Nfreqs>'`.

**`TimeSource_t`** (`init`, `save_params`, `write_txt`, `get_enbw`): `init` takes the orders from the configuration, builds the time axis (printing `dt`, the record length, `df` and the Nyquist frequency) and the series `Ft (Nt, Nsources, 3)`, calls `check_monopole_components` on `Ft` (before the FFT, whose mean removal would hide a constant component) and takes the series to the frequency domain with `compute_rfft` (mean removed, Hann window, zero frequency dropped), keeping the band `fmin`-`fmax` and one bin out of `skipf`. `name_suffix = '_Nt<Nt>'`; `get_enbw` returns `HANN_ENBW`.

### 4.2 `AcousticSolver_t` (`solvers/AcousticSolver.f90`)

Deferred procedures, implemented by every method:

| Procedure | Contract |
|---|---|
| `init(self, sources, config)` | Moves `xyz`, `Force`, `Freqs` and `order` out of `sources` (`MOVE_ALLOC`), stops if a source has an order without a kernel (both methods implement `SOURCE_MONOPOLE` and `SOURCE_DIPOLE`), copies `env`, `run_inputs` and the solver settings (`save_path`, `extra_label`, `block_size`, `debug`, `in_cluster`, `update_existing`, `mask_distance`), reads its own `Config_X_t`, calls `check_sources_in_water` once its environment (seabed included) is set and `set_absorption` (the frequencies are known); `AnalyticalNormalModes` also calls `sources % group_by_axes()` right after `check_sources_in_water`; calls `build_save_path`, prints `get_name()`, and frees `sources` and `config` |
| `get_name(self)` | Returns the human-readable name and stores it in `solver_name` (written as the `Method` dataset) |
| `compute_pressure(self, observers, block_size, total_pressure)` | Pressure `(Nfreqs, Nobs)` at `observers (Nobs, 3)` for every source frequency. The parallel blocks split the observers, never the sources. One kernel per multipole order present is called (`has_dipoles`, `has_monopoles`); each one skips (`cycle`) the sources of the other orders, and its result is accumulated before the next kernel is called. The ANM kernels also receive `sources % axis` and visit the sources by vertical axes: the horizontal distance, the Hankel functions of every mode (with their absorption factor) and the evanescent cutoff are evaluated once per observer, axis and frequency and reused by every source of the axis (an axis without sources of the kernel's order is skipped); the parallel blocks are those of before. With dipoles only, only the dipole kernel runs, so a set of dipoles is computed as before the monopoles (the effect on the bits is in §8). With `env % absorption` the kernels use `alpha(f)` [Np/m]: IM the complex wavenumber k + i*alpha on every path (`exp(-alpha*R)`, and `1/R + alpha - i*k` in the dipole); ANM `exp(-delta_m*r)` on the Hankel functions of every propagating mode, `delta_m = Im sqrt((k + i*alpha)^2 - k_zm^2)`; without it the loops of before run unchanged. Called only through `compute_pressure_masked` |
| `save_parameters(self, emitter)` | In this order: `check_update(emitter)`, its own datasets, `save_shared_parameters(emitter)`, `emitter % save_params`, `write_txt_summary(emitter)` |
| `build_save_path(self)` | Calls `resolve_save_path(method_tag)` with its tag, e.g. `_IM30`, `_ANM24` |
| `write_txt_summary(self, emitter)` | Writes the `.dat` header: title, `emitter % write_txt`, environment (ending with `write_txt_water`) and solver blocks |

Water-column queries (non-deferred, overridable, implemented in the parent):

| Binding | Default | Use |
|---|---|---|
| `z_seabed(x, y)` | `z_seabed_base`: `-depth_base` | Bottom of the water column at `(x, y)` |
| `in_water(points)` | `in_water_base`: `z_seabed(x, y) - tol_z <= z <= tol_z`, `tol_z` = 1e-9 m | Mask of the points `(N, 3)`, observers or sources, inside the water column. It asks `z_seabed`, so a method with its own seabed only overrides `z_seabed` |

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
- `check_sources_in_water()`: stops the run, naming the first offender and the seabed there, if a source is outside the water column (`in_water` on `sources % xyz`). The input readers only check that `z <= 0`, because the seabed belongs to the method. Every `init` must call it (§5.3).
- `save_shared_parameters(emitter)`: writes `Method` (`solver_name`) and, if `emitter % get_enbw() > 0`, `ENBW`; neither goes to the `.dat`. Then the water: `T_wat`, `S_wat`, `pH_wat` when given (not `DEFAULT_R`), and `z_alpha`, `Alpha_dBkm` with `env % absorption`.
- `set_absorption()`: with `env % absorption`, tabulates the volume absorption at every frequency (`alpha_dB` [dB/km], `alpha` [Np/m]) at `z_alpha = depth_base/2` and prints it; without it, does nothing. Every `init` must call it (§5.3).
- `write_txt_water(u)`: the water lines at the end of the `ENVIRONMENT` block of the `.dat` (temperature, salinity, pH when given; formula, depth, alpha at f_min and f_max with absorption); nothing when none is set, so the summaries without absorption keep their content.
- `run()`: executes the enabled steps in a fixed order (spectrums, polar, decay, line, sphere, sliceXY, sliceXZ, sliceVertical, cylinder).
- `run_sphere`: a no-op in the parent (`dummy_run_sphere`, which switches `do_sphere` off); `ImageMethod` overrides it with the free-field sphere.
- `resolve_save_path`, `check_update`, `update_parameters`: see §6.3.

### 4.3 `AcousticSource_t` and `GeometricData_t` (`types/SourceTypes.f90`)

| Field | Shape | Meaning |
|---|---|---|
| `xyz` | `(Nsources, 3)` [m] | Source coordinates |
| `Force` | `(Nfreqs, Nsources, 3)` [N] | Complex source forces (radiation-impedance corrected for turbines); a monopole has its volume velocity Q [m^3/s] in component 1 and 0 in components 2-3. The convention of each order is in `SourceTypes`: peak complex amplitudes, time factor `exp(-i*w*t)` |
| `Freqs` | `(Nfreqs)` [Hz] | Frequency vector |
| `order` | `(Nsources)` [-] | Multipole order label: `SOURCE_MONOPOLE = 0`, `SOURCE_DIPOLE = 1` (1 for every source of a turbine or farm); saved as `Source_order` by the manual sources only |
| `axis` | `(Naxes)` of `AxisGroup_t` | Vertical axes of the sources, not filled in by the emitters: built by `group_by_axes(tolerance)` in the `init` of the method that uses it (ANM). A source joins the first axis whose first source lies within `tolerance` in (x, y) (1 cm by default, `<=`, so 0 asks for the same x and y), or else opens a new axis. Each axis has `x`, `y` (those of its first source, at which every source of the axis is evaluated), `Nsources`, `Nmonopoles`, `Ndipoles` and `idx` (indices into `xyz`, `Force`, `order`, ascending). The axes are numbered in order of first appearance. A DTU 10 MW monopile is one axis of 25 sources; a farm of monopiles, one per turbine; scattered manual sources, one per source. Not saved |
| `geometry % center` | `(2)` [m] | Reference centre: turbine barycentre, or mean of the farm or of the sources |
| `geometry % origin` | `(2)` [m] | Reference origin: turbine axis, or mean of the farm or of the sources |
| `geometry % direction` | [deg] | Reference direction: wind direction, circular mean of the farm, or 0 for point sources |
| `geometry % out_dir`, `name_tag`, `name_suffix` | | Default output name `out_dir // name_tag // <method tag> // name_suffix // ".hdf5"`; folders `./turbine_acoustic_data/`, `./farm_acoustic_data/`, `./source_acoustic_data/` |

What each consumer takes from them:

| Consumer | Uses |
|---|---|
| `compute_pressure` (methods) | `xyz`, `Force`, `Freqs`, `order` (one kernel per order); ANM also `axis` |
| `check_observers_distances` | `xyz` |
| `check_sources_in_water` | `xyz` |
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
1. `src/emitters/TurbineTypes/X.f90`: `X_t EXTENDS(WindTurbine_t)` with its geometry and the deferred procedures, `excitation_spectrum` included. Its `init` sets `path_rpm` (and any data path), `OP_subdir` (its folder in `OP_output_dir`, where its OpenFAST runs go), `is_monopile` (and `pure_tones` if its drivetrain lines are pure tones), and calls `CALL self % init_base(config)`.
2. If its drivetrain is new, add the spectrum to `emitters/DrivetrainExcitation.f90` and call it from `excitation_spectrum`.
3. One `CASE` in `new_turbine` (`emitters/TurbineFactory.f90`).
4. Optionally its short label in `type_alias` (`emitters/WindFarm.f90`) for the default farm file name.
5. Makefile: the object in `OBJS` and its rule, before `TurbineFactory`.

### 5.2 An emitter

**A manual source** (point sources whose forces are given by the user in a new way, e.g. an input HDF5): extend `ManualSources_t`.
1. `src/emitters/X.f90`: `X_t EXTENDS(ManualSources_t)` with:
   - `init(self, config)`: fills `Nsources`, `xyz`, `order` (from its configuration, e.g. `source_order`), `Freqs`, `Force` and sets `kind_name = 'X'`, `title` and `name_suffix`, and calls `check_monopole_components` on what it reads (on the series, before any FFT);
   - optionally `save_params` and `write_txt` overrides that call `save_params_base` / `write_txt_base` first and add what is its own. `write_txt` must list every input value that matters (or checksums of it at full precision): `check_update` compares this block, so an input missing from it would reuse an old output after being changed;
   - `get_enbw`, if its forces come from `compute_rfft`.
   `to_source`, the common datasets and the orders and positions in the `.dat` come from the parent.
2. Its configuration type in `types/ConfigTypes.f90` and its field in `Config_t`; its reader in `input/InputEmitters.f90` (reuse `read_source_file`, which also reads the orders of a fourth column, `check_source_positions`, `check_source_orders`, `check_file_exists`, `read_table`, `compile_expression`); its `CASE` in `read_input_files`; its files `Emitters/<X>.nml` and `ComputationSetups/<X>_CS.nml`, named exactly as its `simulation_type` (§1.1).
3. Its construction branch in `app/main.f90`: a typed variable (the abstract `active_emitter` can neither be allocated nor give access to `init`), `init`, `to_source`, and `MOVE_ALLOC` into `active_emitter` at the end of the branch.
4. Makefile: the object in `OBJS` and its rule, after `ManualSources`, with every module it uses among its prerequisites.
5. Plotting: `readers/hdf5_io.EMITTERS` maps the `Emitter` marker to the script that reads it (`'source'` for the manual sources); add what is specific to it in `readers/source.py` and `wrappers/source.py` (see the plotting README).

**Any other emitter** (a new kind of structure): `X_t EXTENDS(Emitter_t)` with `to_source`, `save_params` (which writes `Emitter = '<X>'` and everything the plotting needs, §6.1), `write_txt` (with the same rule) and, if needed, `get_enbw`; then items 2 to 5, with a reader and a wrapper of its own in `plotting/`.

The planned emitters are in §7.1.

### 5.3 An acoustic method
1. `src/solvers/methods/X.f90`: `X_t EXTENDS(AcousticSolver_t)` with the deferred procedures (§4.2); `save_parameters` calls `save_shared_parameters`.
   - One kernel per multipole order it supports (the convention of each order is in `SourceTypes`); each kernel skips (`cycle`) the sources of the other orders.
   - If part of its kernel depends only on the horizontal source-observer distance (as the Hankel functions of ANM), it can call `CALL self % sources % group_by_axes()` in `init` and evaluate that part once per axis (`sources % axis`, §4.3).
   - The check of `init` lists the orders it supports and stops on any other.
   - `compute_pressure` calls the kernel of each order present and accumulates its result before calling the next one.
   - If it builds images, every node carries the order of its parent source (as `ImageSystem_t % order`).
2. Its `init` calls `CALL self % check_sources_in_water()` once the environment of the method is set (its seabed included) and before anything is saved. The compiler cannot enforce it: without that line the sources are never checked, and one below the seabed gives finite but wrong pressures. It also calls `CALL self % set_absorption()`, and with `env % absorption` its kernels apply `alpha` (wavenumber k + i*alpha) and its `write_txt_summary` calls `write_txt_water`: not enforced either, and without it the file says `Alpha_dBkm` while its pressures lack the absorption (and the plotting no longer applies it).
3. `Config_X_t EXTENDS(Config_Method_t)` in `types/ConfigTypes.f90`.
4. `read_X_config` in `input/InputSolver.f90`, and its `CASE` in `read_input_files` (`input/Args_parser.f90`), with its own `.nml` in `input_files/AcousticMethods/<method>.nml`.
5. Its `CASE` in `new_solver` (`solvers/SolverFactory.f90`).
6. Makefile: the object before `SolverFactory`.
7. Nothing else: the `Method` dataset is written automatically.
8. The free surface is always the plane `z = 0`: the queries, the readers (`z <= 0`) and the steps assume it.

If the method has its own seabed, it overrides `z_seabed` only: `in_water` and `check_sources_in_water` ask it. It never overrides the `run_*` steps.

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

Complex arrays are stored as `(2, N...)` in Fortran order, i.e. `(..., N, 2)` (real, imaginary) when read with h5py; real arrays are also read with their axes reversed (e.g. `Source_Force` `(Nfreqs, Nsources, 3)` in Fortran is `(3, Nsources, Nfreqs, 2)` in h5py). Scalars are 1-element datasets; strings are fixed-length byte strings. Every step also writes `Freqs`.

| Writer | Datasets |
|---|---|
| `ImageMethod % save_parameters` | `N_images`, `c_wat`, `rho_wat`, `Lower_HBC`, `p_ref` |
| `AnalyticalNormalModes % save_parameters` | `Nmodes`, `H`, `c_wat`, `rho_wat`, `Lower_HBC`, `p_ref` |
| `save_shared_parameters` | `Method` (`"Images Method"`, `"Analytical Normal Modes"`), `ENBW` (1.5, only for forces from a windowed FFT: turbines, farms, `TimeSource`); `T_wat` [degC], `S_wat` [ppt], `pH_wat` only when given in `Environment.nml` (or set by `absorption`); `z_alpha` [m] and `Alpha_dBkm` `(Nfreqs)` [dB/km] only with `absorption` (their presence tells the plotting that the pressures include it) |
| turbine (`save_turbine_params`) | `case_type`, `WindSpeed`, `WindDir`, `Depth`, `AxisPos`, `BariPos` |
| turbine (`save_specific_params`) | `Case_type`, `Structure_nodes` plus the type's geometry (`D`, `rho_mat`, `wet_area`; `col_D`, `XsecA`, `XsecB`, `t`, `col_members`, `pon_members`; `Diam_y`, `Diam_z`, `float_members`, `joint_members`) |
| farm (`save_params_WindFarm`) | `Nturbines`, `Source_offset`, and every turbine dataset with the suffix `_i` (e.g. `WindSpeed_2`, `Structure_nodes_2`) |
| manual sources (`save_params_ManualSources`) | `Emitter` (`"FrequencySource"` or `"TimeSource"`), `Nsources`, `NFreqs`, `Source_xyz` `(Nsources, 3)`, `Source_order` `(Nsources)` (0 monopole, 1 dipole), `Source_Force` complex `(Nfreqs, Nsources, 3)` (the force [N] of a dipole; the volume velocity Q [m^3/s] of a monopole in component 1, 0 in 2-3) |
| `TimeSource` (`save_params_TimeSource`) | the above plus `Nt`, `Time_step`, `Source_time` `(Nt)`, `Source_Force_time` `(Nt, Nsources, 3)` (Q(t) of a monopole in component 1) |
| `run_spectrums` | `P_spectrums`, `Obs_spectrums` |
| `run_polar` | `P_polar`, `R_polar`, `Z_polar`, `Theta_deg_polar`, `Obs_polar`, `Center_polar` |
| `run_cylinder` | `P_cylinder`, `R_cylinder`, `Z_cylinder`, `Theta_deg_cylinder`, `Obs_cylinder`, `Center_cylinder`, `dA_cylinder` |
| `run_decay` | `P_decay`, `Z_decay`, `Distance_decay`, `Distances_decay`, `Obs_decay`, `Logspace_decay` |
| `run_line` | `P_line`, `P1_line`, `P2_line`, `Obs_line`, `Distance_line`, `Distances_line`, `Logspace_line` |
| `run_sliceXY` | `P_slicexy`, `Obs_slicexy`, `X_slicexy`, `Y_slicexy`, `Z_slicexy`, `Center_slicexy`, `Nx_slicexy`, `Ny_slicexy` |
| `run_sliceXZ` | `P_slicexz`, `Obs_slicexz`, `X_slicexz`, `Z_slicexz`, `Y_slicexz`, `Nx_slicexz`, `Nz_slicexz` |
| `run_sliceVertical` | `P_sliceVertical`, `Coords_sliceV`, `U_sliceV`, `Z_sliceV`, `Azimuth_sliceV`, `Center_sliceV` |
| `run_sphere` (IM) | `P_sphere`, `Obs_sphere`, `Theta_sphere`, `Z_sphere`, `R_sphere`, `Center_sphere`, `dA_sphere` |

Masked observers (closer than `mask_distance` to a source, or outside the water column) have NaN pressure; the grids keep their shape. The input mode of a manual source (list, grid, file, expressions) is not saved: the same physical input gives bit-identical datasets in every mode. The monopoles add no dataset: they are told apart by `Source_order`.

`Emitter` identifies the emitter that wrote a file and is written only by the emitters added after the refactor (the manual sources); the plotting recognises the older turbine and farm files by their datasets (`Nturbines`, `case_type`, ...). The water depth is not an emitter dataset: both solvers write `Lower_HBC` (`-depth_base`) and the plotting derives the depth from it.

### 6.2 The `.dat` summary

Written next to the HDF5 file (same name, `.dat`), with blocks between separators of 70 `=` or `-`:

1. Title `OFF_COUSTICS SIMULATION SUMMARY` and `HDF5 file: <path>`.
2. Emitter block (`emitter % write_txt`):
   - `TURBINE` and `TURBINE - SPECIFIC PARAMETERS`;
   - or `WIND FARM` (layout table) followed by `>>> Turbine i` and the blocks of every turbine;
   - or `FREQUENCY-DOMAIN SOURCES`: counts, the line `Multipole   : dipole | monopole | mixed`, positions (`F14.6`; only for `mixed`, a last column `order [-]`), the lines `Source file`, `Freq. grid`, `Freq. file`, the expressions or `Force file` when the input was not written in the namelist, then the title of the forces (followed, if there is a monopole, by the line `Monopoles   : volume velocity Q [m^3/s] in the Fx columns (Fy and Fz are 0)`) and the force table (`I4,I4,F12.4,6ES16.8`: index of frequency and source, frequency, Re and Im of Fx, Fy, Fz) if it has at most 200 frequencies and 20 sources, and checksums at full precision (`ES23.15`: sums of `f` and `f**2`; per source and component the sums of Re F, Im F, |F|^2, i Re F, i Im F, with `i` the frequency index) when it is larger or was not written by hand. An input written by hand keeps exactly the table it always had;
   - or `TIME-DOMAIN SOURCES`: counts, the `Multipole` line, positions (with the `order` column for `mixed`), `Source file` if used, the time axis (`N samples`, `First time`, `Time step`, `Skipf`, `Fmin`, `Fmax`, `Time file` if used), the title of the forces and, if there is a monopole, the same `Monopoles` line, the expressions or `Force file`, and per source and component `sum(F)`, `sum(F**2)` and `sum(i*F)` of the series.
3. `ENVIRONMENT` and `SOLVER: <name>` blocks (`write_txt_summary`). The `ENVIRONMENT` block ends with the lines of `write_txt_water`: `Temperature`, `Salinity` and `pH` (`F8.2`) when given and, with absorption, `Absorption       : Ainslie & McColm (1998)`, `Absorption depth` and `Alpha at f_min`, `Alpha at f_max` (`ES10.3` dB/km); so `check_update` rebuilds a file when the absorption or the water changes.
4. `POST-PROCESSING STEPS EXECUTED` (`update_parameters`): one block per step (`Spectrums:`, `Polar:`, ...) with the parameters actually used, defaults resolved.
5. A final `=` separator.

With dipoles only, the block of a manual source is the one it always had (`Multipole   : dipole`, no `order` column, no `Monopoles` line). The orders are part of the header, so `check_update` sees a change of order of any source.

Items 1-3 are the header compared by `check_update`.

### 6.3 `update_existing`

With `update_existing = .true.` in `Solver.nml`, an existing output of the same case is updated instead of rebuilt:
- `resolve_save_path` deletes any stale `.dat.prev`. If the HDF5 file exists, it is kept and its `.dat` is moved to `.dat.prev` (with `.false.` it is deleted).
- `check_update`, the first call of `save_parameters`, keeps it only if it holds the same case:
  - the same `Freqs` (same size, relative tolerance 1e-9);
  - the same `.dat` header (everything before the post-processing block, without the `Block size` and `In cluster` lines; the lines are compared up to 256 characters).

  Otherwise the file and its `.dat.prev` are deleted with a WARNING that gives the reason (`different frequencies`, the first different line of the header quoted old -> new, `different .dat header length`, ...), and the output is rebuilt.
- The steps run now overwrite their datasets; the other datasets stay.
- `update_parameters` writes the blocks of the steps run now and copies the blocks of the other steps from `.dat.prev`, tagged `(previous run)`, then deletes `.dat.prev`.

`Method` and `ENBW` are not in the `.dat`, so they do not take part in the comparison.

**Known gap (verified 08/10/2026).** Settings that change the results but are not in the header, or not at enough precision, are not seen: `UpperBC`, `LowerBC` and `epsilon` of `ImageMethod` (its `SOLVER` block does not print them), `mask_distance` (no `.dat` prints it), and changes of the speed of sound, the density or the depth below the `F8.2` of the `ENVIRONMENT` block. With the shipped `FrequencySource` case and IM (all steps, then `polar` only with one of these changed: `UpperBC` -1 -> 1, `LowerBC` 1 -> -1, `epsilon` 1e-6 -> 900, `mask_distance` 1 -> 900, `depth_base` 30.0 -> 30.004) the second run reported `Same case`, rewrote `P_polar` and kept the other steps of the old setting.

---

## 7. Planned extensions

Decided on 2026-10-02; the status of each item is given.

### 7.1 Manual sources
**Implemented:**
- `ManualSources_t` (ABSTRACT parent, 07/10/2026) with `FrequencySource` (06/10/2026: forces per frequency, no window, `get_enbw = -1`, no radiation-impedance correction) and `TimeSource` (07/10/2026: real series in time through `compute_rfft`, `get_enbw = HANN_ENBW`). `mask_distance` is set by the user in `Solver.nml`.
- Input modes of §1.2 (08/10/2026): positions from `source_file` in both; frequencies as a list, a grid or a file, and forces written, as expressions of `f` or from files, in `FrequencySource`.
- One multipole-order label per source (`order`, saved as `Source_order`).
- The plotting reads both (`plotting/source.py`): `ENBW` (Hann correction only if present), `Method`, and the bin width of uniformly spaced tones.
- Monopoles (13/10/2026):
  - convention: the strength of a monopole is its volume velocity Q [m^3/s], a peak complex amplitude in `Force(:, s, 1)`, with 0 in components 2 and 3 (`SourceTypes`); its free-field pressure is `p = -i*w*rho*Q * exp(i*k*R) / (4*pi*R)`;
  - kernels: `monopole_pressure_images` (IM; the images of a monopole come unchanged from `method_of_images`) and `monopolar_pressure_NM` (ANM), one kernel per order in each method, and a check of both orders in `init`;
  - input: `source_order` in both namelists, or a fourth column of `source_file` (§1.2), and `check_monopole_components` in `init`;
  - plotting: `plots/specific/sources.py` draws |F| [N] and |Q| [m^3/s] on separate axes.

**Planned:**
- An input HDF5 file for many sources with their forces (`xyz`, `order`, `Freqs`, `Force`); `HDF5Utils` only has `read_hdf5_real_1d` today. The text files of §1.2 already remove the limit on the number of sources and frequencies.
- Verification with the dipoles of a turbine, exported and read back as manual sources: they must reproduce its `P_*` datasets. The impedance-corrected forces handed to the solver are not saved in the turbine HDF5 today, so they must first be written out.

### 7.2 Acoustic methods
- `Config_Environment_t` will grow optional arrays: c(z), ρ(z), attenuation, bathymetry H(x, y), sediment. It keeps the base scalars, so the simple methods keep working. Each method checks in `init` that the data it needs is present. The volume absorption of the water is already a scalar option (`absorption`, 14/10/2026); the losses of a penetrable seabed, larger in shallow water, need a sediment.
- The boundaries belong to each method: constants in IM and ANM, arrays in complex methods. A method publishes the seabed by overriding `z_seabed`; `in_water` and `check_sources_in_water` follow it.
- Overriding a whole `run_*` is the exception (precedent: IM's `run_sphere`).

### 7.3 Sphere
The sphere stays IM-only: a no-op in the parent plus the IM override, because `AcousticSolver` cannot `USE ImageMethod` (cycle). A shared free-field kernel would make it generic for every method; that changes the outputs and is a future decision.

---

## 8. Verification protocol

Any change that is meant not to change the results follows `Programming_style_guidelines.txt` §8:

1. Keep a copy of the original sources and build the reference binary from them. Run everything in an isolated copy (never in `FORTRAN/`, `turbine_acoustic_data/`, `farm_acoustic_data/`, `source_acoustic_data/` or `SAITEC_acoustic_data/`): a run deletes and rewrites the outputs of the same name.
2. Builds: `make clean && make -j8`, `fpm build --flag "-O3 -fopenmp"`, and `make clean && make FFLAGS="-O3 -fopenmp -Wall"`, which must give no warning that the reference does not give (keep `-fopenmp` when overriding `FFLAGS`: without it the link fails with undefined `omp_get_wtime_`, `omp_get_max_threads_`, ...). The accepted `-Wall` warnings are the `-Wunused-dummy-argument` of procedures that ignore arguments by design:
   - `z_seabed_base` (`x`, `y`: the default seabed is flat);
   - `to_source_ManualSources` (`config`, which the `Emitter_t` interface requires);
   - with older gfortran versions also the unused `self` of `get_enbw` of `Emitter_t`, `WindTurbine_t`, `WindFarm_t` and `TimeSource_t`, and of `excitation_spectrum` of `DTU10MWMonopile_t` and `DTU10MWFloating_t` (9 warnings in all); gfortran 15.2 does not report an unused passed-object `self`, so it gives only the 3 above (checked 06/11/2026, the same with the sources of before).
3. A debug build, `make FFLAGS="-O0 -g -fopenmp -fcheck=all -fbacktrace"`, must run the cases without a runtime error.
4. Cases, with `OMP_NUM_THREADS=1` and the same `block_size` in both builds. A run is deterministic for given inputs, `block_size` and thread count, and the image method gives the same bits for any of them; the number of threads changes the observer blocks (`block_size` is clamped to `ceil(Nobs / threads)`), and the normal-mode results change in the last bits with the blocks (3.6e-16 of the peak pressure, measured 08/10/2026: with `AnalyticalNormalModes.f90` alone compiled with `-fno-tree-vectorize` the bits no longer depend on the blocks, so the cause is the auto-vectorisation of that module, most likely of the modal shapes tabulated over the observers of each block). The harness of `review/refactor_tools/` ran them (git-ignored and no longer on disk since 10/2026: the cases below are rebuilt by hand in an isolated copy):
   ```bash
   python refactor_cases.py setup   <reference_src_dir> <cases_root>
   python refactor_cases.py run     <binary> <cases_root> <tag>
   python refactor_cases.py compare <cases_root>/<tagA> <cases_root>/<tagB> [extra1,extra2,...]
   ```
   The cases are M1 (DTU 10 MW monopile, IM, all steps), M2 (the same with ANM), M3 (`update_existing`: polar + line, then sphere only, then `N_images = 20`), and L1, L2, L3 (DTU floating, SAITEC, 2-monopile farm, spectrums only). They take about 75 s. The manual sources need no OpenFAST data: their shipped cases (`FrequencySource.nml`, `TimeSource.nml`) with IM and ANM are the quick regression of every change in the input, emitter or solver layers.
5. Acceptance: every HDF5 dataset of the reference exists with the same dtype, shape and bytes; there is no extra dataset except the allowed ones; the `.dat` files are identical; stdout is identical except the `Total elapsed time` line. A difference is never accepted: its cause is found and fixed.
6. After editing a module, rebuild before running (`make` follows the dependencies; a binary older than its sources was the cause of the "segfault" of 07/10/2026, when `main.f90` still had the old `simulation_type` labels).

**Monopoles.** A change in the kernels, the image system or the input of the orders also repeats the monopole checks of 13/10/2026 (`OMP_NUM_THREADS=1`):
- a free-field monopole (IM, `N_images = 0`) against the closed form, 15-250 Hz and 203 observers: 1.2e-15 (relative, maximum; the dipole `Fx` of the same check: 4.8e-16);
- a pair of monopoles +-Q at +-d/2 (d = 0.02 m, `Q = i*F/(w*rho*d)`) against the dipole F, along x and along z, images included: 3e-7 to 1.4e-5 with IM and ANM;
- IM against ANM in the ideal guide (p = 0 above, rigid below, attenuations 1; ANM with 400 modes): a monopole differs by at most 5.0e-2, 2.4e-3, 2.4e-4, 2.4e-5 with `N_images` = 30, 300, 3000, 30000, as 1/N like the dipole `Fz` (1.5e-1, 8.3e-3, 8.3e-4, 8.3e-5);
- an unordered mixed set (dipole, monopole, dipole) against the sum of its sources run one by one: 1.4e-16 (IM) and 2.2e-16 (ANM);
- with dipoles only, the shipped `FrequencySource` and `TimeSource` cases (with the `.nml` files of before) and the DTU 10 MW monopile (`decay`) against the binary of before: identical with IM; with ANM identical except the `P_*` datasets, by at most 2.3e-16 (see below);
- the debug build (`-fcheck=all`) runs all of them without a runtime error, and `-Wall` gives only the 9 accepted warnings.

**Absorption.** A change in the absorption or in the kernels also repeats the checks of 14/10/2026 (`absorption = .true.`; `OMP_NUM_THREADS` 1 to 4):
- `Alpha_dBkm` against Ainslie & McColm in Python, 36 Hz to 20 kHz: at most 5.5e-16 (relative); the plotting's `absorption_coefficient`, the same;
- ANM monopole and dipole (H = 25 m, 46 Hz to 20 kHz, 31 observers to 3 km) against the same model written in Python (`scipy` Hankel functions): 4.5e-13 and 2.3e-12; against the Hankel functions of the complex `k_rm` (the exact solution with k + i*alpha): 1.2e-6 to 4.1e-4;
- one mode (`Nmodes = 1`): `|p_on / p_off| = exp(-delta_1 r)` to 3.3e-16;
- IM without images (`N_images = 0`) against the closed forms with k + i*alpha: 1.5e-16 (monopole), 3.6e-17 (dipole);
- IM (2000 and 4000 images) against ANM, monopole, 5 and 20 kHz, 0.1 to 1 km: 0.001 to 0.017 dB rms (0.012 to 0.055 dB without absorption); the effect of the absorption agrees within 0.07 dB;
- with `absorption = .false.` (and the template `Environment.nml` with the new lines), the shipped `FrequencySource` and `TimeSource` cases with IM and ANM, all steps, against the binary of before: identical HDF5, `.dat` and stdout;
- the debug build runs them without a runtime error, and `-Wall` gives only the 9 accepted warnings.

**Vertical axes.** A change in `group_by_axes` or in the loops of the ANM kernels repeats the checks of 06/11/2026 (`OMP_NUM_THREADS=1`, against the binary of before):
- a `FrequencySource` of 7 sources on 2 axes, both mixed (axis A at (0, 0): dipole, monopole, dipole; axis B at (40, 25): dipole, dipole, monopole, dipole), 3 frequencies, all steps, with and without absorption: `group_by_axes` gives 2 axes of 3 and 4 sources with the right counts, and every HDF5 dataset, the `.dat` and stdout are identical;
- a farm of 4 DTU 10 MW monopiles (2 x 2, 1000 m by 800 m, 11.4, 8.0, 15.0 and 4.5 m/s; 100 sources on 4 axes of 25, 1926 frequencies, 24 modes; spectrums, polar, line): identical;
- the shipped `TimeSource` with ANM, and the shipped `FrequencySource` and `TimeSource` with IM (which does not use the axes): identical;
- the debug build runs them without a runtime error, `-Wall` gives only the accepted warnings, and `fpm build` compiles.

The result is bit-identical because, within an axis, the sources are added in ascending order and every expression keeps its operands in the order of before: `coef_xy * H1 * F_dot_rhat` written as `coef_xy * F_dot_rhat * H1` changed the pressures by up to 3e-13 (relative). With the sources of an axis not contiguous (interleaved with those of another axis) the sum is reordered and the results change in the last bits; with a source within the tolerance but not exactly on its axis, its horizontal distance is that of the first source of the axis.

Time of the same farm (spectrums, polar, line, sliceXY, sliceVertical; 5300 observers): 145 s before and 11 s after with 16 threads; spectrums, polar and line with one thread: 58.6 s and 5.2 s. The gain grows with the number of sources per axis: each Hankel function is evaluated once per axis instead of once per source. The image method has nothing to share along an axis (its `sin`/`cos` depend on the distance to each node), so it is unchanged.

**Normal-mode references before 13/10/2026.** The 1/rho of the two terms of `dipolar_pressure_NM` is applied once to the tabulated source shapes `Psi_m(z_s)` and `Psi_m'(z_s)`, instead of twice inside the sum over the modes (2.2 % less time: 175.9 s against 179.9 s for the DTU 10 MW monopile, `decay`, one thread). The `P_*` datasets of the normal modes with dipoles changed in the last bit (at most 2.3e-16, relative); every other dataset, the `.dat`, stdout and the positions of the NaN are identical, and the image method is unchanged. An ANM reference made before that date is no longer bitwise comparable and must be regenerated; the IM references remain valid. Moving the loop over the sources outside was also tried: identical results and the same time (180.6 s), so the loop order was kept.
