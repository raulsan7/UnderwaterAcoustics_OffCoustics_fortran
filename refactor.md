# Refactor of the off-coustics Fortran solver

This document is the complete work order for one Claude Code session. Follow it literally, phase by phase.

Everything in it was agreed with the user (Raul Sanz) on 2026-10-02. Do not add, drop or "improve" anything that is not written here. If something here is ambiguous, contradicts the code, or cannot be done as written, stop and ask the user. Never improvise a different solution.

Only Phase 4 (§9) changes the outputs, and only in the ways it lists: two new HDF5 datasets and one new `Solver.nml` line. Phase 5 (§10) only adds documentation.

Everything else is a pure reorganisation. It must not change a single calculation, the program flow, or a single byte of the outputs.

---

## 0. Your role and how to talk to the user

- **Language:** talk to the user in Spanish. All code, comments, docstrings and documents stay in English, in his style.
- **Repository:** `/home/hp/simulations/openfast_simulations/off_coustics/FORTRAN`. It is the git root, on branch `main`, with remote `origin` = `github.com/raulsan7/New_fortran`.
- **Pacing:** work in **5 phases** (Phases 1 to 5, §6–§10), after a preparation step (Phase 0, §3).
  - At the end of every phase, produce the patches (§3.5), report, and **stop**.
  - The user applies the patches and commits them himself.
  - Only when he says so do you continue with the next phase.
- **Never** `git commit`, `git push`, or change branches in the user's repository.
- **Never** run `./off_coustics_fortran` in the repository.
- **Never** write into `SAITEC_acoustic_data/`, `turbine_acoustic_data/`, `farm_acoustic_data/` or `figures/`. Every build and every run happens in your scratchpad (§3).
- **Do not touch `plotting/`** (the Python package). Do not edit `input_files/*.nml` except where a phase says so.
- **Precedence:** the style guide `Programming_style_guidelines.txt` (repo root) is mandatory for every line you write or move. Read it completely before starting. Its golden rule applies to Phases 1–3 without exception. Phase 4 has explicitly listed exceptions.

### The confidential turbine file (`SAITEC2MWFloating.f90`)

`src/core/TurbineTypes/SAITEC2MWFloating.f90` contains data of another company.

- **Allowed:** treat it in this refactor exactly like the other turbine types. Move its generic routines out, update its `USE` lines, and move the file to the new folder.
- **Must stay in it:** `compute_force` and `compute_mass` (the routines with the private data).
- **Never:** commit or push it, or copy it out of this machine.
- **Patches:** its changes always go into a **separate** patch (§3.5), so the user decides what he does with them.
- **If an edit is refused:** if the permission system refuses an edit of this file, do not insist. Write the intended change as a diff `refactor_phaseN_SAITEC_manual.diff` in the repository root (next to `refactor.md`) and tell the user.

---

## 1. What the refactor achieves (summary)

1. **New folder and module layout.** Every routine lives in the module its role suggests (§6).
2. **Abstract `Emitter_t`.** Turbines, farms and future manual sources become interchangeable. The solvers no longer know turbines or farms (§7).
3. **Three preparations for new methods** (§8):
   - a factory `new_solver`;
   - a single `filter_frequencies` for every turbine type, through a deferred `excitation_spectrum` and the new `DrivetrainExcitation` module;
   - water-column queries (`z_surface`, `z_seabed`, `in_water`) so the `run_*` steps stay generic when methods with variable bathymetry arrive.
4. **Scaffolding for the next sessions** (§9):
   - new datasets `Method` and `ENBW`;
   - a configurable `mask_distance`;
   - a multipole-order label per source.
5. **`ARCHITECTURE.md`** documenting all of it, including the decisions for future work (§10).

---

## 2. Tools and paths on this machine

| What | Where / how |
|---|---|
| Compiler | `gfortran` 13.3, `h5fc` in `/usr/bin` (HDF5 Fortran), FFTW3 |
| Makefile build | `make -j8` from the source root → binary `off_coustics_fortran`, objects in `build_make/` |
| fpm build | `fpm build --flag "-O3 -fopenmp"` (fpm 0.13 in `~/.local/bin`) |
| Python with h5py + numpy | `~/miniconda3/envs/bempp_env/bin/python` |
| h5diff | `~/miniconda3/envs/bempp_env/bin/h5diff` |
| OpenFAST outputs | `/home/hp/simulations/openfast_simulations/off_coustics/OP_output/` |
| Wind-speed curves | `/home/hp/simulations/openfast_simulations/off_coustics/wind_speed_curves_DTU_10MW/` and `.../wind_speed_curves_SAITEC_2MW/`. The code opens them as `../wind_speed_curves_*` relative to the run folder. |
| Verification harness | `review/refactor_tools/refactor_cases.py`. Its copy is in Appendix A. `review/` is git-ignored. |

---

## 3. Phase 0: preparation (no code change)

Use your session scratchpad as `$S`.

### 3.1 Check that the repository is clean

Run `git status` in the repository. It must be clean; `refactor.md` and the `refactor_phase*.patch` / `*_SAITEC_manual.diff` files of this refactor may be untracked. Note the commit hash: this is the **baseline commit** (`1f2f01d` or later). If anything else is modified or untracked, ask the user before going on.

### 3.2 Create your working copy

Your working copy is a detached git worktree outside the repository:
```bash
git -C /home/hp/simulations/openfast_simulations/off_coustics/FORTRAN worktree add $S/wt --detach HEAD
```
All code changes are made in `$S/wt`, never in the user's repository.

### 3.3 Build the baseline

The baseline comes from the baseline commit and is never rebuilt.
```bash
mkdir -p $S/base && git -C <repo> archive HEAD src app Makefile input_files | tar -x -C $S/base
cd $S/base && make -j8
```
Then:
- **`-Wall` warnings:** record them with `make clean && make FFLAGS="-O3 -fopenmp -Wall" 2>&1 | grep -A3 Warning > $S/base_warnings.txt`, then rebuild normally with `make clean && make -j8`.
- **fpm:** also check `fpm build --flag "-O3 -fopenmp"` in a second copy of the baseline. If fpm fails on the baseline, note it and tell the user: an fpm failure that already exists is not yours to fix.

### 3.4 Create and run the reference cases

```bash
PY=~/miniconda3/envs/bempp_env/bin/python
cp <repo>/review/refactor_tools/refactor_cases.py $S/
$PY $S/refactor_cases.py setup $S/base $S/cases          # case configurations from the baseline input_files
$PY $S/refactor_cases.py run $S/base/off_coustics_fortran $S/cases base
```
This takes about 75 s. The cases (Appendix A) are:

| Case | What it checks |
|---|---|
| M1 | DTU 10 MW monopile, ImageMethod, all steps |
| M2 | DTU 10 MW monopile, AnalyticalNormalModes, all steps |
| M3 | `update_existing` on the monopile: polar + line, then sphere only (update), then `N_images = 20` (incompatible, must be rebuilt) |
| L1, L2, L3 | DTU floating, SAITEC and the 2-monopile farm, spectrums only (force code paths) |

**Do not change the cases.** The monopile is the user's validated fast case.

The case configurations are generated **once**, from the baseline `input_files`. Never regenerate them. In particular, after Phase 4 they still lack `mask_distance`, which tests its default.

### 3.5 Patch protocol, used at the end of every phase N

1. **Build and verify** the phase in `$S/wt` (§4).
2. **Stage everything in your worktree** (not in the user's repo) and write two patches:
   ```bash
   cd $S/wt && git add -A
   git diff --cached -M HEAD -- . ':(exclude,glob)**/SAITEC2MWFloating.f90' ':(exclude)off_coustics_fortran' > <repo>/refactor_phaseN.patch
   git diff --cached -M HEAD -- ':(glob)**/SAITEC2MWFloating.f90'        > <repo>/refactor_phaseN_SAITEC.patch
   ```
   The patches are written in the **repository root**, next to `refactor.md`, as the user asked. They are untracked files there: never put them in a patch, and remind the user not to commit them (he may delete them once applied).

   Delete the SAITEC patch if it is empty. Never put build products in a patch (`build_make/`, `build/`, `*.mod`, `*.smod`, `*.o`). The binary `off_coustics_fortran` is tracked in git, which is why it is excluded explicitly above. This protocol was tested: renames are detected, SAITEC goes to its own patch, and both patches apply in order on HEAD.
3. **Check that both patches apply on the user's HEAD**, in a throw-away worktree:
   ```bash
   git -C <repo> worktree add $S/check --detach HEAD
   git -C $S/check apply <repo>/refactor_phaseN.patch
   git -C $S/check apply <repo>/refactor_phaseN_SAITEC.patch     # if it exists
   ```
   Then build `$S/check` with `make` and run the harness with it (`run ... check`, then `compare $S/cases/base $S/cases/check` with the allowed extras of the phase). It must give `ALL IDENTICAL`, exactly as `$S/wt` did. Remove it afterwards with `git -C <repo> worktree remove --force $S/check`.
4. **Report to the user in Spanish:**
   - what was done;
   - the verification table (§4.3);
   - the exact commands for him:
     ```bash
     git apply refactor_phaseN.patch
     git apply refactor_phaseN_SAITEC.patch
     make
     ```
     then he reviews and commits (he may leave the SAITEC file uncommitted). `git apply` needs his working tree to be clean in the files the patch touches; if it fails, ask him to commit or stash his own changes first, never edit the patch to force it.
5. **Stop.** Wait until he says the phase is committed.
6. **Before starting phase N+1:**
   ```bash
   git -C $S/wt reset --hard && git -C $S/wt clean -fd && git -C $S/wt checkout --detach main
   ```
   If the user did not commit the SAITEC patch, re-apply `refactor_phaseN_SAITEC.patch` in `$S/wt`. Then check that `$S/wt` builds and that the harness still gives `ALL IDENTICAL` against `base`, allowing `Method,ENBW` once Phase 4 is in.

---

## 4. Verification after every phase

Run all of this in `$S/wt`.

### 4.1 Builds

- `make clean && make -j8` must succeed.
- `fpm build --flag "-O3 -fopenmp"` must succeed, if it succeeded on the baseline.
- `make clean && make FFLAGS="-O3 -fopenmp -Wall"` must give no warning that is not in `$S/base_warnings.txt`. Compare by message, not by path or line.
- Rebuild normally before running.

### 4.2 Run and compare

```bash
$PY $S/refactor_cases.py run $S/wt/off_coustics_fortran $S/cases phaseN
$PY $S/refactor_cases.py compare $S/cases/base $S/cases/phaseN            # Phases 1-3
$PY $S/refactor_cases.py compare $S/cases/base $S/cases/phaseN Method,ENBW # Phase 4 and 5
```

### 4.3 Acceptance

`compare` must print `ALL IDENTICAL`. That means:
- every HDF5 dataset of the baseline exists with the same dtype, shape and bytes;
- there are no extra datasets except the allowed ones;
- the `.dat` files are identical;
- stdout is identical except the `Total elapsed time` line.

If anything differs, the phase is not finished: find the cause and fix it. Never "accept" a difference. Report a table with one row per case:

| Case | HDF5 | .dat | stdout | Extra datasets |
|---|---|---|---|---|

---

## 5. General rules for moving code (Phases 1–3)

- **Move procedures verbatim.** Same body, same docstring, same local variable names, same `!$OMP` directives.
  - Only `USE` lines, module names, `PUBLIC`/`PRIVATE` lists, file headers and comments that cite file paths may change.
  - The only body changes allowed are those written in Phases 2 and 3.
- **Visibility:** every new module is `PRIVATE` by default, with an explicit `PUBLIC ::` list of exactly what other modules use. Use aligned `&` for long lists (style 2.6).
- **`USE` statements:**
  - Prefer `USE Module, ONLY: ...`, and keep any `ONLY` list that already exists.
  - Order them by line length (style 6.8).
  - Remove every `USE` that becomes unnecessary; check with `-Wall`.
- **File headers:**
  - **New file:** header per style 7.2. The description starts with "Defines". Author block as in the existing files, `Date : 10/2026`.
  - **Moved file:** keep its header; update the description only if its content changed.
  - **Submodule:** header `! SUBMODULE: AcousticSolver_steps`, and so on.
- **Use `git mv`** for every file that keeps its module, so the history follows it.
- **Stdout:** do not change a single printed string (style 7.7). Keep the order of every `write`, `CALL` and file operation.
- **Comments that cite paths:** update comments in `.f90` files and in `input_files/*.nml` that name a source path or a moved routine, for example the header of `ComputationSetupTurbine.nml`. This is the only edit allowed in `input_files/` in Phases 1–3, and it must touch comments only.
- **Build files:** update the `Makefile` so it lists every new file in dependency order, with submodules after their parent (see §6.4). fpm discovers `src/**` by itself; do not edit `fpm.toml` unless the build needs it.

---

## 6. Phase 1: new layout, pure relocation

No procedure signature changes in this phase.

### 6.1 Target tree

```
src/
├── Kinds.f90                        (unchanged)
├── types/
│   ├── SourceTypes.f90              MODULE SourceTypes
│   └── ConfigTypes.f90              MODULE ConfigTypes
├── utils/
│   ├── FFTW3.f90                    (unchanged)
│   ├── StringUtils.f90              MODULE StringUtils
│   ├── FileUtils.f90                MODULE FileUtils
│   ├── HDF5Utils.f90                MODULE HDF5Utils
│   ├── CSVReader.f90                MODULE CSVReader
│   ├── Sorting.f90                  MODULE Sorting
│   ├── Timing.f90                   MODULE Timing
│   ├── SignalProcessing.f90         MODULE SignalProcessing
│   └── SpecialFunctions.f90         MODULE SpecialFunctions
├── input/
│   ├── Args_parser.f90              MODULE Args_parser   (git mv from src/)
│   ├── InputEmitters.f90            MODULE InputEmitters
│   ├── InputSolver.f90              MODULE InputSolver
│   ├── InputPostProcessing.f90      MODULE InputPostProcessing
│   └── InputChecks.f90              MODULE InputChecks
├── emitters/
│   ├── WindTurbine.f90              (git mv from src/core/)
│   ├── TurbineFactory.f90           MODULE TurbineFactory
│   ├── WindFarm.f90                 (git mv from src/core/)
│   ├── OpenFASTReader.f90           MODULE OpenFASTReader
│   ├── StructureModel.f90           MODULE StructureModel
│   ├── DrivetrainExcitation.f90     MODULE DrivetrainExcitation
│   └── TurbineTypes/                (git mv the three files from src/core/TurbineTypes/)
└── solvers/
    ├── AcousticSolver.f90           (git mv from src/core/)
    ├── AcousticSolver_steps.f90     SUBMODULE (AcousticSolver) AcousticSolver_steps
    ├── AcousticSolver_output.f90    SUBMODULE (AcousticSolver) AcousticSolver_output
    └── methods/
        ├── ImageMethod.f90          (git mv from src/core/AcousticMethods/)
        └── AnalyticalNormalModes.f90 (git mv from src/core/AcousticMethods/)
app/main.f90                         (only USE lines change)
```

**Delete** `src/Types.f90`, `src/utils/IOUtils.f90`, `src/utils/MathUtils.f90` and the empty `src/core/` tree. **No compatibility or "umbrella" modules:** every `USE` in every file, SAITEC included, is updated to the new modules.

### 6.2 Where every entity goes

**`src/Types.f90`:**

| Entity | Goes to |
|---|---|
| `GeometricData_t`, `AcousticSource_t` | `types/SourceTypes.f90` |
| `Config_WindTurbine_t`, `Config_WindFarm_t`, `Config_Environment_t`, `Config_Solver_t`, `Config_Method_t` and its extensions `Config_ImageMethod_t` and `Config_AnalyticalNormalModes_t`, `Config_DoRun_t`, the nine step configurations (`Config_Spectrums_t` … `Config_Sphere_t`), `Config_Computation_t`, `Config_t` | `types/ConfigTypes.f90` (same order as now) |
| `ImageSystem_t` | `solvers/methods/ImageMethod.f90`, defined before `ImageMethod_t` (it is internal to the Image Method) |

If `ConfigTypes` needs a type from `SourceTypes` or the opposite, put the `USE` in the direction that avoids a cycle. Today none of them needs the other.

**`src/utils/IOUtils.f90`:**

| Procedures | Goes to |
|---|---|
| `to_upper`, `split_tokens`, `is_data_line` | `utils/StringUtils.f90` |
| `read_full_line`, `swap_extension`, `ensure_parent_dir`, `move_text_file`, plus `delete_file` from `AcousticSolver.f90` | `utils/FileUtils.f90` |
| `save_variable`, the generic interface `save_to_hdf5` and all 15 `save_hdf5_*`, `read_hdf5_real_1d` | `utils/HDF5Utils.f90` (the only module that does `USE HDF5`) |
| `read_curve`, `parse_csv_row` | `utils/CSVReader.f90` |
| `get_SDsum_variables`, `set_error_status`, `read_input_SD`, `parse_channels_auto`, `parse_channels_binary`, `parse_channels_ascii`, `parse_member_nodes`, `extract_number_nodes` | `emitters/OpenFASTReader.f90` |

**`src/utils/MathUtils.f90`:**

| Procedures | Goes to |
|---|---|
| the C interface `qsort`, `compare_real`, `sort_real_c` | `utils/Sorting.f90` |
| `format_elapsed` | `utils/Timing.f90` |
| `compute_rfft` | `utils/SignalProcessing.f90` |
| `compute_hankel_complex`, `bessel_k0`, `bessel_k1` | `utils/SpecialFunctions.f90` |
| `remove_duplicate_nodes`, `divide_span`, `integral_elipse_z`, `alpha_hankel` | `emitters/StructureModel.f90` (geometry and hydrodynamic corrections of the structural members) |
| `generate_timeseries_banded_sines`, `filter_non_useful_freqs` | `emitters/DrivetrainExcitation.f90` |

**Other files:**

| From | Procedures | Goes to |
|---|---|---|
| `WindTurbine.f90` | `drivetrain10MW_excitation_spectrum` | `emitters/DrivetrainExcitation.f90` |
| `SAITEC2MWFloating.f90` | `drivetrain2MW_excitation_spectrum`, `compute_freqs2MW`, `sorted_amps_by_freq`, `append_family`, `read_bar_csv` (module procedures, not type-bound) | `emitters/DrivetrainExcitation.f90` |
| `WindFarm.f90` | `new_turbine` | `emitters/TurbineFactory.f90`; `main.f90` and `WindFarm.f90` then `USE TurbineFactory` |
| `Args_parser.f90` | `read_turbine_config`, `read_farm_config` | `input/InputEmitters.f90` |
| `Args_parser.f90` | `read_solver_config`, `read_environment_config`, `read_ImageMethod_config`, `read_AnalyticalNormalModes_config` | `input/InputSolver.f90` |
| `Args_parser.f90` | `read_computation_config`, `read_do_config`, and the nine `read_*_config` of the steps (`spectrums` … `sphere`) | `input/InputPostProcessing.f90` |
| `Args_parser.f90` | `check_time_window`, `check_fits`, `check_all_or_none`, `check_count` | `input/InputChecks.f90` |
| `Args_parser.f90` | `read_input_files` (the only public entry) | stays in `input/Args_parser.f90`, which USEs the four modules above |

If a helper is also used from a module that would create a dependency cycle, stop and ask.

`AcousticSolver` currently re-exports `save_variable` in its `PUBLIC` list. After the split it no longer does: every module that writes HDF5 does `USE HDF5Utils, ONLY: save_variable` itself.

### 6.3 Splitting `AcousticSolver.f90` with Fortran 2008 submodules

**`solvers/AcousticSolver.f90`** (the parent module) keeps:
- `MASK_DISTANCE`;
- the type `AcousticSolver_t`, with all its bindings unchanged;
- the `ABSTRACT INTERFACE` block;
- the core procedures `check_observers_distances`, `compute_pressure_masked` and `pressure_norm`;
- **one `INTERFACE` block of separate module procedures** (`MODULE SUBROUTINE` / `MODULE FUNCTION` interface bodies, each with its docstring) for every procedure implemented in a submodule that is type-bound or public.

**`solvers/AcousticSolver_steps.f90`** (`SUBMODULE (AcousticSolver) AcousticSolver_steps`) takes `run_spectrums`, `run_polar`, `run_cylinder`, `run_decay`, `run_line`, `run_sliceXY`, `run_sliceXZ`, `run_sliceVertical`, `dummy_run_sphere`, `run` and `print_step_header`. `print_step_header` is private to the submodule and needs no interface in the parent.

**`solvers/AcousticSolver_output.f90`** (`SUBMODULE (AcousticSolver) AcousticSolver_output`) takes:
- `resolve_save_path`, `check_update`, `update_parameters`;
- `save_emitter_params` and `write_emitter_txt` (they disappear in Phase 2);
- the private helpers `compare_dat_headers`, `read_header_lines` and `write_previous_block`.

In the submodules, implement the procedures with `MODULE SUBROUTINE name(...)` or `MODULE PROCEDURE name`. Use one form consistently, preferably `MODULE SUBROUTINE` with the full argument declarations, so each file reads on its own.

### 6.4 Makefile

- Rewrite `OBJS` and the per-file rules for the new paths in dependency order: Kinds → types → utils → input / emitters → solvers → app.
- A submodule object depends on its parent's object. gfortran writes `acousticsolver@acousticsolver_steps.smod` into `$(BUILD)` through `-J`.
- `HDF5Utils.f90` is compiled with `$(H5FC)` and `$(HDF5_INCLUDE)`, as `IOUtils.f90` was. No other file needs them, unless it does `USE HDF5`.
- Keep the variables, comments style, `run` and `clean` targets of the current Makefile.

### 6.5 Acceptance

§4 with no allowed extras.

---

## 7. Phase 2: abstract `Emitter_t`

### 7.1 New module `emitters/Emitter.f90`

`MODULE Emitter` uses `Kinds`, `ConfigTypes` and `SourceTypes`, and contains:
```fortran
TYPE, ABSTRACT :: Emitter_t
    ! -------------------------------------------------------------------------
    ! Parent of everything that produces acoustic sources (a wind turbine, a
    ! wind farm, and in the future manually defined sources). The solvers only
    ! see this interface: the sources it fills in and the blocks it writes in
    ! the HDF5 and .dat outputs.
    ! -------------------------------------------------------------------------

    CONTAINS

    ! --- Deferred (abstract) procedures --- !
    PROCEDURE(to_source_i)  , DEFERRED :: to_source      ! Fills the AcousticSource_t handed to the solver
    PROCEDURE(save_params_i), DEFERRED :: save_params    ! Emitter parameters -> HDF5 output
    PROCEDURE(write_txt_i)  , DEFERRED :: write_txt      ! Emitter block of the .dat summary

END TYPE Emitter_t
```

The three interfaces take the arguments, intents and attributes of today's procedures, with `CLASS(Emitter_t)` as passed object:
- `to_source(self, config, sources)`:
  - `CLASS(Emitter_t), INTENT(INOUT) :: self`;
  - `TYPE(Config_t), INTENT(IN) :: config`;
  - `TYPE(AcousticSource_t), INTENT(OUT), ALLOCATABLE :: sources`;
  - copy the exact declarations of `turbine_to_source` and `farm_to_source_WindFarm`.
- `save_params(self, save_path)`: `self` `INTENT(IN)`, `CHARACTER(len=*), INTENT(IN) :: save_path`.
- `write_txt(self, unit)`: `self` `INTENT(IN)`, `INTEGER(I32), INTENT(IN) :: unit`.

Overriding bindings require the same dummy names (except the passed object). Write thin wrapper procedures where the existing names differ.

### 7.2 Changes to the emitters and the solvers

**`WindTurbine_t`** becomes `TYPE, ABSTRACT, EXTENDS(Emitter_t) :: WindTurbine_t`:
- `to_source` is bound to the existing `turbine_to_source`, with no body change.
- `save_params` is bound to a new procedure `save_params_turbine(self, save_path)`. Its body is exactly what `save_emitter_params` did for a turbine, in the same order: `CALL save_turbine_params(self, save_path)` then `CALL self % save_specific_params(save_path)`.
- `write_txt` is bound to a wrapper doing exactly what `write_emitter_txt` did for a turbine: `CALL write_turbine_txt(self, unit)`.
- Remove the old `turbine_to_source` binding name; callers use `to_source`.

**`WindFarm_t`** becomes `TYPE, EXTENDS(Emitter_t) :: WindFarm_t`:
- `to_source` is bound to `farm_to_source_WindFarm`. Remove the `farm_to_source` binding name.
- `save_params` and `write_txt` keep their existing implementations.

**Solvers:**
- `save_parameters_i` and `write_txt_summary_i` become `(self, emitter)` with `CLASS(Emitter_t), INTENT(IN) :: emitter` (mandatory). Change the implementations in `ImageMethod` and `AnalyticalNormalModes`, and `check_update(self, emitter)`.
- Each call to `save_emitter_params(path, turbine, farm)` becomes `CALL emitter % save_params(path)`, and each call to `write_emitter_txt(u, turbine, farm)` becomes `CALL emitter % write_txt(u)`, at the same place in the call order. Delete `save_emitter_params` and `write_emitter_txt`.
- `AcousticSolver`, `ImageMethod` and `AnalyticalNormalModes` no longer `USE WindTurbine` or `USE WindFarm`; they `USE Emitter, ONLY: Emitter_t`.

**`main.f90`:**
- Keep the construction of the turbine or the farm exactly as now: the same prints, calls and order.
- Declare `CLASS(Emitter_t), ALLOCATABLE :: active_emitter`.
- Call `active_turbine % to_source(...)` and `active_farm % to_source(...)` where `turbine_to_source` and `farm_to_source` were called.
- Hand the object to the solver with `MOVE_ALLOC` into `active_emitter`, then call `active_solver % save_parameters(active_emitter)` at the place where the old call was.
- Deallocate the emitter where `active_turbine` and `active_farm` are deallocated today, before `run()`.

### 7.3 Acceptance

§4 with no allowed extras.

---

## 8. Phase 3: factory, turbine template, water-column queries

### 8.1 `solvers/SolverFactory.f90`

`MODULE SolverFactory` contains `SUBROUTINE new_solver(method, solver)`:
- `CHARACTER(len=*), INTENT(IN) :: method`;
- `CLASS(AcousticSolver_t), ALLOCATABLE, INTENT(OUT) :: solver`.

Move into it the `SELECT CASE` of `main.f90`, with the same `CASE` strings. Add `CASE DEFAULT` with `error stop "new_solver(): Unknown method " // trim(method)`. The parser already rejects unknown methods, so it never fires.

`main.f90` calls `CALL new_solver(trim(config % solver % method), active_solver)` at the same place, keeps the `>> Launching acoustic propagation solver` line where it is, and no longer `USE`s `ImageMethod` or `AnalyticalNormalModes`.

### 8.2 A single `filter_frequencies` for every turbine type

1. **Check first.** Diff the three `filter_frequencies_*` bodies (DTU10MWMonopile, DTU10MWFloating, SAITEC2MWFloating). They must be identical except for:
   - the call that builds the excitation spectrum;
   - comments and the layout of declarations.

   If they differ in anything else, stop and report to the user.
2. **The common body** becomes a non-deferred `filter_frequencies` of `WindTurbine_t`, in `WindTurbine.f90`. `path_rpm`, `WindSpeed` and `Freqs` are already components of `WindTurbine_t`. Remove `filter_frequencies_i` from the deferred list.
3. **A new deferred procedure** `excitation_spectrum(self, rpm, freqs_amp, keys)` of `WindTurbine_t`, with the declarations the bodies use today:
   - `rpm` `REAL(WP)`;
   - `freqs_amp` `REAL(WP), ALLOCATABLE (:,:)`, `INTENT(OUT)`;
   - `keys` `CHARACTER(len=20), ALLOCATABLE (:)`, `INTENT(OUT)`;
   - copy their exact kinds and lengths.

   Each type implements it with its exact current call:
   - **DTU10MWMonopile** and **DTU10MWFloating:** `CALL drivetrain10MW_excitation_spectrum(rpm=rpm, freqs_amp=freqs_amp, keys=keys)`.
   - **SAITEC2MWFloating:** `CALL drivetrain2MW_excitation_spectrum(rpm=rpm, path_barRB=self % path_barRB, path_barGB=self % path_barGB, fmax=maxval(self % Freqs), Normalised=.false., freqs_amp=freqs_amp, keys=keys)`.
4. **Remove** the three `filter_frequencies_*` procedures and their bindings. Remove the `USE` lines they no longer need.

### 8.3 Water-column queries

The `run_*` steps stay generic. A future method with variable bathymetry overrides these queries; it does not override the steps.

**Add to `AcousticSolver_t`** three non-deferred, overridable type-bound procedures, implemented in the parent module (not in a submodule):

| Binding → implementation | Declarations | Default result |
|---|---|---|
| `z_surface => default_z_surface` | `FUNCTION default_z_surface(self, x, y) RESULT(z)`; `self` `CLASS(AcousticSolver_t), INTENT(IN)`; `x`, `y` `REAL(WP), INTENT(IN)` `[m]`; `z` `REAL(WP)` `[m]` | `0.0_WP` |
| `z_seabed => default_z_seabed` | same arguments as above | `-self % env % depth_base` |
| `in_water => default_in_water` | `FUNCTION default_in_water(self, observers) RESULT(valid)`; `observers(:,:)` `[m]` `(Nobs, 3)`; `valid` `LOGICAL, ALLOCATABLE (:)` | all `.true.` |

Docstrings must explain that a method overrides them to publish its own boundaries (surface, bathymetry), and that the `run_*` steps never read `depth_base` directly.

**In every `run_*`** (the submodule `AcousticSolver_steps`), replace the reads of `self % env % depth_base`:
1. At the step's reference point, set `z_top = self % z_surface(xr, yr)`, `z_bot = self % z_seabed(xr, yr)` and `Depth = z_top - z_bot`.
2. Replace the literal surface `0.0_WP` with `z_top`, and the seabed `-Depth` with `z_bot`.
3. Replace mid-depth `-Depth / 2.0_WP` with `0.5_WP * (z_top + z_bot)`.

With the default queries every expression gives the **same bits** (0 − x and scaling by 2 are exact), and the comparison must prove it. Reference points:

| Step | Reference point (xr, yr) | What uses the queries |
|---|---|---|
| `run_spectrums` | each default observer (10/250/500, 0) | default z = mid-depth at that point |
| `run_polar` | `center_polar` (after its default is resolved) | default `z_polar` |
| `run_cylinder` | `center_cyl` | z levels from `z_bot` to `z_top` |
| `run_decay` | `geometry % origin` | default `z_decay` |
| `run_line` | the default end point (100, 0) | default `p2_line(3)` |
| `run_sliceXY` | `geometry % center` | default `z_xy`, and the "z must be <= 0" check, which becomes `> z_top` (message unchanged) |
| `run_sliceXZ` | the middle of `xlim_xz`, at `y_xz` (after their defaults are resolved) | default `zlim_xz` |
| `run_sliceVertical` | `center_sv` | default and clipping of `zlim_sv` |

**`compute_pressure_masked`:**
- The mask becomes `valid = check_observers_distances(observers) .and. self % in_water(observers)`.
- Count the observers masked by distance and by `in_water` separately.
- The existing WARNING line (distance) keeps its exact text and value.
- Print an extra WARNING only when `in_water` masks something, which never happens with the default: `' WARNING [' // label // '] ', n, ' of ', n_obs, ' observers outside the water column: masked (NaN)'`.

Do not override the queries in `ImageMethod` or `AnalyticalNormalModes` in this refactor; see §10, "Planned extensions".

### 8.4 Acceptance

§4 with no allowed extras.

---

## 9. Phase 4: scaffolding for the next sessions

This is the only phase with allowed output differences.

### 9.1 Dataset `Method`

Each solver's `save_parameters` also writes the solver name returned by `get_name()` (stored in `self % solver_name`) as a character scalar:
```fortran
CALL save_variable(self % save_path, "Method", trim(self % solver_name))
```
If the `CLASS(*)` dispatcher refuses an expression, pass a local `CHARACTER(len=:), ALLOCATABLE` copy of the name.
Write it in a parent procedure `save_shared_parameters(self, emitter)` of `AcousticSolver_t`, called by both solvers right after their own datasets. Every output gets it.

### 9.2 Dataset `ENBW`

- **Constant:** add to `SignalProcessing` `REAL(WP), PARAMETER :: HANN_ENBW = 1.5_WP`, documented as the equivalent noise bandwidth, in bins, of the Hann window of `compute_rfft` (N·Σw²/(Σw)², which tends to 1.5 for large N).
- **Default function:** add to `Emitter_t` a non-deferred function `get_enbw(self) RESULT(enbw)` that returns `-1.0_WP`, meaning "forces not obtained from a windowed FFT".
- **Overrides:** `WindTurbine_t` and `WindFarm_t` override it and return `HANN_ENBW`.
- **Writing it:** `save_shared_parameters` writes `"ENBW"` only when `emitter % get_enbw() > 0.0_WP`.
- **Not in the `.dat`:** neither `Method` nor `ENBW` is printed in the `.dat`. The `.dat` must stay identical, and `check_update` compares its header.

### 9.3 Configurable `mask_distance`

- **Namelist:** add `mask_distance` to `&solver_config` with this line in `input_files/Solver.nml`, after `update_existing`:
  ```
      mask_distance = 5.0                 ! [m] Observers closer than this to a source node get NaN pressure (0 = no masking)
  ```
- **Config and parser:** add a component `mask_distance` to `Config_Solver_t`. In `read_solver_config`, default `5.0_WP`, plus an `error stop` with a clear message if it is `< 0`.
- **Solver:** add a component `mask_distance` to `AcousticSolver_t`, set in both `init` (like `update_existing`).
- **Uses:** replace the uses of the parameter `MASK_DISTANCE` by `self % mask_distance`, then remove the parameter. The WARNING keeps its format and prints `5.0` as now.

### 9.4 Multipole order of every source

- **Constants:** in `SourceTypes`, `INTEGER(I32), PARAMETER :: SOURCE_MONOPOLE = 0_I32, SOURCE_DIPOLE = 1_I32`. They are a label, so that higher orders can be added later without new arrays.
- **Component:** add to `AcousticSource_t` `INTEGER(I32), ALLOCATABLE :: order(:)` `! [-] Multipole order of every source (SOURCE_MONOPOLE, SOURCE_DIPOLE, ...)`.
- **Setting it:** `turbine_to_source` and `farm_to_source_WindFarm` allocate it to the number of sources and set it to `SOURCE_DIPOLE`.
- **Moving it:** both solvers' `init` move it with `MOVE_ALLOC` like `xyz`, `Force` and `Freqs`.
- **Checking it:** right after the move, both `init` stop if any source is not a dipole:
  ```fortran
  if (any(self % sources % order /= SOURCE_DIPOLE)) &
      error stop "<Solver> % init(): only dipole sources (order 1) are implemented"
  ```
  Use the real module name instead of `<Solver>`.
- **Not saved:** it is not written to the HDF5 yet.

### 9.5 Acceptance

`compare ... Method,ENBW` must print `ALL IDENTICAL`. In addition, check with h5py on M1 and on L3 (the farm) that:
- `Method` is `"Images Method"` or `"Analytical Normal Modes"`;
- `ENBW` is `1.5`.

The patch includes the `Solver.nml` line.

---

## 10. Phase 5: `ARCHITECTURE.md` and documentation

### 10.1 `FORTRAN/ARCHITECTURE.md`

Write `FORTRAN/ARCHITECTURE.md` in English with these sections:

1. **Pipeline:** `.nml` → `input/` → emitter (`Emitter_t`) → `AcousticSource_t` → solver (`AcousticSolver_t` + method) → `run_*` steps → `.hdf5` + `.dat` → `plotting/`. Include a small ASCII diagram.
2. **Directory tree:** one line per module saying what it contains and what must not go there.
3. **Dependency layers:** Kinds → types → utils → input / emitters → solvers → app. The rule: a module never `USE`s a higher layer. List which modules each layer contains.
4. **Contracts:**
   - `Emitter_t`: the three deferred procedures and `get_enbw`;
   - `AcousticSolver_t`: each deferred procedure, the water-column queries, `compute_pressure_masked` and `save_shared_parameters`;
   - the fields of `AcousticSource_t` (including `order`) and of `GeometricData_t`, and what each `run_*` takes from them.
5. **"How to add…" recipes:** for each case, the exact list of files to touch:
   - a turbine type: `TurbineTypes/` plus one `CASE` in `TurbineFactory`, the Makefile, and `excitation_spectrum`;
   - an emitter (manual sources);
   - an acoustic method: `solvers/methods/`, `Config_X_t` in `ConfigTypes`, `read_X_config` in `InputSolver` plus the `CASE` in `read_input_files`, the `CASE` in `SolverFactory`, the Makefile, and the `Method` dataset is automatic;
   - a post-processing step: `AcousticSolver_steps`, the configuration type, its reader, its block in `update_parameters`, and the plotting.
6. **Output interface:**
   - the HDF5 datasets written by each step and by `save_parameters` (including `Method` and `ENBW`);
   - the format of the `.dat`;
   - the rules of `update_existing`: same-case check (`Freqs` + `.dat` header), `.dat.prev` merge, `(previous run)` tags;
   - a warning that dataset names, the `.dat` format and namelist names are an interface used by `plotting/`.
7. **Planned extensions** (decided on 2026-10-02, not implemented):
   - **Manual sources:**
     - `simulation_type = 'sources'` and `ManualSources_t EXTENDS(Emitter_t)` in `emitters/ManualSources.f90`;
     - a namelist for a few sources and an input HDF5 file for many (`xyz`, `order`, `Freqs`, `Force`);
     - forces given in the time domain (they go through `compute_rfft`, so `get_enbw` returns `HANN_ENBW`) or in the frequency domain (no window, so `get_enbw` keeps `-1` and no `ENBW` is written);
     - no radiation-impedance correction;
     - one multipole-order label per source; each method implements one kernel per order and drops the dipole-only check when it supports monopoles;
     - `mask_distance` set by the user;
     - the plotting will read `ENBW` (Hann correction only if present) and `Method`.
   - **Acoustic methods:**
     - `Config_Environment_t` will grow optional arrays (c(z), ρ(z), attenuation, bathymetry H(x, y), sediment) and keep the base scalars, so simple methods keep working; each method checks in `init` that the data it needs is present;
     - the boundaries belong to each method (constants in IM/ANM, arrays in complex methods), which publishes them by overriding `z_surface`, `z_seabed` and `in_water`;
     - overriding a whole `run_*` is the exception (precedent: IM's `run_sphere`).
   - **Known inconsistency, pending a decision by the user:** the default queries use the surface at 0 and the seabed at −`depth_base`, while IM and ANM use `Upper_HBC` / `Lower_HBC`. When IM and ANM override the queries with their HBC, results change only for `UpperHBC ≠ 0`.
   - **Sphere:** it stays IM-only (empty version in the parent plus the IM override; an `AcousticSolver` cannot `USE ImageMethod`, because that would be a cycle). A shared free-field kernel would make it generic for every method; that changes outputs and is a future decision.
8. **Verification protocol:** style guide §8, plus the harness of `review/refactor_tools/` and its acceptance rules.

### 10.2 Other documentation

- **Style guide:** update `Programming_style_guidelines.txt` only where it cites paths that no longer exist. For example, `src/utils/MathUtils.f90` becomes a file that exists, such as `src/utils/SignalProcessing.f90`. Change nothing else in it.
- **Acceptance:** §4 with `Method,ENBW` allowed. Nothing in the code should change in this phase; re-run to confirm. The patch contains only documentation.

---

## 11. Final report (after Phase 5)

Give the user, in Spanish:
- one line per phase with its patches (in the repository root);
- the final module tree;
- any decision you had to ask him about;
- the list of "Planned extensions" (§10.1, item 7).

Remove your worktree:
```bash
git -C <repo> worktree remove --force $S/wt
git -C <repo> worktree prune
```

---

## Appendix A: `review/refactor_tools/refactor_cases.py`

This is a copy of the validated harness: two runs of the baseline binary gave `ALL IDENTICAL`, in about 73 s per run. If the file in `review/` is missing, recreate it from this copy.

```python
"""
Reference cases for the bitwise verification of the refactor (see FORTRAN/refactor.md, section 4).

    python refactor_cases.py setup   <repo_dir> <cases_root>             # build the case configurations
    python refactor_cases.py run     <binary>   <cases_root> <tag>        # run every case into <cases_root>/<tag>/
    python refactor_cases.py compare <cases_root>/<tagA> <cases_root>/<tagB> [extra1,extra2,...]

Every case runs with OMP_NUM_THREADS=1 in its own folder, with save_path = 'out/test', so the outputs of two
binaries can be compared byte by byte. "compare" checks: every HDF5 dataset of A exists in B with the same dtype,
shape and bytes; B has no extra dataset except the allowed ones; the .dat files are identical; stdout is identical
except the "Total elapsed time" line.
"""
import os, re, sys, shutil, subprocess
import numpy as np

OP_OUTPUT = "/home/hp/simulations/openfast_simulations/off_coustics/OP_output/"
CURVES    = ["/home/hp/simulations/openfast_simulations/off_coustics/wind_speed_curves_DTU_10MW",
             "/home/hp/simulations/openfast_simulations/off_coustics/wind_speed_curves_SAITEC_2MW"]

SMALL_GRIDS = [(r"^(\s*n_polar = )-?\d+", r"\g<1>36"), (r"^(\s*n_cyl = )-?\d+", r"\g<1>36"),
               (r"^(\s*nz_cyl = )-?\d+", r"\g<1>10"), (r"^(\s*n_decay = )-?\d+", r"\g<1>50"),
               (r"^(\s*n_line = )-?\d+", r"\g<1>40"), (r"^(\s*nx_xy = )-?\d+", r"\g<1>11"),
               (r"^(\s*ny_xy = )-?\d+", r"\g<1>11"), (r"^(\s*nx_xz = )-?\d+", r"\g<1>11"),
               (r"^(\s*nz_xz = )-?\d+", r"\g<1>6"), (r"^(\s*nu_sv = )-?\d+", r"\g<1>20"),
               (r"^(\s*nz_sv = )-?\d+", r"\g<1>10"), (r"^(\s*n_theta_sphere = )-?\d+", r"\g<1>36"),
               (r"^(\s*nz_sphere = )-?\d+", r"\g<1>10")]

def solver(method, sim="turbine"):
    return [(r"simulation_type = '[^']*'", f"simulation_type = '{sim}'"), (r"method = '[^']*'", f"method = '{method}'"),
            (r"save_path = '[^']*'", "save_path = 'out/test'"), (r"extra_label = '[^']*'", "extra_label = 'None'"),
            (r"debug = \.\w+\.", "debug = .false."), (r"update_existing = \.\w+\.", "update_existing = .true.")]

def turbine(case_type, rootname, ws, nm, nn):
    return [(r"case_type = '[^']*'", f"case_type = '{case_type}'"), (r"rootname = '[^']*'", f"rootname = '{rootname}'"),
            (r"OP_output_dir = '[^']*'", f"OP_output_dir = '{OP_OUTPUT}'"), (r"^(\s*WindSpeed = )[-0-9.]+", rf"\g<1>{ws}"),
            (r"^(\s*Nmembers = )\d+", rf"\g<1>{nm}"), (r"^(\s*Nnodes = )\d+", rf"\g<1>{nn}"),
            (r"^(\s*WindDir = )[-0-9.]+", r"\g<1>0.0"), (r"^(\s*AxisPos = )[-0-9., ]+!", r"\g<1>0.0, 0.0     !"),
            (r"^(\s*BariPos = )[-0-9., ]+!", r"\g<1>0.0, 0.0     !"), (r"^(\s*From = )[-0-9.]+", r"\g<1>0.9"),
            (r"^(\s*Upto = )[-0-9.]+", r"\g<1>1.0"), (r"^(\s*skip = )\d+", r"\g<1>1"), (r"^(\s*skipf = )\d+", r"\g<1>1"),
            (r"filter_freqs = \.\w+\.", "filter_freqs = .true.")]

def steps(do_all, on=()):
    rules = [(r"^(\s*do_all\s*=\s*)\.\w+\.", rf"\1.{'true' if do_all else 'false'}.")]
    rules += [(r"^(\s*do_(?!all)\w+\s*=\s*)\.\w+\.", r"\1.false.")]
    rules += [(rf"^(\s*do_{s}\s*=\s*)\.false\.", r"\1.true.") for s in on]
    return rules

CASES = {   # name: (edits per input file, list of runs; a run = (label, edits applied before it))
    # Full checks (every post-processing step): DTU 10 MW monopile, the fastest validated case, with both methods
    "M1_IM_monopile": ({"Solver.nml": solver("ImageMethod"),
                        "WindTurbine.nml": turbine("DTU10MWMonopile", "DTU_DeltaWind_mn_ws11.4", 11.4, 8, 5),
                        "ComputationSetupTurbine.nml": steps(True) + SMALL_GRIDS}, [("run", {})]),
    "M2_ANM_monopile": ({"Solver.nml": solver("AnalyticalNormalModes"),
                         "WindTurbine.nml": turbine("DTU10MWMonopile", "DTU_DeltaWind_mn_ws11.4", 11.4, 8, 5),
                         "ComputationSetupTurbine.nml": steps(True) + SMALL_GRIDS}, [("run", {})]),
    # update_existing: polar + line, then sphere only (same case: update), then N_images changed (rebuild)
    "M3_update_monopile": ({"Solver.nml": solver("ImageMethod"),
                            "WindTurbine.nml": turbine("DTU10MWMonopile", "DTU_DeltaWind_mn_ws11.4", 11.4, 8, 5),
                            "ComputationSetupTurbine.nml": SMALL_GRIDS},
                           [("u1_polar_line", {"ComputationSetupTurbine.nml": steps(False, ("polar", "line"))}),
                            ("u2_sphere",     {"ComputationSetupTurbine.nml": steps(False, ("sphere",))}),
                            ("u3_new_images", {"ComputationSetupTurbine.nml": steps(False, ("line",)),
                                               "ImagesMethod.nml": [(r"^(\s*N_images = )\d+", r"\g<1>20")]})]),
    # Light checks (spectrums only): the other emitters, whose force code the refactor moves
    "L1_floating": ({"Solver.nml": solver("ImageMethod"),
                     "WindTurbine.nml": turbine("DTU10MWFloating", "DTU_DeltaWind_fl_ws11.4", 11.4, 9, 9),
                     "ComputationSetupTurbine.nml": steps(False, ("spectrums",))}, [("run", {})]),
    "L2_SAITEC": ({"Solver.nml": solver("ImageMethod"),
                   "WindTurbine.nml": turbine("SAITEC2MWFloating", "SENVION_2MW_fl_ws15.0", 15.0, 7, 9),
                   "ComputationSetupTurbine.nml": steps(False, ("spectrums",))}, [("run", {})]),
    "L3_farm": ({"Solver.nml": solver("ImageMethod", sim="farm"),
                 "WindFarm.nml": [(r"OP_output_dir = '[^']*'", f"OP_output_dir = '{OP_OUTPUT}'")],
                 "ComputationSetupFarm.nml": steps(False, ("spectrums",))}, [("run", {})]),
}

def edit(path, rules):
    text = open(path).read()
    for pattern, repl in rules:
        text = re.sub(pattern, repl, text, flags=re.M)
    open(path, "w").write(text)

def setup(repo, root):
    tpl = os.path.join(root, "template")
    shutil.rmtree(tpl, ignore_errors=True)
    for name, (edits, _) in CASES.items():
        shutil.copytree(os.path.join(repo, "input_files"), os.path.join(tpl, name, "input_files"))
        for fname, rules in edits.items():
            edit(os.path.join(tpl, name, "input_files", fname), rules)
    print(f"cases written in {tpl}")

def run(binary, root, tag):
    dst = os.path.join(root, tag)
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(os.path.join(root, "template"), dst)
    for c in CURVES:                                    # the code opens ../wind_speed_curves_* from the run folder
        os.symlink(c, os.path.join(dst, os.path.basename(c)))
    env = dict(os.environ, OMP_NUM_THREADS="1")
    for name, (_, runs) in CASES.items():
        cdir = os.path.join(dst, name)
        for label, edits in runs:
            for fname, rules in edits.items():
                edit(os.path.join(cdir, "input_files", fname), rules)
            with open(os.path.join(cdir, f"stdout_{label}.txt"), "w") as out:
                rc = subprocess.run([os.path.abspath(binary)], cwd=cdir, stdout=out, stderr=subprocess.STDOUT, env=env).returncode
            for ext in ("hdf5", "dat"):                      # keep a copy of the outputs of every run
                src = os.path.join(cdir, "out", f"test.{ext}")
                if os.path.exists(src):
                    shutil.copy(src, os.path.join(cdir, f"{label}.{ext}"))
            print(f"{tag}/{name}/{label}: exit {rc}")

def compare(a, b, extra=""):
    import h5py
    allowed, ok = set(filter(None, extra.split(","))), True
    for name in sorted(os.listdir(a)):
        ca, cb = os.path.join(a, name), os.path.join(b, name)
        if not os.path.isdir(ca) or name.startswith("wind_speed"):
            continue
        for f in sorted(os.listdir(ca)):
            pa, pb = os.path.join(ca, f), os.path.join(cb, f)
            if f.endswith(".hdf5"):
                da, db = {}, {}
                h5py.File(pa, "r").visititems(lambda k, o: da.__setitem__(k, o[()]) if isinstance(o, h5py.Dataset) else None)
                h5py.File(pb, "r").visititems(lambda k, o: db.__setitem__(k, o[()]) if isinstance(o, h5py.Dataset) else None)
                diff = [k for k in da if k not in db or np.asarray(da[k]).dtype != np.asarray(db[k]).dtype
                        or np.asarray(da[k]).shape != np.asarray(db[k]).shape
                        or np.asarray(da[k]).tobytes() != np.asarray(db[k]).tobytes()]
                extra_b = sorted(set(db) - set(da) - allowed)
                status = "IDENTICAL" if not diff and not extra_b else f"DIFFERENT {diff} extra {extra_b}"
            elif f.endswith(".dat"):
                status = "IDENTICAL" if open(pa).read() == open(pb).read() else "DIFFERENT"
            elif f.startswith("stdout_"):
                clean = lambda p: [l for l in open(p) if "Total elapsed time" not in l]
                status = "IDENTICAL" if clean(pa) == clean(pb) else "DIFFERENT"
            else:
                continue
            ok &= status == "IDENTICAL"
            print(f"{name}/{f}: {status}")
    print("ALL IDENTICAL" if ok else "DIFFERENCES FOUND")
    return ok

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "setup":
        setup(sys.argv[2], sys.argv[3])
    elif cmd == "run":
        run(sys.argv[2], sys.argv[3], sys.argv[4])
    elif cmd == "compare":
        sys.exit(0 if compare(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else "") else 1)
```
