# Rewrite of the plotting package (`FORTRAN/plotting`)

Specification written on 2026-10-06 by Raul Sanz Ramirez with Claude. It is addressed to a Claude Code session that will carry it out. **Follow it to the letter.** If something in this file contradicts the repository, stop and report it before changing anything (Phase 0). Do not improvise: ask.

Goal: one generic plotting framework for every emitter (single turbine, wind farm, frequency-domain sources, and whatever comes next) with no duplicated plot code, keeping every present functionality.

---

## 1. Hard rules

| # | Rule |
|---|---|
| R1 | **Nothing of the Fortran side is touched.** The calculation logic of the solver is validated and frozen. Do not edit `src/`, `app/`, `input_files/`, `Makefile`, `fpm.toml` or `ARCHITECTURE.md`. Do not run the Fortran binary (it overwrites `.hdf5` and `.dat` files). The Fortran output of the new emitter (§2) is produced by the user. |
| R2 | Work only inside `plotting/` (plus the report of §13). Never modify, import or copy from `../PYTHON` (abandoned). |
| R3 | **Every function, option and flag that exists today stays available and gives the same result** (inventory in §3). The only visible differences allowed are those listed in §10. Any other difference is a bug: find the cause and fix it. If a difference seems desirable, stop and ask. |
| R4 | The plots shown are chosen by **editing the list of calls inside `turbine.py`, `farm.py` and `source.py`**, as today. No command-line switch that selects plots (the user rejected a `--only` option). |
| R5 | Reference data in `turbine_acoustic_data/`, `farm_acoustic_data/`, `SAITEC_acoustic_data/` is **read-only**. Test outputs go to `plotting/tests/_out/` (create `plotting/tests/.gitignore` containing `_out/`; do not edit the root `.gitignore`) or to a scratch directory. `golden_manifest.json` is small and is meant to be committed by the user. |
| R6 | Style: the style of the existing Python (see §12). No new third-party dependency. |
| R7 | No `git commit`, `git push` or branch creation unless the user asks. Report the state with `git status` at the end of every phase. |
| R8 | Files are large (a turbine file is 0.9-1.6 GB). Keep the selective reading that exists (`steps=`): never load a whole file when only line-plot steps are needed. |
| R9 | Python interpreter for every test: `/home/hp/miniconda3/envs/openfast_env/bin/python` (matplotlib 3.9.4, numpy 2.0.2, h5py 3.14.0; LaTeX and dvipng are installed, so `AcademicStyle` uses `usetex`). Use `MPLBACKEND=Agg`. |

---

## 2. What Fortran writes (contract with the new emitter)

The existing turbine and farm outputs are **not changed**. Only the new emitter (`simulation_type = 'FrequencySource'`, `src/emitters/FrequencySource.f90`) changes what it writes in `save_params`. After the user applies it, a source file contains the datasets of the solver (as in every file) plus:

| Dataset | Fortran shape | h5py shape | After `read_ds` | Content |
|---|---|---|---|---|
| `Emitter` | string | `(1,)` bytes | `str` | `"FrequencySource"` (marker) |
| `Nsources` | scalar | `(1,)` | `int` | number of sources |
| `NFreqs` | scalar | `(1,)` | `int` | number of frequencies |
| `Source_xyz` | `(Nsources, 3)` | `(3, Nsources)` | `(Nsources, 3)` | positions [m] |
| `Source_order` | `(Nsources)` | `(Nsources,)` | `(Nsources,)` int | multipole order: 0 monopole, 1 dipole |
| `Source_Force` | complex `(Nfreqs, Nsources, 3)` | `(3, Nsources, Nfreqs, 2)` | see below | force on the fluid [N], peak complex amplitude, time factor `exp(-i w t)` |

- `read_ds` (current `plot_utils.py`) combines the `(real, imaginary)` pair **only for keys starting with `P_`**. `Source_Force` must be combined by the source reader explicitly: `F[..., 0] + 1j * F[..., 1]`, transposed to `(Nfreqs, Nsources, 3)`.
- The frequencies of a source file are **discrete tones** given by the user (e.g. 25, 50, 100 Hz), not the bins of an FFT. This matters for the PSD (§5, `bin_width`).
- No `ENBW` dataset (forces are not windowed), so `read_enbw` returns 1.0. Check what it returns today and keep it.
- No `WindSpeed`, `WindDir`, `AxisPos`, `BariPos`, `case_type`, `Structure_nodes`, `Depth`. The depth comes from `Lower_HBC` as for every file (`read_solver_params`: `Depth = -Lower_HBC`).
- Files written before this change (they have `Structure_nodes` of shape `(2, 3)` and no `Emitter`) are **not supported** by `source.py`: raise a clear error telling the user to regenerate them.
- Detection of the emitter type of a file: `Emitter` dataset if present; otherwise `Nturbines` means farm; otherwise turbine. Each script checks it and gives a precise error for a file of another type (today `load_turbine` already rejects farm files and `load_farm` rejects turbine files).

If `source_acoustic_data/` has no file with `Emitter`, stop and ask the user for one before Phase 5.

---

## 3. Inventory of what exists today (must be preserved)

### 3.1 Files (`plotting/`)

| File | Lines | Role |
|---|---|---|
| `turbine.py`, `farm.py` | 65, 65 | Entry scripts: CLI + the list of plot calls (one branch for one case, one for a comparison) |
| `plot_single_case.py` | 1962 | Plots of one turbine (`get_curve`, `get_ws_curves`, `plot_structure`, `get_spectrum`, `plot_spectrum`, `get_bands`, `plot_polar`, `plot_cylinder`, `plot_distance_decay`, `plot_line`, `plot_sliceXY`, `plot_sliceXZ`, `plot_sliceVertical`, `compute_sphere_metrics`, `compute_cylinder_metrics`, `compute_point_source`) |
| `plot_comparison.py` | 790 | Comparison of turbines: panels that call the functions of `plot_single_case` through `ax_in` (`select_cases`, `make_panels`, `remove_common_title_part`, `wrap_case_titles`, `ws_suffix`, `add_case_legend` and one function per plot) |
| `plot_farm_module.py` | 1249 | Plots of one farm: **a separate implementation** of every step (`get_farm_color`, `get_contour_levels`, `plot_structure`, `plot_spectrums`, `plot_polar`, `plot_cylinder`, `plot_line`, `plot_sliceXY`, `plot_sliceXZ`, `plot_sliceVertical`, `as_single_case`, `compute_cylinder_metrics`, `compute_point_source`) |
| `plot_farm_comparison.py` | 401 | Comparison of farms (`label_titles` and one function per plot) |
| `plot_utils.py` | 1271 | CLI (`parse_*_args`, `add_common_args`, `finish_args`), file discovery (`find_hdf5`, `resolve_paths`, `list_hdf5`, `describe_hdf5`), loaders (`load_turbine(s)`, `load_farm(s)`), `--ws` machinery, labels and colors, HDF5 readers (`read_ds`, `read_array`, `read_float`, `read_solver_params`, `read_enbw`, `read_turbine_params`, `read_steps`, `guess_case_type`), saving (`set_figure_name`, `save_figures`), small helpers (`get_bin_width`, `convert_to`, `is_valid`, `tex_escape`, `draw_direction_arrow`) |
| `plot_style.py` | 93 | `AcademicStyle` (uses LaTeX if available, mathtext otherwise) |
| `utils/AcousticUtils.py`, `utils/EnvironmentalData.py` | 472, 484 | Acoustic conversions (SPL, PSD, OASPL, third octave, absorption) and environment data. Physics helpers: keep them as they are |
| `README.md` | - | Usage. It mentions `plotting/requirements.txt`, which does not exist: report it, do not invent the file |

### 3.2 Command line (both scripts unless noted)

`files...` (paths or bare names, with or without `.hdf5`; one file = single case, several = comparison; none = the most recent file of `--data_dir`), `--data_dir`, `--list`, `--save`, `--output_dir`, `--format png pdf svg eps`, `--no-show`, `--verbose`; **`--ws` only in `turbine.py`**. Figures saved as `figures/<file name>/<plot>_<mode>.png` (comparisons: `figures/comparison_<a>_vs_<b>.../`).

### 3.3 Behaviours that are easy to lose (each one needs a test)

- NaN pressure = observers masked by Fortran (inside a structure, outside the water column): holes/gaps in plots, never a reason to skip a plot; excluded (and reported) in sphere/cylinder metrics.
- `is_valid` checks; `has_<step>` flags: a missing step is skipped with a message, never an exception.
- `read_steps(f, case, steps)` and `WS_STEPS`: other wind speeds of `--ws` read only `spectrums`, `polar`, `decay`, `line`.
- `--ws`: `find_ws_paths` (regex `ws(\d+(?:\.\d+)?)`, whole-name match, so `..._ws11.4_shallow` is not collected), `drop_ws_twins`, `load_ws_cases`, `set_ws_styles` (colormap, line width and z-order by speed rank, shared across panels), `get_ws_cases`, `add_ws_legends`, and the comparison label without wind speed when `--ws`.
- Comparison labels: `set_comparison_labels` / `set_farm_comparison_labels` (type, wind speed, method tag; depth added if two labels match; file name added if they still match), colors by position (`get_case_colors`).
- `guess_case_type` for turbine files without `case_type`; `describe_hdf5` for `--list`.
- `read_steps` warns when `Obs_line` does not match `Distances_line` (files written before the `run_line` fix of 2026-09-29).
- `compute_sphere_metrics` falls back to hard-coded wetted areas when `wet_area` is missing: only for turbine files.
- `plot_distance_decay(turbine_data=True)` draws external data (GE 2025) from `get_turbine_styles()`.
- `plot_farm_module.plot_polar` forces `absorption=False` (a farm has no single axis) and prints a message.
- `AcademicStyle` falls back to mathtext without LaTeX.

---

## 4. Target architecture

Three layers, and plots never see an emitter:

```
 .hdf5 --> READER (one per emitter) --> raw dict --> WRAPPER (one per emitter) --> Case --> generic PLOTS
           what Fortran saved           emitter       translates, labels,         the contract (§5)
                                        vocabulary    families (--ws), colors
```

```
plotting/
├── turbine.py  farm.py  source.py     entry scripts: CLI + the list of calls to edit
├── readers/        hdf5_io.py (read_ds, read_array, read_float, print_keys)
│                   solver.py  (read_solver_params, read_enbw, read_steps: common to every emitter)
│                   turbine.py, farm.py, source.py   (emitter-specific datasets only)
├── wrappers/       case.py (the contract, validate_case, accessors)
│                   turbine.py, farm.py, source.py   (raw -> Case, labels, --ws family, comparison labels)
├── plots/          one module per step: spectrum, polar, cylinder, decay, line, slices, sphere/metrics
│   └── specific/   exclusive functions: structure (turbine, farm), sources (optional), ...
├── compare/        panels, legends and labels of the comparison mode (works with any Case)
├── plot_utils/     cli, files (find/resolve/list), saving, style, colors, tex, conversions
├── utils/          AcousticUtils.py, EnvironmentalData.py (unchanged)
├── legacy/         temporary copy of the current package (§11); deleted only with the user's approval
└── tests/          golden harness, contract tests, CLI tests
```

Rules of dependency: entry script -> wrappers -> readers. `plots/` and `compare/` import only the contract (`wrappers/case.py`), `plot_utils/` and `utils/`; they never import a reader or a wrapper. Readers never import plots. Do not create a folder named like an entry script (`turbine/`, `farm/`, `source/`). `python plotting/turbine.py` must keep working from `FORTRAN/` and from any other directory.

---

## 5. The Case contract (`wrappers/case.py`)

A plain `dict` (the current plot code is dict-based). The wrapper builds it; `validate_case(case)` runs at the end of every wrapper and raises one `ValueError` listing every missing or malformed key. Step arrays keep the names and shapes of `read_steps` today (`P_polar`, `Obs_polar`, `Center_polar`, ..., `has_<step>`): they come from the solver and are identical for every emitter.

| Key | Type / shape | Meaning | Required |
|---|---|---|---|
| `Path`, `Name` | `Path`, `str` | file and stem | yes |
| `Emitter` | `'turbine'` \| `'farm'` \| `'source'` | type of emitter | yes |
| `label` | `str` | text of titles; replaced by a distinct label in comparisons (replaces `Case_type`) | yes |
| `Method`, `Method_name`, `Method_tag` | `str` | as `read_solver_params` today | yes |
| `c_wat`, `rho_wat`, `p_ref`, `Lower_HBC` | `float` | as today | yes |
| `Depth` | `float` | `-Lower_HBC`. **Read it only through `seabed_depth(case, x, y)`** (bathymetry will make it a function of x, y) | yes |
| `Freqs` | `(Nf,)` | frequency vector | yes |
| `ENBW` | `float` | window bandwidth in bins (1.0 without window) | yes |
| `bin_width` | `float` \| `None` | spectral bin width [Hz] for PSD. FFT-based emitters: `get_bin_width(Freqs)` as today. **Tonal sources: `None`** | yes |
| `nodes` | `(N, 3)` | flat list of points drawn over maps; NaN rows removed | yes (may be `(0, 3)`) |
| `node_group` | `(N,)` int \| `None` | group of each node (member, turbine, source) to color/label | no |
| `group_labels` | `list[str]` \| `None` | label of each group | no |
| `node_order` | `(N,)` int \| `None` | multipole order per node (`None` = all dipoles) | no |
| `origin` | `(2,)` \| `None` | reference point for observer distances. `None` = no single reference: absorption is switched off with a message (what the farm does today) | yes |
| `direction` | `float` deg \| `None` | e.g. wind direction. `None` = no arrow is drawn | yes |
| `direction_name` | `str` \| `None` | label of the arrow (`'WindDir'`) | if `direction` |
| `references` | `list` \| `None` | external reference curves for decay (GE 2025 data for turbines) | no |
| `family` | `list[Case]` \| `None` | the case at several values of a parameter (`--ws`): replaces `WS_cases` | no |
| `family_label`, `family_title`, `family_color`, `family_lw`, `family_zorder` | | replace `WindSpeed` label (`'11.4 m/s'`), legend title (`'Wind speed'`), `WS_color`, `WS_lw`, `WS_zorder` | in family members |
| `Color` | color | comparison color | in comparisons |
| `specific` | `dict` | emitter-specific data used only by exclusive functions (turbine: `Structure_nodes (Nm, Nn, 3)`, `wet_area`, `D`, `WindSpeed`, ...; farm: the list of per-turbine dicts, `Nturbines`; source: `Source_Force`, `Nsources`, `NFreqs`) | yes (may be `{}`) |

Accessors (all in `wrappers/case.py`): `seabed_depth(case, x, y)`, `emitter_label(case)`, plus whatever the plots need to read the contract without touching keys directly. Plot code uses accessors or these keys; never `case["WindDir"]`, `case["AxisPos"]`, `case["Structure_nodes"]`, `case["Case_type"]`, `case["WindSpeed"]`, `case["WS_cases"]` (legacy names).

Semantics that the generic functions must implement:

1. `origin is None` -> `absorption` becomes `False` (message), distances to the origin are not computed, and anything labelled "distance" uses what the farm uses today (verify in `plot_farm_module.plot_spectrums`).
2. `bin_width is None` -> modes `PSD` and conversions that need a bin width are not applicable: print one clear message and skip that plot (SPL per tone, OASPL and band sums remain valid).
3. `direction is None` -> no arrow, no `WindDir` text.
4. `nodes` are drawn as **points** over maps (colored by `node_group` when given; monopoles with a different marker when `node_order` says so). They are not joined by lines.
5. `family` present -> line plots draw one curve per member (`family_*` styling); absent -> the case alone.
6. A step with `has_<step>` false is skipped with a message.

---

## 6. Readers (`readers/`)

Return a raw dict with **emitter vocabulary and no interpretation** (dataset names as Fortran wrote them, Fortran order restored by `read_ds`).

- `hdf5_io.py`: `read_ds`, `read_array`, `read_float`, `print_keys` moved as they are, plus a helper that combines a complex `(…, 2)` dataset.
- `solver.py`: `read_solver_params`, `read_enbw`, `read_steps` (with the `steps=` argument), `STEPS`, `WS_STEPS`: **identical behaviour** to today.
- `turbine.py`: generic datasets (`case_type`/`Case_type`, `WindSpeed`, `WindDir`, `Depth`, `AxisPos`, `BariPos`), `Structure_nodes` and `SPECIFIC_KEYS`, `guess_case_type` for old files.
- `farm.py`: `Nturbines`, `Source_offset`, per-turbine datasets with suffix `_i` (`read_turbine_params(f, "_i")`).
- `source.py`: `Emitter`, `Nsources`, `NFreqs`, `Source_xyz`, `Source_order`, `Source_Force` (§2).

Each reader raises a precise error for a file of another emitter type.

## 7. Wrappers (`wrappers/`)

Raw dict -> Case.

- `turbine.py`: `label` = case type (TeX-escaped); `nodes` = `Structure_nodes` flattened, NaN rows dropped, `node_group` = member index; `origin` = `AxisPos`; `direction` = `WindDir`, `direction_name` = `'WindDir'`; `bin_width` = `get_bin_width(Freqs)`; `references` = the GE 2025 data (`get_turbine_styles`) so that `turbine_data=True` keeps working; `specific` = `Structure_nodes` (Nm, Nn, 3) and the type-specific keys (`wet_area`, ...). Comparison labels (`set_comparison_labels`) and `--ws` (`find_ws_paths`, `drop_ws_twins`, `load_ws_cases`, `set_ws_styles`) live here and fill `family`. The other wind speeds are loaded with the same selective `steps=WS_STEPS`.
- `farm.py`: `label` = `get_farm_label`; `nodes` = every turbine's nodes concatenated, `node_group` = turbine index, `group_labels` = `Turbine i`; `origin = None`; `direction`: see how `plot_farm_module` uses the per-turbine `WindDir` and keep it; `bin_width` as turbines; `specific` = `Turbine_parameters`, `Nturbines`, `Source_offset`. No `--ws` (the user declined it for farms). Comparison labels = `set_farm_comparison_labels`.
- `source.py`: `label` = `'FrequencySource'`; `nodes` = `Source_xyz`, `node_group` = `arange(Nsources)`, `group_labels` = `Source i`, `node_order` = `Source_order`; `origin` = `Source_xyz[0, :2]` if `Nsources == 1` else `None`; `direction = None`; **`bin_width = None`**; `references = None`; `specific` = `Source_Force (Nfreqs, Nsources, 3)`, `Nsources`, `NFreqs`. Comparison label: `"<Nsources> sources, <Nfreqs> tones, <Method_tag>"`, with the same disambiguation (depth, then file name) as the others. No `--ws`.

`apply_comparison_labels` and `get_case_colors` are shared.

---

## 8. Generic plots (`plots/`, `compare/`)

- **One function per step for every emitter**, with the **same function name, same keyword arguments and same defaults** as `plot_single_case.py` (e.g. `plot_polar(case, mode='OASPL', target_frequency=0.88, absorption=True, bands=[10., 100.], figsize=(7, 7), ax_in=None)`): copy each signature from the code, do not retype it. `ax_in` is what the comparison layer uses.
- Starting point is `plot_single_case.py` (it already has `--ws`). The code of `plot_farm_module.py` is a fork of it: **for every step, before merging, produce a written diff of the two implementations** (`difflib` on normalised text, then read it) and classify every difference as (a) cosmetic, (b) already covered by the contract (e.g. `origin None`, no `family`), or (c) real farm behaviour. Case (c) becomes a parameter of the generic function or an exclusive function, never silently dropped. Put the table in the report (§13). Measured by the author on 2026-10-06 (similarity of normalised text, single vs farm): polar 51 %, cylinder 75 %, line 41 %, sliceXY 83 %, sliceXZ 77 %, sliceVertical 72 %, spectrum 54 %. Only `plot_polar` was read in detail: its differences are the absence of `--ws`, `absorption` forced off (no single axis) and the title. The others are **not** verified.
- Exclusive functions (they stay emitter-specific, in `plots/specific/`): `plot_structure` of turbine and of farm (they draw member lines from `specific['Structure_nodes']`, with the wind-direction arrow when `direction` is not `None`); `compute_sphere_metrics` (needs the `sphere` step; the per-area lines need `specific['wet_area']`, and the hard-coded fallback applies only to `Emitter == 'turbine'` files); `plot_spectrums` of the farm only if the diff shows real behaviour that `plot_spectrum` cannot express.
- Optional new function (do it last, only if Phases 0-8 are green): `plot_sources(case)` in `plots/specific/` showing the source positions and `|Source_Force|` per source and tone, to show how an exclusive function reads `specific`.
- The comparison layer in `compare/` merges `plot_comparison.py` and `plot_farm_comparison.py`: panels over any list of Cases (`select_cases`, `make_panels`, `add_case_legend`, ...), calling the generic functions with `ax_in`. Preserve every behaviour of both (`--ws` legends in turbine comparisons, common-title removal, etc.).
- Function names may be unified (for instance `plot_spectrum` and the farm's `plot_spectrums`) provided that the entry scripts are updated and the final report gives the table *old call -> new call* of §13.

---

## 9. Entry scripts and CLI

`turbine.py`, `farm.py`, `source.py` keep the current structure: `AcademicStyle.apply()`, parse arguments, load cases, one branch for a single case and one for a comparison with the list of calls (commented alternatives included), then `save_figures` and `plt.show()`. Same flags (§3.2). `source.py`: default `--data_dir` = `FORTRAN/source_acoustic_data`, no `--ws`, `--list` shows `Emitter`, number of sources and tones, method and steps. The default list of calls of `source.py` is that of `turbine.py` minus `plot_structure` and minus anything that needs a turbine, plus `plot_sources` if it exists. `--list` for turbine and farm must print what it prints today.

---

## 10. Approved differences (exhaustive)

| ID | Difference | Where |
|---|---|---|
| D1 | Structure overlay over maps is a set of points instead of member lines | overlays of `plot_sliceXY`, `plot_sliceXZ`, `plot_sliceVertical` (turbine and farm). `plot_structure` itself keeps its lines |
| D2 | Farm calls with default `absorption` print the message of `origin is None` once per plot (the figure is the same) | farm |
| D3 | Titles/labels of source files | `source.py` only (new) |
| D4 | PSD not applicable for tonal sources: message and skip | `source.py` only (new) |
| D5 | Internal key renames of §5 (`Case_type` -> `label`, `WS_*` -> `family_*`, ...) | no visible effect |
| D6 | Messages printed by the loaders may change wording | stdout only; numbers must not change |

Anything else visible (image content, a number in a table, a file name, a flag) needs the user's approval: stop and ask.

---

## 11. Verification protocol

**Reference set** (read-only): `turbine_acoustic_data/plot_mn_IM30_ws8.0.hdf5` (monopile; its `--ws` siblings are ws4.5, ws11.4, ws15.0), `plot_fl_IM30_ws8.0.hdf5` (floating), `SAITEC_acoustic_data/saitec_ANM.hdf5`, `farm_acoustic_data/4.5_kincardine.hdf5` (and `8.0_kincardine.hdf5` for the farm comparison). Comparison cases: mn vs fl, with and without `--ws`; 4.5 vs 8.0 kincardine.

**Golden harness (`plotting/tests/golden.py`, built in Phase 1 on the LEGACY code before anything else moves):**

1. Copy the current package untouched to `plotting/legacy/` (all `.py`, `utils/`, `README.md`). It stays self-contained: the harness runs a tree by setting `sys.path` and importing the modules by their names.
2. `golden_calls.py`: the list of calls = every call of the current entry scripts (commented ones included) + every public plot function with default arguments + the mode variants (`SPL`/`PSD`/`OASPL`, `octave` True/False, `smooth`, `target_frequency`) used in the scripts, for each reference case, in single and comparison form, with and without `--ws` where it applies.
3. For each call record: success/exception, **captured stdout** (paths and timings normalised), and for every figure created: the **md5 of the PNG** (`dpi=60`) and a **structural fingerprint** (titles, axis labels, limits, number of artists by type, a hash of the data of lines, collections and images). Figures are bit-reproducible on this machine (verified: same md5 on two runs).
4. Store `golden_manifest.json` with `matplotlib`/`numpy`/`h5py` versions and whether LaTeX is used; refuse to compare if the environment differs.
5. `python plotting/tests/golden.py --tree legacy|new --record|--check`. The check reports per call: identical / differs (with the structural difference in words) / approved (maps to an ID of §10) / **unapproved**. Unapproved must be zero to pass a phase.

**Other tests:** (a) loaders: for every reference file the new Case must carry the same step arrays as the legacy dictionary (`np.array_equal` with NaN equal) and the same values for the keys that map one to one; (b) contract: `validate_case` on every wrapper; (c) CLI: `--list`, bare names with/without `.hdf5`, default newest file, `--data_dir`, `--save --no-show --output_dir <scratch>` produce the same list of file names as legacy, a file of the wrong emitter gives the precise error; (d) old files: every existing turbine/farm/SAITEC file loads; (e) source files: every generic plot runs on a source file without exception and the figures are looked at once by a human before they are recorded as golden.

---

## 12. Style of the code

Follow the existing Python: file header docstring (`Module`, `Description`, `Author`, `Institution`, `Date`), functions with aligned argument comments (`name : type = default,   # [unit] description`), docstrings with `Returns`, level-2 separators `# ---------- Title ---------- #`, English comments, no dead code or commented-out code in the new modules (except the commented alternative calls in the entry scripts, which are a feature).

---

## 13. Phases, deliverables and final report

Every phase ends with the golden `--check` of the phases before it (zero unapproved), `git status`, and a two-line summary to the user. Keep `plotting/tests/PROGRESS.md` (current phase, what is done, what comes next, open questions) up to date, so that a new session can resume the work from it and from this file.

| Phase | Work | Done when |
|---|---|---|
| 0 | **Verify this document against the repo.** Read the code named here; check every number and claim marked as verified; list the claims marked "not verified" that you did check and what you found; list discrepancies (e.g. the missing `requirements.txt`). Ask about anything blocking | report to the user, no file changed |
| 1 | Copy to `plotting/legacy/`; build the golden harness; record the manifest on legacy | `--tree legacy --check` is green twice |
| 2 | New tree: `readers/`, `wrappers/case.py`, wrappers of turbine and farm, `plot_utils/`; loaders tests of §11(a) | new Cases match legacy dictionaries |
| 3 | Generic plots, step by step, from `plot_single_case.py`, turbine only (entry script `turbine.py` on the new tree, single case) | golden green except D1 |
| 4 | `compare/` and `--ws` as `family` (turbine) | golden green for comparison and `--ws` cases |
| 5 | `readers/source.py`, `wrappers/source.py`, `source.py` and its tests; optional `plot_sources` | tests §11(e) pass, figures reviewed by the user |
| 6 | Farm: per step, written diff single vs farm, then switch the farm to the generic function | farm golden green except D1, D2 |
| 7 | Farm comparison | farm comparison golden green |
| 8 | Update `plotting/README.md` to the new layout and usage. **Do not delete `legacy/` yet** | README matches the code |
| 9 | Final report and demonstration | see below |

**Final report (`plotting/rewrite_report.md`)** and demonstration to the user: (1) the tree as built; (2) for each of `turbine.py`, `farm.py`, `source.py`, how to run it and the list of calls with one example command each, including `--ws` and a comparison; (3) the table *old call -> new call* for every public function; (4) the table of equivalence: every golden call with its result and the approved differences used, and the per-step single-vs-farm diff of §8; (5) the claims of this document that turned out wrong or incomplete; (6) open issues and suggestions. The user decides whether to delete `plotting/legacy/`.

---

## 14. Future changes this design must absorb (do not implement now)

Monopoles and dipoles (`Source_order` -> `node_order`), time-domain input (an emitter whose Freqs come from an FFT sets `bin_width` and `ENBW`), new acoustic methods (`Method` is autodetected in `read_solver_params`), bathymetry (`seabed_depth(case, x, y)` is the only reader of `Depth`), and any new emitter (one reader, one wrapper, one entry script, and its exclusive functions if any). If a requirement of this list cannot be met by adding a reader and a wrapper, the contract is wrong: say so.
