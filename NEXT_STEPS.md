# Next steps

Status on 2026-10-02 and what comes next, in order. The refactor is finished and committed (`036019b`, `d8205fe`); the layout, contracts and decisions are in `ARCHITECTURE.md`. Every output written since then carries the `Method` and `ENBW` datasets, `mask_distance` is set in `Solver.nml`, and every source has a multipole-order label.

The order below goes from small and unblocking to large. Items 2 and 3 need a decision before any code.

---

## 1. Next step: the plotting reads `Method` and `ENBW` from the HDF5

**Why first.** It is small and Python-only. It also removes two hard-coded assumptions that break as soon as manual sources or a third method exist:
- `plot_utils.read_solver_params` detects the method from `Nmodes` / `N_images` and raises `ValueError` for anything else.
- The Hann correction is a hard-coded 1.5 in `AcousticUtils.SPL_to_OASPL` and `AcousticUtils.pressure_to_PSD` (flag `Hanning_correction=True`), and the `enbw=1.5` default of `plot_single_case.compute_point_source`. Forces given in the frequency domain will have no window, so this correction must then be skipped.

**What to do:**
1. In `read_solver_params`, read `Method` when it exists (bytes: `b"Images Method"`, `b"Analytical Normal Modes"`). Keep the current detection only as a fallback for files written before 10/2026. Today `case["Method"]` holds `"ImageMethod"` / `"AnalyticalNormalModes"`, so either map the dataset to those names or store it in a new key (e.g. `case["Method_name"]`). Check every place that compares `case["Method"]`.
2. Store `case["ENBW"]` with this rule (decision needed, this is the proposal):
   - `ENBW` in the file: use it;
   - no `ENBW` and no `Method` (file written before 10/2026, always a turbine or a farm, always Hann): 1.5;
   - `Method` but no `ENBW` (forces without a window): 1.0, i.e. no correction.
3. Use `case["ENBW"]` instead of the hard-coded 1.5 in `SPL_to_OASPL`, `pressure_to_PSD` (and its 1/3 octave PSD) and `compute_point_source`.
4. **Check:** the plots and numbers of an old file and of a new one are exactly the ones obtained now, since both give 1.5.

---

## 2. Decisions pending (no code until decided)

### 2.1 Water-column boundaries of IM and ANM
- The default queries put the surface at `z = 0` and the seabed at `-depth_base`; IM and ANM use `Upper_HBC` and `Lower_HBC = Upper_HBC - depth_base`.
- Options:
  - IM and ANM override `z_surface` / `z_seabed` with their `Upper_HBC` / `Lower_HBC`, so the steps and the method agree;
  - or they stop if `UpperHBC /= 0`.
- Results change only for `UpperHBC /= 0`; both `.nml` files use 0.
- Either option is a few lines in each method.

### 2.2 The `gamma` factor (paper, Appendix B)
- Eq. B.8 has numerator 4, the code has 3 (proposal 14b).
- `gamma` multiplies `m_eff` (steel + water), while Appendix B derives it for the water-only mass.
- Best settled against the BEM case before choosing.

### 2.3 Generic sphere (later)
- `run_sphere` exists only in IM.
- A free-field kernel shared by every method would give ANM a sphere too. It changes outputs, so it is a separate decision.

---

## 3. Manual sources (the next feature)

Design already decided: `ARCHITECTURE.md` §7.1. The recipe for a new emitter is in §5.2.

**Design points still open:**
1. **Monopole strength.** The physical quantity of a monopole (volume velocity, or pressure at 1 m), and where it is stored. The decision was "a label per source, no new arrays", so a convention is needed, e.g. the strength in `Force(:, s, 1)` and components 2-3 unused.
2. **Post-processing file.** Which one `'sources'` uses: `ComputationSetupTurbine.nml`, or a new `ComputationSetupSources.nml`. Also which steps make sense (decay follows `geometry % direction`).
3. **Default output.** The output folder and file-name labels (`geometry % out_dir`, `name_tag`, `name_suffix`).

**Work items, in order:**
1. **Configuration.** `Config_ManualSources_t` in `types/ConfigTypes.f90`, `read_sources_config` in `input/InputEmitters.f90`, and the `'sources'` `CASE` in `read_input_files`. Input:
   - a namelist for a few sources;
   - an input HDF5 file (`xyz`, `order`, `Freqs`, `Force`) for many. `HDF5Utils` only has `read_hdf5_real_1d` today, so it needs readers for real 2D, integer 1D and complex 3D datasets.
2. **`emitters/ManualSources.f90`.** `ManualSources_t EXTENDS(Emitter_t)`:
   - `to_source` fills `xyz`, `Force`, `Freqs`, `order` and `geometry`, with no radiation-impedance correction;
   - `save_params` and `write_txt`;
   - `get_enbw` returns `HANN_ENBW` for time-domain forces (through `compute_rfft`) and keeps -1 for frequency-domain forces.
3. **Program and build.** The construction branch in `app/main.f90`, ending in `MOVE_ALLOC` into `active_emitter`, and the Makefile.
4. **Monopole kernels.**
   - One kernel per order in IM (images of a monopole) and in ANM (modal sum of a monopole).
   - Each method replaces its dipole-only check in `init` by a check of the orders it supports.
5. **Verification:**
   - the reference cases of `review/refactor_tools/` stay `ALL IDENTICAL`, with no extra dataset;
   - the dipoles of a turbine, exported and read back as manual sources, reproduce its `P_*` datasets;
   - a single free-field monopole and dipole (sphere with `N_images = 0`) match the closed form;
   - IM and ANM agree for monopoles as they do for dipoles.

---

## 4. Later: acoustic methods with a non-uniform environment

As decided in `ARCHITECTURE.md` §7.2:
- **Environment data.** `Config_Environment_t` gains optional arrays: c(z), ρ(z), attenuation, bathymetry H(x, y), sediment.
- **Checks.** Each method checks in `init` that the data it needs is present.
- **Boundaries.** The method publishes its surface and seabed by overriding `z_surface`, `z_seabed` and `in_water`, not the `run_*` steps.
- **Recipe.** Adding the method follows `ARCHITECTURE.md` §5.3; the `Method` dataset is automatic.

---

## 5. Small items found on the way

| Item | Where |
|---|---|
| The executable is now `off_coustics` (`9a7431b`), but three places still say `off_coustics_fortran`; the fpm executable takes its name from `fpm.toml` | `ARCHITECTURE.md` §1, the `Makefile` header comment, the header of `app/main.f90` |
| Stale comment: says `D_eff = sqrt(Diam_y*Diam_z)/2`, the code uses `sqrt(Diam_y*Diam_z)` | `SAITEC2MWFloating.f90`, around line 693 |
| SAITEC saves `wet_area = 0`, so the `W/A_wet` metric divides by zero | `plotting/plot_single_case.py`, `compute_sphere_metrics` |
| Rule 1.6 names `Types`, and the example of rule 7.3 uses `USE IOUtils`: modules that no longer exist | `Programming_style_guidelines.txt` |
| The "observers placed over z = 0 m" warning still compares with a literal 0; worth switching to `z_surface` when a method overrides it (item 2.1) | `run_spectrums` |
| Outputs rewritten many times keep dead space; `h5repack` recovers it | `.hdf5` files |
