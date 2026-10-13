# plotting — figures of the Fortran outputs

Plots the `.hdf5` files written by the Fortran code: single turbines (`turbine_acoustic_data/`), wind farms
(`farm_acoustic_data/`) and point sources, with forces given in frequency or in time (`source_acoustic_data/`), alone
or compared. Every
emitter is drawn by the same generic plots: one entry script per emitter reads its files, and the plots never
see which emitter they draw.

## Requirements

Python >= 3.9 with `numpy`, `matplotlib` and `h5py`. LaTeX (`latex` and `dvipng`) is optional: without it (e.g.
on a cluster) the text falls back to Matplotlib's mathtext with a warning.

## Usage (from `FORTRAN/`)

```bash
python plotting/turbine.py                                            # most recent file of turbine_acoustic_data/
python plotting/turbine.py plot_mn_ANM24_ws11.4                       # one file: single case
python plotting/turbine.py plot_mn_ANM24_ws11.4 plot_fl_ANM231_ws11.4 # two or more files: comparison
python plotting/turbine.py plot_mn_IM30_ws11.4 --ws                   # line plots at every wind speed of the case
python plotting/turbine.py --list                                     # files, date, method and saved steps
python plotting/farm.py                                               # most recent file of farm_acoustic_data/
python plotting/farm.py Farm_2_DTU10MN_ANM24_ws8.0-11.4 --save
python plotting/farm.py 4.5_kincardine 15.0_kincardine                # two or more farms: comparison
python plotting/source.py                                             # most recent file of source_acoustic_data/
python plotting/source.py Source_Ns2_Nf3_ANM9 --save --no-show
python plotting/source.py Source_Ns1_Nt2000_ANM8                        # a TimeSource file
```

A file can be given as a path or as a bare name (with or without `.hdf5`), which is searched in `--data_dir`.
Each script reads only the files of its emitter: another file stops it with a message that names the right
script.

| Option | |
|---|---|
| `--ws` | `turbine.py` only: the line plots (spectrum, polar, decay, line) show every wind speed of each file, see below |
| `--save` | saves every figure in `figures/<file name>/<plot>_<mode>.png` (`figures/comparison_<a>_vs_<b>.../` for comparisons) |
| `--output_dir DIR` | other folder for `--save` |
| `--format png pdf ...` | file format(s) of the saved figures (default `png`) |
| `--no-show` | does not open the windows (clusters, batch runs) |
| `--data_dir DIR` | folder of the `.hdf5` files (bare names, `--list`, default file) |
| `--list` | lists the `.hdf5` files of `--data_dir`, newest first, with what they contain and their steps |
| `--verbose` | prints the datasets of every loaded file |

Every file on its own, e.g. on a cluster:

```bash
for f in turbine_acoustic_data/*.hdf5; do python plotting/turbine.py "$f" --save --no-show; done
```

## Choosing the plots

The plots shown are the calls listed in `main()` of `turbine.py`, `farm.py` and `source.py` (one branch for a
single case, one for a comparison): comment, add or change their arguments there. The single-case calls are
`p.<plot>(case, ...)` (package `plots`) and the comparison calls `pc.<plot>(cases, ...)` (package `compare`), with
the same names and arguments for every emitter:

| Plot | |
|---|---|
| `plot_spectrum` | spectrum of the observers: `mode='SPL'` or `'PSD'`, 1/3 octave (`octave=True`) or narrowband (`smooth`, `smooth_fraction`), `do_thresholds`, `filter_under` / `filter_over` |
| `plot_polar` | polar directivity: OASPL per frequency band (`bands`) or one frequency (`mode`, `target_frequency`) |
| `plot_cylinder` | unwrapped cylinder (azimuth x depth) |
| `plot_distance_decay` | decay with distance, 1/r and 1/r² references; `turbine_data=True` adds measured turbine noise (Ge 2025, turbines only) |
| `plot_line` | magnitude along the line P1 -> P2 |
| `plot_sliceXY`, `plot_sliceXZ`, `plot_sliceVertical` | maps of the slices, with the nodes of the emitter drawn as points (`structure`) |
| `plot_structure` | members and wind direction of a turbine or of every turbine of a farm (`mode='3D'`, `'xy'`, `'xz'`) |
| `plot_sources` | positions of point sources (dipoles as circles, monopoles as squares) and their amplitude at every tone, \|F\| [N] of a dipole or \|Q\| [m^3/s] of a monopole (bars up to 20 tones, lines above; with both kinds the monopoles go to a second axis on the right, hatched or dashed); a TimeSource also has its series in time and the amplitude of its spectrum |
| `compute_sphere_metrics` | free-field radiation metrics of the sphere (Image Method) |
| `compute_cylinder_metrics` | power flow through the cylinder |
| `compute_point_source` | equivalent point source (free field and in situ), with its figures |

Every plot takes `mode` (`'OASPL'`, `'SPL'`, `'PSD'`, `'ABS'`, `'REAL'`, `'IMAG'`, `'PHASE'`), `target_frequency`,
`absorption`, `filter_under`, `filter_over` and `ax_in` where it applies. Steps missing in a file (older runs,
interrupted runs) are skipped with a message. Mode `'SWL'` is not implemented: the maps print a message and stay
empty, the decay, line and polar plots fall back to SPL, and the spectrum stops with an error.

## Wind speeds (`--ws`)

With `--ws` every file stands for all the wind speeds of its case: the files of the same folder whose name only
differs in the number after `ws` are loaded too (`plot_mn_IM30_ws11.4` -> `..._ws4.5`, `..._ws8.0`, `..._ws15.0`;
`plot_fl_IM30_ws11.4` does not collect the `_shallow` files). They form the family of the case
(`case["family"]`): the plots that draw curves show one per wind speed; the maps and the sphere/cylinder metrics
keep using the given file. Nothing else changes in `turbine.py`.

| Plot | Color | Line style |
|---|---|---|
| `plot_distance_decay`, `plot_line` | wind speed | - |
| `plot_spectrum` | wind speed | observer (observers at the same distance: energetic mean and range, shaded) |
| `plot_polar` | wind speed | frequency band (OASPL) |

The colors run along `viridis` with the wind speed, and so does the line width (the faster, the thicker and the
further down it is drawn), so curves that coincide, e.g. every speed above rated, stay visible as concentric
lines. A comparison draws one panel per file given, each with its own wind speeds and the same color for the same
speed; a file whose wind speeds were already requested by another one is loaded only once. The other wind
speeds only read the steps they draw (`WS_STEPS` in `readers/solver.py`), a few tens of MB instead of the GB of a
whole file. Without `ws<speed>` in the name the file is plotted alone, and a wind speed whose file lacks a step is
left out of that plot, both with a message. The figures are saved with the suffix `_allws`
(`decay_OASPL_allws.png`), so they do not overwrite the plain ones.

## Wind farms

`farm.py` takes one file (one branch of `main()`) or several (comparison, one panel per farm, labelled with its
turbines, wind speed and method, e.g. `5 x DTU10MWFloating, 11.4 m/s, ANM98`). It has the same plots and options
as the turbine script, except `--ws`, with the same titles and layout (the label of a farm is its turbines, e.g.
`5 x DTU10MWFloating`):

- A farm has no single axis, so the distances of the sound have no origin: absorption is switched off with a
  message. The observers of the spectrum at the same distance to the centre of the farm (the mean turbine axis)
  get their energetic and dB means.
- `plot_structure` draws every turbine with its own wind direction; the maps draw the nodes of every turbine.
- `compute_cylinder_metrics` and `compute_point_source`: the farm as a whole through the cylinder around it, which
  must enclose every turbine. A farm has no sphere, so there are no free-field metrics: the source level is the
  in situ one, with the bias of the cylinder `|p|^2` power (see `compute_point_source` in `plots/metrics.py`).

## Point sources

`source.py` plots the files of the emitters `FrequencySource` and `TimeSource`: monopoles or dipoles at the positions of
`Source_xyz` (`Source_order`: 0 monopole, 1 dipole), with the complex amplitude of `Source_Force` at every frequency:
the force [N] of a dipole, or the volume velocity Q [m^3/s] of a monopole in its first component (the other two 0). Same plots and options as the turbine
script, except `--ws` and the turbine-only ones (`plot_structure`, `turbine_data`). The emitter of a file is its
`Emitter` dataset; a file of another emitter is rejected with the script that reads it.

- **FrequencySource** (forces given per frequency): the frequencies are tones, not the bins of an FFT. They have a
  bin width only when they are uniformly spaced (their spacing, e.g. a `'linear'` grid of `fmin`, `fmax`, `Nfreqs`);
  otherwise (a `'log'` grid, a hand-made list, a file of unevenly spaced frequencies) the PSD, the 1/3 octave bands
  and `compute_point_source` are not applicable and are skipped with a message. SPL per tone, OASPL and the band
  sums of the polar plot remain valid. With a bin width each tone stands for the power of a band `df` wide: that is
  what a sampled continuous force spectrum S_F means (|F_k| = sqrt(2 S_F(f_k) df)), while for isolated tones the PSD
  depends on `df` and has no physical meaning (use the SPL). The way the input was given (list, grid, files,
  expressions of `f`) is not in the file and does not matter to the plotting: the same input gives the same
  datasets in every mode. A comparison labels the file `<N> sources, <Nf> tones, <method>`.
- **TimeSource** (forces given as real series in time): `Source_Force` is the Hann-windowed FFT of the series
  `Source_Force_time` (Nt, Nsources, 3) on the times `Source_time`, so its frequencies are FFT bins: the bin width is
  the frequency resolution `1 / (Nt * Time_step)` (whatever frequencies were kept with `skipf`, `fmin` and `fmax`),
  `ENBW` is 1.5 as for a turbine, and the PSD, the 1/3 octave bands and `compute_point_source` work. A comparison
  labels the file `<N> sources, <Nt> samples, <method>`.
- A single source is the origin of the distances (absorption); several sources have none (absorption off, with a
  message, as for a farm).
- `plot_sources` draws their positions (XZ plane) and the amplitude of their force at every tone (bars, or one line
  per source above 20 tones, on a logarithmic axis if they span two decades or more); for a TimeSource, the force
  series in time and the amplitude of its spectrum (lines). With more than 10 sources (e.g. read from a
  `source_file`) the legends only give their number. The maps draw the sources as points (monopoles as squares).
- Files written before the `Emitter` dataset (06/10/2026) are not supported: regenerate them.

## Layout

```
plotting/
├── turbine.py  farm.py  source.py     entry scripts: command line + the list of calls to edit
├── readers/        one reader per emitter (turbine, farm, source) + hdf5_io, solver (common datasets), describe (--list)
├── wrappers/       case.py: the Case contract and its accessors; one wrapper per emitter (raw dict -> Case, labels, --ws)
├── plots/          generic plots, one module per step (spectrum, polar, cylinder, decay, line, slices, metrics)
│   └── specific/   exclusive functions: structure (turbine, farm), sphere metrics, sources
├── compare/        comparison of any emitter: panels, titles and legends, one function per plot
├── plot_utils/     command line, files, saving, style (AcademicStyle), colors, LaTeX escaping, conversions, drawing
├── utils/          AcousticUtils (SPL, PSD, OASPL, 1/3 octave, absorption), EnvironmentalData (thresholds, Ge 2025)
```

Three layers: `.hdf5` -> reader (what Fortran saved) -> wrapper (the Case) -> plots. The Case is a plain dict
described in `wrappers/case.py` (`validate_case` checks it): solver data and steps as Fortran saved them, plus
`label`, `nodes`, `origin`, `centre`, `direction`, `bin_width`, `family`, `Color` and the emitter-specific data in
`specific`, which only the exclusive functions read. The water depth is read only through
`seabed_depth(case, x, y)`. A new emitter needs a reader, a wrapper and an entry script (and its exclusive
functions, if any); the generic plots do not change.

## Observers without data (NaN)

The Fortran code writes NaN (every frequency) at the observers closer than 5 m to a source node, i.e. inside
the structure. The plots treat them as points without data: holes in the contours, gaps in the curves, left
out of the spectra (with a message) and of the sphere/cylinder metrics (the excluded area is reported). The
OASPL of such an observer stays NaN (frequency sums use `np.sum`, not `np.nansum`).

## Window correction (ENBW) and method

The pressure of every frequency bin is the amplitude of a tone: the Fortran FFT (`compute_rfft`) applies a
Hann window scaled for amplitude. Every sum over bins (OASPL, PSD, band powers of the equivalent point source)
is therefore divided by the equivalent noise bandwidth of the window, read from the file (`case["ENBW"]`,
`readers.solver.read_enbw`):

| File | ENBW [bins] |
|---|---|
| with the `ENBW` dataset (written since 10/2026 when the forces come from a windowed FFT: turbines, farms and TimeSource) | its value, 1.5 |
| without `ENBW` nor `Method` (written before 10/2026: always turbines or farms) | 1.5 |
| with `Method` but without `ENBW` (forces without a window, e.g. sources given in frequency) | 1.0, no correction |

The acoustic method is the `Method` dataset (`case["Method_name"]`, e.g. "Images Method"); older files are
recognised by `N_images` / `Nmodes`.

## Tests

```bash
python plotting/tests/make_synthetic.py      # synthetic inputs of some scenarios (tests/_out, once per machine)
python plotting/tests/golden.py --check      # every plot call against the recorded golden figures
```

`golden_calls.py` lists the calls and the scenarios (`golden.py --list`); `golden_manifest.json` holds their
recorded outputs, figures and fingerprints, from the package before the rewrite of 10/2026. The check reports every
call as identical, approved (a difference accepted with the rewrite, e.g. the nodes of the maps drawn as points)
or unapproved, and fails on the latter. It runs one scenario at a time (a comparison of turbines takes ~3.5 GB)
and needs the same versions of matplotlib, numpy and h5py, and LaTeX, as the manifest. A new plot or scenario is
added to `golden_calls.py` and recorded with `golden.py --record --only <scenario>` once its figures (in
`tests/_out/golden/`) have been looked at.
