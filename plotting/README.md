# plotting — figures of the Fortran outputs

Plots the `.hdf5` files written by the Fortran code: single turbines (`turbine_acoustic_data/`), comparisons
of several turbines, and wind farms (`farm_acoustic_data/`). Same figures, style and functions as the former
`../PYTHON` plotting code, which this folder no longer needs.

## Installation

```bash
pip install -r plotting/requirements.txt
```

LaTeX (`latex` and `dvipng`) is optional: without it (e.g. on a cluster) the text falls back to Matplotlib's
mathtext with a warning.

## Usage (from `FORTRAN/`)

```bash
python plotting/turbine.py                                            # most recent file of turbine_acoustic_data/
python plotting/turbine.py plot_mn_ANM24_ws11.4                       # one file: single case
python plotting/turbine.py plot_mn_ANM24_ws11.4 plot_fl_ANM231_ws11.4 # two or more files: comparison
python plotting/turbine.py plot_mn_IM30_ws11.4 --ws                   # line plots at every wind speed of the case
python plotting/turbine.py --list                                     # files, date, method and saved steps
python plotting/farm.py                                               # most recent file of farm_acoustic_data/
python plotting/farm.py Farm_2_DTU10MN_ANM24_ws8.0-11.4 --save
```

A file can be given as a path or as a bare name (with or without `.hdf5`), which is searched in `--data_dir`.

| Option | |
|---|---|
| `--ws` | turbine script only: the line plots (spectrum, polar, decay, line) show every wind speed of each file, see below |
| `--save` | saves every figure in `figures/<file name>/<plot>_<mode>.png` (`figures/comparison_<a>_vs_<b>.../` for comparisons) |
| `--output_dir DIR` | other folder for `--save` |
| `--format png pdf ...` | file format(s) of the saved figures (default `png`) |
| `--no-show` | does not open the windows (clusters, batch runs) |
| `--data_dir DIR` | folder of the `.hdf5` files (bare names, `--list`, default file) |
| `--verbose` | prints the datasets of every loaded file |

Every file on its own, e.g. on a cluster:

```bash
for f in turbine_acoustic_data/*.hdf5; do python plotting/turbine.py "$f" --save --no-show; done
```

## Choosing the plots

As in `../PYTHON`, the plots shown are the calls listed in `main()` of `turbine.py` (one branch for a
single case, one for a comparison) and of `farm.py`: comment, add or change their arguments there. Every
function keeps the signature of the Python version (`mode`, `filter_under`, `absorption`, `target_frequency`,
...). Steps missing in a file (older runs, interrupted runs) are skipped with a message.

## Wind speeds (`--ws`)

With `--ws` every file stands for all the wind speeds of its case: the files of the same folder whose name only
differs in the number after `ws` are loaded too (`plot_mn_IM30_ws11.4` -> `..._ws4.5`, `..._ws8.0`, `..._ws15.0`;
`plot_fl_IM30_ws11.4` does not collect the `_shallow` files). The plots that draw curves show one per wind
speed; the maps and the sphere/cylinder metrics keep using the given file. Nothing else changes in `turbine.py`.

| Plot | Color | Line style |
|---|---|---|
| `plot_distance_decay`, `plot_line` | wind speed | - |
| `plot_spectrum` | wind speed | observer (observers at the same distance: energetic mean and range, shaded) |
| `plot_polar` | wind speed | frequency band (OASPL) |

The colors run along `viridis` with the wind speed, and so does the line width (the faster, the thicker and the
further down it is drawn), so curves that coincide, e.g. every speed above rated, stay visible as concentric
lines. A comparison draws one panel per file given, each with its own wind speeds and the same color for the same
speed; a file whose wind speeds were already requested by another one is loaded only once. The other wind
speeds only read the steps they draw (`WS_STEPS` in `plot_utils.py`), a few tens of MB instead of the GB of a
whole file. Without `ws<speed>` in the name the file is plotted alone, and a wind speed whose file lacks a step is
left out of that plot, both with a message. The figures are saved with the suffix `_allws`
(`decay_OASPL_allws.png`), so they do not overwrite the plain ones.

## Files

| File | Content |
|---|---|
| `turbine.py`, `farm.py` | entry scripts: CLI, loading and list of plots |
| `plot_utils.py` | CLI, `.hdf5` discovery and loading (`read_ds`, `load_turbine`, `load_farm`), saving, helpers |
| `plot_single_case.py` | single-turbine plots and sphere/cylinder metrics |
| `plot_comparison.py` | several turbines side by side (up to 3 panels per row); the equivalent point source overlays them |
| `plot_farm_module.py` | wind farm plots |
| `plot_style.py` | `AcademicStyle` (fonts, colormaps, colors, line styles) |
| `utils/` | copies of `AcousticUtils` and `EnvironmentalData` (hearing thresholds, Ge 2025 data) |

## Observers without data (NaN)

The Fortran code writes NaN (every frequency) at the observers closer than 5 m to a source node, i.e. inside
the structure. The plots treat them as points without data: holes in the contours, gaps in the curves, left
out of the spectra (with a message) and of the sphere/cylinder metrics (the excluded area is reported). The
OASPL of such an observer stays NaN (frequency sums use `np.sum`, not `np.nansum`).

## Window correction (ENBW) and method

The pressure of every frequency bin is the amplitude of a tone: the Fortran FFT (`compute_rfft`) applies a
Hann window scaled for amplitude. Every sum over bins (OASPL, PSD, band powers of the equivalent point source)
is therefore divided by the equivalent noise bandwidth of the window, read from the file (`case["ENBW"]`,
`plot_utils.read_enbw`):

| File | ENBW [bins] |
|---|---|
| with the `ENBW` dataset (written since 10/2026 when the forces come from a windowed FFT) | its value, 1.5 |
| without `ENBW` nor `Method` (written before 10/2026: always turbines or farms) | 1.5 |
| with `Method` but without `ENBW` (forces without a window, e.g. manual sources given in frequency) | 1.0, no correction |

The acoustic method is the `Method` dataset (`case["Method_name"]`, e.g. "Images Method"); older files are
recognised by `N_images` / `Nmodes`.
