# =====================================================================
# Makefile — standalone build for Off_coustics_fortran
# =====================================================================
# Drop-in replacement for `fpm build` / `fpm run`, for machines (e.g.
# cluster login/compute nodes) where installing fpm is inconvenient.
#
# It reproduces exactly what fpm already does for this project:
#   - compiles every module in the same dependency order fpm resolves
#     (verified against `fpm build`'s own build log and each file's
#     USE statements)
#   - builds/links the one file that needs HDF5 (IOUtils.f90) and the
#     final executable using `h5fc`/`h5pfc` itself as the compiler —
#     that wrapper is HDF5's own drop-in replacement for the plain
#     Fortran compiler and already knows its own include/library paths,
#     so this adapts automatically to whatever HDF5 module is loaded
#     without us having to parse any command it prints. (An earlier
#     version of this Makefile parsed `h5fc -show`'s output instead;
#     that broke on HPC HDF5 wrappers whose `-show` just echoes its own
#     arguments back instead of printing the real underlying command.)
#   - links FFTW3 via a plain -lfftw3, exactly like the
#     `link = ["fftw3"]` line in fpm.toml
#
# Usage:
#   make                        builds ./off_coustics_fortran
#   make FFLAGS="-O2 -g"        override optimization/debug flags
#   make FC=ifort               override the compiler
#   make run                    build, then run from this directory
#   make clean                  remove build_make/ and the executable
#
# Cluster usage: load whatever module provides a Fortran-capable HDF5
# (so that `h5fc` or `h5pfc` is on PATH) and a module providing FFTW3
# (so that -lfftw3 resolves, e.g. via LIBRARY_PATH set by the module)
# BEFORE running make. See the project README/chat notes for the exact
# `module load` commands used on a given cluster.
# =====================================================================

# NOTE: plain `:=` here, not `?=` — GNU Make's built-in rules already
# give FC a default value of "f77" (which may not exist, or may not be
# gfortran, on a cluster), so `?=` alone would never take effect. A
# plain assignment is still safely overridden by `make FC=...` on the
# command line, since command-line variables always win regardless of
# the assignment operator used here.
FC     := gfortran
FFLAGS := -O3 -fopenmp
BUILD  := build_make
EXE    := off_coustics_fortran

# --- HDF5 detection (mirrors fpm's "hdf5" metapackage) ----------------
# h5fc/h5pfc is the wrapper compiler that ships with any Fortran-enabled
# HDF5 install. We call it directly (instead of $(FC)) to compile the
# one file that needs HDF5 and to link the final executable, so it
# injects its own -I/-L/library flags internally — no parsing needed.
H5FC := $(shell command -v h5fc 2>/dev/null || command -v h5pfc 2>/dev/null)
ifeq ($(strip $(H5FC)),)
$(error Could not find h5fc or h5pfc on PATH. Load a Fortran-enabled HDF5 module first, e.g. "module load hdf5")
endif

# Some HDF5 wrappers (seen on HPC/EasyBuild installs) link correctly but
# don't inject the -I flag needed to actually FIND hdf5.mod while
# compiling, since that's not something the C-oriented CPATH env var
# covers for Fortran module lookup. Derive the include dir ourselves
# from where h5fc/h5pfc lives ($PREFIX/bin/h5fc -> $PREFIX/include),
# which is the standard layout for every HDF5 install we've seen.
# Override with `make HDF5_INCLUDE=-I/other/path` if that guess is wrong
# on a given site.
H5_PREFIX    := $(patsubst %/bin/,%,$(dir $(H5FC)))
HDF5_INCLUDE ?= -I$(H5_PREFIX)/include

# --- FFTW3 (mirrors fpm.toml's `link = ["fftw3"]`) --------------------
FFTW3_FLAGS := -lfftw3

# --- Module search/output directory -----------------------------------
MODFLAG := -J$(BUILD) -I$(BUILD)

# --- Source files, in the exact dependency order fpm already builds --
SRC_DIR := src
APP_DIR := app

OBJS := \
    $(BUILD)/Kinds.o \
    $(BUILD)/FFTW3.o \
    $(BUILD)/Types.o \
    $(BUILD)/MathUtils.o \
    $(BUILD)/IOUtils.o \
    $(BUILD)/Args_parser.o \
    $(BUILD)/WindTurbine.o \
    $(BUILD)/DTU10MWFloating.o \
    $(BUILD)/DTU10MWMonopile.o \
    $(BUILD)/SAITEC2MWFloating.o \
    $(BUILD)/WindFarm.o \
    $(BUILD)/AcousticSolver.o \
    $(BUILD)/ImageMethod.o \
    $(BUILD)/AnalyticalNormalModes.o

.PHONY: all run clean
all: $(EXE)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/Kinds.o: $(SRC_DIR)/Kinds.f90 | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/FFTW3.o: $(SRC_DIR)/utils/FFTW3.f90 | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/Types.o: $(SRC_DIR)/Types.f90 $(BUILD)/Kinds.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/MathUtils.o: $(SRC_DIR)/utils/MathUtils.f90 $(BUILD)/FFTW3.o $(BUILD)/Kinds.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/IOUtils.o: $(SRC_DIR)/utils/IOUtils.f90 $(BUILD)/Kinds.o | $(BUILD)
	$(H5FC) $(FFLAGS) $(HDF5_INCLUDE) $(MODFLAG) -c $< -o $@

$(BUILD)/Args_parser.o: $(SRC_DIR)/Args_parser.f90 $(BUILD)/Kinds.o $(BUILD)/Types.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/WindTurbine.o: $(SRC_DIR)/core/WindTurbine.f90 $(BUILD)/Kinds.o $(BUILD)/Types.o $(BUILD)/IOUtils.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/WindFarm.o: $(SRC_DIR)/core/WindFarm.f90 $(BUILD)/Kinds.o $(BUILD)/Types.o $(BUILD)/IOUtils.o $(BUILD)/WindTurbine.o $(BUILD)/DTU10MWMonopile.o $(BUILD)/DTU10MWFloating.o $(BUILD)/SAITEC2MWFloating.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/AcousticSolver.o: $(SRC_DIR)/core/AcousticSolver.f90 $(BUILD)/Kinds.o $(BUILD)/Types.o $(BUILD)/IOUtils.o $(BUILD)/WindTurbine.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/DTU10MWFloating.o: $(SRC_DIR)/core/TurbineTypes/DTU10MWFloating.f90 $(BUILD)/Kinds.o $(BUILD)/Types.o $(BUILD)/WindTurbine.o $(BUILD)/MathUtils.o $(BUILD)/IOUtils.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/DTU10MWMonopile.o: $(SRC_DIR)/core/TurbineTypes/DTU10MWMonopile.f90 $(BUILD)/Kinds.o $(BUILD)/Types.o $(BUILD)/WindTurbine.o $(BUILD)/MathUtils.o $(BUILD)/IOUtils.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/SAITEC2MWFloating.o: $(SRC_DIR)/core/TurbineTypes/SAITEC2MWFloating.f90 $(BUILD)/Kinds.o $(BUILD)/Types.o $(BUILD)/WindTurbine.o $(BUILD)/MathUtils.o $(BUILD)/IOUtils.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/ImageMethod.o: $(SRC_DIR)/core/AcousticMethods/ImageMethod.f90 $(BUILD)/Kinds.o $(BUILD)/Types.o $(BUILD)/IOUtils.o $(BUILD)/WindTurbine.o $(BUILD)/AcousticSolver.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/AnalyticalNormalModes.o: $(SRC_DIR)/core/AcousticMethods/AnalyticalNormalModes.f90 $(BUILD)/Kinds.o $(BUILD)/Types.o $(BUILD)/IOUtils.o $(BUILD)/MathUtils.o $(BUILD)/WindTurbine.o $(BUILD)/AcousticSolver.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/main.o: $(APP_DIR)/main.f90 $(OBJS) | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(EXE): $(OBJS) $(BUILD)/main.o
	$(H5FC) $(FFLAGS) -o $@ $(OBJS) $(BUILD)/main.o $(FFTW3_FLAGS)

run: all
	./$(EXE)

clean:
	rm -rf $(BUILD) $(EXE)
