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
#   - builds/links the one file that needs HDF5 (HDF5Utils.f90) and the
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
#   make run ARGS=./my_inputs   same, reading the .nml files of ./my_inputs instead of ./input_files
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
EXE    := off_coustics

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
    $(BUILD)/SourceTypes.o \
    $(BUILD)/ConfigTypes.o \
    $(BUILD)/StringUtils.o \
    $(BUILD)/FileUtils.o \
    $(BUILD)/HDF5Utils.o \
    $(BUILD)/CSVReader.o \
    $(BUILD)/Sorting.o \
    $(BUILD)/Timing.o \
    $(BUILD)/SignalProcessing.o \
    $(BUILD)/SpecialFunctions.o \
    $(BUILD)/InputChecks.o \
    $(BUILD)/InputEmitters.o \
    $(BUILD)/InputSolver.o \
    $(BUILD)/InputPostProcessing.o \
    $(BUILD)/Args_parser.o \
    $(BUILD)/OpenFASTReader.o \
    $(BUILD)/StructureModel.o \
    $(BUILD)/DrivetrainExcitation.o \
    $(BUILD)/Emitter.o \
    $(BUILD)/WindTurbine.o \
    $(BUILD)/DTU10MWMonopile.o \
    $(BUILD)/DTU10MWFloating.o \
    $(BUILD)/SAITEC2MWFloating.o \
    $(BUILD)/TurbineFactory.o \
    $(BUILD)/WindFarm.o \
    $(BUILD)/AcousticSolver.o \
    $(BUILD)/AcousticSolver_output.o \
    $(BUILD)/AcousticSolver_steps.o \
    $(BUILD)/ImageMethod.o \
    $(BUILD)/AnalyticalNormalModes.o \
    $(BUILD)/SolverFactory.o

.PHONY: all run clean
all: $(EXE)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/Kinds.o: $(SRC_DIR)/Kinds.f90 | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/FFTW3.o: $(SRC_DIR)/utils/FFTW3.f90 | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/SourceTypes.o: $(SRC_DIR)/types/SourceTypes.f90 $(BUILD)/Kinds.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/ConfigTypes.o: $(SRC_DIR)/types/ConfigTypes.f90 $(BUILD)/Kinds.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/StringUtils.o: $(SRC_DIR)/utils/StringUtils.f90 $(BUILD)/Kinds.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/FileUtils.o: $(SRC_DIR)/utils/FileUtils.f90 $(BUILD)/Kinds.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/HDF5Utils.o: $(SRC_DIR)/utils/HDF5Utils.f90 $(BUILD)/Kinds.o | $(BUILD)
	$(H5FC) $(FFLAGS) $(HDF5_INCLUDE) $(MODFLAG) -c $< -o $@

$(BUILD)/CSVReader.o: $(SRC_DIR)/utils/CSVReader.f90 $(BUILD)/Kinds.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/Sorting.o: $(SRC_DIR)/utils/Sorting.f90 $(BUILD)/Kinds.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/Timing.o: $(SRC_DIR)/utils/Timing.f90 $(BUILD)/Kinds.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/SignalProcessing.o: $(SRC_DIR)/utils/SignalProcessing.f90 $(BUILD)/Kinds.o $(BUILD)/FFTW3.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/SpecialFunctions.o: $(SRC_DIR)/utils/SpecialFunctions.f90 $(BUILD)/Kinds.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/InputChecks.o: $(SRC_DIR)/input/InputChecks.f90 $(BUILD)/Kinds.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/InputEmitters.o: $(SRC_DIR)/input/InputEmitters.f90 $(BUILD)/Kinds.o $(BUILD)/ConfigTypes.o $(BUILD)/InputChecks.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/InputSolver.o: $(SRC_DIR)/input/InputSolver.f90 $(BUILD)/Kinds.o $(BUILD)/ConfigTypes.o $(BUILD)/InputChecks.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/InputPostProcessing.o: $(SRC_DIR)/input/InputPostProcessing.f90 $(BUILD)/Kinds.o $(BUILD)/ConfigTypes.o $(BUILD)/InputChecks.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/Args_parser.o: $(SRC_DIR)/input/Args_parser.f90 $(BUILD)/ConfigTypes.o $(BUILD)/InputEmitters.o $(BUILD)/InputSolver.o $(BUILD)/InputPostProcessing.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/OpenFASTReader.o: $(SRC_DIR)/emitters/OpenFASTReader.f90 $(BUILD)/Kinds.o $(BUILD)/StringUtils.o $(BUILD)/FileUtils.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/StructureModel.o: $(SRC_DIR)/emitters/StructureModel.f90 $(BUILD)/Kinds.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/DrivetrainExcitation.o: $(SRC_DIR)/emitters/DrivetrainExcitation.f90 $(BUILD)/Kinds.o $(BUILD)/Sorting.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/Emitter.o: $(SRC_DIR)/emitters/Emitter.f90 $(BUILD)/Kinds.o $(BUILD)/SourceTypes.o $(BUILD)/ConfigTypes.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/WindTurbine.o: $(SRC_DIR)/emitters/WindTurbine.f90 $(BUILD)/Kinds.o $(BUILD)/SourceTypes.o $(BUILD)/ConfigTypes.o $(BUILD)/StringUtils.o $(BUILD)/HDF5Utils.o $(BUILD)/CSVReader.o $(BUILD)/SignalProcessing.o $(BUILD)/OpenFASTReader.o $(BUILD)/DrivetrainExcitation.o $(BUILD)/Emitter.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/DTU10MWMonopile.o: $(SRC_DIR)/emitters/TurbineTypes/DTU10MWMonopile.f90 $(BUILD)/Kinds.o $(BUILD)/ConfigTypes.o $(BUILD)/HDF5Utils.o $(BUILD)/SignalProcessing.o $(BUILD)/StructureModel.o $(BUILD)/DrivetrainExcitation.o $(BUILD)/WindTurbine.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/DTU10MWFloating.o: $(SRC_DIR)/emitters/TurbineTypes/DTU10MWFloating.f90 $(BUILD)/Kinds.o $(BUILD)/ConfigTypes.o $(BUILD)/HDF5Utils.o $(BUILD)/SignalProcessing.o $(BUILD)/StructureModel.o $(BUILD)/DrivetrainExcitation.o $(BUILD)/WindTurbine.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/SAITEC2MWFloating.o: $(SRC_DIR)/emitters/TurbineTypes/SAITEC2MWFloating.f90 $(BUILD)/Kinds.o $(BUILD)/ConfigTypes.o $(BUILD)/HDF5Utils.o $(BUILD)/SignalProcessing.o $(BUILD)/OpenFASTReader.o $(BUILD)/StructureModel.o $(BUILD)/DrivetrainExcitation.o $(BUILD)/WindTurbine.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/TurbineFactory.o: $(SRC_DIR)/emitters/TurbineFactory.f90 $(BUILD)/WindTurbine.o $(BUILD)/DTU10MWMonopile.o $(BUILD)/DTU10MWFloating.o $(BUILD)/SAITEC2MWFloating.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/WindFarm.o: $(SRC_DIR)/emitters/WindFarm.f90 $(BUILD)/Kinds.o $(BUILD)/SourceTypes.o $(BUILD)/ConfigTypes.o $(BUILD)/StringUtils.o $(BUILD)/HDF5Utils.o $(BUILD)/SignalProcessing.o $(BUILD)/Emitter.o $(BUILD)/WindTurbine.o $(BUILD)/TurbineFactory.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/AcousticSolver.o: $(SRC_DIR)/solvers/AcousticSolver.f90 $(BUILD)/Kinds.o $(BUILD)/SourceTypes.o $(BUILD)/ConfigTypes.o $(BUILD)/Emitter.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

# Submodules depend on their parent's object (gfortran writes their .smod
# files into $(BUILD) through -J)
$(BUILD)/AcousticSolver_output.o: $(SRC_DIR)/solvers/AcousticSolver_output.f90 $(BUILD)/StringUtils.o $(BUILD)/FileUtils.o $(BUILD)/HDF5Utils.o $(BUILD)/AcousticSolver.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/AcousticSolver_steps.o: $(SRC_DIR)/solvers/AcousticSolver_steps.f90 $(BUILD)/HDF5Utils.o $(BUILD)/AcousticSolver.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/ImageMethod.o: $(SRC_DIR)/solvers/methods/ImageMethod.f90 $(BUILD)/Kinds.o $(BUILD)/SourceTypes.o $(BUILD)/ConfigTypes.o $(BUILD)/FileUtils.o $(BUILD)/HDF5Utils.o $(BUILD)/Emitter.o $(BUILD)/AcousticSolver.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/AnalyticalNormalModes.o: $(SRC_DIR)/solvers/methods/AnalyticalNormalModes.f90 $(BUILD)/Kinds.o $(BUILD)/SourceTypes.o $(BUILD)/ConfigTypes.o $(BUILD)/FileUtils.o $(BUILD)/HDF5Utils.o $(BUILD)/SpecialFunctions.o $(BUILD)/Emitter.o $(BUILD)/AcousticSolver.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/SolverFactory.o: $(SRC_DIR)/solvers/SolverFactory.f90 $(BUILD)/AcousticSolver.o $(BUILD)/ImageMethod.o $(BUILD)/AnalyticalNormalModes.o | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/main.o: $(APP_DIR)/main.f90 $(OBJS) | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(EXE): $(OBJS) $(BUILD)/main.o
	$(H5FC) $(FFLAGS) -o $@ $(OBJS) $(BUILD)/main.o $(FFTW3_FLAGS)

run: all
	./$(EXE) $(ARGS)

clean:
	rm -rf $(BUILD) $(EXE)
