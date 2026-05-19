# native/get_featureIDs/Makefile — build the get_featureIDs binary.
#
# Usage:
#   make              # build build/get_featureIDs (default target)
#   make clean        # remove build/
#
# Variables you can override on the command line:
#   CXX        - C++ compiler. Defaults to `c++` (the system's default C++
#                front-end: g++ on Linux, Apple Clang on macOS). On macOS
#                with Homebrew GCC you might set `CXX=g++-13`.
#   CXX_STD    - C++ standard. Defaults to c++20. get_featureIDs uses
#                std::span and other C++20 features, so this is required.
#   OPTIMIZE   - Optimization level. Defaults to -O3. Use OPTIMIZE='-O0 -g'
#                for debug builds.
#   CXXFLAGS   - Extra C++ compiler flags appended to ours.
#   CPPFLAGS   - Extra preprocessor flags appended to ours.
#   LDFLAGS    - Extra linker flags appended to ours.
#
# Dependencies (must be installed):
#   - A C++20-capable compiler (gcc 10+, clang 13+, Apple Clang 15+).
#   - zlib development headers and library.
#     Linux  : `apt install zlib1g-dev`  or  `conda install zlib`
#     macOS  : ships with Xcode Command Line Tools, or `brew install zlib`
#
# How the include and library paths are resolved:
#   - The KMC API headers and cxxopts.hpp are vendored under
#     `../external/` (kmc_api/ and cxxopts/). The -I flags below add those
#     to the include search path so the source can use
#     `#include "kmc_api/kmc_file.h"` and `#include "cxxopts.hpp"` without
#     ugly relative paths.
#   - For zlib, we try pkg-config first (works in conda environments and
#     on most Linux/macOS installations). If pkg-config isn't available
#     or can't find zlib, we fall back to plain `-lz`.
#   - If CONDA_PREFIX is set in the environment, we also add its include
#     and lib directories — this catches cases where pkg-config is missing
#     but a conda environment has the headers/libs.

# --- Toolchain ----------------------------------------------------------

CXX      ?= c++
CXX_STD  ?= c++20
OPTIMIZE ?= -O3

# --- Output directories -------------------------------------------------

BUILD_DIR := build
BIN       := $(BUILD_DIR)/get_featureIDs

# --- Source files -------------------------------------------------------

MAIN_SRC    := get_featureIDs.cpp
KMC_API_DIR := ../external/kmc_api
KMC_SOURCES := \
	$(KMC_API_DIR)/kmc_file.cpp \
	$(KMC_API_DIR)/kmer_api.cpp \
	$(KMC_API_DIR)/mmer.cpp

KMC_OBJS := \
	$(BUILD_DIR)/kmc_file.o \
	$(BUILD_DIR)/kmer_api.o \
	$(BUILD_DIR)/mmer.o

MAIN_OBJ := $(BUILD_DIR)/get_featureIDs.o

# --- Include and link configuration ------------------------------------

# Header search paths: vendored KMC API + vendored cxxopts.
BASE_CPPFLAGS := -I../external -I../external/cxxopts

# Compile flags: standard, optimization, warnings, threads.
BASE_CXXFLAGS := -std=$(CXX_STD) $(OPTIMIZE) -pthread -Wall -Wextra

# zlib via pkg-config, falling back to -lz.
# `2>/dev/null || true` so a missing pkg-config doesn't break the makefile.
ZLIB_CFLAGS := $(shell pkg-config --cflags zlib 2>/dev/null || true)
ZLIB_LIBS   := $(shell pkg-config --libs zlib 2>/dev/null || echo -lz)

# Conda-env paths (only added if CONDA_PREFIX is set). Some installs have
# zlib in $CONDA_PREFIX without registering a pkg-config file.
#
# `-Wl,-rpath,$(CONDA_PREFIX)/lib` is essential on macOS: conda-forge dylibs
# are shipped with install_name `@rpath/libfoo.dylib`, which means the
# consuming binary has to embed an `LC_RPATH` entry telling dyld where to
# look. Without it, the binary builds fine but fails at runtime with
# "dyld: Library not loaded: @rpath/libz.1.dylib ... no LC_RPATH's found".
# Conda's compiler-activation scripts inject this automatically when you
# use their compilers; we don't (we use the system `c++`), so we set it
# here. On Linux the same flag adds DT_RUNPATH which works the same way
# — harmless if not needed.
ifdef CONDA_PREFIX
BASE_CPPFLAGS += -I$(CONDA_PREFIX)/include
EXTRA_LDFLAGS := -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib
else
EXTRA_LDFLAGS :=
endif

# --- Combined flag tuples (in compile/link order) ----------------------

# Used by every compile rule.
ALL_CPPFLAGS = $(BASE_CPPFLAGS) $(CPPFLAGS)
ALL_CXXFLAGS = $(BASE_CXXFLAGS) $(ZLIB_CFLAGS) $(CXXFLAGS)
ALL_LDFLAGS  = $(EXTRA_LDFLAGS) $(LDFLAGS)
ALL_LDLIBS   = $(ZLIB_LIBS) -pthread

# --- Targets -----------------------------------------------------------

.PHONY: all clean help

all: $(BIN)

# Final link.
$(BIN): $(MAIN_OBJ) $(KMC_OBJS) | $(BUILD_DIR)
	$(CXX) $(ALL_LDFLAGS) -o $@ $(MAIN_OBJ) $(KMC_OBJS) $(ALL_LDLIBS)

# Compile rule for the main source.
$(MAIN_OBJ): $(MAIN_SRC) | $(BUILD_DIR)
	$(CXX) $(ALL_CXXFLAGS) $(ALL_CPPFLAGS) -c $< -o $@

# Compile rule for each KMC source. Pattern rule means edits to kmc_file.cpp
# rebuild only that .o, not kmer_api.o and mmer.o.
$(BUILD_DIR)/%.o: $(KMC_API_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(ALL_CXXFLAGS) $(ALL_CPPFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)

help:
	@echo "Targets:"
	@echo "  all     (default) Build $(BIN)"
	@echo "  clean             Remove $(BUILD_DIR)/"
	@echo "  help              Show this message"
	@echo ""
	@echo "Variables (override on the command line):"
	@echo "  CXX=$(CXX)"
	@echo "  CXX_STD=$(CXX_STD)"
	@echo "  OPTIMIZE=$(OPTIMIZE)"
	@echo "  CONDA_PREFIX=$(CONDA_PREFIX)"
	@echo ""
	@echo "Examples:"
	@echo "  make                         # release build with -O3"
	@echo "  make OPTIMIZE='-O0 -g'       # debug build"
	@echo "  make CXX=clang++             # force a specific compiler"
	@echo "  make -j4                     # parallel build"
