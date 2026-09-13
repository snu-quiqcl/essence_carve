#
# Copyright 2026 Yunseo Hwang and Giwon Song.
#
# SPDX-License-Identifier: Apache-2.0

# Makefile for OpenMP Offloading using clang++
# This file compiles C++ code with OpenMP directives for GPU offloading.

# --- Compiler and Flags ---
# Use clang++ as the compiler, which has excellent support for OpenMP offloading.
CC := g++

# Target executable name
TARGET_EXEC := code

# Directories
_IDIR := .
_SDIR := .
_ODIR := obj

# --- Common Compiler Flags ---
CXXFLAGS := -std=c++20 -Wall -I$(_IDIR) -O3 -g  -fopenmp

# --- Source and Object Files ---
_MAIN := main
_MODEL := tool/welford tool/timer tool/rng chlib/chstabilizer chlib/core gates ch_runner extended_stabilizer_state file_reader operators
_MODEL_HEADERS_ONLY := tool/mystr

MODELSRC := $(patsubst %, $(_SDIR)/%.cpp, $(_MODEL))
MAINSRC  := $(patsubst %, $(_SDIR)/%.cpp, $(_MAIN))

MODELOBJ := $(patsubst %, $(_ODIR)/%.o, $(_MODEL))
MAINOBJ  := $(patsubst %, $(_ODIR)/%.o, $(_MAIN))

ALL_OBJS := $(MAINOBJ) $(MODELOBJ)
DEPS := $(patsubst %, $(_IDIR)/%.hpp, $(_MODEL) $(_MODEL_HEADERS_ONLY))

# --- Build Rules ---
all: $(TARGET_EXEC)

$(TARGET_EXEC): $(ALL_OBJS)
	@echo "==> Linking target executable: $@"
	$(CC) -o $@ $^ $(CXXFLAGS)

$(_ODIR)/%.o: $(_SDIR)/%.cpp $(DEPS) | $(_ODIR)
	@echo "==> Compiling $< -> $@"
	$(CC) -c -o $@ $< $(CXXFLAGS)

$(_ODIR):
	@echo "==> Creating object directories..."
	@mkdir -p obj/chlib obj/tool

# --- Housekeeping Rules ---
clean:
	@echo "==> Cleaning up object files and executable..."
	@rm -rf $(_ODIR) $(TARGET_EXEC)

rebuild: clean all

.PHONY: all clean rebuild
