#
# Copyright 2026 Yunseo Hwang and Giwon Song.
#
# SPDX-License-Identifier: Apache-2.0

# specify compiler 
CC = g++
# compile option ( -02 : optimization level 2, -Wall : show all warnings )
CFLAGS = -O2 -Wall

# target names
TARGETS = qaoa_gen qke_gen

all: $(TARGETS)

# rule to build each generator
qaoa_gen: E3LIN2.cpp
	$(CC) $(CFLAGS) -o $@ $<

qke_gen: RandomPauliFeatureMap.cpp
	$(CC) $(CFLAGS) -o $@ $<

# make clean: remove generated binaries
clean:
	rm -f $(TARGETS) gen

