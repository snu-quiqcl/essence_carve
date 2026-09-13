/*
 * Copyright 2026 Yunseo Hwang and Giwon Song.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "rng.hpp"

CpuRngEngine::CpuRngEngine() { set_random_seed(); }
CpuRngEngine::CpuRngEngine(size_t seed) { set_seed(seed); }

void CpuRngEngine::set_random_seed() {
	std::random_device rd;
	set_seed(rd());
}

void CpuRngEngine::set_seed(size_t seed) {
	rng.seed(seed);
	initial_seed_ = seed;
}

size_t CpuRngEngine::initial_seed(void) { return initial_seed_; }