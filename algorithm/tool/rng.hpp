/*
 * Copyright 2026 Yunseo Hwang and Giwon Song.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef RNG_HPP
#define RNG_HPP

#include <cstdint>
#include <random>

class CpuRngEngine {
	public:
		CpuRngEngine();
		explicit CpuRngEngine(size_t seed);

		void set_random_seed();
		void set_seed(size_t seed);
		size_t initial_seed(void);

		inline double rand(double a, double b){ return std::uniform_real_distribution<double>(a, b)(rng); }
		inline double rand(double b){ return rand(double(0), b); }
		inline double rand(){ return rand(0, 1); }

		template <typename Integer, typename = std::enable_if_t<std::is_integral<Integer>::value>>
		Integer rand_int(Integer a, Integer b) {
			return std::uniform_int_distribution<Integer>(a, b)(rng);
		}

		template <typename Float, typename = std::enable_if_t<std::is_floating_point<Float>::value>>
		size_t rand_int(const std::vector<Float> &probs) {
			return std::discrete_distribution<size_t>(probs.begin(), probs.end())(rng);
		}
	private:
		std::mt19937_64 rng;
		size_t initial_seed_;
};

#endif