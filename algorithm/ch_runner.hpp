/*
 * This file contains code derived from Qiskit.
 *
 * (C) Copyright IBM 2018, 2019.
 * Modifications Copyright 2026 Yunseo Hwang and Giwon Song.
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * This code is licensed under the Apache License, Version 2.0. You may
 * obtain a copy of this license in the LICENSE.txt file in the root directory
 * of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Any modifications or derivative works of this code must retain this
 * copyright notice, and modified files need to carry a notice indicating
 * that they have been altered from the originals.
 */

#ifndef CH_RUNNER_HPP
#define CH_RUNNER_HPP

#include "chlib/chstabilizer.hpp"
#include "chlib/core.hpp"
#include "tool/rng.hpp"
#include "tool/timer.hpp"
#include "gates.hpp"
#include <unordered_map>
#include <shared_mutex>

inline const uint_t ZERO = 0ULL;

class Runner {
	public:
		uint_t n_qubits_;
		uint_t num_states_;
		stabilizer_t states_[2];
		complex_t coefficients_[2];

	//public:
		Runner(uint_t n_qubits, stabilizer_t state);

		void reset(const stabilizer_t& initial_state);

		uint_t get_num_states() const;
		uint_t get_n_qubits() const;
		
		void apply_cx(uint_t control, uint_t target, uint_t rank);
		void apply_cz(uint_t control, uint_t target, uint_t rank);
		void apply_swap(uint_t qubit_1, uint_t qubit_2, uint_t rank);
		void apply_h(uint_t qubit, uint_t rank);
		void apply_sx(uint_t qubit, uint_t rank);
		void apply_sxdg(uint_t qubit, uint_t rank);
		void apply_s(uint_t qubit, uint_t rank);
		void apply_sdag(uint_t qubit, uint_t rank);
		void apply_x(uint_t qubit, uint_t rank);
		void apply_y(uint_t qubit, uint_t rank);
		void apply_z(uint_t qubit, uint_t rank);

		void apply_t(uint_t qubit, double r, int rank);
		void apply_tdag(uint_t qubit, double r, int rank);
		void apply_u1(uint_t qubit, complex_t lambda, double r, int rank);
		void apply_ccx(uint_t control_1, uint_t control_2, uint_t target,
					   uint_t branch, int rank);
		void apply_ccz(uint_t control_1, uint_t control_2, uint_t target,
					   uint_t branch, int rank);

		void apply_pauli(pauli_t &P);

		complex_t inner_product(uint_t rank1, uint_t rank2);
		complex_t projector_overlap(stabilizer_t projector_state);
		complex_t amplitude(uint_t x_measure);
};

#endif
