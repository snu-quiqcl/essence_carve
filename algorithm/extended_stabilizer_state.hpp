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

#ifndef EXTENDED_STABILIZER_STATE_HPP
#define EXTENDED_STABILIZER_STATE_HPP

#include <random>
#include <string>
#include <omp.h>
#include "ch_runner.hpp"
#include "chlib/chstabilizer.hpp"
#include "chlib/core.hpp"
#include "gates.hpp"
#include "operators.hpp"
#include "tool/timer.hpp"
#include "tool/welford.hpp"

inline constexpr uint_t zero = 0ULL;
inline constexpr uint_t toff_branch_max = 7ULL;

struct ExpvalStats {
	double result = 0.0;
	double xi = 1.0;
	double delta = 0.0;
	double epsilon = 0.0;
	double variance_mean = 0.0;
	double variance_hat = 0.0;
	double std_hat = 0.0;
	double sigma_max = 0.0;
	double variance_bound = 0.0;
	double observable_norm = 0.0;
	bool observable_norm_is_exact = false;
	uint_t requested_shots = 0;
	uint_t hoeffding_shots_initial = 1;
	uint_t variance_shots = 0;
	uint_t bernstein_shots_new = 0;
	uint_t final_shots = 1;
	uint_t n_threads = 0;
	size_t op_count = 0;
	uint_t observable_size = 0;
	bool bernstein_enabled = false;
	bool is_stabilizer = false;
	bool is_projector = false;
};

class State {
	public:
		
		State(uint_t n_qubits): n_qubits(n_qubits) {}
	
		double compute_xi(const Op* first, size_t op_count) const;

		double expval(const Op* first, size_t op_count, Observable obv, uint_t n_threads = omp_get_thread_limit(),
					  uint_t n_shots = 0, double delta = 0.05, double epsilon = 0.1, bool bernstein = true);

		ExpvalStats expval_with_stats(const Op* first, size_t op_count, Observable obv,
									 uint_t n_threads = omp_get_thread_limit(),
									 uint_t n_shots = 0, double delta = 0.05,
									 double epsilon = 0.1, bool bernstein = true);

	protected:		
		void base_apply_pauli(uint_t qubit_size, const int* qubits, const char* pauli, stabilizer_t &base_state);

		// Add the given operation to the extent
		void compute_extent(const Op &op, double &xi) const;

		std::pair<bool, size_t> check_stabilizer_opt(const Op* first, size_t op_count) const;

		stabilizer_t prepare_base_state(const Op* first, size_t op_count);
		stabilizer_t prepare_projector_state(const Observable& obv);

	private:
		uint_t n_qubits;
		double global_phase = 0;
};

template <typename RngEngine>
void expval_variance_kernel(WelfordStats* local_stats_array, Runner& local_runner, const FlatOperators& circuit_ops, const FlatObservables& obv, RngEngine& rng);

template <typename RngEngine>
double expval_mean_kernel(Runner& local_runner, const FlatOperators& circuit_ops, const FlatObservables& obv, RngEngine& rng);

template <typename RngEngine>
double projector_expval_mean_kernel(Runner& local_runner, const FlatOperators& circuit_ops,
									const stabilizer_t& projector_state, double projector_coeff, RngEngine& rng);

template <typename RngEngine, typename FlatOps>
void apply_gate(Runner &runner, uint_t idx, const FlatOps& circuit_ops, RngEngine& rng, uint_t rank);

template <typename FlatOps>
void apply_pauli(Runner &runner, uint_t idx, const FlatOps &circuit_ops, uint_t rank);

#endif
