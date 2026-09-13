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

#include "extended_stabilizer_state.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

double State::expval(const Op* ops, size_t op_count, Observable obv,
					 uint_t n_threads, uint_t n_shots, double delta, double epsilon, bool bernstein) {
	return expval_with_stats(ops, op_count, obv, n_threads, n_shots, delta, epsilon, bernstein).result;
}

ExpvalStats State::expval_with_stats(const Op* ops, size_t op_count, Observable obv,
									 uint_t n_threads, uint_t n_shots, double delta,
									 double epsilon, bool bernstein) {
	
	ExpvalStats stats;
	stats.delta = delta;
	stats.epsilon = epsilon;
	stats.requested_shots = n_shots;
	stats.n_threads = n_threads;
	stats.op_count = op_count;
	stats.bernstein_enabled = bernstein;

	const bool is_projector = obv.type == ObservableType::rank1_projector;
	stats.is_projector = is_projector;
	stats.observable_size = is_projector ? 1 : obv.size;
	stats.observable_norm = obv.norm;
	stats.observable_norm_is_exact = obv.norm_is_exact;

	stabilizer_t projector_state(n_qubits);
	if (is_projector) {
		projector_state = prepare_projector_state(obv);
	}

	std::pair<bool, size_t> stabilizer_opts = check_stabilizer_opt(ops, op_count);
	bool is_stabilizer = stabilizer_opts.first;
	stats.is_stabilizer = is_stabilizer;
	stats.xi = compute_xi(ops, op_count);

	if (is_stabilizer) {
		stats.final_shots = 0;
		stabilizer_t state = prepare_base_state(ops, op_count);
		Runner runner(n_qubits, state);
		if (is_projector) {
			stats.result = obv.projector_coeff * std::real(runner.projector_overlap(projector_state));
			return stats;
		}
		double ret = 0;
		FlatObservables flat_obv = pack_observables(obv);
		CpuRngEngine cpu_rng(0);
		for (uint_t o = 0; o < flat_obv.size; o++) {
			apply_gate(runner, 2*o, flat_obv, cpu_rng, 0);
			ret += flat_obv.coeffs_ptr[o] * std::real(runner.inner_product(0, 1));
			apply_gate(runner, 2*o + 1, flat_obv, cpu_rng, 0);
		}
		stats.result = ret;
		return stats;
	}

	//non-clifford case.
	double xi = stats.xi;
	size_t first_non_clifford = stabilizer_opts.second;
	stabilizer_t base_state(prepare_base_state(ops, first_non_clifford));
	ops += first_non_clifford;
	op_count -= first_non_clifford;

	if (obv.size > MAX_OBV_SIZE) {
		throw std::runtime_error("Observable size exceeds the maximum supported size.");
	}

	FlatOperators flat_ops = pack_operators(ops, op_count);
	FlatObservables flat_obv = {};
	if (!is_projector) {
		flat_obv = pack_observables(obv);
	}
	Runner local_runner(n_qubits, base_state);

	omp_set_num_threads(n_threads);
	if (n_shots == 0) {// defined by epsilon, delta
		n_shots = std::llrint(2 * pow(xi, 2) * pow(epsilon, -2) * log(2.0 / delta));
	}
	stats.hoeffding_shots_initial = n_shots;
	std::cout << "n_shots: " << n_shots << std::endl;
	CpuRngEngine cpu_rng;

	double result = 0;
	uint_t m_shots = 0;

	//Using Bernstein's inequality rather than Hoeffding's.
	if(bernstein){
		m_shots = sqrt(8 * log(4.0 / delta) * n_shots);
		stats.variance_shots = m_shots;
		std::cout << "Number of shots to estimate variance: " << m_shots << std::endl;

		WelfordStats local_stats;
		#pragma omp parallel for if (n_threads > 1) reduction(welford_reduction: local_stats) private (cpu_rng) firstprivate(local_runner)
		for (uint_t shot = 0; shot < m_shots; shot++) {
			local_runner.reset(base_state);
			if (is_projector) {
				local_stats.update(projector_expval_mean_kernel(local_runner, flat_ops, projector_state, obv.projector_coeff, cpu_rng));
			} else {
				local_stats.update(expval_mean_kernel(local_runner, flat_ops, flat_obv, cpu_rng));
			}
		}

		result += local_stats.mean * m_shots;
		double var_hat = local_stats.variance(); //r.v. is estimated in [-|O|, |O|], so the variance is bounded by |O|^2.
		double std_hat = xi * sqrt(var_hat) / obv.norm;

		double sigma_max = std::min(xi, std_hat + xi * sqrt(8 * log(2.0 / delta) / (m_shots - 1)));
		double var_max = pow(sigma_max, 2);
		uint_t new_n_shots = std::llrint((2 * var_max + (xi+1) * epsilon * 2 / 3.0) * pow(epsilon, -2) * log(4.0 / delta)); //for bernstein's inequality
		stats.variance_hat = var_hat;
		stats.std_hat = std_hat;
		stats.sigma_max = sigma_max;
		stats.variance_bound = var_max;
		std::cout << "estimated standard deviation: " << std_hat << ", xi: " << xi << ", sigma_max: " << sigma_max << std::endl;
		std::cout << "n_shots: " << n_shots << " -> " << new_n_shots << std::endl;
		n_shots = new_n_shots;

		if(n_shots < m_shots){ n_shots = m_shots; }
		std::cout << "additional shots for estimate result: " << n_shots - m_shots << std::endl;
	}
	// Variance Estimation end.

	stats.final_shots = n_shots;

	// Estimate result
 	#pragma omp parallel for if (n_threads > 1) reduction(+:result) private(cpu_rng) firstprivate(local_runner)
 	for (uint_t shot = m_shots; shot < n_shots; shot++) {
		local_runner.reset(base_state);
		if (is_projector) {
			result += projector_expval_mean_kernel(local_runner, flat_ops, projector_state, obv.projector_coeff, cpu_rng);
		} else {
			result += expval_mean_kernel(local_runner, flat_ops, flat_obv, cpu_rng);
		}
 	}

	stats.result = result / n_shots * xi;
	return stats;
}

double State::compute_xi(const Op* ops, size_t op_count) const {
	double xi = 1;
	for (size_t i = 0; i < op_count; i++) {
		compute_extent(ops[i], xi);
	}
	return xi;
}

void State::compute_extent(const Op &op, double &xi) const {
	if (op.type == OpType::gate) {
		Gates gate_type = get_gate_enum(op.name);
		switch (gate_type) {
			case Gates::t:
				xi *= t_extent;
				break;
			case Gates::tdg:
				xi *= t_extent;
				break;
			case Gates::ccx:
				xi *= ccx_extent;
				break;
			case Gates::ccz:
				xi *= ccx_extent;
				break;
			case Gates::u1:
				xi *= u1_extent(std::real(op.params));
				break;
			default:
				break;
		}
	}
}

std::pair<bool, size_t> State::check_stabilizer_opt(const Op* ops, size_t op_count) const {
    for (size_t i = 0; i < op_count; ++i) {
        const Op& op = ops[i];

        if (op.type != OpType::gate) {
            continue;
        }
		Gatetypes type = get_gate_type(op.name);

        if (type == Gatetypes::non_clifford) {
            return std::pair<bool, size_t>({false, i});
        }
    }
    return std::pair<bool, size_t>({true, 0});
}

stabilizer_t State::prepare_base_state(const Op* ops, size_t op_count) {
	stabilizer_t base_state(n_qubits);
	for (size_t i = 0; i < op_count; i++) {
		int_t pi2;
    	Gates gate_type = get_gate_enum(ops[i].name);
		switch (gate_type) {
			case Gates::x:
				base_state.X(ops[i].qubits[0]);
				break;
			case Gates::y:
				base_state.Y(ops[i].qubits[0]);
				break;
			case Gates::z:
				base_state.Z(ops[i].qubits[0]);
				break;
			case Gates::s:
				base_state.S(ops[i].qubits[0]);
				break;
			case Gates::sdg:
				base_state.Sdag(ops[i].qubits[0]);
				break;
			case Gates::h:
				base_state.H(ops[i].qubits[0]);
				break;
			case Gates::sx:
				base_state.H(ops[i].qubits[0]);
				base_state.S(ops[i].qubits[0]);
				base_state.H(ops[i].qubits[0]);
				break;
			case Gates::sxdg:
				base_state.H(ops[i].qubits[0]);
				base_state.Sdag(ops[i].qubits[0]);
				base_state.H(ops[i].qubits[0]);
				break;
			case Gates::cx:
				base_state.CX(ops[i].qubits[0], ops[i].qubits[1]);
				break;
			case Gates::cz:
				base_state.CZ(ops[i].qubits[0], ops[i].qubits[1]);
				break;
			case Gates::swap:
				base_state.CX(ops[i].qubits[0], ops[i].qubits[1]);
				base_state.CX(ops[i].qubits[1], ops[i].qubits[0]);
				base_state.CX(ops[i].qubits[0], ops[i].qubits[1]);
				break;
			case Gates::pauli:
				base_apply_pauli(ops[i].qubit_size, ops[i].qubits, ops[i].string_params, base_state);
				break;
			case Gates::ecr:
				base_state.S(ops[i].qubits[0]);
				base_state.Sdag(ops[i].qubits[1]);
				base_state.H(ops[i].qubits[1]);
				base_state.Sdag(ops[i].qubits[1]);
				base_state.CX(ops[i].qubits[0], ops[i].qubits[1]);
				base_state.X(ops[i].qubits[0]);
				break;
			case Gates::rz:
				pi2 = (int_t)std::round(std::real(ops[i].params) * 2.0 / M_PI) & 3;
				if (pi2 == 1) {
					// S
					base_state.S(ops[i].qubits[0]);
				} else if (pi2 == 2) {
					// Z
					base_state.Z(ops[i].qubits[0]);
				} else if (pi2 == 3) {
					// Sdg
					base_state.Sdag(ops[i].qubits[0]);
				}
				break;
			default: // u0 or Identity
				break;
		}
	}
	return base_state;
}

stabilizer_t State::prepare_projector_state(const Observable& obv) {
	if (obv.type != ObservableType::rank1_projector) {
		throw std::logic_error("Observable is not a rank-1 projector.");
	}
	std::pair<bool, size_t> stabilizer_opts = check_stabilizer_opt(obv.projector_ops, obv.projector_op_count);
	if (!stabilizer_opts.first) {
		throw std::logic_error("Rank-1 projector state must be prepared by Clifford gates only.");
	}
	return prepare_base_state(obv.projector_ops, obv.projector_op_count);
}

void State::base_apply_pauli(uint_t qubit_size, const int* qubits, const char* pauli,
                             stabilizer_t &base_state) {
	for (uint_t i = 0; i < qubit_size; ++i) {
		const auto qubit = qubits[i]; //const auto qubit = qubits[qubit_size - 1 - i];
		switch (pauli[i]) {
		case 'I':
			break;
		case 'X':
			base_state.X(qubit);
			break;
		case 'Y':
			base_state.Y(qubit);
			break;
		case 'Z':
			base_state.Z(qubit);
			break;
		default:
			break;
		}
	}
}

template <typename RngEngine>
void expval_variance_kernel(WelfordStats* local_stats_array, Runner& local_runner, const FlatOperators& circuit_ops, const FlatObservables& obv, RngEngine& rng) {
	double value;
	for (uint_t i = 0; i < circuit_ops.num_ops; i++){
		apply_gate(local_runner, i, circuit_ops, rng, 0);
		apply_gate(local_runner, i, circuit_ops, rng, 1);
	}

	for (uint_t o = 0; o < obv.size; o++) { //use Welford's algorithm for memory efficiency.
		apply_gate(local_runner, 2*o, obv, rng, 0);
		value = std::real(local_runner.inner_product(0, 1));
		local_stats_array[o].update(value);
		apply_gate(local_runner, 2*o + 1, obv, rng, 0);
	}
}

template <typename RngEngine>
double expval_mean_kernel(Runner& local_runner, const FlatOperators& circuit_ops, const FlatObservables& obv, RngEngine& rng) {
	double value = 0;
    for (uint_t i = 0; i < circuit_ops.num_ops; i++) {
		apply_gate(local_runner, i, circuit_ops, rng, 0);
		apply_gate(local_runner, i, circuit_ops, rng, 1);
    }
    double ret = 0;
    for (uint_t o = 0; o < obv.size; o++) {
        apply_gate(local_runner, 2*o, obv, rng, 0);
		value = std::real(local_runner.inner_product(0, 1));
        ret += obv.coeffs_ptr[o] * value;
        apply_gate(local_runner, 2*o + 1, obv, rng, 0);
    }
    return ret;
}

template <typename RngEngine>
double projector_expval_mean_kernel(Runner& local_runner, const FlatOperators& circuit_ops,
									const stabilizer_t& projector_state, double projector_coeff, RngEngine& rng) {
	for (uint_t i = 0; i < circuit_ops.num_ops; i++) {
		apply_gate(local_runner, i, circuit_ops, rng, 0);
		apply_gate(local_runner, i, circuit_ops, rng, 1);
	}
	return projector_coeff * std::real(local_runner.projector_overlap(projector_state));
}

template <typename RngEngine, typename FlatOps>
void apply_gate(Runner &runner, uint_t idx, const FlatOps& circuit_ops, RngEngine& rng, uint_t rank) {
	int_t pi2;
	switch (circuit_ops.name_index[idx]) {
		case Gates::x:
			runner.apply_x(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rank);
			break;
		case Gates::y:
			runner.apply_y(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rank);
			break;
		case Gates::z:
			runner.apply_z(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rank);
			break;
		case Gates::s:
			runner.apply_s(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rank);
			break;
		case Gates::sdg:
			runner.apply_sdag(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rank);
			break;
		case Gates::h:
			runner.apply_h(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rank);
			break;
		case Gates::sx:
			runner.apply_sx(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rank);
			break;
		case Gates::sxdg:
			runner.apply_sxdg(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rank);
			break;
		case Gates::cx:
			runner.apply_cx(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]+1], rank);
			break;
		case Gates::cz:
			runner.apply_cz(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]+1], rank);
			break;
		case Gates::swap:
			runner.apply_swap(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx] + 1], rank);
			break;
		case Gates::t:
			runner.apply_t(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rng.rand(), rank);
			break;
		case Gates::tdg:
			runner.apply_tdag(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rng.rand(), rank);
			break;
		case Gates::ccx:
			runner.apply_ccx(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx] + 1], circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx] + 2],
										rng.rand_int(zero, toff_branch_max), rank);
			break;
		case Gates::ccz:
			runner.apply_ccz(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx] + 1], circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx] + 2],
										rng.rand_int(zero, toff_branch_max), rank);
			break;
		case Gates::u1:
			runner.apply_u1(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], circuit_ops.params[idx], rng.rand(), rank);
			break;
		case Gates::pauli:
			apply_pauli(runner, idx, circuit_ops, rank);
			break;
		case Gates::ecr:
			runner.apply_s(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rank);
			runner.apply_sdag(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx] + 1], rank);
			runner.apply_h(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx] + 1], rank);
			runner.apply_sdag(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx] + 1], rank);
			runner.apply_cx(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx] + 1], rank);
			runner.apply_x(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rank);
			break;
		case Gates::rz:
			pi2 = (int_t)std::round(std::real(circuit_ops.params[idx]) * 2.0 / M_PI) & 3;
			if (pi2 == 1) {
				// S
				runner.apply_s(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rank);
			} else if (pi2 == 2) {
				// Z
				runner.apply_z(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rank);
			} else if (pi2 == 3) {
				// Sdg
				runner.apply_sdag(circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx]], rank);
			}
			break;
		default: // u0 or Identity
			break;
	}
}

template <typename FlatOps>
void apply_pauli(Runner &runner, uint_t idx, const FlatOps &circuit_ops, uint_t rank) {
	const uint_t qubit_size = circuit_ops.qubit_size[idx];
	for (uint_t i = 0; i < qubit_size; i++) {
		const int qubit = circuit_ops.all_qubits_buffer[circuit_ops.qubits_offset[idx] + i];
		const char pauli = circuit_ops.all_strings_buffer[circuit_ops.string_params_offset[idx] + i];
		switch (pauli) {
		case 'I':
			break;
		case 'X':
			runner.apply_x(qubit, rank);
			break;
		case 'Y':
			runner.apply_y(qubit, rank);
			break;
		case 'Z':
			runner.apply_z(qubit, rank);
			break;
		case 'H': // Not Pauli, but can be added for observables
			runner.apply_h(qubit, rank);
			break;
		default:
			break;
		}
	}
}
