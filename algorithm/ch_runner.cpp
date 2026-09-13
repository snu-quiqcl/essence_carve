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

#include "ch_runner.hpp"
#include <unordered_map>

Runner::Runner(uint_t n_qubits, stabilizer_t state) : n_qubits_(n_qubits),
    num_states_(2), states_{stabilizer_t(state), stabilizer_t(state)}, 
    coefficients_{complex_t(1., 0.), complex_t(1., 0.)} { }

void Runner::reset(const stabilizer_t& initial_state){
	states_[0] = initial_state;
	states_[1] = initial_state;
	coefficients_[0] = complex_t(1., 0.);
	coefficients_[1] = complex_t(1., 0.);
}

uint_t Runner::get_num_states() const { return num_states_; }

uint_t Runner::get_n_qubits() const { return n_qubits_; }

void Runner::apply_pauli(pauli_t &P) { 
	const int_t END = num_states_;
	for (int_t i = 0; i < END; i++) {
		states_[i].MeasurePauli(P);
	}
}

void Runner::apply_cx(uint_t control, uint_t target, uint_t rank) {
    states_[rank].CX(control, target);
}

void Runner::apply_cz(uint_t control, uint_t target, uint_t rank) {
    states_[rank].CZ(control, target);
}

void Runner::apply_swap(uint_t qubit_1, uint_t qubit_2, uint_t rank) {
    states_[rank].CX(qubit_1, qubit_2);
	states_[rank].CX(qubit_2, qubit_1);
	states_[rank].CX(qubit_1, qubit_2);
}

void Runner::apply_h(uint_t qubit, uint_t rank) {
    states_[rank].H(qubit);
}

void Runner::apply_s(uint_t qubit, uint_t rank) {
    states_[rank].S(qubit);
}

void Runner::apply_sdag(uint_t qubit, uint_t rank) {
    states_[rank].Sdag(qubit);
}

void Runner::apply_sx(uint_t qubit, uint_t rank) {
	states_[rank].H(qubit);
    states_[rank].S(qubit);
    states_[rank].H(qubit);
}

void Runner::apply_sxdg(uint_t qubit, uint_t rank) {
    states_[rank].H(qubit);
	states_[rank].Sdag(qubit);
	states_[rank].H(qubit);
}

void Runner::apply_x(uint_t qubit, uint_t rank) {
    states_[rank].X(qubit);
}

void Runner::apply_y(uint_t qubit, uint_t rank) {
    states_[rank].Y(qubit);
}

void Runner::apply_z(uint_t qubit, uint_t rank) {
    states_[rank].Z(qubit);
}

void Runner::apply_t(uint_t qubit, double r, int rank) {
	thread_local U1Sample t_sample(M_PI/4.);
	sample_branch_t branch = t_sample.sample(r);
    coefficients_[rank] *= branch.first;
    if (branch.second == Gates::s) {
        states_[rank].S(qubit);
    }
}

void Runner::apply_tdag(uint_t qubit, double r, int rank) {
	thread_local U1Sample tdg_sample(-1. * M_PI / 4.);
	sample_branch_t branch = tdg_sample.sample(r);
    coefficients_[rank] *= branch.first;
    if (branch.second == Gates::sdg) {
        states_[rank].Sdag(qubit);
    }
}

void Runner::apply_u1(uint_t qubit, complex_t lambda, double r, int rank) {
	sample_branch_t branch;
	U1Sample rotation(lambda); 
	branch = rotation.sample(r);

	coefficients_[rank] *= branch.first;
	switch (branch.second) {
		case Gates::s:
			states_[rank].S(qubit);
			break;
		case Gates::sdg:
			states_[rank].Sdag(qubit);
			break;
		case Gates::z:
			states_[rank].Z(qubit);
			break;
		default:
			break;
	}
}

void Runner::apply_ccx(uint_t control_1, uint_t control_2, uint_t target,
                       uint_t branch, int rank) {
	switch (branch) // Decomposition of the CCX gate into Cliffords
	{
		case 1:
			states_[rank].CZ(control_1, control_2);
			break;
		case 2:
			states_[rank].CX(control_1, target);
			break;
		case 3:
			states_[rank].CX(control_2, target);
			break;
		case 4:
			states_[rank].CZ(control_1, control_2);
			states_[rank].CX(control_1, target);
			states_[rank].Z(control_1);
			break;
		case 5:
			states_[rank].CZ(control_1, control_2);
			states_[rank].CX(control_2, target);
			states_[rank].Z(control_2);
			break;
		case 6:
			states_[rank].CX(control_1, target);
			states_[rank].CX(control_2, target);
			states_[rank].X(target);
			break;
		case 7:
			states_[rank].CZ(control_1, control_2);
			states_[rank].CX(control_1, target);
			states_[rank].CX(control_2, target);
			states_[rank].Z(control_1);
			states_[rank].Z(control_2);
			states_[rank].X(target);
			coefficients_[rank] *= -1; // Additional phase
			break;
		default: // Identity
			break;
	}
}

void Runner::apply_ccz(uint_t control_1, uint_t control_2, uint_t target,
                       uint_t branch, int rank) {
	switch (branch) // Decomposition of the CCZ gate into Cliffords
	{
		case 1:
			states_[rank].CZ(control_1, control_2);
			break;
		case 2:
			states_[rank].CZ(control_1, target);
			break;
		case 3:
			states_[rank].CZ(control_2, target);
			break;
		case 4:
			states_[rank].CZ(control_1, control_2);
			states_[rank].CZ(control_1, target);
			states_[rank].Z(control_1);
			break;
		case 5:
			states_[rank].CZ(control_1, control_2);
			states_[rank].CZ(control_2, target);
			states_[rank].Z(control_2);
			break;
		case 6:
			states_[rank].CZ(control_1, target);
			states_[rank].CZ(control_2, target);
			states_[rank].Z(target);
			break;
		case 7:
			states_[rank].CZ(control_1, control_2);
			states_[rank].CZ(control_1, target);
			states_[rank].CZ(control_2, target);
			states_[rank].Z(control_1);
			states_[rank].Z(control_2);
			states_[rank].Z(target);
			coefficients_[rank] *= -1; // Additional phase
			break;
		default: // Identity
			break;
	}
}

complex_t Runner::inner_product(uint_t rank1, uint_t rank2) {
	return std::conj(coefficients_[rank1]) * coefficients_[rank2] * states_[rank1].InnerProduct(states_[rank2]).to_complex();  
}

complex_t Runner::projector_overlap(stabilizer_t projector_state) {
	stabilizer_t left_projector(projector_state);
	complex_t left = std::conj(coefficients_[0]) * states_[0].InnerProduct(left_projector).to_complex();
	complex_t right = coefficients_[1] * projector_state.InnerProduct(states_[1]).to_complex();
	return left * right;
}
