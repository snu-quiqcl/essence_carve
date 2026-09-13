/*
 * Copyright 2026 Yunseo Hwang and Giwon Song.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "file_reader.hpp"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <optional>
#include "gates.hpp"

namespace {

bool ops_match(const Op& lhs, const Op& rhs) {
	if (lhs.type != rhs.type || lhs.name != rhs.name || lhs.qubit_size != rhs.qubit_size) {
		return false;
	}
	for (uint_t i = 0; i < lhs.qubit_size; i++) {
		if (lhs.qubits[i] != rhs.qubits[i]) {
			return false;
		}
	}
	if ((lhs.string_params == nullptr) != (rhs.string_params == nullptr)) {
		return false;
	}
	if (lhs.string_params != nullptr && my_strcmp(lhs.string_params, rhs.string_params) != 0) {
		return false;
	}
	return lhs.params == rhs.params;
}

bool can_compute_exact_norm(const Observable& obv) {
	if (obv.type == ObservableType::rank1_projector) {
		return true;
	}
	if (obv.size == 0) {
		return true;
	}
	for (uint_t i = 1; i < obv.size; i++) {
		if (!ops_match(obv.operators[0], obv.operators[i])) {
			return false;
		}
	}
	return true;
}

double compute_frobenius_lower_bound(const Observable& obv) {
	double coeff_sq_sum = 0.0;
	for (uint_t i = 0; i < obv.size; i++) {
		coeff_sq_sum += obv.coeffs[i] * obv.coeffs[i];
	}
	return std::sqrt(coeff_sq_sum);
}

void validate_observable_norm(const Observable& obv) {
	if (!std::isfinite(obv.norm) || obv.norm <= 0.0) {
		throw std::logic_error("Observable norm must be a positive finite value.");
	}
}

void finalize_observable_norm(Observable& obv, std::optional<double> provided_norm) {
	if (obv.type == ObservableType::rank1_projector) {
		const double projector_norm = std::abs(obv.projector_coeff);
		obv.norm = provided_norm.has_value() ? *provided_norm : projector_norm;
		validate_observable_norm(obv);
		obv.norm_is_exact = true;
		return;
	}

	if (provided_norm.has_value()) {
		obv.norm = *provided_norm;
		validate_observable_norm(obv);
		obv.norm_is_exact = true;
		return;
	}

	if (can_compute_exact_norm(obv)) {
		double coeff_sum = 0.0;
		for (uint_t i = 0; i < obv.size; i++) {
			coeff_sum += obv.coeffs[i];
		}
		obv.norm = std::abs(coeff_sum);
		validate_observable_norm(obv);
		obv.norm_is_exact = true;
		return;
	}

	obv.norm = compute_frobenius_lower_bound(obv);
	validate_observable_norm(obv);
	obv.norm_is_exact = false;
}

}

Op CircuitFile::readline(std::fstream& file) {
	std::string name;
	file >> name;
	return readline(file, name);
}

Op CircuitFile::readline(std::fstream& file, const std::string& name) {
	Op op;
	op.type = OpType::gate;

	for (int i = 0; i < NUM_GATE_TYPES; i++){
		if(name == GATE_NAMES[i]){
			op.name = GATE_NAMES[i];
			break;
		}
	}
	if (op.name == nullptr) {
		throw std::logic_error("Unknown gate in input file: " + name);
	}

	if (name == "pauli") {  //ex. pauli ZZZ 3 2 3 8 : op.name = "pauli", op.string_param = "ZZZ", op.qubit_size = 3, op.qubits = {2, 3, 8}
		int N;
		std::string param;
		file >> param >> N;
		op.string_params = new char[param.length() + 1];
		my_strcpy(op.string_params, param.c_str());
		op.qubit_size = N;
		op.qubits = new int[N];
		for (int i = 0; i < N; i++) {
			uint_t q;
			file >> q;
			op.qubits[i] = q;
		}
	}
	else {
		if (name == "cx" || name == "cz" || name == "swap") { // two-qubit gate //ex. cx 9 5 : op.name = "cx", op.qubit_size = 2, op.qubits = {9, 5} 
			int q1, q2;
			file >> q1 >> q2;
			if (q1 == q2) {
				throw std::logic_error("Invalid two-qubit gate specification.");
			}
			op.qubit_size = 2;
			op.qubits = new int[2] {q1, q2};
		}
		else if (name == "rz" || name == "u1") { // one-qubit gate with param  //ex. u1 8 0.157 : op.name = "u1", op.qubits_size = 1, op.qubits = {8}, op.params = {0.157}
			int q;
			complex_t param;
			file >> q >> param;
			op.qubit_size = 1;
			op.qubits = new int[1] {q};
			op.params = param;
		}
		else { // one-qubit gate //ex. h 5 : op.name = "h", op.qubit_size = 1, op.qubits = {5}
			int q;
			file >> q;
			op.qubit_size = 1;
			op.qubits = new int[1] {q};
		}
	}
	return op;
}

void CircuitFile::read_projector(std::fstream& file, Observable& obv, double coeff) {
	uint_t op_count;
	file >> op_count;
	obv.type = ObservableType::rank1_projector;
	obv.size = 1;
	obv.projector_coeff = coeff;
	obv.projector_op_count = op_count;
	obv.projector_ops = new Op[op_count];
	for (uint_t i = 0; i < op_count; i++) {
		obv.projector_ops[i] = readline(file);
	}
}

Observable CircuitFile::get_observable() {
	std::fstream file;
	file.open(filepath + "/observable.txt");
	if (!file.is_open()) {
		throw std::logic_error("Cannot open observable file: " + filepath + "/observable.txt");
	}

	std::string first_token;
	file >> first_token;
	std::optional<double> provided_norm;
	while (first_token == "norm" || first_token == "operator_norm") {
		double value;
		file >> value;
		provided_norm = value;
		file >> first_token;
	}
	Observable obv;
	if (first_token == "projector" || first_token == "rank1_projector") {
		read_projector(file, obv, 1.0);
		finalize_observable_norm(obv, provided_norm);
		return obv;
	}

	uint_t Nop = std::stoull(first_token);
	obv.size = Nop;

	obv.operators = new Op[Nop];
	obv.inv_operators = new Op[Nop];
	obv.coeffs = new double[Nop];

	double coeff;
	std::string name;
	for (uint_t i = 0; i < Nop; i++) {
		file >> coeff >> name;
		if (name == "projector" || name == "rank1_projector") {
			if (Nop != 1) {
				throw std::logic_error("Rank-1 projector observable cannot be mixed with other observable terms.");
			}
			delete[] obv.operators;
			delete[] obv.inv_operators;
			delete[] obv.coeffs;
			obv.operators = nullptr;
			obv.inv_operators = nullptr;
			obv.coeffs = nullptr;
			read_projector(file, obv, coeff);
			return obv;
		}
		obv.coeffs[i] = coeff;
		Op op = readline(file, name);
		obv.operators[i] = op;
		std::string op_name = op.name;
		if (op_name != "s" && op_name != "sdg")
			obv.inv_operators[i] = op;
		else if (op_name == "s") {
			op.name = "sdg";
			obv.inv_operators[i] = op;
		} else {
			op.name = "s";
			obv.inv_operators[i] = op;
		}
	}
	finalize_observable_norm(obv, provided_norm);
	return obv;
}

circuit_list_t CircuitFile::get_circuits() {
	std::fstream file;
	file.open(filepath + "/info.txt");
	int num_circuits;
	file >> num_circuits;
	file.close();
	circuit_list_t circuits;
	for (int index = 0; index < num_circuits; index++) {
		file.open(filepath + "/circuits/" + std::to_string(index) + ".txt");
		uint_t N, G;
		file >> N >> G;
		std::vector<Op> ops;
		for (uint_t i = 0; i < G; i++)
			ops.push_back(readline(file));
		circuits.push_back({N, ops});
		file.close();
	}
	return circuits;
}
