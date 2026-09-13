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

#include "operators.hpp"
#include <stdexcept>

Op::Op() :
    type(OpType::gate),
    name(nullptr),
    qubit_size(0),
    qubits(nullptr),
	params(0),
    string_params(nullptr) { }
    
Op::~Op(){
	delete[] qubits;
	delete[] string_params;
}

Op::Op(const Op& other)
: type(other.type),
  name(other.name), 
  qubit_size(other.qubit_size),
  qubits(nullptr), 
  params(other.params),
  string_params(nullptr) {
	if (other.qubits && qubit_size > 0) {
	    qubits = new int[qubit_size];
	    for (uint_t i = 0; i < qubit_size; ++i) {
	        qubits[i] = other.qubits[i];
	    }
	}
	if (other.string_params) {
	    size_t len = my_strlen(other.string_params);
	    string_params = new char[len + 1];
	    my_strcpy(string_params, other.string_params);
	}
}	

Op& Op::operator=(const Op& other) {
	if (this == &other) {
	    return *this;
	}
	delete[] qubits;
	delete[] string_params;
	type = other.type;
	name = other.name;
	qubit_size = other.qubit_size;
	params = other.params;
	qubits = nullptr;
	if (other.qubits && qubit_size > 0) {
	    qubits = new int[qubit_size];
	    for (uint_t i = 0; i < qubit_size; ++i) {
	        qubits[i] = other.qubits[i];
	    }
	}

	string_params = nullptr;
	if (other.string_params) {
	    size_t len = my_strlen(other.string_params);
	    string_params = new char[len + 1];
	    my_strcpy(string_params, other.string_params);
	}

	return *this;
}

Observable::Observable() :
    type(ObservableType::operator_sum),
    size(0),
    norm(0.0),
    norm_is_exact(false),
    operators(nullptr),
    inv_operators(nullptr),
    coeffs(nullptr),
    projector_coeff(1.0),
    projector_op_count(0),
    projector_ops(nullptr) { }

Observable::~Observable(){
	delete[] operators;
	delete[] inv_operators;
	delete[] coeffs;
	delete[] projector_ops;
}

Observable::Observable(const Observable& other)
: type(other.type),
  size(other.size),
    norm(other.norm),
    norm_is_exact(other.norm_is_exact),
  operators(nullptr),
  inv_operators(nullptr),
  coeffs(nullptr),
  projector_coeff(other.projector_coeff),
  projector_op_count(other.projector_op_count),
  projector_ops(nullptr) {
	if (type == ObservableType::operator_sum && size > 0) {
	    operators = new Op[size];
	    inv_operators = new Op[size];
	    coeffs = new double[size];
	
	    for (uint_t i = 0; i < size; ++i) {
	        operators[i] = other.operators[i];
	        inv_operators[i] = other.inv_operators[i];
	        coeffs[i] = other.coeffs[i];
	    }
	}
	if (projector_op_count > 0) {
	    projector_ops = new Op[projector_op_count];
	    for (uint_t i = 0; i < projector_op_count; ++i) {
	        projector_ops[i] = other.projector_ops[i];
	    }
	}
}

Observable& Observable::operator=(const Observable& other) {
    if (this == &other) {
        return *this;
    }

    delete[] operators;
    delete[] inv_operators;
    delete[] coeffs;
    delete[] projector_ops;

    type = other.type;
    size = other.size;
    norm = other.norm;
    norm_is_exact = other.norm_is_exact;
    operators = nullptr;
    inv_operators = nullptr;
    coeffs = nullptr;
    projector_coeff = other.projector_coeff;
    projector_op_count = other.projector_op_count;
    projector_ops = nullptr;

    if (type == ObservableType::operator_sum && size > 0) {
        operators = new Op[size];
        inv_operators = new Op[size];
        coeffs = new double[size];
	
        for (uint_t i = 0; i < size; ++i) {
            operators[i] = other.operators[i];
            inv_operators[i] = other.inv_operators[i];
            coeffs[i] = other.coeffs[i];
        }
    }
    if (projector_op_count > 0) {
        projector_ops = new Op[projector_op_count];
        for (uint_t i = 0; i < projector_op_count; ++i) {
            projector_ops[i] = other.projector_ops[i];
        }
    }
    return *this;
}

FlatOperators pack_operators(const Op* host_ops, size_t num_ops) {
    FlatOperators gpu_data = {};
    gpu_data.num_ops = num_ops;

    size_t current_qubit_offset = 0;
    size_t current_string_offset = 0;
    for (uint_t i = 0; i < gpu_data.num_ops; ++i) {
        const Op& op = host_ops[i];

        gpu_data.type[i] = op.type;
        gpu_data.name_index[i] = get_gate_enum(op.name);
        gpu_data.params[i] = op.params;

        gpu_data.qubit_size[i] = op.qubit_size;
        gpu_data.qubits_offset[i] = current_qubit_offset;
        if (op.qubit_size > 0) {
            std::copy(op.qubits, op.qubits + op.qubit_size,
                      gpu_data.all_qubits_buffer + current_qubit_offset);
            current_qubit_offset += op.qubit_size;
        }

        gpu_data.string_params_offset[i] = current_string_offset;
        if (op.string_params) {
            size_t len = my_strlen(op.string_params);
            my_strcpy(gpu_data.all_strings_buffer + current_string_offset, op.string_params);
            current_string_offset += len + 1;
        }
    }

    return gpu_data;
}

FlatObservables pack_observables(const Observable& host_obs) {
    FlatObservables gpu_obs = {};
    if (host_obs.type != ObservableType::operator_sum) {
        throw std::runtime_error("Only operator-sum observables can be packed.");
    }
    if (host_obs.size == 0) {
        return gpu_obs;
    }

    gpu_obs.size = host_obs.size;

    size_t qubit_offset = 0;
    size_t string_offset = 0;

    for (uint_t i = 0; i < gpu_obs.size; ++i) {
        gpu_obs.coeffs_ptr[i] = host_obs.coeffs[i];
        
        // operators
        const Op& op = host_obs.operators[i];

        gpu_obs.type[2*i] = op.type;
        gpu_obs.name_index[2*i] = get_gate_enum(op.name);
        gpu_obs.params[2*i] = op.params;

        gpu_obs.qubit_size[2*i] = op.qubit_size;
        gpu_obs.qubits_offset[2*i] = qubit_offset;
        if (op.qubit_size > 0) {
            std::copy(op.qubits, op.qubits + op.qubit_size, gpu_obs.all_qubits_buffer + qubit_offset);
            qubit_offset += op.qubit_size;
        }

        gpu_obs.string_params_offset[2*i] = string_offset;
        if (op.string_params) {
            size_t len = my_strlen(op.string_params);
            my_strcpy(gpu_obs.all_strings_buffer + string_offset, op.string_params);
            string_offset += len + 1;
        }

        // inv_operators
        const Op& inv_op = host_obs.inv_operators[i];
        gpu_obs.type[2*i + 1] = inv_op.type;
        gpu_obs.name_index[2*i + 1] = get_gate_enum(inv_op.name);
        gpu_obs.params[2*i + 1] = inv_op.params;

        gpu_obs.qubit_size[2*i + 1] = inv_op.qubit_size;
        gpu_obs.qubits_offset[2*i + 1] = qubit_offset;
        if (inv_op.qubit_size > 0) {
            std::copy(inv_op.qubits, inv_op.qubits + inv_op.qubit_size, gpu_obs.all_qubits_buffer + qubit_offset);
            qubit_offset += inv_op.qubit_size;
        }

        gpu_obs.string_params_offset[2*i + 1] = string_offset;
        if (inv_op.string_params) {
            size_t len = my_strlen(inv_op.string_params);
            my_strcpy(gpu_obs.all_strings_buffer + string_offset, inv_op.string_params);
            string_offset += len + 1;
        }
    }

    return gpu_obs;
}
