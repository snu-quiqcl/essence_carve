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

#ifndef OPERATORS_HPP
#define OPERATORS_HPP

#include "chlib/core.hpp"
#include "gates.hpp"

enum class OpType {
	gate,
	barrier,
};

enum class ObservableType {
	operator_sum,
	rank1_projector,
};

inline constexpr int MAX_OPS_SIZE = 1000; //we support up to 1000 operations
inline constexpr int MAX_OBV_SIZE = 100; //we support up to 100 observables

struct Op {
	Op();
	~Op();
	Op(const Op& other);

	Op& operator=(const Op& other);

	OpType type;
	const char* name;
	uint_t qubit_size;
	int* qubits;
	complex_t params; 
	char* string_params;
};

struct Observable {
    Observable();
	~Observable();

	Observable(const Observable& other);

	Observable& operator=(const Observable& other);

	ObservableType type;
	uint_t size;
	double norm;
	bool norm_is_exact;
	Op* operators;
	Op* inv_operators;
	double* coeffs;

	double projector_coeff;
	uint_t projector_op_count;
	Op* projector_ops;
};

//offset instead of pointer to qubits and string_params so that we can use it in GPU code.
struct FlatOperators {
    uint_t num_ops; // number of operations

    OpType type[MAX_OPS_SIZE];
    Gates name_index[MAX_OPS_SIZE]; 
    
    uint_t qubit_size[MAX_OPS_SIZE];
    size_t qubits_offset[MAX_OPS_SIZE]; 

    complex_t params[MAX_OPS_SIZE];
    size_t string_params_offset[MAX_OPS_SIZE];

    // all operations' data is stored in a single large buffer
    int all_qubits_buffer[MAX_OPS_SIZE*MAX_OBV_SIZE];
    char all_strings_buffer[MAX_OPS_SIZE*MAX_OBV_SIZE];
};


struct FlatObservables {
    uint_t size; // number of observables

	OpType type[2*MAX_OBV_SIZE];
    Gates name_index[2*MAX_OBV_SIZE]; 
    
    uint_t qubit_size[2*MAX_OBV_SIZE];
    size_t qubits_offset[2*MAX_OBV_SIZE]; 

    complex_t params[2*MAX_OBV_SIZE]; 
    size_t string_params_offset[2*MAX_OBV_SIZE];

    double coeffs_ptr[MAX_OBV_SIZE];

	// all qubits and strings are stored in a single buffer
    int all_qubits_buffer[MAX_OPS_SIZE*MAX_OBV_SIZE];
    char all_strings_buffer[MAX_OPS_SIZE*MAX_OBV_SIZE];
};

FlatOperators pack_operators(const Op* host_ops, size_t num_ops);
FlatObservables pack_observables(const Observable& host_obs);
#endif
