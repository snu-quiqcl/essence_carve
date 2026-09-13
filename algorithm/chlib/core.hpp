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

#ifndef CORE_HPP
#define CORE_HPP

#include <complex>
#include <cstdint>
#include <array>

using complex_t = std::complex<double>;
using uint_t = uint_fast64_t;
using int_t = int_fast64_t;

inline constexpr int MAX_N_SIZE = 64;
inline constexpr std::array<int, 8> RE_PHASE = {1, 1, 0, -1, -1, -1, 0, 1};
inline constexpr std::array<int, 8> IM_PHASE = {0, 1, 1, 1, 0, -1, -1, -1};
inline constexpr float FPE = 1e-8;
inline constexpr uint_t zer = 0U;
inline constexpr uint_t one = 1U;


int popcount(uint_t x);

bool hamming_parity(uint_t x);

bool almost_equal(double x, double y);

struct scalar_t {
	// complex numbers of the form eps * 2^{p/2} * exp(i (pi/4) * e )
	// eps = 0, 1 - if eps = 0 then p, e is arbitrary
	// p = integer
	// e = 0, 1, 2, 3, 4, 5, 6, 7
	int eps = 1;
	int p = 0;
	int e = 0;

	scalar_t() = default;
	scalar_t(const scalar_t& rhs) = default;
	scalar_t(const complex_t& coeff);

	scalar_t &operator*=(const scalar_t &rhs);
	scalar_t operator*(const scalar_t &rhs) const;

	complex_t to_complex() const;
	
	void conjugate();

	void makeOne();
};

struct pauli_t {
	// n-qubit Pauli operators: i^e X(x) Z(z)
	uint_t X; // n-bit string
	uint_t Z; // n-bit string
	unsigned e = 0;  // takes value 0, 1, 2, 3

	pauli_t();
	pauli_t(const pauli_t &p) = default;
	pauli_t(const char* s);

	pauli_t &operator*=(const pauli_t &rhs);
};

std::array<uint_t, MAX_N_SIZE> AtransposeB(const std::array<uint_t, MAX_N_SIZE>& A, const std::array<uint_t, MAX_N_SIZE>& B, unsigned n);

#endif

