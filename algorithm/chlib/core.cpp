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

#include <array>
#include <climits>
#include <complex>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include "core.hpp"
#include "../tool/mystr.hpp"

#include <chrono>

int popcount(uint_t x) {
	int count = 0;
    while (x) {
        x &= (x - 1);
        count++;
    }
    return count;
}

bool hamming_parity(uint_t x) {

	x ^= x >> 32;
	x ^= x >> 16;
	x ^= x >> 8;
	x ^= x >> 4;
	x ^= x >> 2;
	x ^= x >> 1;

	return x & 1;
}

bool almost_equal(double x, double y) {
	return x - y < 1e-8 && y - x < 1e-8;
}

scalar_t::scalar_t(const std::complex<double>& coeff) {
	double abs_val = std::abs(coeff);
	if (std::abs(abs_val - 0.) < FPE) {
		eps = 0;
	} else {
		p = (int)std::log2(my_pow(abs_val, 2));
		bool is_real_zero = (std::abs(coeff.real() - 0) < 1e-8);
		bool is_imag_zero = (std::abs(coeff.imag() - 0) < 1e-8);
		bool is_real_posi = (coeff.real() > 0) && !(is_real_zero);
		bool is_imag_posi = (coeff.imag() > 0) && !(is_imag_zero);
		char switch_var = (is_real_zero * 1) + (is_imag_zero * 2) +
						  (is_real_posi * 4) + (is_imag_posi * 8);
		switch (switch_var) {
			case 0:
				e = 5;
				break;
			case 1:
				e = 6;
				break;
			case 2:
				e = 4;
				break;
			case 4:
				e = 7;
				break;
			case 6:
				e = 0;
				break;
			case 8:
				e = 3;
				break;
			case 9:
				e = 2;
				break;
			case 12:
				e = 1;
				break;
			default:
				//throw std::runtime_error("Unsure what to do here chief.");
				break;
		}
    }
}

scalar_t &scalar_t::operator*=(const scalar_t &rhs) {
	p += (rhs.p);
	e += (rhs.e);
	e %= 8;
	eps *= rhs.eps;
	return *this;
}

scalar_t scalar_t::operator*(const scalar_t &rhs) const {
	scalar_t out;
	out.p = p + rhs.p;
	out.e = (e + rhs.e) % 8;
	out.eps = eps * rhs.eps;
	return out;
}

complex_t scalar_t::to_complex() const {
	if (eps == 0) {
		return {0., 0.};
	}
	complex_t mag(my_pow(2, p / (double)2), 0.);
	complex_t phase(RE_PHASE[e], IM_PHASE[e]);
	if (e % 2) {
		phase /= std::sqrt(2);
	}
	return mag * phase;
}

void scalar_t::conjugate() {
	e %= 8;
    e = (8 - e) % 8;
}

void scalar_t::makeOne() {
	eps = 1;
	p = 0;
	e = 0;
}

pauli_t::pauli_t(): X(zer), Z(zer) {}

pauli_t &pauli_t::operator*=(const pauli_t &rhs) {
	unsigned overlap = popcount(Z & rhs.X);
	X ^= rhs.X;
	Z ^= rhs.Z;
	e = (e + rhs.e + 2 * overlap) % 4;
	return *this;
}

pauli_t::pauli_t(const char* s) {
    for (int i = 0; s[i] != '\0'; ++i) {
        switch (s[i]) {
            case 'I':
                break;
            case 'X':
                X += (1ULL << i);
                break;
            case 'Y':
                X += (1ULL << i);
                Z += (1ULL << i);
                e = (e + 1) % 4;
                break;
            case 'Z':
                Z += (1ULL << i);
                break;
            default:
                return; 
        }
    }
}

std::array<uint_t, MAX_N_SIZE> AtransposeB(const std::array<uint_t, MAX_N_SIZE>& A, const std::array<uint_t, MAX_N_SIZE>& B, unsigned n) {
    std::array<uint_t, MAX_N_SIZE> out{}; 
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            if (hamming_parity(A[i] & B[j])) {
                out[j] ^= (one << i);
            }
        }
    }
    return out;
}