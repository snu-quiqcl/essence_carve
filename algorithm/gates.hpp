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

#ifndef GATES_HPP
#define GATES_HPP

#include <cstdint>
#include <complex>
#include <cmath>
#include <array>

#include "tool/mystr.hpp"

using complex_t = std::complex<double>;
using uint_t = uint_fast64_t;

enum class Gates {
  u0,
  u1,
  id,
  x,
  y,
  z,
  h,
  s,
  sdg,
  sx,
  sxdg,
  t,
  tdg,
  cx,
  cz,
  swap,
  ccx,
  ccz,
  pauli,
  ecr,
  rz,
};

enum class Gatetypes { pauli, clifford, non_clifford };

inline constexpr const char* GATE_NAMES[] = {
    "u0", "delay", "id", "x", "y", "z", "s", "sdg", "sx", "sxdg",
    "h", "t", "tdg", "u1", "p", "rz", "CX", "cx", "cz", "swap",
    "ecr", "ccx", "ccz", "pauli"
};

inline constexpr Gatetypes GATE_TYPES[] = {
    Gatetypes::pauli, Gatetypes::pauli, Gatetypes::pauli, Gatetypes::pauli, Gatetypes::pauli,
    Gatetypes::pauli, Gatetypes::clifford, Gatetypes::clifford, Gatetypes::clifford, Gatetypes::clifford,
    Gatetypes::clifford, Gatetypes::non_clifford, Gatetypes::non_clifford, Gatetypes::non_clifford, Gatetypes::non_clifford,
    Gatetypes::clifford, Gatetypes::clifford, Gatetypes::clifford, Gatetypes::clifford, Gatetypes::clifford,
    Gatetypes::clifford, Gatetypes::non_clifford, Gatetypes::non_clifford, Gatetypes::pauli
};

inline constexpr Gates GATE_ENUMS[] = {
    Gates::u0, Gates::id, Gates::id, Gates::x, Gates::y, Gates::z, Gates::s, Gates::sdg, Gates::sx, Gates::sxdg,
    Gates::h, Gates::t, Gates::tdg, Gates::u1, Gates::u1, Gates::rz, Gates::cx, Gates::cx, Gates::cz, Gates::swap,
    Gates::ecr, Gates::ccx, Gates::ccz, Gates::pauli
};

inline constexpr int NUM_GATE_TYPES = 24;

inline Gatetypes get_gate_type(const char* gate_name) {
    for (int i = 0; i < NUM_GATE_TYPES; ++i) {
        if (my_strcmp(gate_name, GATE_NAMES[i]) == 0) {
            return GATE_TYPES[i];
        }
    }
    return Gatetypes::pauli; 
}

inline Gates get_gate_enum(const char* gate_name) {
    for (int i = 0; i < NUM_GATE_TYPES; ++i) {
        if (my_strcmp(gate_name, GATE_NAMES[i]) == 0) {
            return GATE_ENUMS[i];
        }
    }
    return Gates::id;
}

using sample_branch_t = std::pair<complex_t, Gates>;

inline constexpr double root2 = 1.4142135623730951; // std::sqrt(2);
inline constexpr double root1_2 = 0.70710678118654752; // 1./sqrt(2);
inline constexpr complex_t pi_over_8_phase(0., 3.1415926535897932 / 8); // complex_t(0., M_PI/8);
inline constexpr complex_t omega(root1_2, root1_2); 
inline constexpr complex_t omega_star(root1_2, -1 * root1_2);
inline constexpr complex_t root_omega(0.92387953251128676, 0.38268343236508977); // std::exp(pi_over_8_phase);
inline constexpr complex_t root_omega_star(0.92387953251128676, -0.38268343236508977); // std::conj(root_omega);
inline constexpr double tan_pi_over_8 = 0.41421356237309505; // std::tan(M_PI / 8.);

// Base class for sampling over non-Clifford gates in the Sum Over Cliffords
// routine.
struct Sample {
public:
	Sample() = default;
	std::array<sample_branch_t, 2> branches;
};

// Functor class that defines how to sample branches over a U1 operation
// Used for caching each U1 gate angle we encounter in the circuit
struct U1Sample : public Sample {
	double p_threshold;

	U1Sample() = default;
	U1Sample(complex_t param);
	U1Sample(const U1Sample &other);
	~U1Sample() = default;

	sample_branch_t sample(double r) const;
};

inline const double t_extent = std::pow(1. / (std::cos(M_PI / 8.)), 2);
// Equation 32 in arXiv: 1808.00128
inline const double ccx_extent = 16. / 9.;
inline const double ccx_coeff = 1. / 6.;
// General result for z rotations, Eq. 28 in arXiv 1809.00128
double u1_extent(double lambda);

#endif