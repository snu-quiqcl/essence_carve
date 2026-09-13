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

#include "gates.hpp"

U1Sample::U1Sample(complex_t param) {
	double lambda = std::real(param);
	lambda += M_PI;
	uint_t shift_factor = std::floor(lambda / (2 * M_PI));
	lambda -= shift_factor * (2 * M_PI) + M_PI;
	double angle = std::abs(lambda);
	bool s_z_quadrant = (angle > M_PI / 2);
	if (s_z_quadrant) {
		angle = angle - M_PI / 2;
	}
	angle /= 2;
	complex_t coeff_0 = std::cos(angle) - std::sin(angle);
	complex_t coeff_1 = root2 * std::sin(angle);
	complex_t phase_0, phase_1;
	std::array<Gates, 2> gates;
	if (lambda < 0) {
		coeff_0 *= root_omega_star;
		coeff_1 = coeff_1 * root_omega;
		if (s_z_quadrant) {
			gates[0] = Gates::sdg;
			gates[1] = Gates::z;
		} else {
			gates[0] = Gates::id;
			gates[1] = Gates::sdg;
		}
	} else {
		coeff_0 *= root_omega;
		coeff_1 = coeff_1 * root_omega_star;
		if (s_z_quadrant) {
			gates[0] = Gates::s;
			gates[1] = Gates::z;
		} else {
			gates[0] = Gates::id;
			gates[1] = Gates::s;
		}
	}
	phase_0 = std::polar(1.0, std::arg(coeff_0));
	phase_1 = std::polar(1.0, std::arg(coeff_1));
	branches = {sample_branch_t(phase_0, gates[0]),
				sample_branch_t(phase_1, gates[1])};
	p_threshold = std::abs(coeff_0) / (std::abs(coeff_0) + std::abs(coeff_1));
}

U1Sample::U1Sample(const U1Sample &other) : Sample(other) {
	p_threshold = other.p_threshold;
}

sample_branch_t U1Sample::sample(double r) const {
	if (r < p_threshold) {
		return branches[0];
	}
	else {
		return branches[1];
	}
}

double u1_extent(double lambda) {
	// Shift parameter into +- 2 Pi
	uint_t shift_factor = std::floor(std::abs(lambda) / (2 * M_PI));
	if (shift_factor != 0) {
		if (lambda < 0) {
			lambda += shift_factor * (2 * M_PI);
		} else {
			lambda -= shift_factor * (2 * M_PI);
		}
	}
	// Shift parameter into +- Pi
	if (lambda > M_PI) {
		lambda -= 2 * M_PI;
	} else if (lambda < -1 * M_PI) {
		lambda += 2 * M_PI;
	}
	double angle = std::abs(lambda);
	bool s_z_quadrant = (angle > M_PI / 2);
	if (s_z_quadrant) {
		angle = angle - M_PI / 2;
	}
	angle /= 2;
	return my_pow(std::cos(angle) + tan_pi_over_8 * std::sin(angle), 2.);
}