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

#ifndef CHSTABILIZER_HPP
#define CHSTABILIZER_HPP

#include "core.hpp"

class StabilizerState {
	public:
	    StabilizerState(const unsigned n_qubits);

		uint_t NQubits() const;
		scalar_t Omega() const;
		uint_t Gamma1() const;
		uint_t Gamma2() const;
		uint_t vVector() const;
    	const std::array<uint_t, MAX_N_SIZE> FMatrix() const;
    	const std::array<uint_t, MAX_N_SIZE> GMatrix() const;
    	const std::array<uint_t, MAX_N_SIZE> MMatrix() const;

		void CX(unsigned q, unsigned r);
		void CZ(unsigned q, unsigned r);
		void H(unsigned q);
		void S(unsigned q);
		void Z(unsigned q);
		void X(unsigned q);
		void Y(unsigned q);
		void Sdag(unsigned q);

		void CompBasisVector(uint_t x);
		void HadamardBasisVector(uint_t x);

		inline scalar_t ScalarPart() { return omega; }

		scalar_t Amplitude(uint_t x);

		scalar_t InnerProduct(StabilizerState &rhs);

		void MeasurePauli(const pauli_t P);
	
	private:
	    unsigned n;

		uint_t gamma1;
		uint_t gamma2;
    	std::array<uint_t, MAX_N_SIZE> F;
    	std::array<uint_t, MAX_N_SIZE> G;
		std::array<uint_t, MAX_N_SIZE> M;
		uint_t v;
		uint_t s;
		scalar_t omega;

		void RightCX(unsigned q, unsigned r);
		void RightCZ(unsigned q, unsigned r);
		void RightS(unsigned q);
		void RightZ(unsigned q);
		void RightSdag(unsigned q);

		pauli_t GetPauliX(uint_t x);

		void UpdateSvector(uint_t t, uint_t u, unsigned b);

    	std::array<uint_t, MAX_N_SIZE> FT;
		std::array<uint_t, MAX_N_SIZE> MT;

		void TransposeF();
		void TransposeM();

		bool isReadyFT;
		bool isReadyMT;
};

using stabilizer_t = StabilizerState;

#endif