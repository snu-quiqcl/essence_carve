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

#include "chstabilizer.hpp"
#include <array>

StabilizerState::StabilizerState(const unsigned n_qubits)
	: n(n_qubits),
	  gamma1(zer), gamma2(zer), v(zer), s(zer),
	  omega(), isReadyFT(false), isReadyMT(false) {

    for (unsigned q = 0; q < n; q++) {
        G[q] = (one << q);
        F[q] = (one << q);
        M[q] = zer;
        FT[q] = zer;
        MT[q] = zer;
    }
	
	omega.makeOne();
}

uint_t StabilizerState::NQubits() const { return n; }
scalar_t StabilizerState::Omega() const { return omega; }
uint_t StabilizerState::Gamma1() const { return gamma1; }
uint_t StabilizerState::Gamma2() const { return gamma2; }
uint_t StabilizerState::vVector() const { return v; }
const std::array<uint_t, MAX_N_SIZE> StabilizerState::FMatrix() const { return F; }
const std::array<uint_t, MAX_N_SIZE> StabilizerState::GMatrix() const { return G; }
const std::array<uint_t, MAX_N_SIZE> StabilizerState::MMatrix() const { return M; }

void StabilizerState::CompBasisVector(uint_t x) {
	s = x;
	v = zer;
	gamma1 = zer;
	omega.makeOne();
	for (unsigned q = 0; q < n; q++) {
		M[q] = zer;
		G[q] = (one << q);
		F[q] = (one << q);
	}
	isReadyFT = false;
	isReadyMT = false;
}

void StabilizerState::HadamardBasisVector(uint_t x) {
	s = zer;
	v = x;
	gamma1 = zer;
	omega.makeOne();
	for (unsigned q = 0; q < n; q++) {
		M[q] = zer;
		G[q] = (one << q);
		F[q] = (one << q);
	}
	isReadyFT = false;
	isReadyMT = false;
}

void StabilizerState::RightS(unsigned q) {
	isReadyMT = false;
	M[q] ^= F[q];
	gamma2 ^= F[q] ^ (gamma1 & F[q]);
	gamma1 ^= F[q];
}

void StabilizerState::RightSdag(unsigned q) {
	isReadyMT = false;
	M[q] ^= F[q];
	gamma1 ^= F[q];
	gamma2 ^= F[q] ^ (gamma1 & F[q]);
}

void StabilizerState::RightZ(unsigned q) {
	gamma2 ^= F[q];
}

void StabilizerState::S(unsigned q) {
	isReadyMT = false;
	uint_t C = (one << q);
	for (unsigned p = 0; p < n; p++) {
		M[p] ^= ((G[p] >> q) & one) * C;
	}
	gamma1 ^= C;
	gamma2 ^= ((gamma1 >> q) & one) * C;
}

void StabilizerState::Sdag(unsigned q) {
	isReadyMT = false;
	uint_t C = (one << q);
	for (unsigned p = 0; p < n; p++) {
		M[p] ^= ((G[p] >> q) & one) * C;
	}
	gamma2 ^= ((gamma1 >> q) & one) * C;
	gamma1 ^= C;
}

void StabilizerState::Z(unsigned q) {
	gamma2 ^= (one << q);
}

void StabilizerState::X(unsigned q) {
	if (!isReadyMT) {
		TransposeM();
	}
	if (!isReadyFT) {
		TransposeF();
	}
	uint_t x_string = FT[q];
	uint_t z_string = MT[q];
	int phase = 2 * ((gamma1 >> q) & one) + 4 * ((gamma2 >> q) & one);
	s ^= (z_string & v);
	phase += 4 * (hamming_parity((z_string & ~(v)) & s));
	s ^= (x_string & ~(v));
	phase += 4 * (hamming_parity((x_string & v) & s));
	omega.e = (omega.e + phase) % 8;
}

void StabilizerState::Y(unsigned q) {
	Z(q);
	X(q);
	omega.e = (omega.e + 2) % 8;
}

void StabilizerState::RightCX(unsigned q, unsigned r) {
	G[q] ^= G[r];
	F[r] ^= F[q];
	M[q] ^= M[r];
}

void StabilizerState::CX(unsigned q, unsigned r) {
	isReadyMT = false;
	isReadyFT = false;
	uint_t C = (one << q);
	uint_t T = (one << r);
	bool b = false;
	for (unsigned p = 0; p < n; p++) {
		b = (b != ((M[p] & C) && (F[p] & T)));
		G[p] ^= ((G[p] >> q) & one) * T;
		F[p] ^= ((F[p] >> r) & one) * C;
		M[p] ^= ((M[p] >> r) & one) * C;
	}
	if (b)
		gamma2 ^= C;
	b = ((gamma1 & C) && (gamma1 & T));
	gamma1 ^= ((gamma1 >> r) & one) * C;
	gamma2 ^= ((gamma2 >> r) & one) * C;
	if (b)
		gamma2 ^= C;
}

void StabilizerState::RightCZ(unsigned q, unsigned r) {
	isReadyMT = false;
	M[q] ^= F[r];
	M[r] ^= F[q];
	gamma2 ^= (F[q] & F[r]);
}

void StabilizerState::CZ(unsigned q, unsigned r) {
	isReadyMT = false;
	uint_t C = (one << q);
	uint_t T = (one << r);
	for (unsigned p = 0; p < n; p++) {
		M[p] ^= ((G[p] >> r) & one) * C;
		M[p] ^= ((G[p] >> q) & one) * T;
	}
}

pauli_t StabilizerState::GetPauliX(uint_t x) {
	if (!isReadyMT) {
		TransposeM();
	}
	if (!isReadyFT) {
		TransposeF();
	}
	pauli_t R;

	for (unsigned pos = 0; pos < n; pos++) {
		if (x & (one << pos)) {
			pauli_t P1;
			P1.e = 1 * ((gamma1 >> pos) & one);
			P1.e += 2 * ((gamma2 >> pos) & one);
			P1.X = FT[pos];
			P1.Z = MT[pos];
			R *= P1;
		}
	}

	return R;
}

scalar_t StabilizerState::Amplitude(uint_t x) {
	if (!isReadyMT) {
		TransposeM();
	}
	if (!isReadyFT) {
		TransposeF();
	}
	pauli_t P = GetPauliX(x);

	if (!omega.eps)
		return omega;

	scalar_t amp;
	amp.e = 2 * P.e;
	int p = popcount(v);
	amp.p = -1 * p;
	bool isNonZero = true;

	for (unsigned q = 0; q < n; q++) {
		uint_t pos = (one << q);
		if (v & pos) {
			amp.e += 4 * ((s & pos) && (P.X & pos));
		} else {
			isNonZero = ((P.X & pos) == (s & pos));
		}
		if (!isNonZero)
			break;
	}

	amp.e %= 8;
	if (isNonZero) {
		amp.conjugate();
	} else {
		amp.eps = 0;
		return amp;
	}

	amp.p += omega.p;
	amp.e = (amp.e + omega.e) % 8;
	return amp;
}

void StabilizerState::TransposeF() {
	for (unsigned p = 0; p < n; p++) {
		uint_t row = zer;
		uint_t mask = (one << p);
		for (unsigned q = 0; q < n; q++)
			if (F[q] & mask)
				row ^= (one << q);
		FT[p] = row;
	}
	isReadyFT = true;
}

void StabilizerState::TransposeM() {
	for (unsigned p = 0; p < n; p++) {
		uint_t row = zer;
		uint_t mask = (one << p);
		for (unsigned q = 0; q < n; q++)
			if (M[q] & mask)
				row ^= (one << q);
		MT[p] = row;
	}
	isReadyMT = true;
}

void StabilizerState::UpdateSvector(uint_t t, uint_t u, unsigned b) {
	if (t == u)
		switch (b) {
			case 0:
				omega.p += 1;
				s = t;
				return;
			case 1:
				s = t;
				omega.e = (omega.e + 1) % 8;
				return;
			case 2:
				s = t;
				omega.eps = 0;
				return;
			case 3:
				s = t;
				omega.e = (omega.e + 7) % 8;
				return;
			default:
				break;
		}

	isReadyFT = false;
	isReadyMT = false;
	uint_t ut = u ^ t;
	uint_t nu0 = (~v) & ut;
	uint_t nu1 = v & ut;

	b %= 4;
	unsigned q = 0;
	uint_t qpos = zer;
	if (nu0) {
		q = 0;
		while (!(nu0 & (one << q)))
			q++;

		qpos = (one << q);
		if (!(nu0 & qpos)) {
			//throw std::logic_error(
			//	"Failed to find first bit of nu despite being non-empty.");
			return;
		}

		nu0 ^= qpos; 
		if (nu0)
			for (unsigned q1 = q + 1; q1 < n; q1++)
				if (nu0 & (one << q1))
					RightCX(q, q1);
		if (nu1)
			for (unsigned q1 = 0; q1 < n; q1++)
				if (nu1 & (one << q1))
					RightCZ(q, q1);
	}
	else {
		q = 0;
		while (!(nu1 & (one << q)))
			q++;

		qpos = (one << q);
		if (!(nu1 & qpos)) {
			//throw std::logic_error("Expect at least one non-zero element in nu1.\n");
			return;
		}

		nu1 ^= qpos;
		if (nu1)
			for (unsigned q1 = q + 1; q1 < n; q1++)
				if (nu1 & (one << q1))
					RightCX(q1, q);

	}

	if (t & qpos) {
		s = u;
		omega.e = (omega.e + 2 * b) % 8;
		b = (4 - b) % 4;
		if (!!(u & qpos)) {
		}
	} else
		s = t;

	// change the order of H and S gates
	// H^{a} S^{b} |+> = eta^{e1} S^{e2} H^{e3} |e4>
	// here eta=exp( i (pi/4) )
	// a=0,1
	// b=0,1,2,3
	//
	// H^0 S^0 |+> = eta^0 S^0 H^1 |0>
	// H^0 S^1 |+> = eta^0 S^1 H^1 |0>
	// H^0 S^2 |+> = eta^0 S^0 H^1 |1>
	// H^0 S^3 |+> = eta^0 S^1 H^1 |1>
	//
	// H^1 S^0 |+> = eta^0 S^0 H^0 |0>
	// H^1 S^1 |+> = eta^1 S^1 H^1 |1>
	// H^1 S^2 |+> = eta^0 S^0 H^0 |1>
	// H^1 S^3 |+> = eta^{7} S^1 H^1 |0>
	//
	// "analytic" formula:
	// e1 = a * (b mod 2) * ( 3*b -2 )
	// e2 = b mod 2
	// e3 = not(a) + a * (b mod 2)
	// e4 = not(a)*(b>=2) + a*( (b==1) || (b==2) )

	bool a = ((v & qpos) > 0);
	unsigned e1 = a * (b % 2) * (3 * b - 2);
	unsigned e2 = b % 2;
	bool e3 = ((!a) != (a && ((b % 2) > 0)));
	bool e4 = (((!a) && (b >= 2)) != (a && ((b == 1) || (b == 2))));

	// update CH-form
	// set q-th bit of s to e4
	s &= ~qpos;
	s ^= e4 * qpos;

	// set q-th bit of v to e3
	v &= ~qpos;
	v ^= e3 * qpos;

	// update the scalar factor omega
	omega.e = (omega.e + e1) % 8;

	// multiply the C-layer on the right by S^{e2} on the q-th qubit
	if (e2)
		RightS(q);
}

void StabilizerState::H(unsigned q) {
	isReadyMT = false;
	isReadyFT = false;
	uint_t rowF = zer;
	uint_t rowG = zer;
	uint_t rowM = zer;
	for (unsigned j = 0; j < n; j++) {
		uint_t pos = (one << j);
		rowF ^= pos * ((F[j] >> q) & one);
		rowG ^= pos * ((G[j] >> q) & one);
		rowM ^= pos * ((M[j] >> q) & one);
	}

	// after commuting H through the C and H laters it maps |s> to a state
	// sqrt(0.5)*[  (-1)^alpha |t> + i^{gamma[p]} (-1)^beta |u>  ]
	//
	// compute t,s,alpha,beta
	uint_t t = s ^ (rowG & v);
	uint_t u = s ^ (rowF & (~v)) ^ (rowM & v);

	unsigned alpha = popcount(rowG & (~v) & s);
	unsigned beta = popcount((rowM & (~v) & s) ^ (rowF & v & (rowM ^ s)));

	if (alpha % 2)
		omega.e = (omega.e + 4) % 8;
	// get the phase gamma[q]
	unsigned phase = ((gamma1 >> q) & one) + 2 * ((gamma2 >> q) & one);
	unsigned b = (phase + 2 * alpha + 2 * beta) % 4;

	// now the initial state is sqrt(0.5)*(|t> + i^b |u>)

	// take care of the trivial case
	if (t == u) {
		s = t;

		//std::cout << "b: " << b << ", s: " << s << ", t: " << t << std::endl;

		if (!((b == 1) || (b == 3))) // otherwise the state is not normalized
		{
			//throw std::logic_error(
			//	"State is not properly normalised, b should be 1 or 3.\n");
			return;
		}
		if (b == 1)
			omega.e = (omega.e + 1) % 8;
		else
			omega.e = (omega.e + 7) % 8;
	} else
	UpdateSvector(t, u, b);
}

scalar_t StabilizerState::InnerProduct(StabilizerState &rhs) {
	// compute transposed matrices if needed
	if (!isReadyFT) {
		TransposeF();
	} 
	if (!isReadyMT) {
		TransposeM();
	} 
	if (!rhs.isReadyFT) {
		rhs.TransposeF();
	} 
	if (!rhs.isReadyMT) {
		rhs.TransposeM();
	} 
	// compute U_C1^{-1} U_C2
	StabilizerState temp = StabilizerState(n);
	// F = G_1^T F_2
	temp.F = AtransposeB(G, rhs.F, n);
	// G = F_1^T G_2
	temp.G = AtransposeB(F, rhs.G, n);
	// M = G_1^T M_2 + G_1^T M_1 F_1^T G_2
	std::array<uint_t, MAX_N_SIZE> M1 = AtransposeB(G, rhs.M, n);
	std::array<uint_t, MAX_N_SIZE> M2 = AtransposeB(G, AtransposeB(MT, temp.G, n), n);
	for (size_t i = 0; i < n; i++) {
		temp.M[i] = M1[i] ^ M2[i];
	}
	// gamma = G_1^T (gamma_2 - gamma_1) + 2S
	uint_t sub_gamma1 = gamma1 ^ rhs.gamma1;
	//uint_t sub_gamma2 = (gamma1 ^ gamma2) ^ rhs.gamma2 ^ (gamma1 & rhs.gamma1);
	uint_t sub_gamma2 = (gamma2 ^ rhs.gamma2) ^ ((~rhs.gamma1) & gamma1);
	for (size_t i = 0; i < n; i++) {
		int e = popcount(G[i] & (sub_gamma1));
		if (e & one)
			temp.gamma1 |= (one << i);
		if ((e >> 1) & one)
			temp.gamma2 |= (one << i);
		if (hamming_parity(G[i] & (sub_gamma2)))
			temp.gamma2 ^= (one << i);
	}

	uint_t S = zer;
	uint_t C1, C2;
	for (size_t p = 0; p < n; p++) {
		C1 = zer;
		C2 = zer;
		for (size_t q = 0; q < n; q++) {
			if ((G[p] >> q) & one) {
				if (hamming_parity(C1 & FT[q]) != hamming_parity(C2 & rhs.FT[q]))
					S ^= (one << p);
				C1 ^= MT[q];
				C2 ^= rhs.MT[q];
			}
		}
	}
	temp.gamma2 ^= S;

	temp.s = rhs.s;
	temp.v = rhs.v;
	omega.conjugate();
	temp.omega = omega * rhs.omega;
	omega.conjugate();

	for (size_t i = 0; i < n; i++) {
		if ((v >> i) & one){
			temp.H(i);
		}
	}
	return temp.Amplitude(s);
}

void StabilizerState::MeasurePauli(pauli_t PP) {
	// compute Pauli R = U_C^{-1} P U_C
	pauli_t R;
	R.e = PP.e;

	for (unsigned j = 0; j < n; j++)
		if ((PP.X >> j) & one) {
			// multiply R by U_C^{-1} X_j U_C
			// extract the j-th rows of F and M
			uint_fast64_t rowF = zer;
			uint_fast64_t rowM = zer;
			for (unsigned i = 0; i < n; i++) {
				rowF ^= (one << i) * ((F[i] >> j) & one);
				rowM ^= (one << i) * ((M[i] >> j) & one);
			}
			R.e += 2 * popcount(R.Z & rowF); // extra sign from Pauli commutation
			R.Z ^= rowM;
			R.X ^= rowF;
			R.e += ((gamma1 >> j) & one) + 2 * ((gamma2 >> j) & one);
		}
	for (unsigned q = 0; q < n; q++)
		R.Z ^= (one << q) * (popcount(PP.Z & G[q]) % 2);

	// now R=U_C^{-1} PP U_C
	// next conjugate R by U_H
	uint_fast64_t tempX = ((~v) & R.X) ^ (v & R.Z);
	uint_fast64_t tempZ = ((~v) & R.Z) ^ (v & R.X);
	// the sign flips each time a Hadamard hits Y on some qubit
	R.e = (R.e + 2 * popcount(v & R.X & R.Z)) % 4;
	R.X = tempX;
	R.Z = tempZ;

	// now the initial state |s> becomes 0.5*(|s> + R |s>) = 0.5*(|s> + i^b |s ^
	// R.X>)
	unsigned b = (R.e + 2 * popcount(R.Z & s)) % 4;
	UpdateSvector(s, s ^ R.X, b);
	// account for the extra factor sqrt(1/2)
	omega.p -= 1;

	isReadyMT = false; // we have changed change M and F
	isReadyFT = false;
}