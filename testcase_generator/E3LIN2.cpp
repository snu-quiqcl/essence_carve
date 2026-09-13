/*
 * Copyright 2026 Yunseo Hwang and Giwon Song.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include <vector>
#include <tuple>
#include <random>
#include <algorithm>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <string.h>
#include <math.h>
#include <stdexcept>

namespace {

std::filesystem::path default_testcases_dir() {
	const std::filesystem::path from_repo_root = "testcases";
	if (std::filesystem::exists(from_repo_root) &&
		std::filesystem::is_directory(from_repo_root)) {
		return from_repo_root;
	}

	const std::filesystem::path from_generator_dir = "../testcases";
	if (std::filesystem::exists(from_generator_dir) &&
		std::filesystem::is_directory(from_generator_dir)) {
		return from_generator_dir;
	}

	return from_repo_root;
}

std::filesystem::path resolve_testcase_root(const std::string& output_name) {
	const std::filesystem::path requested(output_name);
	if (requested.is_absolute() || requested.has_parent_path()) {
		return requested;
	}
	return default_testcases_dir() / requested;
}

}  // namespace

class dtuple {
	public:
		bool sign;
		int first, second, third;
		dtuple(bool sign, int first, int second, int third): sign(sign), first(first), second(second), third(third) {}
		friend std::ostream& operator<<(std::ostream& os, const dtuple& T) {
			os << (T.sign ? "-" : "+") << "(Z" << T.first << ")(Z" << T.second << ")(Z" << T.third << ")";
			return os;
		}
};

class dtuple_list {
	public:
		std::vector<dtuple> tuples;
		std::string name;
		dtuple_list(int N, int D): name("N" + std::to_string(N) + "D" + std::to_string(D)), N(N), D(D) {}
		dtuple_list(int N, int D, std::string name): name(name), N(N), D(D) {}
		void generate(int randkey);
		void export_to_file(std::vector<double> gamma_list);
		friend std::ostream& operator<<(std::ostream& os, const dtuple_list& L) {
			for (const dtuple& T: L.tuples) {
				os << T << std::endl;
			}
			return os;
		}
	private:
		int N, D, key;
		std::string build_observable();
		std::string build_circuit(double gamma);
};

void dtuple_list::export_to_file(std::vector<double> gamma_list) {
	std::ofstream file;
	const std::filesystem::path root = resolve_testcase_root(name);
	const std::filesystem::path circuits_dir = root / "circuits";
	std::filesystem::remove_all(root);
	std::filesystem::create_directories(circuits_dir);
	file.open(root / "info.txt");
	file << std::to_string(gamma_list.size());
	file.close();
	file.open(root / "observable.txt");
	file << build_observable();
	file.close();
	file.open(root / "metadata.txt");
	file << "generator,E3LIN2\n";
	file << "num_qubits," << N << "\n";
	file << "degree," << D << "\n";
	file << "num_gammas," << gamma_list.size() << "\n";
	file << "num_circuits," << gamma_list.size() << "\n";
	file << "random_seed," << key << "\n";
	file << "num_clauses," << tuples.size() << "\n";
	file << "gamma_min,0\n";
	file << "gamma_max,pi/2\n";
	file << "observable,Pauli ZZZ clause sum\n";
	file << "ansatz,H^n clause phases then SX^n\n";
	file.close();
	int index = 0;
	for (double& gamma: gamma_list) {
		file.open(circuits_dir / (std::to_string(index++) + ".txt"));
		file << build_circuit(gamma);
		file.close();
	}
	std::cout << "Successfully exported " << gamma_list.size()
	          << " circuits to " << root << std::endl;
}

void dtuple_list::generate(int randkey) { 
	tuples.clear();
	tuples.reserve(N * D / 3);
	std::mt19937 g(randkey);
	key = randkey;
	bool check;
	do {
		check = true;
		std::cout << "Generating " << N * D / 3 << " tuples containing three 0 to " << N-1 << " random numbers that appears " << D << " times each..." << std::endl;

		std::vector<int> selection(N * D);
		std::iota(selection.begin(), selection.end(), 0);
		std::shuffle(selection.begin(), selection.end(), g);

		for (int i = 0; i < N * D / 3; i++) {
			int first = selection[0] / D;
			selection.erase(std::remove(selection.begin(), selection.end(), selection[0]), selection.end());
			int second = -1;
			for (int candidate: selection)
				if (candidate / D != first) {
					second = candidate / D;
					selection.erase(std::remove(selection.begin(), selection.end(), candidate), selection.end());
					break;
				}
			if (second == -1) {
				std::cout << "Generation failed at i = " << i << ". Retrying..." << std::endl;
				tuples.clear();
				check = false;
				break;
			}
			int third = -1;
			for (int candidate: selection)
				if (candidate / D != first && candidate / D != second) {
					third = candidate / D;
					selection.erase(std::remove(selection.begin(), selection.end(), candidate), selection.end());
					break;
				}
			if (third == -1) {
				std::cout << "Generation failed at " << i << "-th tuple. Retrying..." << std::endl;
				tuples.clear();
				check = false;
				break;
			}
			
			double r = std::uniform_real_distribution<double>(0,1)(g);
			tuples.emplace_back(r < 0.5, first, second, third);
		}
	} while (!check);
	return;
}

std::string dtuple_list::build_observable() {
	int G = tuples.size();
	std::string s = std::to_string(G) + "\n";
	for (const dtuple& T: tuples)
		s += std::to_string(T.sign ? -1 : 1) + " pauli ZZZ 3 " + std::to_string(T.first) + " " + std::to_string(T.second) + " " + std::to_string(T.third) + "\n";
	return s;
}

std::string dtuple_list::build_circuit(double gamma) {
	int pi2 = std::round(gamma * 2.0 / M_PI);
	std::string gate = "u1";
	if (std::abs(pi2 * M_PI / 2 - gamma) < 1e-8)
		gate = "rz";
	int G = 5 * tuples.size() + 2*N;
	std::string s = std::to_string(N) + " " + std::to_string(G) + "\n";
	for (int i = 0; i < N; i++)
		s += "h " + std::to_string(i) + "\n";
	for (const dtuple& T: tuples) {
		s += "cx " + std::to_string(T.first) + " " + std::to_string(T.third) + "\n";
		s += "cx " + std::to_string(T.second) + " " + std::to_string(T.third) + "\n";
		s += gate + " " + std::to_string(T.third) + " " + std::to_string((T.sign ? -1 : 1) * gamma) + "\n";
		s += "cx " + std::to_string(T.second) + " " + std::to_string(T.third) + "\n";
		s += "cx " + std::to_string(T.first) + " " + std::to_string(T.third) + "\n";
	}
	for (int i = 0; i < N; i++)
		s += "sx " + std::to_string(i) + "\n";
	return s;
}

int main(int argc, char** argv) {
	try {
		int N = 0;
		int D = 0;
		int Ngamma = 0;
		int randkey = 19;
		std::string output_name;

		if (argc >= 4) {
			N = std::stoi(argv[1]);
			D = std::stoi(argv[2]);
			Ngamma = std::stoi(argv[3]);
			if (argc >= 5) randkey = std::stoi(argv[4]);
			if (argc >= 6) output_name = argv[5];
		} else {
			std::cout << "N: ";
			std::cin >> N;
			std::cout << "D: ";
			std::cin >> D;
			do {
				std::cout << "# of gammas: ";
				std::cin >> Ngamma;
			} while (Ngamma < 2);
		}

		if (N < 3) {
			std::cerr << "N must be at least 3." << std::endl;
			return 1;
		}
		if (D < 1) {
			std::cerr << "D must be at least 1." << std::endl;
			return 1;
		}
		if (Ngamma < 2) {
			std::cerr << "Ngamma must be at least 2." << std::endl;
			return 1;
		}

		dtuple_list L = output_name.empty()
			? dtuple_list(N, D)
			: dtuple_list(N, D, output_name);
		L.generate(randkey);

		std::vector<double> gamma_list;
		for (int i = 0; i < Ngamma; i++) {
			gamma_list.push_back(i / (double)(Ngamma - 1) * M_PI / 2);
		}
		L.export_to_file(gamma_list);
	} catch (const std::exception& ex) {
		std::cerr << "Error: " << ex.what() << std::endl;
		std::cerr << "Usage: <program> <N> <D> <Ngamma> [seed, default 19] "
		          << "[output_name, bare names are written under testcases/]" << std::endl;
		return 1;
	}
}
