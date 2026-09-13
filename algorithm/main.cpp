/*
 * Copyright 2026 Yunseo Hwang and Giwon Song.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "extended_stabilizer_state.hpp"
#include "file_reader.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <ctime>
#include <filesystem>

std::string testfilepath = "../testcases/";
std::string testname;

struct arguments {
	// general
	uint_t N, num_threads = omp_get_num_procs();
	double delta = 0.2, epsilon = 0.2;
	bool bernstein = true;
	uint_t n_shot = 0;
};

std::pair<Observable, circuit_list_t> load(std::string name) {
	CircuitFile file(testfilepath + name);
	Observable obv = file.get_observable();
	circuit_list_t circuits = file.get_circuits();
	return {obv, circuits};
}

void write_detail_header_if_needed(const std::string& filename) {
	std::ofstream fout(filename, std::ios::app);
	if (!fout.is_open()) {
		throw std::logic_error("Cannot open the detail CSV file");
	}
	if (fout.tellp() == 0) {
		fout << "circuit_id,n_qubits,gate_count,observable_size,is_projector,is_stabilizer,"
		     << "epsilon,delta,bernstein,n_threads,xi,h_old,h_new,result,"
		     << "requested_shots,variance_shots,"
		     << "variance_hat,std_hat,sigma_max,variance_bound,observable_norm,observable_norm_is_exact,"
		     << "cpu_time_ns,wall_time_ns\n";
	}
}

void append_detail_row(const std::string& filename, uint_t circuit_id, int n_qubits,
					   size_t gate_count, const ExpvalStats& stats,
					   uint_t cpu_time_ns, uint_t wall_time_ns) {
	std::ofstream fout(filename, std::ios::app);
	if (!fout.is_open()) {
		throw std::logic_error("Cannot open the detail CSV file");
	}
	fout << std::scientific << std::setprecision(12);
	fout << circuit_id << "," << n_qubits << "," << gate_count << ","
	     << stats.observable_size << "," << stats.is_projector << "," << stats.is_stabilizer << "," 
		 << stats.epsilon << "," << stats.delta << "," << stats.bernstein_enabled << "," << stats.n_threads << ","
		 << stats.xi << "," << stats.hoeffding_shots_initial << "," << stats.final_shots << "," << stats.result << "," 
	     << stats.requested_shots << "," << stats.variance_shots << "," 
	     << stats.variance_hat << "," << stats.std_hat << ","
	     << stats.sigma_max << "," << stats.variance_bound << ","
	     << stats.observable_norm << "," << stats.observable_norm_is_exact << ","
	     << cpu_time_ns << "," << wall_time_ns << "\n";
}

void test(std::string name, arguments args) {
	const std::filesystem::path data_dir = "../data";
	std::filesystem::create_directories(data_dir);
	const std::filesystem::path detail_path =
		data_dir / ("details_" + name + "_" + global_run_timestamp + ".csv");
	const std::string detail_filename = detail_path.string();
	write_detail_header_if_needed(detail_filename);

	Observable obv;
	circuit_list_t circuits;
	std::tie(obv, circuits) = load(name);
	uint_t circuit_id = 0;
	for (auto& circuit: circuits) {
		int N;
		std::vector<Op> ops;
		std::tie(N, ops) = circuit;
		State state(N);
		double xi = state.compute_xi(ops.data(), ops.size());
		std::cout << "================================================\nn: " << N << "\nxi: " << xi << "\nobservable: " << obv.size << " different gates\n================================================" << std::endl;
		
		auto start = std::chrono::steady_clock::now();
		std::clock_t cpu_start = std::clock();
		ExpvalStats stats = state.expval_with_stats(ops.data(), ops.size(), obv, args.num_threads, args.n_shot, args.delta, args.epsilon, args.bernstein);
		std::clock_t cpu_end = std::clock();
		auto end = std::chrono::steady_clock::now();
		uint_t time_passed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
		uint_t cpu_time_passed = static_cast<uint_t>(1e9 * (cpu_end - cpu_start) / CLOCKS_PER_SEC);

		std::cout << "Result: " << stats.result << "" << std::endl;	
		std::cout << "Elapsed time: " << format_duration(time_passed) << std::endl;
		append_detail_row(detail_filename, circuit_id, N, ops.size(), stats, cpu_time_passed, time_passed);
		++circuit_id;
	}
}

int main() {
	int num_threads;
	int n_shot;
	double delta, epsilon;
	std::string bernstein;
	arguments args;
	std::string ans;	std::string ans2;	std::string ans3;
	uint_t iter = 1;
	std::cout << "Default setting? (y/n): "; std::cin >> ans;
	if (ans == "n") {
		std::cout << "Thread: "; std::cin >> num_threads; args.num_threads = num_threads;
		std::cout << "Fix n_shot ? (y/n): "; std::cin >> ans3;
			if(ans3 == "y"){
				std::cout << "n_shot: "; std::cin >> n_shot; args.n_shot = n_shot;
				args.bernstein = false;
			}
			else{
				std::cout << "Delta: "; std::cin >> delta; args.delta = delta;
				std::cout << "Epsilon: "; std::cin >> epsilon; args.epsilon = epsilon;
				std::cout << "Use Empirical Variance? (y/n): "; std::cin >> bernstein; args.bernstein = (bernstein == "y");
			}
		std::cout << "Iteration: "; std::cin >> iter;
	}

	std::cout << "test file name: "; std::cin >> testname;
	
	for (uint_t i = 0; i < iter; i++){
		test(testname, args);
	}
}