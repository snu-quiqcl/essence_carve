/*
 * Copyright 2026 Yunseo Hwang and Giwon Song.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef FILE_READER_HPP
#define FILE_READER_HPP
#include "operators.hpp"
#include <fstream>
#include <vector>

using circuit_list_t = std::vector<std::pair<int, std::vector<Op>>>;

class CircuitFile {
	public:
		std::string filepath;
		CircuitFile(std::string filepath): filepath(filepath) {}
		Observable get_observable();
		circuit_list_t get_circuits();
	private:
		Op readline(std::fstream& file);
		Op readline(std::fstream& file, const std::string& name);
		void read_projector(std::fstream& file, Observable& obv, double coeff);
};

#endif
