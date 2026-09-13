/*
 * Copyright 2026 Yunseo Hwang and Giwon Song.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

constexpr double ANGLE_EPS = 1.0e-15;
constexpr const char* PROJECTOR_OBSERVABLE_TOKEN = "projector";

std::string format_double(double value) {
    if (std::abs(value) < ANGLE_EPS) value = 0.0;
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

double canonical_phase_angle(double angle) {
    const double wrapped = std::remainder(angle, 2.0 * M_PI);
    return (std::abs(wrapped) < ANGLE_EPS) ? 0.0 : wrapped;
}

enum class pair_mode_t {
    full_cartesian,
    ordered_excluding_self,
    unordered_excluding_self,
};

std::uint64_t pair_count_for_mode(std::uint64_t n_samples, pair_mode_t mode) {
    switch (mode) {
        case pair_mode_t::full_cartesian:
            return n_samples * n_samples;
        case pair_mode_t::ordered_excluding_self:
            return n_samples * (n_samples - 1);
        case pair_mode_t::unordered_excluding_self:
            return (n_samples * (n_samples - 1)) / 2;
    }
    return 0;
}

std::string pair_mode_name(pair_mode_t mode) {
    switch (mode) {
        case pair_mode_t::full_cartesian:
            return "full_cartesian";
        case pair_mode_t::ordered_excluding_self:
            return "ordered_excluding_self";
        case pair_mode_t::unordered_excluding_self:
            return "unordered_excluding_self";
    }
    return "unknown";
}

pair_mode_t parse_pair_mode(const std::string& mode) {
    if (mode == "full" || mode == "full_cartesian") {
        return pair_mode_t::full_cartesian;
    }
    if (mode == "ordered" || mode == "ordered_excluding_self") {
        return pair_mode_t::ordered_excluding_self;
    }
    if (mode == "unordered" || mode == "unordered_excluding_self") {
        return pair_mode_t::unordered_excluding_self;
    }
    throw std::invalid_argument("pair_mode must be full, ordered, or unordered");
}

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

class random_pauli_feature_dataset {
private:
    int sample_count;
    int num_qubits;
    std::uint64_t seed;
    pair_mode_t pair_mode;

public:
    std::string name;
    std::vector<std::vector<double>> data_points;

    random_pauli_feature_dataset(int n_samples,
                                 int n_qubits,
                                 std::uint64_t random_seed,
                                 pair_mode_t mode = pair_mode_t::full_cartesian)
        : sample_count(n_samples),
          num_qubits(n_qubits),
          seed(random_seed),
          pair_mode(mode),
          name("QKE_RandomPauli_N" + std::to_string(n_samples) +
               "_n" + std::to_string(n_qubits)) {}

    void generate_data();
    void export_to_file() const;

private:
    std::string build_observable() const;
    std::string build_kernel_circuit(const std::vector<double>& x,
                                     const std::vector<double>& xp) const;
};

void random_pauli_feature_dataset::generate_data() {
    data_points.clear();
    data_points.reserve(static_cast<std::size_t>(sample_count));

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, 2.0 * M_PI);

    for (int sample = 0; sample < sample_count; ++sample) {
        std::vector<double> x(static_cast<std::size_t>(num_qubits));
        for (double& value : x) value = dist(rng);
        data_points.push_back(std::move(x));
    }
}

std::string random_pauli_feature_dataset::build_observable() const {
    // Pi_0 = |0^n><0^n| as a rank-one stabilizer projector.
    return std::string(PROJECTOR_OBSERVABLE_TOKEN) + " 0\n";
}

std::string random_pauli_feature_dataset::build_kernel_circuit(
    const std::vector<double>& x,
    const std::vector<double>& xp) const {
    if (static_cast<int>(x.size()) != num_qubits ||
        static_cast<int>(xp.size()) != num_qubits) {
        throw std::invalid_argument("kernel-circuit input vector has wrong dimension");
    }

    std::vector<std::string> gates;
    gates.reserve(static_cast<std::size_t>(3 * num_qubits + 3 * (num_qubits - 1)));

    // Rightmost H^{otimes n}: maps |0^n> to |+^n>.
    for (int i = 0; i < num_qubits; ++i) {
        gates.push_back("h " + std::to_string(i));
    }

    // U^dagger(x) U(x') single-qubit phase differences for exp(i x_i Z_i).
    // u1(lambda) is equivalent to RZ(lambda) up to global phase, so
    // exp(i delta Z) is emitted as u1(-2 delta).
    for (int i = 0; i < num_qubits; ++i) {
        const double dz = xp[static_cast<std::size_t>(i)] - x[static_cast<std::size_t>(i)];
        const double phase_angle = canonical_phase_angle(-2.0 * dz);
        if (std::abs(phase_angle) > ANGLE_EPS) {
            gates.push_back("u1 " + std::to_string(i) + " " + format_double(phase_angle));
        }
    }

    // Nearest-neighbor 1D ZZ terms from the original Pauli feature map:
    // exp(i (pi - x_j)(pi - x_{j+1}) Z_j Z_{j+1}).
    for (int j = 0; j + 1 < num_qubits; ++j) {
        const double theta_x = (M_PI - x[static_cast<std::size_t>(j)]) *
                               (M_PI - x[static_cast<std::size_t>(j + 1)]);
        const double theta_xp = (M_PI - xp[static_cast<std::size_t>(j)]) *
                                (M_PI - xp[static_cast<std::size_t>(j + 1)]);
        const double dzz = theta_xp - theta_x;
        const double phase_angle = canonical_phase_angle(-2.0 * dzz);
        if (std::abs(phase_angle) <= ANGLE_EPS) continue;

        gates.push_back("cx " + std::to_string(j) + " " + std::to_string(j + 1));
        gates.push_back("u1 " + std::to_string(j + 1) + " " + format_double(phase_angle));
        gates.push_back("cx " + std::to_string(j) + " " + std::to_string(j + 1));
    }

    // Leftmost H^{otimes n} in V(x, x') = H^n U^dagger(x) U(x') H^n.
    for (int i = 0; i < num_qubits; ++i) {
        gates.push_back("h " + std::to_string(i));
    }

    std::string s = std::to_string(num_qubits) + " " + std::to_string(gates.size()) + "\n";
    for (const std::string& gate : gates) s += gate + "\n";
    return s;
}

void random_pauli_feature_dataset::export_to_file() const {
    if (data_points.size() != static_cast<std::size_t>(sample_count)) {
        throw std::logic_error("data has not been generated");
    }

    const std::uint64_t pair_count =
        pair_count_for_mode(static_cast<std::uint64_t>(sample_count), pair_mode);

    const std::filesystem::path root = resolve_testcase_root(name);
    const std::filesystem::path circuits_dir = root / "circuits";
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(circuits_dir);

    {
        std::ofstream file(root / "info.txt");
        file << pair_count << "\n";
    }

    {
        std::ofstream file(root / "observable.txt");
        file << build_observable();
    }

    {
        std::ofstream file(root / "metadata.txt");
        file << "num_samples," << sample_count << "\n";
        file << "num_qubits," << num_qubits << "\n";
        file << "num_kernel_circuits," << pair_count << "\n";
        file << "random_seed," << seed << "\n";
        file << "distribution,uniform[0,2*pi]\n";
        file << "pair_mode," << pair_mode_name(pair_mode) << "\n";
        file << "circuit_definition,V(x,x')=H^n U^dagger(x) U(x') H^n\n";
        file << "observable,Pi_0=|0^n><0^n|\n";
        file << "phase_gate,u1(theta) equivalent to RZ(theta) up to global phase\n";
        file << "feature_map,U(x)=exp(i sum_i x_i Z_i + i sum_j (pi-x_j)(pi-x_{j+1}) Z_j Z_{j+1})\n";
        file << "single_phase,x_i\n";
        file << "zz_phase,(pi-x_j)*(pi-x_{j+1})\n";
        file << "zz_connectivity,1d_nearest_neighbor_chain\n";
        file << "zz_pair_count," << std::max(0, num_qubits - 1) << "\n";
    }

    {
        std::ofstream file(root / "samples.csv");
        file << "sample_id";
        for (int i = 0; i < num_qubits; ++i) file << ",x" << i;
        file << "\n";
        for (int sample = 0; sample < sample_count; ++sample) {
            file << sample;
            for (double value : data_points[static_cast<std::size_t>(sample)]) {
                file << "," << format_double(value);
            }
            file << "\n";
        }
    }

    std::ofstream pairs(root / "pairs.csv");
    pairs << "pair_id,circuit_file,sample_x,sample_xprime\n";

    std::uint64_t pair_id = 0;
    for (int a = 0; a < sample_count; ++a) {
        const int b_start =
            (pair_mode == pair_mode_t::unordered_excluding_self) ? a + 1 : 0;
        for (int b = b_start; b < sample_count; ++b) {
            if (a == b && pair_mode != pair_mode_t::full_cartesian) continue;

            const std::string filename = std::to_string(pair_id) + ".txt";
            {
                std::ofstream circuit(circuits_dir / filename);
                circuit << build_kernel_circuit(data_points[static_cast<std::size_t>(a)],
                                                data_points[static_cast<std::size_t>(b)]);
            }
            pairs << pair_id << ",circuits/" << filename << "," << a << "," << b << "\n";
            ++pair_id;
        }
    }

    std::cout << "Successfully exported " << pair_count << " kernel circuits for "
              << sample_count << " random samples to " << root << std::endl;
}

int main(int argc, char** argv) {
    try {
        int sample_count = 0;
        int num_qubits = 0;
        std::uint64_t seed = 5489;
        std::string output_name;
        pair_mode_t pair_mode = pair_mode_t::full_cartesian;

        if (argc >= 3) {
            sample_count = std::stoi(argv[1]);
            num_qubits = std::stoi(argv[2]);
            if (argc >= 4) seed = static_cast<std::uint64_t>(std::stoull(argv[3]));
            if (argc >= 5) output_name = argv[4];
            if (argc >= 6) pair_mode = parse_pair_mode(argv[5]);
        } else {
            std::cout << "--- Random Pauli Feature Map Kernel Generator ---" << std::endl;
            std::cout << "N (number of random vectors): ";
            std::cin >> sample_count;
            std::cout << "n (number of qubits/features): ";
            std::cin >> num_qubits;
            std::cout << "Pair mode (full/ordered/unordered, default full): ";
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::string pair_mode_input;
            std::getline(std::cin, pair_mode_input);
            pair_mode = pair_mode_input.empty()
                ? pair_mode_t::full_cartesian
                : parse_pair_mode(pair_mode_input);
        }

        if (sample_count < 2) {
            std::cerr << "N must be at least 2." << std::endl;
            return 1;
        }
        if (num_qubits < 1) {
            std::cerr << "n must be at least 1." << std::endl;
            return 1;
        }

        random_pauli_feature_dataset dataset(sample_count, num_qubits, seed, pair_mode);
        if (!output_name.empty()) dataset.name = output_name;
        dataset.generate_data();
        dataset.export_to_file();
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        std::cerr << "Usage: <program> <N_random_vectors> <n_qubits> [seed] "
                  << "[output_name, bare names are written under testcases/] "
                  << "[pair_mode:full|ordered|unordered, default full]" << std::endl;
        return 1;
    }

    return 0;
}
