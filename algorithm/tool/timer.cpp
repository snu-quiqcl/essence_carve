/*
 * Copyright 2026 Yunseo Hwang and Giwon Song.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "timer.hpp"
#include <map>
#include <thread>
#include <omp.h>

std::mutex Timer::log_mutex;

Timer::Timer(const std::string& label, bool chron)
    : active_(chron), label_(label)
{
    if (active_) {
        wall_start_ = std::chrono::steady_clock::now();
        cpu_start_ = std::clock();

        extern std::string testname;
        extern const std::string global_run_timestamp;

        log_filename_ = testname.empty() ? "time_analysis_" + global_run_timestamp + ".txt"
                                         : "time_analysis_" + testname + "_" + global_run_timestamp + ".txt";
    }
}

Timer::~Timer() {
    if (!active_) return;

    auto wall_end = std::chrono::steady_clock::now();
    std::clock_t cpu_end = std::clock();
    auto wall_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(wall_end - wall_start_).count();
    auto cpu_ns = static_cast<long long>(1e9 * (cpu_end - cpu_start_) / CLOCKS_PER_SEC);

    auto lock_wait_start = std::chrono::steady_clock::now();  // lock start for mutex
    std::lock_guard<std::mutex> lock(log_mutex);
    auto lock_acquired = std::chrono::steady_clock::now();  // when grabbed mutex
    auto  mutex_wait_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(lock_acquired - lock_wait_start).count();

    std::ofstream& log = getLogFile(log_filename_);

#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = -1;  // single thread
#endif

    log << "[Timer] " << label_ << " | "
        << "Thread: " << thread_id << " | "
        << "Wall: " << format_duration(wall_ns) << " | "
        << "CPU: "  << format_duration(cpu_ns)  << " | "
        << "MutexWait: " << format_duration(mutex_wait_ns)
        << std::endl;
}

std::ofstream& Timer::getLogFile(const std::string& filename) {
    static std::map<std::string, std::ofstream> logs;

    if (logs.find(filename) == logs.end()) {
        logs[filename] = std::ofstream(filename, std::ios::app);
    }

    return logs[filename];
}