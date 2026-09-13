/*
 * Copyright 2026 Yunseo Hwang and Giwon Song.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <string>
#include <fstream>
#include <mutex>
#include <ctime>
#include <sstream>
#include <iomanip>

class Timer {
public:
    Timer(const std::string& label, bool chron = true);
    ~Timer();

private:
    bool active_;
    std::string label_;
    std::string log_filename_;
    std::chrono::steady_clock::time_point wall_start_;
    std::clock_t cpu_start_;

    static std::mutex log_mutex;
    static std::ofstream& getLogFile(const std::string& filename);
};

inline std::string format_duration(long long nanoseconds) {
    using namespace std::chrono;
    std::ostringstream out;
    out << std::fixed << std::setprecision(3);

    if (nanoseconds >= 60'000'000'000LL) {
        out << (nanoseconds / 60'000'000'000.0) << " min";
    } else if (nanoseconds >= 1'000'000'000LL) {
        out << (nanoseconds / 1'000'000'000.0) << " s";
    } else if (nanoseconds >= 1'000'000LL) {
        out << (nanoseconds / 1'000'000.0) << " ms";
    } else if (nanoseconds >= 1'000LL) {
        out << (nanoseconds / 1'000.0) << " μs";
    } else {
        out << nanoseconds << " ns";
    }
    return out.str();
}

inline const std::string global_run_timestamp = [] {
    auto now = std::chrono::system_clock::now();
    auto t = std::chrono::system_clock::to_time_t(now);
    std::tm tm;
#ifdef _WIN32
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
    return oss.str();
}();

#endif