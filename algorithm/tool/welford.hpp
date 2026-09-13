/*
 * Copyright 2026 Yunseo Hwang and Giwon Song.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef WELFORD_HPP
#define WELFORD_HPP

#include <omp.h>

struct WelfordStats {
    double mean = 0;
    double M2 = 0;
    long long count = 0;

    WelfordStats() : mean(0.0), M2(0.0), count(0) {}

    void update(double x);
    void merge(const WelfordStats& other);
    double variance() const;
    bool operator<(const WelfordStats& other) const;
};

inline WelfordStats combine_welford(const WelfordStats& a, const WelfordStats& b) {
    if (b.count == 0) return a;
    if (a.count == 0) return b;

    double delta = b.mean - a.mean;
    long long total_count = a.count + b.count;
    
    WelfordStats result;
    result.count = total_count;
    result.mean = (a.mean * a.count + b.mean * b.count) / total_count;
    result.M2 = a.M2 + b.M2 + delta * delta * a.count * b.count / total_count;
    
    return result;
}

#pragma omp declare reduction(welford_reduction: WelfordStats: \
    omp_out = combine_welford(omp_out, omp_in)) \
    initializer(omp_priv = WelfordStats())

#endif 