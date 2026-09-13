/*
 * Copyright 2026 Yunseo Hwang and Giwon Song.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "welford.hpp"

void WelfordStats::update(double x) {
    count++;
    double delta = x - mean;
    mean += delta / count;
    M2 += delta * (x - mean);
}

void WelfordStats::merge(const WelfordStats& other) {
    if (other.count == 0) return;
    if (count == 0) {
        *this = other;
        return;
    }
    double delta = other.mean - mean;
    int total = count + other.count;
    mean = (mean * count + other.mean * other.count) / total;
    M2 += other.M2 + delta * delta * count * other.count / total;
    count = total;
}

double WelfordStats::variance() const {
    return (count > 1) ? M2 / (count - 1) : 0.0;
}

bool WelfordStats::operator<(const WelfordStats& other) const {
    if(count <= 1 || other.count <= 1){
        return false;
    }
    return M2/(count-1) < other.M2/(other.count-1);
}
