/*
 * Copyright 2026 Yunseo Hwang and Giwon Song.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef MYSTR_HPP
#define MYSTR_HPP

#include <cstddef>

#pragma omp declare target

inline int my_strcmp(const char* s1, const char* s2) {
    while (*s1 && (*s1 == *s2)) {
        s1++;
        s2++;
    }
    return *(const unsigned char*)s1 - *(const unsigned char*)s2;
}

inline char* my_strcpy(char* destination, const char* source) {
    char* ptr = destination;
    while (*source != '\0') {
        *destination = *source;
        destination++;
        source++;
    }
    *destination = '\0';
    return ptr;
}
 
inline size_t my_strlen(const char* s) {
    size_t len = 0;
    if(s != nullptr){
        while (s[len] != '\0') {
            len++;
        }
    }
    return len;
}

inline double my_pow(double base, double exponent) {
    return std::exp(exponent * std::log(base));
}

#pragma omp end declare target

#endif