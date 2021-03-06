/******************************************************************************
 * Copyright (C) 2015-2016 Jamie R. Oaks.
 *
 * This file is part of Ecoevolity.
 *
 * Ecoevolity is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Ecoevolity is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Ecoevolity.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#ifndef ECOEVOLITY_MATH_UTIL_HPP
#define ECOEVOLITY_MATH_UTIL_HPP

#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <numeric>

#include "assert.hpp"
#include "error.hpp"

inline void normalize_log_likelihoods(std::vector<double>& v) {
    double mx = v.at(0);
    for (auto v_iter : v) {
        if (v_iter > mx) {
            mx = v_iter;
        }
    }

    for (unsigned int i = 0; i < v.size(); ++i) {
        v.at(i) -= mx;
    }

    double sum = 0.0;
    for (unsigned int i = 0; i < v.size(); ++i) {
        v.at(i) = std::exp(v.at(i));
        sum += v.at(i);
    }

    double t = 0.0;
    for (unsigned int i = 0; i < v.size(); ++i) {
        v.at(i) /= sum;
        t += v.at(i);
    }
    ECOEVOLITY_ASSERT_APPROX_EQUAL(t, 1.0);
}

/**
 * Calculate the expected number of categories under a Dirichlet process.
 *
 * Calculates and returns the expected number of categories for
 * 'number_of_elements' elements under a Dirichlet process controlled by
 * 'concentration' parameter.
 *
 * @note    Modified from `expNumTables` function of `util.h` from
 *          [`DPPDiv`](http://phylo.bio.ku.edu/content/tracy-heath-dppdiv)
 *          version 1.0b (Copyright Tracy Heath, Mark Holder, and John
 *          Huelsenback; licensed under GPL v3;
 *          <http://phylo.bio.ku.edu/content/tracy-heath-dppdiv>).
 */
inline double get_dpp_expected_number_of_categories(
        double concentration,
        unsigned int number_of_elements) {
    double expected_ncats =  0.0;
    for (unsigned int i = 1; i <= number_of_elements; ++i) {
        expected_ncats += (1.0 / (((double) i) - 1.0 + concentration));
    }
    return expected_ncats * concentration;
}

/**
 * Calculate the Dirichlet-process concentration parameter.
 *
 * Calculates and returns the Dirichlet-process concentration parameter that
 * has an expected number of categories equal to
 * 'expected_number_of_categories' for 'number_of_elements' elements.
 *
 * @note    Modified from `calculateFromPriorMean` function of `util.h` from
 *          [`DPPDiv`](http://phylo.bio.ku.edu/content/tracy-heath-dppdiv)
 *          version 1.0b (Copyright Tracy Heath, Mark Holder, and John
 *          Huelsenback; licensed under GPL v3;
 *          <http://phylo.bio.ku.edu/content/tracy-heath-dppdiv>).
 */
inline double get_dpp_concentration(
        double expected_number_of_categories,
        unsigned int number_of_elements,
        double increment = 0.1,
        double precision = 0.000001,
        double buffer = 0.001) {
    double c = precision;
    double n = (double) number_of_elements;
    double e = expected_number_of_categories;
    double current_e = get_dpp_expected_number_of_categories(c, n);

    ECOEVOLITY_ASSERT(e <= n);
    ECOEVOLITY_ASSERT(e >= 1.0);
    if (e > (n - buffer)) {
        e = n - buffer;
    }
    else if (e < (1.0 + buffer)) {
        e = 1.0 + buffer;
    }
    bool increase = false;
    if (current_e < e) {
        increase = true;
    }
    while (fabs(current_e - e) > precision) {
        if ((current_e < e) && increase) {
            c += increment;
        }
        else if ((current_e > e) && (! increase)) {
            c -= increment;
        }
        else if ((current_e < e) && (! increase)) {
            increment /= 2.0;
            increase = true;
            c += increment;
        }
        else {
            increment /= 2.0;
            increase = false;
            c -= increment;
        }
        current_e = get_dpp_expected_number_of_categories(c, n);
    }
    return c;
}

/**
 * Calculate the scale parameter of a gamma hyper prior on the concentration
 * parameter of the Dirichlet process.
 */
inline double get_dpp_gamma_scale(
        double expected_number_of_categories,
        unsigned int number_of_elements,
        double shape,
        double increment = 0.1,
        double precision = 0.000001,
        double buffer = 0.001) {
    double concentration = get_dpp_concentration(
            expected_number_of_categories,
            number_of_elements,
            increment,
            precision,
            buffer);
    return concentration / shape;
}

template <typename T>
inline double get_dpp_log_prior_probability(
        const std::vector<T>& partition,
        double concentration) {
    ECOEVOLITY_ASSERT(concentration > 0.0);
    double log_concentration = std::log(concentration);
    double log_p = 0.0;
    T current_element;
    std::map<T, unsigned int> subset_counts;
    for (unsigned int i = 0; i < partition.size(); ++i) {
        current_element = partition.at(i);
        if (subset_counts.count(current_element) < 1) {
            log_p += (log_concentration - std::log(concentration + i));
            subset_counts[current_element] = 1;
            continue;
        }
        log_p += (std::log(subset_counts[current_element]) - std::log(concentration + i));
        ++subset_counts[current_element];
    }
    return log_p;
}

inline double get_dpp_log_prior_probability(
        const std::string& partition,
        double concentration) {
    std::vector<char> partition_vector(partition.begin(), partition.end());
    return get_dpp_log_prior_probability<char>(partition_vector, concentration);
}


// Recursion is much slower (~ 1000 times slower!) than dynamic programming
// approach below.
// inline unsigned long long stirling2_recurse(int n, int k) {
//     ECOEVOLITY_ASSERT((n > 0) && (k > 0) && (n >= k));
//     if (n == 0 || k == 1 || k == n) {
//         return 1;
//     }
//     return stirling2_recurse(n - 1, k - 1) + k * stirling2_recurse(n - 1, k);
// }

template <typename T>
inline T stirling2_base(int n, int k) {
    ECOEVOLITY_ASSERT((n > 0) && (k > 0) && (n >= k));
    int maxj = n-k;

    std::vector<T> s_numbers (maxj + 1, 1);

    for (int i = 2; i <= k; ++i) {
        for(int j = 1; j <= maxj; ++j) {
            s_numbers.at(j) += i * s_numbers.at(j-1);
        }
    }

    return s_numbers.at(maxj);
}

inline unsigned long long stirling2(int n, int k) {
    return stirling2_base<unsigned long long>(n, k);
}

inline long double stirling2_float(int n, int k) {
    return stirling2_base<long double>(n, k);
}

template <typename T>
inline T bell_number_base(int n) {
    T b = 0;
    for (int i = 0; i < n; ++i) {
        b += stirling2_base<T>(n, (i + 1));
    }
    return b;
}

inline unsigned long long bell_number(int n) {
    return bell_number_base<unsigned long long>(n);
}

inline long double bell_float(int n) {
    return bell_number_base<long double>(n);
}

template <typename T>
inline double get_uniform_model_log_prior_probability(
        const unsigned int number_of_elements,
        const unsigned int number_of_categories,
        const double split_weight) {
    ECOEVOLITY_ASSERT(split_weight > 0.0);
    ECOEVOLITY_ASSERT(number_of_elements > 0);
    ECOEVOLITY_ASSERT(number_of_categories > 0);
    ECOEVOLITY_ASSERT(number_of_categories <= number_of_elements);
    T denom = 0;
    for (unsigned int k = 1; k <= number_of_elements; ++k) {
        denom += stirling2_base<T>(number_of_elements, k) * std::pow(split_weight, (k - 1));
    }
    return ((number_of_categories - 1) * std::log(split_weight)) - std::log(denom);
}

inline double get_uniform_model_log_prior_probability(
        const unsigned int number_of_elements,
        const unsigned int number_of_categories,
        const double split_weight) {
    return get_uniform_model_log_prior_probability<long double>(
            number_of_elements,
            number_of_categories,
            split_weight);
}

/**
 * Get all possible integer partitions of n into m parts.
 *
 * @note    This is an implementation of Algorithm H on Page 2 of:
 *          Donald E. Knuth. Generating all partitions, 2004. Pre-fascicle 3B
 *          of The Art of Computer Programming, A draft of sections 7.2.1.4–5
 *          http://www-cs-faculty.stanford.edu/~knuth/fasc3b.ps.gz
 */
inline std::vector< std::vector<unsigned int> > get_integer_partitions(
        int n,
        int m) {
    ECOEVOLITY_ASSERT(n >= m);
    ECOEVOLITY_ASSERT(m > 1);
    std::vector< std::vector<unsigned int> > partitions;
    std::vector<int> a(m + 1, 1);
    a.at(0) = n - m + 1;
    a.back() = -1;
    while (true) {
        std::vector<unsigned int> p(a.begin(), a.end()-1);
        partitions.push_back(p);
        while (a.at(1) < (a.at(0) - 1)) {
            --a.at(0);
            ++a.at(1);
            std::vector<unsigned int> p(a.begin(), a.end()-1);
            partitions.push_back(p);
        }
        int j = 3;
        int s = a.at(0) + a.at(1) - 1;
        while (a.at(j-1) >= (a.at(0) - 1)) {
            s = s + a.at(j-1);
            ++j;
        }
        if (j > m) {
            break;
        }
        int x = a.at(j-1) + 1;
        a.at(j-1) = x;
        --j;
        while (j > 1) {
            a.at(j-1) = x;
            s = s - x;
            --j;
        }
        a.at(0) = s;
    }
    return partitions;
}

/**
 * Calculate log factorial.
 *
 * @note    Adapted from Daniel's post on stackoverflow.com
 *          (http://stackoverflow.com/a/30630392)
 */
inline double log_factorial(unsigned int x) {
    if (x == 0 || x == 1) return 0;
    if (x == 2) {
        return std::log(2); // can add more for efficiency
    }
    if (x > 100) {
        return x * std::log(x) - x; // Stirling's approximation
    }
    std::vector<unsigned int> lx(x);
    std::iota(lx.begin(), lx.end(), 1);
    std::vector<double> tx(x);
    std::transform(lx.cbegin(), lx.cend(), tx.begin(),
                   [] (unsigned int a) { return std::log(static_cast<double>(a)); });
    return std::accumulate(tx.cbegin(), tx.cend(), double {});
}

/**
 * Calculate log multinomial coefficient
 *
 * @note    Adapted from Daniel's post on stackoverflow.com
 *          (http://stackoverflow.com/a/30630392)
 */
inline double log_multinomial_coefficient(const std::vector<unsigned int>& il) {
    ECOEVOLITY_ASSERT(il.size() > 1);
    std::vector<double> denoms(il.size());
    std::transform(il.begin(), il.end(), denoms.begin(), log_factorial);
    unsigned int numerator = std::accumulate(il.begin(), il.end(), 0);
    double denominator = std::accumulate(denoms.begin(), denoms.end(), 0.0);
    return log_factorial(numerator) - denominator;
}

#endif
