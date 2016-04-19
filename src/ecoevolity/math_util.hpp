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

#endif
