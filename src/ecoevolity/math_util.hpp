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

#endif
