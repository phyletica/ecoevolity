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

/**
 * Functions and macros for verbose assertions.
 *
 * @note Modified from:
 *       colugo-utilities library licensed under Gnu GPL 3
 *       <https://github.com/jeetsukumaran/colugo-utilities/blob/master/include/colugo/assert.hpp>
 *       2013 Jeet Sukumaran <http://jeetworks.org>
 */

#ifndef ECOEVOLITY_ASSERT_HPP
#define ECOEVOLITY_ASSERT_HPP

#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>

inline void ecoevolity_assertion_failed(char const * expr, char const * function, char const * file, long line) {
    std::cerr << std::endl <<
            "Assertion failed:"     << std::endl <<
            "\texpr: " << expr      << std::endl <<
            "\tfunc: " << function  << std::endl <<
            "\tfile: " << file      << std::endl <<
            "\tline: " << line      << std::endl <<
            std::endl;
#if defined(ASSERT_RAISES_EXCEPTION) && ASSERT_RAISES_EXCEPTION
    throw std::runtime_error("Assertion Error");
#else
    std::exit(EXIT_FAILURE);
#endif
}

inline void ecoevolity_assert_approx_eq_failed(char const * x,
        double val_x,
        const char * y,
        double val_y,
        char const * function,
        char const * file,
        long line) {
    std::cerr << std::endl <<
            "Approximately equal assertion failed:" << std::endl <<
            "\t" << x << " (" << std::fixed << std::setprecision(20) << val_x <<
                    ") approximately equal to " << y << " (" <<
                    std::fixed << std::setprecision(20) << val_y << ")" <<
                    std::endl <<
            "\tfunc: " << function  << std::endl <<
            "\tfile: " << file      << std::endl <<
            "\tline: " << line      << std::endl <<
            std::endl;
#if defined(ASSERT_RAISES_EXCEPTION) && ASSERT_RAISES_EXCEPTION
    throw std::runtime_error("Assertion Error");
#else
    std::exit(EXIT_FAILURE);
#endif
}

#if defined(IGNORE_ECOEVOLITY_ASSERT) || defined(NDEBUG)
#   define ECOEVOLITY_ASSERT(expr)
#   define ECOEVOLITY_ASSERT_APPROX_EQUAL(x, y)
#else
#   define ECOEVOLITY_ASSERT(expr)  if (!(expr)) ecoevolity_assertion_failed((const char *)#expr, (const char *)__FUNCTION__, __FILE__, __LINE__)
#   define ECOEVOLITY_ASSERT_APPROX_EQUAL(x, y)  if (fabs(((x)-(y))/(x)) > 1.0e-6) std::cerr << std::fixed << std::setprecision(20) << (x) << ' ' << std::fixed << std::setprecision(20) << (y) << '\n'; if (fabs(((x)-(y))/(x)) > 1.0e-6) ecoevolity_assert_approx_eq_failed((const char *)#x, x, (const char *)#y, y, (const char *)__FUNCTION__, __FILE__, __LINE__)
#endif

#define ECOEVOLITY_NDEBUG_ASSERT(expr)  if (!(expr)) ecoevolity_assertion_failed((const char *)#expr, (const char *)__FUNCTION__, __FILE__, __LINE__)
#define ECOEVOLITY_NDEBUG_ASSERT_APPROX_EQUAL(x, y)  if (fabs(((x)-(y))/(x)) > 1.0e-6) std::cerr << std::fixed << std::setprecision(20) << (x) << ' ' << std::fixed << std::setprecision(20) << (y) << '\n'; if (fabs(((x)-(y))/(x)) > 1.0e-6) ecoevolity_assert_approx_eq_failed((const char *)#x, x, (const char *)#y, y, (const char *)__FUNCTION__, __FILE__, __LINE__)

#endif
