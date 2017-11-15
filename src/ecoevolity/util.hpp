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

#ifndef ECOEVOLITY_UTIL_HPP
#define ECOEVOLITY_UTIL_HPP

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <stdexcept>

#include "assert.hpp"

/**
 * Function for accessing map elements.
 *
 * If the key does not exist, an out_of_range error is thrown.
 * From:
 *      <http://jeetworks.org/safe-and-const-correct-stdmap-access-in-c-stl/>
 *      2015 Jeet Sukumaran <http://jeetworks.org>
 */
template <typename T>
const typename T::value_type::second_type& map_at(
        const T& container,
        const typename T::value_type::first_type key) {
    typename T::const_iterator it = container.find(key);
    if (it == container.end()) {
        throw std::out_of_range("Key not found");
    }
    return it->second;
}

template <typename T1, typename T2>
bool second_of_pair_is_greater(
        const std::pair<T1, T2> & a,
        const std::pair<T1, T2> & b) {
    return a.second > b.second;
}

template <typename T1, typename T2>
bool second_of_pair_is_lesser(
        const std::pair<T1, T2> & a,
        const std::pair<T1, T2> & b) {
    return a.second < b.second;
}

template <typename T1, typename T2>
bool first_of_pair_is_greater(
        const std::pair<T1, T2> & a,
        const std::pair<T1, T2> & b) {
    return a.first > b.first;
}

template <typename T1, typename T2>
void sort_pairs(std::vector< std::pair<T1, T2> > & pairs,
        bool sort_by_first = true,
        bool reverse = false) {
    if (sort_by_first) {
        if (reverse) {
            std::sort(pairs.begin(), pairs.end(), first_of_pair_is_greater<T1, T2>);
        }
        else {
            std::sort(pairs.begin(), pairs.end());
        }
    }
    else {
        if (reverse) {
            std::sort(pairs.begin(), pairs.end(), second_of_pair_is_greater<T1, T2>);
        }
        else {
            std::sort(pairs.begin(), pairs.end(), second_of_pair_is_lesser<T1, T2>);
        }
    }
}

#endif
