/******************************************************************************
 * Copyright (C) 2016 Jamie R. Oaks.
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


std::vector<std::string> & split(
        const std::string &s,
        char delimiter,
        std::vector<std::string> & elements);

std::vector<std::string> split(
        const std::string &s,
        char delimiter);

#endif

