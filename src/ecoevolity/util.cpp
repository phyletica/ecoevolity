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

#include "util.hpp"

std::vector<std::string> & split(
        const std::string &s,
        char delimiter,
        std::vector<std::string> & elements) {
    std::stringstream str_stream(s);
    std::string item;
    while (std::getline(str_stream, item, delimiter)) {
        elements.push_back(item);
    }
    return elements;
}

std::vector<std::string> split(
        const std::string &s,
        char delimiter) {
    std::vector<std::string> elements;
    split(s, delimiter, elements);
    return elements;
}

