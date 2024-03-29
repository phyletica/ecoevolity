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

#ifndef ECOEVOLITY_STRING_UTIL_HPP
#define ECOEVOLITY_STRING_UTIL_HPP

#include <string>
#include <sstream>
#include <vector>
#include <map>

#include "assert.hpp"

namespace string_util {

inline std::vector<std::string> & split(
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

inline std::vector<std::string> split(
        const std::string &s,
        char delimiter) {
    std::vector<std::string> elements;
    split(s, delimiter, elements);
    return elements;
}

inline std::string join(const std::vector<std::string>& components,
                        const std::string delimiter = "") {
    std::ostringstream ss;
    for (unsigned int i = 0; i < components.size(); ++i) {
        if (i == 0) {
            ss << components.at(i);
        }
        else {
            ss << delimiter << components.at(i);
        }
    }
    return ss.str();
}

inline std::string pad_int(unsigned int n, unsigned int len) {
    ECOEVOLITY_ASSERT(len > 0);
    std::string r = std::to_string(n);
    if (r.size() >= len) {
        return r;
    }
    return std::string(len - r.size(), '0') + r;
}

inline std::string get_indent(unsigned int level = 1) {
    return std::string(4 * level, ' ');
}

inline std::string center(const std::string& s, unsigned int page_width = 70) {
    int w = (int)((page_width - s.length()) / 2);
    std::string indent(w, ' ');
    return indent + s;
}

inline std::string banner(char c, unsigned int page_width = 70) {
    return std::string(page_width, c);
}

inline std::string rstrip(
        const std::string& s,
        const std::string& delimiters = " \f\n\r\t\v" ) {
    std::string::size_type last_nws = s.find_last_not_of(delimiters);
    if (last_nws >= s.size()) {
        return std::string("");
    }
    return s.substr(0, last_nws + 1);
}

inline std::string lstrip(
        const std::string& s,
        const std::string& delimiters = " \f\n\r\t\v" ) {
    std::string::size_type first_nws = s.find_first_not_of(delimiters);
    if (first_nws >= s.size()) {
        return std::string("");
    }
    return s.substr(first_nws);
}

inline std::string strip(
        const std::string& s,
        const std::string& delimiters = " \f\n\r\t\v" ) {
    return lstrip(rstrip(s, delimiters), delimiters);
}

inline bool startswith(
        const std::string& s,
        const std::string& match) {
    return ((s.size() >= match.size()) && std::equal(
                match.begin(),
                match.end(),
                s.begin()));
}

inline void parse_map(
        const std::string &s,
        std::map<std::string, std::string> & key_value_map,
        char item_delimiter = ',',
        char key_value_delimiter = '=') {
    std::vector<std::string> map_items;
    split(s, item_delimiter, map_items);
    for (auto item : map_items) {
        std::vector<std::string> key_value_pair;
        split(item, key_value_delimiter, key_value_pair);
        ECOEVOLITY_ASSERT(key_value_pair.size() == 2);
        key_value_map[strip(key_value_pair.at(0))] = strip(key_value_pair.at(1));
    }
}

inline std::map<std::string, std::string> parse_map(
        const std::string &s,
        char item_delimiter = ',',
        char key_value_delimiter = '=') {
    std::map<std::string, std::string> key_value_map;
    parse_map(s, key_value_map, item_delimiter, key_value_delimiter);
    return key_value_map;
}

} // namespace string_util

#endif
