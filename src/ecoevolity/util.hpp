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

