#ifndef COEVOLITY_UTIL_HPP
#define COEVOLITY_UTIL_HPP

#include <iostream>
#include <map>

/**
 * Function for accessing map elements.
 *
 * If the key does not exist, an out_of_range error is thrown.
 * From:
 *      http://jeetworks.org/safe-and-const-correct-stdmap-access-in-c-stl/
 *      2015 Jeet Sukumaran <http://jeetworks.org
 */
template <typename T>
const typename T::value_type::second_type& map_at(
        const T& container,
        const typename T::value_type::first_type key);

#endif
