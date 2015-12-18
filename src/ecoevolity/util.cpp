#include "util.hpp"

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
