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

