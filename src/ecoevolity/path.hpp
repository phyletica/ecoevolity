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

#ifndef ECOEVOLITY_PATH_HPP
#define ECOEVOLITY_PATH_HPP

#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>

#include "string_util.hpp"

#ifdef WINDOWS
    static const std::string PATH_SEP = "\\";
#else
    static const std::string PATH_SEP = "/";
#endif

namespace path {

inline std::string dirname(const std::string& path) {
    std::string::size_type last_sep = path.find_last_of(PATH_SEP);
    if (last_sep == std::string::npos) {
        return std::string("");
    }
    if (last_sep >= path.size()) {
        return std::string("");
    }
    return path.substr(0, last_sep);
}

inline std::string basename(const std::string& path) {
    std::string::size_type last_sep = path.find_last_of(PATH_SEP);
    if (last_sep == std::string::npos) {
        return path;
    }
    if (last_sep >= path.size()) {
        return std::string("");
    }
    return path.substr(last_sep + 1);
}

inline bool isabs(const std::string& path) {
    if (path.size() < 1) {
        return false;
    }
    std::string::size_type first_sep = path.find_first_of(PATH_SEP);
    if (first_sep == 0) {
        return true;
    }
    return false;
}

inline std::string join(const std::string& parent, const std::string& child) {
    if (isabs(child)) {
        return child;
    }
    std::string p = string_util::rstrip(parent, PATH_SEP);
    return p + PATH_SEP + child;
}

inline bool exists(const std::string& path) {
    struct stat info;
    return (stat(path.c_str(), &info) == 0);
}

inline bool isdir(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        return false;
    }
    if (! (info.st_mode & S_IFDIR)) {
        return false;
    }
    return true;
}

inline bool isfile(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        return false;
    }
    if (! (info.st_mode & S_IFREG)) {
        return false;
    }
    return true;
}

} // namespace path

#endif
