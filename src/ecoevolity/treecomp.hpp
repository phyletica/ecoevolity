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

#ifndef ECOEVOLITY_TREECOMP_HPP
#define ECOEVOLITY_TREECOMP_HPP

#include <cmath>


namespace treecomp {

template<class TreeType>
inline double euclidean_distance(
        const TreeType & tree1,
        const TreeType & tree2,
        const bool resize_splits = false) {
    std::map<Split, double> split_length_map1 = tree1.get_split_length_map(
            resize_splits);
    std::map<Split, double> split_length_map2 = tree2.get_split_length_map(
            resize_splits);
    std::set<Split> all_splits;
    for (auto split_len : split_length_map1) {
        all_splits.insert(split_len.first);
    }
    for (auto split_len : split_length_map2) {
        all_splits.insert(split_len.first);
    }
    double sum_squared_diffs = 0.0;
    for (auto split : all_splits) {
        double l1 = 0.0;
        double l2 = 0.0;
        if (split_length_map1.count(split) > 0) {
            l1 = split_length_map1.at(split);
        }
        if (split_length_map2.count(split) > 0) {
            l2 = split_length_map2.at(split);
        }
        sum_squared_diffs += std::pow((l1 - l2), 2);
    }
    return std::sqrt(sum_squared_diffs);
}

} // treecomp

#endif
