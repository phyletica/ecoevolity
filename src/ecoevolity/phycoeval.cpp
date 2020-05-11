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

#include "phycoeval.hpp"


void write_phy_splash(std::ostream& out) {
    std::string v = "Version ";
    v += PROJECT_DETAILED_VERSION;
    out << string_util::banner('=') << "\n" 
        << string_util::center("Phycoeval") << "\n"
        << string_util::center("Estimating phylogenetic coevality") << "\n\n"
        << string_util::center("Part of:") << "\n"
        << string_util::center(PROJECT_NAME) << "\n"
        << string_util::center(v) << "\n"
        << string_util::banner('=') << "\n";
}

void update_log_paths(
        std::string & tree_log_path,
        std::string & state_log_path,
        std::string & operator_log_path,
        unsigned int max_number_of_attempts) {
    unsigned int ntries = 0;
    while (true) {
        if (
                (! path::exists(tree_log_path))
                && (! path::exists(state_log_path))
                && (! path::exists(operator_log_path))
            ) {
            return;
        }
        increment_log_path(tree_log_path);
        increment_log_path(state_log_path);
        increment_log_path(operator_log_path);
        ++ntries;
        if (ntries > max_number_of_attempts) {
            throw EcoevolityError("Could not generate unique output files");
        }
    }
}

void increment_log_path(
        std::string & log_path) {
    std::vector<std::string> path_elements;
    std::pair<std::string, std::string> prefix_ext;
    std::string new_suffix;
    int run_number;

    prefix_ext = path::splitext(log_path);
    path_elements = string_util::split(prefix_ext.first, '-');
    run_number = std::stoi(path_elements.back());
    path_elements.pop_back();
    ++run_number;
    
    new_suffix = "-" + std::to_string(run_number) + prefix_ext.second;

    log_path = string_util::join(path_elements, "-") + new_suffix;
}
