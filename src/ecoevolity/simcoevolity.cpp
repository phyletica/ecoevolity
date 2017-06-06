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

#include "simcoevolity.hpp"


void write_sim_splash(std::ostream& out) {
    std::string v = "Version ";
    v += PROJECT_DETAILED_VERSION;
    out << string_util::banner('=') << "\n" 
        << string_util::center("Simcoevolity") << "\n"
        << string_util::center("Simulating evolutionary coevality") << "\n"
        << string_util::center("Part of:") << "\n"
        << string_util::center(PROJECT_NAME) << "\n"
        << string_util::center(v) << "\n"
        << string_util::banner('=') << "\n";
}

void check_output_path(const std::string& path) {
    if (path::exists(path)) {
        std::ostringstream message;
        message << "ERROR: The simcoevolity output file \'"
                << path
                << "\' already exists!\n"
                << "Please use a different output file or clear out previous "
                << "simulations\n";
        throw EcoevolityError(message.str());
    }
}
