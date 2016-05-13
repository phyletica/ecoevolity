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

#include "error.hpp"

EcoevolityBaseError::EcoevolityBaseError(
        const std::string & name,
        const std::string & message) {
    std::ostringstream oss;
    oss << name << ": " << message << std::endl;
    this->error_report_ = oss.str();
}

EcoevolityBaseError::EcoevolityBaseError(
        const std::string & name,
        const std::string & message,
        const std::string & file_name) {
    std::ostringstream oss;
    oss << name << ":" << std::endl
        << "    File: " << file_name << std::endl
        << "    Error: " << message << std::endl;
    this->error_report_ = oss.str();
}

EcoevolityBaseError::EcoevolityBaseError(
        const std::string & name,
        const std::string & message,
        unsigned int line_number) {
    std::ostringstream oss;
    oss << name << ":" << std::endl
        << "    Line: " << line_number << std::endl
        << "    Error: " << message << std::endl;
    this->error_report_ = oss.str();
}

EcoevolityBaseError::EcoevolityBaseError(
        const std::string & name,
        const std::string & message,
        const std::string & file_name,
        unsigned int line_number) {
    std::ostringstream oss;
    oss << name << ":" << std::endl
        << "    File: " << file_name << std::endl
        << "    Line: " << line_number << std::endl
        << "    Error: " << message << std::endl;
    this->error_report_ = oss.str();
}

EcoevolityBaseError::EcoevolityBaseError(
        const std::string & name,
        const std::string & message,
        const std::string & file_name,
        const std::string & taxon_label,
        unsigned int character_index) {
    std::ostringstream oss;
    oss << name << ":" << std::endl
        << "    File: " << file_name << std::endl
        << "    Taxon: " << taxon_label << std::endl
        << "    Character: " << character_index + 1 << std::endl
        << "    Error: " << message << std::endl;
    this->error_report_ = oss.str();
}
