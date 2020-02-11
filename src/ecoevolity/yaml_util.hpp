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

#ifndef ECOEVOLITY_YAML_UTIL_HPP
#define ECOEVOLITY_YAML_UTIL_HPP

#include "yaml-cpp/yaml.h"

class YamlCppUtils {

    public:

        static std::string get_node_type(const YAML::Node& node) {
            switch(node.Type()) {
                case YAML::NodeType::Null:
                    return "Null";
                case YAML::NodeType::Scalar:
                    return "Scalar";
                case YAML::NodeType::Sequence:
                    return "Sequence";
                case YAML::NodeType::Map:
                    return "Map";
                case YAML::NodeType::Undefined:
                    return "Undefined";
                default:
                    return "Not recognized";
            }
        }
};

#endif
