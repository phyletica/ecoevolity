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

#include "general_tree.hpp"

void GeneralPopulationTree::init_node_heights() {
    std::vector< std::shared_ptr<PopulationNode> > internal_nodes = this->root_->get_internal_nodes();
    for (unsigned int i = 0; i < internal_nodes.size(); ++i) {
        bool exists = false;
        // Check if we already have a pointer to the same place in memory
        for (unsigned int j = 0; j < this->node_heights_.size(); ++j) {
            if (internal_nodes.at(i)->get_height_parameter() == this->node_heights_.at(j)) {
                exists = true;
            }
        }
        if (! exists) {
            this->node_heights_.push_back(internal_nodes.at(i)->get_height_parameter());
        }
    }
}
