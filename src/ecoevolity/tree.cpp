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

#include "tree.hpp"

PopulationTree::PopulationTree(
        const std::string path, 
        const char population_name_delimiter,
        const bool population_name_is_prefix,
        const bool genotypes_are_diploid,
        const bool markers_are_dominant,
        const bool validate) {
    this->data_ = BiallelicData(
            path,
            population_name_delimiter,
            population_name_is_prefix,
            genotypes_are_diploid,
            markers_are_dominant,
            validate);
    if (this->data_.get_number_of_populations > 2) {
        throw EcoevolityError("PopulationTree(); does not support > 2 tips");
    }

    std::string root_label = ""
    for (unsigned int pop_idx = 0;
            pop_idx < this->data_.get_number_of_populations();
            ++pop_idx) {
        root_label += this->data_.get_population_label(pop_idx);
    }
    this->root_ = new PopulationNode(root_label, 0.01);
    for (pop_idx = 0;
            pop_idx < this->data_.get_number_of_populations();
            ++pop_idx) {
        this->root_.add_child(new PopulationNode(
                this->data_.get_population_label(pop_idx),
                0.0,
                this->data_.get_max_allele_count(pop_idx)));
    }
    this->root_.resize_all();

    this->pattern_probs_.assign(this->data_.get_number_of_patterns(), 0.0);
}

void PopulationTree::compute_pattern_probability(unsigned int pattern_index, PopulationNode * node) {
    if (node->is_leaf()) {
        this->compute_leaf_probability(pattern_index, node);
    }
    else if (node.get_number_of_children() == 1) {
        compute_pattern_probability(pattern_index, node->get_child(0));
        compute_top_of_branch_probability(pattern_index, node->get_child(0));
        compute_internal_probability(pattern_index, node->get_child(0));
    }
    else if (node.get_number_of_children() == 2) {
        compute_pattern_probability(pattern_index, node->get_child(0));
        compute_pattern_probability(pattern_index, node->get_child(1));
        compute_top_of_branch_probability(pattern_index, node->get_child(0));
        compute_top_of_branch_probability(pattern_index, node->get_child(1));
        compute_internal_probability(pattern_index, node->get_child(0), node->get_child(1));
    }
    else {
        throw EcoevolityError(
                "PopulationTree::compute_pattern_probability(); unexpected number of children");
    }
}

void PopulationTree::compute_pattern_probabilities() {
    for (unsigned int pattern_idx = 0;
            pattern_idx < this->data_.get_number_of_patterns();
            ++pattern_idx) {
        this->compute_pattern_probability(pattern_idx, this->root_);
    }
    this->compute_root_likelihood();
}
