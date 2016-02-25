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

void PopulationTree::compute_leaf_partials(
        unsigned int pattern_index,
        PopulationNode * node) {
    unsigned int pop_idx = this->data_.get_population_index(node->get_label());
    unsigned int red_allele_count = this->data_.get_red_allele_count(pattern_index, pop_idx));
    unsigned int allele_count = this->data_.get_allele_count(pattern_index, pop_idx);
    if ((this->data_.markers_are_dominant()) && (red_allele_count > 0)) {
        BiallelicPatternProbabilityMatrix m(allele_count);
        unsigned int n = allele_count / 2;
        unsigned int n_reds = red_allele_count / 2;
        double p_r_k_n = 1.0;
        for (unsigned int r = 1; r <= n_reds; ++r) {
            p_r_k_n = (p_r_k_n * 2.0 * (n - r + 1.0)) / ((2.0 * n) - r + 1.0);
        }
        for (unsigned int k = n_reds; k <= red_allele_count; ++k) {
            if (k > n_reds) {
                p_r_k_n = (p_r_k_n * ((2.0 * n_reds) - k + 1) * k) /
                          (2.0 * ( k - n_reds) * ((2.0 * n) - k + 1.0));
            }
            m.set_pattern_probability(allele_count, k, p_r_k_n);
        }
        return;
    }
    BiallelicPatternProbabilityMatrix m(allele_count, red_allele_count);
    node->copy_bottom_pattern_probs(m);
    return;
}

void PopulationTree::compute_top_of_branch_partials(
        unsigned int pattern_index,
        PopulationNode * node) {
    if (node.get_allele_count() == 0) {
        node.copy_top_pattern_probs(node.get_bottom_pattern_probs());
        return;
    }

    // TODO
    BiallelicPatternProbabilityMatrix m = MatrixExponentiator.expQTtx(...);
    node.copy_pattern_probs(m);
}

void PopulationTree::compute_internal_partials(
        unsigned int pattern_index,
        PopulationNode * node) {
    if (node->get_number_of_children == 1) {
        node.copy_bottom_pattern_probs(node->get_child(0)->get_top_pattern_probs());
        return;
    }
    if (node->get_child(0)->get_allele_count() == 0) {
        node.copy_bottom_pattern_probs(node->get_child(1)->get_top_pattern_probs());
        return;
    }
    if (node->get_child(1)->get_allele_count() == 0) {
        node.copy_bottom_pattern_probs(node->get_child(0)->get_top_pattern_probs());
        return;
    }
    unsigned int allele_count_child1 = node->get_child(0)->get_allele_count();
    unsigned int allele_count_child2 = node->get_child(1)->get_allele_count();

    std::vector<double> pattern_probs_child1 = node->get_child(0)->get_top_pattern_probs()->get_pattern_prob_matrix();
    std::vector<double> pattern_probs_child2 = node->get_child(1)->get_top_pattern_probs()->get_pattern_prob_matrix();

    for (int n = 1; n <= allele_count_child1; ++n) {
        double b_nr = 1.0;
        for (int r = 0; r <= n; ++r) {
            pattern_probs_child1.at(((n*(n+1))/2)-1+r) *= b_nr;
            b_nr *= ((double)n - r)/(r+1);
        }
    }

    for (n = 1; n<= allele_count_child2; ++n) {
        b_nr = 1.0;
        for (r = 0; r <= n; ++r) {
            pattern_probs_child2.at(((n*(n+1))/2)-1+r) *= b_nr;
            b_nr *= ((double)n - r)/(r+1);
        }
    }

    unsigned int allele_count = allele_count_child1 + allele_count_child2;
    std::vector<double> pattern_probs; 
    pattern_probs_.assign(
            ((((allele_count + 1) * (allele_count + 2))/2) - 1),
            0);
    for (int n1 = 1; n1 <= allele_count_child1; ++n1) {
        for (int r1 = 0; r1 <= n1; ++r1) {
            double f11 = pattern_probs_child1.at(n1*(n1+1)/2-1+r1);
            for (int n2 = 1; n2 <= allele_count_child2; ++n2) {
                for (int r2 = 0; r2 <= n2; ++r2) {
                    pattern_probs.at((n1+n2)*(n1+n2+1)/2-1+(r1+r2)) += f11 * pattern_probs_child2.at(n2*(n2+1)/2-1+r2);
                }
            }
        }
    }

    for (n = 1; n <= allele_count; ++n) {
        b_nr = 1.0;
        for (r = 0; r <= n; ++r) {
            double f_nr = pattern_probs.at(n*(n+1)/2-1+r);
            f_nr /= b_nr;
            // TODO: better way to fix this?
            f_nr = std::max(f_nr, 0.0);
            pattern_probs.at(n*(n+1)/2-1+r) = f_nr;
            b_nr *= ((double)n - r)/(r+1);

        }
    }
    BiallelicPatternProbabilityMatrix m(allele_count, pattern_probs);
    node.copy_bottom_pattern_probs(m);
}

void PopulationTree::compute_pattern_partials(
        unsigned int pattern_index,
        PopulationNode * node) {
    if (node->is_leaf()) {
        this->compute_leaf_partials(pattern_index, node);
    }
    else if (node.get_number_of_children() == 1) {
        compute_pattern_partials(pattern_index, node->get_child(0));
        compute_top_of_branch_partials(pattern_index, node->get_child(0));
        compute_internal_partials(pattern_index, node);
    }
    else if (node.get_number_of_children() == 2) {
        compute_pattern_partials(pattern_index, node->get_child(0));
        compute_pattern_partials(pattern_index, node->get_child(1));
        compute_top_of_branch_partials(pattern_index, node->get_child(0));
        compute_top_of_branch_partials(pattern_index, node->get_child(1));
        compute_internal_partials(pattern_index, node);
    }
    else {
        throw EcoevolityError(
            "PopulationTree::compute_pattern_probability(); unexpected number of children");
    }
}

void PopulationTree::compute_pattern_likelihood(unsigned int pattern_index) {
    this->compute_pattern_partials(pattern_idx, this->root_);
    double pattern_likelihood = this->compute_root_likelihood();
    // TODO: store likelihood in pattern_probs_
}

void PopulationTree::compute_pattern_likelihoods() {
    for (unsigned int pattern_idx = 0;
            pattern_idx < this->data_.get_number_of_patterns();
            ++pattern_idx) {
        this->compute_pattern_likelihood(pattern_idx, this->root_);
    }
}
