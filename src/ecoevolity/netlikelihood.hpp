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

#ifndef ECOEVOLITY_NETLIKELIHOOD_HPP
#define ECOEVOLITY_NETLIKELIHOOD_HPP

#include <algorithm>
#include <memory>

#ifdef BUILD_WITH_THREADS
#include <future>
#endif

#include "netnode.hpp"
#include "matrix.hpp"
#include "error.hpp"
#include "assert.hpp"
#include "math_util.hpp"

namespace netlikelihood {

void compute_leaf_partials(
        PopulationNetNode& node,
        const unsigned int red_allele_count,
        const unsigned int allele_count,
        const bool markers_are_dominant
        );

void compute_top_of_branch_partials(
        PopulationNetNode& node,
        const double u,
        const double v,
        const double mutation_rate, 
        const double ploidy
        );

void split_top_of_branch_partials(
        const unsigned int max_num_alleles,
        const BiallelicPatternProbabilityMatrix & top_child_partials,
        const double prob_to_parent1,
        const double prob_to_parent2,
        BiallelicPatternProbabilityMatrix & bottom_parent1_partials,
        BiallelicPatternProbabilityMatrix & bottom_parent2_partials);

void merge_top_of_branch_partials(
        const unsigned int allele_count_child1,
        const unsigned int allele_count_child2,
        const double prob_no_alleles_child1,
        const double prob_no_alleles_child2,
        std::vector<double> & top_partials_child1,
        std::vector<double> & top_partials_child2,
        unsigned int & merged_allele_count,
        std::vector<double> & merged_pattern_probs,
        double & merged_prob_no_alleles,
        const bool do_hypergeom_scaling = true);

void compute_internal_partials(
        PopulationNetNode& node,
        std::set< std::shared_ptr<const PopulationNetNode> > & retic_nodes_visited);

void compute_pattern_partials(
        PopulationNetNode& node,
        const std::vector<unsigned int>& red_allele_counts,
        const std::vector<unsigned int>& allele_counts,
        const double u,
        const double v,
        const double mutation_rate,
        const double ploidy,
        const bool markers_are_dominant,
        std::set< std::shared_ptr<const PopulationNetNode> > & retic_nodes_visited
        );

std::vector< std::vector<double> > compute_root_probabilities(
        const PopulationNetNode& root,
        const double u,
        const double v,
        const double mutation_rate,
        const double ploidy
        );

double compute_root_likelihood(
        const PopulationNetNode& root,
        const double u,
        const double v,
        const double mutation_rate,
        const double ploidy
        );

double compute_pattern_likelihood(
        PopulationNetNode& root,
        const std::vector<unsigned int>& red_allele_counts,
        const std::vector<unsigned int>& allele_counts,
        const double u,
        const double v,
        const double mutation_rate,
        const double ploidy,
        const bool markers_are_dominant
        );

void compute_constant_pattern_log_likelihood_correction(
        PopulationNetNode& root,
        const std::vector< std::vector<unsigned int> > & unique_allele_counts,
        const std::vector<unsigned int> & unique_allele_count_weights,
        const double u,
        const double v,
        const double mutation_rate,
        const double ploidy,
        const bool markers_are_dominant,
        const bool state_frequencies_are_constrained,
        double& constant_log_likelihood_correction
        );

double get_log_likelihood_for_pattern_range(
        PopulationNetNode& root,
        const std::vector< std::vector<unsigned int> >& red_allele_count_matrix,
        const std::vector< std::vector<unsigned int> >& allele_count_matrix,
        const std::vector<unsigned int>& pattern_weights,
        const unsigned int start_index,
        const unsigned int stop_index,
        const double u,
        const double v,
        const double mutation_rate,
        const double ploidy,
        const bool markers_are_dominant
        );

double get_log_likelihood(
        PopulationNetNode& root,
        const std::vector< std::vector<unsigned int> >& red_allele_count_matrix,
        const std::vector< std::vector<unsigned int> >& allele_count_matrix,
        const std::vector<unsigned int>& pattern_weights,
        const std::vector< std::vector<unsigned int> > & unique_allele_counts,
        const std::vector<unsigned int> & unique_allele_count_weights,
        const double u,
        const double v,
        const double mutation_rate,
        const double ploidy,
        const bool markers_are_dominant,
        const bool state_frequencies_are_constrained,
        const bool constant_sites_removed,
        double& constant_log_likelihood_correction,
        unsigned int nthreads = 1
        );
}

#endif
