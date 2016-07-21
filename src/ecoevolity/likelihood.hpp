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

#ifndef ECOEVOLITY_LIKELIHOOD_HPP
#define ECOEVOLITY_LIKELIHOOD_HPP

#include <algorithm>
#include <memory>
#include <future>

#include "data.hpp"
#include "node.hpp"
#include "matrix.hpp"
#include "error.hpp"

void compute_leaf_partials(
        PopulationNode& node,
        const BiallelicData& data,
        int pattern_index
        );

void compute_top_of_branch_partials(
        PopulationNode& node,
        double u,
        double v,
        double mutation_rate, 
        double ploidy
        );

void compute_internal_partials(
        PopulationNode& node);

void compute_pattern_partials(
        PopulationNode& node,
        const BiallelicData& data,
        int pattern_index,
        double u,
        double v,
        double mutation_rate,
        double ploidy
        );

std::vector< std::vector<double> > compute_root_probabilities(
        const PopulationNode& root,
        u,
        v,
        mutation_rate,
        ploidy
        );

double compute_root_likelihood(
        const PopulationNode& root,
        u,
        v,
        mutation_rate,
        ploidy
        );

double compute_pattern_likelihood(
        PopulationNode& root,
        const BiallelicData& data,
        int pattern_index,
        double u,
        double v,
        double mutation_rate,
        double ploidy
        );

void compute_constant_pattern_likelihoods(
        PopulationNode& root,
        const BiallelicData& data,
        double u,
        double v,
        double mutation_rate,
        double ploidy,
        bool state_frequencies_are_constrained,
        double& all_red_pattern_likelihood,
        double& all_green_pattern_likelihood
        );


double get_log_likelihood_for_pattern_range(
        PopulationNode& root,
        const BiallelicData& data,
        unsigned int start_pattern_index,
        unsigned int stop_pattern_index,
        double u,
        double v,
        double mutation_rate,
        double ploidy
        );

double get_log_likelihood(
        PopulationNode& root,
        const BiallelicData& data,
        double u,
        double v,
        double mutation_rate,
        double ploidy,
        bool state_frequencies_are_constrained,
        bool constant_sites_removed,
        double& all_red_pattern_likelihood,
        double& all_green_pattern_likelihood,
        unsigned int nthreads = 1
        );

#endif
