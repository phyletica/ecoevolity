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

#ifndef ECOEVOLITY_COLLECTION_HPP
#define ECOEVOLITY_COLLECTION_HPP

#include "data.hpp"
#include "node.hpp"
#include "tree.hpp"
#include "parameter.hpp"
#include "error.hpp"
#include "assert.hpp"

class ComparisonPopulationTreeCollection {
    private:
        std::vector<ComparisonPopulationTree> trees_;
        std::vector<PositiveRealParameter> node_heights_;
        std::vector<unsigned int> node_height_indices_;
        LogProbabilityDensity log_likelihood_ = LogProbabilityDensity(0.0);
        LogProbabilityDensity log_prior_density_ = LogProbabilityDensity(0.0);
        bool use_multithreading_ = false;

        double compute_log_prior_density_of_node_heights();

    public:
        void store_state();
        void restore_state();
        void compute_log_likelihood_and_prior();

        unsigned int get_number_of_comparisons() const {
            return this->trees_.size();
        }
};

#endif
