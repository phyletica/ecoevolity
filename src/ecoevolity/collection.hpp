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

    friend class DirichletProcessGibbsSampler;

    protected:
        std::vector<ComparisonPopulationTree> trees_;
        std::vector< std::shared_ptr<PositiveRealParameter> > node_heights_;
        std::vector<unsigned int> node_height_indices_;
        LogProbabilityDensity log_likelihood_ = LogProbabilityDensity(0.0);
        LogProbabilityDensity log_prior_density_ = LogProbabilityDensity(0.0);
        std::shared_ptr<ContinuousProbabilityDistribution> concentration_prior_ = std::make_shared<GammaDistribution>(1.0, 1.0);
        std::shared_ptr<PositiveRealParameter> concentration_ = std::make_shared<PositiveRealParameter>(this->concentration_prior_, 1.0);
        std::shared_ptr<ContinuousProbabilityDistribution> node_height_prior_ = std::make_shared<ExponentialDistribution>(100.0);
        bool use_multithreading_ = false;

    public:
        void store_state();
        void restore_state();
        void compute_log_likelihood_and_prior();

        unsigned int get_number_of_trees() const {
            return this->trees_.size();
        }

        unsigned int get_number_of_trees_mapped_to_height(
                unsigned int height_index) const;

        unsigned int get_height_index(unsigned int tree_index) const;
        std::shared_ptr<PositiveRealParameter> get_height_parameter(unsigned int tree_index) const;
        double get_height(unsigned int tree_index) const;

        double get_concentration() const;

        unsigned int get_number_of_auxiliary_heights() const;

        std::vector<unsigned int> get_other_height_indices(
                unsigned int tree_index) const;


        // if the new index is not the same as old
        //      - change index in node_height_indices_ for tree_index
        //      - set_height_parameter for tree accordingly
        // set tree's likelihood value and make clean
        void remap_tree(unsigned int tree_index,
                        unsigned int height_index,
                        double log_likelihood);


        // Add new neight parameter to node_heights_
        // change index in node_height_indices_ for tree
        // set_height_parameter for tree accordingly
        // set tree's likelihood value and make clean
        void map_tree_to_new_height(
                unsigned int tree_index,
                double height,
                double log_likelihood);

};

#endif
