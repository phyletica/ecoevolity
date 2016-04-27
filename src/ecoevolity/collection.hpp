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
#include "settings.hpp"
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
        std::shared_ptr<PositiveRealParameter> concentration_;
        std::shared_ptr<ContinuousProbabilityDistribution> node_height_prior_;
        OperatorSchedule operator_schedule_;
        bool use_multithreading_ = false;

        void compute_tree_partials();
        void compute_tree_partials_threaded();
        void make_trees_clean();

        void remap_tree(unsigned int tree_index,
                        unsigned int height_index,
                        double log_likelihood);

        void map_tree_to_new_height(
                unsigned int tree_index,
                double height,
                double log_likelihood);

        void add_height(
                double height,
                const std::vector<unsigned int>& mapped_tree_indices);

        void remove_height(unsigned int height_index);

    public:
        ComparisonPopulationTreeCollection() { }
        ComparisonPopulationTreeCollection(
                const CollectionSettings & settings,
                RandomNumberGenerator & rng
                );
        virtual ~ComparisonPopulationTreeCollection() { }
        void store_state();
        void restore_state();
        void compute_log_likelihood_and_prior();

        unsigned int get_number_of_trees() const {
            return this->trees_.size();
        }

        unsigned int get_number_of_trees_mapped_to_height(
                unsigned int height_index) const;
        unsigned int get_number_of_partners(
                unsigned int tree_index) const;

        unsigned int get_height_index(unsigned int tree_index) const {
            return this->node_height_indices_.at(tree_index);
        }

        std::shared_ptr<PositiveRealParameter> get_height_parameter(
                unsigned int height_index) const {
            return this->node_heights_.at(height_index);
        }

        double get_height(unsigned int tree_index) const {
            return this->node_heights_.at(tree_index).get_value();
        }

        double get_concentration() const {
            return this->concentration_.get_value();
        }

        unsigned int get_number_of_auxiliary_heights() const {
            return this->number_of_auxiliary_heights_;
        }

        std::vector<unsigned int> get_other_height_indices(
                unsigned int tree_index) const;


        void mcmc();
};

#endif
