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

#include <thread>

#include "data.hpp"
#include "node.hpp"
#include "tree.hpp"
#include "parameter.hpp"
#include "settings.hpp"
#include "operator_schedule.hpp"
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
        std::string state_log_path_ = "ecoevolity-state-run-1.log";
        std::string operator_log_path_ = "ecoevolity-operator-run-1.log";
        bool use_multithreading_ = false;
        unsigned int logging_precision_ = 18;
        std::string logging_delimiter_ = "\t";

        void init_trees(
                const std::vector<ComparisonSettings> & comparison_settings,
                RandomNumberGenerator & rng);

        void compute_tree_partials();
        void compute_tree_partials_threaded();
        void make_trees_clean();
        void make_trees_dirty();

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

        void write_state_log_header(std::ostream& out,
                bool short_summary = false) const;
        void log_state(std::ostream& out,
                unsigned int generation_index,
                bool short_summary = false) const;

        void update_log_paths(unsigned int max_number_of_attempts = 10000);
        void increment_log_paths();

    public:
        ComparisonPopulationTreeCollection() { }
        ComparisonPopulationTreeCollection(
                const CollectionSettings & settings,
                RandomNumberGenerator & rng
                );
        virtual ~ComparisonPopulationTreeCollection() { }
        void store_state();
        void restore_state();
        void compute_log_likelihood_and_prior(bool compute_partials = true);

        void ignore_data() {
            for (unsigned int i = 0; i < this->trees_.size(); ++i) {
                this->trees_.at(i).ignore_data();
            }
        }
        void use_data() {
            for (unsigned int i = 0; i < this->trees_.size();  ++i) {
                this->trees_.at(i).use_data();
            }
        }

        bool using_dpp() const {
            return this->operator_schedule_.using_dpp();
        }

        unsigned int get_logging_precision() const {
            return this->logging_precision_;
        }
        void set_logging_precision(unsigned int precision) {
            this->logging_precision_ = precision;
        }
        const std::string& get_logging_delimiter() const {
            return this->logging_delimiter_;
        }
        void set_logging_delimiter(const std::string& delimiter) {
            this->logging_delimiter_ = delimiter;
        }

        const ComparisonPopulationTree& get_tree(
                unsigned int tree_index) const {
            return this->trees_.at(tree_index);
        }


        unsigned int get_number_of_trees() const {
            return this->trees_.size();
        }

        unsigned int get_number_of_events() const {
            return this->node_heights_.size();
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
            return this->node_heights_.at(tree_index)->get_value();
        }

        double get_concentration() const {
            return this->concentration_->get_value();
        }
        void set_concentration(double value) {
            this->concentration_->set_value(value);
        }

        std::vector<unsigned int> get_other_height_indices(
                unsigned int tree_index) const;

        std::string get_state_log_path() const {
            return this->state_log_path_;
        }
        std::string get_operator_log_path() const {
            return this->operator_log_path_;
        }
        void set_state_log_path(const std::string& path) {
            this->state_log_path_ = path;
        }
        void set_operator_log_path(const std::string& path) {
            this->operator_log_path_ = path;
        }

        void mcmc(RandomNumberGenerator& rng,
                unsigned int chain_length,
                unsigned int sample_frequency);
};

#endif
