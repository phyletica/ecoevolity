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

#include <unordered_set>

#include "data.hpp"
#include "node.hpp"
#include "tree.hpp"
#include "parameter.hpp"
#include "settings.hpp"
#include "operator_schedule.hpp"
#include "error.hpp"
#include "assert.hpp"
#include "rng.hpp"

class BaseComparisonPopulationTreeCollection {

    protected:
        std::vector< std::shared_ptr<PopulationTree> > trees_;
        std::vector< std::shared_ptr<PositiveRealParameter> > node_heights_;
        std::vector<double> stored_node_heights_;
        std::vector<unsigned int> node_height_indices_;
        std::vector<unsigned int> stored_node_height_indices_;
        LogProbabilityDensity log_likelihood_ = LogProbabilityDensity(0.0);
        LogProbabilityDensity log_prior_density_ = LogProbabilityDensity(0.0);
        std::shared_ptr<PositiveRealParameter> concentration_;
        std::shared_ptr<ContinuousProbabilityDistribution> node_height_prior_;
        OperatorSchedule operator_schedule_;
        std::string state_log_path_ = "ecoevolity-state-run-1.log";
        std::string operator_log_path_ = "ecoevolity-operator-run-1.log";
        unsigned int number_of_threads_ = 1;
        unsigned int logging_precision_ = 18;
        std::string logging_delimiter_ = "\t";

        void add_height(
                double height,
                const std::vector<unsigned int>& mapped_tree_indices);

        void remove_height(unsigned int height_index);

    public:
        BaseComparisonPopulationTreeCollection() { }

        void add_log_prefix(const std::string & prefix) {
            this->state_log_path_ = prefix + path::basename(this->state_log_path_);
            this->operator_log_path_ = prefix + path::basename(this->operator_log_path_);
        }

        void store_state();
        void restore_state();
        void store_model_state();
        void restore_model_state();
        void compute_log_likelihood_and_prior(bool compute_partials = true);

        void compute_tree_partials();
        void make_trees_clean();
        void make_trees_dirty();

        void set_number_of_threads(unsigned int n) {
            this->number_of_threads_ = n;
        }
        unsigned int get_number_of_threads() const {
            return this->number_of_threads_;
        }
        void ignore_data() {
            for (unsigned int i = 0; i < this->trees_.size(); ++i) {
                this->trees_.at(i)->ignore_data();
            }
        }
        void use_data() {
            for (unsigned int i = 0; i < this->trees_.size();  ++i) {
                this->trees_.at(i)->use_data();
            }
        }

        bool ignoring_data() const {
            return this->trees_.at(0)->ignoring_data();
        }

        bool using_dpp() const {
            return this->operator_schedule_.using_dpp();
        }

        bool using_reversible_jump() const {
            return this->operator_schedule_.using_reversible_jump();
        }

        bool sampling_models() const {
            return this->operator_schedule_.sampling_models();
        }

        double get_draw_from_node_height_prior(RandomNumberGenerator& rng) {
            return this->node_height_prior_->draw(rng);
        }

        double get_node_height_prior_mean() const {
            return this->node_height_prior_->get_mean();
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

        std::shared_ptr<PopulationTree> get_tree(
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

        std::vector<unsigned int> get_shared_event_indices() const;

        unsigned int get_height_index(unsigned int tree_index) const {
            return this->node_height_indices_.at(tree_index);
        }

        std::vector<unsigned int> get_standardized_height_indices() const;

        unsigned int get_largest_height_index() const;
        std::vector<unsigned int> get_height_indices_sans_largest() const;

        std::shared_ptr<PositiveRealParameter> get_height_parameter(
                unsigned int height_index) const {
            return this->node_heights_.at(height_index);
        }

        double get_log_likelihood() const {
            return this->log_likelihood_.get_value();
        }
        double get_stored_log_likelihood() const {
            return this->log_likelihood_.get_stored_value();
        }
        double get_log_prior_density() const {
            return this->log_prior_density_.get_value();
        }
        double get_stored_log_prior_density() const {
            return this->log_prior_density_.get_stored_value();
        }

        OperatorSchedule & get_operator_schedule() {
            return this->operator_schedule_;
        }
        void set_operator_schedule(OperatorSchedule& os) {
            this->operator_schedule_ = os;
        }

        double get_height(unsigned int height_index) const {
            return this->node_heights_.at(height_index)->get_value();
        }
        void set_height(unsigned int height_index, double height) const {
            this->node_heights_.at(height_index)->set_value(height);
            for (unsigned int i : this->get_indices_of_mapped_trees(height_index)) {
                this->trees_.at(i)->make_dirty();
            }
        }
        double get_height_of_tree(unsigned int tree_index) const {
            return this->get_height(this->get_height_index(tree_index));
        }
        void set_height_of_tree(unsigned int tree_index, double height) const {
            this->set_height(this->get_height_index(tree_index), height);
        }

        std::vector<unsigned int> get_indices_of_mapped_trees(
                unsigned int height_index) const;

        double get_nearest_smaller_height(
                unsigned int height_index) const;
        unsigned int get_nearest_smaller_height_index(
                unsigned int height_index,
                bool allow_smallest_index = false) const;
        unsigned int get_nearest_larger_height_index(
                unsigned int height_index,
                bool allow_largest_index = false) const;

        unsigned int get_distal_height_index_within_move(
                unsigned int starting_height_index,
                double delta_height) const;

        double get_concentration() const {
            return this->concentration_->get_value();
        }
        void set_concentration(double value) {
            this->concentration_->set_value(value);
        }
        bool concentration_is_fixed() const {
            return this->concentration_->is_fixed();
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

        void write_state_log_header(std::ostream& out,
                bool short_summary = false) const;
        void log_state(std::ostream& out,
                unsigned int generation_index,
                bool short_summary = false) const;

        
        std::vector< std::shared_ptr<OperatorInterface> > get_time_operators() const {
            return this->operator_schedule_.get_time_operators();
        }
        std::vector< std::shared_ptr<OperatorInterface> > get_time_operators(
                int tree_index) const {
            return this->operator_schedule_.get_time_operators(tree_index);
        }
        std::vector< std::shared_ptr<OperatorInterface> > get_tree_operators() const {
            return this->operator_schedule_.get_tree_operators();
        }
        std::vector< std::shared_ptr<OperatorInterface> > get_tree_operators(
                int tree_index) const {
            return this->operator_schedule_.get_tree_operators(tree_index);
        }
        std::vector< std::shared_ptr<OperatorInterface> > get_multivariate_time_operators() const {
            return this->operator_schedule_.get_multivariate_time_operators();
        }

        void mcmc(RandomNumberGenerator& rng,
                unsigned int chain_length,
                unsigned int sample_frequency);

        void write_summary(
                std::ostream& out,
                unsigned int indent_level = 0) const;

        void draw_heights_from_prior(RandomNumberGenerator& rng);
        void draw_from_prior(RandomNumberGenerator& rng);

        std::map<std::string, BiallelicData> simulate_biallelic_data_sets(
                RandomNumberGenerator& rng,
                bool validate = true) const;

        std::map<std::string, BiallelicData> simulate_complete_biallelic_data_sets(
                RandomNumberGenerator& rng,
                unsigned int locus_size = 1,
                bool validate = true) const;

        bool all_population_sizes_are_fixed() const {
            for (unsigned int i = 0; i < this->get_number_of_trees(); ++i) {
                if (! this->trees_.at(i)->population_sizes_are_fixed()) {
                    return false;
                }
            }
            return true;
        }
        bool all_mutation_rates_are_fixed() const {
            for (unsigned int i = 0; i < this->get_number_of_trees(); ++i) {
                if (! this->trees_.at(i)->mutation_rate_is_fixed()) {
                    return false;
                }
            }
            return true;
        }

        void remap_tree(unsigned int tree_index,
                        unsigned int height_index);
        void remap_tree(unsigned int tree_index,
                        unsigned int height_index,
                        double log_likelihood);
        unsigned int remap_trees(
                const std::vector<unsigned int>& tree_indices,
                unsigned int height_index);

        unsigned int merge_height(
                unsigned int height_index,
                unsigned int target_height_index);

        void map_tree_to_new_height(
                unsigned int tree_index,
                double height);
        void map_tree_to_new_height(
                unsigned int tree_index,
                double height,
                double log_likelihood);
        unsigned int map_trees_to_new_height(
                const std::vector<unsigned int>& tree_indices,
                double height);

        void set_node_height_indices(
                const std::vector<unsigned int>& indices,
                RandomNumberGenerator & rng);

        void update_log_paths(unsigned int max_number_of_attempts = 10000);
        void increment_log_paths();
};

class ComparisonPopulationTreeCollection: public BaseComparisonPopulationTreeCollection {

    public:
        ComparisonPopulationTreeCollection() : BaseComparisonPopulationTreeCollection() { }
        ComparisonPopulationTreeCollection(
                const CollectionSettings & settings,
                RandomNumberGenerator & rng,
                bool strict_on_constant_sites = true,
                bool strict_on_missing_sites = true,
                bool strict_on_triallelic_sites = true
                );

    protected:
        void init_trees(
                const std::vector<ComparisonSettings> & comparison_settings,
                RandomNumberGenerator & rng,
                bool strict_on_constant_sites = true,
                bool strict_on_missing_sites = true,
                bool strict_on_triallelic_sites = true
                );
};

class ComparisonRelativeRootPopulationTreeCollection: public BaseComparisonPopulationTreeCollection {

    public:
        ComparisonRelativeRootPopulationTreeCollection() : BaseComparisonPopulationTreeCollection() { }
        ComparisonRelativeRootPopulationTreeCollection(
                const RelativeRootCollectionSettings & settings,
                RandomNumberGenerator & rng,
                bool strict_on_constant_sites = true,
                bool strict_on_missing_sites = true,
                bool strict_on_triallelic_sites = true
                );

    protected:
        void init_trees(
                const std::vector<RelativeRootComparisonSettings> & comparison_settings,
                RandomNumberGenerator & rng,
                bool strict_on_constant_sites = true,
                bool strict_on_missing_sites = true,
                bool strict_on_triallelic_sites = true
                );
};

class ComparisonDirichletPopulationTreeCollection: public BaseComparisonPopulationTreeCollection {

    public:
        ComparisonDirichletPopulationTreeCollection() : BaseComparisonPopulationTreeCollection() { }
        ComparisonDirichletPopulationTreeCollection(
                const DirichletCollectionSettings & settings,
                RandomNumberGenerator & rng,
                bool strict_on_constant_sites = true,
                bool strict_on_missing_sites = true,
                bool strict_on_triallelic_sites = true
                );

    protected:
        void init_trees(
                const std::vector<DirichletComparisonSettings> & comparison_settings,
                RandomNumberGenerator & rng,
                bool strict_on_constant_sites = true,
                bool strict_on_missing_sites = true,
                bool strict_on_triallelic_sites = true
                );
};

#endif
