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

#ifndef ECOEVOLITY_TREE_HPP
#define ECOEVOLITY_TREE_HPP

#include "data.hpp"
#include "node.hpp"
#include "likelihood.hpp"
#include "parameter.hpp"
#include "probability.hpp"
#include "error.hpp"
#include "assert.hpp"

template<class NodeType>
class BaseTree {
    protected:
        std::shared_ptr<NodeType> root_;
        std::vector< std::shared_ptr<PositiveRealParameter> > node_heights_;
        LogProbabilityDensity log_likelihood_ = LogProbabilityDensity(0.0);
        LogProbabilityDensity log_prior_density_ = LogProbabilityDensity(0.0);
        bool ignore_data_ = false;
        unsigned int number_of_likelihood_calculations_ = 0;

        void split_singleton_polytomy(
                RandomNumberGenerator & rng,
                std::shared_ptr<NodeType> polytomy_node,
                std::shared_ptr<PositiveRealParameter> new_height_parameter,
                const bool refresh_node_heights = false) {
            unsigned int n_children = polytomy_node->get_number_of_children();
            std::vector< std::vector<unsigned int> > child_subsets;
            // Need to avoid the partitions where all children are assigned to
            // their own subset or they are all assigned to one subset. In
            // these cases, there would be no split.
            while ((child_subsets.size() == 0) ||
                    (child_subsets.size() == 1) ||
                    (child_subsets.size() == n_children)) {
                child_subsets = rng.random_set_partition_as_subsets(
                        n_children);
            }
            // Need to get the pointers to the children to split (can't
            // work with child indices, because these will change as
            // they are split off from polytomy
            std::vector< std::vector< std::shared_ptr<NodeType> > > child_node_subsets;
            for (auto child_subset : child_subsets) {
                // Any singleton children remain at polytomy node
                if (child_subset.size() < 2) {
                    continue;
                }
                // Any groups of children are split from polytomy node
                child_node_subsets.push_back(polytomy_node->get_children(child_subset));
            }
            for (unsigned int subset_index = 0;
                    subset_index < child_node_subsets.size();
                    ++subset_index) {
                polytomy_node->split_children_from_polytomy(
                        child_node_subsets.at(subset_index),
                        new_height_parameter);
            }
            if (refresh_node_heights) {
                this->update_node_heights();
            }
            else {
                this->node_heights_.push_back(new_height_parameter);
                this->sort_node_heights();
            }
        }

    public:
        BaseTree() { }
        BaseTree(std::shared_ptr<NodeType> root) {
            this->set_root(root);
        }

        void update_node_heights() {
            this->node_heights_.clear();
            std::vector< std::shared_ptr<NodeType> > internal_nodes = this->root_->get_internal_nodes();
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
            this->sort_node_heights();
        }

        void sort_node_heights() {
            std::sort(this->node_heights_.begin(), this->node_heights_.end(), PositiveRealParameter::sort_by_value);
        }

        bool root_has_parent() const {
            if (this->root_->has_parent()) {
                return true;
            }
            return false;
        }
        bool root_has_children() const {
            if (this->root_->get_number_of_children() < 1) {
                return false;
            }
            return true;
        }
        bool root_is_valid() const {
            if (this->root_has_parent()) {
                return false;
            }
            if (! this->root_has_children()) {
                return false;
            }
            return true;
        }
        bool node_heights_are_valid() const {
            return this->root_->node_heights_are_valid();
        }
        bool tree_is_valid() const {
            if (! this->root_is_valid()) {
                return false;
            }
            if (! this->node_heights_are_valid()) {
                return false;
            }
            return true;
        }

        void vet_root() const {
            if (this->root_has_parent()) {
                throw EcoevolityError("Root has a parent");
            }
            if (! this->root_has_children()) {
                throw EcoevolityError("Root has fewer than 2 children");
            }
        }
        void vet_tree() const {
            this->vet_root();
            if (! this->node_heights_are_valid()) {
                throw EcoevolityError("Node ages are not valid");
            }
        }

        unsigned int get_node_height_index(const std::shared_ptr<PositiveRealParameter> height) {
            for (unsigned int i = 0; i < this->node_heights_.size(); ++i) {
                if (this->node_heights_.at(i) == height) {
                    return i;
                }
            }
            throw EcoevolityError("Node height does not exist");
        }

        std::vector< std::shared_ptr<NodeType> > get_mapped_nodes(const unsigned int height_index) {
            return this->root_->get_mapped_nodes(this->node_heights_.at(height_index));
        }

        std::vector< std::shared_ptr<NodeType> > get_mapped_polytomy_nodes(const unsigned int height_index) {
            return this->root_->get_mapped_polytomy_nodes(this->node_heights_.at(height_index));
        }

        std::vector< std::shared_ptr<NodeType> > get_polytomy_nodes() {
            return this->root_->get_polytomy_nodes();
        }

        void merge_node_height_up(const unsigned int height_index,
                const bool refresh_node_heights = false) {
            // Make sure we aren't dealing with the root node
            ECOEVOLITY_ASSERT(height_index < (this->get_number_of_node_heights() - 1));

            std::shared_ptr<PositiveRealParameter> new_height = this->node_heights_.at(height_index + 1);
            std::vector< std::shared_ptr<NodeType> > mapped_nodes = this->get_mapped_nodes(height_index);
            for (unsigned int i = 0; i < mapped_nodes.size(); ++i) {
                // If the parent of the node we are moving up is assigned to the next larger
                // node height, we need to add the child to a polytomy
                if (mapped_nodes.at(i)->get_parent()->get_height_parameter() == new_height) {
                    mapped_nodes.at(i)->collapse();
                }
                else {
                    mapped_nodes.at(i)->set_height_parameter(new_height);
                }
            }
            if (refresh_node_heights) {
                this->update_node_heights();
            }
            else {
                this->node_heights_.erase(this->node_heights_.begin() + height_index);
            }
        }

        void split_node_height_down(
                RandomNumberGenerator & rng,
                const unsigned int height_index,
                const bool refresh_node_heights = false) {
            std::vector< std::shared_ptr<NodeType> > mapped_nodes = this->get_mapped_nodes(height_index);
            if (mapped_nodes.size() < 1) {
                return;
            }
            double max_height = this->node_heights_.at(height_index)->get_value();
            double min_height = 0.0;
            if (height_index > 0) {
                min_height = this->node_heights_.at(height_index - 1)->get_value();
            }
            double new_height = rng.uniform_real(min_height, max_height);
            std::shared_ptr<PositiveRealParameter> new_height_parameter = std::make_shared<PositiveRealParameter>(new_height);
            if (mapped_nodes.size() == 1) {
                // If we only have a single polytomy, we need to handle the
                // splitting differently
                ECOEVOLITY_ASSERT(mapped_nodes.at(0)->is_polytomy());
                this->split_singleton_polytomy(rng, mapped_nodes.at(0),
                        new_height_parameter,
                        refresh_node_heights);
                return;
            }

            std::vector< std::vector<unsigned int> > subsets = rng.random_subsets(
                    mapped_nodes.size(), 2);

            unsigned int move_subset_index = rng.uniform_positive_int(0, 1);

            for (auto node_index : subsets.at(move_subset_index)) {
                // If node is a polytomy, we need to deal with possibility of
                // splitting it up
                if (mapped_nodes.at(node_index)->is_polytomy()) {
                    unsigned int n_children = mapped_nodes.at(node_index)->get_number_of_children();
                    std::vector< std::vector<unsigned int> > child_subsets;
                    // Need to avoid the partition where all children are
                    // assigned to their own subset. This would result in the
                    // no nodes being moved to the new height (i.e., the
                    // polytomy remains as is), which complicates Hasting's
                    // ratio, because then there are multiple ways polytomy
                    // nodes are not included in the split (it can either not
                    // be selected in the move pool, or selected, but all the
                    // children are in singleton subsets and thus not split
                    // off).
                    while ((child_subsets.size() == 0) || (child_subsets.size() == n_children)) {
                        child_subsets = rng.random_set_partition_as_subsets(
                                n_children);
                    }
                    // Need to get the pointers to the children to split (can't
                    // work with child indices, because these will change as
                    // they are split off from polytomy
                    std::vector< std::vector< std::shared_ptr<NodeType> > > child_node_subsets;
                    for (auto child_subset : child_subsets) {
                        // Any singleton children remain at polytomy node
                        if (child_subset.size() < 2) {
                            continue;
                        }
                        // Any groups of children are split from polytomy node
                        child_node_subsets.push_back(mapped_nodes.at(node_index)->get_children(child_subset));
                    }
                    for (unsigned int subset_index = 0;
                            subset_index < child_node_subsets.size();
                            ++subset_index) {
                        mapped_nodes.at(node_index)->split_children_from_polytomy(
                                child_node_subsets.at(subset_index),
                                new_height_parameter);
                    }
                }
                // Node is not a polytomy, so we simply assign it to the new height
                else {
                    mapped_nodes.at(node_index)->set_height_parameter(new_height_parameter);
                }
            }
            if (refresh_node_heights) {
                this->update_node_heights();
            }
            else {
                this->node_heights_.push_back(new_height_parameter);
                this->sort_node_heights();
            }
        }

        std::vector<unsigned int> get_indices_of_splittable_heights() const {
            std::vector<unsigned int> splittable_heights;
            for (unsigned int i = 0; i < this->node_heights_.size(); ++i) {
                if (this->height_is_splittable(i)) {
                    splittable_heights.push_back(i);
                }
            }
            return splittable_heights;
        }

        bool height_is_splittable(const unsigned int height_index) const {
            unsigned int mapped_node_count = this->root_->get_mapped_node_count(
                    this->node_heights_.at(height_index));
            if (mapped_node_count > 1) {
                return true;
            }
            unsigned int mapped_polytomy_node_count = this->root_->get_mapped_polytomy_node_count(
                    this->node_heights_.at(height_index));
            if (mapped_polytomy_node_count > 0) {
                return true;
            }
            return false;
        }

        void slide_bump_height(const unsigned int height_index,
                const double new_height) {
            ECOEVOLITY_ASSERT(new_height >= 0.0);
            std::vector<unsigned int> intervening_indices = this->get_intervening_height_indices(
                    height_index,
                    new_height);
            if (intervening_indices.size() < 1) {
                // No intervening nodes to bump
                this->node_heights_.at(height_index)->set_value(new_height);
                return;
            }
            if (height_index < intervening_indices.at(0)) {
                // Older nodes to bump up
                for (auto next_height_idx : intervening_indices) {
                    this->node_heights_.at(next_height_idx - 1)->set_value(
                            this->node_heights_.at(next_height_idx)->get_value()
                            );
                }
                this->node_heights_.at(intervening_indices.back())->set_value(new_height);
                return;
            }
            // Younger nodes to bump down
            for (auto next_height_idx : intervening_indices) {
                this->node_heights_.at(next_height_idx + 1)->set_value(
                        this->node_heights_.at(next_height_idx)->get_value()
                        );
            }
            this->node_heights_.at(intervening_indices.back())->set_value(new_height);
            return;
        }

        void slide_bump_swap_height(
                RandomNumberGenerator & rng,
                const unsigned int height_index,
                const double new_height) {
            ECOEVOLITY_ASSERT(new_height >= 0.0);
            std::vector<unsigned int> intervening_indices = this->get_intervening_height_indices(
                    height_index,
                    new_height);
            // TODO
        }

        std::vector<unsigned int> get_intervening_height_indices(
                const unsigned int height_index,
                const double value) {
            std::vector<unsigned int> indices;
            // value is larger than this height
            if (this->get_height(height_index) < value) {
                if (height_index == (this->get_number_of_node_heights() - 1)) {
                    return indices;
                }
                for (unsigned int i = (height_index + 1);
                        i < this->get_number_of_node_heights();
                        ++i) {
                    if (this->get_height(i) > value) {
                        break;
                    }
                    indices.push_back(i);
                }
                return indices;
            }
            // value is less than this height
            if (height_index == 0) {
                return indices;
            }
            // decrementing an unsigned int to zero is a bit tricky; this stops
            // before we hit zero in the evaluation, but decrements before
            // jumping into loop body
            for (unsigned int i = height_index;
                    i-- > 0;
                ) {
                if (this->get_height(i) < value) {
                    break;
                }
                indices.push_back(i);
            }
            return indices;
        }

        unsigned int get_nearest_height_index(const double value) {
            if (this->get_number_of_node_heights() < 2) {
                return 0;
            }
            if (this->get_height(0) > value) {
                return 0;
            }
            if (this->get_height(this->get_number_of_node_heights() - 1) < value) {
                return this->get_number_of_node_heights() - 1;
            }
            double diff = fabs(this->get_height(0) - value);
            double current_diff = diff;
            for (unsigned int i = 1; i < this->node_heights_.size(); ++i) {
                current_diff = fabs(this->get_height(i) - value);
                if (current_diff > diff) {
                    return i - 1;
                }
                diff = current_diff;
            }
            return this->get_number_of_node_heights() - 1;
        }

        void set_root(std::shared_ptr<NodeType> root) {
            this->root_ = root;
            this->vet_tree();
            this->update_node_heights();
        }

        const NodeType& get_root() const {return *this->root_;}
        NodeType& get_mutable_root() const {return *this->root_;}

        virtual void set_root_node_height_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->root_->set_node_height_prior(prior);
        }
        virtual std::shared_ptr<ContinuousProbabilityDistribution> get_root_node_height_prior() const {
            return this->root_->get_node_height_prior();
        }

        double get_height(const unsigned int height_index) const {
            return this->node_heights_.at(height_index)->get_value();
        }

        std::shared_ptr<PositiveRealParameter> get_height_parameter(const unsigned int height_index) const {
            return this->node_heights_.at(height_index);
        }

        void set_root_height(double height) {
            this->root_->set_height(height);
        }
        double get_root_height() const {
            return this->root_->get_height();
        }

        unsigned int get_degree_of_root() const {
            return this->root_->degree();
        }

        unsigned int get_leaf_node_count() const {
            return this->root_->get_leaf_node_count();
        }
        unsigned int get_node_count() const {
            return this->root_->get_node_count();
        }

        unsigned int get_number_of_node_heights() const {
            return this->node_heights_.size();
        }

        const std::vector< std::shared_ptr<PositiveRealParameter> >& get_node_height_pointers() const {
            return this->node_heights_;
        }

        std::vector<double> get_node_heights() const {
            std::vector<double> heights (this->node_heights_.size());
            for (unsigned int i = 0; i < this->node_heights_.size(); ++i) {
                heights.at(i) = this->node_heights_.at(i)->get_value();
            }
            return heights;
        }

        void ignore_data() {
            this->ignore_data_ = true;
        }
        void use_data() {
            this->ignore_data_ = false;
        }
        bool ignoring_data() const {
            return this->ignore_data_;
        }

        unsigned int get_number_of_likelihood_calculations() {
            return this->number_of_likelihood_calculations_;
        }

        virtual void make_clean() {
            this->root_->make_all_clean();
        }

        virtual void compute_log_likelihood_and_prior(unsigned int nthreads = 1) {
            this->compute_log_likelihood(nthreads);
            this->compute_log_prior_density();
            this->make_clean();
            return;
        }

        virtual double compute_log_likelihood(unsigned int nthreads = 1) {
            ++this->number_of_likelihood_calculations_;
            this->log_likelihood_.set_value(0.0);
            return 0.0;
        }

        void set_log_likelihood_value(double value) {
            this->log_likelihood_.set_value(value);
        }
        double get_log_likelihood_value() const {
            return this->log_likelihood_.get_value();
        }
        double get_stored_log_likelihood_value() const {
            return this->log_likelihood_.get_stored_value();
        }

        virtual double compute_log_prior_density() {
            double d = 0.0;
            d += this->compute_log_prior_density_of_node_heights();
            d += this->compute_relative_log_prior_density_of_toplogy();
            this->log_prior_density_.set_value(d);
            return d;
        }

        virtual double compute_log_prior_density_of_node_heights() const {
            double root_height = this->root_->get_height();
            double d = 0.0;
            d += this->root_->get_height_relative_prior_ln_pdf();
            // prior prob density of non-root internal nodes = (1 / root_height)^n,
            // where n = the number of non-root internal nodes, so on
            // log scale = n * (log(1)-log(root_height)) = n * -log(root_height)
            double internal_node_height_prior_density = -std::log(root_height);
            d += internal_node_height_prior_density * (this->get_number_of_node_heights() - 1);
            return d;
        }
        virtual double compute_relative_log_prior_density_of_toplogy() const {
            return 0.0;
        }

        double get_log_prior_density_value() const {
            return this->log_prior_density_.get_value();
        }
        double get_stored_log_prior_density_value() const {
            return this->log_prior_density_.get_stored_value();
        }

        void store_state() {
            this->store_likelihood();
            this->store_prior_density();
            this->store_parameters();
        }
        void store_likelihood() {
            this->log_likelihood_.store();
        }
        void store_prior_density() {
            this->log_prior_density_.store();
        }
        virtual void store_parameters() {
            this->store_all_heights();
            this->store_all_height_pointers();
        }
        virtual void store_all_heights() {
            this->root_->store_all_heights();
        }
        virtual void store_all_height_pointers() {
            this->root_->store_all_height_pointers();
        }
        void restore_state() {
            this->restore_likelihood();
            this->restore_prior_density();
            this->restore_parameters();
        }
        void restore_likelihood() {
            this->log_likelihood_.restore();
        }
        void restore_prior_density() {
            this->log_prior_density_.restore();
        }
        virtual void restore_parameters() {
            this->restore_all_height_pointers();
            this->restore_all_heights();
        }
        virtual void restore_all_heights() {
            this->root_->restore_all_heights();
        }
        virtual void restore_all_height_pointers() {
            this->root_->restore_all_height_pointers();
            this->update_node_heights();
        }

        virtual void write_state_log_header(std::ostream& out,
                bool include_event_index,
                const std::string& delimiter = "\t") const {
            throw EcoevolityError("write_state_log_header called from base BaseTree class");
        }
        virtual void log_state(std::ostream& out,
                unsigned int event_index,
                const std::string& delimiter = "\t") const {
            throw EcoevolityError("log_state called from base BaseTree class");
        }
        virtual void log_state(std::ostream& out,
                const std::string& delimiter = "\t") const {
            throw EcoevolityError("log_state called from base BaseTree class");
        }

        virtual void draw_from_prior(RandomNumberGenerator& rng) {
            throw EcoevolityError("draw_from_prior called from base BaseTree class");
        }
};


class BasePopulationTree : public BaseTree<PopulationNode> {
    protected:
        BiallelicData data_;
        std::shared_ptr<ContinuousProbabilityDistribution> population_size_prior_ = std::make_shared<GammaDistribution>(1.0, 0.001);
        double ploidy_ = 2.0;
        std::shared_ptr<PositiveRealParameter> freq_1_ = std::make_shared<PositiveRealParameter>(
                std::make_shared<BetaDistribution>(1.0, 1.0),
                0.5);
        std::shared_ptr<PositiveRealParameter> mutation_rate_ = std::make_shared<PositiveRealParameter>(
                1.0,
                true);
        LogProbabilityDensity log_likelihood_correction_ = LogProbabilityDensity(0.0);
        bool likelihood_correction_was_calculated_ = false;
        bool constant_sites_removed_ = true;
        // int provided_number_of_constant_red_sites_ = -1;
        // int provided_number_of_constant_green_sites_ = -1;
        // bool use_removed_constant_site_counts_ = false;
        bool population_sizes_are_constrained_ = false;
        bool state_frequencies_are_constrained_ = false;
        bool is_dirty_ = true;

        // Vectors for storing unique allele counts and associated weights.
        // These are used for calculating the likelihood correction term for
        // constant site patterns. It is a bit weird to store data here, but
        // it's cheaper than calling 'this->data_.get_unique_allele_counts()'
        // every time the likelihood needs to be calculated.
        // These vectors are populated in 'init' method.
        std::vector< std::vector<unsigned int> > unique_allele_counts_;
        std::vector<unsigned int> unique_allele_count_weights_;

        // methods
        void init_tree();
        // bool constant_site_counts_were_provided();
        void calculate_likelihood_correction();

        double calculate_log_binomial(
                unsigned int red_allele_count,
                unsigned int allele_count) const;

        void set_population_sizes(
                std::shared_ptr<PopulationNode> node,
                const std::vector<double> & sizes);
        
        void get_population_sizes(
                std::shared_ptr<PopulationNode> node,
                std::vector<double> & sizes) const;

        void update_unique_allele_counts();

    public:
        BasePopulationTree() { }
        BasePopulationTree(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                );

        BasePopulationTree(
                std::shared_ptr<PopulationNode> root,
                unsigned int number_of_loci = 10000,
                unsigned int length_of_loci = 1,
                bool validate_data = false);

        void init(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                );

        bool has_seq_loci_info() const {
            return this->data_.has_seq_loci_info();
        }

        void fold_patterns();

        bool constant_sites_removed() const {
            return this->constant_sites_removed_;
        }

        // int get_provided_number_of_constant_red_sites() const {
        //     return this->provided_number_of_constant_red_sites_;
        // }
        // int get_provided_number_of_constant_green_sites() const {
        //     return this->provided_number_of_constant_green_sites_;
        // }

        bool initialized() const {return (bool)this->root_;}

        const std::vector<std::string>& get_population_labels() const {
            return this->data_.get_population_labels();
        }

        const BiallelicData& get_data() const {
            return this->data_;
        }

        void set_data(const BiallelicData & data, bool constant_sites_removed);

        void set_ploidy(double ploidy) {
            this->ploidy_ = ploidy;
        }
        double get_ploidy() const {
            return this->ploidy_;
        }

        virtual double get_node_theta(const PopulationNode& node) const {
            return (2 * this->get_ploidy() *
                    node.get_population_size() *
                    this->get_mutation_rate());
        }
        double get_node_length_in_subs_per_site(const PopulationNode& node) const {
            return (node.get_length() * this->get_mutation_rate());
        }
        double get_node_height_in_subs_per_site(const PopulationNode& node) const {
            return (node.get_height() * this->get_mutation_rate());
        }

        void set_freq_1(double p);
        double get_freq_1() const;
        double get_freq_0() const;
        void store_freq_1();
        void restore_freq_1();
        double get_u() const;
        double get_v() const;

        void set_mutation_rate(double m);
        double get_mutation_rate() const;
        void store_mutation_rate();
        void restore_mutation_rate();

        bool is_dirty() const;
        void make_dirty();
        void make_clean();

        // void provide_number_of_constant_sites(
        //         unsigned int number_all_red,
        //         unsigned int number_all_green);

        std::shared_ptr<PositiveRealParameter> get_freq_1_parameter() const;

        void set_mutation_rate_parameter(std::shared_ptr<PositiveRealParameter> h);
        std::shared_ptr<PositiveRealParameter> get_mutation_rate_parameter() const;

        virtual void set_root_population_size(double size);
        virtual void set_all_population_sizes(double size);
        virtual unsigned int scale_all_population_sizes(double scale);
        virtual unsigned int scale_root_population_size(double scale);
        virtual double get_root_population_size() const;
        virtual std::shared_ptr<PositiveRealParameter> get_root_population_size_parameter() const;

        std::vector<double> get_population_sizes() const;

        double get_mean_population_size() const;
        virtual void set_mean_population_size(double size);

        double get_leaf_mean_population_size() const;

        std::vector<double> get_population_sizes_as_proportions() const;
        std::vector<double> get_population_sizes_as_multipliers() const;
        virtual void set_population_sizes_as_proportions(const std::vector<double> & proportions);
        virtual void set_population_sizes(const std::vector<double> & sizes);

        double get_likelihood_correction(bool force = false);

        void compute_log_likelihood_and_prior(unsigned int nthreads = 1) {
            if (this->is_dirty()) {
                this->compute_log_likelihood(nthreads);
                ++this->number_of_likelihood_calculations_;
                this->compute_log_prior_density();
                this->make_clean();
            }
            return;
        }

        double compute_log_likelihood(unsigned int nthreads = 1);

        double compute_log_prior_density();
        double compute_log_prior_density_of_state_frequencies() const;
        double compute_log_prior_density_of_mutation_rate() const;
        virtual double compute_log_prior_density_of_population_sizes() const;

        void store_parameters();
        virtual void store_all_population_sizes();
        void restore_parameters();
        virtual void restore_all_population_sizes();

        virtual void set_population_size_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior);
        virtual std::shared_ptr<ContinuousProbabilityDistribution> get_population_size_prior() const {
            return this->population_size_prior_;
        }

        void set_freq_1_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior);
        std::shared_ptr<ContinuousProbabilityDistribution> get_freq_1_prior() const {
            return this->freq_1_->prior;
        }

        void set_mutation_rate_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior);
        std::shared_ptr<ContinuousProbabilityDistribution> get_mutation_rate_prior() const {
            return this->mutation_rate_->prior;
        }

        virtual void fix_population_sizes() {
            this->root_->fix_all_population_sizes();
        }
        virtual void estimate_population_sizes() {
            this->root_->estimate_all_population_sizes();
        }
        virtual bool population_sizes_are_fixed() const {
            return this->root_->all_population_sizes_are_fixed();
        }

        virtual bool root_population_size_is_fixed() const {
            return this->root_->population_size_is_fixed();
        }

        virtual void constrain_population_sizes() {
            this->population_sizes_are_constrained_ = true;
            this->root_->set_all_population_size_parameters();
        }
        virtual bool population_sizes_are_constrained() const {
            return this->population_sizes_are_constrained_;
        }

        void fix_state_frequencies() {
            this->freq_1_->fix();
        }
        void estimate_state_frequencies() {
            if (this->state_frequencies_are_constrained_) {
                throw EcoevolityError("Cannot estimate constrained state frequencies");
            }
            this->freq_1_->estimate();
        }
        bool state_frequencies_are_fixed() const {
            return this->freq_1_->is_fixed();
        }

        void fix_mutation_rate() {
            this->mutation_rate_->fix();
        }
        void estimate_mutation_rate() {
            this->mutation_rate_->estimate();
        }
        bool mutation_rate_is_fixed() const {
            return this->mutation_rate_->is_fixed();
        }

        void constrain_state_frequencies() {
            this->state_frequencies_are_constrained_ = true;
            this->freq_1_->set_value(0.5);
            this->freq_1_->fix();
            this->make_dirty();
        }
        bool state_frequencies_are_constrained() const {
            return this->state_frequencies_are_constrained_;
        }

        void simulate_gene_tree(
                const std::shared_ptr<PopulationNode> node,
                std::unordered_map<unsigned int, std::vector< std::shared_ptr<GeneTreeSimNode> > > & branch_lineages,
                const unsigned int pattern_index,
                RandomNumberGenerator & rng,
                const bool use_max_allele_counts = false) const;

        std::shared_ptr<GeneTreeSimNode> simulate_gene_tree(
                const unsigned int pattern_index,
                RandomNumberGenerator& rng,
                const bool use_max_allele_counts = false) const;

        static double coalesce_in_branch(
                std::vector< std::shared_ptr<GeneTreeSimNode> >& lineages,
                double population_size,
                RandomNumberGenerator& rng,
                double bottom_of_branch_height = 0.0,
                double top_of_branch_height = std::numeric_limits<double>::infinity(),
                unsigned int branch_index = 0
                );

        bool sample_pattern(
                RandomNumberGenerator& rng,
                const float singleton_sample_probability,
                const std::vector<unsigned int>& red_allele_counts,
                const std::vector<unsigned int>& allele_counts) const;

        BiallelicData simulate_biallelic_data_set(
                RandomNumberGenerator& rng,
                float singleton_sample_probability = 1.0,
                bool validate = true) const;

        BiallelicData simulate_linked_biallelic_data_set(
                RandomNumberGenerator& rng,
                float singleton_sample_probability,
                bool max_one_variable_site_per_locus = false,
                bool validate = true) const;

        std::pair<BiallelicData, unsigned int>
        simulate_complete_biallelic_data_set(
                RandomNumberGenerator& rng,
                unsigned int locus_size = 1,
                float singleton_sample_probability = 1.0,
                bool validate = true) const;

        std::pair<BiallelicData, unsigned int>
        simulate_data_set_max_one_variable_site_per_locus(
                RandomNumberGenerator& rng,
                unsigned int locus_size,
                float singleton_sample_probability,
                bool validate = true) const;

        std::pair<
                std::pair<std::vector<unsigned int>, std::vector<unsigned int> >,
                std::shared_ptr<GeneTreeSimNode> >
        simulate_biallelic_site(
                const unsigned int pattern_idx,
                RandomNumberGenerator& rng,
                const bool use_max_allele_counts = false) const;

        std::pair<std::vector<unsigned int>, std::vector<unsigned int> >
        simulate_biallelic_site(
                std::shared_ptr<GeneTreeSimNode> gene_tree,
                RandomNumberGenerator& rng) const;

        std::pair<std::vector<unsigned int>, std::vector<unsigned int> >
        simulate_biallelic_site_sans_missing(
                std::shared_ptr<GeneTreeSimNode> gene_tree,
                const std::vector<unsigned int> & site_allele_counts,
                RandomNumberGenerator& rng) const;

        void write_data_summary(
                std::ostream& out,
                unsigned int indent_level = 0) const {
            this->data_.write_summary(out, indent_level);
        }

};



// TODO: PopulationTree is misnomer for this class and should be changed.
// This name made sense when this class used to be at the very base of the tree
// class hierarchy. Now it is an intermediate class in the hierarchy, an
// intermediate leading to the comparison tree classes (only one or two tips).
class PopulationTree : public BasePopulationTree {

    protected:
        std::shared_ptr<ContinuousProbabilityDistribution> root_node_height_prior_ = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<DirichletDistribution> population_size_multiplier_prior_;
        bool population_size_multipliers_are_fixed_ = false;
        bool mean_population_size_is_fixed_ = false;

        void update_root_population_size() { return; }
        void update_relative_root_population_size() { return; }

    public:

        PopulationTree() { }
        PopulationTree(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                ) : BasePopulationTree(
                    path,
                    population_name_delimiter,
                    population_name_is_prefix,
                    genotypes_are_diploid,
                    markers_are_dominant,
                    constant_sites_removed,
                    validate,
                    strict_on_constant_sites,
                    strict_on_missing_sites,
                    strict_on_triallelic_sites,
                    ploidy,
                    store_seq_loci_info) { }

        PopulationTree(
                std::shared_ptr<PopulationNode> root,
                unsigned int number_of_loci = 10000,
                unsigned int length_of_loci = 1,
                bool validate_data = false
                ) : BasePopulationTree(
                    root,
                    number_of_loci,
                    length_of_loci,
                    validate_data) { }

        virtual bool using_population_size_multipliers() const {
            return false;
        }

        virtual bool using_relative_root_population_size() const {
            return false;
        }

        // These are overloaded by RelativeRootPopulationTree
        bool relative_root_population_size_is_fixed() const { return false; }
        void fix_relative_root_population_size() { return; }
        void estimate_relative_root_population_size() { return; }


        void set_root_node_height_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->root_node_height_prior_ = prior;
            this->root_->set_all_node_height_priors(prior);
        }
        std::shared_ptr<ContinuousProbabilityDistribution> get_root_node_height_prior() const {
            return this->root_node_height_prior_;
        }
        void set_node_height_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->set_root_node_height_prior(prior);
        }
        std::shared_ptr<ContinuousProbabilityDistribution> get_node_height_prior() const {
            return this->get_root_node_height_prior();
        }

        void set_root_height_parameter(std::shared_ptr<PositiveRealParameter> h);
        std::shared_ptr<PositiveRealParameter> get_root_height_parameter() const;


        void store_root_height();
        void restore_root_height();

        virtual void set_relative_root_population_size_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) { return; }

        virtual double get_relative_root_population_size() const {
            return this->get_root_population_size() / this->get_leaf_mean_population_size();
        }

        virtual void fix_population_sizes() {
            this->fix_mean_population_size();
            this->fix_population_size_multipliers();
            this->root_->fix_all_population_sizes();
        }
        virtual void estimate_population_sizes() {
            this->estimate_mean_population_size();
            this->estimate_population_size_multipliers();
            this->root_->estimate_all_population_sizes();
        }

        void fix_population_size_multipliers() {
            this->population_size_multipliers_are_fixed_ = true;
            this->make_dirty();
        }
        void estimate_population_size_multipliers() {
            this->population_size_multipliers_are_fixed_ = false;
            this->make_dirty();
        }
        bool population_size_multipliers_are_fixed() const {
            return this->population_size_multipliers_are_fixed_;
        }

        void fix_mean_population_size() {
            this->mean_population_size_is_fixed_ = true;
            this->make_dirty();
        }
        void estimate_mean_population_size() {
            this->mean_population_size_is_fixed_ = false;
            this->make_dirty();
        }
        bool mean_population_size_is_fixed() const {
            return this->mean_population_size_is_fixed_;
        }

        void set_population_size_multiplier_prior(std::shared_ptr<DirichletDistribution> prior);
        std::shared_ptr<DirichletDistribution> get_population_size_multiplier_prior() const{
            return this->population_size_multiplier_prior_;
        }

        void set_population_sizes(
                const std::vector<double> & sizes);
        
        void set_population_sizes_as_proportions(const std::vector<double> & proportions);
        
        void set_mean_population_size(double size);

        void store_parameters();
        void restore_parameters();

        // Override this method to old node height prior behavior
        double compute_log_prior_density_of_node_heights() const;

        // TODO: This PopulationTree hierarchy of classes is messy. The problem
        // is that each derived class has its own subset of methods in addition
        // to the base class methods. Thus, I can't simply use "virtual ... =
        // 0;" here, because some of the derived classes would be left with
        // invalid methods. Declaring methods here that throw errors if not
        // overloaded works for now.
        // Methods to be overloaded
        virtual void set_child_population_size(unsigned int child_index, double size) {
            throw EcoevolityError("set_child_population_size called from PopulationTree");
        }
        virtual double get_child_population_size(unsigned int child_index) const {
            throw EcoevolityError("get_child_population_size called from PopulationTree");
        }
        virtual std::shared_ptr<PositiveRealParameter> get_child_population_size_parameter(
                unsigned int child_index) const {
            throw EcoevolityError("get_child_population_size_parameter called from PopulationTree");
        }
};


class RelativeRootPopulationTree : public PopulationTree {

    protected:
        std::shared_ptr<PositiveRealParameter> relative_root_population_size_ = std::make_shared<PositiveRealParameter>(
                std::make_shared<GammaDistribution>(10.0, 0.1),
                1.0);
        bool leaf_population_sizes_are_fixed_ = false;

        void update_root_population_size();
        void update_relative_root_population_size();

    public:
        RelativeRootPopulationTree() { }
        RelativeRootPopulationTree(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                ) : PopulationTree(
                    path,
                    population_name_delimiter,
                    population_name_is_prefix,
                    genotypes_are_diploid,
                    markers_are_dominant,
                    constant_sites_removed,
                    validate,
                    strict_on_constant_sites,
                    strict_on_missing_sites,
                    strict_on_triallelic_sites,
                    ploidy,
                    store_seq_loci_info) { }

        bool using_relative_root_population_size() const {
            return true;
        }

        bool relative_root_population_size_is_fixed() const {
            return this->relative_root_population_size_->is_fixed();
        }
        void fix_relative_root_population_size() {
            this->relative_root_population_size_->fix();
        }
        void estimate_relative_root_population_size() {
            this->relative_root_population_size_->estimate();
            this->root_->estimate_population_size();
        }

        void set_root_population_size(double size);
        double get_root_population_size() const;
        void set_all_population_sizes(double size);
        unsigned int scale_all_population_sizes(double scale);
        unsigned int scale_root_population_size(double scale);
        void set_mean_population_size(double size);

        double get_relative_root_population_size() const {
            return this->relative_root_population_size_->get_value();
        }

        void set_population_sizes_as_proportions(
                const std::vector<double> & proportions);
        void set_population_sizes(const std::vector<double> & sizes);

        void set_relative_root_population_size_prior(
                std::shared_ptr<ContinuousProbabilityDistribution> prior);

        void fix_population_sizes() {
            PopulationTree::fix_population_sizes();
            if (! this->relative_root_population_size_is_fixed()) {
                this->root_->estimate_population_size();
            }
            this->leaf_population_sizes_are_fixed_ = true;
        }
        void estimate_population_sizes() {
            PopulationTree::estimate_population_sizes();
            this->leaf_population_sizes_are_fixed_ = false;
        }
        virtual void constrain_population_sizes() {
            PopulationTree::constrain_population_sizes();
            this->relative_root_population_size_->set_value(1.0);
            this->relative_root_population_size_->fix();
        }

        void store_all_population_sizes();
        void restore_all_population_sizes();

        double compute_log_prior_density_of_population_sizes() const;
};


class DirichletPopulationTree: public PopulationTree {

    public:
        DirichletPopulationTree() { }
        DirichletPopulationTree(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                );

        bool using_population_size_multipliers() const {
            return true;
        }

        void constrain_population_sizes() {
            throw EcoevolityError("This method is not supported; multipliers should be set and fixed");
        }
        bool population_sizes_are_constrained() const {
            throw EcoevolityError("This method is not supported; check if multipliers are fixed");
        }

        double compute_log_prior_density_of_population_sizes() const;
};


class ComparisonPopulationTree: public PopulationTree {

    public:
        ComparisonPopulationTree() { }
        ComparisonPopulationTree(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                );
        ComparisonPopulationTree(
                const ComparisonSettings& settings,
                RandomNumberGenerator& rng,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                bool store_seq_loci_info = false
                );
        void comparison_init(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                );

        void set_child_population_size(unsigned int child_index, double size);
        double get_child_population_size(unsigned int child_index) const;
        std::shared_ptr<PositiveRealParameter> get_child_population_size_parameter(
                unsigned int child_index) const;

        void store_all_heights() {
            this->store_root_height();
        }
        void restore_all_heights() {
            this->restore_root_height();
        }
        void store_parameters();
        void restore_parameters();

        double compute_log_prior_density();

        void write_state_log_header(std::ostream& out,
                bool include_event_index,
                const std::string& delimiter = "\t") const;
        void log_state(std::ostream& out,
                unsigned int event_index,
                const std::string& delimiter = "\t") const;
        void log_state(std::ostream& out,
                const std::string& delimiter = "\t") const;

        void draw_from_prior(RandomNumberGenerator& rng);

};

class ComparisonRelativeRootPopulationTree: public RelativeRootPopulationTree {

    public:
        ComparisonRelativeRootPopulationTree() { }
        ComparisonRelativeRootPopulationTree(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                );
        ComparisonRelativeRootPopulationTree(
                const RelativeRootComparisonSettings& settings,
                RandomNumberGenerator& rng,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                bool store_seq_loci_info = false
                );
        void comparison_init(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                );

        void set_child_population_size(unsigned int child_index, double size);
        double get_child_population_size(unsigned int child_index) const;
        std::shared_ptr<PositiveRealParameter> get_child_population_size_parameter(
                unsigned int child_index) const;

        void store_all_heights() {
            this->store_root_height();
        }
        void restore_all_heights() {
            this->restore_root_height();
        }
        void store_parameters();
        void restore_parameters();

        double compute_log_prior_density();

        void write_state_log_header(std::ostream& out,
                bool include_event_index,
                const std::string& delimiter = "\t") const;
        void log_state(std::ostream& out,
                unsigned int event_index,
                const std::string& delimiter = "\t") const;
        void log_state(std::ostream& out,
                const std::string& delimiter = "\t") const;

        void draw_from_prior(RandomNumberGenerator& rng);

};

class ComparisonDirichletPopulationTree: public ComparisonPopulationTree {

    public:
        ComparisonDirichletPopulationTree() { }
        ComparisonDirichletPopulationTree(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                );
        ComparisonDirichletPopulationTree(
                const DirichletComparisonSettings& settings,
                RandomNumberGenerator& rng,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                bool store_seq_loci_info = false
                );
        void comparison_init(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                );

        bool using_population_size_multipliers() const {
            return true;
        }

        void constrain_population_sizes() {
            throw EcoevolityError("This method is not supported; multipliers should be set and fixed");
        }
        bool population_sizes_are_constrained() const {
            throw EcoevolityError("This method is not supported; check if multipliers are fixed");
        }

        double compute_log_prior_density_of_population_sizes() const;

        void write_state_log_header(std::ostream& out,
                bool include_event_index,
                const std::string& delimiter = "\t") const;
        void log_state(std::ostream& out,
                unsigned int event_index,
                const std::string& delimiter = "\t") const;
        void log_state(std::ostream& out,
                const std::string& delimiter = "\t") const;

        void draw_from_prior(RandomNumberGenerator& rng);
};

#endif
