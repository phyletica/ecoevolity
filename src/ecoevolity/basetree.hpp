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

#ifndef ECOEVOLITY_BASETREE_HPP
#define ECOEVOLITY_BASETREE_HPP

#include "node.hpp"
#include "parameter.hpp"
#include "probability.hpp"
#include "error.hpp"
#include "assert.hpp"

template<class NodeType>
class BaseTree {
    protected:
        std::shared_ptr<NodeType> root_;
        std::shared_ptr<NodeType> stored_root_;
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

        std::vector< std::shared_ptr<NodeType> > get_collision_parents(
                const unsigned int older_height_index,
                const unsigned int younger_height_index) {
            ECOEVOLITY_ASSERT(older_height_index == younger_height_index + 1);
            std::vector< std::shared_ptr<NodeType> > collision_parents;
            std::vector< std::shared_ptr<NodeType> > older_nodes = this->get_mapped_nodes(older_height_index);
            std::vector< std::shared_ptr<NodeType> > younger_nodes = this->get_mapped_nodes(younger_height_index);
            for (auto older_nd : older_nodes) {
                for (auto younger_nd : younger_nodes) {
                    if (older_nd->is_child(younger_nd)) {
                        collision_parents.push_back(older_nd);
                        break;
                    }
                }
            }
            return collision_parents;
        }

        void collision_node_swap(
                RandomNumberGenerator & rng,
                const unsigned int older_height_index,
                const unsigned int younger_height_index) {
            std::vector< std::shared_ptr<NodeType> > collision_parents = this->get_collision_parents(
                    older_height_index,
                    younger_height_index);
            for (auto parent_nd : collision_parents) {
                std::vector< std::shared_ptr<NodeType> > swap_node_pool;
                std::vector< std::shared_ptr<NodeType> > child_colliding_nodes;
                std::vector< std::shared_ptr<NodeType> > child_non_colliding_nodes;
                for (auto child_nd : parent_nd->get_all_children()) {
                    // Each colliding child contributes one random grandchild to
                    // the swap pool
                    if ((! child_nd->is_leaf()) &&
                            (this->get_node_height_index(child_nd->get_height_parameter()) ==
                             younger_height_index)) {
                        ECOEVOLITY_ASSERT(child_nd->has_children());
                        child_colliding_nodes.push_back(child_nd);
                        unsigned int random_child_index = rng.uniform_positive_int(child_nd->get_number_of_children() - 1);
                        swap_node_pool.push_back(child_nd->get_child(random_child_index));
                    }
                    else {
                        child_non_colliding_nodes.push_back(child_nd);
                    }
                }
                // If parent node in collision has non-colliding descendants,
                // it contributes one of these descendants (randomly) to the
                // swap pool
                if (child_non_colliding_nodes.size() > 0) {
                    unsigned int random_child_index = rng.uniform_positive_int(child_non_colliding_nodes.size() - 1);
                    swap_node_pool.push_back(child_non_colliding_nodes.at(random_child_index));
                }

                // Randomize the swap pool
                std::shuffle(std::begin(swap_node_pool), std::end(swap_node_pool), rng.engine_);

                std::shared_ptr<NodeType> random_swap_node;
                // If the parent contributed a non-colliding child to the swap
                // pool, randomly assign it a new child from the swap pool
                if (child_non_colliding_nodes.size() > 0) {
                    std::shared_ptr<NodeType> random_swap_node = swap_node_pool.back();
                    swap_node_pool.pop_back();
                    random_swap_node->remove_parent();
                    random_swap_node->add_parent(parent_nd);
                }
                for (auto child_colliding_nd : child_colliding_nodes) {
                    random_swap_node = swap_node_pool.back();
                    swap_node_pool.pop_back();
                    random_swap_node->remove_parent();
                    random_swap_node->add_parent(child_colliding_nd);
                }
                ECOEVOLITY_ASSERT(swap_node_pool.size() == 0);
            }
        }

        void slide_bump_height(
                RandomNumberGenerator & rng,
                const unsigned int height_index,
                const double new_height,
                bool collisions_swap_nodes = false) {
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
                    if (collisions_swap_nodes) {
                        this->collision_node_swap(rng, next_height_idx, next_height_idx - 1);
                    }
                }
                this->node_heights_.at(intervening_indices.back())->set_value(new_height);
                return;
            }
            // Younger nodes to bump down
            for (auto next_height_idx : intervening_indices) {
                this->node_heights_.at(next_height_idx + 1)->set_value(
                        this->node_heights_.at(next_height_idx)->get_value()
                        );
                if (collisions_swap_nodes) {
                    this->collision_node_swap(rng, next_height_idx + 1, next_height_idx);
                }
            }
            this->node_heights_.at(intervening_indices.back())->set_value(new_height);
            return;
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
            this->store_topology();
        }
        virtual void store_all_heights() {
            this->root_->store_all_heights();
        }
        virtual void store_all_height_pointers() {
            this->root_->store_all_height_pointers();
        }
        virtual void store_topology() {
            this->stored_root_ = this->root_->get_copy();
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
            this->restore_topology();
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
        virtual void restore_topology() {
            this->root_ = this->stored_root_;
        }

        std::string to_parentheses() const {
            return this->root_->to_parentheses();
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

#endif
