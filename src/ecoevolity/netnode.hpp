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

#ifndef ECOEVOLITY_NETNODE_HPP
#define ECOEVOLITY_NETNODE_HPP

#include <memory>

#include "basenetnode.hpp"
#include "parameter.hpp"
#include "rng.hpp"

/**
 * Base class for a node of a phylogenetic network.
 *
 * @note    Many of the class' methods modified from:
 *              BasicTNode class of Bio++ Library
 *              <http://biopp.univ-montp2.fr/wiki/index.php/Main_Page>
 *              License:    CeCILL <http://www.cecill.info>
 *              Author:     Sylvain Gaillard
 *              Copyright:  CNRS, (January 12, 2011)
 */
class NetNode: public BaseNetNode<Node>{
    private:
        typedef BaseNetNode<NetNode> BaseClass;

    public:
        NetNode() { }
        NetNode(int index) : BaseClass(index) { }
        NetNode(std::string label) : BaseClass(label) { }
        NetNode(int index, std::string label) : BaseClass(index, label) { }
        NetNode(double height) : BaseClass(height) { }
        NetNode(int index, double height) : BaseClass(index, height) { }
        NetNode(std::shared_ptr<PositiveRealParameter> height) : BaseClass(height) { }
        NetNode(int index, std::shared_ptr<PositiveRealParameter> height) : BaseClass(index, height) { }
        NetNode(std::string label, double height) : BaseClass(label, height) { }
        NetNode(int index, std::string label, double height) : BaseClass(index, label, height) { }
};


class PopulationNetNode: public BaseNetNode<PopulationNetNode>{
    protected:
        typedef BaseNetNode<PopulationNetNode> BaseClass;
        int population_index_ = -1;
        BiallelicPatternProbabilityMatrix bottom_pattern_probs_;
        BiallelicPatternProbabilityMatrix top_pattern_probs_;
        std::shared_ptr<PositiveRealParameter> population_size_ = std::make_shared<PositiveRealParameter>(0.001);
        std::shared_ptr<PositiveRealParameter> stored_population_size_ = std::make_shared<PositiveRealParameter>(0.001);

        void add_ln_relative_population_size_prior_density(
                double& density,
                std::vector< std::shared_ptr<PositiveRealParameter> >& parameters,
                std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes) const {
            visited_nodes.insert(this->shared_from_this());
            bool parameter_found = false;
            for (auto parameter_iter : parameters) {
                if (parameter_iter == this->population_size_) {
                    parameter_found = true;
                    break;
                }
            }
            if (! parameter_found) {
                density += this->population_size_->relative_prior_ln_pdf();
                parameters.push_back(this->population_size_);
            }
            for (unsigned int i = 0; i < this->children_.size(); ++i) {
                if (visited_nodes.count(this->children_.at(i)) < 1) {
                    this->children_.at(i)->add_ln_relative_population_size_prior_density(density, parameters, visited_nodes);
                }
            }
        }

        void get_all_population_size_parameters(
                std::vector< std::shared_ptr<PositiveRealParameter> >& parameters,
                std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes) const {
            visited_nodes.insert(this->shared_from_this());
            bool parameter_found = false;
            for (auto parameter_iter : parameters) {
                if (parameter_iter == this->population_size_) {
                    parameter_found = true;
                    break;
                }
            }
            if (! parameter_found) {
                parameters.push_back(this->population_size_);
            }
            for (unsigned int i = 0; i < this->children_.size(); ++i) {
                if (visited_nodes.count(this->children_.at(i)) < 1) {
                    this->children_.at(i)->get_all_population_size_parameters(
                            parameters, visited_nodes);
                }
            }
        }

        void scale_all_population_sizes(
                double scale,
                std::vector< std::shared_ptr<PositiveRealParameter> >& parameters,
                std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes,
                unsigned int & number_of_free_parameters_scaled) {
            visited_nodes.insert(this->shared_from_this());
            bool parameter_found = false;
            for (auto parameter_iter : parameters) {
                if (parameter_iter == this->population_size_) {
                    parameter_found = true;
                    break;
                }
            }
            if (! parameter_found) {
                parameters.push_back(this->population_size_);
                if (! this->population_size_is_fixed()) {
                    this->set_population_size(this->get_population_size() * scale);
                    ++number_of_free_parameters_scaled;
                }
            }
            for (unsigned int i = 0; i < this->children_.size(); ++i) {
                if (visited_nodes.count(this->children_.at(i)) < 1) {
                    this->children_.at(i)->scale_all_population_sizes(
                            scale,
                            parameters,
                            visited_nodes
                            number_of_free_parameters_scaled);
                }
            }
        }

    public:
        PopulationNetNode() { }
        PopulationNetNode(std::string label) : BaseClass(label) { }
        PopulationNetNode(double height) : BaseClass(height) { }
        PopulationNetNode(std::shared_ptr<PositiveRealParameter> height) : BaseClass(height) { }
        PopulationNetNode(std::string label, double height) :
            BaseClass(label, height)
            { }
        PopulationNetNode(unsigned int allele_count) : BaseClass() {
            this->bottom_pattern_probs_.resize(allele_count);
            this->top_pattern_probs_.resize(allele_count);
        }
        PopulationNetNode(std::string label, unsigned int allele_count) :
            BaseClass(label)
        {
            this->bottom_pattern_probs_.resize(allele_count);
            this->top_pattern_probs_.resize(allele_count);
        }
        PopulationNetNode(double height, unsigned int allele_count) :
            BaseClass(height)
        {
            this->bottom_pattern_probs_.resize(allele_count);
            this->top_pattern_probs_.resize(allele_count);
        }
        PopulationNetNode(std::string label,
                double height,
                unsigned int allele_count) :
            BaseClass(label, height)
        {
            this->bottom_pattern_probs_.resize(allele_count);
            this->top_pattern_probs_.resize(allele_count);
        }
        PopulationNetNode(
                int population_index,
                double height) :
            BaseClass(population_index, height)
        { }
        PopulationNetNode(
                int population_index,
                std::string label,
                double height) :
            BaseClass(population_index, label, height)
        { }
        PopulationNetNode(
                int population_index,
                std::string label,
                double height,
                unsigned int allele_count) :
            BaseClass(population_index, label, height)
        {
            this->bottom_pattern_probs_.resize(allele_count);
            this->top_pattern_probs_.resize(allele_count);
        }
        PopulationNetNode(const PopulationNetNode& node) :
            BaseClass(node.index_, node.label_, node.height_)
        {
            this->population_size_ = node.population_size_;
            this->bottom_pattern_probs_ = node.bottom_pattern_probs_;
            this->top_pattern_probs_ = node.top_pattern_probs_;
        }

        // overload copy operator
        // PopulationNetNode& operator=(const PopulationNetNode& node) {
        //     this->children_ = node.children_;
        //     this->parent_ = node.parent_;
        //     this->height_->set_value(node.height_->get_value());
        //     this->label_ = node.label_;
        //     this->is_dirty_ = node.is_dirty_;
        //     this->bottom_pattern_probs_ = node.bottom_pattern_probs_;
        //     this->top_pattern_probs_ = node.top_pattern_probs_;
        //     return * this;
        // }
        void copy_node_type_specific_members(std::shared_ptr<PopulationNetNode> copy) const {
            copy->population_size_ = this->population_size_;
            copy->stored_population_size_ = this->stored_population_size_;
        }
        void deep_copy_node_type_specific_members(std::shared_ptr<PopulationNetNode> copy) const {
            copy->population_size_ = std::make_shared<PositiveRealParameter>(*this->population_size_);
            copy->stored_population_size_ = std::make_shared<PositiveRealParameter>(*this->stored_population_size_);
        }

        void get_parameter_map(
                std::map<std::string, double> & parameter_map) const {
            BaseNetNode::get_parameter_map(parameter_map);
            parameter_map["pop_size"] = this->get_population_size();
        }

        /**
         * Method to populate non-height related data (e.g., pop size) from
         * node comments in the form of a map. Overridding from BaseNetNode.
         */
        void extract_data_from_node_comments(
                const std::map<std::string, std::string> & comment_map) {
            if (comment_map.count("pop_size") > 0) {
                double pop_size;
                std::stringstream s_converter(comment_map.at("pop_size"));
                if (! (s_converter >> pop_size)) {
                    throw EcoevolityError("could not convert pop_size \'" +
                            s_converter.str() + "\'");
                }
                this->set_population_size(pop_size);
            }
        }

        // methods for accessing/changing pattern probabilities
        unsigned int get_allele_count() const {
            return this->bottom_pattern_probs_.get_allele_count();
        }

        unsigned int get_leaf_allele_count(std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes) const {
            visited_nodes.insert(this->shared_from_this());
            if (this->is_leaf()) {
                return this->get_allele_count();
            }
            unsigned int n = 0;
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    n += child_iter->get_leaf_allele_count(visited_nodes);
                }
            }
            return n;
        }
        unsigned int get_leaf_allele_count() const {
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            return this->get_leaf_allele_count(visited_nodes);
        }

        void resize(unsigned int allele_count) {
            this->bottom_pattern_probs_.resize(allele_count);
            this->top_pattern_probs_.resize(allele_count);
        }
        void reset(unsigned int allele_count) {
            this->bottom_pattern_probs_.reset(allele_count);
            this->top_pattern_probs_.reset(allele_count);
        }
        void resize_all(std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes) {
            visited_nodes.insert(this->shared_from_this());
            this->resize(this->get_leaf_allele_count());
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->resize_all(visited_nodes);
                }
            }
        }
        void resize_all() {
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            this->resize_all(visited_nodes);
        }

        std::shared_ptr<PopulationNetNode> get_clone() const {
            return std::make_shared<PopulationNetNode>(*this);
        }

        // std::shared_ptr<PopulationNetNode> get_clade_clone() const {
        //     if (this->is_leaf()) {
        //         return this->get_clone();
        //     }
        //     std::shared_ptr<PopulationNetNode> n = this->get_clone();
        //     for (auto child_iter: this->children_) {
        //         n->add_child(child_iter->get_clade_clone());
        //     }
        //     return n;
        // }

        const BiallelicPatternProbabilityMatrix& get_bottom_pattern_probs() const{
            return this->bottom_pattern_probs_;
        }
        const BiallelicPatternProbabilityMatrix& get_top_pattern_probs() const{
            return this->top_pattern_probs_;
        }
        BiallelicPatternProbabilityMatrix* clone_bottom_pattern_probs() const{
            return this->bottom_pattern_probs_.clone();
        }
        BiallelicPatternProbabilityMatrix* clone_top_pattern_probs() const{
            return this->top_pattern_probs_.clone();
        }

        void copy_bottom_pattern_probs(const BiallelicPatternProbabilityMatrix& m) {
            // No check here; bottom probs will be updated first and can differ
            // in size until the top is also updated.
            // if (m.get_allele_count() != this->top_pattern_probs_.get_allele_count()) {
            //     throw EcoevolityError(
            //             "PopulationNetNode:copy_bottom_pattern_probs(); allele counts must be the same between top and bottom of branch");
            // }
            this->bottom_pattern_probs_.copy(m);
        }
        void copy_top_pattern_probs(const BiallelicPatternProbabilityMatrix& m) {
            if (m.get_allele_count() != this->bottom_pattern_probs_.get_allele_count()) {
                throw EcoevolityError(
                        "PopulationNetNode:copy_top_pattern_probs(); allele counts must be the same between top and bottom of branch");
            }
            this->top_pattern_probs_.copy(m);
        }
        void copy_pattern_probs(
                const BiallelicPatternProbabilityMatrix& bottom_probs,
                const BiallelicPatternProbabilityMatrix& top_probs) {
            if (bottom_probs.get_allele_count() != top_probs.get_allele_count()) {
                throw EcoevolityError(
                        "PopulationNetNode:copy_pattern_probs(); allele counts must be the same between top and bottom of branch");
            }
            this->bottom_pattern_probs_.copy(bottom_probs);
            this->top_pattern_probs_.copy(top_probs);
        }

        double get_bottom_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count) const {
            return this->bottom_pattern_probs_.get_pattern_probability(
                    allele_count,
                    red_allele_count);
        }
        double get_top_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count) const {
            return this->top_pattern_probs_.get_pattern_probability(
                    allele_count,
                    red_allele_count);
        }
        void set_bottom_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count,
                double probability) {
            this->bottom_pattern_probs_.set_pattern_probability(
                    allele_count,
                    red_allele_count,
                    probability);
        }
        void set_top_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count,
                double probability) {
            this->top_pattern_probs_.set_pattern_probability(
                    allele_count,
                    red_allele_count,
                    probability);
        }

        double get_population_size() const {
            return this->population_size_->get_value();
        }
        std::shared_ptr<PositiveRealParameter> get_population_size_parameter() const {
            return this->population_size_;
        }
        std::vector< std::shared_ptr<PositiveRealParameter> > get_all_population_size_parameters() const {
            std::vector< std::shared_ptr<PositiveRealParameter> > parameters;
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            // parameters.reserve(this->get_node_count());
            this->get_all_population_size_parameters(
                    parameters, visited_nodes);
            return parameters;
        }

        void set_population_size_parameter(std::shared_ptr<PositiveRealParameter> size) {
            this->population_size_ = size;
            this->make_all_dirty();
        }
        void set_all_population_size_parameters(std::shared_ptr<PositiveRealParameter> size,
                std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes) {
            visited_nodes.insert(this->shared_from_this());
            this->population_size_ = size;
            this->make_dirty();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->set_all_population_size_parameters(size, visited_nodes);
                }
            }
        }
        void set_all_population_size_parameters(std::shared_ptr<PositiveRealParameter> size) {
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            this->set_all_population_size_parameters(size, visited_nodes);
        }
        void set_all_population_size_parameters() {
            std::shared_ptr<PositiveRealParameter> size = this->population_size_;
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            this->set_all_population_size_parameters(size, visited_nodes);
        }

        void set_population_size(double size) {
            this->population_size_->set_value(size);
            this->make_all_dirty();
        }
        void set_all_population_sizes(double size,
                std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes) {
            visited_nodes.insert(this->shared_from_this());
            this->population_size_->set_value(size);
            this->make_dirty();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->set_all_population_sizes(size, visited_nodes);
                }
            }
        }
        void set_all_population_sizes(double size) {
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            this->set_all_population_sizes(size, visited_nodes);
        }
        void update_population_size(double size) {
            this->population_size_->update_value(size);
            this->make_all_dirty();
        }
        void update_all_population_sizes(double size,
                std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes) {
            visited_nodes.insert(this->shared_from_this());
            this->population_size_->update_value(size);
            this->make_dirty();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->update_all_population_sizes(size, visited_nodes);
                }
            }
        }
        void update_all_population_sizes(double size) {
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            this->update_all_population_sizes(size, visited_nodes);
        }

        void store_population_size() {
            this->population_size_->store();
        }
        void restore_population_size() {
            this->population_size_->restore();
            this->make_all_dirty();
        }
        void store_all_population_sizes(std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes) {
            visited_nodes.insert(this->shared_from_this());
            this->population_size_->store();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->store_all_population_sizes(visited_nodes);
                }
            }
        }
        void restore_all_population_sizes(std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes) {
            visited_nodes.insert(this->shared_from_this());
            this->population_size_->restore();
            this->make_dirty();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->restore_all_population_sizes(visited_nodes);
                }
            }
        }
        void store_all_population_sizes() {
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            this->store_all_population_sizes(visited_nodes);
        }
        void restore_all_population_sizes() {
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            this->restore_all_population_sizes(visited_nodes);
        }

        void store_population_size_pointer() {
            this->stored_population_size_ = this->population_size_;
        }
        void store_all_population_size_pointers(std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes) {
            visited_nodes.insert(this->shared_from_this());
            this->stored_population_size_ = this->population_size_;
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->store_all_population_size_pointers(visited_nodes);
                }
            }
        }
        void store_all_population_size_pointers() {
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            this->store_all_population_size_pointers(visited_nodes);
        }

        void restore_population_size_pointer() {
            this->population_size_ = this->stored_population_size_;
            this->make_all_dirty();
        }
        void restore_all_population_size_pointers(std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes) {
            visited_nodes.insert(this->shared_from_this());
            this->population_size_ = this->stored_population_size_;
            this->make_dirty();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->restore_all_population_size_pointers(visited_nodes);
                }
            }
        }
        void restore_all_population_size_pointers() {
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            this->restore_all_population_size_pointers(visited_nodes);
        }

        void store_all_parameter_values() {
            this->store_all_heights();
            this->store_all_population_sizes();
        }
        void restore_all_parameter_values() {
            this->restore_all_heights();
            this->restore_all_population_sizes();
        }
        void store_all_parameter_pointers() {
            this->store_all_height_pointers();
            this->store_all_population_size_pointers();
        }
        void restore_all_parameter_pointers() {
            this->restore_all_height_pointers();
            this->restore_all_population_size_pointers();
        }

        void set_population_size_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->population_size_->set_prior(prior);
            this->make_all_dirty();
        }
        void set_all_population_size_priors(std::shared_ptr<ContinuousProbabilityDistribution> prior,
                std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes) {
            visited_nodes.insert(this->shared_from_this());
            this->population_size_->set_prior(prior);
            this->make_dirty();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->set_all_population_size_priors(prior, visited_nodes);
                }
            }
        }
        void set_all_population_size_priors(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            this->set_all_population_size_priors(prior, visited_nodes);
        }

        void fix_population_size() {
            this->population_size_->fix();
            this->make_dirty();
        }
        void fix_all_population_sizes(
                std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes) {
            visited_nodes.insert(this->shared_from_this());
            this->population_size_->fix();
            this->make_dirty();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->fix_all_population_sizes(visited_nodes);
                }
            }
        }
        void fix_all_population_sizes() {
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            this->fix_all_population_sizes(visited_nodes);
        }

        void estimate_population_size() {
            this->population_size_->estimate();
            this->make_dirty();
        }
        void estimate_all_population_sizes(
                std::set< std::shared_ptr<PopulationNetNode> > & visited_nodes) {
            visited_nodes.insert(this->shared_from_this());
            this->population_size_->estimate();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->estimate_all_population_sizes(visited_nodes);
                }
            }
        }
        void estimate_all_population_sizes() {
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            this->estimate_all_population_sizes(visited_nodes);
        }

        ////////////////////////////////////////

        double calculate_ln_relative_population_size_prior_density() const {
            double d = 0.0;
            std::vector< std::shared_ptr<PositiveRealParameter> > parameters;
            parameters.reserve(this->get_node_count());
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            this->add_ln_relative_population_size_prior_density(d, parameters, visited_nodes);
            return d;
        }

        unsigned int scale_all_population_sizes(double scale) {
            std::vector< std::shared_ptr<PositiveRealParameter> > parameters;
            parameters.reserve(this->get_node_count());
            unsigned int number_of_free_parameters_scaled = 0;
            std::set< std::shared_ptr<PopulationNetNode> > visited_nodes;
            this->scale_all_population_sizes(
                    scale,
                    parameters,
                    visited_nodes,
                    number_of_free_parameters_scaled);
            return number_of_free_parameters_scaled;
        }

        bool population_size_is_fixed() const {
            return this->population_size_->is_fixed();
        }

        bool all_population_sizes_are_fixed() const {
            if (! this->population_size_is_fixed()) {
                return false;
            }
            for (auto child_iter: this->children_) {
                if (! child_iter->all_population_sizes_are_fixed()) {
                    return false;
                }
            }
            return true;
        }

        double get_population_size_relative_prior_ln_pdf() const {
            return this->population_size_->relative_prior_ln_pdf();
        }

        void get_node_indices(std::vector<unsigned int> & internal_indices,
                std::vector<unsigned int> & leaf_indices) {
            if (this->is_leaf()) {
                leaf_indices.push_back(this->get_index());
            } else {
                internal_indices.push_back(this->get_index());
            }
            for (auto child_iter: this->children_) {
                child_iter->get_node_indices(internal_indices, leaf_indices);
            }
        }

        // Overriding this method from BaseNetNode to return pop size
        std::string get_additional_comment_data_string(
                unsigned int precision = 12) const {
            std::ostringstream s;
            s.precision(precision);
            s << "pop_size="
              << this->get_population_size();
            return s.str();
        }

        // Overriding this method from BaseNetNode to make sure newly inserted
        // nodes are fully initialized
        void finish_initializing_inserted_internal_node(
                RandomNumberGenerator & rng) {
            ECOEVOLITY_ASSERT(this->has_parent() && this->has_children());
            if (this->get_parent()->get_population_size_parameter() == this->children_.at(0)->population_size_) {
                this->population_size_ = this->children_.at(0)->population_size_;
                return;
            }
            this->population_size_->estimate();
            this->population_size_->set_prior(this->children_.at(0)->population_size_->get_prior());
            if (this->children_.at(0)->population_size_->is_fixed()) {
                double child_pop_size_sum = 0.0;
                for (auto child : this->children_) {
                    child_pop_size_sum += child->get_population_size();
                }
                double new_pop_size = child_pop_size_sum / this->children_.size();
                this->population_size_->set_value(new_pop_size);
                this->population_size_->fix();
                return;
            }
            // double child_max_pop_size = -1.0;
            // double child_min_pop_size = std::numeric_limits<double>::max();
            // for (auto child : this->children_) {
            //     double pop_size = child->get_population_size();
            //     if (pop_size < child_min_pop_size) {
            //         child_min_pop_size = pop_size;
            //     }
            //     if (pop_size > child_max_pop_size) {
            //         child_max_pop_size = pop_size;
            //     }
            // }
            // double pop_size_diff = child_max_pop_size - child_min_pop_size;
            // if (pop_size_diff == 0.0) {
            //     this->population_size_->set_value(child_max_pop_size);
            //     return;
            // }
            // double new_pop_size = rng.uniform_real(
            //         child_min_pop_size,
            //         child_max_pop_size);
            // this->population_size_->set_value(new_pop_size);
            this->population_size_->set_value_from_prior(rng);
        }

        virtual double get_ln_prob_of_drawing_state() {
            ECOEVOLITY_ASSERT(this->has_parent() && this->has_children());
            if (this->population_size_->is_fixed()) {
                return 0.0;
            }
            if (this->children_.at(0)->population_size_ == this->population_size_) {
                return 0.0;
            }
            // double child_max_pop_size = -1.0;
            // double child_min_pop_size = std::numeric_limits<double>::max();
            // for (auto child : this->children_) {
            //     double pop_size = child->get_population_size();
            //     if (pop_size < child_min_pop_size) {
            //         child_min_pop_size = pop_size;
            //     }
            //     if (pop_size > child_max_pop_size) {
            //         child_max_pop_size = pop_size;
            //     }
            // }
            // double pop_size_diff = child_max_pop_size - child_min_pop_size;
            // if (pop_size_diff == 0.0) {
            //     return std::numeric_limits<double>::infinity();
            // }
            // return -std::log(pop_size_diff);
            return this->population_size_->prior_ln_pdf();
        }
};

#endif
