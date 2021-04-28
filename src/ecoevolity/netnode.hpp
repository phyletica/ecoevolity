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
class NetNode: public BaseNetNode<NetNode>{
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
        Probability parent_inheritance_proportion_ = 1.0;
        std::vector<BiallelicPatternProbabilityMatrix> bottom_pattern_probs_ = {BiallelicPatternProbabilityMatrix()};
        std::vector<BiallelicPatternProbabilityMatrix> top_pattern_probs_ = {BiallelicPatternProbabilityMatrix()};
        std::vector<std::shared_ptr<PositiveRealParameter> > population_sizes_ = { std::make_shared<PositiveRealParameter>(0.001) };
        std::vector<std::shared_ptr<PositiveRealParameter> > stored_population_sizes_;

        void add_ln_relative_population_size_prior_density(
                double& density,
                std::vector< std::shared_ptr<PositiveRealParameter> >& parameters,
                std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) const {
            this->visit(visited_nodes);
            for (auto pop_size : this->population_sizes_) {
                bool parameter_found = false;
                for (auto parameter_iter : parameters) {
                    if (parameter_iter == pop_size) {
                        parameter_found = true;
                        break;
                    }
                }
                if (! parameter_found) {
                    density += pop_size->relative_prior_ln_pdf();
                    parameters.push_back(pop_size);
                }
            }
            for (unsigned int i = 0; i < this->children_.size(); ++i) {
                if (visited_nodes.count(this->children_.at(i)) < 1) {
                    this->children_.at(i)->add_ln_relative_population_size_prior_density(density, parameters, visited_nodes);
                }
            }
        }

        void get_all_population_size_parameters(
                std::vector< std::shared_ptr<PositiveRealParameter> >& parameters,
                std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) const {
            this->visit(visited_nodes);
            for (auto pop_size : this->population_sizes_) {
                bool parameter_found = false;
                for (auto parameter_iter : parameters) {
                    if (parameter_iter == pop_size) {
                        parameter_found = true;
                        break;
                    }
                }
                if (! parameter_found) {
                    parameters.push_back(pop_size);
                }
                for (unsigned int i = 0; i < this->children_.size(); ++i) {
                    if (visited_nodes.count(this->children_.at(i)) < 1) {
                        this->children_.at(i)->get_all_population_size_parameters(
                                parameters, visited_nodes);
                    }
                }
            }
        }

        void scale_all_population_sizes(
                double scale,
                std::vector< std::shared_ptr<PositiveRealParameter> >& parameters,
                std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes,
                unsigned int & number_of_free_parameters_scaled) {
            this->visit(visited_nodes);
            for (unsigned int branch_idx = 0;
                    branch_idx < this->population_sizes_.size();
                    ++branch_idx) {
                bool parameter_found = false;
                for (auto parameter_iter : parameters) {
                    if (parameter_iter == this->population_sizes_.at(branch_idx)) {
                        parameter_found = true;
                        break;
                    }
                }
                if (! parameter_found) {
                    parameters.push_back(this->population_sizes_.at(branch_idx));
                    if (! this->population_size_is_fixed()) {
                        this->set_population_size(branch_idx, this->get_population_size(branch_idx) * scale);
                        ++number_of_free_parameters_scaled;
                    }
                }
            }
            for (unsigned int i = 0; i < this->children_.size(); ++i) {
                if (visited_nodes.count(this->children_.at(i)) < 1) {
                    this->children_.at(i)->scale_all_population_sizes(
                            scale,
                            parameters,
                            visited_nodes,
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
            for (auto bpp : this->bottom_pattern_probs_) {
                bpp.resize(allele_count);
            }
            for (auto tpp : this->top_pattern_probs_) {
                tpp.resize(allele_count);
            }
        }
        PopulationNetNode(std::string label, unsigned int allele_count) :
            BaseClass(label)
        {
            for (auto bpp : this->bottom_pattern_probs_) {
                bpp.resize(allele_count);
            }
            for (auto tpp : this->top_pattern_probs_) {
                tpp.resize(allele_count);
            }
        }
        PopulationNetNode(double height, unsigned int allele_count) :
            BaseClass(height)
        {
            for (auto bpp : this->bottom_pattern_probs_) {
                bpp.resize(allele_count);
            }
            for (auto tpp : this->top_pattern_probs_) {
                tpp.resize(allele_count);
            }
        }
        PopulationNetNode(std::string label,
                double height,
                unsigned int allele_count) :
            BaseClass(label, height)
        {
            for (auto bpp : this->bottom_pattern_probs_) {
                bpp.resize(allele_count);
            }
            for (auto tpp : this->top_pattern_probs_) {
                tpp.resize(allele_count);
            }
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
            for (auto bpp : this->bottom_pattern_probs_) {
                bpp.resize(allele_count);
            }
            for (auto tpp : this->top_pattern_probs_) {
                tpp.resize(allele_count);
            }
        }
        PopulationNetNode(const PopulationNetNode& node) :
            BaseClass(node.index_, node.label_, node.height_)
        {
            this->population_sizes_ = node.population_sizes_;
            this->stored_population_sizes_ = node.stored_population_sizes_;
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
            copy->population_sizes_ = this->population_sizes_;
            copy->stored_population_sizes_ = this->stored_population_sizes_;
        }
        void deep_copy_node_type_specific_members(std::shared_ptr<PopulationNetNode> copy) const {
            copy->population_sizes_.clear();
            copy->stored_population_sizes_.clear();
            for (auto pop_size : this->population_sizes_) {
                copy->population_sizes_.push_back(std::make_shared<PositiveRealParameter>(*pop_size));
            }
            for (auto stored_pop_size : this->stored_population_sizes_) {
                copy->stored_population_sizes_.push_back(std::make_shared<PositiveRealParameter>(*stored_pop_size));
            }
        }

        void get_parameter_map(
                std::map<std::string, double> & parameter_map) const {
            BaseNetNode::get_parameter_map(parameter_map);
            for (unsigned int i = 0; i < this->population_sizes_.size(); ++i) {
                parameter_map["pop_size" + std::to_string(i)] = this->get_population_size(i);
            }
        }

        /**
         * Method to populate non-height related data (e.g., pop size) from
         * node comments in the form of a map. Overridding from BaseNetNode.
         */
        void extract_data_from_node_comments(
                const std::map<std::string, std::string> & comment_map) {
            for (unsigned int i = 0; i < this->population_sizes_.size(); ++i) {
                std::string key = "pop_size" + std::to_string(i);
                if (comment_map.count(key) > 0) {
                    double pop_size;
                    std::stringstream s_converter(comment_map.at(key));
                    if (! (s_converter >> pop_size)) {
                        throw EcoevolityError("could not convert pop_size \'" +
                                s_converter.str() + "\'");
                    }
                    this->set_population_size(i, pop_size);
                }
            }
        }

        // Overriding adding/removing parents to make sure we can't have more
        // than two parents and that the parent_inheritance_proportion gets
        // updated
        void add_parent(std::shared_ptr<PopulationNetNode> node) {
            if (this->get_number_of_parents() > 1) {
                throw EcoevolityError("PopulationNetNode::add_parent(), tried to add 3rd parent");
            }
            if (this->has_parent()) {
                this->parent_inheritance_proportion_.set_value(0.5);
                ECOEVOLITY_ASSERT(this->population_sizes_.size() == 1);
                ECOEVOLITY_ASSERT(this->bottom_pattern_probs_.size() == 1);
                ECOEVOLITY_ASSERT(this->top_pattern_probs_.size() == 1);
                this->population_sizes_.push_back(this->get_population_size_parameter(0));
                this->bottom_pattern_probs_.push_back(BiallelicPatternProbabilityMatrix(this->get_allele_count()));
                this->top_pattern_probs_.push_back(BiallelicPatternProbabilityMatrix(this->get_allele_count()));
            }
            BaseNetNode::add_parent(node);
        }
        void add_parent(std::shared_ptr<PopulationNetNode> node, double inheritance_proportion) {
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            if (! this->has_parent()) {
                ECOEVOLITY_ASSERT(inheritance_proportion == 1.0);
                this->parent_inheritance_proportion_.set_value(1.0);
                this->add_parent(node);
                return;
            }
            this->parent_inheritance_proportion_.set_value(1.0 - inheritance_proportion);
            this->add_parent(node);
        }

        std::shared_ptr<PopulationNetNode> remove_parent() {
            this->parent_inheritance_proportion_.set_value(1.0);
            // Always keep one pop size parameter and bottom/top pattern probs
            // so no need to remove any of these. The call to remove_parent
            // asserts that there is currently less than 2 parents
            return BaseNetNode::remove_parent();
        }
        int remove_parent(std::shared_ptr<PopulationNetNode> node) {
            this->parent_inheritance_proportion_.set_value(1.0);
            unsigned int parent_idx = BaseNetNode::remove_parent(node);
            this->population_sizes_.at(parent_idx).reset();
            this->population_sizes_.erase(this->population_sizes_.begin() + parent_idx);
            this->bottom_pattern_probs_.erase(this->bottom_pattern_probs_.begin() + parent_idx);
            this->top_pattern_probs_.erase(this->top_pattern_probs_.begin() + parent_idx);
        }
        std::shared_ptr<PopulationNetNode> remove_parent(const unsigned int index) {
            this->parent_inheritance_proportion_.set_value(1.0);
            return BaseNetNode::remove_parent(index);
            this->population_sizes_.at(index).reset();
            this->population_sizes_.erase(this->population_sizes_.begin() + index);
            this->bottom_pattern_probs_.erase(this->bottom_pattern_probs_.begin() + index);
            this->top_pattern_probs_.erase(this->top_pattern_probs_.begin() + index);
        }

        double get_inheritance_proportion(unsigned int parent_index) const {
            if (parent_index == 0) {
                return this->parent_inheritance_proportion_.get_value();
            }
            if (parent_index == 1) {
                return (1.0 - this->parent_inheritance_proportion_.get_value());
            }
            throw EcoevolityError("PopulationNetNode::get_inheritance_proportion() called with index greater than 1");
        }
        void set_inheritance_proportion(unsigned int parent_index, double proportion) {
            if (parent_index == 0) {
                this->parent_inheritance_proportion_.set_value(proportion);
                return;
            }
            if (parent_index == 1) {
                this->parent_inheritance_proportion_.set_value(1.0 - proportion);
                return;
            }
            throw EcoevolityError("PopulationNetNode::set_inheritance_proportion() called with index greater than 1");
        }

        // methods for accessing/changing pattern probabilities
        unsigned int get_allele_count() const {
            return this->bottom_pattern_probs_.at(0).get_allele_count();
        }

        unsigned int get_leaf_allele_count(std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) const {
            this->visit(visited_nodes);
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
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            return this->get_leaf_allele_count(visited_nodes);
        }

        void resize(unsigned int allele_count) {
            for (auto bpp : this->bottom_pattern_probs_) {
                bpp.resize(allele_count);
            }
            for (auto tpp : this->top_pattern_probs_) {
                tpp.resize(allele_count);
            }
        }
        void reset(unsigned int allele_count) {
            for (auto bpp : this->bottom_pattern_probs_) {
                bpp.reset(allele_count);
            }
            for (auto tpp : this->top_pattern_probs_) {
                tpp.reset(allele_count);
            }
        }
        void resize_all(std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) {
            this->visit(visited_nodes);
            this->resize(this->get_leaf_allele_count());
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->resize_all(visited_nodes);
                }
            }
        }
        void resize_all() {
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            this->resize_all(visited_nodes);
        }

        std::shared_ptr<PopulationNetNode> get_clone() const {
            return std::make_shared<PopulationNetNode>(*this);
        }

        std::shared_ptr<PopulationNetNode> get_clade_clone(
                std::map< std::shared_ptr<const PopulationNetNode>, std::shared_ptr<PopulationNetNode> > & visited_nodes
                ) const {
            std::shared_ptr<PopulationNetNode> n = this->get_clone();
            if (this->has_multiple_parents()) {
                ECOEVOLITY_ASSERT(visited_nodes.count(this->shared_from_this()) < 1);
                visited_nodes[this->shared_from_this()] = n;
            }
            if (this->is_leaf()) {
                return n;
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) > 0) {
                    n->add_child(visited_nodes.at(child_iter));
                }
                else {
                    n->add_child(child_iter->get_clade_clone());
                }
            }
            return n;
        }
        std::shared_ptr<PopulationNetNode> get_clade_clone() const {
            std::map< std::shared_ptr<const PopulationNetNode>, std::shared_ptr<PopulationNetNode> > visited_nodes;
            return this->get_clade_clone(visited_nodes);
        }

        const BiallelicPatternProbabilityMatrix& get_bottom_pattern_probs(unsigned int branch_idx) const{
            return this->bottom_pattern_probs_.at(branch_idx);
        }
        const BiallelicPatternProbabilityMatrix& get_bottom_pattern_probs() const{
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            return this->get_bottom_pattern_probs(0);
        }
        const BiallelicPatternProbabilityMatrix& get_top_pattern_probs(unsigned int branch_idx) const{
            return this->top_pattern_probs_.at(branch_idx);
        }
        const BiallelicPatternProbabilityMatrix& get_top_pattern_probs() const{
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            return this->get_top_pattern_probs(0);
        }
        BiallelicPatternProbabilityMatrix* clone_bottom_pattern_probs(unsigned int branch_idx) const{
            return this->bottom_pattern_probs_.at(branch_idx).clone();
        }
        BiallelicPatternProbabilityMatrix* clone_bottom_pattern_probs() const{
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            return clone_bottom_pattern_probs(0);
        }
        BiallelicPatternProbabilityMatrix* clone_top_pattern_probs(unsigned int branch_idx) const{
            return this->top_pattern_probs_.at(branch_idx).clone();
        }
        BiallelicPatternProbabilityMatrix* clone_top_pattern_probs() const{
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            return clone_top_pattern_probs(0);
        }

        void copy_bottom_pattern_probs(unsigned int branch_idx, const BiallelicPatternProbabilityMatrix& m) {
            this->bottom_pattern_probs_.at(branch_idx)..copy(m);
        }
        void copy_bottom_pattern_probs(const BiallelicPatternProbabilityMatrix& m) {
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            this->copy_bottom_pattern_probs(0, m);
        }
        void copy_top_pattern_probs(unsigned int branch_idx, const BiallelicPatternProbabilityMatrix& m) {
            if (m.get_allele_count() != this->bottom_pattern_probs_.at(branch_idx).get_allele_count()) {
                throw EcoevolityError(
                        "PopulationNetNode:copy_top_pattern_probs(); allele counts must be the same between top and bottom of branch");
            }
            this->top_pattern_probs_.at(branch_idx).copy(m);
        }
        void copy_top_pattern_probs(const BiallelicPatternProbabilityMatrix& m) {
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            this->copy_top_pattern_probs(0, m);
        }
        void copy_pattern_probs(unsigned int branch_idx,
                const BiallelicPatternProbabilityMatrix& bottom_probs,
                const BiallelicPatternProbabilityMatrix& top_probs) {
            if (bottom_probs.get_allele_count() != top_probs.get_allele_count()) {
                throw EcoevolityError(
                        "PopulationNetNode:copy_pattern_probs(); allele counts must be the same between top and bottom of branch");
            }
            this->bottom_pattern_probs_.at(branch_idx).copy(bottom_probs);
            this->top_pattern_probs_.at(branch_idx).copy(top_probs);
        }
        void copy_pattern_probs(
                const BiallelicPatternProbabilityMatrix& bottom_probs,
                const BiallelicPatternProbabilityMatrix& top_probs) {
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            this->copy_pattern_probs(0, bottom_probs, top_probs);
        }

        double get_bottom_pattern_probability(
                unsigned int branch_idx,
                unsigned int allele_count,
                unsigned int red_allele_count) const {
            return this->bottom_pattern_probs_.at(branch_idx).get_pattern_probability(
                    allele_count,
                    red_allele_count);
        }
        double get_bottom_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count) const {
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            return this->get_bottom_pattern_probability(0, allele_count, red_allele_count);
        }
        double get_top_pattern_probability(
                unsigned int branch_idx,
                unsigned int allele_count,
                unsigned int red_allele_count) const {
            return this->top_pattern_probs_.at(branch_idx).get_pattern_probability(
                    allele_count,
                    red_allele_count);
        }
        double get_top_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count) const {
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            return this->get_top_pattern_probability(0, allele_count, red_allele_count);
        }
        void set_bottom_pattern_probability(
                unsigned int branch_idx,
                unsigned int allele_count,
                unsigned int red_allele_count,
                double probability) {
            this->bottom_pattern_probs_.at(branch_idx).set_pattern_probability(
                    allele_count,
                    red_allele_count,
                    probability);
        }
        void set_bottom_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count,
                double probability) {
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            this->set_bottom_pattern_probability(0, allele_count, red_allele_count, probability);
        }
        void set_top_pattern_probability(
                unsigned int branch_idx,
                unsigned int allele_count,
                unsigned int red_allele_count,
                double probability) {
            this->top_pattern_probs_.at(branch_idx).set_pattern_probability(
                    allele_count,
                    red_allele_count,
                    probability);
        }
        void set_top_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count,
                double probability) {
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            this->set_top_pattern_probability(0, allele_count, red_allele_count, probability);
        }

        double get_population_size(unsigned int branch_idx) const {
            return this->population_sizes_.at(branch_idx)->get_value();
        }
        double get_population_size() const {
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            return this->get_population_size(0);
        }
        double get_population_size(std::shared_ptr<PopulationNetNode> parent_node) const {
            bool parent_found = false;
            for (unsigned int i = 0; i < this->parents_.size(); ++i) {
                if (this->parents_.at(i).lock() == parent_node) {
                    parent_found = true;
                    break;
                }
            }
            if (! parent_found) {
                throw EcoevolityError("PopulationNetNode::get_population_size(parent_node); not a parent!");
            }
            return this->get_population_size(i);
        }
        std::shared_ptr<PositiveRealParameter> get_population_size_parameter(
                unsigned int branch_idx) const {
            return this->population_sizes_.at(branch_idx);
        }
        std::shared_ptr<PositiveRealParameter> get_population_size_parameter() const {
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            return this->get_population_size_parameter(0);
        }
        std::vector< std::shared_ptr<PositiveRealParameter> > get_all_population_size_parameters() const {
            std::vector< std::shared_ptr<PositiveRealParameter> > parameters;
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            // parameters.reserve(this->get_node_count());
            this->get_all_population_size_parameters(
                    parameters, visited_nodes);
            return parameters;
        }

        void set_population_size_parameter(unsigned int branch_idx,
                std::shared_ptr<PositiveRealParameter> size) {
            this->population_size_.at(branch_idx) = size;
            this->make_all_dirty();
        }
        void set_population_size_parameter(std::shared_ptr<PositiveRealParameter> size) {
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            this->set_population_size_parameter(0, size);
        }
        void set_all_population_size_parameters(std::shared_ptr<PositiveRealParameter> size,
                std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) {
            this->visit(visited_nodes);
            for (auto pop_size : this->population_sizes_) {
                pop_size = size;
            }
            this->make_dirty();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->set_all_population_size_parameters(size, visited_nodes);
                }
            }
        }
        void set_all_population_size_parameters(std::shared_ptr<PositiveRealParameter> size) {
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            this->set_all_population_size_parameters(size, visited_nodes);
        }
        void set_all_population_size_parameters() {
            std::shared_ptr<PositiveRealParameter> size = this->population_size_.at(0);
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            this->set_all_population_size_parameters(size, visited_nodes);
        }

        void set_population_size(unsigned int branch_index, double size) {
            this->population_sizes_.at(branch_index)->set_value(size);
            this->make_all_dirty();
        }
        void set_population_size(double size) {
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            this->set_population_size(0, size);
        }
        void set_all_population_sizes(double size,
                std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) {
            this->visit(visited_nodes);
            for (auto pop_size : this->population_sizes_) {
                pop_size->set_value(size);
            }
            this->make_dirty();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->set_all_population_sizes(size, visited_nodes);
                }
            }
        }
        void set_all_population_sizes(double size) {
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            this->set_all_population_sizes(size, visited_nodes);
        }

        void update_population_size(unsigned int branch_index, double size) {
            this->population_sizes_.at(branch_index)->update_value(size);
            this->make_all_dirty();
        }
        void update_population_size(double size) {
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            this->update_population_size(0, size);
        }
        void update_all_population_sizes(double size,
                std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) {
            this->visit(visited_nodes);
            for (auto pop_size : this->population_sizes_) {
                pop_size->update_value(size);
            }
            this->make_dirty();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->update_all_population_sizes(size, visited_nodes);
                }
            }
        }
        void update_all_population_sizes(double size) {
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            this->update_all_population_sizes(size, visited_nodes);
        }

        void store_population_size() {
            for (auto pop_size : this->population_sizes_) {
                pop_size->store();
            }
        }
        void restore_population_size() {
            for (auto pop_size : this->population_sizes_) {
                pop_size->restore();
            }
            this->make_all_dirty();
        }
        void store_all_population_sizes(std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) {
            this->visit(visited_nodes);
            this->store_population_size();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->store_all_population_sizes(visited_nodes);
                }
            }
        }
        void restore_all_population_sizes(std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) {
            this->visit(visited_nodes);
            this->restore_population_size();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->restore_all_population_sizes(visited_nodes);
                }
            }
        }
        void store_all_population_sizes() {
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            this->store_all_population_sizes(visited_nodes);
        }
        void restore_all_population_sizes() {
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            this->restore_all_population_sizes(visited_nodes);
        }

        void store_population_size_pointer() {
            this->stored_population_sizes_ = this->population_sizes_;
        }
        void store_all_population_size_pointers(std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) {
            this->visit(visited_nodes);
            this->store_population_size_pointer();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->store_all_population_size_pointers(visited_nodes);
                }
            }
        }
        void store_all_population_size_pointers() {
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            this->store_all_population_size_pointers(visited_nodes);
        }

        void restore_population_size_pointer() {
            this->population_sizes_ = this->stored_population_sizes_;
            this->make_all_dirty();
        }
        void restore_all_population_size_pointers(std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) {
            this->visit(visited_nodes);
            this->restore_population_size_pointer();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->restore_all_population_size_pointers(visited_nodes);
                }
            }
        }
        void restore_all_population_size_pointers() {
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
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

        std::shared_ptr<ContinuousProbabilityDistribution> get_population_size_prior() {
            return this->population_sizes_.at(0)->get_prior();
        }
        void set_population_size_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            for (auto pop_size : this->population_sizes_) {
                pop_size->set_prior(prior);
            }
            this->make_all_dirty();
        }
        void set_all_population_size_priors(std::shared_ptr<ContinuousProbabilityDistribution> prior,
                std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) {
            this->visit(visited_nodes);
            this->set_population_size_prior(prior);
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->set_all_population_size_priors(prior, visited_nodes);
                }
            }
        }
        void set_all_population_size_priors(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            this->set_all_population_size_priors(prior, visited_nodes);
        }

        void fix_population_size() {
            for (auto pop_size : this->population_sizes_) {
                pop_size->fix();
            }
            this->make_dirty();
        }
        void fix_all_population_sizes(
                std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) {
            this->visit(visited_nodes);
            this->fix_population_size();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->fix_all_population_sizes(visited_nodes);
                }
            }
        }
        void fix_all_population_sizes() {
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            this->fix_all_population_sizes(visited_nodes);
        }

        void estimate_population_size() {
            for (auto pop_size : this->population_sizes_) {
                pop_size->estimate();
            }
            this->make_dirty();
        }
        void estimate_all_population_sizes(
                std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) {
            this->visit(visited_nodes);
            this->estimate_population_size();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->estimate_all_population_sizes(visited_nodes);
                }
            }
        }
        void estimate_all_population_sizes() {
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            this->estimate_all_population_sizes(visited_nodes);
        }

        double calculate_ln_relative_population_size_prior_density() const {
            double d = 0.0;
            std::vector< std::shared_ptr<PositiveRealParameter> > parameters;
            parameters.reserve(this->get_node_count());
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            this->add_ln_relative_population_size_prior_density(d, parameters, visited_nodes);
            return d;
        }

        unsigned int scale_all_population_sizes(double scale) {
            std::vector< std::shared_ptr<PositiveRealParameter> > parameters;
            parameters.reserve(this->get_node_count());
            unsigned int number_of_free_parameters_scaled = 0;
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            this->scale_all_population_sizes(
                    scale,
                    parameters,
                    visited_nodes,
                    number_of_free_parameters_scaled);
            return number_of_free_parameters_scaled;
        }

        bool population_size_is_fixed() const {
            bool f = this->population_sizes_.at(0)->is_fixed();
            for (unsigned int i = 1; i < this->population_sizes_.size(); ++i) {
                ECOEVOLITY_ASSERT(this->population_sizes_.at(i)->is_fixed() == f);
            }
            return f;
        }

        bool all_population_sizes_are_fixed(
                std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) const {
            this->visit(visited_nodes);
            if (! this->population_size_is_fixed()) {
                return false;
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    if (! child_iter->all_population_sizes_are_fixed(visited_nodes)) {
                        return false;
                    }
                }
            }
            return true;
        }
        bool all_population_sizes_are_fixed() const {
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            return this->all_population_sizes_are_fixed(visited_nodes);
        }

        double get_population_size_relative_prior_ln_pdf(unsigned int branch_index) const {
            return this->population_size_.at(branch_index)->relative_prior_ln_pdf();
        }
        double get_population_size_relative_prior_ln_pdf() const {
            ECOEVOLITY_ASSERT(this->get_number_of_parents() < 2);
            return this->get_population_size_relative_prior_ln_pdf(0);
        }

        void get_node_indices(std::vector<unsigned int> & internal_indices,
                std::vector<unsigned int> & leaf_indices,
                std::set< std::shared_ptr<const PopulationNetNode> > & visited_nodes) {
            this->visit(visited_nodes);
            if (this->is_leaf()) {
                leaf_indices.push_back(this->get_index());
            } else {
                internal_indices.push_back(this->get_index());
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->get_node_indices(internal_indices, leaf_indices, visited_nodes);
                }
            }
        }
        void get_node_indices(std::vector<unsigned int> & internal_indices,
                std::vector<unsigned int> & leaf_indices) {
            std::set< std::shared_ptr<const PopulationNetNode> > visited_nodes;
            this->get_node_indices(internal_indices, leaf_indices, visited_nodes);
        }

        // Overriding this method from BaseNetNode to return pop size
        std::string get_additional_comment_data_string(
                unsigned int precision = 12) const {
            std::ostringstream s;
            s.precision(precision);
            for (unsigned int i = 0; i < this->population_sizes_.size(); ++i) {
                s << "pop_size" << i << "="
                  << this->get_population_size(i);
            }
            return s.str();
        }

        // Overriding this method from BaseNetNode to make sure newly inserted
        // nodes are fully initialized
        void finish_initializing_inserted_internal_node(
                RandomNumberGenerator & rng) {
            ECOEVOLITY_ASSERT(this->has_parent() && this->has_children());
            if (this->get_parent()->get_population_size_parameter(0) == this->children_.at(0)->get_population_size_parameter(0)) {
                for (auto pop_size : this->population_sizes_) {
                    pop_size = this->children_.at(0)->get_population_size_parameter(0);
                }
                return;
            }
            this->estimate_population_size();
            this->set_population_size_prior(this->children_.at(0)->get_population_size_prior());
            if (this->children_.at(0)->population_size_is_fixed()) {
                double child_pop_size_sum = 0.0;
                for (auto child : this->children_) {
                    child_pop_size_sum += child->get_population_size(this->shared_from_this());
                }
                double new_pop_size = child_pop_size_sum / this->children_.size();
                for (auto pop_size : this->population_sizes_) {
                    pop_size->set_value(new_pop_size);
                    pop_size->fix();
                }
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
            for (auto pop_size : this->population_sizes_) {
                pop_size->set_value_from_prior(rng);
            }
        }

        virtual double get_ln_prob_of_drawing_state() {
            ECOEVOLITY_ASSERT(this->has_parent() && this->has_children());
            if (this->population_size_is_fixed()) {
                return 0.0;
            }
            if (this->children_.at(0)->population_sizes_.at(0) == this->population_sizes_.at(0)) {
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
            double ln_prob = this->population_sizes_.at(0)->prior_ln_pdf();
            for (unsigned int i = 1; i < this->population_sizes_.size(); ++i) {
                ln_prob += this->population_sizes_.at(i)->prior_ln_pdf()
            }
            return ln_prob;
        }
};

#endif
