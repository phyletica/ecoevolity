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

#ifndef ECOEVOLITY_BASENODE_HPP
#define ECOEVOLITY_BASENODE_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>

#include "matrix.hpp"
#include "parameter.hpp"
#include "probability.hpp"
#include "debug.hpp"
#include "assert.hpp"
#include "error.hpp"

/**
 * Base class for a node of a phylogenetic tree.
 *
 * @note    Many of the class' methods modified from:
 *              BasicTNode class of Bio++ Library
 *              <http://biopp.univ-montp2.fr/wiki/index.php/Main_Page>
 *              License:    CeCILL <http://www.cecill.info>
 *              Author:     Sylvain Gaillard
 *              Copyright:  CNRS, (January 12, 2011)
 */
template<class DerivedNodeT>
class BaseNode : public std::enable_shared_from_this<DerivedNodeT> {
    protected:
        std::vector< std::shared_ptr<DerivedNodeT> > children_;
        std::shared_ptr<DerivedNodeT> parent_ = nullptr;
        std::string label_ = "";
        std::shared_ptr<PositiveRealParameter> height_ = std::make_shared<PositiveRealParameter>(0.0);
        bool is_dirty_ = true;

        void add_ln_relative_node_height_prior_density(
                double& density,
                std::vector< std::shared_ptr<PositiveRealParameter> >& parameters) const {
            bool parameter_found = false;
            for (auto parameter_iter : parameters) {
                if (parameter_iter == this->height_) {
                    parameter_found = true;
                    break;
                }
            }
            if (! parameter_found) {
                density += this->height_->relative_prior_ln_pdf();
                parameters.push_back(this->height_);
            }
            for (auto child_iter: this->children_) {
                child_iter->add_ln_relative_node_height_prior_density(density, parameters);
            }
        }

    public:
        // Constructors
        BaseNode() { }
        BaseNode(std::string label) {
            this->label_ = label;
        }
        BaseNode(double height) {
            this->height_->set_value(height);
        }
        BaseNode(std::string label, double height) {
            this->label_ = label;
            this->height_->set_value(height);
        }

        // Destructor
        // virtual ~BaseNode() { }

        /* DerivedNodeT* clone() const { */
        /*     return new DerivedNodeT(static_cast<DerivedNodeT const &>(* this)); */
        /* } */
        
        //Methods
        unsigned int degree() const {
            unsigned int d = children_.size();
            if (this->has_parent()) {
                d += 1;
            }
            return d;
        }

        bool has_parent() const { return parent_ ? true : false; }
        bool is_root() const { return parent_ ? false : true; }

        unsigned int get_number_of_parents() const { return parent_ ? 1 : 0; }

        const std::shared_ptr<DerivedNodeT>& get_parent() const {
            return this->parent_;
        }
        std::shared_ptr<DerivedNodeT> get_parent() {
            return this->parent_;
        }

        bool is_parent(const std::shared_ptr<DerivedNodeT>& node) const {
            if (this->parent_ == node) {
                return true;
            }
            return false;
        }

        void add_parent(std::shared_ptr<DerivedNodeT> node) {
            if (!node) {
                throw EcoevolityNullPointerError("BaseNode::add_parent(), empty node given");
            }
            if (this->parent_) {
                throw EcoevolityError("BaseNode::add_parent(), this node already has a parent");
            }
            this->parent_ = node;
            if (! node->is_child(this->shared_from_this())) {
                node->add_child(this->shared_from_this());
            }
        }

        std::shared_ptr<DerivedNodeT> remove_parent() {
            if (this->has_parent()) {
                std::shared_ptr<DerivedNodeT> p = this->parent_;
                this->parent_ = nullptr;
                p->remove_child(this->shared_from_this());
                return p;
            }
            return nullptr;
        }

        bool has_children() const { return !this->children_.empty(); }
        bool is_leaf() const { return this->children_.empty(); }

        unsigned int get_number_of_children() const { return this->children_.size(); }
        const std::shared_ptr<DerivedNodeT>& get_child(unsigned int index) const {
            if (index >= this->children_.size()) {
                throw std::out_of_range("BaseNode::get_child() index out of range");
            }
            return this->children_.at(index);
        }

        std::shared_ptr<DerivedNodeT> get_child(unsigned int index) {
            if (index >= this->children_.size()) {
                throw std::out_of_range("BaseNode::get_child() index out of range");
            }
            return this->children_.at(index);
        }

        bool is_child(const std::shared_ptr<DerivedNodeT>& node) const {
            for (auto child_iter: this->children_) {
                if (child_iter == node) {
                    return true;
                }
            }
            return false;
        }

        void add_child(std::shared_ptr<DerivedNodeT> node) {
            if (!node) {
                throw EcoevolityNullPointerError("BaseNode::add_child(), empty node given");
            }
            if (! this->is_child(node)) {
                this->children_.push_back(node);
            }
            if (! node->is_parent(this->shared_from_this())) {
                node->add_parent(this->shared_from_this());
            }
        }

        void remove_child(std::shared_ptr<DerivedNodeT> node) {
            if (!node) {
                throw EcoevolityNullPointerError("BaseNode::remove_child(), empty node given");
            }
            for (unsigned int i = 0; i < this->children_.size(); ++i) {
                if (this->children_.at(i) == node) {
                    this->children_.erase(this->children_.begin() + i);
                    node->remove_parent();
                }
            }
        }

        std::shared_ptr<DerivedNodeT> remove_child(unsigned int index) {
            if (index >= this->children_.size()) {
                throw std::out_of_range("BaseNode::remove_child() index out of range");
            }
            std::shared_ptr<DerivedNodeT> c = this->children_.at(index);
            this->children_.erase(this->children_.begin() + index);
            c->remove_parent();
            return c;
        }

        double get_height() const {
            return this->height_->get_value();
        }

        void set_height(double height) {
            this->height_->set_value(height);
            this->make_all_dirty();
        }

        void set_height_parameter(std::shared_ptr<PositiveRealParameter> height_parameter) {
            this->height_ = height_parameter;
            this->make_all_dirty();
        }

        std::shared_ptr<PositiveRealParameter> get_height_parameter() const {
            return this->height_;
        }

        void update_height(double height) {
            this->height_->update_value(height);
            this->make_all_dirty();
        }

        void store_height() {
            this->height_->store();
        }
        void store_all_heights() {
            this->height_->store();
            for (auto child_iter: this->children_) {
                child_iter->store_all_heights();
            }
        }

        void restore_height() {
            this->height_->restore();
            this->make_all_dirty();
        }
        void restore_all_heights() {
            this->height_->restore();
            this->make_dirty();
            for (auto child_iter: this->children_) {
                child_iter->restore_all_heights();
            }
        }

        double get_length() const {
            if (this->has_parent()) {
                return this->get_parent()->get_height() - this->get_height();
            }
            return 0.0;
        }

        const std::string& get_label() const {
            return this->label_;
        }
        void set_label(std::string label) {
            this->label_ = label;
        }

        bool is_dirty() const {
            return this->is_dirty_;
        }

        bool clade_has_dirt() const {
            if (this->is_dirty()) {
                return true;
            }
            for (auto child_iter: this->children_) {
                if (child_iter->clade_has_dirt()) {
                    return true;
                }
            }
            return false;
        }

        void make_dirty() {
            this->is_dirty_ = true;
        }
        void make_all_dirty() {
            this->is_dirty_ = true;
            for (auto child_iter: this->children_) {
                child_iter->make_all_dirty();
            }
        }
        void make_clean() {
            this->is_dirty_ = false;
        }
        void make_all_clean() {
            this->is_dirty_ = false;
            for (auto child_iter: this->children_) {
                child_iter->make_all_clean();
            }
        }

        unsigned int get_node_count() const {
            unsigned int n = 1;
            for (auto child_iter: this->children_) {
                n += child_iter->get_node_count();
            }
            return n;
        }
        unsigned int get_leaf_node_count() const {
            if (this->is_leaf()) {
                return 1;
            }
            unsigned int n = 0;
            for (auto child_iter: this->children_) {
                n += child_iter->get_leaf_node_count();
            }
            return n;
        }
        unsigned int get_internal_node_count() const {
            if (this->is_leaf()) {
                return 0;
            }
            unsigned int n = 1;
            for (auto child_iter: this->children_) {
                n += child_iter->get_internal_node_count();
            }
            return n;
        }

        void set_node_height_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->height_->set_prior(prior);
            this->make_all_dirty();
        }
        void set_all_node_height_priors(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->height_->set_prior(prior);
            this->make_dirty();
            for (auto child_iter: this->children_) {
                child_iter->set_all_node_height_priors(prior);
            }
        }

        void fix_node_height() {
            this->height_->fix();
        }
        void fix_all_node_heights() {
            this->height_->fix();
            for (auto child_iter: this->children_) {
                child_iter->fix_all_node_heights();
            }
        }
        void estimate_node_height() {
            this->height_->estimate();
        }
        void estimate_all_node_heights() {
            this->height_->estimate();
            for (auto child_iter: this->children_) {
                child_iter->estimate_all_node_heights();
            }
        }

        bool node_height_is_fixed() const {
            return this->height_->is_fixed();
        }

        bool all_node_heights_are_fixed() const {
            if (! this->node_height_is_fixed()) {
                return false;
            }
            for (auto child_iter: this->children_) {
                if (! child_iter->all_node_heights_are_fixed()) {
                    return false;
                }
            }
            return true;
        }

        double calculate_ln_relative_node_height_prior_density() const {
            double d = 0.0;
            std::vector< std::shared_ptr<PositiveRealParameter> > parameters;
            parameters.reserve(this->get_node_count());
            this->add_ln_relative_node_height_prior_density(d, parameters);
            return d;
        }

};

#endif
